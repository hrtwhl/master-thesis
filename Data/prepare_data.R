# Load necessary libraries
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(zoo)
library(ggplot2)

# 1. Define the file names
files <- c(
  "Data/raw/bonds.csv", 
  "Data/raw/stocks.csv", 
  "Data/raw/oil.csv", 
  "Data/raw/gold.csv", 
  "Data/raw/usd.csv",
  "Data/raw/copper.csv", 
  "Data/raw/vix.csv", 
  "Data/raw/us3mo.csv", 
  "Data/raw/stock_bond_corr.csv", 
  "Data/raw/yield_curve.csv"
)

# 2. Helper function to read files and ensure date formatting
read_asset_csv <- function(file_name) {
  read_csv(file_name, show_col_types = FALSE) %>%
    mutate(date = as.Date(date))
}

# 3. Read all files into a list of dataframes
list_of_dfs <- map(files, read_asset_csv)

# 4. Merge all dataframes together 
# We use full_join to ensure no rows are dropped if a date exists in one file but not another
df_master <- list_of_dfs %>%
  reduce(full_join, by = "date") %>%
  arrange(date) %>%
  # 5. Fill forward missing values (LOCF)
  # .direction = "down" carries the last valid observation forward
  fill(everything(), .direction = "down") %>%
  # 6. Filter for your specific date range
  filter(date >= as.Date("1990-01-02") & date <= as.Date("2025-12-31"))

# 7. Create df_assets
df_assets <- df_master %>%
  select(date, stocks, bonds, oil, gold, usd)

# 8. Create df_macro
df_macro <- df_master %>%
  select(date, stocks, oil, copper, yield_curve, stock_bond_corr, vix, us3mo)

# Check the results
head(df_assets)
tail(df_assets)

head(df_macro)
tail(df_macro)

# 9. Export dataframes as csv files

write.csv(df_assets, "Data/asset_data.csv", row.names = FALSE)
write.csv(df_macro, "Data/macro_data.csv", row.names = FALSE)


# 10. Plot data for a visual sanity check
plot_time_series_grid <- function(df, title) {
  df %>%
    # Convert to long format for ggplot
    pivot_longer(cols = -date, names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = date, y = value)) +
    geom_line(color = "darkblue", linewidth =0.2) +
    # Create the grid: scales = "free_y" lets each chart have its own Y-axis range
    facet_wrap(~ variable, scales = "free_y", ncol = 3) + 
    labs(
      title = title,
      x = "Year",
      y = "Value"
    ) +
    theme_minimal() 
}

# plot for df_assets
plot_assets <- plot_time_series_grid(df_assets, "Financial Assets Time Series (1990-2025)")
print(plot_assets)

# plot for df_macro
plot_macro <- plot_time_series_grid(df_macro, "Macroeconomic Indicators Time Series (1990-2025)")
print(plot_macro)
