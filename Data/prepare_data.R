# Load necessary libraries
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(zoo)

# 1. Define the file names
files <- c(
  "Data/bonds.csv", "Data/stocks.csv", "Data/oil.csv", "Data/gold.csv", "Data/usd.csv",
  "Data/copper.csv", "Data/vix.csv", "Data/us3mo.csv", "Data/stock_bond_corr.csv", "Data/yield_curve.csv"
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

View(df_assets)
View(df_macro)

library(ggplot2)
library(tidyr)

# Helper function to create the grid plot
plot_time_series_grid <- function(df, title) {
  df %>%
    # Convert to long format for ggplot
    pivot_longer(cols = -date, names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = date, y = value)) +
    geom_line(color = "steelblue") +
    # Create the grid: scales = "free_y" lets each chart have its own Y-axis range
    facet_wrap(~ variable, scales = "free_y", ncol = 3) + 
    labs(
      title = title,
      x = "Year",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold")
    )
}

# 1. Plot for df_assets
plot_assets <- plot_time_series_grid(df_assets, "Financial Assets Time Series (1990-2025)")
print(plot_assets)

# 2. Plot for df_macro
plot_macro <- plot_time_series_grid(df_macro, "Macroeconomic Indicators Time Series (1990-2025)")
print(plot_macro)
