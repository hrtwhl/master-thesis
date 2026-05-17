# Load necessary libraries
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(zoo)
library(ggplot2)

# 1. Define the file names (Added oil.csv to the list)
files <- c(
  "Data/raw/bonds.csv", 
  "Data/raw/stocks.csv", 
  "Data/raw/brent.csv", 
  "Data/raw/oil.csv", # WTI Crude
  "Data/raw/gold.csv", 
  "Data/raw/usd.csv",
  "Data/raw/copper.csv", 
  "Data/raw/vix.csv", 
  "Data/raw/us3mo.csv", 
  "Data/raw/stock_bond_corr.csv", 
  "Data/raw/yield_curve.csv"
)

# 2. Helper function to read files, ensure date formatting, and prevent naming conflicts
read_asset_csv <- function(file_name) {
  df <- read_csv(file_name, show_col_types = FALSE) %>%
    mutate(date = as.Date(date))
  
  # Check if the file is brent.csv or oil.csv and rename the 'oil' column accordingly 
  # to prevent duplicate column names (oil.x, oil.y) during the full_join
  if (grepl("brent.csv", file_name) && "oil" %in% names(df)) {
    df <- df %>% rename(brent = oil)
  }
  if (grepl("oil.csv", file_name) && "oil" %in% names(df)) {
    df <- df %>% rename(wti = oil)
  }
  
  return(df)
}

# 3. Read all files into a list of dataframes
list_of_dfs <- map(files, read_asset_csv)

# 4. Merge all dataframes together 
df_master <- list_of_dfs %>%
  reduce(full_join, by = "date") %>%
  arrange(date) %>%
  # 5. Fill forward missing values (LOCF)
  fill(everything(), .direction = "down") %>%
  # 6. Filter for your specific date range AND keep only weekdays
  filter(date >= as.Date("1990-01-02") & date <= as.Date("2025-12-31")) %>%
  # format(date, "%u") returns 1 for Monday through 7 for Sunday. 1:5 filters out weekends.
  filter(format(date, "%u") %in% 1:5)

# 7. Create df_assets (Using Brent)
df_assets <- df_master %>%
  # Select brent and rename it to 'oil' on the fly for the final dataframe
  select(date, stocks, bonds, oil = brent, gold, usd)

# 8. Create df_macro (Using WTI)
df_macro <- df_master %>%
  # Select wti and rename it to 'oil' on the fly for the final dataframe
  select(date, stocks, oil = wti, copper, yield_curve, stock_bond_corr, vix, us3mo)

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

#----------------------------------------------------------------
# Code snippet to get data from FRED, clean and export as a csv
#----------------------------------------------------------------

library(fredr)
library(tidyverse)
fredr_set_key("263fff36e95136e168c2b9128597195d")

# Your list of FRED IDs
fred_ids <- c(
  oil = "DCOILBRENTEU"
)

# 1. Download data and combine into a "long" format dataframe
df_long <- imap_dfr(fred_ids, function(id, name) {
  # Fetch data for each ID starting from 1900
  fredr(
    series_id = id,
    observation_start = as.Date("1960-01-01")
  ) %>%
    # Select only the relevant columns and assign your custom variable name
    select(date, value) %>%
    mutate(variable = name)
})

# 2. Pivot the data into a "wide" format to get your final df_fred
df_fred <- df_long %>%
  pivot_wider(
    names_from = variable, 
    values_from = value
  ) %>%
  arrange(date) # Ensure dates are chronologically ordered

# View the first few rows
View(df_fred)

write.csv(df_fred, "Data/raw/brent.csv", row.names = FALSE)
