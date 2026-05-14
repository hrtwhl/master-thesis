# Install required packages if you don't have them already:
# install.packages(c("quantmod", "zoo"))

library(quantmod)
library(zoo)

# Define the tickers and their desired column names
ticker_mapping <- c(
  stocks   = "SPY",
  bonds = "AGG",
  gold  = "GLD",
  oil   = "USO",
  usd   = "UUP"
)

# Create an empty environment to store the downloaded data
data_env <- new.env()

# Download daily data from Yahoo Finance
# Setting 'from' to a very early date ensures we get the earliest available data for each ticker
getSymbols(ticker_mapping, env = data_env, src = "yahoo", from = "1900-01-01", auto.assign = TRUE)

# Extract the Adjusted Close prices for each ticker into a list
# (Adjusted close is best practice for daily data as it accounts for dividends and splits)
price_list <- lapply(ticker_mapping, function(sym) {
  Ad(get(sym, envir = data_env))
})

# Merge all the individual xts time-series objects into one combined xts object
# This automatically aligns the dates and introduces NAs where days don't match
merged_xts <- do.call(merge, price_list)

# Rename the columns to your specified names (spx, bonds, gold, oil, usd)
colnames(merged_xts) <- names(ticker_mapping)

# Forward fill the missing values
# na.rm = FALSE ensures leading NAs (before a specific ETF existed) are kept as NA
merged_xts_filled <- na.locf(merged_xts, na.rm = FALSE)

# Convert the final xts object into a standard R dataframe
# Extract the date index into its own 'date' column
final_df <- data.frame(
  date = index(merged_xts_filled), 
  coredata(merged_xts_filled)
)

# Reset row names for a clean dataframe
rownames(final_df) <- NULL

# View the first and last few rows to verify
head(final_df)
tail(final_df)

final_df <- na.omit(final_df)

write.csv(final_df, "Data/asset_data_yf.csv", row.names = FALSE)







#------


# Load required library
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# 1. Load the data (Make sure the file paths match your working directory)
data_bbg <- read.csv("Data/asset_data.csv", stringsAsFactors = FALSE)
data_yf  <- read.csv("Data/asset_data_yf.csv", stringsAsFactors = FALSE)

# Convert dates to proper Date objects
data_bbg$date <- as.Date(data_bbg$date)
data_yf$date  <- as.Date(data_yf$date)

# 2. Find the overlapping date range
start_date <- max(min(data_bbg$date, na.rm=TRUE), min(data_yf$date, na.rm=TRUE))
end_date   <- min(max(data_bbg$date, na.rm=TRUE), max(data_yf$date, na.rm=TRUE))

cat("========================================================\n")
cat("DATA COMPARISON: OVERLAPPING PERIOD\n")
cat("From:", as.character(start_date), "to", as.character(end_date), "\n")
cat("========================================================\n\n")

# Filter both datasets to the overlapping period and arrange chronologically
bbg_overlap <- data_bbg %>% filter(date >= start_date & date <= end_date) %>% arrange(date)
yf_overlap  <- data_yf  %>% filter(date >= start_date & date <= end_date) %>% arrange(date)

# Merge datasets for a direct day-to-day comparison
merged_data <- inner_join(bbg_overlap, yf_overlap, by = "date", suffix = c("_bbg", "_yf"))

assets <- c("stocks", "bonds", "oil", "gold", "usd")

# 3. Summary Statistics (Prices)
cat("--- 1. RAW PRICE COMPARISON ---\n")
cat("(Note: Large differences here are fine if one is an Index and the other is an ETF)\n\n")

for (asset in assets) {
  col_bbg <- paste0(asset, "_bbg")
  col_yf  <- paste0(asset, "_yf")
  
  cor_price <- cor(merged_data[[col_bbg]], merged_data[[col_yf]], use = "complete.obs")
  
  cat(toupper(asset), "\n")
  cat(sprintf("  Price Correlation: %.4f\n", cor_price))
  cat(sprintf("  BBG Source - Mean: %8.2f | Min: %8.2f | Max: %8.2f\n",
              mean(merged_data[[col_bbg]], na.rm=TRUE), 
              min(merged_data[[col_bbg]], na.rm=TRUE), 
              max(merged_data[[col_bbg]], na.rm=TRUE)))
  cat(sprintf("  YF  Source - Mean: %8.2f | Min: %8.2f | Max: %8.2f\n\n",
              mean(merged_data[[col_yf]], na.rm=TRUE), 
              min(merged_data[[col_yf]], na.rm=TRUE), 
              max(merged_data[[col_yf]], na.rm=TRUE)))
}

# 4. Compare Log Returns (The most critical metric for the HMM model)
cat("--- 2. DAILY LOG RETURN COMPARISON ---\n")
cat("(We want correlations here to be as close to 1.000 as possible)\n\n")

for (asset in assets) {
  col_bbg <- paste0(asset, "_bbg")
  col_yf  <- paste0(asset, "_yf")
  
  # Calculate daily log returns
  ret_bbg <- c(NA, diff(log(merged_data[[col_bbg]])))
  ret_yf  <- c(NA, diff(log(merged_data[[col_yf]])))
  
  # Return metrics
  cor_ret <- cor(ret_bbg, ret_yf, use = "complete.obs")
  mae_ret <- mean(abs(ret_bbg - ret_yf), na.rm=TRUE)
  
  vol_bbg <- sd(ret_bbg, na.rm=TRUE) * sqrt(252) * 100
  vol_yf  <- sd(ret_yf, na.rm=TRUE) * sqrt(252) * 100
  
  cat(toupper(asset), "RETURNS\n")
  cat(sprintf("  Return Correlation: %.4f\n", cor_ret))
  cat(sprintf("  Mean Abs Difference: %.6f (Daily)\n", mae_ret))
  cat(sprintf("  Annualized Volatility - BBG: %5.2f%% | YF: %5.2f%%\n\n", vol_bbg, vol_yf))
}

# 5. Missing Data Check
cat("--- 3. MISSING VALUES IN OVERLAPPING PERIOD ---\n")
missing_counts <- merged_data %>% summarise(across(everything(), ~sum(is.na(.))))
missing_df <- data.frame(
  Asset = assets,
  Missing_BBG = as.numeric(missing_counts[paste0(assets, "_bbg")]),
  Missing_YF  = as.numeric(missing_counts[paste0(assets, "_yf")])
)
print(missing_df, row.names = FALSE)
cat("\nDone.\n")

