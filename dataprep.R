df_bbg <- read.csv("market_data.csv")
df_bbg <- df_bbg[, c("date", "sp500", "copper", "vix")]
colnames(df_bbg)

# Load required libraries
library(fredr)
library(dplyr)
library(tidyr)
library(purrr)

# Your API key script (already run by you)
source("Data/APIKey.R")

# Your list of FRED IDs
fred_ids <- c(
  us10y = "DGS10",
  us3m  = "DGS3MO",
  yield_curve = "T10Y3M",
  oil = "WTISPLC"
)

# 1. Download data and combine into a "long" format dataframe
df_long <- imap_dfr(fred_ids, function(id, name) {
  # Fetch data for each ID starting from 1900
  fredr(
    series_id = id,
    observation_start = as.Date("1900-01-01")
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
colnames(df_fred)


# Load required libraries
library(dplyr)
library(tidyr)
library(zoo)

# 1. Clean and format the Bloomberg data (df_bbg)
# Convert 'date' to actual Date format, and the rest to numeric
df_bbg_clean <- df_bbg %>%
  mutate(
    date = as.Date(date),
    across(c(sp500, copper, vix), as.numeric)
  )

# 2. Merge and chronologically sort the dataframes
# Now both have 'date' as a Date object, so they will join perfectly
df_merged <- full_join(df_bbg_clean, df_fred, by = "date") %>%
  arrange(date)

# 3. Perform calculations
df_calc <- df_merged %>%
  mutate(
    # Calculate daily log returns for S&P 500 and daily delta for US10Y
    sp500_ret = log(sp500 / lag(sp500)),
    us10y_delta = us10y - lag(us10y)
  ) %>%
  mutate(
    # Calculate 22-day rolling realized volatility 
    # Annualized (* sqrt(252)) and converted to percentage (* 100) to match VIX scale
    sp500_rv = rollapplyr(
      sp500_ret, 
      width = 22, 
      FUN = function(x) sd(x, na.rm = TRUE) * sqrt(252) * 100, 
      fill = NA
    ),
    # Extend VIX using coalesce: keeps actual VIX where available, fills with RV otherwise
    vix = coalesce(vix, sp500_rv)
  ) %>%
  mutate(
    # Calculate 3-year (approx. 756 trading days) rolling correlation
    stock_bond_corr = rollapplyr(
      data = cbind(sp500_ret, us10y_delta),
      width = 756,
      FUN = function(x) {
        # Only calculate correlation if we have a reasonable amount of overlapping valid data
        if(sum(complete.cases(x)) > 30) { 
          cor(x[,1], x[,2], use = "pairwise.complete.obs")
        } else {
          NA
        }
      },
      by.column = FALSE,
      fill = NA
    )
  ) %>%
  # 4. Select the final variables 
  select(date, sp500, yield_curve, oil, copper, us3m, vix, stock_bond_corr)

# 5. Find the common starting date
# Find the first non-NA date for each variable, then take the maximum of those dates
start_dates <- sapply(df_calc %>% select(-date), function(x) min(df_calc$date[!is.na(x)], na.rm = TRUE))
common_start_date <- max(as.Date(start_dates, origin = "1970-01-01"))

# 6. Filter to the start date and forward-fill
df_final <- df_calc %>%
  filter(date >= common_start_date) %>%
  # Forward fill any subsequent missing values (e.g., daily gaps in monthly oil data)
  fill(everything(), .direction = "down")

# View the result
head(df_final)


library(ggplot2)

df_final |> ggplot(aes(x=date)) +
  geom_line(aes(y=stock_bond_corr), color = "black", linewidth=0.5) +
  theme_bw() +
  labs(title = "Stock Bond Correlation", x = NULL, y=NULL) +
  theme(panel.grid.major = element_blank(), 
panel.grid.minor = element_blank())


out_path <- "etf_prices.csv"
write.csv(df_final, "data_final.csv", row.names = FALSE)
