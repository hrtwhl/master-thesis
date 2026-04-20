# Download Bloomberg historical time series and export to CSV
# Requires: active Bloomberg terminal session + Rblpapi

library(Rblpapi)

blpConnect()

# name -> ticker, requested start date
# Alternatives noted in comments; swap if a ticker doesn't pull as expected.
series <- list(
  sp500  = list(ticker = "SPX Index",      start = "1927-01-01"),  # S&P 500 (BBG history usually starts ~1928)
  copper = list(ticker = "HG1 Comdty",     start = "1959-01-01"),  # COMEX copper generic front month; alt: "LP1 Comdty" (LME 3M)
  vix    = list(ticker = "VIX Index",      start = "1990-01-01"),  # CBOE VIX
  us10y  = list(ticker = "USGG10YR Index", start = "1962-01-01")   # BBG generic 10Y yield; alt: "H15T10Y Index" (Fed H.15 CMT)
)

opts <- c(
  periodicitySelection    = "DAILY",
  nonTradingDayFillOption = "ACTIVE_DAYS_ONLY"
)

fetch_one <- function(name, spec) {
  message("Fetching ", spec$ticker)
  df <- bdh(
    securities = spec$ticker,
    fields     = "PX_LAST",
    start.date = as.Date(spec$start),
    end.date   = Sys.Date(),
    options    = opts
  )
  names(df)[names(df) == "PX_LAST"] <- name
  df
}

data_list <- Map(fetch_one, names(series), series)

# Full outer join on date, sorted ascending
merged <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), data_list)
merged <- merged[order(merged$date), ]

out_path <- "bloomberg_series.csv"
write.csv(merged, out_path, row.names = FALSE)

message(sprintf("Saved %d rows to %s. Range: %s to %s",
                nrow(merged), out_path, min(merged$date), max(merged$date)))