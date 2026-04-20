# Download SPDR sector ETF historical prices from Bloomberg and export to CSV.
# Requires: active Bloomberg terminal, Rblpapi.

library(Rblpapi)

blpConnect()

etfs    <- c("SPY", "XLB", "XLE", "XLF", "XLI", "XLK", "XLP", "XLU", "XLV", "XLY")
tickers <- paste(etfs, "US Equity")

opts <- c(
  periodicitySelection    = "DAILY",
  nonTradingDayFillOption = "ACTIVE_DAYS_ONLY"
)

# bdh with multiple securities returns a named list, one frame per ticker.
# Start well before any ETF inception (SPY: 1993, XL* sector SPDRs: Dec 1998)
# -- Bloomberg returns from actual inception onwards.
raw <- bdh(
  securities = tickers,
  fields     = "PX_LAST",
  start.date = as.Date("1990-01-01"),
  end.date   = Sys.Date(),
  options    = opts
)

# Rename PX_LAST column in each frame to the short ETF ticker
df_list <- Map(function(d, short) {
  names(d)[names(d) == "PX_LAST"] <- short
  d
}, raw, etfs)

df <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), df_list)
df <- df[order(df$date), ]

out_path <- "etf_prices.csv"
write.csv(df, out_path, row.names = FALSE)

message(sprintf("Saved %d rows to %s. Range: %s to %s",
                nrow(df), out_path, min(df$date), max(df$date)))