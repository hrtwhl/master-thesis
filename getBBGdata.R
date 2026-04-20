# Download Bloomberg + FRED historical series and build derived regime variables.
# Requires: active Bloomberg terminal, Rblpapi, fredr (FRED API key), zoo.

library(Rblpapi)
library(fredr)
library(zoo)

# fredr_set_key("YOUR_FRED_KEY")   # or set FRED_API_KEY in ~/.Renviron

blpConnect()

# --- 1. Bloomberg series ----------------------------------------------------

bbg_series <- list(
  sp500  = list(ticker = "SPX Index",          start = "1927-01-01"),
  copper = list(ticker = "LMCADY LME Comdty",  start = "1986-01-01"),
  vix    = list(ticker = "VIX Index",          start = "1990-01-01"),
  us10y  = list(ticker = "USGG10YR Index",     start = "1962-01-01")
)

bbg_opts <- c(
  periodicitySelection    = "DAILY",
  nonTradingDayFillOption = "ACTIVE_DAYS_ONLY"
)

fetch_bbg <- function(name, spec) {
  message("Bloomberg: ", spec$ticker)
  df <- bdh(
    securities = spec$ticker,
    fields     = "PX_LAST",
    start.date = as.Date(spec$start),
    end.date   = Sys.Date(),
    options    = bbg_opts
  )
  names(df)[names(df) == "PX_LAST"] <- name
  df
}

# --- 2. FRED series ---------------------------------------------------------

fred_series <- list(
  oil  = list(ticker = "DCOILWTICO", start = "1986-01-01"),  # WTI Cushing spot
  us3m = list(ticker = "DTB3",       start = "1954-01-01")   # 3M T-bill secondary market
)

fetch_fred <- function(name, spec) {
  message("FRED: ", spec$ticker)
  df <- fredr(series_id         = spec$ticker,
              observation_start = as.Date(spec$start),
              observation_end   = Sys.Date())
  out <- data.frame(date = df$date, value = df$value)
  names(out)[2] <- name
  out
}

# --- 3. Fetch and merge -----------------------------------------------------

bbg_data  <- Map(fetch_bbg,  names(bbg_series),  bbg_series)
fred_data <- Map(fetch_fred, names(fred_series), fred_series)

df <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE),
             c(bbg_data, fred_data))
df <- df[order(df$date), ]

# --- 4. Derived regime variables -------------------------------------------

# 4a. Yield curve slope (10Y - 3M)
df$yield_curve <- df$us10y - df$us3m

# 4b. Extended VIX: trailing 21d annualised realised vol of SPX pre-1990,
#     actual VIX from 1990 onwards. Trailing window -> no look-ahead.
spx_logret <- c(NA, diff(log(df$sp500)))
rv_window  <- 21
df$rv_spx <- rollapplyr(spx_logret, width = rv_window,
                        FUN = function(x) sd(x, na.rm = FALSE) * sqrt(252) * 100,
                        fill = NA)

vix_start <- as.Date("1990-01-02")
df$vix_extended <- ifelse(df$date < vix_start, df$rv_spx, df$vix)

# 4c. Stock-bond correlation: rolling 63d corr of SPX returns and bond-return
#     proxy (= -dY on the 10Y). Sign convention matches the textbook stock-bond
#     correlation: positive means stocks and bonds move together.
bond_ret_proxy <- -c(NA, diff(df$us10y))
corr_window    <- 63

pair <- cbind(spx_logret, bond_ret_proxy)
df$stock_bond_corr <- rollapplyr(
  pair, width = corr_window, by.column = FALSE,
  FUN = function(z) {
    ok <- complete.cases(z)
    if (sum(ok) >= 0.8 * corr_window) cor(z[ok, 1], z[ok, 2]) else NA_real_
  },
  fill = NA
)

# --- 5. Export --------------------------------------------------------------

out_path <- "market_data.csv"
write.csv(df, out_path, row.names = FALSE)

message(sprintf("Saved %d rows to %s. Range: %s to %s",
                nrow(df), out_path, min(df$date), max(df$date)))