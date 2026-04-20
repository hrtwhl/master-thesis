# Download Bloomberg + FRED historical series and build derived regime variables.
# Requires: active Bloomberg terminal, Rblpapi, fredr (FRED API key), zoo.

library(Rblpapi)
library(fredr)
library(zoo)

fredr_set_key("263fff36e95136e168c2b9128597195d")   # or set FRED_API_KEY in ~/.Renviron

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
  oil         = list(ticker = "DCOILWTICO", start = "1986-01-01"),  # WTI Cushing spot
  us3m        = list(ticker = "DTB3",       start = "1954-01-01"),  # 3M T-bill secondary market
  yield_curve = list(ticker = "T10Y3M",     start = "1982-01-01")   # 10Y CMT minus 3M CMT
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

# --- 4. Derived variables ---------------------------------------------------

rv_window   <- 21
vix_start   <- as.Date("1990-01-02")
corr_window <- 63

# 4a. Extended VIX: realised vol computed on the clean SPX calendar
#     (avoids merge-gap NAs), switching to actual VIX from 1990-01-02.
spx_raw <- bbg_data$sp500[order(bbg_data$sp500$date), ]
spx_raw$logret <- c(NA, diff(log(spx_raw$sp500)))
spx_raw$rv_spx <- rollapplyr(
  spx_raw$logret, width = rv_window,
  FUN = function(x) sd(x) * sqrt(252) * 100,
  fill = NA
)

ext <- merge(spx_raw[, c("date", "rv_spx")], bbg_data$vix, by = "date", all = TRUE)
ext$vix_extended <- ifelse(ext$date < vix_start, ext$rv_spx, ext$vix)
df <- merge(df, ext[, c("date", "vix_extended")], by = "date", all.x = TRUE)

# 4b. Stock-bond correlation: rolling 63d corr(SPX returns, bond return proxy = -dY).
#     Computed on inner-joined SPX + US10Y calendar for the same reason.
sb <- merge(
  bbg_data$sp500[order(bbg_data$sp500$date), ],
  bbg_data$us10y[order(bbg_data$us10y$date), ],
  by = "date"
)
sb$spx_logret     <- c(NA, diff(log(sb$sp500)))
sb$bond_ret_proxy <- -c(NA, diff(sb$us10y))
sb$stock_bond_corr <- rollapplyr(
  cbind(sb$spx_logret, sb$bond_ret_proxy),
  width = corr_window, by.column = FALSE,
  FUN = function(z) {
    ok <- complete.cases(z)
    if (sum(ok) >= 0.8 * corr_window) cor(z[ok, 1], z[ok, 2]) else NA_real_
  },
  fill = NA
)
df <- merge(df, sb[, c("date", "stock_bond_corr")], by = "date", all.x = TRUE)

# --- 5. Final cleanup -------------------------------------------------------

df$us10y <- NULL
df <- df[order(df$date), ]
first_full <- min(df$date[complete.cases(df)])
df <- df[df$date >= first_full, ]

# --- 6. Export --------------------------------------------------------------

out_path <- "market_data.csv"
write.csv(df, out_path, row.names = FALSE)

message(sprintf("Saved %d rows to %s. Range: %s to %s",
                nrow(df), out_path, min(df$date), max(df$date)))