# =========================================================
# PACKAGES
# =========================================================
required <- c("fredr", "dplyr", "tidyr", "purrr", "tidyquant", "stringr", "tibble")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(required, library, character.only = TRUE))

# =========================================================
# FRED API KEY
# =========================================================
# Your file should call: fredr_set_key("<YOUR_KEY>")
source("APIKey.R")

# =========================================================
# FRED SERIES (THE "CORE 6" MACRO LIST)
# =========================================================
# This list is based on our analysis of Heiden/AQR/Man
# and designed for low correlation and high EM-relevance.
fred_series <- c(

  monetary_policy  = "DGS3MO",
  yield_curve      = "T10Y3M",
  copper           = "PCOPPUSDM",
  vix_index        = "VIXCLS",
  oil              = "DCOILWTICO",
  yield_10y        = "DGS10"
)


# Helper: download FRED as long table
get_fred_data <- function(named_series, start_date = "2000-01-01", freq = NULL) {
  purrr::imap_dfr(named_series, function(ticker, nice_name) {
    fredr(
      series_id = ticker,
      observation_start = as.Date(start_date),
      frequency = freq
    ) |>
      transmute(date, series = nice_name, value)
  }) |>
    mutate(source = "FRED")
}

# =========================================================
# YAHOO FUTURES (THE "CORE 6" MACRO LIST)
# =========================================================
# 5. Global Growth (Man Group Theme)
# 6. Inflation Shock (Man Group Theme)
yf_futures <- tibble::tibble(
  yf = "^SPX",
  series = "SPX"  # Clean names for the series
)

get_yahoo_close <- function(tickers_tbl) {
  tq_get(
    tickers_tbl$yf,
    get  = "stock.prices",
    from = as.Date("1900-01-01"),
    to   = Sys.Date()
  ) |>
    transmute(date, symbol, value = close) |>
    left_join(tickers_tbl, by = c("symbol" = "yf")) |>
    transmute(date, series, value, source = "YF") |>
    arrange(series, date)
}

# ----------------------------
# Additional packages (if needed)
# ----------------------------
required_extra <- c("zoo", "lubridate", "readr")    # readr used for saving; replace with saveRDS if you prefer
to_install2 <- setdiff(required_extra, rownames(installed.packages()))
if (length(to_install2)) install.packages(to_install2)
invisible(lapply(required_extra, library, character.only = TRUE))

# ----------------------------
# 1) Download the series (if not already)
# ----------------------------
# fred_df: long table from get_fred_data(...) with columns date, series, value, source
fred_df <- get_fred_data(fred_series, start_date = "1990-01-01")  # you can change start_date
spx_df  <- get_yahoo_close(yf_futures)  # from your function; columns date, series (SPX), value, source

# ----------------------------
# 2) Extract the 10y yield (DGS10) and SPX and align daily
# ----------------------------
# isolate series names used in your fred_series vector
# fred_series has named element "yield_10y" -> "DGS10"
yield_name <- "yield_10y"  # this is how you named it in fred_series
spx_name   <- "SPX"

# wide merging: keep raw values, one row per date for either source (outer join)
wide_raw <- bind_rows(fred_df, spx_df) %>%
  select(date, series, value) %>%
  # keep only SPX & yield_10y for correlation calculation, but retain other series for macro_wide later
  pivot_wider(names_from = series, values_from = value)

# check we have the columns expected
if (!all(c(spx_name, yield_name) %in% colnames(wide_raw))) {
  stop("Either SPX or yield_10y not present in the combined dataset. Check series names.")
}

# ----------------------------
# 3) Create transformed series for correlation:
#    - spx_log_return: daily log returns of SPX
#    - d10y: daily first difference of 10y yield (absolute change)
# ----------------------------
corr_df <- wide_raw %>%
  arrange(date) %>%
  mutate(
    spx_log = !!rlang::sym(spx_name),
    yield10 = !!rlang::sym(yield_name)
  ) %>%
  # Compute transforms. Use log returns for SPX, simple diff for yield.
  mutate(
    spx_log_ret = ifelse(!is.na(spx_log), c(NA, diff(log(spx_log))), NA_real_),
    d10y = ifelse(!is.na(yield10), c(NA, diff(yield10)), NA_real_)
  ) %>%
  select(date, spx_log_ret, d10y)

# ----------------------------
# 4) Rolling 3-year correlation on daily data
#    - use window = 756 (3 * 252 trading days). Adjust if you prefer exact calendar days.
# ----------------------------
window_days <- 3 * 252  # approximate trading days in 3 years (756)
# prepare matrix for rollapply: need complete cases inside function
# Make a zoo object with date as the explicit time index
mat <- zoo(
  corr_df %>% select(spx_log_ret, d10y),
  order.by = corr_df$date
)

# define rolling correlation function with pairwise complete.obs
rolling_cor <- zoo::rollapplyr(
  data = mat,
  width = window_days,
  FUN = function(x) {
    # x is a matrix window_days x 2
    # return NA if too many NAs (optional)
    if (sum(complete.cases(x)) < floor(window_days * 0.5)) return(NA_real_)
    cor(x[,1], x[,2], use = "pairwise.complete.obs")
  },
  by.column = FALSE,
  fill = NA_real_
)

corr_series <- tibble(
  date = index(rolling_cor),
  spx_10y_corr = coredata(rolling_cor)
)

# ----------------------------
# 5) Merge correlation back into the full macro dataset (wide)
# ----------------------------
# create macro_wide: pivot the full fred + spx combined long dataset into wide,
# then left_join the corr_series; then drop the raw 10y yield column.
macro_wide <- bind_rows(fred_df, spx_df) %>%
  select(date, series, value) %>%
  pivot_wider(names_from = series, values_from = value) %>%
  # attach correlation series
  left_join(corr_series, by = "date") %>%
  # drop the raw 10y yield column now that correlation is added
  select(-all_of(yield_name))

# ----------------------------
# 6) Create long version
# ----------------------------
macro_long <- macro_wide %>%
  pivot_longer(-date, names_to = "series", values_to = "value")






