# =========================================================
# 01_data.R — extended history version
# =========================================================
suppressPackageStartupMessages({
  library(fredr); library(tidyquant); library(frenchdata)
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(lubridate)
})
source("Data/APIKey.R")

# ---- DAILY FRED (used for stock-bond corr, VIX post-1990) ----
fred_daily_ids <- c(
  yield_10y_d = "DGS10",      # 1962+
  vix_d       = "VIXCLS"      # 1990+
)

# ---- MONTHLY FRED (longer history than daily counterparts) ----
fred_monthly_ids <- c(
  tbill_3m    = "TB3MS",      # 1934+  (monthly avg, 3m T-bill secondary market)
  yield_10y_m = "GS10",       # 1953+  (monthly avg, 10y Treasury)
  oil         = "WTISPLC",    # 1946+  (monthly WTI spot)
  copper      = "PCOPPUSDM"   # 1990+  (still the binding constraint)
)

pull_fred <- function(ids, start = "1925-01-01") {
  imap_dfr(ids, function(id, nm) {
    fredr(series_id = id, observation_start = as.Date(start)) %>%
      transmute(date, series = nm, value)
  })
}

fred_daily   <- pull_fred(fred_daily_ids)
fred_monthly <- pull_fred(fred_monthly_ids)

# ---- S&P 500 daily (1927+) ----
spx_daily <- tq_get("^GSPC", from = "1927-01-01", to = Sys.Date()) %>%
  transmute(date, series = "sp500", value = adjusted)

# ---- Fama-French 5 + Momentum ----
ff5_raw  <- download_french_data("Fama/French 5 Factors (2x3)")
mom_raw  <- download_french_data("Momentum Factor (Mom)")
parse_ff <- function(tbl, cols) {
  tbl %>% as_tibble() %>%
    mutate(date = ceiling_date(ymd(paste0(date, "01")), "month") - days(1)) %>%
    mutate(across(all_of(cols), ~ as.numeric(.) / 100)) %>%
    select(date, all_of(cols))
}
ff5 <- parse_ff(ff5_raw$subsets$data[[1]],
                c("Mkt-RF","SMB","HML","RMW","CMA","RF"))
mom <- parse_ff(mom_raw$subsets$data[[1]], "Mom")
ff_factors <- ff5 %>% inner_join(mom, by = "date") %>%
  rename(market=`Mkt-RF`, size=SMB, value=HML,
         profitability=RMW, investment=CMA, momentum=Mom, rf=RF) %>%
  select(date, market, size, value, profitability, investment, momentum, rf)

# ---- OPTIONAL: long copper from Refinitiv ----
# If you pull LME copper from Refinitiv Workspace as a monthly CSV with
# columns `date` (YYYY-MM-DD, month-end) and `copper`, save it to
# data/copper_long.csv and uncomment the block below. It will override
# the short FRED copper series in fred_monthly.
#
# if (file.exists("data/copper_long.csv")) {
#   long_copper <- readr::read_csv("data/copper_long.csv", show_col_types = FALSE) %>%
#     transmute(date = as.Date(date), series = "copper", value = copper)
#   fred_monthly <- fred_monthly %>% filter(series != "copper") %>%
#     bind_rows(long_copper)
# }

dir.create("data", showWarnings = FALSE)
saveRDS(
  list(fred_daily   = fred_daily,
       fred_monthly = fred_monthly,
       spx_daily    = spx_daily,
       ff_factors   = ff_factors),
  "data/regimes_raw.rds"
)
message("Saved data/regimes_raw.rds")