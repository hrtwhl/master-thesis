# ---------------------------------------------------------------------------
# Replication of Mulliner, Harvey, Xia, Fang, van Hemert (2025), "Regimes"
#
# Entry point. Sources modules in order and runs the full pipeline.
# Run from the project root: source("main.R")
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(zoo)
  library(slider)
  library(patchwork)
  library(scales)
})

# ---- Configuration --------------------------------------------------------

CFG <- list(
  paths = list(
    macro_csv   = "data_final.csv",
    ff5_url     = "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_5_Factors_2x3_CSV.zip",
    mom_url     = "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Momentum_Factor_CSV.zip",
    output_dir  = "Mulliner/output"
  ),

  # End of sample (binding constraint from the user)
  sample_end = as.Date("2024-12-31"),

  # Z-score construction
  change_horizon  = 12L,    # 12-month change
  rolling_window  = 120L,   # 10 years of monthly 12m-diffs for rolling std
  winsorize_at    = 3,

  # Similarity
  mask_months   = 36L,     # exclude last 3 years when measuring similarity
  paper_top_pct = 0.15,    # 15% most similar for Exhibits 6-8

  # Strategy
  quantile_default = 5L,           # quintiles
  quantile_robust  = c(2, 3, 4, 5, 10, 20),
  lookback_robust  = c(12L, 36L, 60L),  # 1y, 3y, 5y z-score lookbacks
  vol_target       = 0.15,
  vol_window       = 36L,          # realised-vol window for vol targeting

  # EWMA half-life lookbacks (in months)
  ewma_lookbacks = c(12L, 24L, 36L, 48L),

  # Variable metadata. diff_type controls how the 12-month change is formed.
  # Use log differences for multiplicative / trending series (S&P, oil, copper,
  # VIX) so that the resulting z-scores are ~stationary and every variable
  # carries roughly equal weight in the Euclidean distance. Level differences
  # for rate- or correlation-valued series that are already stationary and can
  # cross zero.
  variables = tibble::tribble(
    ~var,              ~label,            ~diff_type,
    "sp500",           "Market",          "log",
    "yield_curve",     "Yield curve",     "level",
    "oil",             "Oil",             "log",
    "copper",          "Copper",          "log",
    "us3m",            "Monetary policy", "level",
    "vix",             "Volatility",      "log",
    "stock_bond_corr", "Stock-bond",      "level"
  )
)

dir.create(CFG$paths$output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Pipeline -------------------------------------------------------------

source("Mulliner/utils.R")
source("Mulliner/01_data.R")
source("Mulliner/02_similarity.R")
source("Mulliner/03_strategy.R")
source("Mulliner/04_exhibits.R")

message("\n[1/5] Loading macro data and building state variables ...")
macro <- load_macro(CFG$paths$macro_csv, CFG$sample_end)
state <- build_state_variables(macro, CFG$variables,
                               horizon   = CFG$change_horizon,
                               window    = CFG$rolling_window,
                               winsorize = CFG$winsorize_at)

message("[2/5] Downloading Fama-French factors and momentum ...")
factors <- load_factors(CFG$paths$ff5_url, CFG$paths$mom_url, CFG$sample_end)

message("[3/5] Computing similarity (Euclidean distance matrix) ...")
dist_mat <- compute_distance_matrix(state$transformed_winsorized)

message("[4/5] Running factor-timing strategies ...")
strat <- run_strategies(
  distance_matrix = dist_mat,
  factor_returns  = factors,
  mask_months     = CFG$mask_months,
  q_default       = CFG$quantile_default,
  q_robust        = CFG$quantile_robust,
  lookback_robust = CFG$lookback_robust,
  macro_monthly   = state$monthly,
  variables       = CFG$variables,
  horizon         = CFG$change_horizon,
  winsorize       = CFG$winsorize_at
)

message("[5/5] Producing exhibits ...")
make_all_exhibits(
  macro_monthly = state$monthly,
  transformed   = state$transformed,
  winsorized    = state$transformed_winsorized,
  distance_mat  = dist_mat,
  factors       = factors,
  strategies    = strat,
  cfg           = CFG
)

message("\nDone. Output written to '", CFG$paths$output_dir, "'.")


source("Mulliner/diagnostics.R")
run_diagnostics(state, dist_mat, factors, strat, CFG)
