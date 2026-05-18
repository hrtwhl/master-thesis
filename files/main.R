# =============================================================================
#  main.R
#  Wasserstein HMM + MVO daily cross-asset allocation
#  Source: Boukardagha (2026), "Explainable Regime-Aware Investing"
#
#  Pipeline:
#    1. Load asset_data_yf.csv, compute log returns, build features
#    2. Train/test split at 2023-06-03
#    3. Compile the C++ kernel (one-time, ~30 seconds)
#    4. Run Wasserstein HMM + MVO backtest (strictly causal expanding window)
#    5. Run Equal-Weight, 60/40, SPX buy-and-hold benchmarks on the same dates
#    6. Compute diagnostics, save CSVs and PNGs into output/
#
#  --- Discrepancies between paper and notebook (resolved in this implementation):
#
#  (a) MONOTONE K. The paper says K is selected predictively each period and
#      "can go up/down". The notebook hard-codes a non-decreasing ratchet
#      (K_curr = max(K_curr, K_candidate)). This implementation follows the
#      PAPER. Reason: a non-decreasing ratchet ossifies K at the first big
#      validation slice; the predictive-selection machinery loses its purpose.
#
#  (b) REGIME MOMENTS ON FORWARD RETURNS. The paper text is silent on the
#      timing of r in the regime-conditional moments. The notebook uses
#      r_{s+1} given z_s, which is the right object for one-step-ahead
#      forecasting under MVO. This implementation follows the NOTEBOOK.
#
#  (c) KNN BENCHMARK. Dropped per instructions. Replaced with a daily-
#      rebalanced 60/40 stocks/bonds benchmark, evaluated on the same OOS
#      dates as the strategy.
#
#  (d) HMMLEARN VS CUSTOM RCPP. The C++ Gaussian HMM uses Baum-Welch with
#      K-means++ init. Exact numerical equivalence with hmmlearn is not
#      possible (different init RNG paths). Sharpe / drawdown should agree
#      with the paper within ~5%.
# =============================================================================

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(Matrix)
  library(osqp)
  library(caTools)
  library(abind)
})

# ---------------- Configuration ----------------

PROJ_DIR <- if (interactive()) getwd() else commandArgs(trailingOnly = FALSE) |>
  (\(x) dirname(sub("^--file=", "", x[grepl("^--file=", x)])))()
if (length(PROJ_DIR) == 0 || !nzchar(PROJ_DIR)) PROJ_DIR <- getwd()

DATA_PATH <- file.path(PROJ_DIR, "Data",  "asset_data.csv")
OUT_DIR   <- file.path(PROJ_DIR, "files/output_2")
CPP_FILE  <- file.path(PROJ_DIR, "files",   "wasserstein_hmm.cpp")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CFG <- list(
  asset_names   = c("SPX","BOND","GOLD","OIL","USD"),
  split_date    = "2005-01-01",
  # feature windows
  vol_window    = 60,
  mom_window    = 20,
  # HMM and K selection
  K_min         = 5,
  K_max         = 6,
  F_K           = 5,        # re-select K every 5 trading days
  L_val         = 252,      # validation slice for predictive K selection
  lambda_K      = 1.0,      # complexity penalty
  em_iter       = 300,
  seed          = 42,
  lookback      = 60,
  # templates
  G_max         = 8,
  eta_tpl       = 0.05,
  spawn_thresh  = 2.5,
  # MVO
  lambda        = 5.0,
  tc            = 0.0002,
  w_max         = 0.6
)

# ---------------- Compile C++ + source R modules ----------------

cat("Compiling C++ kernels...\n")
Rcpp::sourceCpp(CPP_FILE)

cat("Sourcing R modules...\n")
source(file.path(PROJ_DIR, "files", "01_data.R"))
source(file.path(PROJ_DIR, "files", "02_hmm_engine.R"))
source(file.path(PROJ_DIR, "files", "03_strategy.R"))
source(file.path(PROJ_DIR, "files", "04_diagnostics.R"))

# ---------------- Load data ----------------

cat("Loading data from ", DATA_PATH, "\n", sep = "")
prices  <- load_asset_data(DATA_PATH, asset_order = CFG$asset_names)
cat(sprintf("  rows=%d  dates=%s -> %s\n",
            nrow(prices), format(min(prices$date)), format(max(prices$date))))

returns_df  <- compute_log_returns(prices)
features_df <- build_features(returns_df, CFG$vol_window, CFG$mom_window)

# Realign returns_df to start at the first feature date so train/test are aligned
returns_df <- returns_df[returns_df$date >= min(features_df$date), ]

splits <- train_test_split(returns_df, features_df, CFG$split_date)
cat(sprintf("Train: %s -> %s  (n=%d)\n",
            format(min(splits$returns_train$date)),
            format(max(splits$returns_train$date)),
            nrow(splits$returns_train)))
cat(sprintf("Test : %s -> %s  (n=%d)\n",
            format(min(splits$returns_test$date)),
            format(max(splits$returns_test$date)),
            nrow(splits$returns_test)))

# ---------------- Run strategy ----------------

cat("\nRunning Wasserstein HMM + MVO backtest...\n")
t0 <- Sys.time()
results <- run_strategy(
  returns_train = splits$returns_train,
  returns_test  = splits$returns_test,
  features_all  = features_df,
  asset_names   = CFG$asset_names,
  cfg           = CFG,
  verbose       = TRUE
)
cat(sprintf("Done. Elapsed: %.1f sec\n", as.numeric(Sys.time() - t0, units = "secs")))

# ---------------- Benchmarks ----------------

test_dates <- results$pnl$date
benchmarks <- list(
  ew          = benchmark_equal_weight(returns_df, CFG$asset_names, test_dates),
  sixty_forty = benchmark_60_40      (returns_df, test_dates),
  spx         = benchmark_spx        (returns_df, test_dates),
  returns_test = splits$returns_test
)

# ---------------- Diagnostics: export ----------------

cat("\nExporting CSV tables...\n")
perf <- export_tables(results, benchmarks, CFG$asset_names, OUT_DIR)
print(perf, row.names = FALSE)

cat("\nTurnover stats (Wasserstein HMM):\n")
print(turnover_stats(results$weights, CFG$asset_names), row.names = FALSE)

cat("\nAllocation summary (Wasserstein HMM):\n")
print(allocation_summary(results$weights, CFG$asset_names), row.names = FALSE)

cat("\nConcentration (Wasserstein HMM):\n")
print(concentration_stats(results$weights, CFG$asset_names), row.names = FALSE)

cat("\nPer-regime portfolio performance:\n")
print(regime_performance(results$pnl, results$tpl_label), row.names = FALSE)

cat("\nPer-regime asset Sharpe:\n")
print(asset_sharpe_by_regime(splits$returns_test, results$tpl_label, CFG$asset_names),
      row.names = FALSE)

cat("\nExporting PNG plots...\n")
export_all_plots(results, benchmarks, CFG$asset_names, OUT_DIR)

cat(sprintf("\nAll outputs written to: %s\n", normalizePath(OUT_DIR)))

