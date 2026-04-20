# main.R -- Oliveira, Sandfelder, Fujita, Dong, Cucuringu (2025) replication
#           on the paper's ETF universe (Table 2):
#
#     SPY + XLB, XLE, XLF, XLI, XLK, XLP, XLU, XLV, XLY   (10 assets)
#
# Run:
#   Rscript main.R
#   source("main.R")
#
# Outputs in ./output/ :
#   - PNG plots (figs 1-13 equivalent)
#   - performance_table.csv
#   - backtest_results.rds (full object for further analysis)

# ---- locate script directory -------------------------------------------------

script_dir <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- cmd_args[grep("^--file=", cmd_args)]
    if (length(file_arg))
      normalizePath(dirname(sub("^--file=", "", file_arg)))
    else
      normalizePath(getwd())
  }
)

source(file.path(script_dir, "01_data.R"))
source(file.path(script_dir, "02_regimes.R"))
source(file.path(script_dir, "03_models.R"))
source(file.path(script_dir, "04_backtest.R"))
source(file.path(script_dir, "05_plots.R"))

out_dir <- file.path(script_dir, "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- config (paper-consistent) ----------------------------------------------

config <- list(
  macro_start    = "1960-01-01",     # FRED-MD start
  etf_start      = "1998-12-01",     # 9 sector ETFs inception (SPY from 1993)
  backtest_start = as.Date("2001-01-01"),
  train_min      = 24L,              # months of prior returns for forecast fit
  k_inner        = 5L,               # paper r=5 inner regimes (+1 outlier = 6)
  pca_var        = 0.95,
  exclude_groups = 6L,               # FRED-MD group 6 = Interest & Exchange
  l_grid         = c(2, 3, 4),
  vol_target     = 0.10,
  refit_regime_every = 12L           # annual walk-forward regime refit with Hungarian matching
)

# ---- 1. Load data ------------------------------------------------------------

message("[1/5] Loading FRED-MD...")
fredmd <- load_fredmd_data(date_start = config$macro_start,
                           exclude_groups = config$exclude_groups)
message(sprintf("  %d variables x %d months",
                ncol(fredmd$data) - 1, nrow(fredmd$data)))

message("[2/5] Loading ETF returns from Yahoo Finance...")
etf <- load_etf_returns(tickers    = default_etf_tickers(),
                        date_start = config$etf_start)

common_start <- max(min(fredmd$data$date), min(etf$returns$date))
common_end   <- min(max(fredmd$data$date), max(etf$returns$date))
message(sprintf("  common window: %s to %s", common_start, common_end))

# ---- 2. In-sample regime model (for descriptive plots) ----------------------

message("[3/5] Fitting in-sample regime model...")
macro_for_plots <- fredmd$data |>
  filter(date >= common_start, date <= common_end)
pca_is <- fit_pca(macro_for_plots, var_thresh = config$pca_var)
message(sprintf("  PCA kept %d components (>= %.0f%% variance)",
                pca_is$k, config$pca_var * 100))

regime_is <- fit_regime_model(pca_is$scores, k_inner = config$k_inner)
class_is  <- classify_regimes(regime_is, pca_is$scores)
E_is      <- regime_transition_matrix(class_is$labels, config$k_inner + 1L)
E_is_cond <- conditional_transition_matrix(E_is)

message("  Fitting GMM for comparison...")
gmm_is <- fit_gmm_regimes(pca_is$scores, G = config$k_inner + 1L)

# ---- 3. Walk-forward backtest ------------------------------------------------

message("[4/5] Running walk-forward backtest...")
bt <- run_backtest(
  macro_raw          = fredmd$data,
  factor_rets        = etf$returns,       # parameter keeps name, now holds ETFs
  backtest_start     = config$backtest_start,
  train_min          = config$train_min,
  k_inner            = config$k_inner,
  pca_var            = config$pca_var,
  l_grid             = config$l_grid,
  vol_target         = config$vol_target,
  refit_regime_every = config$refit_regime_every,
  verbose            = TRUE
)

# ---- 4. Performance table ----------------------------------------------------

perf <- compute_performance_table(bt$portfolio_returns_scaled)
perf_sorted <- perf |> arrange(desc(Sharpe))
print(perf_sorted, n = nrow(perf_sorted))
write.csv(perf_sorted, file.path(out_dir, "performance_table.csv"),
          row.names = FALSE)

# ---- 5. Plots ----------------------------------------------------------------

message("[5/5] Generating plots...")

regime_label_map <- c(
  "0" = "Economic Difficulty",
  "1" = "Economic Recovery",
  "2" = "Expansionary Growth",
  "3" = "Stagflationary Pressure",
  "4" = "Pre-Recession Transition",
  "5" = "Reflationary Boom"
)

plots <- generate_all_plots(
  pca_fit = pca_is,
  dates   = as.Date(rownames(pca_is$scores)),
  km_labels  = class_is$labels,
  gmm_labels = gmm_is$labels,
  km_probs_R0  = class_is$probs[, "R0"],
  gmm_probs_R0 = if (!is.null(gmm_is$probs)) gmm_is$probs[, 1] else NULL,
  macro_raw_aligned = macro_for_plots,
  transition_matrix      = E_is,
  transition_matrix_cond = E_is_cond,
  regime_labels_table    = regime_label_map,
  portfolio_returns_scaled = bt$portfolio_returns_scaled,
  nber_recessions = NULL
)

save_plots(plots, out_dir = out_dir)

# ---- persist --------------------------------------------------------------

saveRDS(bt, file.path(out_dir, "backtest_results.rds"))
saveRDS(list(
  pca_fit = pca_is,
  regime_model = regime_is,
  labels = class_is$labels,
  probs  = class_is$probs,
  transition_matrix = E_is,
  gmm = gmm_is
), file.path(out_dir, "regime_insample.rds"))

message(sprintf("Done. Outputs in %s/", out_dir))
