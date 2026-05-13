###############################################################################
# main.R — Orchestrator for Boukardagha (2026) Replication
# Wasserstein HMM Regime-Aware Portfolio Construction
###############################################################################
#
# USAGE:
#   Rscript main.R                                # full run (uses cache if valid)
#   Rscript main.R --skip-sensitivity             # skip slow sensitivity loops
#   Rscript main.R --force-rerun                  # ignore cache, recompute everything
#   Rscript main.R --force-rerun --skip-sensitivity
#
# CACHING:
#   After the first run, backtest results are saved to cache/.
#   Subsequent runs with identical parameters + data skip the backtest
#   and go straight to diagnostics and plots (~seconds).
#   Change any parameter in config.R or update asset_data.csv → auto-invalidates.
#   Set CACHE_ENABLED <- FALSE in config.R to disable entirely.
#
# OUTPUT:
#   output/      — PNG charts
#   diagnostics/ — CSV diagnostic tables
#   cache/       — cached backtest results (.rds)
#
###############################################################################

cat("
\u2554\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2557
\u2551  Boukardagha (2026) Replication                                \u2551
\u2551  Wasserstein HMM Regime-Aware Portfolio Construction           \u2551
\u255a\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u255d
")

# ── Parse CLI args ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
SKIP_SENSITIVITY <- "--skip-sensitivity" %in% args
FORCE_RERUN      <- "--force-rerun" %in% args

# ── Set working directory to script location ─────────────────────────────────
if (interactive()) {
  # When sourced interactively, assume we're already in the project dir
} else {
  script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  if (script_dir != ".") setwd(script_dir)
}

# ── Source modules ───────────────────────────────────────────────────────────
t_total <- proc.time()

cat("\n[main] Loading modules...\n")
source("config.R")
source("data_processor.R")
source("model_functions.R")
source("backtest_engine.R")
source("diagnostics_and_robustness.R")
source("visualizer.R")

# Ensure output dirs exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DIAG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)


# ═══════════════════════════════════════════════════════════════════════════════
# CACHE INFRASTRUCTURE
# ═══════════════════════════════════════════════════════════════════════════════

#' Build a fingerprint list from all computation-affecting parameters + data.
#' Any change in this list invalidates the cache.
build_fingerprint <- function(data_path) {
  # Data identity: md5 of the CSV file
  data_md5 <- tools::md5sum(data_path)
  names(data_md5) <- NULL
  
  list(
    data_md5          = data_md5,
    # Feature params
    VOL_WINDOW        = VOL_WINDOW,
    MEAN_WINDOW       = MEAN_WINDOW,
    LOG_RETURN_CAP    = LOG_RETURN_CAP,
    ASSET_NAMES       = ASSET_NAMES,
    # Backtest params
    OOS_START_DATE    = as.character(OOS_START_DATE),
    K_MIN             = K_MIN,
    K_MAX             = K_MAX,
    EM_MAX_ITER       = EM_MAX_ITER,
    EM_TOL            = EM_TOL,
    COV_REG_FACTOR    = COV_REG_FACTOR,
    MAX_HMM_WINDOW    = MAX_HMM_WINDOW,
    ORDER_SELECT_FREQ = ORDER_SELECT_FREQ,
    VALIDATION_LEN    = VALIDATION_LEN,
    LAMBDA_K          = LAMBDA_K,
    G_TEMPLATES       = G_TEMPLATES,
    ETA_SMOOTH        = ETA_SMOOTH,
    GAMMA_RISK        = GAMMA_RISK,
    TAU_TCOST         = TAU_TCOST,
    W_MAX             = W_MAX,
    REFIT_FREQ        = REFIT_FREQ,
    # Benchmark params
    EW_REBAL_FREQ     = EW_REBAL_FREQ,
    BENCH_6040_ALLOC  = BENCH_6040_ALLOC
  )
}

#' Build a separate fingerprint for sensitivity analysis (extends backtest fp)
build_sensitivity_fingerprint <- function(backtest_fp) {
  c(backtest_fp, list(
    SENS_G_GRID     = SENS_G_GRID,
    SENS_ETA_GRID   = SENS_ETA_GRID,
    SENS_GAMMA_GRID = SENS_GAMMA_GRID,
    SENS_TAU_GRID   = SENS_TAU_GRID
  ))
}

#' Check if a cached fingerprint matches the current one
cache_is_valid <- function(cache_dir, current_fp, label = "backtest") {
  fp_path <- file.path(cache_dir, paste0(label, "_fingerprint.rds"))
  if (!file.exists(fp_path)) return(FALSE)
  
  saved_fp <- tryCatch(readRDS(fp_path), error = function(e) NULL)
  if (is.null(saved_fp)) return(FALSE)
  
  identical(saved_fp, current_fp)
}

#' Save results + fingerprint to cache
cache_save <- function(cache_dir, results_list, fingerprint, label = "backtest") {
  fp_path   <- file.path(cache_dir, paste0(label, "_fingerprint.rds"))
  data_path <- file.path(cache_dir, paste0(label, "_results.rds"))
  
  saveRDS(fingerprint, fp_path)
  saveRDS(results_list, data_path)
  
  size_mb <- file.info(data_path)$size / 1024^2
  cat(sprintf("[cache] Saved %s results (%.1f MB) to %s/\n", label, size_mb, cache_dir))
}

#' Load cached results
cache_load <- function(cache_dir, label = "backtest") {
  data_path <- file.path(cache_dir, paste0(label, "_results.rds"))
  if (!file.exists(data_path)) return(NULL)
  
  results <- tryCatch(readRDS(data_path), error = function(e) NULL)
  if (!is.null(results)) {
    size_mb <- file.info(data_path)$size / 1024^2
    cat(sprintf("[cache] Loaded %s results (%.1f MB) from %s/\n", label, size_mb, cache_dir))
  }
  results
}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: DATA LOADING AND FEATURE ENGINEERING
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("  STEP 1: DATA & FEATURES\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

price_data <- load_prices(DATA_PATH)
feat <- build_features(price_data$dates, price_data$prices)

cat(sprintf("\n  Date range: %s to %s\n", min(feat$dates), max(feat$dates)))
cat(sprintf("  Feature dim: %d (3 x %d assets)\n", ncol(feat$features), N_ASSETS))
cat(sprintf("  Valid feature rows: %d / %d (%.1f%%)\n",
            sum(feat$valid), length(feat$valid),
            mean(feat$valid) * 100))


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 + 3: BACKTEST + BENCHMARKS (with cache)
# ═══════════════════════════════════════════════════════════════════════════════

current_fp <- build_fingerprint(DATA_PATH)
use_cache  <- CACHE_ENABLED && !FORCE_RERUN && cache_is_valid(CACHE_DIR, current_fp)

if (use_cache) {
  cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  cat("  STEPS 2-3: LOADING CACHED RESULTS (parameters + data unchanged)\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")
  
  cached <- cache_load(CACHE_DIR, "backtest")
  hmm_result <- cached$hmm_result
  sf_result  <- cached$sf_result
  spx_result <- cached$spx_result
  
  cat(sprintf("  Cached OOS period: %s to %s (%d days)\n",
              min(hmm_result$dates), max(hmm_result$dates), length(hmm_result$dates)))
  
} else {
  if (FORCE_RERUN) {
    cat("\n[main] --force-rerun flag set, ignoring cache.\n")
  } else if (!CACHE_ENABLED) {
    cat("\n[main] Cache disabled (CACHE_ENABLED = FALSE).\n")
  } else {
    cat("\n[main] Cache miss (parameters or data changed). Running full backtest.\n")
  }
  
  # ── STEP 2: HMM BACKTEST ──
  cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  cat("  STEP 2: WASSERSTEIN HMM BACKTEST\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")
  
  hmm_result <- run_wasserstein_hmm(feat, OOS_START_DATE)
  
  # ── STEP 3: BENCHMARKS ──
  cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  cat("  STEP 3: BENCHMARKS\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")
  
  sf_result  <- run_sixty_forty(feat, OOS_START_DATE)
  spx_result <- run_spx_buyhold(feat, OOS_START_DATE)
  
  cat("[benchmarks] 60/40 and SPX Buy & Hold computed.\n")
  
  # ── Save to cache ──
  if (CACHE_ENABLED) {
    cache_save(CACHE_DIR,
               list(hmm_result = hmm_result,
                    sf_result  = sf_result,
                    spx_result = spx_result),
               current_fp, "backtest")
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: PERFORMANCE COMPARISON TABLE
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("  STEP 4: PERFORMANCE COMPARISON\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

metrics_hmm <- compute_metrics(hmm_result)
metrics_sf  <- compute_metrics(sf_result)
metrics_spx <- compute_metrics(spx_result)

comp_table <- rbind(
  data.frame(Strategy = "Wasserstein HMM", metrics_hmm),
  data.frame(Strategy = "60/40 Stock/Bond", metrics_sf),
  data.frame(Strategy = "SPX Buy & Hold", metrics_spx)
)

print(comp_table, row.names = FALSE, digits = 4)

write.csv(comp_table,
          file.path(DIAG_DIR, "performance_comparison.csv"),
          row.names = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: DIAGNOSTICS (always re-run — cheap)
# ═══════════════════════════════════════════════════════════════════════════════

diag_results <- run_diagnostics(hmm_result, feat)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6: VISUALIZATION (always re-run — cheap, and the point of caching)
# ═══════════════════════════════════════════════════════════════════════════════

reg_metrics <- compute_regime_metrics(hmm_result)

generate_all_plots(
  hmm_result, sf_result, spx_result,
  reg_metrics = reg_metrics,
  sens_df = NULL  # filled below if sensitivity runs
)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7: SENSITIVITY ANALYSIS (optional, slow — separately cached)
# ═══════════════════════════════════════════════════════════════════════════════

if (!SKIP_SENSITIVITY) {
  sens_fp      <- build_sensitivity_fingerprint(current_fp)
  sens_cached  <- CACHE_ENABLED && !FORCE_RERUN &&
                  cache_is_valid(CACHE_DIR, sens_fp, "sensitivity")
  
  if (sens_cached) {
    cat("\n[main] Loading cached sensitivity results...\n")
    sens_df <- cache_load(CACHE_DIR, "sensitivity")$sens_df
  } else {
    cat("\n[main] Running sensitivity analysis (this may take a while)...\n")
    cat("[main] To skip, run with: Rscript main.R --skip-sensitivity\n\n")
    
    sens_df <- run_all_sensitivity(feat)
    
    if (CACHE_ENABLED) {
      cache_save(CACHE_DIR, list(sens_df = sens_df), sens_fp, "sensitivity")
    }
  }
  
  plot_sensitivity(sens_df)
} else {
  cat("\n[main] Sensitivity analysis SKIPPED (--skip-sensitivity flag).\n")
  sens_df <- NULL
}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: EXPORT SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("  SUMMARY\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

elapsed <- (proc.time() - t_total)[3]
cat(sprintf("  Total runtime: %.1f seconds (%.1f minutes)\n", elapsed, elapsed / 60))
if (use_cache) {
  cat("  (backtest loaded from cache — computation was skipped)\n")
}
cat(sprintf("  OOS period: %s to %s (%d trading days)\n",
            min(hmm_result$dates), max(hmm_result$dates), length(hmm_result$dates)))
cat(sprintf("\n  -- Key Results --\n"))
cat(sprintf("  Wasserstein HMM Sharpe:  %.3f\n", metrics_hmm$Sharpe))
cat(sprintf("  Wasserstein HMM Max DD:  %.2f%%\n", metrics_hmm$Max_DD * 100))
cat(sprintf("  60/40 Sharpe:            %.3f\n", metrics_sf$Sharpe))
cat(sprintf("  SPX B&H Sharpe:          %.3f\n", metrics_spx$Sharpe))

cat(sprintf("\n  Output directory: %s/\n", OUTPUT_DIR))
cat(sprintf("  Diagnostics:     %s/\n", DIAG_DIR))
if (CACHE_ENABLED) {
  cat(sprintf("  Cache directory:  %s/\n", CACHE_DIR))
  cat("  (re-run with --force-rerun to recompute from scratch)\n")
}

# Save workspace for interactive inspection
save.image(file.path(OUTPUT_DIR, "workspace.RData"))
cat(sprintf("  Workspace saved: %s/workspace.RData\n", OUTPUT_DIR))

cat("\n[main] Done.\n")
