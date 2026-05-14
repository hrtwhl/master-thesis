###############################################################################
# diagnostics_and_robustness.R — Sensitivity Analysis & Regime Diagnostics
# Boukardagha (2026) Replication
###############################################################################

library(parallel)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. SENSITIVITY ANALYSIS LOOPS
# ═══════════════════════════════════════════════════════════════════════════════

#' Run sensitivity analysis for a single parameter
#' @param feat        output of build_features()
#' @param param_name  one of "G", "eta", "gamma", "tau"
#' @param grid        numeric vector of values to test
#' @param n_cores     parallel cores
#' @return data.frame with metrics for each grid value
run_sensitivity <- function(feat, param_name, grid, n_cores = N_CORES) {
  
  cat(sprintf("[sensitivity] Running %s grid: %s\n",
              param_name, paste(grid, collapse = ", ")))
  
  run_one <- function(val) {
    cfg <- list(verbose = FALSE)
    cfg[[param_name]] <- val
    
    result <- tryCatch(
      run_wasserstein_hmm(feat, OOS_START_DATE, cfg),
      error = function(e) {
        cat(sprintf("  [WARN] %s=%s failed: %s\n", param_name, val, e$message))
        NULL
      }
    )
    
    if (is.null(result)) {
      return(data.frame(param = param_name, value = val,
                        Sharpe = NA, Max_DD = NA, Avg_Turnover = NA,
                        Sortino = NA, Ann_Return = NA, Ann_Vol = NA))
    }
    
    m <- compute_metrics(result)
    data.frame(param = param_name, value = val,
               Sharpe = m$Sharpe, Max_DD = m$Max_DD,
               Avg_Turnover = m$Avg_Turnover,
               Sortino = m$Sortino,
               Ann_Return = m$Ann_Return, Ann_Vol = m$Ann_Vol)
  }
  
  if (n_cores > 1 && length(grid) > 1) {
    cl <- makeCluster(min(n_cores, length(grid)))
    on.exit(stopCluster(cl))
    
    # Export required objects and functions to cluster
    clusterExport(cl, envir = globalenv(),
                  varlist = c("feat", "run_wasserstein_hmm", "compute_metrics",
                              "OOS_START_DATE", "ASSET_NAMES", "N_ASSETS",
                              "ASSET_LABELS", "G_TEMPLATES", "ETA_SMOOTH",
                              "GAMMA_RISK", "TAU_TCOST", "LAMBDA_K", "W_MAX",
                              "K_MIN", "K_MAX", "EM_MAX_ITER", "EM_TOL",
                              "COV_REG_FACTOR", "ORDER_SELECT_FREQ",
                              "VAL_SLICE_DAYS", "HMM_FIT_FREQ", "N_CORES",
                              "HMM_INIT_WINDOW", "SEED",
                              "fit_gaussian_hmm", "hmm_init", "hmm_filter_step",
                              "dmvnorm_log",
                              "log_sum_exp", "log_sum_exp_cols",
                              "select_model_order",
                              "wasserstein2_gaussian", "init_templates",
                              "assign_to_templates", "update_templates",
                              "compute_conditional_moments",
                              "ledoit_wolf_shrink", "shrink_cov",
                              "psd_project", "solve_mvo",
                              "LEDOIT_WOLF"))
    # Compile C++ in each worker
    cpp_path <- normalizePath("hmm_rcpp.cpp")
    clusterExport(cl, "cpp_path", envir = environment())
    clusterEvalQ(cl, Rcpp::sourceCpp(cpp_path))
    clusterEvalExpr <- function(cl) {
      clusterEvalQ(cl, { library(MASS); library(Matrix) })
    }
    tryCatch(clusterEvalExpr(cl), error = function(e) NULL)
    
    results <- parLapply(cl, grid, run_one)
  } else {
    results <- lapply(grid, run_one)
  }
  
  do.call(rbind, results)
}

#' Run all sensitivity analyses
run_all_sensitivity <- function(feat) {
  sens_results <- list()
  
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("  SENSITIVITY ANALYSIS\n")
  cat("══════════════════════════════════════════════════════════════\n\n")
  
  # G (templates)
  sens_results$G <- run_sensitivity(feat, "G", SENS_G_GRID, n_cores = 1)
  
  # eta (smoothing)
  sens_results$eta <- run_sensitivity(feat, "eta", SENS_ETA_GRID, n_cores = 1)
  
  # gamma (risk aversion)
  sens_results$gamma <- run_sensitivity(feat, "gamma", SENS_GAMMA_GRID, n_cores = 1)
  
  # tau (transaction costs)
  sens_results$tau <- run_sensitivity(feat, "tau", SENS_TAU_GRID, n_cores = 1)
  
  all_sens <- do.call(rbind, sens_results)
  
  # Save to CSV
  write.csv(all_sens, file.path(DIAG_DIR, "sensitivity_results.csv"),
            row.names = FALSE)
  cat(sprintf("[sensitivity] Results saved to %s/sensitivity_results.csv\n", DIAG_DIR))
  
  all_sens
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. REGIME FAILURE EPISODE DIAGNOSTIC
# ═══════════════════════════════════════════════════════════════════════════════

#' Identify "Regime Failure Episodes" — periods where the model assigns high
#' probability to a template but realized portfolio returns are significantly
#' negative. This is the "Regime 4 diagnostic" from the user's prior work.
#'
#' @param result    output of run_wasserstein_hmm()
#' @param prob_threshold  minimum template probability to consider "high confidence"
#' @param dd_window       rolling window for drawdown assessment (days)
#' @param dd_threshold    max-drawdown threshold to flag as "failure" (negative)
#' @return data.frame of failure episodes
identify_regime_failures <- function(result,
                                      prob_threshold = 0.60,
                                      dd_window = 21L,
                                      dd_threshold = -0.02) {
  
  dates <- result$dates
  ret   <- result$returns
  dom   <- result$dominant_regime
  tprob <- result$template_probs
  n     <- length(dates)
  
  episodes <- data.frame()
  
  for (t in (dd_window + 1):n) {
    if (is.na(dom[t]) || is.na(ret[t])) next
    
    # Check if dominant regime has high confidence
    max_prob <- max(tprob[t, ], na.rm = TRUE)
    if (max_prob < prob_threshold) next
    
    # Rolling drawdown over past dd_window days
    window_ret <- ret[(t - dd_window + 1):t]
    cum <- cumsum(window_ret)
    dd <- min(cum - cummax(cum))
    
    if (dd < dd_threshold) {
      episodes <- rbind(episodes, data.frame(
        date              = dates[t],
        dominant_regime   = dom[t],
        regime_prob       = max_prob,
        rolling_dd        = dd,
        daily_return      = ret[t],
        cum_return_window = sum(window_ret),
        stringsAsFactors  = FALSE
      ))
    }
  }
  
  if (nrow(episodes) > 0) {
    # Deduplicate overlapping episodes: keep only the worst day per regime-cluster
    episodes <- episodes[order(episodes$rolling_dd), ]
    
    write.csv(episodes,
              file.path(DIAG_DIR, "regime_failure_episodes.csv"),
              row.names = FALSE)
    cat(sprintf("[diagnostics] %d regime failure episodes identified and saved.\n",
                nrow(episodes)))
  } else {
    cat("[diagnostics] No regime failure episodes found.\n")
  }
  
  episodes
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. MACRO-INFORMATION GAP
# ═══════════════════════════════════════════════════════════════════════════════

#' Compute the "Macro-Information Gap": correlation between "wrong template"
#' assignments and a proxy for economic stress (realized volatility).
#'
#' A "wrong template" event is defined as a day where:
#'   1. The dominant regime has high probability (confident assignment)
#'   2. Realized returns are in the opposite direction of the template's
#'      conditional expected return sign
#'
#' Economic stress proxy: 21-day realized volatility of the stock return.
#'
#' @param result   output of run_wasserstein_hmm()
#' @return list with correlation, time series, and summary
compute_macro_info_gap <- function(result) {
  
  dates <- result$dates
  ret   <- result$returns
  dom   <- result$dominant_regime
  tprob <- result$template_probs
  asset_ret <- result$asset_returns
  templates <- result$templates
  n <- length(dates)
  
  if (is.null(templates) || is.null(dom)) {
    cat("[diagnostics] Cannot compute macro-info gap: no template data.\n")
    return(NULL)
  }
  
  n_assets <- ncol(asset_ret)
  
  # Compute 21-day realized vol of stocks (economic stress proxy)
  stock_col <- which(ASSET_NAMES == "stocks")
  rv <- rep(NA_real_, n)
  for (t in 22:n) {
    rv[t] <- sd(asset_ret[(t - 20):t, stock_col]) * sqrt(252)
  }
  
  # Identify "wrong template" days
  wrong_template <- rep(NA_real_, n)
  
  for (t in seq_len(n)) {
    if (is.na(dom[t]) || is.na(ret[t])) next
    g <- dom[t]
    if (g < 1 || g > length(templates)) next
    
    # Template's expected return sign (from return sub-block)
    tmpl_mu <- templates[[g]]$mu[1:n_assets]
    expected_port_ret <- mean(tmpl_mu)  # simple average across assets
    
    # "Wrong" if signs disagree and confidence is high
    max_prob <- max(tprob[t, ], na.rm = TRUE)
    if (max_prob > 0.5) {
      wrong_template[t] <- as.numeric(sign(ret[t]) != sign(expected_port_ret))
    }
  }
  
  # Rolling "wrong template" rate (21-day)
  wrong_rate <- rep(NA_real_, n)
  for (t in 22:n) {
    window <- wrong_template[(t - 20):t]
    if (sum(!is.na(window)) > 10) {
      wrong_rate[t] <- mean(window, na.rm = TRUE)
    }
  }
  
  # Correlation between wrong-template rate and realized vol
  valid <- !is.na(wrong_rate) & !is.na(rv)
  if (sum(valid) < 30) {
    cat("[diagnostics] Insufficient data for macro-info gap.\n")
    return(NULL)
  }
  
  corr <- cor(wrong_rate[valid], rv[valid])
  rank_corr <- cor(wrong_rate[valid], rv[valid], method = "spearman")
  
  cat(sprintf("[diagnostics] Macro-Information Gap:\n"))
  cat(sprintf("  Pearson correlation (wrong-template rate vs. realized vol): %.4f\n", corr))
  cat(sprintf("  Spearman rank correlation: %.4f\n", rank_corr))
  cat(sprintf("  Interpretation: %.0f%% of variance in template errors is\n", corr^2 * 100))
  cat(sprintf("    linearly associated with equity market stress.\n"))
  
  gap_df <- data.frame(
    date = dates[valid],
    wrong_rate = wrong_rate[valid],
    realized_vol = rv[valid]
  )
  
  write.csv(gap_df, file.path(DIAG_DIR, "macro_info_gap.csv"), row.names = FALSE)
  
  list(
    pearson  = corr,
    spearman = rank_corr,
    data     = gap_df
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. WEIGHT STABILITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

compute_weight_diagnostics <- function(result) {
  w <- result$weights
  if (is.null(w)) return(NULL)
  
  n <- nrow(w)
  n_assets <- ncol(w)
  
  # Average weights
  avg_wt <- colMeans(w, na.rm = TRUE)
  
  # Weight volatility
  wt_vol <- apply(w, 2, sd, na.rm = TRUE)
  
  # Fraction of time > 10%
  time_above_10 <- colMeans(w > 0.10, na.rm = TRUE)
  
  # Average absolute daily weight change
  dw <- abs(diff(w))
  avg_dw <- colMeans(dw, na.rm = TRUE)
  
  # Effective N (HHI inverse)
  neff <- 1 / rowSums(w^2, na.rm = TRUE)
  
  weight_table <- data.frame(
    Asset       = ASSET_LABELS,
    Avg_Weight  = avg_wt,
    Weight_Vol  = wt_vol,
    Time_Above_10pct = time_above_10,
    Avg_Abs_dW  = avg_dw,
    stringsAsFactors = FALSE
  )
  
  list(
    table     = weight_table,
    avg_neff  = mean(neff, na.rm = TRUE),
    med_neff  = median(neff, na.rm = TRUE),
    neff_ts   = neff
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. COMPREHENSIVE DIAGNOSTIC REPORT
# ═══════════════════════════════════════════════════════════════════════════════

run_diagnostics <- function(result, feat) {
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("  DIAGNOSTIC ANALYSIS\n")
  cat("══════════════════════════════════════════════════════════════\n\n")
  
  # Performance metrics
  metrics <- compute_metrics(result)
  cat("─── Overall Performance ───\n")
  print(metrics, row.names = FALSE)
  
  # Regime-conditional performance
  reg_met <- compute_regime_metrics(result)
  if (!is.null(reg_met)) {
    cat("\n─── Portfolio Performance by Regime ───\n")
    print(reg_met$portfolio, row.names = FALSE)
    cat("\n─── Asset Sharpe by Regime ───\n")
    print(reg_met$asset_sharpe, row.names = FALSE, digits = 3)
  }
  
  # Weight diagnostics
  wt_diag <- compute_weight_diagnostics(result)
  if (!is.null(wt_diag)) {
    cat("\n─── Weight Diagnostics ───\n")
    print(wt_diag$table, row.names = FALSE, digits = 4)
    cat(sprintf("  Avg N_eff: %.2f | Median N_eff: %.2f\n",
                wt_diag$avg_neff, wt_diag$med_neff))
  }
  
  # Regime failures
  cat("\n─── Regime Failure Episodes ───\n")
  failures <- identify_regime_failures(result)
  
  # Macro-information gap
  cat("\n─── Macro-Information Gap ───\n")
  gap <- compute_macro_info_gap(result)
  
  list(
    metrics   = metrics,
    regime    = reg_met,
    weights   = wt_diag,
    failures  = failures,
    macro_gap = gap
  )
}


cat("[diagnostics_and_robustness.R] Functions loaded.\n")
