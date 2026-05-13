###############################################################################
# backtest_engine.R — Backtest Loop for Wasserstein HMM + Benchmarks
# Boukardagha (2026) Replication
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 1. MAIN WASSERSTEIN HMM STRATEGY
# ═══════════════════════════════════════════════════════════════════════════════

run_wasserstein_hmm <- function(feat, oos_start = OOS_START_DATE, cfg = list()) {
  
  # Override config with any sensitivity params
  g_tmpl    <- if (!is.null(cfg$G))     cfg$G     else G_TEMPLATES
  eta       <- if (!is.null(cfg$eta))   cfg$eta   else ETA_SMOOTH
  gamma_ra  <- if (!is.null(cfg$gamma)) cfg$gamma else GAMMA_RISK
  tau_tc    <- if (!is.null(cfg$tau))   cfg$tau   else TAU_TCOST
  lambda_k  <- if (!is.null(cfg$lambda_k)) cfg$lambda_k else LAMBDA_K
  refit_f   <- if (!is.null(cfg$refit_freq)) cfg$refit_freq else REFIT_FREQ
  k_min     <- if (!is.null(cfg$k_min)) cfg$k_min else K_MIN
  k_max     <- if (!is.null(cfg$k_max)) cfg$k_max else K_MAX
  w_max     <- if (!is.null(cfg$w_max)) cfg$w_max else W_MAX
  verbose   <- if (!is.null(cfg$verbose)) cfg$verbose else TRUE
  
  dates    <- feat$dates
  returns  <- feat$returns
  features <- feat$features
  valid    <- feat$valid
  n_days   <- length(dates)
  n_assets <- ncol(returns)
  
  oos_idx <- which(dates >= oos_start)
  if (length(oos_idx) == 0) stop("OOS start date is beyond data range.")
  t0 <- oos_idx[1]
  
  if (verbose) cat(sprintf("[backtest] OOS starts at index %d (%s), %d OOS days\n",
                           t0, dates[t0], n_days - t0 + 1))
  
  # ── Storage ────────────────────────────────────────────────────────────────
  weights_mat   <- matrix(NA_real_, n_days, n_assets)
  colnames(weights_mat) <- ASSET_NAMES
  port_returns  <- rep(NA_real_, n_days)
  dominant_reg  <- rep(NA_integer_, n_days)
  template_prob_mat <- matrix(NA_real_, n_days, g_tmpl)
  current_K     <- rep(NA_integer_, n_days)
  turnover      <- rep(NA_real_, n_days)
  w2_distances  <- rep(NA_real_, n_days)
  
  # ── State ──────────────────────────────────────────────────────────────────
  # Templates live in ORIGINAL scale (persist across refits where normalization changes)
  templates      <- NULL
  # HMM params in STANDARDIZED (z) scale — used for daily forward filtering
  hmm_params_z   <- NULL
  alpha_filt     <- NULL
  # Standardization params from last refit (needed to z-score new observations)
  cur_feat_means <- NULL
  cur_feat_sds   <- NULL
  
  current_Kt     <- k_min
  w_prev         <- rep(1 / n_assets, n_assets)
  days_since_refit   <- refit_f   # force refit on first OOS day
  days_since_order   <- ORDER_SELECT_FREQ
  initialized    <- FALSE
  
  # ── Main loop ──────────────────────────────────────────────────────────────
  t_start <- proc.time()
  
  for (t in t0:n_days) {
    if (!valid[t]) {
      weights_mat[t, ] <- w_prev
      port_returns[t] <- sum(w_prev * returns[t, ])
      next
    }
    
    # Expanding history up to t-1 (strict causality)
    hist_idx <- which(valid[1:(t - 1)])
    if (length(hist_idx) < 100) {
      weights_mat[t, ] <- w_prev
      port_returns[t] <- sum(w_prev * returns[t, ])
      next
    }
    
    X_hist <- features[hist_idx, , drop = FALSE]
    
    # Standardize features for HMM numerical stability
    feat_means <- colMeans(X_hist)
    feat_sds   <- apply(X_hist, 2, sd)
    feat_sds[feat_sds < 1e-10] <- 1e-10
    X_hist_z <- scale(X_hist, center = feat_means, scale = feat_sds)
    
    # ── Model-order selection (periodic) ───────────────────────────────────
    days_since_order <- days_since_order + 1
    if (days_since_order >= ORDER_SELECT_FREQ && nrow(X_hist_z) > VALIDATION_LEN + 100) {
      order_result <- select_model_order(X_hist_z, VALIDATION_LEN,
                                          k_min, k_max, lambda_k)
      current_Kt <- order_result$K
      days_since_order <- 0
    }
    
    # ── HMM fitting (periodic) or filtering (daily) ───────────────────────
    days_since_refit <- days_since_refit + 1
    
    if (days_since_refit >= refit_f || !initialized) {
      # ==== FULL HMM RE-FIT ====
      max_window <- min(nrow(X_hist_z), MAX_HMM_WINDOW)
      X_fit <- X_hist_z[(nrow(X_hist_z) - max_window + 1):nrow(X_hist_z), , drop = FALSE]
      
      hmm_fit <- tryCatch(
        fit_gaussian_hmm(X_fit, current_Kt, params = hmm_params_z),
        error = function(e) NULL
      )
      
      if (!is.null(hmm_fit)) {
        # Store z-space params for daily filtering
        hmm_params_z <- list(
          A = hmm_fit$A, mu = hmm_fit$mu, Sigma = hmm_fit$Sigma,
          K = hmm_fit$K, d = hmm_fit$d
        )
        cur_feat_means <- feat_means
        cur_feat_sds   <- feat_sds
        alpha_filt <- hmm_fit$filtered[nrow(hmm_fit$filtered), ]
        
        # De-standardize for template tracking (original scale)
        # mu_orig = mu_z * sd + mean; Sigma_orig = diag(sd) Sigma_z diag(sd)
        hmm_fit_orig <- hmm_fit
        sd_diag <- diag(cur_feat_sds)
        for (kk in seq_len(hmm_fit$K)) {
          hmm_fit_orig$mu[[kk]] <- hmm_fit$mu[[kk]] * cur_feat_sds + cur_feat_means
          hmm_fit_orig$Sigma[[kk]] <- sd_diag %*% hmm_fit$Sigma[[kk]] %*% sd_diag
        }
        
        # Template init / mapping (all in original scale)
        if (!initialized) {
          templates <- init_templates(hmm_fit_orig, g_tmpl)
          initialized <- TRUE
          if (verbose) cat(sprintf("  [%s] Templates initialized (K=%d, G=%d)\n",
                                   dates[t], current_Kt, g_tmpl))
        }
        
        map_result <- map_components_to_templates(hmm_fit_orig, templates)
        templates <- update_templates(templates, map_result, eta)
        
        template_prob_mat[t, ] <- map_result$template_probs
        dominant_reg[t] <- which.max(map_result$template_probs)
        w2_distances[t] <- min(map_result$dist_mat)
        
        days_since_refit <- 0
      }
      
    } else if (initialized && !is.null(hmm_params_z)) {
      # ==== DAILY FILTER UPDATE (no re-fit) ====
      # Standardize today's feature using last refit's normalization
      x_t_z <- (features[t, ] - cur_feat_means) / cur_feat_sds
      
      # Forward-filter one step in z-space
      alpha_filt <- hmm_filter_step(alpha_filt, x_t_z, hmm_params_z)
      
      # Map filtered probs to templates
      # Need original-scale component params for W2 distances
      K_cur <- hmm_params_z$K
      G_cur <- length(templates)
      sd_diag <- diag(cur_feat_sds)
      
      dist_mat <- matrix(Inf, K_cur, G_cur)
      for (k in seq_len(K_cur)) {
        mu_k_orig  <- hmm_params_z$mu[[k]] * cur_feat_sds + cur_feat_means
        Sig_k_orig <- sd_diag %*% hmm_params_z$Sigma[[k]] %*% sd_diag
        for (g in seq_len(G_cur)) {
          dist_mat[k, g] <- wasserstein2_gaussian(
            mu_k_orig, Sig_k_orig,
            templates[[g]]$mu, templates[[g]]$Sigma
          )
        }
      }
      mapping <- apply(dist_mat, 1, which.min)
      
      tprobs <- rep(0, G_cur)
      for (k in seq_len(K_cur)) {
        tprobs[mapping[k]] <- tprobs[mapping[k]] + alpha_filt[k]
      }
      template_prob_mat[t, ] <- tprobs
      dominant_reg[t] <- which.max(tprobs)
    }
    
    # ── Conditional moments & MVO ──────────────────────────────────────────
    if (initialized) {
      t_probs <- template_prob_mat[t, ]
      if (any(is.na(t_probs))) t_probs <- rep(1 / g_tmpl, g_tmpl)
      
      cond_mom <- compute_conditional_moments(templates, t_probs, n_assets)
      cond_mom$Sigma <- ledoit_wolf_shrink(cond_mom$Sigma, 252)
      
      w_new <- solve_mvo(cond_mom$mu, cond_mom$Sigma, w_prev,
                          gamma_ra, tau_tc, w_max)
    } else {
      w_new <- rep(1 / n_assets, n_assets)
    }
    
    # ── Record ─────────────────────────────────────────────────────────────
    weights_mat[t, ] <- w_new
    port_returns[t] <- sum(w_new * returns[t, ])
    turnover[t] <- 0.5 * sum(abs(w_new - w_prev))
    current_K[t] <- current_Kt
    w_prev <- w_new
    
    if (verbose && (t - t0) %% 500 == 0) {
      elapsed <- (proc.time() - t_start)[3]
      pct <- (t - t0) / (n_days - t0) * 100
      cat(sprintf("  [%s] %.1f%% done | K=%d | dom_regime=%s | elapsed=%.0fs\n",
                  dates[t], pct, current_Kt,
                  ifelse(is.na(dominant_reg[t]), "NA", dominant_reg[t]),
                  elapsed))
    }
  }
  
  elapsed_total <- (proc.time() - t_start)[3]
  if (verbose) cat(sprintf("[backtest] Completed in %.1f seconds.\n", elapsed_total))
  
  oos <- t0:n_days
  
  list(
    dates           = dates[oos],
    weights         = weights_mat[oos, , drop = FALSE],
    returns         = port_returns[oos],
    asset_returns   = returns[oos, , drop = FALSE],
    turnover        = turnover[oos],
    dominant_regime = dominant_reg[oos],
    template_probs  = template_prob_mat[oos, , drop = FALSE],
    current_K       = current_K[oos],
    w2_distances    = w2_distances[oos],
    templates       = templates,
    config = list(G = g_tmpl, eta = eta, gamma = gamma_ra, tau = tau_tc,
                  lambda_k = lambda_k, w_max = w_max)
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. BENCHMARK: EQUAL WEIGHT
# ═══════════════════════════════════════════════════════════════════════════════

run_equal_weight <- function(feat, oos_start = OOS_START_DATE,
                              rebal_freq = EW_REBAL_FREQ) {
  dates   <- feat$dates
  returns <- feat$returns
  n_days  <- length(dates)
  n_assets <- ncol(returns)
  
  oos_idx <- which(dates >= oos_start)
  t0 <- oos_idx[1]
  
  target_w <- rep(1 / n_assets, n_assets)
  w_prev <- target_w
  
  weights_mat  <- matrix(NA_real_, length(oos_idx), n_assets)
  port_returns <- rep(NA_real_, length(oos_idx))
  turnover_vec <- rep(NA_real_, length(oos_idx))
  
  days_since_rebal <- rebal_freq
  
  for (i in seq_along(oos_idx)) {
    t <- oos_idx[i]
    days_since_rebal <- days_since_rebal + 1
    
    if (days_since_rebal >= rebal_freq) {
      w_new <- target_w
      days_since_rebal <- 0
    } else {
      w_drifted <- w_prev * (1 + returns[t, ])
      w_new <- w_drifted / sum(w_drifted)
    }
    
    weights_mat[i, ] <- w_new
    port_returns[i] <- sum(w_new * returns[t, ])
    turnover_vec[i] <- 0.5 * sum(abs(w_new - w_prev))
    w_prev <- w_new
  }
  
  colnames(weights_mat) <- ASSET_NAMES
  
  list(
    dates    = dates[oos_idx],
    weights  = weights_mat,
    returns  = port_returns,
    asset_returns = returns[oos_idx, , drop = FALSE],
    turnover = turnover_vec,
    label    = "Equal Weight (20%)"
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. BENCHMARK: 60/40 STOCK/BOND
# ═══════════════════════════════════════════════════════════════════════════════

run_sixty_forty <- function(feat, oos_start = OOS_START_DATE,
                             rebal_freq = EW_REBAL_FREQ) {
  dates   <- feat$dates
  returns <- feat$returns
  n_days  <- length(dates)
  n_assets <- ncol(returns)
  
  oos_idx <- which(dates >= oos_start)
  t0 <- oos_idx[1]
  
  target_w <- BENCH_6040_ALLOC[ASSET_NAMES]
  w_prev <- target_w
  
  weights_mat  <- matrix(NA_real_, length(oos_idx), n_assets)
  port_returns <- rep(NA_real_, length(oos_idx))
  turnover_vec <- rep(NA_real_, length(oos_idx))
  
  days_since_rebal <- rebal_freq
  
  for (i in seq_along(oos_idx)) {
    t <- oos_idx[i]
    days_since_rebal <- days_since_rebal + 1
    
    if (days_since_rebal >= rebal_freq) {
      w_new <- target_w
      days_since_rebal <- 0
    } else {
      w_drifted <- w_prev * (1 + returns[t, ])
      w_new <- w_drifted / sum(w_drifted)
    }
    
    weights_mat[i, ] <- w_new
    port_returns[i] <- sum(w_new * returns[t, ])
    turnover_vec[i] <- 0.5 * sum(abs(w_new - w_prev))
    w_prev <- w_new
  }
  
  colnames(weights_mat) <- ASSET_NAMES
  
  list(
    dates    = dates[oos_idx],
    weights  = weights_mat,
    returns  = port_returns,
    asset_returns = returns[oos_idx, , drop = FALSE],
    turnover = turnover_vec,
    label    = "60/40 Stock/Bond"
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. BUY & HOLD SPX
# ═══════════════════════════════════════════════════════════════════════════════

run_spx_buyhold <- function(feat, oos_start = OOS_START_DATE) {
  dates   <- feat$dates
  returns <- feat$returns
  n_assets <- ncol(returns)
  
  oos_idx <- which(dates >= oos_start)
  spx_col <- which(ASSET_NAMES == "stocks")
  
  port_returns <- returns[oos_idx, spx_col]
  
  weights_mat <- matrix(0, length(oos_idx), n_assets)
  weights_mat[, spx_col] <- 1
  colnames(weights_mat) <- ASSET_NAMES
  
  list(
    dates    = dates[oos_idx],
    weights  = weights_mat,
    returns  = port_returns,
    asset_returns = returns[oos_idx, , drop = FALSE],
    turnover = rep(0, length(oos_idx)),
    label    = "SPX Buy & Hold"
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. PERFORMANCE METRICS
# ═══════════════════════════════════════════════════════════════════════════════

compute_metrics <- function(result, annualize = 252) {
  r <- result$returns
  r <- r[!is.na(r)]
  
  n <- length(r)
  cum_ret <- cumsum(r)
  
  ann_ret <- mean(r) * annualize
  ann_vol <- sd(r) * sqrt(annualize)
  sharpe  <- ann_ret / ann_vol
  
  down_r <- r[r < 0]
  down_vol <- sqrt(mean(down_r^2)) * sqrt(annualize)
  sortino  <- ann_ret / down_vol
  
  running_max <- cummax(cum_ret)
  drawdowns <- cum_ret - running_max
  max_dd <- min(drawdowns)
  
  hit_rate <- mean(r > 0)
  calmar <- ann_ret / abs(max_dd)
  avg_turnover <- mean(result$turnover, na.rm = TRUE)
  
  data.frame(
    Ann_Return   = ann_ret,
    Ann_Vol      = ann_vol,
    Sharpe       = sharpe,
    Sortino      = sortino,
    Max_DD       = max_dd,
    Hit_Rate     = hit_rate,
    Calmar       = calmar,
    Avg_Turnover = avg_turnover,
    N_Days       = n,
    stringsAsFactors = FALSE
  )
}

compute_regime_metrics <- function(result, annualize = 252) {
  if (is.null(result$dominant_regime)) return(NULL)
  
  r <- result$returns
  reg <- result$dominant_regime
  valid <- !is.na(r) & !is.na(reg)
  r <- r[valid]
  reg <- reg[valid]
  asset_r <- result$asset_returns[valid, , drop = FALSE]
  
  regimes <- sort(unique(reg))
  
  regime_stats <- do.call(rbind, lapply(regimes, function(g) {
    idx <- which(reg == g)
    rg <- r[idx]
    data.frame(
      Regime    = g,
      Days      = length(idx),
      Ann_Mean  = mean(rg) * annualize,
      Ann_Vol   = sd(rg) * sqrt(annualize),
      Sharpe    = (mean(rg) * annualize) / (sd(rg) * sqrt(annualize)),
      Hit_Rate  = mean(rg > 0),
      Max_DD    = min(cumsum(rg) - cummax(cumsum(rg)))
    )
  }))
  
  asset_sharpe <- do.call(rbind, lapply(regimes, function(g) {
    idx <- which(reg == g)
    sapply(seq_len(ncol(asset_r)), function(j) {
      ar <- asset_r[idx, j]
      (mean(ar) * annualize) / (sd(ar) * sqrt(annualize))
    })
  }))
  colnames(asset_sharpe) <- ASSET_LABELS
  asset_sharpe <- data.frame(Regime = regimes, asset_sharpe)
  
  list(portfolio = regime_stats, asset_sharpe = asset_sharpe)
}


cat("[backtest_engine.R] Functions loaded.\n")
