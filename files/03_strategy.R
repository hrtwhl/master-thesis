# =============================================================================
#  03_strategy.R
#  - Template-based regime tracking helpers
#  - Main Wasserstein-HMM + MVO backtest (strictly causal expanding window)
#  - Benchmarks: equal-weight, 60/40 stocks/bonds, SPX buy-and-hold
# =============================================================================

# ---------------- Regime moments (forward returns r_{s+1} | z_s = k) ----------

#' Conditional moments per regime k using forward returns r_{s+1} grouped by z_s.
#' This is what the author does in the notebook (and what makes economic sense
#' for one-step-ahead forecasting) although the paper text is silent on it.
compute_regime_moments_forward <- function(ret_align, z_s, K, N, min_obs = 60) {
  # z_s is integer regime label at time s (0-indexed length T-1)
  R <- as.matrix(ret_align)
  r_fwd <- R[-1, , drop = FALSE]                       # r_{s+1}
  if (length(z_s) > nrow(r_fwd)) z_s <- z_s[seq_len(nrow(r_fwd))]
  if (length(z_s) < nrow(r_fwd)) r_fwd <- r_fwd[seq_along(z_s), , drop = FALSE]
  
  MU  <- matrix(0, K, N)
  SIG <- array(0, dim = c(N, N, K))
  for (k in seq_len(K)) {
    idx <- which(z_s == (k - 1L))                       # 0-indexed in cpp output
    if (length(idx) < min_obs) {
      MU[k, ]    <- 0
      SIG[,,k]   <- diag(N)
    } else {
      Rk <- r_fwd[idx, , drop = FALSE]
      MU[k, ]    <- colMeans(Rk)
      SIG[,,k]   <- shrunk_cov(Rk)
    }
  }
  list(MU = MU, SIG = SIG)
}

# ---------------- Template tracking ----------------

#' Map each HMM component k to the nearest template g by Wasserstein-2 distance.
#' Returns 1-indexed integer vector of length K (template assignment for each k).
map_components_to_templates <- function(MU_K, SIG_K, MU_tpl, SIG_tpl) {
  K <- nrow(MU_K); G <- nrow(MU_tpl)
  map <- integer(K); dist <- numeric(K)
  for (k in seq_len(K)) {
    d <- numeric(G)
    for (g in seq_len(G)) {
      d[g] <- w2_gaussian(MU_tpl[g, ], SIG_tpl[,,g],
                          MU_K[k, ],   SIG_K[,,k])
    }
    g_best <- which.min(d)
    map[k]  <- g_best
    dist[k] <- d[g_best]
  }
  list(map = map, dist = dist)
}

#' Spawn a new template if a component is far from all templates and we have room.
spawn_template_if_needed <- function(MU_K, SIG_K, pK, MU_tpl, SIG_tpl,
                                     spawn_thresh, G_max) {
  if (is.null(MU_tpl)) return(list(MU = MU_K, SIG = SIG_K))
  G <- nrow(MU_tpl)
  if (G >= G_max) return(list(MU = MU_tpl, SIG = SIG_tpl))
  
  mp <- map_components_to_templates(MU_K, SIG_K, MU_tpl, SIG_tpl)
  cand <- which(mp$dist > spawn_thresh)
  if (length(cand) == 0) return(list(MU = MU_tpl, SIG = SIG_tpl))
  
  k_new <- cand[which.max(pK[cand])]
  MU_new  <- rbind(MU_tpl, MU_K[k_new, , drop = FALSE])
  SIG_new <- abind::abind(SIG_tpl, SIG_K[,,k_new, drop = FALSE], along = 3)
  list(MU = MU_new, SIG = SIG_new)
}

#' EMA update of templates using posterior-weighted component averages.
update_templates_ema <- function(MU_tpl, SIG_tpl, MU_K, SIG_K, pK, map_k_to_g,
                                 eta = 0.05) {
  G <- nrow(MU_tpl); K <- nrow(MU_K); N <- ncol(MU_tpl)
  for (g in seq_len(G)) {
    ks <- which(map_k_to_g == g)
    if (length(ks) == 0) next
    w <- pK[ks]
    ws <- sum(w)
    if (ws <= 1e-12) next
    w <- w / ws
    
    mu_bar <- as.numeric(crossprod(w, MU_K[ks, , drop = FALSE]))
    S_bar  <- matrix(0, N, N)
    for (i in seq_along(ks)) S_bar <- S_bar + w[i] * SIG_K[,,ks[i]]
    
    MU_tpl[g, ]  <- (1 - eta) * MU_tpl[g, ] + eta * mu_bar
    Snew         <- (1 - eta) * SIG_tpl[,,g] + eta * S_bar
    SIG_tpl[,,g] <- 0.5 * (Snew + t(Snew)) + 1e-8 * diag(N)
  }
  list(MU = MU_tpl, SIG = SIG_tpl)
}

#' Aggregate component posteriors p_K into template posteriors p_G.
aggregate_template_posteriors <- function(pK, map_k_to_g, G) {
  pG <- numeric(G)
  for (k in seq_along(pK)) pG[map_k_to_g[k]] <- pG[map_k_to_g[k]] + pK[k]
  s <- sum(pG)
  if (s > 1e-12) pG / s else rep(1 / G, G)
}

# ---------------- Main backtest ----------------

#' Strictly causal expanding-window Wasserstein HMM + MVO backtest.
#'
#' On each OOS day d:
#'   1. Build features on history ending at d-1.
#'   2. Every F_K days, run predictive K selection on the current history.
#'   3. Fit a K-state Gaussian HMM, get smoothed posteriors gamma_t(k),
#'      compute forward-return regime moments.
#'   4. Map each component to a persistent template via W2 distance, spawn
#'      new templates if necessary, EMA-update templates.
#'   5. Compute template-mixture moments and solve the L1-MVO QP.
#'   6. Realize P&L on day d as w_d . r_d (weights set at d-1 close).
run_strategy <- function(returns_train, returns_test,
                         features_all, asset_names, cfg, verbose = TRUE) {
  
  N <- length(asset_names)
  
  # Cache the full features and returns indexed by date for O(1) lookups
  feat_dates <- features_all$date
  feat_mat   <- as.matrix(features_all[, setdiff(names(features_all), "date")])
  
  all_ret <- rbind(returns_train, returns_test)
  all_ret <- all_ret[!duplicated(all_ret$date), ]
  all_ret <- all_ret[order(all_ret$date), ]
  ret_dates <- all_ret$date
  ret_mat   <- as.matrix(all_ret[, asset_names])
  
  test_dates <- returns_test$date
  Ntest <- length(test_dates)
  
  weights_oos <- matrix(NA_real_, Ntest, N)
  pnl_oos     <- numeric(Ntest)
  K_history   <- integer(Ntest)
  tpl_label   <- integer(Ntest); tpl_label[] <- NA_integer_
  tpl_prob    <- numeric(Ntest);  tpl_prob[]  <- NA_real_
  tpl_count   <- integer(Ntest)
  
  w_prev <- rep(0, N)
  K_curr <- cfg$K_min
  MU_tpl <- NULL; SIG_tpl <- NULL
  
  pb_step <- max(1L, floor(Ntest / 50L))
  
  for (i in seq_len(Ntest)) {
    d <- test_dates[i]
    
    # ---- history of features available at d (use up to d-1)
    feat_idx_end <- findInterval(d - 1L, feat_dates)
    if (feat_idx_end < cfg$lookback) {
      weights_oos[i, ] <- w_prev
      pnl_oos[i] <- as.numeric(w_prev %*% ret_mat[ret_dates == d, ])
      K_history[i] <- NA_integer_
      next
    }
    X <- feat_mat[seq_len(feat_idx_end), , drop = FALSE]
    feat_hist_dates <- feat_dates[seq_len(feat_idx_end)]
    
    # Align return history to feature dates (for forward-return regime moments)
    ret_align <- ret_mat[match(feat_hist_dates, ret_dates), , drop = FALSE]
    
    # ---- (A) Predictive K selection (every F_K days)
    if (((i - 1L) %% cfg$F_K) == 0L) {
      K_curr <- select_K_predictive(
        X,
        K_candidates = seq.int(cfg$K_min, cfg$K_max),
        L_val   = cfg$L_val,
        n_iter  = cfg$em_iter,
        lambda_K = cfg$lambda_K,
        seed     = cfg$seed
      )
    }
    
    # ---- (B) Fit HMM at K_curr on full history up to d-1
    fit <- try(fit_hmm(X, K = K_curr, n_iter = cfg$em_iter, seed = cfg$seed),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      weights_oos[i, ] <- w_prev
      pnl_oos[i] <- as.numeric(w_prev %*% ret_mat[ret_dates == d, ])
      K_history[i] <- K_curr
      next
    }
    probs <- predict_proba_hmm(X, fit)
    pK    <- probs[nrow(probs), ]
    zK    <- max.col(probs, ties.method = "first") - 1L  # 0-indexed
    
    # Regime moments using forward returns r_{s+1} | z_s = k
    rm <- compute_regime_moments_forward(ret_align, z_s = zK[-length(zK)],
                                         K = K_curr, N = N, min_obs = 60)
    MU_K <- rm$MU; SIG_K <- rm$SIG
    
    # ---- (C) Templates: init / spawn / map / EMA update
    if (is.null(MU_tpl)) {
      MU_tpl  <- MU_K
      SIG_tpl <- SIG_K
    } else {
      sp <- spawn_template_if_needed(MU_K, SIG_K, pK,
                                     MU_tpl, SIG_tpl,
                                     spawn_thresh = cfg$spawn_thresh,
                                     G_max = cfg$G_max)
      MU_tpl <- sp$MU; SIG_tpl <- sp$SIG
    }
    
    mp <- map_components_to_templates(MU_K, SIG_K, MU_tpl, SIG_tpl)
    G  <- nrow(MU_tpl)
    pG <- aggregate_template_posteriors(pK, mp$map, G)
    
    upd <- update_templates_ema(MU_tpl, SIG_tpl, MU_K, SIG_K, pK, mp$map,
                                eta = cfg$eta_tpl)
    MU_tpl <- upd$MU; SIG_tpl <- upd$SIG
    
    # ---- (D) Template mixture moments + MVO
    mu_t    <- as.numeric(crossprod(pG, MU_tpl))                       # N-vec
    Sigma_t <- matrix(0, N, N)
    for (g in seq_len(G)) Sigma_t <- Sigma_t + pG[g] * SIG_tpl[,,g]
    
    w_t <- solve_mvo(mu_t, Sigma_t, w_prev,
                     lambda = cfg$lambda, tc = cfg$tc, w_max = cfg$w_max)
    
    weights_oos[i, ] <- w_t
    pnl_oos[i]       <- as.numeric(w_t %*% ret_mat[ret_dates == d, ])
    K_history[i]     <- K_curr
    tpl_label[i]     <- which.max(pG)                                  # 1-indexed
    tpl_prob[i]      <- max(pG)
    tpl_count[i]     <- G
    
    w_prev <- w_t
    
    if (verbose && (i %% pb_step == 0L)) {
      cat(sprintf("  [%4d / %4d]  %s  K=%d  G=%d\n",
                  i, Ntest, format(d), K_curr, G))
    }
  }
  
  colnames(weights_oos) <- asset_names
  list(
    weights    = data.frame(date = test_dates, weights_oos, check.names = FALSE),
    pnl        = data.frame(date = test_dates, pnl = pnl_oos),
    K_history  = data.frame(date = test_dates, K = K_history),
    tpl_label  = data.frame(date = test_dates, label = tpl_label, prob = tpl_prob,
                            G = tpl_count)
  )
}

# ---------------- Benchmarks ----------------

#' Equal-weight portfolio (daily-rebalanced), aligned to the same dates as the
#' strategy. Daily log return = w' r_t with w fixed at 1/N.
benchmark_equal_weight <- function(returns_df, asset_names, dates) {
  R <- as.matrix(returns_df[match(dates, returns_df$date), asset_names])
  N <- length(asset_names)
  w <- rep(1 / N, N)
  data.frame(date = dates, pnl = as.numeric(R %*% w))
}

#' 60/40 stocks/bonds (daily-rebalanced).
benchmark_60_40 <- function(returns_df, dates, stocks_col = "SPX", bonds_col = "BOND") {
  R <- returns_df[match(dates, returns_df$date), c(stocks_col, bonds_col)]
  data.frame(date = dates, pnl = 0.6 * R[[stocks_col]] + 0.4 * R[[bonds_col]])
}

#' SPX buy-and-hold (no rebalancing needed since single asset).
benchmark_spx <- function(returns_df, dates, stocks_col = "SPX") {
  data.frame(date = dates,
             pnl = returns_df[match(dates, returns_df$date), stocks_col])
}
