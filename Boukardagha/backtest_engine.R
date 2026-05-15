###############################################################################
# backtest_engine.R — Matching Python replication.py
###############################################################################

run_wasserstein_hmm <- function(feat, oos_start = OOS_START_DATE, cfg = list()) {
  
  g_tmpl    <- if (!is.null(cfg$G))     cfg$G     else G_TEMPLATES
  eta       <- if (!is.null(cfg$eta))   cfg$eta   else ETA_SMOOTH
  gamma_ra  <- if (!is.null(cfg$gamma)) cfg$gamma else GAMMA_RISK
  tau_tc    <- if (!is.null(cfg$tau))   cfg$tau   else TAU_TCOST
  lambda_k  <- if (!is.null(cfg$lambda_k)) cfg$lambda_k else LAMBDA_K
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
  
  if (verbose) cat(sprintf("[backtest] OOS: index %d (%s), %d days\n",
                           t0, dates[t0], n_days - t0 + 1))
  
  weights_mat       <- matrix(NA_real_, n_days, n_assets); colnames(weights_mat) <- ASSET_NAMES
  port_returns      <- rep(NA_real_, n_days)
  dominant_reg      <- rep(NA_integer_, n_days)
  template_prob_mat <- matrix(NA_real_, n_days, g_tmpl)
  current_K_vec     <- rep(NA_integer_, n_days)
  turnover_vec      <- rep(NA_real_, n_days)
  w2_distances      <- rep(NA_real_, n_days)
  
  templates      <- NULL
  last_model     <- NULL
  current_K      <- K_MIN
  w_prev         <- rep(1 / n_assets, n_assets)
  initialized    <- FALSE
  steps_since_refit     <- 0L
  steps_since_order_sel <- 0L
  
  # Initialize templates
  init_end <- t0 - 1
  init_start <- max(which(valid)[1], init_end - HMM_INIT_WINDOW + 1)
  init_idx <- which(valid[init_start:init_end]) + init_start - 1
  if (length(init_idx) < 100) stop("Insufficient data for template initialization")
  X_init <- features[init_idx, , drop = FALSE]
  
  if (verbose) cat(sprintf("  Initializing templates on %d rows...\n", nrow(X_init)))
  tmpl_result <- init_templates(X_init, g_tmpl)
  templates <- tmpl_result$templates
  last_model <- tmpl_result$model
  current_K <- g_tmpl
  initialized <- TRUE
  if (verbose) cat(sprintf("  Templates initialized (K=%d)\n", last_model$K))
  
  t_start <- proc.time()
  
  for (i in seq_along(oos_idx)) {
    t <- oos_idx[i]
    
    if (!valid[t]) {
      weights_mat[t, ] <- w_prev
      port_returns[t] <- sum(w_prev * returns[t, ])
      next
    }
    
    hist_idx <- which(valid[1:t])
    if (length(hist_idx) < 50) {
      weights_mat[t, ] <- w_prev
      port_returns[t] <- sum(w_prev * returns[t, ])
      next
    }
    X_hist <- features[hist_idx, , drop = FALSE]
    
    # Periodic K selection
    if (steps_since_order_sel == 0 || i == 1) {
      order_result <- select_model_order(X_hist, VAL_SLICE_DAYS,
                                          K_MIN, K_MAX, lambda_k)
      if (length(order_result$K) == 1) current_K <- order_result$K
    }
    steps_since_order_sel <- (steps_since_order_sel + 1L) %% ORDER_SELECT_FREQ
    
    # HMM refit
    do_refit <- (steps_since_refit == 0) ||
                is.null(last_model) ||
                (last_model$K != current_K)
    
    if (do_refit) {
      new_fit <- fit_gaussian_hmm(X_hist, current_K)
      if (!is.null(new_fit)) last_model <- new_fit
    }
    steps_since_refit <- (steps_since_refit + 1L) %% HMM_FIT_FREQ
    
    model <- last_model
    if (is.null(model)) {
      weights_mat[t, ] <- w_prev
      port_returns[t] <- sum(w_prev * returns[t, ])
      next
    }
    
    # Filtered posterior
    posterior <- model$filtered[nrow(model$filtered), ]
    if (any(!is.finite(posterior))) posterior <- rep(1 / model$K, model$K)
    
    # Assign to templates + track W2
    g_of_k <- assign_to_templates(model$mu, model$Sigma, templates)
    min_w2 <- Inf
    for (k in seq_len(model$K)) {
      d_w2 <- wasserstein2_gaussian(model$mu[[k]], model$Sigma[[k]],
                                     templates[[g_of_k[k]]]$mu, templates[[g_of_k[k]]]$Sigma)
      if (d_w2 < min_w2) min_w2 <- d_w2
    }
    w2_distances[t] <- min_w2
    
    # Template probabilities
    p_template <- rep(0, g_tmpl)
    for (k in seq_len(model$K)) {
      p_template[g_of_k[k]] <- p_template[g_of_k[k]] + posterior[k]
    }
    p_sum <- sum(p_template)
    if (p_sum > 1e-12) p_template <- p_template / p_sum
    
    template_prob_mat[t, ] <- p_template
    dominant_reg[t] <- which.max(p_template)
    current_K_vec[t] <- current_K
    
    # Update templates
    templates <- update_templates(templates, model$mu, model$Sigma,
                                  posterior, g_of_k, eta)
    
    # Conditional moments
    cond_mom <- compute_conditional_moments(templates, p_template, n_assets)
    
    # Shrinkage on raw return sample
    R_hist_idx <- which(valid[1:(t - 1)])
    if (length(R_hist_idx) > 50) {
      R_sample <- returns[R_hist_idx, , drop = FALSE]
      cond_mom$Sigma <- shrink_cov(cond_mom$Sigma, R_sample)
    }
    
    # MVO
    w_new <- solve_mvo(cond_mom$mu, cond_mom$Sigma, w_prev,
                        gamma_ra, tau_tc, w_max)
    
    weights_mat[t, ] <- w_new
    port_returns[t] <- sum(w_new * returns[t, ])
    turnover_vec[t] <- 0.5 * sum(abs(w_new - w_prev))
    w_prev <- w_new
    
    if (verbose && (i - 1) %% 500 == 0) {
      elapsed <- (proc.time() - t_start)[3]
      pct <- i / length(oos_idx) * 100
      cat(sprintf("  [%s] %.1f%% | K=%d | dom=%d | elapsed=%.0fs\n",
                  dates[t], pct, current_K,
                  ifelse(is.na(dominant_reg[t]), 0, dominant_reg[t]), elapsed))
    }
  }
  
  elapsed_total <- (proc.time() - t_start)[3]
  if (verbose) cat(sprintf("[backtest] Done in %.1f seconds.\n", elapsed_total))
  
  oos <- t0:n_days
  list(
    dates           = dates[oos],
    weights         = weights_mat[oos, , drop = FALSE],
    returns         = port_returns[oos],
    asset_returns   = returns[oos, , drop = FALSE],
    turnover        = turnover_vec[oos],
    dominant_regime = dominant_reg[oos],
    template_probs  = template_prob_mat[oos, , drop = FALSE],
    current_K       = current_K_vec[oos],
    w2_distances    = w2_distances[oos],
    templates       = templates,
    config = list(G = g_tmpl, eta = eta, gamma = gamma_ra, tau = tau_tc,
                  lambda_k = lambda_k, w_max = w_max)
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# BENCHMARKS + METRICS
# ═══════════════════════════════════════════════════════════════════════════════

run_sixty_forty <- function(feat, oos_start = OOS_START_DATE, rebal_freq = EW_REBAL_FREQ) {
  dates <- feat$dates; returns <- feat$returns; n_assets <- ncol(returns)
  oos_idx <- which(dates >= oos_start)
  target_w <- BENCH_6040_ALLOC[ASSET_NAMES]; w_prev <- target_w
  weights_mat <- matrix(NA_real_, length(oos_idx), n_assets)
  port_returns <- turnover_vec <- rep(NA_real_, length(oos_idx))
  dsr <- rebal_freq
  for (i in seq_along(oos_idx)) {
    t <- oos_idx[i]; dsr <- dsr + 1
    if (dsr >= rebal_freq) { w_new <- target_w; dsr <- 0
    } else { w_d <- w_prev * (1 + returns[t, ]); w_new <- w_d / sum(w_d) }
    weights_mat[i, ] <- w_new; port_returns[i] <- sum(w_new * returns[t, ])
    turnover_vec[i] <- 0.5 * sum(abs(w_new - w_prev)); w_prev <- w_new
  }
  colnames(weights_mat) <- ASSET_NAMES
  list(dates = dates[oos_idx], weights = weights_mat, returns = port_returns,
       asset_returns = returns[oos_idx, , drop = FALSE], turnover = turnover_vec)
}

run_spx_buyhold <- function(feat, oos_start = OOS_START_DATE) {
  dates <- feat$dates; returns <- feat$returns; n_assets <- ncol(returns)
  oos_idx <- which(dates >= oos_start); spx_col <- which(ASSET_NAMES == "stocks")
  weights_mat <- matrix(0, length(oos_idx), n_assets)
  weights_mat[, spx_col] <- 1; colnames(weights_mat) <- ASSET_NAMES
  list(dates = dates[oos_idx], weights = weights_mat, returns = returns[oos_idx, spx_col],
       asset_returns = returns[oos_idx, , drop = FALSE], turnover = rep(0, length(oos_idx)))
}

compute_metrics <- function(result, annualize = 252) {
  r <- result$returns[!is.na(result$returns)]; n <- length(r); cum_ret <- cumsum(r)
  ann_ret <- mean(r)*annualize; ann_vol <- sd(r)*sqrt(annualize)
  sharpe <- ann_ret/ann_vol
  down_r <- r[r<0]; down_vol <- sqrt(mean(down_r^2))*sqrt(annualize)
  sortino <- ann_ret/down_vol; max_dd <- min(cum_ret - cummax(cum_ret))
  data.frame(Ann_Return=ann_ret, Ann_Vol=ann_vol, Sharpe=sharpe, Sortino=sortino,
             Max_DD=max_dd, Hit_Rate=mean(r>0), Calmar=ann_ret/abs(max_dd),
             Avg_Turnover=mean(result$turnover,na.rm=TRUE), N_Days=n, stringsAsFactors=FALSE)
}

compute_regime_metrics <- function(result, annualize = 252) {
  if (is.null(result$dominant_regime)) return(NULL)
  r <- result$returns; reg <- result$dominant_regime
  v <- !is.na(r) & !is.na(reg); r <- r[v]; reg <- reg[v]
  asset_r <- result$asset_returns[v, , drop = FALSE]; regimes <- sort(unique(reg))
  regime_stats <- do.call(rbind, lapply(regimes, function(g) {
    rg <- r[reg==g]
    data.frame(Regime=g, Days=length(rg), Ann_Mean=mean(rg)*annualize,
               Ann_Vol=sd(rg)*sqrt(annualize),
               Sharpe=(mean(rg)*annualize)/(sd(rg)*sqrt(annualize)),
               Hit_Rate=mean(rg>0), Max_DD=min(cumsum(rg)-cummax(cumsum(rg))))
  }))
  asset_sharpe <- do.call(rbind, lapply(regimes, function(g) {
    idx <- which(reg==g)
    sapply(seq_len(ncol(asset_r)), function(j) {
      ar <- asset_r[idx,j]; (mean(ar)*annualize)/(sd(ar)*sqrt(annualize))
    })
  }))
  colnames(asset_sharpe) <- ASSET_LABELS
  asset_sharpe <- data.frame(Regime=regimes, asset_sharpe)
  list(portfolio=regime_stats, asset_sharpe=asset_sharpe)
}

cat("[backtest_engine.R] Loaded.\n")
