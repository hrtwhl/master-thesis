# 04_backtest.R -- walk-forward backtest engine and performance metrics.
#
# At each month t in the OOS window:
#   1. Use macro up to t-1 to classify current regime distribution.
#   2. Apply transition matrix to get next-period regime distribution.
#   3. Refit forecast models on returns 1:(t-1) conditioned on labels 1:(t-1).
#   4. Generate forecasts for t+1; apply sizing to get portfolio weights.
#   5. Realise r_t = w . r[t].
#
# Optional walk-forward regime refit every `refit_regime_every` months
# (Hungarian label matching).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---- performance metrics ----------------------------------------------------

sharpe_annual <- function(r, rf = 0, freq = 12) {
  r <- na.omit(r - rf)
  if (length(r) < 2) return(NA_real_)
  (mean(r) / stats::sd(r)) * sqrt(freq)
}

sortino_annual <- function(r, rf = 0, freq = 12) {
  r <- na.omit(r - rf)
  if (length(r) < 2) return(NA_real_)
  downside <- r[r < 0]
  if (length(downside) < 2) return(NA_real_)
  (mean(r) / sqrt(mean(downside^2))) * sqrt(freq)
}

max_drawdown <- function(r) {
  r <- na.omit(r)
  if (!length(r)) return(NA_real_)
  equity <- cumprod(1 + r); peaks <- cummax(equity)
  min(equity / peaks - 1) * 100
}

avg_drawdown <- function(r) {
  r <- na.omit(r)
  if (!length(r)) return(NA_real_)
  equity <- cumprod(1 + r); peaks <- cummax(equity)
  dd <- equity / peaks - 1
  mean(dd[dd < 0]) * 100
}

pct_positive <- function(r) mean(na.omit(r) > 0)

performance_row <- function(name, r) {
  tibble(
    Model         = name,
    Sharpe        = sharpe_annual(r),
    Sortino       = sortino_annual(r),
    AvgDD         = avg_drawdown(r),
    MaxDD         = max_drawdown(r),
    PctPositive   = pct_positive(r),
    AnnReturn     = mean(na.omit(r)) * 12,
    AnnVol        = stats::sd(na.omit(r)) * sqrt(12)
  )
}

compute_performance_table <- function(returns_df) {
  strat_cols <- setdiff(colnames(returns_df), "date")
  dplyr::bind_rows(lapply(strat_cols, function(nm) {
    performance_row(nm, returns_df[[nm]])
  }))
}

# ---- walk-forward backtest --------------------------------------------------

run_backtest <- function(macro_raw, factor_rets,
                         backtest_start = as.Date("2001-01-01"),
                         train_min = 24L,
                         k_inner = 5L,
                         pca_var = 0.95,
                         l_grid = c(2, 3, 4),
                         vol_target = 0.10,
                         refit_regime_every = NULL,
                         verbose = TRUE) {

  all_dates <- sort(as.Date(intersect(macro_raw$date, factor_rets$date),
                            origin = "1970-01-01"))
  macro_raw   <- macro_raw   |> filter(date %in% all_dates) |> arrange(date)
  factor_rets <- factor_rets |> filter(date %in% all_dates) |> arrange(date)

  oos_idx <- which(all_dates >= backtest_start & all_dates >= all_dates[train_min])
  if (length(oos_idx) == 0) stop("Backtest window empty after alignment.")

  # initial regime model fit on all pre-OOS macro data
  fit_end_ix <- oos_idx[1] - 1L
  macro_fit <- macro_raw[seq_len(fit_end_ix), ]
  pca_fit   <- fit_pca(macro_fit, var_thresh = pca_var)
  X_fit     <- pca_fit$scores
  regime_model <- fit_regime_model(X_fit, k_inner = k_inner)

  hist_class <- classify_regimes(regime_model, X_fit)
  labels_hist <- hist_class$labels
  n_reg <- k_inner + 1L
  E_hist <- regime_transition_matrix(labels_hist, n_reg)

  factor_cols <- setdiff(colnames(factor_rets), "date")

  # ---- strategy space ----
  base_with_l <- c(
    paste0("naive_", c("lo", "lns", "los", "mx")),
    paste0("ridge_", c("lo", "lns", "los", "mx")),
    paste0("bl_",    c("lo", "lns"))
  )
  strategies <- c(
    as.vector(outer(base_with_l, l_grid, function(x, l) paste(x, l, sep = "_"))),
    "mvo", "ew"
  )
  has_spy <- "SPY" %in% factor_cols
  if (has_spy) strategies <- c(strategies, "spy")

  regime_log    <- vector("list", length(oos_idx))
  portfolio_ret <- tibble(date = all_dates[oos_idx])
  for (s in strategies) portfolio_ret[[s]] <- NA_real_

  # ---- main walk-forward loop ----
  for (k in seq_along(oos_idx)) {
    t <- oos_idx[k]
    if (verbose && k %% 24 == 0) {
      message(sprintf("  [%d/%d] %s", k, length(oos_idx), all_dates[t]))
    }

    # optional walk-forward regime refit
    if (!is.null(refit_regime_every) && k > 1 &&
        (k - 1) %% refit_regime_every == 0) {
      macro_refit <- macro_raw[seq_len(t - 1), ]
      pca_fit <- fit_pca(macro_refit, var_thresh = pca_var)
      new_model <- fit_regime_model(pca_fit$scores, k_inner = k_inner)
      new_model <- match_regime_labels(new_model, regime_model)
      regime_model <- new_model
      hist_class <- classify_regimes(regime_model, pca_fit$scores)
      labels_hist <- hist_class$labels
      E_hist <- regime_transition_matrix(labels_hist, n_reg)
    }

    # current regime distribution (use macro at t-1 to stay out-of-sample)
    x_now <- project_pca(pca_fit, macro_raw[t - 1, ])
    class_now <- classify_regimes(regime_model, x_now)
    p_now  <- class_now$probs[1, ]
    p_next <- next_regime_probs(p_now, E_hist)

    past_ix <- seq_len(t - 1L)
    X_past <- project_pca(pca_fit, macro_raw[past_ix, ])
    labels_past <- classify_regimes(regime_model, X_past)$labels
    rets_past <- factor_rets[past_ix, ]

    # ---- forecasts ----
    y_naive <- forecast_naive(rets_past, labels_past, p_next)
    y_bl    <- regime_conditional_mean(rets_past, labels_past, p_next)
    y_ridge <- forecast_ridge(rets_past, X_past, labels_past, p_next,
                              lag_features = TRUE)

    # realised return at t
    r_t <- as.numeric(factor_rets[t, factor_cols])
    names(r_t) <- factor_cols

    forecast_reg_is_zero <- (which.max(p_next) - 1L) == 0L

    apply_and_log <- function(w, strat_name) {
      w_ordered <- setNames(rep(0, length(factor_cols)), factor_cols)
      common_nm <- intersect(names(w), factor_cols)
      w_ordered[common_nm] <- w[common_nm]
      portfolio_ret[[strat_name]][k] <<- sum(w_ordered * r_t, na.rm = TRUE)
    }

    for (l in l_grid) {
      apply_and_log(sizing_long_only(y_naive, l),     paste0("naive_lo_",  l))
      apply_and_log(sizing_long_short(y_naive, l),    paste0("naive_lns_", l))
      apply_and_log(sizing_long_or_short(y_naive, l), paste0("naive_los_", l))
      apply_and_log(sizing_mixed(y_naive, l, forecast_reg_is_zero),
                    paste0("naive_mx_",  l))

      apply_and_log(sizing_long_only(y_ridge, l),     paste0("ridge_lo_",  l))
      apply_and_log(sizing_long_short(y_ridge, l),    paste0("ridge_lns_", l))
      apply_and_log(sizing_long_or_short(y_ridge, l), paste0("ridge_los_", l))
      apply_and_log(sizing_mixed(y_ridge, l, forecast_reg_is_zero),
                    paste0("ridge_mx_",  l))
    }

    # BL with LW-shrunk prior covariance
    prior_R <- as.matrix(rets_past[, factor_cols])
    mu_prior <- colMeans(prior_R)
    S_prior  <- ledoit_wolf_cov(prior_R)
    w_bl_full <- bl_weights(mu_prior, S_prior, q = y_bl)

    for (l in l_grid) {
      trim_top <- function(w, l, sign_allowed = c("long", "both")) {
        sign_allowed <- match.arg(sign_allowed)
        if (sign_allowed == "long") {
          ord <- order(w, decreasing = TRUE)
          keep <- ord[seq_len(min(l, sum(w > 0)))]
        } else {
          ord <- order(abs(w), decreasing = TRUE)
          keep <- ord[seq_len(min(l, length(w)))]
        }
        w2 <- numeric(length(w)); names(w2) <- names(w)
        w2[keep] <- w[keep]
        w2 / max(sum(abs(w2)), 1e-12)
      }
      apply_and_log(trim_top(w_bl_full, l, "long"), paste0("bl_lo_",  l))
      apply_and_log(trim_top(w_bl_full, l, "both"), paste0("bl_lns_", l))
    }

    # MVO benchmark: long-only, LW-shrunk, uses the full universe (no l)
    w_mvo <- mvo_weights(rets_past, long_only = TRUE)
    apply_and_log(w_mvo, "mvo")

    # EW and SPY benchmarks
    w_ew <- sizing_equal_weight(factor_cols)
    portfolio_ret[["ew"]][k] <- sum(w_ew * r_t)
    if (has_spy) portfolio_ret[["spy"]][k] <- r_t["SPY"]

    regime_log[[k]] <- dplyr::bind_cols(
      tibble(date = all_dates[t], current_label = class_now$labels[1]),
      as_tibble(rbind(p_now)),
      as_tibble(rbind(setNames(p_next, paste0("pnext_R", 0:(n_reg - 1)))))
    )
  }

  # apply vol target to all strategy columns
  portfolio_ret_vs <- portfolio_ret
  for (s in setdiff(colnames(portfolio_ret_vs), "date")) {
    portfolio_ret_vs[[s]] <- vol_scale(portfolio_ret_vs[[s]],
                                       vol_target = vol_target)
  }

  list(
    portfolio_returns        = portfolio_ret,
    portfolio_returns_scaled = portfolio_ret_vs,
    regime_log               = dplyr::bind_rows(regime_log),
    regime_model             = regime_model,
    pca_fit                  = pca_fit,
    transition_matrix        = E_hist,
    config = list(backtest_start = backtest_start, k_inner = k_inner,
                  vol_target = vol_target, l_grid = l_grid,
                  refit_regime_every = refit_regime_every)
  )
}

# ---- statistical tests (for Nemenyi / t-tests) -----------------------------

paired_t_test_one_sided <- function(control, treatment) {
  stats::t.test(treatment, control, paired = TRUE, alternative = "greater")
}

nemenyi_rank <- function(control, treatment) {
  mat <- cbind(ctrl = control, trt = treatment)
  ranks <- t(apply(mat, 1, rank))
  colMeans(ranks)
}
