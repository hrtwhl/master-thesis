# 03_models.R -- forecasting models and position sizing.
#
# Forecasts:
#   - Naive: Sharpe conditional on argmax next regime (eq. 9)
#   - Black-Litterman: LW-shrunk sample cov, regime-conditional views (eq. 19)
#   - Ridge: per-regime linear ridge on LAGGED macro features (x_{s-1} -> r_s),
#     probability-weighted aggregation across regimes (eq. 12-14)
#
# Position sizing (sec. 5.3):
#   lo  : long-only top-l                los : l largest |forecast|
#   lns : top-l long + bottom-l short    mx  : lo normally, lns in Regime 0
#
# MVO benchmark is long-only with Ledoit-Wolf covariance shrinkage (no l).

suppressPackageStartupMessages({
  library(dplyr)
})

# =============================================================================
# Covariance shrinkage (Ledoit-Wolf 2004 constant-correlation target)
# =============================================================================

#' Ledoit-Wolf shrinkage estimator of the covariance matrix with the
#' constant-correlation target (Ledoit & Wolf 2003, JEF).
#'
#' Pure-R implementation; no external package required. Degrades gracefully
#' to sample cov + small diagonal if shrinkage intensity cannot be computed.
ledoit_wolf_cov <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (n < 2L) return(diag(p))

  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- crossprod(Xc) / n      # MLE (Ledoit-Wolf uses this, not n-1)

  # target F: constant-correlation
  v   <- diag(S)
  sds <- sqrt(pmax(v, 1e-16))
  corr <- S / outer(sds, sds)
  rbar <- (sum(corr) - p) / max(p * (p - 1), 1)
  F <- rbar * outer(sds, sds)
  diag(F) <- v

  # pi: sum of asymptotic variances of S_ij
  # Var(sqrt(n) * S_ij) estimated by 1/n sum_t (Xc_ti Xc_tj - S_ij)^2
  Y <- Xc
  pi_mat <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) {
    pi_mat[i, j] <- mean((Y[, i] * Y[, j] - S[i, j])^2)
  }
  pi_sum <- sum(pi_mat)

  # rho: LW's "covariance between S and F" term
  # formula from LW(2003) appendix -- for constant-corr target
  theta_ii_ij <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) {
    if (i == j) next
    theta_ii_ij[i, j] <-
      mean((Y[, i]^2 - S[i, i]) * (Y[, i] * Y[, j] - S[i, j]))
  }
  rho_diag <- sum(diag(pi_mat))
  rho_off  <- 0
  for (i in seq_len(p)) for (j in seq_len(p)) {
    if (i == j) next
    rho_off <- rho_off +
      (rbar / 2) * (sds[j] / sds[i] * theta_ii_ij[i, j] +
                    sds[i] / sds[j] * theta_ii_ij[j, i])
  }
  rho_sum <- rho_diag + rho_off

  # gamma: Frobenius distance (F - S)
  gamma_sum <- sum((F - S)^2)

  # kappa = (pi - rho) / gamma; shrinkage = max(0, min(1, kappa / n))
  kappa <- (pi_sum - rho_sum) / max(gamma_sum, 1e-12)
  delta <- max(0, min(1, kappa / n))

  S_shrunk <- delta * F + (1 - delta) * S
  # ensure symmetric PSD
  S_shrunk <- 0.5 * (S_shrunk + t(S_shrunk))
  S_shrunk + diag(1e-10, p)
}

# =============================================================================
# Forecasting
# =============================================================================

#' Naive regime-conditional Sharpe forecast (eq. 9).
forecast_naive <- function(returns_past, regime_labels_past, next_regime_prob,
                           mode = c("argmax", "weighted")) {
  mode <- match.arg(mode)
  fac <- setdiff(colnames(returns_past), "date")
  R <- as.matrix(returns_past[, fac])
  regimes <- sort(unique(regime_labels_past))

  sr_by_regime <- t(sapply(regimes, function(r) {
    sel <- regime_labels_past == r
    if (sum(sel) < 3) return(rep(NA_real_, length(fac)))
    mu <- colMeans(R[sel, , drop = FALSE])
    sd <- apply(R[sel, , drop = FALSE], 2, stats::sd)
    sd[sd == 0] <- NA
    mu / sd
  }))
  rownames(sr_by_regime) <- as.character(regimes)
  colnames(sr_by_regime) <- fac

  overall <- colMeans(R) / apply(R, 2, stats::sd)
  for (j in seq_len(ncol(sr_by_regime))) {
    sr_by_regime[is.na(sr_by_regime[, j]), j] <- overall[j]
  }

  if (mode == "argmax") {
    i_star <- which.max(next_regime_prob)
    key <- as.character(regimes[min(i_star, length(regimes))])
    return(sr_by_regime[key, ])
  }
  p <- next_regime_prob[seq_len(nrow(sr_by_regime))]
  p <- p / sum(p)
  as.numeric(t(sr_by_regime) %*% p) |> setNames(fac)
}

#' Regime-conditional mean return (used as BL views).
regime_conditional_mean <- function(returns_past, regime_labels_past,
                                    next_regime_prob,
                                    mode = c("argmax", "weighted")) {
  mode <- match.arg(mode)
  fac <- setdiff(colnames(returns_past), "date")
  R <- as.matrix(returns_past[, fac])
  regimes <- sort(unique(regime_labels_past))

  mu_by_regime <- t(sapply(regimes, function(r) {
    sel <- regime_labels_past == r
    if (sum(sel) < 3) return(rep(NA_real_, length(fac)))
    colMeans(R[sel, , drop = FALSE])
  }))
  rownames(mu_by_regime) <- as.character(regimes)
  colnames(mu_by_regime) <- fac

  overall <- colMeans(R)
  for (j in seq_len(ncol(mu_by_regime))) {
    mu_by_regime[is.na(mu_by_regime[, j]), j] <- overall[j]
  }

  if (mode == "argmax") {
    i_star <- which.max(next_regime_prob)
    key <- as.character(regimes[min(i_star, length(regimes))])
    return(mu_by_regime[key, ])
  }
  p <- next_regime_prob[seq_len(nrow(mu_by_regime))]
  p <- p / sum(p)
  as.numeric(t(mu_by_regime) %*% p) |> setNames(fac)
}

#' Black-Litterman weights (eq. 19) with LW-shrunk prior covariance.
#'
#' Uses shrunk sample mean/cov as prior, pick matrix P = I. `q` is the
#' full-vector view of expected returns per asset (regime-conditional).
bl_weights <- function(prior_mean, prior_cov, q, tau = 0.05,
                       risk_aversion = 1.0, omega_scale = NULL) {
  if (is.null(omega_scale)) omega_scale <- tau
  n <- length(prior_mean)
  P <- diag(n)
  Omega <- diag(diag(P %*% prior_cov %*% t(P)) * omega_scale, n)

  tauSigma_inv <- solve(tau * prior_cov + diag(1e-8, n))
  Omega_inv    <- solve(Omega + diag(1e-8, n))

  post_prec <- tauSigma_inv + t(P) %*% Omega_inv %*% P
  post_mean <- solve(post_prec,
                     tauSigma_inv %*% prior_mean + t(P) %*% Omega_inv %*% q)
  post_cov  <- solve(post_prec) + prior_cov

  w <- solve(risk_aversion * post_cov, post_mean)
  w <- as.numeric(w) / max(sum(abs(w)), 1e-12)
  setNames(w, names(prior_mean))
}

#' Per-regime ridge on LAGGED macro features (predictive version of paper
#' eq. 12-14):
#'
#'   Fit beta_i on pairs (x_{s-1}, r_s) where regime at s == i.
#'   Forecast r_{t+1} via sum_i P(R_{t+1}=i) * beta_i * x_t.
#'
#' This avoids look-ahead (x_t is known at time t, used to predict r_{t+1}).
forecast_ridge <- function(returns_past, macro_past, regime_labels_past,
                           next_regime_prob, min_obs = 30L,
                           lambda = "cv", lag_features = TRUE) {
  if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
  fac <- setdiff(colnames(returns_past), "date")
  R <- as.matrix(returns_past[, fac])
  X <- as.matrix(macro_past)
  stopifnot(nrow(R) == nrow(X))
  n <- nrow(X)
  if (n < 4L) return(setNames(rep(0, length(fac)), fac))

  if (lag_features) {
    # fit pairs (x_{s-1}, r_s); label conditions on regime(s)
    X_fit <- X[-n, , drop = FALSE]
    R_fit <- R[-1, , drop = FALSE]
    labels_fit <- regime_labels_past[-1]
    x_pred <- matrix(X[n, ], nrow = 1)   # known at time t; forecast r_{t+1}
  } else {
    X_fit <- X; R_fit <- R; labels_fit <- regime_labels_past
    x_pred <- matrix(X[n, ], nrow = 1)
  }

  regimes <- sort(unique(labels_fit))
  preds <- matrix(NA_real_, nrow = length(regimes), ncol = length(fac),
                  dimnames = list(as.character(regimes), fac))

  for (r in regimes) {
    sel <- labels_fit == r
    if (sum(sel) < min_obs) {
      Xr <- X_fit; yr <- R_fit    # global fallback
    } else {
      Xr <- X_fit[sel, , drop = FALSE]; yr <- R_fit[sel, , drop = FALSE]
    }
    for (j in seq_along(fac)) {
      y <- yr[, j]
      if (length(unique(y)) < 2) { preds[as.character(r), j] <- mean(y); next }
      if (identical(lambda, "cv")) {
        fit <- try(glmnet::cv.glmnet(Xr, y, alpha = 0, nfolds = 5,
                                     standardize = FALSE), silent = TRUE)
        if (inherits(fit, "try-error")) {
          preds[as.character(r), j] <- mean(y); next
        }
        preds[as.character(r), j] <- as.numeric(
          predict(fit, newx = x_pred, s = "lambda.min"))
      } else {
        fit <- glmnet::glmnet(Xr, y, alpha = 0, lambda = lambda,
                              standardize = FALSE)
        preds[as.character(r), j] <- as.numeric(predict(fit, newx = x_pred))
      }
    }
  }

  p <- next_regime_prob[seq_len(nrow(preds))]
  p <- p / sum(p)
  as.numeric(t(preds) %*% p) |> setNames(fac)
}

# =============================================================================
# Position sizing (sec 5.3)
# =============================================================================

sizing_long_short <- function(forecast, l) {
  fac <- names(forecast)
  ord <- order(forecast, decreasing = TRUE)
  n <- length(forecast)
  l <- min(l, floor(n / 2))
  long_idx  <- ord[seq_len(l)]
  short_idx <- ord[(n - l + 1):n]
  w <- numeric(n); names(w) <- fac
  selected <- c(long_idx, short_idx)
  w[selected] <- forecast[selected]
  w / max(sum(abs(w[selected])), 1e-12)
}

sizing_long_or_short <- function(forecast, l) {
  fac <- names(forecast)
  ord <- order(abs(forecast), decreasing = TRUE)
  keep <- ord[seq_len(min(l, length(forecast)))]
  w <- numeric(length(forecast)); names(w) <- fac
  w[keep] <- forecast[keep]
  w / max(sum(abs(w)), 1e-12)
}

sizing_long_only <- function(forecast, l) {
  fac <- names(forecast)
  pos <- which(forecast > 0)
  if (length(pos) == 0) return(setNames(numeric(length(forecast)), fac))
  ord <- order(forecast, decreasing = TRUE)
  keep <- intersect(ord, pos)[seq_len(min(l, length(pos)))]
  w <- numeric(length(forecast)); names(w) <- fac
  w[keep] <- forecast[keep]
  w / max(sum(abs(w)), 1e-12)
}

sizing_mixed <- function(forecast, l, forecast_regime_is_zero) {
  if (isTRUE(forecast_regime_is_zero)) sizing_long_short(forecast, l)
  else sizing_long_only(forecast, l)
}

# ---- benchmark weights ------------------------------------------------------

sizing_equal_weight <- function(factor_names) {
  w <- rep(1 / length(factor_names), length(factor_names))
  setNames(w, factor_names)
}

#' Long-only mean-variance optimal weights with Ledoit-Wolf covariance
#' shrinkage. Solves:
#'   max_w  mu' w - (lambda/2) w' Sigma w   s.t.  sum(w) = 1, w >= 0
#' Falls back to equal-weight if quadprog fails.
mvo_weights <- function(returns_past, risk_aversion = 1.0, long_only = TRUE) {
  fac <- setdiff(colnames(returns_past), "date")
  R <- as.matrix(returns_past[, fac])
  mu <- colMeans(R)
  S  <- ledoit_wolf_cov(R)
  n  <- length(mu)

  if (long_only) {
    if (!requireNamespace("quadprog", quietly = TRUE)) install.packages("quadprog")
    # solve.QP minimizes 0.5 w' D w - d' w with A' w >= b, meq equalities
    Dmat <- risk_aversion * S + diag(1e-6, n)
    dvec <- mu
    Amat <- cbind(rep(1, n), diag(n))
    bvec <- c(1, rep(0, n))
    sol <- tryCatch(
      quadprog::solve.QP(Dmat = Dmat, dvec = dvec,
                         Amat = Amat, bvec = bvec, meq = 1),
      error = function(e) NULL)
    if (is.null(sol)) {
      w <- rep(1 / n, n)
    } else {
      w <- sol$solution
      w <- pmax(w, 0)
      w <- w / max(sum(w), 1e-12)
    }
  } else {
    # unconstrained (legacy)
    w <- solve(risk_aversion * S + diag(1e-6, n), mu)
    w <- w / max(sum(abs(w)), 1e-12)
  }
  setNames(w, fac)
}

# =============================================================================
# Volatility scaling
# =============================================================================

#' Scale a return series to a target annualised volatility.
#' 36-month rolling sd (min 6 by default -- shorter warmup than before).
vol_scale <- function(returns, vol_target = 0.10, window = 36, min_obs = 6) {
  n <- length(returns)
  sd_roll <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    lo <- max(1, t - window + 1)
    if (t - lo + 1 >= min_obs) {
      sd_roll[t] <- stats::sd(returns[lo:t], na.rm = TRUE)
    }
  }
  ann_sd <- sd_roll * sqrt(12)
  scale <- vol_target / pmax(ann_sd, 1e-8)
  scale <- pmin(scale, 5)           # leverage cap
  scale_l <- c(NA, scale[-n])       # lag to avoid look-ahead
  returns * scale_l
}
