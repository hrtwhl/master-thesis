###############################################################################
# model_functions.R — Matching Python replication.py exactly
# Key differences from previous R versions:
#   1. select_K: penalty is LAMBDA_K * K (per state), not per free parameter
#   2. select_K: fit on FULL history, score(full) - score(train) = pred LL
#   3. Conditional moments: simple weighted sum (no mixture cross-term)
#   4. Shrinkage: LedoitWolf on raw return sample, not conditional covariance
#   5. No feature standardization, no history window cap
###############################################################################

library(parallel)

# ── Compile C++ ──────────────────────────────────────────────────────────────
if (!requireNamespace("Rcpp", quietly = TRUE) ||
    !requireNamespace("RcppArmadillo", quietly = TRUE))
  stop("Required: install.packages(c('Rcpp', 'RcppArmadillo'))", call. = FALSE)

cpp_file <- file.path(getwd(), "Boukardagha/hmm_rcpp.cpp")
if (!file.exists(cpp_file)) stop("hmm_rcpp.cpp not found in ", getwd(), call. = FALSE)

cat("[model_functions] Compiling C++... ")
t_c <- proc.time(); Rcpp::sourceCpp(cpp_file)
cat(sprintf("done (%.1fs)\n", (proc.time() - t_c)[3]))


# ═══════════════════════════════════════════════════════════════════════════════
# HMM
# ═══════════════════════════════════════════════════════════════════════════════

hmm_init <- function(X, K, seed = SEED) {
  n <- nrow(X); d <- ncol(X)
  set.seed(seed)
  idx <- sample.int(n, 1); centroids <- list(X[idx, ])
  for (k in 2:K) {
    dists <- sapply(centroids, function(mu) rowSums(sweep(X, 2, mu)^2))
    if (is.matrix(dists)) dists <- apply(dists, 1, min)
    probs <- dists / sum(dists); probs[!is.finite(probs)] <- 1/n
    idx <- sample.int(n, 1, prob = probs); centroids[[k]] <- X[idx, ]
  }
  assign <- integer(n)
  for (i in seq_len(n)) {
    dists <- sapply(centroids, function(mu) sum((X[i, ] - mu)^2))
    assign[i] <- which.min(dists)
  }
  mu <- Sigma <- vector("list", K)
  for (k in seq_len(K)) {
    idx_k <- which(assign == k)
    if (length(idx_k) < d + 1) idx_k <- sample.int(n, max(d + 1, 20))
    mu[[k]] <- colMeans(X[idx_k, , drop = FALSE])
    S <- cov(X[idx_k, , drop = FALSE]); S <- (S + t(S)) / 2
    Sigma[[k]] <- S + COV_REG_FACTOR * diag(d)
  }
  list(pi0 = rep(1/K, K), A = matrix(1/K, K, K), mu = mu, Sigma = Sigma)
}

fit_gaussian_hmm <- function(X, K, params = NULL,
                              max_iter = EM_MAX_ITER, tol = EM_TOL) {
  if (is.null(params) || length(params$mu) != K) params <- hmm_init(X, K)
  pi0 <- if (!is.null(params$pi0)) params$pi0 else rep(1/K, K)
  result <- fit_hmm_cpp(X, params$mu, params$Sigma, params$A, pi0,
                         as.integer(max_iter), tol, COV_REG_FACTOR)
  result$mu <- lapply(seq_len(result$K), function(k) result$mu[[k]])
  result$Sigma <- lapply(seq_len(result$K), function(k) result$Sigma[[k]])
  result
}

hmm_filter_step <- function(alpha_prev, x_new, hmm_params) {
  as.numeric(hmm_filter_step_cpp(
    alpha_prev, x_new, hmm_params$mu, hmm_params$Sigma, hmm_params$A))
}

dmvnorm_log <- function(X, mu, Sigma) as.numeric(dmvnorm_log_cpp(X, mu, Sigma))

log_sum_exp <- function(x) {
  x <- x[is.finite(x)]; if (length(x) == 0L) return(-Inf)
  mx <- max(x); mx + log(sum(exp(x - mx)))
}
log_sum_exp_cols <- function(M) {
  mx <- apply(M, 2, max); mx + log(colSums(exp(sweep(M, 2, mx))))
}


# ═══════════════════════════════════════════════════════════════════════════════
# MODEL-ORDER SELECTION — matching Python select_K() exactly
#
# Python approach:
#   1. For each K: fit HMM on FULL X_hist
#   2. pred_ll = model.score(X_hist) - model.score(X_hist[:-val_slice])
#   3. criterion = pred_ll - LAMBDA_K * K   (penalty per STATE)
#   4. Pick K with highest criterion
# ═══════════════════════════════════════════════════════════════════════════════

select_model_order <- function(X_history, val_len = VAL_SLICE_DAYS,
                                k_min = K_MIN, k_max = K_MAX,
                                lambda_k = LAMBDA_K, n_cores = N_CORES) {
  n <- nrow(X_history)
  if (n <= val_len + 50) return(list(K = k_min, scores = NULL))
  
  X_pre <- X_history[1:(n - val_len), , drop = FALSE]  # training portion
  
  score_one_K <- function(K) {
    fit <- tryCatch(fit_gaussian_hmm(X_history, K), error = function(e) NULL)
    if (is.null(fit)) return(-Inf)
    
    # Score full sequence and training-only under the fitted model
    ll_full <- score_sequence_cpp(X_history, fit$mu, fit$Sigma, fit$A, fit$pi0)
    ll_pre  <- score_sequence_cpp(X_pre,     fit$mu, fit$Sigma, fit$A, fit$pi0)
    
    pred_ll <- ll_full - ll_pre
    pred_ll - lambda_k * K   # penalty PER STATE (matching Python)
  }
  
  k_candidates <- k_min:k_max
  use_par <- (.Platform$OS.type == "unix") && (n_cores > 1L) && (length(k_candidates) > 1L)
  
  if (use_par) {
    scores <- unlist(parallel::mclapply(k_candidates, score_one_K,
                                         mc.cores = min(n_cores, length(k_candidates))))
  } else {
    scores <- sapply(k_candidates, score_one_K)
  }
  names(scores) <- k_candidates
  list(K = as.integer(names(which.max(scores))), scores = scores)
}


# ═══════════════════════════════════════════════════════════════════════════════
# WASSERSTEIN TEMPLATES
# ═══════════════════════════════════════════════════════════════════════════════

wasserstein2_gaussian <- function(mu1, Sigma1, mu2, Sigma2) {
  d <- length(mu1); mean_sq <- sum((mu1 - mu2)^2)
  eig2 <- eigen(Sigma2, symmetric = TRUE)
  sqrt_S2 <- eig2$vectors %*% diag(sqrt(pmax(eig2$values, 0)), d) %*% t(eig2$vectors)
  M <- sqrt_S2 %*% Sigma1 %*% sqrt_S2; M <- (M + t(M)) / 2
  eig_M <- eigen(M, symmetric = TRUE)
  sqrt_M <- eig_M$vectors %*% diag(sqrt(pmax(eig_M$values, 0)), d) %*% t(eig_M$vectors)
  sqrt(mean_sq + max(sum(diag(Sigma1)) + sum(diag(Sigma2)) - 2 * sum(diag(sqrt_M)), 0))
}

#' Initialize templates by fitting K=G on the init window
#' (matching Python: model = _fit_hmm(X_init, G_TEMPLATES))
init_templates <- function(X_init, G = G_TEMPLATES) {
  fit <- fit_gaussian_hmm(X_init, G)
  templates <- vector("list", G)
  for (g in seq_len(G)) {
    templates[[g]] <- list(mu = fit$mu[[g]], Sigma = fit$Sigma[[g]], active = TRUE)
  }
  list(templates = templates, model = fit)
}

assign_to_templates <- function(comp_means, comp_covs, templates) {
  K <- length(comp_means); G <- length(templates)
  mapping <- integer(K)
  for (k in seq_len(K)) {
    dists <- sapply(seq_len(G), function(g)
      wasserstein2_gaussian(comp_means[[k]], comp_covs[[k]],
                             templates[[g]]$mu, templates[[g]]$Sigma))
    mapping[k] <- which.min(dists)
  }
  mapping
}

update_templates <- function(templates, comp_means, comp_covs, comp_probs,
                              g_of_k, eta = ETA_SMOOTH) {
  G <- length(templates)
  for (g in seq_len(G)) {
    mask <- which(g_of_k == g)
    if (length(mask) == 0) next
    w <- comp_probs[mask]; tot <- sum(w)
    if (tot < 1e-12) next
    w <- w / tot
    mean_bar <- rep(0, length(templates[[g]]$mu))
    cov_bar  <- matrix(0, nrow(templates[[g]]$Sigma), ncol(templates[[g]]$Sigma))
    for (i in seq_along(mask)) {
      mean_bar <- mean_bar + w[i] * comp_means[[mask[i]]]
      cov_bar  <- cov_bar  + w[i] * comp_covs[[mask[i]]]
    }
    templates[[g]]$mu    <- (1 - eta) * templates[[g]]$mu    + eta * mean_bar
    templates[[g]]$Sigma <- (1 - eta) * templates[[g]]$Sigma + eta * cov_bar
  }
  templates
}


# ═══════════════════════════════════════════════════════════════════════════════
# CONDITIONAL MOMENTS — matching Python: simple weighted sum, no cross-term
# ═══════════════════════════════════════════════════════════════════════════════

compute_conditional_moments <- function(templates, p_template, n_assets) {
  d <- length(templates[[1]]$mu)
  mu_full <- rep(0, d)
  cov_full <- matrix(0, d, d)
  for (g in seq_along(templates)) {
    mu_full  <- mu_full  + p_template[g] * templates[[g]]$mu
    cov_full <- cov_full + p_template[g] * templates[[g]]$Sigma
  }
  idx <- seq_len(n_assets)
  list(mu = mu_full[idx],
       Sigma = cov_full[idx, idx, drop = FALSE])
}


# ═══════════════════════════════════════════════════════════════════════════════
# SHRINKAGE — matching Python: LedoitWolf on raw return SAMPLE
# ═══════════════════════════════════════════════════════════════════════════════

#' Analytical Ledoit-Wolf shrinkage (toward scaled identity)
#' Applied to the RAW return sample, not the conditional covariance
ledoit_wolf_shrink <- function(sample_matrix) {
  n <- nrow(sample_matrix); p <- ncol(sample_matrix)
  if (n < 30 || p <= 1) {
    S <- cov(sample_matrix)
    return(0.9 * S + 0.1 * diag(diag(S)))
  }
  S <- cov(sample_matrix)
  tr_S <- sum(diag(S))
  target <- (tr_S / p) * diag(p)
  
  # Optimal shrinkage intensity (Ledoit-Wolf 2004)
  X_c <- scale(sample_matrix, scale = FALSE)
  S2 <- crossprod(X_c) / n
  delta <- S2 - S
  d2 <- sum(delta^2) / p
  b2 <- 0
  for (i in seq_len(n)) {
    xi <- tcrossprod(X_c[i, ])
    b2 <- b2 + sum((xi - S2)^2)
  }
  b2 <- b2 / (n^2 * p)
  intensity <- min(1, max(0, b2 / d2))
  
  (1 - intensity) * S + intensity * target
}

shrink_cov <- function(cond_cov, R_hist, n_sample = 750) {
  if (!LEDOIT_WOLF) return(cond_cov)
  n <- nrow(R_hist)
  if (n >= 30) {
    sample <- R_hist[max(1, n - n_sample + 1):n, , drop = FALSE]
    return(ledoit_wolf_shrink(sample))
  }
  0.9 * cond_cov + 0.1 * diag(diag(cond_cov))
}


# ═══════════════════════════════════════════════════════════════════════════════
# MVO — exact L1 via SLSQP-style slack variables (matching Python)
# ═══════════════════════════════════════════════════════════════════════════════

#' PSD projection
psd_project <- function(S, eps = 1e-10) {
  S <- (S + t(S)) / 2
  eig <- eigen(S, symmetric = TRUE)
  eig$vectors %*% diag(pmax(eig$values, eps)) %*% t(eig$vectors)
}

#' Solve MVO with exact L1 TC penalty
#' Uses augmented-variable reformulation:
#'   w - w_prev = d_pos - d_neg, ||w-w_prev||_1 = sum(d_pos + d_neg)
#' Objective becomes smooth QP, solved with L-BFGS-B + equality via penalty
solve_mvo <- function(mu, Sigma, w_prev, gamma = GAMMA_RISK,
                       tau = TAU_TCOST, w_max = W_MAX) {
  N <- length(mu)
  Sigma <- psd_project(Sigma)
  
  # Variable: x = [w(N), d_pos(N), d_neg(N)]
  # Objective: -mu'w + gamma*w'Sigma*w + tau*sum(d_pos + d_neg)
  # Equality: w - w_prev = d_pos - d_neg; sum(w) = 1
  # Bounds: 0 <= w <= w_max; d_pos, d_neg >= 0
  
  objective <- function(x) {
    w <- x[1:N]
    dp <- x[(N+1):(2*N)]; dn <- x[(2*N+1):(3*N)]
    # Core objective
    obj <- -sum(mu * w) + gamma * as.numeric(t(w) %*% Sigma %*% w) + tau * sum(dp + dn)
    # Penalty for equality constraints
    eq_sum <- sum(w) - 1
    eq_slack <- w - w_prev - (dp - dn)
    penalty <- 1e4 * (eq_sum^2 + sum(eq_slack^2))
    obj + penalty
  }
  
  gradient <- function(x) {
    w <- x[1:N]
    dp <- x[(N+1):(2*N)]; dn <- x[(2*N+1):(3*N)]
    eq_sum <- sum(w) - 1
    eq_slack <- w - w_prev - (dp - dn)
    
    grad_w  <- -mu + 2 * gamma * as.numeric(Sigma %*% w) +
                1e4 * (2 * eq_sum + 2 * eq_slack)
    grad_dp <- rep(tau, N) + 1e4 * (-2 * eq_slack)
    grad_dn <- rep(tau, N) + 1e4 * ( 2 * eq_slack)
    c(grad_w, grad_dp, grad_dn)
  }
  
  x0 <- c(w_prev, pmax(w_prev - w_prev, 0), pmax(w_prev - w_prev, 0))  # all zeros for slacks
  lower <- c(rep(0, N), rep(0, 2*N))
  upper <- c(rep(w_max, N), rep(Inf, 2*N))
  
  result <- optim(par = x0, fn = objective, gr = gradient,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   control = list(maxit = 300, factr = 1e7))
  
  w <- result$par[1:N]
  w <- pmax(w, 0)
  s <- sum(w)
  if (s < 1e-10 || any(!is.finite(w))) return(w_prev)
  w / s
}


# ═══════════════════════════════════════════════════════════════════════════════
# SPEED CHECK
# ═══════════════════════════════════════════════════════════════════════════════

.run_speed_check <- function() {
  d <- 15; n <- 500; K <- 3
  X_test <- matrix(rnorm(n * d), n, d)
  t0 <- proc.time()
  fit <- fit_gaussian_hmm(X_test, K, max_iter = 20L, tol = 1e-3)
  el <- (proc.time() - t0)[3]
  cat(sprintf("[model_functions] Speed: HMM(n=%d,d=%d,K=%d,%d iters)=%.3fs\n",
              n, d, K, fit$n_iter, el))
  if (el > 1.0) cat("  WARNING: >1s, expected <0.1s with C++\n")
}
.run_speed_check()
cat("[model_functions.R] Loaded (Rcpp/C++, Python-matched).\n")
