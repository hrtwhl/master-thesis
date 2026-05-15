###############################################################################
# model_functions.R — Optimized: dual param storage, C++ MVO
###############################################################################

library(parallel)

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

#' Fit Gaussian HMM with internal standardization.
#' Returns BOTH original-scale and z-scale params to avoid re-standardization.
fit_gaussian_hmm <- function(X, K, params = NULL,
                              max_iter = EM_MAX_ITER, tol = EM_TOL,
                              n_restarts = HMM_N_RESTARTS) {
  n <- nrow(X); d <- ncol(X)
  
  center <- colMeans(X)
  sds <- apply(X, 2, sd)
  sds[sds < 1e-10] <- 1e-10
  X_z <- scale(X, center = center, scale = sds)
  
  best <- NULL; best_ll <- -Inf
  for (r in seq_len(n_restarts)) {
    init <- hmm_init(X_z, K, seed = SEED + r - 1L)
    result <- tryCatch(
      fit_hmm_cpp(X_z, init$mu, init$Sigma, init$A, init$pi0,
                   as.integer(max_iter), tol, COV_REG_FACTOR),
      error = function(e) NULL)
    if (is.null(result)) next
    if (is.finite(result$log_lik) && result$log_lik > best_ll) {
      best_ll <- result$log_lik; best <- result
    }
  }
  if (is.null(best)) return(NULL)
  
  best$mu <- lapply(seq_len(best$K), function(k) best$mu[[k]])
  best$Sigma <- lapply(seq_len(best$K), function(k) best$Sigma[[k]])
  
  # Store z-space params (for scoring + filtering — avoids re-standardization)
  best$mu_z <- best$mu
  best$Sigma_z <- best$Sigma
  best$z_center <- center
  best$z_sds <- sds
  
  # De-standardize to original scale (for template tracking)
  S_diag <- diag(sds)
  for (k in seq_len(best$K)) {
    best$mu[[k]] <- best$mu_z[[k]] * sds + center
    best$Sigma[[k]] <- S_diag %*% best$Sigma_z[[k]] %*% S_diag
  }
  
  best
}

#' Score sequence using stored z-space params (no re-standardization)
score_hmm <- function(X, fit) {
  X_z <- scale(X, center = fit$z_center, scale = fit$z_sds)
  score_sequence_cpp(X_z, fit$mu_z, fit$Sigma_z, fit$A, fit$pi0)
}

#' Filter step using stored z-space params (no re-standardization)
hmm_filter_step <- function(alpha_prev, x_new, hmm_params) {
  x_z <- (x_new - hmm_params$z_center) / hmm_params$z_sds
  as.numeric(hmm_filter_step_cpp(alpha_prev, x_z,
                                  hmm_params$mu_z, hmm_params$Sigma_z, hmm_params$A))
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
# MODEL-ORDER SELECTION
# ═══════════════════════════════════════════════════════════════════════════════

select_model_order <- function(X_history, val_len = VAL_SLICE_DAYS,
                                k_min = K_MIN, k_max = K_MAX,
                                lambda_k = LAMBDA_K, n_cores = N_CORES) {
  n <- nrow(X_history)
  if (n <= val_len + 50) return(list(K = k_min, scores = NULL))
  X_pre <- X_history[1:(n - val_len), , drop = FALSE]
  
  score_one_K <- function(K) {
    fit <- fit_gaussian_hmm(X_history, K, n_restarts = 1L)  # fewer restarts for speed
    if (is.null(fit)) return(-Inf)
    ll_full <- score_hmm(X_history, fit)
    ll_pre  <- score_hmm(X_pre, fit)
    (ll_full - ll_pre) - lambda_k * K
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
  best <- which.max(scores)
  if (length(best) == 0) return(list(K = k_min, scores = scores))
  list(K = as.integer(names(best)), scores = scores)
}


# ═══════════════════════════════════════════════════════════════════════════════
# WASSERSTEIN TEMPLATES (original scale)
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

init_templates <- function(X_init, G = G_TEMPLATES) {
  fit <- fit_gaussian_hmm(X_init, G)
  if (is.null(fit)) stop("Template initialization failed")
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
    d <- length(templates[[g]]$mu)
    mean_bar <- rep(0, d); cov_bar <- matrix(0, d, d)
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
# CONDITIONAL MOMENTS + SHRINKAGE + MVO (C++)
# ═══════════════════════════════════════════════════════════════════════════════

compute_conditional_moments <- function(templates, p_template, n_assets) {
  d <- length(templates[[1]]$mu)
  mu_full <- rep(0, d); cov_full <- matrix(0, d, d)
  for (g in seq_along(templates)) {
    mu_full  <- mu_full  + p_template[g] * templates[[g]]$mu
    cov_full <- cov_full + p_template[g] * templates[[g]]$Sigma
  }
  idx <- seq_len(n_assets)
  list(mu = mu_full[idx], Sigma = cov_full[idx, idx, drop = FALSE])
}

ledoit_wolf_shrink <- function(sample_matrix) {
  n <- nrow(sample_matrix); p <- ncol(sample_matrix)
  if (n < 30 || p <= 1) {
    S <- cov(sample_matrix)
    return(0.9 * S + 0.1 * diag(diag(S)))
  }
  S <- cov(sample_matrix)
  tr_S <- sum(diag(S)); target <- (tr_S / p) * diag(p)
  X_c <- scale(sample_matrix, scale = FALSE)
  S2 <- crossprod(X_c) / n
  d2 <- sum((S2 - S)^2) / p
  b2 <- 0
  for (i in seq_len(n)) { xi <- tcrossprod(X_c[i, ]); b2 <- b2 + sum((xi - S2)^2) }
  b2 <- b2 / (n^2 * p)
  intensity <- min(1, max(0, b2 / max(d2, 1e-20)))
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

#' MVO — calls C++ solver
solve_mvo <- function(mu, Sigma, w_prev, gamma = GAMMA_RISK,
                       tau = TAU_TCOST, w_max = W_MAX) {
  as.numeric(solve_mvo_cpp(mu, Sigma, w_prev, gamma, tau, w_max))
}


# ═══════════════════════════════════════════════════════════════════════════════
# SPEED CHECK
# ═══════════════════════════════════════════════════════════════════════════════
.run_speed_check <- function() {
  X_test <- matrix(rnorm(500 * 15), 500, 15)
  t0 <- proc.time()
  fit <- fit_gaussian_hmm(X_test, 4L, max_iter = 20L, tol = 1e-3)
  el <- (proc.time() - t0)[3]
  cat(sprintf("[speed] HMM(n=500,d=15,K=4,5 restarts)=%.3fs\n", el))
  
  mu <- rnorm(5) * 0.001; Sigma <- cov(matrix(rnorm(250*5), 250, 5))
  w0 <- rep(0.2, 5)
  t0 <- proc.time()
  for (i in 1:1000) solve_mvo(mu, Sigma, w0)
  el <- (proc.time() - t0)[3]
  cat(sprintf("[speed] MVO x1000 = %.3fs (%.0f µs/call)\n", el, el * 1000))
}
.run_speed_check()
cat("[model_functions.R] Loaded (Rcpp/C++, C++ MVO).\n")
