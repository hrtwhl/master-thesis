###############################################################################
# model_functions.R — HMM Engine, Wasserstein Templates, MVO
# Boukardagha (2026) Replication
#
# Performance optimizations (no C++, no parameter changes):
#   1. parallel::mclapply for model-order selection (K candidates in parallel)
#   2. compiler::cmpfun byte-compilation on all hot-path functions
#   3. Pre-cached Cholesky factors inside EM to avoid redundant decomposition
#   4. Vectorized forward-backward over K dimension
###############################################################################

library(compiler)
library(parallel)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. GAUSSIAN HMM (BAUM-WELCH EM)
# ═══════════════════════════════════════════════════════════════════════════════

#' K-means++ initialization for HMM
hmm_init <- function(X, K, seed = NULL) {
  n <- nrow(X); d <- ncol(X)
  if (!is.null(seed)) set.seed(seed)
  
  idx <- sample.int(n, 1)
  centroids <- list(X[idx, ])
  for (k in 2:K) {
    dists <- sapply(centroids, function(mu) rowSums(sweep(X, 2, mu)^2))
    if (is.matrix(dists)) dists <- apply(dists, 1, min)
    probs <- dists / sum(dists)
    probs[!is.finite(probs)] <- 1 / n
    idx <- sample.int(n, 1, prob = probs)
    centroids[[k]] <- X[idx, ]
  }
  
  assign <- integer(n)
  for (i in seq_len(n)) {
    dists <- sapply(centroids, function(mu) sum((X[i, ] - mu)^2))
    assign[i] <- which.min(dists)
  }
  
  mu <- vector("list", K)
  Sigma <- vector("list", K)
  for (k in seq_len(K)) {
    idx_k <- which(assign == k)
    if (length(idx_k) < d + 1) idx_k <- sample.int(n, max(d + 1, 20))
    mu[[k]] <- colMeans(X[idx_k, , drop = FALSE])
    Sigma[[k]] <- cov(X[idx_k, , drop = FALSE])
    Sigma[[k]] <- regularize_cov(Sigma[[k]], d)
  }
  
  list(pi0 = rep(1 / K, K), A = matrix(1 / K, K, K), mu = mu, Sigma = Sigma)
}

#' Regularize covariance toward scaled identity
regularize_cov <- function(S, d, eps = COV_REG_FACTOR) {
  S <- (S + t(S)) / 2
  S + eps * diag(d)
}

#' Log-density of MVN using pre-computed Cholesky factor L (upper triangular)
#' Avoids redundant chol() calls when the same Sigma is used many times.
dmvnorm_log_chol <- function(X, mu, L, log_det) {
  d <- length(mu)
  X_c <- sweep(X, 2, mu)
  Z <- forwardsolve(t(L), t(X_c))
  mahal <- colSums(Z^2)
  -0.5 * (d * log(2 * pi) + log_det + mahal)
}

#' Log-density of MVN (self-contained, for external callers)
dmvnorm_log <- function(X, mu, Sigma) {
  d <- length(mu)
  L <- tryCatch(chol(Sigma), error = function(e) chol(Sigma + diag(1e-4, d)))
  log_det <- 2 * sum(log(diag(L)))
  dmvnorm_log_chol(X, mu, L, log_det)
}

#' Log-sum-exp (numerically stable)
log_sum_exp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(-Inf)
  mx <- max(x)
  mx + log(sum(exp(x - mx)))
}

#' Vectorized log-sum-exp over columns
log_sum_exp_cols <- function(M) {
  mx <- apply(M, 2, max)
  mx + log(colSums(exp(sweep(M, 2, mx))))
}

#' Fit Gaussian HMM via Baum-Welch EM
#' Hot path — byte-compiled after definition.
fit_gaussian_hmm <- function(X, K, params = NULL,
                              max_iter = EM_MAX_ITER, tol = EM_TOL) {
  n <- nrow(X)
  d <- ncol(X)
  
  if (is.null(params) || length(params$mu) != K) {
    params <- hmm_init(X, K)
  }
  
  pi0   <- if (!is.null(params$pi0)) params$pi0 else rep(1 / K, K)
  A     <- params$A
  mu    <- params$mu
  Sigma <- params$Sigma
  
  log_A <- log(A + 1e-300)
  prev_ll <- -Inf
  ll <- -Inf
  iter <- 1L
  
  for (iter in seq_len(max_iter)) {
    # ── Pre-compute Cholesky factors (reused in emission + pred_log_score) ──
    L_list      <- vector("list", K)
    log_det_vec <- numeric(K)
    for (k in seq_len(K)) {
      L_list[[k]] <- tryCatch(chol(Sigma[[k]]),
                               error = function(e) chol(Sigma[[k]] + diag(1e-4, d)))
      log_det_vec[k] <- 2 * sum(log(diag(L_list[[k]])))
    }
    
    # ── E-step: emission log-probs using cached Cholesky ──
    B <- matrix(NA_real_, n, K)
    for (k in seq_len(K)) {
      B[, k] <- dmvnorm_log_chol(X, mu[[k]], L_list[[k]], log_det_vec[k])
    }
    B[!is.finite(B)] <- -1e10
    
    # ── Forward pass (vectorized over K) ──
    log_alpha <- matrix(-Inf, n, K)
    log_alpha[1L, ] <- log(pi0 + 1e-300) + B[1L, ]
    
    for (t in 2:n) {
      M <- matrix(log_alpha[t - 1L, ], K, K) + log_A
      log_alpha[t, ] <- log_sum_exp_cols(M) + B[t, ]
    }
    
    # ── Backward pass (vectorized over K) ──
    log_beta <- matrix(-Inf, n, K)
    log_beta[n, ] <- 0
    
    for (t in (n - 1L):1L) {
      vals <- B[t + 1L, ] + log_beta[t + 1L, ]
      M <- log_A + matrix(vals, K, K, byrow = TRUE)
      log_beta[t, ] <- apply(M, 1, log_sum_exp)
    }
    
    # ── Posterior gamma ──
    log_gamma <- log_alpha + log_beta
    log_normalizer <- apply(log_gamma, 1, log_sum_exp)
    log_gamma <- log_gamma - log_normalizer
    gamma <- exp(log_gamma)
    gamma[!is.finite(gamma)] <- 1 / K
    gamma <- pmax(gamma, 1e-300)
    gamma <- gamma / rowSums(gamma)
    
    ll <- log_sum_exp(log_alpha[n, ])
    if (!is.finite(ll)) break
    if (abs(ll - prev_ll) < tol && iter > 2L) break
    prev_ll <- ll
    
    # ── M-step ──
    pi0 <- pmax(gamma[1L, ], 1e-6)
    pi0 <- pi0 / sum(pi0)
    
    A_new <- matrix(1e-10, K, K)
    for (j in seq_len(K)) {
      for (k in seq_len(K)) {
        log_xi_sum <- log_sum_exp(
          log_alpha[1:(n - 1L), j] + log_A[j, k] +
            B[2:n, k] + log_beta[2:n, k]
        )
        if (is.finite(log_xi_sum)) A_new[j, k] <- exp(log_xi_sum)
      }
    }
    A <- pmax(A_new, 1e-10)
    A <- A / rowSums(A)
    log_A <- log(A + 1e-300)
    
    for (k in seq_len(K)) {
      w <- gamma[, k]
      w_sum <- sum(w)
      if (!is.finite(w_sum) || w_sum < 1e-6) next
      mu[[k]] <- colSums(w * X) / w_sum
      if (any(!is.finite(mu[[k]]))) { mu[[k]] <- colMeans(X) }
      X_c <- sweep(X, 2, mu[[k]])
      Sigma[[k]] <- crossprod(X_c * sqrt(w)) / w_sum
      Sigma[[k]] <- regularize_cov(Sigma[[k]], d)
    }
  }
  
  # Filtered probabilities (forward only)
  log_norm <- apply(log_alpha, 1, log_sum_exp)
  filtered <- exp(log_alpha - log_norm)
  filtered[!is.finite(filtered)] <- 1 / K
  filtered <- pmax(filtered, 0)
  filtered <- filtered / rowSums(filtered)
  
  list(
    pi0 = pi0, A = A, mu = mu, Sigma = Sigma,
    K = K, d = d, log_lik = ll,
    gamma = gamma, filtered = filtered,
    n_iter = min(iter, max_iter)
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. MODEL-ORDER SELECTION (PARALLEL)
# ═══════════════════════════════════════════════════════════════════════════════

#' Predictive log-score for a single K
pred_log_score <- function(X_train, X_valid, K, lambda_k = LAMBDA_K) {
  d <- ncol(X_train)
  
  fit <- tryCatch(
    fit_gaussian_hmm(X_train, K, max_iter = 30L, tol = 1e-3),
    error = function(e) NULL
  )
  if (is.null(fit)) return(-Inf)
  
  pred_ll <- 0
  alpha_prev <- fit$filtered[nrow(X_train), ]
  
  # Pre-compute Cholesky factors for validation scoring
  L_list <- vector("list", K)
  log_det_vec <- numeric(K)
  for (k in seq_len(K)) {
    L_list[[k]] <- tryCatch(chol(fit$Sigma[[k]]),
                             error = function(e) chol(fit$Sigma[[k]] + diag(1e-4, d)))
    log_det_vec[k] <- 2 * sum(log(diag(L_list[[k]])))
  }
  
  for (s in seq_len(nrow(X_valid))) {
    pred_state <- as.numeric(alpha_prev %*% fit$A)
    pred_state <- pmax(pred_state, 1e-300)
    
    log_emission <- sapply(seq_len(K), function(k) {
      dmvnorm_log_chol(matrix(X_valid[s, ], nrow = 1),
                        fit$mu[[k]], L_list[[k]], log_det_vec[k])
    })
    
    log_pred <- log_sum_exp(log(pred_state) + log_emission)
    pred_ll <- pred_ll + log_pred
    
    log_alpha_new <- log(pred_state) + log_emission
    alpha_prev <- exp(log_alpha_new - log_sum_exp(log_alpha_new))
    alpha_prev[!is.finite(alpha_prev)] <- 1 / K
  }
  
  n_params <- (K - 1) + K * (K - 1) + K * d + K * d * (d + 1) / 2
  pred_ll - lambda_k * n_params
}

#' Select optimal K via predictive scoring — PARALLEL over K candidates
#' Uses mclapply (fork-based, Unix/macOS) which avoids the overhead of
#' socket clusters and variable export. Falls back to sequential on Windows.
select_model_order <- function(X_history, val_len = VALIDATION_LEN,
                                k_min = K_MIN, k_max = K_MAX,
                                lambda_k = LAMBDA_K,
                                n_cores = N_CORES) {
  n <- nrow(X_history)
  if (n <= val_len + 50) return(list(K = k_min, scores = NULL))
  
  X_train <- X_history[1:(n - val_len), , drop = FALSE]
  X_valid <- X_history[(n - val_len + 1):n, , drop = FALSE]
  
  k_candidates <- k_min:k_max
  n_cand <- length(k_candidates)
  use_parallel <- (.Platform$OS.type == "unix") && (n_cores > 1L) && (n_cand > 1L)
  
  if (use_parallel) {
    # mclapply: fork-based, inherits all objects, no export needed
    scores_list <- mclapply(k_candidates, function(K) {
      pred_log_score(X_train, X_valid, K, lambda_k)
    }, mc.cores = min(n_cores, n_cand))
    scores <- unlist(scores_list)
  } else {
    scores <- sapply(k_candidates, function(K) {
      pred_log_score(X_train, X_valid, K, lambda_k)
    })
  }
  names(scores) <- k_candidates
  
  list(K = as.integer(names(which.max(scores))), scores = scores)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. WASSERSTEIN DISTANCE AND TEMPLATE TRACKING
# ═══════════════════════════════════════════════════════════════════════════════

#' Closed-form 2-Wasserstein distance between two Gaussians
wasserstein2_gaussian <- function(mu1, Sigma1, mu2, Sigma2) {
  d <- length(mu1)
  mean_sq <- sum((mu1 - mu2)^2)
  
  eig2 <- eigen(Sigma2, symmetric = TRUE)
  vals2 <- pmax(eig2$values, 0)
  sqrt_S2 <- eig2$vectors %*% diag(sqrt(vals2), d) %*% t(eig2$vectors)
  
  M <- sqrt_S2 %*% Sigma1 %*% sqrt_S2
  M <- (M + t(M)) / 2
  eig_M <- eigen(M, symmetric = TRUE)
  sqrt_M <- eig_M$vectors %*% diag(sqrt(pmax(eig_M$values, 0)), d) %*% t(eig_M$vectors)
  
  cov_term <- sum(diag(Sigma1)) + sum(diag(Sigma2)) - 2 * sum(diag(sqrt_M))
  cov_term <- max(cov_term, 0)
  
  sqrt(mean_sq + cov_term)
}

#' Initialize G template slots from first HMM fit
init_templates <- function(hmm_fit, G = G_TEMPLATES) {
  K_fit <- hmm_fit$K; d <- hmm_fit$d
  templates <- vector("list", G)
  for (g in seq_len(G)) {
    if (g <= K_fit) {
      templates[[g]] <- list(mu = hmm_fit$mu[[g]], Sigma = hmm_fit$Sigma[[g]], active = TRUE)
    } else {
      templates[[g]] <- list(
        mu = colMeans(do.call(rbind, hmm_fit$mu)),
        Sigma = Reduce("+", hmm_fit$Sigma) / K_fit * 2,
        active = FALSE
      )
    }
  }
  templates
}

#' Map HMM components to templates via W2 distance
map_components_to_templates <- function(hmm_fit, templates) {
  K <- hmm_fit$K; G <- length(templates)
  
  dist_mat <- matrix(Inf, K, G)
  for (k in seq_len(K)) {
    for (g in seq_len(G)) {
      dist_mat[k, g] <- wasserstein2_gaussian(
        hmm_fit$mu[[k]], hmm_fit$Sigma[[k]],
        templates[[g]]$mu, templates[[g]]$Sigma
      )
    }
  }
  
  mapping <- apply(dist_mat, 1, which.min)
  
  filtered_last <- hmm_fit$filtered[nrow(hmm_fit$filtered), ]
  template_probs <- rep(0, G)
  for (k in seq_len(K)) {
    template_probs[mapping[k]] <- template_probs[mapping[k]] + filtered_last[k]
  }
  
  d <- hmm_fit$d
  template_mu_avg <- vector("list", G)
  template_Sigma_avg <- vector("list", G)
  
  for (g in seq_len(G)) {
    assigned <- which(mapping == g)
    if (length(assigned) == 0) {
      template_mu_avg[[g]] <- templates[[g]]$mu
      template_Sigma_avg[[g]] <- templates[[g]]$Sigma
      next
    }
    w <- filtered_last[assigned]; w_sum <- sum(w)
    if (w_sum < 1e-10) {
      template_mu_avg[[g]] <- templates[[g]]$mu
      template_Sigma_avg[[g]] <- templates[[g]]$Sigma
      next
    }
    w <- w / w_sum
    mu_avg <- rep(0, d); Sigma_avg <- matrix(0, d, d)
    for (i in seq_along(assigned)) {
      mu_avg <- mu_avg + w[i] * hmm_fit$mu[[assigned[i]]]
      Sigma_avg <- Sigma_avg + w[i] * hmm_fit$Sigma[[assigned[i]]]
    }
    template_mu_avg[[g]] <- mu_avg
    template_Sigma_avg[[g]] <- Sigma_avg
  }
  
  list(mapping = mapping, template_probs = template_probs, dist_mat = dist_mat,
       template_mu_avg = template_mu_avg, template_Sigma_avg = template_Sigma_avg)
}

#' Update templates via exponential smoothing
update_templates <- function(templates, map_result, eta = ETA_SMOOTH) {
  G <- length(templates)
  for (g in seq_len(G)) {
    if (map_result$template_probs[g] > 1e-6) {
      templates[[g]]$mu <- (1 - eta) * templates[[g]]$mu +
        eta * map_result$template_mu_avg[[g]]
      templates[[g]]$Sigma <- (1 - eta) * templates[[g]]$Sigma +
        eta * map_result$template_Sigma_avg[[g]]
      templates[[g]]$Sigma <- (templates[[g]]$Sigma + t(templates[[g]]$Sigma)) / 2
      templates[[g]]$active <- TRUE
    }
  }
  templates
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. CONDITIONAL MOMENTS AND PORTFOLIO OPTIMIZATION
# ═══════════════════════════════════════════════════════════════════════════════

compute_conditional_moments <- function(templates, template_probs, n_assets) {
  G <- length(templates); d <- length(templates[[1]]$mu)
  
  mu_full <- rep(0, d); Sigma_full <- matrix(0, d, d)
  for (g in seq_len(G)) {
    if (template_probs[g] < 1e-10) next
    mu_full <- mu_full + template_probs[g] * templates[[g]]$mu
    Sigma_full <- Sigma_full + template_probs[g] * (
      templates[[g]]$Sigma + tcrossprod(templates[[g]]$mu - mu_full)
    )
  }
  
  ret_idx <- seq_len(n_assets)
  mu_ret <- mu_full[ret_idx]
  Sigma_ret <- Sigma_full[ret_idx, ret_idx, drop = FALSE]
  Sigma_ret <- (Sigma_ret + t(Sigma_ret)) / 2
  
  list(mu = mu_ret, Sigma = Sigma_ret)
}

ledoit_wolf_shrink <- function(S, n) {
  p <- nrow(S)
  if (p <= 1) return(S)
  tr_S <- sum(diag(S))
  target <- (tr_S / p) * diag(p)
  delta <- S - target
  d2 <- sum(delta^2) / p
  intensity <- min(1, max(0, d2 / (n * sum(diag(S %*% S)) / p)))
  (1 - intensity) * S + intensity * target
}

solve_mvo <- function(mu, Sigma, w_prev,
                       gamma = GAMMA_RISK, tau = TAU_TCOST,
                       w_max = W_MAX) {
  n <- length(mu)
  eps_smooth <- 1e-6
  
  objective <- function(w) {
    port_ret  <- sum(mu * w)
    port_risk <- as.numeric(t(w) %*% Sigma %*% w)
    tc_cost   <- sum(sqrt((w - w_prev)^2 + eps_smooth))
    -(port_ret - gamma * port_risk - tau * tc_cost)
  }
  
  gradient <- function(w) {
    grad_ret  <- mu
    grad_risk <- 2 * as.numeric(Sigma %*% w)
    diff_w    <- w - w_prev
    grad_tc   <- diff_w / sqrt(diff_w^2 + eps_smooth)
    -(grad_ret - gamma * grad_risk - tau * grad_tc)
  }
  
  result <- optim(
    par = w_prev, fn = objective, gr = gradient,
    method = "L-BFGS-B",
    lower = rep(0, n), upper = rep(w_max, n),
    control = list(maxit = 200, factr = 1e7)
  )
  
  w_opt <- pmax(result$par, 0)
  s <- sum(w_opt)
  if (s < 1e-10) return(rep(1 / n, n))
  w_opt <- w_opt / s
  w_opt <- pmin(w_opt, w_max)
  w_opt / sum(w_opt)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. DAILY FILTER STEP
# ═══════════════════════════════════════════════════════════════════════════════

hmm_filter_step <- function(alpha_prev, x_new, hmm_params) {
  K <- hmm_params$K
  pred <- as.numeric(alpha_prev %*% hmm_params$A)
  pred <- pmax(pred, 1e-300)
  
  log_emission <- sapply(seq_len(K), function(k) {
    dmvnorm_log(matrix(x_new, nrow = 1), hmm_params$mu[[k]], hmm_params$Sigma[[k]])
  })
  
  log_alpha_new <- log(pred) + log_emission
  mx <- max(log_alpha_new[is.finite(log_alpha_new)])
  if (!is.finite(mx)) return(rep(1 / K, K))
  
  alpha_new <- exp(log_alpha_new - mx)
  alpha_new <- alpha_new / sum(alpha_new)
  alpha_new[!is.finite(alpha_new)] <- 1 / K
  alpha_new
}


# ═══════════════════════════════════════════════════════════════════════════════
# 6. BYTE-COMPILE HOT FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════
# compiler::cmpfun JIT-compiles R bytecode. Typical 2-4× on tight numeric loops.
# No external dependencies, no C++.

dmvnorm_log_chol    <- cmpfun(dmvnorm_log_chol)
dmvnorm_log         <- cmpfun(dmvnorm_log)
log_sum_exp         <- cmpfun(log_sum_exp)
log_sum_exp_cols    <- cmpfun(log_sum_exp_cols)
fit_gaussian_hmm    <- cmpfun(fit_gaussian_hmm)
pred_log_score      <- cmpfun(pred_log_score)
hmm_filter_step     <- cmpfun(hmm_filter_step)
wasserstein2_gaussian <- cmpfun(wasserstein2_gaussian)
solve_mvo           <- cmpfun(solve_mvo)

cat("[model_functions.R] Engine loaded (byte-compiled, parallel model selection).\n")
