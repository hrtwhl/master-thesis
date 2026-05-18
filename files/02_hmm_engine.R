# =============================================================================
#  02_hmm_engine.R
#  R wrappers around the C++ HMM kernels + predictive K selection + MVO QP
# =============================================================================

# ---------------- 1. HMM wrappers ----------------

#' Fit a Gaussian HMM with full covariance (Baum-Welch via Rcpp).
#' Returns a list with means/covars/transmat/startprob and the final log-likelihood.
fit_hmm <- function(X, K, n_iter = 100, tol = 1e-2, seed = 42) {
  storage.mode(X) <- "double"
  gaussian_hmm_em(X, as.integer(K), as.integer(n_iter), as.numeric(tol),
                  as.integer(seed))
}

#' Log-likelihood of a fresh sequence given a fitted model. Used for predictive
#' K selection (validation slice scoring).
score_hmm <- function(X, model) {
  storage.mode(X) <- "double"
  gaussian_hmm_score(X, model$means, model$covars, model$transmat, model$startprob)
}

#' Smoothed posteriors gamma_t(k). Mirrors hmmlearn's predict_proba.
predict_proba_hmm <- function(X, model) {
  storage.mode(X) <- "double"
  gaussian_hmm_predict_proba(X, model$means, model$covars, model$transmat, model$startprob)
}

# ---------------- 2. Predictive K selection ----------------

#' Predictive (out-of-sample) model-order selection.
#'
#' Splits X into X_fit = first (T - L_val) rows and X_val = last L_val rows.
#' For each K in K_candidates: fit on X_fit, score on X_val (joint log-lik,
#' which equals sum_s log p(x_s | x_{1:s-1}) by chain rule). The complexity-
#' penalized score is  predLL(K) - lambda_K * K.
#'
#' Returns the argmax-K (paper formulation: K can go up OR down within the
#' candidate set). The notebook adds a monotone-K ratchet that I intentionally
#' do not implement — see explanation in main.R.
select_K_predictive <- function(X, K_candidates, L_val = 252,
                                n_iter = 300, lambda_K = 1.0, seed = 42) {
  Tn <- nrow(X)
  if (Tn <= (L_val + 50)) return(as.integer(min(K_candidates)))
  
  X_fit <- X[seq_len(Tn - L_val), , drop = FALSE]
  X_val <- X[seq.int(Tn - L_val + 1, Tn), , drop = FALSE]
  
  best_score <- -Inf
  best_K     <- as.integer(min(K_candidates))
  
  for (K in K_candidates) {
    fit <- try(fit_hmm(X_fit, K = K, n_iter = n_iter, seed = seed), silent = TRUE)
    if (inherits(fit, "try-error")) next
    s <- try(score_hmm(X_val, fit), silent = TRUE)
    if (inherits(s, "try-error") || !is.finite(s)) next
    s_pen <- s - lambda_K * K
    if (s_pen > best_score) { best_score <- s_pen; best_K <- as.integer(K) }
  }
  best_K
}

# ---------------- 3. Wasserstein-2 wrapper ----------------

#' W2 distance between two multivariate Gaussians (closed form).
w2_gaussian <- function(mu1, S1, mu2, S2) {
  wasserstein2_gaussian(as.numeric(mu1), as.matrix(S1),
                        as.numeric(mu2), as.matrix(S2))
}

# ---------------- 4. Ledoit-Wolf wrapper ----------------

shrunk_cov <- function(X) {
  if (nrow(X) < 2) return(diag(ncol(X)))
  ledoit_wolf_cov(X)
}

# ---------------- 5. MVO with L1 turnover cost ----------------

#' Transaction-cost-aware MVO solved as a QP via OSQP.
#'
#'   maximize   mu' w - lambda w' Sigma w - tc * |w - w_prev|_1
#'   s.t.       sum(w) = 1,  0 <= w <= w_max
#'
#' Reformulated with auxiliary u = |w - w_prev|, u >= 0:
#'   variables  z = [w (N) ; u (N)]
#'   OSQP form: minimize  (1/2) z' P z + q' z   s.t.  l <= A z <= u
#'   with  P  = blockdiag(2 lambda Sigma, 0_{NxN}),
#'         q  = [-mu ; tc * 1]
#'   constraints:
#'     sum(w) = 1
#'     0 <= w <= w_max
#'     0 <= u <= +Inf
#'     w - u <=  w_prev      (i.e.  u >= w - w_prev)
#'     w + u >= w_prev       (i.e.  u >= -(w - w_prev))
solve_mvo <- function(mu, Sigma, w_prev,
                      lambda = 5.0, tc = 0.0002, w_max = 0.6) {
  if (!requireNamespace("osqp", quietly = TRUE))
    stop("Install osqp: install.packages('osqp')")
  
  N <- length(mu)
  Sigma <- 0.5 * (Sigma + t(Sigma)) + 1e-8 * diag(N)
  
  # ---- P, q
  P <- Matrix::bdiag(2 * lambda * Sigma, Matrix::Matrix(0, N, N, sparse = TRUE))
  q <- c(-as.numeric(mu), rep(tc, N))
  
  # ---- A and bounds
  I  <- Matrix::Diagonal(N)
  Z  <- Matrix::Matrix(0, N, N, sparse = TRUE)
  
  A_eq      <- cbind(matrix(1, 1, N), matrix(0, 1, N))       # sum w = 1
  A_w_box   <- cbind(I, Z)                                    # 0 <= w <= w_max
  A_u_box   <- cbind(Z, I)                                    # 0 <= u <= +Inf
  A_pos_dif <- cbind(I,  -I)                                  # w - u <= w_prev
  A_neg_dif <- cbind(I,   I)                                  # w + u >= w_prev
  
  A <- rbind(A_eq, A_w_box, A_u_box, A_pos_dif, A_neg_dif)
  
  l <- c(1,
         rep(0, N),                # w >= 0
         rep(0, N),                # u >= 0
         rep(-Inf, N),             # w - u <= w_prev
         w_prev)                   # w + u >= w_prev
  u <- c(1,
         rep(w_max, N),
         rep(Inf, N),
         w_prev,
         rep(Inf, N))
  
  settings <- osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-7, eps_rel = 1e-7,
                                 max_iter = 20000, polish = TRUE)
  model <- try(osqp::osqp(P, q, A, l, u, pars = settings), silent = TRUE)
  if (inherits(model, "try-error")) return(w_prev)
  
  sol <- model$Solve()
  if (sol$info$status_val != 1 && sol$info$status_val != 2) return(w_prev)
  
  w <- sol$x[seq_len(N)]
  w <- pmax(w, 0)
  s <- sum(w)
  if (s <= 1e-12) return(w_prev)
  w / s
}
