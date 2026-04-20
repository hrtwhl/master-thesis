# 02_regimes.R -- modified k-means regime detection (Oliveira et al. 2025).
#
# Pipeline (Algorithm 1 in the paper):
#   1. k-means k=2 under L2 -> smaller cluster = Regime 0 (outliers / crisis).
#   2. k-means k=r under cosine distance on the larger cluster -> Regimes 1..r.
#   3. Per-month regime probability via fuzzy distances to centroids, combined
#      across the two stages via eq. (4).
#
# Also: transition matrix (eq. 5), Hungarian regime-label matching for
# walk-forward consistency, GMM baseline.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---- distance helpers --------------------------------------------------------

l2_dist_to_centroids <- function(X, centroids) {
  n <- nrow(X); k <- nrow(centroids)
  d <- matrix(NA_real_, n, k)
  for (j in seq_len(k)) {
    d[, j] <- sqrt(rowSums(sweep(X, 2, centroids[j, ], "-")^2))
  }
  d
}

cosine_dist_to_centroids <- function(X, centroids) {
  Xn <- X / pmax(sqrt(rowSums(X^2)), 1e-12)
  Cn <- centroids / pmax(sqrt(rowSums(centroids^2)), 1e-12)
  1 - Xn %*% t(Cn)
}

# ---- fuzzy probabilities (eq. 1) ---------------------------------------------

fuzzy_probs <- function(dist_mat) {
  dist_mat <- pmax(dist_mat, 1e-12)
  row_sum  <- rowSums(dist_mat)
  w <- 1 - dist_mat / row_sum
  w <- pmax(w, 0)
  w / pmax(rowSums(w), 1e-12)
}

# ---- cosine k-means (Lloyd on unit sphere) ----------------------------------

kmeans_cosine <- function(X, k, n_init = 20, max_iter = 100, seed = 1L) {
  set.seed(seed)
  Xn <- X / pmax(sqrt(rowSums(X^2)), 1e-12)

  best <- NULL
  for (trial in seq_len(n_init)) {
    init_ix <- sample.int(nrow(Xn), k)
    cent <- Xn[init_ix, , drop = FALSE]
    for (it in seq_len(max_iter)) {
      sim <- Xn %*% t(cent)
      lab <- max.col(sim, ties.method = "first")
      new_cent <- t(vapply(seq_len(k), function(j) {
        idx <- which(lab == j)
        if (length(idx) == 0) return(cent[j, ])
        m <- colMeans(Xn[idx, , drop = FALSE])
        m / pmax(sqrt(sum(m^2)), 1e-12)
      }, numeric(ncol(Xn))))
      if (max(abs(new_cent - cent)) < 1e-8) { cent <- new_cent; break }
      cent <- new_cent
    }
    inertia <- sum(1 - rowSums(Xn * cent[lab, , drop = FALSE]))
    if (is.null(best) || inertia < best$inertia) {
      best <- list(centroids = cent, labels = lab, inertia = inertia)
    }
  }
  best
}

kmeans_elbow_cosine <- function(X, k_min = 2, k_max = 10, ...) {
  inertias <- vapply(k_min:k_max, function(k) kmeans_cosine(X, k, ...)$inertia,
                     numeric(1))
  if (length(inertias) < 3) return(k_min)
  d2 <- diff(inertias, differences = 2)
  (k_min:k_max)[which.max(-d2) + 1]
}

# ---- main regime fitter ------------------------------------------------------

fit_regime_model <- function(X, k_inner = 5, seed = 1L) {
  set.seed(seed)

  # stage 1: L2 k=2 outlier split
  km2 <- stats::kmeans(X, centers = 2, nstart = 25, iter.max = 100)
  sizes <- km2$size
  outlier_id <- which.min(sizes)
  typical_id <- setdiff(1:2, outlier_id)
  outer_centroids <- km2$centers[c(outlier_id, typical_id), , drop = FALSE]
  rownames(outer_centroids) <- c("outlier", "typical")

  is_typical <- km2$cluster == typical_id

  # stage 2: cosine k-means on typical months
  X_typ <- X[is_typical, , drop = FALSE]
  km_inner <- kmeans_cosine(X_typ, k_inner, seed = seed)

  labels <- integer(nrow(X))
  labels[!is_typical] <- 0L
  labels[is_typical]  <- km_inner$labels

  list(
    outer_centroids = outer_centroids,
    inner_centroids = km_inner$centroids,
    k_inner = k_inner,
    labels  = labels,
    dates   = as.Date(rownames(X))
  )
}

# ---- probability assignment for new samples ---------------------------------

classify_regimes <- function(model, X_new) {
  d_outer <- l2_dist_to_centroids(X_new, model$outer_centroids)
  p_outer <- fuzzy_probs(d_outer)
  p_reg0  <- p_outer[, 1]

  d_inner <- cosine_dist_to_centroids(X_new, model$inner_centroids)
  p_inner <- fuzzy_probs(d_inner)

  # eq. (4): scale P(Regime 0) by -P_max * log2(1 - P(Regime 0))
  p_max <- apply(p_inner, 1, max)
  safe <- pmin(p_reg0, 1 - 1e-9)
  p_reg0_scaled <- -p_max * log2(1 - safe)
  p_reg0_scaled[p_reg0 == 0] <- 0

  P <- cbind(p_reg0_scaled, p_inner)
  colnames(P) <- c("R0", paste0("R", seq_len(ncol(p_inner))))
  P <- P / pmax(rowSums(P), 1e-12)
  rownames(P) <- rownames(X_new)

  hard <- apply(P, 1, which.max) - 1L
  list(probs = P, labels = hard)
}

# ---- transition matrix (eq. 5) ----------------------------------------------

regime_transition_matrix <- function(labels, n_regimes) {
  states <- 0:(n_regimes - 1L)
  E <- matrix(0, nrow = n_regimes, ncol = n_regimes,
              dimnames = list(from = states, to = states))
  for (t in seq_len(length(labels) - 1L)) {
    E[as.character(labels[t]), as.character(labels[t + 1L])] <-
      E[as.character(labels[t]), as.character(labels[t + 1L])] + 1
  }
  row_tot <- rowSums(E)
  row_tot[row_tot == 0] <- 1
  E / row_tot
}

conditional_transition_matrix <- function(E) {
  diag_e <- diag(E)
  off <- E; diag(off) <- 0
  denom <- 1 - diag_e
  denom[denom == 0] <- 1
  sweep(off, 1, denom, "/")
}

# ---- regime label matching for walk-forward ---------------------------------
#
# Cleaner rewrite: use Hungarian algorithm on cosine-distance cost matrix
# between new centroids (rows of A) and reference centroids (rows of B).
# assignment[i] = j means "new cluster i corresponds to reference cluster j".
# We remap so that centroids come out in reference order and labels get
# translated through the assignment table.

match_regime_labels <- function(new_model, ref_model) {
  if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")

  A <- new_model$inner_centroids     # new (current refit)
  B <- ref_model$inner_centroids     # canonical reference
  k <- nrow(A)

  An <- A / pmax(sqrt(rowSums(A^2)), 1e-12)
  Bn <- B / pmax(sqrt(rowSums(B^2)), 1e-12)
  cost <- 1 - An %*% t(Bn)            # k x k cosine-distance matrix
  assignment <- as.integer(clue::solve_LSAP(cost, maximum = FALSE))

  # Reorder new_model centroids so row j matches reference row j:
  # want centroid_out[j] = A[i] where assignment[i] == j
  perm_inv <- integer(k)
  for (i in seq_len(k)) perm_inv[assignment[i]] <- i
  new_model$inner_centroids <- A[perm_inv, , drop = FALSE]

  # Translate labels: each old inner label i (1..k) becomes assignment[i].
  # Regime 0 (outlier) stays 0.
  lbl <- new_model$labels
  new_lbl <- lbl
  inner_mask <- lbl >= 1L & lbl <= k
  new_lbl[inner_mask] <- assignment[lbl[inner_mask]]
  new_model$labels <- new_lbl
  new_model
}

# ---- forecasted next-regime distribution (eq. 7) ----------------------------

next_regime_probs <- function(p_t, transition_matrix) {
  p_t <- p_t / sum(p_t)
  as.numeric(p_t %*% transition_matrix)
}

# ---- GMM baseline for comparison --------------------------------------------

fit_gmm_regimes <- function(X, G = 6, seed = 1L,
                            model_names = c("VVI", "EEE", "VVV")) {
  if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
  suppressPackageStartupMessages(library(mclust))
  set.seed(seed)

  fit <- NULL
  for (mn in model_names) {
    fit_try <- tryCatch(
      mclust::Mclust(X, G = G, modelNames = mn, verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit_try) && !is.null(fit_try$classification) &&
        length(fit_try$classification) == nrow(X)) {
      fit <- fit_try
      if (mn != model_names[1]) {
        message(sprintf("  GMM: '%s' failed, using '%s' instead.",
                        model_names[1], mn))
      }
      break
    }
  }
  if (is.null(fit)) {
    warning("Mclust failed for all candidate models; returning NA labels.")
    return(list(
      labels = rep(NA_integer_, nrow(X)),
      probs  = matrix(NA_real_, nrow(X), G),
      fit    = NULL,
      dates  = as.Date(rownames(X))
    ))
  }

  list(
    labels = fit$classification - 1L,
    probs  = fit$z,
    fit    = fit,
    dates  = as.Date(rownames(X))
  )
}
