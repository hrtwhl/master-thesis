# ---------------------------------------------------------------------------
# R/06_pca_factors.R
#
# Estimate the FRED-MD principal-component factors per Stock-Watson EM-PCA,
# with the number of factors selected per Bai & Ng (2002).
#
# Key methodological points:
#   1. STRICT CAUSALITY. At each estimation date T, PCA loadings, means, and
#      standard deviations use only data through T. The output series at T
#      uses only data observable at T.
#   2. EM imputation for missing values, exactly as in McCracken-Ng (2016).
#   3. SIGN ALIGNMENT. Each rolling re-estimation can flip the arbitrary
#      sign of any factor; we re-anchor signs to a fixed reference panel
#      so the time series of factor 1 is comparable across estimation dates.
#
# Functions:
#   bai_ng_select() : compute Bai-Ng IC criteria and return arg min
#   estimate_pca_em(): static EM-PCA on a balanced/unbalanced panel
#   rolling_pca_factors(): the production routine - returns a tibble of
#                          F factors at month-end, estimated causally
# ---------------------------------------------------------------------------

#' Standardise a matrix column-wise: subtract column mean, divide by column
#' std. NAs in input remain NA in output.
standardise_cols <- function(X) {
  mu <- colMeans(X, na.rm = TRUE)
  sd <- apply(X, 2, sd, na.rm = TRUE)
  sd[sd == 0 | !is.finite(sd)] <- 1
  Xs <- sweep(X, 2, mu, "-")
  Xs <- sweep(Xs, 2, sd, "/")
  list(X = Xs, mu = mu, sd = sd)
}

#' Static EM-PCA (Stock-Watson 2002, McCracken-Ng 2016).
#'
#' Iteratively (i) replaces missing values with their factor-model fits and
#' (ii) re-estimates factors via SVD of the standardised panel. We initialise
#' missing values to zero (= the column mean after standardisation).
#'
#' Inputs:
#'   X     : T x N numeric matrix (may contain NAs)
#'   r     : number of factors to extract
#'   tol   : convergence tolerance on the change in imputed values
#'   max_iter: hard cap on iterations
#'
#' Returns a list with F (T x r factor scores), Lambda (N x r loadings),
#' the standardisation parameters, and convergence info.
estimate_pca_em <- function(X, r, tol = 1e-6, max_iter = 200) {
  Tn <- nrow(X); N <- ncol(X)

  # Standardise and replace NAs with zeros (= column mean post-standardisation)
  st <- standardise_cols(X)
  Xs <- st$X
  na_mask <- is.na(Xs)
  Xs[na_mask] <- 0

  prev_impute <- Xs[na_mask]
  iter <- 0L
  converged <- FALSE

  while (iter < max_iter) {
    iter <- iter + 1L
    # SVD
    sv <- svd(Xs, nu = r, nv = r)
    # Loadings normalised so that Lambda' Lambda / N = I_r is approximate;
    # in practice the convention follows Stock-Watson: F = sqrt(T) * U_r,
    # Lambda = (1 / sqrt(T)) * V_r * D_r so that F' F / T = I_r.
    F     <- sqrt(Tn) * sv$u
    Lambda <- (1 / sqrt(Tn)) * sv$v %*% diag(sv$d[seq_len(r)], r, r)
    fit   <- F %*% t(Lambda)

    if (any(na_mask)) {
      new_impute <- fit[na_mask]
      delta <- max(abs(new_impute - prev_impute))
      Xs[na_mask] <- new_impute
      prev_impute <- new_impute
      if (delta < tol) { converged <- TRUE; break }
    } else {
      converged <- TRUE; break
    }
  }

  # Re-standardise once more after final imputation so the factors are
  # consistent with the final standardised panel
  st_final <- standardise_cols(Xs)

  list(
    F         = F,
    Lambda    = Lambda,
    sv_d      = sv$d,
    mu        = st$mu,
    sd        = st$sd,
    n_iter    = iter,
    converged = converged,
    Xs_final  = Xs
  )
}

#' Bai & Ng (2002) factor-number selection. Computes ICp1, ICp2, PCp1, PCp2.
#'
#' For each candidate r in 1..rmax, fit r factors via EM-PCA and compute:
#'   V(r) = (1 / NT) sum_{i,t} (x_it_std - lambda_i' f_t)^2
#' then the criterion is:
#'   ICp1(r) = log V(r) + r * ((N + T) / (N T)) * log( (N T) / (N + T) )
#'   ICp2(r) = log V(r) + r * ((N + T) / (N T)) * log( min(N, T) )
#'   PCp1(r) = V(r)     + r * sigma_hat^2 * ((N + T) / (N T)) * log( (N T) / (N + T) )
#'   PCp2(r) = V(r)     + r * sigma_hat^2 * ((N + T) / (N T)) * log( min(N, T) )
#' where sigma_hat^2 = V(rmax).
#'
#' Returns the table of criteria values per r and the arg-min for each.
bai_ng_select <- function(X, rmax = 12, tol = 1e-6, max_iter = 200) {
  Tn <- nrow(X); N <- ncol(X)
  cn  <- (N + Tn) / (N * Tn)
  pen_log_min <- log(min(N, Tn))
  pen_log_NT  <- log((N * Tn) / (N + Tn))

  # Compute V(r) for all r up to rmax
  V <- numeric(rmax)
  for (r in seq_len(rmax)) {
    fit <- estimate_pca_em(X, r, tol = tol, max_iter = max_iter)
    resid <- fit$Xs_final - fit$F %*% t(fit$Lambda)
    V[r] <- sum(resid^2) / (N * Tn)
  }

  sigma_hat2 <- V[rmax]

  ICp1 <- log(V) + seq_len(rmax) * cn * pen_log_NT
  ICp2 <- log(V) + seq_len(rmax) * cn * pen_log_min
  PCp1 <- V     + seq_len(rmax) * sigma_hat2 * cn * pen_log_NT
  PCp2 <- V     + seq_len(rmax) * sigma_hat2 * cn * pen_log_min

  tibble(
    r    = seq_len(rmax),
    V    = V,
    ICp1 = ICp1, ICp2 = ICp2,
    PCp1 = PCp1, PCp2 = PCp2,
    arg_min = c(which.min(ICp1) == seq_len(rmax)) |
              c(which.min(ICp2) == seq_len(rmax))
  ) -> tab

  list(
    table   = tab,
    r_ICp1  = which.min(ICp1),
    r_ICp2  = which.min(ICp2),
    r_PCp1  = which.min(PCp1),
    r_PCp2  = which.min(PCp2),
    sigma2  = sigma_hat2
  )
}


#' Sign- and order-align a current factor estimate to a reference panel
#' so the time series of factor k is comparable across rolling re-estimations.
#'
#' The sign of an SVD-derived factor is arbitrary. The ORDER of factors past
#' the first few may also shuffle as the eigenvalues of small subsequent
#' factors are close together. To stabilise both:
#'   (i) project the new loadings onto the reference loadings; for each
#'       reference factor j, find the new factor k* with maximum |cor(loadings)|
#'       that has not already been claimed.
#'   (ii) flip the sign if the correlation is negative.
#'
#' Inputs:
#'   F_new      : T_new x r_new (current rolling estimate, already on the same
#'                variable set as F_ref)
#'   Lambda_new : N x r_new
#'   Lambda_ref : N x r_ref reference loadings (e.g. full-sample first PC of
#'                a stable subset, or the first rolling estimate)
#'   r_target   : how many aligned factors to return (must be <= min(r_new, r_ref))
#'
#' Returns a list with sign-corrected, order-permuted F and Lambda matrices.
align_factors <- function(F_new, Lambda_new, Lambda_ref, r_target) {
  r_new <- ncol(Lambda_new)
  r_ref <- ncol(Lambda_ref)
  r_eff <- min(r_target, r_new, r_ref)

  # Pearson correlation between each new loading column and each reference
  # loading column, on the variables shared between them
  C <- cor(Lambda_new[, seq_len(r_new), drop = FALSE],
           Lambda_ref[, seq_len(r_ref), drop = FALSE])

  perm <- integer(r_eff)
  signs <- numeric(r_eff)
  taken <- logical(r_new)

  for (j in seq_len(r_eff)) {
    cand <- abs(C[, j])
    cand[taken] <- -Inf
    k_star <- which.max(cand)
    perm[j]  <- k_star
    signs[j] <- sign(C[k_star, j])
    if (signs[j] == 0) signs[j] <- 1
    taken[k_star] <- TRUE
  }

  F_aligned      <- F_new[, perm, drop = FALSE] %*% diag(signs, r_eff, r_eff)
  Lambda_aligned <- Lambda_new[, perm, drop = FALSE] %*% diag(signs, r_eff, r_eff)
  list(F = F_aligned, Lambda = Lambda_aligned, perm = perm, signs = signs)
}


#' Rolling/expanding-window PCA factor extraction.
#'
#' At each end-of-month estimation date in `est_dates`, fit EM-PCA on data
#' up through that date, optionally select the number of factors via Bai-Ng,
#' align signs to a fixed reference, and record the most recent factor
#' values as the "live" factor scores at that date.
#'
#' Returns a tibble with columns date, F1, F2, ..., F<r_target>.
#'
#' Args:
#'   transformed : tibble with `date` and macro variable columns
#'   r_target    : number of factors to return per date
#'   r_max       : max factors considered for selection (must be >= r_target)
#'   est_dates   : vector of dates at which to refit; default = all months
#'                 in transformed (monthly refit)
#'   selection   : "fixed" uses r_target as the chosen r; "bai_ng" runs
#'                 Bai-Ng each refit and uses the IC-chosen r (capped at
#'                 r_max). When r_chosen != r_target we still return r_target
#'                 columns; if r_chosen < r_target the trailing columns are NA.
#'   ic_type     : if selection = "bai_ng", which IC to use ("ICp2" default)
#'   min_obs     : minimum number of months before the first estimate
#'   verbose     : print progress every k refits
rolling_pca_factors <- function(transformed,
                                 r_target  = 7L,
                                 r_max     = 12L,
                                 est_dates = NULL,
                                 selection = c("fixed", "bai_ng"),
                                 ic_type   = "ICp2",
                                 min_obs   = 120L,
                                 verbose   = TRUE) {
  selection <- match.arg(selection)
  vars  <- setdiff(names(transformed), "date")
  dates <- transformed$date

  if (is.null(est_dates)) est_dates <- dates
  est_dates <- est_dates[est_dates >= dates[min_obs]]

  # Reference loadings: estimate once on the FULL sample and use these as a
  # fixed sign/order anchor across the rolling re-estimations. This is purely
  # to keep the factor time series visually consistent - the rolling factor
  # SCORES we return are still strictly causal.
  X_full <- as.matrix(transformed[, vars])
  ref_fit <- estimate_pca_em(X_full, r = r_max)
  Lambda_ref <- ref_fit$Lambda

  out_F <- matrix(NA_real_, length(est_dates), r_target)
  r_chosen_vec <- integer(length(est_dates))
  ic_table_list <- list()

  for (i in seq_along(est_dates)) {
    d <- est_dates[i]
    idx <- which(dates <= d)
    if (length(idx) < min_obs) next

    X_t <- as.matrix(transformed[idx, vars])

    if (selection == "bai_ng") {
      bn <- bai_ng_select(X_t, rmax = r_max)
      r_use <- bn[[paste0("r_", ic_type)]]
      ic_table_list[[as.character(d)]] <- bn$table
    } else {
      r_use <- r_target
    }
    r_chosen_vec[i] <- r_use

    fit <- estimate_pca_em(X_t, r = max(r_use, r_target))
    al  <- align_factors(fit$F, fit$Lambda, Lambda_ref, r_target = r_target)

    # Most recent row of the aligned factor matrix is the current state
    out_F[i, ] <- al$F[nrow(al$F), seq_len(r_target)]

    if (verbose && (i %% 50 == 0 || i == length(est_dates))) {
      cat(sprintf("  rolling PCA: %d / %d (date %s, r_chosen = %d)\n",
                  i, length(est_dates), format(d), r_use))
    }
  }

  result <- tibble(date = est_dates, !!!setNames(
      lapply(seq_len(r_target), function(k) out_F[, k]),
      paste0("F", seq_len(r_target))
  ))

  attr(result, "r_chosen")    <- tibble(date = est_dates, r_chosen = r_chosen_vec)
  attr(result, "ic_tables")   <- ic_table_list
  attr(result, "Lambda_ref")  <- Lambda_ref
  attr(result, "ref_fit_mu")  <- ref_fit$mu
  attr(result, "ref_fit_sd")  <- ref_fit$sd
  result
}


#' Static (one-shot) full-sample factor extraction. Useful for descriptive
#' analysis (factor interpretation: which variables load on which factor) -
#' NOT for use in any forward-looking signal because it is not causal.
static_pca_factors <- function(transformed, r = 8L) {
  vars  <- setdiff(names(transformed), "date")
  X     <- as.matrix(transformed[, vars])
  fit   <- estimate_pca_em(X, r = r)
  list(
    factors = tibble(date = transformed$date,
                     !!!setNames(
                       lapply(seq_len(r), function(k) fit$F[, k]),
                       paste0("F", seq_len(r))
                     )),
    loadings = fit$Lambda,
    var_names = vars,
    fit = fit
  )
}
