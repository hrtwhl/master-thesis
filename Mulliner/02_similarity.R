# ---- Similarity metric ----------------------------------------------------

#' Compute the full NxN Euclidean distance matrix across the transformed
#' state variables. Entry [i, j] is the sqrt of the sum of squared
#' differences across the V variables between month i and month j.
#' Rows with any missing variable are NA (we drop them when selecting).
compute_distance_matrix <- function(transformed_w) {
  dates <- transformed_w$date
  X <- as.matrix(transformed_w[ , -1])
  n <- nrow(X)

  # Use the standard fast route: ||a-b||^2 = ||a||^2 + ||b||^2 - 2 a.b
  sq <- rowSums(X * X)
  cross <- X %*% t(X)
  D2 <- outer(sq, sq, "+") - 2 * cross
  D2[D2 < 0] <- 0
  D <- sqrt(D2)

  rownames(D) <- as.character(dates)
  colnames(D) <- as.character(dates)

  # propagate NAs where any variable was missing
  row_has_na <- apply(X, 1, function(r) any(is.na(r)))
  D[row_has_na, ] <- NA
  D[ , row_has_na] <- NA
  diag(D) <- 0

  list(D = D, dates = dates, X = X,
       valid = !row_has_na)
}

#' For each target date T, return the indices of the historical months that
#' are the most similar (lowest distance), excluding the last `mask_months`
#' before T.
#' Returns a list indexed by target date, each element an integer vector of
#' row indices in `D`.
select_similar <- function(D, dates, mask_months, top_pct = 0.15,
                           from_date = NULL) {
  n <- nrow(D)
  out <- vector("list", n)
  names(out) <- as.character(dates)
  for (i in seq_len(n)) {
    if (all(is.na(D[i, ]))) next
    # eligible historical months: j < i - mask_months and valid
    upper <- i - mask_months
    if (upper < 2) next
    d <- D[i, 1:upper]
    valid_idx <- which(!is.na(d))
    if (length(valid_idx) < 5) next
    k <- max(1L, floor(length(valid_idx) * top_pct))
    ord <- valid_idx[order(d[valid_idx])]
    out[[i]] <- head(ord, k)
  }
  out
}

#' EWMA of the distance-from-T series, used to flag regime shifts.
#' At each time T, the weights favour t close to T (recent past) so a large
#' C_T means recent months look unlike T (the environment is changing).
#' beta = 1 - 1/n, with n in months.
compute_ewma_scores <- function(D, dates, lookback_months) {
  n <- nrow(D)
  beta <- 1 - 1 / lookback_months
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (all(is.na(D[i, 1:i]))) next
    if (i < 2) next
    j <- seq_len(i - 1)
    d <- D[i, j]
    keep <- !is.na(d)
    if (!any(keep)) next
    # weights: most recent (largest j) gets weight ~1, decay going back
    w <- beta ^ (i - j)
    out[i] <- sum(w[keep] * d[keep]) / sum(w[keep])
  }
  tibble(date = dates, ewma = out)
}

#' Convenience: EWMA half-life given beta = 1 - 1/n.
ewma_half_life <- function(n) -log(2) / log(1 - 1 / n)
