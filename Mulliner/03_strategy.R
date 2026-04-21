# ---- Factor-timing strategies ---------------------------------------------

#' Given a distance matrix, for each date T split historical months (eligible
#' months before T - mask) into `q` equal-sized buckets ranked by similarity.
#' Bucket 1 is the most similar; bucket q is the most dissimilar.
#' Returns an integer matrix of row indices per bucket for each T (as a list).
bucket_by_similarity <- function(D, mask_months, q) {
  n <- nrow(D)
  buckets <- vector("list", n)
  for (i in seq_len(n)) {
    upper <- i - mask_months
    if (upper < 5) next
    d <- D[i, 1:upper]
    valid_idx <- which(!is.na(d))
    if (length(valid_idx) < q * 2) next
    ranks <- rank(d[valid_idx], ties.method = "first")
    b_assign <- ceiling(ranks / length(valid_idx) * q)
    # positional list so buckets[[i]][[k]] is always bucket k
    by_bucket <- vector("list", q)
    for (k in seq_len(q)) by_bucket[[k]] <- valid_idx[b_assign == k]
    buckets[[i]] <- by_bucket
  }
  buckets
}

#' Given a bucket membership list, compute the factor-timing strategy return
#' for each month T:
#'   - signal is computed from the state at T-1 (buckets[[T-1]]) so that the
#'     trade placed at T uses only information through T-1
#'   - for each factor f, look up returns in month (t+1) for every t in the
#'     selected bucket (t is strictly older than T-1 - mask_months)
#'   - sign = +1 if the average is > 0, else -1
#'   - strategy return at T = factor return at T multiplied by sign, averaged
#'     across factors
timing_strategy <- function(buckets, factor_mat, dates, factor_names) {
  n <- length(dates)
  R <- as.matrix(factor_mat[, factor_names, drop = FALSE])  # n x F
  strat <- matrix(NA_real_, n, length(factor_names),
                  dimnames = list(NULL, factor_names))
  signs <- matrix(NA_real_, n, length(factor_names),
                  dimnames = list(NULL, factor_names))

  for (i in seq_len(n - 1)) {
    b <- buckets[[i]]   # signal from state at month i
    if (is.null(b) || length(b) == 0) next
    if (is.list(b)) stop("timing_strategy expects a single bucket per date")
    idx_next <- b + 1
    idx_next <- idx_next[idx_next <= n]
    if (length(idx_next) == 0) next
    sub <- R[idx_next, , drop = FALSE]
    avg <- colMeans(sub, na.rm = TRUE)
    sg  <- sign(avg)
    sg[!is.finite(sg) | sg == 0] <- 0   # undefined -> flat
    # realised at month i+1
    strat[i + 1, ] <- sg * R[i + 1, ]
    signs[i + 1, ] <- sg
  }

  list(
    returns = bind_cols(tibble(date = dates), as_tibble(strat)),
    signs   = bind_cols(tibble(date = dates), as_tibble(signs))
  )
}

#' Equally-weighted combination of per-factor timed returns.
equal_weight <- function(returns_tbl, factor_names) {
  returns_tbl %>%
    mutate(avg = rowMeans(across(all_of(factor_names)), na.rm = TRUE)) %>%
    select(date, avg)
}

#' Vol-target a monthly return series to `target_ann` annual vol, using a
#' trailing window of realised vol (ex-ante).
vol_target <- function(r_tbl, target_ann = 0.15, window = 36) {
  r <- r_tbl$avg
  rv <- slider::slide_dbl(r, ~ sd(.x, na.rm = TRUE),
                          .before = window, .after = -1,
                          .complete = FALSE)
  scale <- (target_ann / sqrt(12)) / rv
  scale[!is.finite(scale)] <- NA
  tibble(date = r_tbl$date, r = r, r_vt = r * scale, scale = scale)
}

#' Run all strategies (quintiles, difference portfolio, robustness sweeps).
run_strategies <- function(distance_matrix, factor_returns, mask_months,
                           q_default, q_robust, lookback_robust,
                           macro_monthly, variables, horizon, winsorize) {
  D <- distance_matrix$D
  dates_d <- distance_matrix$dates

  # Align factor returns to the macro-state dates
  fr <- tibble(date = dates_d) %>%
    left_join(factor_returns, by = "date")
  factor_names <- setdiff(names(factor_returns), "date")

  # ---- Default (quintile) run ---------------------------------------------
  bks_q <- bucket_by_similarity(D, mask_months, q_default)

  per_bucket <- vector("list", q_default)
  names(per_bucket) <- paste0("Q", seq_len(q_default))
  for (q in seq_len(q_default)) {
    this_buckets <- lapply(bks_q, function(b) if (is.null(b)) NULL else b[[q]])
    per_bucket[[q]] <- timing_strategy(this_buckets, fr, dates_d, factor_names)
  }

  # Equally-weighted across factors for each quintile
  ew_q <- map_dfr(seq_len(q_default), function(q) {
    equal_weight(per_bucket[[q]]$returns, factor_names) %>%
      mutate(quintile = paste0("Q", q))
  })

  # Long-only baseline = average across factors
  lo <- fr %>%
    mutate(avg = rowMeans(across(all_of(factor_names)), na.rm = TRUE)) %>%
    select(date, avg)

  # Difference portfolio: Q1 - Q5
  diff_qs <- ew_q %>%
    filter(quintile %in% c("Q1", paste0("Q", q_default))) %>%
    pivot_wider(names_from = quintile, values_from = avg) %>%
    mutate(diff = .data[["Q1"]] - .data[[paste0("Q", q_default)]]) %>%
    select(date, diff)

  # Vol-targeted version of the difference portfolio (for Exhibit 1)
  diff_vt <- vol_target(diff_qs %>% rename(avg = diff),
                        target_ann = CFG$vol_target,
                        window     = CFG$vol_window)

  # ---- Robustness: different quantile choices -----------------------------
  rb_quantile <- map_dfr(q_robust, function(q) {
    bks <- bucket_by_similarity(D, mask_months, q)
    top_buckets <- lapply(bks, function(b) if (is.null(b)) NULL else b[[1]])
    bot_buckets <- lapply(bks, function(b) if (is.null(b)) NULL else b[[q]])
    top <- equal_weight(timing_strategy(top_buckets, fr, dates_d, factor_names)$returns,
                        factor_names)
    bot <- equal_weight(timing_strategy(bot_buckets, fr, dates_d, factor_names)$returns,
                        factor_names)
    tibble(date    = top$date,
           top_ret = top$avg,
           bot_ret = bot$avg,
           diff    = top$avg - bot$avg,
           q       = q)
  })

  # ---- Robustness: different z-score lookbacks ----------------------------
  rb_lookback <- map_dfr(lookback_robust, function(lk_months) {
    state_alt <- build_state_variables(macro_monthly, variables,
                                       horizon   = horizon,
                                       window    = lk_months,
                                       winsorize = winsorize)
    dm_alt <- compute_distance_matrix(state_alt$transformed_winsorized)
    bks <- bucket_by_similarity(dm_alt$D, mask_months, q_default)
    top_buckets <- lapply(bks, function(b) if (is.null(b)) NULL else b[[1]])
    bot_buckets <- lapply(bks, function(b) if (is.null(b)) NULL else b[[q_default]])
    fr_alt <- tibble(date = dm_alt$dates) %>%
      left_join(factor_returns, by = "date")
    top <- equal_weight(timing_strategy(top_buckets, fr_alt, dm_alt$dates, factor_names)$returns,
                        factor_names)
    bot <- equal_weight(timing_strategy(bot_buckets, fr_alt, dm_alt$dates, factor_names)$returns,
                        factor_names)
    tibble(date = top$date, diff = top$avg - bot$avg,
           lookback_years = lk_months / 12)
  })

  # ---- Individual factor analysis (Exhibits A1, A2) -----------------------
  per_factor <- map_dfr(seq_len(q_default), function(q) {
    this_buckets <- lapply(bks_q, function(b) if (is.null(b)) NULL else b[[q]])
    ts <- timing_strategy(this_buckets, fr, dates_d, factor_names)
    ts$returns %>%
      pivot_longer(-date, names_to = "factor", values_to = "r") %>%
      mutate(quintile = paste0("Q", q))
  })

  list(
    factor_names  = factor_names,
    per_bucket    = ew_q,                   # monthly returns per quintile
    diff_default  = diff_qs,                # Q1 - Q5
    diff_vol_targeted = diff_vt,            # for Exhibit 1
    long_only     = lo,
    robust_qtl    = rb_quantile,
    robust_lkback = rb_lookback,
    per_factor    = per_factor
  )
}
