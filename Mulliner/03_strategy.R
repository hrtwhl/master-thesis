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
#' Accepts either a tibble with a `r` / `avg` / `diff` column (auto-detected)
#' or a raw numeric vector. Returns a tibble with date, r (raw), r_vt
#' (vol-targeted), and scale (the leverage applied at each month).
vol_target <- function(r_tbl, target_ann = 0.15, window = 36) {
  if (is.numeric(r_tbl)) {
    r <- r_tbl
    d <- NULL
  } else {
    # auto-detect the return column
    rc <- intersect(c("r", "avg", "diff", "ret"), names(r_tbl))[1]
    if (is.na(rc)) stop("vol_target: no return column found in tibble")
    r <- r_tbl[[rc]]
    d <- r_tbl$date
  }
  rv <- slider::slide_dbl(r, ~ sd(.x, na.rm = TRUE),
                          .before = window, .after = -1,
                          .complete = FALSE)
  scale <- (target_ann / sqrt(12)) / rv
  scale[!is.finite(scale)] <- NA
  if (is.null(d)) {
    tibble(r = r, r_vt = r * scale, scale = scale)
  } else {
    tibble(date = d, r = r, r_vt = r * scale, scale = scale)
  }
}

#' Apply vol targeting to a long-format per-quintile tibble with columns
#' (date, avg, quintile). Each quintile is scaled using its own trailing
#' realised vol, so all quintiles end up at ~target_ann vol individually.
vol_target_by_group <- function(tbl, group_col, ret_col = "avg",
                                target_ann = 0.15, window = 36) {
  tbl %>%
    group_by(.data[[group_col]]) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      rv    = slider::slide_dbl(.data[[ret_col]], ~ sd(.x, na.rm = TRUE),
                                .before = window, .after = -1,
                                .complete = FALSE),
      scale = (target_ann / sqrt(12)) / rv,
      scale = if_else(is.finite(scale), scale, NA_real_),
      avg_vt = .data[[ret_col]] * scale
    ) %>%
    ungroup() %>%
    select(-rv)
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

  # ---- Vol-targeted variants (paper convention, ~15% ann vol) -------------
  # Apply a 15% portfolio vol target to each series individually using
  # trailing 36m realised vol. This matches the paper's Exhibit-11 scale
  # (LO reaches ~50% drawdowns only after ~3-4x leverage).
  ew_q_vt <- vol_target_by_group(ew_q, "quintile",
                                 target_ann = CFG$vol_target,
                                 window     = CFG$vol_window) %>%
    select(date, quintile, avg_vt) %>%
    rename(avg = avg_vt)

  lo_vt_full <- vol_target(lo %>% rename(r = avg),
                           target_ann = CFG$vol_target,
                           window     = CFG$vol_window)
  lo_vt <- tibble(date = lo_vt_full$date, avg = lo_vt_full$r_vt)

  diff_vt_full <- vol_target(diff_qs %>% rename(r = diff),
                             target_ann = CFG$vol_target,
                             window     = CFG$vol_window)
  diff_vt <- tibble(date = diff_vt_full$date, diff = diff_vt_full$r_vt)

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

  # Vol-target the per-factor series too (grouped by factor x quintile).
  per_factor_vt <- per_factor %>%
    group_by(factor, quintile) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      rv    = slider::slide_dbl(r, ~ sd(.x, na.rm = TRUE),
                                .before = CFG$vol_window, .after = -1,
                                .complete = FALSE),
      scale = (CFG$vol_target / sqrt(12)) / rv,
      scale = if_else(is.finite(scale), scale, NA_real_),
      r     = r * scale
    ) %>%
    ungroup() %>%
    select(date, factor, r, quintile)

  list(
    factor_names      = factor_names,
    per_bucket        = ew_q,                   # monthly returns per quintile (unlevered)
    per_bucket_vt     = ew_q_vt,                # vol-targeted version
    diff_default      = diff_qs,                # Q1 - Q5 (unlevered)
    diff_default_vt   = diff_vt,                # Q1 - Q5 (vol-targeted)
    diff_vol_targeted = diff_vt_full,           # kept for Exhibit 1 (legacy interface)
    long_only         = lo,                     # unlevered
    long_only_vt      = lo_vt,                  # vol-targeted (paper scale)
    robust_qtl        = rb_quantile,
    robust_lkback     = rb_lookback,
    per_factor        = per_factor,             # unlevered per-factor x quintile
    per_factor_vt     = per_factor_vt           # vol-targeted
  )
}
