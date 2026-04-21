# ---------------------------------------------------------------------------
# Diagnostic checks. Run AFTER main.R has populated `state`, `dist_mat`,
# `factors`, and `strat` in the workspace (or re-source main.R first).
#
# Usage:
#   source("main.R")
#   source("R/diagnostics.R")
#   run_diagnostics(state, dist_mat, factors, strat, CFG)
# ---------------------------------------------------------------------------

run_diagnostics <- function(state, dist_mat, factors, strat, cfg) {

  cat("\n================ DIAGNOSTIC 1: transformed z-score health ================\n")
  cat("Paper Exhibit 4 reports std ~ 0.89 - 0.98 for every variable.\n")
  cat("Deviations above ~1.1 mean a variable is over-weighted in distance.\n\n")
  z <- state$transformed_winsorized
  stats <- z %>%
    pivot_longer(-date, names_to = "var", values_to = "z") %>%
    group_by(var) %>%
    summarise(
      n_valid = sum(!is.na(z)),
      mean    = mean(z, na.rm = TRUE),
      std     = sd(z, na.rm = TRUE),
      q01     = quantile(z, 0.01, na.rm = TRUE),
      q99     = quantile(z, 0.99, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(match(var, cfg$variables$var))
  print(stats %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

  cat("\n================ DIAGNOSTIC 2: distance-matrix influence =================\n")
  cat("Contribution of each variable to the *average* squared distance.\n")
  cat("If one variable is > ~25% the rest are under-weighted in the metric.\n\n")
  X <- as.matrix(z[, cfg$variables$var])
  # For each variable v, avg over all i != j of (X[i,v] - X[j,v])^2
  contrib <- sapply(seq_along(cfg$variables$var), function(k) {
    v <- X[, k]
    v <- v[!is.na(v)]
    # mean squared pairwise difference = 2 * var(v)
    2 * var(v)
  })
  names(contrib) <- cfg$variables$label
  share <- contrib / sum(contrib)
  print(tibble(variable = names(share),
               squared_contribution = round(contrib, 3),
               share = round(share, 3)) %>%
        arrange(desc(share)))

  cat("\n================ DIAGNOSTIC 3: Jan 2009 similar dates ====================\n")
  cat("Top 15 most similar historical months to Jan 2009 (paper Exhibit 6).\n")
  cat("With our sample start (~1997), we cannot reproduce the paper's 1980s\n")
  cat("hits, but 2001-02 recession or late-1998 LTCM should rank highly.\n\n")
  target <- as.Date("2009-01-31")
  i <- which(dist_mat$dates == target)
  if (length(i) == 0) {
    i <- which.min(abs(as.numeric(dist_mat$dates - target)))
    cat("(using", as.character(dist_mat$dates[i]), "as closest month)\n")
  }
  d <- dist_mat$D[i, ]
  upper <- i - cfg$mask_months
  elig <- seq_len(max(upper, 1))
  top <- tibble(date = dist_mat$dates[elig], d = d[elig]) %>%
    arrange(d) %>%
    slice_head(n = 15)
  print(top %>% mutate(d = round(d, 3)))

  cat("\n================ DIAGNOSTIC 4: signal distribution =======================\n")
  cat("For the default quintile strategy (Q1 = most similar), frequency of\n")
  cat("long (+1), flat (0), short (-1) signs per factor. In a healthy setup\n")
  cat("the +1 share should dominate (factor premia are positive on average).\n\n")
  # Reconstruct Q1 signals by rerunning the pieces we need
  bks_q <- bucket_by_similarity(dist_mat$D, cfg$mask_months, cfg$quantile_default)
  top_buckets <- lapply(bks_q, function(b) if (is.null(b)) NULL else b[[1]])
  fr <- tibble(date = dist_mat$dates) %>% left_join(factors, by = "date")
  ts_q1 <- timing_strategy(top_buckets, fr, dist_mat$dates,
                           setdiff(names(factors), "date"))
  sig_q1 <- ts_q1$signs %>%
    pivot_longer(-date, names_to = "factor", values_to = "sg") %>%
    filter(!is.na(sg)) %>%
    group_by(factor) %>%
    summarise(
      n      = n(),
      pct_long  = mean(sg > 0),
      pct_flat  = mean(sg == 0),
      pct_short = mean(sg < 0),
      .groups = "drop"
    )
  cat("Q1 (most similar):\n")
  print(sig_q1 %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

  bot_buckets <- lapply(bks_q, function(b) if (is.null(b)) NULL else b[[cfg$quantile_default]])
  ts_q5 <- timing_strategy(bot_buckets, fr, dist_mat$dates,
                           setdiff(names(factors), "date"))
  sig_q5 <- ts_q5$signs %>%
    pivot_longer(-date, names_to = "factor", values_to = "sg") %>%
    filter(!is.na(sg)) %>%
    group_by(factor) %>%
    summarise(
      n      = n(),
      pct_long  = mean(sg > 0),
      pct_flat  = mean(sg == 0),
      pct_short = mean(sg < 0),
      .groups = "drop"
    )
  cat("\nQ5 (least similar):\n")
  print(sig_q5 %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

  cat("\n================ DIAGNOSTIC 5: date alignment ===========================\n")
  cat("First and last overlap of macro dates and factor dates, and NA counts.\n\n")
  d_state  <- range(state$monthly$date)
  d_factor <- range(factors$date)
  d_common <- intersect(as.character(state$monthly$date), as.character(factors$date))
  cat("  macro monthly range :", format(d_state[1]), "to", format(d_state[2]), "\n")
  cat("  factor monthly range:", format(d_factor[1]), "to", format(d_factor[2]), "\n")
  cat("  common month-ends   :", length(d_common), "\n")
  # join quality
  fr_full <- tibble(date = state$monthly$date) %>% left_join(factors, by = "date")
  na_pct <- fr_full %>%
    summarise(across(-date, ~ mean(is.na(.x)))) %>%
    pivot_longer(everything(), names_to = "factor", values_to = "pct_NA")
  cat("\nShare of macro-date rows with NA factor return after join:\n")
  print(na_pct %>% mutate(pct_NA = round(pct_NA, 3)))

  cat("\n================ DIAGNOSTIC 6: effective strategy sample =================\n")
  strat_returns <- strat$diff_default %>% filter(is.finite(diff))
  cat("  Q1-Q5 diff strategy: ", nrow(strat_returns), "months,",
      as.character(min(strat_returns$date)), "to",
      as.character(max(strat_returns$date)), "\n\n")

  # alpha regressions
  alpha_diagnostics(strat, factors, cfg)

  invisible(list(stats = stats, top = top, sig_q1 = sig_q1, sig_q5 = sig_q5))
}


# ---------------------------------------------------------------------------
# Alpha regressions against the six Fama-French factors. Called at the end
# of run_diagnostics() but also callable standalone.
#
# For each of (Q1-Q5, Q1, Q5) we fit:
#    r_t = alpha + beta_M * Market_t + beta_S * Size_t + beta_V * Value_t
#          + beta_P * Profitability_t + beta_I * Investment_t
#          + beta_M * Momentum_t + eps_t
# and report the annualised alpha, its t-stat, regression R^2, and the Sharpe
# of the residuals (factor-neutral risk-adjusted performance).
#
# The paper reports "alpha three standard errors from zero" for the
# Q1-Q5 difference portfolio. That corresponds to the t-stat on `alpha` here.
# ---------------------------------------------------------------------------
alpha_diagnostics <- function(strat, factors, cfg) {
  cat("\n================ DIAGNOSTIC 7: factor-model alpha regressions ===========\n")
  cat("Fit strategy returns on the six FF factors; report the intercept\n")
  cat("(alpha) annualised, its t-stat, and the residual (factor-neutral)\n")
  cat("Sharpe. Paper reports Q1-Q5 alpha ~3 standard errors from zero.\n\n")

  fc <- c("Market", "Size", "Value", "Profitability",
          "Investment", "Momentum")

  run_reg <- function(df_y, label) {
    df <- df_y %>%
      inner_join(factors, by = "date") %>%
      filter(if_all(everything(), is.finite))
    if (nrow(df) < 24) {
      cat(sprintf("  %s: too few observations (n = %d), skipping\n\n",
                  label, nrow(df)))
      return(invisible(NULL))
    }
    m <- lm(reformulate(fc, response = "y"), data = df)
    s <- summary(m)

    a_month <- coef(m)[1]
    a_t     <- s$coefficients[1, "t value"]
    a_p     <- s$coefficients[1, "Pr(>|t|)"]
    resids  <- residuals(m)
    resid_sr <- mean(resids) / sd(resids) * sqrt(12)
    raw_sr   <- mean(df$y, na.rm = TRUE) / sd(df$y, na.rm = TRUE) * sqrt(12)
    raw_t    <- raw_sr * sqrt(nrow(df) / 12)

    cat(sprintf("  %s (n = %d months)\n", label, nrow(df)))
    cat(sprintf("    raw SR     = %+.2f       (naive t-stat on SR = %+.2f)\n",
                raw_sr, raw_t))
    cat(sprintf("    alpha      = %+.2f%% p.a.  t-stat = %+.2f   p = %.3f\n",
                a_month * 12 * 100, a_t, a_p))
    cat(sprintf("    R^2        = %.2f        residual SR = %+.2f\n",
                s$r.squared, resid_sr))

    betas <- s$coefficients[-1, , drop = FALSE]
    sig   <- ifelse(betas[, 4] < 0.05, " *", "  ")
    cat("    loadings:\n")
    for (i in seq_len(nrow(betas))) {
      cat(sprintf("      %-14s  beta = %+.3f  (t = %+.2f)%s\n",
                  rownames(betas)[i], betas[i, 1], betas[i, 3], sig[i]))
    }
    cat("\n")
    invisible(m)
  }

  run_reg(strat$diff_default %>% select(date, y = diff),
          "Q1 minus Q5 (difference portfolio)")
  run_reg(strat$per_bucket %>% filter(quintile == "Q1") %>%
            select(date, y = avg),
          "Q1 only (long similarity)")
  run_reg(strat$per_bucket %>%
            filter(quintile == paste0("Q", cfg$quantile_default)) %>%
            select(date, y = avg),
          "Q5 only (long dissimilarity)")

  invisible(NULL)
}
