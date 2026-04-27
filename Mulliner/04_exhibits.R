# ---- Exhibits: plots and tables -------------------------------------------

# Color palette close to the paper's aesthetic
pal <- list(
  blue  = "#1f4e79",
  red   = "#c0504d",
  green = "#2e9a46",
  grey  = "#7f7f7f",
  light_blue = "#4f81bd"
)

# -------------------- Exhibit 2: raw state variables -----------------------

exhibit_2 <- function(macro_m, variables) {
  df <- macro_m %>%
    pivot_longer(-date, names_to = "var", values_to = "x") %>%
    left_join(variables, by = "var") %>%
    mutate(label = factor(label, levels = variables$label))

  p <- ggplot(df, aes(date, x)) +
    geom_line(colour = pal$blue, linewidth = 0.4) +
    facet_wrap(~ label, scales = "free_y", ncol = 2) +
    labs(title = "Exhibit 2. Raw economic state variables",
         x = NULL, y = NULL) +
    theme_regimes()

  save_plot(p, "exhibit_02_raw_state_variables", width = 10, height = 10)
}

# -------------------- Exhibit 3: z-score + winsorized ----------------------

exhibit_3 <- function(transformed, transformed_w, variables) {
  long <- transformed %>%
    pivot_longer(-date, names_to = "var", values_to = "z_raw") %>%
    left_join(
      transformed_w %>% pivot_longer(-date, names_to = "var", values_to = "z_wins"),
      by = c("date", "var")
    ) %>%
    left_join(variables, by = "var") %>%
    mutate(label = factor(label, levels = variables$label))

  p <- ggplot(long, aes(date)) +
    geom_line(aes(y = z_raw, colour = "raw z-score"), linewidth = 0.4) +
    geom_line(aes(y = z_wins, colour = "winsorized"), linewidth = 0.4) +
    facet_wrap(~ label, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = c("raw z-score" = pal$blue,
                                   "winsorized"  = pal$red),
                        name = NULL) +
    labs(title = "Exhibit 3. Transformed economic state variables",
         x = NULL, y = "z-score") +
    theme_regimes()

  save_plot(p, "exhibit_03_transformed_state_variables", width = 10, height = 10)
}

# -------------------- Exhibit 4: autocorrelations --------------------------

exhibit_4 <- function(transformed_w, variables) {
  acf_at <- function(x, k) {
    x <- x[!is.na(x)]
    if (length(x) < k + 5) return(NA_real_)
    a <- acf(x, lag.max = k, plot = FALSE, na.action = na.pass)
    as.numeric(a$acf[k + 1])
  }

  tbl <- map_dfr(variables$var, function(v) {
    x <- transformed_w[[v]]
    tibble(
      variable = variables$label[variables$var == v],
      `1-month`  = acf_at(x, 1),
      `3-month`  = acf_at(x, 3),
      `12-month` = acf_at(x, 12),
      `3-year`   = acf_at(x, 36),
      `10-year`  = acf_at(x, 120),
      monthly_mean = mean(x, na.rm = TRUE),
      std          = sd(x, na.rm = TRUE),
      frequency    = "monthly"
    )
  })

  save_table(tbl, "exhibit_04_autocorrelations")
  print(tbl %>% mutate(across(where(is.numeric), ~ round(.x, 2))))
}

# -------------------- Exhibit 5: correlation heatmap -----------------------

exhibit_5 <- function(transformed_w, variables) {
  X <- as.matrix(transformed_w[, variables$var])
  colnames(X) <- variables$label
  C <- cor(X, use = "pairwise.complete.obs")
  C_long <- as_tibble(C, rownames = "row") %>%
    pivot_longer(-row, names_to = "col", values_to = "rho") %>%
    mutate(row = factor(row, levels = rev(variables$label)),
           col = factor(col, levels = variables$label))

  p <- ggplot(C_long, aes(col, row, fill = rho)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3.2) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c53030",
                         midpoint = 0, limits = c(-1, 1), name = NULL) +
    labs(title = "Exhibit 5. Correlation of the economic state variables",
         x = NULL, y = NULL) +
    theme_regimes() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())

  save_plot(p, "exhibit_05_correlation_heatmap", width = 7, height = 6)
  save_table(as_tibble(C, rownames = "variable"),
             "exhibit_05_correlation_matrix")
}

# -------------------- Exhibits 6-8: crisis similarity plots ----------------

similarity_strip_plot <- function(dates, D, target_date, top_pct, mask_months,
                                  title, highlight_bands = NULL) {
  i <- which(dates == target_date)
  if (length(i) == 0) {
    i <- which.min(abs(as.numeric(dates - target_date)))
    warning("Target date not exact; using ", dates[i])
  }
  d <- D[i, ]
  upper <- i - mask_months
  eligible <- seq_len(max(upper, 1))

  df <- tibble(date = dates, d = d) %>%
    mutate(eligible = row_number() <= upper,
           masked   = row_number() >  upper & row_number() <= i)

  # Restrict the plot to months where the distance is actually computable.
  # Without this, ggplot auto-expands the x-axis back to the first macro
  # date even though no line is drawn there.
  plot_df <- df %>% filter(!is.na(d))
  x_min <- min(plot_df$date); x_max <- max(plot_df$date)

  # selected similar = lowest 'top_pct' within eligible
  d_elig <- df$d[df$eligible]
  valid <- !is.na(d_elig)
  k <- max(1L, floor(sum(valid) * top_pct))
  thresh <- sort(d_elig[valid])[k]
  sim_dates <- df %>%
    filter(eligible, !is.na(d), d <= thresh) %>%
    pull(date)

  # Build rectangle spans for similar months as month-wide bars. This avoids
  # geom_col's reliance on the current y-range and always covers the full
  # vertical extent of the panel.
  sim_rects <- if (length(sim_dates) > 0) {
    tibble(
      xmin = sim_dates - 15,   # ~half a month
      xmax = sim_dates + 15
    )
  } else tibble(xmin = as.Date(character()), xmax = as.Date(character()))

  masked_df <- df %>% filter(masked) %>%
    summarise(xmin = min(date), xmax = max(date))

  p <- ggplot(plot_df, aes(date, d))

  if (nrow(masked_df) && is.finite(masked_df$xmin)) {
    p <- p + geom_rect(data = masked_df,
                       aes(xmin = xmin, xmax = xmax,
                           ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE,
                       fill = "grey80", alpha = 0.6)
  }
  if (!is.null(highlight_bands)) {
    p <- p + geom_rect(data = highlight_bands,
                       aes(xmin = start, xmax = end,
                           ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE,
                       fill = pal$red, alpha = 0.10, colour = NA)
  }
  if (nrow(sim_rects) > 0) {
    p <- p + geom_rect(data = sim_rects,
                       aes(xmin = xmin, xmax = xmax,
                           ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE,
                       fill = pal$green, alpha = 0.30, colour = NA)
  }

  p +
    geom_line(colour = pal$blue, linewidth = 0.5, linetype = "dashed") +
    scale_x_date(limits = c(x_min, x_max), expand = expansion(mult = 0.01)) +
    labs(title = title, x = NULL, y = "Global score (distance)") +
    theme_regimes()
}

exhibit_6 <- function(distance_matrix, mask_months, top_pct) {
  p <- similarity_strip_plot(
    distance_matrix$dates, distance_matrix$D,
    target_date = as.Date("2009-01-31"),
    top_pct     = top_pct,
    mask_months = mask_months,
    title       = "Exhibit 6. Global similarity score calculated for January 2009"
  )
  save_plot(p, "exhibit_06_similarity_2009_01", width = 10, height = 4.5)
}

exhibit_7 <- function(distance_matrix, mask_months, top_pct) {
  p1 <- similarity_strip_plot(
    distance_matrix$dates, distance_matrix$D,
    target_date = as.Date("2020-02-29"),
    top_pct = top_pct, mask_months = mask_months,
    title = "Global similarity score calculated for February 2020"
  )
  p2 <- similarity_strip_plot(
    distance_matrix$dates, distance_matrix$D,
    target_date = as.Date("2020-04-30"),
    top_pct = top_pct, mask_months = mask_months,
    title = "Global similarity score calculated for April 2020"
  )
  p <- p1 / p2 + plot_annotation(title = "Exhibit 7. Similarity during COVID-19")
  save_plot(p, "exhibit_07_similarity_covid", width = 10, height = 8)
}

exhibit_8 <- function(distance_matrix, mask_months, top_pct) {
  # historical inflation periods (Table in Exhibit 8) -- we can only
  # highlight the ones within our sample
  bands <- tibble::tribble(
    ~start,         ~end,
    "1941-04-01",   "1942-05-31",
    "1946-03-01",   "1947-03-31",
    "1950-08-01",   "1951-02-28",
    "1966-02-01",   "1970-01-31",
    "1972-07-01",   "1974-12-31",
    "1977-02-01",   "1980-03-31",
    "1987-02-01",   "1990-11-30",
    "2007-09-01",   "2008-07-31",
    "2022-01-01",   "2022-10-31"
  ) %>% mutate(start = as.Date(start), end = as.Date(end))

  p <- similarity_strip_plot(
    distance_matrix$dates, distance_matrix$D,
    target_date = as.Date("2022-08-31"),
    top_pct = top_pct, mask_months = mask_months,
    title = "Exhibit 8. Similarity of August 2022 inflation surge",
    highlight_bands = bands
  )
  save_plot(p, "exhibit_08_similarity_2022_08", width = 10, height = 5)
}

# -------------------- Exhibit 9: EWMA of global score ----------------------

exhibit_9 <- function(distance_matrix, lookbacks) {
  df <- map_dfr(lookbacks, function(n_mo) {
    out <- compute_ewma_scores(distance_matrix$D, distance_matrix$dates, n_mo)
    out %>% mutate(lookback = paste0(n_mo / 12, "-year"))
  })
  mean_tbl <- df %>%
    group_by(date) %>%
    summarise(ewma = mean(ewma, na.rm = TRUE), .groups = "drop") %>%
    mutate(lookback = "Mean")

  plot_df <- bind_rows(df, mean_tbl) %>%
    filter(is.finite(ewma)) %>%   # drop pre-data rows so x-axis doesn't extend
    mutate(lookback = factor(lookback,
              levels = c(paste0(lookbacks / 12, "-year"), "Mean")))

  p <- ggplot(plot_df, aes(date, ewma, colour = lookback, linewidth = lookback)) +
    geom_line() +
    scale_linewidth_manual(values = c(rep(0.4, length(lookbacks)), 0.9),
                           guide = "none") +
    scale_colour_manual(values = c(pal$blue, pal$light_blue,
                                   pal$green, pal$grey, "black")) +
    scale_x_date(expand = expansion(mult = 0.01)) +
    labs(title = "Exhibit 9. EWMA of global score across lookbacks",
         x = NULL, y = "EWMA of global score", colour = NULL) +
    theme_regimes()

  save_plot(p, "exhibit_09_ewma_regime_shifts", width = 10, height = 5)
}

# -------------------- Exhibit 1: yearly returns (vol-targeted) -------------

exhibit_1 <- function(diff_vt) {
  yearly <- diff_vt %>%
    filter(is.finite(r_vt)) %>%
    mutate(year = year(date)) %>%
    group_by(year) %>%
    summarise(ret = sum(r_vt, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n == 12)

  pos <- yearly$ret[yearly$ret > 0]
  neg <- yearly$ret[yearly$ret <= 0]
  stats <- tibble(
    years = paste0(min(yearly$year), "-", max(yearly$year)),
    pct_positive = length(pos) / nrow(yearly),
    mean_when_pos = mean(pos),
    mean_when_neg = mean(neg)
  )
  save_table(stats, "exhibit_01_yearly_return_stats")
  save_table(yearly, "exhibit_01_yearly_returns")
  print(stats)

  p <- ggplot(yearly, aes(year, ret * 100, fill = ret > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = pal$blue, "FALSE" = pal$red),
                      guide = "none") +
    labs(title = "Exhibit 1. Yearly returns: long similarity, short dissimilarity",
         subtitle = "Q1 minus Q5, scaled to 15% annualised vol",
         x = NULL, y = "% return") +
    theme_regimes()

  save_plot(p, "exhibit_01_yearly_returns", width = 10, height = 4.5)
}

# -------------------- Exhibit 10: quintile cumulative returns + diff -------

cum_pct <- function(r) cumsum(ifelse(is.na(r), 0, r)) * 100

# Internal: renders one Exhibit-10 variant (unlevered OR vol-targeted) given
# the relevant per-bucket, long-only, and diff inputs.
exhibit_10_one <- function(per_bucket, long_only, diff_default,
                           suffix, subtitle_tail) {
  qs <- per_bucket %>% filter(is.finite(avg))
  lo <- long_only  %>% filter(is.finite(avg))

  start_date <- qs %>% filter(quintile == "Q1") %>%
    summarise(d = min(date)) %>% pull(d)
  qs <- qs %>% filter(date >= start_date)
  lo <- lo %>% filter(date >= start_date)

  qs <- qs %>% group_by(quintile) %>% mutate(cum = cum_pct(avg)) %>% ungroup()
  lo <- lo %>% mutate(cum = cum_pct(avg))

  stats_q <- qs %>%
    group_by(quintile) %>%
    summarise(SR = sharpe(avg),
              corr = cor(avg, lo$avg[match(date, lo$date)], use = "pairwise"),
              .groups = "drop")
  save_table(stats_q, paste0("exhibit_10_quintile_stats", suffix))
  cat("  quintile stats:\n")
  print(stats_q)

  p1 <- ggplot() +
    geom_line(data = qs, aes(date, cum, colour = quintile), linewidth = 0.6) +
    geom_line(data = lo, aes(date, cum), colour = "black",
              linetype = "dashed", linewidth = 0.6) +
    labs(title = "Exhibit 10 (left). Similarity quintiles vs long-only",
         subtitle = paste0("Quintile 1 = most similar, Quintile 5 = least similar. ",
                           "Black dashed: LO. ", subtitle_tail),
         x = NULL, y = "Cumulative return (%)", colour = NULL) +
    theme_regimes()

  diff_df <- diff_default %>% filter(is.finite(diff)) %>%
    mutate(cum = cum_pct(diff))
  sr_diff <- sharpe(diff_default$diff)
  corr_diff <- cor(diff_default$diff,
                   long_only$avg[match(diff_default$date, long_only$date)],
                   use = "pairwise")

  p2 <- ggplot(diff_df, aes(date, cum)) +
    geom_line(colour = pal$blue, linewidth = 0.7) +
    labs(title = sprintf("Exhibit 10 (right). Q1 minus Q5 (SR = %.2f, corr LO = %.2f)",
                         sr_diff, corr_diff),
         subtitle = subtitle_tail,
         x = NULL, y = "Cumulative return (%)") +
    theme_regimes()

  save_plot(p1 / p2, paste0("exhibit_10_quintiles_and_diff", suffix),
            width = 11, height = 9)
}

exhibit_10 <- function(strategies) {
  # Vol-targeted at 15% ann (paper convention)
  exhibit_10_one(strategies$per_bucket_vt,
                 strategies$long_only_vt,
                 strategies$diff_default_vt,
                 suffix        = "",
                 subtitle_tail = "Each series scaled to 15% annualised vol (paper convention).")
}

# -------------------- Exhibit 11: drawdowns --------------------------------

exhibit_11_one <- function(long_only, diff_default, suffix, subtitle_tail) {
  lo_dd   <- tibble(date = long_only$date,
                    dd   = drawdown(long_only$avg),
                    what = "Long-only factor model")

  diff_dd <- tibble(date = diff_default$date,
                    dd   = drawdown(diff_default$diff),
                    what = "Similarity model (Q1 - Q5)")

  df <- bind_rows(lo_dd, diff_dd) %>% filter(is.finite(dd))

  p <- ggplot(df, aes(date, dd * 100, colour = what)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = c("Long-only factor model"        = pal$grey,
                                   "Similarity model (Q1 - Q5)"    = pal$blue),
                        name = NULL) +
    labs(title = "Exhibit 11. Drawdown profile",
         subtitle = subtitle_tail,
         x = NULL, y = "Drawdown (%)") +
    theme_regimes()

  save_plot(p, paste0("exhibit_11_drawdowns", suffix), width = 10, height = 4.5)
}

exhibit_11 <- function(strategies) {
  exhibit_11_one(strategies$long_only_vt,
                 strategies$diff_default_vt,
                 suffix        = "",
                 subtitle_tail = "Each series scaled to 15% annualised vol (paper convention).")
}

# -------------------- Exhibit 12: quantile robustness ----------------------

exhibit_12 <- function(strategies) {
  # Vol-target the diff series per quantile-count group to 15% ann
  rb <- strategies$robust_qtl %>%
    filter(is.finite(diff)) %>%
    group_by(q) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      rv    = slider::slide_dbl(diff, ~ sd(.x, na.rm = TRUE),
                                .before = CFG$vol_window, .after = -1,
                                .complete = FALSE),
      scale = (CFG$vol_target / sqrt(12)) / rv,
      scale = if_else(is.finite(scale), scale, NA_real_),
      diff  = diff * scale
    ) %>%
    filter(is.finite(diff)) %>%
    mutate(cum = cum_pct(diff)) %>%
    ungroup()

  stats <- rb %>% group_by(q) %>%
    summarise(SR = sharpe(diff), .groups = "drop") %>%
    arrange(q)
  save_table(stats, "exhibit_12_quantile_robustness_stats")
  print(stats)

  rb <- rb %>% mutate(label = sprintf("%d quantiles (SR %.2f)", q,
                                      stats$SR[match(q, stats$q)]))
  p <- ggplot(rb, aes(date, cum, colour = label)) +
    geom_line(linewidth = 0.55) +
    labs(title = "Exhibit 12. Robustness to the number of quantiles",
         subtitle = "Each series scaled to 15% annualised vol (paper convention).",
         x = NULL, y = "Cumulative return (%)", colour = NULL) +
    theme_regimes()
  save_plot(p, "exhibit_12_quantile_robustness", width = 10, height = 5)
}

# -------------------- Exhibit 13: lookback robustness ----------------------

exhibit_13 <- function(strategies) {
  # Vol-target the diff series per lookback group to 15% ann
  rb <- strategies$robust_lkback %>%
    filter(is.finite(diff)) %>%
    group_by(lookback_years) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      rv    = slider::slide_dbl(diff, ~ sd(.x, na.rm = TRUE),
                                .before = CFG$vol_window, .after = -1,
                                .complete = FALSE),
      scale = (CFG$vol_target / sqrt(12)) / rv,
      scale = if_else(is.finite(scale), scale, NA_real_),
      diff  = diff * scale
    ) %>%
    filter(is.finite(diff)) %>%
    mutate(cum = cum_pct(diff)) %>%
    ungroup()

  stats <- rb %>% group_by(lookback_years) %>%
    summarise(SR = sharpe(diff), .groups = "drop") %>%
    arrange(lookback_years)
  save_table(stats, "exhibit_13_lookback_robustness_stats")
  print(stats)

  rb <- rb %>%
    mutate(label = sprintf("%d-year lookback (SR %.2f)",
                           as.integer(lookback_years),
                           stats$SR[match(lookback_years, stats$lookback_years)]))
  p <- ggplot(rb, aes(date, cum, colour = label)) +
    geom_line(linewidth = 0.55) +
    labs(title = "Exhibit 13. Sensitivity to z-score lookback window",
         subtitle = "Each series scaled to 15% annualised vol (paper convention).",
         x = NULL, y = "Cumulative return (%)", colour = NULL) +
    theme_regimes()
  save_plot(p, "exhibit_13_lookback_robustness", width = 10, height = 5)
}

# -------------------- Appendix A1 & A2: individual factors -----------------

exhibit_a_one <- function(per_factor_data, suffix, subtitle_tail) {
  per_fac <- per_factor_data %>%
    filter(is.finite(r)) %>%
    group_by(factor, quintile) %>%
    arrange(date) %>%
    mutate(cum = cum_pct(r)) %>%
    ungroup()

  # A1: quintile cumulative returns per factor
  p_a1 <- ggplot(per_fac, aes(date, cum, colour = quintile)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ factor, scales = "free_y") +
    labs(title = "Exhibit A1. Similarity quintiles per factor",
         subtitle = subtitle_tail,
         x = NULL, y = "Cumulative return (%)", colour = NULL) +
    theme_regimes()
  save_plot(p_a1, paste0("exhibit_A1_per_factor_quintiles", suffix),
            width = 11, height = 7)

  # A2: Q1 - Q5 per factor
  diff_fac <- per_fac %>%
    filter(quintile %in% c("Q1", "Q5")) %>%
    select(date, factor, quintile, r) %>%
    pivot_wider(names_from = quintile, values_from = r) %>%
    mutate(diff = .data[["Q1"]] - .data[["Q5"]]) %>%
    filter(is.finite(diff)) %>%
    group_by(factor) %>% arrange(date) %>%
    mutate(cum = cum_pct(diff)) %>% ungroup()

  stats_a2 <- diff_fac %>% group_by(factor) %>%
    summarise(SR = sharpe(diff), .groups = "drop")
  save_table(stats_a2, paste0("exhibit_A2_per_factor_stats", suffix))
  cat("  per-factor Q1-Q5 Sharpes:\n")
  print(stats_a2)

  p_a2 <- ggplot(diff_fac, aes(date, cum)) +
    geom_line(colour = pal$blue, linewidth = 0.55) +
    facet_wrap(~ factor, scales = "free_y") +
    labs(title = "Exhibit A2. Q1 minus Q5 per factor",
         subtitle = paste0(subtitle_tail, " Horizontal stretches are months where Q1 and Q5 signals agree for that factor (zero contribution to the diff)."),
         x = NULL, y = "Cumulative return (%)") +
    theme_regimes()
  save_plot(p_a2, paste0("exhibit_A2_per_factor_diff", suffix),
            width = 11, height = 7)
}

exhibit_a <- function(strategies) {
  exhibit_a_one(strategies$per_factor_vt,
                suffix        = "",
                subtitle_tail = "Each factor-quintile scaled to 15% annualised vol (paper convention).")
}

# -------------------- Orchestrator -----------------------------------------

make_all_exhibits <- function(macro_monthly, transformed, winsorized,
                              distance_mat, factors, strategies, cfg) {
  message("  - Exhibit 2  (raw variables)")
  exhibit_2(macro_monthly, cfg$variables)
  message("  - Exhibit 3  (transformed variables)")
  exhibit_3(transformed, winsorized, cfg$variables)
  message("  - Exhibit 4  (autocorrelations)")
  exhibit_4(winsorized, cfg$variables)
  message("  - Exhibit 5  (correlations)")
  exhibit_5(winsorized, cfg$variables)
  message("  - Exhibit 6  (Jan 2009)")
  exhibit_6(distance_mat, cfg$mask_months, cfg$paper_top_pct)
  message("  - Exhibit 7  (COVID 2020)")
  exhibit_7(distance_mat, cfg$mask_months, cfg$paper_top_pct)
  message("  - Exhibit 8  (Aug 2022)")
  exhibit_8(distance_mat, cfg$mask_months, cfg$paper_top_pct)
  message("  - Exhibit 9  (EWMA)")
  exhibit_9(distance_mat, cfg$ewma_lookbacks)
  message("  - Exhibit 10 (quintile strategy)")
  exhibit_10(strategies)
  message("  - Exhibit 11 (drawdowns)")
  exhibit_11(strategies)
  message("  - Exhibit 12 (quantile robustness)")
  exhibit_12(strategies)
  message("  - Exhibit 13 (lookback robustness)")
  exhibit_13(strategies)
  message("  - Exhibit 1  (yearly returns, vol-targeted)")
  exhibit_1(strategies$diff_vol_targeted)
  message("  - Appendix A1 & A2 (per-factor)")
  exhibit_a(strategies)
}
