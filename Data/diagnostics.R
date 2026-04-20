# =========================================================
# Diagnostics for Regimes replication
# Run after 02_analysis.R (objects state_full2, D2, positions,
# F2, dates2, rf_aligned, monthly_wide, transformed_long, state_z
# must be in the global env)
# =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(lubridate); library(ggplot2)
})

sep <- function(s) cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")

# =========================================================
# CHECK 1 — Data coverage per state variable
# What's the binding constraint on sample start?
# =========================================================
sep("CHECK 1: Data coverage per state variable")

coverage <- monthly_wide %>%
  pivot_longer(-date, names_to = "series", values_to = "value") %>%
  filter(series %in% VAR_NAMES, !is.na(value)) %>%
  group_by(series) %>%
  summarise(
    raw_start = min(date),
    raw_end   = max(date),
    n_raw     = n(),
    .groups = "drop"
  )

z_coverage <- transformed_long %>%
  filter(series %in% VAR_NAMES, !is.na(z)) %>%
  group_by(series) %>%
  summarise(
    z_start = min(date),
    z_end   = max(date),
    n_z     = n(),
    .groups = "drop"
  )

cov_tbl <- coverage %>% left_join(z_coverage, by = "series")
print(cov_tbl)

cat(sprintf("\nFirst date with ALL 7 z-scores non-NA: %s\n",
            min(state_full2$date)))
cat(sprintf("Effective factor-timing sample start (after 36m mask + 1m warmup): %s\n",
            dates2[37 + 20]))  # need 20+ eligible historical obs

# =========================================================
# CHECK 2 — State variable descriptive stats vs paper Exhibit 4
# Paper reports monthly mean ~0 and std ~1 for all variables
# =========================================================
sep("CHECK 2: State variable stats (compare to paper Exhibit 4)")
stats_check <- transformed_long %>%
  filter(series %in% VAR_NAMES, !is.na(z)) %>%
  group_by(series) %>%
  summarise(
    mean_z = mean(z),
    std_z  = sd(z),
    min_z  = min(z),
    max_z  = max(z),
    frac_winsorized = mean(abs(z_raw) > 3, na.rm = TRUE),
    .groups = "drop"
  )
print(stats_check)
cat("\nPaper's values (mean, std):\n",
    "  market: 0.49, 0.89 | yield curve: -0.07, 0.98 | oil: 0.24, 0.97\n",
    "  copper: 0.15, 0.94 | monetary policy: 0.13, 0.98\n",
    "  volatility: 0.05, 0.97 | stock-bond: -0.10, 0.98\n")

# =========================================================
# CHECK 3 — Cross-correlations vs paper Exhibit 5
# =========================================================
sep("CHECK 3: Cross-correlations (compare to paper Exhibit 5)")
cor_check <- cor(state_full2 %>% select(all_of(VAR_NAMES)))
print(round(cor_check, 2))
cat("\nPaper's key correlations to check:\n",
    "  oil-copper: 0.33 | slope-tbill_3m: -0.67\n",
    "  slope-monetary: -0.36 | market-volatility: -0.25\n")

# =========================================================
# CHECK 4 — Fama-French factor sanity
# Long-only SR per factor should roughly match published values
# =========================================================
sep("CHECK 4: FF factor Sharpe ratios (sanity)")
ff_check <- as_tibble(F2) %>%
  summarise(across(everything(),
                   list(mean_ann = ~ mean(.x) * 12,
                        sd_ann   = ~ sd(.x) * sqrt(12),
                        SR       = ~ mean(.x)/sd(.x) * sqrt(12)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(),
               names_to = c("factor", "stat"), names_sep = "__",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value)
print(ff_check)
cat("\nRough expected SRs (full-sample published): market~0.5, value~0.3, momentum~0.5\n")
cat("If momentum SR is near 0 or negative in your sample, you're in the 'momentum drought' post-2001 era — that's a real thing, not a bug.\n")

# =========================================================
# CHECK 5 — Date alignment between state and factor returns
# Pick a random month; verify state_at_T matches FF_at_T
# =========================================================
sep("CHECK 5: Date alignment state vs factors")
align_check <- state_full2 %>%
  select(date, sp500, volatility) %>%
  inner_join(rf_aligned %>% select(date, market, momentum), by = "date") %>%
  tail(10)
print(align_check)
cat("\nDates should all be end-of-month, no gaps. If `market` column has many NAs but sp500 doesn't, FF dates are off by one month.\n")

# Explicit verification: FF factor dates
ff_date_sample <- rf_aligned %>%
  filter(year(date) %in% c(2008, 2009, 2020)) %>%
  head(15)
cat("\nFF factor date sample (should be month-ends):\n")
print(ff_date_sample$date)

# =========================================================
# CHECK 6 — Position distribution per factor
# For each factor, what fraction of months is the position +1?
# Long-only SR ~ 0.7 in your sample means factors mostly pay,
# so positions should be +1 >50% of the time. If near 50/50 the
# signal has no directional conviction.
# =========================================================
sep("CHECK 6: Position distribution (Q1) per factor")

pos_q1 <- positions[, 1, ]   # shape: (N, K)
pos_dist <- apply(pos_q1, 2, function(x) {
  x <- x[!is.na(x)]
  c(n = length(x),
    frac_long  = mean(x > 0),
    frac_short = mean(x < 0),
    frac_zero  = mean(x == 0))
}) %>% t()
print(round(pos_dist, 3))

cat("\nRed flag: if frac_long is near 0.5 for every factor, the similar-date mean is just noise.\n")
cat("Paper's strategy is long 80%+ of the time for strong factors (hence 0.76 corr to LO).\n")

# Same per quintile — are all quintiles similar? Would indicate signal collapse
pos_dist_all <- map_dfr(1:5, function(q) {
  mat <- positions[, q, ]
  tibble(
    quintile = paste0("Q", q),
    factor   = FACTOR_NAMES,
    frac_long = apply(mat, 2, function(x) mean(x > 0, na.rm = TRUE)),
    n_obs     = apply(mat, 2, function(x) sum(!is.na(x)))
  )
})
cat("\nFrac long by quintile (wide):\n")
print(pos_dist_all %>%
        select(-n_obs) %>%
        pivot_wider(names_from = quintile, values_from = frac_long) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2))))

# =========================================================
# CHECK 7 — Look-ahead integrity
# Position at T must only use state info observable by T and
# returns from periods strictly before T (minus 36m mask).
# =========================================================
sep("CHECK 7: Look-ahead check")

test_T <- floor(nrow(state_full2) * 0.7)   # some month in the middle
anchor_date <- dates2[test_T]
cat(sprintf("Testing anchor T = %s (index %d)\n", anchor_date, test_T))

# Recompute manually
max_i <- test_T - 36 - 1
eligible <- seq_len(max_i)
ord <- eligible[order(D2[test_T, eligible])]
q1_sel <- ord[1:ceiling(length(ord)/5)]

cat(sprintf("  Eligible dates run from %s to %s\n",
            dates2[1], dates2[max_i]))
cat(sprintf("  Gap between last eligible and T: %d months (should be ≥36)\n",
            interval(dates2[max_i], anchor_date) %/% months(1)))
cat(sprintf("  Q1 selected dates range: %s to %s\n",
            min(dates2[q1_sel]), max(dates2[q1_sel])))
cat(sprintf("  Most recent eligible 'subsequent' return used: %s (must be ≤ T=%s)\n",
            dates2[max(q1_sel) + 1], anchor_date))
stopifnot(all(dates2[q1_sel + 1] < anchor_date))
cat("  ✓ No look-ahead: all subsequent returns used are strictly before T\n")

# =========================================================
# CHECK 8 — Replicate paper Exhibit 6 (Jan 2009 similar months)
# Paper says result should include double-dip recessions of early 1980s
# =========================================================
sep("CHECK 8: Jan 2009 most-similar dates (paper Exhibit 6)")

jan09_idx <- which(dates2 == as.Date("2009-01-31"))
if (length(jan09_idx) == 1) {
  eligible09 <- seq_len(jan09_idx - 37)
  ord09 <- eligible09[order(D2[jan09_idx, eligible09])]
  top15 <- ord09[1:ceiling(length(ord09) * 0.15)]
  cat("Top 15% most similar months to Jan 2009, grouped by year:\n")
  print(table(year(dates2[top15])))
  cat("\nPaper expects hits in: 1969-70, 1973-74, 1980-82, 1990, 2001-02\n")
  cat("If your sample starts ~2001, you cannot reproduce this — only 2001-02 is in your history.\n")
} else {
  cat("Jan 2009 not in sample — check start date\n")
}

# =========================================================
# CHECK 9 — Manual return reconstruction for one month
# Take a specific T, compute Q1 return manually, compare to pipeline
# =========================================================
sep("CHECK 9: Manual vs pipeline Q1 return for a specific month")

T_test <- floor(nrow(state_full2) * 0.6)
max_i  <- T_test - 37
eligible_t <- seq_len(max_i)
ord_t <- eligible_t[order(D2[T_test, eligible_t])]
q1_t  <- ord_t[1:ceiling(length(ord_t) / 5)]

subs_rets <- F2[q1_t + 1, ]
mean_subs <- colMeans(subs_rets, na.rm = TRUE)
pos_manual <- sign(mean_subs)
fwd_ret <- F2[T_test + 1, ]
manual_q1_ret <- mean(pos_manual * fwd_ret)

# Pipeline Q1 return for T_test+1
pipeline_q1 <- qrets_agg %>%
  filter(quintile == "Q1", date == dates2[T_test + 1]) %>%
  pull(ew)

cat(sprintf("Anchor T = %s, trading month = %s\n",
            dates2[T_test], dates2[T_test + 1]))
cat("Per-factor manual positions:\n")
print(round(pos_manual, 0))
cat(sprintf("Manual Q1 EW return:   %.5f\n", manual_q1_ret))
cat(sprintf("Pipeline Q1 EW return: %.5f\n", pipeline_q1))
cat(sprintf("Match: %s\n", isTRUE(all.equal(manual_q1_ret, pipeline_q1, tolerance = 1e-8))))

# =========================================================
# CHECK 10 — Is the distance metric doing what we think?
# For very similar months, Euclidean distance should be small.
# Print distances for months near T and verify not-too-small due
# to the 12m diff autocorrelation leaking in.
# =========================================================
sep("CHECK 10: Distance structure sanity")

T_spot <- which(dates2 == as.Date("2020-03-31"))
if (length(T_spot) == 1) {
  nearby <- D2[T_spot, (T_spot-6):(T_spot-1)]
  names(nearby) <- format(dates2[(T_spot-6):(T_spot-1)])
  cat("Distances from 2020-03-31 to 6 preceding months (should be small, 12m diff makes them autocorrelated):\n")
  print(round(nearby, 3))

  far <- D2[T_spot, c(1, 50, 100)]
  names(far) <- format(dates2[c(1, 50, 100)])
  cat("\nDistances to early-sample dates (should be larger on average):\n")
  print(round(far, 3))
}

# =========================================================
# CHECK 11 — Stock-bond correlation sanity
# This is the most-likely-to-be-wrong variable.
# 3y rolling corr of SPX logret vs d(10y) should be mostly negative
# post-2000, positive in 1970s-90s.
# =========================================================
sep("CHECK 11: Stock-bond correlation regimes")

sb_plot_df <- monthly_wide %>%
  select(date, stock_bond) %>%
  filter(!is.na(stock_bond)) %>%
  mutate(decade = floor(year(date) / 10) * 10) %>%
  group_by(decade) %>%
  summarise(mean_corr = mean(stock_bond),
            n = n(), .groups = "drop")
print(sb_plot_df)
cat("\nExpected: positive in 1970s-80s-90s (inflation dominant),\n")
cat("          negative post-2000 (flight-to-quality dominant).\n")

# =========================================================
# CHECK 12 — VIX prepend continuity
# Check the 1989-1990 boundary where we switch from realized vol to VIX
# =========================================================
sep("CHECK 12: VIX/realized-vol prepend continuity")

vix_boundary <- monthly_wide %>%
  select(date, volatility) %>%
  filter(year(date) %in% 1988:1991, !is.na(volatility))
print(vix_boundary)
cat("\nIf there's a big jump between Dec-1989 and Jan-1990, the realized-vol scaling is off.\n")

# =========================================================
# CHECK 13 — Positions and returns time-series plot per factor
# Visual check of whether positions make sense
# =========================================================
sep("CHECK 13: Saving position diagnostic plot")

pos_long <- map_dfr(seq_along(FACTOR_NAMES), function(k) {
  tibble(date = dates2,
         factor = FACTOR_NAMES[k],
         pos   = positions[, 1, k])
}) %>% filter(!is.na(pos))

p_diag <- ggplot(pos_long, aes(date, pos)) +
  geom_step(colour = "#1f4e79") +
  facet_wrap(~factor, ncol = 2) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  labs(title = "Q1 positions over time (1 = long, -1 = short)",
       x = NULL, y = "position")
ggsave("plots/diag_positions.png", p_diag, width = 10, height = 7, dpi = 140)
cat("Saved plots/diag_positions.png\n")

cat("\n", strrep("=", 70), "\nDONE — scan CHECK 2, 3, 5, 6, 8 first; those flag the most common issues.\n", sep = "")