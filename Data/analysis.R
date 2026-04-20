# =========================================================
# Regimes (Mulliner et al., 2025) — Analysis & Exhibits
# =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(lubridate); library(zoo); library(ggplot2); library(scales)
  library(patchwork)
})

# source("01_data.R")   # uncomment to refresh data
raw <- readRDS("data/regimes_raw.rds")

dir.create("plots", showWarnings = FALSE)
theme_set(theme_minimal(base_size = 11))

# =========================================================
# 1) Build daily & monthly feeds for the 7 state variables
# =========================================================

# Wide daily frame on trading days (use SPX calendar as master)
daily_wide <- raw$spx_daily %>%
  transmute(date, sp500 = value) %>%
  left_join(raw$fred_daily %>% pivot_wider(names_from = series, values_from = value),
            by = "date") %>%
  arrange(date)

# Stock-bond correlation: rolling 3Y of daily logret(SPX) vs diff(10y yield)
sb_corr_daily <- daily_wide %>%
  arrange(date) %>%
  transmute(date,
            spx_ret  = c(NA, diff(log(sp500))),
            d_y10    = c(NA, diff(yield_10y))) %>%
  # Inner-filter to rows where both are observable
  mutate(pair_ok = !is.na(spx_ret) & !is.na(d_y10))

win_corr <- 756  # ~3 * 252 trading days
sb_corr_daily$stock_bond <- NA_real_
ok_idx <- which(sb_corr_daily$pair_ok)
# Compact series and roll-correlate, then map back
x <- sb_corr_daily$spx_ret[ok_idx]
y <- sb_corr_daily$d_y10[ok_idx]
rc <- rollapplyr(cbind(x, y), width = win_corr,
                 FUN = function(m) cor(m[,1], m[,2], use = "pairwise.complete.obs"),
                 by.column = FALSE, fill = NA_real_)
sb_corr_daily$stock_bond[ok_idx] <- rc

# Realized vol on SPX: 21-day trailing annualized in %, used pre-1990
vol_daily <- daily_wide %>%
  transmute(date, spx_ret = c(NA, diff(log(sp500)))) %>%
  mutate(realized_vol = rollapplyr(spx_ret, 21, sd, fill = NA) * sqrt(252) * 100)

# Combined volatility: VIX when available, realized vol before
vol_combined <- daily_wide %>%
  select(date, vix) %>%
  left_join(vol_daily %>% select(date, realized_vol), by = "date") %>%
  mutate(volatility = coalesce(vix, realized_vol)) %>%
  select(date, volatility)

# Build daily_wide with stock-bond + combined volatility
daily_full <- daily_wide %>%
  left_join(sb_corr_daily %>% select(date, stock_bond), by = "date") %>%
  left_join(vol_combined, by = "date")

# Monthly end-of-month values (last trading day of each month)
to_monthly_eom <- function(df) {
  df %>%
    mutate(ym = floor_date(date, "month")) %>%
    group_by(ym) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(date = ceiling_date(ym, "month") - days(1)) %>%
    select(-ym)
}

monthly_daily_src <- daily_full %>%
  select(date, sp500, slope, tbill_3m, oil, volatility, stock_bond) %>%
  to_monthly_eom()

# Copper is already monthly; align to end-of-month
monthly_copper <- raw$fred_monthly %>%
  filter(series == "copper") %>%
  transmute(date = ceiling_date(floor_date(date, "month"), "month") - days(1),
            copper = value)

monthly_wide <- monthly_daily_src %>%
  full_join(monthly_copper, by = "date") %>%
  arrange(date)

# =========================================================
# 2) Transformation: 12m change / 10y rolling σ of 12m changes, winsorize ±3
# =========================================================
VAR_NAMES <- c("sp500", "slope", "oil", "copper",
               "tbill_3m", "volatility", "stock_bond")
LABELS <- c(sp500 = "S&P 500", slope = "Yield curve (10y–3m)",
            oil = "Oil (WTI)", copper = "Copper",
            tbill_3m = "US 3m yield", volatility = "Volatility",
            stock_bond = "Stock–bond corr")

transform_series <- function(x, diff_lag = 12, roll_win = 120, wins = 3) {
  d <- c(rep(NA_real_, diff_lag), diff(x, lag = diff_lag))
  s <- rollapplyr(d, width = roll_win,
                  FUN = function(v) if (sum(!is.na(v)) >= 0.8 * roll_win) sd(v, na.rm = TRUE) else NA_real_,
                  fill = NA_real_)
  z_raw <- d / s
  z_win <- pmin(pmax(z_raw, -wins), wins)
  tibble(z_raw = z_raw, z = z_win)
}

# Long frame of transformed variables (z_raw and z winsorized)
transformed_long <- monthly_wide %>%
  pivot_longer(-date, names_to = "series", values_to = "value") %>%
  filter(series %in% VAR_NAMES) %>%
  arrange(series, date) %>%
  group_by(series) %>%
  mutate(transform_series(value)) %>%
  ungroup()

# Wide frame of winsorized z-scores only
state_z <- transformed_long %>%
  select(date, series, z) %>%
  pivot_wider(names_from = series, values_from = z) %>%
  arrange(date) %>%
  # Keep columns in our preferred order
  select(date, all_of(VAR_NAMES))

# Keep only complete-case rows for distance calculation
state_full <- state_z %>% drop_na()
cat(sprintf("Full (all 7 vars) z-score data: %s to %s  (%d months)\n",
            format(min(state_full$date)), format(max(state_full$date)),
            nrow(state_full)))

# =========================================================
# 3) EXHIBIT 2 — raw state variables
# =========================================================
raw_long <- monthly_wide %>%
  pivot_longer(-date, names_to = "series", values_to = "value") %>%
  filter(series %in% VAR_NAMES, !is.na(value)) %>%
  mutate(series = factor(series, levels = VAR_NAMES, labels = LABELS[VAR_NAMES]))

p_ex2 <- ggplot(raw_long, aes(date, value)) +
  geom_line(colour = "#1f4e79") +
  facet_wrap(~series, scales = "free_y", ncol = 2) +
  labs(title = "Exhibit 2 — Raw economic state variables",
       x = NULL, y = NULL)
ggsave("plots/exhibit_02_raw.png", p_ex2, width = 10, height = 9, dpi = 160)

# =========================================================
# 4) EXHIBIT 3 — transformed (raw z vs winsorized)
# =========================================================
trans_long_plot <- transformed_long %>%
  filter(series %in% VAR_NAMES) %>%
  pivot_longer(c(z_raw, z), names_to = "type", values_to = "val") %>%
  filter(!is.na(val)) %>%
  mutate(series = factor(series, levels = VAR_NAMES, labels = LABELS[VAR_NAMES]),
         type   = factor(type, levels = c("z_raw", "z"),
                         labels = c("raw z-score", "winsorized")))

p_ex3 <- ggplot(trans_long_plot, aes(date, val, colour = type)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~series, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("raw z-score" = "#1f4e79",
                                 "winsorized"  = "#e07b00")) +
  labs(title = "Exhibit 3 — Transformed state variables (raw z vs winsorized)",
       x = NULL, y = "z-score", colour = NULL) +
  theme(legend.position = "bottom")
ggsave("plots/exhibit_03_transformed.png", p_ex3, width = 10, height = 9, dpi = 160)

# =========================================================
# 5) EXHIBIT 4 — autocorrelations + descriptive stats
# =========================================================
auto_at <- function(x, lag_months) {
  x <- x[!is.na(x)]
  if (length(x) <= lag_months + 5) return(NA_real_)
  cor(x[-seq_len(lag_months)], head(x, -lag_months))
}

stats_tbl <- transformed_long %>%
  filter(series %in% VAR_NAMES) %>%
  group_by(series) %>%
  summarise(
    `1-month`    = auto_at(z, 1),
    `3-month`    = auto_at(z, 3),
    `12-month`   = auto_at(z, 12),
    `3-year`     = auto_at(z, 36),
    `10-year`    = auto_at(z, 120),
    `monthly mean` = mean(z, na.rm = TRUE),
    `std`        = sd(z, na.rm = TRUE),
    frequency    = "monthly",
    .groups = "drop"
  ) %>%
  mutate(series = LABELS[series]) %>%
  arrange(factor(series, levels = LABELS[VAR_NAMES]))

print(stats_tbl, n = Inf)
write.csv(stats_tbl, "plots/exhibit_04_autocorrelations.csv", row.names = FALSE)

# =========================================================
# 6) EXHIBIT 5 — cross-correlation heatmap
# =========================================================
cor_mat <- cor(state_full %>% select(-date), use = "pairwise.complete.obs")
cor_df  <- as_tibble(cor_mat, rownames = "v1") %>%
  pivot_longer(-v1, names_to = "v2", values_to = "rho") %>%
  mutate(across(c(v1, v2), ~ factor(LABELS[.], levels = rev(LABELS[VAR_NAMES]))))

p_ex5 <- ggplot(cor_df, aes(v2, v1, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.2) +
  scale_fill_gradient2(low = "#3a7dca", mid = "#f6f6f6", high = "#c0392b",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Exhibit 5 — Cross-correlations of state variables",
       x = NULL, y = NULL, fill = "Corr") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/exhibit_05_cross_corr.png", p_ex5, width = 7.5, height = 6.5, dpi = 160)

# =========================================================
# 7) Distance matrix
# =========================================================
X <- as.matrix(state_full %>% select(-date))
dates <- state_full$date
N <- nrow(X)
D <- as.matrix(dist(X, method = "euclidean"))
rownames(D) <- colnames(D) <- as.character(dates)

# =========================================================
# 8) EXHIBIT 6–8 — similarity plots for specific anchor dates
# =========================================================
plot_similarity <- function(anchor_date, title_str,
                            top_pct = 0.15, mask_months = 36,
                            extra_shading = NULL) {
  anchor <- as.Date(anchor_date)
  if (!(anchor %in% dates)) {
    # take nearest month-end
    anchor <- dates[which.min(abs(dates - anchor))]
  }
  T_idx <- which(dates == anchor)
  dvec  <- D[T_idx, ]
  eligible <- seq_len(T_idx - mask_months - 1)
  if (length(eligible) < 20) stop("Not enough history before ", anchor)

  n_sel <- ceiling(length(eligible) * top_pct)
  sel_idx <- eligible[order(dvec[eligible])][seq_len(n_sel)]

  plot_df <- tibble(date = dates, d = dvec) %>%
    mutate(idx = row_number(),
           masked   = idx > (T_idx - mask_months - 1) & idx <= T_idx,
           selected = idx %in% sel_idx)

  ggplot(plot_df, aes(date, d)) +
    {if (!is.null(extra_shading))
      lapply(extra_shading, function(s)
        geom_rect(data = tibble(), aes(xmin = s$start, xmax = s$end,
                                       ymin = -Inf, ymax = Inf),
                  fill = "#f4d58d", alpha = 0.35, inherit.aes = FALSE))} +
    geom_rect(data = filter(plot_df, selected),
              aes(xmin = date - 15, xmax = date + 15, ymin = -Inf, ymax = Inf),
              fill = "#7bc96f", alpha = 0.6, inherit.aes = FALSE) +
    geom_rect(data = filter(plot_df, masked),
              aes(xmin = date - 15, xmax = date + 15, ymin = -Inf, ymax = Inf),
              fill = "grey80", alpha = 0.5, inherit.aes = FALSE) +
    geom_line(linetype = "dashed", colour = "#1f4e79") +
    labs(title = title_str, x = NULL, y = "Global score")
}

p_ex6 <- plot_similarity("2009-01-31",
                         "Exhibit 6 — Global similarity to January 2009 (GFC)")
ggsave("plots/exhibit_06_jan2009.png", p_ex6, width = 10, height = 4.5, dpi = 160)

p_ex7a <- plot_similarity("2020-02-29",
                          "Exhibit 7a — Similarity to February 2020 (COVID)")
p_ex7b <- plot_similarity("2020-04-30",
                          "Exhibit 7b — Similarity to April 2020 (COVID)")
ggsave("plots/exhibit_07_covid.png", p_ex7a / p_ex7b, width = 10, height = 8, dpi = 160)

# Inflation shading for Exhibit 8
inflation_periods <- list(
  list(start = as.Date("1941-04-01"), end = as.Date("1942-05-31")),
  list(start = as.Date("1946-03-01"), end = as.Date("1947-03-31")),
  list(start = as.Date("1950-08-01"), end = as.Date("1951-02-28")),
  list(start = as.Date("1966-02-01"), end = as.Date("1970-01-31")),
  list(start = as.Date("1972-07-01"), end = as.Date("1974-12-31")),
  list(start = as.Date("1977-02-01"), end = as.Date("1980-03-31")),
  list(start = as.Date("1987-02-01"), end = as.Date("1990-11-30")),
  list(start = as.Date("2007-09-01"), end = as.Date("2008-07-31")),
  list(start = as.Date("2022-01-01"), end = as.Date("2022-10-31"))
)
p_ex8 <- plot_similarity("2022-08-31",
                         "Exhibit 8 — Similarity to August 2022 (inflation surge)",
                         extra_shading = inflation_periods)
ggsave("plots/exhibit_08_aug2022.png", p_ex8, width = 10, height = 4.5, dpi = 160)

# =========================================================
# 9) EXHIBIT 9 — EWMA of global scores, multiple lookbacks
# =========================================================
# For each T, C_T = EWMA_{t<=T} of d_{T,t}, higher weight on recent t.
ewma_of_row <- function(d_row, n_months) {
  beta <- 1 - 1/n_months
  # weights: newest date (index = T) gets weight 1, previous gets β, etc.
  # here d_row is distances to T; the "time index" of those distances is position i
  # We treat the distance's own position-in-time as the weighting variable,
  # so the most recent entry gets the highest weight.
  T_idx <- length(d_row)
  w <- beta ^ ((T_idx - 1):0)
  sum(d_row * w, na.rm = TRUE) / sum(w[!is.na(d_row)])
}

ewma_series <- function(n_months) {
  sapply(seq_len(N), function(T_idx) ewma_of_row(D[T_idx, 1:T_idx], n_months))
}

ewma_df <- tibble(
  date        = dates,
  `1-year`    = ewma_series(12),
  `2-year`    = ewma_series(24),
  `3-year`    = ewma_series(36),
  `4-year`    = ewma_series(48)
) %>%
  mutate(Mean = rowMeans(select(., -date), na.rm = TRUE)) %>%
  pivot_longer(-date, names_to = "lookback", values_to = "ewma")

p_ex9 <- ggplot(ewma_df, aes(date, ewma, colour = lookback)) +
  geom_line(linewidth = 0.45) +
  scale_colour_manual(values = c(`1-year` = "#a2c5f2", `2-year` = "#5b9bd5",
                                 `3-year` = "#1f4e79", `4-year` = "#0b2d4a",
                                 Mean = "#c0392b")) +
  labs(title = "Exhibit 9 — EWMA of global similarity score",
       x = NULL, y = "EWMA", colour = "Lookback")
ggsave("plots/exhibit_09_ewma.png", p_ex9, width = 10, height = 5, dpi = 160)

# =========================================================
# 10) Factor-timing strategy
# =========================================================
# Align Fama-French factors to state dates
FACTOR_NAMES <- c("market", "size", "value", "profitability", "investment", "momentum")
rf_aligned <- raw$ff_factors %>%
  mutate(date = ceiling_date(floor_date(date, "month"), "month") - days(1))

# Join factor rets to state dates
state_with_rets <- state_full %>%
  left_join(rf_aligned %>% select(date, all_of(FACTOR_NAMES)), by = "date")

# Keep only months where both state and all factors are available
state_full2 <- state_with_rets %>% drop_na()
X2 <- as.matrix(state_full2 %>% select(all_of(VAR_NAMES)))
F2 <- as.matrix(state_full2 %>% select(all_of(FACTOR_NAMES)))
dates2 <- state_full2$date
N2 <- nrow(X2)
D2 <- as.matrix(dist(X2))

cat(sprintf("Factor-timing sample: %s to %s (%d months)\n",
            format(min(dates2)), format(max(dates2)), N2))

# ---------------------------------------------------------
# Core: for each month T, for each factor, quintile position
# Position at T+1 = sign(mean of r[i+1] over selected i)
# i eligible: 1:(T - exclude - 1); selected by distance quintile
# ---------------------------------------------------------
compute_quintile_positions <- function(D, F, dates,
                                       n_quantiles = 5,
                                       exclude_months = 36) {
  N <- nrow(D); K <- ncol(F)
  # positions[T, q, k] -> sign for month T+1
  pos <- array(NA_real_, dim = c(N, n_quantiles, K))
  for (T_idx in seq_len(N - 1)) {
    max_i <- T_idx - exclude_months - 1
    if (max_i < 20) next
    eligible <- seq_len(max_i)
    ord <- eligible[order(D[T_idx, eligible])]
    per_q <- ceiling(length(ord) / n_quantiles)
    for (q in seq_len(n_quantiles)) {
      start <- (q - 1) * per_q + 1
      end   <- min(q * per_q, length(ord))
      if (start > length(ord)) next
      sel <- ord[start:end]
      subs <- F[sel + 1, , drop = FALSE]    # subsequent 1m returns
      mu <- colMeans(subs, na.rm = TRUE)
      pos[T_idx, q, ] <- sign(mu)
    }
  }
  dimnames(pos) <- list(as.character(dates),
                        paste0("Q", seq_len(n_quantiles)),
                        colnames(F))
  pos
}

positions <- compute_quintile_positions(D2, F2, dates2, n_quantiles = 5)

# Strategy returns: pos[T,q,k] * F[T+1, k] => return realised at T+1
# Build a long tibble of returns per quintile/factor
quintile_returns <- function(pos, F) {
  N <- dim(pos)[1]; Q <- dim(pos)[2]; K <- dim(pos)[3]
  out <- list()
  for (q in seq_len(Q)) {
    sgn <- pos[, q, , drop = FALSE]
    # r[T+1] * sign[T]; align by shifting
    sig_mat <- sgn[1:(N-1), 1, ]
    fwd_ret <- F[2:N, ]
    r <- sig_mat * fwd_ret
    out[[q]] <- tibble(date = dimnames(pos)[[1]][2:N]) %>%
      mutate(date = as.Date(date)) %>%
      bind_cols(as_tibble(r)) %>%
      mutate(quintile = paste0("Q", q))
  }
  bind_rows(out)
}

qrets <- quintile_returns(positions, F2)

# Aggregate: equal-weight across 6 factors per quintile
qrets_agg <- qrets %>%
  rowwise() %>%
  mutate(ew = mean(c_across(all_of(FACTOR_NAMES)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(date, quintile, ew)

# Long-only benchmark: +1 on each factor
lo_ret <- as_tibble(F2) %>%
  mutate(date = dates2) %>%
  rowwise() %>%
  mutate(ew = mean(c_across(all_of(FACTOR_NAMES)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(date, ew) %>%
  mutate(quintile = "LO")

all_rets <- bind_rows(qrets_agg, lo_ret)

# Cumulative (log-scale-friendly) paths
cum_df <- all_rets %>%
  filter(!is.na(ew)) %>%
  arrange(quintile, date) %>%
  group_by(quintile) %>%
  mutate(cum = cumsum(ew) * 100) %>%   # in %
  ungroup()

# Sharpe helper
ann_sharpe <- function(r) {
  r <- r[!is.na(r)]
  mean(r) / sd(r) * sqrt(12)
}
sr_by_q <- all_rets %>%
  filter(!is.na(ew)) %>%
  group_by(quintile) %>%
  summarise(SR = ann_sharpe(ew), .groups = "drop")
print(sr_by_q)

# =========================================================
# 11) EXHIBIT 10 — quintile cumulative + 1st-5th spread
# =========================================================
q_labels <- sr_by_q %>%
  mutate(label = sprintf("%s (SR %.2f)", quintile, SR))

p_ex10a <- cum_df %>%
  left_join(q_labels, by = "quintile") %>%
  ggplot(aes(date, cum, colour = label)) +
  geom_line() +
  labs(title = "Exhibit 10 — Quintile cumulative returns (equal-weighted across 6 factors)",
       x = NULL, y = "Cumulative return (%)", colour = NULL)

spread <- all_rets %>%
  filter(quintile %in% c("Q1", "Q5")) %>%
  select(date, quintile, ew) %>%
  pivot_wider(names_from = quintile, values_from = ew) %>%
  mutate(spread = Q1 - Q5) %>%
  arrange(date) %>%
  mutate(cum = cumsum(replace_na(spread, 0)) * 100)

sr_spread <- ann_sharpe(spread$spread)

p_ex10b <- ggplot(spread, aes(date, cum)) +
  geom_line(colour = "#1f4e79") +
  labs(title = sprintf("Exhibit 10 — Q1 minus Q5 spread (SR %.2f)", sr_spread),
       x = NULL, y = "Cumulative return (%)")

ggsave("plots/exhibit_10_quintiles.png", p_ex10a / p_ex10b,
       width = 10, height = 8, dpi = 160)

# =========================================================
# 12) EXHIBIT 11 — drawdowns: LO vs Q1-Q5 model
# =========================================================
dd <- function(r) {
  r[is.na(r)] <- 0
  eq <- cumprod(1 + r)
  eq / cummax(eq) - 1
}

dd_df <- bind_rows(
  tibble(date = dates2, dd = dd(lo_ret$ew), series = "Long-only"),
  tibble(date = spread$date, dd = dd(spread$spread), series = "Q1 − Q5 model")
)

p_ex11 <- ggplot(dd_df, aes(date, dd * 100, colour = series)) +
  geom_line() +
  scale_colour_manual(values = c("Long-only" = "#888", "Q1 − Q5 model" = "#1f4e79")) +
  labs(title = "Exhibit 11 — Drawdown comparison",
       x = NULL, y = "% of capital", colour = NULL)
ggsave("plots/exhibit_11_drawdown.png", p_ex11, width = 10, height = 4.5, dpi = 160)

# =========================================================
# 13) EXHIBIT 12 — robustness to number of quantiles (top - bottom)
# =========================================================
top_minus_bottom <- function(n_q) {
  pos_q <- compute_quintile_positions(D2, F2, dates2, n_quantiles = n_q)
  qr <- quintile_returns(pos_q, F2)
  qr_ew <- qr %>%
    rowwise() %>% mutate(ew = mean(c_across(all_of(FACTOR_NAMES)), na.rm = TRUE)) %>%
    ungroup() %>% select(date, quintile, ew)
  top <- qr_ew %>% filter(quintile == "Q1") %>% select(date, top = ew)
  bot <- qr_ew %>% filter(quintile == paste0("Q", n_q)) %>% select(date, bot = ew)
  inner_join(top, bot, by = "date") %>%
    mutate(r = top - bot,
           n_q = paste0(n_q, " quantiles"),
           cum = cumsum(replace_na(r, 0)) * 100)
}

ex12_df <- map_dfr(c(2, 3, 4, 5, 10, 20), top_minus_bottom)
ex12_sr <- ex12_df %>% group_by(n_q) %>% summarise(SR = ann_sharpe(r), .groups = "drop")
ex12_df <- ex12_df %>% left_join(ex12_sr, by = "n_q") %>%
  mutate(label = sprintf("%s  (SR %.2f)", n_q, SR))

p_ex12 <- ggplot(ex12_df, aes(date, cum, colour = label)) +
  geom_line() +
  labs(title = "Exhibit 12 — Top−bottom spread for various quantile choices",
       x = NULL, y = "Cumulative return (%)", colour = NULL)
ggsave("plots/exhibit_12_quantile_robust.png", p_ex12, width = 10, height = 5, dpi = 160)

# =========================================================
# 14) EXHIBIT 13 — robustness to z-score lookback window
# =========================================================
build_state_for_lookback <- function(years) {
  wide <- monthly_wide %>%
    pivot_longer(-date, names_to = "series", values_to = "value") %>%
    filter(series %in% VAR_NAMES) %>%
    arrange(series, date) %>%
    group_by(series) %>%
    mutate(transform_series(value, roll_win = years * 12)) %>%
    ungroup() %>%
    select(date, series, z) %>%
    pivot_wider(names_from = series, values_from = z) %>%
    select(date, all_of(VAR_NAMES)) %>%
    drop_na()
  wide
}

run_pipeline_for_lookback <- function(years) {
  sf <- build_state_for_lookback(years) %>%
    left_join(rf_aligned %>% select(date, all_of(FACTOR_NAMES)), by = "date") %>%
    drop_na()
  if (nrow(sf) < 60) return(NULL)
  Xl <- as.matrix(sf %>% select(all_of(VAR_NAMES)))
  Fl <- as.matrix(sf %>% select(all_of(FACTOR_NAMES)))
  Dl <- as.matrix(dist(Xl))
  posl <- compute_quintile_positions(Dl, Fl, sf$date, n_quantiles = 5)
  qrl <- quintile_returns(posl, Fl) %>%
    rowwise() %>% mutate(ew = mean(c_across(all_of(FACTOR_NAMES)), na.rm = TRUE)) %>%
    ungroup() %>% select(date, quintile, ew)
  top <- qrl %>% filter(quintile == "Q1") %>% select(date, top = ew)
  bot <- qrl %>% filter(quintile == "Q5") %>% select(date, bot = ew)
  inner_join(top, bot, by = "date") %>%
    mutate(r = top - bot,
           lookback = paste0(years, "-year"),
           cum = cumsum(replace_na(r, 0)) * 100)
}

ex13_df <- map_dfr(c(1, 3, 5), run_pipeline_for_lookback)
p_ex13 <- ggplot(ex13_df, aes(date, cum, colour = lookback)) +
  geom_line() +
  labs(title = "Exhibit 13 — Robustness to z-score lookback",
       x = NULL, y = "Cumulative return (%)", colour = NULL)
ggsave("plots/exhibit_13_lookback.png", p_ex13, width = 10, height = 5, dpi = 160)

# =========================================================
# 15) EXHIBIT 1 — yearly returns of Q1−Q5 model at 15% vol target
# =========================================================
r_model <- spread$spread
target_vol <- 0.15
ann_vol <- sd(r_model, na.rm = TRUE) * sqrt(12)
scale_k <- target_vol / ann_vol
r_scaled <- r_model * scale_k

yearly <- tibble(date = spread$date, r = r_scaled) %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(yr = prod(1 + replace_na(r, 0)) - 1, .groups = "drop")

p_ex1 <- ggplot(yearly, aes(year, yr * 100)) +
  geom_col(fill = "#1f77b4") +
  labs(title = "Exhibit 1 — Yearly returns at 15% vol target",
       x = "Year", y = "% return")
ggsave("plots/exhibit_01_yearly.png", p_ex1, width = 10, height = 4.5, dpi = 160)

# Descriptive stats that match the paper's text
pos_yrs <- mean(yearly$yr > 0)
mu_pos  <- mean(yearly$yr[yearly$yr > 0]) * 100
mu_neg  <- mean(yearly$yr[yearly$yr < 0]) * 100
cat(sprintf("Positive years: %.0f%%  |  mean(up): %.1f%%  mean(down): %.1f%%\n",
            pos_yrs * 100, mu_pos, mu_neg))

# =========================================================
# 16) APPENDIX A1 — per-factor quintile analysis
# =========================================================
per_factor_q <- qrets %>%
  pivot_longer(all_of(FACTOR_NAMES), names_to = "factor", values_to = "ret") %>%
  arrange(factor, quintile, date) %>%
  group_by(factor, quintile) %>%
  mutate(cum = cumsum(replace_na(ret, 0)) * 100) %>%
  ungroup()

# Add LO per factor: +1 on that factor always
lo_per_factor <- as_tibble(F2) %>%
  mutate(date = dates2) %>%
  pivot_longer(-date, names_to = "factor", values_to = "ret") %>%
  group_by(factor) %>%
  mutate(cum = cumsum(replace_na(ret, 0)) * 100,
         quintile = "LO") %>%
  ungroup()

p_A1 <- bind_rows(per_factor_q %>% select(date, factor, quintile, cum),
                  lo_per_factor %>% select(date, factor, quintile, cum)) %>%
  ggplot(aes(date, cum, colour = quintile)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~factor, scales = "free_y", ncol = 3) +
  labs(title = "Appendix A1 — Individual factor analysis by quintile",
       x = NULL, y = "Cumulative return (%)", colour = NULL)
ggsave("plots/appendix_A1_per_factor_quintiles.png", p_A1,
       width = 11, height = 7, dpi = 160)

# Appendix A2: per-factor Q1 - Q5 vs per-factor LO
p_A2_df <- qrets %>%
  pivot_longer(all_of(FACTOR_NAMES), names_to = "factor", values_to = "ret") %>%
  select(date, factor, quintile, ret) %>%
  pivot_wider(names_from = quintile, values_from = ret) %>%
  mutate(spread = Q1 - Q5) %>%
  group_by(factor) %>%
  arrange(date) %>%
  mutate(cum_spread = cumsum(replace_na(spread, 0)) * 100) %>%
  ungroup()

p_A2 <- ggplot() +
  geom_line(data = p_A2_df, aes(date, cum_spread, colour = "Q1 − Q5"), linewidth = 0.4) +
  geom_line(data = lo_per_factor, aes(date, cum, colour = "Long-only"),
            linewidth = 0.4, linetype = "dashed") +
  facet_wrap(~factor, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("Q1 − Q5" = "#1f4e79", "Long-only" = "grey40")) +
  labs(title = "Appendix A2 — Per-factor similarity minus dissimilarity",
       x = NULL, y = "Cumulative return (%)", colour = NULL)
ggsave("plots/appendix_A2_per_factor_spread.png", p_A2,
       width = 11, height = 7, dpi = 160)

message("Done. See ./plots/ for all exhibits.")