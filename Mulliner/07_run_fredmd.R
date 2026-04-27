# ---------------------------------------------------------------------------
# R/07_run_fredmd.R
#
# FRED-MD extension: runs the Mulliner regime strategy with the seven hand-
# picked financial state variables replaced by k principal components of the
# FRED-MD panel. Everything else (distance metric, mask, quintiles, factor
# universe, vol targeting) is identical to the baseline.
#
# Important methodological note on the z-score step:
#
# In the baseline, each raw state variable goes through:
#   level -> 12-month difference -> divide by 10-yr rolling std -> winsorize
#
# That transformation is what turns the seven non-stationary financial series
# into something near-unit-variance and zero-mean. The PCA factor scores from
# EM-PCA on FRED-MD are ALREADY standardised (the algorithm standardises every
# series before SVD), and they are extracted FROM ALREADY-STATIONARY series
# (the McCracken-Ng tcodes deliver stationarity per series). They have, by
# construction, mean ~ 0 and variance ~ 1 over the estimation window.
#
# So the PCA factors are plugged directly into the distance metric without
# the 12m-diff / rolling-std / winsorize step. They are still winsorized at
# +-3 just to be safe against extreme tail draws.
#
# This is the cleanest interpretation of "all else equal": each side of the
# comparison delivers state variables that satisfy the same input conditions
# (~stationary, ~unit variance, near-zero mean) to the distance metric.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(slider)
  library(patchwork)
})

# ---- Configuration: identical to main.R except for the variables block ----

CFG_FREDMD <- list(
  paths = list(
    fredmd_csv = "data/2025-12-md.csv",
    ff5_url    = "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_5_Factors_2x3_CSV.zip",
    mom_url    = "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Momentum_Factor_CSV.zip",
    output_dir = "output/fredmd"
  ),

  sample_end = as.Date("2025-12-31"),

  # PCA configuration
  r_target  = 8L,            # # of PCA factors fed into distance metric
                             # (Bai-Ng on FRED-MD typically selects 8)
  r_max     = 12L,           # max for Bai-Ng selection (only used when
                             # selection = "bai_ng")
  selection = "fixed",       # "fixed" (12x faster) or "bai_ng"
  ic_type   = "ICp2",        # which IC to use under "bai_ng"
  min_obs   = 120L,          # months of history before first PCA estimate

  # If TRUE, apply Mulliner's 12m-diff + 10y-rolling-std + winsorise to
  # the PCA factor SCORES before they enter the distance metric. This makes
  # the FRED-MD path methodologically identical to the baseline (just with
  # different inputs). If FALSE, we use the PCA scores directly (already
  # standardised by EM-PCA), as in our initial run.
  apply_rolling_zscore = TRUE,
  zscore_horizon       = 12L,
  zscore_window        = 120L,

  # Apply +-3 winsorisation to the (post-z-score, if applicable) PCA scores
  winsorize_at = 3,

  # Caching: skip rolling PCA recomputation if cache exists. Set to FALSE
  # to force a fresh recomputation (e.g., after changing r_target).
  use_pca_cache = TRUE,

  # Strategy parameters - identical to baseline
  mask_months      = 36L,
  paper_top_pct    = 0.15,
  quantile_default = 5L,
  quantile_robust  = c(2, 3, 4, 5, 10, 20),
  vol_target       = 0.15,
  vol_window       = 36L,
  ewma_lookbacks   = c(12L, 24L, 36L, 48L)
)

dir.create(CFG_FREDMD$paths$output_dir, showWarnings = FALSE, recursive = TRUE)

# CFG must be set so utility functions in the existing pipeline find their
# output path. We point CFG to the FRED-MD output dir for this run.
CFG <- CFG_FREDMD

# ---- Source the existing pipeline modules (we reuse most of them) --------
source("R/utils.R")
source("R/01_data.R")          # for load_factors() - the FF download
source("R/02_similarity.R")    # for compute_distance_matrix() etc.
source("R/03_strategy.R")      # for run_strategies(), bucket_by_similarity()
source("R/04_exhibits.R")      # for the exhibit functions

# ---- Source the FRED-MD specific modules ---------------------------------
source("R/05_fredmd_data.R")
source("R/06_pca_factors.R")

# ---- Pipeline ------------------------------------------------------------

message("\n[1/7] Loading FRED-MD CSV ...")
fmd <- load_fredmd(CFG_FREDMD$paths$fredmd_csv)

message("[2/7] Applying McCracken-Ng transformations ...")
xform <- transform_fredmd(fmd, sample_end = CFG_FREDMD$sample_end)
describe_fredmd(fmd, xform)

# Optional: outlier removal per McCracken-Ng (10x IQR from median).
# We apply this because BN factor selection is sensitive to extreme draws.
message("[3/7] Removing outliers (>10x IQR from median, set to NA) ...")
xform_clean <- remove_outliers_mn(xform, factor = 10)
cat("  Replaced ", attr(xform_clean, "n_outliers"), " values with NA.\n")

message("[4/7] Rolling EM-PCA factor extraction ...")
pca_cache_path <- file.path(CFG$paths$output_dir,
  sprintf("pca_rolling_r%d_%s.rds",
          CFG_FREDMD$r_target, CFG_FREDMD$selection))

if (CFG_FREDMD$use_pca_cache && file.exists(pca_cache_path)) {
  message("  loading cached rolling PCA from ", pca_cache_path)
  pca_out <- readRDS(pca_cache_path)
} else {
  message("  no cache found at ", pca_cache_path, " - computing now")
  message("  (this takes ~5 min with selection='fixed', ~1 hour with 'bai_ng')")
  pca_out <- rolling_pca_factors(
    transformed = xform_clean,
    r_target    = CFG_FREDMD$r_target,
    r_max       = CFG_FREDMD$r_max,
    selection   = CFG_FREDMD$selection,
    ic_type     = CFG_FREDMD$ic_type,
    min_obs     = CFG_FREDMD$min_obs,
    verbose     = TRUE
  )
  saveRDS(pca_out, pca_cache_path)
  message("  cached to ", pca_cache_path)
}

# Save the time series of chosen r (one diagnostic Mulliner doesn't have)
r_chosen <- attr(pca_out, "r_chosen")
if (!is.null(r_chosen)) {
  write_csv(r_chosen, file.path(CFG$paths$output_dir, "fredmd_r_chosen.csv"))
}

# Also save a static full-sample factor extraction for interpretation only
message("[5/7] Static full-sample PCA for factor interpretation ...")
static_out <- static_pca_factors(xform_clean, r = CFG_FREDMD$r_target)

# Compute mR^2 per factor for each variable - which series load on which factor
loadings <- static_out$loadings
F_full   <- static_out$factors
var_names <- static_out$var_names
F_mat <- as.matrix(F_full[, paste0("F", seq_len(CFG_FREDMD$r_target))])
X_full <- as.matrix(xform_clean[, var_names])
X_full_std <- standardise_cols(X_full)$X
X_full_std[is.na(X_full_std)] <- 0

# For each (variable, factor) pair, R^2 of regressing variable on factor
# (factors are already orthogonal in the ref fit, so each marginal R^2 is
# additive)
R2_table <- matrix(0, length(var_names), CFG_FREDMD$r_target,
                   dimnames = list(var_names,
                                   paste0("F", seq_len(CFG_FREDMD$r_target))))
for (j in seq_len(CFG_FREDMD$r_target)) {
  for (i in seq_along(var_names)) {
    fit <- lm(X_full_std[, i] ~ F_mat[, j] - 1)
    R2_table[i, j] <- 1 - sum(residuals(fit)^2) / sum(X_full_std[, i]^2)
  }
}
mR2 <- colMeans(R2_table)
write_csv(as_tibble(R2_table, rownames = "variable"),
          file.path(CFG$paths$output_dir, "fredmd_factor_r2_per_variable.csv"))
write_csv(tibble(factor = paste0("F", seq_along(mR2)), mR2 = mR2),
          file.path(CFG$paths$output_dir, "fredmd_factor_mR2.csv"))

cat("\nMean R^2 per factor (avg across all variables):\n")
print(round(mR2, 3))
cat("\nTop 10 variables by |loading| on factor 1:\n")
print(head(R2_table[order(-R2_table[, 1]), 1, drop = FALSE], 10))

# ---- Wire PCA factors into the existing strategy pipeline ----------------
#
# The existing pipeline expects two objects:
#   state$transformed_winsorized : tibble with `date` + variable columns
#                                  (z-scored, winsorized, ready for distance)
#   state$monthly                : raw monthly state for re-running variants
#
# We don't have the latter for the FRED-MD path (and don't need it because
# the lookback robustness sweep doesn't apply: the PCA standardises by
# construction, no rolling-std denominator to vary). We pass through only
# what's needed.

message("[6/7] Building distance matrix from PCA factors ...")
factor_cols <- paste0("F", seq_len(CFG_FREDMD$r_target))

if (CFG_FREDMD$apply_rolling_zscore) {
  message("  Applying Mulliner-style rolling z-score to PCA factor scores")
  message("  (12-month diff / ", CFG_FREDMD$zscore_window,
          "-month rolling std / +-",
          CFG_FREDMD$winsorize_at, " winsorize)")

  # Reuse the existing build_state_variables() machinery. We treat each
  # PCA factor F_k as a "level" series and apply the same transformation
  # Mulliner applies to his 7 raw variables.
  fake_variables <- tibble(var = factor_cols, label = factor_cols,
                           diff_type = "level")
  pca_state <- build_state_variables(
    macro_m   = pca_out,
    variables = fake_variables,
    horizon   = CFG_FREDMD$zscore_horizon,
    window    = CFG_FREDMD$zscore_window,
    winsorize = CFG_FREDMD$winsorize_at
  )
  state_fmd <- list(
    transformed            = pca_state$transformed,
    transformed_winsorized = pca_state$transformed_winsorized,
    raw_pca                = pca_out
  )
} else {
  message("  Skipping rolling z-score: PCA scores already standardised by EM-PCA")
  state_fmd <- list(
    transformed            = pca_out,
    transformed_winsorized = pca_out %>%
      mutate(across(all_of(factor_cols),
                    ~ winsorize(., bound = CFG_FREDMD$winsorize_at))),
    raw_pca                = pca_out
  )
}

dist_mat_fmd <- compute_distance_matrix(state_fmd$transformed_winsorized)

# ---- Load factor returns -------------------------------------------------
message("[7/7] Loading FF factors and running strategies ...")
factors <- load_factors(CFG_FREDMD$paths$ff5_url,
                        CFG_FREDMD$paths$mom_url,
                        CFG_FREDMD$sample_end)

# run_strategies expects `macro_monthly` and `variables` for the lookback
# robustness rerun. The lookback sweep doesn't apply on the PCA path
# (no rolling-std denominator to vary), so we pass NULL and skip that
# robustness in this run. We achieve this by passing an empty
# `lookback_robust`.

# We need a dummy `variables` tibble for the robustness rerun - won't be used
# because lookback_robust is empty.
dummy_variables <- tibble(var = factor_cols, label = factor_cols,
                          diff_type = "level")

strat_fmd <- run_strategies(
  distance_matrix = dist_mat_fmd,
  factor_returns  = factors,
  mask_months     = CFG_FREDMD$mask_months,
  q_default       = CFG_FREDMD$quantile_default,
  q_robust        = CFG_FREDMD$quantile_robust,
  lookback_robust = integer(0),                # skip lookback sweep for PCA
  macro_monthly   = state_fmd$transformed_winsorized,
  variables       = dummy_variables,
  horizon         = 12L,
  winsorize       = CFG_FREDMD$winsorize_at
)

# ---- Exhibits ------------------------------------------------------------

message("Producing FRED-MD exhibits ...")

# Exhibit 4 equivalent: autocorrelations and cross-correlations of factors
cat("\nFRED-MD factor autocorrelations:\n")
acf_at <- function(x, k) {
  x <- x[!is.na(x)]
  if (length(x) < k + 5) return(NA_real_)
  a <- acf(x, lag.max = k, plot = FALSE, na.action = na.pass)
  as.numeric(a$acf[k + 1])
}
acf_table <- map_dfr(factor_cols, function(v) {
  x <- state_fmd$transformed_winsorized[[v]]
  tibble(
    factor = v,
    `1m`  = acf_at(x, 1),
    `3m`  = acf_at(x, 3),
    `12m` = acf_at(x, 12),
    `3y`  = acf_at(x, 36),
    mean  = mean(x, na.rm = TRUE),
    std   = sd(x, na.rm = TRUE)
  )
})
print(acf_table %>% mutate(across(where(is.numeric), ~ round(.x, 2))))
write_csv(acf_table, file.path(CFG$paths$output_dir,
                                "fredmd_factor_autocorrelations.csv"))

# --- Exhibit 10: quintile + LO cumulative return -------------------------
# Generate both unlevered and vol-targeted versions, mirroring the baseline.

# Helper: vol-target a return series (ex-ante, trailing window) to 15% ann.
vol_target_simple <- function(r, target = CFG_FREDMD$vol_target,
                               win = CFG_FREDMD$vol_window) {
  rv <- slider::slide_dbl(r, ~ sd(.x, na.rm = TRUE),
                          .before = win, .after = -1, .complete = FALSE)
  scale <- (target / sqrt(12)) / rv
  scale[!is.finite(scale)] <- NA
  r * scale
}

# Determine the start date once
start_date <- strat_fmd$per_bucket %>%
  filter(quintile == "Q1", is.finite(avg)) %>%
  summarise(d = min(date)) %>% pull(d)

# --- Unlevered ---
qs <- strat_fmd$per_bucket %>%
  filter(is.finite(avg), date >= start_date) %>%
  group_by(quintile) %>% mutate(cum = cum_pct(avg)) %>% ungroup()
lo <- strat_fmd$long_only %>%
  filter(is.finite(avg), date >= start_date) %>%
  mutate(cum = cum_pct(avg))

stats_q <- qs %>%
  group_by(quintile) %>%
  summarise(SR = sharpe(avg),
            corr = cor(avg, lo$avg[match(date, lo$date)],
                       use = "pairwise"),
            .groups = "drop")
cat("\n[unlevered] FRED-MD quintile stats:\n")
print(stats_q)
write_csv(stats_q, file.path(CFG$paths$output_dir,
                              "fredmd_exhibit_10_quintile_stats.csv"))

p10 <- ggplot() +
  geom_line(data = qs, aes(date, cum, colour = quintile), linewidth = 0.6) +
  geom_line(data = lo, aes(date, cum), colour = "black",
            linetype = "dashed", linewidth = 0.6) +
  labs(title = "FRED-MD: similarity quintiles vs long-only (unlevered)",
       subtitle = paste0(CFG_FREDMD$r_target,
                         " PCA factors of FRED-MD (",
                         CFG_FREDMD$selection, " selection",
                         if (CFG_FREDMD$apply_rolling_zscore)
                           "; with rolling z-score" else
                           "; PCA scores direct",
                         ")"),
       x = NULL, y = "Cumulative return (%)", colour = NULL) +
  theme_regimes()
ggsave(file.path(CFG$paths$output_dir, "fredmd_exhibit_10_quintiles.png"),
       p10, width = 11, height = 5, dpi = 150)

# --- Vol-targeted ---
qs_vt <- strat_fmd$per_bucket %>%
  filter(is.finite(avg), date >= start_date) %>%
  group_by(quintile) %>%
  mutate(r_vt = vol_target_simple(avg)) %>%
  filter(is.finite(r_vt)) %>%
  mutate(cum = cum_pct(r_vt)) %>%
  ungroup()
lo_vt <- strat_fmd$long_only %>%
  filter(is.finite(avg), date >= start_date) %>%
  mutate(r_vt = vol_target_simple(avg)) %>%
  filter(is.finite(r_vt)) %>%
  mutate(cum = cum_pct(r_vt))

stats_q_vt <- qs_vt %>%
  group_by(quintile) %>%
  summarise(SR = sharpe(r_vt),
            corr = cor(r_vt, lo_vt$r_vt[match(date, lo_vt$date)],
                       use = "pairwise"),
            .groups = "drop")
cat("\n[vol-targeted] FRED-MD quintile stats:\n")
print(stats_q_vt)
write_csv(stats_q_vt,
          file.path(CFG$paths$output_dir,
                    "fredmd_exhibit_10_quintile_stats_voltgt.csv"))

p10_vt <- ggplot() +
  geom_line(data = qs_vt, aes(date, cum, colour = quintile), linewidth = 0.6) +
  geom_line(data = lo_vt, aes(date, cum), colour = "black",
            linetype = "dashed", linewidth = 0.6) +
  labs(title = "FRED-MD: similarity quintiles vs long-only (15% vol target)",
       subtitle = paste0(CFG_FREDMD$r_target,
                         " PCA factors of FRED-MD"),
       x = NULL, y = "Cumulative return (%)", colour = NULL) +
  theme_regimes()
ggsave(file.path(CFG$paths$output_dir,
                 "fredmd_exhibit_10_quintiles_voltgt.png"),
       p10_vt, width = 11, height = 5, dpi = 150)

# --- Exhibit 11: drawdowns, both versions ---
# unlevered
diff_dd <- tibble(date = strat_fmd$diff_default$date,
                  dd   = drawdown(strat_fmd$diff_default$diff),
                  what = "Similarity model (Q1 - Q5)")
lo_dd <- tibble(date = strat_fmd$long_only$date,
                dd   = drawdown(strat_fmd$long_only$avg),
                what = "Long-only factor model")
df_dd <- bind_rows(lo_dd, diff_dd) %>% filter(is.finite(dd))
p11 <- ggplot(df_dd, aes(date, dd * 100, colour = what)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Long-only factor model" = "grey50",
                                  "Similarity model (Q1 - Q5)" = "#1f4e79"),
                      name = NULL) +
  labs(title = "FRED-MD: drawdown profile (unlevered)",
       x = NULL, y = "Drawdown (%)") +
  theme_regimes()
ggsave(file.path(CFG$paths$output_dir, "fredmd_exhibit_11_drawdowns.png"),
       p11, width = 10, height = 4.5, dpi = 150)

# vol-targeted
diff_r_vt <- vol_target_simple(strat_fmd$diff_default$diff)
lo_r_vt   <- vol_target_simple(strat_fmd$long_only$avg)
diff_dd_vt <- tibble(date = strat_fmd$diff_default$date,
                     dd   = drawdown(diff_r_vt),
                     what = "Similarity model (Q1 - Q5)")
lo_dd_vt   <- tibble(date = strat_fmd$long_only$date,
                     dd   = drawdown(lo_r_vt),
                     what = "Long-only factor model")
df_dd_vt <- bind_rows(lo_dd_vt, diff_dd_vt) %>% filter(is.finite(dd))
p11_vt <- ggplot(df_dd_vt, aes(date, dd * 100, colour = what)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("Long-only factor model" = "grey50",
                                  "Similarity model (Q1 - Q5)" = "#1f4e79"),
                      name = NULL) +
  labs(title = "FRED-MD: drawdown profile (15% vol target)",
       x = NULL, y = "Drawdown (%)") +
  theme_regimes()
ggsave(file.path(CFG$paths$output_dir,
                 "fredmd_exhibit_11_drawdowns_voltgt.png"),
       p11_vt, width = 10, height = 4.5, dpi = 150)

# Per-factor Q1-Q5 Sharpes
per_fac <- strat_fmd$per_factor %>%
  filter(is.finite(r), quintile %in% c("Q1", "Q5")) %>%
  pivot_wider(names_from = quintile, values_from = r) %>%
  mutate(diff = .data[["Q1"]] - .data[["Q5"]]) %>%
  filter(is.finite(diff))
stats_pf <- per_fac %>%
  group_by(factor) %>%
  summarise(SR = sharpe(diff), .groups = "drop")
cat("\nFRED-MD per-factor Q1-Q5 Sharpes:\n")
print(stats_pf)
write_csv(stats_pf, file.path(CFG$paths$output_dir,
                               "fredmd_per_factor_Q1minusQ5_SR.csv"))

# Save the strategy returns themselves so we can re-load and analyse later
write_csv(strat_fmd$diff_default,
          file.path(CFG$paths$output_dir, "fredmd_diff_returns.csv"))
write_csv(strat_fmd$per_bucket %>%
            pivot_wider(names_from = quintile, values_from = avg),
          file.path(CFG$paths$output_dir, "fredmd_quintile_returns.csv"))

# Persist the in-memory objects for downstream diagnostics
saveRDS(list(state_fmd = state_fmd,
             dist_mat_fmd = dist_mat_fmd,
             factors = factors,
             strat_fmd = strat_fmd,
             pca_out = pca_out,
             static_out = static_out,
             R2_table = R2_table,
             mR2 = mR2,
             CFG_FREDMD = CFG_FREDMD),
        file.path(CFG$paths$output_dir, "fredmd_results.rds"))

message("\nDone. FRED-MD output written to '", CFG$paths$output_dir, "'.")
message("To run the diagnostic suite on the FRED-MD strategies:")
message("  source('R/diagnostics.R')")
message("  run_diagnostics(state_fmd, dist_mat_fmd, factors, strat_fmd, CFG_FREDMD)")
