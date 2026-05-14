###############################################################################
# visualizer.R — Publication-Quality Charts (Base R Graphics)
# Boukardagha (2026) Replication
#
# Original: 10 charts (cumulative PnL, comparison, weights, turnover,
#           regime heatmap, asset Sharpe, W2 distances, N_eff, drawdowns, K path)
# New:      6 charts (regime return densities, rolling Sharpe, regime transitions,
#           regime durations, per-asset weight panels, correlation by regime)
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

save_png <- function(filename, width = PLOT_WIDTH, height = PLOT_HEIGHT,
                     dpi = PLOT_DPI) {
  filepath <- file.path(OUTPUT_DIR, filename)
  png(filepath, width = width, height = height, units = "in", res = dpi)
  filepath
}

date_axis <- function(dates, ...) {
  at_dates <- pretty(dates, n = 8)
  axis(1, at = at_dates, labels = format(at_dates, "%Y-%m"), ...)
}

add_grid <- function() {
  grid(col = "gray90", lty = 1)
}

# Polished par defaults for multi-panel layouts
par_clean <- function(mar = c(4, 4, 3, 1)) {
  par(mar = mar, family = "", cex.main = 1.1, cex.lab = 0.95,
      cex.axis = 0.85, las = 1)
}

# Compute kernel density safely (handles edge cases)
safe_density <- function(x, bw = "SJ", n = 256, ...) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(NULL)
  tryCatch(density(x, bw = bw, n = n, ...), error = function(e) {
    tryCatch(density(x, bw = "nrd0", n = n, ...), error = function(e2) NULL)
  })
}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. CUMULATIVE PnL — COLORED BY REGIME
# ═══════════════════════════════════════════════════════════════════════════════

plot_cumulative_pnl_by_regime <- function(result) {
  fp <- save_png("01_cumulative_pnl_regime.png")
  par_clean()

  dates <- result$dates
  cum_ret <- cumsum(result$returns)
  dom_reg <- result$dominant_regime
  n <- length(dates)

  plot(dates, cum_ret, type = "n",
       main = "Cumulative PnL Colored by Template Regime",
       xlab = "Date", ylab = "Cumulative Log Return", xaxt = "n")
  add_grid()
  date_axis(dates)

  for (i in 2:n) {
    g <- dom_reg[i]
    col <- if (!is.na(g) && g >= 1 && g <= length(REGIME_COLORS)) {
      REGIME_COLORS[g]
    } else { "gray50" }
    segments(dates[i - 1], cum_ret[i - 1], dates[i], cum_ret[i],
             col = col, lwd = 1.5)
  }

  reg_used <- sort(unique(dom_reg[!is.na(dom_reg)]))
  legend("topleft", legend = paste0("Regime ", LETTERS[reg_used]),
         col = REGIME_COLORS[reg_used], lwd = 2, bg = "white", cex = 0.8)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. CUMULATIVE PnL COMPARISON — HMM vs 60/40 vs SPX
# ═══════════════════════════════════════════════════════════════════════════════

plot_cumulative_comparison <- function(hmm_result, sf_result, spx_result) {
  fp <- save_png("02_cumulative_comparison.png")
  par_clean()

  dates <- hmm_result$dates
  cum_hmm <- cumsum(hmm_result$returns)
  cum_sf  <- cumsum(sf_result$returns)
  cum_spx <- cumsum(spx_result$returns)

  ylim <- range(c(cum_hmm, cum_sf, cum_spx), na.rm = TRUE)

  plot(dates, cum_hmm, type = "l", col = REGIME_COLORS[1], lwd = 2.5,
       ylim = ylim, main = "Cumulative PnL Comparison (OOS)",
       xlab = "Date", ylab = "Cumulative Log Return", xaxt = "n")
  add_grid()
  date_axis(dates)
  lines(dates, cum_sf,  col = REGIME_COLORS[3], lwd = 1.8, lty = 5)
  lines(dates, cum_spx, col = REGIME_COLORS[4], lwd = 1.8, lty = 4)

  legend("topleft",
         legend = c("Wasserstein HMM + MVO", "60/40 Stock/Bond", "SPX Buy & Hold"),
         col = c(REGIME_COLORS[1], REGIME_COLORS[3], REGIME_COLORS[4]),
         lwd = c(2.5, 1.8, 1.8), lty = c(1, 5, 4), bg = "white", cex = 0.85)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. STACKED PORTFOLIO WEIGHTS
# ═══════════════════════════════════════════════════════════════════════════════

plot_stacked_weights <- function(result, title = "Portfolio Weights (Stacked)") {
  fp <- save_png("03_stacked_weights.png")
  par(mar = c(4, 4, 3, 8), xpd = TRUE, cex.main = 1.1, las = 1)

  dates <- result$dates; w <- result$weights; n <- nrow(w)
  cum_w <- t(apply(w, 1, cumsum))

  plot(dates, rep(0, n), type = "n", ylim = c(0, 1),
       main = title, xlab = "Date", ylab = "Weight", xaxt = "n")
  add_grid()
  date_axis(dates)

  for (j in N_ASSETS:1) {
    upper <- cum_w[, j]
    lower <- if (j > 1) cum_w[, j - 1] else rep(0, n)
    polygon(c(dates, rev(dates)), c(upper, rev(lower)),
            col = adjustcolor(ASSET_COLORS[j], alpha.f = 0.8), border = NA)
  }

  legend("right", inset = c(-0.12, 0), legend = ASSET_LABELS,
         fill = ASSET_COLORS, bg = "white", cex = 0.75)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. TURNOVER TIME SERIES
# ═══════════════════════════════════════════════════════════════════════════════

plot_turnover <- function(result) {
  fp <- save_png("04_turnover.png")
  par_clean()

  dates <- result$dates; to <- result$turnover

  plot(dates, to, type = "h", col = adjustcolor(REGIME_COLORS[1], 0.5),
       main = "Daily Turnover — Wasserstein HMM", xlab = "Date",
       ylab = "Daily Turnover (0.5 × L1)", xaxt = "n")
  add_grid()
  date_axis(dates)

  avg_to <- mean(to, na.rm = TRUE)
  abline(h = avg_to, col = REGIME_COLORS[4], lty = 2, lwd = 1.5)
  text(min(dates), avg_to * 1.2, sprintf("Avg: %.4f", avg_to),
       pos = 4, col = REGIME_COLORS[4], cex = 0.85)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. REGIME PROBABILITY TIMELINE (stacked area)
# ═══════════════════════════════════════════════════════════════════════════════

plot_regime_heatmap <- function(result) {
  fp <- save_png("05_regime_heatmap.png", height = 4)
  par(mar = c(4, 5, 3, 1), cex.main = 1.1, las = 1)

  dates <- result$dates; tprob <- result$template_probs
  G <- ncol(tprob); n <- nrow(tprob)

  cum_p <- t(apply(tprob, 1, function(x) { x[is.na(x)] <- 0; cumsum(x) }))

  plot(dates, rep(0, n), type = "n", ylim = c(0, 1),
       main = "Template Regime Probabilities Over Time",
       xlab = "Date", ylab = "Probability", xaxt = "n")
  date_axis(dates)

  for (g in G:1) {
    upper <- cum_p[, g]
    lower <- if (g > 1) cum_p[, g - 1] else rep(0, n)
    polygon(c(dates, rev(dates)), c(upper, rev(lower)),
            col = adjustcolor(REGIME_COLORS[g], alpha.f = 0.7), border = NA)
  }

  legend("topright", legend = paste0("Tmpl ", seq_len(G)),
         fill = REGIME_COLORS[seq_len(G)], bg = "white", cex = 0.7)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 6. ASSET SHARPE BY REGIME (grouped bar chart)
# ═══════════════════════════════════════════════════════════════════════════════

plot_asset_sharpe_by_regime <- function(reg_metrics) {
  fp <- save_png("06_asset_sharpe_by_regime.png")

  if (is.null(reg_metrics) || is.null(reg_metrics$asset_sharpe)) {
    plot.new(); text(0.5, 0.5, "No regime data"); dev.off(); return()
  }

  df <- reg_metrics$asset_sharpe
  regimes <- df$Regime; sharpe_mat <- as.matrix(df[, -1])

  par(mar = c(5, 4, 3, 8), xpd = TRUE, cex.main = 1.1, las = 1)
  bp <- barplot(t(sharpe_mat), beside = TRUE,
                names.arg = paste0("Regime ", LETTERS[regimes]),
                col = ASSET_COLORS, ylim = range(sharpe_mat) * 1.2,
                main = "Asset Sharpe Ratios by Regime",
                ylab = "Sharpe Ratio (Annualized)", border = NA)
  abline(h = 0, col = "black", lwd = 0.5)
  add_grid()

  legend("right", inset = c(-0.14, 0), legend = ASSET_LABELS,
         fill = ASSET_COLORS, bg = "white", cex = 0.75)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 7. WASSERSTEIN DISTANCE OVER TIME
# ═══════════════════════════════════════════════════════════════════════════════

plot_wasserstein_distances <- function(result) {
  fp <- save_png("07_wasserstein_distances.png")
  par_clean()

  dates <- result$dates; w2d <- result$w2_distances
  valid <- !is.na(w2d)
  if (sum(valid) < 10) { plot.new(); text(0.5, 0.5, "Insufficient W2 data"); dev.off(); return() }

  plot(dates[valid], w2d[valid], type = "l",
       col = adjustcolor(REGIME_COLORS[1], 0.5), lwd = 0.8,
       main = "Min Wasserstein Distance (Component → Template)",
       xlab = "Date", ylab = "W2 Distance", xaxt = "n")
  add_grid()
  date_axis(dates[valid])

  if (sum(valid) > 50) {
    w2_smooth <- filter(w2d[valid], rep(1 / 21, 21), sides = 1)
    lines(dates[valid], w2_smooth, col = REGIME_COLORS[4], lwd = 2)
    legend("topright", legend = c("Raw", "21d MA"),
           col = c(REGIME_COLORS[1], REGIME_COLORS[4]), lwd = c(1, 2), bg = "white")
  }

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 8. EFFECTIVE DIVERSIFICATION (N_eff)
# ═══════════════════════════════════════════════════════════════════════════════

plot_neff <- function(result) {
  fp <- save_png("08_neff.png", height = 5)
  par_clean()

  dates <- result$dates; w <- result$weights
  neff <- 1 / rowSums(w^2, na.rm = TRUE)

  plot(dates, neff, type = "l", col = REGIME_COLORS[1], lwd = 1,
       main = "Portfolio Concentration Over Time",
       xlab = "Date", ylab = expression("Effective # Positions ("*N[eff]*")"),
       xaxt = "n", ylim = c(1, N_ASSETS + 0.5))
  add_grid()
  date_axis(dates)
  abline(h = mean(neff, na.rm = TRUE), col = REGIME_COLORS[4], lty = 2)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 9. DRAWDOWN CHART — HMM vs 60/40 vs SPX
# ═══════════════════════════════════════════════════════════════════════════════

plot_drawdowns <- function(hmm_result, sf_result, spx_result) {
  fp <- save_png("09_drawdowns.png")
  par_clean()

  dd_fn <- function(r) { cum <- cumsum(r); cum - cummax(cum) }

  dates  <- hmm_result$dates
  dd_hmm <- dd_fn(hmm_result$returns)
  dd_sf  <- dd_fn(sf_result$returns)
  dd_spx <- dd_fn(spx_result$returns)

  ylim <- range(c(dd_hmm, dd_sf, dd_spx), na.rm = TRUE)

  plot(dates, dd_hmm, type = "l", col = REGIME_COLORS[1], lwd = 2, ylim = ylim,
       main = "Drawdown Comparison (OOS)",
       xlab = "Date", ylab = "Drawdown (log)", xaxt = "n")
  add_grid()
  date_axis(dates)
  lines(dates, dd_sf,  col = REGIME_COLORS[3], lwd = 1.5, lty = 5)
  lines(dates, dd_spx, col = REGIME_COLORS[4], lwd = 1.5, lty = 4)
  abline(h = 0, col = "black", lwd = 0.5)

  legend("bottomleft",
         legend = c("Wasserstein HMM", "60/40 Stock/Bond", "SPX Buy & Hold"),
         col = c(REGIME_COLORS[1], REGIME_COLORS[3], REGIME_COLORS[4]),
         lwd = c(2, 1.5, 1.5), lty = c(1, 5, 4), bg = "white", cex = 0.8)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 10. SENSITIVITY GRIDS
# ═══════════════════════════════════════════════════════════════════════════════

plot_sensitivity <- function(sens_df) {
  if (is.null(sens_df) || nrow(sens_df) == 0) return()

  params <- unique(sens_df$param)

  for (p in params) {
    fp <- save_png(sprintf("10_sensitivity_%s.png", p), height = 5)
    sub <- sens_df[sens_df$param == p, ]
    sub <- sub[!is.na(sub$Sharpe), ]
    if (nrow(sub) == 0) { dev.off(); next }

    par(mfrow = c(1, 3), mar = c(4, 4, 2, 1), cex.main = 1.0)

    plot(sub$value, sub$Sharpe, type = "b", pch = 19, col = REGIME_COLORS[1],
         main = paste0(p, " → Sharpe"), xlab = p, ylab = "Sharpe")
    add_grid(); abline(v = sub$value[which.max(sub$Sharpe)], col = REGIME_COLORS[4], lty = 2)

    plot(sub$value, sub$Max_DD, type = "b", pch = 19, col = REGIME_COLORS[4],
         main = paste0(p, " → Max DD"), xlab = p, ylab = "Max Drawdown")
    add_grid()

    plot(sub$value, sub$Avg_Turnover, type = "b", pch = 19, col = REGIME_COLORS[3],
         main = paste0(p, " → Turnover"), xlab = p, ylab = "Avg Daily Turnover")
    add_grid()

    par(mfrow = c(1, 1))
    dev.off()
    cat(sprintf("[plot] Saved: sensitivity_%s.png\n", p))
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# 11. MODEL ORDER (K) PATH
# ═══════════════════════════════════════════════════════════════════════════════

plot_model_order <- function(result) {
  fp <- save_png("11_model_order_K.png", height = 4)
  par_clean()

  dates <- result$dates; k_path <- result$current_K
  valid <- !is.na(k_path)

  plot(dates[valid], k_path[valid], type = "s", col = REGIME_COLORS[1], lwd = 1.5,
       main = "Predictive Model Order (K) Over Time",
       xlab = "Date", ylab = "Number of HMM States",
       ylim = c(K_MIN - 0.5, K_MAX + 0.5), xaxt = "n")
  add_grid()
  date_axis(dates[valid])

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NEW CHART 12: REGIME RETURN DENSITY PLOTS
#  Overlaid KDE per regime for each asset — inspired by paper Fig. 7 style
# ═══════════════════════════════════════════════════════════════════════════════

plot_regime_return_densities <- function(result) {
  fp <- save_png("12_regime_return_densities.png", width = 14, height = 10)

  r <- result$asset_returns
  reg <- result$dominant_regime
  valid <- !is.na(reg)
  r <- r[valid, ]; reg <- reg[valid]
  regimes <- sort(unique(reg))

  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), cex.main = 1.1, las = 1)

  for (j in seq_len(N_ASSETS)) {
    # Compute densities per regime
    dens_list <- list()
    max_y <- 0
    x_range <- range(r[, j], na.rm = TRUE)

    for (g in regimes) {
      d <- safe_density(r[reg == g, j])
      if (!is.null(d)) {
        dens_list[[as.character(g)]] <- d
        max_y <- max(max_y, max(d$y))
      }
    }
    if (length(dens_list) == 0) { plot.new(); next }

    plot(NULL, xlim = x_range, ylim = c(0, max_y * 1.1),
         main = paste0(ASSET_LABELS[j], " — Daily Return Distributions"),
         xlab = "Daily Log Return", ylab = "Density")
    add_grid()

    for (gi in seq_along(dens_list)) {
      g_name <- names(dens_list)[gi]
      d <- dens_list[[g_name]]
      g_idx <- as.integer(g_name)
      col_fill <- adjustcolor(REGIME_COLORS[g_idx], alpha.f = 0.20)
      col_line <- REGIME_COLORS[g_idx]
      polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))),
              col = col_fill, border = NA)
      lines(d$x, d$y, col = col_line, lwd = 2)
    }
    abline(v = 0, col = "gray40", lty = 3)
  }

  # Legend panel in the 6th slot
  plot.new()
  legend("center", legend = paste0("Regime ", LETTERS[regimes]),
         fill = adjustcolor(REGIME_COLORS[regimes], alpha.f = 0.4),
         border = REGIME_COLORS[regimes],
         title = "Template Regimes", cex = 1.2, bty = "n")

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NEW CHART 13: PORTFOLIO RETURN DENSITY BY REGIME
#  Single panel — regime distributions of the W-HMM portfolio itself
# ═══════════════════════════════════════════════════════════════════════════════

plot_portfolio_return_density <- function(result) {
  fp <- save_png("13_portfolio_return_by_regime.png", width = 10, height = 6)
  par_clean(mar = c(4, 4, 3, 2))

  r <- result$returns
  reg <- result$dominant_regime
  valid <- !is.na(reg) & !is.na(r)
  r <- r[valid]; reg <- reg[valid]
  regimes <- sort(unique(reg))

  dens_list <- list(); max_y <- 0
  for (g in regimes) {
    d <- safe_density(r[reg == g])
    if (!is.null(d)) { dens_list[[as.character(g)]] <- d; max_y <- max(max_y, max(d$y)) }
  }
  d_all <- safe_density(r)

  x_range <- range(r)
  plot(NULL, xlim = x_range, ylim = c(0, max(max_y, max(d_all$y)) * 1.1),
       main = "W-HMM Portfolio Return Density by Regime",
       xlab = "Daily Log Return", ylab = "Density")
  add_grid()

  for (gi in seq_along(dens_list)) {
    g_name <- names(dens_list)[gi]; g_idx <- as.integer(g_name)
    d <- dens_list[[g_name]]
    polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))),
            col = adjustcolor(REGIME_COLORS[g_idx], alpha.f = 0.18), border = NA)
    lines(d$x, d$y, col = REGIME_COLORS[g_idx], lwd = 2)
  }
  # Overall in gray
  lines(d_all$x, d_all$y, col = "gray30", lwd = 2, lty = 2)
  abline(v = 0, col = "gray40", lty = 3)

  legend("topright",
         legend = c(paste0("Regime ", LETTERS[regimes]), "Overall"),
         col = c(REGIME_COLORS[regimes], "gray30"),
         lwd = 2, lty = c(rep(1, length(regimes)), 2),
         bg = "white", cex = 0.8)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NEW CHART 14: ROLLING SHARPE RATIO COMPARISON
# ═══════════════════════════════════════════════════════════════════════════════

plot_rolling_sharpe <- function(hmm_result, sf_result, spx_result, window = 252) {
  fp <- save_png("14_rolling_sharpe.png")
  par_clean()

  rolling_sharpe <- function(r, w) {
    n <- length(r); out <- rep(NA_real_, n)
    if (n < w) return(out)
    cs  <- cumsum(r); cs2 <- cumsum(r^2)
    for (i in w:n) {
      s  <- cs[i] - ifelse(i - w > 0, cs[i - w], 0)
      s2 <- cs2[i] - ifelse(i - w > 0, cs2[i - w], 0)
      mu <- s / w; var_est <- s2 / w - mu^2
      if (var_est > 0) out[i] <- (mu * 252) / (sqrt(var_est * 252 / (w - 1)))
    }
    out
  }

  dates <- hmm_result$dates
  rs_hmm <- rolling_sharpe(hmm_result$returns, window)
  rs_sf  <- rolling_sharpe(sf_result$returns, window)
  rs_spx <- rolling_sharpe(spx_result$returns, window)

  ylim <- range(c(rs_hmm, rs_sf, rs_spx), na.rm = TRUE)
  ylim <- c(max(ylim[1], -4), min(ylim[2], 6))

  plot(dates, rs_hmm, type = "l", col = REGIME_COLORS[1], lwd = 2,
       ylim = ylim, main = sprintf("Rolling %d-Day Sharpe Ratio", window),
       xlab = "Date", ylab = "Sharpe Ratio (Annualized)", xaxt = "n")
  add_grid()
  date_axis(dates)
  lines(dates, rs_sf,  col = REGIME_COLORS[3], lwd = 1.5, lty = 5)
  lines(dates, rs_spx, col = REGIME_COLORS[4], lwd = 1.5, lty = 4)
  abline(h = 0, col = "black", lwd = 0.5)

  legend("topright",
         legend = c("W-HMM", "60/40", "SPX"),
         col = c(REGIME_COLORS[1], REGIME_COLORS[3], REGIME_COLORS[4]),
         lwd = c(2, 1.5, 1.5), lty = c(1, 5, 4), bg = "white", cex = 0.85)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NEW CHART 15: REGIME TRANSITION HEATMAP
#  Empirical transition frequencies between dominant regimes
# ═══════════════════════════════════════════════════════════════════════════════

plot_regime_transitions <- function(result) {
  fp <- save_png("15_regime_transitions.png", width = 8, height = 7)

  reg <- result$dominant_regime
  valid <- !is.na(reg)
  reg <- reg[valid]
  regimes <- sort(unique(reg))
  n_reg <- length(regimes)

  # Build empirical transition count matrix
  trans_mat <- matrix(0, n_reg, n_reg)
  rownames(trans_mat) <- colnames(trans_mat) <- paste0("", LETTERS[regimes])
  for (i in 2:length(reg)) {
    from <- which(regimes == reg[i - 1])
    to   <- which(regimes == reg[i])
    trans_mat[from, to] <- trans_mat[from, to] + 1
  }

  # Normalize to probabilities
  row_sums <- rowSums(trans_mat)
  row_sums[row_sums == 0] <- 1
  trans_prob <- trans_mat / row_sums

  par(mar = c(5, 5, 4, 6), xpd = TRUE)

  # Draw heatmap with rect()
  n <- n_reg
  color_ramp <- colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100)

  plot(NULL, xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5),
       xlab = "To Regime", ylab = "From Regime",
       main = "Empirical Regime Transition Probabilities",
       xaxt = "n", yaxt = "n", asp = 1)

  axis(1, at = 1:n, labels = paste0("", LETTERS[regimes]), tick = FALSE)
  axis(2, at = 1:n, labels = paste0("", LETTERS[regimes]), tick = FALSE)

  for (i in 1:n) {
    for (j in 1:n) {
      p <- trans_prob[i, j]
      col_idx <- max(1, min(100, ceiling(p * 100)))
      rect(j - 0.45, (n - i + 1) - 0.45, j + 0.45, (n - i + 1) + 0.45,
           col = color_ramp[col_idx], border = "white", lwd = 0.5)
      # Label
      txt_col <- if (p > 0.5) "white" else "gray20"
      text(j, n - i + 1, sprintf("%.0f%%", p * 100), col = txt_col, cex = 0.95)
    }
  }

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NEW CHART 16: REGIME DURATION DISTRIBUTION
#  How long does each regime persist before switching?
# ═══════════════════════════════════════════════════════════════════════════════

plot_regime_durations <- function(result) {
  fp <- save_png("16_regime_durations.png", width = 12, height = 6)

  reg <- result$dominant_regime
  valid <- !is.na(reg)
  reg <- reg[valid]
  regimes <- sort(unique(reg))

  # Extract spell durations
  rle_result <- rle(reg)
  spell_regime <- rle_result$values
  spell_len    <- rle_result$lengths

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), cex.main = 1.05, las = 1)

  # Panel 1: boxplot of duration by regime
  dur_by_regime <- split(spell_len, spell_regime)
  bp_data <- dur_by_regime[as.character(regimes)]

  boxplot(bp_data, names = paste0("", LETTERS[regimes]),
          col = adjustcolor(REGIME_COLORS[regimes], 0.5),
          border = REGIME_COLORS[regimes], outline = TRUE, pch = 16, cex = 0.4,
          main = "Regime Spell Duration (days)",
          xlab = "Regime", ylab = "Duration (trading days)")
  add_grid()

  # Panel 2: histogram of all durations (stacked by regime)
  max_dur <- quantile(spell_len, 0.98)
  breaks <- seq(0, max(spell_len) + 5, by = 5)
  plot(NULL, xlim = c(0, max_dur), ylim = c(0, 1),
       main = "Duration Distribution (all regimes)",
       xlab = "Duration (trading days)", ylab = "Density")
  add_grid()

  for (g in regimes) {
    d <- safe_density(spell_len[spell_regime == g])
    if (!is.null(d)) {
      max_d <- max(d$y)
      d$y <- d$y / max(max_d, 1e-6)  # normalize for overlay
      polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))),
              col = adjustcolor(REGIME_COLORS[g], 0.25), border = NA)
      lines(d$x, d$y, col = REGIME_COLORS[g], lwd = 2)
    }
  }

  legend("topright", legend = paste0("Regime ", LETTERS[regimes]),
         col = REGIME_COLORS[regimes], lwd = 2, bg = "white", cex = 0.75)

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
#  NEW CHART 17: PER-ASSET WEIGHT EVOLUTION (individual panels)
# ═══════════════════════════════════════════════════════════════════════════════

plot_individual_weights <- function(result) {
  fp <- save_png("17_individual_weights.png", width = 14, height = 10)

  dates <- result$dates; w <- result$weights
  reg <- result$dominant_regime

  par(mfrow = c(3, 2), mar = c(3, 4, 2.5, 1), cex.main = 1.0, las = 1)

  for (j in seq_len(N_ASSETS)) {
    plot(dates, w[, j], type = "n", ylim = c(0, max(w[, j], na.rm = TRUE) * 1.1),
         main = paste0(ASSET_LABELS[j], " Weight"),
         xlab = "", ylab = "Weight", xaxt = "n")
    add_grid()
    date_axis(dates)

    # Draw weight colored by regime
    n <- length(dates)
    for (i in 2:n) {
      g <- reg[i]
      col <- if (!is.na(g)) adjustcolor(REGIME_COLORS[g], 0.6) else "gray70"
      segments(dates[i - 1], w[i - 1, j], dates[i], w[i, j], col = col, lwd = 1.2)
    }

    # Rolling mean overlay
    if (n > 63) {
      w_smooth <- filter(w[, j], rep(1 / 63, 63), sides = 1)
      lines(dates, w_smooth, col = "black", lwd = 1.8, lty = 1)
    }
  }

  # Legend in 6th panel
  plot.new()
  reg_used <- sort(unique(reg[!is.na(reg)]))
  legend("center",
         legend = c(paste0("Regime ", LETTERS[reg_used]), "63d MA"),
         col = c(REGIME_COLORS[reg_used], "black"),
         lwd = c(rep(2, length(reg_used)), 2),
         lty = c(rep(1, length(reg_used)), 1),
         title = "Coloring", cex = 1.1, bty = "n")

  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# MASTER PLOT FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

generate_all_plots <- function(hmm_result, sf_result, spx_result,
                                reg_metrics = NULL, sens_df = NULL) {
  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("  GENERATING PLOTS\n")
  cat("══════════════════════════════════════════════════════════════\n\n")

  # ── Original charts ──
  plot_cumulative_pnl_by_regime(hmm_result)
  plot_cumulative_comparison(hmm_result, sf_result, spx_result)
  plot_stacked_weights(hmm_result)
  plot_turnover(hmm_result)
  plot_regime_heatmap(hmm_result)
  if (!is.null(reg_metrics)) plot_asset_sharpe_by_regime(reg_metrics)
  plot_wasserstein_distances(hmm_result)
  plot_neff(hmm_result)
  plot_drawdowns(hmm_result, sf_result, spx_result)
  plot_model_order(hmm_result)

  # ── New charts ──
  plot_regime_return_densities(hmm_result)
  plot_portfolio_return_density(hmm_result)
  plot_rolling_sharpe(hmm_result, sf_result, spx_result)
  plot_regime_transitions(hmm_result)
  plot_regime_durations(hmm_result)
  plot_individual_weights(hmm_result)

  if (!is.null(sens_df)) plot_sensitivity(sens_df)

  cat(sprintf("\n[visualizer] All plots saved to %s/\n", OUTPUT_DIR))
}


cat("[visualizer.R] Functions loaded.\n")
