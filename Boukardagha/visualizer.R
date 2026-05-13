###############################################################################
# visualizer.R — Publication-Quality Charts (Base R Graphics)
# Boukardagha (2026) Replication
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


# ═══════════════════════════════════════════════════════════════════════════════
# 1. CUMULATIVE PnL — COLORED BY REGIME (Paper Fig. 1)
# ═══════════════════════════════════════════════════════════════════════════════

plot_cumulative_pnl_by_regime <- function(result) {
  fp <- save_png("01_cumulative_pnl_regime.png")
  
  dates <- result$dates
  cum_ret <- cumsum(result$returns)
  dom_reg <- result$dominant_regime
  n <- length(dates)
  
  par(mar = c(4, 4, 3, 1))
  plot(dates, cum_ret, type = "n",
       main = "Cumulative PnL Colored by Template Regime",
       xlab = "Date", ylab = "Cumulative Log Return",
       xaxt = "n")
  add_grid()
  date_axis(dates)
  
  # Plot segments colored by regime
  for (i in 2:n) {
    g <- dom_reg[i]
    col <- if (!is.na(g) && g >= 1 && g <= length(REGIME_COLORS)) {
      REGIME_COLORS[g]
    } else {
      "gray50"
    }
    segments(dates[i - 1], cum_ret[i - 1], dates[i], cum_ret[i],
             col = col, lwd = 1.5)
  }
  
  # Legend
  reg_used <- sort(unique(dom_reg[!is.na(dom_reg)]))
  reg_labels <- paste0("Regime ", LETTERS[reg_used])
  legend("topleft", legend = reg_labels,
         col = REGIME_COLORS[reg_used], lwd = 2,
         bg = "white", cex = 0.8)
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. CUMULATIVE PnL COMPARISON — HMM vs EW vs 60/40 vs SPX (Paper Fig. 11)
# ═══════════════════════════════════════════════════════════════════════════════

plot_cumulative_comparison <- function(hmm_result, sf_result, spx_result) {
  fp <- save_png("02_cumulative_comparison.png")
  
  dates <- hmm_result$dates
  cum_hmm <- cumsum(hmm_result$returns)
  cum_sf  <- cumsum(sf_result$returns)
  cum_spx <- cumsum(spx_result$returns)
  
  ylim <- range(c(cum_hmm, cum_sf, cum_spx), na.rm = TRUE)
  
  par(mar = c(4, 4, 3, 1))
  plot(dates, cum_hmm, type = "l", col = "#1f77b4", lwd = 2,
       ylim = ylim,
       main = "Cumulative PnL Comparison (OOS)",
       xlab = "Date", ylab = "Cumulative Log Return",
       xaxt = "n")
  add_grid()
  date_axis(dates)
  lines(dates, cum_sf,  col = "#2ca02c", lwd = 1.5, lty = 3)
  lines(dates, cum_spx, col = "#d62728", lwd = 1.5, lty = 4)
  
  legend("topleft",
         legend = c("Wasserstein HMM + MVO",
                    "60/40 Stock/Bond", "SPX Buy & Hold"),
         col = c("#1f77b4", "#2ca02c", "#d62728"),
         lwd = c(2, 1.5, 1.5), lty = c(1, 3, 4),
         bg = "white", cex = 0.8)
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. STACKED PORTFOLIO WEIGHTS (Paper Fig. 6)
# ═══════════════════════════════════════════════════════════════════════════════

plot_stacked_weights <- function(result, title = "Portfolio Weights (Stacked)") {
  fp <- save_png("03_stacked_weights.png")
  
  dates <- result$dates
  w <- result$weights
  n <- nrow(w)
  n_assets <- ncol(w)
  
  # Colors for assets
  asset_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  
  par(mar = c(4, 4, 3, 8), xpd = TRUE)
  
  # Compute cumulative areas
  cum_w <- t(apply(w, 1, cumsum))
  
  plot(dates, rep(0, n), type = "n", ylim = c(0, 1),
       main = title, xlab = "Date", ylab = "Weight",
       xaxt = "n")
  add_grid()
  date_axis(dates)
  
  for (j in n_assets:1) {
    upper <- cum_w[, j]
    lower <- if (j > 1) cum_w[, j - 1] else rep(0, n)
    polygon(c(dates, rev(dates)), c(upper, rev(lower)),
            col = asset_colors[j], border = NA)
  }
  
  legend("right", inset = c(-0.12, 0),
         legend = ASSET_LABELS, fill = asset_colors,
         bg = "white", cex = 0.75)
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 4. TURNOVER TIME SERIES (Paper Fig. 4)
# ═══════════════════════════════════════════════════════════════════════════════

plot_turnover <- function(result, title = "Daily Turnover — Wasserstein HMM") {
  fp <- save_png("04_turnover.png")
  
  dates <- result$dates
  to <- result$turnover
  
  par(mar = c(4, 4, 3, 1))
  plot(dates, to, type = "h", col = "#1f77b4",
       main = title, xlab = "Date", ylab = "Daily Turnover (0.5 * L1)",
       xaxt = "n")
  add_grid()
  date_axis(dates)
  
  # Average line
  avg_to <- mean(to, na.rm = TRUE)
  abline(h = avg_to, col = "red", lty = 2, lwd = 1.5)
  text(min(dates), avg_to, sprintf("Avg: %.4f", avg_to),
       pos = 4, col = "red", cex = 0.8)
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. REGIME HEATMAP / TIMELINE (Paper Fig. 10 variant)
# ═══════════════════════════════════════════════════════════════════════════════

plot_regime_heatmap <- function(result) {
  fp <- save_png("05_regime_heatmap.png", height = 4)
  
  dates <- result$dates
  tprob <- result$template_probs
  G <- ncol(tprob)
  n <- nrow(tprob)
  
  par(mar = c(4, 5, 3, 1))
  
  # Plot template probability time series (stacked area)
  cum_p <- t(apply(tprob, 1, function(x) {
    x[is.na(x)] <- 0
    cumsum(x)
  }))
  
  plot(dates, rep(0, n), type = "n", ylim = c(0, 1),
       main = "Template Regime Probabilities Over Time",
       xlab = "Date", ylab = "Probability",
       xaxt = "n")
  date_axis(dates)
  
  for (g in G:1) {
    upper <- cum_p[, g]
    lower <- if (g > 1) cum_p[, g - 1] else rep(0, n)
    polygon(c(dates, rev(dates)), c(upper, rev(lower)),
            col = adjustcolor(REGIME_COLORS[g], alpha.f = 0.7), border = NA)
  }
  
  reg_labels <- paste0("Tmpl ", seq_len(G))
  legend("topright", legend = reg_labels,
         fill = REGIME_COLORS[seq_len(G)], bg = "white", cex = 0.7)
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 6. ASSET SHARPE BY REGIME (Paper Fig. 9)
# ═══════════════════════════════════════════════════════════════════════════════

plot_asset_sharpe_by_regime <- function(reg_metrics) {
  fp <- save_png("06_asset_sharpe_by_regime.png")
  
  if (is.null(reg_metrics) || is.null(reg_metrics$asset_sharpe)) {
    plot.new(); text(0.5, 0.5, "No regime data"); dev.off(); return()
  }
  
  df <- reg_metrics$asset_sharpe
  regimes <- df$Regime
  n_reg <- nrow(df)
  n_assets <- ncol(df) - 1  # exclude Regime column
  
  sharpe_mat <- as.matrix(df[, -1])
  
  # Grouped bar chart
  asset_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  
  par(mar = c(5, 4, 3, 8), xpd = TRUE)
  bp <- barplot(t(sharpe_mat), beside = TRUE,
                names.arg = paste0("Regime ", LETTERS[regimes]),
                col = asset_colors, ylim = range(sharpe_mat) * 1.2,
                main = "Asset Sharpe Ratios by Regime",
                ylab = "Sharpe Ratio (Annualized)",
                border = NA)
  abline(h = 0, col = "black", lwd = 0.5)
  add_grid()
  
  legend("right", inset = c(-0.14, 0),
         legend = ASSET_LABELS, fill = asset_colors,
         bg = "white", cex = 0.75)
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 7. WASSERSTEIN DISTANCE OVER TIME
# ═══════════════════════════════════════════════════════════════════════════════

plot_wasserstein_distances <- function(result) {
  fp <- save_png("07_wasserstein_distances.png")
  
  dates <- result$dates
  w2d <- result$w2_distances
  valid <- !is.na(w2d)
  
  if (sum(valid) < 10) {
    plot.new(); text(0.5, 0.5, "Insufficient W2 data"); dev.off(); return()
  }
  
  par(mar = c(4, 4, 3, 1))
  plot(dates[valid], w2d[valid], type = "l", col = "#1f77b4",
       main = "Min Wasserstein Distance (Component → Template)",
       xlab = "Date", ylab = "W2 Distance",
       xaxt = "n")
  add_grid()
  date_axis(dates[valid])
  
  # Smoothed trend
  if (sum(valid) > 50) {
    w2_smooth <- filter(w2d[valid], rep(1 / 21, 21), sides = 1)
    lines(dates[valid], w2_smooth, col = "red", lwd = 2)
    legend("topright", legend = c("Raw", "21d MA"),
           col = c("#1f77b4", "red"), lwd = c(1, 2), bg = "white")
  }
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 8. EFFECTIVE DIVERSIFICATION (N_eff) OVER TIME (Paper Fig. 8)
# ═══════════════════════════════════════════════════════════════════════════════

plot_neff <- function(result) {
  fp <- save_png("08_neff.png", height = 5)
  
  dates <- result$dates
  w <- result$weights
  neff <- 1 / rowSums(w^2, na.rm = TRUE)
  
  par(mar = c(4, 4, 3, 1))
  plot(dates, neff, type = "l", col = "#1f77b4", lwd = 1,
       main = "Portfolio Concentration Over Time",
       xlab = "Date", ylab = "Effective # Positions (N_eff)",
       xaxt = "n", ylim = c(1, N_ASSETS + 0.5))
  add_grid()
  date_axis(dates)
  abline(h = mean(neff, na.rm = TRUE), col = "red", lty = 2)
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 9. DRAWDOWN CHART
# ═══════════════════════════════════════════════════════════════════════════════

plot_drawdowns <- function(hmm_result, sf_result, spx_result) {
  fp <- save_png("09_drawdowns.png")
  
  dd_fn <- function(r) {
    cum <- cumsum(r)
    cum - cummax(cum)
  }
  
  dates  <- hmm_result$dates
  dd_hmm <- dd_fn(hmm_result$returns)
  dd_sf  <- dd_fn(sf_result$returns)
  dd_spx <- dd_fn(spx_result$returns)
  
  ylim <- range(c(dd_hmm, dd_sf, dd_spx), na.rm = TRUE)
  
  par(mar = c(4, 4, 3, 1))
  plot(dates, dd_hmm, type = "l", col = "#1f77b4", lwd = 2,
       ylim = ylim,
       main = "Drawdown Comparison (OOS)",
       xlab = "Date", ylab = "Drawdown (log)",
       xaxt = "n")
  add_grid()
  date_axis(dates)
  lines(dates, dd_sf,  col = "#2ca02c", lwd = 1.5, lty = 3)
  lines(dates, dd_spx, col = "#d62728", lwd = 1.5, lty = 4)
  abline(h = 0, col = "black", lwd = 0.5)
  
  legend("bottomleft",
         legend = c("Wasserstein HMM", "60/40 Stock/Bond", "SPX Buy & Hold"),
         col = c("#1f77b4", "#2ca02c", "#d62728"),
         lwd = c(2, 1.5, 1.5), lty = c(1, 3, 4),
         bg = "white", cex = 0.8)
  
  dev.off()
  cat(sprintf("[plot] Saved: %s\n", fp))
}


# ═══════════════════════════════════════════════════════════════════════════════
# 10. SENSITIVITY ANALYSIS GRIDS
# ═══════════════════════════════════════════════════════════════════════════════

plot_sensitivity <- function(sens_df) {
  if (is.null(sens_df) || nrow(sens_df) == 0) return()
  
  params <- unique(sens_df$param)
  
  for (p in params) {
    fp <- save_png(sprintf("10_sensitivity_%s.png", p), height = 5)
    
    sub <- sens_df[sens_df$param == p, ]
    sub <- sub[!is.na(sub$Sharpe), ]
    if (nrow(sub) == 0) { dev.off(); next }
    
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    
    # Sharpe
    plot(sub$value, sub$Sharpe, type = "b", pch = 19, col = "#1f77b4",
         main = paste0(p, " → Sharpe"), xlab = p, ylab = "Sharpe")
    add_grid()
    abline(v = sub$value[which.max(sub$Sharpe)], col = "red", lty = 2)
    
    # Max DD
    plot(sub$value, sub$Max_DD, type = "b", pch = 19, col = "#d62728",
         main = paste0(p, " → Max DD"), xlab = p, ylab = "Max Drawdown")
    add_grid()
    
    # Turnover
    plot(sub$value, sub$Avg_Turnover, type = "b", pch = 19, col = "#2ca02c",
         main = paste0(p, " → Turnover"), xlab = p, ylab = "Avg Daily Turnover")
    add_grid()
    
    par(mfrow = c(1, 1))
    dev.off()
    cat(sprintf("[plot] Saved: sensitivity_%s.png\n", p))
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# 11. MODEL ORDER (K) PATH OVER TIME
# ═══════════════════════════════════════════════════════════════════════════════

plot_model_order <- function(result) {
  fp <- save_png("11_model_order_K.png", height = 4)
  
  dates <- result$dates
  k_path <- result$current_K
  valid <- !is.na(k_path)
  
  par(mar = c(4, 4, 3, 1))
  plot(dates[valid], k_path[valid], type = "s", col = "#1f77b4", lwd = 1.5,
       main = "Predictive Model Order (K) Over Time",
       xlab = "Date", ylab = "Number of HMM States",
       ylim = c(K_MIN - 0.5, K_MAX + 0.5),
       xaxt = "n")
  add_grid()
  date_axis(dates[valid])
  
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
  
  plot_cumulative_pnl_by_regime(hmm_result)
  plot_cumulative_comparison(hmm_result, sf_result, spx_result)
  plot_stacked_weights(hmm_result)
  plot_turnover(hmm_result)
  plot_regime_heatmap(hmm_result)
  
  if (!is.null(reg_metrics)) {
    plot_asset_sharpe_by_regime(reg_metrics)
  }
  
  plot_wasserstein_distances(hmm_result)
  plot_neff(hmm_result)
  plot_drawdowns(hmm_result, sf_result, spx_result)
  plot_model_order(hmm_result)
  
  if (!is.null(sens_df)) {
    plot_sensitivity(sens_df)
  }
  
  cat(sprintf("\n[visualizer] All plots saved to %s/\n", OUTPUT_DIR))
}


cat("[visualizer.R] Functions loaded.\n")
