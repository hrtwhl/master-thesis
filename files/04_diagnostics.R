# =============================================================================
#  04_diagnostics.R
#  Performance metrics, diagnostic tables, plots, CSV/PNG exports
# =============================================================================

TRADING_DAYS_PER_YEAR <- 252

# ---------------- Metrics ----------------

#' Sharpe, Sortino, MaxDD, ann mean, ann vol from a daily log-return series.
performance_metrics <- function(pnl) {
  pnl <- as.numeric(pnl)
  pnl <- pnl[is.finite(pnl)]
  if (length(pnl) == 0) {
    return(data.frame(ann_mean = NA, ann_vol = NA, sharpe = NA,
                      sortino = NA, max_dd = NA))
  }
  ann_mean <- mean(pnl) * TRADING_DAYS_PER_YEAR
  ann_vol  <- sd(pnl)   * sqrt(TRADING_DAYS_PER_YEAR)
  sharpe   <- if (ann_vol > 0) ann_mean / ann_vol else NA_real_
  
  downside <- pmin(pnl, 0)
  d_vol    <- sqrt(mean(downside^2)) * sqrt(TRADING_DAYS_PER_YEAR)
  sortino  <- if (d_vol > 0) ann_mean / d_vol else NA_real_
  
  cum <- cumsum(pnl)
  max_dd <- min(cum - cummax(cum))
  data.frame(ann_mean = ann_mean, ann_vol = ann_vol,
             sharpe = sharpe, sortino = sortino, max_dd = max_dd)
}

#' Compare performance across an arbitrary set of strategies.
compare_performance <- function(pnl_list) {
  rows <- lapply(names(pnl_list), function(nm) {
    cbind(strategy = nm, performance_metrics(pnl_list[[nm]]))
  })
  do.call(rbind, rows)
}

#' Turnover stats. Expects a data.frame with date + asset weight columns.
turnover_stats <- function(weights_df, asset_names) {
  W  <- as.matrix(weights_df[, asset_names])
  dW <- apply(W, 2, function(v) c(NA, diff(v)))
  TO <- 0.5 * rowSums(abs(dW), na.rm = FALSE)
  TO <- TO[is.finite(TO)]
  data.frame(
    avg_turnover  = mean(TO),
    median_turnover = median(TO),
    q95_turnover  = unname(quantile(TO, 0.95)),
    frac_above_1pct = mean(TO > 0.01),
    frac_above_5pct = mean(TO > 0.05)
  )
}

#' Average allocation summary (avg weight, weight vol, time>10%, avg |Δw|).
allocation_summary <- function(weights_df, asset_names) {
  W  <- as.matrix(weights_df[, asset_names])
  dW <- apply(W, 2, function(v) c(NA, diff(v)))
  data.frame(
    asset      = asset_names,
    avg_weight = colMeans(W),
    weight_vol = apply(W, 2, sd),
    time_gt_10 = colMeans(W > 0.10),
    avg_abs_dw = colMeans(abs(dW), na.rm = TRUE)
  )
}

#' Concentration metric: N_eff = 1 / sum(w_i^2). Average and median.
concentration_stats <- function(weights_df, asset_names) {
  W <- as.matrix(weights_df[, asset_names])
  n_eff <- 1 / pmax(rowSums(W^2), 1e-12)
  data.frame(avg_neff = mean(n_eff), median_neff = median(n_eff))
}

#' Per-regime portfolio performance.
regime_performance <- function(pnl_df, label_df) {
  m <- merge(pnl_df, label_df[, c("date","label")], by = "date")
  m <- m[!is.na(m$label), ]
  do.call(rbind, lapply(sort(unique(m$label)), function(g) {
    rg <- m$pnl[m$label == g]
    if (length(rg) < 2) return(NULL)
    cum <- cumsum(rg); dd <- min(cum - cummax(cum))
    data.frame(
      regime   = g,
      days     = length(rg),
      ann_mean = mean(rg) * TRADING_DAYS_PER_YEAR,
      ann_vol  = sd(rg)   * sqrt(TRADING_DAYS_PER_YEAR),
      sharpe   = if (sd(rg) > 0) mean(rg) * TRADING_DAYS_PER_YEAR /
                                 (sd(rg) * sqrt(TRADING_DAYS_PER_YEAR)) else NA_real_,
      hit_rate = mean(rg > 0),
      max_dd_within = dd
    )
  }))
}

#' Per-regime, per-asset Sharpe (annualized).
asset_sharpe_by_regime <- function(returns_df, label_df, asset_names) {
  m <- merge(returns_df[, c("date", asset_names)],
             label_df[, c("date","label")], by = "date")
  m <- m[!is.na(m$label), ]
  do.call(rbind, lapply(sort(unique(m$label)), function(g) {
    sub <- m[m$label == g, asset_names, drop = FALSE]
    sr  <- sapply(sub, function(x) {
      if (length(x) < 2 || sd(x) == 0) return(NA_real_)
      mean(x) * TRADING_DAYS_PER_YEAR / (sd(x) * sqrt(TRADING_DAYS_PER_YEAR))
    })
    data.frame(regime = g, t(sr))
  }))
}

# ---------------- CSV exports ----------------

export_tables <- function(results, benchmarks, asset_names, out_dir) {
  pnl_strat <- results$pnl$pnl
  pnl_ew    <- benchmarks$ew$pnl
  pnl_6040  <- benchmarks$sixty_forty$pnl
  pnl_spx   <- benchmarks$spx$pnl
  
  perf <- compare_performance(list(
    "Wasserstein HMM + MVO" = pnl_strat,
    "Equal Weight (20%)"    = pnl_ew,
    "60/40 Stocks/Bonds"    = pnl_6040,
    "SPX Buy & Hold"        = pnl_spx
  ))
  write.csv(perf, file.path(out_dir, "performance_metrics.csv"), row.names = FALSE)
  
  write.csv(turnover_stats(results$weights, asset_names),
            file.path(out_dir, "turnover_stats.csv"), row.names = FALSE)
  
  write.csv(allocation_summary(results$weights, asset_names),
            file.path(out_dir, "allocation_summary.csv"), row.names = FALSE)
  
  write.csv(concentration_stats(results$weights, asset_names),
            file.path(out_dir, "concentration.csv"), row.names = FALSE)
  
  write.csv(regime_performance(results$pnl, results$tpl_label),
            file.path(out_dir, "regime_performance.csv"), row.names = FALSE)
  
  write.csv(asset_sharpe_by_regime(benchmarks$returns_test, results$tpl_label, asset_names),
            file.path(out_dir, "asset_sharpe_by_regime.csv"), row.names = FALSE)
  
  write.csv(results$weights, file.path(out_dir, "daily_weights.csv"), row.names = FALSE)
  
  daily_pnl <- data.frame(
    date = results$pnl$date,
    wasserstein_hmm = pnl_strat,
    equal_weight    = pnl_ew,
    sixty_forty     = pnl_6040,
    spx_bh          = pnl_spx
  )
  write.csv(daily_pnl, file.path(out_dir, "daily_pnl.csv"), row.names = FALSE)
  
  write.csv(results$K_history, file.path(out_dir, "K_history.csv"), row.names = FALSE)
  write.csv(results$tpl_label, file.path(out_dir, "tpl_label_history.csv"), row.names = FALSE)
  
  perf
}

# ---------------- Plots ----------------

# Color palette tuned for clarity (colorblind-friendly-ish).
.PALETTE <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
              "#17becf", "#bcbd22", "#e377c2")

#' Save a PNG. Closes the device. width/height in inches.
.save_png <- function(file, expr, width = 10, height = 5, res = 150) {
  png(file, width = width, height = height, units = "in", res = res)
  on.exit(dev.off(), add = TRUE)
  expr
}

plot_cumulative_comparison <- function(pnl_list, file) {
  dates <- pnl_list[[1]]$date
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    cumvals <- lapply(pnl_list, function(df) cumsum(df$pnl))
    ymax <- max(sapply(cumvals, max)); ymin <- min(sapply(cumvals, min))
    plot(dates, cumvals[[1]], type = "l", ylim = c(ymin, ymax),
         col = .PALETTE[1], lwd = 2,
         main = "Cumulative OOS log return — strategy vs benchmarks",
         xlab = "", ylab = "Cumulative log return")
    grid(col = "gray85")
    for (i in seq_along(cumvals)[-1]) {
      lines(dates, cumvals[[i]], col = .PALETTE[i], lwd = 2, lty = i)
    }
    legend("topleft", names(pnl_list), col = .PALETTE[seq_along(cumvals)],
           lwd = 2, lty = seq_along(cumvals), bty = "n", cex = 0.9)
  })
}

plot_cumulative_by_regime <- function(pnl_df, label_df, file) {
  m <- merge(pnl_df, label_df[, c("date","label")], by = "date")
  m <- m[order(m$date), ]
  m$cum <- cumsum(m$pnl)
  regs  <- sort(unique(m$label[!is.na(m$label)]))
  cols  <- .PALETTE[seq_along(regs)]; names(cols) <- as.character(regs)
  
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    plot(m$date, m$cum, type = "n",
         main = "Wasserstein HMM cumulative P&L, colored by template regime",
         xlab = "", ylab = "Cumulative log return")
    grid(col = "gray85")
    # plot points colored by regime, plus connecting line in gray
    lines(m$date, m$cum, col = "gray60", lwd = 0.6)
    for (g in regs) {
      ix <- which(m$label == g)
      points(m$date[ix], m$cum[ix], col = cols[as.character(g)], pch = 16, cex = 0.6)
    }
    legend("topleft", paste("Regime", regs), col = cols, pch = 16, bty = "n",
           cex = 0.9)
  })
}

plot_stacked_weights <- function(weights_df, asset_names, file, title) {
  W <- as.matrix(weights_df[, asset_names])
  dates <- weights_df$date
  cols <- .PALETTE[seq_along(asset_names)]; names(cols) <- asset_names
  
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    plot(range(dates), c(0, 1), type = "n",
         main = title, xlab = "", ylab = "Portfolio weight")
    # build cumulative top-of-band per asset
    cum <- matrix(0, nrow(W), length(asset_names) + 1)
    for (j in seq_along(asset_names)) cum[, j + 1] <- cum[, j] + W[, j]
    for (j in seq_along(asset_names)) {
      polygon(c(dates, rev(dates)),
              c(cum[, j], rev(cum[, j + 1])),
              col = cols[j], border = NA)
    }
    legend("topleft", asset_names, fill = cols, bty = "n",
           horiz = TRUE, cex = 0.85)
  })
}

plot_turnover <- function(weights_df, asset_names, file, title) {
  W <- as.matrix(weights_df[, asset_names])
  dW <- apply(W, 2, function(v) c(NA, diff(v)))
  TO <- 0.5 * rowSums(abs(dW))
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    plot(weights_df$date, TO, type = "l", col = .PALETTE[1], lwd = 1,
         main = title, xlab = "", ylab = "Daily turnover (0.5 * L1)")
    grid(col = "gray85")
  })
}

plot_neff <- function(weights_df, asset_names, file, title) {
  W <- as.matrix(weights_df[, asset_names])
  ne <- 1 / pmax(rowSums(W^2), 1e-12)
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    plot(weights_df$date, ne, type = "l", col = .PALETTE[3], lwd = 1.2,
         main = title, xlab = "", ylab = expression(N[eff]))
    grid(col = "gray85")
  })
}

plot_K_history <- function(K_df, file) {
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    plot(K_df$date, K_df$K, type = "s", col = .PALETTE[4], lwd = 1.2,
         main = "Selected K over time (predictive selection)",
         xlab = "", ylab = "K", ylim = range(K_df$K, na.rm = TRUE) + c(-0.5, 0.5))
    grid(col = "gray85")
  })
}

plot_drawdown <- function(pnl_list, file) {
  dates <- pnl_list[[1]]$date
  dds <- lapply(pnl_list, function(df) { c <- cumsum(df$pnl); c - cummax(c) })
  
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    ymin <- min(sapply(dds, min)); ymax <- 0
    plot(dates, dds[[1]], type = "l", col = .PALETTE[1], lwd = 2,
         ylim = c(ymin, ymax),
         main = "Drawdown comparison (cum log return - running max)",
         xlab = "", ylab = "Drawdown (log return)")
    grid(col = "gray85"); abline(h = 0, col = "black")
    for (i in seq_along(dds)[-1])
      lines(dates, dds[[i]], col = .PALETTE[i], lwd = 2, lty = i)
    legend("bottomleft", names(pnl_list),
           col = .PALETTE[seq_along(dds)], lwd = 2, lty = seq_along(dds),
           bty = "n", cex = 0.9)
  })
}

plot_rolling_sharpe <- function(pnl_list, file, window = 60) {
  dates <- pnl_list[[1]]$date
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    plot(range(dates), c(-2, 6), type = "n",
         main = sprintf("Rolling %d-day annualized Sharpe", window),
         xlab = "", ylab = "Rolling Sharpe")
    grid(col = "gray85"); abline(h = 0, col = "black", lty = 2)
    for (i in seq_along(pnl_list)) {
      pnl <- pnl_list[[i]]$pnl
      m   <- caTools::runmean(pnl, window, endrule = "NA", align = "right")
      s   <- caTools::runsd(pnl, window, endrule = "NA", align = "right")
      rs  <- (m / s) * sqrt(TRADING_DAYS_PER_YEAR)
      lines(dates, rs, col = .PALETTE[i], lwd = 1.5, lty = i)
    }
    legend("topleft", names(pnl_list), col = .PALETTE[seq_along(pnl_list)],
           lwd = 1.5, lty = seq_along(pnl_list), bty = "n", cex = 0.85)
  })
}

plot_asset_sharpe_by_regime <- function(reg_tbl, asset_names, file) {
  vals <- as.matrix(reg_tbl[, asset_names])
  regs <- reg_tbl$regime
  rownames(vals) <- paste("Regime", regs)
  cols <- .PALETTE[seq_along(asset_names)]
  
  .save_png(file, width = 10, height = 5, expr = {
    par(mar = c(4, 4, 3, 1), las = 1)
    barplot(t(vals), beside = TRUE, col = cols,
            main = "Per-regime annualized Sharpe by asset",
            ylab = "Sharpe ratio", legend.text = asset_names,
            args.legend = list(x = "topright", bty = "n", cex = 0.85))
    abline(h = 0)
  })
}

plot_stacked_cum_by_regime <- function(pnl_df, label_df, file) {
  m <- merge(pnl_df, label_df[, c("date","label")], by = "date")
  m <- m[order(m$date), ]
  regs <- sort(unique(m$label[!is.na(m$label)]))
  
  # cumulative PnL contribution by regime
  M <- matrix(0, nrow(m), length(regs))
  colnames(M) <- as.character(regs)
  for (g in regs) {
    M[, as.character(g)] <- ifelse(m$label == g & !is.na(m$label), m$pnl, 0)
  }
  cumM <- apply(M, 2, cumsum)
  
  cols <- .PALETTE[seq_along(regs)]
  .save_png(file, {
    par(mar = c(4, 4, 3, 1), las = 1)
    plot(range(m$date), c(0, sum(m$pnl)), type = "n",
         main = "Stacked cumulative P&L attribution by template regime",
         xlab = "", ylab = "Cumulative log return")
    cum_top <- matrix(0, nrow(m), length(regs) + 1)
    for (j in seq_along(regs)) cum_top[, j + 1] <- cum_top[, j] + cumM[, j]
    for (j in seq_along(regs)) {
      polygon(c(m$date, rev(m$date)),
              c(cum_top[, j], rev(cum_top[, j + 1])),
              col = cols[j], border = NA)
    }
    legend("topleft", paste("Regime", regs), fill = cols, bty = "n",
           horiz = TRUE, cex = 0.85)
  })
}

# ---------------- Aggregate export ----------------

export_all_plots <- function(results, benchmarks, asset_names, out_dir) {
  pnl_list <- list(
    "Wasserstein HMM + MVO" = results$pnl,
    "Equal Weight (20%)"    = benchmarks$ew,
    "60/40 Stocks/Bonds"    = benchmarks$sixty_forty,
    "SPX Buy & Hold"        = benchmarks$spx
  )
  
  plot_cumulative_comparison(pnl_list,
    file.path(out_dir, "01_cumulative_comparison.png"))
  
  plot_cumulative_by_regime(results$pnl, results$tpl_label,
    file.path(out_dir, "02_cumulative_by_regime.png"))
  
  plot_stacked_weights(results$weights, asset_names,
    file.path(out_dir, "03_weights_stacked.png"),
    "Wasserstein HMM — daily portfolio weights")
  
  plot_turnover(results$weights, asset_names,
    file.path(out_dir, "04_turnover.png"),
    "Wasserstein HMM — daily turnover")
  
  plot_neff(results$weights, asset_names,
    file.path(out_dir, "05_neff.png"),
    "Wasserstein HMM — effective number of positions")
  
  plot_K_history(results$K_history,
    file.path(out_dir, "06_K_history.png"))
  
  plot_drawdown(pnl_list,
    file.path(out_dir, "07_drawdown_comparison.png"))
  
  plot_rolling_sharpe(pnl_list,
    file.path(out_dir, "08_rolling_sharpe.png"))
  
  reg_asset <- asset_sharpe_by_regime(benchmarks$returns_test,
                                      results$tpl_label, asset_names)
  plot_asset_sharpe_by_regime(reg_asset, asset_names,
    file.path(out_dir, "09_asset_sharpe_by_regime.png"))
  
  plot_stacked_cum_by_regime(results$pnl, results$tpl_label,
    file.path(out_dir, "10_stacked_cum_by_regime.png"))
}
