# 05_plots.R -- reproduces figures 1-13 from the paper (adapted to FF factors).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ---- Fig 1: PCA cumulative variance -----------------------------------------

plot_pca_variance <- function(pca_fit, threshold = 0.95) {
  df <- tibble(n = seq_along(pca_fit$cum_var), cum_var = pca_fit$cum_var)
  ggplot(df, aes(n, cum_var)) +
    geom_line(colour = "steelblue", linewidth = 0.7) +
    geom_hline(yintercept = threshold, linetype = "dashed", colour = "firebrick") +
    annotate("text", x = max(df$n) * 0.7, y = threshold + 0.02,
             label = paste0(threshold * 100, "% variance"), colour = "firebrick") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "# components", y = "cumulative explained variance",
         title = "PCA explained variance") +
    theme_minimal(base_size = 11)
}

# ---- Fig 2: modified k-means vs GMM regime labels over time -----------------

plot_regime_timeline <- function(dates, km_labels, gmm_labels = NULL,
                                 nber_recessions = NULL) {
  df <- tibble(date = dates, `Modified k-means` = km_labels)
  # only add GMM column if it's the right length and not all NA
  if (!is.null(gmm_labels) && length(gmm_labels) == nrow(df) &&
      !all(is.na(gmm_labels))) {
    df$GMM <- gmm_labels
  }
  long <- df |>
    pivot_longer(-date, names_to = "method", values_to = "regime")

  p <- ggplot(long, aes(date, regime, colour = factor(regime))) +
    facet_wrap(~ method, ncol = 1) +
    geom_point(size = 1.2) +
    scale_colour_brewer(palette = "Set1", name = "regime") +
    labs(x = NULL, y = "regime", title = "Regime classifications over time") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")

  if (!is.null(nber_recessions)) {
    p <- p + geom_rect(
      data = nber_recessions, inherit.aes = FALSE,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      fill = "grey60", alpha = 0.3)
  }
  p
}

# ---- Fig 3: crisis probability (k-means vs GMM) -----------------------------

plot_crisis_probabilities <- function(dates, km_probs_R0, gmm_probs_R0,
                                      nber_recessions = NULL) {
  has_gmm <- !is.null(gmm_probs_R0) && length(gmm_probs_R0) == length(dates) &&
             !all(is.na(gmm_probs_R0))
  df <- tibble(date = dates, `Modified k-means` = km_probs_R0)
  if (has_gmm) df$GMM <- gmm_probs_R0
  df <- df |> pivot_longer(-date, names_to = "method", values_to = "p_crisis")

  colour_vals <- c("Modified k-means" = "firebrick", GMM = "steelblue")
  p <- ggplot(df, aes(date, p_crisis, colour = method)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ method, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = colour_vals) +
    labs(x = NULL, y = "P(crisis / Regime 0)",
         title = if (has_gmm)
                   "Crisis probability: modified k-means vs GMM"
                 else
                   "Crisis probability: modified k-means") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  if (!is.null(nber_recessions)) {
    p <- p + geom_rect(data = nber_recessions, inherit.aes = FALSE,
                       aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                       fill = "grey60", alpha = 0.3)
  }
  p
}

# ---- Fig 4: FRED-MD per-regime averages heatmap -----------------------------

plot_regime_heatmap <- function(macro_raw, labels, key_vars = NULL) {
  if (is.null(key_vars)) {
    candidates <- c(
      # activity
      "RPI", "W875RX1", "INDPRO", "IPFINAL", "IPCONGD", "IPMANSICS",
      "DPCERA3M086SBEA", "CMRMTSPLx", "RETAILx",
      # labour
      "UNRATE", "PAYEMS", "CES0600000007", "CES0600000008",
      "CLAIMSx", "HWIx",
      # prices & sentiment
      "CPIAUCSL", "CPIULFSL", "PCEPI", "WPSFD49207",
      "UMCSENTx",
      # housing
      "HOUST", "HOUSTS", "PERMIT",
      # credit / money
      "M2SL", "M2REAL", "BUSLOANS", "CONSPI",
      # financial / market
      "S&P 500", "S&P: indust", "S&P div yield", "S&P PE ratio",
      "VIXCLSx"
    )
    key_vars <- intersect(candidates, colnames(macro_raw))
    # cap at 12 for readability
    if (length(key_vars) > 12) key_vars <- key_vars[seq_len(12)]
  }
  if (!length(key_vars)) {
    warning("No known macro variables in the data for heatmap.")
    return(ggplot() + theme_void() +
           labs(title = "Regime heatmap: no target vars available"))
  }
  df <- macro_raw |>
    mutate(regime = labels) |>
    select(regime, all_of(key_vars)) |>
    pivot_longer(-regime, names_to = "variable", values_to = "value") |>
    group_by(regime, variable) |>
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") |>
    group_by(variable) |>
    mutate(mean_norm = (mean_val - min(mean_val)) /
                       pmax(max(mean_val) - min(mean_val), 1e-12)) |>
    ungroup()

  ggplot(df, aes(factor(regime), variable, fill = mean_norm)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", mean_norm)), size = 3) +
    scale_fill_viridis_c(option = "D", limits = c(0, 1),
                         name = "min-max\nnormalised") +
    labs(x = "regime", y = NULL,
         title = "Min-max normalised regime statistics (tcode-transformed)") +
    theme_minimal(base_size = 11)
}

# ---- Fig 5: transition matrix (raw + conditional) ---------------------------

plot_transition_matrix <- function(E, title = "Regime transition matrix") {
  df <- as.data.frame(as.table(E))
  names(df) <- c("from", "to", "prob")
  ggplot(df, aes(from, to, fill = prob)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", prob)), size = 3) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1)) +
    labs(x = "from regime", y = "to regime", title = title) +
    theme_minimal(base_size = 11) +
    coord_equal()
}

# ---- Fig 6: transition network graph ----------------------------------------

plot_transition_network <- function(E_cond, regime_labels = NULL) {
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
  if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
  n <- nrow(E_cond)
  regimes <- rownames(E_cond)
  if (is.null(regime_labels)) {
    regime_labels <- setNames(paste("Regime", regimes), regimes)
  }
  edges <- expand.grid(from = regimes, to = regimes,
                       stringsAsFactors = FALSE) |>
    mutate(weight = as.numeric(E_cond)) |>
    filter(from != to, weight > 0.01)

  g <- igraph::graph_from_data_frame(edges,
                                     vertices = data.frame(
                                       name  = regimes,
                                       label = regime_labels))

  ggraph::ggraph(g, layout = "circle") +
    ggraph::geom_edge_link(aes(alpha = weight, width = weight),
                           arrow = arrow(length = unit(2, "mm")),
                           end_cap = ggraph::circle(4, "mm")) +
    ggraph::scale_edge_alpha(range = c(0.15, 0.9), guide = "none") +
    ggraph::scale_edge_width(range = c(0.3, 2), name = "P(j | i, transition)") +
    ggraph::geom_node_point(size = 14, colour = "seagreen") +
    ggraph::geom_node_text(aes(label = name), colour = "white",
                           fontface = "bold") +
    ggraph::geom_node_text(aes(label = label), nudge_y = 0.15,
                           size = 3.2) +
    theme_void() +
    labs(title = "Regime transition network (conditional on transition)")
}

# ---- Figs 7-9: boxplots for random vs non-random regime comparisons ---------

plot_metric_boxplot <- function(control, treatment, metric_name,
                                 title = NULL, p_value = NULL) {
  df <- tibble(
    group = rep(c("Control (random)", "Treatment"),
                c(length(control), length(treatment))),
    value = c(control, treatment)
  )
  p <- ggplot(df, aes(group, value)) +
    geom_boxplot(width = 0.5, outlier.size = 1) +
    geom_jitter(width = 0.1, size = 1.2, alpha = 0.7) +
    labs(x = NULL, y = metric_name,
         title = title %||% metric_name) +
    theme_minimal(base_size = 11)
  if (!is.null(p_value)) {
    p <- p + annotate("text", x = 1.5, y = min(df$value, na.rm = TRUE),
                      label = sprintf("p = %.3f", p_value),
                      vjust = -0.5, size = 3.5)
  }
  p
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Figs 10-13: cumulative log returns -------------------------------------

plot_cumulative_returns <- function(portfolio_returns_scaled,
                                    groups = NULL, title = NULL) {
  df <- portfolio_returns_scaled |>
    pivot_longer(-date, names_to = "strategy", values_to = "ret") |>
    group_by(strategy) |>
    arrange(date) |>
    mutate(cum_log = cumsum(replace_na(log(1 + ret), 0))) |>
    ungroup()

  if (!is.null(groups)) {
    df <- df |> mutate(group = groups[strategy])
    df <- df |> filter(!is.na(group))
    p <- ggplot(df, aes(date, cum_log, colour = group, group = strategy)) +
      geom_line(alpha = 0.55, linewidth = 0.5)
  } else {
    p <- ggplot(df, aes(date, cum_log, colour = strategy)) +
      geom_line(alpha = 0.7, linewidth = 0.5)
  }
  p +
    labs(x = NULL, y = "cumulative log return",
         title = title %||% "Cumulative log returns (10% vol target)") +
    theme_minimal(base_size = 11)
}

strategy_groups <- function(strategy_names) {
  prefix <- function(s) {
    # mvo, ew, spy are single strategies; others are prefix_suffix_l
    if (s %in% c("mvo", "ew", "spy")) return(s)
    sub("_.*", "", s)
  }
  setNames(vapply(strategy_names, prefix, character(1)), strategy_names)
}

# ---- generator: all plots ---------------------------------------------------

generate_all_plots <- function(pca_fit,
                               dates, km_labels, gmm_labels, km_probs_R0, gmm_probs_R0,
                               macro_raw_aligned,
                               transition_matrix, transition_matrix_cond,
                               regime_labels_table,
                               portfolio_returns_scaled,
                               nber_recessions = NULL) {

  plots <- list()
  plots$fig01_pca_variance <- plot_pca_variance(pca_fit)
  plots$fig02_regime_timeline <- plot_regime_timeline(
    dates, km_labels, gmm_labels, nber_recessions)
  plots$fig03_crisis_probabilities <- plot_crisis_probabilities(
    dates, km_probs_R0, gmm_probs_R0, nber_recessions)
  plots$fig04_regime_heatmap <- plot_regime_heatmap(
    macro_raw_aligned, km_labels)
  plots$fig05a_transition_matrix <- plot_transition_matrix(
    transition_matrix, "Regime transition matrix")
  plots$fig05b_transition_matrix_cond <- plot_transition_matrix(
    transition_matrix_cond, "Conditional transition matrix (given transition)")
  plots$fig06_transition_network <- plot_transition_network(
    transition_matrix_cond, regime_labels_table)

  gmap <- strategy_groups(setdiff(colnames(portfolio_returns_scaled), "date"))
  plots$fig10_cumulative_all <- plot_cumulative_returns(
    portfolio_returns_scaled, groups = gmap,
    title = "Cumulative log returns, all strategies (10% vol target)")

  naive_cols <- c("date", grep("^(naive|ew|mvo|spy)", colnames(portfolio_returns_scaled), value = TRUE))
  plots$fig11_naive_vs_benchmarks <- plot_cumulative_returns(
    portfolio_returns_scaled |> select(any_of(naive_cols)),
    groups = gmap, title = "Naive regime portfolio vs benchmarks")

  bl_cols <- c("date", grep("^(bl|mvo|ew|spy)", colnames(portfolio_returns_scaled), value = TRUE))
  plots$fig12_bl_vs_benchmarks <- plot_cumulative_returns(
    portfolio_returns_scaled |> select(any_of(bl_cols)),
    groups = gmap, title = "Black-Litterman regime portfolio vs benchmarks")

  ridge_cols <- c("date", grep("^(ridge|mvo|ew|spy)", colnames(portfolio_returns_scaled), value = TRUE))
  plots$fig13_ridge_vs_benchmarks <- plot_cumulative_returns(
    portfolio_returns_scaled |> select(any_of(ridge_cols)),
    groups = gmap, title = "Ridge regime portfolio vs benchmarks")

  plots
}

# ---- convenience: save all --------------------------------------------------

save_plots <- function(plots, out_dir = "output",
                       width = 10, height = 5.5, dpi = 140) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  for (nm in names(plots)) {
    ggplot2::ggsave(file.path(out_dir, paste0(nm, ".png")),
                    plots[[nm]], width = width, height = height, dpi = dpi)
  }
  invisible(plots)
}
