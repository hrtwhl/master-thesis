# ---- Small helpers --------------------------------------------------------

#' End-of-month sample of a daily time series
#' Takes the last non-NA observation within each calendar month.
to_monthly_eom <- function(df, date_col = "date") {
  df %>%
    mutate(ym = floor_date(.data[[date_col]], "month")) %>%
    group_by(ym) %>%
    summarise(across(-all_of(date_col), ~ {
      v <- .[!is.na(.)]
      if (length(v) == 0) NA_real_ else tail(v, 1)
    }), .groups = "drop") %>%
    rename(date = ym) %>%
    mutate(date = ceiling_date(date, "month") - days(1))  # anchor at month-end
}

#' Winsorize a numeric vector at symmetric bounds
winsorize <- function(x, bound = 3) pmin(pmax(x, -bound), bound)

#' Annualised volatility from monthly log-returns (or simple returns)
ann_vol <- function(r) sd(r, na.rm = TRUE) * sqrt(12)

#' Annualised Sharpe ratio for monthly returns
sharpe <- function(r) {
  r <- r[is.finite(r)]
  if (length(r) < 6 || sd(r) == 0) return(NA_real_)
  (mean(r) / sd(r)) * sqrt(12)
}

#' Relative drawdown on a compounded equity curve.
#' Input: monthly simple returns (e.g. 0.02 for +2%). NAs are treated as 0.
#' Output: a numeric vector in [-1, 0] where -0.30 means "down 30% from the
#' equity-curve high-water mark." This is the convention the paper uses
#' (Exhibit 11 reads -50% as "capital dropped 50% from peak").
#'
#' For guardrails, if compounded equity goes non-positive (extremely negative
#' returns - shouldn't happen with our vol-targeted monthly series), we floor
#' at -1 rather than let NaN propagate.
drawdown <- function(r) {
  r <- ifelse(is.na(r), 0, r)
  eq <- cumprod(1 + r)
  peak <- cummax(eq)
  dd <- eq / peak - 1
  dd[!is.finite(dd)] <- -1
  dd
}

#' Absolute (P&L-sum) drawdown - the difference between cumulative P&L and
#' its running maximum. Useful for long-short strategies expressed as
#' cumulative return in percentage points, but NOT the conventional
#' "drawdown" for an equity curve.
drawdown_abs <- function(r) {
  eq <- cumsum(ifelse(is.na(r), 0, r))
  eq - cummax(eq)
}

#' Simple ggplot theme used for all exhibits
theme_regimes <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey30"),
      strip.text    = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_plot <- function(p, name, width = 9, height = 5, dpi = 150) {
  path <- file.path(CFG$paths$output_dir, paste0(name, ".png"))
  ggsave(path, p, width = width, height = height, dpi = dpi)
  invisible(path)
}

save_table <- function(tbl, name) {
  path <- file.path(CFG$paths$output_dir, paste0(name, ".csv"))
  write_csv(tbl, path)
  invisible(path)
}
