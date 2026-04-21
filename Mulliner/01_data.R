# ---- Macro data pipeline and factor loading -------------------------------

#' Load the user-provided daily macro CSV, clip to sample end, and resample to
#' month-end.
load_macro <- function(path, sample_end) {
  raw <- read_csv(path, show_col_types = FALSE) %>%
    mutate(date = as.Date(date)) %>%
    filter(date <= sample_end) %>%
    arrange(date)

  # End-of-month snapshot
  to_monthly_eom(raw, "date")
}

#' Apply the z-score-like transformation described in the paper:
#'   - take 12-month change
#'   - divide by std of the last 10 years of 12-month changes
#'   - winsorize at +/- winsorize
#' diff_type = "level" uses x_t - x_{t-horizon};
#' diff_type = "log"   uses log(x_t) - log(x_{t-horizon}).
build_state_variables <- function(macro_m, variables, horizon, window, winsorize) {
  vars <- variables$var

  # Safe log: returns NA for non-positive input rather than NaN / -Inf.
  safe_log <- function(x) {
    out <- rep(NA_real_, length(x))
    ok  <- is.finite(x) & x > 0
    out[ok] <- log(x[ok])
    out
  }

  out_long <- macro_m %>%
    pivot_longer(-date, names_to = "var", values_to = "x") %>%
    left_join(variables, by = "var") %>%
    group_by(var) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      x_used = if (unique(diff_type) == "log") safe_log(x) else x,
      diff   = x_used - lag(x_used, horizon),
      # rolling std of the 12m difference over the last `window` months
      roll_sd = slider::slide_dbl(diff, ~ sd(.x, na.rm = TRUE),
                                  .before = window - 1, .after = 0,
                                  .complete = FALSE),
      # need a full window of diffs to declare a valid z-score
      diff_nonna_count = slider::slide_dbl(diff, ~ sum(!is.na(.x)),
                                           .before = window - 1, .after = 0,
                                           .complete = FALSE),
      roll_sd = if_else(diff_nonna_count >= window, roll_sd, NA_real_),
      z_raw   = diff / roll_sd,
      z_wins  = winsorize(z_raw, winsorize)
    ) %>%
    ungroup()

  # Wide matrices keyed by date
  transformed <- out_long %>%
    select(date, var, z_raw) %>%
    pivot_wider(names_from = var, values_from = z_raw) %>%
    arrange(date) %>%
    select(date, all_of(vars))

  transformed_w <- out_long %>%
    select(date, var, z_wins) %>%
    pivot_wider(names_from = var, values_from = z_wins) %>%
    arrange(date) %>%
    select(date, all_of(vars))

  list(
    monthly                 = macro_m,
    long                    = out_long,
    transformed             = transformed,
    transformed_winsorized  = transformed_w
  )
}

# ---- Fama-French factors --------------------------------------------------

#' Download, unzip, parse a Ken French monthly factor CSV.
#' The Dartmouth CSVs have a text preamble and a second (annual) block after a
#' blank row; this parses the monthly block only.
fetch_french_csv <- function(url) {
  zip_tmp <- tempfile(fileext = ".zip")
  utils::download.file(url, zip_tmp, mode = "wb", quiet = TRUE)
  ex_dir <- tempfile(); dir.create(ex_dir)
  utils::unzip(zip_tmp, exdir = ex_dir)
  csv_file <- list.files(ex_dir, pattern = "\\.csv$", full.names = TRUE)[1]
  lines <- readLines(csv_file)

  # find the first row that looks like a header (contains the first factor)
  header_idx <- which(grepl("^\\s*,.*Mkt-?RF|,Mom|,WML", lines))[1]
  # data starts the line after the header
  # first blank line after the header marks end of monthly block
  tail_lines <- lines[(header_idx + 1):length(lines)]
  blank_idx  <- which(trimws(tail_lines) == "" |
                      grepl("Annual", tail_lines, ignore.case = TRUE))
  end_idx <- if (length(blank_idx) > 0) blank_idx[1] - 1 else length(tail_lines)
  monthly_block <- c(lines[header_idx], tail_lines[1:end_idx])

  df <- read_csv(paste(monthly_block, collapse = "\n"),
                 show_col_types = FALSE)
  names(df)[1] <- "ym"
  df %>%
    filter(!is.na(ym), nchar(as.character(ym)) == 6) %>%
    mutate(
      ym_str = sprintf("%06d", as.integer(ym)),
      date   = as.Date(paste0(ym_str, "01"), "%Y%m%d"),
      date   = ceiling_date(date, "month") - days(1)
    ) %>%
    select(-ym, -ym_str) %>%
    mutate(across(-date, ~ as.numeric(.x) / 100))   # FF reports percent
}

#' Load six long-short factors: Mkt-RF, SMB, HML, RMW, CMA, Mom
load_factors <- function(ff5_url, mom_url, sample_end) {
  ff5 <- fetch_french_csv(ff5_url)
  mom <- fetch_french_csv(mom_url)
  # momentum column is named "Mom" with surrounding spaces in some files
  mom_col <- setdiff(names(mom), "date")[1]
  names(mom)[names(mom) == mom_col] <- "Mom"

  df <- ff5 %>%
    select(date, `Mkt-RF`, SMB, HML, RMW, CMA) %>%
    inner_join(mom %>% select(date, Mom), by = "date") %>%
    filter(date <= sample_end) %>%
    arrange(date)

  names(df) <- c("date", "Market", "Size", "Value", "Profitability",
                 "Investment", "Momentum")
  df
}
