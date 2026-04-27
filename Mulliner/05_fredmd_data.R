# ---------------------------------------------------------------------------
# R/05_fredmd_data.R
#
# Load and transform the FRED-MD database following McCracken & Ng (2016).
#
# The CSV downloaded from
#   https://www.stlouisfed.org/research/economists/mccracken/fred-databases
# has two header rows:
#   row 1: column names (first column is "sasdate")
#   row 2: integer transformation codes per series (1-7), labelled "Transform:"
# Followed by monthly observations.
#
# Transformation codes (McCracken-Ng):
#   1: x_t                          (no transformation)
#   2: x_t - x_{t-1}                (first difference)
#   3: x_t - 2 x_{t-1} + x_{t-2}    (second difference)
#   4: log(x_t)
#   5: log(x_t) - log(x_{t-1})      (log first difference)
#   6: log(x_t) - 2 log(x_{t-1}) + log(x_{t-2})   (log second difference)
#   7: (x_t / x_{t-1}) - 1          (rate of change)
#
# Outputs:
#   - load_fredmd(path)            : raw data + tcodes + dates
#   - transform_fredmd(raw)        : applies tcodes, returns stationary panel
# ---------------------------------------------------------------------------

#' Read FRED-MD CSV and split header rows.
#'
#' Returns a list with:
#'   raw    : tibble of date + raw levels (untransformed)
#'   tcodes : named integer vector of transformation codes per variable
#'   dates  : end-of-month dates aligned with rows of `raw`
load_fredmd <- function(path) {
  # Read header rows separately so we can capture the transformation row.
  header_lines <- readLines(path, n = 2)
  col_names <- strsplit(header_lines[1], ",")[[1]]
  col_names <- gsub("\r|\"", "", col_names)
  tcode_tokens <- strsplit(header_lines[2], ",")[[1]]
  tcode_tokens <- gsub("\r|\"", "", tcode_tokens)
  # First column of row 2 is the literal string "Transform:"; the rest are codes
  tcodes <- as.integer(tcode_tokens[-1])
  names(tcodes) <- col_names[-1]

  # Read the data rows (skip the two header lines)
  data_rows <- read_csv(path, skip = 2, col_names = col_names,
                        show_col_types = FALSE,
                        col_types = cols(.default = col_double(),
                                         sasdate = col_character()))

  # FRED-MD dates use M/D/YYYY format (e.g., "1/1/1959") and indicate the
  # first day of each month. We anchor them at end-of-month for compatibility
  # with the rest of the pipeline.
  dates <- as.Date(data_rows$sasdate, format = "%m/%d/%Y")
  dates <- ceiling_date(dates, "month") - days(1)

  raw <- data_rows %>%
    select(-sasdate) %>%
    mutate(date = dates, .before = 1) %>%
    arrange(date)

  list(
    raw    = raw,
    tcodes = tcodes,
    dates  = sort(dates),
    n_vars = length(tcodes)
  )
}

#' Apply McCracken-Ng transformation code `tc` to a numeric vector.
#'
#' For codes 4-6 we require positive values; non-positive entries are returned
#' as NA, matching how McCracken-Ng treat such cases (they avoid log() on
#' non-positive series in their published transformations).
apply_tcode <- function(x, tc) {
  # Safe log: NA on non-positive values
  safe_log <- function(v) {
    out <- rep(NA_real_, length(v))
    ok  <- is.finite(v) & v > 0
    out[ok] <- log(v[ok])
    out
  }

  switch(as.character(tc),
    "1" = x,
    "2" = c(NA_real_, diff(x, lag = 1)),
    "3" = c(NA_real_, NA_real_, diff(x, lag = 1, differences = 2)),
    "4" = safe_log(x),
    "5" = c(NA_real_, diff(safe_log(x), lag = 1)),
    "6" = c(NA_real_, NA_real_, diff(safe_log(x), lag = 1, differences = 2)),
    "7" = c(NA_real_, x[-1] / x[-length(x)] - 1),
    stop("Unknown tcode: ", tc)
  )
}

#' Apply the transformation codes to every series and clip to a sample
#' window. Returns a tibble keyed by `date` with one column per variable
#' containing the stationary transformed series.
#'
#' By default, drops the rows lost to the longest transformation (2 obs lost
#' for tcodes 3 and 6) so the panel is well-defined from the start.
transform_fredmd <- function(fmd, sample_end = NULL, drop_initial_NAs = TRUE) {
  vars <- setdiff(names(fmd$raw), "date")
  out  <- fmd$raw %>% select(date)

  for (v in vars) {
    tc <- fmd$tcodes[[v]]
    out[[v]] <- apply_tcode(fmd$raw[[v]], tc)
  }

  if (!is.null(sample_end)) {
    out <- out %>% filter(date <= as.Date(sample_end))
  }

  if (drop_initial_NAs) {
    # Drop the first 2 rows (where 2nd-difference codes 3 and 6 are NA by
    # construction). We do NOT drop later NAs - those are real missing data
    # for individual series, handled by the EM-PCA in the next step.
    out <- out[-c(1, 2), , drop = FALSE]
  }

  out
}

#' McCracken-Ng outlier removal (optional, not applied by default).
#'
#' Definition: a point in the transformed series is an outlier if its
#' deviation from the median exceeds 10 times the interquartile range.
#' Outliers are set to NA so the EM-PCA imputes them.
remove_outliers_mn <- function(transformed, factor = 10) {
  out <- transformed
  vars <- setdiff(names(out), "date")
  n_outlier <- 0L
  for (v in vars) {
    x <- out[[v]]
    if (sum(!is.na(x)) < 12) next
    med <- median(x, na.rm = TRUE)
    iqr <- IQR(x, na.rm = TRUE)
    is_out <- !is.na(x) & abs(x - med) > factor * iqr
    n_outlier <- n_outlier + sum(is_out)
    out[[v]][is_out] <- NA_real_
  }
  attr(out, "n_outliers") <- n_outlier
  out
}

#' Quick descriptive summary of the FRED-MD panel: how many variables,
#' how much missingness, the date range, and a per-tcode count.
describe_fredmd <- function(fmd, transformed) {
  vars <- setdiff(names(transformed), "date")
  na_pct <- sapply(vars, function(v) mean(is.na(transformed[[v]])))
  tcode_counts <- table(fmd$tcodes[vars])

  cat("FRED-MD panel summary\n")
  cat("  Variables          : ", length(vars), "\n")
  cat("  Date range (raw)   : ", format(min(fmd$dates)), "to",
      format(max(fmd$dates)), "\n")
  cat("  Date range (xform) : ", format(min(transformed$date)), "to",
      format(max(transformed$date)), "\n")
  cat("  Variables with > 5%  NA: ", sum(na_pct > 0.05),  "\n")
  cat("  Variables with > 25% NA: ", sum(na_pct > 0.25), "\n")
  cat("  Variables with > 50% NA: ", sum(na_pct > 0.50), "\n")
  cat("  Transformation code counts:\n")
  for (tc in names(tcode_counts))
    cat(sprintf("    tcode %s : %d series\n", tc, tcode_counts[tc]))
  invisible(NULL)
}
