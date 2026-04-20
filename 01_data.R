# 01_data.R -- FRED-MD macro data + sector ETF returns (Oliveira et al. 2025,
# Table 2).
#
# Universe: SPY + 9 SPDR sector ETFs (XLB, XLE, XLF, XLI, XLK, XLP, XLU, XLV,
# XLY). All have clean monthly data from 1998-12 (SPY from 1993-01).
#
# Notes on look-ahead:
#   - FRED-MD uses the CURRENT vintage (not ALFRED real-time).
#   - Yahoo adjusted close is the current split/dividend-adjusted series, so
#     survivorship and adjustment revisions are implicit.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
  library(tibble)
})

# ---- dependency bootstrap ----------------------------------------------------

ensure_pkg <- function(pkg, github = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible())
  if (!is.null(github)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github(github, upgrade = "never")
  } else {
    install.packages(pkg)
  }
}

# =============================================================================
# FRED-MD
# =============================================================================

#' Quick structural check that a candidate file really is a FRED-MD CSV.
is_valid_fredmd_csv <- function(path) {
  ok <- tryCatch({
    first_two <- readLines(path, n = 2L, warn = FALSE)
    if (length(first_two) < 2L) return(FALSE)
    header <- strsplit(first_two[1], ",", fixed = TRUE)[[1]]
    tcodes <- strsplit(first_two[2], ",", fixed = TRUE)[[1]]
    if (length(header) < 100L || length(tcodes) < 100L) return(FALSE)
    h1 <- tolower(gsub("[^a-z]", "", header[1]))
    if (!h1 %in% c("sasdate", "date")) return(FALSE)
    tc_num <- suppressWarnings(as.integer(tcodes[-1]))
    frac_valid <- mean(tc_num >= 1L & tc_num <= 7L, na.rm = TRUE)
    isTRUE(frac_valid > 0.9)
  }, error = function(e) FALSE)
  isTRUE(ok)
}

#' Walk back through recent FRED-MD monthly vintages until one validates.
download_fredmd_csv <- function(url = NULL, max_lookback = 6L, quiet = FALSE) {
  if (!is.null(url)) urls <- url else {
    base <- "https://www.stlouisfed.org/-/media/project/frbstl/stlouisfed/research/fred-md/monthly"
    months <- seq(Sys.Date(), by = "-1 month", length.out = max_lookback + 1L)
    urls <- sprintf("%s/%s-md.csv", base, format(months, "%Y-%m"))
  }
  for (u in urls) {
    if (!quiet) message("  trying ", u)
    dest <- tempfile(fileext = ".csv")
    ok <- tryCatch({
      utils::download.file(u, destfile = dest, mode = "wb", quiet = TRUE)
      file.exists(dest) && file.size(dest) > 500000L
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!isTRUE(ok)) {
      if (!quiet) message("    skipped (download failed or file too small)")
      next
    }
    if (!is_valid_fredmd_csv(dest)) {
      if (!quiet) message("    skipped (not a valid FRED-MD; ",
                          round(file.size(dest) / 1024), " KB)")
      next
    }
    if (!quiet) message("  -> OK (", round(file.size(dest) / 1024), " KB)")
    return(dest)
  }
  stop("Could not download a valid FRED-MD CSV. Pass url = '...' manually.")
}

#' Apply FRED-MD tcode transformations to a single series.
apply_tcode <- function(x, tcode) {
  tcode <- as.integer(tcode)
  safe_log <- function(v) { v[!is.finite(v) | v <= 0] <- NA; log(v) }
  switch(as.character(tcode),
         "1" = x,
         "2" = c(NA, diff(x)),
         "3" = c(NA, NA, diff(x, differences = 2)),
         "4" = safe_log(x),
         "5" = c(NA, diff(safe_log(x))),
         "6" = c(NA, NA, diff(safe_log(x), differences = 2)),
         "7" = c(NA, diff(x / dplyr::lag(x) - 1)),
         x)
}

#' Parse a FRED-MD CSV directly (header + tcode row + data rows).
parse_fredmd_csv <- function(path, transform = TRUE) {
  all_chr <- suppressWarnings(
    readr::read_csv(path, col_names = FALSE,
                    col_types = readr::cols(.default = "c"))
  )
  if (nrow(all_chr) < 5L) stop("FRED-MD file has fewer than 5 rows.")

  header <- as.character(unlist(all_chr[1, ], use.names = FALSE))
  tcodes <- suppressWarnings(
    as.integer(unlist(all_chr[2, ], use.names = FALSE)))
  body   <- as.data.frame(all_chr[-c(1, 2), , drop = FALSE],
                          stringsAsFactors = FALSE)
  colnames(body) <- header

  date_col <- header[1]
  body <- body[!is.na(body[[date_col]]) & nzchar(trimws(body[[date_col]])), ,
               drop = FALSE]

  d <- suppressWarnings(lubridate::mdy(body[[date_col]]))
  if (mean(is.na(d)) > 0.5) d <- suppressWarnings(lubridate::ymd(body[[date_col]]))
  body[[date_col]] <- d
  body <- body[!is.na(d), , drop = FALSE]
  body[[date_col]] <- as.Date(format(body[[date_col]], "%Y-%m-01"))

  series_cols <- setdiff(header, date_col)
  for (j in series_cols) {
    body[[j]] <- suppressWarnings(as.numeric(body[[j]]))
  }

  if (transform) {
    tcode_map <- setNames(tcodes[-1], series_cols)
    for (j in series_cols) {
      body[[j]] <- apply_tcode(body[[j]], tcode_map[[j]])
    }
  }

  out <- tibble::as_tibble(body)
  out <- dplyr::rename(out, date = !!date_col)
  dplyr::arrange(out, date)
}

#' Load FRED-MD, apply tcodes, return cleaned wide panel.
#'
#' Two-stage NA handling (fixes earlier over-aggressive drop):
#'   1. Skip the first 2 rows to wipe tcode-induced NAs from second diffs.
#'   2. Drop columns with residual NAs (= genuine missing data).
#'   3. Drop remaining rows with any NA (should be zero after step 2).
load_fredmd_data <- function(date_start = "1960-01-01",
                             date_end   = Sys.Date(),
                             exclude_groups = 6L,
                             url = NULL,
                             max_lookback = 6L,
                             verbose = TRUE) {
  first_of_month <- function(d) as.Date(format(as.Date(d), "%Y-%m-01"))
  ds <- first_of_month(date_start)
  de <- first_of_month(date_end)

  if (!is.null(url) && file.exists(url)) {
    csv_path <- url
  } else {
    csv_path <- download_fredmd_csv(url = url, max_lookback = max_lookback,
                                    quiet = !verbose)
  }

  raw <- parse_fredmd_csv(csv_path, transform = TRUE)
  raw <- raw |> filter(date >= ds, date <= de)
  if (verbose) message(sprintf("  parsed %d rows x %d cols",
                               nrow(raw), ncol(raw)))

  # optional group filter using fbi::fredmd_description
  desc <- tryCatch({
    if (requireNamespace("fbi", quietly = TRUE)) {
      data("fredmd_description", package = "fbi", envir = environment())
      as_tibble(fredmd_description)
    } else NULL
  }, error = function(e) NULL)

  avail_cols <- setdiff(colnames(raw), "date")
  if (!is.null(desc) && all(c("fred", "group") %in% colnames(desc))) {
    keep_vars <- desc |> filter(!group %in% exclude_groups) |>
                        pull(fred) |> intersect(avail_cols)
    if (verbose) message(sprintf("  keeping %d / %d variables (excluded groups: %s)",
                                 length(keep_vars), length(avail_cols),
                                 paste(exclude_groups, collapse = ", ")))
  } else {
    keep_vars <- avail_cols
    if (verbose) message("  fbi::fredmd_description unavailable; keeping all.")
  }

  if (length(keep_vars) == 0L) stop("No variables survived the filter.")

  wide <- raw |> select(date, all_of(keep_vars)) |> arrange(date)

  # ---- two-stage NA handling ----
  # Stage 1: skip first 2 rows (tcode 3/6 create 2 leading NAs).
  if (nrow(wide) > 2L) wide <- wide[-seq_len(2L), , drop = FALSE]
  # Stage 2: drop columns with residual NAs (genuine missingness).
  has_na <- sapply(wide[, setdiff(colnames(wide), "date"), drop = FALSE],
                   function(v) any(is.na(v)))
  bad <- names(has_na)[has_na]
  if (length(bad)) wide <- wide |> select(-all_of(bad))
  if (verbose) message(sprintf("  after tcode-NA scrub: %d variables remain",
                               ncol(wide) - 1))
  wide <- wide |> drop_na()

  list(
    data    = wide,
    groups  = if (!is.null(desc)) desc |> filter(fred %in% colnames(wide)) else NULL,
    dropped_groups = exclude_groups
  )
}

# =============================================================================
# PCA
# =============================================================================

fit_pca <- function(macro_wide, var_thresh = 0.95) {
  X <- as.matrix(macro_wide[, -1])
  mu <- colMeans(X)
  sd <- apply(X, 2, stats::sd)
  sd[sd == 0] <- 1
  Z <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")

  pr <- stats::prcomp(Z, center = FALSE, scale. = FALSE)
  var_explained <- pr$sdev^2 / sum(pr$sdev^2)
  cum_var <- cumsum(var_explained)
  above <- which(cum_var >= var_thresh)
  k <- if (length(above)) above[1] else length(cum_var)

  scores <- pr$x[, seq_len(k), drop = FALSE]
  rownames(scores) <- as.character(macro_wide$date)

  list(
    scores  = scores,
    loadings = pr$rotation[, seq_len(k), drop = FALSE],
    mean    = mu, sd = sd, k = k,
    var_explained = var_explained,
    cum_var = cum_var,
    dates   = macro_wide$date
  )
}

project_pca <- function(pca_fit, macro_wide_new) {
  X_new <- as.matrix(macro_wide_new[, -1])
  common <- intersect(colnames(X_new), names(pca_fit$mean))
  X_new <- X_new[, common, drop = FALSE]
  mu <- pca_fit$mean[common]
  sd <- pca_fit$sd[common]
  Z <- sweep(sweep(X_new, 2, mu, "-"), 2, sd, "/")
  scores <- Z %*% pca_fit$loadings[common, , drop = FALSE]
  rownames(scores) <- as.character(macro_wide_new$date)
  scores
}

# =============================================================================
# Sector ETFs via Yahoo Finance (quantmod)
# =============================================================================

#' Default universe per Oliveira et al. (2025) Table 2.
default_etf_tickers <- function() {
  c("SPY", "XLB", "XLE", "XLF", "XLI", "XLK", "XLP", "XLU", "XLV", "XLY")
}

#' Download monthly adjusted-close returns for a set of tickers from Yahoo
#' Finance. Returns a tibble(date, ticker1, ticker2, ...) of simple monthly
#' returns with date = first of each month.
load_etf_returns <- function(tickers = default_etf_tickers(),
                             date_start = "1998-12-01",
                             date_end   = Sys.Date(),
                             verbose = TRUE) {
  ensure_pkg("quantmod")
  ensure_pkg("xts")

  if (verbose) message(sprintf("  downloading %d tickers from Yahoo...",
                               length(tickers)))

  fetch_one <- function(tkr) {
    px <- tryCatch(
      quantmod::getSymbols(tkr, src = "yahoo",
                           from = as.Date(date_start),
                           to   = as.Date(date_end),
                           auto.assign = FALSE,
                           warnings = FALSE),
      error = function(e) {
        warning("Failed to fetch ", tkr, ": ", conditionMessage(e))
        NULL
      })
    if (is.null(px)) return(NULL)
    adj <- quantmod::Ad(px)
    # month-end adjusted close, then simple monthly returns
    monthly_last <- xts::apply.monthly(adj, function(xs) as.numeric(xs[nrow(xs)]))
    ret <- as.numeric(monthly_last) / as.numeric(stats::lag(monthly_last, 1)) - 1
    out <- tibble::tibble(
      date = as.Date(format(zoo::index(monthly_last), "%Y-%m-01")),
      !!tkr := ret
    )
    out[!is.na(out[[tkr]]), , drop = FALSE]
  }

  per_tkr <- lapply(tickers, fetch_one)
  names(per_tkr) <- tickers
  missing <- tickers[sapply(per_tkr, is.null)]
  if (length(missing)) {
    warning("Missing tickers: ", paste(missing, collapse = ", "))
  }
  per_tkr <- per_tkr[!sapply(per_tkr, is.null)]

  # full outer join by date
  rets <- Reduce(function(a, b) full_join(a, b, by = "date"), per_tkr)

  # align to first-of-month and sort
  rets <- rets |>
    mutate(date = as.Date(format(date, "%Y-%m-01"))) |>
    arrange(date)

  # drop leading rows before ALL tickers have data
  avail_cols <- setdiff(colnames(rets), "date")
  has_all <- complete.cases(rets[, avail_cols, drop = FALSE])
  first_full <- which(has_all)[1]
  if (is.na(first_full)) stop("No month has all tickers.")
  rets <- rets[first_full:nrow(rets), , drop = FALSE]

  if (verbose) {
    message(sprintf("  ETF returns: %s to %s (%d months x %d assets)",
                    min(rets$date), max(rets$date),
                    nrow(rets), ncol(rets) - 1))
  }

  list(
    returns = rets,
    tickers = avail_cols
  )
}

# =============================================================================
# Alignment helper
# =============================================================================

align_macro_factors <- function(pca_scores_mat, factor_df) {
  macro_dates <- as.Date(rownames(pca_scores_mat))
  fac_dates   <- factor_df$date
  common <- sort(as.Date(intersect(macro_dates, fac_dates),
                         origin = "1970-01-01"))
  list(
    macro   = pca_scores_mat[as.character(common), , drop = FALSE],
    factors = factor_df |> filter(date %in% common) |> arrange(date),
    dates   = common
  )
}
