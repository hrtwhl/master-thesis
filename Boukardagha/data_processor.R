###############################################################################
# data_processor.R — Data Loading, Cleaning, and Feature Engineering
# Boukardagha (2026) Replication
#
# Paper-faithful: no return clipping, no price flooring.
# Negative/zero prices are treated as missing and forward-filled, since the
# paper uses ETF prices (USO, GLD, etc.) which cannot go negative.
###############################################################################

load_prices <- function(path = DATA_PATH) {
  raw <- read.csv(path, stringsAsFactors = FALSE)
  raw$date <- as.Date(raw$date)
  stopifnot(all(ASSET_NAMES %in% names(raw)))
  
  raw <- raw[order(raw$date), ]
  raw <- raw[!duplicated(raw$date), ]
  
  prices <- as.matrix(raw[, ASSET_NAMES])
  storage.mode(prices) <- "double"
  dates <- raw$date
  
  # Treat negative/zero prices as missing (paper uses ETFs which can't go
  # negative; our raw WTI data has negative prices on 2020-04-20).
  # Forward-fill replaces them with the last valid price → return = 0 on
  # those days, which is the minimal-intervention approach.
  n_invalid <- sum(prices <= 0, na.rm = TRUE)
  if (n_invalid > 0) {
    cat(sprintf("  [INFO] %d negative/zero prices set to NA and forward-filled\n", n_invalid))
    prices[prices <= 0] <- NA
  }
  
  na_pct <- colMeans(is.na(prices)) * 100
  if (any(na_pct > 5)) {
    warning(sprintf("High NA rate: %s", paste(names(which(na_pct > 5)), collapse = ", ")))
  }
  
  # Forward-fill NAs
  for (j in seq_len(ncol(prices))) {
    v <- prices[, j]
    idx <- which(!is.na(v))
    if (length(idx) > 0 && length(idx) < length(v)) {
      prices[, j] <- approx(idx, v[idx], xout = seq_along(v),
                             method = "constant", rule = 2)$y
    }
  }
  
  cat(sprintf("[data_processor] Loaded %d rows, %s to %s\n",
              nrow(prices), dates[1], dates[length(dates)]))
  
  list(dates = dates, prices = prices)
}

# ── Compute log returns (no clipping — paper does not mention it) ────────────
compute_log_returns <- function(prices) {
  n <- nrow(prices)
  log_p <- log(prices)
  log_p[-1, , drop = FALSE] - log_p[-n, , drop = FALSE]
}

# ── Rolling statistics (strictly causal) ─────────────────────────────────────
rolling_sd_fast <- function(x, window) {
  n <- length(x); out <- rep(NA_real_, n)
  if (n <= window) return(out)
  cs <- cumsum(x); cs2 <- cumsum(x^2)
  for (i in (window + 1):n) {
    s  <- cs[i - 1] - ifelse(i - window - 1 > 0, cs[i - window - 1], 0)
    s2 <- cs2[i - 1] - ifelse(i - window - 1 > 0, cs2[i - window - 1], 0)
    mu <- s / window; var_est <- s2 / window - mu^2
    out[i] <- sqrt(max(0, var_est * window / (window - 1)))
  }
  out
}

rolling_mean_fast <- function(x, window) {
  n <- length(x); out <- rep(NA_real_, n)
  if (n <= window) return(out)
  cs <- cumsum(x)
  for (i in (window + 1):n) {
    s <- cs[i - 1] - ifelse(i - window - 1 > 0, cs[i - window - 1], 0)
    out[i] <- s / window
  }
  out
}

# ── Build full feature matrix ────────────────────────────────────────────────
build_features <- function(dates, prices,
                           vol_window = VOL_WINDOW, mean_window = MEAN_WINDOW) {
  log_ret <- compute_log_returns(prices)
  ret_dates <- dates[-1]
  n <- nrow(log_ret); p <- ncol(log_ret)
  
  vol_mat  <- matrix(NA_real_, n, p)
  mean_mat <- matrix(NA_real_, n, p)
  for (j in seq_len(p)) {
    vol_mat[, j]  <- rolling_sd_fast(log_ret[, j], vol_window)
    mean_mat[, j] <- rolling_mean_fast(log_ret[, j], mean_window)
  }
  
  # Feature vector: x_t = [r_{t-1}; sigma_t; m_t] (lagged return for causality)
  lagged_ret <- rbind(rep(NA_real_, p), log_ret[-n, , drop = FALSE])
  features <- cbind(lagged_ret, vol_mat, mean_mat)
  colnames(features) <- c(paste0("ret_", ASSET_NAMES),
                           paste0("vol_", ASSET_NAMES),
                           paste0("mean_", ASSET_NAMES))
  
  valid <- complete.cases(features)
  first_valid <- which(valid)[1]
  
  cat(sprintf("[data_processor] Features: %d x %d | first valid: %d (%s)\n",
              n, ncol(features), first_valid, ret_dates[first_valid]))
  
  list(dates = ret_dates, returns = log_ret, features = features,
       valid = valid, first_valid = first_valid)
}

cat("[data_processor.R] Functions loaded.\n")
