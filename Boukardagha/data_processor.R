###############################################################################
# data_processor.R — Data Loading, Cleaning, and Feature Engineering
# Boukardagha (2026) Replication
###############################################################################

# ── Load prices ──────────────────────────────────────────────────────────────
load_prices <- function(path = DATA_PATH) {
  raw <- read.csv(path, stringsAsFactors = FALSE)
  raw$date <- as.Date(raw$date)
  
  # Ensure columns match expected assets
  stopifnot(all(ASSET_NAMES %in% names(raw)))
  
  # Sort by date, remove duplicates
  raw <- raw[order(raw$date), ]
  raw <- raw[!duplicated(raw$date), ]
  
  # Convert to numeric matrix for speed
  prices <- as.matrix(raw[, ASSET_NAMES])
  storage.mode(prices) <- "double"
  
  dates <- raw$date
  
  # Handle negative prices (e.g., WTI oil April 2020)
  # Floor at 0.01 to avoid log(negative). This affects ~1 day in 35 years
  # and is standard practice for price-index-based backtests.
  neg_count <- sum(prices <= 0, na.rm = TRUE)
  if (neg_count > 0) {
    cat(sprintf("  [WARN] %d negative/zero prices floored to 0.01\n", neg_count))
    prices[prices <= 0] <- 0.01
  }
  
  # Quality checks
  na_pct <- colMeans(is.na(prices)) * 100
  if (any(na_pct > 5)) {
    warning(sprintf("High NA rate: %s", paste(names(which(na_pct > 5)),
                                                collapse = ", ")))
  }
  
  # Forward-fill NAs (standard for daily price data)
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
  cat(sprintf("  NA forward-filled. Max NA%% per col: %.1f%%\n", max(na_pct)))
  
  list(dates = dates, prices = prices)
}

# ── Compute log returns ──────────────────────────────────────────────────────
compute_log_returns <- function(prices, cap = LOG_RETURN_CAP) {
  # r_t = log(P_t) - log(P_{t-1}), clamped to [-cap, cap]
  n <- nrow(prices)
  log_p <- log(prices)
  returns <- log_p[-1, , drop = FALSE] - log_p[-n, , drop = FALSE]
  
  # Clamp extreme returns (e.g. WTI April 2020: $18 → -$37 → floor $0.01
  # produces ln(0.01/18) = -7.5, which is artificial)
  n_capped <- sum(abs(returns) > cap, na.rm = TRUE)
  if (n_capped > 0) {
    cat(sprintf("  [WARN] %d daily log returns capped at +/-%.1f\n", n_capped, cap))
    returns <- pmin(pmax(returns, -cap), cap)
  }
  
  returns
}

# ── Rolling statistics (strictly causal: uses data up to t-1) ────────────────
rolling_sd <- function(x, window) {
  # Returns rolling SD using previous `window` observations
  # Output[t] = sd(x[(t-window):(t-1)])
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in (window + 1):n) {
    out[i] <- sd(x[(i - window):(i - 1)])
  }
  out
}

rolling_mean <- function(x, window) {
  # Output[t] = mean(x[(t-window):(t-1)])
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in (window + 1):n) {
    out[i] <- mean(x[(i - window):(i - 1)])
  }
  out
}

# Vectorized versions for speed
rolling_sd_fast <- function(x, window) {
  n <- length(x)
  out <- rep(NA_real_, n)
  if (n <= window) return(out)
  
  # Cumulative sums for online variance
  cs  <- cumsum(x)
  cs2 <- cumsum(x^2)
  
  for (i in (window + 1):n) {
    s  <- cs[i - 1] - ifelse(i - window - 1 > 0, cs[i - window - 1], 0)
    s2 <- cs2[i - 1] - ifelse(i - window - 1 > 0, cs2[i - window - 1], 0)
    mu <- s / window
    var_est <- s2 / window - mu^2
    # Bessel correction
    out[i] <- sqrt(max(0, var_est * window / (window - 1)))
  }
  out
}

rolling_mean_fast <- function(x, window) {
  n <- length(x)
  out <- rep(NA_real_, n)
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
                           vol_window  = VOL_WINDOW,
                           mean_window = MEAN_WINDOW) {
  
  log_ret <- compute_log_returns(prices)
  # Align dates: returns start from row 2
  ret_dates <- dates[-1]
  
  n <- nrow(log_ret)
  p <- ncol(log_ret)
  
  # Rolling vol and mean for each asset (strictly causal)
  vol_mat  <- matrix(NA_real_, n, p)
  mean_mat <- matrix(NA_real_, n, p)
  
  for (j in seq_len(p)) {
    vol_mat[, j]  <- rolling_sd_fast(log_ret[, j], vol_window)
    mean_mat[, j] <- rolling_mean_fast(log_ret[, j], mean_window)
  }
  
  # Feature vector: x_t = [r_t; sigma_t; m_t] ∈ R^{3N}
  # Strict causality: r_t uses returns up to t-1, but the return itself is
  # at time t. Per the paper, features at time t use info up to t-1.
  # So the "return" component is the lagged return: r_{t-1}.
  # We shift: feature row i corresponds to day i, using return from day i-1.
  
  # Lag returns by 1 for strict causality
  lagged_ret <- rbind(rep(NA_real_, p), log_ret[-n, , drop = FALSE])
  
  features <- cbind(lagged_ret, vol_mat, mean_mat)
  colnames(features) <- c(
    paste0("ret_", ASSET_NAMES),
    paste0("vol_", ASSET_NAMES),
    paste0("mean_", ASSET_NAMES)
  )
  
  # Identify rows with complete features
  valid <- complete.cases(features)
  first_valid <- which(valid)[1]
  
  cat(sprintf("[data_processor] Features built: %d x %d (3N=%d)\n",
              n, ncol(features), 3 * p))
  cat(sprintf("  First valid feature row: %d (%s)\n",
              first_valid, ret_dates[first_valid]))
  
  list(
    dates       = ret_dates,
    returns     = log_ret,      # actual daily log returns (for portfolio P&L)
    features    = features,     # 3N-dimensional feature vectors
    valid       = valid,        # logical: rows with complete features
    first_valid = first_valid
  )
}

cat("[data_processor.R] Functions loaded.\n")
