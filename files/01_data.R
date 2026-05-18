# =============================================================================
#  01_data.R
#  Data loading, log returns, feature construction, train/test split
# =============================================================================

#' Load asset price CSV. Expects 'date' + asset columns.
load_asset_data <- function(path, asset_order = c("SPX","BOND","GOLD","OIL","USD")) {
  raw <- read.csv(path, stringsAsFactors = FALSE)
  raw$date <- as.Date(raw$date)
  
  # Map CSV lowercase column names to canonical asset names.
  cmap <- c(stocks = "SPX", bonds = "BOND", gold = "GOLD",
            oil = "OIL", usd = "USD")
  hits <- intersect(names(raw), names(cmap))
  if (length(hits) < length(cmap)) {
    stop("CSV missing expected columns. Found: ",
         paste(names(raw), collapse = ", "))
  }
  names(raw)[match(hits, names(raw))] <- cmap[hits]
  
  raw <- raw[, c("date", asset_order)]
  raw <- raw[order(raw$date), ]
  raw <- raw[complete.cases(raw), ]
  raw
}

#' Daily log returns. Returns data.frame with date + asset log-returns.
compute_log_returns <- function(prices) {
  assets <- setdiff(names(prices), "date")
  P <- as.matrix(prices[, assets])
  R <- diff(log(P))
  data.frame(date = prices$date[-1], R, check.names = FALSE)
}

#' Build features x_t = [r_t, sigma_t, m_t]
#'   sigma_t = trailing std of r over vol_window
#'   m_t     = trailing mean of r over mom_window
#' Right-aligned trailing windows: x_t uses returns up to and including t.
#' Strict causality is enforced at the strategy level (x_{t-1} -> w_t).
build_features <- function(returns_df, vol_window = 60, mom_window = 20) {
  if (!requireNamespace("caTools", quietly = TRUE))
    stop("Install caTools: install.packages('caTools')")
  
  assets <- setdiff(names(returns_df), "date")
  R <- as.matrix(returns_df[, assets])
  Tn <- nrow(R); N <- ncol(R)
  
  vol <- matrix(NA_real_, Tn, N); mom <- matrix(NA_real_, Tn, N)
  for (j in seq_len(N)) {
    vol[, j] <- caTools::runsd(R[, j],  vol_window, endrule = "NA", align = "right")
    mom[, j] <- caTools::runmean(R[, j], mom_window, endrule = "NA", align = "right")
  }
  
  feats <- cbind(R, vol, mom)
  colnames(feats) <- c(
    paste0("ret_", assets),
    paste0("vol_", assets),
    paste0("mom_", assets)
  )
  out <- data.frame(date = returns_df$date, feats, check.names = FALSE)
  out[complete.cases(out), ]
}

#' Split into train/test by date. Returns aligned returns and features.
train_test_split <- function(returns_df, features_df, split_date) {
  split_date <- as.Date(split_date)
  list(
    returns_train  = returns_df[returns_df$date <  split_date, ],
    returns_test   = returns_df[returns_df$date >= split_date, ],
    features_train = features_df[features_df$date <  split_date, ],
    features_test  = features_df[features_df$date >= split_date, ]
  )
}
