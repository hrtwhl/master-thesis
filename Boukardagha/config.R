###############################################################################
# config.R — Parameters matching Python replication.py
###############################################################################

FAST_MODE <- FALSE

DATA_PATH   <- "Data/asset_data.csv"
OUTPUT_DIR  <- "Boukardagha/output"
DIAG_DIR    <- "Boukardagha/diagnostics"
CACHE_DIR   <- "Boukardagha/cache"
CACHE_ENABLED <- TRUE

ASSET_NAMES  <- c("stocks", "bonds", "gold", "oil", "usd")
N_ASSETS     <- length(ASSET_NAMES)
ASSET_LABELS <- c("SPX", "BOND", "GOLD", "OIL", "USD")

OOS_START_DATE <- as.Date("2000-01-03")

VOL_WINDOW  <- 60L
MEAN_WINDOW <- 20L

K_MIN       <- 2L
K_MAX       <- 6L
G_TEMPLATES <- 6L
ETA_SMOOTH  <- 0.05
LAMBDA_K    <- 5.0
VAL_SLICE_DAYS <- 252L
ORDER_SELECT_FREQ <- 21L
HMM_FIT_FREQ <- 5L
HMM_INIT_WINDOW <- 1000L
HMM_N_RESTARTS <- 1L
EM_MAX_ITER <- 50L
EM_TOL      <- 1e-3
COV_REG_FACTOR <- 1e-3   # matches hmmlearn default min_covar
SEED        <- 0L

GAMMA_RISK  <- 5.0
TAU_TCOST   <- 0.001
W_MAX       <- 1.0
LEDOIT_WOLF <- TRUE

EW_REBAL_FREQ    <- 5L
BENCH_6040_ALLOC <- c(stocks = 0.60, bonds = 0.40, gold = 0, oil = 0, usd = 0)

SENS_G_GRID      <- c(4L, 6L, 8L)
SENS_ETA_GRID    <- c(0.02, 0.05, 0.10)
SENS_GAMMA_GRID  <- c(2, 5, 10)
SENS_TAU_GRID    <- c(0.0005, 0.001, 0.005)
SENS_LAMBDA_GRID <- c(1.0, 5.0, 10.0)
SENS_WMAX_GRID   <- c(0.40, 0.60, 1.0)

N_CORES <- max(1L, parallel::detectCores() - 1L)

PLOT_WIDTH <- 12; PLOT_HEIGHT <- 7; PLOT_DPI <- 150
REGIME_COLORS <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860")
ASSET_COLORS  <- c("#2171B5", "#FD8D3C", "#238B45", "#CB181D", "#6A51A3")

cat("[config.R] Loaded.\n")
cat(sprintf("  K=[%d,%d] | G=%d | eta=%.2f | lambda_K=%.1f | restarts=%d\n",
            K_MIN, K_MAX, G_TEMPLATES, ETA_SMOOTH, LAMBDA_K, HMM_N_RESTARTS))
