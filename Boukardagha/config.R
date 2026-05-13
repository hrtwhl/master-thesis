###############################################################################
# config.R — Centralized Parameters for Boukardagha (2026) Replication
# Wasserstein HMM Regime-Aware Portfolio Construction
###############################################################################

# ── Speed mode ──────────────────────────────────────────────────────────────
# Set FAST_MODE <- TRUE for quick testing (~5 min); FALSE for full fidelity
# (~15-30 min depending on hardware).
FAST_MODE <- FALSE

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_PATH   <- "Data/asset_data.csv"
OUTPUT_DIR  <- "Boukardagha/output"
DIAG_DIR    <- "Boukardagha/diagnostics"
CACHE_DIR   <- "Boukardagha/cache"

# ── Caching ──────────────────────────────────────────────────────────────────
# When TRUE, backtest results are saved after computation and loaded on
# subsequent runs if all parameters + data are unchanged.
# Override with --force-rerun CLI flag to ignore cache.
CACHE_ENABLED <- TRUE

# ── Asset universe ───────────────────────────────────────────────────────────
ASSET_NAMES <- c("stocks", "bonds", "gold", "oil", "usd")
N_ASSETS    <- length(ASSET_NAMES)
ASSET_LABELS <- c("SPX", "BOND", "GOLD", "OIL", "USD")

# ── Date boundaries ──────────────────────────────────────────────────────────
# 10 years of burn-in from 1990 → OOS starts 2000-01-03
OOS_START_DATE <- as.Date("2000-01-03")

# ── Feature engineering ──────────────────────────────────────────────────────
VOL_WINDOW  <- 60L   # wσ: rolling standard deviation lookback
MEAN_WINDOW <- 20L   # wm: rolling mean return lookback

# Cap any single-day log return at ±100% (±ln(2.72)).
# Prevents the April 2020 WTI negative price event from producing
# artificial ±750% log returns that contaminate features and portfolios.
LOG_RETURN_CAP <- 1.0

# ── HMM estimation ───────────────────────────────────────────────────────────
K_MIN           <- 2L
K_MAX           <- if (FAST_MODE) 4L else 6L
EM_MAX_ITER     <- if (FAST_MODE) 30L else 75L
EM_TOL          <- 1e-4
COV_REG_FACTOR  <- 1e-5
MAX_HMM_WINDOW  <- if (FAST_MODE) 500L else 1000L  # cap expanding window

# ── Model-order selection ────────────────────────────────────────────────────
ORDER_SELECT_FREQ <- if (FAST_MODE) 63L else 21L   # trading days between K re-selection
VALIDATION_LEN    <- 63L   # |V|: ~3 months for predictive LL scoring
LAMBDA_K          <- 0.5   # complexity penalty per free parameter

# ── Wasserstein template tracking ────────────────────────────────────────────
G_TEMPLATES <- 6L
ETA_SMOOTH  <- 0.10

TEMPLATE_INIT_METHOD <- "first_fit"

# ── MVO optimization ─────────────────────────────────────────────────────────
GAMMA_RISK   <- 5.0
TAU_TCOST    <- 0.002
W_MAX        <- 0.40
OPTIM_METHOD <- "L-BFGS-B"

# ── Rebalancing ──────────────────────────────────────────────────────────────
REFIT_FREQ <- if (FAST_MODE) 10L else 5L  # days between full HMM refits

# ── Benchmark parameters ─────────────────────────────────────────────────────
EW_REBAL_FREQ   <- 5L
BENCH_6040_ALLOC <- c(stocks = 0.60, bonds = 0.40, gold = 0, oil = 0, usd = 0)

# ── Sensitivity analysis grids ───────────────────────────────────────────────
SENS_G_GRID     <- c(3L, 4L, 5L, 6L, 7L, 8L)
SENS_ETA_GRID   <- c(0.01, 0.05, 0.10, 0.20, 0.30)
SENS_GAMMA_GRID <- c(1, 2, 5, 10, 20)
SENS_TAU_GRID   <- c(0.0005, 0.001, 0.002, 0.005, 0.01)

# ── Parallelization ──────────────────────────────────────────────────────────
N_CORES <- max(1L, parallel::detectCores() - 1L)

# ── Plotting ─────────────────────────────────────────────────────────────────
PLOT_WIDTH  <- 12
PLOT_HEIGHT <- 7
PLOT_DPI    <- 150

REGIME_COLORS <- c(
  "#1f77b4",  # Template 1 — blue
  "#ff7f0e",  # Template 2 — orange
  "#2ca02c",  # Template 3 — green
  "#d62728",  # Template 4 — red
  "#9467bd",  # Template 5 — purple
  "#8c564b"   # Template 6 — brown
)

cat("[config.R] Configuration loaded.\n")
cat(sprintf("  Mode: %s | Cache: %s | Assets: %s | OOS from: %s\n",
            ifelse(FAST_MODE, "FAST", "FULL"),
            ifelse(CACHE_ENABLED, "ON", "OFF"),
            paste(ASSET_LABELS, collapse=","), OOS_START_DATE))
cat(sprintf("  HMM: K=[%d,%d] | EM_iter=%d | window=%d | refit=%dd | order_select=%dd\n",
            K_MIN, K_MAX, EM_MAX_ITER, MAX_HMM_WINDOW, REFIT_FREQ, ORDER_SELECT_FREQ))
