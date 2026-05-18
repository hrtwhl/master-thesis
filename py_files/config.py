"""
config.py
---------
Central configuration for the Wasserstein-HMM regime-aware allocation backtest.

All hyperparameters that govern the experiment live here so that `main.py`
remains a thin orchestration layer. Settings match Boukardagha (2026)
exactly unless explicitly noted.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


# --------------------------------------------------------------------- #
# 1. Asset universe
# --------------------------------------------------------------------- #
# Yahoo-Finance ticker mapping. Keys are the friendly names used everywhere
# downstream; values are the Yahoo symbols pulled from the API.
TICKERS: dict[str, str] = {
    "SPX":  "^GSPC",   # S&P 500 Index
    "BOND": "IEF",     # iShares 7-10y Treasury  (broad bond proxy)
    "GOLD": "GLD",     # SPDR Gold Shares
    "OIL":  "USO",     # US Oil Fund
    "USD":  "UUP",     # Invesco USD Index Bullish
}
ASSET_NAMES: list[str] = list(TICKERS.keys())
N_ASSETS: int = len(ASSET_NAMES)


# --------------------------------------------------------------------- #
# 2. Sample period & train/test split
# --------------------------------------------------------------------- #
START_DATE = "2005-01-01"
END_DATE   = None         # None  ⇒  use today (set by data.py)
SPLIT_DATE = "2023-06-03" # First strictly OOS day (matches paper code)


# --------------------------------------------------------------------- #
# 3. Feature construction
# --------------------------------------------------------------------- #
VOL_WINDOW = 60   # 60-day rolling volatility per asset
MOM_WINDOW = 20   # 20-day rolling mean return per asset


# --------------------------------------------------------------------- #
# 4. Wasserstein-HMM hyperparameters
# --------------------------------------------------------------------- #
@dataclass(frozen=True)
class HMMConfig:
    # Predictive model-order selection
    min_regimes: int = 5
    max_regimes: int = 6
    f_k:         int = 5    # Re-select K every F_K days (weekly)
    l_val:       int = 252  # Validation window inside history (~1 trading year)
    lam_k:       float = 1.0  # Complexity penalty in PredLL - lam_k * K
    n_iter:      int = 300    # EM iterations for hmmlearn.GaussianHMM
    random_state: int = 42

    # Template tracking (persistent regime identities)
    g_max:        int = 8     # Max number of templates allowed
    eta_tpl:      float = 0.05  # EMA learning rate for template updates
    spawn_thresh: float = 2.5   # W2 distance threshold to spawn a new template
    min_obs:      int = 60      # Min obs per regime to estimate its moments


HMM = HMMConfig()


# --------------------------------------------------------------------- #
# 5. Portfolio optimization
# --------------------------------------------------------------------- #
@dataclass(frozen=True)
class MVOConfig:
    lam:    float = 5.0      # Risk aversion γ
    tc:     float = 0.0002   # L1 transaction-cost penalty τ
    w_max:  float = 0.6      # Max single-asset weight
    solver: str   = "OSQP"   # CVXPY solver (falls back to SCS on failure)


MVO = MVOConfig()


# --------------------------------------------------------------------- #
# 6. Reporting / output paths
# --------------------------------------------------------------------- #
OUTPUT_DIR = Path(__file__).resolve().parent / "output"
FIG_DIR    = OUTPUT_DIR / "figures"
TBL_DIR    = OUTPUT_DIR / "tables"
CACHE_DIR  = OUTPUT_DIR / "cache"


# Performance + speed knobs (do NOT affect methodology)
@dataclass(frozen=True)
class RunConfig:
    verbose:      bool = True   # Print K each iteration?
    use_cache:    bool = True   # Cache raw Yahoo download to disk
    parallel_k:   bool = True   # Joblib-parallel K candidate fitting
    n_jobs:       int  = -1     # -1 = all cores
    tqdm:         bool = True   # Show progress bar


RUN = RunConfig()


# Trading-day convention used for annualization
TRADING_DAYS = 252


def ensure_dirs() -> None:
    """Create output folders if they don't exist."""
    for d in (OUTPUT_DIR, FIG_DIR, TBL_DIR, CACHE_DIR):
        d.mkdir(parents=True, exist_ok=True)
