"""
config.py
---------
Central configuration for the Wasserstein-HMM regime-aware allocation backtest.

All hyperparameters that govern the experiment live here so that `main.py`
remains a thin orchestration layer.

Settings match Boukardagha (2026) exactly **except for three deliberate
deviations** required to handle the extended 1990–2025 sample. Each is
labelled below with «EXTENSION:» and is justified in `README.md`:

    EXTENSION 1.  max_regimes      :  6  → 8
    EXTENSION 2.  g_max            :  8  → 10
    EXTENSION 3.  f_hmm            :  1  → 5   (HMM refit cadence)

Restoring `max_regimes=6`, `g_max=8`, `f_hmm=1` here recovers the exact
paper specification.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


# --------------------------------------------------------------------- #
# 1. Asset universe
# --------------------------------------------------------------------- #
# Friendly names used everywhere downstream.
ASSET_NAMES: list[str] = ["SPX", "BOND", "GOLD", "OIL", "USD"]
N_ASSETS:    int       = len(ASSET_NAMES)

# CSV column → friendly name. Lets us point at an external file with
# arbitrary lower-case asset labels (as in `asset_data.csv`).
CSV_COLUMN_MAP: dict[str, str] = {
    "stocks": "SPX",
    "bonds":  "BOND",
    "gold":   "GOLD",
    "oil":    "OIL",
    "usd":    "USD",
}


# --------------------------------------------------------------------- #
# 2. Data source
# --------------------------------------------------------------------- #
# Path to the long-history CSV. Resolved relative to this file unless absolute.
CSV_PATH: Path = Path(__file__).resolve().parent / "asset_data.csv"


# --------------------------------------------------------------------- #
# 3. Sample period & train/test split
# --------------------------------------------------------------------- #
# Use the full long history. END_DATE = None means "go to the last row
# available in the source".
START_DATE = "1990-01-02"
END_DATE   = None

# *** OOS start date — see README §"OOS start choice" for the reasoning. ***
# 15 years of training history (1990–2004) precede this; 21 years of OOS follow.
# The window deliberately covers the GFC, COVID, and the 2022 inflation regime.
SPLIT_DATE = "2005-01-03"


# --------------------------------------------------------------------- #
# 4. Feature construction
# --------------------------------------------------------------------- #
VOL_WINDOW = 60   # 60-day rolling volatility per asset
MOM_WINDOW = 20   # 20-day rolling mean return per asset


# --------------------------------------------------------------------- #
# 5. Wasserstein-HMM hyperparameters
# --------------------------------------------------------------------- #
@dataclass(frozen=True)
class HMMConfig:
    # ---- Predictive model-order selection ----
    min_regimes: int = 5
    # EXTENSION 1: raised 6 → 8 to give the K-selector enough headroom to
    # capture the more diverse regime set observed across 1990–2025
    # (90s bull, dot-com bust, GFC, ZIRP, COVID, 2022 inflation, etc.).
    # The PredLL − λ_K·K rule still penalises complexity, so this only
    # widens the candidate set, it does not force higher K.
    max_regimes: int = 8

    f_k:          int = 5    # Re-select K every F_K days (weekly)

    # EXTENSION 3: HMM EM refit cadence. The paper refits daily (f_hmm=1)
    # which is prohibitive over a 35-year expanding window. Refitting only
    # on K-selection days (f_hmm=5) is consistent with the existing K-selection
    # schedule and reduces backtest wall-time by ~5×. Regime probabilities,
    # labels, conditional moments, and templates *still update daily* via
    # `predict_proba` on the extended feature matrix; only the HMM transition/
    # emission *parameters* stay frozen for up to 4 days between refits.
    # Set f_hmm = 1 here to recover the paper's exact daily-refit behaviour.
    f_hmm:        int = 5

    l_val:        int = 252  # Validation window inside history (~1 trading year)
    lam_k:        float = 1.0  # Complexity penalty in PredLL - lam_k * K
    n_iter:       int = 300    # EM iterations for hmmlearn.GaussianHMM
    random_state: int = 42

    # ---- Template tracking (persistent regime identities) ----
    # EXTENSION 2: raised 8 → 10 to accommodate the wider max_regimes and
    # the larger set of distinct macro regimes plausibly emerging over a
    # 21-year OOS window.
    g_max:        int = 10
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


# Performance + speed knobs (do NOT affect methodology)
@dataclass(frozen=True)
class RunConfig:
    verbose:    bool = False  # Print K-state transitions (off by default — they're rare and the 5%-interval progress log already shows current K)
    parallel_k: bool = True   # Joblib-parallel K candidate fitting
    n_jobs:     int  = -1     # -1 = all cores


RUN = RunConfig()


# Trading-day convention used for annualization
TRADING_DAYS = 252


def ensure_dirs() -> None:
    """Create output folders if they don't exist."""
    for d in (OUTPUT_DIR, FIG_DIR, TBL_DIR):
        d.mkdir(parents=True, exist_ok=True)
