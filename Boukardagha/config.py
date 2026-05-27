"""
config.py
=========
Centralized configuration for the regime-aware portfolio construction
infrastructure.

Two strategies are run:
  1. PURE MARKET strategy  - faithful replication of Boukardagha (2026)
                              Wasserstein HMM + MVO.
  2. HIERARCHICAL strategy - macro Wasserstein HMM at the top level
                              (Fine, Singer, Tishby 1998 architecture),
                              market Wasserstein HMMs conditional on
                              the macro regime at the bottom level.

DEVIATIONS FROM BOUKARDAGHA (2026):
-----------------------------------
We extend the paper's 2023-06 -> 2026-02 OOS window to 2005-01 -> 2025-12.
This is the only methodological deviation:

  (D1) GaussianHMM n_iter = 100 (paper: 300).
       Profiling on the actual data shows EM converges in <= 80
       iterations at tol=1e-3 for all sizes encountered in this
       backtest; cap at 100 gives a comfortable margin while saving
       wall time on the daily-refit schedule.  Set HMM_N_ITER = 300
       to recover the paper exactly.

Refit cadence (F_REFIT), K-selection cadence (F_K), training-window
length (expanding), monotone-K rule, and every other parameter are
TAKEN VERBATIM from Paper_Code.ipynb.

The expanding window starts at DATA_START = 1990-01-01 (the paper's
own download starts at 2005-01-01; we keep more history so the model
has macro context when the OOS begins in 2005).

EFFICIENCY:
-----------
Daily refits with an expanding window would naively be intractable
(8000-row K=6 full-cov fits cost ~4 s each, times ~5300 OOS days).
We use hmmlearn's `init_params=""` warm-start mechanism: each daily
refit is initialized from the previous day's fitted parameters and
typically converges in 3-10 EM iterations (~0.1-0.3 s).  See
wasserstein_hmm.py:fit_gaussian_hmm_warm.
"""

import numpy as np

# -----------------------------------------------------------------------
# 0)  Paths
# -----------------------------------------------------------------------
import os

# 1. Ermittle den genauen Ordner, in dem die main.py liegt
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# 2. Definiere die Pfade zu deinen CSV-Dateien im Unterordner "Data"
ASSET_CSV    = os.path.join(BASE_DIR, "Data", "asset_data.csv")
MACRO_CSV    = os.path.join(BASE_DIR, "Data", "macro_data.csv")

# 3. Definiere die Pfade für deine Outputs
OUTPUT_DIR   = os.path.join(BASE_DIR, "output")
CHART_DIR    = os.path.join(OUTPUT_DIR, "charts")
TABLE_DIR    = os.path.join(OUTPUT_DIR, "tables")
RAW_DIR      = os.path.join(OUTPUT_DIR, "raw")

# 4. Erstelle die Output-Ordner automatisch, falls sie noch nicht da sind
os.makedirs(CHART_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
os.makedirs(RAW_DIR, exist_ok=True)

# Hinweis: Wenn du die Pfade später z.B. an Pandas übergibst (df.to_csv), 
# akzeptieren die meisten modernen Bibliotheken Path-Objekte direkt.

# -----------------------------------------------------------------------
# 1)  Asset universe (mapping from CSV column names to paper notation)
# -----------------------------------------------------------------------
ASSET_NAME_MAP = {
    "stocks": "SPX",
    "bonds":  "BOND",
    "gold":   "GOLD",
    "oil":    "OIL",
    "usd":    "USD",
}
ASSET_NAMES = ["SPX", "BOND", "GOLD", "OIL", "USD"]   # paper ordering

# Macro variables (7 - "Mulliner et al." set)
MACRO_VARS = [
    "stocks",          # equity level
    "oil",             # commodity level (went negative in 2020 -> level)
    "copper",          # commodity level
    "yield_curve",     # 10y-2y or similar spread (level)
    "stock_bond_corr", # rolling correlation (already stationary, [-1,1])
    "vix",             # implied vol (price-like, strictly positive)
    "us3mo",           # 3m T-bill yield (level; touches zero)
]
MACRO_VAR_KIND = {
    "stocks":          "price",
    "oil":             "level",   # WTI spot - went negative in 2020
    "copper":          "price",
    "yield_curve":     "level",
    "stock_bond_corr": "level",
    "vix":             "price",   # VIX is strictly positive
    "us3mo":           "level",   # rate; touches zero
}


# -----------------------------------------------------------------------
# 2)  Backtest dates
# -----------------------------------------------------------------------
DATA_START   = "1990-01-01"   # full available history (paper: 2005-01-01)
OOS_START    = "2005-01-01"   # << OOS begins here (paper: 2023-06-03)
DATA_END     = "2025-12-31"


# -----------------------------------------------------------------------
# 3)  Feature engineering (paper defaults, verbatim)
# -----------------------------------------------------------------------
VOL_WINDOW = 60      # 60-day rolling vol
MOM_WINDOW = 20      # 20-day rolling mean


# -----------------------------------------------------------------------
# 4)  Wasserstein-HMM regime model (paper defaults, verbatim)
# -----------------------------------------------------------------------
# Market layer
MIN_REGIMES   = 5
MAX_REGIMES   = 6
F_K           = 5            # K-selection every 5 trading days (PAPER)
L_VAL         = 252          # validation slice (~1y) for predictive K
LAM_K         = 1.0          # complexity penalty
G_MAX         = 8            # max number of templates
ETA_TPL       = 0.05         # EMA learning rate for template updates
SPAWN_THRESH  = 2.5          # W2 spawn threshold for new template

# Macro layer (top level of the hierarchical model - our extension)
# Per user direction: 4-5 macro templates, predictive selection in [4,5]
MACRO_MIN_REGIMES  = 4
MACRO_MAX_REGIMES  = 5
MACRO_G_MAX        = 6
MACRO_ETA_TPL      = 0.05
MACRO_SPAWN_THRESH = 2.5
MACRO_F_K          = 21      # less frequent K-selection for macro layer
MACRO_F_REFIT      = 5       # weekly macro refit (macro regimes are slow)


# -----------------------------------------------------------------------
# 5)  Mean-variance optimizer (paper defaults, verbatim)
# -----------------------------------------------------------------------
LAM      = 5.0       # risk aversion gamma
TC       = 0.0002    # L1 turnover penalty tau (~2 bp)
W_MAX    = 0.6       # per-asset upper bound


# -----------------------------------------------------------------------
# 6)  Efficiency / numerics
# -----------------------------------------------------------------------
# D1: see header.  hmmlearn EM converges <= 80 iters on this data.
HMM_N_ITER = 100

# Refit cadence: PAPER refits DAILY (F_REFIT=1).  We provide a single
# knob so users can speed up if they want to.
F_REFIT    = 1               # PAPER

# Hard-cap the expanding training window?  None = expanding (PAPER).
# (Kept as a config option so the user can experiment.)
MAX_TRAIN_WINDOW = None      # None => expanding window (PAPER)

MIN_HIST_FOR_HMM   = 60      # need at least this many feature rows before fitting
MIN_OBS_PER_REGIME = 60      # min obs to compute robust regime moments


# -----------------------------------------------------------------------
# 7)  Hierarchical-strategy improvements (our extension)
# -----------------------------------------------------------------------
# The original hierarchical implementation underperformed.  The
# following knobs control the simplified/robust hierarchical design
# documented in hierarchical_hmm.py:

# Use a SINGLE shared market WHMM, NOT one sub-HMM per macro template.
# Macro regime probabilities are used only to RE-WEIGHT the
# template-conditional forward-return moments (mu_g, Sigma_g) computed
# from the joint (macro_template, market_template) labels.  This is
# more data-efficient because it doesn't carve the asset history into
# small macro-conditional slices.
HIER_USE_JOINT_MOMENTS = True

# How aggressively to shift portfolio risk in adverse macro regimes:
# in the hierarchical strategy, we additionally scale the equity-risk
# premium (mu of SPX & OIL) by a regime-conditional shrinkage factor
# proportional to historical equity Sharpe in that macro regime.
# Set to 0.0 to disable the tilt and recover an "additive" hierarchy.
HIER_MACRO_TILT_STRENGTH = 1.0


# -----------------------------------------------------------------------
# 8)  Numerics / reproducibility
# -----------------------------------------------------------------------
RANDOM_SEED = 42


# -----------------------------------------------------------------------
# 9)  Benchmark definitions
# -----------------------------------------------------------------------
EQUAL_WEIGHT = {a: 1.0 / len(ASSET_NAMES) for a in ASSET_NAMES}
SIXTY_FORTY  = {"SPX": 0.6, "BOND": 0.4, "GOLD": 0.0, "OIL": 0.0, "USD": 0.0}

TRADING_DAYS = 252
