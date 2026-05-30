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
# 4R)  Robustness analyses
# -----------------------------------------------------------------------
# R1 - Loosened-K macro sweep (robustness_ksweep.py).
#      Tests whether the 2-effective-regime collapse is data-driven or
#      imposed by the tight, penalised default macro config.  We re-run
#      ONLY the macro Wasserstein-HMM (no MVO) under each config below
#      and count how many DURABLE templates emerge (templates whose
#      dominant-share over the OOS exceeds R1_DURABLE_SHARE).
R1_KSWEEP_CONFIGS = [
    # label                  Kmin Kmax Gmax lamK spawn monotone
    ("default_4_5",            4,   5,   6,  1.0,  2.5,  True),
    ("loose_4_8",              4,   8,  10,  1.0,  2.5,  True),
    ("loose_4_8_nopenalty",    4,   8,  10,  0.0,  2.5,  True),
    ("loose_4_8_lowspawn",     4,   8,  10,  0.0,  1.0,  True),
    ("free_3_8_nonmono",       3,   8,  10,  0.0,  1.0,  False),
]
R1_DURABLE_SHARE = 0.02   # a template is "durable" if dominant on >=2% of OOS days

# R2 - No-overlap macro panel (robustness_nooverlap.py).
#      'stocks' (=SPX) and 'oil' are also tradeable assets; the market
#      WHMM already sees them.  Including them in the macro panel makes
#      the macro and market layers informationally dependent.  R2 re-runs
#      Hierarchical B and C with a macro panel that EXCLUDES the
#      overlapping variables, leaving only genuinely exogenous
#      macro-financial series.
R2_MACRO_VARS = ["copper", "yield_curve", "stock_bond_corr", "vix", "us3mo"]


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
# 7)  Hierarchical strategy B  (macro x market joint mixture, no tilt)
# -----------------------------------------------------------------------
# This is the original hierarchical extension WITHOUT any expected-return
# tilt (the earlier 'off'/'asymmetric' modes, which were numerically
# identical because the macro posterior is one-hot ~99% of the time and
# the asymmetric tilt self-cancels).
#
# Mechanism (see methodology.md SS4-B):
#   - Macro WHMM on 21-d macro features, templates tracked in MACRO
#     FEATURE space.
#   - Market WHMM (shared) on asset features.
#   - Joint conditional moments on (g_macro, h_market) cells, with
#     market-only fallback for sparse cells.
#   - Composite probability p(t,g,h) = p_macro(t,g) * p_market(t,h).
#   - mu_t, Sigma_t = sum_{g,h} p(t,g,h) (mu_{g,h}, Sigma_{g,h}).
# No expected-return tilt is applied: mu_t feeds MVO unchanged.
#
# (No tunable parameters here beyond the macro-layer settings in
#  block 4 above.)


# -----------------------------------------------------------------------
# 7C)  Hierarchical strategy C  (macro = risk modulator, our Fix-C)
# -----------------------------------------------------------------------
# Hierarchical C addresses the three problems diagnosed in
# analysis_results.md by changing HOW the macro layer is used, while
# KEEPING the full 21-feature Mulliner macro set (the 2-effective-regime
# degeneracy is itself a result we report).
#
# Three design changes vs. Hierarchical B:
#
#   (C1) Macro templates are tracked in ASSET-OUTCOME space.
#        Instead of tracking macro templates by what macro FEATURES look
#        like, we track them by how ASSETS behave (forward-return mean
#        and covariance) in each macro regime.  This makes the macro
#        regimes maximally discriminative for allocation by construction.
#
#   (C2) No joint-cell fragmentation; macro modulates RISK, not direction.
#        The market layer keeps FULL control of expected returns mu_t
#        (exactly the pure-market mu_t).  The macro layer modulates only
#        the RISK side:
#          - the effective risk aversion gamma_t, and
#          - a multiplicative scale on Sigma_t,
#        both driven by a macro "stress" score.  This preserves the
#        clean directional signal of the market layer (which was being
#        diluted by joint-cell averaging in Hierarchical B).
#
#   (C3) Softened (tempered) macro posterior.
#        The macro posterior is one-hot ~99% of the time, carrying no
#        usable uncertainty.  We temper it with exponent
#        HIER_C_MACRO_TEMPERATURE and blend with a uniform prior of
#        weight HIER_C_PRIOR_BLEND, so the optimizer responds to graded
#        macro uncertainty rather than a hard switch.

# Macro stress -> risk multiplier.  gamma_t = gamma * (1 + kappa * stress_t)
# and Sigma_t *= (1 + sigma_scale * stress_t), where stress_t in [0, 1]
# is the tempered-posterior-weighted, cross-sectionally normalized
# historical turbulence (annualized vol of the portfolio-relevant
# forward returns) of the active macro regime.
HIER_C_KAPPA_GAMMA   = 4.0    # risk-aversion sensitivity to macro stress
HIER_C_SIGMA_SCALE   = 1.0    # covariance inflation sensitivity to macro stress
HIER_C_MACRO_TEMPERATURE = 4.0  # >1 softens (flattens) the macro posterior
HIER_C_PRIOR_BLEND   = 0.10   # blend weight on a uniform macro prior in [0,1]
# Stress definition: 'vol' uses regime forward-return vol; 'drawdown'
# uses regime within-spell max drawdown; 'sharpe' uses negative Sharpe.
HIER_C_STRESS_METRIC = "vol"


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
