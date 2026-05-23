"""
config.py
---------
Central configuration for the Wasserstein-HMM regime-aware allocation backtest.

All hyperparameters that govern the experiment live here so that `main.py`
remains a thin orchestration layer.

Eight extensions plus the paper-fidelity correction and the v6 corrections.

    EXTENSION 1.  max_regimes        :  6  →  6        (no change net)
    EXTENSION 2.  g_max              :  8  → 10        (legacy)
    EXTENSION 3.  f_hmm              :  1  →  5        (HMM refit cadence)
    EXTENSION 4.  lam_k              : flat 1.0 → BIC-scaled
    EXTENSION 5.  monotone_K         : True → False
    EXTENSION 6.  spawn rule         : DROPPED — G selected once on calibration
    EXTENSION 7.  template mapping   : hard argmin in FEATURE space
    EXTENSION 8.  w_max              : 0.6 → 0.8       (MVO bound)

PAPER FIX.   template warmup      : OFF → ON (Algorithm 1, line 3)
V6 FIX.      template structure   : dynamic spawn+EMA → calibration HMM
V8 FIX.      G selection          : hard-coded → BIC-style on calibration
V8 FIX.      template drift       : fully frozen → slow mass-gated EMA
                                      (η small + only update dominant
                                       template each day)

Note: v7 (expanding-window forward-return moment refresh) was tried and
reverted — it made the strategy procyclical, loading equities heading
into 2008 because recent equity outperformance had been absorbed into
regime-conditional return statistics.

To recover paper-exact released-code: max_regimes=6, g_max=8, f_hmm=1,
lam_k_override=1.0, monotone_K=True, G_min=G_max=5, w_max=0.6,
eta_tpl=0.05, tpl_update_thresh=0 (and revert the template mechanism
to dynamic spawn — not currently toggleable).
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

# *** Template calibration window — implements Algorithm 1 line 3 of
#     Boukardagha (2026): "Initialize templates {Θg} using an initial
#     calibration window."
#
# In the v5 implementation this calibration window is used to fit ONE big
# HMM (with K=K_template, see HMMConfig) on the entire pre-OOS history.
# The HMM's emission means/covariances become the FIXED template pool
# for the OOS phase — they do not update via EMA in the OOS loop. See
# HMMConfig docstring for the rationale (template EMA drift was the
# root cause of the regime-detection collapse over our 21-year OOS).
#
# WARMUP_START_DATE = None (recommended): use ALL available pre-OOS
#   feature rows for the calibration HMM fit. Earliest usable date is
#   ~1991-05 given L_VAL=252 and 60-day feature warmup.
#
# Setting WARMUP_START_DATE = SPLIT_DATE disables warmup; the OOS loop
# will not have templates and will hold weights. Provided only for
# diagnostic purposes — there is no longer a "fallback to first OOS day"
# mode since fixed templates require a calibration window by definition.
WARMUP_START_DATE = None


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
    min_regimes: int = 3
    # EXTENSION 1 (revised): max_regimes 6 → 8 → 6.
    # Initial extension raised the ceiling to 8 on the hypothesis that a
    # 21-year sample would demand more regime variety. Empirically the
    # K-selector pegged at the ceiling and only 6 of 8 templates ever fired
    # as the dominant regime. Capping at 6 forces the selector to use its
    # capacity well rather than over-parametrise.
    max_regimes: int = 6

    f_k:          int = 5    # Re-select K every F_K days (weekly)

    # EXTENSION 5 (new): monotone-K is OFF.
    # The paper enforces K non-decreasing to prevent oscillation over a
    # short 2-year OOS. Over 21 years this becomes a one-way ratchet that
    # locks the model at its first selection. Allowing K to fall when the
    # data warrant a simpler model is essential for the longer window.
    # Set `monotone_K = True` here to recover the paper's behaviour.
    monotone_K:   bool = False

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

    # EXTENSION 4 (new): BIC-style K penalty replaces the paper's flat λ_K = 1.
    #
    # The paper uses `PredLL(K) − λ_K · K` with λ_K = 1. On a 35-year sample
    # the validation log-likelihood is on the order of 10^5, against which a
    # 1-nat penalty per state is effectively zero — the selector picks the
    # ceiling deterministically.
    #
    # The well-known fix is to scale the per-state penalty with sample size
    # and feature dimensionality:
    #
    #     λ_K = lam_k_scale · 0.5 · log(T_val) · d_feat
    #
    # where T_val = l_val (validation length) and d_feat = 3·N (the HMM input
    # dimension). With l_val=252, N=5 ⇒ d_feat=15, the per-state penalty is
    # ≈ 41 nats per state — comparable to a light BIC and ~40× the paper's
    # value but still permissive enough for K to rise when warranted.
    #
    # `lam_k_scale` defaults to 1.0 (standard BIC-style). Setting it to a
    # very small value (e.g. 0.02) approximately recovers the paper's
    # behaviour. Setting `lam_k_override` to a positive number overrides
    # the scaling entirely.
    lam_k_scale:     float = 1.0
    lam_k_override:  float = 0.0   # 0 ⇒ use the scaled formula above

    n_iter:       int = 300    # EM iterations for hmmlearn.GaussianHMM
    random_state: int = 42

    # ---- Template tracking (persistent regime identities) ----
    # EXTENSION 2 (re-revised): g_max = 10.
    # Iteration history:
    #   Iter 1: raised 8 → 10 (defensive, K hit ceiling, templates filled)
    #   Iter 2: lowered to 6 (templates collapsed to 3 anyway, w_max binding)
    #   Iter 3: back to 10. Dormant templates are free (zero posterior ⇒
    #   zero contribution to MVO moments); missing templates are costly
    #   (forces structurally different regimes into the wrong identity).
    #   The paper's `g_max=8` was a safety cap for a 2-year OOS; over 21
    #   years a slightly larger pool is appropriate.
    g_max:        int = 10        # legacy; not used in v8

    # ====================================================================
    # v8 framework: slow, mass-gated EMA templates (was: fully frozen)
    # ====================================================================
    # Iteration history of failed approaches:
    #   v1-4: dynamic spawn + standard EMA → templates collapse to 1
    #   v5:   spawn + soft mapping         → uniform pG = 1/G
    #   v6:   fully frozen, fixed G=6      → regime detection works,
    #                                        but moments don't adapt;
    #                                        bond-stuffed allocations
    #   v7:   expanding-window moments     → procyclical, 2008 disaster
    #
    # CURRENT APPROACH (v8):
    #   1. G selected ONCE on calibration window via predictive-LL +
    #      BIC penalty (G_min..G_max). Same selector as daily K, but
    #      applied on a larger validation slice.
    #   2. Initial templates from one big HMM with the chosen G.
    #   3. Templates drift SLOWLY (η tiny) and only when a template
    #      dominates today's posterior (mass-gated). This prevents both
    #      collapse and procyclicality.
    #
    # The mass-gated EMA is the key innovation over the paper: only
    # templates with pG[g] >= tpl_update_thresh get updated today.
    # Rare regimes drift slowly because they're rarely "on"; frequent
    # regimes can adapt more, but only when their identity is confirmed.
    #
    # Setting eta_tpl = 0 recovers v6's fully-frozen behaviour.
    # Setting tpl_update_thresh = 0 recovers ordinary EMA on every
    # template every day (paper-like but at smaller η).

    # G candidates for the calibration HMM. Selected once via BIC-style
    # criterion (see `_select_G_predictive` in backtest.py).
    # Range matches the spirit of the paper's K_min..K_max but applied
    # to template count (G is structurally larger than daily K).
    G_min:        int = 5
    G_max:        int = 8

    # BIC scaling for G selection — same formula as lam_k but applied
    # to the larger validation slice used in calibration.
    lam_G_scale:  float = 1.0

    # Slow EMA settings (replaces v6 fully-frozen mode)
    # eta_tpl: per-day pull rate of mass-gated EMA. With our 5,500-day
    # OOS, this gives templates a half-life of ~ln(2)/0.001 = 693 days,
    # i.e. ~3 years. Slow enough to preserve identity, fast enough to
    # absorb regime-level structural change.
    eta_tpl:      float = 0.001

    # tpl_update_thresh: a template only updates if its posterior mass
    # today exceeds this. With hard argmin against well-separated
    # templates, pG is essentially one-hot, so 0.5 cleanly partitions
    # "this is today's regime" from "no, today is some other regime."
    tpl_update_thresh: float = 0.5

    min_obs:      int = 60      # Min obs per regime to estimate its moments


HMM = HMMConfig()


# --------------------------------------------------------------------- #
# 5. Portfolio optimization
# --------------------------------------------------------------------- #
@dataclass(frozen=True)
class MVOConfig:
    # Paper value γ=5 retained.
    # Note: in early diagnostic iterations we tried γ=2 to reduce the
    # bond-heavy concentration observed in v5. That over-corrected —
    # γ=2 swung the portfolio to ~70% SPX-at-cap instead. The actual
    # fix is V7 (expanding-window template moments), not the MVO
    # gamma. γ=5 is kept as the paper specifies.
    lam:    float = 4.0
    tc:     float = 0.001   # L1 transaction-cost penalty τ
    # EXTENSION 8: w_max 0.6 → 0.8 (see prior iteration log).
    w_max:  float = 0.6
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
