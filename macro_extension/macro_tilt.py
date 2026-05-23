"""
macro_tilt.py
-------------
Macro-conditional return tilts and volatility scaling for the MVO
allocation step.

Architecture
~~~~~~~~~~~~
v10.2 separates concerns cleanly:

  1. Market regime detection: unchanged from baseline (constant A,
     k-means init, K-selection via BIC, W2 template-mapping). This
     preserves baseline's natural regime-label diversity over long
     OOS windows.

  2. Macro layer: frozen-after-calibration Gaussian HMM on the 7-dim
     macro state vector. Viterbi decode each day to get m_t.

  3. MVO input adjustment: the macro layer's signal is folded into
     the allocation step as an additive return tilt:

         mu_t^adj    = mu_t      + alpha_mu  * macro_mu_tilt[m_t]
         Sigma_t^adj = Sigma_t   * scale_sig(m_t, alpha_sig)

     where `macro_mu_tilt[m]` is the macro-regime-conditional excess
     return (vs unconditional mean) measured on the **calibration
     sample only** (strictly causal, no lookahead), and `alpha_mu`
     is a tilt-strength hyperparameter in [0, 1].

Why this architecture
~~~~~~~~~~~~~~~~~~~~~
The earlier v10/v10.1 architectures gated the market HMM's transition
matrix on the macro regime (A_t = A_{m_t}). This interacted badly
with the W2 template-mapping under long-horizon EMA template drift:
over 21 years of OOS data, the slow EMA pulls all regimes' MU_ret
toward the OOS unconditional mean, causing MVO outputs to collapse
into a single constant allocation regardless of which regime fires.

v10.2 avoids this by NOT touching the market HMM at all. The macro
contribution is purely additive at the MVO step, so the regime
detection diversity that baseline already produces is preserved. The
macro layer's value-add is measured directly as the MVO performance
delta from `alpha_mu = 0` (baseline) to `alpha_mu > 0`.

Causality
~~~~~~~~~
`macro_mu_tilt[m]` is calibrated ONCE on the pre-OOS sample (returns
< OOS start, paired with macro Viterbi sequence < OOS start). It does
NOT update during OOS. This means at every OOS day, the tilt applied
is identical to the day-1 tilt — no information from after t is used.

The macro Viterbi sequence m_t IS updated daily (using the frozen
calibration macro HMM applied to the macro state up to t). But the
*tilt values themselves* are frozen at calibration.

Outputs of `calibrate_macro_tilts`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    macro_mu_tilt   : (M, N_assets) — per-regime mean return DEVIATION
                                       from unconditional calibration mean
    macro_vol_scale : (M,)           — per-regime volatility multiplier
                                       (sqrt of ratio: regime cov vs uncond)
    macro_counts    : (M,)           — # calibration days per macro regime
"""

from __future__ import annotations

import _paths  # noqa: F401

from typing import Optional

import numpy as np
import pandas as pd

from config import N_ASSETS, ASSET_NAMES


def calibrate_macro_tilts(
    macro_seq_cal:  np.ndarray,    # (T_cal,) — Viterbi macro regimes on calibration
    returns_cal:    np.ndarray,    # (T_cal, N) — asset returns on calibration
    M_macro:        int,
    min_obs:        int = 50,
    shrink:         float = 0.5,
) -> dict:
    """Compute the per-macro-regime mean-return tilt and volatility scale
    from the calibration sample.

    Parameters
    ----------
    macro_seq_cal : (T_cal,) — Viterbi-decoded macro regime labels.
    returns_cal   : (T_cal, N) — daily asset returns.
    M_macro       : int — number of macro regimes.
    min_obs       : int — minimum days for a regime to get its own tilt;
                          if fewer, the regime's tilt is set to zero
                          (i.e., it falls back to unconditional mean).
    shrink        : float in [0, 1] — Bayesian shrinkage of the regime
                          mean toward the unconditional mean. With s=0
                          we use raw regime means (high variance, prone
                          to overfitting in small regimes). With s=1 we
                          use zero tilt (recovers baseline). Default 0.5
                          is a compromise.

    Returns
    -------
    dict with keys:
        macro_mu_tilt     : (M, N) — per-regime return deviation
                                       from unconditional mean
        macro_vol_scale   : (M,)   — per-regime volatility multiplier
        macro_counts      : (M,)   — calibration days per regime
        unconditional_mu  : (N,)   — calibration unconditional mean
        unconditional_vol : (N,)   — calibration unconditional vol
    """
    if returns_cal.shape[0] != macro_seq_cal.shape[0]:
        raise ValueError(
            f"Length mismatch: returns has {returns_cal.shape[0]}, "
            f"macro_seq has {macro_seq_cal.shape[0]}."
        )

    N = returns_cal.shape[1]
    unconditional_mu  = returns_cal.mean(axis=0)
    unconditional_vol = returns_cal.std(axis=0, ddof=1)

    macro_mu_tilt   = np.zeros((M_macro, N))
    macro_vol_scale = np.ones(M_macro)
    macro_counts    = np.zeros(M_macro, dtype=int)

    for m in range(M_macro):
        mask = (macro_seq_cal == m)
        n_m  = int(mask.sum())
        macro_counts[m] = n_m
        if n_m < min_obs:
            # Too few obs — leave at 0 tilt (= unconditional)
            continue

        regime_mu  = returns_cal[mask].mean(axis=0)
        regime_vol = returns_cal[mask].std(axis=0, ddof=1)

        # Bayesian shrinkage of the regime mean toward unconditional
        raw_tilt        = regime_mu - unconditional_mu
        macro_mu_tilt[m] = (1.0 - shrink) * raw_tilt

        # Vol scale: average of per-asset (regime_vol / unconditional_vol)
        # Anchored: shrink also pulls scale toward 1.0
        raw_scale = float(np.mean(regime_vol / np.clip(unconditional_vol, 1e-10, None)))
        macro_vol_scale[m] = (1.0 - shrink) * raw_scale + shrink * 1.0

    return {
        "macro_mu_tilt":     macro_mu_tilt,
        "macro_vol_scale":   macro_vol_scale,
        "macro_counts":      macro_counts,
        "unconditional_mu":  unconditional_mu,
        "unconditional_vol": unconditional_vol,
    }


def apply_macro_tilt(
    mu_baseline:    np.ndarray,    # (N,) — baseline MVO mean
    Sigma_baseline: np.ndarray,    # (N, N) — baseline MVO cov
    m_today:        int,
    macro_tilts:    dict,
    alpha_mu:       float = 1.0,
    alpha_sig:      float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the macro tilt to baseline MVO inputs.

    Parameters
    ----------
    mu_baseline    : (N,) — baseline expected return vector for today.
    Sigma_baseline : (N, N) — baseline covariance for today.
    m_today        : int — macro regime label for today (>= 0).
    macro_tilts    : dict from `calibrate_macro_tilts`.
    alpha_mu       : float — multiplier on the mean tilt. 0 disables
                              the macro contribution entirely (=baseline);
                              1 applies the full calibration-time tilt.
    alpha_sig      : float — multiplier on volatility scaling adjustment.
                              0 leaves Sigma alone; 1 fully applies the
                              calibration-time vol scale.

    Returns
    -------
    mu_adj    : (N,) — tilted mean.
    Sigma_adj : (N, N) — possibly scaled covariance.
    """
    if m_today < 0 or m_today >= macro_tilts["macro_mu_tilt"].shape[0]:
        # Unknown macro regime — no tilt
        return mu_baseline.copy(), Sigma_baseline.copy()

    tilt = macro_tilts["macro_mu_tilt"][m_today]
    mu_adj = mu_baseline + alpha_mu * tilt

    if alpha_sig > 0:
        scale = macro_tilts["macro_vol_scale"][m_today]
        scale_blended = (1.0 - alpha_sig) * 1.0 + alpha_sig * scale
        Sigma_adj = Sigma_baseline * (scale_blended ** 2)
    else:
        Sigma_adj = Sigma_baseline.copy()

    return mu_adj, Sigma_adj


def format_tilts_for_report(
    macro_tilts: dict,
    M_macro:     int,
) -> pd.DataFrame:
    """Format calibration-time tilts into a readable DataFrame
    (annualized %) for thesis tables."""
    rows = []
    mu_tilt = macro_tilts["macro_mu_tilt"] * 252 * 100  # annualized %
    vol_scale = macro_tilts["macro_vol_scale"]
    counts = macro_tilts["macro_counts"]
    uncond_mu = macro_tilts["unconditional_mu"] * 252 * 100

    for m in range(M_macro):
        row = {"macro_regime": m, "n_cal_days": int(counts[m]), "vol_scale": float(vol_scale[m])}
        for i, a in enumerate(ASSET_NAMES):
            row[f"{a}_tilt"]    = float(mu_tilt[m, i])
            row[f"{a}_implied"] = float(uncond_mu[i] + mu_tilt[m, i])
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------- #
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    T, N, M = 1000, 5, 3
    # Synthetic: 3 macro regimes with different means
    macro = rng.integers(0, M, size=T)
    R = np.zeros((T, N))
    base_mu = [
        np.array([0.0005, 0.0002, 0.0001, 0.0003, 0.0000]),   # bull
        np.array([-0.0003, 0.0004, 0.0002, -0.0001, 0.0001]), # bear
        np.array([0.0001, 0.0001, 0.0001, 0.0001, 0.0001]),   # neutral
    ]
    for t in range(T):
        R[t] = rng.multivariate_normal(base_mu[macro[t]], 0.0001 * np.eye(N))

    out = calibrate_macro_tilts(macro, R, M_macro=M)
    print("Synthetic test:")
    print(f"Calibration unconditional mean (annualized %): {out['unconditional_mu']*252*100}")
    print()
    print("Per-regime tilts (annualized %):")
    for m in range(M):
        print(f"  M{m}: tilt={out['macro_mu_tilt'][m]*252*100}, vol_scale={out['macro_vol_scale'][m]:.3f}, n={out['macro_counts'][m]}")
    print()
    print("Pretty-format DataFrame:")
    print(format_tilts_for_report(out, M).round(2))
