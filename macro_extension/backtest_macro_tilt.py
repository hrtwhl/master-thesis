"""
backtest_macro_tilt.py
----------------------
v10.2 backtest: baseline market HMM (unchanged) + macro-tilted MVO.

Architecture
~~~~~~~~~~~~
Step 1. Calibration (once, on pre-OOS data):
    a. Fit macro HMM (BIC-selected M with persistence prior)
    b. Decode macro Viterbi sequence m_t on calibration period
    c. Compute macro_tilt[m] = E[r | m] - E[r] on calibration sample

Step 2. Warmup (templates): unchanged from baseline.
    Builds MU_ret_tpl, SIG_ret_tpl, MU_tpl_feat, SIG_tpl_feat.

Step 3. OOS loop (each day t):
    a. Standard baseline regime step: K-selection, HMM refit, W2 mapping
       -> pG, mu_baseline, Sigma_baseline (as in baseline)
    b. Macro decode: Viterbi on macro_state[0:t] using FROZEN macro HMM
       -> m_t
    c. Tilt application:
         mu_t    = mu_baseline + alpha_mu * macro_mu_tilt[m_t]
         Sigma_t = Sigma_baseline * scale(m_t, alpha_sig)
    d. MVO: solve_mvo(mu_t, Sigma_t, w_prev, ...) -> w_t
    e. Slow EMA template update (unchanged from baseline)

The macro tilts are computed once at calibration and never updated
during OOS. This is strictly causal — the tilt at time t uses ONLY
information available before OOS started.

Compared to v10/v10.1
~~~~~~~~~~~~~~~~~~~~~
v10 and v10.1 modified the market HMM's transition matrix to depend
on the macro regime. This required custom EM with macro-gated A_m
matrices. Over long horizons it produced flat weights because:
    - Template-anchored emissions kept K-components pinned to templates
    - Slow EMA pulled per-template MU_ret toward OOS conditional means
    - These conditional means were close enough that MVO outputs converged
      to a single static allocation

v10.2 leaves the market HMM unchanged. The macro contribution is
purely additive at the MVO step. This guarantees that:
    - Regime detection diversity is preserved (baseline behavior unchanged)
    - The macro layer's value-add is cleanly attributable as the
      MVO performance delta from alpha_mu=0 (= baseline) to alpha_mu>0
    - Long-horizon stability is automatic — tilts don't drift
"""

from __future__ import annotations

import _paths  # noqa: F401

import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from config import HMM, MVO, N_ASSETS, RUN, ASSET_NAMES, WARMUP_START_DATE
from data import build_features

# Baseline machinery (unchanged)
from backtest import (
    BacktestResult, RegimeState,
    solve_mvo,
    _step_regime_fixed_templates,
    _build_fixed_templates,
    _slow_ema_update_templates,
)

# v10.2-specific
from macro_regime_hmm import MacroRegimeHMM
from macro_tilt import calibrate_macro_tilts, apply_macro_tilt
from data_macro import load_macro_state


def run_backtest_macro_tilt(
    returns:        pd.DataFrame,
    returns_train:  pd.DataFrame,
    returns_test:   pd.DataFrame,
    features_all:   Optional[pd.DataFrame] = None,
    macro_state:    Optional[pd.DataFrame] = None,
    verbose:        bool  = RUN.verbose,
    progress_pct:   float = 5.0,
    warmup_start:   Optional[str] = None,
    # Macro layer config
    macro_M_candidates:        tuple[int, ...] = (3, 4, 5),
    macro_transmat_prior_diag: float = 10.0,
    # Tilt config
    alpha_mu:       float = 1.0,
    alpha_sig:      float = 0.0,
    tilt_shrink:    float = 0.5,
    tilt_min_obs:   int   = 50,
) -> tuple[BacktestResult, dict]:
    """Run the macro-tilted backtest.

    Parameters
    ----------
    returns, returns_train, returns_test, features_all, macro_state:
        Same conventions as baseline.
    alpha_mu : float, default 1.0
        Strength of macro mean tilt. 0 = baseline (no tilt). 1 = full
        calibration-time tilt magnitude.
    alpha_sig : float, default 0.0
        Strength of macro volatility-scale adjustment. 0 = leave Sigma
        alone. Recommended to start at 0 (tilt only on the mean).
    tilt_shrink : float in [0, 1], default 0.5
        Bayesian shrinkage of regime-conditional means toward the
        unconditional mean. 0 = raw (high variance); 1 = full shrink
        (no tilt). 0.5 is a sensible default.
    tilt_min_obs : int, default 50
        Minimum calibration days for a regime to get its own tilt.
        Smaller regimes default to zero tilt.

    Returns
    -------
    (BacktestResult, diagnostics dict)
    """
    if features_all is None:
        features_all = build_features(returns)

    if macro_state is None:
        macro_state = load_macro_state(feature_index=features_all.index)
    else:
        macro_state = macro_state.reindex(features_all.index, method="ffill")

    ret_hist_full = pd.concat([returns_train, returns_test]).sort_index()
    ret_hist_full = ret_hist_full[~ret_hist_full.index.duplicated(keep="last")]
    hist_idx_np   = ret_hist_full.index.values

    feat_idx_np     = features_all.index.values
    feat_values_np  = features_all.to_numpy()
    macro_values_np = macro_state.to_numpy()

    test_dates = returns_test.index
    n_total    = len(test_dates)
    split_date = test_dates[0]

    if warmup_start is None:
        warmup_start = WARMUP_START_DATE

    # ===================================================================
    # WARMUP — baseline template construction (unchanged)
    # ===================================================================
    state = RegimeState(K_curr=HMM.min_regimes, K_candidate=HMM.min_regimes)

    if warmup_start is None:
        min_feat_needed = HMM.l_val + 50
        if min_feat_needed < len(features_all):
            auto_start = features_all.index[min_feat_needed]
            warmup_dates = returns_train.index[returns_train.index >= auto_start]
        else:
            warmup_dates = returns_train.index[:0]
    else:
        ws = pd.Timestamp(warmup_start)
        if ws >= split_date:
            warmup_dates = returns_train.index[:0]
        else:
            warmup_dates = returns_train.index[returns_train.index >= ws]

    if len(warmup_dates) > 0:
        state = _build_fixed_templates(
            state=state,
            warmup_dates=warmup_dates,
            feat_idx_np=feat_idx_np,
            feat_values_np=feat_values_np,
            features_index=features_all.index,
            returns_train=returns_train,
            verbose=verbose,
        )
    else:
        raise RuntimeError(
            "v10.2 requires warmup; set WARMUP_START_DATE in config or "
            "pass warmup_start explicitly."
        )

    # ===================================================================
    # MACRO HMM CALIBRATION + TILT COMPUTATION (once, pre-OOS)
    # ===================================================================
    print("\n  Calibrating macro layer HMM...", flush=True)
    cal_cutoff_feat = np.searchsorted(feat_idx_np, np.datetime64(split_date), side="left")
    s_cal = macro_values_np[:cal_cutoff_feat]
    valid = ~np.isnan(s_cal).any(axis=1)
    s_cal_clean = s_cal[valid]
    cal_clean_dates = features_all.index[:cal_cutoff_feat][valid]

    print(f"    Macro calibration sample: {s_cal_clean.shape[0]} days "
          f"(of {s_cal.shape[0]} total)", flush=True)

    macro_hmm = MacroRegimeHMM(
        M_candidates=macro_M_candidates,
        transmat_prior_diag=macro_transmat_prior_diag,
        n_iter=HMM.n_iter,
        random_state=HMM.random_state,
    )
    macro_hmm.fit_initial(s_cal_clean, verbose=verbose)
    M_macro = macro_hmm.M_
    print(f"    Macro HMM fitted: M = {M_macro}, "
          f"BIC scores = {macro_hmm.bic_scores_}", flush=True)

    # Viterbi decode on calibration period
    m_cal_seq = macro_hmm.predict_viterbi(s_cal_clean)
    print(f"    Calibration macro occupancy: "
          f"{pd.Series(m_cal_seq).value_counts(normalize=True).sort_index().round(3).to_dict()}",
          flush=True)

    # Build calibration returns aligned to the clean macro dates
    # (the macro state has NaN warmup at the start, so we drop those
    # days from BOTH the macro sequence and the returns)
    returns_cal_full = returns_train.reindex(cal_clean_dates)
    valid_ret = ~returns_cal_full.isna().any(axis=1)
    R_cal      = returns_cal_full[valid_ret].to_numpy()
    m_cal_used = m_cal_seq[valid_ret.to_numpy()]
    print(f"    Calibration sample for tilts: {R_cal.shape[0]} days", flush=True)

    macro_tilts = calibrate_macro_tilts(
        macro_seq_cal=m_cal_used,
        returns_cal=R_cal,
        M_macro=M_macro,
        min_obs=tilt_min_obs,
        shrink=tilt_shrink,
    )

    print(f"    Calibration unconditional mean (ann %): "
          f"{(macro_tilts['unconditional_mu']*252*100).round(2)}", flush=True)
    print(f"    Per-regime mean tilts (ann %, after shrinkage={tilt_shrink}):", flush=True)
    for m in range(M_macro):
        tilt_pct = macro_tilts['macro_mu_tilt'][m] * 252 * 100
        print(f"      M{m} (n={macro_tilts['macro_counts'][m]:4d}): tilt = {tilt_pct.round(2)}", flush=True)

    # ===================================================================
    # OOS LOOP
    # ===================================================================
    weights_oos:       list = []
    pnl_oos:           list = []
    dates_oos:         list = []
    K_history:         list = []
    tpl_label_history: list = []
    tpl_prob_history:  list = []
    tpl_count_history: list = []
    macro_history:     list = []

    t_start  = time.time()
    last_log = -progress_pct
    w_prev   = np.zeros(N_ASSETS)

    print(f"\nRunning v10.2 backtest "
          f"(macro M={M_macro}, alpha_mu={alpha_mu}, alpha_sig={alpha_sig}, "
          f"OOS days = {n_total})", flush=True)

    for t_i, date in enumerate(test_dates):
        cutoff_feat = np.searchsorted(feat_idx_np, np.datetime64(date), side="left")

        pct_done = (t_i + 1) / n_total * 100.0
        if progress_pct > 0 and (pct_done - last_log >= progress_pct
                                 or t_i == n_total - 1):
            elapsed = time.time() - t_start
            rate = (t_i + 1) / max(elapsed, 1e-6)
            eta_min = (n_total - t_i - 1) / rate / 60.0
            arr = np.asarray(pnl_oos, dtype=float)
            if arr.size > 5 and arr.std(ddof=1) > 0:
                sh = np.sqrt(252) * arr.mean() / arr.std(ddof=1)
                sh_str = f"Sharpe={sh:+.2f}"
            else:
                sh_str = "Sharpe=  n/a"
            m_today_print = macro_history[-1] if macro_history else -1
            print(f"  [{pct_done:5.1f}%]  {t_i+1:>5}/{n_total}  "
                  f"date {pd.Timestamp(date).date()}  "
                  f"elapsed {elapsed/60:5.1f}m  ETA {eta_min:5.1f}m  "
                  f"K={state.K_curr}  m={m_today_print}  {sh_str}",
                  flush=True)
            last_log = pct_done

        if cutoff_feat < 2:
            r_today = returns.loc[date].to_numpy()
            pnl     = float(w_prev @ r_today)
            weights_oos.append(w_prev.copy())
            pnl_oos.append(pnl); dates_oos.append(date)
            K_history.append(np.nan)
            tpl_label_history.append(np.nan)
            tpl_prob_history.append(np.nan)
            tpl_count_history.append(state.G)
            macro_history.append(-1)
            continue

        X = feat_values_np[:cutoff_feat]

        # --- (A) Standard baseline regime step (unchanged) ---
        pK, map_k_to_g, pG, K_active = _step_regime_fixed_templates(
            state=state, t_i=t_i, X=X,
            f_hmm=HMM.f_hmm, verbose=verbose, date_for_print=date,
        )

        if pK is None:
            r_today = returns.loc[date].to_numpy()
            pnl     = float(w_prev @ r_today)
            weights_oos.append(w_prev.copy())
            pnl_oos.append(pnl); dates_oos.append(date)
            K_history.append(np.nan)
            tpl_label_history.append(np.nan)
            tpl_prob_history.append(np.nan)
            tpl_count_history.append(state.G)
            macro_history.append(-1)
            continue

        # --- (B) Macro decode (Viterbi on frozen macro HMM) ---
        s_macro_today = macro_values_np[:cutoff_feat]
        # Forward-fill any NaN
        if np.isnan(s_macro_today).any():
            s_macro_today = pd.DataFrame(s_macro_today).ffill().fillna(0).to_numpy()
        try:
            m_seq = macro_hmm.predict_viterbi(s_macro_today)
            m_today = int(m_seq[-1])
        except Exception:
            m_today = -1

        # --- (C) Compute baseline MVO inputs ---
        mu_baseline    = state.MU_ret_tpl.T @ pG
        Sigma_baseline = np.tensordot(pG, state.SIG_ret_tpl, axes=(0, 0))

        # --- (D) Apply macro tilt ---
        mu_t, Sigma_t = apply_macro_tilt(
            mu_baseline=mu_baseline,
            Sigma_baseline=Sigma_baseline,
            m_today=m_today,
            macro_tilts=macro_tilts,
            alpha_mu=alpha_mu,
            alpha_sig=alpha_sig,
        )

        # --- (E) MVO ---
        w_t = solve_mvo(mu_t, Sigma_t, w_prev,
                        lam=MVO.lam, tc=MVO.tc, w_max=MVO.w_max,
                        solver=MVO.solver)

        r_today = returns.loc[date].to_numpy()
        pnl = float(w_t @ r_today)

        # --- (F) Slow EMA template update (unchanged) ---
        _slow_ema_update_templates(
            state=state, pG=pG, map_k_to_g=map_k_to_g, pK=pK,
            X_today=X[-1] if X.shape[0] > 0 else None,
            cached_hmm=state.hmm_cached,
            r_today=r_today,
            r_prev=state.r_prev_for_ema,
            pG_prev=state.pG_prev_for_ema,
        )
        state.pG_prev_for_ema = pG.copy()
        state.r_prev_for_ema  = r_today.copy()

        weights_oos.append(w_t)
        pnl_oos.append(pnl); dates_oos.append(date)
        K_history.append(K_active)
        g_hat = int(np.argmax(pG)) if np.all(np.isfinite(pG)) else np.nan
        tpl_label_history.append(g_hat)
        tpl_prob_history.append(float(np.max(pG))
                                if np.all(np.isfinite(pG)) else np.nan)
        tpl_count_history.append(state.G)
        macro_history.append(m_today)
        w_prev = w_t

    if state.fail_count > 0:
        print(f"  [note] {state.fail_count} HMM refit(s) fell back to cached "
              f"model.", flush=True)

    pnl = pd.Series(pnl_oos, index=pd.to_datetime(dates_oos)).sort_index()
    W   = pd.DataFrame(weights_oos, index=pnl.index, columns=ASSET_NAMES)
    turnover = 0.5 * W.diff().abs().sum(axis=1)

    result = BacktestResult(
        pnl=pnl,
        weights=W,
        cum_pnl=pnl.cumsum(),
        K_history    =pd.Series(K_history,         index=pnl.index, name="K"),
        tpl_label    =pd.Series(tpl_label_history, index=pnl.index, name="regime"),
        tpl_max_prob =pd.Series(tpl_prob_history,  index=pnl.index, name="max_p"),
        tpl_count    =pd.Series(tpl_count_history, index=pnl.index, name="G"),
        turnover     =turnover.fillna(0.0).rename("turnover"),
    )

    diagnostics = {
        "macro_hmm":       macro_hmm,
        "macro_tilts":     macro_tilts,
        "macro_state":     macro_state,
        "macro_seq_oos":   pd.Series(macro_history, index=pnl.index, name="macro"),
        "M_macro":         M_macro,
        "alpha_mu":        alpha_mu,
        "alpha_sig":       alpha_sig,
        "tilt_shrink":     tilt_shrink,
        "templates": {
            "MU_tpl_feat":  state.MU_tpl_feat,
            "SIG_tpl_feat": state.SIG_tpl_feat,
            "MU_ret_tpl":   state.MU_ret_tpl,
            "SIG_ret_tpl":  state.SIG_ret_tpl,
            "is_valid_tpl": state.is_valid_tpl,
        },
    }
    return result, diagnostics
