"""
backtest_hierarchical.py
------------------------
Backtest harness for the hierarchical (two-layer HMM) extension.

Architecture
~~~~~~~~~~~~
Macro layer (slow, calibrated once)
    `MacroRegimeHMM`   — Gaussian HMM on the 7-dim macro state s_t
                          with BIC-selected M and persistence prior.
                          Calibrated on data < OOS start, then frozen.
                          Viterbi-decoded each day to produce m_t.

Market layer (fast, refit weekly like baseline)
    `HierarchicalMarketHMM` — custom HMM on the 15-dim market feature
                              space X_t with macro-regime-gated
                              transitions A_t = A_{m_t}.

Templates and MVO
    Identical to baseline:
        * fixed G templates from calibration
        * W2 mapping of market HMM components to templates
        * slow EMA on template forward-return moments
        * MVO with same gamma, tc, w_max

Strict causality
~~~~~~~~~~~~~~~~
1. Macro HMM is calibrated on macro state up to (and including) the
   day before OOS start.
2. At OOS day t, the macro Viterbi sequence m_{1:t} is decoded using
   the frozen macro HMM applied to the macro state up to t.
3. The market HMM is fit on X_{1:t} and m_{1:t}, with refit every
   f_hmm market days (default 5, same as baseline).

No lookahead at any point.

Outputs
~~~~~~~
Same `BacktestResult` dataclass as baseline + a diagnostics dict
containing the macro HMM, the cached market HMM, and the macro
regime sequence — all needed for the comparison module.
"""

from __future__ import annotations

import _paths  # noqa: F401  (injects v9 + baseline into sys.path)

import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from config import HMM, MVO, N_ASSETS, RUN, ASSET_NAMES, WARMUP_START_DATE
from data import build_features

# Reuse baseline machinery for the market layer
from backtest import (
    BacktestResult, RegimeState,
    solve_mvo,
    _build_fixed_templates,
    _slow_ema_update_templates,
)
from wasserstein_hmm import (
    aggregate_template_posteriors_hard,
    map_components_to_fixed_templates,
    select_K_predictive,
)

# v10-specific pieces
from macro_regime_hmm import MacroRegimeHMM
from market_hmm_hierarchical import (
    HierarchicalMarketHMM,
    fit_hierarchical_market_hmm,
)
from data_macro import load_macro_state


# --------------------------------------------------------------------- #
# 1. State container — extends baseline RegimeState
# --------------------------------------------------------------------- #
@dataclass
class HierarchicalRegimeState(RegimeState):
    """State for the hierarchical backtest. Adds macro HMM, the
    decoded macro regime sequence, and the cached market HMM.
    """
    macro_hmm:      Optional[MacroRegimeHMM] = None
    M_macro:        Optional[int] = None
    macro_seq:      Optional[np.ndarray] = None       # Viterbi decoded m_{1:T}
    market_hmm:     Optional[HierarchicalMarketHMM] = None
    market_K:       Optional[int] = None
    refit_count:    int = 0
    refit_failures: int = 0


# --------------------------------------------------------------------- #
# 1b. Template-to-K initialization helper
# --------------------------------------------------------------------- #
def _build_K_init_from_templates(
    K: int,
    MU_tpl: np.ndarray,
    SIG_tpl: np.ndarray,
    is_valid_tpl: np.ndarray,
    random_state: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Produce K initial (means, covs) for the market HMM by selecting
    or aggregating from the G fixed templates.

    Why this exists
    ~~~~~~~~~~~~~~~
    When `frozen_emissions=False` in the hierarchical backtest, the
    market HMM must learn K-dimensional emissions. Naive k-means init
    (on X) under macro-gated EM collapses all K components toward the
    global mean — diagnosed empirically in v10.0 unfrozen test runs.

    The fix is to **initialize EM from the baseline templates**: K
    components start at well-separated regime locations, then EM is
    free to update them within their attraction basins. This produces
    HMM components that map cleanly back to templates via W2.

    Strategy by K vs G
    ~~~~~~~~~~~~~~~~~~
      * K == G: use templates 1-to-1
      * K < G : k-means-cluster the G template means into K groups
                and use the per-cluster mean (mean of grouped templates)
                + mean of grouped covariances (mixture-style)
      * K > G : use all G templates plus (K-G) duplicates with
                slight noise on the means (EM will separate them)

    All cases skip invalid templates (those with insufficient
    calibration data) and use the global empirical mean/cov as
    fallback for any unfilled component.
    """
    from sklearn.cluster import KMeans

    d = MU_tpl.shape[1]
    eye_d = np.eye(d)
    valid_idx = np.where(is_valid_tpl)[0]
    G_valid = len(valid_idx)

    if G_valid == 0:
        raise ValueError("No valid templates to initialize K-components.")

    MU_valid  = MU_tpl[valid_idx]
    SIG_valid = SIG_tpl[valid_idx]

    if K == G_valid:
        return MU_valid.copy(), SIG_valid.copy()

    if K < G_valid:
        # k-means cluster the G template means into K groups
        km = KMeans(n_clusters=K, n_init=10, random_state=random_state)
        labels = km.fit_predict(MU_valid)
        means_init = np.zeros((K, d))
        covs_init  = np.zeros((K, d, d))
        for k in range(K):
            grp = np.where(labels == k)[0]
            if len(grp) == 0:
                # Empty cluster — extremely unlikely with k-means
                means_init[k] = MU_valid.mean(axis=0)
                covs_init[k]  = SIG_valid.mean(axis=0) + 1e-6 * eye_d
            else:
                # Mean of the templates in this cluster
                means_init[k] = MU_valid[grp].mean(axis=0)
                # Mixture covariance: avg cov + cross-mean spread
                cov_avg = SIG_valid[grp].mean(axis=0)
                spread  = np.cov(MU_valid[grp].T) if len(grp) > 1 else np.zeros((d, d))
                covs_init[k] = cov_avg + spread + 1e-6 * eye_d
                covs_init[k] = 0.5 * (covs_init[k] + covs_init[k].T)
        return means_init, covs_init

    # K > G: use all G templates + duplicates with small perturbation
    means_init = np.zeros((K, d))
    covs_init  = np.zeros((K, d, d))
    means_init[:G_valid] = MU_valid
    covs_init[:G_valid]  = SIG_valid
    rng = np.random.default_rng(random_state)
    for k in range(G_valid, K):
        # Duplicate a random template with small noise so EM can
        # separate them. Noise scale = 0.1 * sqrt(diag(template cov))
        src = k % G_valid
        noise_scale = 0.1 * np.sqrt(np.diag(SIG_valid[src]))
        means_init[k] = MU_valid[src] + rng.normal(0, noise_scale)
        covs_init[k]  = SIG_valid[src].copy() + 1e-6 * eye_d
    return means_init, covs_init


# --------------------------------------------------------------------- #
# 2. Per-day step: macro decode + market fit/predict
# --------------------------------------------------------------------- #
def _step_regime_hierarchical(
    state:              HierarchicalRegimeState,
    t_i:                int,
    X:                  np.ndarray,
    s_macro:            np.ndarray,
    f_hmm:              int,
    macro_f_hmm_mult:   int = 1,
    dirichlet_alpha:    float = 1.0,
    market_n_iter:      int = 15,
    frozen_emissions:   bool = False,
    mean_anchor_strength: float = 0.0,
    freeze_covars:        bool  = False,
    verbose:            bool = False,
    date_for_print:     Optional[pd.Timestamp] = None,
) -> tuple[Optional[np.ndarray], Optional[list[int]], Optional[np.ndarray], int]:
    """One OOS step:
        (A)  Viterbi-decode macro sequence m_{1:t} using FROZEN macro HMM.
        (B)  Possibly refit the market HMM on (X, m_{1:t}).
        (C)  predict_proba on (X, m_{1:t}); map to templates; return pG.
    """
    K_prev = state.K_curr

    # ---- (A) Macro Viterbi decode ----------------------------- #
    # Macro HMM is frozen from calibration; we just run inference
    # on the macro history up to today.
    if state.macro_hmm is None:
        # No macro HMM yet — fall back to a single macro regime
        # (i.e., baseline behavior).
        m_seq = np.zeros(s_macro.shape[0], dtype=np.int64)
        M = 1
    else:
        try:
            m_seq = state.macro_hmm.predict_viterbi(s_macro).astype(np.int64)
            M = state.macro_hmm.M_
        except Exception:
            # If macro decode fails (very rare), use cached or zeros
            if state.macro_seq is not None and state.macro_seq.shape[0] >= s_macro.shape[0]:
                m_seq = state.macro_seq[:s_macro.shape[0]]
            else:
                m_seq = np.zeros(s_macro.shape[0], dtype=np.int64)
            M = state.M_macro or 1
    state.macro_seq = m_seq.copy()

    # ---- (B) K selection + market HMM fit/reuse --------------- #
    if frozen_emissions:
        # K=G enforced when emissions are frozen to templates
        G = state.G
        if G == 0:
            return None, None, None, 0
        state.K_curr = G
        state.K_candidate = G
    else:
        is_K_day = (t_i % HMM.f_k == 0) or (t_i == 0)
        if is_K_day:
            K_candidates = list(range(HMM.min_regimes, HMM.max_regimes + 1))
            state.K_candidate = select_K_predictive(
                X_all=X, K_candidates=K_candidates,
                l_val=HMM.l_val, n_iter=HMM.n_iter,
                random_state=HMM.random_state, lam_k=None,
            )
        if HMM.monotone_K:
            state.K_curr = max(state.K_curr, state.K_candidate)
        else:
            state.K_curr = int(state.K_candidate)

    if verbose and state.K_curr != K_prev and date_for_print is not None:
        print(f"  {pd.Timestamp(date_for_print).date()}  "
              f"K: {K_prev} -> {state.K_curr}", flush=True)

    refit_cadence = f_hmm * max(1, macro_f_hmm_mult)
    need_refit = (
        state.market_hmm is None
        or state.market_K != state.K_curr
        or (t_i % refit_cadence == 0)
    )

    if need_refit:
        if frozen_emissions:
            fee_means = state.MU_tpl_feat
            fee_covs  = state.SIG_tpl_feat
            ti_means  = None
            ti_covs   = None
        else:
            fee_means = None
            fee_covs  = None
            # Template-anchored initialization (v10.1 fix): start EM
            # from K templates (clustered if K!=G) so emissions stay
            # near interpretable regime locations instead of collapsing
            # toward the global mean under macro-gated EM.
            ti_means, ti_covs = _build_K_init_from_templates(
                K=state.K_curr,
                MU_tpl=state.MU_tpl_feat,
                SIG_tpl=state.SIG_tpl_feat,
                is_valid_tpl=state.is_valid_tpl,
                random_state=HMM.random_state,
            )

        new_hmm = fit_hierarchical_market_hmm(
            X=X, m_seq=m_seq, K=state.K_curr, M=M,
            n_iter=market_n_iter,
            dirichlet_alpha=dirichlet_alpha,
            random_state=HMM.random_state,
            frozen_emissions=frozen_emissions,
            fixed_emission_means=fee_means,
            fixed_emission_covs=fee_covs,
            template_init_means=ti_means,
            template_init_covs=ti_covs,
            mean_anchor_strength=mean_anchor_strength,
            freeze_covars=freeze_covars,
        )
        if new_hmm is not None:
            state.market_hmm = new_hmm
            state.market_K   = state.K_curr
            state.refit_count += 1
        else:
            state.refit_failures += 1

    if state.market_hmm is None:
        return None, None, None, state.G

    # ---- (C) Inference + template mapping --------------------- #
    K_active = int(state.market_K)
    probsK   = state.market_hmm.predict_proba(X, m_seq)
    pK       = probsK[-1].copy()

    if frozen_emissions:
        G = state.MU_tpl_feat.shape[0]
        map_k_to_g = list(range(G))
        pG = pK.copy()
    else:
        MU_K_feat  = state.market_hmm.means_
        SIG_K_feat = state.market_hmm.covars_
        map_k_to_g, _ = map_components_to_fixed_templates(
            MU_K_feat, SIG_K_feat,
            state.MU_tpl_feat, state.SIG_tpl_feat,
            state.is_valid_tpl,
        )
        G = state.MU_tpl_feat.shape[0]
        pG = aggregate_template_posteriors_hard(pK, map_k_to_g, G)

    return pK, map_k_to_g, pG, K_active


# --------------------------------------------------------------------- #
# 3. Hierarchical backtest entry point
# --------------------------------------------------------------------- #
def run_backtest_hierarchical(
    returns:           pd.DataFrame,
    returns_train:     pd.DataFrame,
    returns_test:      pd.DataFrame,
    features_all:      Optional[pd.DataFrame] = None,
    macro_state:       Optional[pd.DataFrame] = None,
    verbose:           bool = RUN.verbose,
    progress_pct:      float = 5.0,
    warmup_start:      Optional[str] = None,
    macro_f_hmm_mult:  int = 1,
    macro_M_candidates: tuple[int, ...] = (3, 4, 5),
    macro_transmat_prior_diag: float = 10.0,
    dirichlet_alpha:   float = 1.0,
    market_n_iter:     int = 15,
    frozen_emissions:  bool = False,
    mean_anchor_strength: float = 100.0,
    freeze_covars:     bool  = True,
) -> tuple[BacktestResult, dict]:
    """Run the hierarchical (two-layer) Wasserstein-HMM + MVO backtest.

    Macro layer is calibrated once on data < OOS start, then frozen.
    Market layer is refit every f_hmm * macro_f_hmm_mult days with
    macro-regime-conditional transition matrices.
    """
    if features_all is None:
        features_all = build_features(returns)

    if macro_state is None:
        macro_state = load_macro_state(feature_index=features_all.index)
    else:
        macro_state = macro_state.reindex(features_all.index, method="ffill")

    # Combine train+test returns sorted
    ret_hist_full = pd.concat([returns_train, returns_test]).sort_index()
    ret_hist_full = ret_hist_full[~ret_hist_full.index.duplicated(keep="last")]
    hist_idx_np   = ret_hist_full.index.values

    feat_idx_np      = features_all.index.values
    feat_values_np   = features_all.to_numpy()
    macro_values_np  = macro_state.to_numpy()

    test_dates = returns_test.index
    n_total = len(test_dates)
    split_date = test_dates[0]

    if warmup_start is None:
        warmup_start = WARMUP_START_DATE

    state = HierarchicalRegimeState(
        K_curr=HMM.min_regimes, K_candidate=HMM.min_regimes,
    )

    # =================================================================
    # WARMUP — same template-building as baseline
    # =================================================================
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
        print("  WARNING: hierarchical extension requires warmup.")
        return BacktestResult(
            pnl=pd.Series(dtype=float, index=test_dates),
            weights=pd.DataFrame(0.0, index=test_dates, columns=ASSET_NAMES),
            cum_pnl=pd.Series(dtype=float, index=test_dates),
            K_history=pd.Series(dtype=float, index=test_dates, name="K"),
            tpl_label=pd.Series(dtype=float, index=test_dates, name="regime"),
            tpl_max_prob=pd.Series(dtype=float, index=test_dates, name="max_p"),
            tpl_count=pd.Series(0, index=test_dates, name="G"),
            turnover=pd.Series(dtype=float, index=test_dates, name="turnover"),
        ), {}

    # =================================================================
    # MACRO HMM CALIBRATION — fit once, then frozen
    # =================================================================
    print("\n  Calibrating macro layer HMM...", flush=True)
    # Macro state is z-scored from data_macro.py and aligned to feature dates.
    # Use ALL macro observations up to the last warmup date (strictly causal).
    cal_cutoff_feat = np.searchsorted(feat_idx_np, np.datetime64(split_date), side="left")
    s_cal = macro_values_np[:cal_cutoff_feat]
    # Drop rows with NaN (the Mulliner warmup period)
    valid = ~np.isnan(s_cal).any(axis=1)
    s_cal_clean = s_cal[valid]
    print(f"    Macro calibration sample: {s_cal_clean.shape[0]} days "
          f"(of {s_cal.shape[0]} total)", flush=True)

    macro_hmm = MacroRegimeHMM(
        M_candidates=macro_M_candidates,
        transmat_prior_diag=macro_transmat_prior_diag,
        n_iter=HMM.n_iter,
        random_state=HMM.random_state,
    )
    macro_hmm.fit_initial(s_cal_clean, verbose=verbose)
    state.macro_hmm = macro_hmm
    state.M_macro = macro_hmm.M_
    print(f"    Macro HMM fitted: M = {macro_hmm.M_}, "
          f"BIC scores = {macro_hmm.bic_scores_}", flush=True)

    # Decode macro regime for the full calibration period for diagnostics
    m_cal_seq = macro_hmm.predict_viterbi(s_cal_clean)
    cal_occ = pd.Series(m_cal_seq).value_counts(normalize=True).sort_index().round(3).to_dict()
    print(f"    Calibration macro occupancy: {cal_occ}", flush=True)

    # =================================================================
    # OOS LOOP
    # =================================================================
    weights_oos:       list = []
    pnl_oos:           list = []
    dates_oos:         list = []
    K_history:         list = []
    tpl_label_history: list = []
    tpl_prob_history:  list = []
    tpl_count_history: list = []
    macro_history:     list = []      # the m_t value at each OOS day

    t_start  = time.time()
    last_log = -progress_pct
    w_prev   = np.zeros(N_ASSETS)

    print(f"\nRunning HIERARCHICAL Wasserstein-HMM + MVO backtest "
          f"(M={state.M_macro} macro regimes, "
          f"K variable, OOS days = {n_total})", flush=True)

    for t_i, date in enumerate(test_dates):
        cutoff_feat = np.searchsorted(feat_idx_np, np.datetime64(date), side="left")
        cutoff_ret  = np.searchsorted(hist_idx_np, np.datetime64(date), side="left")

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
            macro_today = macro_history[-1] if macro_history else -1
            print(f"  [{pct_done:5.1f}%]  {t_i+1:>5}/{n_total}  "
                  f"date {pd.Timestamp(date).date()}  "
                  f"elapsed {elapsed/60:5.1f}m  ETA {eta_min:5.1f}m  "
                  f"K={state.K_curr}  G={state.G}  "
                  f"m={macro_today}  "
                  f"refits={state.refit_count}  {sh_str}", flush=True)
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
        s_macro_today = macro_values_np[:cutoff_feat]
        # Forward-fill any NaN in macro (Mulliner warmup)
        if np.isnan(s_macro_today).any():
            s_macro_today = pd.DataFrame(s_macro_today).ffill().fillna(0).to_numpy()

        pK, map_k_to_g, pG, K_active = _step_regime_hierarchical(
            state=state, t_i=t_i, X=X, s_macro=s_macro_today,
            f_hmm=HMM.f_hmm, macro_f_hmm_mult=macro_f_hmm_mult,
            dirichlet_alpha=dirichlet_alpha,
            market_n_iter=market_n_iter,
            frozen_emissions=frozen_emissions,
            mean_anchor_strength=mean_anchor_strength,
            freeze_covars=freeze_covars,
            verbose=verbose, date_for_print=date,
        )

        m_today = int(state.macro_seq[-1]) if state.macro_seq is not None else -1

        if pK is None:
            r_today = returns.loc[date].to_numpy()
            pnl     = float(w_prev @ r_today)
            weights_oos.append(w_prev.copy())
            pnl_oos.append(pnl); dates_oos.append(date)
            K_history.append(np.nan)
            tpl_label_history.append(np.nan)
            tpl_prob_history.append(np.nan)
            tpl_count_history.append(state.G)
            macro_history.append(m_today)
            continue

        # MVO using the same template-mixture moments as baseline
        mu_t    = state.MU_ret_tpl.T @ pG
        Sigma_t = np.tensordot(pG, state.SIG_ret_tpl, axes=(0, 0))

        w_t = solve_mvo(mu_t, Sigma_t, w_prev,
                        lam=MVO.lam, tc=MVO.tc, w_max=MVO.w_max,
                        solver=MVO.solver)

        r_today = returns.loc[date].to_numpy()
        pnl = float(w_t @ r_today)

        # Slow EMA template update (same as baseline/v9)
        _slow_ema_update_templates(
            state=state, pG=pG, map_k_to_g=map_k_to_g, pK=pK,
            X_today=X[-1] if X.shape[0] > 0 else None,
            cached_hmm=state.market_hmm,
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

    if state.refit_failures > 0:
        print(f"  [note] {state.refit_failures} market HMM refit(s) failed "
              f"({state.refit_failures/(state.refit_count+1):.1%} of attempts); "
              f"cached HMM was reused.", flush=True)

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
        "final_market_hmm":   state.market_hmm,
        "macro_hmm":          state.macro_hmm,
        "macro_state":        macro_state,
        "macro_seq_oos":      pd.Series(macro_history, index=pnl.index, name="macro"),
        "refit_count":        state.refit_count,
        "refit_failures":     state.refit_failures,
        "templates": {
            "MU_tpl_feat":   state.MU_tpl_feat,
            "SIG_tpl_feat":  state.SIG_tpl_feat,
            "MU_ret_tpl":    state.MU_ret_tpl,
            "SIG_ret_tpl":   state.SIG_ret_tpl,
            "is_valid_tpl":  state.is_valid_tpl,
        },
    }
    return result, diagnostics
