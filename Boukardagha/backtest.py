"""
backtest.py
-----------
* `solve_mvo`        : transaction-cost-aware mean-variance optimizer (CVXPY/OSQP).
* `_step_regime`     : single-day regime-state update shared by warmup and OOS.
* `_warmup_templates`: runs the daily template-update logic across the pre-OOS
                       calibration window. Implements Algorithm 1, line 3 of
                       Boukardagha (2026) ("Initialize templates using an
                       initial calibration window") which the released code
                       skips. See `config.WARMUP_START_DATE`.
* `run_backtest`     : the strictly causal expanding-window backtest of the
                       Wasserstein-HMM strategy described in Algorithm 1.

Speed-relevant deviation from the reference notebook
----------------------------------------------------
The reference code recomputes the full feature matrix (rolling vol + rolling
momentum over the entire history) **inside** the daily loop. We instead build
features once on the full return series and slice them, which is identical
(rolling windows are causal) but ~700× faster across a 700-day OOS window.
Methodology — HMM hyperparameters, K selection rule, monotone-K toggle,
template tracking, EMA rate, MVO solve — is unchanged.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import cvxpy as cp
import numpy as np
import pandas as pd

from config import HMM, MVO, N_ASSETS, RUN, ASSET_NAMES
from data import build_features
from wasserstein_hmm import (aggregate_template_posteriors_hard,
                             fit_fixed_templates_from_history,
                             fit_hmm,
                             map_components_to_fixed_templates,
                             select_K_predictive)


# --------------------------------------------------------------------- #
# 1.  Mean–variance optimizer
# --------------------------------------------------------------------- #
def solve_mvo(
    mu:     np.ndarray,
    Sigma:  np.ndarray,
    w_prev: np.ndarray,
    lam:    float = MVO.lam,
    tc:     float = MVO.tc,
    w_max:  float = MVO.w_max,
    solver: str   = MVO.solver,
) -> np.ndarray:
    """Solve
        max   µᵀw − γ wᵀΣw − τ ‖w − w_prev‖₁
        s.t.  1ᵀw = 1,  w ≥ 0,  ‖w‖∞ ≤ w_max

    Falls back to SCS if OSQP fails. Always returns a feasible nonneg-summing
    vector; if every solve fails, returns the previous weights unchanged.
    """
    mu = np.asarray(mu).ravel()
    n  = mu.shape[0]
    Sigma = 0.5 * (Sigma + Sigma.T) + 1e-8 * np.eye(n)

    w = cp.Variable(n)
    objective = cp.Maximize(
        mu @ w
        - lam * cp.quad_form(w, cp.psd_wrap(Sigma))
        - tc * cp.norm1(w - w_prev)
    )
    constraints = [cp.sum(w) == 1, w >= 0, cp.norm_inf(w) <= w_max]
    prob = cp.Problem(objective, constraints)

    try:
        prob.solve(solver=getattr(cp, solver), verbose=False)
        if w.value is None:
            prob.solve(solver=cp.SCS, verbose=False)
    except Exception:
        try:
            prob.solve(solver=cp.SCS, verbose=False)
        except Exception:
            return w_prev.copy()

    if w.value is None:
        return w_prev.copy()

    w_sol = np.maximum(np.asarray(w.value).ravel(), 0.0)
    s = w_sol.sum()
    return w_prev.copy() if s <= 1e-12 else w_sol / s


# --------------------------------------------------------------------- #
# 2.  Backtest result container
# --------------------------------------------------------------------- #
@dataclass
class BacktestResult:
    """Bundles every per-day diagnostic the reporting layer needs."""
    pnl:            pd.Series           # daily portfolio log return
    weights:        pd.DataFrame        # (T x N) weights
    cum_pnl:        pd.Series           # cumulative sum of `pnl`
    K_history:      pd.Series           # selected K per day
    tpl_label:      pd.Series           # hard template label per day
    tpl_max_prob:   pd.Series           # max template posterior per day
    tpl_count:      pd.Series           # # templates active per day
    turnover:       pd.Series           # ½‖Δw‖₁ per day

    # --- Convenience -------------------------------------------------- #
    @property
    def daily_returns(self) -> pd.Series:
        """Alias: same as `pnl` (logs are additive ⇒ ≈ simple returns daily)."""
        return self.pnl


# --------------------------------------------------------------------- #
# 2b.  Per-day regime state, shared by warmup and OOS loops
# --------------------------------------------------------------------- #
@dataclass
class RegimeState:
    """Persistent state carried day-to-day in the OOS loop.

    Templates (MU_tpl_feat, SIG_tpl_feat, MU_ret_tpl, SIG_ret_tpl,
    is_valid_tpl) are initialised from the calibration HMM and then
    *slowly* drift via a mass-gated EMA (see `_slow_ema_update_templates`
    and `HMM.eta_tpl`, `HMM.tpl_update_thresh`).

    The mass-gated EMA prevents the "dominant-template absorbs everything"
    collapse pathology observed in earlier iterations: a template only
    updates when its posterior pG[g] exceeds the threshold, so rare
    regimes drift slowly while frequent regimes can still adapt.
    """
    K_curr:           int                      # most-recently-chosen K
    K_candidate:      int                      # last K candidate from selector
    hmm_cached:       Optional[object] = None
    hmm_K:            Optional[int]    = None
    fail_count:       int              = 0

    # --- Templates (initialised from warmup, slowly EMA-updated) ------ #
    MU_tpl_feat:  Optional[np.ndarray] = None  # (G, d_feat) feature-space means
    SIG_tpl_feat: Optional[np.ndarray] = None  # (G, d_feat, d_feat) feature-space covs
    MU_ret_tpl:   Optional[np.ndarray] = None  # (G, N) forward-return means
    SIG_ret_tpl:  Optional[np.ndarray] = None  # (G, N, N) forward-return covs
    is_valid_tpl: Optional[np.ndarray] = None  # (G,) bool

    # --- Carry-state for the strictly-causal forward-return EMA ------- #
    # We need yesterday's regime posterior and yesterday's realized return
    # (= today's realized return at MVO time; tomorrow's pre-OOS step uses
    # both to update each template's forward-return moments).
    pG_prev_for_ema:  Optional[np.ndarray] = None   # (G,) yesterday's pG
    r_prev_for_ema:   Optional[np.ndarray] = None   # (N,) yesterday's realized return

    @property
    def G(self) -> int:
        return 0 if self.MU_tpl_feat is None else int(self.MU_tpl_feat.shape[0])


def _step_regime_fixed_templates(
    state:          RegimeState,
    t_i:            int,
    X:              np.ndarray,
    f_hmm:          int,
    verbose:        bool = False,
    date_for_print: Optional[pd.Timestamp] = None,
) -> tuple[Optional[np.ndarray], Optional[list[int]], Optional[np.ndarray], int]:
    """Single OOS-day update with FIXED templates.

    Compared to the previous dynamic-template variant, this:
      * fits the HMM each refit day,
      * computes HMM emission parameters in FEATURE space (the HMM's
        own state distributions, not derived forward-return moments),
      * W2-maps each emission to the nearest FIXED template,
      * aggregates posteriors via hard argmin (paper-style).

    No template spawn, no EMA — templates are immutable after warmup.
    This is what fixes the "all components map to one template" collapse.

    Returns
    -------
    pK         : (K,) HMM posterior
    map_k_to_g : list of length K
    pG         : (G,) template posterior over the fixed pool
    K_active   : int, the HMM dim actually used today
    """
    # ---- (A) K-selection ---------------------------------------- #
    is_K_day = (t_i % HMM.f_k == 0) or (t_i == 0)
    if is_K_day:
        K_candidates = list(range(HMM.min_regimes, HMM.max_regimes + 1))
        state.K_candidate = select_K_predictive(
            X_all=X,
            K_candidates=K_candidates,
            l_val=HMM.l_val,
            n_iter=HMM.n_iter,
            random_state=HMM.random_state,
            lam_k=None,
        )
    K_prev = state.K_curr
    if HMM.monotone_K:
        state.K_curr = max(state.K_curr, state.K_candidate)
    else:
        state.K_curr = int(state.K_candidate)
    if verbose and state.K_curr != K_prev and date_for_print is not None:
        print(f"  {pd.Timestamp(date_for_print).date()}  "
              f"K: {K_prev} → {state.K_curr}", flush=True)

    # ---- (B) Fit OR reuse HMM ----------------------------------- #
    need_refit = (
        state.hmm_cached is None
        or state.hmm_K != state.K_curr
        or (t_i % f_hmm == 0)
    )
    if need_refit:
        new_hmm = fit_hmm(X, K=state.K_curr, n_iter=HMM.n_iter,
                          random_state=HMM.random_state)
        if new_hmm is not None:
            state.hmm_cached = new_hmm
            state.hmm_K      = state.K_curr
        else:
            state.fail_count += 1

    if state.hmm_cached is None:
        return None, None, None, state.G

    K_active   = int(state.hmm_K)
    probsK     = state.hmm_cached.predict_proba(X)
    pK         = probsK[-1].copy()

    # ---- (C) HMM emission parameters in FEATURE space ----------- #
    MU_K_feat  = state.hmm_cached.means_     # (K, d_feat)
    SIG_K_feat = state.hmm_cached.covars_    # (K, d_feat, d_feat)

    # ---- (D) Map to FIXED templates via W2 in feature space ----- #
    map_k_to_g, _ = map_components_to_fixed_templates(
        MU_K_feat, SIG_K_feat,
        state.MU_tpl_feat, state.SIG_tpl_feat,
        state.is_valid_tpl,
    )
    G = state.MU_tpl_feat.shape[0]
    pG = aggregate_template_posteriors_hard(pK, map_k_to_g, G)

    return pK, map_k_to_g, pG, K_active


# --------------------------------------------------------------------- #
# 2b-bis.  Slow, mass-gated EMA template update
# --------------------------------------------------------------------- #
def _select_G_predictive(X_warm: np.ndarray) -> int:
    """BIC-style selection of the calibration-HMM order G.

    Same predictive-LL + BIC penalty as `select_K_predictive` is used
    for the daily K-selector, applied once on the full calibration
    window. Validation slice is the last `HMM.l_val * 4` rows of
    `X_warm` — a longer window than the daily K-selector uses because
    we have the luxury of doing this fit only once.

    Returns G in [HMM.G_min, HMM.G_max] that maximises
        PredLL(G) − λ_G · G
    with λ_G = lam_G_scale · 0.5 · log(n_val) · d_feat (matches the
    daily lam_k formula).
    """
    from wasserstein_hmm import _fit_and_score, _compute_lam_k

    n = len(X_warm)
    d_feat = X_warm.shape[1]
    # Use a generous validation slice — we only do this once.
    l_val_G = min(HMM.l_val * 4, max(HMM.l_val, n // 5))
    X_fit = X_warm[:-l_val_G]
    X_val = X_warm[-l_val_G:]

    lam_G = HMM.lam_G_scale * 0.5 * float(np.log(max(l_val_G, 2))) * float(d_feat)
    print(f"  G-selection: n_fit={len(X_fit)}, n_val={len(X_val)}, "
          f"lam_G={lam_G:.1f} (per state)", flush=True)

    G_candidates = list(range(HMM.G_min, HMM.G_max + 1))
    scored = []
    for G in G_candidates:
        # Use the standard helper that scores PredLL - lam · G
        _, score = _fit_and_score(
            K=G, X_fit=X_fit, X_val=X_val,
            lam_k=lam_G, n_iter=HMM.n_iter,
            random_state=HMM.random_state,
        )
        # _fit_and_score returns PredLL - lam_k * G, so score directly comparable
        scored.append((G, score))
        print(f"    G={G}: score={score:.1f}")
    best_G, best_score = max(scored, key=lambda gs: gs[1])
    return int(best_G)


def _slow_ema_update_templates(
    state:        RegimeState,
    pG:           np.ndarray,
    map_k_to_g:   list[int],
    pK:           np.ndarray,
    X_today:      Optional[np.ndarray],
    cached_hmm:   object,
    r_today:      np.ndarray,
    r_prev:       Optional[np.ndarray],
    pG_prev:      Optional[np.ndarray],
) -> None:
    """Slow, mass-gated EMA update of templates.

    Two updates happen here, both gated on template posterior mass to
    prevent the "dominant template absorbs everything" collapse:

    (1) FEATURE-space templates (MU_tpl_feat, SIG_tpl_feat):
        Each template g whose CURRENT posterior pG[g] >= tpl_update_thresh
        is pulled (with rate eta_tpl) toward the weighted average of HMM
        components mapped to g today. Templates not assigned today don't
        move at all. Rare regimes therefore drift much slower than
        frequent ones.

    (2) FORWARD-RETURN moments (MU_ret_tpl, SIG_ret_tpl):
        Strictly causal: uses YESTERDAY's regime posterior pG_prev
        and TODAY's realized return r_today (since r_today is what
        "forward return given yesterday's regime" actually means).
        Each template g with pG_prev[g] >= tpl_update_thresh is pulled
        toward r_today (mean) and toward outer(r_today - mu_old, ...) for
        covariance. Standard online-mean / online-covariance updates with
        the SAME small eta_tpl.

    Setting `eta_tpl = 0` disables this update entirely (recovers v6
    fully-frozen behaviour).
    """
    eta = HMM.eta_tpl
    thr = HMM.tpl_update_thresh
    if eta <= 0:
        return  # frozen templates

    # ---------- (1) Feature-space template update ---------- #
    if X_today is not None and cached_hmm is not None:
        MU_K_feat = cached_hmm.means_      # (K, d_feat)
        SIG_K_feat = cached_hmm.covars_    # (K, d_feat, d_feat)
        N_feat = MU_K_feat.shape[1]
        eye_d = np.eye(N_feat)

        for g in range(state.G):
            if not state.is_valid_tpl[g]:
                continue
            if pG[g] < thr:
                continue
            # Which HMM components mapped to template g?
            ks = [k for k, gk in enumerate(map_k_to_g) if gk == g]
            if not ks:
                continue
            w = np.asarray([pK[k] for k in ks], dtype=float)
            ws = w.sum()
            if ws <= 1e-12:
                continue
            w = w / ws

            mu_bar = (w[:, None] * MU_K_feat[ks]).sum(axis=0)
            S_bar  = np.einsum("k,kij->ij", w, SIG_K_feat[ks])

            state.MU_tpl_feat[g]  = (1 - eta) * state.MU_tpl_feat[g]  + eta * mu_bar
            state.SIG_tpl_feat[g] = (1 - eta) * state.SIG_tpl_feat[g] + eta * S_bar
            # Numerical hygiene
            state.SIG_tpl_feat[g] = 0.5 * (state.SIG_tpl_feat[g] + state.SIG_tpl_feat[g].T) + 1e-10 * eye_d

    # ---------- (2) Forward-return moment update ---------- #
    # Causality: use yesterday's regime posterior + today's realized return.
    # This is exactly "what forward return given yesterday's regime"; we
    # update each regime's mean and covariance estimate online.
    if r_prev is None or pG_prev is None:
        return
    # Note: by convention we treat "r_today" as the realised return
    # whose precursor regime distribution was "pG_prev".  This naming
    # follows the standard online-update for r_{s+1} | z_s = g.
    r_realised = np.asarray(r_today, dtype=float).ravel()
    N_ret = r_realised.shape[0]
    eye_n = np.eye(N_ret)

    for g in range(state.G):
        if not state.is_valid_tpl[g]:
            continue
        if pG_prev[g] < thr:
            continue
        # Online EMA of mean
        mu_old = state.MU_ret_tpl[g].copy()
        state.MU_ret_tpl[g]  = (1 - eta) * mu_old + eta * r_realised
        # Online EMA of covariance using deviation from the OLD mean
        # (Welford-style; the deviation from old mean is what gets
        # absorbed into the running second moment in standard online
        # covariance updates).
        dev = r_realised - mu_old
        S_new = np.outer(dev, dev)
        state.SIG_ret_tpl[g] = (1 - eta) * state.SIG_ret_tpl[g] + eta * S_new
        state.SIG_ret_tpl[g] = 0.5 * (state.SIG_ret_tpl[g] + state.SIG_ret_tpl[g].T) + 1e-10 * eye_n


# --------------------------------------------------------------------- #
# 2c.  Template warmup pass — Algorithm 1 line 3 of Boukardagha (2026)
# --------------------------------------------------------------------- #
def _build_fixed_templates(
    state:           RegimeState,
    warmup_dates:    pd.DatetimeIndex,
    feat_idx_np:     np.ndarray,
    feat_values_np:  np.ndarray,
    features_index:  pd.DatetimeIndex,
    returns_train:   pd.DataFrame,
    verbose:         bool = False,
) -> RegimeState:
    """Fit the calibration HMM ONCE on the full pre-OOS data and store
    its (G, d_feat) emission parameters as the template pool.

    G is selected on the calibration window via a predictive-LL +
    BIC-style penalty, using the same criterion as the daily K-selector
    (see HMMConfig.G_min, G_max, lam_G_scale). This avoids hard-coding
    G and makes the choice defensible (same selection rule applied at
    a different scale).

    Forward-return moments per template are also pre-computed once and
    stored on the state. They subsequently drift slowly via mass-gated
    EMA in the OOS loop (see `_slow_ema_update_templates`,
    `HMM.eta_tpl`, `HMM.tpl_update_thresh`).

    This is the implementation of Algorithm 1 line 3 of
    Boukardagha (2026): "Initialize templates using an initial
    calibration window." We use the entire pre-OOS history as that
    calibration window.
    """
    n = len(warmup_dates)
    if n == 0:
        print("  Warmup disabled (no calibration window).")
        return state

    print(f"  Warmup pass — building templates on the calibration window "
          f"({warmup_dates[0].date()} → {warmup_dates[-1].date()}, {n} days)…",
          flush=True)
    t0 = time.time()

    # Use ALL features up to the day before OOS begins (strictly causal).
    last_warmup_date = warmup_dates[-1]
    cutoff_feat = np.searchsorted(feat_idx_np, np.datetime64(last_warmup_date), side="right")
    X_warm = feat_values_np[:cutoff_feat]
    feat_dates_warm = features_index[:cutoff_feat]

    # ----- BIC-style G selection on the calibration window ----- #
    G_chosen = _select_G_predictive(X_warm)
    print(f"  G selected via predictive-LL + BIC penalty: G = {G_chosen} "
          f"(from candidates {HMM.G_min}..{HMM.G_max})", flush=True)

    print(f"  Fitting calibration HMM with G={G_chosen} on {len(X_warm)} pre-OOS feature rows...")
    MU_feat, SIG_feat, MU_ret, SIG_ret, is_valid = fit_fixed_templates_from_history(
        X=X_warm,
        returns_df=returns_train,
        feat_dates=feat_dates_warm,
        K_template=G_chosen,
        n_iter=HMM.n_iter,
        random_state=HMM.random_state,
        min_obs=HMM.min_obs,
        n_assets=N_ASSETS,
    )

    state.MU_tpl_feat  = MU_feat
    state.SIG_tpl_feat = SIG_feat
    state.MU_ret_tpl   = MU_ret
    state.SIG_ret_tpl  = SIG_ret
    state.is_valid_tpl = is_valid

    print(f"  Warmup complete: G={state.G} templates initialised "
          f"({is_valid.sum()} valid, {(~is_valid).sum()} dropped due to insufficient data). "
          f"({(time.time()-t0)/60:.1f}m)",
          flush=True)
    return state


# --------------------------------------------------------------------- #
# 3.  Main backtest loop
# --------------------------------------------------------------------- #
def run_backtest(
    returns:       pd.DataFrame,
    returns_train: pd.DataFrame,
    returns_test:  pd.DataFrame,
    features_all:  Optional[pd.DataFrame] = None,
    verbose:       bool = RUN.verbose,
    progress_pct:  float = 5.0,
    warmup_start:  Optional[str] = None,
    f_hmm_warmup:  Optional[int] = None,
) -> BacktestResult:
    """Run the strictly causal expanding-window Wasserstein-HMM + MVO backtest.

    Parameters
    ----------
    returns
        Full daily log-return DataFrame (used to look up the realized r_t).
    returns_train, returns_test
        The train/test split. Test index drives the OOS loop.
    features_all
        Optional precomputed feature matrix on the full `returns`. If `None`
        it is built here. **Strongly recommended** for speed.
    verbose
        If True, print K-state transitions (rare).
    progress_pct
        Print a status line every this-many percent done (default 5.0).
        Set to 0 or negative to silence progress entirely.
    warmup_start
        Optional ISO date string. If provided (and < SPLIT_DATE), the
        template-update logic runs across the pre-OOS window from this
        date to SPLIT_DATE-1, implementing Algorithm 1 line 3 of the
        paper. If None, the function reads `config.WARMUP_START_DATE`.
        Pass `warmup_start = str(SPLIT_DATE)` or any later date to skip
        warmup and recover the released-code behaviour.
    f_hmm_warmup
        HMM refit cadence during warmup (default = config.WARMUP_F_HMM).

    Returns
    -------
    `BacktestResult` — see the dataclass above.
    """
    if features_all is None:
        features_all = build_features(returns)

    # Pre-build the concatenated history once. Inside the loop we use np.searchsorted
    # to find "all dates strictly before today" — O(log n) instead of pd.concat+sort.
    ret_hist_full   = pd.concat([returns_train, returns_test]).sort_index()
    ret_hist_full   = ret_hist_full[~ret_hist_full.index.duplicated(keep="last")]
    hist_idx_np     = ret_hist_full.index.values

    feat_idx_np     = features_all.index.values
    feat_values_np  = features_all.to_numpy()

    test_dates = returns_test.index
    n_total = len(test_dates)
    split_date = test_dates[0]

    # ===================================================================
    # WARMUP PASS — Algorithm 1, line 3 of Boukardagha (2026)
    # ===================================================================
    import config as _cfg
    if warmup_start is None:
        warmup_start = _cfg.WARMUP_START_DATE

    # Initialise the persistent state.
    state = RegimeState(K_curr=HMM.min_regimes, K_candidate=HMM.min_regimes)

    # Determine warmup date range
    warmup_dates: pd.DatetimeIndex
    if warmup_start is None:
        # Auto-mode: start as early as the K-selector can run.
        min_feat_needed = HMM.l_val + 50
        if min_feat_needed < len(features_all):
            auto_start = features_all.index[min_feat_needed]
            warmup_dates = returns_train.index[returns_train.index >= auto_start]
        else:
            warmup_dates = returns_train.index[:0]
    else:
        ws = pd.Timestamp(warmup_start)
        if ws >= split_date:
            warmup_dates = returns_train.index[:0]   # disabled
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
        print("  WARNING: warmup disabled — fixed-template mode requires "
              "a calibration window. Falling back to no-template mode "
              "(weights will hold).", flush=True)
        return BacktestResult(
            pnl=pd.Series(dtype=float, index=test_dates),
            weights=pd.DataFrame(0.0, index=test_dates, columns=ASSET_NAMES),
            cum_pnl=pd.Series(dtype=float, index=test_dates),
            K_history=pd.Series(dtype=float, index=test_dates, name="K"),
            tpl_label=pd.Series(dtype=float, index=test_dates, name="regime"),
            tpl_max_prob=pd.Series(dtype=float, index=test_dates, name="max_p"),
            tpl_count=pd.Series(0, index=test_dates, name="G"),
            turnover=pd.Series(dtype=float, index=test_dates, name="turnover"),
        )

    # ===================================================================
    # OOS LOOP
    # ===================================================================
    weights_oos:        list[np.ndarray] = []
    pnl_oos:            list[float]      = []
    dates_oos:          list[pd.Timestamp] = []
    K_history:          list[float]      = []
    tpl_label_history:  list[float]      = []
    tpl_prob_history:   list[float]      = []
    tpl_count_history:  list[int]        = []

    t_start   = time.time()
    last_log  = -progress_pct
    w_prev    = np.zeros(N_ASSETS)

    if verbose:
        print(f"\nRunning OOS backtest "
              f"(strictly causal, OOS days = {len(test_dates)})\n")

    for t_i, date in enumerate(test_dates):
        cutoff_feat = np.searchsorted(feat_idx_np, np.datetime64(date), side="left")
        cutoff_ret  = np.searchsorted(hist_idx_np, np.datetime64(date), side="left")

        # Periodic progress log
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
            print(f"  [{pct_done:5.1f}%]  {t_i+1:>5}/{n_total}  "
                  f"date {pd.Timestamp(date).date()}  "
                  f"elapsed {elapsed/60:5.1f}m  ETA {eta_min:5.1f}m  "
                  f"K={state.K_curr}  G={state.G}  {sh_str}", flush=True)
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
            continue

        X = feat_values_np[:cutoff_feat]

        # Step the regime state forward
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
            continue

        # ---- (D) Template-mixture moments + MVO ---------------------- #
        # Forward-return moments per template come from the calibration
        # window. The slow EMA below adapts BOTH the feature-space templates
        # (for W2 distance — preserves identity but allows gentle drift)
        # AND the forward-return moments (for MVO inputs — they should
        # absorb new realized return data).
        #
        # Strict causality: forward-return moments are updated using
        # YESTERDAY's regime posterior and TODAY's realized return; we
        # apply this update AFTER computing today's weights, not before.
        mu_t    = state.MU_ret_tpl.T @ pG
        Sigma_t = np.tensordot(pG, state.SIG_ret_tpl, axes=(0, 0))

        w_t = solve_mvo(mu_t, Sigma_t, w_prev,
                        lam=MVO.lam, tc=MVO.tc, w_max=MVO.w_max,
                        solver=MVO.solver)

        r_today = returns.loc[date].to_numpy()
        pnl = float(w_t @ r_today)

        # Slow EMA update of templates — mass-gated to prevent the
        # "dominant template absorbs everything" collapse.
        _slow_ema_update_templates(
            state=state,
            pG=pG,
            map_k_to_g=map_k_to_g,
            pK=pK,
            X_today=X[-1] if X.shape[0] > 0 else None,
            cached_hmm=state.hmm_cached,
            r_today=r_today,
            r_prev=state.r_prev_for_ema,
            pG_prev=state.pG_prev_for_ema,
        )
        # Stash today's regime info for tomorrow's forward-return update
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
        w_prev = w_t

    # End-of-run note about fit fallbacks, if any.
    if state.fail_count > 0:
        print(f"  [note] {state.fail_count} HMM refit(s) fell back to the cached "
              f"model due to non-PSD covariance "
              f"(~{state.fail_count/len(test_dates):.1%} of OOS days).",
              flush=True)

    # ------- Assemble outputs ------- #
    pnl = pd.Series(pnl_oos, index=pd.to_datetime(dates_oos)).sort_index()
    W   = pd.DataFrame(weights_oos, index=pnl.index, columns=ASSET_NAMES)
    turnover = 0.5 * W.diff().abs().sum(axis=1)

    return BacktestResult(
        pnl=pnl,
        weights=W,
        cum_pnl=pnl.cumsum(),
        K_history    =pd.Series(K_history,         index=pnl.index, name="K"),
        tpl_label    =pd.Series(tpl_label_history, index=pnl.index, name="regime"),
        tpl_max_prob =pd.Series(tpl_prob_history,  index=pnl.index, name="max_p"),
        tpl_count    =pd.Series(tpl_count_history, index=pnl.index, name="G"),
        turnover     =turnover.fillna(0.0).rename("turnover"),
    )
