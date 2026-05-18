"""
backtest.py
-----------
* `solve_mvo`        : transaction-cost-aware mean-variance optimizer (CVXPY/OSQP).
* `run_backtest`     : the strictly causal expanding-window backtest of the
                       Wasserstein-HMM strategy described in Algorithm 1 of
                       Boukardagha (2026).

Speed-relevant deviation from the reference notebook
----------------------------------------------------
The reference code recomputes the full feature matrix (rolling vol + rolling
momentum over the entire history) **inside** the daily loop. We instead build
features once on the full return series and slice them, which is identical
(rolling windows are causal) but ~700× faster across a 700-day OOS window.
Methodology — HMM hyperparameters, K selection rule, monotone-K safety,
template tracking, EMA rate, MVO solve — is unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import cvxpy as cp
import numpy as np
import pandas as pd
from tqdm import tqdm

from config import HMM, MVO, N_ASSETS, RUN, ASSET_NAMES
from data import build_features
from wasserstein_hmm import (aggregate_template_posteriors,
                             compute_regime_moments_forward, fit_hmm,
                             map_components_to_templates,
                             select_K_predictive, spawn_template_if_needed,
                             update_templates_ema)


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
# 3.  Main backtest loop
# --------------------------------------------------------------------- #
def run_backtest(
    returns:       pd.DataFrame,
    returns_train: pd.DataFrame,
    returns_test:  pd.DataFrame,
    features_all:  Optional[pd.DataFrame] = None,
    verbose:       bool = RUN.verbose,
    use_tqdm:      bool = RUN.tqdm,
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
    iterator = tqdm(test_dates, desc="Backtest", ncols=100) if use_tqdm else test_dates

    # ------- Loop state ------- #
    weights_oos:        list[np.ndarray] = []
    pnl_oos:            list[float]      = []
    dates_oos:          list[pd.Timestamp] = []
    K_history:          list[float]      = []
    tpl_label_history:  list[float]      = []
    tpl_prob_history:   list[float]      = []
    tpl_count_history:  list[int]        = []

    w_prev = np.zeros(N_ASSETS)
    K_curr = HMM.min_regimes
    K_candidate = K_curr
    MU_tpl: Optional[np.ndarray]  = None
    SIG_tpl: Optional[np.ndarray] = None

    if verbose:
        print(f"Running Wasserstein-HMM + MVO backtest "
              f"(strictly causal, OOS days = {len(test_dates)})\n")

    for t_i, date in enumerate(iterator):
        # Index of "first date >= today" in the feature matrix.
        # Everything strictly before is allowed.
        cutoff_feat = np.searchsorted(feat_idx_np, np.datetime64(date), side="left")
        cutoff_ret  = np.searchsorted(hist_idx_np, np.datetime64(date), side="left")

        if cutoff_feat < 2:
            # Not enough history yet → hold previous weights.
            r_today = returns.loc[date].to_numpy()
            pnl     = float(w_prev @ r_today)
            weights_oos.append(w_prev.copy())
            pnl_oos.append(pnl); dates_oos.append(date)
            K_history.append(np.nan)
            tpl_label_history.append(np.nan)
            tpl_prob_history.append(np.nan)
            tpl_count_history.append(0 if MU_tpl is None else MU_tpl.shape[0])
            continue

        X         = feat_values_np[:cutoff_feat]
        feat_idx  = feat_idx_np[:cutoff_feat]
        # Align return rows to feature rows (1-to-1 by date)
        ret_align = ret_hist_full.iloc[:cutoff_ret].loc[feat_idx]

        # ---- (A) Periodic predictive K selection -------------------- #
        if (t_i % HMM.f_k == 0) or (t_i == 0):
            K_candidates = list(range(K_curr, HMM.max_regimes + 1))
            K_candidate = select_K_predictive(
                X_all=X,
                K_candidates=K_candidates,
                l_val=HMM.l_val,
                n_iter=HMM.n_iter,
                random_state=HMM.random_state,
                lam_k=HMM.lam_k,
            )
        K_curr = max(K_curr, K_candidate)     # enforce non-decreasing K

        if verbose:
            print(f"{pd.Timestamp(date).date()}  K = {K_curr}")

        # ---- (B) Fit HMM at chosen K -------------------------------- #
        hmm = fit_hmm(X, K=K_curr, n_iter=HMM.n_iter,
                      random_state=HMM.random_state)
        probsK = hmm.predict_proba(X)
        pK     = probsK[-1].copy()
        zK_raw = np.argmax(probsK, axis=1)

        # Forward-return regime moments
        z_s = zK_raw[:-1]
        MU_K, SIG_K = compute_regime_moments_forward(
            ret_align, z_s, K=K_curr, n_assets=N_ASSETS, min_obs=HMM.min_obs
        )

        # ---- (C) Initialize / spawn / map / update templates -------- #
        if MU_tpl is None or SIG_tpl is None:
            MU_tpl  = MU_K.copy()
            SIG_tpl = SIG_K.copy()
        else:
            MU_tpl, SIG_tpl = spawn_template_if_needed(
                MU_K, SIG_K, pK, MU_tpl, SIG_tpl,
                spawn_thresh=HMM.spawn_thresh, g_max=HMM.g_max,
            )

        map_k_to_g, _ = map_components_to_templates(MU_K, SIG_K, MU_tpl, SIG_tpl)
        G = MU_tpl.shape[0]
        pG = aggregate_template_posteriors(pK, map_k_to_g, G)

        MU_tpl, SIG_tpl = update_templates_ema(
            MU_tpl, SIG_tpl, MU_K, SIG_K, pK, map_k_to_g,
            eta=HMM.eta_tpl,
        )

        # ---- (D) Template-mixture moments + MVO ---------------------- #
        mu_t    = MU_tpl.T @ pG
        Sigma_t = np.tensordot(pG, SIG_tpl, axes=(0, 0))

        w_t = solve_mvo(mu_t, Sigma_t, w_prev,
                        lam=MVO.lam, tc=MVO.tc, w_max=MVO.w_max,
                        solver=MVO.solver)

        # ---- Realize today's PnL ------------------------------------ #
        r_today = returns.loc[date].to_numpy()
        pnl = float(w_t @ r_today)

        weights_oos.append(w_t)
        pnl_oos.append(pnl); dates_oos.append(date)
        K_history.append(K_curr)
        g_hat = int(np.argmax(pG)) if np.all(np.isfinite(pG)) else np.nan
        tpl_label_history.append(g_hat)
        tpl_prob_history.append(float(np.max(pG))
                                if np.all(np.isfinite(pG)) else np.nan)
        tpl_count_history.append(G)
        w_prev = w_t

    # ------- Assemble outputs ------- #
    idx = pd.DatetimeIndex(dates_oos).sort_values()
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
