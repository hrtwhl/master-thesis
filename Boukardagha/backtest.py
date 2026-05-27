"""
backtest.py
===========
Strictly-causal daily backtest engines.

Two main engines:

  run_pure_market_strategy()    - faithful Boukardagha (2026) replication
                                  (expanding window, daily refit, weekly K).
  run_hierarchical_strategy()   - macro-conditional extension
                                  (see hierarchical_hmm.py).

Plus a static-weight engine for passive benchmarks:

  run_static_weight()  - executes a fixed weight vector daily.

All engines produce, for every OOS date d:
    weights    : vector of asset weights chosen using info up to d-1
    pnl        : realised log-return of the portfolio that day
    side_info  : K, G, dominant template, max template posterior, ...

The pure-market engine also returns auxiliary diagnostics needed for
the unified `daily_backtest_output.csv` file (regime, K, max_p, G,
turnover, cum_pnl).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from tqdm import tqdm

from config import (
    ASSET_NAMES, OOS_START,
    LAM, TC, W_MAX, MIN_HIST_FOR_HMM,
    F_REFIT, MAX_TRAIN_WINDOW,
    MIN_REGIMES, MAX_REGIMES, F_K, L_VAL, LAM_K, G_MAX, ETA_TPL, SPAWN_THRESH,
)
from features import build_asset_features
from wasserstein_hmm import (
    WassersteinHMMState, step_wasserstein_hmm,
)
from mvo import solve_mvo


def _cap_window(df: pd.DataFrame, max_rows) -> pd.DataFrame:
    """Optional rolling-window cap.  Returns df unchanged when max_rows
    is None (paper's expanding window)."""
    if max_rows is None or len(df) <= max_rows:
        return df
    return df.iloc[-max_rows:]


# =======================================================================
#  Pure market Wasserstein HMM strategy (Boukardagha replication)
# =======================================================================
def run_pure_market_strategy(returns: pd.DataFrame,
                             oos_start: str = OOS_START,
                             verbose: bool = True) -> dict:
    """
    Faithful replication of Paper_Code.ipynb 'Commercial V2.0' loop.

    Returns
    -------
    dict with keys:
        pnl, weights, K_history, tpl_label, tpl_count, tpl_prob,
        turnover, cum_pnl
    """
    asset_names = list(returns.columns)
    N = len(asset_names)

    oos_start_ts = pd.Timestamp(oos_start)
    oos_dates = returns.index[returns.index >= oos_start_ts]

    state = WassersteinHMMState(
        n_features_returns=N,
        min_regimes=MIN_REGIMES, max_regimes=MAX_REGIMES,
        g_max=G_MAX, eta_tpl=ETA_TPL, spawn_thresh=SPAWN_THRESH,
        f_k=F_K, l_val=L_VAL, lam_k=LAM_K,
        monotone_K=True,
    )

    weights_oos, pnl_oos, dates_oos = [], [], []
    K_history, tpl_label_history = [], []
    tpl_prob_history, tpl_count_history = [], []

    w_prev = np.zeros(N)

    iterator = enumerate(oos_dates)
    if verbose:
        iterator = enumerate(tqdm(oos_dates, desc="Pure market", ncols=100))

    for t_i, date in iterator:
        # ---- 1) Strictly causal history up to date-1 ----
        ret_hist = returns.loc[returns.index < date]
        X_df = build_asset_features(ret_hist)
        X_df = _cap_window(X_df, MAX_TRAIN_WINDOW)

        if len(X_df) < max(MIN_HIST_FOR_HMM, 2):
            w_t = w_prev.copy()
            pnl_t = float(np.dot(w_t, returns.loc[date].values))
            weights_oos.append(w_t)
            pnl_oos.append(pnl_t)
            dates_oos.append(date)
            K_history.append(np.nan)
            tpl_label_history.append(np.nan)
            tpl_prob_history.append(np.nan)
            tpl_count_history.append(0)
            continue

        ret_align = ret_hist.loc[X_df.index]
        X_full    = X_df.values

        # ---- 2) Wasserstein-HMM step ----
        out = step_wasserstein_hmm(
            state, X_full, ret_align,
            step_index=t_i, refit_every=F_REFIT,
        )

        # ---- 3) MVO ----
        w_t = solve_mvo(out["mu_t"], out["Sigma_t"], w_prev,
                        lam=LAM, tc=TC, w_max=W_MAX)

        # ---- 4) Realised PnL on date `date` ----
        pnl_t = float(np.dot(w_t, returns.loc[date].values))

        weights_oos.append(w_t)
        pnl_oos.append(pnl_t)
        dates_oos.append(date)
        K_history.append(out["K"])
        tpl_label_history.append(out["dominant_template"])
        tpl_prob_history.append(out["p_max"])
        tpl_count_history.append(out["G"])
        w_prev = w_t

    return _package_results(
        dates_oos, pnl_oos, weights_oos, asset_names,
        K_history, tpl_label_history, tpl_prob_history, tpl_count_history,
    )


def _package_results(dates, pnl_list, weights_list, asset_names,
                     K_list, lbl_list, prob_list, cnt_list) -> dict:
    """Assemble the diagnostics DataFrame/Series bundle returned by
    the pure-market and hierarchical engines."""
    idx = pd.to_datetime(dates)
    pnl_series = pd.Series(pnl_list, index=idx, name="pnl").sort_index()
    W_df       = pd.DataFrame(weights_list, index=pnl_series.index,
                              columns=asset_names)
    cum_pnl    = pnl_series.cumsum().rename("cum_pnl")
    turnover   = (0.5 * W_df.diff().abs().sum(axis=1)).rename("turnover")
    turnover.iloc[0] = 0.0

    return dict(
        pnl       = pnl_series,
        weights   = W_df,
        cum_pnl   = cum_pnl,
        turnover  = turnover,
        K_history = pd.Series(K_list, index=pnl_series.index, name="K"),
        tpl_label = pd.Series(lbl_list, index=pnl_series.index,
                              name="dominant_template"),
        tpl_prob  = pd.Series(prob_list, index=pnl_series.index,
                              name="max_template_posterior"),
        tpl_count = pd.Series(cnt_list, index=pnl_series.index,
                              name="template_count"),
    )


def make_daily_backtest_csv(result: dict, path: str | None = None) -> pd.DataFrame:
    """
    Build the unified daily backtest output:
        date, w_<asset>, pnl, cum_pnl, K, regime, max_p, G, turnover
    Schema matches the user request.
    """
    W = result["weights"]
    df = pd.DataFrame(index=W.index)
    df.index.name = "date"
    for c in W.columns:
        df[f"w_{c}"] = W[c].values
    df["pnl"]      = result["pnl"].values
    df["cum_pnl"]  = result["cum_pnl"].values
    df["K"]        = result["K_history"].values
    df["regime"]   = result["tpl_label"].values
    df["max_p"]    = result["tpl_prob"].values
    df["G"]        = result["tpl_count"].values
    df["turnover"] = result["turnover"].values
    if path is not None:
        df.to_csv(path)
    return df


# =======================================================================
#  Static-weight passive benchmarks (frictionless, daily rebalance)
# =======================================================================
def run_static_weight(returns: pd.DataFrame,
                      weight_dict: dict,
                      oos_start: str = OOS_START,
                      name: str = "static") -> dict:
    """
    Execute a constant target weight vector daily.
    Boukardagha's paper benchmarks (Equal-Weight, SPX B&H) are also
    treated frictionlessly, so we mirror that.
    """
    asset_names = list(returns.columns)
    w = np.array([weight_dict.get(a, 0.0) for a in asset_names])
    s = w.sum()
    if s > 0:
        w = w / s

    oos_dates = returns.index[returns.index >= pd.Timestamp(oos_start)]
    R = returns.loc[oos_dates, asset_names]

    pnl = R.values @ w
    pnl_series = pd.Series(pnl, index=R.index, name=f"pnl_{name}")
    W_df = pd.DataFrame(
        np.tile(w, (len(R), 1)),
        index=R.index, columns=asset_names,
    )
    return dict(pnl=pnl_series, weights=W_df, target=w)
