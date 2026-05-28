"""
hierarchical_hmm.py
===================
Hierarchical Hidden Markov Model for portfolio construction, following
the architecture of Fine, Singer & Tishby (1998).

Design (REVISED v2)
-------------------
The original v1 fitted a SEPARATE market Wasserstein HMM PER macro
template, which carved the asset history into small subsets and led to
unstable conditional moment estimates.  v2 corrects this with a
data-efficient design:

  Root
   |
   +-- Macro Wasserstein HMM (operates on macro features)
   |       templates g = 1..G_macro
   |
   +-- Market Wasserstein HMM (operates on asset features, SHARED)
           templates h = 1..G_market

At each OOS day t:

  1. Run the macro WHMM step on macro features through t-1.
     -> p_macro(t, g)  for g = 1..G_macro
     -> macro template label for every historical row: g_s

  2. Run the market WHMM step on asset features through t-1.
     -> p_market(t, h)  for h = 1..G_market
     -> market template label for every historical row: h_s

  3. Conditional moments on the JOINT (g, h) cells.
     For each cell (g, h) with N_{g,h} >= MIN_OBS_PER_REGIME, compute
        mu_{g,h}    = mean(r_{s+1} | g_s = g, h_s = h)
        Sigma_{g,h} = LedoitWolf(r_{s+1} | g_s = g, h_s = h)
     Cells with too few obs fall back to the market-only moments
     mu_h, Sigma_h, which use ALL data with market label h.

  4. Composite probability under conditional independence
        p(t, g, h) = p_macro(t, g) * p_market(t, h)

  5. Template-mixture moments fed to MVO
        mu_t    = sum_{g,h} p(t,g,h) mu_{g,h}
        Sigma_t = sum_{g,h} p(t,g,h) Sigma_{g,h}

  6. (Optional, controlled by HIER_MACRO_TILT_STRENGTH)  Apply a
     macro-regime tilt on expected equity returns: scale the SPX/OIL
     entries of mu_t by a factor proportional to historical SPX
     Sharpe in the dominant macro regime, smoothly interpolated by the
     macro posterior.

Both HMMs use the same WassersteinHMMState/step_wasserstein_hmm
infrastructure as the pure-market strategy.  Daily refits are
warm-started so the full backtest stays tractable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf
from tqdm import tqdm

from config import (
    ASSET_NAMES, OOS_START, LAM, TC, W_MAX, F_REFIT, MAX_TRAIN_WINDOW,
    MIN_HIST_FOR_HMM, MIN_OBS_PER_REGIME,
    # Market layer
    MIN_REGIMES, MAX_REGIMES, F_K, L_VAL, LAM_K,
    G_MAX, ETA_TPL, SPAWN_THRESH,
    # Macro layer
    MACRO_MIN_REGIMES, MACRO_MAX_REGIMES, MACRO_G_MAX,
    MACRO_ETA_TPL, MACRO_SPAWN_THRESH, MACRO_F_K, MACRO_F_REFIT,
    # Tilt
    HIER_MACRO_TILT_STRENGTH, HIER_MACRO_TILT_MODE,
    HIER_MACRO_TILT_ASSETS,
)
from features import build_asset_features, build_macro_features
from wasserstein_hmm import (
    WassersteinHMMState, step_wasserstein_hmm,
    fit_gaussian_hmm_warm, select_K_predictive,
    map_components_to_templates, aggregate_template_posteriors,
    spawn_template_if_needed, update_templates_ema,
)
from mvo import solve_mvo
from backtest import _package_results, _cap_window


# =======================================================================
#  Macro Wasserstein HMM step
#  Similar to market step but tracks templates in MACRO FEATURE space
#  (not in forward-return space), since the macro layer's role is to
#  identify slow-moving environment regimes, not to forecast assets.
# =======================================================================
def step_macro_whmm(state: WassersteinHMMState,
                    X_macro_full: np.ndarray,
                    step_index: int,
                    refit_every: int):
    """
    Returns
    -------
    dict with:
        K            : selected number of components
        G            : current macro template count
        pG           : (G,) macro template posteriors at current t
        z_argmax     : (T,) argmax HMM-component label per historical row
        comp_to_tpl  : list len K mapping components to templates
        history_g    : (T,) dominant TEMPLATE label per historical row

    The macro features have wildly different scales (e.g. oil
    increments at O(1) vs stock_bond_corr increments at O(0.001)),
    which makes a full-covariance Gaussian HMM ill-conditioned.  We
    therefore z-score the features INSIDE this step using statistics
    computed over the full current causal history -- which is still
    strictly causal (we only ever see data up to t-1) and reproducible.
    """
    # Causal standardisation (paper does not need this for asset
    # features because returns/vol/momentum are already on similar
    # scales; macro variables span many orders of magnitude).
    mu = X_macro_full.mean(axis=0, keepdims=True)
    sd = X_macro_full.std(axis=0, keepdims=True) + 1e-8
    X_macro_full = (X_macro_full - mu) / sd

    # -- K-selection (uses macro F_K cadence) ---------------------------
    if (step_index % state.f_k) == 0 or step_index == 0:
        if state.monotone_K:
            Ks = list(range(state.K_curr, state.max_regimes + 1))
        else:
            Ks = list(range(state.min_regimes, state.max_regimes + 1))
        K_cand, new_models = select_K_predictive(
            X_macro_full, Ks, L_val=state.l_val,
            lamK=state.lam_k, warm_pool=state.warm_pool,
            return_models=True,
        )
        for k_, mdl in new_models.items():
            state.warm_pool[k_] = mdl
        if state.monotone_K:
            state.K_curr = max(state.K_curr, K_cand)
        else:
            state.K_curr = int(K_cand)

    # -- Refit cadence (warm-started) ----------------------------------
    must_refit = ((state.last_hmm is None)
                  or (step_index % refit_every == 0))
    if must_refit or state.last_hmm.n_components != state.K_curr:
        init = state.last_hmm if (
            state.last_hmm is not None
            and state.last_hmm.n_components == state.K_curr
        ) else state.warm_pool.get(state.K_curr)
        state.last_hmm = fit_gaussian_hmm_warm(
            X_macro_full, K=state.K_curr, init_hmm=init,
        )
    hmm = state.last_hmm

    # -- Filtered posteriors over components ---------------------------
    probsK = hmm.predict_proba(X_macro_full)
    pK = probsK[-1].copy()
    zK_raw = np.argmax(probsK, axis=1)

    # -- Component moments in MACRO FEATURE SPACE for W2 distance ------
    N_feat = X_macro_full.shape[1]
    MU_K  = np.asarray(hmm.means_)                       # (K, N_feat)
    SIG_K = np.asarray(hmm.covars_).reshape(state.K_curr, N_feat, N_feat)

    # -- Template management -----------------------------------------
    if state.MU_tpl is None or state.SIG_tpl is None:
        state.MU_tpl  = MU_K.copy()
        state.SIG_tpl = SIG_K.copy()
    else:
        state.MU_tpl, state.SIG_tpl = spawn_template_if_needed(
            MU_K, SIG_K, pK,
            state.MU_tpl, state.SIG_tpl,
            spawn_thresh=state.spawn_thresh, G_max=state.g_max,
        )

    map_k_to_g, _ = map_components_to_templates(
        MU_K, SIG_K, state.MU_tpl, state.SIG_tpl
    )
    G  = state.MU_tpl.shape[0]
    pG = aggregate_template_posteriors(pK, map_k_to_g, G)

    state.MU_tpl, state.SIG_tpl = update_templates_ema(
        state.MU_tpl, state.SIG_tpl,
        MU_K, SIG_K, pK, map_k_to_g,
        eta=state.eta_tpl, N=N_feat,
    )

    history_g = np.array([map_k_to_g[int(k)] for k in zK_raw], dtype=int)

    return dict(
        K=int(state.K_curr), G=int(G), pG=pG,
        z_argmax=zK_raw, comp_to_tpl=map_k_to_g,
        history_g=history_g,
    )


# =======================================================================
#  Joint conditional moments (g, h) over forward returns
# =======================================================================
def _joint_moments(ret_align: pd.DataFrame,
                   macro_g_s: np.ndarray,
                   market_h_s: np.ndarray,
                   G_macro: int, G_market: int, N: int,
                   min_obs: int = MIN_OBS_PER_REGIME):
    """
    Compute (mu_{g,h}, Sigma_{g,h}) on FORWARD returns r_{s+1} grouped
    by joint label (g_s, h_s).  Falls back to market-only moments
    (mu_h, Sigma_h) for cells with N_{g,h} < min_obs.

    Returns
    -------
    MU    : (G_macro, G_market, N)
    SIG   : (G_macro, G_market, N, N)
    """
    r_fwd = ret_align.shift(-1).dropna()
    # macro_g_s and market_h_s are aligned to ret_align; we need to
    # truncate one row off the front to match r_fwd's length
    L = min(len(macro_g_s), len(market_h_s), len(r_fwd))
    g_arr = macro_g_s[:L]
    h_arr = market_h_s[:L]
    r_fwd = r_fwd.iloc[:L]

    # Market-only moments per h (fallback)
    MU_h  = np.zeros((G_market, N))
    SIG_h = np.tile(np.eye(N), (G_market, 1, 1))
    for h in range(G_market):
        mask = (h_arr == h)
        if mask.sum() >= min_obs:
            rh = r_fwd.iloc[mask]
            MU_h[h]  = rh.mean().values
            SIG_h[h] = LedoitWolf().fit(rh.values).covariance_

    # Joint moments per (g, h)
    MU  = np.zeros((G_macro, G_market, N))
    SIG = np.zeros((G_macro, G_market, N, N))
    for g in range(G_macro):
        for h in range(G_market):
            mask = (g_arr == g) & (h_arr == h)
            if mask.sum() >= min_obs:
                rgh = r_fwd.iloc[mask]
                MU[g, h]  = rgh.mean().values
                SIG[g, h] = LedoitWolf().fit(rgh.values).covariance_
            else:
                MU[g, h]  = MU_h[h]
                SIG[g, h] = SIG_h[h]
    return MU, SIG, MU_h, SIG_h


# =======================================================================
#  Macro tilt on equity-risk expected returns
# =======================================================================
def _macro_sharpe_per_g(ret_align: pd.DataFrame,
                        macro_g_s: np.ndarray,
                        G_macro: int, asset_names: list[str],
                        min_obs: int = MIN_OBS_PER_REGIME) -> np.ndarray:
    """
    For each macro template g, compute the historical Sharpe of every
    asset across all days with that macro template.  Returns an
    (G_macro, N) array.
    """
    r_fwd = ret_align.shift(-1).dropna()
    L = min(len(macro_g_s), len(r_fwd))
    g_arr = macro_g_s[:L]
    r_fwd = r_fwd.iloc[:L]

    S = np.zeros((G_macro, len(asset_names)))
    for g in range(G_macro):
        mask = (g_arr == g)
        if mask.sum() >= min_obs:
            sub = r_fwd.iloc[mask]
            mu  = sub.mean().values
            sd  = sub.std(ddof=1).values + 1e-12
            S[g] = np.sqrt(252) * mu / sd
    return S


def _apply_macro_tilt(mu_t: np.ndarray,
                      pG_macro: np.ndarray,
                      sharpe_per_g: np.ndarray,
                      asset_names: list[str],
                      strength: float,
                      mode: str = "symmetric",
                      tilt_assets: tuple = ("SPX", "OIL")) -> np.ndarray:
    """
    Re-scale risky-asset expected returns by a smooth macro-conditional
    factor.  Two modes:

      'off':         no tilt (identity).
      'symmetric':   factor_i = 1 + alpha * sum_g pG[g] * tanh(S_{g,i})
                     - always positive if S_{g,i} > 0 in all dominant
                       regimes -> can only push mu UP.  Original design.
      'asymmetric':  factor_i = 1 + alpha * sum_g pG[g] *
                                tanh(S_{g,i} - S_bar_i)
                     where S_bar_i = mean_g pG[g] * S_{g,i} is the
                     current macro-mixed unconditional Sharpe.  Now
                     tanh can fire negative when the current macro
                     regime has worse-than-average Sharpe for asset i,
                     enabling DOWNWARD tilt in adverse regimes.

    Defensive assets (BOND, GOLD, USD) are left untouched.

    tanh(.) saturates around +/- 1 to avoid wild scaling on regimes
    with extreme historical Sharpe.
    """
    if mode == "off" or strength <= 0.0:
        return mu_t

    tilt_set = set(tilt_assets)
    factor = np.ones(len(asset_names))

    if mode == "asymmetric":
        # Per-asset unconditional Sharpe mixed by current macro
        # posterior.  Equivalent to sum_g pG[g] * S[g, i].
        S_bar = pG_macro @ sharpe_per_g          # shape (N,)

    for i, a in enumerate(asset_names):
        if a not in tilt_set:
            continue
        s_a = sharpe_per_g[:, i]                 # (G_macro,)
        if mode == "symmetric":
            tilt = float(np.dot(pG_macro, np.tanh(s_a)))
        elif mode == "asymmetric":
            tilt = float(np.dot(pG_macro, np.tanh(s_a - S_bar[i])))
        else:
            raise ValueError(f"Unknown HIER_MACRO_TILT_MODE: {mode}")
        factor[i] = 1.0 + strength * tilt

    return mu_t * factor


# =======================================================================
#  Backtest engine for the hierarchical strategy
# =======================================================================
def run_hierarchical_strategy(returns: pd.DataFrame,
                              macro_levels: pd.DataFrame,
                              oos_start: str = OOS_START,
                              verbose: bool = True) -> dict:
    """
    Hierarchical (macro X market) Wasserstein-HMM allocation.

    Returns a result dict in the same shape as the pure-market engine,
    plus extra diagnostics: macro_label, macro_prob, G_macro_history,
    K_macro_history.
    """
    asset_names = list(returns.columns)
    N = len(asset_names)

    macro_features_full = build_macro_features(macro_levels)

    oos_start_ts = pd.Timestamp(oos_start)
    common_idx = returns.index.intersection(macro_features_full.index)
    common_idx = common_idx[common_idx >= oos_start_ts]
    oos_dates = common_idx

    # Two independent WHMM states
    market_state = WassersteinHMMState(
        n_features_returns=N,
        min_regimes=MIN_REGIMES, max_regimes=MAX_REGIMES, g_max=G_MAX,
        eta_tpl=ETA_TPL, spawn_thresh=SPAWN_THRESH,
        f_k=F_K, l_val=L_VAL, lam_k=LAM_K, monotone_K=True,
    )
    macro_state = WassersteinHMMState(
        n_features_returns=N,        # unused in macro step API
        min_regimes=MACRO_MIN_REGIMES, max_regimes=MACRO_MAX_REGIMES,
        g_max=MACRO_G_MAX, eta_tpl=MACRO_ETA_TPL,
        spawn_thresh=MACRO_SPAWN_THRESH,
        f_k=MACRO_F_K, l_val=L_VAL, lam_k=LAM_K, monotone_K=True,
    )

    # Outputs
    weights_oos, pnl_oos, dates_oos = [], [], []
    macro_label_history, macro_prob_history = [], []
    K_macro_history, G_macro_history = [], []
    market_label_history, market_prob_history = [], []
    K_market_history, G_market_history = [], []

    w_prev = np.zeros(N)

    iterator = enumerate(oos_dates)
    if verbose:
        iterator = enumerate(tqdm(oos_dates, desc="Hierarchical", ncols=100))

    for t_i, date in iterator:
        # ---- 1) Strictly causal histories -------------------------
        ret_hist   = returns.loc[returns.index < date]
        macro_hist = macro_levels.loc[macro_levels.index < date]

        X_asset_df = build_asset_features(ret_hist)
        X_macro_df = build_macro_features(macro_hist)

        common = X_asset_df.index.intersection(X_macro_df.index)
        common = _cap_window(common.to_frame(index=False), MAX_TRAIN_WINDOW).iloc[:, 0]
        common = pd.DatetimeIndex(common)

        if len(common) < max(MIN_HIST_FOR_HMM, 2):
            w_t = w_prev.copy()
            pnl_t = float(np.dot(w_t, returns.loc[date].values))
            weights_oos.append(w_t)
            pnl_oos.append(pnl_t)
            dates_oos.append(date)
            macro_label_history.append(np.nan)
            macro_prob_history.append(np.nan)
            K_macro_history.append(np.nan)
            G_macro_history.append(0)
            market_label_history.append(np.nan)
            market_prob_history.append(np.nan)
            K_market_history.append(np.nan)
            G_market_history.append(0)
            continue

        X_asset   = X_asset_df.loc[common].values
        X_macro   = X_macro_df.loc[common].values
        ret_align = ret_hist.loc[common]

        # ---- 2) Macro step ---------------------------------------
        macro_out = step_macro_whmm(
            macro_state, X_macro,
            step_index=t_i, refit_every=MACRO_F_REFIT,
        )
        pG_macro    = macro_out["pG"]            # (G_macro,)
        G_macro     = macro_out["G"]
        history_g_s = macro_out["history_g"]

        # ---- 3) Market step (shared) -----------------------------
        market_out = step_wasserstein_hmm(
            market_state, X_asset, ret_align,
            step_index=t_i, refit_every=F_REFIT,
        )
        pH_market   = market_out["pG"]           # (G_market,)
        G_market    = market_out["G"]
        zK_market   = market_out["zK_raw"]
        # market template label per historical row
        map_k_to_g_mkt = market_out["map_k_to_g"]
        history_h_s = np.array([map_k_to_g_mkt[int(k)] for k in zK_market],
                               dtype=int)

        # Align macro labels to market labels (both on `common` index)
        # `history_g_s` has length = T (full common history), same as
        # `history_h_s`.  We use entries [:-1] (forward-return alignment).
        g_s = history_g_s[:-1]
        h_s = history_h_s[:-1]

        # ---- 4) Joint conditional moments ------------------------
        MU_gh, SIG_gh, MU_h, SIG_h = _joint_moments(
            ret_align, g_s, h_s,
            G_macro=G_macro, G_market=G_market, N=N,
            min_obs=MIN_OBS_PER_REGIME,
        )

        # ---- 5) Composite probability & template-mixture moments -
        # p(t, g, h) = pG_macro[g] * pH_market[h]
        # mu_t    = sum_{g,h} p(t,g,h) mu_{g,h}
        # Sigma_t = sum_{g,h} p(t,g,h) Sigma_{g,h}
        # Vectorised:
        # p_gh: (G_macro, G_market)
        p_gh = np.outer(pG_macro, pH_market)
        mu_t = np.einsum("gh,ghn->n", p_gh, MU_gh)
        Sigma_t = np.einsum("gh,ghnm->nm", p_gh, SIG_gh)
        Sigma_t = 0.5 * (Sigma_t + Sigma_t.T)

        # ---- 6) Optional macro tilt on risky-asset returns -------
        if HIER_MACRO_TILT_MODE != "off" and HIER_MACRO_TILT_STRENGTH > 0.0:
            S_g = _macro_sharpe_per_g(ret_align, g_s,
                                      G_macro=G_macro,
                                      asset_names=asset_names)
            mu_t = _apply_macro_tilt(
                mu_t, pG_macro, S_g, asset_names,
                strength=HIER_MACRO_TILT_STRENGTH,
                mode=HIER_MACRO_TILT_MODE,
                tilt_assets=HIER_MACRO_TILT_ASSETS,
            )

        # ---- 7) MVO -----------------------------------------------
        w_t = solve_mvo(mu_t, Sigma_t, w_prev,
                        lam=LAM, tc=TC, w_max=W_MAX)

        # ---- 8) Realised PnL -------------------------------------
        pnl_t = float(np.dot(w_t, returns.loc[date].values))

        weights_oos.append(w_t)
        pnl_oos.append(pnl_t)
        dates_oos.append(date)
        macro_label_history.append(int(np.argmax(pG_macro)))
        macro_prob_history.append(float(np.max(pG_macro)))
        K_macro_history.append(macro_out["K"])
        G_macro_history.append(G_macro)
        market_label_history.append(market_out["dominant_template"])
        market_prob_history.append(market_out["p_max"])
        K_market_history.append(market_out["K"])
        G_market_history.append(G_market)

        w_prev = w_t

    # ---- Package -----------------------------------------------------
    # For the universal "regime" column we use the MARKET template
    # (so the daily_backtest_output.csv format is consistent with
    # pure-market).  Macro-specific diagnostics are returned separately.
    base = _package_results(
        dates_oos, pnl_oos, weights_oos, asset_names,
        K_market_history, market_label_history,
        market_prob_history, G_market_history,
    )
    idx = base["pnl"].index
    base["macro_label"]     = pd.Series(macro_label_history, index=idx,
                                        name="dominant_macro")
    base["macro_prob"]      = pd.Series(macro_prob_history, index=idx,
                                        name="max_macro_posterior")
    base["K_macro_history"] = pd.Series(K_macro_history, index=idx,
                                        name="K_macro")
    base["G_macro_history"] = pd.Series(G_macro_history, index=idx,
                                        name="G_macro")
    return base


def make_daily_hier_csv(result: dict, path: str | None = None) -> pd.DataFrame:
    """
    Build the hierarchical daily output, mirroring the pure-market
    schema but including extra macro columns:
        date, w_<asset>, pnl, cum_pnl,
        K_market, regime_market, max_p_market, G_market,
        K_macro, regime_macro, max_p_macro, G_macro,
        turnover
    """
    W = result["weights"]
    df = pd.DataFrame(index=W.index)
    df.index.name = "date"
    for c in W.columns:
        df[f"w_{c}"] = W[c].values
    df["pnl"]            = result["pnl"].values
    df["cum_pnl"]        = result["cum_pnl"].values
    df["K_market"]       = result["K_history"].values
    df["regime_market"]  = result["tpl_label"].values
    df["max_p_market"]   = result["tpl_prob"].values
    df["G_market"]       = result["tpl_count"].values
    df["K_macro"]        = result["K_macro_history"].values
    df["regime_macro"]   = result["macro_label"].values
    df["max_p_macro"]    = result["macro_prob"].values
    df["G_macro"]        = result["G_macro_history"].values
    df["turnover"]       = result["turnover"].values
    if path is not None:
        df.to_csv(path)
    return df
