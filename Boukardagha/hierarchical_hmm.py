"""
hierarchical_hmm.py
===================
Hierarchical Hidden Markov Models for portfolio construction, following
the architecture of Fine, Singer & Tishby (1998).

This module provides TWO hierarchical strategies that sit on top of the
pure-market Wasserstein-HMM (Boukardagha 2026):

    run_hierarchical_strategy_B()   -- "Hierarchical B": macro x market
                                       joint-mixture moments, NO tilt.
    run_hierarchical_strategy_C()   -- "Hierarchical C": macro layer is a
                                       RISK MODULATOR; market layer owns
                                       direction (mu).  Implements the
                                       Fix-C recommendations.

----------------------------------------------------------------------
HIERARCHICAL B  (macro x market joint mixture)
----------------------------------------------------------------------
  Root
   |-- Macro Wasserstein HMM   (macro features; templates in MACRO
   |                            feature space)
   |-- Market Wasserstein HMM  (asset features; SHARED with pure-market)

At each OOS day t:
  1. Macro step -> p_macro(t, g), per-row macro label g_s.
  2. Market step -> p_market(t, h), per-row market label h_s.
  3. Joint conditional moments on (g, h) cells from FORWARD returns,
     with market-only fallback for sparse cells.
  4. p(t, g, h) = p_macro(t, g) * p_market(t, h).
  5. mu_t, Sigma_t = sum_{g,h} p(t,g,h) (mu_{g,h}, Sigma_{g,h}).
  6. MVO with these moments.  (NO expected-return tilt.)

----------------------------------------------------------------------
HIERARCHICAL C  (macro as risk modulator -- Fix C)
----------------------------------------------------------------------
Three design changes vs. B (see analysis_results.md / methodology.md):

  (C1) Macro templates tracked in ASSET-OUTCOME space.  The macro WHMM
       still infers latent states from the 21 macro features, but the
       persistent TEMPLATES are matched and updated using the Gaussian
       of FORWARD ASSET RETURNS conditional on each macro state.  Macro
       regimes are therefore discriminative for allocation by
       construction.

  (C2) No joint-cell fragmentation; macro modulates RISK, not mu.
       The market layer produces exactly the pure-market mu_t and
       Sigma_t.  The macro layer produces a scalar "stress" score in
       [0, 1].  We then:
          gamma_t = gamma * (1 + kappa * stress_t)
          Sigma_t <- Sigma_t * (1 + sigma_scale * stress_t)
       and solve MVO with (mu_t, Sigma_t, gamma_t).  The market layer's
       clean directional signal is preserved.

  (C3) Tempered + prior-blended macro posterior.  The raw macro
       posterior is one-hot ~99% of the time.  We temper it
       (exponent 1/T) and blend with a uniform prior so the stress
       score varies smoothly rather than switching hard.

Both engines reuse the WassersteinHMMState / step_wasserstein_hmm
infrastructure.  Daily refits are warm-started.
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
    # Hierarchical C
    HIER_C_KAPPA_GAMMA, HIER_C_SIGMA_SCALE,
    HIER_C_MACRO_TEMPERATURE, HIER_C_PRIOR_BLEND, HIER_C_STRESS_METRIC,
    TRADING_DAYS,
)
from features import build_asset_features, build_macro_features
from wasserstein_hmm import (
    WassersteinHMMState, step_wasserstein_hmm,
    fit_gaussian_hmm_warm, select_K_predictive,
    map_components_to_templates, aggregate_template_posteriors,
    spawn_template_if_needed, update_templates_ema,
    compute_regime_moments_forward,
)
from mvo import solve_mvo
from backtest import _package_results, _cap_window


# =======================================================================
#  Macro Wasserstein HMM step
# =======================================================================
def step_macro_whmm(state: WassersteinHMMState,
                    X_macro_full: np.ndarray,
                    step_index: int,
                    refit_every: int,
                    ret_align: pd.DataFrame | None = None,
                    track_space: str = "feature"):
    """
    One macro-layer step.

    The latent states are always inferred from the macro features
    X_macro_full.  How the persistent TEMPLATES are tracked depends on
    `track_space`:

      'feature' (Hierarchical B):
          templates are matched/updated in MACRO FEATURE space using the
          HMM's own emission means/covariances.

      'outcome' (Hierarchical C, change C1):
          templates are matched/updated in ASSET FORWARD-RETURN space:
          for each macro component k we compute the Gaussian of forward
          asset returns r_{s+1} on days assigned to component k, and use
          those (mu, Sigma) for the Wasserstein template tracking.  This
          requires `ret_align` (asset returns aligned to the macro
          feature rows) and uses its column count as N.

    Returns dict with: K, G, pG (template posteriors at t), history_g
    (per-row dominant template label), and -- for 'outcome' tracking --
    the per-template forward-asset-return moments MU_tpl, SIG_tpl.

    The macro features have wildly different scales, so we z-score them
    inside this step (strictly causal: only data up to t-1 is seen).
    """
    # Causal standardisation of macro features.
    mu_f = X_macro_full.mean(axis=0, keepdims=True)
    sd_f = X_macro_full.std(axis=0, keepdims=True) + 1e-8
    Xz = (X_macro_full - mu_f) / sd_f

    # -- K-selection (macro F_K cadence) --------------------------------
    if (step_index % state.f_k) == 0 or step_index == 0:
        if state.monotone_K:
            Ks = list(range(state.K_curr, state.max_regimes + 1))
        else:
            Ks = list(range(state.min_regimes, state.max_regimes + 1))
        K_cand, new_models = select_K_predictive(
            Xz, Ks, L_val=state.l_val, lamK=state.lam_k,
            warm_pool=state.warm_pool, return_models=True,
        )
        for k_, mdl in new_models.items():
            state.warm_pool[k_] = mdl
        if state.monotone_K:
            state.K_curr = max(state.K_curr, K_cand)
        else:
            state.K_curr = int(K_cand)

    # -- Refit cadence (warm-started) -----------------------------------
    must_refit = ((state.last_hmm is None)
                  or (step_index % refit_every == 0))
    if must_refit or state.last_hmm.n_components != state.K_curr:
        init = state.last_hmm if (
            state.last_hmm is not None
            and state.last_hmm.n_components == state.K_curr
        ) else state.warm_pool.get(state.K_curr)
        state.last_hmm = fit_gaussian_hmm_warm(Xz, K=state.K_curr, init_hmm=init)
    hmm = state.last_hmm

    # -- Filtered posteriors over components ---------------------------
    probsK = hmm.predict_proba(Xz)
    pK = probsK[-1].copy()
    zK_raw = np.argmax(probsK, axis=1)

    # -- Component moments for template tracking ------------------------
    if track_space == "feature":
        N_track = Xz.shape[1]
        MU_K  = np.asarray(hmm.means_)
        SIG_K = np.asarray(hmm.covars_).reshape(state.K_curr, N_track, N_track)
    elif track_space == "outcome":
        if ret_align is None:
            raise ValueError("track_space='outcome' requires ret_align")
        N_track = ret_align.shape[1]
        z_s = zK_raw[:-1]
        MU_K, SIG_K = compute_regime_moments_forward(
            ret_align, z_s, K=state.K_curr, N=N_track,
            min_obs=MIN_OBS_PER_REGIME,
        )
    else:
        raise ValueError(f"Unknown track_space: {track_space}")

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
        eta=state.eta_tpl, N=N_track,
    )

    history_g = np.array([map_k_to_g[int(k)] for k in zK_raw], dtype=int)

    return dict(
        K=int(state.K_curr), G=int(G), pG=pG,
        z_argmax=zK_raw, comp_to_tpl=map_k_to_g,
        history_g=history_g,
        MU_tpl=state.MU_tpl.copy(), SIG_tpl=state.SIG_tpl.copy(),
    )


# =======================================================================
#  Joint conditional moments (g, h) over forward returns  [Hier B]
# =======================================================================
def _joint_moments(ret_align: pd.DataFrame,
                   macro_g_s: np.ndarray,
                   market_h_s: np.ndarray,
                   G_macro: int, G_market: int, N: int,
                   min_obs: int = MIN_OBS_PER_REGIME):
    """(mu_{g,h}, Sigma_{g,h}) on forward returns grouped by joint label,
    with market-only (mu_h, Sigma_h) fallback for sparse cells."""
    r_fwd = ret_align.shift(-1).dropna()
    L = min(len(macro_g_s), len(market_h_s), len(r_fwd))
    g_arr = macro_g_s[:L]
    h_arr = market_h_s[:L]
    r_fwd = r_fwd.iloc[:L]

    MU_h  = np.zeros((G_market, N))
    SIG_h = np.tile(np.eye(N), (G_market, 1, 1))
    for h in range(G_market):
        mask = (h_arr == h)
        if mask.sum() >= min_obs:
            rh = r_fwd.iloc[mask]
            MU_h[h]  = rh.mean().values
            SIG_h[h] = LedoitWolf().fit(rh.values).covariance_

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
#  Macro stress score per template  [Hier C]
# =======================================================================
def _macro_stress_per_g(ret_align: pd.DataFrame,
                        macro_g_s: np.ndarray,
                        G_macro: int,
                        metric: str = "vol",
                        min_obs: int = MIN_OBS_PER_REGIME) -> np.ndarray:
    """
    Per macro template g, compute a raw turbulence score from forward
    EQUAL-WEIGHT portfolio returns on days in that regime, then
    cross-sectionally normalise to [0, 1].

      'vol'      : annualised vol of the equal-weight forward return.
      'drawdown' : within-regime max drawdown magnitude.
      'sharpe'   : max(0, -annualised Sharpe)  (only downside counts).

    We use the equal-weight portfolio as a regime-agnostic proxy for
    "how turbulent is this macro regime for the cross-section", so the
    stress score does not depend on the market layer's allocation.
    """
    r_fwd = ret_align.shift(-1).dropna()
    L = min(len(macro_g_s), len(r_fwd))
    g_arr = macro_g_s[:L]
    r_fwd = r_fwd.iloc[:L]
    ew = r_fwd.mean(axis=1)   # equal-weight daily forward return

    raw = np.full(G_macro, np.nan)
    for g in range(G_macro):
        mask = (g_arr == g)
        if mask.sum() < min_obs:
            continue
        s = ew.iloc[mask.values] if hasattr(mask, "values") else ew.iloc[mask]
        if metric == "vol":
            raw[g] = float(s.std(ddof=1) * np.sqrt(TRADING_DAYS))
        elif metric == "drawdown":
            wealth = (1.0 + s).cumprod()
            raw[g] = float(-(wealth / wealth.cummax() - 1.0).min())
        elif metric == "sharpe":
            sd = s.std(ddof=1) + 1e-12
            sharpe = np.sqrt(TRADING_DAYS) * s.mean() / sd
            raw[g] = float(max(0.0, -sharpe))
        else:
            raise ValueError(f"Unknown stress metric: {metric}")

    # Cross-sectional min-max normalisation to [0, 1].
    finite = raw[np.isfinite(raw)]
    if len(finite) <= 1:
        out = np.zeros(G_macro)
    else:
        lo, hi = np.nanmin(raw), np.nanmax(raw)
        if hi - lo < 1e-12:
            out = np.zeros(G_macro)
        else:
            out = (raw - lo) / (hi - lo)
    # Undefined (sparse) regimes -> neutral 0 stress.
    out = np.where(np.isfinite(out), out, 0.0)
    return out


def _temper_posterior(pG: np.ndarray,
                      temperature: float,
                      prior_blend: float) -> np.ndarray:
    """
    Soften a (possibly one-hot) posterior:
        p_temp[g] proportional to pG[g] ** (1/temperature)
        p_out = (1 - prior_blend) * p_temp + prior_blend * uniform
    """
    pG = np.asarray(pG, dtype=float)
    G = len(pG)
    if G == 0:
        return pG
    inv_t = 1.0 / max(temperature, 1e-6)
    pt = np.power(np.clip(pG, 1e-12, 1.0), inv_t)
    s = pt.sum()
    pt = pt / s if s > 1e-12 else np.ones(G) / G
    uniform = np.ones(G) / G
    out = (1.0 - prior_blend) * pt + prior_blend * uniform
    out = out / out.sum()
    return out


# =======================================================================
#  Shared state container for the two-layer model
# =======================================================================
def _make_market_state(N: int) -> WassersteinHMMState:
    return WassersteinHMMState(
        n_features_returns=N,
        min_regimes=MIN_REGIMES, max_regimes=MAX_REGIMES, g_max=G_MAX,
        eta_tpl=ETA_TPL, spawn_thresh=SPAWN_THRESH,
        f_k=F_K, l_val=L_VAL, lam_k=LAM_K, monotone_K=True,
    )


def _make_macro_state(N: int) -> WassersteinHMMState:
    return WassersteinHMMState(
        n_features_returns=N,
        min_regimes=MACRO_MIN_REGIMES, max_regimes=MACRO_MAX_REGIMES,
        g_max=MACRO_G_MAX, eta_tpl=MACRO_ETA_TPL,
        spawn_thresh=MACRO_SPAWN_THRESH,
        f_k=MACRO_F_K, l_val=L_VAL, lam_k=LAM_K, monotone_K=True,
    )


def _prep_step(returns, macro_levels, date, t_i):
    """Build causal, calendar-aligned, window-capped feature panels for
    one OOS day.  Returns (X_asset, X_macro, ret_align, common) or None
    if insufficient history."""
    ret_hist   = returns.loc[returns.index < date]
    macro_hist = macro_levels.loc[macro_levels.index < date]

    X_asset_df = build_asset_features(ret_hist)
    X_macro_df = build_macro_features(macro_hist)

    common = X_asset_df.index.intersection(X_macro_df.index)
    if MAX_TRAIN_WINDOW is not None and len(common) > MAX_TRAIN_WINDOW:
        common = common[-MAX_TRAIN_WINDOW:]

    if len(common) < max(MIN_HIST_FOR_HMM, 2):
        return None

    X_asset   = X_asset_df.loc[common].values
    X_macro   = X_macro_df.loc[common].values
    ret_align = ret_hist.loc[common]
    return X_asset, X_macro, ret_align, common


# =======================================================================
#  HIERARCHICAL B  (joint mixture, no tilt)
# =======================================================================
def run_hierarchical_strategy_B(returns: pd.DataFrame,
                                macro_levels: pd.DataFrame,
                                oos_start: str = OOS_START,
                                verbose: bool = True) -> dict:
    asset_names = list(returns.columns)
    N = len(asset_names)

    macro_features_full = build_macro_features(macro_levels)
    common_idx = returns.index.intersection(macro_features_full.index)
    oos_dates = common_idx[common_idx >= pd.Timestamp(oos_start)]

    market_state = _make_market_state(N)
    macro_state  = _make_macro_state(N)

    weights_oos, pnl_oos, dates_oos = [], [], []
    macro_label_history, macro_prob_history = [], []
    K_macro_history, G_macro_history = [], []
    market_label_history, market_prob_history = [], []
    K_market_history, G_market_history = [], []

    w_prev = np.zeros(N)

    iterator = enumerate(oos_dates)
    if verbose:
        iterator = enumerate(tqdm(oos_dates, desc="Hierarchical B", ncols=100))

    for t_i, date in iterator:
        prep = _prep_step(returns, macro_levels, date, t_i)
        if prep is None:
            w_t = w_prev.copy()
            pnl_t = float(np.dot(w_t, returns.loc[date].values))
            weights_oos.append(w_t); pnl_oos.append(pnl_t); dates_oos.append(date)
            macro_label_history.append(np.nan); macro_prob_history.append(np.nan)
            K_macro_history.append(np.nan); G_macro_history.append(0)
            market_label_history.append(np.nan); market_prob_history.append(np.nan)
            K_market_history.append(np.nan); G_market_history.append(0)
            continue
        X_asset, X_macro, ret_align, common = prep

        # Macro step (feature-space tracking)
        macro_out = step_macro_whmm(macro_state, X_macro,
                                    step_index=t_i, refit_every=MACRO_F_REFIT,
                                    track_space="feature")
        pG_macro = macro_out["pG"]; G_macro = macro_out["G"]
        history_g_s = macro_out["history_g"]

        # Market step (shared)
        market_out = step_wasserstein_hmm(market_state, X_asset, ret_align,
                                          step_index=t_i, refit_every=F_REFIT)
        pH_market = market_out["pG"]; G_market = market_out["G"]
        zK_market = market_out["zK_raw"]; map_mkt = market_out["map_k_to_g"]
        history_h_s = np.array([map_mkt[int(k)] for k in zK_market], dtype=int)

        g_s = history_g_s[:-1]
        h_s = history_h_s[:-1]

        # Joint conditional moments
        MU_gh, SIG_gh, MU_h, SIG_h = _joint_moments(
            ret_align, g_s, h_s, G_macro=G_macro, G_market=G_market,
            N=N, min_obs=MIN_OBS_PER_REGIME,
        )

        # Composite probability & mixture moments
        p_gh = np.outer(pG_macro, pH_market)
        mu_t = np.einsum("gh,ghn->n", p_gh, MU_gh)
        Sigma_t = np.einsum("gh,ghnm->nm", p_gh, SIG_gh)
        Sigma_t = 0.5 * (Sigma_t + Sigma_t.T)

        # MVO (NO tilt)
        w_t = solve_mvo(mu_t, Sigma_t, w_prev, lam=LAM, tc=TC, w_max=W_MAX)
        pnl_t = float(np.dot(w_t, returns.loc[date].values))

        weights_oos.append(w_t); pnl_oos.append(pnl_t); dates_oos.append(date)
        macro_label_history.append(int(np.argmax(pG_macro)))
        macro_prob_history.append(float(np.max(pG_macro)))
        K_macro_history.append(macro_out["K"]); G_macro_history.append(G_macro)
        market_label_history.append(market_out["dominant_template"])
        market_prob_history.append(market_out["p_max"])
        K_market_history.append(market_out["K"]); G_market_history.append(G_market)
        w_prev = w_t

    base = _package_results(
        dates_oos, pnl_oos, weights_oos, asset_names,
        K_market_history, market_label_history,
        market_prob_history, G_market_history,
    )
    idx = base["pnl"].index
    base["macro_label"]     = pd.Series(macro_label_history, index=idx, name="dominant_macro")
    base["macro_prob"]      = pd.Series(macro_prob_history, index=idx, name="max_macro_posterior")
    base["K_macro_history"] = pd.Series(K_macro_history, index=idx, name="K_macro")
    base["G_macro_history"] = pd.Series(G_macro_history, index=idx, name="G_macro")
    return base


# =======================================================================
#  HIERARCHICAL C  (macro as risk modulator -- Fix C)
# =======================================================================
def run_hierarchical_strategy_C(returns: pd.DataFrame,
                                macro_levels: pd.DataFrame,
                                oos_start: str = OOS_START,
                                verbose: bool = True) -> dict:
    asset_names = list(returns.columns)
    N = len(asset_names)

    macro_features_full = build_macro_features(macro_levels)
    common_idx = returns.index.intersection(macro_features_full.index)
    oos_dates = common_idx[common_idx >= pd.Timestamp(oos_start)]

    market_state = _make_market_state(N)
    macro_state  = _make_macro_state(N)

    weights_oos, pnl_oos, dates_oos = [], [], []
    macro_label_history, macro_prob_history = [], []
    K_macro_history, G_macro_history = [], []
    market_label_history, market_prob_history = [], []
    K_market_history, G_market_history = [], []
    stress_history, gamma_history = [], []

    w_prev = np.zeros(N)

    iterator = enumerate(oos_dates)
    if verbose:
        iterator = enumerate(tqdm(oos_dates, desc="Hierarchical C", ncols=100))

    for t_i, date in iterator:
        prep = _prep_step(returns, macro_levels, date, t_i)
        if prep is None:
            w_t = w_prev.copy()
            pnl_t = float(np.dot(w_t, returns.loc[date].values))
            weights_oos.append(w_t); pnl_oos.append(pnl_t); dates_oos.append(date)
            macro_label_history.append(np.nan); macro_prob_history.append(np.nan)
            K_macro_history.append(np.nan); G_macro_history.append(0)
            market_label_history.append(np.nan); market_prob_history.append(np.nan)
            K_market_history.append(np.nan); G_market_history.append(0)
            stress_history.append(np.nan); gamma_history.append(LAM)
            continue
        X_asset, X_macro, ret_align, common = prep

        # --- Macro step with ASSET-OUTCOME-SPACE template tracking (C1) ---
        macro_out = step_macro_whmm(macro_state, X_macro,
                                    step_index=t_i, refit_every=MACRO_F_REFIT,
                                    ret_align=ret_align, track_space="outcome")
        pG_macro = macro_out["pG"]; G_macro = macro_out["G"]
        history_g_s = macro_out["history_g"]

        # --- Market step: owns mu_t and Sigma_t exactly as pure-market ---
        market_out = step_wasserstein_hmm(market_state, X_asset, ret_align,
                                          step_index=t_i, refit_every=F_REFIT)
        mu_t    = market_out["mu_t"]
        Sigma_t = market_out["Sigma_t"]
        G_market = market_out["G"]

        # --- Macro stress score, tempered posterior (C3) ---
        g_s = history_g_s[:-1]
        stress_g = _macro_stress_per_g(
            ret_align, g_s, G_macro=G_macro,
            metric=HIER_C_STRESS_METRIC, min_obs=MIN_OBS_PER_REGIME,
        )
        pG_temp = _temper_posterior(
            pG_macro, temperature=HIER_C_MACRO_TEMPERATURE,
            prior_blend=HIER_C_PRIOR_BLEND,
        )
        stress_t = float(np.dot(pG_temp, stress_g))   # in [0, 1]

        # --- Macro modulates RISK only (C2) ---
        gamma_t = LAM * (1.0 + HIER_C_KAPPA_GAMMA * stress_t)
        Sigma_mod = Sigma_t * (1.0 + HIER_C_SIGMA_SCALE * stress_t)
        Sigma_mod = 0.5 * (Sigma_mod + Sigma_mod.T)

        # --- MVO with macro-modulated risk; market owns direction ---
        w_t = solve_mvo(mu_t, Sigma_mod, w_prev,
                        lam=gamma_t, tc=TC, w_max=W_MAX)
        pnl_t = float(np.dot(w_t, returns.loc[date].values))

        weights_oos.append(w_t); pnl_oos.append(pnl_t); dates_oos.append(date)
        macro_label_history.append(int(np.argmax(pG_macro)))
        macro_prob_history.append(float(np.max(pG_macro)))
        K_macro_history.append(macro_out["K"]); G_macro_history.append(G_macro)
        market_label_history.append(market_out["dominant_template"])
        market_prob_history.append(market_out["p_max"])
        K_market_history.append(market_out["K"]); G_market_history.append(G_market)
        stress_history.append(stress_t); gamma_history.append(gamma_t)
        w_prev = w_t

    base = _package_results(
        dates_oos, pnl_oos, weights_oos, asset_names,
        K_market_history, market_label_history,
        market_prob_history, G_market_history,
    )
    idx = base["pnl"].index
    base["macro_label"]     = pd.Series(macro_label_history, index=idx, name="dominant_macro")
    base["macro_prob"]      = pd.Series(macro_prob_history, index=idx, name="max_macro_posterior")
    base["K_macro_history"] = pd.Series(K_macro_history, index=idx, name="K_macro")
    base["G_macro_history"] = pd.Series(G_macro_history, index=idx, name="G_macro")
    base["macro_stress"]    = pd.Series(stress_history, index=idx, name="macro_stress")
    base["gamma_eff"]       = pd.Series(gamma_history, index=idx, name="gamma_eff")
    return base


# =======================================================================
#  Unified daily-output CSV writers
# =======================================================================
def make_daily_hier_csv(result: dict, path: str | None = None) -> pd.DataFrame:
    """Daily output for Hierarchical B (market + macro columns)."""
    W = result["weights"]
    df = pd.DataFrame(index=W.index); df.index.name = "date"
    for c in W.columns:
        df[f"w_{c}"] = W[c].values
    df["pnl"]           = result["pnl"].values
    df["cum_pnl"]       = result["cum_pnl"].values
    df["K_market"]      = result["K_history"].values
    df["regime_market"] = result["tpl_label"].values
    df["max_p_market"]  = result["tpl_prob"].values
    df["G_market"]      = result["tpl_count"].values
    df["K_macro"]       = result["K_macro_history"].values
    df["regime_macro"]  = result["macro_label"].values
    df["max_p_macro"]   = result["macro_prob"].values
    df["G_macro"]       = result["G_macro_history"].values
    df["turnover"]      = result["turnover"].values
    if path is not None:
        df.to_csv(path)
    return df


def make_daily_hierC_csv(result: dict, path: str | None = None) -> pd.DataFrame:
    """Daily output for Hierarchical C (adds macro_stress and gamma_eff)."""
    df = make_daily_hier_csv(result, path=None)
    df["macro_stress"] = result["macro_stress"].values
    df["gamma_eff"]    = result["gamma_eff"].values
    if path is not None:
        df.to_csv(path)
    return df


# Backwards-compatible alias: the previous single hierarchical entry
# point now maps to Hierarchical B (no tilt).
run_hierarchical_strategy = run_hierarchical_strategy_B
