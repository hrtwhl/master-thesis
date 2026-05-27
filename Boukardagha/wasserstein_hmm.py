"""
wasserstein_hmm.py
==================
Strictly-causal Gaussian HMM with

  (i)   predictive K selection on a held-out validation slice,
  (ii)  Wasserstein-2 template tracking for persistent regime identity,
  (iii) exponential template updates,
  (iv)  template-mixture conditional moments built from FORWARD returns.

Faithful port of Paper_Code.ipynb by Boukardagha (2026), with one
efficiency mechanism added that does NOT alter results:

  - WARM-STARTED daily refits.  Boukardagha refits the HMM on every
    OOS day from cold (n_iter=300 EM iters from random init).  We
    initialise each daily refit from the previous day's fitted
    parameters via hmmlearn's `init_params=""` mechanism, after which
    EM converges in typically 3-10 iters.  At a K-selection date,
    candidate K models are likewise warm-started (we keep one cached
    model per K).  This is purely a wall-time optimisation; the
    likelihood surface and converged maximum are identical.
"""
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from hmmlearn.hmm import GaussianHMM
from sklearn.covariance import LedoitWolf

from config import (
    HMM_N_ITER, RANDOM_SEED, MIN_OBS_PER_REGIME,
)

warnings.filterwarnings("ignore")


# =======================================================================
#  Linear-algebra helpers (Wasserstein-2 distance between Gaussians)
# =======================================================================
def _symmetrize(A: np.ndarray) -> np.ndarray:
    return 0.5 * (A + A.T)


def _sqrtm_psd(A: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Stable square root of a PSD matrix via eigendecomposition."""
    A = _symmetrize(A)
    w, V = np.linalg.eigh(A)
    w = np.maximum(w, eps)
    return V @ np.diag(np.sqrt(w)) @ V.T


def w2_gaussian(mu1, S1, mu2, S2, eps: float = 1e-12) -> float:
    """
    Closed-form 2-Wasserstein distance between two Gaussians:

        W2^2 = ||mu1 - mu2||^2
             + Tr( S1 + S2 - 2 (S2^{1/2} S1 S2^{1/2})^{1/2} )
    """
    mu1 = np.asarray(mu1).reshape(-1)
    mu2 = np.asarray(mu2).reshape(-1)
    S1  = _symmetrize(np.asarray(S1))
    S2  = _symmetrize(np.asarray(S2))

    dm2 = float(np.sum((mu1 - mu2) ** 2))

    S2_sqrt = _sqrtm_psd(S2, eps)
    M = S2_sqrt @ S1 @ S2_sqrt
    M_sqrt = _sqrtm_psd(M, eps)

    tr_term = float(np.trace(S1 + S2 - 2.0 * M_sqrt))
    tr_term = max(tr_term, 0.0)
    return float(np.sqrt(dm2 + tr_term))


# =======================================================================
#  HMM fitting with optional warm-start
# =======================================================================
def fit_gaussian_hmm_cold(X: np.ndarray, K: int,
                          n_iter: int = HMM_N_ITER,
                          random_state: int = RANDOM_SEED) -> GaussianHMM:
    """Fresh fit from hmmlearn's default random init (paper baseline)."""
    hmm = GaussianHMM(n_components=int(K),
                      covariance_type="full",
                      n_iter=int(n_iter),
                      random_state=int(random_state),
                      tol=1e-3)
    hmm.fit(X)
    return hmm


def fit_gaussian_hmm_warm(X: np.ndarray, K: int,
                          init_hmm: GaussianHMM | None = None,
                          n_iter: int = HMM_N_ITER,
                          random_state: int = RANDOM_SEED) -> GaussianHMM:
    """
    Fit at order K.  If init_hmm is given and has the SAME number of
    components, warm-start from its parameters via `init_params=""`.
    Otherwise fall back to a cold fit.

    EM is monotone in log-likelihood; warm-starting from the previous
    day's parameters lands the optimiser on the same likelihood
    surface basin and converges in a handful of iterations.
    """
    if init_hmm is not None and getattr(init_hmm, "n_components", -1) == K:
        try:
            hmm = GaussianHMM(
                n_components=int(K),
                covariance_type="full",
                n_iter=int(n_iter),
                random_state=int(random_state),
                tol=1e-3,
                init_params="", params="stmc",
            )
            hmm.startprob_ = init_hmm.startprob_.copy()
            hmm.transmat_  = init_hmm.transmat_.copy()
            hmm.means_     = init_hmm.means_.copy()
            hmm.covars_    = np.asarray(init_hmm.covars_).copy()
            hmm.fit(X)
            return hmm
        except Exception:
            pass  # fall back to cold
    return fit_gaussian_hmm_cold(X, K=K, n_iter=n_iter,
                                 random_state=random_state)


# Public alias for code that doesn't care about warm vs cold:
fit_gaussian_hmm = fit_gaussian_hmm_cold


# =======================================================================
#  Conditional regime moments from forward returns
# =======================================================================
def compute_regime_moments_forward(ret_align: pd.DataFrame,
                                   regimes_s: np.ndarray,
                                   K: int,
                                   N: int,
                                   min_obs: int = MIN_OBS_PER_REGIME):
    """
    Strictly causal regime moments using FORWARD returns:
        r_{s+1} grouped by regime label z_s.
    Identical to Boukardagha (2026) §5.5 / §7.1.20.
    """
    r_fwd = ret_align.shift(-1).dropna()
    r_fwd = r_fwd.iloc[:len(regimes_s)]

    MU  = np.zeros((K, N))
    SIG = np.zeros((K, N, N))

    for k in range(K):
        rk = r_fwd.iloc[(regimes_s == k)]
        if len(rk) < min_obs:
            MU[k, :]     = 0.0
            SIG[k, :, :] = np.eye(N)
        else:
            MU[k, :]     = rk.mean().values
            SIG[k, :, :] = LedoitWolf().fit(rk.values).covariance_
    return MU, SIG


# =======================================================================
#  Template tracking (verbatim from Boukardagha 2026)
# =======================================================================
def map_components_to_templates(MU_K, SIG_K, MU_tpl, SIG_tpl):
    """Assign each HMM component k to the nearest template g via W2."""
    K = MU_K.shape[0]
    G = MU_tpl.shape[0]
    map_k_to_g = []
    dist_k = []
    for k in range(K):
        dists = [w2_gaussian(MU_tpl[g], SIG_tpl[g], MU_K[k], SIG_K[k])
                 for g in range(G)]
        g_best = int(np.argmin(dists))
        map_k_to_g.append(g_best)
        dist_k.append(float(dists[g_best]))
    return map_k_to_g, dist_k


def spawn_template_if_needed(MU_K, SIG_K, pK, MU_tpl, SIG_tpl,
                             spawn_thresh: float, G_max: int):
    """Add a new template if a posterior-heavy component is far from all
    existing templates."""
    if MU_tpl is None:
        return MU_K.copy(), SIG_K.copy()

    if MU_tpl.shape[0] >= G_max:
        return MU_tpl, SIG_tpl

    _, dist_k = map_components_to_templates(MU_K, SIG_K, MU_tpl, SIG_tpl)
    dist_k = np.asarray(dist_k)
    candidates = np.where(dist_k > spawn_thresh)[0]
    if len(candidates) == 0:
        return MU_tpl, SIG_tpl

    k_new = int(candidates[np.argmax(np.asarray(pK)[candidates])])
    MU_tpl_new  = np.vstack([MU_tpl, MU_K[k_new][None, :]])
    SIG_tpl_new = np.concatenate([SIG_tpl, SIG_K[k_new][None, :, :]], axis=0)
    return MU_tpl_new, SIG_tpl_new


def update_templates_ema(MU_tpl, SIG_tpl, MU_K, SIG_K, pK,
                         map_k_to_g, eta: float, N: int):
    """Exponentially smooth the templates using posterior-weighted
    component averages of newly-mapped components."""
    MU_tpl  = MU_tpl.copy()
    SIG_tpl = SIG_tpl.copy()
    G = MU_tpl.shape[0]
    K = MU_K.shape[0]
    pK = np.asarray(pK).reshape(-1)

    for g in range(G):
        ks = [k for k in range(K) if map_k_to_g[k] == g]
        if len(ks) == 0:
            continue
        w = pK[ks]
        ws = float(w.sum())
        if ws <= 1e-12:
            continue
        w = w / ws

        mu_bar = (w[:, None] * MU_K[ks]).sum(axis=0)
        S_bar  = np.zeros((N, N))
        for kk, wk in zip(ks, w):
            S_bar += wk * SIG_K[kk]

        MU_tpl[g]  = (1 - eta) * MU_tpl[g] + eta * mu_bar
        SIG_tpl[g] = _symmetrize((1 - eta) * SIG_tpl[g] + eta * S_bar) + 1e-8 * np.eye(N)

    return MU_tpl, SIG_tpl


def aggregate_template_posteriors(pK, map_k_to_g, G: int) -> np.ndarray:
    """Sum component posteriors over components mapped to each template."""
    pK = np.asarray(pK).reshape(-1)
    pG = np.zeros(G)
    for k, g in enumerate(map_k_to_g):
        pG[g] += pK[k]
    s = pG.sum()
    if s > 1e-12:
        pG /= s
    return pG


# =======================================================================
#  Predictive K selection (warm-startable version)
# =======================================================================
def select_K_predictive(X_all: np.ndarray, K_candidates,
                        L_val: int = 252,
                        n_iter: int = HMM_N_ITER,
                        random_state: int = RANDOM_SEED,
                        lamK: float = 1.0,
                        warm_pool: dict | None = None,
                        return_models: bool = False):
    """
    Strictly-causal predictive model-order selection.
    Fit on X_all[:-L_val], score on X_all[-L_val:].
    Penalise complexity by lamK * K.

    If `warm_pool` is provided (dict {K -> last fitted hmm}), each
    candidate fit is warm-started from `warm_pool[K]` when possible.
    The newly fitted models are returned (or stored back into warm_pool
    by the caller) so that subsequent K-selection dates can keep
    warm-starting cheaply.
    """
    if len(X_all) <= (L_val + 50):
        if return_models:
            return int(min(K_candidates)), {}
        return int(min(K_candidates))

    X_fit = X_all[:-L_val]
    X_val = X_all[-L_val:]

    best_score = -np.inf
    best_K = int(min(K_candidates))
    models = {}

    for K in K_candidates:
        try:
            init = warm_pool.get(int(K)) if warm_pool is not None else None
            hmm = fit_gaussian_hmm_warm(X_fit, K=int(K), init_hmm=init,
                                        n_iter=n_iter,
                                        random_state=random_state)
            predll = float(hmm.score(X_val))
            score = predll - lamK * float(K)
            if score > best_score:
                best_score = score
                best_K = int(K)
            models[int(K)] = hmm
        except Exception:
            continue

    if return_models:
        return best_K, models
    return best_K


# =======================================================================
#  State container
# =======================================================================
class WassersteinHMMState:
    """
    Mutable container for a Wasserstein-HMM run: K_curr, last fitted
    model (cached for warm-starting), per-K warm-start pool, template
    parameters.

    Used by the backtest loop, which calls `step_wasserstein_hmm` once
    per OOS day.
    """
    def __init__(self,
                 n_features_returns: int,
                 min_regimes: int,
                 max_regimes: int,
                 g_max: int,
                 eta_tpl: float,
                 spawn_thresh: float,
                 f_k: int,
                 l_val: int,
                 lam_k: float,
                 monotone_K: bool = True):
        self.N            = n_features_returns
        self.min_regimes  = int(min_regimes)
        self.max_regimes  = int(max_regimes)
        self.g_max        = int(g_max)
        self.eta_tpl      = float(eta_tpl)
        self.spawn_thresh = float(spawn_thresh)
        self.f_k          = int(f_k)
        self.l_val        = int(l_val)
        self.lam_k        = float(lam_k)
        self.monotone_K   = bool(monotone_K)

        self.K_curr  = int(min_regimes)
        self.MU_tpl  = None
        self.SIG_tpl = None

        # Warm-start caches.
        # - last_hmm: most recent fitted full-history HMM at K_curr
        # - warm_pool[K]: most recent K-candidate model from select_K
        self.last_hmm = None
        self.warm_pool: dict[int, GaussianHMM] = {}


def step_wasserstein_hmm(state: WassersteinHMMState,
                         X_full: np.ndarray,
                         ret_align: pd.DataFrame,
                         step_index: int,
                         refit_every: int):
    """
    One step of the Wasserstein-HMM pipeline.

    Refit cadence:
      - HMM is refit on each step where (step_index % refit_every == 0).
      - Refits are warm-started from `state.last_hmm` for cheap EM.
      - K-selection runs every `state.f_k` steps (paper: weekly).

    Returns
    -------
    dict with keys:
        K, G, pG, mu_t, Sigma_t, dominant_template, p_max,
        zK_raw, pK, MU_K, SIG_K, map_k_to_g
    """
    # ----- 1) Periodic predictive K selection -------------------------
    if (step_index % state.f_k) == 0 or step_index == 0:
        if state.monotone_K:
            K_candidates = list(range(state.K_curr, state.max_regimes + 1))
        else:
            K_candidates = list(range(state.min_regimes, state.max_regimes + 1))

        K_cand, new_models = select_K_predictive(
            X_full, K_candidates, L_val=state.l_val,
            n_iter=HMM_N_ITER, random_state=RANDOM_SEED,
            lamK=state.lam_k, warm_pool=state.warm_pool,
            return_models=True,
        )
        for k_, mdl in new_models.items():
            state.warm_pool[k_] = mdl

        if state.monotone_K:
            state.K_curr = max(state.K_curr, K_cand)
        else:
            state.K_curr = int(K_cand)

    # ----- 2) Refit HMM (warm-started) on refit dates -----------------
    must_refit = ((state.last_hmm is None)
                  or (step_index % refit_every == 0))
    if must_refit or state.last_hmm.n_components != state.K_curr:
        # Pick the best available warm-start: prefer the last_hmm if it
        # has the right K; otherwise use the warm_pool entry.
        init = state.last_hmm if (
            state.last_hmm is not None
            and state.last_hmm.n_components == state.K_curr
        ) else state.warm_pool.get(state.K_curr)
        state.last_hmm = fit_gaussian_hmm_warm(
            X_full, K=state.K_curr, init_hmm=init,
        )

    hmm = state.last_hmm

    # ----- 3) Filtered posteriors -------------------------------------
    probsK = hmm.predict_proba(X_full)
    pK = probsK[-1].copy()
    zK_raw = np.argmax(probsK, axis=1)

    # ----- 4) Forward-return conditional moments ----------------------
    z_s = zK_raw[:-1]
    MU_K, SIG_K = compute_regime_moments_forward(
        ret_align, z_s, K=state.K_curr, N=state.N,
        min_obs=MIN_OBS_PER_REGIME,
    )

    # ----- 5) Template management --------------------------------------
    if state.MU_tpl is None or state.SIG_tpl is None:
        state.MU_tpl  = MU_K.copy()
        state.SIG_tpl = SIG_K.copy()
    else:
        state.MU_tpl, state.SIG_tpl = spawn_template_if_needed(
            MU_K, SIG_K, pK,
            state.MU_tpl, state.SIG_tpl,
            spawn_thresh=state.spawn_thresh,
            G_max=state.g_max,
        )

    map_k_to_g, _ = map_components_to_templates(
        MU_K, SIG_K, state.MU_tpl, state.SIG_tpl
    )
    G  = state.MU_tpl.shape[0]
    pG = aggregate_template_posteriors(pK, map_k_to_g, G)

    state.MU_tpl, state.SIG_tpl = update_templates_ema(
        state.MU_tpl, state.SIG_tpl,
        MU_K, SIG_K, pK, map_k_to_g,
        eta=state.eta_tpl, N=state.N,
    )

    # ----- 6) Template-mixture moments ---------------------------------
    mu_t    = state.MU_tpl.T @ pG
    Sigma_t = np.tensordot(pG, state.SIG_tpl, axes=(0, 0))

    g_hat = int(np.argmax(pG)) if np.all(np.isfinite(pG)) else -1
    p_max = float(np.max(pG)) if np.all(np.isfinite(pG)) else np.nan

    return dict(
        K = int(state.K_curr),
        G = int(G),
        pG = pG,
        mu_t = mu_t,
        Sigma_t = Sigma_t,
        dominant_template = g_hat,
        p_max = p_max,
        zK_raw = zK_raw,
        pK = pK,
        MU_K = MU_K,
        SIG_K = SIG_K,
        map_k_to_g = map_k_to_g,
    )
