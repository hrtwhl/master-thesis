"""
wasserstein_hmm.py
------------------
Wasserstein-Hidden-Markov-Model machinery for daily regime inference.

Section references match Boukardagha (2026):
    §5.2 — Predictive model-order selection      (`select_K_predictive`)
    §5.3 — Rolling Gaussian HMM inference        (`fit_hmm`)
    §5.4 — 2-Wasserstein template tracking       (`w2_gaussian`,
                                                  `map_components_to_templates`,
                                                  `spawn_template_if_needed`,
                                                  `update_templates_ema`)
    §7.1 — Strictly-causal regime moments via
            forward returns                      (`compute_regime_moments_forward`)
"""

from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd
from hmmlearn.hmm import GaussianHMM
from sklearn.covariance import LedoitWolf

from config import HMM, RUN

# hmmlearn emits its "Model is not converging" (and a few other ill-conditioning
# notices) at WARNING level via `logging.warning`. Suppressing them at ERROR
# level is clean; we already retry the fit with stronger regularization and
# fall back to the cached HMM on persistent failure, so these warnings carry
# no actionable information for the end user.
logging.getLogger("hmmlearn").setLevel(logging.ERROR)
logging.getLogger("hmmlearn.base").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=RuntimeWarning, module="hmmlearn")


# --------------------------------------------------------------------- #
# 1.  Numerical helpers for the W2 distance between Gaussians
# --------------------------------------------------------------------- #
def _symmetrize(A: np.ndarray) -> np.ndarray:
    """½(A + Aᵀ) — kills any numerical asymmetry from float arithmetic."""
    return 0.5 * (A + A.T)


def _sqrtm_psd(A: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Symmetric positive-(semi)definite matrix square root via eigendecomp.

    Cheaper and more stable than `scipy.linalg.sqrtm` for the small (5x5)
    covariances we use.
    """
    A = _symmetrize(A)
    w, V = np.linalg.eigh(A)
    w = np.maximum(w, eps)
    return (V * np.sqrt(w)) @ V.T


def w2_gaussian(
    mu1: np.ndarray, S1: np.ndarray,
    mu2: np.ndarray, S2: np.ndarray,
    S2_sqrt: np.ndarray | None = None,
    eps: float = 1e-12,
) -> float:
    """2-Wasserstein distance between two N-variate Gaussians.

    W2² = ‖µ1-µ2‖² + Tr( S1 + S2 − 2 (S2^½ S1 S2^½)^½ )

    Pass `S2_sqrt` precomputed for speed when distances to one fixed reference
    distribution (e.g. a template) are evaluated in a hot loop.
    """
    mu1 = np.asarray(mu1).ravel()
    mu2 = np.asarray(mu2).ravel()
    S1  = _symmetrize(np.asarray(S1))
    S2  = _symmetrize(np.asarray(S2))

    dm2 = float(((mu1 - mu2) ** 2).sum())

    if S2_sqrt is None:
        S2_sqrt = _sqrtm_psd(S2, eps=eps)

    M = S2_sqrt @ S1 @ S2_sqrt
    tr_term = float(np.trace(S1 + S2 - 2.0 * _sqrtm_psd(M, eps=eps)))
    tr_term = max(tr_term, 0.0)  # numerical floor
    return float(np.sqrt(dm2 + tr_term))


# --------------------------------------------------------------------- #
# 2.  Predictive model-order selection
# --------------------------------------------------------------------- #
# Progressive `min_covar` ladder used by both `_fit_and_score` (K-selection)
# and `fit_hmm` (main fit) to recover from rare non-PSD covariance updates
# during EM. Each rung adds more diagonal regularization.
_MIN_COVAR_LADDER: tuple[float, ...] = (1e-3, 1e-2, 1e-1)

# Seed offsets tried at each `min_covar` rung in `fit_hmm`. Different EM
# initializations very often resolve the rare degenerate-regime case that
# produces non-PSD covariances.
_SEED_OFFSETS: tuple[int, ...] = (0, 1, 7, 23)


def _fit_and_score(
    K: int, X_fit: np.ndarray, X_val: np.ndarray,
    lam_k: float, n_iter: int, random_state: int,
) -> tuple[int, float]:
    """Helper for joblib: fit a K-state HMM on X_fit, score on X_val.

    Uses the same `min_covar` regularization ladder as `fit_hmm` so a
    candidate K isn't unfairly knocked out of the selection by a transient
    numerical hiccup.
    """
    for mc in _MIN_COVAR_LADDER:
        try:
            hmm = GaussianHMM(
                n_components=int(K),
                covariance_type="full",
                n_iter=n_iter,
                random_state=random_state,
                min_covar=float(mc),
            )
            hmm.fit(X_fit)
            predll = float(hmm.score(X_val))
            return int(K), predll - lam_k * float(K)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
    return int(K), -np.inf


def select_K_predictive(
    X_all: np.ndarray,
    K_candidates: list[int],
    l_val: int = HMM.l_val,
    n_iter: int = HMM.n_iter,
    random_state: int = HMM.random_state,
    lam_k: float = HMM.lam_k,
    parallel: bool = RUN.parallel_k,
    n_jobs: int = RUN.n_jobs,
) -> int:
    """Pick K maximizing  PredLL(K) - λ_K · K  on a validation slice.

    Fully strictly-causal: validation is the tail of the history *within*
    `X_all`, so no information from after the OOS date leaks in.
    """
    if len(X_all) <= (l_val + 50):
        return int(min(K_candidates))

    X_fit = X_all[:-l_val]
    X_val = X_all[-l_val:]

    if parallel and len(K_candidates) > 1:
        from joblib import Parallel, delayed
        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_fit_and_score)(K, X_fit, X_val, lam_k, n_iter, random_state)
            for K in K_candidates
        )
    else:
        results = [_fit_and_score(K, X_fit, X_val, lam_k, n_iter, random_state)
                   for K in K_candidates]

    best_K, best_score = min(K_candidates), -np.inf
    for K, score in results:
        if score > best_score:
            best_score = score
            best_K = K
    return int(best_K)


# --------------------------------------------------------------------- #
# 3.  HMM fitting + regime moments (forward returns)
# --------------------------------------------------------------------- #
def fit_hmm(
    X: np.ndarray,
    K: int,
    n_iter: int = HMM.n_iter,
    random_state: int = HMM.random_state,
) -> "GaussianHMM | None":
    """Fit a K-state Gaussian HMM (full covariance) with robust retry.

    On rare occasions hmmlearn's EM produces a non-positive-definite
    emission covariance — typically when one regime collects too few
    observations during an iteration. The next `predict_proba` call then
    raises ``ValueError: 'covars' must be symmetric, positive-definite``,
    crashing the backtest mid-run.

    This function tries the fit with progressively stronger diagonal
    regularization (`min_covar`) and a small ladder of alternative random
    seeds. On every attempt it also validates that `predict_proba` works
    on a small slice of `X` before returning. If every attempt fails it
    returns ``None``; the caller is expected to fall back to a cached HMM.
    """
    last_err: Exception | None = None
    for mc in _MIN_COVAR_LADDER:
        for seed_off in _SEED_OFFSETS:
            try:
                hmm = GaussianHMM(
                    n_components=K,
                    covariance_type="full",
                    n_iter=n_iter,
                    random_state=int(random_state) + int(seed_off),
                    min_covar=float(mc),
                )
                hmm.fit(X)
                # Validate: predict_proba is the actual failure surface.
                _ = hmm.predict_proba(X[-min(50, len(X)):])
                return hmm
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as e:
                last_err = e
                continue

    # Total failure — caller will keep using the previously cached HMM.
    return None


def compute_regime_moments_forward(
    ret_align: pd.DataFrame,
    regimes_s: np.ndarray,
    K: int,
    n_assets: int,
    min_obs: int = HMM.min_obs,
) -> tuple[np.ndarray, np.ndarray]:
    """Conditional mean and Ledoit–Wolf-shrunk covariance per regime,
    using **forward** returns r_{s+1} grouped by the regime label z_s.

    Strictly causal: forward labels at the last feature row are dropped.

    Returns
    -------
    MU  : (K, N) array, row k = E[r_{s+1} | z_s = k]
    SIG : (K, N, N) array, slab k = Cov(r_{s+1} | z_s = k)
    """
    r_fwd = ret_align.shift(-1).dropna()
    r_fwd = r_fwd.iloc[: len(regimes_s)]

    MU  = np.zeros((K, n_assets))
    SIG = np.zeros((K, n_assets, n_assets))

    r_arr = r_fwd.to_numpy()
    for k in range(K):
        mask = (regimes_s == k)
        if mask.sum() < min_obs:
            MU[k]      = 0.0
            SIG[k]     = np.eye(n_assets)
        else:
            rk = r_arr[mask]
            MU[k]  = rk.mean(axis=0)
            SIG[k] = LedoitWolf().fit(rk).covariance_
    return MU, SIG


# --------------------------------------------------------------------- #
# 4.  Template tracking
# --------------------------------------------------------------------- #
def map_components_to_templates(
    MU_K: np.ndarray, SIG_K: np.ndarray,
    MU_tpl: np.ndarray, SIG_tpl: np.ndarray,
) -> tuple[list[int], list[float]]:
    """For each HMM component k, find the closest template g under W2.

    The square roots of template covariances are computed *once* per call
    (rather than once per component-template pair), which roughly halves
    the cost of the inner loop.
    """
    K = MU_K.shape[0]
    G = MU_tpl.shape[0]

    # Pre-compute template Σ^{1/2} once
    SIG_tpl_sqrt = [_sqrtm_psd(SIG_tpl[g]) for g in range(G)]

    map_k_to_g: list[int] = []
    dist_k:     list[float] = []
    for k in range(K):
        dists = [
            w2_gaussian(MU_tpl[g], SIG_tpl[g],
                        MU_K[k],   SIG_K[k],
                        S2_sqrt=SIG_tpl_sqrt[g])
            for g in range(G)
        ]
        g_best = int(np.argmin(dists))
        map_k_to_g.append(g_best)
        dist_k.append(float(dists[g_best]))
    return map_k_to_g, dist_k


def spawn_template_if_needed(
    MU_K: np.ndarray, SIG_K: np.ndarray, pK: np.ndarray,
    MU_tpl: np.ndarray | None, SIG_tpl: np.ndarray | None,
    spawn_thresh: float = HMM.spawn_thresh,
    g_max:        int   = HMM.g_max,
) -> tuple[np.ndarray, np.ndarray]:
    """If any component is far (W2 > thresh) from all templates and we still
    have room (G < g_max), create a new template from the highest-mass
    far component.
    """
    if MU_tpl is None:
        return MU_K.copy(), SIG_K.copy()

    G = MU_tpl.shape[0]
    if G >= g_max:
        return MU_tpl, SIG_tpl

    _, dist_k = map_components_to_templates(MU_K, SIG_K, MU_tpl, SIG_tpl)
    dist_k = np.asarray(dist_k)
    candidates = np.where(dist_k > spawn_thresh)[0]
    if candidates.size == 0:
        return MU_tpl, SIG_tpl

    k_new = int(candidates[np.argmax(np.asarray(pK)[candidates])])
    return (np.vstack([MU_tpl, MU_K[k_new][None, :]]),
            np.concatenate([SIG_tpl, SIG_K[k_new][None, :, :]], axis=0))


def update_templates_ema(
    MU_tpl: np.ndarray, SIG_tpl: np.ndarray,
    MU_K:   np.ndarray, SIG_K:   np.ndarray,
    pK:     np.ndarray, map_k_to_g: list[int],
    eta:    float = HMM.eta_tpl,
) -> tuple[np.ndarray, np.ndarray]:
    """Exponential-moving-average update of templates with the posterior-weighted
    average of the HMM components that map to each one.
    """
    MU_tpl  = MU_tpl.copy()
    SIG_tpl = SIG_tpl.copy()
    G, N    = MU_tpl.shape
    pK      = np.asarray(pK).ravel()

    eye_N = np.eye(N)
    for g in range(G):
        ks = [k for k, gk in enumerate(map_k_to_g) if gk == g]
        if not ks:
            continue
        w = pK[ks]
        ws = float(w.sum())
        if ws <= 1e-12:
            continue
        w = w / ws

        mu_bar = (w[:, None] * MU_K[ks]).sum(axis=0)
        # einsum is fast and clear: weighted sum of (N,N) slabs
        S_bar  = np.einsum("k,kij->ij", w, SIG_K[ks])

        MU_tpl[g]  = (1 - eta) * MU_tpl[g]  + eta * mu_bar
        SIG_tpl[g] = _symmetrize((1 - eta) * SIG_tpl[g] + eta * S_bar) + 1e-8 * eye_N
    return MU_tpl, SIG_tpl


def aggregate_template_posteriors(
    pK: np.ndarray, map_k_to_g: list[int], G: int,
) -> np.ndarray:
    """p_t,g = Σ_{k : g(k)=g} p_t,k    (then normalized)."""
    pK = np.asarray(pK).ravel()
    pG = np.zeros(G)
    for k, g in enumerate(map_k_to_g):
        pG[g] += pK[k]
    s = pG.sum()
    if s > 1e-12:
        pG /= s
    return pG
