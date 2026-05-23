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


def _compute_lam_k(
    n_val: int,
    d_feat: int,
    lam_k_scale: float = HMM.lam_k_scale,
    lam_k_override: float = HMM.lam_k_override,
) -> float:
    """Compute the per-state complexity penalty λ_K used in K-selection.

    Two modes:
      * If `lam_k_override > 0`, return it directly (paper-style flat penalty;
        set to 1.0 for the paper's exact value).
      * Else: ``λ_K = lam_k_scale · 0.5 · log(n_val) · d_feat``  (BIC-style).

    Rationale: with a 35-year history the validation log-likelihood is on the
    order of 10^5, against which a flat λ_K=1 (the paper's value) cannot
    discriminate between candidate K. Scaling with `0.5·log(T)` (BIC) and
    the feature dimensionality gives a penalty that grows correctly with
    sample size and model dimension.
    """
    if lam_k_override > 0:
        return float(lam_k_override)
    return float(lam_k_scale) * 0.5 * float(np.log(max(n_val, 2))) * float(d_feat)


def select_K_predictive(
    X_all: np.ndarray,
    K_candidates: list[int],
    l_val: int = HMM.l_val,
    n_iter: int = HMM.n_iter,
    random_state: int = HMM.random_state,
    lam_k: float | None = None,
    parallel: bool = RUN.parallel_k,
    n_jobs: int = RUN.n_jobs,
) -> int:
    """Pick K maximizing  PredLL(K) - λ_K · K  on a validation slice.

    Fully strictly-causal: validation is the tail of the history *within*
    `X_all`, so no information from after the OOS date leaks in.

    The penalty `λ_K` is computed from the BIC-style formula in
    `_compute_lam_k` unless an explicit value is passed in. The default
    formula adapts to both validation length and feature dimension, which
    matters on multi-decade samples where a flat penalty becomes invisible.
    """
    if len(X_all) <= (l_val + 50):
        return int(min(K_candidates))

    X_fit = X_all[:-l_val]
    X_val = X_all[-l_val:]

    if lam_k is None:
        d_feat = int(X_all.shape[1])
        lam_k = _compute_lam_k(n_val=len(X_val), d_feat=d_feat)

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


def fit_fixed_templates_from_history(
    X: np.ndarray,
    returns_df: pd.DataFrame,
    feat_dates: pd.DatetimeIndex,
    K_template: int = 6,
    n_iter: int = HMM.n_iter,
    random_state: int = HMM.random_state,
    min_obs: int = HMM.min_obs,
    n_assets: int = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build a FIXED template pool by fitting one HMM on full pre-OOS history.

    This is the corrected template initialization strategy after diagnosing
    that the paper's spawn + EMA-update mechanism collapses templates onto
    each other over a 21-year OOS (template drift over 9,000 EMA updates
    is incompatible with stable identity).

    Strategy:
      1. Fit ONE big HMM with `K_template` components on all pre-OOS data.
      2. Take its emission parameters (means in feature space, covs in
         feature space) as the fixed templates for W2 distance.
      3. Pre-compute the forward-return mean and covariance for each
         template by grouping pre-OOS returns by HMM state assignment.
         These are the MVO inputs that will be plugged into the
         portfolio optimizer at each OOS day, weighted by template
         posterior.

    Returns
    -------
    MU_feat   : (G, d_feat) template means in 15-dim feature space —
                used for the W2 distance to HMM components in the loop.
    SIG_feat  : (G, d_feat, d_feat) template covariances in feature space.
    MU_ret    : (G, N) forward-return mean per template — used in MVO.
    SIG_ret   : (G, N, N) forward-return Ledoit-Wolf covariance per template.
    is_valid  : (G,) bool — True if the template has >= min_obs samples
                in the calibration window (and so MU_ret/SIG_ret are
                meaningful estimates). The OOS loop must skip templates
                where is_valid is False — i.e. drop their posterior mass.
    """
    if n_assets is None:
        n_assets = returns_df.shape[1]

    print(f"  Fitting fixed-template HMM with K_template={K_template} on {len(X)} pre-OOS rows...")
    hmm = fit_hmm(X, K=K_template, n_iter=n_iter, random_state=random_state)
    if hmm is None:
        raise RuntimeError("Fixed-template HMM fit failed.")

    # Hard-assign states to ALL training rows
    z_all = hmm.predict(X)
    print(f"  HMM trained. State occupancy: {np.bincount(z_all)}")

    # Feature-space emission parameters (templates for W2 distance)
    MU_feat  = hmm.means_.copy()             # (G, d_feat)
    SIG_feat = hmm.covars_.copy()            # (G, d_feat, d_feat)

    # Forward-return moments per template (for MVO inputs)
    # Align z_all with returns: z_all[i] is the regime ASSIGNED at feature row i,
    # which uses data up to feat_dates[i]. The corresponding forward return is
    # returns[feat_dates[i] + 1 trading day], i.e. the next-day return.
    feat_idx = pd.Index(feat_dates)
    ret_idx  = returns_df.index

    MU_ret  = np.zeros((K_template, n_assets))
    SIG_ret = np.zeros((K_template, n_assets, n_assets))
    is_valid = np.zeros(K_template, dtype=bool)

    # For each state, collect the next-day returns of feature rows assigned to it
    for k in range(K_template):
        mask = (z_all == k)
        if mask.sum() < min_obs:
            MU_ret[k]  = 0.0
            SIG_ret[k] = np.eye(n_assets)
            is_valid[k] = False
            continue

        # For each feature row in this state, find next-day return
        rows_in_state = np.where(mask)[0]
        next_rets = []
        for i in rows_in_state:
            # next trading day after feat_dates[i]
            try:
                pos = ret_idx.get_loc(feat_idx[i])
                if pos + 1 < len(ret_idx):
                    next_rets.append(returns_df.iloc[pos + 1].to_numpy())
            except KeyError:
                continue

        if len(next_rets) < min_obs:
            is_valid[k] = False
            MU_ret[k]  = 0.0
            SIG_ret[k] = np.eye(n_assets)
            continue

        next_rets_arr = np.asarray(next_rets)
        MU_ret[k]  = next_rets_arr.mean(axis=0)
        SIG_ret[k] = LedoitWolf().fit(next_rets_arr).covariance_
        is_valid[k] = True

    print(f"  Templates valid: {is_valid.sum()}/{K_template} (need >= {min_obs} samples each)")
    return MU_feat, SIG_feat, MU_ret, SIG_ret, is_valid



def map_components_to_fixed_templates(
    MU_K_feat: np.ndarray, SIG_K_feat: np.ndarray,
    MU_tpl_feat: np.ndarray, SIG_tpl_feat: np.ndarray,
    is_valid_tpl: np.ndarray,
) -> tuple[list[int], np.ndarray]:
    """Hard argmin assignment of each HMM component to its nearest fixed
    template, using W2 distance over FEATURE-space emission parameters
    (not forward returns — see diagnostic finding).

    Skips invalid templates (those with insufficient calibration data).

    Returns
    -------
    map_k_to_g : list of length K, integer assignments.
    D          : (K, G) full W2 distance matrix.
    """
    K = MU_K_feat.shape[0]
    G = MU_tpl_feat.shape[0]
    if not is_valid_tpl.any():
        raise RuntimeError("No valid templates available.")

    SIG_tpl_sqrt = [_sqrtm_psd(SIG_tpl_feat[g]) for g in range(G)]
    D = np.zeros((K, G))
    for k in range(K):
        for g in range(G):
            if not is_valid_tpl[g]:
                D[k, g] = np.inf
            else:
                D[k, g] = w2_gaussian(
                    MU_tpl_feat[g], SIG_tpl_feat[g],
                    MU_K_feat[k],   SIG_K_feat[k],
                    S2_sqrt=SIG_tpl_sqrt[g] if is_valid_tpl[g] else None,
                )
    map_k_to_g = [int(np.argmin(D[k])) for k in range(K)]
    return map_k_to_g, D


def aggregate_template_posteriors_hard(
    pK: np.ndarray, map_k_to_g: list[int], G: int,
) -> np.ndarray:
    """Hard argmin aggregation (paper's original): sum component
    posteriors over all components mapped to each template.

    Used with the new fixed-template framework where templates are
    genuinely distinct (because they come from a single large-K HMM fit
    on the full pre-OOS history), so hard argmin gives clean,
    interpretable regime labels.
    """
    pK = np.asarray(pK).ravel()
    pG = np.zeros(G)
    for k, g in enumerate(map_k_to_g):
        pG[g] += pK[k]
    s = pG.sum()
    if s > 1e-12:
        pG /= s
    return pG



