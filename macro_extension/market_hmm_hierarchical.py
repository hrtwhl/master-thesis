"""
market_hmm_hierarchical.py
--------------------------
Market-layer Hidden Markov Model with macro-regime-conditional
transition matrices. This is the "fast layer" of the hierarchical
extension.

Model
~~~~~
Standard Gaussian HMM on the 15-dim market feature space X, except
the transition matrix at time t depends on the categorical macro
regime m_t (computed beforehand by the macro layer in
`macro_regime_hmm.py`):

    P(Z_t = j | Z_{t-1} = i, m_t)  =  A_{m_t}[i, j]

We learn M separate K-by-K transition matrices, one per macro regime,
all jointly via EM:

    E-step
        Forward-backward with time-varying A_t = A_{m_t}. We reuse
        the Numba-JIT'd forward-backward kernel from
        wasserstein_hmm_macro.py (v9), since that's already a
        time-varying-A solver. Only the *log_A tensor* changes here:
        instead of softmax(beta . s_t), we build it by indexing.

    M-step
        Emissions:        weighted MLE (closed-form Gaussian).
        Transition A_m:   Dirichlet-smoothed counting,
            A_m[i, j] = (sum_{t: m_t=m} xi_t[i, j] + alpha)
                        / (sum_{t: m_t=m} sum_k xi_t[i, k] + K*alpha).
                This is the standard MAP estimator for a multinomial
                with Dirichlet(alpha+1) prior. With alpha=1 (default),
                a macro regime with no transition counts at all gets
                a uniform 1/K transition row (sensible default).
        Initial dist:     gamma_1 normalized.

Why this is dramatically simpler than v9
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
v9 modeled A_t = softmax(beta_0 + beta . s_t), which required
gradient-based optimization (L-BFGS), L2 regularization,
warm-starts, and convergence drama. Over 1,096 refits these
ingredients compounded into a mode-collapse pathology where one
state absorbed everything.

v10 replaces this with categorical gating: M discrete macro regimes,
each with its own K-by-K transition matrix, estimated by closed-form
counting with a Dirichlet smoother. Closed-form M-step means no
gradient optimizer, no regularization hyperparameter for transitions,
no warm-start chains, no mode collapse. Just count-and-normalize.

K (the number of market regimes) is still variable: BIC-selected from
{3..6} like the baseline. K varies independently of M (the number of
macro regimes).

Strict causality
~~~~~~~~~~~~~~~~
At OOS time t, the macro layer has produced m_{1:t} using only
information up to t. The market HMM fits / predicts using X_{1:t}
and m_{1:t} only. No lookahead.
"""

from __future__ import annotations

import _paths  # noqa: F401  (injects v9 + baseline into sys.path)

import logging
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from sklearn.cluster import KMeans

# Reuse the JIT forward-backward + log-Gaussian-emit from v9
from wasserstein_hmm_macro import (
    _forward_backward,
    _log_gaussian_emit,
    _gaussian_M_step,
    _ffill_2d,
)

logging.getLogger("hmmlearn").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# --------------------------------------------------------------------- #
# 1. Build log-A tensor from a categorical m_t and a set of A_m matrices
# --------------------------------------------------------------------- #
def build_log_A_tensor_categorical(
    A_per_macro: np.ndarray,   # (M, K, K)
    m_seq: np.ndarray,         # (T,) int  — macro regime label per market day
) -> np.ndarray:
    """Construct the (T-1, K, K) log transition tensor by indexing the
    macro regime at each step.

        log_A[t, i, j] = log A_{m_{t+1}}[i, j]

    (transition from t to t+1 is conditioned on m_{t+1}, matching the
    convention used in v9: A_t uses s_{t+1}.)
    """
    T = m_seq.shape[0]
    # Use m_{t+1} to gate the transition t -> t+1
    m_use = m_seq[1:]                                # (T-1,)
    # Clip to valid range (defensive)
    M = A_per_macro.shape[0]
    m_use = np.clip(m_use, 0, M - 1)

    # Compute log A_per_macro once
    A_safe = np.clip(A_per_macro, 1e-12, None)
    log_A_per_macro = np.log(A_safe)                 # (M, K, K)

    # Index into log_A_per_macro by m_use
    log_A_tensor = log_A_per_macro[m_use]            # (T-1, K, K)
    return log_A_tensor


# --------------------------------------------------------------------- #
# 2. M-step: estimate A_m by Dirichlet-smoothed counting
# --------------------------------------------------------------------- #
def estimate_A_per_macro(
    xi:        np.ndarray,       # (T-1, K, K) soft transition counts
    m_seq:     np.ndarray,       # (T,) macro regime labels
    M:         int,              # number of macro regimes
    alpha:     float = 1.0,      # Dirichlet pseudo-count
) -> np.ndarray:
    """For each macro regime m, sum xi over days where m_{t+1} = m,
    then normalize with Dirichlet smoothing.

    Returns A_per_macro of shape (M, K, K).

    The smoothing parameter `alpha` is a per-cell pseudo-count. With
    K cells per row, the effective prior strength per row is K*alpha.
    For our K~5, alpha=1.0 means a row with no observations defaults
    to uniform 1/K (which is what we want for unseen origin states).
    """
    T_minus_1, K, _ = xi.shape
    m_use = m_seq[1:]                                # (T-1,)
    A_per_macro = np.zeros((M, K, K))

    for m in range(M):
        mask = (m_use == m)
        if mask.any():
            N_m = xi[mask].sum(axis=0)                # (K, K)
        else:
            N_m = np.zeros((K, K))
        # Dirichlet-smoothed posterior
        N_m_smoothed = N_m + alpha
        row_sums = N_m_smoothed.sum(axis=1, keepdims=True)
        A_per_macro[m] = N_m_smoothed / row_sums

    return A_per_macro


# --------------------------------------------------------------------- #
# 3. The HierarchicalMarketHMM class
# --------------------------------------------------------------------- #
@dataclass
class HierarchicalMarketHMM:
    """Market HMM with macro-regime-gated transitions, fitted via custom EM.

    Same API surface as v9's MacroHMM and hmmlearn's GaussianHMM so it
    slots into the backtest with minimal change.

    Parameters
    ----------
    n_components : int
        Number of market (latent) states K.
    n_macro_regimes : int
        Number of macro regimes M (from the macro layer).
    n_iter : int
        Max EM iterations.
    tol : float
        Convergence tolerance on log-likelihood improvement.
    dirichlet_alpha : float
        Pseudo-count for transition smoothing. Larger -> rows pulled
        toward uniform 1/K. Default 1.0.
    random_state : int
        Seed for the k-means initialization.
    frozen_emissions : bool
        If True, emissions are held constant at fixed_emission_means /
        _covs and only A_m matrices are learned. Forces K = G when
        used with the hierarchical backtest. Diagnosed in v10.0 as
        producing pathologically one-hot pG and constant MVO weights
        — not recommended for production runs, kept for ablation.
    fixed_emission_means : (K, d) array or None
        Required when frozen_emissions=True. Held constant during EM.
    fixed_emission_covs : (K, d, d) array or None
        Required when frozen_emissions=True. Held constant during EM.
    template_init_means : (K, d) array or None  ***NEW in v10.1***
        Optional initialization for emission means when
        frozen_emissions=False. If provided, EM starts from these
        values instead of k-means. Used by the hierarchical backtest
        to anchor the unfrozen HMM near the baseline template
        structure (which prevents the components from collapsing
        toward the global mean under macro-gated EM, which we
        diagnosed in v10.0 unfrozen test runs).
    template_init_covs : (K, d, d) array or None
        Companion to template_init_means.
    """
    n_components:        int            = 5
    n_macro_regimes:     int            = 5
    n_iter:              int            = 30
    tol:                 float          = 1e-3
    dirichlet_alpha:     float          = 1.0
    random_state:        int            = 42
    frozen_emissions:    bool           = False
    fixed_emission_means: Optional[np.ndarray] = None
    fixed_emission_covs:  Optional[np.ndarray] = None
    template_init_means:  Optional[np.ndarray] = None
    template_init_covs:   Optional[np.ndarray] = None
    mean_anchor_strength: float         = 0.0   # effective # of "anchor" observations
    freeze_covars:        bool          = False  # if True, covars stay at init values

    # Fitted attributes
    means_:              Optional[np.ndarray] = field(default=None, init=False, repr=False)
    covars_:             Optional[np.ndarray] = field(default=None, init=False, repr=False)
    A_per_macro_:        Optional[np.ndarray] = field(default=None, init=False, repr=False)
    startprob_:          Optional[np.ndarray] = field(default=None, init=False, repr=False)
    log_lik_history_:    list = field(default_factory=list, init=False, repr=False)
    converged_:          bool = field(default=False, init=False, repr=False)

    # =================================================================
    def _init_emissions(self, X: np.ndarray):
        """Initialize emissions in priority order:
            1. If frozen_emissions: use fixed_emission_means/_covs (held constant)
            2. If template_init_means provided: use as init but allow EM updates
            3. Otherwise: k-means on X (legacy behavior)
        """
        K = self.n_components

        if self.frozen_emissions:
            if self.fixed_emission_means is None or self.fixed_emission_covs is None:
                raise ValueError("frozen_emissions=True requires "
                                 "fixed_emission_means and fixed_emission_covs.")
            means = self.fixed_emission_means.copy()
            covs  = self.fixed_emission_covs.copy()
            startprob = np.full(K, 1.0 / K)
            return means, covs, startprob

        if self.template_init_means is not None and self.template_init_covs is not None:
            # Template-anchored init: start EM from these means/covs,
            # but allow them to update. The hierarchical backtest passes
            # baseline template emissions here (with K-template clustering
            # if K != G) so the K-components stay near interpretable
            # regime locations rather than collapsing to the global mean.
            if self.template_init_means.shape[0] != K:
                raise ValueError(
                    f"template_init_means has K={self.template_init_means.shape[0]}, "
                    f"but n_components={K}. Caller should pre-cluster templates to K."
                )
            means = self.template_init_means.copy()
            covs  = self.template_init_covs.copy()
            startprob = np.full(K, 1.0 / K)
            return means, covs, startprob

        # Legacy k-means init (used only when no template init provided)
        km = KMeans(n_clusters=K, n_init=5, random_state=self.random_state)
        labels = km.fit_predict(X)
        means = km.cluster_centers_
        covs = np.zeros((K, X.shape[1], X.shape[1]))
        eye_d = np.eye(X.shape[1])
        for k in range(K):
            mask = labels == k
            if mask.sum() > 10:
                covs[k] = np.cov(X[mask].T) + 1e-3 * eye_d
            else:
                covs[k] = np.cov(X.T) + 1e-3 * eye_d
        startprob = np.bincount(labels, minlength=K).astype(float)
        startprob = (startprob + 1.0) / (startprob.sum() + K)
        return means, covs, startprob

    # =================================================================
    def fit(self, X: np.ndarray, m_seq: np.ndarray, verbose: bool = False) -> "HierarchicalMarketHMM":
        """Run EM on (X, m_seq).

        Parameters
        ----------
        X     : (T, d_feat) market feature matrix.
        m_seq : (T,) int — macro regime label per day, from the macro layer.
        """
        T, d = X.shape
        if m_seq.shape[0] != T:
            raise ValueError(f"X and m_seq must have same length; got {T} vs {m_seq.shape[0]}")
        m_seq = np.asarray(m_seq, dtype=np.int64)

        K = self.n_components
        M = self.n_macro_regimes

        means, covs, startprob = self._init_emissions(X)
        log_pi = np.log(startprob + 1e-12)

        # Initialize A_per_macro to uniform 1/K (Dirichlet prior alone)
        A_per_macro = np.full((M, K, K), 1.0 / K)

        prev_log_lik = -np.inf
        for it in range(self.n_iter):
            # E-step
            log_emit = _log_gaussian_emit(X, means, covs)             # (T, K)
            log_A = build_log_A_tensor_categorical(A_per_macro, m_seq) # (T-1, K, K)
            log_alpha, log_beta, log_gamma, log_lik = _forward_backward(
                log_pi, log_A, log_emit
            )
            gamma = np.exp(log_gamma)
            log_xi = (log_alpha[:-1, :, None]
                      + log_A
                      + log_emit[1:, None, :]
                      + log_beta[1:, None, :]
                      - log_lik)
            xi = np.exp(log_xi)                                       # (T-1, K, K)

            self.log_lik_history_.append(log_lik)
            if verbose:
                print(f"  EM iter {it+1:2d}: log_lik = {log_lik:.2f}", flush=True)

            if it > 0 and log_lik - prev_log_lik < self.tol:
                self.converged_ = True
                break
            prev_log_lik = log_lik

            # M-step: emissions (skipped if frozen)
            if not self.frozen_emissions:
                means_mle, covs_mle = _gaussian_M_step(gamma, X)
                # Anchor means toward template_init if requested.
                # μ_k_new = (Σ γ x + λ * μ_anchor) / (Σ γ + λ)
                # = (N_k * μ_mle + λ * μ_anchor) / (N_k + λ)
                # where N_k = Σ_t γ_t,k is effective sample size.
                if (self.mean_anchor_strength > 0
                        and self.template_init_means is not None):
                    N_k = gamma.sum(axis=0)            # (K,)
                    lam = self.mean_anchor_strength
                    means_anchored = np.zeros_like(means_mle)
                    for k in range(K):
                        nk = N_k[k]
                        # If a component has zero mass, snap to anchor
                        if nk < 1e-6:
                            means_anchored[k] = self.template_init_means[k]
                        else:
                            means_anchored[k] = (
                                (nk * means_mle[k] + lam * self.template_init_means[k])
                                / (nk + lam)
                            )
                    means = means_anchored
                else:
                    means = means_mle
                # Covariances: either update via MLE or freeze at init
                if self.freeze_covars and self.template_init_covs is not None:
                    pass  # keep current covs (which are init since not updated above)
                else:
                    covs = covs_mle

            # M-step: initial distribution
            startprob = gamma[0] / gamma[0].sum()
            log_pi = np.log(startprob + 1e-12)

            # M-step: A_per_macro by Dirichlet-smoothed counting
            A_per_macro = estimate_A_per_macro(
                xi=xi, m_seq=m_seq, M=M, alpha=self.dirichlet_alpha,
            )

        self.means_       = means
        self.covars_      = covs
        self.A_per_macro_ = A_per_macro
        self.startprob_   = startprob
        return self

    # =================================================================
    def predict_proba(self, X: np.ndarray, m_seq: np.ndarray) -> np.ndarray:
        """Posterior P(z_t | x_{1:T}, m_{1:T}). Returns (T, K)."""
        if self.means_ is None:
            raise RuntimeError("Call .fit before .predict_proba")
        m_seq = np.asarray(m_seq, dtype=np.int64)
        log_emit = _log_gaussian_emit(X, self.means_, self.covars_)
        log_A    = build_log_A_tensor_categorical(self.A_per_macro_, m_seq)
        log_pi   = np.log(self.startprob_ + 1e-12)
        _, _, log_gamma, _ = _forward_backward(log_pi, log_A, log_emit)
        return np.exp(log_gamma)

    @property
    def n_components_(self):
        return self.n_components


# --------------------------------------------------------------------- #
# 4. Convenience wrapper with try/except fallback
# --------------------------------------------------------------------- #
def fit_hierarchical_market_hmm(
    X: np.ndarray,
    m_seq: np.ndarray,
    K: int,
    M: int,
    n_iter: int = 20,
    dirichlet_alpha: float = 1.0,
    random_state: int = 42,
    frozen_emissions: bool = False,
    fixed_emission_means: Optional[np.ndarray] = None,
    fixed_emission_covs:  Optional[np.ndarray] = None,
    template_init_means:  Optional[np.ndarray] = None,
    template_init_covs:   Optional[np.ndarray] = None,
    mean_anchor_strength: float = 0.0,
    freeze_covars:        bool  = False,
) -> Optional[HierarchicalMarketHMM]:
    """Fit a HierarchicalMarketHMM with fallback to None on irrecoverable failure."""
    try:
        model = HierarchicalMarketHMM(
            n_components=K, n_macro_regimes=M,
            n_iter=n_iter, dirichlet_alpha=dirichlet_alpha,
            random_state=random_state,
            frozen_emissions=frozen_emissions,
            fixed_emission_means=fixed_emission_means,
            fixed_emission_covs=fixed_emission_covs,
            template_init_means=template_init_means,
            template_init_covs=template_init_covs,
            mean_anchor_strength=mean_anchor_strength,
            freeze_covars=freeze_covars,
        )
        model.fit(X, m_seq, verbose=False)
        return model
    except Exception:
        return None


# --------------------------------------------------------------------- #
# 5. Self-test
# --------------------------------------------------------------------- #
if __name__ == "__main__":
    import time
    rng = np.random.default_rng(0)
    T, d, K, M = 600, 5, 3, 2
    X = rng.normal(size=(T, d))
    # Two macro regimes alternating in 200-day blocks
    m_seq = np.concatenate([
        np.zeros(200, dtype=int),
        np.ones(200, dtype=int),
        np.zeros(200, dtype=int),
    ])

    print("=== HierarchicalMarketHMM self-test ===")
    t0 = time.time()
    mdl = HierarchicalMarketHMM(n_components=K, n_macro_regimes=M, n_iter=15)
    mdl.fit(X, m_seq, verbose=True)
    print(f"\nFit took {time.time() - t0:.2f}s, converged={mdl.converged_}")
    print(f"A_per_macro shape: {mdl.A_per_macro_.shape}  (expected ({M}, {K}, {K}))")
    print(f"\nA_0 (macro regime 0):\n{mdl.A_per_macro_[0].round(3)}")
    print(f"\nA_1 (macro regime 1):\n{mdl.A_per_macro_[1].round(3)}")
    print(f"Row sums check (should be 1):")
    print(f"  A_0: {mdl.A_per_macro_[0].sum(axis=1).round(3)}")
    print(f"  A_1: {mdl.A_per_macro_[1].sum(axis=1).round(3)}")

    pp = mdl.predict_proba(X, m_seq)
    print(f"\npredict_proba shape: {pp.shape}, row sums in [{pp.sum(1).min():.3f}, {pp.sum(1).max():.3f}]")
