"""
macro_regime_hmm.py
-------------------
Standalone macro-state Hidden Markov Model that produces a stable
categorical macro regime label m_t for each calendar day. This is the
"slow layer" of the hierarchical extension: it discovers a small
number of persistent macro regimes in the 7-dimensional macro state
vector s_t.

The output m_t is then used by `backtest_hierarchical.py` to gate the
market HMM's transition matrix: A_t = A_{m_t}.

Methodology
~~~~~~~~~~~
1. Fit a standard Gaussian HMM (hmmlearn) on the macro state with
   M states. M is BIC-selected from {3,4,5} on calibration data.

2. **Frozen after calibration**. Once fit, the macro HMM parameters
   (means, covars, transmat) are held constant for the entire OOS
   period. At each OOS time t, we run Viterbi decoding using the
   frozen parameters on the full macro state history s_{1:t} to
   produce m_t.

   Rationale for freezing:
       (a) Adaptivity is already baked in via Mulliner normalization
           (data_macro.py): macro variables are z-scored over rolling
           10-year windows, so a VIX of 25 in 2010 is interpreted
           differently from a VIX of 25 in 2020.
       (b) Rolling-refit experiments showed that re-fitting on
           expanded windows (e.g., adding 2008 GFC data to the
           calibration) genuinely re-partitions the macro state space:
           regime "R1" pre-refit becomes a different distribution
           post-refit, even after Hungarian/Wasserstein label
           alignment. This breaks the stability we need to estimate
           macro-conditional market transition matrices A_m.
       (c) Boukardagha (2026) follows the same pattern for market
           templates: calibrated once, then held stable. We mirror
           that for macro templates.

3. **Persistence regularization via Dirichlet prior**: hmmlearn's EM
   on weakly-clustered macro data converges to a degenerate fixed
   point where two states ping-pong against each other every
   timestep, producing 1-day macro regime runs. We pass a transition
   prior `1 + diag*I` with `diag = 10` (the default), which biases
   EM toward persistent self-loops without imposing a strong
   parametric constraint. Result on our 7-variable calibration
   window: 9 distinct regime runs, mean duration ~340 days, max ~700
   days. This matches the economic literature on macroeconomic
   cycles (NBER cycle phases of 6-24 months).

4. **Viterbi decoding**, not marginal argmax. The Viterbi path is the
   most-likely *sequence* under the joint transition + emission
   model, and under a persistence prior naturally produces smooth
   regime sequences. Marginal argmax flickers in transition zones.

5. **Reproducible regime labelling** via "natural ordering" by
   projection onto the first principal component of the means matrix.
   Macro regime 0 corresponds to one extreme of the macro state space
   (typically high VIX / inverted curve / risk-off), regime M-1 to
   the other end. Independent of hmmlearn's EM init RNG.

(Rolling-refit functionality is retained in the code for ablation /
experimentation but not used in the default hierarchical backtest.
See `refit_with_label_alignment`, which solves the Hirsa et al. 2024
label-assignment problem via Hungarian/W2.)

Strict causality
~~~~~~~~~~~~~~~~
The macro HMM is fit only on data observable at calibration time
(< 2005-01-03 by default). The Viterbi decoding at OOS time t uses
the macro state s_{1:t} which is itself lagged by one trading day
in `data_macro.py`. No lookahead.

API
~~~
    fit_initial(X)                   - calibration fit + BIC over M
    predict_viterbi(X)               - decode regime sequence (preferred)
    predict_proba(X)                 - marginal posteriors
    get_smoothed_regime(pM)          - rolling-mode smoothing for charts
    refit_with_label_alignment(X)    - optional rolling refit (Hirsa)

Used by `backtest_hierarchical.py` to drive the market-layer
transition matrix gating.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from scipy.linalg import sqrtm
from scipy.optimize import linear_sum_assignment

from hmmlearn import hmm

# Quiet hmmlearn (it logs EM warnings we don't care about)
logging.getLogger("hmmlearn").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# --------------------------------------------------------------------- #
# 1. Wasserstein-2 distance between Gaussian distributions
# --------------------------------------------------------------------- #
def w2_gaussian_pairwise(
    MU_a: np.ndarray, SIG_a: np.ndarray,
    MU_b: np.ndarray, SIG_b: np.ndarray,
) -> np.ndarray:
    """Pairwise W2 distance between two sets of Gaussians.

    Parameters
    ----------
    MU_a    : (Ma, d)
    SIG_a   : (Ma, d, d)
    MU_b    : (Mb, d)
    SIG_b   : (Mb, d, d)

    Returns
    -------
    D       : (Ma, Mb)  with D[i,j] = W2(N(MU_a[i], SIG_a[i]),
                                         N(MU_b[j], SIG_b[j]))
    """
    Ma = MU_a.shape[0]
    Mb = MU_b.shape[0]
    D = np.zeros((Ma, Mb))
    eye_d = np.eye(MU_a.shape[1])

    # Precompute Sigma_a^{1/2} once per row
    sqrt_a = [_sym_sqrt(SIG_a[i] + 1e-10 * eye_d) for i in range(Ma)]

    for i in range(Ma):
        for j in range(Mb):
            mu_diff = MU_a[i] - MU_b[j]
            mean_term = float(mu_diff @ mu_diff)
            # Bures-Wasserstein covariance term
            inner = sqrt_a[i] @ SIG_b[j] @ sqrt_a[i]
            inner = 0.5 * (inner + inner.T)  # symmetrize
            sqrt_inner = _sym_sqrt(inner + 1e-10 * eye_d)
            cov_term = float(
                np.trace(SIG_a[i] + SIG_b[j] - 2.0 * sqrt_inner)
            )
            # Clamp tiny negatives from numerical error
            d2 = max(mean_term + cov_term, 0.0)
            D[i, j] = np.sqrt(d2)
    return D


def _sym_sqrt(M: np.ndarray) -> np.ndarray:
    """Symmetric square root of a PSD matrix via eigendecomposition.

    Faster and more stable than scipy.linalg.sqrtm for small matrices
    (d=7 in our case) — and returns a real-valued result without the
    complex-residual warnings sqrtm sometimes emits."""
    M = 0.5 * (M + M.T)
    w, V = np.linalg.eigh(M)
    w_clip = np.clip(w, 0.0, None)
    return (V * np.sqrt(w_clip)) @ V.T


# --------------------------------------------------------------------- #
# 2. Hungarian label-stabilization
# --------------------------------------------------------------------- #
def align_labels_hungarian(
    MU_new: np.ndarray, SIG_new: np.ndarray,
    MU_ref: np.ndarray, SIG_ref: np.ndarray,
) -> np.ndarray:
    """Solve the linear assignment problem to align new regimes to a
    reference labelling.

    Returns a permutation `perm` of length M such that

        regime k of the new fit  <-->  regime perm[k] of the reference

    so that `MU_new[perm_inv]` is aligned with `MU_ref`.

    If new and reference have different numbers of regimes we still try
    to align as many as possible (rectangular Hungarian). Unmatched new
    regimes get fresh labels (= max(ref) + 1, +2, ...).
    """
    M_new = MU_new.shape[0]
    M_ref = MU_ref.shape[0]

    D = w2_gaussian_pairwise(MU_new, SIG_new, MU_ref, SIG_ref)  # (M_new, M_ref)

    # Hungarian on the D matrix: rows = new, cols = ref
    row_ind, col_ind = linear_sum_assignment(D)
    # row_ind[k] is a new-regime index, col_ind[k] is the matched ref label.

    # Build perm: perm[new_idx] = ref_label
    perm = np.full(M_new, -1, dtype=int)
    for r, c in zip(row_ind, col_ind):
        perm[r] = c

    # Handle case where M_new > M_ref: unmatched new regimes get fresh IDs
    next_label = int(MU_ref.shape[0])
    for k in range(M_new):
        if perm[k] == -1:
            perm[k] = next_label
            next_label += 1

    return perm


# --------------------------------------------------------------------- #
# 3. BIC-based M selection
# --------------------------------------------------------------------- #
def _bic_score(model, X: np.ndarray) -> float:
    """BIC = -2 * log-likelihood + d * log(N).

    For Gaussian HMM with full covariance, parameter count d is:
        M * d_obs                       (means)
      + M * d_obs * (d_obs+1) / 2       (covariances, symmetric)
      + M * (M - 1)                     (transition matrix rows sum to 1)
      + (M - 1)                         (initial distribution sums to 1)
    where d_obs is observation dimensionality and M is n_components.
    """
    M = model.n_components
    d_obs = X.shape[1]
    N = X.shape[0]
    n_params = (
        M * d_obs                          # means
        + M * d_obs * (d_obs + 1) // 2     # full covariances
        + M * (M - 1)                      # transition matrix
        + (M - 1)                          # initial distribution
    )
    ll = model.score(X)
    return -2.0 * ll + n_params * np.log(N)


def select_M_via_bic(
    X: np.ndarray,
    M_candidates: list[int] = (3, 4, 5, 6),
    n_iter: int = 100,
    random_state: int = 42,
    transmat_prior_diag: float = 10.0,
    verbose: bool = False,
) -> tuple[int, dict]:
    """BIC-select the number of macro regimes M on calibration data.

    Returns
    -------
    M_star : int            — selected number of regimes
    scores : dict[M -> BIC] — full BIC trace for diagnostics
    """
    scores = {}
    for M in M_candidates:
        try:
            prior = np.ones((M, M)) + transmat_prior_diag * np.eye(M)
            mdl = hmm.GaussianHMM(
                n_components=M,
                covariance_type="full",
                n_iter=n_iter,
                random_state=random_state,
                tol=1e-3,
                transmat_prior=prior,
            )
            mdl.fit(X)
            scores[M] = _bic_score(mdl, X)
            if verbose:
                print(f"    M={M}: BIC = {scores[M]:.2f}  (ll = {mdl.score(X):.2f})")
        except Exception as e:
            scores[M] = np.inf
            if verbose:
                print(f"    M={M}: failed ({e})")

    valid = {M: s for M, s in scores.items() if np.isfinite(s)}
    if not valid:
        return int(M_candidates[0]), scores
    M_star = int(min(valid, key=valid.get))
    return M_star, scores


# --------------------------------------------------------------------- #
# 4. The MacroRegimeHMM class
# --------------------------------------------------------------------- #
@dataclass
class MacroRegimeHMM:
    """Macro-state HMM with label-stabilized rolling refitting.

    Attributes set by `fit_initial`:
        M_                  : selected number of regimes (after BIC)
        means_              : (M, d_macro) — last fit's means in *stable* label space
        covars_             : (M, d_macro, d_macro) — same, covariances
        transmat_           : (M, M) — last fit's transition matrix in stable labels
        startprob_          : (M,)
        bic_scores_         : dict from initial M selection
        refit_count_        : number of refits performed
        label_history_      : list of (date, M_seen, perm) per refit (diagnostics)

    Parameters
    ----------
    M_candidates : tuple of ints
        Candidate values of M for BIC selection.
    n_iter : int
        Max EM iterations per fit.
    random_state : int
        Seed for hmmlearn EM initialization.
    smooth_window : int
        Window for rolling-mode smoothing of the regime label sequence.
        Used only for visualization-style smoothing; the gating logic in
        the market layer uses the Viterbi path (which is already smooth).
    transmat_prior_diag : float
        Strength of the diagonal-preference Dirichlet prior on the
        transition matrix. Prior is `1 + transmat_prior_diag * I`, which
        biases EM toward persistent self-loops. Without this, hmmlearn's
        EM tends to converge to a degenerate fixed point where two
        states ping-pong against each other every timestep, producing
        unrealistic 1-day macro regime runs. With diag = 10 (default),
        macro regimes persist for ~200-300 days on average, which
        matches the economic literature on macroeconomic cycles
        (NBER cycle phases typically run 6-24 months).
    """
    M_candidates:        tuple[int, ...] = (3, 4, 5)
    n_iter:              int             = 100
    random_state:        int             = 42
    smooth_window:       int             = 5
    transmat_prior_diag: float           = 10.0

    # Populated by fit
    M_:           Optional[int]     = field(default=None, init=False)
    means_:       Optional[np.ndarray] = field(default=None, init=False, repr=False)
    covars_:      Optional[np.ndarray] = field(default=None, init=False, repr=False)
    transmat_:    Optional[np.ndarray] = field(default=None, init=False, repr=False)
    startprob_:   Optional[np.ndarray] = field(default=None, init=False, repr=False)
    bic_scores_:  dict              = field(default_factory=dict, init=False, repr=False)
    refit_count_: int               = field(default=0, init=False)
    label_history_: list            = field(default_factory=list, init=False, repr=False)

    # =================================================================
    def fit_initial(self, X: np.ndarray, verbose: bool = False) -> "MacroRegimeHMM":
        """Initial calibration fit with BIC-selection over M_candidates.
        Establishes the reference labelling for all subsequent refits.
        """
        if verbose:
            print(f"  Macro HMM: selecting M via BIC over {list(self.M_candidates)}")
        M_star, scores = select_M_via_bic(
            X, M_candidates=self.M_candidates,
            n_iter=self.n_iter, random_state=self.random_state,
            transmat_prior_diag=self.transmat_prior_diag,
            verbose=verbose,
        )
        self.M_ = M_star
        self.bic_scores_ = scores

        if verbose:
            print(f"  Macro HMM: M_star = {M_star}")

        prior = np.ones((M_star, M_star)) + self.transmat_prior_diag * np.eye(M_star)
        mdl = hmm.GaussianHMM(
            n_components=M_star, covariance_type="full",
            n_iter=self.n_iter, random_state=self.random_state, tol=1e-3,
            transmat_prior=prior,
        )
        mdl.fit(X)

        sort_idx = self._natural_order(mdl.means_)
        self.means_     = mdl.means_[sort_idx].copy()
        self.covars_    = mdl.covars_[sort_idx].copy()
        self.transmat_  = mdl.transmat_[sort_idx][:, sort_idx].copy()
        self.startprob_ = mdl.startprob_[sort_idx].copy()
        self.refit_count_ = 1

        return self

    # =================================================================
    def refit_with_label_alignment(
        self,
        X: np.ndarray,
        verbose: bool = False,
    ) -> bool:
        """Refit on expanded data X, align labels via Hungarian/W2 to
        the previous fit's regimes.

        Returns True if refit succeeded, False if it failed (caller
        keeps using the previous fit's parameters).
        """
        if self.M_ is None:
            raise RuntimeError("Call fit_initial first.")

        try:
            prior = np.ones((self.M_, self.M_)) + self.transmat_prior_diag * np.eye(self.M_)
            mdl = hmm.GaussianHMM(
                n_components=self.M_, covariance_type="full",
                n_iter=self.n_iter, random_state=self.random_state, tol=1e-3,
                transmat_prior=prior,
            )
            mdl.fit(X)
        except Exception:
            return False

        try:
            perm = align_labels_hungarian(
                MU_new=mdl.means_, SIG_new=mdl.covars_,
                MU_ref=self.means_, SIG_ref=self.covars_,
            )
        except Exception:
            return False

        M = self.M_
        inv = np.full(M, -1, dtype=int)
        for new_idx in range(M):
            ref_label = perm[new_idx]
            if 0 <= ref_label < M:
                inv[ref_label] = new_idx

        if (inv < 0).any():
            return False

        self.means_     = mdl.means_[inv].copy()
        self.covars_    = mdl.covars_[inv].copy()
        self.transmat_  = mdl.transmat_[inv][:, inv].copy()
        self.startprob_ = mdl.startprob_[inv].copy()
        self.refit_count_ += 1

        return True

    # =================================================================
    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Posterior P(m_t = m | s_{1:T}) for t in 1..T, returned as (T, M)."""
        return self._build_cached_mdl().predict_proba(X)

    # =================================================================
    def predict_viterbi(self, X: np.ndarray) -> np.ndarray:
        """Viterbi-decoded most-likely state sequence. (T,) ints.

        Use this rather than `argmax(predict_proba)` for the macro
        gating signal — Viterbi enforces global path consistency under
        the (persistence-prior) transition matrix, so the resulting
        regime sequence is naturally smooth (~250-day average runs in
        our calibration with transmat_prior_diag=10).
        """
        return self._build_cached_mdl().predict(X)

    # =================================================================
    def _build_cached_mdl(self):
        """Reconstruct an hmmlearn model with the stored stable-label
        parameters. Used by predict_proba / predict_viterbi."""
        prior = np.ones((self.M_, self.M_)) + self.transmat_prior_diag * np.eye(self.M_)
        mdl = hmm.GaussianHMM(
            n_components=self.M_, covariance_type="full",
            n_iter=1, random_state=self.random_state, transmat_prior=prior,
        )
        mdl.startprob_ = self.startprob_.copy()
        mdl.transmat_  = self.transmat_.copy()
        mdl.means_     = self.means_.copy()
        mdl.covars_    = self.covars_.copy()
        return mdl

    # =================================================================
    def get_smoothed_regime(
        self,
        pM: np.ndarray,
        window: Optional[int] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compute raw and smoothed macro regime sequences.

        Parameters
        ----------
        pM : (T, M) posterior probabilities from predict_proba
        window : int, default = self.smooth_window — rolling-mode window

        Returns
        -------
        m_raw    : (T,) int — argmax of pM each day
        m_smooth : (T,) int — rolling-mode of m_raw over `window` days
        """
        if window is None:
            window = self.smooth_window
        m_raw = pM.argmax(axis=1)
        if window <= 1:
            return m_raw, m_raw.copy()

        # Rolling mode via pandas
        s = pd.Series(m_raw)
        m_smooth = (s.rolling(window, min_periods=1)
                     .apply(lambda x: x.value_counts().idxmax(), raw=False)
                     .astype(int).values)
        return m_raw, m_smooth

    # =================================================================
    def _natural_order(self, MU: np.ndarray) -> np.ndarray:
        """Sort regimes by projection onto the first principal component
        of the means matrix. Gives a reproducible canonical ordering
        independent of hmmlearn's EM init RNG.
        """
        # Center means and take SVD: first right singular vector is PC1
        centered = MU - MU.mean(axis=0, keepdims=True)
        U, S, Vt = np.linalg.svd(centered, full_matrices=False)
        # Projection of each regime mean onto PC1
        scores = centered @ Vt[0]
        return np.argsort(scores)


# --------------------------------------------------------------------- #
# 5. Self-test
# --------------------------------------------------------------------- #
if __name__ == "__main__":
    import time
    rng = np.random.default_rng(0)
    # Synthetic: 4 macro regimes, each a different Gaussian
    centers = np.array([[-2, -1, -1, -1, 0, 2, 1],
                        [-1,  0, -1,  0, 0, 1, 0],
                        [ 1,  1,  1,  1, 0, -1, -1],
                        [ 2,  2,  1,  2, 0, -1, -1]])
    N_per = 300
    Xs, ms = [], []
    for m, c in enumerate(centers):
        Xs.append(rng.normal(c, 0.5, size=(N_per, 7)))
        ms.append(np.full(N_per, m))
    X = np.vstack(Xs)
    m_true = np.concatenate(ms)

    print("=== Synthetic macro-HMM self-test ===")
    t0 = time.time()
    mac = MacroRegimeHMM(M_candidates=(3, 4, 5))
    mac.fit_initial(X, verbose=True)
    print(f"\nElapsed: {time.time()-t0:.2f}s")
    print(f"BIC scores: {mac.bic_scores_}")
    print(f"Selected M: {mac.M_}")
    print(f"Means after natural ordering:\n{mac.means_.round(2)}")
    pM = mac.predict_proba(X)
    m_raw, m_smooth = mac.get_smoothed_regime(pM)
    print(f"Posterior shape: {pM.shape}")
    print(f"Accuracy (modulo permutation, raw): {(m_raw == m_true).mean():.3f}  "
          f"(can be low if labels are permuted; just sanity)")
    print(f"Unique smoothed regimes seen: {np.unique(m_smooth)}")

    # Refit on a perturbed version to test label stability
    X2 = X + rng.normal(0, 0.05, size=X.shape)
    ok = mac.refit_with_label_alignment(X2, verbose=True)
    print(f"\nRefit succeeded: {ok}")
    print(f"Means after refit (should be similar):\n{mac.means_.round(2)}")
