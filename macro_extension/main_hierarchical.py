"""
main_hierarchical.py
--------------------
Entry point for the HIERARCHICAL (two-layer HMM) macro extension.

Usage
~~~~~
    python main_hierarchical.py

What it does
~~~~~~~~~~~~
1. Loads daily asset prices + the 7-variable macro state vector.
2. Loads the baseline backtest **from disk** if it has been run already
   (cached as `Boukardagha/output/tables/daily_backtest_output.csv`).
   If the cache is missing, runs the baseline and writes it to disk.
3. Runs the HIERARCHICAL backtest:
       * Macro layer: Gaussian HMM with M BIC-selected from {3,4,5},
         calibrated once on macro state up to OOS start, then frozen.
       * Market layer: HMM with macro-regime-gated transitions
         A_t = A_{m_t}, with K = G (frozen emissions to baseline templates).
4. Runs the side-by-side comparison (charts + tables) and writes
   everything to `macro_extension_v10/output/compare/`.

To force a baseline rerun:
    main(force_baseline_rerun=True)

Wall time expectations
~~~~~~~~~~~~~~~~~~~~~~
With Numba installed and HMM.f_hmm = 5 (default):
    Baseline backtest        ~10-20  minutes (cached after first run)
    Hierarchical backtest    ~15-30 minutes
    Comparison + charts      <30 sec

v10's market layer is faster than v9 because:
  - Closed-form Dirichlet-counting M-step (no L-BFGS for transitions)
  - Frozen emissions skip the Gaussian M-step
The macro layer adds ~1 min for the one-time calibration fit.
"""

from __future__ import annotations

import _paths  # noqa: F401  (injects baseline + v9 into sys.path)

import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from config import RUN, OUTPUT_DIR as BASELINE_OUTPUT_DIR, ASSET_NAMES, ensure_dirs
from data import load_all
from data_macro import load_macro_state
from backtest import run_backtest, BacktestResult
from backtest_hierarchical import run_backtest_hierarchical
from compare_baseline_vs_hierarchical import run_full_comparison
from reporting import export_all


# --------------------------------------------------------------------- #
# Baseline cache loader (identical pattern to v9's main_macro.py)
# --------------------------------------------------------------------- #
def load_baseline_from_disk(
    tables_dir: Path | None = None,
) -> Optional[BacktestResult]:
    """Reconstruct a BacktestResult from the cached baseline daily CSV.

    Looks for `daily_backtest_output.csv` in `tables_dir` (default
    `Boukardagha/output/tables/`). Returns None if missing.
    """
    if tables_dir is None:
        tables_dir = BASELINE_OUTPUT_DIR / "tables"
    csv = Path(tables_dir) / "daily_backtest_output.csv"
    if not csv.exists():
        return None

    df = pd.read_csv(csv, parse_dates=["date"]).set_index("date").sort_index()
    weight_cols = [f"w_{a}" for a in ASSET_NAMES]
    missing = [c for c in (weight_cols + ["pnl", "cum_pnl", "K", "regime",
                                           "max_p", "G", "turnover"])
               if c not in df.columns]
    if missing:
        raise ValueError(
            f"Cached daily_backtest_output.csv at {csv} is missing "
            f"expected columns: {missing}. Re-run the baseline."
        )

    weights = df[weight_cols].copy()
    weights.columns = ASSET_NAMES

    return BacktestResult(
        pnl          = df["pnl"].copy(),
        weights      = weights,
        cum_pnl      = df["cum_pnl"].copy(),
        K_history    = df["K"].rename("K"),
        tpl_label    = df["regime"].rename("regime"),
        tpl_max_prob = df["max_p"].rename("max_p"),
        tpl_count    = df["G"].rename("G"),
        turnover     = df["turnover"].rename("turnover"),
    )


def _banner(msg: str) -> None:
    print("\n" + "=" * 78)
    print(f"  {msg}")
    print("=" * 78)


def _headline(label: str, pnl: pd.Series, elapsed_min: float | None = None) -> None:
    if len(pnl) < 2 or pnl.std(ddof=1) == 0:
        print(f"  {label}:  Sharpe = n/a    MaxDD = n/a")
        return
    sharpe = (pnl.mean() / (pnl.std(ddof=1) + 1e-12)) * (252 ** 0.5)
    cum    = pnl.cumsum()
    mdd    = float((cum - cum.cummax()).min())
    extra  = f"    elapsed = {elapsed_min:.1f}m" if elapsed_min is not None else ""
    print(f"\n  {label}:  Sharpe = {sharpe:.2f}    MaxDD = {mdd:.2%}{extra}")


# --------------------------------------------------------------------- #
# Main entry
# --------------------------------------------------------------------- #
def main(
    macro_f_hmm_mult: int = 1,
    macro_M_candidates: tuple[int, ...] = (3, 4, 5),
    macro_transmat_prior_diag: float = 10.0,
    dirichlet_alpha: float = 1.0,
    market_n_iter: int = 15,
    frozen_emissions: bool = False,
    mean_anchor_strength: float = 100.0,
    freeze_covars: bool = True,
    force_baseline_rerun: bool = False,
) -> None:
    """Run (or load) baseline, then run hierarchical, then compare.

    Parameters
    ----------
    macro_f_hmm_mult : int, default 1
        Refit market HMM every `macro_f_hmm_mult * HMM.f_hmm` OOS days.
    macro_M_candidates : tuple, default (3, 4, 5)
        Candidate values of M for BIC selection.
    macro_transmat_prior_diag : float, default 10.0
        Diagonal-preference Dirichlet prior strength on macro HMM transitions.
    dirichlet_alpha : float, default 1.0
        Pseudo-count for market layer's transition smoothing.
    market_n_iter : int, default 15
        Max EM iterations per market HMM refit.
    frozen_emissions : bool, default False  ***v10.1 default***
        v10.0 had this True (K=G=5, template emissions held constant),
        which produced one-hot pG and constant MVO weights due to the
        L1 transaction cost. v10.1 default unfreezes emissions so K
        varies via BIC and pG is naturally mixed.
    mean_anchor_strength : float, default 100.0  ***NEW in v10.1***
        L2 anchor strength on the market HMM's emission means toward
        the baseline templates. Effective sample size of the prior:
        lam=100 with N_k~700 per regime means the prior is worth ~14%
        of the data. Without this anchor, macro-gated EM lets the K
        components drift toward the global mean (we diagnosed this in
        v10.0 unfrozen test runs), collapsing the W2 mapping onto 1-2
        templates. The anchor keeps components near their initial
        template locations so the W2 mapping produces diverse
        regime labels.
    freeze_covars : bool, default True  ***NEW in v10.1***
        Keep emission covariances at their template-init values rather
        than re-estimating each refit. Reduces parameter count and
        further stabilizes the regime structure. Companion to
        `mean_anchor_strength`.
    force_baseline_rerun : bool, default False
        Re-run baseline backtest even if cached.
    """
    ensure_dirs()
    t0_total = time.time()

    # ---------------------------------------------------------------- #
    # 1. Data
    # ---------------------------------------------------------------- #
    _banner("STEP 1 - Load data, compute returns, build features + macro")
    bundle = load_all()
    print(f"  prices     : {bundle['prices'].shape}")
    print(f"  returns    : {bundle['returns'].shape}")
    print(f"  features   : {bundle['features'].shape}")
    print(f"  train/test : {len(bundle['returns_train'])} / "
          f"{len(bundle['returns_test'])}  "
          f"(split @ {bundle['split_date'].date()})")

    macro_state = load_macro_state(feature_index=bundle["features"].index)
    print(f"  macro      : {macro_state.shape} "
          f"(NaN-free from {macro_state.dropna().index.min().date()})")

    # ---------------------------------------------------------------- #
    # 2. Baseline - load cache if present, else compute
    # ---------------------------------------------------------------- #
    cached = None if force_baseline_rerun else load_baseline_from_disk()
    if cached is not None:
        _banner("STEP 2 - BASELINE loaded from cache "
                "(Boukardagha/output/tables/daily_backtest_output.csv)")
        baseline_result = cached
        print(f"  Loaded {len(baseline_result.pnl)} OOS days "
              f"from {baseline_result.pnl.index[0].date()} "
              f"-> {baseline_result.pnl.index[-1].date()}")
        _headline("Baseline headline", baseline_result.pnl)
        print("\n  (To force re-run: set force_baseline_rerun=True.)")
    else:
        _banner("STEP 2 - BASELINE backtest (constant transition matrix)")
        t0 = time.time()
        baseline_result = run_backtest(
            returns       = bundle["returns"],
            returns_train = bundle["returns_train"],
            returns_test  = bundle["returns_test"],
            features_all  = bundle["features"],
            verbose       = RUN.verbose,
            progress_pct  = 10.0,
        )
        _headline("Baseline headline", baseline_result.pnl,
                  elapsed_min=(time.time() - t0) / 60)
        _banner("STEP 2b - Export baseline figures + tables (caches CSV)")
        export_all(baseline_result, bundle["returns"])

    # ---------------------------------------------------------------- #
    # 3. Hierarchical backtest
    # ---------------------------------------------------------------- #
    _banner("STEP 3 - HIERARCHICAL backtest (macro-regime-gated A_t)")
    print(f"  Settings: M_candidates={list(macro_M_candidates)}, "
          f"transmat_prior_diag={macro_transmat_prior_diag}, "
          f"dirichlet_alpha={dirichlet_alpha}")
    print(f"            market_n_iter={market_n_iter}, "
          f"frozen_emissions={frozen_emissions}, "
          f"mean_anchor_strength={mean_anchor_strength}, "
          f"freeze_covars={freeze_covars}, "
          f"f_hmm_mult={macro_f_hmm_mult}")
    t0 = time.time()
    hier_result, hier_diagnostics = run_backtest_hierarchical(
        returns                   = bundle["returns"],
        returns_train             = bundle["returns_train"],
        returns_test              = bundle["returns_test"],
        features_all              = bundle["features"],
        macro_state               = macro_state,
        verbose                   = RUN.verbose,
        progress_pct              = 5.0,
        macro_f_hmm_mult          = macro_f_hmm_mult,
        macro_M_candidates        = macro_M_candidates,
        macro_transmat_prior_diag = macro_transmat_prior_diag,
        dirichlet_alpha           = dirichlet_alpha,
        market_n_iter             = market_n_iter,
        frozen_emissions          = frozen_emissions,
        mean_anchor_strength      = mean_anchor_strength,
        freeze_covars             = freeze_covars,
    )
    _headline("Hierarchical headline", hier_result.pnl,
              elapsed_min=(time.time() - t0) / 60)
    print(f"  Refits: {hier_diagnostics['refit_count']}  "
          f"failures: {hier_diagnostics['refit_failures']}")
    if hier_diagnostics["macro_hmm"] is not None:
        print(f"  Macro M selected: {hier_diagnostics['macro_hmm'].M_}")

    # ---------------------------------------------------------------- #
    # 4. Comparison + charts
    # ---------------------------------------------------------------- #
    _banner("STEP 4 - A/B comparison: Baseline vs Hierarchical")

    spx_pnl = bundle["returns"]["SPX"].reindex(
        hier_result.pnl.index).fillna(0.0)

    run_full_comparison(
        baseline_result  = baseline_result,
        hier_result      = hier_result,
        hier_diagnostics = hier_diagnostics,
        spx_pnl          = spx_pnl,
    )

    _banner(f"DONE in {(time.time() - t0_total) / 60:.1f} minutes")


if __name__ == "__main__":
    main()
