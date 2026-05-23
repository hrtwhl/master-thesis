"""
main_macro_tilt.py
------------------
Entry point for v10.2 — baseline market HMM + macro-tilted MVO.

Usage
~~~~~
    python main_macro_tilt.py

What it does
~~~~~~~~~~~~
1. Loads daily asset prices + 7-variable macro state vector.
2. Loads baseline backtest from cache if available, else runs it.
3. Runs v10.2 backtest:
       * Calibration: fit macro HMM (BIC over {3,4,5}), then compute
         per-macro-regime return-tilt vectors on calibration sample.
       * OOS: baseline market HMM (unchanged) for regime detection;
         macro Viterbi for m_t; tilt MVO inputs by alpha_mu * tilt[m_t].
4. Comparison vs baseline (charts + CSVs in output/compare/).

Architecture
~~~~~~~~~~~~
v10.2 abandons the macro-gated transition matrix approach (v10/v10.1)
which interacted badly with long-horizon EMA template drift and
produced flat weights. Instead, the macro signal enters cleanly at the
MVO step as a frozen-after-calibration additive return tilt.

This separation has several advantages:
  * Baseline regime detection is preserved unchanged — diverse weights
    are guaranteed (baseline already produces them).
  * Macro contribution is a single hyperparameter `alpha_mu` whose
    value-add is cleanly measurable as the performance delta from 0.
  * No long-horizon stability issues: tilts are frozen, no drift.
  * Strictly causal: tilts use only pre-OOS calibration data.

Hyperparameters
~~~~~~~~~~~~~~~
alpha_mu     (default 1.0)  : strength of mean tilt
alpha_sig    (default 0.0)  : strength of vol-scale tilt (off by default)
tilt_shrink  (default 0.5)  : Bayesian shrink of regime mean → uncond mean
"""

from __future__ import annotations

import _paths  # noqa: F401

import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from config import RUN, OUTPUT_DIR as BASELINE_OUTPUT_DIR, ASSET_NAMES, ensure_dirs
from data import load_all
from data_macro import load_macro_state
from backtest import run_backtest, BacktestResult
from backtest_macro_tilt import run_backtest_macro_tilt
from compare_baseline_vs_v10p2 import run_full_comparison
from reporting import export_all


def load_baseline_from_disk(tables_dir: Path | None = None) -> Optional[BacktestResult]:
    """Reconstruct a BacktestResult from the cached baseline daily CSV."""
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
            f"Cached daily_backtest_output.csv at {csv} is missing columns: {missing}."
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


def main(
    macro_M_candidates: tuple[int, ...] = (3, 4, 5),
    macro_transmat_prior_diag: float = 10.0,
    alpha_mu: float = 1.0,
    alpha_sig: float = 0.0,
    tilt_shrink: float = 0.5,
    tilt_min_obs: int = 50,
    force_baseline_rerun: bool = False,
) -> None:
    """Run (or load) baseline, then run v10.2 macro-tilted, then compare.

    Parameters
    ----------
    macro_M_candidates : tuple, default (3, 4, 5)
        Candidate values of M (macro regimes) for BIC selection.
    macro_transmat_prior_diag : float, default 10.0
        Diagonal-preference Dirichlet prior strength on macro HMM transitions.
    alpha_mu : float, default 1.0
        Strength of macro mean tilt at MVO step. 0 = baseline (no tilt).
        1 = full calibration-time tilt magnitude.
    alpha_sig : float, default 0.0
        Strength of macro vol-scale adjustment. 0 = leave Sigma alone.
    tilt_shrink : float in [0, 1], default 0.5
        Bayesian shrinkage of regime means toward unconditional mean.
        0 = raw (high variance, overfits small regimes); 1 = no tilt.
    tilt_min_obs : int, default 50
        Minimum calibration days for a regime to get its own tilt.
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
    # 2. Baseline (cached if available)
    # ---------------------------------------------------------------- #
    cached = None if force_baseline_rerun else load_baseline_from_disk()
    if cached is not None:
        _banner("STEP 2 - BASELINE loaded from cache")
        baseline_result = cached
        print(f"  Loaded {len(baseline_result.pnl)} OOS days "
              f"from {baseline_result.pnl.index[0].date()} "
              f"-> {baseline_result.pnl.index[-1].date()}")
        _headline("Baseline headline", baseline_result.pnl)
    else:
        _banner("STEP 2 - BASELINE backtest (constant A, no macro tilt)")
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
    # 3. v10.2 backtest
    # ---------------------------------------------------------------- #
    _banner("STEP 3 - v10.2 backtest (baseline market HMM + macro tilt at MVO)")
    print(f"  Macro config: M_candidates={list(macro_M_candidates)}, "
          f"transmat_prior_diag={macro_transmat_prior_diag}")
    print(f"  Tilt config:  alpha_mu={alpha_mu}, alpha_sig={alpha_sig}, "
          f"tilt_shrink={tilt_shrink}, tilt_min_obs={tilt_min_obs}")
    t0 = time.time()
    hier_result, hier_diagnostics = run_backtest_macro_tilt(
        returns                   = bundle["returns"],
        returns_train             = bundle["returns_train"],
        returns_test              = bundle["returns_test"],
        features_all              = bundle["features"],
        macro_state               = macro_state,
        verbose                   = RUN.verbose,
        progress_pct              = 5.0,
        macro_M_candidates        = macro_M_candidates,
        macro_transmat_prior_diag = macro_transmat_prior_diag,
        alpha_mu                  = alpha_mu,
        alpha_sig                 = alpha_sig,
        tilt_shrink               = tilt_shrink,
        tilt_min_obs              = tilt_min_obs,
    )
    _headline("v10.2 headline", hier_result.pnl,
              elapsed_min=(time.time() - t0) / 60)
    if hier_diagnostics["macro_hmm"] is not None:
        print(f"  Macro M selected: {hier_diagnostics['macro_hmm'].M_}")

    # ---------------------------------------------------------------- #
    # 4. Comparison
    # ---------------------------------------------------------------- #
    _banner("STEP 4 - A/B comparison: Baseline vs v10.2")
    spx_pnl = bundle["returns"]["SPX"].reindex(hier_result.pnl.index).fillna(0.0)
    run_full_comparison(
        baseline_result  = baseline_result,
        hier_result      = hier_result,
        hier_diagnostics = hier_diagnostics,
        spx_pnl          = spx_pnl,
    )

    _banner(f"DONE in {(time.time() - t0_total) / 60:.1f} minutes")


if __name__ == "__main__":
    main()
