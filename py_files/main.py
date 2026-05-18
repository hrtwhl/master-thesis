"""
main.py
-------
Entry point for the Wasserstein-HMM regime-aware allocation pipeline.

Usage
~~~~~
    python main.py

What it does
~~~~~~~~~~~~
1. Load (or download) daily prices, compute log returns, precompute features.
2. Run the strictly causal expanding-window Wasserstein-HMM + MVO backtest.
3. Build SPX Buy & Hold, Equal-Weight, and 60/40 benchmarks aligned to the
   OOS window.
4. Export every table (CSV) and every figure (PNG) to `./output/`.

Everything is parameterized in `config.py` — modify there, not here.
"""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from backtest import run_backtest
from config import RUN, OUTPUT_DIR, ensure_dirs
from data import load_all
from reporting import export_all


def _banner(msg: str) -> None:
    print("\n" + "=" * 78)
    print(f"  {msg}")
    print("=" * 78)


def main() -> None:
    ensure_dirs()
    t0 = time.time()

    # ----------------------------------------------------- #
    # 1. Data
    # ----------------------------------------------------- #
    _banner("STEP 1 — Load data, compute returns, build features")
    bundle = load_all()
    print(f"  prices     : {bundle['prices'].shape}")
    print(f"  returns    : {bundle['returns'].shape}")
    print(f"  features   : {bundle['features'].shape}")
    print(f"  train/test : {len(bundle['returns_train'])} / "
          f"{len(bundle['returns_test'])}  "
          f"(split @ {bundle['split_date'].date()})")

    # ----------------------------------------------------- #
    # 2. Backtest
    # ----------------------------------------------------- #
    _banner("STEP 2 — Wasserstein-HMM + MVO backtest (strictly causal)")
    result = run_backtest(
        returns       = bundle["returns"],
        returns_train = bundle["returns_train"],
        returns_test  = bundle["returns_test"],
        features_all  = bundle["features"],
        verbose       = RUN.verbose,
        progress_pct  = 5.0,
    )

    # Quick summary so the operator sees the headline numbers immediately
    sharpe = (result.pnl.mean() / (result.pnl.std(ddof=1) + 1e-12)) * (252 ** 0.5)
    cum    = result.pnl.cumsum()
    mdd    = float((cum - cum.cummax()).min())
    print(f"\n  Headline:  Sharpe = {sharpe:.2f}    Max DD = {mdd:.2%}    "
          f"OOS days = {len(result.pnl)}")

    # ----------------------------------------------------- #
    # 3. Export everything (CSVs + PNGs)
    # ----------------------------------------------------- #
    _banner("STEP 3 — Export figures (PNG) and tables (CSV)")
    artefacts = export_all(result, bundle["returns"])

    print(f"\n  Wrote {len(artefacts['figures'])} figures and "
          f"{len(artefacts['tables']) + 1} tables to:")
    print(f"  {OUTPUT_DIR}")

    # Headline performance table to console
    print("\n--- Performance vs benchmarks ----------------------------------")
    print(artefacts["tables"]["performance"].round(3).to_string())

    print("\n--- Regime-conditioned portfolio performance --------------------")
    print(artefacts["tables"]["regime_perf"].round(3).to_string())

    _banner(f"DONE in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
