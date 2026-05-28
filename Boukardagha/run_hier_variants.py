"""
run_hier_variants.py
====================
Runs the hierarchical strategy under all three macro-tilt modes
('off', 'symmetric', 'asymmetric') and saves them as three separate
result files for direct comparison.

This is what we recommend after the first OOS run revealed that the
*symmetric* tilt was actively hurting performance in 2022 by holding
too much SPX during the rates-driven equity sell-off.

Output CSVs (in output/):
    daily_backtest_output_hier_off.csv         (FIX A - tilt disabled)
    daily_backtest_output_hier_symmetric.csv   (original design)
    daily_backtest_output_hier_asymmetric.csv  (FIX B - new design)
"""
import os, sys, time, warnings, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

import pandas as pd
import config
from data import load_asset_returns, load_macro_panel, align_calendars
import hierarchical_hmm   # imported as module so we can reload it

os.makedirs(config.RAW_DIR, exist_ok=True)
os.makedirs(config.OUTPUT_DIR, exist_ok=True)


def main():
    print("[hier-var] loading data...", flush=True)
    ar = load_asset_returns()
    mp = load_macro_panel()
    ar, mp = align_calendars(ar, mp)
    print(f"[hier-var] {len(ar)} rows, "
          f"{ar.index[0].date()} -> {ar.index[-1].date()}", flush=True)

    for mode in ("off", "symmetric", "asymmetric"):
        print(f"\n--- HIER variant: tilt mode = '{mode}' ---", flush=True)
        config.HIER_MACRO_TILT_MODE = mode
        importlib.reload(hierarchical_hmm)   # rebind the module-level closure
        from hierarchical_hmm import (
            run_hierarchical_strategy, make_daily_hier_csv,
        )

        t0 = time.time()
        out = run_hierarchical_strategy(ar, mp,
                                        oos_start=config.OOS_START,
                                        verbose=True)
        print(f"[hier-var] '{mode}' done in {time.time()-t0:.1f}s",
              flush=True)

        suffix = mode
        out["pnl"].to_csv(os.path.join(
            config.RAW_DIR, f"hier_{suffix}_pnl.csv"))
        out["weights"].to_csv(os.path.join(
            config.RAW_DIR, f"hier_{suffix}_weights.csv"))
        out["macro_label"].to_csv(os.path.join(
            config.RAW_DIR, f"hier_{suffix}_macro_label.csv"))

        csv_path = os.path.join(
            config.OUTPUT_DIR,
            f"daily_backtest_output_hier_{suffix}.csv",
        )
        make_daily_hier_csv(out, csv_path)
        print(f"[hier-var] wrote {csv_path}", flush=True)


if __name__ == "__main__":
    main()
