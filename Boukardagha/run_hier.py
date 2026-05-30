"""
run_hier.py
===========
Runs BOTH hierarchical strategies and saves raw outputs + unified
daily CSVs:

  Hierarchical B : macro x market joint-mixture moments, no tilt.
  Hierarchical C : macro as risk modulator (Fix C), market owns mu.

Output CSVs (in output/):
    daily_backtest_output_hierB.csv
    daily_backtest_output_hierC.csv
"""
import os, sys, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

import pandas as pd

import config
from data import load_asset_returns, load_macro_panel, align_calendars
from hierarchical_hmm import (
    run_hierarchical_strategy_B, run_hierarchical_strategy_C,
    make_daily_hier_csv, make_daily_hierC_csv,
)

os.makedirs(config.RAW_DIR, exist_ok=True)
os.makedirs(config.OUTPUT_DIR, exist_ok=True)


def _save_common(out, tag):
    """Save the raw series shared by both hierarchical variants."""
    R = config.RAW_DIR
    out["pnl"].to_csv(os.path.join(R, f"{tag}_pnl.csv"))
    out["weights"].to_csv(os.path.join(R, f"{tag}_weights.csv"))
    out["cum_pnl"].to_csv(os.path.join(R, f"{tag}_cum_pnl.csv"))
    out["turnover"].to_csv(os.path.join(R, f"{tag}_turnover.csv"))
    out["K_history"].to_csv(os.path.join(R, f"{tag}_K_market.csv"))
    out["tpl_label"].to_csv(os.path.join(R, f"{tag}_market_label.csv"))
    out["tpl_count"].to_csv(os.path.join(R, f"{tag}_G_market.csv"))
    out["tpl_prob"].to_csv(os.path.join(R, f"{tag}_market_prob.csv"))
    out["macro_label"].to_csv(os.path.join(R, f"{tag}_macro_label.csv"))
    out["macro_prob"].to_csv(os.path.join(R, f"{tag}_macro_prob.csv"))
    out["K_macro_history"].to_csv(os.path.join(R, f"{tag}_K_macro.csv"))
    out["G_macro_history"].to_csv(os.path.join(R, f"{tag}_G_macro.csv"))


def main():
    t0 = time.time()
    print("[hier] loading data...", flush=True)
    ar = load_asset_returns()
    mp = load_macro_panel()
    ar, mp = align_calendars(ar, mp)
    print(f"[hier] {len(ar)} rows, {ar.index[0].date()} -> {ar.index[-1].date()}",
          flush=True)

    # ---- Hierarchical B -------------------------------------------------
    print("\n[hier] running Hierarchical B (joint mixture, no tilt)...", flush=True)
    t1 = time.time()
    outB = run_hierarchical_strategy_B(ar, mp, oos_start=config.OOS_START,
                                       verbose=True)
    _save_common(outB, "hierB")
    csvB = os.path.join(config.OUTPUT_DIR, "daily_backtest_output_hierB.csv")
    make_daily_hier_csv(outB, csvB)
    print(f"[hier] B done in {time.time()-t1:.1f}s -> {csvB}", flush=True)

    # ---- Hierarchical C -------------------------------------------------
    print("\n[hier] running Hierarchical C (macro risk modulator, Fix C)...",
          flush=True)
    t1 = time.time()
    outC = run_hierarchical_strategy_C(ar, mp, oos_start=config.OOS_START,
                                       verbose=True)
    _save_common(outC, "hierC")
    # C-specific extras
    outC["macro_stress"].to_csv(os.path.join(config.RAW_DIR, "hierC_macro_stress.csv"))
    outC["gamma_eff"].to_csv(os.path.join(config.RAW_DIR, "hierC_gamma_eff.csv"))
    csvC = os.path.join(config.OUTPUT_DIR, "daily_backtest_output_hierC.csv")
    make_daily_hierC_csv(outC, csvC)
    print(f"[hier] C done in {time.time()-t1:.1f}s -> {csvC}", flush=True)

    print(f"\n[hier] all done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
