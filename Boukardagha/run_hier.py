"""
run_hier.py
===========
Runs ONLY the hierarchical (macro + market) Wasserstein-HMM + MVO
backtest, saving raw outputs and the unified daily CSV.
"""
import os, sys, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

import pandas as pd

import config
from data import load_asset_returns, load_macro_panel, align_calendars
from hierarchical_hmm import run_hierarchical_strategy, make_daily_hier_csv

os.makedirs(config.RAW_DIR, exist_ok=True)
os.makedirs(config.OUTPUT_DIR, exist_ok=True)


def main():
    t0 = time.time()
    print("[hier] loading data...", flush=True)
    ar = load_asset_returns()
    mp = load_macro_panel()
    ar, mp = align_calendars(ar, mp)
    print(f"[hier] {len(ar)} rows, {ar.index[0].date()} -> {ar.index[-1].date()}",
          flush=True)

    print("[hier] running hierarchical Wasserstein HMM + MVO...", flush=True)
    out = run_hierarchical_strategy(ar, mp,
                                    oos_start=config.OOS_START, verbose=True)

    # Raw outputs (one CSV per series for incremental inspection)
    out["pnl"].to_csv(os.path.join(config.RAW_DIR, "hier_pnl.csv"))
    out["weights"].to_csv(os.path.join(config.RAW_DIR, "hier_weights.csv"))
    out["cum_pnl"].to_csv(os.path.join(config.RAW_DIR, "hier_cum_pnl.csv"))
    out["turnover"].to_csv(os.path.join(config.RAW_DIR, "hier_turnover.csv"))
    # Market layer
    out["K_history"].to_csv(os.path.join(config.RAW_DIR, "hier_K_market.csv"))
    out["tpl_label"].to_csv(os.path.join(config.RAW_DIR, "hier_market_label.csv"))
    out["tpl_count"].to_csv(os.path.join(config.RAW_DIR, "hier_G_market.csv"))
    out["tpl_prob"].to_csv(os.path.join(config.RAW_DIR, "hier_market_prob.csv"))
    # Macro layer
    out["macro_label"].to_csv(os.path.join(config.RAW_DIR, "hier_macro_label.csv"))
    out["macro_prob"].to_csv(os.path.join(config.RAW_DIR, "hier_macro_prob.csv"))
    out["K_macro_history"].to_csv(os.path.join(config.RAW_DIR, "hier_K_macro.csv"))
    out["G_macro_history"].to_csv(os.path.join(config.RAW_DIR, "hier_G_macro.csv"))

    csv_path = os.path.join(config.OUTPUT_DIR, "daily_backtest_output_hier.csv")
    make_daily_hier_csv(out, csv_path)
    print(f"[hier] wrote {csv_path}", flush=True)

    print(f"[hier] done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
