"""
run_pure.py
===========
Runs ONLY the pure-market Wasserstein-HMM + MVO backtest, saving:
  - raw daily PnL / weights / regime time series  (raw/)
  - the unified daily_backtest_output.csv         (output/)
"""
import os, sys, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

import pandas as pd

import config
from data import load_asset_returns, load_macro_panel, align_calendars
from backtest import run_pure_market_strategy, make_daily_backtest_csv

os.makedirs(config.RAW_DIR, exist_ok=True)
os.makedirs(config.OUTPUT_DIR, exist_ok=True)


def main():
    t0 = time.time()
    print("[pure] loading data...", flush=True)
    ar = load_asset_returns()
    mp = load_macro_panel()
    ar, mp = align_calendars(ar, mp)
    print(f"[pure] {len(ar)} rows, {ar.index[0].date()} -> {ar.index[-1].date()}",
          flush=True)
    ar.to_csv(os.path.join(config.RAW_DIR, "asset_log_returns.csv"))
    mp.to_csv(os.path.join(config.RAW_DIR, "macro_levels.csv"))

    print("[pure] running pure-market Wasserstein HMM + MVO...", flush=True)
    out = run_pure_market_strategy(ar, oos_start=config.OOS_START, verbose=True)

    out["pnl"].to_csv(os.path.join(config.RAW_DIR, "pure_pnl.csv"))
    out["weights"].to_csv(os.path.join(config.RAW_DIR, "pure_weights.csv"))
    out["cum_pnl"].to_csv(os.path.join(config.RAW_DIR, "pure_cum_pnl.csv"))
    out["turnover"].to_csv(os.path.join(config.RAW_DIR, "pure_turnover.csv"))
    out["K_history"].to_csv(os.path.join(config.RAW_DIR, "pure_K_history.csv"))
    out["tpl_label"].to_csv(os.path.join(config.RAW_DIR, "pure_tpl_label.csv"))
    out["tpl_count"].to_csv(os.path.join(config.RAW_DIR, "pure_tpl_count.csv"))
    out["tpl_prob"].to_csv(os.path.join(config.RAW_DIR, "pure_tpl_prob.csv"))

    # Unified daily output CSV (user-requested schema)
    csv_path = os.path.join(config.OUTPUT_DIR, "daily_backtest_output_pure.csv")
    make_daily_backtest_csv(out, csv_path)
    print(f"[pure] wrote {csv_path}", flush=True)

    print(f"[pure] done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
