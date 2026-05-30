"""
robustness_nooverlap.py  (Robustness R2)
========================================
'stocks' (=SPX) and 'oil' are in BOTH the asset universe and the macro
panel.  Including tradeable assets in the macro panel makes the macro
and market layers informationally dependent: a skeptic can argue that
Hierarchical C's "macro stress" is partly just "SPX/oil are volatile
right now" -- endogenous market information the market WHMM already sees.

R2 re-runs Hierarchical B and C with a macro panel that EXCLUDES the
overlapping variables (config.R2_MACRO_VARS = copper, yield_curve,
stock_bond_corr, vix, us3mo), leaving only genuinely exogenous
macro-financial series.  If C still beats Pure with this non-overlapping
panel, the result is not driven by re-using tradeable-asset information.

Outputs (in output/):
    daily_backtest_output_hierB_nooverlap.csv
    daily_backtest_output_hierC_nooverlap.csv
    tables/R2_nooverlap_performance.csv   (vs pure & full-panel B/C)
    charts/R2_nooverlap_comparison.png
"""
import os, sys, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config
from data import load_asset_returns, load_macro_panel, align_calendars
from hierarchical_hmm import (
    run_hierarchical_strategy_B, run_hierarchical_strategy_C,
    make_daily_hier_csv, make_daily_hierC_csv,
)
from diagnostics import performance_summary, turnover_from_weights

os.makedirs(config.RAW_DIR, exist_ok=True)
os.makedirs(config.TABLE_DIR, exist_ok=True)
os.makedirs(config.CHART_DIR, exist_ok=True)


def _metrics(pnl, W):
    from diagnostics import (ann_sharpe, ann_sortino, ann_mean, ann_vol,
                             max_drawdown_from_returns, calmar_ratio, hit_rate)
    to = turnover_from_weights(W).mean()
    return {
        "Sharpe":  ann_sharpe(pnl), "Sortino": ann_sortino(pnl),
        "Ann Mean": ann_mean(pnl), "Ann Vol": ann_vol(pnl),
        "Max DD":  max_drawdown_from_returns(pnl),
        "Calmar":  calmar_ratio(pnl), "Hit Rate": hit_rate(pnl),
        "Turnover": float(to),
    }


def main():
    t0 = time.time()
    print("[R2] loading data...", flush=True)
    ar = load_asset_returns()
    mp_full = load_macro_panel()
    ar, mp_full = align_calendars(ar, mp_full)

    # No-overlap macro panel
    mp_no = mp_full[config.R2_MACRO_VARS].copy()
    print(f"[R2] macro panel reduced from {list(mp_full.columns)} "
          f"to {list(mp_no.columns)}", flush=True)

    # ---- Hierarchical B (no-overlap) -----------------------------------
    print("\n[R2] Hierarchical B (no-overlap macro panel)...", flush=True)
    t1 = time.time()
    outB = run_hierarchical_strategy_B(ar, mp_no, oos_start=config.OOS_START,
                                       verbose=True)
    make_daily_hier_csv(outB, os.path.join(
        config.OUTPUT_DIR, "daily_backtest_output_hierB_nooverlap.csv"))
    print(f"[R2] B done in {time.time()-t1:.1f}s", flush=True)

    # ---- Hierarchical C (no-overlap) -----------------------------------
    print("\n[R2] Hierarchical C (no-overlap macro panel)...", flush=True)
    t1 = time.time()
    outC = run_hierarchical_strategy_C(ar, mp_no, oos_start=config.OOS_START,
                                       verbose=True)
    make_daily_hierC_csv(outC, os.path.join(
        config.OUTPUT_DIR, "daily_backtest_output_hierC_nooverlap.csv"))
    print(f"[R2] C done in {time.time()-t1:.1f}s", flush=True)

    # ---- Performance comparison ----------------------------------------
    # Pull in the full-panel results if their raw PnL exists, plus pure.
    R = config.RAW_DIR
    pnls, weights = {}, {}

    def _try_pnl(path):
        if os.path.exists(path):
            return pd.read_csv(path, index_col=0, parse_dates=True).iloc[:, 0]
        return None

    pure_pnl = _try_pnl(os.path.join(R, "pure_pnl.csv"))
    if pure_pnl is not None:
        pnls["Pure"] = pure_pnl
        weights["Pure"] = pd.read_csv(os.path.join(R, "pure_weights.csv"),
                                      index_col=0, parse_dates=True)
    fullB = _try_pnl(os.path.join(R, "hierB_pnl.csv"))
    fullC = _try_pnl(os.path.join(R, "hierC_pnl.csv"))
    if fullB is not None:
        pnls["HierB_full"] = fullB
        weights["HierB_full"] = pd.read_csv(os.path.join(R, "hierB_weights.csv"),
                                            index_col=0, parse_dates=True)
    if fullC is not None:
        pnls["HierC_full"] = fullC
        weights["HierC_full"] = pd.read_csv(os.path.join(R, "hierC_weights.csv"),
                                            index_col=0, parse_dates=True)

    pnls["HierB_nooverlap"] = outB["pnl"]
    pnls["HierC_nooverlap"] = outC["pnl"]
    weights["HierB_nooverlap"] = outB["weights"]
    weights["HierC_nooverlap"] = outC["weights"]

    # Common window
    common = None
    for s in pnls.values():
        common = s.index if common is None else common.intersection(s.index)
    rows = {k: _metrics(v.loc[common], weights[k].loc[common])
            for k, v in pnls.items()}
    perf = pd.DataFrame(rows).T
    perf.to_csv(os.path.join(config.TABLE_DIR, "R2_nooverlap_performance.csv"))
    print("\n=========== R2 NO-OVERLAP PERFORMANCE ===========")
    print(perf.round(4).to_string(), flush=True)

    # ---- Comparison chart ---------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 6))
    for label, s in pnls.items():
        cum = s.loc[common].cumsum()
        cum = cum - cum.iloc[0]
        style = "--" if "nooverlap" in label else "-"
        ax.plot(cum.index, cum.values, style, linewidth=1.4, label=label)
    ax.set_title("R2: Cumulative log return — full vs no-overlap macro panel")
    ax.set_xlabel("Date"); ax.set_ylabel("Cumulative log return")
    ax.grid(True); ax.legend()
    fig.savefig(os.path.join(config.CHART_DIR, "R2_nooverlap_comparison.png"),
                bbox_inches="tight")
    plt.close(fig)

    print(f"\n[R2] done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
