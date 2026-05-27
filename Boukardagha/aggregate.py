"""
aggregate.py
============
After run_pure.py and run_hier.py have populated output/raw/, this
script reads those CSVs and produces all paper-ready CSV tables and
PNG charts under output/tables/ and output/charts/.

Separated from main.py so the long backtest runs can be done once,
then aggregation/plotting can be re-run quickly while iterating.
"""
import os, sys, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import config
from data import load_asset_returns, load_macro_panel, align_calendars
from backtest import run_static_weight
from diagnostics import (
    performance_summary, turnover_summary, allocation_summary,
    concentration_summary, performance_by_regime,
    asset_performance_by_regime, regime_transition_matrix,
    regime_persistence,
)
from plots import (
    plot_strategy_bundle, plot_strategy_comparison,
    plot_drawdown_comparison, plot_rolling_sharpe,
    plot_annual_returns_bar, plot_underwater,
    plot_macro_diagnostics,
)

os.makedirs(config.CHART_DIR, exist_ok=True)
os.makedirs(config.TABLE_DIR, exist_ok=True)


def _csv(df, name, subdir=config.TABLE_DIR):
    path = os.path.join(subdir, f"{name}.csv")
    df.to_csv(path)
    return path


def _read_series(path: str, name: str) -> pd.Series:
    df = pd.read_csv(path, index_col=0, parse_dates=True)
    s = df.iloc[:, 0]
    s.name = name
    return s


def _read_df(path: str) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0, parse_dates=True)


def main():
    print("[agg] loading raw backtest outputs...", flush=True)
    R = config.RAW_DIR

    # Pure-market raw outputs
    pure_pnl   = _read_series(os.path.join(R, "pure_pnl.csv"),       "pnl")
    pure_W     = _read_df    (os.path.join(R, "pure_weights.csv"))
    pure_K     = _read_series(os.path.join(R, "pure_K_history.csv"), "K")
    pure_lbl   = _read_series(os.path.join(R, "pure_tpl_label.csv"), "tpl_label")
    pure_cnt   = _read_series(os.path.join(R, "pure_tpl_count.csv"), "tpl_count")
    pure_prob  = _read_series(os.path.join(R, "pure_tpl_prob.csv"),  "tpl_prob")

    # Hierarchical raw outputs
    hier_pnl     = _read_series(os.path.join(R, "hier_pnl.csv"),           "pnl")
    hier_W       = _read_df    (os.path.join(R, "hier_weights.csv"))
    hier_Km      = _read_series(os.path.join(R, "hier_K_market.csv"),      "K")
    hier_mlbl    = _read_series(os.path.join(R, "hier_market_label.csv"),  "market_lbl")
    hier_mcnt    = _read_series(os.path.join(R, "hier_G_market.csv"),      "market_cnt")
    hier_mprob   = _read_series(os.path.join(R, "hier_market_prob.csv"),   "market_prob")
    hier_macro_lbl  = _read_series(os.path.join(R, "hier_macro_label.csv"), "macro_lbl")
    hier_macro_prob = _read_series(os.path.join(R, "hier_macro_prob.csv"),  "macro_prob")
    hier_Kmac    = _read_series(os.path.join(R, "hier_K_macro.csv"),       "K_macro")
    hier_Gmac    = _read_series(os.path.join(R, "hier_G_macro.csv"),       "G_macro")

    # Asset + macro data
    print("[agg] loading data for benchmarks...", flush=True)
    asset_rets  = load_asset_returns()
    macro_panel = load_macro_panel()
    asset_rets, macro_panel = align_calendars(asset_rets, macro_panel)

    # Benchmarks
    print("[agg] running passive benchmarks...", flush=True)
    eqw = run_static_weight(asset_rets, config.EQUAL_WEIGHT,
                            oos_start=config.OOS_START, name="equal_weight")
    s60 = run_static_weight(asset_rets, config.SIXTY_FORTY,
                            oos_start=config.OOS_START, name="sixty_forty")
    eqw["pnl"].to_csv(os.path.join(R, "eqw_pnl.csv"))
    s60["pnl"].to_csv(os.path.join(R, "s60_pnl.csv"))

    # Common OOS window
    pnls = {
        "PureMarket_WHMM":   pure_pnl,
        "Hierarchical_WHMM": hier_pnl,
        "EqualWeight":       eqw["pnl"],
        "SixtyForty":        s60["pnl"],
    }
    common = pnls["PureMarket_WHMM"].index
    for s in pnls.values():
        common = common.intersection(s.index)
    pnls = {k: v.loc[common] for k, v in pnls.items()}
    print(f"[agg] common OOS index: {len(common)} days "
          f"({common[0].date()} -> {common[-1].date()})", flush=True)

    # ---------------- 1) Performance summary ----------------------
    perf = performance_summary(pnls)
    _csv(perf, "T01_performance_summary")
    print("\n=========== OOS PERFORMANCE SUMMARY ===========")
    print(perf.round(4).to_string(), flush=True)

    # ---------------- 2) Turnover / allocation / concentration ----
    weights_dict = {
        "PureMarket_WHMM":   pure_W.loc[common],
        "Hierarchical_WHMM": hier_W.loc[common],
        "EqualWeight":       eqw["weights"].loc[common],
        "SixtyForty":        s60["weights"].loc[common],
    }
    to_df = pd.DataFrame({k: turnover_summary(W) for k, W in weights_dict.items()}).T
    _csv(to_df, "T02_turnover_summary")

    _csv(allocation_summary(weights_dict["PureMarket_WHMM"]),
         "T03a_allocation_pure_market")
    _csv(allocation_summary(weights_dict["Hierarchical_WHMM"]),
         "T03b_allocation_hierarchical")

    conc = pd.DataFrame({k: concentration_summary(W) for k, W in weights_dict.items()}).T
    _csv(conc, "T04_concentration_summary")

    # ---------------- 3) Performance by regime ---------------------
    pure_by_reg = performance_by_regime(pure_pnl.loc[common], pure_lbl.loc[common])
    _csv(pure_by_reg, "T05a_pure_portfolio_by_regime")
    _csv(asset_performance_by_regime(asset_rets.loc[common], pure_lbl.loc[common]),
         "T06a_pure_asset_by_regime")

    hier_by_macro = performance_by_regime(hier_pnl.loc[common],
                                          hier_macro_lbl.loc[common])
    _csv(hier_by_macro, "T05b_hier_portfolio_by_macro_regime")
    _csv(asset_performance_by_regime(asset_rets.loc[common],
                                     hier_macro_lbl.loc[common]),
         "T06b_hier_asset_by_macro_regime")

    hier_by_market = performance_by_regime(hier_pnl.loc[common],
                                           hier_mlbl.loc[common])
    _csv(hier_by_market, "T05c_hier_portfolio_by_market_regime")

    # ---------------- 4) Regime transitions / persistence ----------
    _csv(regime_transition_matrix(pure_lbl.loc[common]),
         "T08a_pure_regime_transitions")
    _csv(regime_persistence(pure_lbl.loc[common]),
         "T09a_pure_regime_persistence")
    _csv(regime_transition_matrix(hier_macro_lbl.loc[common]),
         "T08b_hier_macro_regime_transitions")
    _csv(regime_persistence(hier_macro_lbl.loc[common]),
         "T09b_hier_macro_regime_persistence")
    _csv(regime_transition_matrix(hier_mlbl.loc[common]),
         "T08c_hier_market_regime_transitions")
    _csv(regime_persistence(hier_mlbl.loc[common]),
         "T09c_hier_market_regime_persistence")

    # ---------------- 5) Crisis window stress test -----------------
    crises = {
        "GFC_2008":        ("2008-01-01", "2009-06-30"),
        "Euro_2011":       ("2011-06-01", "2012-06-30"),
        "Taper_2013":      ("2013-04-01", "2013-12-31"),
        "China_2015":      ("2015-06-01", "2016-02-29"),
        "Q4_2018":         ("2018-09-01", "2019-01-31"),
        "Covid_2020":      ("2020-02-01", "2020-06-30"),
        "Rates_2022":      ("2022-01-01", "2022-12-31"),
        "Liberation_2025": ("2025-01-01", "2025-09-30"),
    }
    rows = []
    for label, (s, e) in crises.items():
        s_ts, e_ts = pd.Timestamp(s), pd.Timestamp(e)
        for strat, pnl in pnls.items():
            p = pnl.loc[(pnl.index >= s_ts) & (pnl.index <= e_ts)]
            if len(p) == 0:
                continue
            wealth = (1.0 + p.fillna(0.0)).cumprod()
            dd = (wealth / wealth.cummax() - 1.0).min()
            rows.append({
                "Crisis":        label,
                "Start":         s,
                "End":           e,
                "Days":          int(len(p)),
                "Strategy":      strat,
                "Total log ret": float(p.sum()),
                "Ann Sharpe":    float(np.sqrt(252) * p.mean() / (p.std() + 1e-12)),
                "Max DD":        float(dd),
            })
    crisis_df = pd.DataFrame(rows)
    if not crisis_df.empty:
        _csv(crisis_df.set_index(["Crisis", "Strategy"]),
             "T07_crisis_performance")

    # ---------------- 6) Per-strategy chart bundles ----------------
    print("[agg] plotting pure-market bundle...", flush=True)
    plot_strategy_bundle(
        pnl=pure_pnl.loc[common],
        weights=pure_W.loc[common],
        regime=pure_lbl.loc[common],
        K_history=pure_K.loc[common],
        G_history=pure_cnt.loc[common],
        tpl_prob=pure_prob.loc[common],
        asset_returns=asset_rets.loc[common],
        prefix="pure_",
    )

    print("[agg] plotting hierarchical (market view) bundle...", flush=True)
    plot_strategy_bundle(
        pnl=hier_pnl.loc[common],
        weights=hier_W.loc[common],
        regime=hier_mlbl.loc[common],
        K_history=hier_Km.loc[common],
        G_history=hier_mcnt.loc[common],
        tpl_prob=hier_mprob.loc[common],
        asset_returns=asset_rets.loc[common],
        prefix="hier_",
    )

    print("[agg] plotting hierarchical macro-layer diagnostics...", flush=True)
    plot_macro_diagnostics(
        macro_label=hier_macro_lbl.loc[common],
        macro_prob=hier_macro_prob.loc[common],
        market_label=hier_mlbl.loc[common],
        K_macro=hier_Kmac.loc[common],
        G_macro=hier_Gmac.loc[common],
        asset_returns=asset_rets.loc[common],
    )

    # ---------------- 7) Cross-strategy comparison charts ---------
    print("[agg] cross-strategy comparison charts...", flush=True)
    plot_strategy_comparison(pnls, name="11_strategy_comparison")
    plot_drawdown_comparison(pnls, name="12_drawdown_comparison")
    plot_rolling_sharpe(pnls, window=252, name="13_rolling_sharpe_1y")
    plot_annual_returns_bar(pnls, name="14_annual_returns")
    plot_underwater(pnls, name="15_underwater")

    # ---------------- 8) Unified PnL panel ------------------------
    pnl_panel = pd.DataFrame(pnls).sort_index()
    pnl_panel.to_csv(os.path.join(R, "pnl_panel.csv"))
    pnl_panel.cumsum().to_csv(os.path.join(R, "cum_pnl_panel.csv"))

    # ---------------- 9) Benchmark daily CSVs --------------------
    bench_panel = pd.concat([
        eqw["weights"].add_prefix("w_").loc[common],
        eqw["pnl"].loc[common].rename("pnl_equal_weight"),
        s60["weights"].add_prefix("w_60_40_").loc[common],
        s60["pnl"].loc[common].rename("pnl_sixty_forty"),
    ], axis=1)
    bench_panel.to_csv(os.path.join(config.OUTPUT_DIR,
                                    "daily_backtest_output_benchmarks.csv"))

    print("[agg] done.", flush=True)
    print(f"  Charts: {config.CHART_DIR}")
    print(f"  Tables: {config.TABLE_DIR}")
    print(f"  Raw:    {config.RAW_DIR}")


if __name__ == "__main__":
    main()
