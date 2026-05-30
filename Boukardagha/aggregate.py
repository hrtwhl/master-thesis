"""
aggregate.py
============
Reads the raw CSVs produced by run_pure.py and run_hier.py and builds
all paper-ready CSV tables and PNG charts.

Strategies compared:
    PureMarket_WHMM     - Boukardagha (2026) replication
    Hierarchical_B      - macro x market joint mixture, no tilt
    Hierarchical_C      - macro as risk modulator (Fix C)
    EqualWeight, SixtyForty - passive benchmarks
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
    plot_macro_diagnostics, plot_macro_stress,
)

os.makedirs(config.CHART_DIR, exist_ok=True)
os.makedirs(config.TABLE_DIR, exist_ok=True)


def _csv(df, name, subdir=config.TABLE_DIR):
    path = os.path.join(subdir, f"{name}.csv")
    df.to_csv(path)
    return path


def _rs(path, name):
    df = pd.read_csv(path, index_col=0, parse_dates=True)
    s = df.iloc[:, 0]; s.name = name
    return s


def _rd(path):
    return pd.read_csv(path, index_col=0, parse_dates=True)


def _load_hier(tag):
    """Load the raw bundle for a hierarchical variant (tag = hierB/hierC)."""
    R = config.RAW_DIR
    d = dict(
        pnl       = _rs(os.path.join(R, f"{tag}_pnl.csv"),          "pnl"),
        weights   = _rd(os.path.join(R, f"{tag}_weights.csv")),
        K_market  = _rs(os.path.join(R, f"{tag}_K_market.csv"),     "K"),
        mkt_lbl   = _rs(os.path.join(R, f"{tag}_market_label.csv"), "mkt_lbl"),
        mkt_cnt   = _rs(os.path.join(R, f"{tag}_G_market.csv"),     "mkt_cnt"),
        mkt_prob  = _rs(os.path.join(R, f"{tag}_market_prob.csv"),  "mkt_prob"),
        macro_lbl = _rs(os.path.join(R, f"{tag}_macro_label.csv"),  "macro_lbl"),
        macro_prob= _rs(os.path.join(R, f"{tag}_macro_prob.csv"),   "macro_prob"),
        K_macro   = _rs(os.path.join(R, f"{tag}_K_macro.csv"),      "K_macro"),
        G_macro   = _rs(os.path.join(R, f"{tag}_G_macro.csv"),      "G_macro"),
    )
    return d


def main():
    print("[agg] loading raw backtest outputs...", flush=True)
    R = config.RAW_DIR

    # Pure-market
    pure = dict(
        pnl   = _rs(os.path.join(R, "pure_pnl.csv"),       "pnl"),
        W     = _rd(os.path.join(R, "pure_weights.csv")),
        K     = _rs(os.path.join(R, "pure_K_history.csv"), "K"),
        lbl   = _rs(os.path.join(R, "pure_tpl_label.csv"), "tpl_label"),
        cnt   = _rs(os.path.join(R, "pure_tpl_count.csv"), "tpl_count"),
        prob  = _rs(os.path.join(R, "pure_tpl_prob.csv"),  "tpl_prob"),
    )

    hierB = _load_hier("hierB")
    hierC = _load_hier("hierC")
    # C-specific extras
    hierC_stress = _rs(os.path.join(R, "hierC_macro_stress.csv"), "macro_stress")
    hierC_gamma  = _rs(os.path.join(R, "hierC_gamma_eff.csv"),    "gamma_eff")

    print("[agg] loading data for benchmarks...", flush=True)
    asset_rets  = load_asset_returns()
    macro_panel = load_macro_panel()
    asset_rets, macro_panel = align_calendars(asset_rets, macro_panel)

    print("[agg] running passive benchmarks...", flush=True)
    eqw = run_static_weight(asset_rets, config.EQUAL_WEIGHT,
                            oos_start=config.OOS_START, name="equal_weight")
    s60 = run_static_weight(asset_rets, config.SIXTY_FORTY,
                            oos_start=config.OOS_START, name="sixty_forty")
    eqw["pnl"].to_csv(os.path.join(R, "eqw_pnl.csv"))
    s60["pnl"].to_csv(os.path.join(R, "s60_pnl.csv"))

    # Common OOS window across all strategies
    pnls = {
        "PureMarket_WHMM":  pure["pnl"],
        "Hierarchical_B":   hierB["pnl"],
        "Hierarchical_C":   hierC["pnl"],
        "EqualWeight":      eqw["pnl"],
        "SixtyForty":       s60["pnl"],
    }
    common = pnls["PureMarket_WHMM"].index
    for s in pnls.values():
        common = common.intersection(s.index)
    pnls = {k: v.loc[common] for k, v in pnls.items()}
    print(f"[agg] common OOS index: {len(common)} days "
          f"({common[0].date()} -> {common[-1].date()})", flush=True)

    # ---- 1) Performance summary --------------------------------------
    perf = performance_summary(pnls)
    _csv(perf, "T01_performance_summary")
    print("\n=========== OOS PERFORMANCE SUMMARY ===========")
    print(perf.round(4).to_string(), flush=True)

    # ---- 2) Turnover / allocation / concentration --------------------
    weights_dict = {
        "PureMarket_WHMM": pure["W"].loc[common],
        "Hierarchical_B":  hierB["weights"].loc[common],
        "Hierarchical_C":  hierC["weights"].loc[common],
        "EqualWeight":     eqw["weights"].loc[common],
        "SixtyForty":      s60["weights"].loc[common],
    }
    _csv(pd.DataFrame({k: turnover_summary(W) for k, W in weights_dict.items()}).T,
         "T02_turnover_summary")
    _csv(allocation_summary(weights_dict["PureMarket_WHMM"]), "T03a_allocation_pure")
    _csv(allocation_summary(weights_dict["Hierarchical_B"]),  "T03b_allocation_hierB")
    _csv(allocation_summary(weights_dict["Hierarchical_C"]),  "T03c_allocation_hierC")
    _csv(pd.DataFrame({k: concentration_summary(W) for k, W in weights_dict.items()}).T,
         "T04_concentration_summary")

    # ---- 3) Performance by regime ------------------------------------
    _csv(performance_by_regime(pure["pnl"].loc[common], pure["lbl"].loc[common]),
         "T05a_pure_by_regime")
    _csv(asset_performance_by_regime(asset_rets.loc[common], pure["lbl"].loc[common]),
         "T06a_pure_asset_by_regime")

    _csv(performance_by_regime(hierB["pnl"].loc[common], hierB["macro_lbl"].loc[common]),
         "T05b_hierB_by_macro_regime")
    _csv(performance_by_regime(hierB["pnl"].loc[common], hierB["mkt_lbl"].loc[common]),
         "T05b_hierB_by_market_regime")
    _csv(asset_performance_by_regime(asset_rets.loc[common], hierB["macro_lbl"].loc[common]),
         "T06b_hierB_asset_by_macro_regime")

    _csv(performance_by_regime(hierC["pnl"].loc[common], hierC["macro_lbl"].loc[common]),
         "T05c_hierC_by_macro_regime")
    _csv(performance_by_regime(hierC["pnl"].loc[common], hierC["mkt_lbl"].loc[common]),
         "T05c_hierC_by_market_regime")
    _csv(asset_performance_by_regime(asset_rets.loc[common], hierC["macro_lbl"].loc[common]),
         "T06c_hierC_asset_by_macro_regime")

    # ---- 4) Regime transitions / persistence -------------------------
    _csv(regime_transition_matrix(pure["lbl"].loc[common]),    "T08a_pure_transitions")
    _csv(regime_persistence(pure["lbl"].loc[common]),          "T09a_pure_persistence")
    _csv(regime_transition_matrix(hierB["macro_lbl"].loc[common]), "T08b_hierB_macro_transitions")
    _csv(regime_persistence(hierB["macro_lbl"].loc[common]),       "T09b_hierB_macro_persistence")
    _csv(regime_transition_matrix(hierC["macro_lbl"].loc[common]), "T08c_hierC_macro_transitions")
    _csv(regime_persistence(hierC["macro_lbl"].loc[common]),       "T09c_hierC_macro_persistence")

    # ---- 5) Crisis-window stress test --------------------------------
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
                "Crisis": label, "Start": s, "End": e, "Days": int(len(p)),
                "Strategy": strat, "Total log ret": float(p.sum()),
                "Ann Sharpe": float(np.sqrt(252) * p.mean() / (p.std() + 1e-12)),
                "Max DD": float(dd),
            })
    crisis_df = pd.DataFrame(rows)
    if not crisis_df.empty:
        _csv(crisis_df.set_index(["Crisis", "Strategy"]), "T07_crisis_performance")

    # ---- 6) Hierarchical C stress diagnostics ------------------------
    sC = hierC_stress.loc[common]
    gC = hierC_gamma.loc[common]
    stress_tab = pd.DataFrame({
        "macro_stress": [sC.mean(), sC.median(), sC.std(), sC.min(), sC.max()],
        "gamma_eff":    [gC.mean(), gC.median(), gC.std(), gC.min(), gC.max()],
    }, index=["mean", "median", "std", "min", "max"])
    _csv(stress_tab, "T10_hierC_stress_summary")

    # ---- 7) Per-strategy chart bundles -------------------------------
    print("[agg] plotting pure-market bundle...", flush=True)
    plot_strategy_bundle(
        pnl=pure["pnl"].loc[common], weights=pure["W"].loc[common],
        regime=pure["lbl"].loc[common], K_history=pure["K"].loc[common],
        G_history=pure["cnt"].loc[common], tpl_prob=pure["prob"].loc[common],
        asset_returns=asset_rets.loc[common], prefix="pure_",
    )

    for tag, hd in [("hierB", hierB), ("hierC", hierC)]:
        print(f"[agg] plotting {tag} bundle...", flush=True)
        plot_strategy_bundle(
            pnl=hd["pnl"].loc[common], weights=hd["weights"].loc[common],
            regime=hd["mkt_lbl"].loc[common], K_history=hd["K_market"].loc[common],
            G_history=hd["mkt_cnt"].loc[common], tpl_prob=hd["mkt_prob"].loc[common],
            asset_returns=asset_rets.loc[common], prefix=f"{tag}_",
        )
        plot_macro_diagnostics(
            macro_label=hd["macro_lbl"].loc[common],
            macro_prob=hd["macro_prob"].loc[common],
            market_label=hd["mkt_lbl"].loc[common],
            K_macro=hd["K_macro"].loc[common], G_macro=hd["G_macro"].loc[common],
            asset_returns=asset_rets.loc[common], prefix=f"{tag}_",
        )

    # Hier-C stress / gamma chart
    plot_macro_stress(sC, gC, name="hierC_macro_stress")

    # ---- 8) Cross-strategy comparison --------------------------------
    print("[agg] cross-strategy comparison charts...", flush=True)
    plot_strategy_comparison(pnls, name="11_strategy_comparison")
    plot_drawdown_comparison(pnls, name="12_drawdown_comparison")
    plot_rolling_sharpe(pnls, window=252, name="13_rolling_sharpe_1y")
    plot_annual_returns_bar(pnls, name="14_annual_returns")
    plot_underwater(pnls, name="15_underwater")

    # ---- 9) Unified PnL panel ----------------------------------------
    pnl_panel = pd.DataFrame(pnls).sort_index()
    pnl_panel.to_csv(os.path.join(R, "pnl_panel.csv"))
    pnl_panel.cumsum().to_csv(os.path.join(R, "cum_pnl_panel.csv"))

    bench_panel = pd.concat([
        eqw["weights"].add_prefix("w_eqw_").loc[common],
        eqw["pnl"].loc[common].rename("pnl_equal_weight"),
        s60["weights"].add_prefix("w_6040_").loc[common],
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
