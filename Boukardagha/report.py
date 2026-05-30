"""
report.py
=========
Single-page HTML report knitting together charts, tables, and findings
for the three regime-aware strategies and two benchmarks.
Run AFTER aggregate.py.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
import config

REPORT_PATH = os.path.join(config.OUTPUT_DIR, "report.html")


def _img(path, caption="", width=900):
    if not os.path.exists(path):
        return f"<p><em>(missing {os.path.basename(path)})</em></p>"
    rel = os.path.relpath(path, config.OUTPUT_DIR)
    return (f'<figure style="margin:1em 0"><img src="{rel}" '
            f'style="width:{width}px;max-width:100%"/>'
            f'<figcaption style="font-size:90%;color:#555">{caption}'
            f'</figcaption></figure>')


def _tbl(csv_path, caption="", index_col=0, max_rows=60):
    if not os.path.exists(csv_path):
        return f"<p><em>(missing {os.path.basename(csv_path)})</em></p>"
    df = pd.read_csv(csv_path, index_col=index_col)
    if len(df) > max_rows:
        df = df.head(max_rows)
    return f"<h3>{caption}</h3>" + df.to_html(
        border=0, classes="data", float_format=lambda x: f"{x:.4f}")


def main():
    print("[report] building HTML report...", flush=True)
    T, C = config.TABLE_DIR, config.CHART_DIR

    h = ["<!DOCTYPE html><html><head><meta charset='utf-8'>",
         "<title>Regime-Aware Investing — Replication & Two Hierarchical Extensions</title>",
         """<style>
         body{font-family:-apple-system,BlinkMacSystemFont,sans-serif;
              max-width:1100px;margin:2em auto;padding:0 1em;color:#222;line-height:1.45;}
         h1{border-bottom:2px solid #999;padding-bottom:6px;}
         h2{border-bottom:1px solid #ccc;padding-bottom:4px;margin-top:2em;}
         table.data{border-collapse:collapse;margin:1em 0;font-size:90%;}
         table.data th,table.data td{padding:4px 10px;border-bottom:1px solid #ddd;}
         table.data th{background:#f0f0f0;text-align:left;}
         figure{text-align:center;}
         .note{background:#eef6ff;padding:.5em 1em;border-left:4px solid #3b82f6;margin:.5em 0;}
         .dev{background:#fff8e0;padding:.5em 1em;border-left:4px solid #d4a017;margin:.5em 0;}
         </style></head><body>"""]

    h.append("<h1>Explainable Regime-Aware Investing</h1>")
    h.append("<p><b>Boukardagha (2026) replication</b> + "
             "<b>two hierarchical (FST98) extensions</b></p>")
    h.append(f"<p>OOS window: <b>{config.OOS_START}</b> &rarr; <b>{config.DATA_END}</b>.</p>")
    h.append("""<div class="note">
        <b>Strategies.</b>
        <b>PureMarket_WHMM</b>: Boukardagha market-regime WHMM + MVO.
        <b>Hierarchical_B</b>: macro × market joint-mixture moments, no
        expected-return tilt. <b>Hierarchical_C</b>: macro layer used as
        a <i>risk modulator</i> — macro templates tracked in
        asset-outcome space, tempered macro posterior, stress-scaled
        risk aversion &gamma;<sub>t</sub> and covariance; the market
        layer retains full control of the expected-return direction.
        </div>""")

    # 1. Headline
    h.append("<h2>1. Headline performance</h2>")
    h.append(_tbl(os.path.join(T, "T01_performance_summary.csv"),
                  "Risk-adjusted performance"))
    h.append(_img(os.path.join(C, "11_strategy_comparison.png"),
                  "Cumulative log return — all strategies"))
    h.append(_img(os.path.join(C, "12_drawdown_comparison.png"),
                  "Drawdown comparison (on wealth)"))
    h.append(_img(os.path.join(C, "15_underwater.png"), "Underwater drawdown"))
    h.append(_img(os.path.join(C, "13_rolling_sharpe_1y.png"),
                  "Rolling 1y annualised Sharpe"))
    h.append(_img(os.path.join(C, "14_annual_returns.png"),
                  "Annual log returns per strategy"))

    # 2. Frictions
    h.append("<h2>2. Trading frictions and concentration</h2>")
    h.append(_tbl(os.path.join(T, "T02_turnover_summary.csv"), "Turnover summary"))
    h.append(_tbl(os.path.join(T, "T04_concentration_summary.csv"), "Effective # positions"))
    h.append(_tbl(os.path.join(T, "T03a_allocation_pure.csv"), "Pure-market allocation"))
    h.append(_tbl(os.path.join(T, "T03b_allocation_hierB.csv"), "Hierarchical B allocation"))
    h.append(_tbl(os.path.join(T, "T03c_allocation_hierC.csv"), "Hierarchical C allocation"))

    # 3. Pure-market diagnostics
    h.append("<h2>3. Pure-Market diagnostics</h2>")
    for fn, cap in [
        ("pure_01_cum_pnl_scatter.png", "Cum PnL coloured by template"),
        ("pure_03_weights_stacked.png", "Portfolio weights"),
        ("pure_05_asset_sharpe_by_regime.png", "Asset Sharpe by template"),
        ("pure_06_stacked_cum_pnl_by_regime.png", "Stacked cumulative PnL by template"),
        ("pure_07_K_history.png", "Selected K over time"),
        ("pure_09_template_label.png", "Dominant template over time"),
    ]:
        h.append(_img(os.path.join(C, fn), cap))
    h.append(_tbl(os.path.join(T, "T05a_pure_by_regime.csv"), "Portfolio performance by template"))
    h.append(_tbl(os.path.join(T, "T06a_pure_asset_by_regime.csv"), "Asset Sharpe/Mean/Vol by template"))
    h.append(_tbl(os.path.join(T, "T09a_pure_persistence.csv"), "Regime persistence"))

    # 4. Hierarchical B
    h.append("<h2>4. Hierarchical B — macro × market joint mixture</h2>")
    for fn, cap in [
        ("hierB_macro_label.png", "Dominant macro template over time"),
        ("hierB_macro_prob.png", "Max macro posterior over time"),
        ("hierB_macro_vs_market.png", "Joint frequency of (macro, market) templates"),
        ("hierB_macro_macro_05_asset_sharpe_by_regime.png", "Asset Sharpe by macro template"),
        ("hierB_03_weights_stacked.png", "Portfolio weights"),
    ]:
        h.append(_img(os.path.join(C, fn), cap))
    h.append(_tbl(os.path.join(T, "T05b_hierB_by_macro_regime.csv"), "Performance by macro template"))
    h.append(_tbl(os.path.join(T, "T06b_hierB_asset_by_macro_regime.csv"), "Asset stats by macro template"))
    h.append(_tbl(os.path.join(T, "T09b_hierB_macro_persistence.csv"), "Macro regime persistence (feature-space templates)"))

    # 5. Hierarchical C
    h.append("<h2>5. Hierarchical C — macro as risk modulator (Fix C)</h2>")
    h.append("""<div class="note">
        Macro templates are tracked in <b>asset-outcome space</b>, the
        macro posterior is <b>tempered</b> and blended with a uniform
        prior, and the resulting stress score scales the effective risk
        aversion &gamma;<sub>t</sub> = &gamma;·(1 + &kappa;·stress) and
        inflates &Sigma;<sub>t</sub>. The market layer owns &mu;.
        </div>""")
    for fn, cap in [
        ("hierC_macro_stress.png", "Macro stress and effective risk aversion"),
        ("hierC_macro_label.png", "Dominant macro template over time"),
        ("hierC_macro_prob.png", "Max macro posterior over time"),
        ("hierC_macro_vs_market.png", "Joint frequency of (macro, market) templates"),
        ("hierC_macro_macro_05_asset_sharpe_by_regime.png", "Asset Sharpe by macro template"),
        ("hierC_03_weights_stacked.png", "Portfolio weights"),
    ]:
        h.append(_img(os.path.join(C, fn), cap))
    h.append(_tbl(os.path.join(T, "T05c_hierC_by_macro_regime.csv"), "Performance by macro template"))
    h.append(_tbl(os.path.join(T, "T06c_hierC_asset_by_macro_regime.csv"), "Asset stats by macro template"))
    h.append(_tbl(os.path.join(T, "T09c_hierC_macro_persistence.csv"), "Macro regime persistence (outcome-space templates)"))
    h.append(_tbl(os.path.join(T, "T10_hierC_stress_summary.csv"), "Macro stress / effective gamma summary"))

    # 6. Crisis
    cpath = os.path.join(T, "T07_crisis_performance.csv")
    if os.path.exists(cpath):
        h.append("<h2>6. Crisis-window stress test</h2>")
        h.append(_tbl(cpath, "Per-strategy performance over named crises", index_col=[0, 1]))

    # 7. Methodology
    h.append("<h2>7. Methodology notes</h2>")
    h.append("<p>Full mathematical specification in <code>methodology.md</code>; "
             "analysis of the macro-overlay behaviour in "
             "<code>analysis_results.md</code>.</p>")
    h.append('<div class="dev">Only deviation from Boukardagha (2026): '
             '<code>HMM_N_ITER = 100</code> (paper: 300). EM converges in '
             '&le; 80 iters at tol=1e-3 on this data; set to 300 to recover '
             'the paper exactly.</div>')
    h.append('<div class="dev">Daily refits are warm-started from the '
             'previous day&apos;s fitted parameters (EM monotonicity); the '
             'converged maximum is identical to a cold start.</div>')

    h.append("</body></html>")

    with open(REPORT_PATH, "w") as f:
        f.write("\n".join(h))
    print(f"[report] wrote {REPORT_PATH}", flush=True)


if __name__ == "__main__":
    main()
