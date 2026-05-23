"""
compare_baseline_vs_hierarchical.py
-----------------------------------
A/B comparison module for the hierarchical (v10) macro extension.

Reuses most of v9's chart machinery but replaces the macro-specific
diagnostics (beta heatmap, transition matrix snapshots) with v10's:
    * A_per_macro matrices (one per macro regime)
    * Macro regime timeline (over the full OOS)
    * Joint macro+market regime heatmap (which combinations fire?)
    * Macro regime distinctness (W2 distances between A_m matrices)

Outputs (saved under output/compare/):
    Charts
    ------
    01_cumulative_pnl.png              — Cumulative log PnL, baseline vs hierarchical
    02_drawdown_overlay.png            — Underwater curves
    03_yearly_sharpe.png               — Per-year Sharpe ratio bar chart
    04_market_regime_occupancy.png     — Market template occupancy
    05_stress_period_zoom.png          — 2008-09, 2020, 2022 close-ups
    06_macro_regime_means.png          — Per-macro-regime mean of 7 macro variables
    07_A_per_macro_heatmaps.png        — One K×K transition matrix per macro regime
    08_macro_timeline.png              — Macro regime over time
    09_weight_evolution.png            — Weights over time, side by side
    10_market_regime_timeline.png      — Market template regime over time
    11_joint_regime_heatmap.png        — Co-occurrence of macro × market regimes

    Tables
    ------
    performance_comparison.csv         — Headline metrics side by side
    yearly_metrics.csv                 — Year-by-year metrics
    market_regime_occupancy.csv        — Days per market regime
    macro_regime_occupancy.csv         — Days per macro regime
    stress_period_metrics.csv          — Per stress window
    A_per_macro.csv                    — Long-format transition matrices
    macro_regime_means.csv             — Mean of each macro variable per regime
    joint_regime_occupancy.csv         — Co-occurrence counts
    hierarchical_daily_backtest_output.csv — Full daily output
"""

from __future__ import annotations

import _paths  # noqa: F401

import warnings
from pathlib import Path
from typing import Optional

import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from config import ASSET_NAMES
from data_macro import MACRO_COLS

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
warnings.filterwarnings("ignore", category=FutureWarning)

# Visual style
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 180, "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.titleweight": "bold", "axes.titlesize": 12,
})

COL_BASE = "#2C5282"
COL_HIER = "#22543D"      # dark green (distinguish from v9's red)
COL_SPX = "#718096"

PALETTE_REGIME = ["#3182CE", "#DD6B20", "#38A169", "#E53E3E", "#805AD5",
                  "#D69E2E", "#319795", "#B83280"]


# --------------------------------------------------------------------- #
# 1. Metrics helpers
# --------------------------------------------------------------------- #
def _metrics(pnl: pd.Series) -> dict:
    if pnl.std() == 0 or len(pnl) < 2:
        return {"Sharpe": np.nan, "AnnRet": np.nan, "AnnVol": np.nan, "MaxDD": np.nan}
    cum = pnl.cumsum()
    mdd = (cum - cum.cummax()).min()
    return {
        "Sharpe": pnl.mean() / pnl.std() * np.sqrt(252),
        "AnnRet": pnl.mean() * 252,
        "AnnVol": pnl.std() * np.sqrt(252),
        "MaxDD": float(mdd),
    }


def compare_performance(baseline_pnl, hier_pnl, spx_pnl) -> pd.DataFrame:
    aligned = pd.concat({"baseline": baseline_pnl, "hierarchical": hier_pnl, "spx_bh": spx_pnl}, axis=1).dropna()
    out = pd.DataFrame({
        "Baseline":      _metrics(aligned["baseline"]),
        "Hierarchical":  _metrics(aligned["hierarchical"]),
        "SPX B&H":       _metrics(aligned["spx_bh"]),
    })
    out.loc["CumLogRet"] = [aligned["baseline"].sum(), aligned["hierarchical"].sum(), aligned["spx_bh"].sum()]
    return out


def yearly_comparison(baseline_pnl, hier_pnl) -> pd.DataFrame:
    aligned = pd.concat({"baseline": baseline_pnl, "hierarchical": hier_pnl}, axis=1).dropna()
    rows = []
    for y in sorted(aligned.index.year.unique()):
        sub = aligned[aligned.index.year == y]
        if len(sub) < 50:
            continue
        row = {"Year": y}
        for label in ("baseline", "hierarchical"):
            p = sub[label]
            sh = p.mean() / p.std() * np.sqrt(252) if p.std() > 0 else np.nan
            cum = p.cumsum()
            mdd = (cum - cum.cummax()).min()
            row[f"{label}_Sharpe"] = sh
            row[f"{label}_Ret"]    = p.sum()
            row[f"{label}_MaxDD"]  = float(mdd)
        rows.append(row)
    return pd.DataFrame(rows).set_index("Year")


# --------------------------------------------------------------------- #
# 2. Charts shared with v9 (style adjusted for "hierarchical")
# --------------------------------------------------------------------- #
def chart_cumulative(out_dir, baseline, hier, spx_pnl):
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(baseline["pnl"].cumsum(), label="Baseline (constant A)", color=COL_BASE, linewidth=1.6)
    ax.plot(hier["pnl"].cumsum(), label="Baseline + Macro Tilt (v10.2)", color=COL_HIER, linewidth=1.6)
    ax.plot(spx_pnl.cumsum(), label="SPX Buy & Hold", color=COL_SPX, linewidth=1.2, linestyle="--", alpha=0.8)
    ax.set_ylabel("Cumulative Log Return")
    ax.set_title("Out-of-Sample Cumulative PnL")
    ax.legend(loc="upper left", frameon=False)
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    fig.tight_layout()
    fig.savefig(out_dir / "01_cumulative_pnl.png")
    plt.close(fig)


def chart_drawdown(out_dir, baseline, hier, spx_pnl):
    fig, ax = plt.subplots(figsize=(11, 4.5))
    for name, pnl, col in [("Baseline", baseline["pnl"], COL_BASE),
                            ("Hierarchical", hier["pnl"], COL_HIER),
                            ("SPX B&H", spx_pnl, COL_SPX)]:
        cum = pnl.cumsum()
        dd = cum - cum.cummax()
        ax.fill_between(dd.index, dd, 0, color=col, alpha=0.3, label=name)
    ax.set_ylabel("Drawdown (log)")
    ax.set_title("Underwater Curves")
    ax.legend(loc="lower left", frameon=False)
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    fig.tight_layout()
    fig.savefig(out_dir / "02_drawdown_overlay.png")
    plt.close(fig)


def chart_yearly_sharpe(out_dir, yearly):
    fig, ax = plt.subplots(figsize=(13, 4.5))
    years = yearly.index.values
    x = np.arange(len(years))
    w = 0.35
    ax.bar(x - w/2, yearly["baseline_Sharpe"], w, label="Baseline", color=COL_BASE)
    ax.bar(x + w/2, yearly["hierarchical_Sharpe"], w, label="Hierarchical", color=COL_HIER)
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_xticks(x); ax.set_xticklabels(years, rotation=45)
    ax.set_ylabel("Annualized Sharpe ratio")
    ax.set_title("Year-by-Year Sharpe Ratio")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "03_yearly_sharpe.png")
    plt.close(fig)


def chart_market_regime_occupancy(out_dir, baseline, hier):
    G_base = baseline["tpl_label"].dropna().astype(int)
    G_hier = hier["tpl_label"].dropna().astype(int)
    all_regs = sorted(set(G_base.unique()) | set(G_hier.unique()))
    occ_base = G_base.value_counts(normalize=True).reindex(all_regs, fill_value=0)
    occ_hier = G_hier.value_counts(normalize=True).reindex(all_regs, fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(all_regs)); w = 0.35
    ax.bar(x - w/2, occ_base.values * 100, w, label="Baseline", color=COL_BASE)
    ax.bar(x + w/2, occ_hier.values * 100, w, label="Hierarchical", color=COL_HIER)
    ax.set_xticks(x); ax.set_xticklabels([f"R{r}" for r in all_regs])
    ax.set_ylabel("Fraction of OOS Days (%)")
    ax.set_title("Market Template Regime Occupancy")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "04_market_regime_occupancy.png")
    plt.close(fig)


def chart_stress_zoom(out_dir, baseline, hier, spx_pnl):
    windows = [
        ("2008-09 GFC",      "2008-08-01", "2009-06-30"),
        ("2020 COVID",       "2020-02-01", "2020-06-30"),
        ("2022 Inflation",   "2022-01-01", "2022-12-31"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, (label, lo, hi) in zip(axes, windows):
        b = baseline["pnl"].loc[lo:hi].cumsum()
        h = hier["pnl"].loc[lo:hi].cumsum()
        s = spx_pnl.loc[lo:hi].cumsum()
        ax.plot(b, color=COL_BASE,  linewidth=1.7, label="Baseline")
        ax.plot(h, color=COL_HIER, linewidth=1.7, label="Hierarchical")
        ax.plot(s, color=COL_SPX,   linewidth=1.2, linestyle="--", alpha=0.8, label="SPX B&H")
        ax.set_title(label); ax.axhline(0, color="k", linewidth=0.4)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b-%y"))
        for tick in ax.get_xticklabels(): tick.set_rotation(45)
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].set_ylabel("Cumulative Log Return (window)")
    fig.suptitle("Stress-Period Zoom", fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_dir / "05_stress_period_zoom.png")
    plt.close(fig)


# --------------------------------------------------------------------- #
# 3. v10-specific charts
# --------------------------------------------------------------------- #
def chart_macro_regime_means(out_dir, macro_hmm):
    """Bar chart showing how each macro variable behaves under each
    macro regime. Most direct interpretability chart for thesis."""
    if macro_hmm is None or macro_hmm.means_ is None:
        return
    means = macro_hmm.means_                # (M, d_macro)
    M, d = means.shape
    fig, ax = plt.subplots(figsize=(13, 5))
    x = np.arange(d)
    w = 0.8 / M
    for m in range(M):
        ax.bar(x + (m - (M-1)/2) * w, means[m], w,
               label=f"Macro R{m}",
               color=PALETTE_REGIME[m % len(PALETTE_REGIME)])
    ax.set_xticks(x); ax.set_xticklabels(MACRO_COLS, rotation=20)
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_ylabel("Mean z-score of macro variable")
    ax.set_title("Macro Regime Profiles  (mean of each macro variable per regime)")
    ax.legend(ncol=M, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout()
    fig.savefig(out_dir / "06_macro_regime_means.png", bbox_inches="tight")
    plt.close(fig)


def chart_macro_tilt_bars(out_dir, macro_tilts, M_macro):
    """Bar chart of per-macro-regime mean return tilts (annualized %).

    This is v10.2's equivalent of the v10 A_per_macro heatmap: it shows
    the per-regime contribution that the macro layer adds to the MVO.
    """
    if macro_tilts is None:
        return
    tilt = macro_tilts["macro_mu_tilt"] * 252 * 100  # annualized %
    counts = macro_tilts["macro_counts"]
    uncond = macro_tilts["unconditional_mu"] * 252 * 100

    fig, ax = plt.subplots(figsize=(13, 5.5))
    n_assets = tilt.shape[1]
    x = np.arange(n_assets)
    w = 0.8 / M_macro

    for m in range(M_macro):
        ax.bar(x + (m - (M_macro - 1) / 2) * w,
               tilt[m], w,
               label=f"Macro R{m}  (n={counts[m]})",
               color=PALETTE_REGIME[m % len(PALETTE_REGIME)])

    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_xticks(x); ax.set_xticklabels(ASSET_NAMES)
    ax.set_ylabel("Mean return tilt vs unconditional (annualized %)")
    ax.set_title("v10.2 Macro Tilts: Calibration-Sample Return Deviation Per Macro Regime")
    ax.legend(ncol=M_macro, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.10))

    # Annotate unconditional means below
    txt = "Unconditional calibration mean (ann %): " + \
          ", ".join(f"{a}={uncond[i]:+.1f}" for i, a in enumerate(ASSET_NAMES))
    ax.text(0.5, -0.22, txt, transform=ax.transAxes,
            ha="center", fontsize=9, style="italic")

    fig.tight_layout()
    fig.savefig(out_dir / "07_macro_tilt_bars.png", bbox_inches="tight")
    plt.close(fig)


def chart_macro_timeline(out_dir, hier, macro_seq_oos, n_regimes):
    """Macro regime over time with SPX cum log-return for context."""
    fig, ax = plt.subplots(figsize=(14, 4))
    spx_cum = hier["pnl"].cumsum() * 0 + 1  # placeholder (will overlay raw SPX below)

    lbl = macro_seq_oos.dropna().astype(int)
    change_idx = (lbl != lbl.shift()).cumsum()
    for _, grp in lbl.groupby(change_idx):
        m = int(grp.iloc[0])
        if m < 0:
            continue
        ax.axvspan(grp.index[0], grp.index[-1],
                   color=PALETTE_REGIME[m % len(PALETTE_REGIME)], alpha=0.65, lw=0)

    ax2 = ax.twinx()
    hier_cum = hier["pnl"].cumsum()
    ax2.plot(hier_cum.index, hier_cum.values, color="#222", linewidth=0.9, alpha=0.8,
             label="Hierarchical cum. log-ret")
    ax2.set_ylabel("Hierarchical cum. log-ret")
    ax.set_ylim(0, 1); ax.set_yticks([])
    ax.set_title("Macro Regime Over Time  (color = macro regime per day)")

    handles = [mpatches.Patch(color=PALETTE_REGIME[m], label=f"Macro R{m}", alpha=0.65)
               for m in range(n_regimes)]
    ax.legend(handles=handles, loc="upper left", ncol=n_regimes,
              bbox_to_anchor=(0, 1.15), frameon=False, fontsize=9)
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    fig.tight_layout()
    fig.savefig(out_dir / "08_macro_timeline.png", bbox_inches="tight")
    plt.close(fig)


def chart_weight_evolution(out_dir, baseline, hier):
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    cols = ["#3182CE", "#DD6B20", "#D69E2E", "#E53E3E", "#805AD5"]
    for ax, (label, dct) in zip(axes, [("Baseline (constant A)", baseline),
                                       ("Baseline + Macro Tilt (v10.2)", hier)]):
        wts = dct["weights"]
        ax.stackplot(wts.index, wts.T.values, labels=wts.columns,
                     colors=cols, alpha=0.85)
        ax.set_title(label)
        ax.set_ylabel("Portfolio weight")
        ax.set_ylim(0, 1)
    axes[0].legend(loc="upper right", ncol=5, frameon=False, fontsize=9)
    axes[1].xaxis.set_major_locator(mdates.YearLocator(2))
    fig.tight_layout()
    fig.savefig(out_dir / "09_weight_evolution.png")
    plt.close(fig)


def chart_market_regime_timeline(out_dir, hier, n_regimes_market):
    fig, ax = plt.subplots(figsize=(14, 3))
    lbl = hier["tpl_label"].dropna().astype(int)
    change_idx = (lbl != lbl.shift()).cumsum()
    for _, grp in lbl.groupby(change_idx):
        r = int(grp.iloc[0])
        ax.axvspan(grp.index[0], grp.index[-1],
                   color=PALETTE_REGIME[r % len(PALETTE_REGIME)], alpha=0.65, lw=0)
    ax.set_ylim(0, 1); ax.set_yticks([])
    ax.set_title("Market Regime Over Time  (color = template regime per day)")
    handles = [mpatches.Patch(color=PALETTE_REGIME[r], label=f"Market R{r}", alpha=0.65)
               for r in range(n_regimes_market)]
    ax.legend(handles=handles, loc="upper left", ncol=n_regimes_market,
              bbox_to_anchor=(0, 1.4), frameon=False, fontsize=9)
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    fig.tight_layout()
    fig.savefig(out_dir / "10_market_regime_timeline.png", bbox_inches="tight")
    plt.close(fig)


def chart_joint_regime_heatmap(out_dir, hier, macro_seq_oos):
    """Co-occurrence of macro × market regimes (which combinations fire?)"""
    market_lbl = hier["tpl_label"].dropna().astype(int)
    aligned = pd.concat({"macro": macro_seq_oos, "market": market_lbl}, axis=1).dropna()
    aligned = aligned[aligned["macro"] >= 0]  # drop -1 sentinel
    co = pd.crosstab(aligned["macro"].astype(int), aligned["market"].astype(int))
    co_pct = co / co.values.sum() * 100.0

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(co_pct.values, cmap="Blues", aspect="auto")
    ax.set_xticks(range(co_pct.shape[1]))
    ax.set_xticklabels([f"Market R{c}" for c in co_pct.columns])
    ax.set_yticks(range(co_pct.shape[0]))
    ax.set_yticklabels([f"Macro R{r}" for r in co_pct.index])
    for i in range(co_pct.shape[0]):
        for j in range(co_pct.shape[1]):
            v = co_pct.values[i, j]
            color = "white" if v > co_pct.values.max() * 0.5 else "black"
            ax.text(j, i, f"{v:.1f}%", ha="center", va="center", fontsize=9, color=color)
    plt.colorbar(im, ax=ax, label="% of OOS days")
    ax.set_title("Joint Regime Co-occurrence  (macro × market)")
    fig.tight_layout()
    fig.savefig(out_dir / "11_joint_regime_heatmap.png")
    plt.close(fig)


# --------------------------------------------------------------------- #
# 4. CSV exports
# --------------------------------------------------------------------- #
def export_market_occupancy_csv(out_dir, baseline, hier):
    G_base = baseline["tpl_label"].dropna().astype(int)
    G_hier = hier["tpl_label"].dropna().astype(int)
    all_regs = sorted(set(G_base.unique()) | set(G_hier.unique()))
    df = pd.DataFrame({
        "regime":             all_regs,
        "baseline_days":      [(G_base == r).sum() for r in all_regs],
        "hierarchical_days":  [(G_hier == r).sum() for r in all_regs],
        "baseline_pct":       [(G_base == r).mean() for r in all_regs],
        "hierarchical_pct":   [(G_hier == r).mean() for r in all_regs],
    })
    df.to_csv(out_dir / "market_regime_occupancy.csv", index=False)
    return df


def export_macro_occupancy_csv(out_dir, macro_seq_oos):
    s = macro_seq_oos.dropna()
    s = s[s >= 0]
    df = pd.DataFrame({
        "macro_regime": sorted(s.unique()),
        "days":         [(s == m).sum() for m in sorted(s.unique())],
        "pct":          [(s == m).mean() for m in sorted(s.unique())],
    })
    df.to_csv(out_dir / "macro_regime_occupancy.csv", index=False)
    return df


def export_macro_tilts_csv(out_dir, macro_tilts, M_macro):
    """Export the per-macro-regime calibration-time tilts as CSV."""
    if macro_tilts is None:
        return None
    rows = []
    tilt = macro_tilts["macro_mu_tilt"] * 252 * 100
    uncond = macro_tilts["unconditional_mu"] * 252 * 100
    for m in range(M_macro):
        row = {
            "macro_regime": m,
            "n_cal_days": int(macro_tilts["macro_counts"][m]),
            "vol_scale": float(macro_tilts["macro_vol_scale"][m]),
        }
        for i, a in enumerate(ASSET_NAMES):
            row[f"{a}_tilt_ann_pct"]    = float(tilt[m, i])
            row[f"{a}_uncond_ann_pct"]  = float(uncond[i])
            row[f"{a}_implied_ann_pct"] = float(uncond[i] + tilt[m, i])
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "macro_tilts.csv", index=False)
    return df


def export_macro_means_csv(out_dir, macro_hmm):
    if macro_hmm is None or macro_hmm.means_ is None:
        return None
    df = pd.DataFrame(macro_hmm.means_, columns=MACRO_COLS)
    df.index.name = "macro_regime"
    df.to_csv(out_dir / "macro_regime_means.csv")
    return df


def export_joint_occupancy_csv(out_dir, hier, macro_seq_oos):
    market_lbl = hier["tpl_label"].dropna().astype(int)
    aligned = pd.concat({"macro": macro_seq_oos, "market": market_lbl}, axis=1).dropna()
    aligned = aligned[aligned["macro"] >= 0]
    co = pd.crosstab(aligned["macro"].astype(int), aligned["market"].astype(int))
    co.index.name = "macro_regime"; co.columns.name = "market_regime"
    co.to_csv(out_dir / "joint_regime_occupancy.csv")
    return co


def export_stress_period_csv(out_dir, baseline, hier, spx_pnl):
    windows = [
        ("2008-09 GFC",          "2008-08-01", "2009-06-30"),
        ("2011 EU debt crisis",  "2011-07-01", "2011-12-31"),
        ("2015-16 China/oil",    "2015-08-01", "2016-02-29"),
        ("2018 Q4 sell-off",     "2018-09-15", "2018-12-31"),
        ("2020 COVID",           "2020-02-01", "2020-06-30"),
        ("2022 Inflation",       "2022-01-01", "2022-12-31"),
        ("2025 Liberation Day",  "2025-04-01", "2025-06-30"),
    ]
    rows = []
    for name, lo, hi in windows:
        b = baseline["pnl"].loc[lo:hi]
        h = hier["pnl"].loc[lo:hi]
        s = spx_pnl.loc[lo:hi]
        if b.empty or h.empty:
            continue
        def _met(p):
            if p.std() == 0 or len(p) < 2:
                return np.nan, np.nan, np.nan
            sh = p.mean() / p.std() * np.sqrt(252)
            cum = p.cumsum()
            mdd = float((cum - cum.cummax()).min())
            return sh, p.sum(), mdd
        b_sh, b_ret, b_mdd = _met(b)
        h_sh, h_ret, h_mdd = _met(h)
        s_sh, s_ret, s_mdd = _met(s)
        rows.append({
            "period": name, "start": lo, "end": hi, "n_days": len(b),
            "baseline_Sharpe": b_sh,     "baseline_Ret": b_ret,     "baseline_MaxDD": b_mdd,
            "hier_Sharpe":     h_sh,     "hier_Ret":     h_ret,     "hier_MaxDD":     h_mdd,
            "spx_Sharpe":      s_sh,     "spx_Ret":      s_ret,     "spx_MaxDD":      s_mdd,
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "stress_period_metrics.csv", index=False)
    return df


def export_daily_output_csv(out_dir, hier_result, macro_seq_oos):
    df = pd.DataFrame(index=hier_result.pnl.index)
    for c in hier_result.weights.columns:
        df[f"w_{c}"] = hier_result.weights[c]
    df["pnl"]      = hier_result.pnl
    df["cum_pnl"]  = hier_result.cum_pnl
    df["K"]        = hier_result.K_history
    df["regime"]   = hier_result.tpl_label
    df["max_p"]    = hier_result.tpl_max_prob
    df["G"]        = hier_result.tpl_count
    df["turnover"] = hier_result.turnover
    df["macro"]    = macro_seq_oos.reindex(df.index)
    df.index.name = "date"
    df.to_csv(out_dir / "hierarchical_daily_backtest_output.csv")
    return df


# --------------------------------------------------------------------- #
# 5. Master comparison
# --------------------------------------------------------------------- #
def run_full_comparison(
    baseline_result,
    hier_result,
    hier_diagnostics: dict,
    spx_pnl: pd.Series,
    output_dir: Path | str = None,
) -> dict:
    out_dir = Path(output_dir) if output_dir is not None else (OUTPUT_DIR / "compare")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[compare] writing outputs to {out_dir}", flush=True)

    baseline = {
        "pnl":          baseline_result.pnl,
        "weights":      baseline_result.weights,
        "tpl_label":    baseline_result.tpl_label,
        "tpl_max_prob": baseline_result.tpl_max_prob,
    }
    hier = {
        "pnl":          hier_result.pnl,
        "weights":      hier_result.weights,
        "tpl_label":    hier_result.tpl_label,
        "tpl_max_prob": hier_result.tpl_max_prob,
    }
    macro_seq_oos = hier_diagnostics.get("macro_seq_oos")
    macro_hmm     = hier_diagnostics.get("macro_hmm")
    macro_tilts   = hier_diagnostics.get("macro_tilts")
    M_macro       = hier_diagnostics.get("M_macro", macro_hmm.M_ if macro_hmm else 1)

    n_market_regimes = int(pd.concat([baseline["tpl_label"], hier["tpl_label"]]).max() + 1)
    n_macro_regimes  = M_macro

    # Charts
    chart_cumulative(out_dir, baseline, hier, spx_pnl)
    chart_drawdown(out_dir, baseline, hier, spx_pnl)
    yearly = yearly_comparison(baseline_result.pnl, hier_result.pnl)
    chart_yearly_sharpe(out_dir, yearly)
    chart_market_regime_occupancy(out_dir, baseline, hier)
    chart_stress_zoom(out_dir, baseline, hier, spx_pnl)
    chart_macro_regime_means(out_dir, macro_hmm)
    chart_macro_tilt_bars(out_dir, macro_tilts, n_macro_regimes)
    if macro_seq_oos is not None:
        chart_macro_timeline(out_dir, hier, macro_seq_oos, n_macro_regimes)
    chart_weight_evolution(out_dir, baseline, hier)
    chart_market_regime_timeline(out_dir, hier, n_market_regimes)
    if macro_seq_oos is not None:
        chart_joint_regime_heatmap(out_dir, hier, macro_seq_oos)

    # CSVs
    perf = compare_performance(baseline_result.pnl, hier_result.pnl, spx_pnl)
    perf.to_csv(out_dir / "performance_comparison.csv")
    yearly.to_csv(out_dir / "yearly_metrics.csv")
    export_market_occupancy_csv(out_dir, baseline, hier)
    if macro_seq_oos is not None:
        export_macro_occupancy_csv(out_dir, macro_seq_oos)
        export_joint_occupancy_csv(out_dir, hier, macro_seq_oos)
    export_macro_tilts_csv(out_dir, macro_tilts, n_macro_regimes)
    export_macro_means_csv(out_dir, macro_hmm)
    export_stress_period_csv(out_dir, baseline, hier, spx_pnl)
    export_daily_output_csv(out_dir, hier_result, macro_seq_oos)

    print(f"\n{'='*60}")
    print(f"  PERFORMANCE COMPARISON")
    print(f"{'='*60}")
    print(perf.round(4).to_string())
    print()

    return {
        "performance":  perf,
        "yearly":       yearly,
        "output_dir":   out_dir,
    }
