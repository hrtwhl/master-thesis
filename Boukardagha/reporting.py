"""
reporting.py
------------
Everything the paper's “Empirical Results”, “Economic Interpretation”, and
“Benchmarking” sections do, packaged as deterministic export routines:

    * performance metrics (Sharpe, Sortino, max drawdown …)
    * benchmark strategies  — SPX Buy & Hold, Equal-Weight, **60/40** stocks/bonds
    * per-regime portfolio & asset analytics
    * PNG figures   (FIG_DIR)
    * CSV tables    (TBL_DIR)

The KNN benchmark from the paper is deliberately **not** reproduced here, per
the project brief.

Visual style
~~~~~~~~~~~~
A consistent palette and grid style are applied via `_apply_style`. Every
figure is saved at 150 dpi with `bbox_inches="tight"` and never `plt.show()`d
so this runs cleanly in CI / headless environments.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from backtest import BacktestResult
from config import (ASSET_NAMES, FIG_DIR, TBL_DIR, TRADING_DAYS, ensure_dirs)


# --------------------------------------------------------------------- #
# 0.  Style
# --------------------------------------------------------------------- #
_PALETTE = {
    "parametric": "#283d3b",   # deep navy
    "spx":        "#edddd4",   # green
    "ew":         "#197278",   # amber
    "sixty40":    "#c44536",   # purple
    "neutral":    "#374151",
    "accent":     "#dc2626",
    "gridcolor":  "#d4d4d8",
}
_ASSET_COLORS = {
    "SPX":  "#283D3B",
    "BOND": "#197278",
    "GOLD": "#EDDDD4",
    "OIL":  "#D99185",
    "USD":  "#C44536",
}
_REGIME_COLORS = ["#1f4e79", "#dc2626", "#16a34a", "#d97706", "#7c3aed",
                  "#0891b2", "#be185d", "#65a30d"]


def _apply_style() -> None:
    """Light, paper-ready matplotlib defaults."""
    plt.rcParams.update({
        "figure.dpi":        110,
        "savefig.dpi":       150,
        "savefig.bbox":      "tight",
        "font.size":         11,
        "axes.titlesize":    12,
        "axes.titleweight":  "bold",
        "axes.labelsize":    11,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.grid":         False,
        "grid.color":        _PALETTE["gridcolor"],
        "grid.alpha":        0.6,
        "grid.linewidth":    0.6,
        "legend.frameon":    False,
    })


def _save(fig: plt.Figure, name: str) -> Path:
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


# --------------------------------------------------------------------- #
# 1.  Core metrics
# --------------------------------------------------------------------- #
def ann_sharpe(r: pd.Series) -> float:
    r = r.dropna()
    return float(np.sqrt(TRADING_DAYS) * r.mean() / (r.std(ddof=1) + 1e-12))


def ann_sortino(r: pd.Series) -> float:
    r = r.dropna()
    downside = r[r < 0]
    if downside.empty:
        return float("inf")
    return float(np.sqrt(TRADING_DAYS) * r.mean() / (downside.std(ddof=1) + 1e-12))


def ann_mean(r: pd.Series) -> float:
    return float(r.mean() * TRADING_DAYS)


def ann_vol(r: pd.Series) -> float:
    return float(r.std(ddof=1) * np.sqrt(TRADING_DAYS))


def hit_rate(r: pd.Series) -> float:
    return float((r.dropna() > 0).mean())


def max_drawdown_log(cum_log: pd.Series) -> float:
    """Max drawdown computed on cumulative LOG returns (as in the paper)."""
    return float((cum_log - cum_log.cummax()).min())


def max_drawdown_simple(r: pd.Series) -> float:
    """Max drawdown of simple-return wealth path (used inside regimes)."""
    wealth = (1.0 + r.fillna(0.0)).cumprod()
    return float((wealth / wealth.cummax() - 1.0).min())


def drawdown_series(cum_log: pd.Series) -> pd.Series:
    return cum_log - cum_log.cummax()


# --------------------------------------------------------------------- #
# 2.  Benchmarks
# --------------------------------------------------------------------- #
def build_benchmarks(
    returns: pd.DataFrame,
    oos_index: pd.DatetimeIndex,
) -> dict[str, pd.Series]:
    """Return daily log-return series for each benchmark, aligned to `oos_index`.

    Benchmarks
    ----------
    * SPX Buy & Hold           : 100% SPX
    * Equal-Weight             : 20% each of the 5 assets, daily rebalanced
    * 60/40                    : 60% SPX, 40% BOND, daily rebalanced

    (Daily rebalancing is the standard convention for these passive baselines
    in the academic literature and matches the parametric strategy's daily
    decision frequency.)
    """
    R = returns.loc[oos_index, ASSET_NAMES]

    w_eq = np.ones(len(ASSET_NAMES)) / len(ASSET_NAMES)
    eq   = (R.to_numpy() @ w_eq).astype(float)

    w_6040 = np.zeros(len(ASSET_NAMES))
    w_6040[ASSET_NAMES.index("SPX")]  = 0.60
    w_6040[ASSET_NAMES.index("BOND")] = 0.40
    sixty40 = (R.to_numpy() @ w_6040).astype(float)

    spx = R["SPX"].to_numpy().astype(float)

    return {
        "SPX Buy & Hold":        pd.Series(spx,     index=oos_index, name="SPX B&H"),
        "Equal-Weight (20%)":    pd.Series(eq,      index=oos_index, name="EW"),
        "60/40 SPX/BOND":        pd.Series(sixty40, index=oos_index, name="60/40"),
    }


# --------------------------------------------------------------------- #
# 3.  Strategy-level tables (CSV)
# --------------------------------------------------------------------- #
def write_performance_table(
    result: BacktestResult,
    benchmarks: dict[str, pd.Series],
    path: Path = TBL_DIR / "performance_summary.csv",
) -> pd.DataFrame:
    """Sharpe, Sortino, Max DD vs. benchmarks (Table 7 of the paper)."""
    rows = {
        "Wasserstein HMM + MVO": result.pnl,
        **benchmarks,
    }
    df = pd.DataFrame({
        "Sharpe":       [ann_sharpe(r)                       for r in rows.values()],
        "Sortino":      [ann_sortino(r)                      for r in rows.values()],
        "Ann Mean":     [ann_mean(r)                         for r in rows.values()],
        "Ann Vol":      [ann_vol(r)                          for r in rows.values()],
        "Max Drawdown": [max_drawdown_log(r.cumsum())        for r in rows.values()],
        "Hit Rate":     [hit_rate(r)                         for r in rows.values()],
        "Total Days":   [int(r.dropna().shape[0])            for r in rows.values()],
    }, index=list(rows.keys()))
    df.to_csv(path, float_format="%.6f")
    return df


def write_turnover_table(
    result: BacktestResult,
    path: Path = TBL_DIR / "turnover_stats.csv",
) -> pd.DataFrame:
    """Table 2 of the paper (parametric column only)."""
    t = result.turnover.dropna()
    df = pd.DataFrame({
        "Metric": [
            "Average Daily Turnover",
            "Median Daily Turnover",
            "95% Turnover Quantile",
            "Days with > 1% Turnover",
            "Days with > 5% Turnover",
        ],
        "Value": [
            t.mean(),
            t.median(),
            t.quantile(0.95),
            float((t > 0.01).mean()),
            float((t > 0.05).mean()),
        ],
    })
    df.to_csv(path, index=False, float_format="%.6f")
    return df


def write_allocation_table(
    result: BacktestResult,
    path: Path = TBL_DIR / "allocation_summary.csv",
) -> pd.DataFrame:
    """Table 3 of the paper (parametric column only)."""
    W  = result.weights
    dW = W.diff().dropna()
    df = pd.DataFrame({
        "Avg Wt":      W.mean(),
        "Wt Vol":      W.std(ddof=1),
        "Time > 10%":  (W > 0.10).mean(),
        "Avg |Δw|":    dW.abs().mean().reindex(W.columns),
    })
    df.to_csv(path, float_format="%.6f")
    return df


def write_concentration_table(
    result: BacktestResult,
    path: Path = TBL_DIR / "concentration.csv",
) -> pd.DataFrame:
    """Table 4 of the paper (parametric column only)."""
    n_eff = 1.0 / (result.weights.pow(2).sum(axis=1) + 1e-12)
    df = pd.DataFrame({
        "Metric": ["Average N_eff", "Median N_eff",
                   "Min N_eff", "Max N_eff"],
        "Value":  [n_eff.mean(), n_eff.median(), n_eff.min(), n_eff.max()],
    })
    df.to_csv(path, index=False, float_format="%.6f")
    return df


def write_daily_outputs(
    result: BacktestResult,
    path: Path = TBL_DIR / "daily_backtest_output.csv",
) -> pd.DataFrame:
    """Per-day weights + pnl + regime label + selected K + turnover."""
    df = pd.concat([
        result.weights.add_prefix("w_"),
        result.pnl.rename("pnl"),
        result.cum_pnl.rename("cum_pnl"),
        result.K_history,
        result.tpl_label,
        result.tpl_max_prob,
        result.tpl_count,
        result.turnover,
    ], axis=1)
    df.index.name = "date"
    df.to_csv(path, float_format="%.8f")
    return df


# --------------------------------------------------------------------- #
# 4.  Regime analytics (Table 5 + Table 6 of the paper)
# --------------------------------------------------------------------- #
def regime_portfolio_table(
    result: BacktestResult,
    path: Path = TBL_DIR / "regime_portfolio_performance.csv",
) -> pd.DataFrame:
    """Portfolio performance broken down by template regime."""
    df = pd.DataFrame({
        "pnl":    result.pnl.astype(float),
        "regime": result.tpl_label.astype(float),
    }).dropna()
    grp = df.groupby("regime")["pnl"]

    out = pd.DataFrame({
        "Days":      grp.size().astype(int),
        "Frac":      grp.size() / grp.size().sum(),
        "Ann Mean":  grp.apply(ann_mean),
        "Ann Vol":   grp.apply(ann_vol),
        "Sharpe":    grp.apply(ann_sharpe),
        "Hit Rate":  grp.apply(hit_rate),
        "Max DD (within)": grp.apply(max_drawdown_simple),
    }).sort_index()
    out.index.name = "Regime"
    out.to_csv(path, float_format="%.6f")
    return out


def regime_asset_sharpe(
    result: BacktestResult,
    returns: pd.DataFrame,
    path: Path = TBL_DIR / "regime_asset_sharpe.csv",
) -> pd.DataFrame:
    """Per-asset annualized Sharpe inside each template regime."""
    R = returns.loc[result.pnl.index, ASSET_NAMES]
    reg = result.tpl_label.astype(float)
    df = pd.concat([R, reg.rename("regime")], axis=1).dropna()

    out = (df.groupby("regime")[ASSET_NAMES]
             .apply(lambda g: g.apply(ann_sharpe))
             .reindex(columns=ASSET_NAMES))
    out.index.name = "Regime"
    out.to_csv(path, float_format="%.6f")
    return out


# --------------------------------------------------------------------- #
# 5.  Plots — strategy diagnostics
# --------------------------------------------------------------------- #
def plot_cumulative_pnl_colored_by_regime(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 5))
    cum = result.pnl.cumsum()
    reg = result.tpl_label.fillna(-1).astype(int)
    unique_regs = sorted([r for r in reg.unique() if r >= 0])

    for i, r in enumerate(unique_regs):
        mask = reg == r
        ax.scatter(cum.index[mask], cum.values[mask],
                   s=10, color=_REGIME_COLORS[i % len(_REGIME_COLORS)],
                   label=f"Regime {chr(ord('A') + i)}", alpha=0.85)
    ax.set_title("Cumulative PnL Colored by Template Regime")
    ax.set_xlabel("Date"); ax.set_ylabel("Cumulative log return")
    ax.legend(loc="upper left", ncol=2)
    return _save(fig, "01_cumulative_pnl_colored_by_regime")


def plot_turnover(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(result.turnover.index, result.turnover.values,
            color=_PALETTE["parametric"], lw=0.9)
    ax.fill_between(result.turnover.index, 0, result.turnover.values,
                    color=_PALETTE["parametric"], alpha=0.25)
    ax.set_title("Daily Turnover — Wasserstein HMM")
    ax.set_ylabel("0.5 · ‖Δw‖₁")
    return _save(fig, "02_turnover_timeseries")


def plot_weights_stack(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 5))
    W = result.weights
    colors = [_ASSET_COLORS[c] for c in W.columns]
    ax.stackplot(W.index, [W[c].values for c in W.columns],
                 labels=W.columns, colors=colors, alpha=0.95)
    ax.set_title("Daily Portfolio Weights (Stacked)")
    ax.set_ylim(0, 1.001)
    ax.legend(loc="upper left", ncol=5)
    return _save(fig, "03_portfolio_weights_stacked")


def plot_neff(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 4))
    n_eff = 1.0 / (result.weights.pow(2).sum(axis=1) + 1e-12)
    ax.plot(n_eff.index, n_eff.values, color=_PALETTE["parametric"], lw=1.1)
    ax.axhline(n_eff.mean(), color=_PALETTE["accent"], ls="--", lw=0.9,
               label=f"mean = {n_eff.mean():.2f}")
    ax.set_title("Portfolio Effective # Positions (N_eff)")
    ax.set_ylabel("N_eff = 1 / Σ wᵢ²"); ax.legend()
    return _save(fig, "04_neff_timeseries")


def plot_K_history(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.step(result.K_history.index, result.K_history.values,
            color=_PALETTE["parametric"], lw=1.1, where="post")
    ax.set_title("Selected K Over Time  (Predictive Model-Order Selection)")
    ax.set_ylabel("K"); ax.set_ylim(bottom=0)
    return _save(fig, "05_K_history")


def plot_tpl_count(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.step(result.tpl_count.index, result.tpl_count.values,
            color=_PALETTE["parametric"], lw=1.1, where="post")
    ax.set_title("Number of Active Templates Over Time")
    ax.set_ylabel("# Templates")
    return _save(fig, "06_template_count")


def plot_tpl_label_timeline(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.plot(result.tpl_label.index, result.tpl_label.values,
            color=_PALETTE["parametric"], lw=1.0)
    ax.set_title("Dominant Template Label Over Time")
    ax.set_ylabel("Regime (hard label)")
    return _save(fig, "07_regime_label_timeline")


def plot_tpl_confidence(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.plot(result.tpl_max_prob.index, result.tpl_max_prob.values,
            color=_PALETTE["parametric"], lw=0.9)
    ax.fill_between(result.tpl_max_prob.index, 0,
                    result.tpl_max_prob.values,
                    color=_PALETTE["parametric"], alpha=0.2)
    ax.axhline(0.5, color="grey", ls="--", lw=0.8)
    ax.set_title("Max Template Posterior (Confidence Proxy)")
    ax.set_ylim(0, 1.01); ax.set_ylabel("max_g p_t,g")
    return _save(fig, "08_regime_confidence")


# --------------------------------------------------------------------- #
# 6.  Plots — regime analytics
# --------------------------------------------------------------------- #
def plot_regime_occupancy(reg_perf: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(9, 4))
    labels = [f"Regime {chr(ord('A') + i)}" for i in range(len(reg_perf))]
    ax.bar(labels, reg_perf["Frac"].values,
           color=[_REGIME_COLORS[i % len(_REGIME_COLORS)]
                  for i in range(len(reg_perf))], alpha=0.9)
    ax.set_ylabel("Fraction of OOS Days")
    ax.set_title("Template Regime Occupancy")
    for i, v in enumerate(reg_perf["Frac"].values):
        ax.text(i, v, f"{v:.0%}", ha="center", va="bottom", fontsize=10)
    return _save(fig, "09_regime_occupancy")


def plot_regime_sharpe(reg_perf: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(9, 4))
    labels = [f"Regime {chr(ord('A') + i)}" for i in range(len(reg_perf))]
    ax.bar(labels, reg_perf["Sharpe"].values,
           color=[_REGIME_COLORS[i % len(_REGIME_COLORS)]
                  for i in range(len(reg_perf))], alpha=0.9)
    ax.set_ylabel("Annualized Sharpe")
    ax.set_title("Portfolio Sharpe Ratio by Regime")
    for i, v in enumerate(reg_perf["Sharpe"].values):
        ax.text(i, v, f"{v:.2f}", ha="center",
                va="bottom" if v >= 0 else "top", fontsize=10)
    return _save(fig, "10_regime_portfolio_sharpe")


def plot_asset_sharpe_by_regime(asset_sharpe: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(11, 5))
    reg_labels = [f"Regime {chr(ord('A') + i)}" for i in range(len(asset_sharpe))]
    x = np.arange(len(asset_sharpe))
    m = len(ASSET_NAMES)
    bar_w = 0.8 / m
    offsets = (np.arange(m) - (m - 1) / 2.0) * bar_w

    for i, a in enumerate(ASSET_NAMES):
        ax.bar(x + offsets[i], asset_sharpe[a].values,
               width=bar_w, label=a, color=_ASSET_COLORS[a], alpha=0.9)

    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(reg_labels)
    ax.set_ylabel("Sharpe (Annualized)")
    ax.set_title("Asset Sharpe Ratios by Regime")
    ax.legend(ncol=len(ASSET_NAMES), loc="upper center",
              bbox_to_anchor=(0.5, 1.18))
    return _save(fig, "11_asset_sharpe_by_regime_bars")


def plot_asset_sharpe_heatmap(asset_sharpe: pd.DataFrame) -> Path:
    """A heatmap version of Table 6 — easier to read at a glance.

    Robust to two pathologies that can hit small-sample / synthetic data:
    * Regimes with too few observations → ``NaN`` Sharpes.
    * A single extreme value washing out the colour scale → clipped via
      the 95th absolute percentile (rare in real data but worth doing).
    """
    fig, ax = plt.subplots(figsize=(8, 4.5))
    reg_labels = [f"Regime {chr(ord('A') + i)}"
                  for i in range(len(asset_sharpe))]
    data = asset_sharpe.values.astype(float)
    mask_nan = ~np.isfinite(data)

    # Robust color scale: clip at 95th percentile of |data|, floored at 0.5
    abs_finite = np.abs(data[~mask_nan]) if (~mask_nan).any() else np.array([1.0])
    vmax = max(float(np.nanpercentile(abs_finite, 95)) if abs_finite.size else 1.0, 0.5)

    cmap = LinearSegmentedColormap.from_list(
        "rbg", ["#b91c1c", "#f1f5f9", "#1e3a8a"], N=256)
    plot_data = np.where(mask_nan, 0.0, data)
    im = ax.imshow(plot_data, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")

    # Grey overlay for NaN cells
    if mask_nan.any():
        nan_overlay = np.where(mask_nan, 1.0, np.nan)
        ax.imshow(nan_overlay, cmap="Greys", vmin=0, vmax=1.4, alpha=0.55,
                  aspect="auto")

    ax.set_xticks(range(len(ASSET_NAMES))); ax.set_xticklabels(ASSET_NAMES)
    ax.set_yticks(range(len(reg_labels)));  ax.set_yticklabels(reg_labels)
    ax.set_title("Annualized Asset Sharpe by Regime  (Heatmap)")

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if not np.isfinite(v):
                ax.text(j, i, "n/a", ha="center", va="center",
                        color="#374151", fontsize=10, style="italic")
                continue
            color = "white" if abs(v) > vmax * 0.6 else "#111827"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color=color, fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.04)
    return _save(fig, "12_asset_sharpe_by_regime_heatmap")


def plot_stacked_pnl_by_regime(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 5))
    df = pd.DataFrame({
        "pnl":    result.pnl.astype(float),
        "regime": result.tpl_label.astype(float),
    }).dropna()
    regs = sorted(df["regime"].unique())
    cums = {}
    for i, g in enumerate(regs):
        pnl_g = df["pnl"].where(df["regime"] == g, 0.0)
        cums[f"Regime {chr(ord('A') + i)}"] = pnl_g.cumsum()

    cum_df = pd.DataFrame(cums, index=df.index).ffill().fillna(0.0)
    ax.stackplot(cum_df.index, [cum_df[c].values for c in cum_df.columns],
                 labels=cum_df.columns,
                 colors=_REGIME_COLORS[:len(cum_df.columns)], alpha=0.92)
    ax.set_title("Stacked Cumulative PnL by Template Regime")
    ax.set_ylabel("Cumulative PnL contribution")
    ax.legend(loc="upper left")
    return _save(fig, "13_stacked_cumulative_pnl_by_regime")


# --------------------------------------------------------------------- #
# 7.  Plots — benchmark comparison & extras
# --------------------------------------------------------------------- #
def plot_cumulative_vs_benchmarks(
    result: BacktestResult, benchmarks: dict[str, pd.Series],
) -> Path:
    fig, ax = plt.subplots(figsize=(12, 6))
    cum_param = result.pnl.cumsum()
    cum_param -= cum_param.iloc[0]
    ax.plot(cum_param.index, cum_param.values,
            label="Wasserstein HMM + MVO", lw=2.2, color=_PALETTE["parametric"])

    styles = {"SPX Buy & Hold":  ("dotted", _PALETTE["spx"]),
              "Equal-Weight (20%)": ("dashed", _PALETTE["ew"]),
              "60/40 SPX/BOND":  ("dashdot", _PALETTE["sixty40"])}
    for name, series in benchmarks.items():
        c = series.cumsum(); c -= c.iloc[0]
        ls, col = styles.get(name, ("--", "#000"))
        ax.plot(c.index, c.values, label=name, lw=1.8, linestyle=ls, color=col)

    ax.set_title("Cumulative OOS Log Return — Strategy vs Benchmarks")
    ax.set_ylabel("Cumulative log return")
    ax.legend(loc="upper left")
    return _save(fig, "14_cumulative_vs_benchmarks")


def plot_drawdown_vs_benchmarks(
    result: BacktestResult, benchmarks: dict[str, pd.Series],
) -> Path:
    fig, ax = plt.subplots(figsize=(12, 5))
    series = {"Wasserstein HMM + MVO": result.pnl, **benchmarks}
    style_map = {"Wasserstein HMM + MVO": ("solid",  _PALETTE["parametric"], 2.2),
                 "SPX Buy & Hold":        ("dotted", _PALETTE["spx"], 1.8),
                 "Equal-Weight (20%)":    ("dashed", _PALETTE["ew"], 1.8),
                 "60/40 SPX/BOND":        ("dashdot",_PALETTE["sixty40"], 1.8)}
    for name, r in series.items():
        cum = r.cumsum()
        dd  = drawdown_series(cum)
        ls, col, lw = style_map.get(name, ("--", "#000", 1.5))
        ax.fill_between(dd.index, 0, dd.values,
                        alpha=0.15 if name == "Wasserstein HMM + MVO" else 0.0,
                        color=col)
        ax.plot(dd.index, dd.values, label=name, lw=lw, linestyle=ls, color=col)
    ax.set_title("Drawdown (Underwater) — Strategy vs Benchmarks")
    ax.set_ylabel("Drawdown (log-return units)")
    ax.legend(loc="lower left")
    return _save(fig, "15_drawdown_vs_benchmarks")


def plot_rolling_sharpe(
    result: BacktestResult, benchmarks: dict[str, pd.Series], window: int = 60,
) -> Path:
    fig, ax = plt.subplots(figsize=(12, 4.5))
    series = {"Wasserstein HMM + MVO": result.pnl, **benchmarks}
    style_map = {"Wasserstein HMM + MVO": ("solid",  _PALETTE["parametric"], 2.0),
                 "SPX Buy & Hold":        ("dotted", _PALETTE["spx"], 1.4),
                 "Equal-Weight (20%)":    ("dashed", _PALETTE["ew"], 1.4),
                 "60/40 SPX/BOND":        ("dashdot",_PALETTE["sixty40"], 1.4)}
    for name, r in series.items():
        mu  = r.rolling(window).mean()
        sd  = r.rolling(window).std(ddof=1)
        rs  = np.sqrt(TRADING_DAYS) * mu / (sd + 1e-12)
        ls, col, lw = style_map.get(name, ("--", "#000", 1.4))
        ax.plot(rs.index, rs.values, label=name, lw=lw, linestyle=ls, color=col)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_title(f"{window}-Day Rolling Sharpe")
    ax.set_ylabel("Annualized Sharpe")
    ax.legend(loc="upper left", ncol=2)
    return _save(fig, "16_rolling_sharpe")


def plot_rolling_vol(
    result: BacktestResult, benchmarks: dict[str, pd.Series], window: int = 60,
) -> Path:
    fig, ax = plt.subplots(figsize=(12, 4.5))
    series = {"Wasserstein HMM + MVO": result.pnl, **benchmarks}
    style_map = {"Wasserstein HMM + MVO": ("solid",  _PALETTE["parametric"], 2.0),
                 "SPX Buy & Hold":        ("dotted", _PALETTE["spx"], 1.4),
                 "Equal-Weight (20%)":    ("dashed", _PALETTE["ew"], 1.4),
                 "60/40 SPX/BOND":        ("dashdot",_PALETTE["sixty40"], 1.4)}
    for name, r in series.items():
        v = r.rolling(window).std(ddof=1) * np.sqrt(TRADING_DAYS)
        ls, col, lw = style_map.get(name, ("--", "#000", 1.4))
        ax.plot(v.index, v.values, label=name, lw=lw, linestyle=ls, color=col)
    ax.set_title(f"{window}-Day Rolling Volatility (Annualized)")
    ax.set_ylabel("Volatility")
    ax.legend(loc="upper left", ncol=2)
    return _save(fig, "17_rolling_volatility")


def plot_monthly_returns_heatmap(result: BacktestResult) -> Path:
    """A calendar heatmap of monthly returns for the parametric strategy."""
    fig, ax = plt.subplots(figsize=(10, 4.2))
    monthly = result.pnl.resample("ME").sum()
    if monthly.empty:
        return _save(fig, "18_monthly_returns_heatmap")
    df = monthly.to_frame("ret")
    df["Year"]  = df.index.year
    df["Month"] = df.index.month
    pivot = df.pivot(index="Year", columns="Month", values="ret")
    pivot = pivot.reindex(columns=range(1, 13))

    data = pivot.values
    finite = np.abs(data[np.isfinite(data)])
    vmax = max(float(np.nanpercentile(finite, 95)) if finite.size else 0.01, 0.005)
    cmap = LinearSegmentedColormap.from_list(
        "rbg", ["#b91c1c", "#f1f5f9", "#1e3a8a"], N=256)
    plot_data = np.where(np.isfinite(data), data, 0.0)
    im = ax.imshow(plot_data, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    # Grey out NaN months (e.g. partial start/end years)
    nan_mask = ~np.isfinite(data)
    if nan_mask.any():
        ax.imshow(np.where(nan_mask, 1.0, np.nan), cmap="Greys",
                  vmin=0, vmax=1.4, alpha=0.55, aspect="auto")
    ax.set_xticks(range(12))
    ax.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Monthly Log Returns — Wasserstein HMM Strategy")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if not np.isfinite(v):
                continue
            color = "white" if abs(v) > vmax * 0.55 else "black"
            ax.text(j, i, f"{v:.1%}", ha="center", va="center",
                    color=color, fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04)
    return _save(fig, "18_monthly_returns_heatmap")


def plot_strategy_drawdown(result: BacktestResult) -> Path:
    fig, ax = plt.subplots(figsize=(11, 4))
    dd = drawdown_series(result.pnl.cumsum())
    ax.fill_between(dd.index, 0, dd.values,
                    color=_PALETTE["accent"], alpha=0.3)
    ax.plot(dd.index, dd.values, color=_PALETTE["accent"], lw=1.0)
    ax.set_title("Wasserstein HMM Strategy — Drawdown Path")
    ax.set_ylabel("Drawdown")
    return _save(fig, "19_strategy_drawdown")


# --------------------------------------------------------------------- #
# 8.  Top-level orchestrator
# --------------------------------------------------------------------- #
def export_all(
    result:  BacktestResult,
    returns: pd.DataFrame,
) -> dict[str, pd.DataFrame | Path | dict]:
    """Run every CSV + PNG export. Returns a dict of generated artefacts."""
    ensure_dirs()
    _apply_style()

    benchmarks = build_benchmarks(returns, result.pnl.index)

    # ---- Tables ---- #
    tables = {
        "performance":   write_performance_table(result, benchmarks),
        "turnover":      write_turnover_table(result),
        "allocation":    write_allocation_table(result),
        "concentration": write_concentration_table(result),
        "daily":         write_daily_outputs(result),
        "regime_perf":   regime_portfolio_table(result),
        "regime_assets": regime_asset_sharpe(result, returns),
    }

    # Also dump the benchmark daily returns for downstream use
    bench_df = pd.DataFrame(benchmarks)
    bench_df.index.name = "date"
    bench_df.to_csv(TBL_DIR / "benchmark_daily_returns.csv",
                    float_format="%.8f")

    # ---- Figures ---- #
    figures = {
        "cum_pnl_regime":     plot_cumulative_pnl_colored_by_regime(result),
        "turnover":           plot_turnover(result),
        "weights":            plot_weights_stack(result),
        "neff":               plot_neff(result),
        "K_history":          plot_K_history(result),
        "tpl_count":          plot_tpl_count(result),
        "tpl_timeline":       plot_tpl_label_timeline(result),
        "tpl_confidence":     plot_tpl_confidence(result),
        "regime_occupancy":   plot_regime_occupancy(tables["regime_perf"]),
        "regime_sharpe":      plot_regime_sharpe(tables["regime_perf"]),
        "asset_sharpe_bars":  plot_asset_sharpe_by_regime(tables["regime_assets"]),
        "asset_sharpe_heat":  plot_asset_sharpe_heatmap(tables["regime_assets"]),
        "stacked_pnl_regime": plot_stacked_pnl_by_regime(result),
        "vs_bench_cum":       plot_cumulative_vs_benchmarks(result, benchmarks),
        "vs_bench_dd":        plot_drawdown_vs_benchmarks(result, benchmarks),
        "rolling_sharpe":     plot_rolling_sharpe(result, benchmarks),
        "rolling_vol":        plot_rolling_vol(result, benchmarks),
        "monthly_heat":       plot_monthly_returns_heatmap(result),
        "strategy_dd":        plot_strategy_drawdown(result),
    }

    return {"tables": tables, "figures": figures, "benchmarks": benchmarks}
