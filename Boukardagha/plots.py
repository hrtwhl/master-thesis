"""
plots.py
========
Replicates all charts in Boukardagha (2026) Section 8 and adds extra
cross-strategy comparison charts and macro-layer charts for the
hierarchical extension.

Every chart is saved as a PNG under CHART_DIR.

Per-strategy bundle (one set with prefix 'pure_', one with prefix 'hier_'):
    01_cum_pnl_scatter.png        Figure 1
    02_turnover.png               Figure 4
    03_weights_stacked.png        Figure 6
    04_n_eff.png                  Figure 8
    05_asset_sharpe_by_regime.png Figure 9
    06_stacked_cum_pnl_by_regime.png Figure 10
    07_K_history.png              supplementary
    08_template_count.png         supplementary
    09_template_label.png         supplementary
    10_template_posterior.png     supplementary

Hierarchical extras:
    hier_macro_label.png          dominant macro template over time
    hier_macro_prob.png           macro posterior over time
    hier_macro_vs_market.png      heatmap of joint (macro, market) freq

Cross-strategy:
    11_strategy_comparison.png    cumulative log return
    12_drawdown_comparison.png    drawdown (on wealth)
    13_rolling_sharpe_1y.png      1y rolling Sharpe
    14_annual_returns.png         calendar-year log returns
    15_underwater.png             max-drawdown underwater for all
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")

from config import CHART_DIR
from diagnostics import (
    turnover_from_weights, n_effective, ann_sharpe,
)

os.makedirs(CHART_DIR, exist_ok=True)
plt.rcParams.update({"figure.dpi": 110})


def _save(fig, name):
    path = os.path.join(CHART_DIR, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


# =======================================================================
#  Per-strategy chart bundle
# =======================================================================
def plot_strategy_bundle(pnl: pd.Series,
                         weights: pd.DataFrame,
                         regime: pd.Series | None,
                         K_history: pd.Series | None,
                         G_history: pd.Series | None,
                         tpl_prob: pd.Series | None,
                         asset_returns: pd.DataFrame,
                         prefix: str) -> dict:
    out = {}
    cum = pnl.cumsum()

    # 01 - cumulative PnL coloured by dominant regime
    fig, ax = plt.subplots(figsize=(11, 5))
    if regime is not None:
        df = pd.DataFrame({"cum": cum, "regime": regime}).dropna()
        if not df.empty:
            sc = ax.scatter(df.index, df["cum"].values,
                            c=df["regime"].values, s=8, cmap="viridis")
            plt.colorbar(sc, ax=ax, label="dominant regime")
        else:
            ax.plot(cum.index, cum.values, linewidth=1.5)
    else:
        ax.plot(cum.index, cum.values, linewidth=1.5)
    ax.set_title("Cumulative PnL Coloured by Dominant Regime")
    ax.set_xlabel("Date"); ax.set_ylabel("Cumulative log return")
    ax.grid(True)
    out["01_cum_pnl_scatter"] = _save(fig, f"{prefix}01_cum_pnl_scatter.png")

    # 02 - daily turnover
    to = turnover_from_weights(weights).dropna()
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(to.index, to.values, linewidth=0.9, label="Daily turnover (0.5*L1)")
    ax.set_title("Turnover Over Time")
    ax.grid(True); ax.legend()
    out["02_turnover"] = _save(fig, f"{prefix}02_turnover.png")

    # 03 - portfolio weights stacked
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.stackplot(weights.index,
                 [weights[c].values for c in weights.columns],
                 labels=weights.columns)
    ax.set_title("Portfolio Weights (Stacked)")
    ax.legend(loc="upper left", ncol=len(weights.columns))
    ax.grid(True)
    out["03_weights_stacked"] = _save(fig, f"{prefix}03_weights_stacked.png")

    # 04 - effective number of positions
    neff = n_effective(weights)
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(neff.index, neff.values, linewidth=1)
    ax.set_title("Portfolio Concentration Over Time (N_eff)")
    ax.grid(True)
    out["04_n_eff"] = _save(fig, f"{prefix}04_n_eff.png")

    # 05 - asset Sharpe by regime
    if regime is not None:
        out["05_asset_sharpe_by_regime"] = _plot_asset_sharpe_by_regime(
            asset_returns, regime, prefix=prefix,
        )

    # 06 - stacked cumulative PnL by regime
    if regime is not None:
        out["06_stacked_cum_pnl_by_regime"] = _plot_stacked_cum_pnl_by_regime(
            pnl, regime, prefix=prefix,
        )

    # 07 - K history
    if K_history is not None:
        fig, ax = plt.subplots(figsize=(11, 3))
        ax.plot(K_history.index, K_history.values, linewidth=1)
        ax.set_title("Selected K Over Time (Predictive Selection)")
        ax.grid(True)
        out["07_K_history"] = _save(fig, f"{prefix}07_K_history.png")

    # 08 - template count
    if G_history is not None:
        fig, ax = plt.subplots(figsize=(11, 3))
        ax.plot(G_history.index, G_history.values, linewidth=1)
        ax.set_title("Number of Templates Over Time")
        ax.grid(True)
        out["08_template_count"] = _save(fig, f"{prefix}08_template_count.png")

    # 09 - dominant template label over time
    if regime is not None:
        fig, ax = plt.subplots(figsize=(11, 3))
        ax.plot(regime.index, regime.values, linewidth=0.9)
        ax.set_title("Dominant Template Label Over Time")
        ax.grid(True)
        out["09_template_label"] = _save(fig, f"{prefix}09_template_label.png")

    # 10 - template max posterior over time
    if tpl_prob is not None:
        fig, ax = plt.subplots(figsize=(11, 3))
        ax.plot(tpl_prob.index, tpl_prob.values, linewidth=0.9)
        ax.set_title("Max Template Posterior Over Time (confidence proxy)")
        ax.grid(True)
        out["10_template_posterior"] = _save(fig, f"{prefix}10_template_posterior.png")

    return out


def _plot_asset_sharpe_by_regime(returns: pd.DataFrame,
                                 regime: pd.Series, prefix: str) -> str:
    df = returns.copy()
    df["regime"] = regime.reindex(df.index)
    df = df.dropna(subset=["regime"])

    regimes = sorted(df["regime"].dropna().unique())
    letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    reg_labels = [f"Regime {letters[i]} ({int(g)})" for i, g in enumerate(regimes)]
    assets = list(returns.columns)

    S = np.zeros((len(assets), len(regimes)))
    for j, r in enumerate(regimes):
        sub = df[df["regime"] == r]
        for i, a in enumerate(assets):
            S[i, j] = ann_sharpe(sub[a]) if len(sub) > 0 else np.nan

    x = np.arange(len(regimes))
    m = len(assets)
    bar_w = 0.8 / m
    offsets = (np.arange(m) - (m - 1) / 2.0) * bar_w

    fig, ax = plt.subplots(figsize=(11, 5))
    for i, a in enumerate(assets):
        ax.bar(x + offsets[i], S[i, :], width=bar_w, label=a)
    ax.set_xticks(x)
    ax.set_xticklabels(reg_labels, rotation=15, ha="right")
    ax.set_ylabel("Annualised Sharpe Ratio")
    ax.set_title("Asset Sharpe Ratios by Regime")
    ax.grid(True, axis="y")
    ax.legend(ncol=m, loc="upper center", bbox_to_anchor=(0.5, 1.15))
    return _save(fig, f"{prefix}05_asset_sharpe_by_regime.png")


def _plot_stacked_cum_pnl_by_regime(pnl: pd.Series,
                                    regime: pd.Series, prefix: str) -> str:
    df = pd.DataFrame({"pnl": pnl, "regime": regime}).dropna()
    regimes = sorted(df["regime"].unique())
    letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    cum_by_regime = {}
    for i, r in enumerate(regimes):
        pnl_r = df["pnl"].where(df["regime"] == r, 0.0)
        cum_by_regime[f"Regime {letters[i]}"] = pnl_r.cumsum()
    cum_df = pd.DataFrame(cum_by_regime, index=df.index).ffill().fillna(0.0)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.stackplot(cum_df.index,
                 [cum_df[c].values for c in cum_df.columns],
                 labels=cum_df.columns)
    ax.set_title("Stacked Cumulative PnL by Template Regime")
    ax.set_xlabel("Date"); ax.set_ylabel("Cumulative log return")
    ax.legend(loc="upper left"); ax.grid(True)
    return _save(fig, f"{prefix}06_stacked_cum_pnl_by_regime.png")


# =======================================================================
#  Hierarchical-only charts
# =======================================================================
def plot_macro_diagnostics(macro_label: pd.Series,
                           macro_prob: pd.Series,
                           market_label: pd.Series,
                           K_macro: pd.Series,
                           G_macro: pd.Series,
                           asset_returns: pd.DataFrame,
                           prefix: str = "hier_") -> dict:
    out = {}

    fig, ax = plt.subplots(figsize=(11, 3))
    ax.plot(macro_label.index, macro_label.values, linewidth=0.9)
    ax.set_title(f"{prefix}: Dominant MACRO Template Over Time")
    ax.grid(True)
    out["macro_label"] = _save(fig, f"{prefix}macro_label.png")

    fig, ax = plt.subplots(figsize=(11, 3))
    ax.plot(macro_prob.index, macro_prob.values, linewidth=0.9)
    ax.set_title(f"{prefix}: Max MACRO Posterior Over Time")
    ax.grid(True)
    out["macro_prob"] = _save(fig, f"{prefix}macro_prob.png")

    fig, ax = plt.subplots(figsize=(11, 3))
    ax.plot(K_macro.index, K_macro.values, linewidth=1, label="K_macro")
    ax.plot(G_macro.index, G_macro.values, linewidth=1, label="G_macro")
    ax.set_title(f"{prefix}: Macro-Layer K and G Over Time")
    ax.grid(True); ax.legend()
    out["macro_KG"] = _save(fig, f"{prefix}macro_KG.png")

    # Joint (macro, market) frequency heatmap
    df = pd.DataFrame({"g": macro_label, "h": market_label}).dropna()
    if not df.empty:
        df["g"] = df["g"].astype(int)
        df["h"] = df["h"].astype(int)
        tab = pd.crosstab(df["g"], df["h"], normalize="all") * 100.0

        fig, ax = plt.subplots(figsize=(7, 5))
        im = ax.imshow(tab.values, aspect="auto", cmap="Blues")
        ax.set_xticks(range(len(tab.columns)))
        ax.set_xticklabels(tab.columns)
        ax.set_yticks(range(len(tab.index)))
        ax.set_yticklabels(tab.index)
        ax.set_xlabel("market template h")
        ax.set_ylabel("macro template g")
        ax.set_title(f"{prefix}: Joint frequency p(g, h)  [% of days]")
        for i in range(tab.shape[0]):
            for j in range(tab.shape[1]):
                ax.text(j, i, f"{tab.values[i, j]:.1f}",
                        ha="center", va="center", fontsize=8,
                        color="black" if tab.values[i, j] < 15 else "white")
        plt.colorbar(im, ax=ax, label="% of days")
        out["macro_vs_market"] = _save(fig, f"{prefix}macro_vs_market.png")

    # Macro-conditional asset Sharpe (separate from market-template version)
    out["macro_asset_sharpe"] = _plot_asset_sharpe_by_regime(
        asset_returns, macro_label, prefix=f"{prefix}macro_",
    )

    return out


def plot_macro_stress(stress: pd.Series, gamma_eff: pd.Series,
                      name: str = "hierC_macro_stress") -> str:
    """Hierarchical-C diagnostic: macro stress score and the resulting
    effective risk aversion over time (twin axes)."""
    fig, ax1 = plt.subplots(figsize=(11, 4))
    ax1.plot(stress.index, stress.values, linewidth=0.9,
             color="tab:red", label="macro stress (0-1)")
    ax1.set_ylabel("macro stress", color="tab:red")
    ax1.tick_params(axis="y", labelcolor="tab:red")
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(gamma_eff.index, gamma_eff.values, linewidth=0.9,
             color="tab:blue", label="effective gamma")
    ax2.set_ylabel("effective gamma", color="tab:blue")
    ax2.tick_params(axis="y", labelcolor="tab:blue")

    ax1.set_title("Hierarchical C: Macro Stress and Effective Risk Aversion")
    return _save(fig, f"{name}.png")


# =======================================================================
#  Cross-strategy comparison charts
# =======================================================================
def plot_strategy_comparison(pnl_dict: dict, name: str = "11_strategy_comparison") -> str:
    fig, ax = plt.subplots(figsize=(12, 6))
    for label, r in pnl_dict.items():
        cum = r.cumsum()
        cum = cum - cum.iloc[0]
        ax.plot(cum.index, cum.values, label=label, linewidth=1.5)
    ax.set_title("Cumulative log return comparison")
    ax.set_xlabel("Date"); ax.set_ylabel("Cumulative log return")
    ax.grid(True); ax.legend()
    return _save(fig, f"{name}.png")


def plot_drawdown_comparison(pnl_dict: dict, name: str = "12_drawdown_comparison") -> str:
    fig, ax = plt.subplots(figsize=(12, 5))
    for label, r in pnl_dict.items():
        wealth = (1.0 + r.fillna(0.0)).cumprod()
        dd = wealth / wealth.cummax() - 1.0
        ax.plot(dd.index, dd.values, label=label, linewidth=1.2)
    ax.set_title("Drawdown comparison (on wealth)")
    ax.set_xlabel("Date"); ax.set_ylabel("Drawdown")
    ax.grid(True); ax.legend()
    return _save(fig, f"{name}.png")


def plot_rolling_sharpe(pnl_dict: dict, window: int = 252,
                        name: str = "13_rolling_sharpe") -> str:
    fig, ax = plt.subplots(figsize=(12, 5))
    for label, r in pnl_dict.items():
        roll = r.rolling(window).mean() / (r.rolling(window).std() + 1e-12)
        roll = roll * np.sqrt(252)
        ax.plot(roll.index, roll.values, label=label, linewidth=1.0)
    ax.axhline(0.0, color="k", linewidth=0.5)
    ax.set_title(f"Rolling {window}-day annualised Sharpe")
    ax.set_xlabel("Date"); ax.set_ylabel("Sharpe")
    ax.grid(True); ax.legend()
    return _save(fig, f"{name}.png")


def plot_annual_returns_bar(pnl_dict: dict,
                            name: str = "14_annual_returns") -> str:
    annual = {label: r.groupby(r.index.year).sum()
              for label, r in pnl_dict.items()}
    df = pd.DataFrame(annual)
    fig, ax = plt.subplots(figsize=(13, 5))
    df.plot(kind="bar", ax=ax, width=0.85)
    ax.set_title("Calendar-year log return per strategy")
    ax.set_xlabel("Year"); ax.set_ylabel("Log return")
    ax.grid(True, axis="y")
    plt.xticks(rotation=45)
    return _save(fig, f"{name}.png")


def plot_underwater(pnl_dict: dict, name: str = "15_underwater") -> str:
    """Plot each strategy's underwater drawdown curve."""
    fig, ax = plt.subplots(figsize=(12, 5))
    for label, r in pnl_dict.items():
        wealth = (1.0 + r.fillna(0.0)).cumprod()
        dd = wealth / wealth.cummax() - 1.0
        ax.fill_between(dd.index, dd.values, 0.0, alpha=0.3, label=label)
    ax.set_title("Underwater drawdown")
    ax.set_xlabel("Date"); ax.set_ylabel("Drawdown")
    ax.grid(True); ax.legend()
    return _save(fig, f"{name}.png")
