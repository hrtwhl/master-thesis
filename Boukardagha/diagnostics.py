"""
diagnostics.py
==============
Performance and risk diagnostics.

All metrics are computed on daily log returns.  We follow Boukardagha's
notebook conventions:
    - Sharpe / Sortino annualised by sqrt(252).
    - Ann mean = mean * 252, Ann vol = std * sqrt(252).
    - "Max DD (wealth)" computed on cumulative WEALTH = cumprod(1 + r).
    - "Max DD (cum log)" matches Boukardagha's Table 1 cumulative-log
      drawdown convention.
"""
import numpy as np
import pandas as pd

from config import TRADING_DAYS


def ann_mean(r: pd.Series) -> float:
    return float(r.mean() * TRADING_DAYS)


def ann_vol(r: pd.Series) -> float:
    return float(r.std(ddof=1) * np.sqrt(TRADING_DAYS))


def ann_sharpe(r: pd.Series) -> float:
    s = r.std(ddof=1)
    if s <= 1e-12 or np.isnan(s):
        return np.nan
    return float(np.sqrt(TRADING_DAYS) * r.mean() / s)


def ann_sortino(r: pd.Series) -> float:
    down = r[r < 0]
    if len(down) == 0 or down.std(ddof=1) <= 1e-12:
        return np.nan
    return float(np.sqrt(TRADING_DAYS) * r.mean() / down.std(ddof=1))


def hit_rate(r: pd.Series) -> float:
    return float((r > 0).mean())


def calmar_ratio(r: pd.Series) -> float:
    dd = max_drawdown_from_returns(r)
    if dd >= 0 or np.isnan(dd):
        return np.nan
    return float(ann_mean(r) / abs(dd))


def max_drawdown_from_returns(r: pd.Series) -> float:
    """Max drawdown on cumulative WEALTH = (1+r).cumprod()."""
    wealth = (1.0 + r.fillna(0.0)).cumprod()
    dd = wealth / wealth.cummax() - 1.0
    return float(dd.min())


def max_drawdown_from_logsum(r: pd.Series) -> float:
    """Drawdown on cum_pnl = r.cumsum() (Boukardagha Table 1 convention)."""
    cum = r.cumsum()
    dd = cum - cum.cummax()
    return float(dd.min())


def turnover_from_weights(W: pd.DataFrame) -> pd.Series:
    """Daily L1/2 turnover."""
    dW = W.diff()
    return 0.5 * dW.abs().sum(axis=1)


def n_effective(W: pd.DataFrame) -> pd.Series:
    """Effective number of positions = 1 / sum(w^2)."""
    return 1.0 / (W.pow(2).sum(axis=1) + 1e-12)


# -----------------------------------------------------------------------
#  Aggregated tables
# -----------------------------------------------------------------------
def performance_summary(pnl_dict: dict) -> pd.DataFrame:
    """
    pnl_dict : {name: pd.Series of daily log returns}
    """
    rows = {}
    for name, r in pnl_dict.items():
        rows[name] = {
            "Ann Mean":         ann_mean(r),
            "Ann Vol":          ann_vol(r),
            "Sharpe":           ann_sharpe(r),
            "Sortino":          ann_sortino(r),
            "Calmar":           calmar_ratio(r),
            "Max DD (wealth)":  max_drawdown_from_returns(r),
            "Max DD (cum log)": max_drawdown_from_logsum(r),
            "Hit Rate":         hit_rate(r),
            "Days":             int(len(r)),
        }
    return pd.DataFrame(rows).T


def turnover_summary(W: pd.DataFrame) -> dict:
    to = turnover_from_weights(W).dropna()
    return {
        "Mean":         float(to.mean()),
        "Median":       float(to.median()),
        "95% Quantile": float(to.quantile(0.95)),
        "Days > 1%":    float((to > 0.01).mean()),
        "Days > 5%":    float((to > 0.05).mean()),
    }


def allocation_summary(W: pd.DataFrame) -> pd.DataFrame:
    dW = W.diff().dropna()
    return pd.DataFrame({
        "Average Weight":    W.mean(),
        "Weight Volatility": W.std(),
        "Time > 10%":        (W > 0.10).mean(),
        "Avg |dW|":          dW.abs().mean().reindex(W.columns),
    })


def concentration_summary(W: pd.DataFrame) -> dict:
    neff = n_effective(W)
    return {
        "Average N_eff": float(neff.mean()),
        "Median N_eff":  float(neff.median()),
    }


# -----------------------------------------------------------------------
#  Regime breakdown tables
# -----------------------------------------------------------------------
def performance_by_regime(pnl: pd.Series, regime: pd.Series) -> pd.DataFrame:
    df = pd.DataFrame({"pnl": pnl, "regime": regime}).dropna()
    if df.empty:
        return pd.DataFrame()
    grp = df.groupby("regime")["pnl"]
    out = pd.DataFrame({
        "Days":            grp.size().astype(int),
        "Frac of Sample":  grp.size() / grp.size().sum(),
        "Ann Mean":        grp.apply(ann_mean),
        "Ann Vol":         grp.apply(ann_vol),
        "Sharpe":          grp.apply(ann_sharpe),
        "Hit Rate":        grp.apply(hit_rate),
        "Max DD (within)": grp.apply(max_drawdown_from_returns),
    }).sort_index()
    return out


def asset_performance_by_regime(returns: pd.DataFrame,
                                regime: pd.Series) -> pd.DataFrame:
    """Annualised mean / vol / Sharpe of each asset within each regime."""
    df = returns.copy()
    df["regime"] = regime.reindex(df.index)
    df = df.dropna(subset=["regime"])
    rows = []
    for r_val, sub in df.groupby("regime"):
        for col in returns.columns:
            s = sub[col].dropna()
            rows.append({
                "regime":   r_val,
                "asset":    col,
                "Ann Mean": ann_mean(s),
                "Ann Vol":  ann_vol(s),
                "Sharpe":   ann_sharpe(s),
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # Use set_index + unstack instead of pivot to avoid a pandas
    # version compatibility issue with 'performance_warnings'.
    out = out.set_index(["regime", "asset"])
    out = out.unstack("asset")
    return out


def regime_transition_matrix(regime: pd.Series) -> pd.DataFrame:
    """Empirical day-on-day regime transition matrix (rows->cols)."""
    r = regime.dropna().astype(int)
    if len(r) < 2:
        return pd.DataFrame()
    pairs = pd.DataFrame({"from": r.iloc[:-1].values, "to": r.iloc[1:].values})
    tab = pd.crosstab(pairs["from"], pairs["to"], normalize="index")
    return tab


def regime_persistence(regime: pd.Series) -> pd.DataFrame:
    """Mean spell length and # spells per regime label."""
    r = regime.dropna().astype(int).values
    if len(r) < 2:
        return pd.DataFrame()
    spells = []
    cur, n = r[0], 1
    for x in r[1:]:
        if x == cur:
            n += 1
        else:
            spells.append((cur, n))
            cur, n = x, 1
    spells.append((cur, n))
    sp = pd.DataFrame(spells, columns=["regime", "length"])
    out = sp.groupby("regime")["length"].agg(["count", "mean", "median", "max"])
    out = out.rename(columns={"count": "n_spells",
                              "mean": "mean_length",
                              "median": "median_length",
                              "max": "max_length"})
    return out
