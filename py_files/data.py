"""
data.py
-------
Data acquisition and feature engineering.

Prices are loaded from the local CSV at `config.CSV_PATH`. The CSV is
expected to have a `date` column plus one column per asset, named per
`config.CSV_COLUMN_MAP`. Any extra columns are silently ignored.

* `load_prices`     : reads and aligns the price CSV.
* `compute_returns` : daily log returns.
* `build_features`  : the 3N-dimensional feature vector x_t = [r_t, σ_t, m_t]
                      used by the Gaussian HMM. Strict causality is preserved
                      because all rolling windows look strictly backwards.

Speed note
~~~~~~~~~~
The backtest loop in `backtest.py` calls `build_features(returns)` **once**
on the entire return history and slices the result inside the loop, instead
of re-computing features at every OOS date as in the original notebook.
Because every rolling window is causal, slicing the precomputed feature
matrix up to `t-1` is mathematically identical to recomputing features on
`returns[:t-1]`, but ~1000× cheaper for a 5000-day OOS period.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from config import (ASSET_NAMES, CSV_COLUMN_MAP, CSV_PATH, END_DATE,
                    MOM_WINDOW, START_DATE, VOL_WINDOW)


# --------------------------------------------------------------------- #
# 1. Prices
# --------------------------------------------------------------------- #
def load_prices(
    start: str | None = None,
    end:   str | None = None,
    path:  Path | None = None,
) -> pd.DataFrame:
    """Read prices from `config.CSV_PATH` (or an override path).

    The CSV must have a `date` column and the source columns named in
    `CSV_COLUMN_MAP` (extras are ignored). Dates are parsed, the frame
    is sorted, and any rows with missing data are dropped.

    Parameters
    ----------
    start, end
        Optional date filters. Default to `config.START_DATE` and
        `config.END_DATE` (None ⇒ last available row).
    path
        Optional CSV path override.

    Returns
    -------
    DataFrame indexed by date, columns = `ASSET_NAMES`, no NaNs.
    """
    start = start or START_DATE
    end   = end   or END_DATE
    path  = path  or CSV_PATH

    if not path.exists():
        raise FileNotFoundError(
            f"CSV data file not found at {path}. Update config.CSV_PATH."
        )

    df = pd.read_csv(path, parse_dates=["date"]).set_index("date").sort_index()

    missing = [src for src in CSV_COLUMN_MAP if src not in df.columns]
    if missing:
        raise ValueError(
            f"CSV at {path} is missing expected columns: {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    prices = df.rename(columns=CSV_COLUMN_MAP)[ASSET_NAMES]

    if start:
        prices = prices.loc[prices.index >= pd.Timestamp(start)]
    if end:
        prices = prices.loc[prices.index <= pd.Timestamp(end)]

    prices = prices.dropna(how="any")
    print(f"[data] Loaded {len(prices):,} daily observations "
          f"({prices.index.min().date()} → {prices.index.max().date()}).")
    return prices


# --------------------------------------------------------------------- #
# 2. Returns
# --------------------------------------------------------------------- #
def compute_returns(prices: pd.DataFrame) -> pd.DataFrame:
    """Daily log returns r_t = log(P_t) - log(P_{t-1})."""
    return np.log(prices).diff().dropna()


# --------------------------------------------------------------------- #
# 3. Features
# --------------------------------------------------------------------- #
def build_features(
    returns: pd.DataFrame,
    vol_window: int = VOL_WINDOW,
    mom_window: int = MOM_WINDOW,
) -> pd.DataFrame:
    """Concatenate raw returns, rolling vol, and rolling momentum.

    x_t = [r_t ;  σ_t ;  m_t]   ∈ ℝ^{3N}

    All windows are strictly backwards-looking, so the row indexed by date
    `s` uses information only up to `s` inclusive.
    """
    vol = returns.rolling(vol_window).std()
    mom = returns.rolling(mom_window).mean()
    feats = pd.concat([returns, vol, mom], axis=1)
    feats.columns = (
        [f"ret_{c}" for c in returns.columns]
        + [f"vol_{c}" for c in returns.columns]
        + [f"mom_{c}" for c in returns.columns]
    )
    return feats.dropna()


# --------------------------------------------------------------------- #
# 4. One-shot loader
# --------------------------------------------------------------------- #
def load_all() -> dict:
    """Convenience: prices, returns, features, and the train/test split.

    Returns
    -------
    dict with keys:
        prices, returns, features, returns_train, returns_test, split_date
    """
    from config import SPLIT_DATE   # local import keeps import-cost low

    prices   = load_prices()
    returns  = compute_returns(prices)
    features = build_features(returns)

    split = pd.Timestamp(SPLIT_DATE)
    returns_train = returns.loc[returns.index <  split]
    returns_test  = returns.loc[returns.index >= split]

    return dict(
        prices=prices,
        returns=returns,
        features=features,
        returns_train=returns_train,
        returns_test=returns_test,
        split_date=split,
    )
