"""
data.py
=======
Load and pre-process the asset and macro daily time series.

Boukardagha (2026) used yfinance to download adjusted close prices.
We instead read pre-supplied CSV files which extend the universe back
to 1990-01-02, giving us a much larger training window before the
2005-01 OOS start.

Outputs
-------
load_asset_returns() -> pd.DataFrame
    Daily log returns for the 5 paper assets (SPX, BOND, GOLD, OIL, USD)
    on the calendar intersection of the available data.

load_macro_panel()   -> pd.DataFrame
    Daily levels of the 7 Mulliner et al. macro variables, intersected
    with the asset-return calendar.
"""
import numpy as np
import pandas as pd

from config import (
    ASSET_CSV, MACRO_CSV, ASSET_NAME_MAP, ASSET_NAMES,
    MACRO_VARS, DATA_START, DATA_END,
)


def load_asset_prices() -> pd.DataFrame:
    """Load asset prices and rename columns to paper notation."""
    df = pd.read_csv(ASSET_CSV, parse_dates=["date"]).set_index("date").sort_index()
    df = df.rename(columns=ASSET_NAME_MAP)
    df = df[ASSET_NAMES]   # enforce ordering
    df = df.loc[(df.index >= pd.Timestamp(DATA_START)) &
                (df.index <= pd.Timestamp(DATA_END))]
    df = df.dropna(how="any")
    return df


def load_asset_returns() -> pd.DataFrame:
    """Daily log returns of the five assets."""
    prices = load_asset_prices()
    rets   = np.log(prices).diff().dropna()
    return rets


def load_macro_panel() -> pd.DataFrame:
    """Load the 7 macro state variables (levels)."""
    df = pd.read_csv(MACRO_CSV, parse_dates=["date"]).set_index("date").sort_index()
    df = df[MACRO_VARS]
    df = df.loc[(df.index >= pd.Timestamp(DATA_START)) &
                (df.index <= pd.Timestamp(DATA_END))]
    df = df.dropna(how="any")
    return df


def align_calendars(*dfs):
    """Intersect the indices of all supplied frames and reindex."""
    idx = dfs[0].index
    for d in dfs[1:]:
        idx = idx.intersection(d.index)
    return tuple(d.loc[idx].copy() for d in dfs)
