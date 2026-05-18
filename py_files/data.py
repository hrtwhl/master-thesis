"""
data.py
-------
Data acquisition and feature engineering.

* `load_prices`     : pulls auto-adjusted close prices from Yahoo Finance.
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
`returns[:t-1]`, but ~700x cheaper for a 700-day OOS period.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from config import (ASSET_NAMES, CACHE_DIR, END_DATE, MOM_WINDOW, RUN,
                    START_DATE, TICKERS, VOL_WINDOW)


# --------------------------------------------------------------------- #
# 1. Prices
# --------------------------------------------------------------------- #
def load_prices(
    start: str | None = None,
    end:   str | None = None,
    use_cache: bool   = True,
) -> pd.DataFrame:
    """Download (or load cached) daily adjusted close prices.

    Returns
    -------
    DataFrame with one column per asset (using friendly names from
    `config.ASSET_NAMES`) and DatetimeIndex.
    """
    start = start or START_DATE
    end   = end   or END_DATE or pd.Timestamp.today().strftime("%Y-%m-%d")

    cache_file: Path = CACHE_DIR / f"prices_{start}_{end}.pkl"
    if use_cache and cache_file.exists():
        with cache_file.open("rb") as fh:
            return pickle.load(fh)

    # Lazy import: yfinance is slow to import and we want fast CLI help text.
    import yfinance as yf

    print(f"[data] Downloading prices for {ASSET_NAMES}  ({start} → {end}) …")
    raw = yf.download(
        list(TICKERS.values()),
        start=start,
        end=end,
        auto_adjust=True,
        progress=False,
    )["Close"]

    # Rename Yahoo symbol → friendly names and enforce column order
    inverse = {v: k for k, v in TICKERS.items()}
    raw = raw.rename(columns=inverse)[ASSET_NAMES]

    prices = raw.dropna(how="any").sort_index()
    print(f"[data] Got {len(prices):,} clean daily observations.")

    if use_cache:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with cache_file.open("wb") as fh:
            pickle.dump(prices, fh)

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
def load_all(use_cache: bool = True) -> dict:
    """Convenience: prices, returns, features, and the train/test split.

    Returns
    -------
    dict with keys:
        prices, returns, features, returns_train, returns_test, split_date
    """
    from config import SPLIT_DATE   # local import keeps import-cost low

    prices   = load_prices(use_cache=use_cache)
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
