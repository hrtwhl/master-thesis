"""
features.py
===========
Feature construction.

Asset features (Boukardagha §4):
    x_t = [r_t ; sigma_t ; m_t]
where r_t is the vector of daily log returns, sigma_t is the 60-day
rolling std of returns, and m_t is the 20-day rolling mean of returns.

Macro features (our extension):
    Same recipe applied to all 7 macro variables.
    Price-like variables  -> log returns
    Level-like variables  -> first differences
    Then 60-day rolling std and 20-day rolling mean of those increments
    are appended, giving 7 * 3 = 21 macro features.

All features are STRICTLY CAUSAL: features at time t only use data
available up to t-1.
"""
import numpy as np
import pandas as pd

from config import VOL_WINDOW, MOM_WINDOW, MACRO_VAR_KIND


def build_asset_features(returns: pd.DataFrame,
                         vol_window: int = VOL_WINDOW,
                         mom_window: int = MOM_WINDOW) -> pd.DataFrame:
    """
    Build the 15-dimensional (N*3) asset feature matrix used by the
    market-regime HMM.

    Parameters
    ----------
    returns : pd.DataFrame
        Daily log returns of the N assets indexed by date.
    """
    vol = returns.rolling(vol_window).std()
    mom = returns.rolling(mom_window).mean()
    feats = pd.concat([returns, vol, mom], axis=1)
    feats.columns = (
        [f"ret_{c}" for c in returns.columns] +
        [f"vol_{c}" for c in returns.columns] +
        [f"mom_{c}" for c in returns.columns]
    )
    return feats.dropna()


def macro_increments(macro_levels: pd.DataFrame) -> pd.DataFrame:
    """
    Transform raw macro levels into stationary increments.

    price-like variables  -> log returns
    level-like variables  -> first differences

    Returns a DataFrame with one column per macro variable, indexed by
    date.  NaNs (from the first observation per series) are dropped.
    """
    out = {}
    for col in macro_levels.columns:
        kind = MACRO_VAR_KIND[col]
        s = macro_levels[col]
        if kind == "price":
            out[col] = np.log(s).diff()
        elif kind == "level":
            out[col] = s.diff()
        else:
            raise ValueError(f"Unknown macro variable kind for {col}: {kind}")
    df = pd.DataFrame(out).dropna()
    return df


def build_macro_features(macro_levels: pd.DataFrame,
                         vol_window: int = VOL_WINDOW,
                         mom_window: int = MOM_WINDOW) -> pd.DataFrame:
    """
    Build the 21-dimensional macro feature matrix used by the macro
    layer of the hierarchical model.

    Following the user's directive, the Boukardagha asset-feature
    recipe (raw increment, 60d rolling std, 20d rolling mean) is
    applied to all 7 macro variables AFTER the stationary increment
    transform in `macro_increments`.
    """
    incs = macro_increments(macro_levels)
    vol  = incs.rolling(vol_window).std()
    mom  = incs.rolling(mom_window).mean()
    feats = pd.concat([incs, vol, mom], axis=1)
    feats.columns = (
        [f"inc_{c}" for c in incs.columns] +
        [f"vol_{c}" for c in incs.columns] +
        [f"mom_{c}" for c in incs.columns]
    )
    return feats.dropna()
