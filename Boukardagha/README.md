# Wasserstein-HMM Regime-Aware Investing — Extended OOS

A modular Python reproduction of **Boukardagha (2026), *Explainable
Regime-Aware Investing***, extended to a **21-year out-of-sample window**
(Jan 2005 – Dec 2025) using a 1990-onwards price history.

The project keeps the paper's methodology (strictly-causal expanding-window
Gaussian HMM, predictive K-selection with monotone-K safety, 2-Wasserstein
template tracking with EMA updates, forward-return regime moments,
transaction-cost-aware long-only MVO) and replaces the paper's
non-parametric KNN benchmark with the three baselines requested in the
brief: **SPX Buy & Hold**, **Equal-Weight (20% each)**, and
**60/40 SPX/BOND**.

## Run

```bash
pip install -r requirements.txt
python main.py
```

The 1990–2025 price history is loaded from `asset_data.csv` (path
configurable in `config.py`). Outputs land in `output/`:

```
output/
├── figures/        # 19 PNG figures
└── tables/         # 8 CSV tables
```

## OOS start choice — *why January 3, 2005*

Three OOS-start candidates were considered against a 1990-01-02 → 2025-12-31
sample. The decision criteria were: (i) enough training history for the HMM
to learn reasonably stable regime moments, (ii) the OOS window must cover
the **2008 Global Financial Crisis** as a true out-of-sample test, and
(iii) the model should not enter the GFC having seen only a single benign
decade.

| Start | Train years | OOS years | GFC in OOS? | Decision |
|---|---|---|---|---|
| **2005-01-03** | 15 (1990–2004) | 21 (2005–2025) | ✅ | **chosen** |
| 2008-01-02 | 18 | 18 | ✅ (but model has only seen 90s+early-00s) | rejected |
| 2010-01-04 | 20 | 16 | ❌ | rejected |

**Why not 2008**: with the HMM's 252-day internal validation window plus
the 60-day feature warm-up, the model needs at least ~1.5 years of
in-sample history before it has anything to validate against. Starting
OOS in 2005 means the model enters the GFC having *already* experienced
the dot-com bust and recovery (2001–2003) in-sample — a fairer test of
its regime-detection ability than dropping it cold into the GFC.

**Why not 2010**: throws away the GFC, which is one of the highest-value
regime-detection observations in modern history.

**Why not earlier (1996–2000)**: the HMM would be trained on a small
sample of mostly benign 1990s data and asked to detect regimes it has
never seen. Starting 2005 ensures at least the 2001 recession is in-sample.

## Parameter and methodology deviations from the paper

Three deliberate changes were made for the longer sample. All are flagged
in `config.py` with `EXTENSION 1/2/3` and can be reverted to paper values
to recover the exact original specification.

### `max_regimes : 6 → 8`  (parameter)
The paper calibrates 5–6 regimes for a 2005–2023 training window with a
2-year OOS. A 1990–2025 sample plausibly contains more distinct macro
regimes (90s bull, dot-com bust, GFC, ZIRP, taper tantrum, COVID,
post-COVID inflation, 2022 rate-hike cycle). The K-selection rule
`PredLL − λ_K · K` is unchanged — this only widens the candidate ceiling,
it does not force higher K. If the data don't support K=8, K=8 is never
chosen.

### `g_max : 8 → 10`  (parameter)
Templates accumulate over time and a higher `max_regimes` raises the
number of HMM components a new template might spawn from. Increasing
`g_max` proportionally avoids artificially saturating the template pool
over a 21-year window. The W2 spawn threshold and EMA rate are unchanged.

### `f_hmm : 1 → 5`  (methodology — for tractability)
The paper refits the HMM EM **every day** on the full expanding window.
Over a 35-year history this is prohibitive (~10–15 hours on a laptop for
the full OOS). I add an `f_hmm` parameter to control the HMM-refit cadence
and default it to 5 (same as `f_k`, weekly).

**What still updates daily**: regime probabilities `p_K`, hard labels
`z_s`, conditional moments `μ_K, Σ_K`, the W2 component→template mapping,
the EMA template update, the aggregate template posterior `p_G`, and the
MVO weight `w_t`. The only thing held fixed for up to 4 days between
refits is the HMM's *parameters themselves* (transition matrix, emission
means/covariances). Those are estimated on tens of thousands of daily
observations — letting them age by ≤4 days is negligible compared to
estimation noise.

Setting `f_hmm = 1` in `config.py` recovers the paper's exact daily-refit
behaviour for an apples-to-apples cross-check.

### Everything else preserved exactly
`LAM=5, TC=0.0002, W_MAX=0.6`,  `VOL_WINDOW=60, MOM_WINDOW=20`,
`L_VAL=252, LAM_K=1, ETA_TPL=0.05, SPAWN_THRESH=2.5, MIN_REGIMES=5`,
the monotone-K safety, Ledoit–Wolf shrinkage, forward-return regime
moments, strict-causality everywhere, random seed = 42, full covariance,
`n_iter=300`.

## Wall-time

The 5,478-day OOS run takes roughly **2–3 hours** on a modern laptop with
the default settings (`f_hmm=5`, `n_iter=300`, joblib parallelism on
across cores). The paper-spec `f_hmm=1` would take **~10–15 hours**.

If you need to iterate faster while developing thesis extensions: drop
`HMM.n_iter` to 100 in `config.py` while prototyping — the EM almost
always converges well before 300 iterations.

## Files

| File                  | Role                                                                                       |
|-----------------------|--------------------------------------------------------------------------------------------|
| `config.py`           | Every hyperparameter, all three EXTENSIONS flagged.                                        |
| `data.py`             | CSV loader, log returns, features.                                                         |
| `wasserstein_hmm.py`  | Predictive K-selection, HMM fit, regime moments, 2-Wasserstein, template tracking.         |
| `backtest.py`         | CVXPY/OSQP MVO solver + strictly-causal expanding-window loop with periodic HMM refit.    |
| `reporting.py`        | Metrics, three benchmarks, regime analytics, 19 PNG figures, 8 CSV tables.                |
| `main.py`             | Orchestrator.                                                                              |
| `asset_data.csv`      | Price history, 1990-01-02 → 2025-12-31 (stocks/bonds/oil/gold/usd).                       |

## Reproducibility

All randomness is controlled by `HMM.random_state = 42`. The CSV is a
fixed local file, so the entire pipeline is deterministic.

## Outputs

### Tables (`output/tables/`)
- `performance_summary.csv` — Sharpe / Sortino / Ann mean & vol / Max DD / Hit rate vs all 3 benchmarks
- `turnover_stats.csv` — paper Table 2 (parametric column)
- `allocation_summary.csv` — paper Table 3
- `concentration.csv` — paper Table 4
- `regime_portfolio_performance.csv` — paper Table 5
- `regime_asset_sharpe.csv` — paper Table 6
- `daily_backtest_output.csv` — full per-day weights, PnL, regime label, K, turnover
- `benchmark_daily_returns.csv` — SPX B&H, EW, 60/40 daily returns aligned to OOS

### Figures (`output/figures/`)
Paper figures plus five extras: drawdown overlay vs all benchmarks,
60-day rolling Sharpe and volatility, asset-Sharpe-by-regime heatmap,
monthly-returns calendar heatmap, and a strategy-only underwater plot.
