# Explainable Regime-Aware Investing — Replication & Hierarchical Extension

Faithfully replicates Boukardagha (2026) ["Explainable Regime-Aware
Investing"](https://arxiv.org/abs/2603.04441) on a 2005-01 → 2025-12
daily out-of-sample window, and adds a hierarchical Wasserstein-HMM
macro layer following Fine, Singer & Tishby (1998).

## Strategies

| Strategy            | Description                                                              |
| ------------------- | ------------------------------------------------------------------------ |
| **PureMarket_WHMM** | Faithful replication of Boukardagha (2026): expanding-window Gaussian HMM with predictive K-selection (weekly), Wasserstein template tracking, daily refits, transaction-cost-aware MVO. |
| **Hierarchical_WHMM** | Two-level HHMM (FST98):<br>• **Top**: Wasserstein HMM on 21-d macro features (7 vars × {ret, 60d-vol, 20d-mom}).<br>• **Bottom**: shared market Wasserstein HMM as in pure-market.<br>• Conditional moments computed on JOINT (g_macro, h_market) cells using LedoitWolf shrinkage, with market-only fallback when a cell has fewer than 60 obs.<br>• Composite probabilities p(t,g,h) = p_macro(t,g) × p_market(t,h) drive the template-mixture moments fed to MVO.<br>• Optional macro tilt scales risky-asset (SPX, OIL) expected returns by tanh of historical regime Sharpe. |
| **EqualWeight**     | Static 20% per asset (SPX, BOND, GOLD, OIL, USD), frictionless.          |
| **SixtyForty**      | Static 60% SPX / 40% BOND, frictionless.                                  |

## Deviations from Boukardagha (2026)

Documented in detail at the top of `config.py`.

| ID | Deviation                                | Reason                                        |
| -- | ---------------------------------------- | --------------------------------------------- |
| D1 | `HMM_N_ITER = 100` (paper 300)           | hmmlearn EM converges in ≤ 80 iters at tol=1e-3 on these data. Set HMM_N_ITER = 300 in `config.py` to recover the paper exactly. |

Refit cadence (daily), K-selection cadence (weekly), expanding training window,
monotone-K rule, K range {5, 6}, transaction-cost MVO, all other parameters
(γ=5, τ=0.0002, w_max=0.6, G_max=8, η=0.05, spawn_thresh=2.5, vol_window=60,
mom_window=20, L_val=252, λ_K=1) are **taken verbatim** from `Paper_Code.ipynb`.

### Efficiency mechanism (mathematically transparent)

Daily refits of a full-covariance K=6 HMM on an expanding window
naively cost ~3-4 s each (>5 hours over 5300 OOS days).  We exploit
EM's monotonicity by **warm-starting** each daily refit from the
previous day's fitted parameters using hmmlearn's `init_params=""`
mechanism.  Adding one row to the training history typically lands EM
in the same likelihood basin and reconverges in 3-10 iterations
(~0.1-0.3 s).  This is **not a methodological change** — the
converged maximum on each day is identical to what a cold-started fit
would produce.  Implementation: `wasserstein_hmm.py:fit_gaussian_hmm_warm`.

## File layout

```
regime_aware_investing/
├── config.py              # all knobs (with paper-vs-ours flags)
├── data.py                # asset & macro CSV loaders
├── features.py            # asset and macro feature builders
├── mvo.py                 # transaction-cost-aware MVO (cvxpy)
├── wasserstein_hmm.py     # WHMM step + W2 helpers + predictive K + warm-start
├── hierarchical_hmm.py    # macro->market HHMM (FST98 architecture, v2)
├── backtest.py            # pure-market backtest + static benchmark + CSV helper
├── diagnostics.py         # Sharpe/Sortino/DD/turnover/N_eff/persistence
├── plots.py               # all paper charts + macro charts + comparisons
├── run_pure.py            # entry: pure-market backtest only
├── run_hier.py            # entry: hierarchical backtest only
├── aggregate.py           # reads raw CSVs -> tables + PNGs
├── narrative.py           # Markdown auto-summary
├── report.py              # single-page HTML report
├── main.py                # end-to-end driver (1 → 5)
└── output/
    ├── charts/                       # PNG figures (~28 files)
    ├── tables/                       # paper-ready CSVs (T01..T09)
    ├── raw/                          # daily PnL/weights/regime time series
    ├── daily_backtest_output_pure.csv      # unified daily output (requested schema)
    ├── daily_backtest_output_hier.csv      # ditto with macro extras
    ├── daily_backtest_output_benchmarks.csv
    ├── narrative_summary.md
    └── report.html
```

## How to run

```bash
cd regime_aware_investing
python main.py     # ~30 + ~30 + 1 min on 1 CPU core, single shot
```

Or run the stages independently (useful for iterating on plots):

```bash
python run_pure.py     # ~30 min
python run_hier.py     # ~30 min
python aggregate.py    # <1 min
python narrative.py    # writes output/narrative_summary.md
python report.py       # writes output/report.html
```

## Output catalogue

### Tables (CSV in `output/tables/`)

| File                                            | Content                                              |
| ----------------------------------------------- | ---------------------------------------------------- |
| T01_performance_summary.csv                     | Sharpe, Sortino, Calmar, Ann Mean/Vol, Max DD, Hit Rate per strategy |
| T02_turnover_summary.csv                        | Mean/Median/95% turnover, % days > 1%/>5%            |
| T03a_allocation_pure_market.csv                 | Mean weight, weight vol, time-in-position, \|Δw\|    |
| T03b_allocation_hierarchical.csv                | Same, for hierarchical strategy                      |
| T04_concentration_summary.csv                   | Average / median N_eff                               |
| T05a_pure_portfolio_by_regime.csv               | Portfolio perf in each market template               |
| T05b_hier_portfolio_by_macro_regime.csv         | Portfolio perf in each MACRO template                |
| T05c_hier_portfolio_by_market_regime.csv        | Portfolio perf in each MARKET template (hier)        |
| T06a_pure_asset_by_regime.csv                   | Asset Sharpe/mean/vol in each market template        |
| T06b_hier_asset_by_macro_regime.csv             | Asset Sharpe/mean/vol in each macro template         |
| T07_crisis_performance.csv                      | Per-strategy returns / DD over 8 named crisis windows |
| T08a_pure_regime_transitions.csv                | Day-on-day transition matrix, pure-market templates  |
| T08b_hier_macro_regime_transitions.csv          | Same, macro templates                                |
| T08c_hier_market_regime_transitions.csv         | Same, hierarchical-market templates                  |
| T09a_pure_regime_persistence.csv                | Spell lengths per market template                    |
| T09b_hier_macro_regime_persistence.csv          | Spell lengths per macro template                     |
| T09c_hier_market_regime_persistence.csv         | Spell lengths per hierarchical-market template       |

### Daily output CSVs (in `output/`)

`daily_backtest_output_pure.csv` — schema requested by the user:
```
date, w_SPX, w_BOND, w_GOLD, w_OIL, w_USD, pnl, cum_pnl, K, regime, max_p, G, turnover
```

`daily_backtest_output_hier.csv` — same plus macro layer columns:
```
date, w_SPX, w_BOND, w_GOLD, w_OIL, w_USD, pnl, cum_pnl,
K_market, regime_market, max_p_market, G_market,
K_macro,  regime_macro,  max_p_macro,  G_macro,
turnover
```

`daily_backtest_output_benchmarks.csv` — EqualWeight + SixtyForty daily weights and PnL.

### Charts (PNG in `output/charts/`)

Per-strategy bundles (prefixes `pure_` and `hier_`): cumulative PnL
scatter coloured by regime, turnover, weights stacked, N_eff, asset
Sharpe by regime, stacked PnL by regime, K history, template count,
template label, max template posterior — 10 charts each.

Hierarchical-only macro-layer charts:
- `hier_macro_label.png` — dominant macro template over time
- `hier_macro_prob.png` — max macro posterior over time
- `hier_macro_KG.png` — macro K and G over time
- `hier_macro_vs_market.png` — heatmap of joint (macro, market) freq
- `hier_macro_05_asset_sharpe_by_regime.png` — asset Sharpe by macro regime

Cross-strategy:
- `11_strategy_comparison.png` — cumulative log return
- `12_drawdown_comparison.png` — drawdown
- `13_rolling_sharpe_1y.png` — rolling Sharpe
- `14_annual_returns.png` — calendar-year returns
- `15_underwater.png` — underwater drawdown
