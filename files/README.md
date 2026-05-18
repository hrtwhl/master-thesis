# Wasserstein-HMM Portfolio Strategy — R + C++ port

Implementation of Boukardagha (2026), *Explainable Regime-Aware Investing*,
ported from the author's Python notebook to R with C++ kernels via
RcppArmadillo. Designed to reproduce the paper's June 2023 → end-2025 results
as a first milestone, and to be extended to a longer OOS window thereafter.

## Project layout

```
proj/
├── main.R                      # Orchestrator
├── R/
│   ├── 01_data.R               # Loading, log returns, features, split
│   ├── 02_hmm_engine.R         # HMM wrappers + predictive K + MVO QP
│   ├── 03_strategy.R           # Strategy backtest + benchmarks
│   └── 04_diagnostics.R        # Metrics, tables, plots
├── src/
│   └── wasserstein_hmm.cpp     # HMM EM, W2, Ledoit-Wolf (Rcpp + RcppArmadillo)
├── data/
│   └── asset_data_yf.csv       # ^GSPC, IEF, GLD, USO, UUP (2007 → 2025)
└── output/                     # CSVs and PNGs land here
```

## Dependencies

R packages: `Rcpp`, `RcppArmadillo`, `Matrix`, `osqp`, `caTools`, `abind`.

```r
install.packages(c("Rcpp", "RcppArmadillo", "Matrix",
                   "osqp", "caTools", "abind"))
```

You also need a C++17 toolchain (gcc/clang) accessible to R — i.e. a working
`Rcpp::sourceCpp()`. Test with:

```r
Rcpp::cppFunction('int test() { return 42; }')
test()    # should print 42
```

## Run

From the project root:

```bash
Rscript main.R
```

or interactively, `setwd("proj/")` then `source("main.R")`. The first run
takes 30–60 s to compile the C++ kernel; subsequent runs cache.

Expected runtime for the paper's OOS window (~650 days, K=5–6, 5 assets):
**~1–3 minutes** depending on machine.

## Configuration

All knobs live in the `CFG` list at the top of `main.R`:

| Parameter        | Value      | Meaning                                      |
|------------------|------------|----------------------------------------------|
| `split_date`     | 2023-06-03 | OOS start (matches paper)                    |
| `vol_window`     | 60         | Trailing-volatility feature window           |
| `mom_window`     | 20         | Trailing-momentum feature window             |
| `K_min, K_max`   | 5, 6       | Predictive-K search range (paper default)    |
| `F_K`            | 5          | Re-select K every 5 trading days             |
| `L_val`          | 252        | Validation slice for predictive K            |
| `lambda_K`       | 1.0        | Complexity penalty in K selection            |
| `em_iter`        | 300        | Max EM iterations                            |
| `G_max`          | 8          | Max number of persistent templates           |
| `eta_tpl`        | 0.05       | Template EMA learning rate                   |
| `spawn_thresh`   | 2.5        | W2 distance threshold to spawn new template  |
| `lambda`         | 5.0        | MVO risk-aversion                            |
| `tc`             | 0.0002     | L1 turnover penalty (2bp per unit reb.)      |
| `w_max`          | 0.6        | Per-asset weight cap                         |

## Outputs

`output/` will contain after a run:

- **CSV tables** — `performance_metrics.csv`, `turnover_stats.csv`,
  `allocation_summary.csv`, `concentration.csv`, `regime_performance.csv`,
  `asset_sharpe_by_regime.csv`, `daily_weights.csv`, `daily_pnl.csv`,
  `K_history.csv`, `tpl_label_history.csv`
- **PNG plots** — cumulative comparison, cumulative colored by regime,
  stacked weights, turnover, N_eff, K-history, drawdown comparison,
  rolling Sharpe, per-regime asset Sharpe bars, stacked regime attribution.

## Implementation notes — where the code differs from the paper / notebook

These are flagged inline as well; centralized here for the thesis writeup.

1. **K can decrease.** The notebook hard-codes a non-decreasing ratchet
   (`K_curr = max(K_curr, K_candidate)`). The paper says K is selected each
   period and "can go up or down." I follow the paper. A ratchet ossifies K
   at the first noisy validation slice and defeats predictive selection.

2. **Regime moments use forward returns** `r_{s+1} | z_s = k`. The paper text
   is silent on the timing; the notebook does this; it is also the right
   object for one-step-ahead forecasting. Followed.

3. **KNN benchmark dropped**, replaced by daily-rebalanced 60/40
   stocks/bonds. SPX buy-and-hold and equal-weight kept per paper.

4. **HMM implementation.** The Python uses `hmmlearn`; we re-implement
   Baum-Welch with full covariance in C++. K-means++ init is seeded but
   does not exactly trace `hmmlearn`'s sklearn-K-means RNG path, so
   numerical equivalence to the digit is not expected. Target tolerance:
   Sharpe within ~5% of the paper's 2.18.

## Extending to the long sample

When you move to `asset_data.csv` (1990-, pre-cleaned to drop weekends):

- Push `split_date` to ~2005–2010 so 15–20 years of feature/template
  burn-in are available.
- Consider lowering `eta_tpl` from 0.05 toward 0.01 — slower template
  drift over a multi-decade window keeps regime identities economically
  stable.
- `K_min, K_max` stays at 5, 6 — agree that 10 is too many; the dynamic
  templates already absorb evolving regime character without needing
  extra HMM components.

These changes are configuration only, no code changes required.
