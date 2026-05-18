# Wasserstein-HMM Regime-Aware Investing — Reproduction

A modular Python reproduction of **Boukardagha (2026), *Explainable
Regime-Aware Investing***. The project keeps the paper's methodology
(strictly-causal rolling Gaussian HMM with predictive K-selection,
2-Wasserstein template tracking, transaction-cost-aware MVO) but drops the
non-parametric KNN benchmark and adds a **60/40 SPX/BOND** benchmark
alongside the paper's SPX Buy & Hold and Equal-Weight baselines.

## Run

```bash
pip install -r requirements.txt
python main.py
```

Outputs land in `output/`:

```
output/
├── cache/          # cached raw Yahoo download
├── figures/        # 19 PNG figures
└── tables/         # 8 CSV tables
```

## Files

| File                  | Role                                                                                       |
|-----------------------|--------------------------------------------------------------------------------------------|
| `config.py`           | Every hyperparameter. Edit here.                                                           |
| `data.py`             | Yahoo download, log returns, feature matrix.                                               |
| `wasserstein_hmm.py`  | Predictive K-selection, HMM fit, regime moments (forward returns), 2-Wasserstein, templates.|
| `backtest.py`         | CVXPY/OSQP MVO solver + strictly-causal expanding-window backtest loop.                    |
| `reporting.py`        | Metrics, benchmarks, regime analytics, PNG figures, CSV tables.                            |
| `main.py`             | Orchestrator.                                                                              |

## Speed notes

The reproduction is materially faster than the reference notebook without
changing methodology. The main wins:

1. **Features computed once.** Rolling vol/momentum are strictly causal, so
   the full feature matrix is precomputed on the entire return history and
   sliced inside the daily loop. This alone eliminates a redundant
   ~700 × O(T) computation.
2. **Concatenated history precomputed.** Train/test returns are merged
   once; daily history is taken via `np.searchsorted` rather than
   `pd.concat`+sort on every iteration.
3. **Parallel K-selection.** Candidate-K HMM fits inside `select_K_predictive`
   run in parallel via `joblib` (`RUN.parallel_k=True`, `RUN.n_jobs=-1`).
4. **Cached template square roots.** `Σ^{1/2}` of each template is computed
   once per call to `map_components_to_templates`, not once per
   (template, component) pair.
5. **Vectorised template EMA update** with `np.einsum`.

The HMM EM (`hmmlearn`), the Ledoit–Wolf shrinkage (`scikit-learn`), the
matrix square roots (`numpy.linalg.eigh`), and OSQP are all already backed by
optimised C / Fortran, so further speed-ups from a C++ extension would be
marginal and not worth the maintenance cost.

## Reproducibility

All randomness is controlled by `HMM.random_state` (default `42`). Yahoo
prices are cached after first download in `output/cache/`, so re-runs are
fully deterministic.

## Outputs at a glance

### Tables (`output/tables/`)
- `performance_summary.csv` — Sharpe / Sortino / Max DD vs all benchmarks (Table 7 + extra columns)
- `turnover_stats.csv` — paper Table 2 (parametric column)
- `allocation_summary.csv` — paper Table 3
- `concentration.csv` — paper Table 4
- `regime_portfolio_performance.csv` — paper Table 5
- `regime_asset_sharpe.csv` — paper Table 6
- `daily_backtest_output.csv` — full per-day weights, PnL, regime label, K
- `benchmark_daily_returns.csv` — SPX B&H, EW, 60/40 daily returns aligned to OOS

### Figures (`output/figures/`)
Paper figures plus extras (drawdown overlay, rolling Sharpe/vol, asset-Sharpe
heatmap, monthly-return calendar heatmap, strategy underwater plot).

## Methodology fidelity

Every algorithmic choice exactly matches the paper:

- Strictly causal expanding-window estimation
- Predictive log-likelihood K-selection with monotone-K safety
- 2-Wasserstein distance between full Gaussians for template mapping
- EMA template updates with `η = 0.05`
- Forward-return regime moments + Ledoit–Wolf shrinkage
- Long-only MVO with γ=5, τ=0.0002, w_max=0.6
- Random seed = 42, HMM `n_iter = 300`, `covariance_type = "full"`

If you change defaults to test extensions for your thesis, do it in
`config.py` — every module reads from there.
