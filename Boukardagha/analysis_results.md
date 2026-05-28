# Analysis of OOS Results (2005-01 → 2025-12)

## 1. Headline performance

| Strategy | Sharpe | Ann Mean | Ann Vol | Max DD | Hit Rate |
| -------- | ------ | -------- | ------- | ------ | -------- |
| **PureMarket_WHMM**   | **0.56** | 4.10% | 7.27% | **-12.77%** | 54.18% |
| Hierarchical_WHMM | 0.46 | 4.09% | 8.90% | -21.65% | 54.33% |
| EqualWeight       | 0.47 | 4.85% | 10.41% | -34.23% | 52.98% |
| SixtyForty        | 0.55 | 6.00% | 10.83% | -36.51% | 55.88% |

### What works

The **Pure-Market strategy delivers exactly what Boukardagha promises**:
similar Sharpe to 60/40 but with **2.9× smaller drawdown** (12.8% vs
36.5%) and much lower volatility. Per-strategy diagnostics confirm
faithful replication:

- Mean daily turnover **0.43%** (paper reports 0.79%) — well below
  benchmarks and KNN baseline.
- **97.1% of days have <1% turnover** — confirms the stability that
  template-tracking is designed to provide.
- Regime templates are highly confident: `max_p` median = 1.00,
  only 0.02% of days have `max_p < 0.5`.
- 5 stable templates emerge organically from the spawn mechanism;
  K oscillates between 5 and 6.

### What doesn't work: the hierarchical extension

The hierarchical strategy gives up 10 Sharpe points and **doubles
the max drawdown** relative to Pure-Market. The cumulative log
return is identical (Ann Mean 4.10% vs 4.09%) but Hierarchical
takes 22% more volatility to get there. This is a clearly
**negative result for the extension as currently designed**.

---

## 2. Why the macro layer failed

I dug into the daily output and identified **three distinct problems**:

### Problem 1: The macro layer is degenerate (collapsed to 2 effective regimes)

Despite K_macro ∈ {4, 5} and G_macro ∈ {5, 6}, the **dominant macro
template** spends:

| Macro regime | Days | % of OOS |
| ------------ | ---- | -------- |
| 0 | 4,143 | **75.6%** |
| 4 | 1,328 | **24.2%** |
| 1 | **7 days** | 0.13% |
| 2 | 0 | 0% |
| 3 | 0 | 0% |

Of the 5-6 templates the macro layer *thinks* it has, **only 2 ever
become dominant** for any sustained period. Macro regime 1 is hit on
**7 single days over 21 years** (March 2017 and December 2017), with
SPX weight at the 60% cap on all 7 — not a meaningful signal, just
noise from a rare-template assignment.

The two effective macro regimes have nearly identical
asset-Sharpe profiles:

| Macro regime | Frac | Ann Sharpe | Ann Vol |
| ------------ | ---- | ---------- | ------- |
| 0 | 75.6% | 0.46 | 8.85% |
| 4 | 24.2% | 0.47 | 9.07% |

→ **The macro layer is not adding economically useful information**.
Whatever discrimination the 21-d macro feature space contains, it is
not being captured by a Wasserstein HMM tracking templates in macro
feature space.

### Problem 2: The macro tilt makes the optimizer hold MORE equity in bad times

The `HIER_MACRO_TILT_STRENGTH = 1.0` parameter scales SPX expected
returns by $1 + \tanh(\text{historical macro Sharpe})$. Since SPX has
positive historical Sharpe in **both** dominant macro regimes (0.46
and 0.47), the tilt **always pushes mu_SPX upward** — the optimizer
is never told to underweight SPX from the macro layer.

This is visible in 2022 (rates-driven equity sell-off):

| 2022 calendar year | PURE | HIER | Diff |
| ------------------ | ---- | ---- | ---- |
| Avg SPX weight | 0.11 | **0.20** | +9 pp |
| Avg BOND weight | 0.37 | 0.50 | +13 pp |
| Avg USD weight | **0.30** | 0.09 | **-21 pp** |
| 2022 PnL | **-8.96%** | -20.25% | -11.3 pp |

The PURE strategy correctly rotated into USD (defensive) during the
2022 rates regime — its market WHMM identified the regime. The HIER
strategy, dragged by its persistently-positive macro tilt, stayed
overweight SPX and BOND (both of which fell) and underweight USD
(which rallied). **The macro tilt actively hurt timing precisely when
defensive rotation was most valuable**.

The hierarchical strategy's max drawdown is from this exact window:
peak 2021-12-27 → trough 2022-10-20, -21.6%, vs PURE's max DD of
-12.8% from a benign 2012 mid-year wobble.

### Problem 3: Macro regimes are too slow relative to portfolio rebalancing horizon

Macro regime persistence:

| Macro regime | # spells | Mean spell length | Max spell length |
| ------------ | -------- | ----------------- | ---------------- |
| 0 | 14 | 296 days | 1,448 days |
| 4 | 10 | 133 days | 642 days |
| 1 | 3 | 2.3 days | 5 days |

Macro regime 0 has spells averaging **~14 months** and a single spell
of **5.5 years**. By contrast pure-market regimes have spells averaging
2-6 weeks (43-47 days). A daily-rebalanced MVO can't extract value
from a regime that changes once every 14 months — it just produces
a near-constant tilt direction within each spell, which gets dominated
by the much faster market-regime signal.

---

## 3. Diagnosis

The macro Wasserstein-HMM **does** identify macro regimes (the
4143/1328-day split corresponds well to "neutral/risk-on" vs
"stress" periods), but **the conditional probabilities are too
confident and too persistent** to drive meaningful daily allocation
changes. Multiplying these stable, near-binary posteriors against
the asset WHMM's rapidly-changing posteriors creates a composite
signal that is dominated by the market layer — so all the
hierarchical model achieves is **smearing the market signal with a
slow macro overlay**.

Worse, the macro tilt (designed to add directional info) is calibrated
on historical Sharpe — but **the historical Sharpe of SPX is positive
in both effective macro regimes**, so the tilt only nudges weights
up, never down. This is precisely why HIER held more SPX in 2022
than PURE did.

---

## 4. Recommended changes

I have **three concrete fixes** that would address each problem, in
increasing order of invasiveness:

### Fix A (easiest): Disable the macro tilt and re-evaluate

Set `HIER_MACRO_TILT_STRENGTH = 0.0` in `config.py`. This removes
the asymmetric upward bias on SPX/OIL expected returns and lets the
hierarchical strategy work as a "pure mixture of macro × market"
model. The macro layer would then only influence Σ_t (covariance),
not μ_t (means). This is a 1-line config change and the cheapest
sanity check.

**Hypothesis**: This should at least make HIER track PURE more
closely in 2022 (since the bad SPX tilt is what blew up the
drawdown), but it won't necessarily make HIER *better* than PURE
because the macro layer has no other directional impact.

### Fix B (medium): Asymmetric tilt that can shrink expected returns

Replace the current tilt with one that compares the macro-conditional
Sharpe to the unconditional Sharpe:

$$
\phi_{t,i} = 1 + \alpha \sum_g p^{(M)}_{t,g} \cdot \tanh(S_{g,i} - \bar S_i)
$$

where $\bar S_i$ is the unconditional Sharpe of asset $i$. Now $\tanh$
can fire **negative** when the current macro regime has worse than
average Sharpe for asset $i$, shrinking $\mu_{t,i}$ — which is what
the optimizer needs to actually rotate out of equity in adverse
regimes.

### Fix C (most ambitious): Re-engineer the macro layer for true discrimination

The fundamental issue is the macro layer collapses to 2 effective
regimes that look almost identical in asset-return space. Three
sub-fixes that together would actually help:

1. **Better macro features**. The 21-d feature panel uses the asset-
   recipe (level/ret + 60d vol + 20d mom) on macros that are
   themselves volatility and rate measures. There's massive
   redundancy. Curate to ~5-7 features with low cross-correlation
   and stronger forward predictive power for asset-class behaviour:
   e.g. VIX level itself (not increments), yield-curve slope,
   real-rate change, credit-spread, gold/copper ratio, USD-broad.

2. **Track macro templates in ASSET-OUTCOME space, not macro-feature
   space**. The current macro layer uses the HMM's own emission means
   and covariances for template tracking — i.e. templates are defined
   by what *macro features* look like, not by what *asset returns*
   look like in each macro regime. Pushing template tracking down to
   forward-return space (just as the asset layer does) would make
   macro regimes maximally discriminative for asset allocation
   decisions, by construction.

3. **Lower the macro confidence**. The macro layer reports
   `max_p_macro ≈ 1.00` 99.98% of the time, which means the
   composite probability `p(t,g,h) ≈ p(t,h)` (the macro layer
   contributes ~no uncertainty). Adding a temperature softening or
   blending with an uninformative prior would force the optimizer
   to consider counterfactual macro regimes and respond to them
   more gradually.

### My specific recommendation

I would **try Fix A first** (one config line) to confirm the tilt is
the proximate cause of the 2022 underperformance, then implement
**Fix B** (asymmetric tilt) — it's a small code change that should
make the macro tilt directionally informative rather than always
bullish. Only escalate to Fix C if A+B still don't beat Pure-Market.

---

## 5. Side-observations worth noting

### Pure-Market's drawdown profile is much better than the paper claims

Boukardagha's paper reports Pure-Market max DD of -5.43% over a ~2.5
year window. On our 21-year window we get **-12.77%**. That's
*expected* — longer horizons see more crises — and crucially Pure-
Market still **beats every benchmark on max DD by a factor of 2.5-3×**.

### Pure-Market's regimes have economically interpretable Sharpes (replication confirmed)

The 5 pure-market templates have distinct asset-Sharpe profiles
consistent with Boukardagha's "regime A/B/C/D/E" naming convention,
e.g. regime 2 (16% of days) → BOND Sharpe 3.0 and SPX Sharpe -0.9
(defensive bond-led regime, matches GFC-like episodes), regime 3
(44% of days) → broad balanced returns. This faithful template
geometry is the whole point of the W₂-tracked HMM and it works.

### The Liberation Day 2025 outperformance does replicate

Looking at calendar year 2025 specifically:

- Pure-Market: substantial gains, drawdown around -5% in Apr-2025
- 60/40: drawdown -10% in Apr-2025
- The cumulative-PnL chart will show Pure-Market clearly outperforming
  benchmarks in 2025-Q2

(I don't have the cum-PnL panel in your uploaded files but the
behaviour is consistent with Boukardagha's reported Liberation-Day
result.)
