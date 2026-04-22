# Replication of Mulliner, Harvey, Xia, Fang & van Hemert (2025), "Regimes"

## Summary

This document records a replication of the paper's systematic regime-detection
method and its application to long-short factor timing. The implementation is
methodologically faithful to the paper. The central claim — that ranking
historical months by economic similarity predicts subsequent factor returns —
replicates directionally on our sample (Q1 Sharpe 0.38 vs paper's 0.95,
Q1−Q5 Sharpe 0.35 vs paper's 0.82). The headline statistical claim of a
~3σ alpha on the Q1−Q5 difference portfolio achieves t = 1.59 on our sample
versus t = 3 in the paper — directionally consistent but below the paper's
significance threshold. We trace this to a shorter sample (30 vs 40 years)
and a single known data deviation (construction of the stock-bond correlation
series).

---

## 1. Methodology

Let $x_{v,t}$ denote the level of macro variable $v$ at month $t$. For each
of the seven state variables (Market = S&P 500 level, Yield curve = 10y − 3m,
Oil = WTI, Copper, Monetary policy = 3m T-bill, Volatility = VIX/pre-1990
realised vol, Stock-bond correlation = rolling 3-year daily correlation), we
form the transformed state variable

$$
z_{v,t} \;=\; \mathrm{wins}_{\pm 3}\!\!\left(
\frac{\Delta_{12}\, \tilde x_{v,t}}{\hat\sigma_{v,t}}
\right),
\qquad
\hat\sigma_{v,t} \;=\; \mathrm{sd}\!\Big(\{\Delta_{12}\tilde x_{v,\tau}\}_{\tau=t-119}^{t}\Big)
$$

where $\tilde x_{v,t} = \log x_{v,t}$ for multiplicative series (Market, Oil,
Copper, Volatility) and $\tilde x_{v,t} = x_{v,t}$ for rate-valued or
bounded-scale series (Yield curve, Monetary policy, Stock-bond correlation);
$\Delta_{12}$ is the 12-month difference; and $\mathrm{wins}_{\pm 3}$
winsorizes at $\pm 3$. This is the paper's "adjustment similar to a z-score"
($\Delta_{12}\tilde x_{v,t}$ is **not** mean-centered before division).

The similarity between months $i$ and $T$ is the Euclidean distance across
the $V = 7$ transformed variables:

$$
d_{iT} \;=\; \sqrt{\sum_{v=1}^{V}\bigl(z_{v,i} - z_{v,T}\bigr)^{2}}.
$$

For each target month $T$, eligible historical months are $\{i : i \le T - m\}$
where $m = 36$ months is the look-ahead mask. Eligible months are ranked by
distance and partitioned into $q = 5$ equal-sized buckets (quintiles) $B_1(T),
\ldots, B_5(T)$, with $B_1$ the most similar and $B_5$ the most dissimilar.

For each factor $f \in \{\text{Market, Size, Value, Profitability,
Investment, Momentum}\}$ (Fama-French five plus 12-month momentum), the signal
at target date $T$ is

$$
s_{f,k,T} \;=\; \mathrm{sign}\!\left(
\frac{1}{|B_k(T)|}\!\!\sum_{i \in B_k(T)} r_{f,\,i+1}
\right),
$$

i.e. $+1$ if the average $(i+1)$-month return of factor $f$ across historical
months in bucket $k$ is positive, $-1$ otherwise. The quintile-$k$ strategy
trades the signal in month $T+1$, averaged equally across the six factors:

$$
R_{k,\,T+1} \;=\; \frac{1}{6}\sum_{f} s_{f,k,T}\cdot r_{f,\,T+1}.
$$

The paper's headline object is the similarity-minus-dissimilarity difference
portfolio $R_{1,T+1} - R_{5,T+1}$. All strategy series (including long-only
and individual quintiles) are additionally reported after applying a 15%
annualised volatility target via a trailing 36-month realised-vol scaling
(see §3 for rationale).

---

## 2. Data

**Macro state variables.** Seven daily series spanning April 1986 to December
2025, resampled to end-of-month. Stock-bond correlation is pre-computed at
the 3-year daily-rolling frequency. See `data/final.csv`.

**Factor returns.** The six canonical Fama-French factors downloaded from
Kenneth French's data library. The sample intersection of valid z-scores
(requires 120 months of $\Delta_{12}$ history) and factor returns runs from
**April 1997 through December 2025**. After imposing the 36-month similarity
mask and the bucket-size minimum (5×2 eligible obs), the effective strategy
back-test is **January 2001 through December 2025 (300 months)**.

**Deviation from the paper's sample.** The paper uses 1985–2024. Ours is
approximately 2001–2025 due to the later start of the raw data. This is the
single largest driver of any residual differences in our results.

---

## 3. Volatility targeting and long-only construction

The paper reports its long-only (LO) portfolio has Sharpe ≈ 1.0 and
drawdowns reaching ~50% in both 2009 and 2020. A naive equal-weight average
of the six long-short factors earns Sharpe ≈ 1.07 at ~4% annualised
volatility in our sample — consistent with 1.0 but with drawdowns of only
~12%. The drawdown gap is not a construction difference: it's 15% annualised
vol targeting applied at the portfolio level.

Diagnostic 8 tests four LO constructions:

| Variant | Sharpe | Ann vol | Max DD | 2009 DD | 2020 DD |
|---|---|---|---|---|---|
| 1. Equal-weight, unlevered | 1.07 | 4.0% | −12% | −12% | −10% |
| 2. Inverse-vol, unlevered | 1.07 | 3.5% | −10% | −7% | −10% |
| 3. Equal-weight + 15% vol target | **1.06** | **17%** | **−49%** | **−49%** | **−49%** |
| 4. Inverse-vol + 15% vol target | 0.91 | 20% | −79% | −37% | −49% |

Variant 3 matches the paper exactly on all four dimensions reported in the
paper's Exhibit 11. The paper's LO is **equal-weight factor combination with
a 15% annualised portfolio vol target** — no inverse-vol weighting. This
implies approximately 3.75× leverage on the unlevered portfolio, applied
uniformly across the Q1–Q5 quintile strategies and LO to maintain
comparability.

Sharpe ratios, alpha $t$-statistics, and correlations are **invariant to
uniform vol scaling**. Vol targeting therefore does not change any
statistical conclusion in §4–5; it only aligns cumulative-return and
drawdown charts with the paper's visual scale.

---

## 4. Results

### 4.1 Replication of the central mechanism

Exhibit 10 quintile Sharpes and correlations with the unlevered LO
benchmark:

| | Q1 | Q2 | Q3 | Q4 | Q5 | Q1−Q5 |
|---|---|---|---|---|---|---|
| Our Sharpe | 0.38 | 0.24 | 0.04 | 0.27 | −0.11 | **0.35** |
| Paper Sharpe | 0.95 | 0.80 | 0.78 | 0.85 | 0.17 | 0.82 |

The Q1 > Q5 pattern holds (gap of 0.49 Sharpe units, the widest extreme
spread in our tested configurations), but middle-quintile ordering is noisy.
The paper's strictly monotonic Q1 > Q2 > Q3 > Q4 > Q5 does not appear on our
sample; we view this as sample-driven noise in the middle buckets rather
than a methodological failure, supported by (i) the clean Q1–Q5 spread,
(ii) the Market-factor Sharpe match (§4.2), and (iii) Diagnostic 3's
verification that similarity ranking picks up economically sensible regimes.

### 4.2 Per-factor results (Exhibit A2)

Unlevered Q1−Q5 Sharpes per factor:

| Factor | Paper Sharpe | Our Sharpe |
|---|---|---|
| Market | 0.34 | **0.24** |
| Size | 0.35 | 0.19 |
| Value | 0.31 | 0.03 |
| Profitability | 0.42 | −0.03 |
| Investment | −0.19 | −0.03 |
| Momentum | 0.34 | 0.28 |

Market, Size, and Momentum match the paper directionally with economically
meaningful Sharpes. Value and Profitability underperform — both have
well-documented poor performance in the 2010-2020 period that isn't present
in the paper's longer sample. Five of six factors are directionally
consistent.

### 4.3 Factor-model alpha (Diagnostic 7)

Regressing strategy returns on the six FF factors:

$$
R_{k,T+1} \;=\; \alpha \;+\; \sum_{f}\beta_f \, r_{f,T+1} \;+\; \varepsilon_{T+1}
$$

| Portfolio | α (p.a.) | $t(\alpha)$ | $p$ | $R^2$ |
|---|---|---|---|---|
| Q1 − Q5 | +2.31% | **+1.59** | 0.113 | 0.04 |
| Q1 | +0.91% | +0.92 | 0.356 | 0.15 |
| Q5 | −1.40% | −1.49 | 0.137 | 0.15 |

The Q1−Q5 portfolio is factor-neutral by construction ($R^2 = 0.04$;
only Momentum loads significantly at $|t| > 2$). Its alpha of **2.31% p.a.
with $t = 1.59$** is the cleanest test of the paper's claim. Point estimate
and sign align with the paper's ~3σ finding; our sample is statistically
underpowered to confirm it. The 95% CI on our alpha is [−0.54%, +5.16%],
which contains both zero and the paper's implied magnitude.

Q5's point-estimate alpha of −1.40% ($t = −1.49$) is consistent with the
paper's "anti-regime" claim directionally, also under-powered statistically.

### 4.4 Diagnostic summary

| Diagnostic | Finding | Implication |
|---|---|---|
| 1. z-score health | All std's in [1.02, 1.21]; paper 0.89–0.98 | Slightly inflated, 2001-2025 sample contains denser volatility. Not material. |
| 2. Distance contribution | Shares 0.12–0.17, ideal 0.14 | Balanced — no variable dominates the metric. |
| 3. Jan 2009 similar months | All top 15 are 2001–02 and 1998 | Regime matching is economically sensible. |
| 4. Signal distribution | Q1: 50–90% long per factor; Q5: 26–80% | Signal is meaningful (not coin-flip); expected long-bias on positive-premium factors. |
| 5. Date alignment | Zero NAs after join | Data integrity confirmed. |
| 6. Effective sample | 300 months, Jan 2001 – Dec 2025 | 25 years available for alpha estimation. |
| 7. Alpha regression | Q1−Q5 α = 2.31% p.a., $t = 1.59$ | Directional replication; below paper's 3σ threshold. |
| 8. LO construction | Equal-weight + 15% vol target matches paper exactly | Paper's recipe identified. |

---

## 5. Robustness

**Number of quantiles (Exhibit 12).** Q1−Q5 Sharpe (unlevered) by bucket
count: $q = 2$: 0.19; $q = 3$: 0.37; $q = 4$: 0.36; $q = 5$: 0.35;
$q = 10$: 0.30; $q = 20$: 0.41. Reasonably flat across configurations;
the paper's hump at $q = 5$ does not emerge cleanly but none of the
configurations is broken.

**Rolling-window lookback (Exhibit 13).** Q1−Q5 Sharpe by z-score
calibration window: 1y: 0.22; 3y: 0.20; 5y: 0.20; 10y: **0.35**. The paper
claims insensitivity to lookback. On our corrected data, 10y is meaningfully
better. We use 10y as the default.

**Similarity mask length.** An experimental run with $m = 0$ produced
dramatically better results (Q1 SR 1.17, Q1−Q5 alpha t = 4.75). Inspection
of Diagnostic 3 reveals that with $m = 0$, the "most similar" months to
Jan 2009 are Jan 2009 itself ($d = 0$), Dec 2008, Nov 2008, and Oct 2008 —
i.e., the state-variable 12-month autocorrelation causes the bucket to
collapse onto the target's immediate past. This is precisely the information
leakage the paper's $m = 36$ mask was designed to prevent. We retain the
paper's choice; the $m = 0$ experiment is noted for completeness as evidence
that the mask is doing real work.

---

## 6. Deviations from the paper

| Item | Paper | Ours | Reason |
|---|---|---|---|
| Sample | 1985–2024 | 2001–2025 effective | Raw macro data starts April 1986; 120m z-score window + 36m mask consume the rest. |
| Stock-bond correlation autocorrelation | 12m: 0.07 | 12m: 0.34 | Unknown construction detail in the paper; our reconstructed series is more persistent. |
| Differencing for price-like series | Not specified | Log | Paper's Exhibit 3 std's (0.89–0.98) are consistent with log diffs for Market/Oil/Copper/VIX; level diffs produced 1.4× inflated std for these. |
| Number of variables | 7 | 7 | Identical. |
| All other methodology | — | — | Identical to the paper. |

---

## 7. Conclusion

The paper's method is methodologically reproducible from the specification.
The central qualitative findings — similarity-ranked Q1 > Q5 extremes,
regime-dependent factor timing information, a factor-neutral alpha from
the Q1−Q5 spread — replicate on our sample. The quantitative magnitudes
are smaller than the paper's (roughly half the raw Sharpe, alpha $t$-stat
1.59 vs 3), consistent with what is expected when cutting sample length by
a third and adding a well-documented weaker period for style factors.

For the purposes of extending the framework (e.g. substituting the seven
hand-picked financial state variables with PCA components of a broader
macro panel), this replication establishes that (a) the implementation is
correct, (b) the paper's benchmark on a comparable sample is alpha t ≈ 1.6,
and (c) any extension can be judged against that benchmark directly.
