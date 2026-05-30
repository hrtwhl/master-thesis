# Methodology

Detailed mathematical specification of the regime-aware portfolio
construction infrastructure.  Two strategies are defined here:

1. **Pure-Market WHMM** — faithful reproduction of Boukardagha (2026).
2. **Hierarchical WHMM** — our extension that adds a macro layer
   on top of Boukardagha's market layer, following the
   Hierarchical-HMM architecture of Fine, Singer & Tishby (1998).

Both strategies feed their template-mixture conditional moments into
the same transaction-cost-aware mean-variance optimiser.

---

## 1. Notation and Universe

Let $t$ index daily trading dates and $N = 5$ the number of risky
assets in the cross-asset universe:

$$
\mathcal{A} = \{ \text{SPX}, \text{BOND}, \text{GOLD}, \text{OIL}, \text{USD} \}.
$$

Adjusted-close prices $P_t \in \mathbb{R}^N$ are converted to daily
log returns

$$
r_t = \log P_t - \log P_{t-1} \in \mathbb{R}^N.
$$

For the hierarchical strategy we additionally use a macro panel of
$M = 7$ state variables (Mulliner et al.):

$$
\mathcal{M} = \{\text{stocks, oil, copper, yield\_curve, stock\_bond\_corr, vix, us3mo}\}.
$$

The out-of-sample window is $t \in [t_0, T]$ with
$t_0 = $ 2005-01-03 and $T = $ 2025-12-31.  All quantities indexed by
$t$ use only data observed up to $t-1$ (strict causality).

---

## 2. Feature Construction

### 2.1 Asset features (Boukardagha §4)

For each asset $i \in \mathcal{A}$ define the 60-day rolling
volatility and 20-day rolling mean of returns:

$$
\sigma_{t,i} = \sqrt{\frac{1}{59}\sum_{s=t-59}^{t}(r_{s,i}-\bar r_{t,i})^2},
\qquad
m_{t,i} = \frac{1}{20}\sum_{s=t-19}^{t} r_{s,i}.
$$

The asset feature vector at time $t$ is the stacked

$$
\boxed{\; x_t = [\,r_t;\; \sigma_t;\; m_t\,] \in \mathbb{R}^{3N=15}. \;}
$$

### 2.2 Macro features (hierarchical strategy only)

Each macro variable $j \in \mathcal{M}$ is classified as
**price-like** (strictly positive) or **level-like**.  Stationary
increments are computed as

$$
\Delta_{t,j} =
\begin{cases}
\log m_{t,j} - \log m_{t-1,j} & j \text{ price-like (stocks, copper, vix)}\\
m_{t,j} - m_{t-1,j}            & j \text{ level-like (oil, yield\_curve, stock\_bond\_corr, us3mo)}
\end{cases}
$$

(Oil is treated as level-like because WTI spot printed negative in
April 2020.)  Rolling 60-day std $\sigma^{(M)}_{t,j}$ and 20-day mean
$m^{(M)}_{t,j}$ of $\Delta_{\cdot,j}$ are stacked into the
21-dimensional macro feature vector:

$$
\boxed{\; x^{(M)}_t = [\,\Delta_t;\; \sigma^{(M)}_t;\; m^{(M)}_t\,] \in \mathbb{R}^{3M=21}. \;}
$$

Because raw macro features span four orders of magnitude (oil
increments at $\mathcal{O}(1)$ vs stock-bond correlation increments
at $\mathcal{O}(10^{-3})$), the macro feature vector is **causally
z-scored** at each step $t$ using sample mean and standard deviation
computed on the in-sample training slice $\{x^{(M)}_s\}_{s < t}$.

---

## 3. Wasserstein Hidden Markov Model (Pure-Market Strategy)

### 3.1 Gaussian HMM with predictive model-order selection

At each trading day $t$, an expanding-window training set
$\mathcal{H}_t = \{ x_s : s < t \}$ is built.  A Gaussian HMM with
$K$ latent states is fitted to $\mathcal{H}_t$:

$$
z_t \in \{1, \ldots, K\}, \qquad
x_t \mid z_t = k \sim \mathcal{N}(\mu_{t,k}, \Sigma_{t,k}),
$$

with state-transition matrix $A_t \in [0,1]^{K \times K}$.  Parameter
estimation is by Baum-Welch (EM) at `tol=1e-3` for at most 100
iterations.

#### Predictive K-selection (every $F_K = 5$ days)

The validation slice is the most recent $L_\text{val} = 252$ rows.
For each candidate $K \in \{K_{\min}, \ldots, K_{\max}\} = \{5, 6\}$:
fit a $K$-state HMM on $\mathcal{H}_t \setminus \mathcal{V}_t$,
score on $\mathcal{V}_t$, and select

$$
\boxed{\; K_t \;=\; \arg\max_{K \in \{K_\text{curr},\, \ldots, K_{\max}\}} \;
\sum_{s \in \mathcal{V}_t} \log p(x_s \mid \mathcal{F}_{s-1}, K) \;-\;\lambda_K \cdot K, \;}
$$

with the monotone-K constraint $K_t \geq K_{\text{curr}}$ (Boukardagha
2026 §7.1.14).  Complexity penalty $\lambda_K = 1$.

#### Filtered posteriors

Given $K_t$ and the fitted HMM, smoothed component posteriors over
the full history are obtained via the forward-backward algorithm:

$$
p_{t,k} = P(z_t = k \mid \mathcal{F}_{t-1}), \qquad k = 1, \ldots, K_t.
$$

The most likely component at $t$ is $\hat z_t = \arg\max_k p_{t,k}$.

### 3.2 Conditional moments from forward returns

To preserve strict causality, regime-conditional moments are computed
on **forward** returns $r_{s+1}$ grouped by $\hat z_s$:

$$
\mu_{t,k} = \frac{1}{N_{t,k}}\sum_{s : \hat z_s = k, s < t-1} r_{s+1},
\qquad
\Sigma_{t,k} = \text{LedoitWolf}\!\left( \{ r_{s+1} : \hat z_s = k \} \right),
$$

where $N_{t,k}$ is the number of obs in component $k$ and
LedoitWolf$(\cdot)$ is the Ledoit-Wolf shrinkage covariance estimator.
Components with $N_{t,k} < 60$ obs default to $\mu = 0$ and
$\Sigma = I_N$.

### 3.3 Template-based regime identity (Wasserstein tracking)

Repeated HMM re-fits can permute the latent component labels.
Boukardagha replaces discrete assignment matching with persistent
**templates** $\{\Theta_g = (\mu_g, \Sigma_g)\}_{g=1}^G$ tracked over
time using the closed-form 2-Wasserstein distance between Gaussian
distributions:

$$
\boxed{\;
W_2^2\big(\mathcal{N}(\mu_1, \Sigma_1), \mathcal{N}(\mu_2, \Sigma_2)\big)
= \lVert \mu_1 - \mu_2 \rVert_2^2
+ \operatorname{Tr}\!\left( \Sigma_1 + \Sigma_2 - 2\,\big(\Sigma_2^{1/2}\Sigma_1\Sigma_2^{1/2}\big)^{1/2} \right). \;}
$$

#### Mapping

Each HMM component $k$ at time $t$ is mapped to its nearest template:

$$
g(k) = \arg\min_{g \leq G} W_2\big(\mathcal{N}(\mu_g, \Sigma_g),\, \mathcal{N}(\mu_{t,k}, \Sigma_{t,k})\big).
$$

Aggregated template posteriors:

$$
p_{t,g} = \sum_{k:\, g(k)=g} p_{t,k}.
$$

#### Template spawning

If some component $k^\star$ is too far from every existing template
($\min_g W_2 > \tau_\text{spawn} = 2.5$) and $G < G_\text{max} = 8$,
a new template is spawned from the unmatched component with the
highest posterior:

$$
k^\star = \arg\max_{\{k : \min_g W_2(\Theta_g, \cdot_k) > \tau_\text{spawn}\}} p_{t,k},
\qquad
\Theta_{G+1} \leftarrow (\mu_{t,k^\star}, \Sigma_{t,k^\star}).
$$

#### Template EMA update

Templates are exponentially smoothed using posterior-weighted
component averages.  For each template $g$ with components mapped
to it:

$$
\bar\mu_{t,g} = \frac{\sum_{k:\, g(k)=g} p_{t,k}\, \mu_{t,k}}{\sum_{k:\, g(k)=g} p_{t,k}},
\qquad
\bar\Sigma_{t,g} = \frac{\sum_{k:\, g(k)=g} p_{t,k}\, \Sigma_{t,k}}{\sum_{k:\, g(k)=g} p_{t,k}},
$$

$$
\boxed{\;
\mu_g \leftarrow (1-\eta)\mu_g + \eta\, \bar\mu_{t,g},
\qquad
\Sigma_g \leftarrow (1-\eta)\Sigma_g + \eta\, \bar\Sigma_{t,g}. \;}
$$

Learning rate $\eta = 0.05$.

### 3.4 Template-mixture moments

The conditional moments fed to the optimizer are

$$
\boxed{\;
\mu_t = \sum_{g=1}^G p_{t,g}\, \mu_g,
\qquad
\Sigma_t = \sum_{g=1}^G p_{t,g}\, \Sigma_g. \;}
$$

### 3.5 Transaction-cost-aware mean-variance optimization

Daily portfolio weights $w_t \in \mathbb{R}^N$ solve

$$
\boxed{\;
\max_{w_t}\;\; \mu_t^\top w_t \;-\; \gamma\, w_t^\top \Sigma_t w_t \;-\; \tau\, \lVert w_t - w_{t-1} \rVert_1
\quad \text{s.t.} \quad
\mathbf{1}^\top w_t = 1,\;\; w_t \geq 0,\;\; \lVert w_t \rVert_\infty \leq w_\text{max}. \;}
$$

Parameter values (Boukardagha verbatim):
$\gamma = 5$, $\tau = 0.0002$, $w_\text{max} = 0.6$.  Solved by
OSQP via cvxpy with SCS fallback.

---

## 4. Hierarchical Wasserstein HMM (two strategies)

The hierarchical strategies add a top-level macro Wasserstein HMM
above the market WHMM.  Following Fine, Singer & Tishby (1998), the
two levels are organised as a tree.  We define **two** hierarchical
strategies that differ in how the macro layer is used:

* **Hierarchical B** (§4.1-4.5) — *joint mixture*.  Macro
  templates are tracked in macro-feature space; conditional moments
  are computed on the joint $(g, h)$ cells and mixed by
  $p(t,g,h) = p^{(M)}_{t,g}\, p^{(A)}_{t,h}$.  No expected-return
  tilt.

* **Hierarchical C** (§4.6) — *macro as risk modulator*.  Macro
  templates are tracked in **asset-outcome space**; the macro
  posterior is tempered; the market layer alone produces $\mu_t$ and
  $\Sigma_t$ (exactly as in the pure-market strategy); the macro
  layer produces a scalar **stress** score that scales the effective
  risk aversion $\gamma_t$ and inflates $\Sigma_t$.

```
Hierarchical B:
Root
├── Macro template g = 1, …, G_macro     (macro-feature-space templates)
└── Market template h = 1, …, G_market   (asset features)
        ↓
    Joint cell (g, h) provides (μ_{g,h}, Σ_{g,h}); mix by p(g,h)

Hierarchical C:
Root
├── Macro template g = 1, …, G_macro     (asset-OUTCOME-space templates)
│        ↓ tempered posterior → stress_t ∈ [0,1]
│        ↓ γ_t = γ·(1+κ·stress_t),  Σ_t ← Σ_t·(1+s·stress_t)
└── Market template h = 1, …, G_market   (asset features) → μ_t, Σ_t
```

### 4.1 Macro Wasserstein HMM (top layer)

Identical machinery to Section 3 applied to the macro feature panel
$X^{(M)} = \{x^{(M)}_s\}_{s < t}$ (causally z-scored).  Parameters:
$K_{\min}^{(M)} = 4$, $K_{\max}^{(M)} = 5$, $G_{\max}^{(M)} = 6$,
$\eta^{(M)} = 0.05$, $\tau^{(M)}_\text{spawn} = 2.5$.  K-selection
cadence $F_K^{(M)} = 21$ days, refit cadence $F_\text{refit}^{(M)} = 5$
days (macro regimes are slow).

One key difference: macro template moments are tracked in **macro
feature space** (the HMM's own emission means and covariances),
because the macro layer's role is to identify slow-moving environment
regimes, not to forecast asset returns directly.

Output at time $t$: macro template posterior $p^{(M)}_{t,g}$ for
$g = 1, \ldots, G^{(M)}_t$, and per-historical-row macro template
label $g_s$.

### 4.2 Market Wasserstein HMM (bottom layer, shared)

Identical to Section 3.  Output at time $t$: market template
posterior $p^{(A)}_{t,h}$ for $h = 1, \ldots, G^{(A)}_t$, and
per-historical-row market template label $h_s$.

### 4.3 Joint conditional moments

For each joint cell $(g, h)$:

$$
\mathcal{S}_{g,h} \;=\; \{ s : g_s = g,\; h_s = h,\; s < t-1 \},
\qquad
N_{g,h} = |\mathcal{S}_{g,h}|.
$$

If $N_{g,h} \geq N_\text{min} = 60$:

$$
\mu_{g,h} = \frac{1}{N_{g,h}}\sum_{s \in \mathcal{S}_{g,h}} r_{s+1},
\qquad
\Sigma_{g,h} = \text{LedoitWolf}\!\left(\{r_{s+1} : s \in \mathcal{S}_{g,h}\}\right).
$$

Otherwise fall back to the market-only marginal:

$$
\mathcal{S}_h \;=\; \{ s : h_s = h \},
\qquad
\mu_{g,h} = \mu_h,\;\; \Sigma_{g,h} = \Sigma_h.
$$

### 4.4 Composite probability and mixture moments

Under conditional independence of the macro and market layers given
their joint history:

$$
\boxed{\;
p_t(g, h) \;=\; p^{(M)}_{t,g} \cdot p^{(A)}_{t,h}. \;}
$$

Template-mixture moments:

$$
\boxed{\;
\mu_t = \sum_{g,h} p_t(g, h)\, \mu_{g,h},
\qquad
\Sigma_t = \sum_{g,h} p_t(g, h)\, \Sigma_{g,h}. \;}
$$

### 4.5 Optimization (Hierarchical B)

Hierarchical B passes the joint-mixture moments $(\mu_t, \Sigma_t)$
of §4.4 directly into the same optimiser as the pure-market strategy
(§3.5), with **no expected-return tilt**:

$$
\max_{w_t}\;\; \mu_t^\top w_t - \gamma\, w_t^\top \Sigma_t w_t - \tau \lVert w_t - w_{t-1} \rVert_1,
$$

subject to the same budget, long-only, and box constraints.

> **Design history.**  An earlier version of Hierarchical B applied a
> macro-conditional tilt to the expected returns of risky assets,
> $\widetilde\mu_{t,i} = \phi_{t,i}\,\mu_{t,i}$ with
> $\phi_{t,i} = 1 + \alpha\sum_g p^{(M)}_{t,g}\tanh(S_{g,i})$.  This
> *symmetric* tilt only ever pushed $\mu_\text{SPX}$ **upward** (SPX
> Sharpe is positive in both dominant macro regimes), worsening the
> 2022 drawdown.  An *asymmetric* variant centred on the macro-mixed
> mean Sharpe, $\tanh(S_{g,i}-\bar S_{t,i})$, turned out to be a
> numerical no-op because the macro posterior is one-hot $\approx99\%$
> of days (so $\bar S_{t,i}=S_{g,i}$ and the tilt vanishes).  Both are
> documented in `analysis_results.md`; the tilt has been **removed**,
> and the lessons motivate Hierarchical C below.

### 4.6 Hierarchical C — macro layer as risk modulator

Hierarchical C addresses the three failure modes diagnosed for the
joint-mixture design (`analysis_results.md`): (1) macro regimes that
are not discriminative in asset-return space, (2) dilution of the
market layer's clean directional signal by joint-cell averaging, and
(3) an over-confident, one-hot macro posterior.  The 21-feature
Mulliner macro set is retained unchanged.

**(C1) Macro templates tracked in asset-outcome space.**
The macro WHMM still infers latent states $z^{(M)}_t$ from the
z-scored 21-d macro features, but the persistent **templates** are
matched and EMA-updated (via the same $W_2$ machinery, §3.3) using
the Gaussian of **forward asset returns** conditional on each macro
component:

$$
\mu^{(C)}_{g} = \mathbb{E}[\,r_{s+1} \mid z^{(M)}_s = g\,], \qquad
\Sigma^{(C)}_{g} = \text{LedoitWolf}\big(\{ r_{s+1} : z^{(M)}_s = g \}\big).
$$

Macro regimes are thereby discriminative for allocation **by
construction** — two macro states that imply the same asset-return
distribution collapse to the same template.

**(C2) Macro modulates risk, not direction.**
The market WHMM (§3.1-3.4) alone produces the expected-return vector
$\mu_t$ and covariance $\Sigma_t$ — identical to the pure-market
strategy, preserving its clean directional signal.  The macro layer
contributes only a scalar **stress** score $\text{stress}_t \in [0,1]$.

For each macro template $g$ with $\geq N_\text{min}$ observations,
let $r^\text{EW}_{s+1} = \tfrac1N\mathbf 1^\top r_{s+1}$ be the
equal-weight forward return (a regime-agnostic turbulence proxy).
Define a raw turbulence score (default `HIER_C_STRESS_METRIC='vol'`):

$$
\rho_g = \sqrt{252}\;\operatorname{std}\big(\{ r^\text{EW}_{s+1} : z^{(M)}_s = g \}\big),
$$

(or, optionally, the within-regime max drawdown, or
$\max(0,-\text{Sharpe})$).  Cross-sectionally min-max normalise to
$[0,1]$:

$$
\widehat\rho_g = \frac{\rho_g - \min_{g'}\rho_{g'}}{\max_{g'}\rho_{g'} - \min_{g'}\rho_{g'}}.
$$

**(C3) Tempered, prior-blended macro posterior.**
Because the raw macro posterior $p^{(M)}_{t,\cdot}$ is one-hot
$\approx 99\%$ of the time, we soften it with temperature
$T = $ `HIER_C_MACRO_TEMPERATURE` and blend with a uniform prior of
weight $\beta = $ `HIER_C_PRIOR_BLEND`:

$$
\tilde p^{(M)}_{t,g} \propto \big(p^{(M)}_{t,g}\big)^{1/T},
\qquad
\bar p^{(M)}_{t,g} = (1-\beta)\,\tilde p^{(M)}_{t,g} + \beta \cdot \tfrac1{G^{(M)}}.
$$

The composite stress score is the tempered-posterior expectation of
the normalised turbulence:

$$
\boxed{\; \text{stress}_t = \sum_{g} \bar p^{(M)}_{t,g}\,\widehat\rho_g \;\in\;[0,1]. \;}
$$

**Risk modulation and optimisation.**
The macro stress raises the effective risk aversion and inflates the
covariance:

$$
\boxed{\;
\gamma_t = \gamma\,\big(1 + \kappa\,\text{stress}_t\big),
\qquad
\Sigma^{(C)}_t = \Sigma_t\,\big(1 + s\,\text{stress}_t\big), \;}
$$

with $\kappa = $ `HIER_C_KAPPA_GAMMA` $= 4$ and $s = $
`HIER_C_SIGMA_SCALE` $= 1$.  The day-$t$ weights solve

$$
\max_{w_t}\;\; \mu_t^\top w_t - \gamma_t\, w_t^\top \Sigma^{(C)}_t w_t - \tau \lVert w_t - w_{t-1} \rVert_1,
$$

subject to the same constraints as §3.5.  In benign macro regimes
($\text{stress}_t \to 0$) Hierarchical C reduces exactly to the
pure-market strategy; in turbulent macro regimes it de-risks
(higher $\gamma_t$, inflated $\Sigma_t$) toward a flatter, more
diversified allocation, while letting the market layer continue to
choose the direction of any active tilt.

---

## 5. Benchmarks

* **Equal-Weight**: static $w_i = 1/N$ for all assets, frictionless,
  rebalanced daily.
* **60/40**: static $w_\text{SPX} = 0.6$, $w_\text{BOND} = 0.4$,
  frictionless, rebalanced daily.

---

## 6. Implementation notes

### 6.1 Strict causality

Every quantity indexed by $t$ uses **only** observations available
strictly before $t$:

$$
\mathcal{H}_t \;=\; \{ x_s : s < t \}.
$$

Forward-return moments in Sections 3.2 and 4.3 use $r_{s+1}$ for
$s < t-1$ to avoid look-ahead.  The portfolio realised on day $t$
uses weights $w_t$ chosen from $\mathcal{H}_t$ multiplied by the
actual day-$t$ return $r_t$.

### 6.2 Warm-started EM

Daily refits over a 21-year OOS with full-covariance K=6 HMMs on
8000-row training sets would naively take 5+ hours of compute.  We
exploit EM's monotonicity by **warm-starting** every daily refit
from the previous day's fitted parameters via hmmlearn's
`init_params=""` mechanism:

$$
\theta_t^{(0)} = \theta_{t-1}^{(\text{converged})}.
$$

Adding a single row to the training history typically lands EM in
the same likelihood-surface basin and reconverges in 3-10 iterations
(~0.1-0.3 s).  This is **not a methodological change**: the converged
maximum on each day is identical to what a cold-started fit would
produce; only the optimisation path differs.

### 6.3 Deviations from Boukardagha (2026)

The only methodological deviation is `HMM_N_ITER = 100` (paper: 300).
EM converges in $\leq 80$ iterations at `tol=1e-3` for every fit
size encountered in this backtest, so the cap is non-binding.  Set
`HMM_N_ITER = 300` in `config.py` to recover the paper exactly.

---

## 6R. Robustness analyses

### 6R.1 Loosened-K macro sweep (R1)

The macro layer collapses to $\approx 2$ effective regimes under the
default configuration ($K^{(M)}\in\{4,5\}$, $\lambda_K=1$,
$\tau_\text{spawn}=2.5$, monotone-$K$).  To establish that this is a
property of the *data* rather than an artifact of a tight, penalised
configuration, we re-run the macro Wasserstein-HMM alone (no market
layer, no MVO) under a ladder of progressively freer configurations:

| Config | $K_{\max}$ | $\lambda_K$ | $\tau_\text{spawn}$ | monotone-$K$ |
| ------ | ---------- | ----------- | ------------------- | ------------ |
| `default_4_5`          | 5 | 1.0 | 2.5 | yes |
| `loose_4_8`            | 8 | 1.0 | 2.5 | yes |
| `loose_4_8_nopenalty`  | 8 | 0.0 | 2.5 | yes |
| `loose_4_8_lowspawn`   | 8 | 0.0 | 1.0 | yes |
| `free_3_8_nonmono`     | 8 | 0.0 | 1.0 | no  |

For each config we count the number of **durable** templates — those
that are the dominant macro regime on at least $\rho_\text{min}$
(`R1_DURABLE_SHARE` $= 2\%$) of OOS days.  A template that the HMM
nominally creates but never sustains as a dominant state is not a
genuine regime.  Macro templates are tracked in asset-outcome space
(as in Hierarchical C).  If the number of durable regimes stays
$\approx 2$ even when $K$ may rise to 8 with no penalty and a low
spawn threshold, the collapse is data-driven.

### 6R.2 No-overlap macro panel (R2)

Two of the seven Mulliner macro variables — `stocks` ($\equiv$ SPX)
and `oil` — are *also* tradeable assets in the universe
$\mathcal{A}$.  Their inclusion in the macro panel makes the macro and
market layers informationally dependent, and weakens the conditional
independence assumption $p(t,g,h) = p^{(M)}_{t,g} p^{(A)}_{t,h}$ used
in Hierarchical B.  R2 re-runs both hierarchical strategies with the
non-overlapping macro panel

$$
\mathcal{M}_\text{R2} = \{\text{copper, yield\_curve, stock\_bond\_corr, vix, us3mo}\},
$$

containing only genuinely exogenous macro-financial series.  If
Hierarchical C still outperforms the pure-market strategy under
$\mathcal{M}_\text{R2}$, its edge cannot be attributed to the macro
layer re-using tradeable-asset information that the market layer
already observes.

---

## 7. Strategy parameters at a glance

| Parameter | Pure-Market | Hierarchical (top) | Hierarchical (bottom) |
| --------- | ----------- | ------------------- | --------------------- |
| Feature dim | 15 (asset) | 21 (macro, z-scored) | 15 (asset) |
| $K_{\min}, K_{\max}$ | 5, 6 | 4, 5 | 5, 6 |
| K-selection cadence | 5 days | 21 days | 5 days |
| Refit cadence | 1 day (warm) | 5 days | 1 day (warm) |
| $L_\text{val}$ | 252 | 252 | 252 |
| $\lambda_K$ | 1.0 | 1.0 | 1.0 |
| $G_\text{max}$ | 8 | 6 | 8 |
| $\eta$ (EMA) | 0.05 | 0.05 | 0.05 |
| $\tau_\text{spawn}$ | 2.5 | 2.5 | 2.5 |
| Monotone-K | yes | yes | yes |

| MVO parameter | Value |
| ------------- | ----- |
| $\gamma$ (risk aversion) | 5 |
| $\tau$ (L1 turnover penalty) | 0.0002 |
| $w_\text{max}$ (per-asset cap) | 0.6 |
| $\sum w_i$ (budget) | 1 |
| $w_i$ (long-only) | $\geq 0$ |
| Solver | OSQP, SCS fallback |

| Hierarchical C parameter | Value |
| ------------------------ | ----- |
| $\kappa$ (`HIER_C_KAPPA_GAMMA`) — risk-aversion sensitivity | 4.0 |
| $s$ (`HIER_C_SIGMA_SCALE`) — covariance inflation | 1.0 |
| $T$ (`HIER_C_MACRO_TEMPERATURE`) — posterior temper | 4.0 |
| $\beta$ (`HIER_C_PRIOR_BLEND`) — uniform-prior blend | 0.10 |
| `HIER_C_STRESS_METRIC` | `'vol'` (alt: `'drawdown'`, `'sharpe'`) |

(Hierarchical B has no extra parameters beyond the macro-layer
settings in §4.1; the expected-return tilt has been removed.)

---

## References

[1] Boukardagha, A. (2026). *Explainable Regime-Aware Investing*. arXiv:2603.04441.
[2] Fine, S., Singer, Y., Tishby, N. (1998). *The Hierarchical Hidden Markov Model: Analysis and Applications*. Machine Learning, 32, 41-62.
[3] Ledoit, O., Wolf, M. (2004). *A well-conditioned estimator for large-dimensional covariance matrices*. J. Multivariate Analysis, 88(2), 365-411.
[4] Takatsu, A. (2011). *Wasserstein geometry of Gaussian measures*. Osaka Journal of Mathematics, 48, 1005-1026.
[5] Hamilton, J. D. (1989). *A new approach to the economic analysis of nonstationary time series and the business cycle*. Econometrica, 57(2), 357-384.
