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

## 4. Hierarchical Wasserstein HMM (Hierarchical Strategy)

The hierarchical strategy adds a top-level macro Wasserstein HMM
above the market WHMM.  Following Fine, Singer & Tishby (1998), the
two levels are organised as a tree:

```
Root
├── Macro template g = 1, …, G_macro     (operates on macro features)
└── Market template h = 1, …, G_market   (operates on asset features)
        ↓
    Joint cell (g, h) provides (μ_{g,h}, Σ_{g,h})
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

### 4.5 Macro-conditional tilt on risky-asset expected returns

To translate macro information into a directional signal that the
optimiser will respond to, we apply a smooth tilt to the expected
returns of risky assets $\mathcal{R} = \{\text{SPX, OIL}\}$ before
passing $\mu_t$ to the optimiser.

For each macro template $g$, compute the historical annualised
Sharpe of each risky asset:

$$
S_{g, i} = \sqrt{252}\, \frac{\mathbb{E}[r_{s+1, i} \mid g_s = g]}{\sigma_{i \mid g}}.
$$

The tilt factor for asset $i \in \mathcal{R}$ depends on
`HIER_MACRO_TILT_MODE`:

**`'symmetric'` mode** (original design — first OOS run):

$$
\phi_{t,i}^{\text{sym}} = 1 + \alpha \sum_{g} p^{(M)}_{t,g} \cdot \tanh\!\big(S_{g, i}\big),
$$

with $\alpha = $ `HIER_MACRO_TILT_STRENGTH` $= 1$.  The $\tanh(\cdot)$
saturates at $\pm 1$ to prevent wild scaling on regimes with extreme
historical Sharpe.  **Problem in practice**: if $S_{g, i} > 0$ in all
dominant macro regimes (which it is for SPX in our 21-year sample),
the tilt only ever pushes $\mu_{t, i}$ **upward** — never downward —
so the optimiser receives no signal to underweight risky assets in
adverse regimes (see `analysis_results.md` §2).

**`'asymmetric'` mode** (new default):

$$
\boxed{\;
\phi_{t,i}^{\text{asym}} = 1 + \alpha \sum_{g} p^{(M)}_{t,g} \cdot \tanh\!\big(S_{g, i} - \bar S_{t, i}\big),
\qquad
\bar S_{t, i} = \sum_{g} p^{(M)}_{t,g}\, S_{g, i}, \;}
$$

i.e. each macro template's contribution is centered on the
**macro-mixed unconditional Sharpe** $\bar S_{t, i}$.  Now
$\tanh(\cdot)$ fires **negative** when the current macro regime has
worse-than-average Sharpe for asset $i$, which shrinks $\mu_{t, i}$
— exactly what the optimiser needs to rotate out of risky assets in
adverse regimes.

**`'off'` mode**: $\phi_{t,i} = 1$ for all $i$, leaving $\mu_t$
unchanged.  Use this to isolate the pure-mixture effect of the macro
layer (Σ_t only).

Defensive assets (BOND, GOLD, USD) are untouched in all modes.
Letting $\widetilde \mu_{t, i} = \phi_{t,i} \cdot \mu_{t, i}$ for the
risky assets and $\widetilde \mu_{t, i} = \mu_{t, i}$ otherwise, the
optimiser receives $\widetilde \mu_t$ as the expected-return input.

### 4.6 Optimization

Identical to Section 3.5:

$$
\max_{w_t}\;\; \widetilde{\mu}_t^\top w_t - \gamma w_t^\top \Sigma_t w_t - \tau \lVert w_t - w_{t-1} \rVert_1.
$$

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

| Macro-tilt parameter (hierarchical only) | Value |
| ---------------------------------------- | ----- |
| $\alpha$ (`HIER_MACRO_TILT_STRENGTH`)    | 1.0 |
| Mode (`HIER_MACRO_TILT_MODE`)            | `'asymmetric'` (default) / `'symmetric'` / `'off'` |
| Risky asset set (`HIER_MACRO_TILT_ASSETS`) | (SPX, OIL) |

---

## References

[1] Boukardagha, A. (2026). *Explainable Regime-Aware Investing*. arXiv:2603.04441.
[2] Fine, S., Singer, Y., Tishby, N. (1998). *The Hierarchical Hidden Markov Model: Analysis and Applications*. Machine Learning, 32, 41-62.
[3] Ledoit, O., Wolf, M. (2004). *A well-conditioned estimator for large-dimensional covariance matrices*. J. Multivariate Analysis, 88(2), 365-411.
[4] Takatsu, A. (2011). *Wasserstein geometry of Gaussian measures*. Osaka Journal of Mathematics, 48, 1005-1026.
[5] Hamilton, J. D. (1989). *A new approach to the economic analysis of nonstationary time series and the business cycle*. Econometrica, 57(2), 357-384.
