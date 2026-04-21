# Regimes replication (Mulliner, Harvey, Xia, Fang, van Hemert, 2025)

## Setup

```
regimes/
├── main.R
├── R/
│   ├── utils.R
│   ├── 01_data.R         # load CSV, resample to monthly, download FF factors
│   ├── 02_similarity.R   # z-scores, Euclidean distance matrix, EWMA
│   ├── 03_strategy.R     # quintile buckets, factor-timing strategy
│   └── 04_exhibits.R     # all plots and tables
├── data/
│   └── final.csv         # put your macro CSV here
└── output/               # all exhibits land here (created on first run)
```

Packages: `tidyverse`, `lubridate`, `zoo`, `slider`, `patchwork`, `scales`.
Install once with `install.packages(c("tidyverse","zoo","slider","patchwork"))`
(lubridate and scales come with tidyverse).

Run:

```r
setwd(".../regimes")
source("main.R")
```

`main.R` sources all modules in order and writes exhibits to `output/`.

## What the scripts do

1. **`01_data.R`** — resamples `final.csv` to month-end (last daily observation
   per calendar month), filters to `sample_end` (2025-12-31), and downloads the
   six Fama-French factors (Market, Size, Value, Profitability, Investment,
   Momentum) from Ken French's website. Factor returns are converted from the
   percent units used on the website to decimals.

2. **`02_similarity.R`** — builds the transformed state series. For each
   variable at month *t*, takes the 12-month change, divides by the std of the
   last 120 monthly 12-month changes, and winsorizes at ±3. Then builds the
   full N×N Euclidean distance matrix across the seven transformed variables
   and computes the EWMA regime-shift series for the four lookbacks.

3. **`03_strategy.R`** — groups the distance-matrix rows into quantile buckets
   (quintiles by default) after masking the last 36 months before each target
   date. For each target month T and each factor, it averages the returns in
   month (t+1) across all historical months t in the bucket. The sign of that
   average is the signal applied to the factor's return in month T+1 (so the
   trade placed at T uses only info through T-1). The equally-weighted average
   across the six factors is the quintile's monthly return. Also runs the
   robustness sweeps for different quantile counts (2, 3, 4, 5, 10, 20) and
   different z-score lookbacks (1-, 3-, 5-year).

4. **`04_exhibits.R`** — builds every numbered exhibit plus the appendix
   panels. Outputs are PNGs under `output/`; stats tables are CSVs.

## Known differences from the paper

- **Sample start**: the paper uses a ~1963 start; our macro CSV begins April
  1986. The first valid z-score is around April 1997, and because the strategy
  needs a 36-month similarity mask, the factor-timing results effectively run
  from ~mid-2000 through 2025-12-31. The FF factor history extends back to
  the 1920s but the strategy back-test is bounded by when we have usable
  similarity scores.

- **12-month change**: I've used level differences (not log differences) for
  all seven variables, as this matches the paper's Exhibit 3 behaviour (in
  particular the raw z-score for oil reaching ~30 in the 1970s is only
  reproducible with level differences against very small rolling std in
  pre-OPEC-era; our sample doesn't cover that but we keep the spec). If you
  want log changes for sp500/oil/copper/vix, flip their `diff_type` to `"log"`
  in `CFG$variables` in `main.R`.

- **Momentum factor**: Ken French labels this column "Mom"; the code handles
  renaming automatically.

- **Vol-targeting for Exhibit 1**: I use a trailing 36-month realised vol to
  ex-ante scale monthly returns to 15% annualised. The paper doesn't specify
  the exact vol-targeting recipe so this is a reasonable default. Yearly
  returns are computed by summing monthly vol-targeted returns within the
  calendar year and only complete calendar years are shown.

- **Inflation bands (Exhibit 8)**: the paper's table lists nine historical
  inflation episodes going back to 1941. Only the three within our sample
  (1987-1990 Reagan boom, 2007-2008, 2022) will be visible on the plot.

## Sanity checks to run after the pipeline

The CSVs under `output/` include:

- `exhibit_04_autocorrelations.csv` — compare the 1m/3m/12m autocorrs and
  std's to the paper's Exhibit 4. Should be close if not identical.
- `exhibit_05_correlation_matrix.csv` — eyeball against Exhibit 5; yield-curve
  vs monetary-policy should be strongly negative (it's mechanical), copper-oil
  should be materially positive.
- `exhibit_10_quintile_stats.csv` — paper reports SR 0.95 for Q1 and 0.17 for
  Q5 with Q1-Q5 Sharpe of 0.82. Our shorter sample should show a similar
  monotonic pattern though exact Sharpes will differ.
- `exhibit_01_yearly_return_stats.csv` — paper reports 80% positive years,
  +13.3% conditional-on-pos, -5.1% conditional-on-neg. Should be directionally
  similar on our shorter sample.

If the FF download fails (corporate proxy, SSL, etc.), the cleanest workaround
is to download the two zip files manually from
`https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html`,
unzip them into `data/`, and point `fetch_french_csv()` at the local `.csv`.
