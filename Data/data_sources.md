# Data Inventory & Source Documentation

The project utilizes a mix of raw market data and derived indicators. All timeseries have been standardized to a daily frequency and merged into a master dataset covering the period from **January 1990 to December 2025**.

---

### 1. Market Data Overview

| Asset / Indicator | Ticker | Source | Description | Date Range |
| :--- | :--- | :--- | :--- | :--- |
| **Gold** | `XAU=` | LSEG Workspace | Gold Spot Price (USD/oz) | 1968 – 2026 |
| **Bonds** | `.TRXVUSGOV10U` | LSEG Workspace | US 10Y Treasuries Total Return Index | 1980 – 2026 |
| **Stocks** | `SPXT INDEX` | Bloomberg | S&P 500 Total Return Index | 1927 – 2026 |
| **Oil** | `DCOILWTICO` | FRED | Crude Oil Prices: West Texas Intermediate (WTI) | 1986 – 2026 |
| **USD** | `D.N.N.US` | BIS | Nominal Effective Exchange Rate (Narrow Basket) | 1983 – 2026 |
| **VIX** | `VIXCLS` | FRED | CBOE Volatility Index | 1990 – 2026 |
| **Copper** | `LMCADY` | Bloomberg | LME Copper Cash Official Price | ~1990 – 2026 |
| **US3MO** | `DGS3MO` | FRED | Market Yield on U.S. Treasury Securities at 3-Month Constant Maturity, Quoted on an Investment Basis | ~1954 – 2026 |
| **Yield Curve** | `T10Y3M` | FRED | 10-Year Minus 3-Month Treasury Constant Maturity | ~1982 – 2026 |
| **Stock-Bond Corr** | *Derived* | FRED / BBG | 3Y Rolling Correlation: 10Y Yields vs. SPX Returns | 1990 – 2025 |

---

### 2. Methodology & Derived Indicators

To ensure analytical consistency, several indicators were calculated or processed as follows:

*   **Stock-Bond Correlation:** 
    *   **Calculation:** 3-year rolling Pearson correlation.
    *   **Inputs:** Daily log-returns of the S&P 500 (`SPXT`) and the daily change ($\Delta$) in the 10-Year Treasury Yield. 
*   **Missing Values:** All series were merged using a `full_join` on the `date` index. Missing observations (holidays, weekends) were handled via **Last Observation Carried Forward (LOCF)**.
*   **Time Horizon:** While raw data ranges vary significantly, the analytical dataframes (`df_assets` and `df_macro`) are truncated to start on **January 2nd, 1990** and end on **December 31st, 2025**.

---

### 3. Folder Structure
The data is organized within the repository as follows:
*   `/Data/raw/`: Contains the original, individual `.csv` files for each ticker.
*   `/Data/`: Contains the processed, merged, and cleaned `df_assets.csv` and `df_macro.csv` files used for final analysis.

---

> **Note on Data Retrieval:** 
> Data retrieved via API and LSEG/Bloomberg terminals. The "Stock-Bond Correlation" is a custom-calculated metric and should be interpreted as a proxy for cross-asset hedging effectiveness rather than a direct market price.




### Citations: 

#### Oil
U.S. Energy Information Administration, Crude Oil Prices: West Texas Intermediate (WTI) - Cushing, Oklahoma [DCOILWTICO], retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/DCOILWTICO, May 11, 2026.

#### USD Index
Bank for International Settlements (2026), https://data.bis.org/topics/EER/BIS,WS_EER,1.0/D.N.N.US (accessed on 11 May 2026).

#### VIX
Chicago Board Options Exchange, CBOE Volatility Index: VIX [VIXCLS], retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/VIXCLS, May 11, 2026.

#### Yield Spread
Federal Reserve Bank of St. Louis, 10-Year Treasury Constant Maturity Minus 3-Month Treasury Constant Maturity [T10Y3M], retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/T10Y3M, May 11, 2026.