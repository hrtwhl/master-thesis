"""
robustness_ksweep.py  (Robustness R1)
=====================================
Tests whether the 2-effective-regime collapse of the macro layer is
DATA-DRIVEN or merely IMPOSED by the tight, penalised default macro
config (Kmax=5, lambda_K=1, spawn=2.5, monotone-K).

For each config in config.R1_KSWEEP_CONFIGS we run ONLY the macro
Wasserstein-HMM (no market layer, no MVO) over the full OOS and count
how many DURABLE templates emerge -- templates that are the dominant
macro regime on at least R1_DURABLE_SHARE of OOS days.

If the number of durable regimes stays ~2 even when K is allowed up to
8 with no complexity penalty and a low spawn threshold, the collapse is
a genuine property of daily macro data, not a config artifact.

Outputs:
    output/tables/R1_ksweep_summary.csv         (one row per config)
    output/tables/R1_ksweep_<label>_shares.csv  (regime dominant-shares)
    output/raw/R1_ksweep_<label>_path.csv       (daily regime path)
    output/charts/R1_ksweep_durable_regimes.png (bar chart)
"""
import os, sys, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config
from data import load_asset_returns, load_macro_panel, align_calendars
from hierarchical_hmm import run_macro_regime_path

os.makedirs(config.RAW_DIR, exist_ok=True)
os.makedirs(config.TABLE_DIR, exist_ok=True)
os.makedirs(config.CHART_DIR, exist_ok=True)


def _cfg_to_dict(row):
    label, kmin, kmax, gmax, lamk, spawn, mono = row
    return label, dict(
        min_regimes=kmin, max_regimes=kmax, g_max=gmax,
        lam_k=lamk, spawn_thresh=spawn, monotone_K=mono,
        eta_tpl=config.MACRO_ETA_TPL, f_k=config.MACRO_F_K,
    )


def main():
    t0 = time.time()
    print("[R1] loading data...", flush=True)
    ar = load_asset_returns()
    mp = load_macro_panel()
    ar, mp = align_calendars(ar, mp)

    summary_rows = []
    durable_counts = {}

    for row in config.R1_KSWEEP_CONFIGS:
        label, cfg = _cfg_to_dict(row)
        print(f"\n[R1] config '{label}': {cfg}", flush=True)
        t1 = time.time()
        path = run_macro_regime_path(
            ar, mp, macro_cfg=cfg, oos_start=config.OOS_START,
            track_space="outcome", verbose=True,
        )
        path.to_csv(os.path.join(config.RAW_DIR, f"R1_ksweep_{label}_path.csv"))

        lbl = path["regime_macro"].dropna().astype(int)
        shares = lbl.value_counts(normalize=True).sort_index()
        shares.to_csv(os.path.join(config.TABLE_DIR,
                                   f"R1_ksweep_{label}_shares.csv"))

        durable = int((shares >= config.R1_DURABLE_SHARE).sum())
        durable_counts[label] = durable

        summary_rows.append({
            "config":            label,
            "K_max_allowed":     cfg["max_regimes"],
            "lambda_K":          cfg["lam_k"],
            "spawn_thresh":      cfg["spawn_thresh"],
            "monotone_K":        cfg["monotone_K"],
            "K_reached":         int(path["K_macro"].dropna().max()),
            "G_reached":         int(path["G_macro"].dropna().max()),
            "n_durable_regimes": durable,
            "top1_share":        float(shares.iloc[shares.values.argmax()]) if len(shares) else np.nan,
            "n_distinct_labels": int(lbl.nunique()),
            "elapsed_s":         round(time.time() - t1, 1),
        })
        print(f"[R1] '{label}': K_reached={summary_rows[-1]['K_reached']}, "
              f"durable_regimes={durable} "
              f"(shares>{config.R1_DURABLE_SHARE:.0%}), "
              f"{summary_rows[-1]['elapsed_s']}s", flush=True)

    summary = pd.DataFrame(summary_rows).set_index("config")
    summary.to_csv(os.path.join(config.TABLE_DIR, "R1_ksweep_summary.csv"))
    print("\n=========== R1 K-SWEEP SUMMARY ===========")
    print(summary.to_string(), flush=True)

    # Bar chart of durable regimes per config
    fig, ax = plt.subplots(figsize=(10, 4))
    labels = list(durable_counts.keys())
    vals = [durable_counts[k] for k in labels]
    ax.bar(range(len(labels)), vals, color="steelblue")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel(f"# durable macro regimes (>{config.R1_DURABLE_SHARE:.0%} of days)")
    ax.set_title("R1: Durable macro regimes vs model freedom")
    ax.grid(True, axis="y")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.05, str(v), ha="center", fontsize=10)
    fig.savefig(os.path.join(config.CHART_DIR, "R1_ksweep_durable_regimes.png"),
                bbox_inches="tight")
    plt.close(fig)

    print(f"\n[R1] done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
