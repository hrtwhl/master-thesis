"""
main.py
=======
End-to-end driver:
    1) run pure-market Wasserstein-HMM backtest        (run_pure.py)
    2) run BOTH hierarchical backtests (B and C)        (run_hier.py)
    3) aggregate -> tables + charts + benchmarks        (aggregate.py)
    4) build narrative summary                          (narrative.py)
    5) build HTML report                                (report.py)
    6) R1: loosened-K macro sweep robustness            (robustness_ksweep.py)
    7) R2: no-overlap macro panel robustness            (robustness_nooverlap.py)

Pass --no-robustness to run only stages 1-5.

Strategies compared:
    PureMarket_WHMM  - Boukardagha (2026) replication
    Hierarchical_B   - macro x market joint mixture, no tilt
    Hierarchical_C   - macro as risk modulator (Fix C)
    EqualWeight, SixtyForty - passive benchmarks

Robustness:
    R1 - does the 2-regime macro collapse survive a much looser model
         configuration (K up to 8, no penalty, low spawn threshold)?
    R2 - does Hierarchical C still beat Pure when the macro panel
         excludes the tradeable assets (stocks, oil)?

Each stage is also a standalone script:

    python run_pure.py
    python run_hier.py
    python aggregate.py
    python narrative.py
    python report.py
    python robustness_ksweep.py
    python robustness_nooverlap.py
"""
import os, sys, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")


def main(run_robustness: bool = True):
    t0 = time.time()
    n = 7 if run_robustness else 5

    print(f"===========  STAGE 1 / {n} : PURE MARKET  ===========", flush=True)
    import run_pure
    run_pure.main()

    print(f"\n===========  STAGE 2 / {n} : HIERARCHICAL (B & C)  ===========", flush=True)
    import run_hier
    run_hier.main()

    print(f"\n===========  STAGE 3 / {n} : AGGREGATE  ===========", flush=True)
    import aggregate
    aggregate.main()

    print(f"\n===========  STAGE 4 / {n} : NARRATIVE  ===========", flush=True)
    import narrative
    narrative.main()

    print(f"\n===========  STAGE 5 / {n} : HTML REPORT  ===========", flush=True)
    import report
    report.main()

    if run_robustness:
        print(f"\n===========  STAGE 6 / {n} : R1 K-SWEEP  ===========", flush=True)
        import robustness_ksweep
        robustness_ksweep.main()

        print(f"\n===========  STAGE 7 / {n} : R2 NO-OVERLAP  ===========", flush=True)
        import robustness_nooverlap
        robustness_nooverlap.main()

    print(f"\n[main] total wall time: {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    # Pass --no-robustness to skip the two robustness stages.
    run_rob = "--no-robustness" not in sys.argv
    main(run_robustness=run_rob)
