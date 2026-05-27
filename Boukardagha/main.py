"""
main.py
=======
End-to-end driver:
    1) run pure-market Wasserstein-HMM backtest      (run_pure.py)
    2) run hierarchical macro+market backtest        (run_hier.py)
    3) aggregate -> tables + charts + benchmarks     (aggregate.py)
    4) build narrative summary                       (narrative.py)
    5) build HTML report                             (report.py)

Each stage is also a standalone script:

    python run_pure.py        # pure-market backtest only
    python run_hier.py        # hierarchical backtest only
    python aggregate.py       # reads raw CSVs and plots
    python narrative.py       # builds Markdown summary
    python report.py          # builds HTML report
"""
import os, sys, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")


def main():
    t0 = time.time()

    print("===========  STAGE 1 / 5 : PURE MARKET  ===========", flush=True)
    import run_pure
    run_pure.main()

    print("\n===========  STAGE 2 / 5 : HIERARCHICAL  ===========", flush=True)
    import run_hier
    run_hier.main()

    print("\n===========  STAGE 3 / 5 : AGGREGATE  ===========", flush=True)
    import aggregate
    aggregate.main()

    print("\n===========  STAGE 4 / 5 : NARRATIVE  ===========", flush=True)
    import narrative
    narrative.main()

    print("\n===========  STAGE 5 / 5 : HTML REPORT  ===========", flush=True)
    import report
    report.main()

    print(f"\n[main] total wall time: {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
