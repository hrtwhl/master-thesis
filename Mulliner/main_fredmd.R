# ---------------------------------------------------------------------------
# main_fredmd.R
#
# Top-level entry point for the FRED-MD extension to Mulliner et al. (2025).
# Replaces the seven hand-picked financial state variables with k principal-
# component factors of the McCracken-Ng (2016) FRED-MD panel of ~120 monthly
# US macroeconomic series. Everything else is identical to the baseline.
#
# Run from the project root:
#   source("main_fredmd.R")
#
# This script does NOT touch any of:
#   - main.R
#   - R/01_data.R, R/02_similarity.R, R/03_strategy.R, R/04_exhibits.R,
#     R/utils.R, R/diagnostics.R
#   - output/  (the baseline replication output)
#
# Output goes to:
#   output/fredmd/
#
# Required input file:
#   data/2025-12-md.csv  (the FRED-MD vintage from the St Louis Fed,
#                         downloaded once and never updated thereafter)
#
# After running, an .rds with the in-memory objects is saved to
# output/fredmd/fredmd_results.rds. To re-run diagnostics on the FRED-MD
# strategy, load that file and pass to run_diagnostics:
#
#   src <- readRDS("output/fredmd/fredmd_results.rds")
#   list2env(src, envir = .GlobalEnv)
#   source("R/diagnostics.R")
#   run_diagnostics(state_fmd, dist_mat_fmd, factors, strat_fmd, CFG_FREDMD)
# ---------------------------------------------------------------------------

source("Mulliner/07_run_fredmd.R")


src <- readRDS("output/fredmd/fredmd_results.rds")
list2env(src, envir = .GlobalEnv)
source("Mulliner/diagnostics_v2.R")
run_diagnostics(state_fmd, dist_mat_fmd, factors, strat_fmd, CFG_FREDMD)

alpha_diagnostics_windowed(strat_fmd, factors, CFG_FREDMD,
                           start_date = "2001-01-31",
                           label_suffix = " [FRED-MD, matched window]")
