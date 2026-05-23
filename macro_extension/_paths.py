"""
_paths.py
---------
Makes the sibling `Boukardagha/` folder importable from within
`macro_extension_v10p2/`.

v10.2 is a cleaner architecture than v10/v10.1: instead of macro-gated
HMM transitions (which interact badly with the W2 template-mapping
when combined with long-horizon EMA drift), v10.2 keeps the baseline
market HMM intact and applies the macro signal as a *MVO input tilt*
at the allocation step. This preserves baseline's regime-detection
diversity while still incorporating macro-conditional return information.

Layout assumed
~~~~~~~~~~~~~~
    repo_root/
        Boukardagha/             <- baseline scripts
        macro_extension/         <- v9 (kept for ablation comparison)
        macro_extension_v10/     <- v10/v10.1 (kept for ablation)
        macro_extension_v10p2/   <- this folder (v10.2 — clean macro-tilt)
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_here = Path(__file__).resolve().parent
_parent = _here.parent

# Locate the Boukardagha folder
_baseline_default = (_parent / "Boukardagha").resolve()
_baseline_dir = Path(os.environ.get("BOUKARDAGHA_DIR", _baseline_default)).resolve()
if not _baseline_dir.exists():
    raise FileNotFoundError(
        f"Could not find baseline folder at {_baseline_dir}. "
        f"Expected layout: repo_root/Boukardagha sibling to "
        f"repo_root/macro_extension_v10p2, or set env var "
        f"BOUKARDAGHA_DIR to the correct absolute path."
    )

if str(_baseline_dir) not in sys.path:
    sys.path.insert(0, str(_baseline_dir))

BASELINE_DIR = _baseline_dir
