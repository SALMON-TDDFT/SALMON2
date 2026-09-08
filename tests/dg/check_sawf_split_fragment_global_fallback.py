#!/usr/bin/env python3
from pathlib import Path

source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()
required = [
    "SAWF-GLOBAL-SPLIT",
    "split_fragment_global_mode",
    ".not.split_fragment_global_mode",
    "build_sawf_dmn_operation_fragment_local",
]
missing = [token for token in required if token not in source]
if missing:
    raise SystemExit("SAWF split-fragment global fallback missing: " + ", ".join(missing))
print("PASS SAWF split-fragment global D_band fallback wiring")
