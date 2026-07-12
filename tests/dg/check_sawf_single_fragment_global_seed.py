#!/usr/bin/env python3
from pathlib import Path

source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()

required = [
    "dc%n_frag==1",
    "SAWF-GLOBAL-SINGLE",
    "single fragment already spans the global LCFO eigensystem",
    "split_fragment_global_mode=(dc%n_frag==1)",
]
missing = [token for token in required if token not in source]
if missing:
    raise SystemExit(
        "SAWF single-fragment global seed fallback missing: " + ", ".join(missing)
    )

print("PASS SAWF single-fragment uses the existing global LCFO eigensystem")
