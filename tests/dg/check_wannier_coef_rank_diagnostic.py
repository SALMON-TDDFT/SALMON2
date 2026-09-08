#!/usr/bin/env python3
from pathlib import Path


root = Path(__file__).resolve().parents[2]
src = (root / "src" / "gs" / "dc" / "lcfo_flux.f90").read_text()

markers = [
    "subroutine diagnose_wannier_coef_rank",
    "call diagnose_wannier_coef_rank(nband_wann)",
    "[DC-LCFO-WANNIER-COEF]",
    "gram_min=",
    "near_null=",
]

for marker in markers:
    if marker not in src:
        raise SystemExit(f"missing Wannier coef-rank diagnostic marker: {marker}")

print("Wannier export logs coef_wf Gram-rank diagnostics before AMN projection")
