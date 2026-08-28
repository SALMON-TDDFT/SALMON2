#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/main_dft.f90").read_text()

flag = "yn_dg_hybrid_continuation_scf == 'y'"
assert flag in source, "missing explicit DG continuation production branch"
branch = source[source.index(flag):]
branch = branch[:branch.index("endif")]

required = [
    "dc%rho_tot",
    "close_dg_hybrid_selection",
    "divided_lcfo_srows",
    "dg_hybrid_interface_rows",
    "run_dg_hybrid_continuation_scf",
]
for token in required:
    assert token in branch, f"continuation branch does not use {token}"

for forbidden in [
    "dg_hybrid_production_continuation_adapter",
    "build_dg_hybrid_production_catalog",
    "s_dg_hybrid_continuation_callbacks",
    "solve_dg_hybrid_generalized_once_and_publish",
    "write_overlapping_wannier_occupied_checkpoint",
]:
    assert forbidden not in branch, f"continuation branch uses forbidden {forbidden}"

print("PASS concrete DG continuation production route contract")
