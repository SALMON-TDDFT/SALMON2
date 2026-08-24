#!/usr/bin/env python3
"""Contract for reusing authoritative DC controls in divided Hybrid SCF."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GLOBAL = (ROOT / "src/io/salmon_global.f90").read_text(errors="replace").lower()
INPUT = (ROOT / "src/io/inputoutput.f90").read_text(errors="replace").lower()
DCDTF = (ROOT / "src/gs/dc/dcdft.f90").read_text(errors="replace").lower()
MAIN = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()

assert "character(1)   :: yn_dg_hybrid_divided_scf" in GLOBAL
assert "yn_dg_hybrid_divided_scf = 'n'" in INPUT, "divided route must default off"
assert "call comm_bcast(yn_dg_hybrid_divided_scf" in INPUT
assert "call yn_argument_check(yn_dg_hybrid_divided_scf)" in INPUT

for token in (
    "subroutine prepare_dg_hybrid_divided_dc_controls",
    "convergence_mode=convergence",
    "density_threshold=threshold",
    "initial_total_density=dc%rho_tot%f",
):
    assert token in DCDTF, f"DC adapter is missing authoritative control: {token}"

assert "yn_dg_hybrid_divided_scf" in MAIN
branch_start = "if(yn_dg_hybrid_divided_scf=='y')then"
assert branch_start in MAIN, "missing default-off divided DC preparation branch"
branch = MAIN[MAIN.index(branch_start) :].split("endif", 1)[0]
assert "prepare_dg_hybrid_divided_dc_controls" in branch
assert "dc%rho_tot" in DCDTF
for forbidden in (
    "dg_hybrid_divided_density_tolerance",
    "dg_dc_gs_final_density_tolerance",
    "post_lcfo_density",
):
    assert forbidden not in branch, f"divided route introduced a new density gate: {forbidden}"

print("divided Hybrid DC controls contract: PASS")
