#!/usr/bin/env python3
"""Production route contract for divided WF+PW SCF and one-shot LCFO."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()
BRANCH_START = "if(yn_dg_hybrid_divided_scf=='y')"

assert BRANCH_START in SOURCE, "missing divided WF+PW SCF production branch"
branch = SOURCE[SOURCE.index(BRANCH_START) :]
branch = branch.split("endif", 1)[0]

scf_position = branch.find("call run_dg_hybrid_divided_scf")
lcfo_position = branch.find("call assemble_dg_hybrid_lcfo_rows")
assert scf_position >= 0, "divided branch must run fragment WF+PW SCF"
assert lcfo_position > scf_position, "LCFO assembly must follow divided SCF"
assert branch.count("solve_dg_hybrid_generalized_") == 1, (
    "divided branch must perform exactly one generalized eigensolve"
)
assert "run_dg_hybrid_self_consistent_ground_state" not in branch, (
    "divided branch must not use the repeated full-cell Hybrid SCF"
)

for forbidden in (
    "allocate(ow_hybrid_hrows(ntarget,ntarget)",
    "allocate(ow_hybrid_srows(ntarget,ntarget)",
    "complex(8)::hybrid_h(ntarget,ntarget)",
    "complex(8)::hybrid_s(ntarget,ntarget)",
):
    assert forbidden not in branch.replace(" ", ""), (
        f"divided branch must not replicate a dense LCFO matrix: {forbidden}"
    )

print("divided WF+PW LCFO route contract: PASS")
