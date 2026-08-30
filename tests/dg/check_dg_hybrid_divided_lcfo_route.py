#!/usr/bin/env python3
"""Production route contract for divided WF+PW SCF and one-shot LCFO."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text(errors="replace").lower()
BRANCH_START = "if(yn_dg_hybrid_divided_scf=='y'.or.yn_dg_hybrid_continuation_scf=='y')then"

assert BRANCH_START in SOURCE, "missing divided WF+PW SCF production branch"
branch = SOURCE[SOURCE.index(BRANCH_START) :]
continuation_branch = "if(yn_dg_hybrid_continuation_scf=='y')then"
assert continuation_branch in branch, "missing continuation/divided route selection"
branch = branch[branch.index(continuation_branch) :]
assert "else" in branch, "missing divided-only route"
branch = branch.split("else", 1)[1]
branch = branch.split("if(yn_dg_hybrid_scf=='y')then", 1)[0]

scf_position = branch.find("call run_dg_hybrid_divided_scf")
lcfo_position = branch.find("call assemble_dg_hybrid_lcfo_rows")
assert scf_position >= 0, "divided branch must run fragment WF+PW SCF"
assert lcfo_position > scf_position, "LCFO assembly must follow divided SCF"
assert branch.count("solve_dg_hybrid_generalized_") == 1, (
    "divided branch must perform exactly one generalized eigensolve"
)
checkpoint_position = branch.find("call write_rt_dg_hybrid_occupied_checkpoint")
assert checkpoint_position > lcfo_position, "occupied checkpoint must follow the one-shot LCFO solve"
assert branch.count("write_rt_dg_hybrid_occupied_checkpoint") == 1
for forbidden in (
    "reconstruct_dg_hybrid_density",
    "post_lcfo_density",
    "run_dg_hybrid_divided_scf",  # checked separately below for exactly the pre-LCFO call
):
    if forbidden == "run_dg_hybrid_divided_scf":
        assert branch.count(forbidden) == 1, "divided SCF must not repeat after LCFO"
    else:
        assert forbidden not in branch, f"divided route added post-LCFO density work: {forbidden}"
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
