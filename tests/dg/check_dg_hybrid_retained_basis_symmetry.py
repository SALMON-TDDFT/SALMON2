#!/usr/bin/env python3
"""Contract for the projected WF+PW retained-basis symmetry action."""

from pathlib import Path


root = Path(__file__).resolve().parents[2]
main = (root / "src/gs/main_dft.f90").read_text(errors="replace").lower()
source = (
    root / "src/gs/dc/dg_hybrid_retained_basis_symmetry.f90"
).read_text(errors="replace").lower()
acceptance_runner = (
    root / "tests/dg/run_dg_hybrid_si64_continuation_rt.py"
).read_text(errors="replace").lower()

name = "subroutine build_dg_hybrid_retained_basis_representation"
assert name in source, "missing retained WF+PW symmetry representation builder"
assert "call build_dg_hybrid_retained_basis_representation" in main
for token in (
    "assemble_dg_distributed_basis_symmetry_overlap",
    "fragment_basis_arg%buffer_values",
    "s_rows_arg",
    "call zgesv",
    "matmul(metric,representation_arg",
    "matmul(conjg(transpose(representation_arg",
):
    assert token in source, f"retained-basis action is missing {token}"
assert "representation_arg=identity" not in source, (
    "nontrivial retained-basis action must not fall back to the identity"
)
packed_source = "".join(source.split())
assert "if(noperation==0)then" in packed_source
assert "allocate(representation_arg(nbasis,nbasis,0))" in packed_source
assert "closure_defect_arg=0d0;closure_ok_arg=.true.;callback_ok=.true." in packed_source, (
    "identity-only systems must expose a vacuous nonidentity-generator action"
)
for token in ("closure_defect_arg", "closure_ok_arg"):
    assert token in source, f"retained-basis diagnostic omits {token}"
assert "callback_ok=ierr==mpi_success.and.ieee_is_finite(global_defect)" in source.replace(" ", ""), (
    "finite full-basis closure failure must remain structurally usable"
)
assert "closure_ok_arg=callback_ok.and.global_defect<=tolerance_arg" in source.replace(" ", ""), (
    "full-basis closure result is not reported separately"
)
for token in (
    "divided_full_basis_closure_defect",
    "divided_full_basis_closure_ok",
    "[hybrid-retained-basis-symmetry]",
):
    assert token in main, f"production route omits retained-basis diagnostic {token}"
assert "if(.not.divided_full_basis_closure_ok)error stop" not in main.replace(" ", ""), (
    "finite full-basis nonclosure still stops LCFO"
)
assert "retained_basis_receipt" in acceptance_runner, (
    "Si64 acceptance does not require the retained-basis diagnostic"
)
assert "closed=0" in acceptance_runner, (
    "Si64 acceptance self-test does not permit finite full-basis nonclosure"
)

print("retained WF+PW symmetry representation contract: PASS")
