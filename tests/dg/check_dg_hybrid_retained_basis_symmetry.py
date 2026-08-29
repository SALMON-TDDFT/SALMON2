#!/usr/bin/env python3
"""Contract for the projected WF+PW retained-basis symmetry action."""

from pathlib import Path


root = Path(__file__).resolve().parents[2]
main = (root / "src/gs/main_dft.f90").read_text(errors="replace").lower()
source = (
    root / "src/gs/dc/dg_hybrid_retained_basis_symmetry.f90"
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

print("retained WF+PW symmetry representation contract: PASS")
