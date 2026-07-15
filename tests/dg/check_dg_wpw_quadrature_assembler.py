#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
src = ROOT / "src/rt/dg/rt_dg_wpw_quadrature_assembler.f90"
point = ROOT / "src/rt/dg/rt_dg_wpw_point_evaluator.f90"
assert src.exists(), "missing WPW quadrature assembler"
assert point.exists(), "missing WPW point evaluator"
text = src.read_text().lower()
for token in ("assemble_wpw_volume_point", "assemble_wpw_canonical_face_point",
              "pack_wpw_point_candidates", "s_dg_wpw_sparse_candidates",
              "wpw_candidate_volume_face", "wpw_volume_weak_pair", "wpw_sipg_face_pair", "k_minus", "k_plus"):
    assert token in text, f"missing quadrature assembler contract: {token}"
assert "pp_face" not in text, "quadrature assembler invents a PP face contribution"
assert "k_minus >= k_plus" in text, "face is not canonicalized by fragment id"
point_text = point.read_text().lower()
for token in ("evaluate_windowed_kg_point", "evaluate_wannier_point", "sqrt(omega_cell)", "grad_chi"):
    assert token in point_text, f"missing normalized point-evaluator contract: {token}"
assert "select rank" not in point_text and "values(..)" not in point_text, "point evaluator exceeds Fortran 2008"
print("PASS support-local WP/PP quadrature assembler source contract")
