#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
provider = ROOT / "src/rt/dg/rt_dg_wpw_face_trace_provider.f90"
scanner = ROOT / "src/rt/dg/rt_dg_wpw_face_trace_scanner.f90"
fixture = ROOT / "tests/dg/test_dg_wpw_face_trace_scanner.f90"
runner = ROOT / "tests/dg/run_dg_wpw_face_trace_fixture.py"

assert provider.exists(), "missing WPW face trace provider"
text = provider.read_text().lower()
for token in (
    "module rt_dg_wpw_face_trace_provider",
    "type, public :: s_wpw_face_trace_provider",
    "abstract interface",
    "class(*), pointer",
    "procedure(wpw_face_trace_callback), pointer, nopass",
    "bind_wpw_face_trace_provider",
    "unbind_wpw_face_trace_provider",
    "evaluate_wpw_face_traces",
    "k_minus", "k_plus", "axis", "side_from_k_minus", "unwrapped_grid",
    "w_minus", "w_plus", "grad_w_minus", "grad_w_plus",
    "p_minus", "p_plus", "grad_p_minus", "grad_p_plus", "info",
    "associated(provider%callback)", "associated(provider%user_context)",
):
    assert token in text, f"missing face trace provider contract: {token}"
for forbidden in ("mpi_", "h_mat", "s_mat", "dense_h", "dense_s", "pp_face"):
    assert forbidden not in text, f"provider crosses forbidden boundary: {forbidden}"

assert scanner.exists(), "missing WPW canonical face trace scanner"
scan = scanner.read_text().lower()
for token in (
    "module rt_dg_wpw_face_trace_scanner",
    "assemble_wpw_canonical_face_grid",
    "k_minus", "k_plus", "axis", "side_from_k_minus",
    "normal(axis)=dble(side_from_k_minus)",
    "h_normal=hgs(axis)",
    "face_weight=hgs(tangent(1))*hgs(tangent(2))",
    "evaluate_wpw_face_traces",
    "assemble_wpw_canonical_face_point",
    "strictly_increasing(w_row_ids)",
    "strictly_increasing(p_column_ids)",
    "temporary_block", "wp_face_h=temporary_block",
):
    assert token in scan, f"missing canonical face scanner contract: {token}"
for forbidden in ("do k=1,n_frag", "do ifrag=1,n_frag", "mpi_", "pp_face", "h_mat", "s_mat", "dense_h", "dense_s"):
    assert forbidden not in scan, f"scanner crosses forbidden boundary: {forbidden}"

assert fixture.exists() and runner.exists(), "missing linked canonical face fixture"
fixture_text = fixture.read_text().lower()
for token in ("point_count", "expected_grids", "deterministic_traces", "assemble_wpw_canonical_face_point"):
    assert token in fixture_text, f"numerical fixture misses contract: {token}"
assert "--build-dir" in runner.read_text(), "fixture runner must accept the configured build directory"

print("PASS WPW face trace provider/scanner source contract")
