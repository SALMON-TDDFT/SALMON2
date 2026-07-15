#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
provider = ROOT / "src/rt/dg/rt_dg_wpw_face_trace_provider.f90"

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

print("PASS WPW face trace provider source contract")
