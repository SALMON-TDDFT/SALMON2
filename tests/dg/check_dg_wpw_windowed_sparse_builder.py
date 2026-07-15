#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
layout = (ROOT / "src/rt/dg/rt_dg_wpw_column_layout.f90").read_text().lower()
builder_path = ROOT / "src/rt/dg/rt_dg_wpw_sparse_builder.f90"

assert "pw_fragment_ids" in layout, "PW columns still lack the K index in P_(K,G)"
assert "pw_g_ids" in layout, "PW columns still lack a stable G index in P_(K,G)"
assert builder_path.exists(), "missing direct windowed sparse WPW builder"
builder = builder_path.read_text().lower()
for token in (
    "build_windowed_sparse_wpw_operators",
    "wpw_normalized_window_at_grid",
    "wpw_grad_chi",
    "s_dg_wpw_sparse_blocks",
    "wpw_column_owner",
    "pw_fragment_ids",
    "pw_g_ids",
    "wpw_candidate_volume",
    "wpw_candidate_face",
    "pp_origin",
):
    assert token in builder, f"missing windowed sparse-builder token: {token}"
for forbidden in (
    "allocate(s_frag_pw(n_frag, n_pw",
    "allocate(h_frag_pw(n_frag, n_pw",
    "allocate(h_pw(n_pw, n_pw",
    "allocate(s_pw(n_pw, n_pw",
):
    assert forbidden not in builder, f"sparse builder allocates global dense operator: {forbidden}"
assert "info = 12" in builder, "non-owned candidates are not rejected"
assert "info = 13" in builder, "PP face candidates are not rejected"
print("PASS direct windowed sparse WPW builder contract")
