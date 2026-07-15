#!/usr/bin/env python3
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
path=ROOT/"src/rt/dg/rt_dg_wpw_grid_point_adapter.f90"
faces=ROOT/"src/rt/dg/rt_dg_wpw_canonical_faces.f90"
assert path.exists(), "missing real-grid WPW point adapter"
assert faces.exists(), "missing canonical fragment-face enumerator"
src=path.read_text().lower()
for token in ("evaluate_local_wannier_grid_point", "evaluate_windowed_kg_grid_point",
              "global_wannier_local_ids", "global_wannier_local_coef",
              "gradient_basis_cache", "wpw_chi", "wpw_grad_chi",
              "wpw_window_box_lo", "wpw_window_box_hi",
              "evaluate_windowed_kg_point", "evaluate_wannier_point"):
    assert token in src, f"missing grid adapter contract: {token}"
for token in ("size(dg_frag%k_pw,1)/=3", "grid_context_extents_ready",
              "size(dg_frag%gradient_basis_cache,5)", "size(dg_frag%gradient_basis_cache,6)"):
    assert token in src, f"grid adapter does not fail closed on cache extents: {token}"
for token in ("dg_frag%ifrag_start<1", "dg_frag%ifrag_end>dg_frag%n_frag",
              "size(dg_frag%wpw_chi,4)/=nlocal",
              "size(dg_frag%wpw_grad_chi,5)/=nlocal",
              "window_extent(1)>size(dg_frag%wpw_chi,1)",
              "window_extent(3)>size(dg_frag%wpw_grad_chi,4)"):
    assert token in src, f"window cache shape is not fail closed: {token}"
for forbidden in ("h_mat_pw", "s_mat_pw", "allocate(h_", "allocate(s_"):
    assert forbidden not in src, f"grid adapter creates a dense operator: {forbidden}"
assert "call wpw_normalized_window_at_grid" not in src, \
    "production adapter must not scan every fragment at each grid point"
face_src=faces.read_text().lower()
for token in ("s_wpw_canonical_face_list", "build_wpw_canonical_face_list",
              "wpw_face_neighbor_fragment", "k_minus", "k_plus", "side_from_k_minus"):
    assert token in face_src, f"missing canonical-face contract: {token}"
assert "if (ifrag >= jfrag) cycle" in face_src
assert "do ifrag=dg_frag%ifrag_start,dg_frag%ifrag_end" in face_src
assert "do ifrag=1,dg_frag%n_frag" not in face_src
print("PASS rank-local real-grid W/(K,G) point adapter contract")
