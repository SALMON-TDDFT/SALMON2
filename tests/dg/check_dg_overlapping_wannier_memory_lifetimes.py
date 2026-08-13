#!/usr/bin/env python3
"""Protect the production overlapping-Wannier peak-memory lifetimes."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (ROOT / "src/gs/main_dft.f90").read_text()


def require_between(token: str, start: str, end: str, message: str) -> None:
    first = SOURCE.index(start)
    last = SOURCE.index(end, first)
    if token not in SOURCE[first:last]:
        raise AssertionError(message)


require_between(
    "deallocate(ow_box_values,ow_box_gradients)",
    "ow_core_gradients(:,:,core_index)=ow_box_gradients(:,:,p)",
    "call invert_ow_lattice",
    "the initial full buffer must be released immediately after core extraction",
)

if "w90_anchors=global_seed_values" in SOURCE:
    raise AssertionError("the retained seed frame must not be copied into a duplicate anchor frame")
if "call move_alloc(global_seed_values,w90_anchors)" not in SOURCE:
    raise AssertionError("the retained seed allocation must be moved into the W90 anchor owner")

require_between(
    "deallocate(w90_anchors,w90_fractional)",
    "call fingerprint_ow_w90_matrices",
    "call run_dg_w90_gamma_library",
    "assembly-only W90 anchor/fractional arrays must be released before localization",
)

require_between(
    "deallocate(w90_m_matrix,w90_a_matrix,w90_eigenvalues)",
    "call run_dg_w90_gamma_library",
    "call apply_dg_w90_gamma_transform",
    "coordinator matrices must be released immediately after Wannier90 returns",
)

require_between(
    "deallocate(lcfo_occupied_core)",
    "occupied_density_before=sum(abs(lcfo_occupied_core",
    "global_seed_values(1:nstate,:)=adapted_occupied_candidates",
    "the occupied LCFO source must be released after its final density receipt",
)

print("PASS overlapping-Wannier production array lifetimes are bounded")
