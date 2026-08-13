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
    "deallocate(ow_box_values)",
    "ow_core_values(:,core_index)=ow_box_values(:,p)",
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

w90_source = (ROOT / "src/gs/dc/dg_overlapping_wannier_w90.f90").read_text()
gamma_transform = w90_source[
    w90_source.index("subroutine apply_dg_w90_gamma_transform") :
    w90_source.index("end subroutine apply_dg_w90_gamma_transform")
]
if "new_values(:,:)" in gamma_transform or "new_gradients(:,:,:)" in gamma_transform:
    raise AssertionError("the Gamma transform must not allocate full output-sized temporaries")
if "point_values(:)" not in gamma_transform or "point_gradients(:,:)" not in gamma_transform:
    raise AssertionError("the Gamma transform must use bounded one-point work vectors")

adapter_start = SOURCE.index("call materialize_ow_distributed_core_to_buffer")
adapter_end = SOURCE.index("call validate_dg_factored_point_cogroup_gauge", adapter_start)
pre_final_gauge = SOURCE[adapter_start:adapter_end]
if "call periodic_box_gradients" in pre_final_gauge:
    raise AssertionError("production must not form gradients before the final character gauge")
if "ow_core_gradients(3,ntarget,ncore)" in pre_final_gauge:
    raise AssertionError("production must not retain a pre-final-gauge core gradient tensor")
if "transform=w90_transform" not in pre_final_gauge:
    raise AssertionError("production must invoke the Gamma transform through its values-only contract")
if "optional::gradients(:,:,:)" not in gamma_transform.replace(" ", ""):
    raise AssertionError("the Gamma transform gradient payload must be optional")

print("PASS overlapping-Wannier production array lifetimes are bounded")
