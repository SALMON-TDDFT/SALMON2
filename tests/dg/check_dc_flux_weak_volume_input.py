#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
global_src = (root / "src/io/salmon_global.f90").read_text()
input_src = (root / "src/io/inputoutput.f90").read_text()
flux_src = (root / "src/gs/dc/lcfo_flux.f90").read_text()
name = "yn_dc_lcfo_flux_weak_volume"

assert f"character(1)   :: {name}" in global_src
for token in (
    f"& {name}, &",
    f"{name} = 'y'",
    f"call comm_bcast({name}, nproc_group_global)",
    f'"{name}",{name}',
    f"call yn_argument_check({name})",
):
    assert token in input_src, token
assert f"enabled = ({name} == 'y')" in flux_src
weak_start = flux_src.index("subroutine calc_hamiltonian_matrix")
weak_end = flux_src.index("end subroutine calc_hamiltonian_matrix", weak_start)
weak_body = flux_src[weak_start:weak_end]
assert "build_fragment_nonlocal_basis_action" in weak_body
assert "V_local(ispin)%f" in weak_body
assert "mat_H_volume_local(io,jo,ispin) &\n        & -" not in weak_body
assert "[DC-LCFO-SAWF-H-COVARIANCE]" in flux_src
assert "label='volume'" in flux_src
assert "label='volume_kinetic'" in flux_src
assert "label='volume_local'" in flux_src
assert "label='volume_nonlocal'" in flux_src
assert "label='surface_self'" in flux_src
assert "label='surface_cross'" in flux_src
assert "private(io,jo,idir,term_sum,term_local,term_nonlocal)" in flux_src
assert "[DC-LCFO-SAWF-VLOCAL-SYMMETRY]" in flux_src
assert "[DC-LCFO-SAWF-H-COVARIANCE-BLOCK]" in flux_src
print("PASS DC Flux weak-volume namelist wiring")
