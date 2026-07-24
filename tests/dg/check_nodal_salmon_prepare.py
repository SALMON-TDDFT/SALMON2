#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
src=root/'src/rt/dg/rt_dg_nodal_salmon_prepare.f90'
assert src.exists(),'missing nodal SALMON preparation entry point'
body=src.read_text()
assert 'subroutine prepare_nodal_salmon_ground_state' in body
assert 'call initialize_nodal_from_dg_coefficients' in body
assert 'call build_nodal_local_potential' in body
assert 'call get_nodal_stencil_coefficients' in body
assert 'call relax_nodal_salmon_ground_state_mpi' in body
assert 'if (.not. state%ground_state_ready)' in body
print('PASS SALMON nodal preparation ends with a verified complete-H ground state')
