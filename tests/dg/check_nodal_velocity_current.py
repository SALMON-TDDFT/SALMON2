#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = root / 'src/rt/dg/rt_dg_nodal_current.f90'
assert source.exists(), 'missing nodal velocity-gauge current'
body = source.read_text()
assert 'subroutine calculate_nodal_velocity_current_mpi' in body
assert 'call exchange_nodal_face_halos' in body
assert 'gradient_offset' in body
assert 'system%vec_k' in body
assert 'system%vec_Ac' in body
assert 'ppg%rinv_uvu' in body
assert 'ppg%rxyz' in body
assert 'call MPI_Allreduce' in body
compact = body.replace(' ', '')
assert '/(dble(system%ngrid)*system%hvol)' in compact
assert 'current=current_global/dble(system%ngrid)' not in compact

main = (root / 'src/rt/main_tddft.f90').read_text()
assert 'use rt_dg_nodal_current' in main
assert 'call calculate_nodal_velocity_current_mpi' in main
assert "' current='" in main
print('PASS nodal velocity current includes local derivative, A, and nonlocal terms')
