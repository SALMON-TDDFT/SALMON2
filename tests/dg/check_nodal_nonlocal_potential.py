#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
src=root/'src/rt/dg/rt_dg_nodal_nonlocal.f90'
assert src.exists(),'missing nodal nonlocal pseudopotential action'
body=src.read_text()
assert 'subroutine apply_nodal_nonlocal_potential_mpi' in body
assert 'conjg(ppg%zekr_uV(j,ilma,ik)) * state%psi_core' in body
assert 'ppg%rinv_uvu(ilma)' in body
assert 'MPI_Allreduce' in body
assert 'ppg%zekr_uV(j,ilma,ik) * projector_global' in body
assert 'global_to_local_core_index' in body
wrapper=(root/'src/rt/dg/rt_dg_nodal_salmon_hamiltonian.f90').read_text()
assert 'subroutine apply_nodal_salmon_hamiltonian_mpi' in wrapper
assert 'call exchange_nodal_face_halos' in wrapper
assert 'call apply_nodal_local_hamiltonian' in wrapper
assert 'call apply_nodal_nonlocal_potential_mpi' in wrapper
print('PASS SALMON nodal Hamiltonian combines face, local, and nonlocal actions')
