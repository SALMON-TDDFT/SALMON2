#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
taylor=(root/'src/rt/dg/rt_dg_nodal_taylor.f90').read_text()
salmon=root/'src/rt/dg/rt_dg_nodal_salmon_taylor.f90'
assert 'subroutine propagate_nodal_taylor_action' in taylor
assert 'procedure(nodal_taylor_hamiltonian_action) :: apply_hamiltonian' in taylor
assert 'call apply_hamiltonian(state,hterm)' in taylor
assert salmon.exists(),'missing SALMON complete-H Taylor driver'
body=salmon.read_text()
assert 'subroutine propagate_nodal_salmon_taylor' in body
assert 'call propagate_nodal_taylor_action' in body
assert 'call apply_nodal_salmon_hamiltonian_mpi' in body
print('PASS nodal GS and Taylor share the complete SALMON Hamiltonian callback')
