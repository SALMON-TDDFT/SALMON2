#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
solver=(root/'src/rt/dg/rt_dg_nodal_ground_state_solver.f90').read_text()
interfaces=(root/'src/common/dg_nodal_interfaces.f90').read_text()
salmon=(root/'src/rt/dg/rt_dg_nodal_salmon_ground_state.f90')
assert 'use dg_nodal_interfaces' in solver
assert 'abstract interface' in interfaces
assert 'subroutine nodal_complete_h_action' in interfaces
assert 'subroutine relax_nodal_ground_state_action_mpi' in solver
assert 'procedure(nodal_hamiltonian_action) :: apply_hamiltonian' in solver
assert 'call apply_hamiltonian(state,hpsi)' in solver
assert 'call global_complex_array_sum' in solver
assert '[DG-NODAL-GS-ITER]' in solver
assert salmon.exists(),'missing SALMON complete-H ground-state driver'
body=salmon.read_text()
assert 'subroutine relax_nodal_salmon_ground_state_mpi' in body
assert 'call solve_nodal_ground_state_cg_mpi' in body
assert 'call apply_nodal_salmon_hamiltonian_mpi' in body
print('PASS nodal GS solver accepts the complete SALMON Hamiltonian action')
