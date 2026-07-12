#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = root / 'src/rt/dg/rt_dg_nodal_cg.f90'
assert source.exists(), 'missing nodal callback CG eigensolver'
body = source.read_text()
assert 'subroutine solve_nodal_ground_state_cg_mpi' in body
assert 'procedure(nodal_hamiltonian_action) :: apply_hamiltonian' in body
assert 'search_direction' in body
assert 'previous_residual_norm' in body
assert 'call zhegv' in body
assert 'call orthonormalize_nodal_states_mpi' in body
assert 'call rotate_subspace(state,hpsi,eigenvalues,communicator)' in body
assert 'max_residual <= tolerance .and. mod(iteration-1,cg_block_length) == 0' in body
accepted = body.split('niteration=iteration', 1)[1].split('deallocate(', 1)[0]
assert 'orthonormalize_nodal_states_mpi' not in accepted
assert 'apply_hamiltonian' not in accepted

wrapper = (root / 'src/rt/dg/rt_dg_nodal_salmon_ground_state.f90').read_text()
assert 'use rt_dg_nodal_cg' in wrapper
assert 'call solve_nodal_ground_state_cg_mpi' in wrapper
assert 'relax_nodal_ground_state_action_mpi' not in wrapper

cmake = (root / 'src/rt/CMakeLists.txt').read_text()
assert 'dg/rt_dg_nodal_cg.f90' in cmake
print('PASS nodal complete-H ground state uses callback conjugate gradient')
