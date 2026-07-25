#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = root / 'src/rt/dg/rt_dg_nodal_cg.f90'
assert source.exists(), 'missing nodal callback CG eigensolver'
body = source.read_text()
assert 'subroutine solve_nodal_ground_state_cg_mpi' in body
assert 'procedure(nodal_hamiltonian_action) :: apply_hamiltonian' in body
assert 'logical, intent(in), optional :: require_convergence' in body
assert 'must_converge=.' in body
assert 'search_direction' in body
assert 'previous_residual_norm' in body
assert 'call zhegv' in body
assert 'call project_search_outside_occupied' not in body
assert 'hpsi(:,:,:,istate,ispin)=h2(1,1)*hpsi(:,:,:,istate,ispin)' in body
assert 'call orthonormalize_nodal_states_mpi' in body
assert 'call rotate_subspace(state,hpsi,eigenvalues,communicator)' in body
assert 'if(must_converge) stop' in body
accepted = body.split('niteration=min(iteration,max_iteration)', 1)[1].split('deallocate(', 1)[0]
assert accepted.index('orthonormalize_nodal_states_mpi') < accepted.index('apply_hamiltonian')
assert accepted.index('apply_hamiltonian') < accepted.index('rotate_subspace')
assert accepted.index('rotate_subspace') < accepted.index('build_residuals')

wrapper = (root / 'src/rt/dg/rt_dg_nodal_salmon_ground_state.f90').read_text()
assert 'use rt_dg_nodal_cg' in wrapper
assert 'call solve_nodal_ground_state_cg_mpi' in wrapper
assert 'relax_nodal_ground_state_action_mpi' not in wrapper

cmake = (root / 'src/rt/CMakeLists.txt').read_text()
assert 'dg/rt_dg_nodal_cg.f90' in cmake
print('PASS nodal complete-H ground state uses callback conjugate gradient')
