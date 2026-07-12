#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
rr=root/'src/rt/dg/rt_dg_nodal_rayleigh_ritz.f90'
assert rr.exists(),'missing nodal Rayleigh-Ritz rotation'
body=rr.read_text()
assert 'subroutine rayleigh_ritz_nodal_subspace_mpi' in body
assert 'conjg(state%psi_core' in body
assert 'call MPI_Allreduce' in body
assert '[DG-NODAL-HERMITICITY]' in body
assert 'antiherm_rel=' in body
assert 'call eigen_zheev' in body
assert 'psi_matrix=matmul(psi_matrix,evec)' in body
assert 'hpsi_matrix=matmul(hpsi_matrix,evec)' in body
solver=(root/'src/rt/dg/rt_dg_nodal_ground_state_solver.f90').read_text()
assert 'procedure(nodal_subspace_rotation), optional :: rotate_subspace' in solver
assert 'call rotate_subspace(state,hpsi,eigenvalues,communicator)' in solver
salmon=(root/'src/rt/dg/rt_dg_nodal_salmon_ground_state.f90').read_text()
assert 'rayleigh_ritz_nodal_subspace_mpi' in salmon
assert 'eigenvalues,max_residual,niteration,rotate_complete_h_subspace' in salmon
print('PASS nodal GS removes occupied-subspace gauge rotations before residual evaluation')
