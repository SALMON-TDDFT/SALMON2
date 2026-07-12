#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
adapter=root/'src/rt/dg/rt_dg_nodal_salmon_adapter.f90'
assert adapter.exists(),'missing SALMON-to-nodal adapter'
body=adapter.read_text()
assert 'subroutine initialize_nodal_from_full_orbital' in body
assert 'subroutine initialize_nodal_from_dg_coefficients' in body
assert 'dg_frag%coef_global_to_local' in body
assert 'dg_frag%phi_frag' in body
assert 'psi_matrix=matmul(phi_matrix,coef_matrix)' in body
assert "dg_frag%nproc_frag /= 1" in body
assert 'info%imap(fragment_coords(1),fragment_coords(2),fragment_coords(3)' in body
assert 'spsi%zwf(ixg,iyg,izg,ispin,io,ik,im)' in body
assert 'spsi%rwf(ixg,iyg,izg,ispin,io,ik,im)' in body
assert 'call initialize_nodal_face_topology' in body
assert 'count(system%rocc(:,1,ispin) > 1.0d-12)' in body
assert 'subroutine build_nodal_local_potential' in body
assert 'Vh%f(ixg,iyg,izg) + Vxc(ispin)%f(ixg,iyg,izg) + Vpsl%f(ixg,iyg,izg)' in body
assert 'kinetic_center = stencil%coef_lap0' in body
assert 'kinetic_offset = -0.5d0 * stencil%coef_lap' in body
assert 'gradient_offset = stencil%coef_nab' in body
print('PASS SALMON orbitals, potential, and stencil map into the nodal core route')
