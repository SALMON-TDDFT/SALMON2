#!/usr/bin/env python3
import os, shutil, subprocess, sys, tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[2]
sources=[root/'src/common/dg_nodal_state.f90',root/'src/common/dg_nodal_interfaces.f90',
 root/'src/rt/dg/rt_dg_nodal_types.f90',root/'src/rt/dg/rt_dg_nodal_halo.f90',
 root/'src/rt/dg/rt_dg_nodal_mpi.f90',root/'src/rt/dg/rt_dg_nodal_hamiltonian.f90',
 root/'src/rt/dg/rt_dg_nodal_ground_state.f90',root/'src/rt/dg/rt_dg_nodal_ground_state_solver.f90']
assert sources[-1].exists(),'missing matrix-free nodal ground-state solver'
body=sources[-1].read_text()
assert 'subroutine relax_nodal_dg_ground_state_mpi' in body
assert 'call exchange_nodal_face_halos' in body
assert 'hpsi(:,:,:,istate,ispin) - eigenvalues(istate,ispin)' in body
assert 'call orthonormalize_nodal_states_mpi' in body
orthonormalize = body.split('subroutine orthonormalize_nodal_states_mpi', 1)[1].split(
 'end subroutine orthonormalize_nodal_states_mpi', 1
)[0]
assert "call zgemm('C','N'" in orthonormalize
assert "call zpotrf('U'" in orthonormalize
assert "call ztrsm('R','U','N','N'" in orthonormalize
assert 'do jstate = 1, istate - 1' not in orthonormalize
assert 'call verify_nodal_dg_eigenstate_mpi' in body
driver=r"""
program check_nodal_gs
 use mpi
 use rt_dg_nodal_types
 use rt_dg_nodal_halo
 use rt_dg_nodal_ground_state_solver
 implicit none
 type(s_dg_nodal_state)::state
 integer::ierr,rank,nproc,ix,iface,niter
 integer::ranks(2)
 real(8)::v(4,1,1,1),kin(1,3),grad(1,3),eval(1,1),res
 call MPI_Init(ierr); call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr); call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
 if(nproc/=2) error stop 'requires two ranks'
 ranks=[0,1]; call allocate_nodal_state(state,rank+1,[4,1,1],1,1,1)
 call initialize_nodal_face_topology(state%faces,[rank,0,0],[2,1,1],ranks,1)
 do iface=1,6
  call allocate_nodal_face_buffers(state%faces(iface),state%core_size,1,1,1)
 end do
 do ix=1,4
  state%psi_core(ix,1,1,1,1)=cmplx(1d0+0.17d0*dble(4*rank+ix),0.11d0*dble(ix-rank),8)
 end do
 v=0d0; kin=0d0; grad=0d0; kin(1,1)=-0.5d0
 call relax_nodal_dg_ground_state_mpi(state,v,1d0,kin,grad,0.2d0,2000,1d-10,MPI_COMM_WORLD,eval,res,niter)
 if(.not.state%ground_state_ready .or. abs(eval(1,1))>1d-9 .or. res>1d-10) error stop 'nodal GS failed'
 if(rank==0) print *,'PASS nodal matrix-free GS',eval(1,1),res,niter
 call MPI_Finalize(ierr)
end program
"""
fc=shutil.which('mpifort'); run=shutil.which('mpiexec'); assert fc and run
with tempfile.TemporaryDirectory() as tmp:
 tmp=Path(tmp); src=tmp/'check.f90'; src.write_text(driver); exe=tmp/'check'
 (tmp/'config.h').write_text('')
 linear_algebra=['-framework','Accelerate'] if sys.platform=='darwin' else ['-llapack','-lblas']
 subprocess.run([fc,'-cpp','-DUSE_MPI','-I',str(tmp),*map(str,sources),str(src),
                 *linear_algebra,'-o',str(exe)],check=True,cwd=tmp)
 env=os.environ.copy(); env.setdefault('OMPI_MCA_rmaps_base_oversubscribe','1')
 out=subprocess.run([run,'-n','2',str(exe)],check=True,text=True,capture_output=True,env=env).stdout
 assert 'PASS nodal matrix-free GS' in out
print('PASS matrix-free fragment Hamiltonian generates its own stationary initial state')
