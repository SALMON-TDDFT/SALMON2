#!/usr/bin/env python3
import os, shutil, subprocess, tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[2]
sources=[root/'src/rt/dg/rt_dg_nodal_types.f90',root/'src/rt/dg/rt_dg_nodal_halo.f90',
 root/'src/rt/dg/rt_dg_nodal_mpi.f90',root/'src/rt/dg/rt_dg_nodal_hamiltonian.f90',
 root/'src/rt/dg/rt_dg_nodal_ground_state.f90']
assert sources[-1].exists(),'missing nodal DG ground-state verifier'
body=sources[-1].read_text()
assert 'subroutine verify_nodal_dg_eigenstate_mpi' in body
assert 'hpsi - eigenvalues' in body
assert 'call accept_nodal_dg_ground_state' in body
driver=r"""
program check_ground_residual
 use mpi
 use rt_dg_nodal_types
 use rt_dg_nodal_halo
 use rt_dg_nodal_ground_state
 implicit none
 type(s_dg_nodal_state)::state
 integer::ierr,rank,nproc,ix,ig,iface
 integer::ranks(2)
 real(8)::pi,k,eig,res,v(4,1,1,1),kin(1,3),grad(1,3),eval(1,1)
 call MPI_Init(ierr); call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr); call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
 if(nproc/=2) error stop 'requires two ranks'
 ranks=[0,1]; call allocate_nodal_state(state,rank+1,[4,1,1],1,1,1)
 call initialize_nodal_face_topology(state%faces,[rank,0,0],[2,1,1],ranks,1)
 do iface=1,6
  call allocate_nodal_face_buffers(state%faces(iface),state%core_size,1,1,1)
 end do
 pi=acos(-1d0); k=2d0*pi/8d0; eig=1d0-cos(k); eval(1,1)=eig
 do ix=1,4
  ig=4*rank+ix; state%psi_core(ix,1,1,1,1)=exp((0d0,1d0)*k*dble(ig-1))
 end do
 v=0d0; kin=0d0; grad=0d0; kin(1,1)=-0.5d0
 call verify_nodal_dg_eigenstate_mpi(state,eval,v,1d0,kin,grad,[0d0,0d0,0d0],MPI_COMM_WORLD,1d-12,res)
 if(.not.state%dg_ground_state_ready .or. res>1d-12) error stop 'ground-state verification failed'
 if(rank==0) print *,'PASS nodal DG eigen residual',res
 call MPI_Finalize(ierr)
end program
"""
fc=shutil.which('mpifort'); run=shutil.which('mpiexec'); assert fc and run
with tempfile.TemporaryDirectory() as tmp:
 tmp=Path(tmp); src=tmp/'check.f90'; src.write_text(driver); exe=tmp/'check'
 subprocess.run([fc,'-cpp','-I',str(root/'cmakefiles/build-mpi'),*map(str,sources),str(src),'-o',str(exe)],check=True)
 env=os.environ.copy(); env.setdefault('OMPI_MCA_rmaps_base_oversubscribe','1')
 out=subprocess.run([run,'-n','2',str(exe)],check=True,text=True,capture_output=True,env=env).stdout
 assert 'PASS nodal DG eigen residual' in out
print('PASS MPI eigen-residual verification gates nodal Taylor initialization')
