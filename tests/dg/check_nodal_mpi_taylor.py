#!/usr/bin/env python3
import os, shutil, subprocess, tempfile
from pathlib import Path

root=Path(__file__).resolve().parents[2]
sources=[root/'src/rt/dg/rt_dg_nodal_types.f90',root/'src/rt/dg/rt_dg_nodal_halo.f90',
         root/'src/rt/dg/rt_dg_nodal_mpi.f90',root/'src/rt/dg/rt_dg_nodal_hamiltonian.f90',
         root/'src/rt/dg/rt_dg_nodal_taylor.f90']
body=sources[-1].read_text()
assert 'subroutine propagate_nodal_taylor_mpi' in body
assert 'call exchange_nodal_face_halos(state, communicator)' in body

driver=r"""
program check_nodal_mpi_taylor
  use mpi
  use rt_dg_nodal_types
  use rt_dg_nodal_halo
  use rt_dg_nodal_taylor
  implicit none
  type(s_dg_nodal_state) :: state
  integer :: ierr,rank,nproc,ix,ig,iface
  integer :: ranks(2)
  real(8) :: pi,k,eig,dt,err,global_err,v(4,1,1,1),kin(1,3),grad(1,3)
  complex(8) :: initial(4),expected
  call MPI_Init(ierr); call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) error stop 'test requires two ranks'
  ranks=[0,1]
  call allocate_nodal_state(state,rank+1,[4,1,1],1,1,1)
  call initialize_nodal_face_topology(state%faces,[rank,0,0],[2,1,1],ranks,1)
  do iface=1,6
    call allocate_nodal_face_buffers(state%faces(iface),state%core_size,1,1,1)
  end do
  pi=acos(-1.0d0); k=2.0d0*pi/8.0d0; eig=1.0d0-cos(k); dt=0.02d0
  do ix=1,4
    ig=4*rank+ix
    initial(ix)=exp((0.0d0,1.0d0)*k*dble(ig-1))
    state%psi_core(ix,1,1,1,1)=initial(ix)
  end do
    call accept_nodal_dg_ground_state(state,0.0d0,1.0d-10)
  v=0.0d0; kin=0.0d0; grad=0.0d0; kin(1,1)=-0.5d0
  call propagate_nodal_taylor_mpi(state,v,1.0d0,kin,grad,[0.0d0,0.0d0,0.0d0],dt,8,MPI_COMM_WORLD)
  err=0.0d0
  do ix=1,4
    expected=initial(ix)*exp(-(0.0d0,1.0d0)*eig*dt)
    err=max(err,abs(state%psi_core(ix,1,1,1,1)-expected))
  end do
  call MPI_Allreduce(err,global_err,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(global_err>1.0d-13) error stop 'MPI Taylor mismatch'
  if(rank==0) print *,'PASS nodal MPI Taylor',global_err
  call MPI_Finalize(ierr)
end program
"""
mpifort=shutil.which('mpifort'); mpiexec=shutil.which('mpiexec'); assert mpifort and mpiexec
with tempfile.TemporaryDirectory() as tmp:
    tmp=Path(tmp); src=tmp/'check.f90'; src.write_text(driver); exe=tmp/'check'
    subprocess.run([mpifort,'-cpp','-I',str(root/'cmakefiles/build-mpi'),*map(str,sources),str(src),'-o',str(exe)],check=True)
    env=os.environ.copy(); env.setdefault('OMPI_MCA_rmaps_base_oversubscribe','1')
    out=subprocess.run([mpiexec,'-n','2',str(exe)],check=True,text=True,capture_output=True,env=env).stdout
    assert 'PASS nodal MPI Taylor' in out
print('PASS every MPI Taylor Hamiltonian action refreshes fragment face halos')
