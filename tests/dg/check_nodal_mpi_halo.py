#!/usr/bin/env python3
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sources = [
    root / "src/rt/dg/rt_dg_nodal_types.f90",
    root / "src/rt/dg/rt_dg_nodal_halo.f90",
    root / "src/rt/dg/rt_dg_nodal_mpi.f90",
    root / "src/rt/dg/rt_dg_nodal_hamiltonian.f90",
]
assert sources[2].exists(), "missing nodal MPI face exchange"
assert '#include "config.h"' in sources[2].read_text()

driver = r"""
program check_nodal_mpi
  use mpi
  use rt_dg_nodal_types
  use rt_dg_nodal_halo
  use rt_dg_nodal_mpi
  use rt_dg_nodal_hamiltonian
  implicit none
  type(s_dg_nodal_state) :: state
  integer :: ierr,rank,nproc,ix,ig,iface
  integer :: ranks(2)
  real(8) :: pi,err,global_err,v(4,1,1,1),kin(1,3),grad(1,3)
  complex(8) :: hpsi(4,1,1,1,1),ref,phase
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if (nproc /= 2) error stop 'test requires two ranks'
  ranks=[0,1]
  call allocate_nodal_state(state,rank+1,[4,1,1],1,1,1)
  call initialize_nodal_face_topology(state%faces,[rank,0,0],[2,1,1],ranks,1)
  do iface=1,6
    call allocate_nodal_face_buffers(state%faces(iface),state%core_size,1,1,1)
  end do
  pi=acos(-1.0d0)
  do ix=1,4
    ig=4*rank+ix
    state%psi_core(ix,1,1,1,1)=exp((0.0d0,1.0d0)*2.0d0*pi*dble(ig-1)/8.0d0)
  end do
  call exchange_nodal_face_halos(state,MPI_COMM_WORLD)
  v=0.0d0; kin=0.0d0; grad=0.0d0; kin(1,1)=-0.5d0
  call apply_nodal_local_hamiltonian(state,v,1.0d0,kin,grad,[0.0d0,0.0d0,0.0d0],hpsi)
  err=0.0d0
  do ix=1,4
    ig=4*rank+ix
    phase=state%psi_core(ix,1,1,1,1)
    ref=phase-0.5d0*(exp((0.0d0,1.0d0)*2.0d0*pi*dble(modulo(ig,8))/8.0d0)+ &
      exp((0.0d0,1.0d0)*2.0d0*pi*dble(modulo(ig-2,8))/8.0d0))
    err=max(err,abs(hpsi(ix,1,1,1,1)-ref))
  end do
  call MPI_Allreduce(err,global_err,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if (global_err > 1.0d-12) error stop 'two-fragment halo mismatch'
  if (rank == 0) print *, 'PASS nodal MPI halo',global_err
  call MPI_Finalize(ierr)
end program check_nodal_mpi
"""

mpifort = shutil.which("mpifort")
mpiexec = shutil.which("mpiexec")
assert mpifort and mpiexec, "MPI compiler and launcher are required"
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    source = tmp / "check.f90"
    source.write_text(driver)
    exe = tmp / "check"
    subprocess.run([mpifort, "-cpp", "-I", str(root/"cmakefiles/build-mpi"), *map(str, sources), str(source), "-o", str(exe)], check=True)
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    out = subprocess.run([mpiexec, "-n", "2", str(exe)], check=True, text=True, capture_output=True, env=env).stdout
    assert "PASS nodal MPI halo" in out
print("PASS two-fragment MPI halo preserves the periodic Hamiltonian action")
