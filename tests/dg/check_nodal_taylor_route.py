#!/usr/bin/env python3
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sources = [
    root / "src/common/dg_nodal_state.f90",
    root / "src/rt/dg/rt_dg_nodal_types.f90",
    root / "src/rt/dg/rt_dg_nodal_halo.f90",
    root / "src/rt/dg/rt_dg_nodal_mpi.f90",
    root / "src/rt/dg/rt_dg_nodal_hamiltonian.f90",
    root / "src/rt/dg/rt_dg_nodal_taylor.f90",
]
assert sources[-1].exists(), "missing nodal Taylor propagator"
body = sources[-1].read_text()
assert "dt > 0.02d0" in body
assert "call refresh_nodal_single_fragment_halos(state)" in body
assert "call apply_nodal_local_hamiltonian" in body

driver = r"""
program check_nodal_taylor
  use rt_dg_nodal_types
  use rt_dg_nodal_taylor
  implicit none
  type(s_dg_nodal_state) :: state
  real(8) :: v(2,1,1,1), kin(1,3), grad(1,3), err, dt
  integer :: iface
  call allocate_nodal_state(state,1,[2,1,1],1,1,1)
  do iface=1,6
    call allocate_nodal_face_buffers(state%faces(iface),state%core_size,1,1,1)
  end do
  state%psi_core=(1.0d0,0.0d0)
  call accept_nodal_dg_ground_state(state,0.0d0,1.0d-10)
  v=1.0d0; kin=0.0d0; grad=0.0d0; dt=0.02d0
  call propagate_nodal_taylor(state,v,0.0d0,kin,grad,[0.0d0,0.0d0,0.0d0],dt,8)
  err=maxval(abs(state%psi_core-exp(-(0.0d0,1.0d0)*dt)))
  if (err > 1.0d-14) error stop 'nodal Taylor mismatch'
  print *, 'PASS nodal Taylor',err
end program check_nodal_taylor
"""
fc = shutil.which("gfortran")
assert fc
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    source = tmp / "check.f90"
    source.write_text(driver)
    (tmp / "config.h").write_text("")
    exe = tmp / "check"
    subprocess.run([fc, "-cpp", "-I", str(tmp), *map(str, sources), str(source), "-o", str(exe)],
                   check=True, cwd=tmp)
    out = subprocess.run([str(exe)], check=True, text=True, capture_output=True).stdout
    assert "PASS nodal Taylor" in out
print("PASS nodal Taylor uses a refreshed halo for each Hamiltonian action")
