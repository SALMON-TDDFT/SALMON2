#!/usr/bin/env python3
import shutil
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[2]
types = root / "src/rt/dg/rt_dg_nodal_types.f90"
common = root / "src/common/dg_nodal_state.f90"
operator = root / "src/rt/dg/rt_dg_nodal_hamiltonian.f90"
assert operator.exists(), "missing matrix-free nodal Hamiltonian"

driver = r"""
program check_nodal_h
  use rt_dg_nodal_types
  use rt_dg_nodal_hamiltonian
  implicit none
  type(s_dg_nodal_state) :: state
  complex(8) :: hpsi(4,1,1,1,1), ref, phase
  real(8) :: v(4,1,1,1), kin(1,3), grad(1,3), pi, err
  integer :: ix, iface

  pi = acos(-1.0d0)
  call allocate_nodal_state(state, 1, [4,1,1], 1, 1, 1)
  do ix=1,4
    state%psi_core(ix,1,1,1,1)=exp((0.0d0,1.0d0)*2.0d0*pi*dble(ix-1)/4.0d0)
  end do
  v=0.0d0; kin=0.0d0; grad=0.0d0
  kin(1,1)=-0.5d0
  do iface=1,6
    call allocate_nodal_face_buffers(state%faces(iface), state%core_size, 1, 1, 1)
  end do
  state%faces(nodal_face_slot(1,-1))%recv_value(1,1,1,1,1)=state%psi_core(4,1,1,1,1)
  state%faces(nodal_face_slot(1,+1))%recv_value(1,1,1,1,1)=state%psi_core(1,1,1,1,1)
  call apply_nodal_local_hamiltonian(state,v,1.0d0,kin,grad,[0.0d0,0.0d0,0.0d0],hpsi)
  err=0.0d0
  do ix=1,4
    phase=state%psi_core(ix,1,1,1,1)
    ref=phase-0.5d0*(state%psi_core(modulo(ix,4)+1,1,1,1,1)+ &
      state%psi_core(modulo(ix-2,4)+1,1,1,1,1))
    err=max(err,abs(hpsi(ix,1,1,1,1)-ref))
  end do
  if (err > 1.0d-12) error stop 'decomposed nodal kinetic mismatch'
  print *, 'PASS nodal Hamiltonian action', err
end program check_nodal_h
"""

fc = shutil.which("gfortran")
assert fc, "gfortran is required"
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    source = tmp / "check.f90"
    source.write_text(driver)
    (tmp / "config.h").write_text("")
    exe = tmp / "check"
    subprocess.run([fc, "-cpp", "-I", str(tmp), str(common), str(types), str(operator), str(source),
                    "-o", str(exe)], check=True, cwd=tmp)
    result = subprocess.run([str(exe)], check=True, text=True, capture_output=True)
    assert "PASS nodal Hamiltonian action" in result.stdout
print("PASS decomposed periodic nodal Hamiltonian matches reference")
