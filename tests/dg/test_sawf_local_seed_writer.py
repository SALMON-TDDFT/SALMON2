#!/usr/bin/env python3
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FC = shutil.which("gfortran")
if not FC:
    raise SystemExit("gfortran is required")

driver = r'''program check_local_seed_writer
  use lcfo_wannier_sawf_seed, only: write_sawf_local_eig_amn_mmn,build_sawf_local_seed_matrices
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  real(8) :: energy(2)
  complex(8) :: amn(2,2),mmn(2,2,3)
  complex(8) :: states(2,2),projections(2,2),phase(2,3)
  logical :: ok
  character(256) :: message
  integer :: i
  energy=[-1d0,2d0]
  states=(0d0,0d0);states(1,1)=1;states(2,2)=1
  projections=states;phase=(1d0,0d0)
  call build_sawf_local_seed_matrices(states,projections,phase,1d0,amn,mmn,ok,message)
  if(.not.ok.or.maxval(abs(amn-states))>1d-14)error stop 3
  if(maxval(abs(mmn(:,:,1)-states))>1d-14)error stop 4
  amn=(0d0,0d0);amn(1,1)=1;amn(2,2)=1
  mmn=(0d0,0d0)
  do i=1,3;mmn(1,1,i)=1;mmn(2,2,i)=1;end do
  call write_sawf_local_eig_amn_mmn('.', 'local',energy,amn,mmn,ok,message)
  if(.not.ok)then;write(*,'(a)')trim(message);error stop 1;end if
  mmn(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,kind=8)
  call write_sawf_local_eig_amn_mmn('.', 'bad',energy,amn,mmn,ok,message)
  if(ok)error stop 2
  write(*,'(a)')'PASS local SAWF eig/amn/mmn writer'
end program'''

with tempfile.TemporaryDirectory(prefix="sawf-local-seed-") as td:
    td = Path(td)
    (td / "driver.f90").write_text(driver)
    result = subprocess.run(
        [FC, str(ROOT / "src/gs/dc/lcfo_wannier_sawf_seed.f90"),
         str(td / "driver.f90"), "-o", str(td / "check")],
        cwd=td, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    if result.returncode:
        raise SystemExit(result.stdout)
    result = subprocess.run([str(td / "check")], cwd=td, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(result.stdout)
    eig = (td / "local.eig").read_text().splitlines()
    amn = (td / "local.amn").read_text().splitlines()
    mmn = (td / "local.mmn").read_text().splitlines()
    assert len(eig) == 2 and eig[0].split()[:2] == ["1", "1"]
    assert [int(x) for x in amn[1].split()] == [2, 1, 2]
    assert [int(x) for x in mmn[1].split()] == [2, 1, 3]
    assert "PASS local SAWF" in result.stdout
    print(result.stdout.strip())
