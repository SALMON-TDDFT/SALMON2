#!/usr/bin/env python3
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FC = shutil.which("gfortran")
if not FC:
    raise SystemExit("gfortran is required")

driver = r'''program check_local_win
  use lcfo_wannier_sawf_win, only: write_sawf_local_preprocess_win
  implicit none
  real(8) :: lattice(3,3),atoms(3,2)
  logical :: ok
  character(256) :: message
  lattice=0d0;lattice(1,1)=4d0;lattice(2,2)=5d0;lattice(3,3)=6d0
  atoms=reshape([0.1d0,0.2d0,0.3d0,0.6d0,0.7d0,0.8d0],[3,2])
  call write_sawf_local_preprocess_win('local.win',4,2,lattice,atoms,ok,message)
  if(.not.ok)then;write(*,'(a)')trim(message);error stop 1;end if
  write(*,'(a)')'PASS atomic local SAWF preprocess WIN writer'
end program'''

with tempfile.TemporaryDirectory(prefix="sawf-local-win-") as td:
    td = Path(td);(td / "driver.f90").write_text(driver)
    result = subprocess.run(
        [FC, str(ROOT / "src/gs/dc/lcfo_wannier_sawf_win.f90"),str(td / "driver.f90"),
         "-o",str(td / "check")],cwd=td,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if result.returncode: raise SystemExit(result.stdout)
    result = subprocess.run([str(td / "check")],cwd=td,text=True,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if result.returncode: raise SystemExit(result.stdout)
    text=(td / "local.win").read_text().lower()
    for token in ("num_bands = 4","num_wann = 2","mp_grid = 1 1 1","gamma_only = true",
                  "begin unit_cell_cart","begin atoms_frac","begin kpoints"):
        assert token in text, token
    assert not list(td.glob("local.win.tmp.*"))
    print(result.stdout.strip())
