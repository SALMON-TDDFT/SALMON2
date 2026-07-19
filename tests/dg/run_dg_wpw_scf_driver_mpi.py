#!/usr/bin/env python3
from pathlib import Path
import subprocess,tempfile
ROOT=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='dg_wpw_scf_driver_') as tmp:
  exe=Path(tmp)/'test'
  subprocess.run(['mpifort','-J',tmp,str(ROOT/'src/gs/dc/dg_wpw_scf_driver.f90'),
    str(ROOT/'tests/dg/test_dg_wpw_scf_driver_mpi.f90'),'-o',str(exe)],check=True,cwd=tmp)
  subprocess.run(['mpirun','-np','2',str(exe)],check=True,cwd=ROOT,timeout=20)
