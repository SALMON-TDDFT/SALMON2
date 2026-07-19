#!/usr/bin/env python3
from pathlib import Path
import subprocess,tempfile
ROOT=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='dg_wpw_map_') as tmp:
  exe=Path(tmp)/'test'
  subprocess.run(['mpifort','-J',tmp,str(ROOT/'src/gs/dc/dg_wpw_lcfo_potential_map.f90'),
    str(ROOT/'tests/dg/test_dg_wpw_lcfo_potential_map_mpi.f90'),'-o',str(exe)],cwd=tmp,check=True)
  subprocess.run(['mpirun','-np','2',str(exe)],cwd=ROOT,check=True,timeout=20)
