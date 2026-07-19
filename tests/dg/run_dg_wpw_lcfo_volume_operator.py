#!/usr/bin/env python3
from pathlib import Path
import subprocess,tempfile
ROOT=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='dg_wpw_volume_') as tmp:
  exe=Path(tmp)/'test'
  subprocess.run(['gfortran',str(ROOT/'src/gs/dc/dg_wpw_lcfo_volume_operator.f90'),
    str(ROOT/'tests/dg/test_dg_wpw_lcfo_volume_operator.f90'),'-o',str(exe)],check=True,cwd=tmp)
  subprocess.run([str(exe)],check=True,cwd=tmp)
