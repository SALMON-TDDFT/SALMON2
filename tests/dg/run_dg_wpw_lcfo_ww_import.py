#!/usr/bin/env python3
from pathlib import Path
import subprocess,tempfile
ROOT=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='dg_wpw_lcfo_ww_') as tmp:
    exe=Path(tmp)/'test'
    subprocess.run(['mpifort','-J',tmp,str(ROOT/'src/common/dg_wpw_owner_exchange.f90'),
      str(ROOT/'src/common/dg_wpw_bounded_operator.f90'),str(ROOT/'src/common/dg_wpw_production_context.f90'),
      str(ROOT/'src/gs/dc/dg_wpw_lcfo_operator_adapter.f90'),
      str(ROOT/'tests/dg/test_dg_wpw_lcfo_ww_import.f90'),'-o',str(exe)],check=True,cwd=tmp)
    subprocess.run([str(exe)],check=True)
