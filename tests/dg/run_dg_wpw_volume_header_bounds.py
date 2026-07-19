#!/usr/bin/env python3
from pathlib import Path
import os,subprocess
ROOT=Path(__file__).resolve().parents[2]
mod=Path('/tmp/test_dg_wpw_volume_header_bounds_mod');mod.mkdir(exist_ok=True)
exe=Path('/tmp/test_dg_wpw_volume_header_bounds')
subprocess.run([os.environ.get('MPIFC','mpifort'),'-J',str(mod),'-I',str(mod),
  str(ROOT/'src/rt/dg/rt_dg_wpw_volume_halo_provider.f90'),
  str(ROOT/'tests/dg/test_dg_wpw_volume_header_bounds.f90'),'-o',str(exe)],check=True,cwd='/tmp')
subprocess.run([str(exe)],check=True,cwd=ROOT)
