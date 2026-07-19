#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT=Path(__file__).resolve().parents[2]
mod=Path('/tmp/test_dg_wpw_volume_halo_schedule_mod');mod.mkdir(exist_ok=True)
exe=Path('/tmp/test_dg_wpw_volume_halo_schedule_mpi')
subprocess.run([os.environ.get('MPIFC','mpifort'),'-J',str(mod),'-I',str(mod),
  str(ROOT/'src/rt/dg/rt_dg_wpw_volume_halo_provider.f90'),
  str(ROOT/'tests/dg/test_dg_wpw_volume_halo_schedule_mpi.f90'),'-o',str(exe)],check=True,cwd='/tmp')
subprocess.run(['mpirun','-np','4',str(exe)],check=True,cwd=ROOT,timeout=20)
