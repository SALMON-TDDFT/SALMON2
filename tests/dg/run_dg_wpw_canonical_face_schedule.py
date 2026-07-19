#!/usr/bin/env python3
from pathlib import Path
import os,subprocess
ROOT=Path(__file__).resolve().parents[2]
mod=Path('/tmp/test_dg_wpw_canonical_face_schedule_mod');mod.mkdir(exist_ok=True)
exe=Path('/tmp/test_dg_wpw_canonical_face_schedule')
subprocess.run([os.environ.get('FC','gfortran'),'-J',str(mod),'-I',str(mod),
  str(ROOT/'src/common/dg_wpw_canonical_face_schedule.f90'),
  str(ROOT/'tests/dg/test_dg_wpw_canonical_face_schedule.f90'),'-o',str(exe)],check=True,cwd='/tmp')
subprocess.run([str(exe)],check=True,cwd=ROOT)
