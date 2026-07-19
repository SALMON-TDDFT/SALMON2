#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT=Path(__file__).resolve().parents[2]
MOD=Path('/tmp/test_dg_wpw_nonlocal_projector_mod');MOD.mkdir(exist_ok=True)
EXE=Path('/tmp/test_dg_wpw_nonlocal_projector')
fc=os.environ.get('MPIFC','mpifort')
subprocess.run([fc,'-J',str(MOD),'-I',str(MOD),
    str(ROOT/'src/common/dg_wpw_nonlocal_projector.f90'),
    str(ROOT/'tests/dg/test_dg_wpw_nonlocal_projector.f90'),
    '-fcheck=all','-o',str(EXE)],check=True,cwd='/tmp')
subprocess.run([str(EXE)],check=True,cwd=ROOT,timeout=10)
