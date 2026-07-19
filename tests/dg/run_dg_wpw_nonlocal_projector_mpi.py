#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT=Path(__file__).resolve().parents[2]
MOD=Path('/tmp/test_dg_wpw_nonlocal_projector_mpi_mod');MOD.mkdir(exist_ok=True)
EXE=Path('/tmp/test_dg_wpw_nonlocal_projector_mpi')
fc=os.environ.get('MPIFC','mpifort')
subprocess.run([fc,'-fcheck=all','-J',str(MOD),'-I',str(MOD),
    str(ROOT/'src/common/dg_wpw_nonlocal_projector.f90'),
    str(ROOT/'tests/dg/test_dg_wpw_nonlocal_projector_mpi.f90'),
    '-o',str(EXE)],check=True,cwd='/tmp')
subprocess.run(['mpirun','-np','2',str(EXE)],check=True,cwd=ROOT,timeout=20)
