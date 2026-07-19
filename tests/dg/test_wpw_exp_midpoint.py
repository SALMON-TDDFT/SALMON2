#!/usr/bin/env python3
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='wpw_exp_midpoint_') as tmp:
    tmp=Path(tmp);exe=tmp/'test_wpw_exp_midpoint'
    subprocess.run([
        'gfortran','-std=f2008','-J',str(tmp),
        str(ROOT/'src/rt/dg/rt_dg_wpw_exp_production.f90'),
        str(ROOT/'tests/dg/test_wpw_exp_midpoint.f90'),
        '-llapack','-lblas','-o',str(exe)
    ],check=True,cwd=tmp)
    subprocess.run([str(exe)],check=True,cwd=tmp)
