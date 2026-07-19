#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
mod = Path('/tmp/test_dg_wpw_fragment_root_layout_mod')
mod.mkdir(exist_ok=True)
exe = Path('/tmp/test_dg_wpw_fragment_root_layout_mpi')
subprocess.run([
    os.environ.get('MPIFC', 'mpifort'), '-J', str(mod), '-I', str(mod),
    str(ROOT / 'src/rt/dg/rt_dg_wpw_column_layout.f90'),
    str(ROOT / 'src/common/dg_wpw_owner_exchange.f90'),
    str(ROOT / 'src/common/dg_wpw_bounded_operator.f90'),
    str(ROOT / 'src/common/dg_wpw_production_context.f90'),
    str(ROOT / 'tests/dg/test_dg_wpw_fragment_root_layout_mpi.f90'),
    '-o', str(exe),
], check=True, cwd='/tmp')
subprocess.run(['mpirun', '-np', '4', str(exe)], check=True, cwd=ROOT)
