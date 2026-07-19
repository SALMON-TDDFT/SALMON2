#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
MOD = Path('/tmp/test_dg_wpw_full_rank_seed_mod')
MOD.mkdir(exist_ok=True)
EXE = Path('/tmp/test_dg_wpw_full_rank_seed_mpi')
fc = os.environ.get('MPIFC', 'mpifort')
subprocess.run([
    fc, '-J', str(MOD), '-I', str(MOD),
    str(ROOT / 'src/common/dg_generalized_algebra.f90'),
    str(ROOT / 'src/common/dg_wpw_matrix_free_operator.f90'),
    str(ROOT / 'src/gs/dc/dg_wpw_matrix_free_scf.f90'),
    str(ROOT / 'tests/dg/test_dg_wpw_full_rank_seed.f90'),
    '-llapack', '-lblas', '-o', str(EXE),
], check=True, cwd='/tmp')
subprocess.run(['mpirun', '-np', '8', str(EXE)], check=True, cwd=ROOT, timeout=30)
