#!/usr/bin/env python3
from pathlib import Path
import os,subprocess
ROOT=Path(__file__).resolve().parents[2]
mod=Path('/tmp/test_dg_wpw_rank_local_quadrature_mod');mod.mkdir(exist_ok=True)
exe=Path('/tmp/test_dg_wpw_rank_local_quadrature')
subprocess.run([os.environ.get('MPIFC','mpifort'),'-J',str(mod),'-I',str(mod),
  str(ROOT/'src/rt/dg/rt_dg_wpw_weak_form_evaluator.f90'),
  str(ROOT/'src/rt/dg/rt_dg_wpw_sparse_blocks.f90'),
  str(ROOT/'src/rt/dg/rt_dg_wpw_column_layout.f90'),
  str(ROOT/'src/rt/dg/rt_dg_wpw_sparse_builder.f90'),
  str(ROOT/'src/rt/dg/rt_dg_wpw_quadrature_assembler.f90'),
  str(ROOT/'src/gs/dc/dg_wpw_rank_local_quadrature.f90'),
  str(ROOT/'tests/dg/test_dg_wpw_rank_local_quadrature.f90'),'-o',str(exe)],check=True,cwd='/tmp')
subprocess.run([str(exe)],check=True,cwd=ROOT)
