#!/usr/bin/env python3
from pathlib import Path
import os,subprocess
ROOT=Path(__file__).resolve().parents[2]
mod=Path('/tmp/test_dg_wpw_trace_halo_periodic_faces_mod');mod.mkdir(exist_ok=True)
exe=Path('/tmp/test_dg_wpw_trace_halo_periodic_faces')
sources=['src/common/dg_wpw_owner_exchange.f90','src/common/dg_wpw_bounded_operator.f90',
 'src/common/dg_wpw_production_context.f90','src/rt/dg/rt_dg_wpw_face_trace_provider.f90',
 'src/rt/dg/rt_dg_wpw_weak_form_evaluator.f90','src/rt/dg/rt_dg_wpw_sparse_blocks.f90',
 'src/rt/dg/rt_dg_wpw_column_layout.f90','src/rt/dg/rt_dg_wpw_sparse_builder.f90',
 'src/rt/dg/rt_dg_wpw_quadrature_assembler.f90','src/rt/dg/rt_dg_wpw_face_trace_scanner.f90',
 'src/rt/dg/rt_dg_wpw_trace_halo_provider.f90','tests/dg/test_dg_wpw_trace_halo_periodic_faces.f90']
subprocess.run([os.environ.get('MPIFC','mpifort'),'-J',str(mod),'-I',str(mod),
  *[str(ROOT/s) for s in sources],'-o',str(exe)],check=True,cwd='/tmp')
subprocess.run([str(exe)],check=True,cwd=ROOT)
