#!/usr/bin/env python3
"""Build and run the small kernel verification, retaining both mandatory audits."""
import os
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
checker = Path(os.environ.get('SALMON_FORBIDDEN_CHECKER', root / 'tools/check_forbidden_eigensolvers.py'))
compiler = os.environ.get('FC', 'gfortran')

def run(args):
    subprocess.run([str(x) for x in args], check=True, cwd=root)

run(['python3', checker, '--source-root', root / 'src'])
with tempfile.TemporaryDirectory(prefix='nc_mgga23_') as scratch:
    build = Path(scratch)
    flags = ['-std=f2008', '-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
             '-I', build, '-J', build]
    run([compiler, *flags, '-c', root / 'src/xc/nc_mgga23_kernel.f90', '-o', build / 'kernel.o'])
    run([compiler, *flags, '-c', root / 'testsuites/unit_nc_mgga23/test_kernel.f90', '-o', build / 'test.o'])
    run([compiler, build / 'kernel.o', build / 'test.o', '-o', build / 'test_kernel'])
    run(['python3', checker, '--object-tree', build])
    run([build / 'test_kernel'])
