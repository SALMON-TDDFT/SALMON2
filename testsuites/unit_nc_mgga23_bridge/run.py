#!/usr/bin/env python3
"""Link a synthetic bridge test against one freshly audited SALMON build."""
from pathlib import Path
import argparse
import os
import shlex
import subprocess
import tempfile
import math

p = argparse.ArgumentParser()
p.add_argument('--build-dir', type=Path, required=True)
p.add_argument('--checker', type=Path, required=True,
               help='workspace tools/check_forbidden_eigensolvers.py')
p.add_argument('--fc', default='gfortran')
p.add_argument('--mpi', action='store_true')
p.add_argument('--rspace', action='store_true', help='two-rank real-space regression (implies --mpi)')
p.add_argument('--launcher', default='mpiexec')
p.add_argument('--libs', default='-lopenblas')
a = p.parse_args()
if a.rspace:
    a.mpi = True
root = Path(__file__).resolve().parents[2]
build = a.build_dir.resolve()
objroot = build / 'src/CMakeFiles/salmon.dir'
subprocess.run(['python3', str(a.checker), '--source-root', str(root / 'src')], check=True)
subprocess.run(['python3', str(a.checker), '--object-tree', str(objroot)], check=True)
objects = sorted(x for x in objroot.rglob('*.o') if x.name != 'main.f90.o')
if not objects:
    raise SystemExit('no application objects found')
with tempfile.TemporaryDirectory(prefix='nc-mgga23-bridge-') as tmp:
    tmp = Path(tmp)
    test_obj = tmp / 'test.o'
    exe = tmp / 'test_bridge'
    subprocess.run([a.fc, '-cpp', *(['-DTEST_MPI'] if a.mpi else []), '-fopenmp', '-ffree-line-length-none', '-fcheck=all',
                    '-I' + str(build), '-c', str(Path(__file__).with_name('test_rspace.f90' if a.rspace else 'test_bridge.f90')),
                    '-o', str(test_obj)], check=True)
    subprocess.run(['python3', str(a.checker), '--object-tree', str(tmp)], check=True)
    subprocess.run([a.fc, '-fopenmp', str(test_obj), *map(str, objects),
                    *shlex.split(a.libs), '-o', str(exe)], check=True)
    launch = [a.launcher, '-n', '2'] if a.mpi else []
    subprocess.run([*launch, str(exe)], check=True,
                   env=dict(os.environ, OMP_NUM_THREADS='1'), cwd=tmp)

    if a.rspace:
        raise SystemExit(0)
    norm = (tmp / 'synthetic_rt_magnetization_norm.data').read_text()
    if 'NOT total XC torque' not in norm or 'integral_abs_local_spin' not in norm:
        raise SystemExit('spin norm diagnostic header missing semantics')
    rows = [list(map(float, line.split())) for line in norm.splitlines() if line and not line.startswith('#')]
    if len(rows) != 1 or len(rows[0]) != 6 or not all(map(math.isfinite, rows[0])):
        raise SystemExit('invalid spin norm diagnostic sample')
    if rows[0][2] + 1e-12 < rows[0][1]:
        raise SystemExit('integral local spin smaller than total spin')
    print('PASS actual write_rt_spin diagnostic:', rows[0])
    for mode in ['--nonzero-core', '--nonzero-core-tau']:
        rejected = subprocess.run([*launch, str(exe), mode], capture_output=True, text=True,
                                  env=dict(os.environ, OMP_NUM_THREADS='1'), cwd=tmp)
        if rejected.returncode == 0 or 'consistent core tau/current model' not in rejected.stderr:
            raise SystemExit('expected explicit nonzero NLCC rejection: ' + mode)
        print('PASS explicit NLCC rejection:', mode)
