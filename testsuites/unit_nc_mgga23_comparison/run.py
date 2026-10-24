#!/usr/bin/env python3
"""Audit, link and run deterministic synthetic LSDA/MGGA trajectories and compare."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--build-dir', type=Path, required=True)
p.add_argument('--checker', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--fc', default='gfortran')
p.add_argument('--libs', default='-lopenblas')
p.add_argument('--plot', action='store_true')
a = p.parse_args()
root = Path(__file__).resolve().parents[2]
objroot = a.build_dir.resolve() / 'src/CMakeFiles/salmon.dir'
objects = sorted(x for x in objroot.rglob('*.o') if x.name != 'main.f90.o')
if not objects:
    raise SystemExit('no SALMON application objects found')
if a.output.exists():
    raise SystemExit('output must be a new directory')


def run(args, **kwargs):
    return subprocess.run([str(x) for x in args], check=True, **kwargs)

run(['python3', a.checker, '--source-root', root / 'src'])
run(['python3', a.checker, '--object-tree', objroot])
with tempfile.TemporaryDirectory(prefix='nc-mgga23-comparison-') as scratch:
    scratch = Path(scratch)
    obj = scratch / 'comparison.o'
    exe = scratch / 'comparison'
    run([a.fc, '-fopenmp', '-ffree-line-length-none', '-fcheck=all',
         '-I' + str(a.build_dir.resolve()), '-c', Path(__file__).with_name('test_comparison.f90'), '-o', obj])
    run(['python3', a.checker, '--object-tree', scratch])
    run([a.fc, '-fopenmp', obj, *objects, *shlex.split(a.libs), '-o', exe])
    result = run([exe], cwd=scratch, env=dict(os.environ, OMP_NUM_THREADS='1'), capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print(result.stderr)
    for label in ['LSDA', 'MGGA23']:
        data = (scratch / f'{label}_rt_spin.data').read_text()
        rows = [line for line in data.splitlines() if line and not line.startswith('#')]
        if len(rows) != 21 or float(rows[0].split()[0]) != 0.0:
            raise SystemExit('expected initial-inclusive 21-row native spin stream')
    a.output.mkdir(parents=True, exist_ok=False)
    for filename in ['LSDA_rt_spin.data', 'MGGA23_rt_spin.data',
                     'LSDA_rt_magnetization_norm.data', 'MGGA23_rt_magnetization_norm.data', 'validation.txt']:
        shutil.copy2(scratch / filename, a.output / filename)
    (a.output / 'run.log').write_text(result.stdout + result.stderr)
    metadata = {'scope': 'synthetic deterministic non-ground-state spinors; no material, no GS collection',
                'executable_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
                'sources': str(root), 'build': str(a.build_dir.resolve()),
                'source_audit': 'PASS', 'full_application_object_audit': 'PASS',
                'test_object_audit': 'PASS', 'samples_per_run': 21, 'dt_au': .002, 'steps': 20,
                'same_initial_state': True, 'own_functional_ground_states': False,
                'ionic_soc': False, 'scalar_field': '.02*sin(5*t)*cos(2*pi*x/8)',
                'compiler': a.fc, 'taylor_order': 6, 'omp_threads': 1}
    (a.output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
run(['python3', root / 'tools/compare_spin.py', '--run', f'LSDA={a.output}/LSDA_rt_spin.data',
     '--run', f'MGGA23={a.output}/MGGA23_rt_spin.data', '--output', a.output / 'comparison',
     *(['--plot'] if a.plot else [])])
print((a.output / 'validation.txt').read_text())
