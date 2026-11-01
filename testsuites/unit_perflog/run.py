#!/usr/bin/env python3
"""Exercise production timer/perflog output with deterministic rank timings."""
from pathlib import Path
import argparse
import csv
import math
import os
import re
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('--compiler', default=os.environ.get('MPIFC', 'mpifort'))
parser.add_argument('--mpiexec', default=os.environ.get('MPIEXEC', 'mpiexec'))
args = parser.parse_args()
here = Path(__file__).resolve().parent
root = here.parents[1]
checker = root / 'tools/check_forbidden_eigensolvers.py'
subprocess.run(['python3', str(checker), '--source-root', str(root / 'src')], check=True)
with tempfile.TemporaryDirectory(prefix='perflog-regression-') as directory:
    build = Path(directory)
    flags = ['-cpp', '-fopenmp', '-ffpe-trap=invalid,zero,overflow', '-fcheck=all',
             '-O0', '-g', '-J', str(build), '-I', str(build)]
    sources = [here / 'stubs.f90', root / 'src/math/math_constants.f90',
               root / 'src/io/salmon_global.f90', root / 'src/misc/timer.f90',
               root / 'src/io/perflog.f90', here / 'probe.f90']
    objects = []
    for source in sources:
        obj = build / (source.stem + '.o')
        subprocess.run([args.compiler, *flags, '-c', str(source), '-o', str(obj)],
                       check=True, cwd=build)
        objects.append(str(obj))
    subprocess.run(['python3', str(checker), '--object-tree', str(build)], check=True)
    executable = build / 'probe'
    subprocess.run([args.compiler, *flags, *objects, '-o', str(executable)],
                   check=True, cwd=build)
    for ranks in (1, 2, 4):
        case = build / str(ranks)
        case.mkdir()
        subprocess.run([args.mpiexec, '-n', str(ranks), str(executable)], check=True,
                       cwd=case, env=dict(os.environ, OMP_NUM_THREADS='1'), timeout=60)
        text = (case / 'readable.rank00000000.txt').read_text()
        csv_text = (case / 'csv.rank00000000.txt').read_text()
        assert not re.search(r'\b(?:NaN|Infinity)\b', text + csv_text, re.I)
        rows = {row[0].strip(): [value.strip() for value in row[1:]]
                for row in csv.reader(csv_text.splitlines()) if len(row) == 7
                and row[0].strip() != 'scf calculation'
                and row[1].strip() != 'min'}
        # All-zero timers remain omitted; ordinary finite ratios are unchanged.
        assert 'reading data' not in rows
        assert float(rows['scf iterations'][-1]) == 1
        values = [float(x) for x in rows['writing scf data']]
        # E12.5E2 output has five significant digits.
        expected_mean = (ranks + 1) / 2
        expected_std = math.sqrt((ranks * ranks - 1) / 12)
        assert math.isclose(values[2], expected_mean, rel_tol=5e-5)
        assert math.isclose(values[3], expected_std, rel_tol=5e-5, abs_tol=1e-12)
        assert values[-1] == ranks
        # Centered variance retains differences between nearly equal timers.
        close_std = float(rows['gs initialization'][3])
        assert math.isclose(close_std, expected_std * 1e-9, rel_tol=5e-5, abs_tol=1e-15)
        if ranks == 1:
            assert 'init. scf iter.' not in rows
        else:
            assert rows['init. scf iter.'][-1] == 'N/A'
            assert rows['writing LDA data'][-1] == 'N/A'  # overflow, positive minimum
            assert 'N/A' in text
            z = [float(x) for x in rows['init. scf iter.'][:-1]]
            assert z[0] == 0 and z[1] == 1e-6
            assert math.isclose(z[2], (ranks - 1) / ranks * 1e-6, rel_tol=5e-5)
            assert math.isclose(z[3], math.sqrt(ranks - 1) / ranks * 1e-6, rel_tol=5e-5)
        print(f'PERFLOG_REGRESSION_PASS ranks={ranks}', flush=True)
