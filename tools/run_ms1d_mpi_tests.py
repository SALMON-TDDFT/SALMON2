"""Bounded 2-point SALMON/serial/MPI comparison and communication measurements.
python tools/run_ms1d_mpi_tests.py EMBEDDED_BUILD_JSON PREPARED_SOURCE NEW_OUTPUT
"""
import argparse
import json
import os
from pathlib import Path
import sys
import numpy as np
from run_parallel_mpi_tests import run_bounded

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('build'); p.add_argument('source'); p.add_argument('output')
a = p.parse_args()
out = Path(a.output).resolve(); out.mkdir(parents=True, exist_ok=False)
root = Path(__file__).resolve().parents[1]
record = Path(a.build).resolve()
exe = record.parent/'build'/'salmon'
if not exe.is_file(): raise ValueError('embedded build must retain its standalone SALMON executable')
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
           SALMONPY_RUNTIME_SOURCE=str(Path(a.source).resolve()))
reference = out/'reference'
cases = [('reference', [sys.executable, str(root/'tools/test_longitudinal_salmon.py'),
                       str(exe), str(reference), '--trans-longi', 'ms1d', '--cells', '2'])]
for ranks in (3, 5):
    cases.append((f'mpi{ranks}', ['mpiexec', '-n', str(ranks), sys.executable,
        str(root/'tools/test_ms1d_mpi.py'), str(record), str(reference), str(out/f'mpi{ranks}')]))
cases.append(('contracts', ['mpiexec', '-n', '3', sys.executable, str(root/'tools/test_ms1d_contract_mpi.py')]))
for name, command in cases:
    with (out/(name+'.log')).open('w') as log:
        run_bounded(command, env=env, log=log, timeout=180)
    print(name, 'PASS', flush=True)
with np.load(out/'mpi3/history.npz') as one, np.load(out/'mpi5/history.npz') as two:
    errors = {}
    for key in ('A', 'current', 'grid_A'):
        np.testing.assert_allclose(one[key], two[key], rtol=2e-12, atol=2e-14)
        errors[key] = float(np.max(np.abs(one[key]-two[key])))
summary = {'scope': 'synthetic interface parity, not converged material physics',
           'material_rank_parity_max_absolute': errors,
           'reference': json.loads((reference/'comparison.json').read_text()),
           'communication': {str(n): json.loads((out/f'mpi{n}/communication.json').read_text()) for n in (3, 5)}}
(out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print('rank parity PASS', errors)
