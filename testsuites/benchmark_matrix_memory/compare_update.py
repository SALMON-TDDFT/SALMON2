"""Compare one updated ACE run with a saved, matched ACE baseline."""
import argparse
import json
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--before', type=Path, required=True)
p.add_argument('--after', type=Path, required=True)
p.add_argument('--states', type=int, default=128)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
datasets = [json.loads(path.read_text()) for path in [a.before, a.after]]
for d in datasets:
    if not d['complete']:
        raise RuntimeError('Batch is incomplete')
for key in ['mpi', 'omp', 'blas_threads', 'steps', 'dt']:
    if datasets[0][key] != datasets[1][key]:
        raise RuntimeError('Configuration mismatch: '+key)
rows = []
for d in datasets:
    matches = [r for r in d['runs'] if r['states'] == a.states and r['mode'] == 'ace']
    if len(matches) != 1:
        raise RuntimeError('Expected exactly one ACE sample per batch')
    rows.append(matches[0])
before, after = rows
for key in ['input_sha256', 'seed']:
    if before[key] != after[key]:
        raise RuntimeError('Workload mismatch: '+key)
if any(r['returncode'] for r in rows):
    raise RuntimeError('Failed run')
dj = float(np.max(np.abs(np.asarray(before['observables'])[:, 13:16] -
                         np.asarray(after['observables'])[:, 13:16])))
de = float(np.max(np.abs(np.asarray(before['energies'])[:, 1] -
                         np.asarray(after['energies'])[:, 1])))
if max(dj, de) >= 2e-8:
    raise RuntimeError(f'Observable mismatch: current {dj}, energy {de}')
keys = ['directory', 'wall_seconds', 'rt_seconds_per_step', 'rt_max_seconds', 'peak_rank_bytes']
result = dict(states=a.states, mpi=datasets[0]['mpi'], omp=datasets[0]['omp'],
              before={k: before[k] for k in keys}, after={k: after[k] for k in keys},
              input_sha256=before['input_sha256'], seed=before['seed'],
              current_difference=dj, energy_difference_ha=de,
              memory_reduction_percent=100*(1-after['peak_rank_bytes']/before['peak_rank_bytes']),
              time_change_percent=100*(after['rt_seconds_per_step']/before['rt_seconds_per_step']-1),
              before_results=str(a.before.resolve()), after_results=str(a.after.resolve()))
a.output.parent.mkdir(parents=True, exist_ok=True)
a.output.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
