"""Aggregate completed matched routing runs and plot measured time/memory."""
import argparse
import json
from pathlib import Path
import statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('results', type=Path, nargs='+')
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
runs = []
provenance = None
for path in a.results:
    data = json.loads(path.read_text())
    if not data['complete']:
        raise RuntimeError(f'Incomplete batch: {path}')
    identity = (data['mpi'], data['omp'], data['blas_threads'], data['steps'],
                {k: v['binary_sha256'] for k, v in data['variants']['variants'].items()})
    if provenance is not None and identity != provenance:
        raise RuntimeError('Different execution configuration/binaries')
    provenance = identity
    runs.extend(data['runs'])
labels = {'replicated': 'Replicated', 'ace': 'Distributed ACE', 'all': 'Distributed ACE + MLWF'}
summary = []
for states in sorted({r['states'] for r in runs}):
    subset = [r for r in runs if r['states'] == states]
    reference = next(r for r in subset if r['mode'] == 'replicated')
    for r in subset:
        if r['input_sha256'] != reference['input_sha256'] or r['seed'] != reference['seed']:
            raise RuntimeError('Different input/seed')
        r['current_diff'] = float(np.max(np.abs(np.asarray(r['observables'])[:, 13:16] -
                                               np.asarray(reference['observables'])[:, 13:16])))
        r['energy_diff'] = float(np.max(np.abs(np.asarray(r['energies'])[:, 1] -
                                              np.asarray(reference['energies'])[:, 1])))
        if r['returncode'] or max(r['current_diff'], r['energy_diff']) >= 2e-8:
            raise RuntimeError(f'Numerical mismatch: {r["directory"]}')
    for mode in labels:
        group = [r for r in subset if r['mode'] == mode]
        if not group:
            raise RuntimeError(f'Missing {states} {mode}')
        row = dict(states=states, mode=mode, count=len(group))
        for key in ['rt_seconds_per_step', 'rt_max_seconds', 'wall_seconds', 'peak_rank_bytes']:
            values = [r[key] for r in group]
            row[key] = dict(median=statistics.median(values), minimum=min(values), maximum=max(values))
        row['current_difference_max'] = max(r['current_diff'] for r in group)
        row['energy_difference_max_ha'] = max(r['energy_diff'] for r in group)
        row['directories'] = [r['directory'] for r in group]
        summary.append(row)
a.output.mkdir(parents=True, exist_ok=True)
(a.output / 'summary.json').write_text(json.dumps(dict(
    mpi=provenance[0], omp=provenance[1], blas_threads=provenance[2], steps=provenance[3],
    binaries=provenance[4], source_results=[str(p.resolve()) for p in a.results],
    summary=summary), indent=2) + '\n')
fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
states_list = sorted({r['states'] for r in summary})
colors = ['#555555', '#0072B2', '#D55E00']
for ax, key, scale, ylabel in zip(axes, ['rt_seconds_per_step', 'peak_rank_bytes'],
                                [1, 2**20], ['RT time (s/step)', 'Peak RSS per rank (MiB)']):
    for i, (mode, label) in enumerate(labels.items()):
        group = [next(r for r in summary if r['states'] == n and r['mode'] == mode) for n in states_list]
        y = np.array([r[key]['median'] / scale for r in group])
        lo = y - np.array([r[key]['minimum'] / scale for r in group])
        hi = np.array([r[key]['maximum'] / scale for r in group]) - y
        x = np.arange(len(group)) + (i - 1) * .25
        errors = dict(yerr=[lo, hi], capsize=3) if any(r['count'] > 1 for r in group) else {}
        ax.bar(x, y, width=.24, label=label, color=colors[i], **errors)
        for px, py in zip(x, y):
            ax.text(px, py, f'{py:.2f}' if scale == 1 else f'{py:.0f}', ha='center', va='bottom', fontsize=8)
    ax.set_xticks(range(len(states_list)), [f'{n} H₂' for n in states_list])
    ax.set_ylabel(ylabel)
    ax.grid(axis='y', alpha=.2)
    ax.set_axisbelow(True)
    ax.margins(y=.18)
axes[0].legend(fontsize=8)
sample_note = 'One run per condition' if all(r['count'] == 1 for r in summary) else 'Median and observed range'
fig.suptitle(f'PBEh(40) saved DC → Γ RT, 16 steps; MPI {provenance[0]} / OMP {provenance[1]}\n{sample_note}; RSS includes initialization', fontsize=11)
fig.savefig(a.output / 'time-memory.png', dpi=180)
for r in summary:
    print(r['states'], r['mode'], 'n=', r['count'], 's/step=', r['rt_seconds_per_step'],
          'MiB=', r['peak_rank_bytes']['median']/2**20, 'dE=', r['energy_difference_max_ha'])
