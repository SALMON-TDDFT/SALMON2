"""Run an explicitly sampled A/c waveform from a pinned saved GS catalog.

Waveform NPZ: A shape (steps+1, 3), scalar dt, all in atomic units; A[0]=0.
Output is a new directory. No downloading, compilation or implicit SCF.
"""
import argparse
import json
from pathlib import Path
import numpy as np

from salmonpy.builtin_materials import builtin_material_registry
from salmonpy.config import load_materials
from salmonpy.errors import InputError
from salmonpy.initial_states import tddft_from_catalog
from salmonpy.runtime import load_runtime


def run(catalog_path, initial_state_id, runtime, waveform, directory):
    path = Path(catalog_path).resolve()
    with np.load(waveform, allow_pickle=False) as data:
        raw_A = np.asarray(data['A'])
        dt_array = np.asarray(data['dt'])
        if raw_A.dtype.kind not in 'iuf' or dt_array.dtype.kind not in 'iuf':
            raise InputError('A and dt must have real numeric dtypes, excluding bool')
        A = np.array(raw_A, dtype=float, copy=True)
        if dt_array.shape != ():
            raise InputError('dt must be scalar')
        dt = float(dt_array)
    if (A.ndim != 2 or A.shape[1] != 3 or len(A) < 2 or
            not np.isfinite(A).all() or not np.isfinite(dt) or dt <= 0 or np.any(A[0])):
        raise InputError('require finite A (steps+1,3), A[0]=0, and positive finite dt')
    engine = load_runtime(runtime, required=('dft', 'tddft', 'single-rank'))
    catalog = load_materials(path, registry=builtin_material_registry())
    with tddft_from_catalog(catalog, initial_state_id, base_directory=path.parent,
                           executable=engine, max_steps=len(A)-1, threads=1) as initial:
        output = Path(directory).resolve()
        output.mkdir(parents=True, exist_ok=False)
        report = {**initial.record(), 'status': 'running', 'executable_sha256': engine.sha256,
                  'drive': 'A/c in atomic units', 'current': 'charge-current density j = -J_SALMON'}
        try:
            with initial.new_state(dt=dt) as state:
                current = [state.current.copy()]
                for field in A[1:]:
                    current.append(state.step(field).copy())
            np.savez(output/'current.npz', time=np.arange(len(A))*dt, A=A, current=current)
            report['status'] = 'passed'
        except BaseException as exc:
            report.update(status='failed', error=str(exc))
            raise
        finally:
            (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('catalog', 'initial_state_id', 'runtime', 'waveform', 'directory'):
        parser.add_argument(name)
    args = parser.parse_args()
    print(json.dumps(run(args.catalog, args.initial_state_id, args.runtime,
                         args.waveform, args.directory), indent=2))


if __name__ == '__main__':
    main()
