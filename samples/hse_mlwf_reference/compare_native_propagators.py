"""Gauge-independent native RT comparison at a common physical endpoint."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from model import NativeModel


def gauge_invariant_wavefunction_error(u, reference):
    """Relative distance after the best occupied unitary rotation at each k."""
    distance = 0.0
    for k in range(u.shape[2]):
        left, _, right = np.linalg.svd(u[:, :, k].conj().T @ reference[:, :, k])
        aligned = u[:, :, k] @ (left @ right)
        distance += np.linalg.norm(aligned - reference[:, :, k])**2
    return float(np.sqrt(distance) / np.linalg.norm(reference))


def load_run(path, model):
    path = Path(path)
    log = (path / 'run.log').read_text()
    if 'end SALMON' not in log:
        raise ValueError(f'Incomplete run: {path}')
    current = np.atleast_2d(np.loadtxt(path / 'Si_rt.data'))
    energy = np.atleast_2d(np.loadtxt(path / 'Si_rt_energy.data'))
    checkpoint = sorted(path.glob('checkpoint_rt_*/wfn.bin'))[-1]
    dt = float(re.search(r'^\s*dt\s*=\s*([\d.eEdD+-]+)',
                         (path / 'inputfile').read_text(), re.M).group(1).replace('d', 'e').replace('D', 'E'))
    checkpoint_step = int(checkpoint.parent.name.rsplit('_', 1)[1])
    if abs(checkpoint_step * dt - current[-1, 0]) > 1e-8:
        raise ValueError(f'Last checkpoint is not the current endpoint: {path}')
    u = np.fromfile(checkpoint, np.complex128).reshape(
        (int(np.prod(model.shape)), model.no, model.nk), order='F')
    rho = 2 * np.sum(abs(u)**2, axis=(1, 2)) / model.nk
    gram = np.einsum('gik,gjk->kij', u.conj(), u) * model.dv
    return dict(path=path, current=current, energy=energy, rho=rho,
                gram=gram, u=u)


def compare(export, reference, directories):
    model = NativeModel(export)
    ref = load_run(reference, model)
    rho0 = model.density(model.psi).reshape(-1, order='F')
    response_scale = np.linalg.norm(ref['rho'] - rho0)
    runs = [load_run(path, model) for path in directories]
    common_times = ref['current'][:, 0]
    for run in runs:
        common_times = np.intersect1d(common_times, run['current'][:, 0])
    if not len(common_times):
        raise ValueError('No common current samples')
    jref = ref['current'][np.searchsorted(ref['current'][:, 0], common_times), 13:16]
    results = {}
    for run in runs:
        t = run['current'][:, 0]
        if abs(t[-1] - ref['current'][-1, 0]) > 1e-8:
            raise ValueError('Compare runs at the same physical final time')
        # Use matching time samples, never interpolate away time-step errors.
        dj = run['current'][np.searchsorted(t, common_times), 13:16] - jref
        drho = run['rho'] - ref['rho']
        e = run['energy'][run['energy'][:, 0] > 0, 1]
        results[run['path'].name] = dict(
            final_time_au=float(t[-1]), common_current_samples=len(common_times),
            current_max_abs_error_au=float(np.max(abs(dj))),
            current_relative_rms_error=float(np.linalg.norm(dj) / np.linalg.norm(jref)),
            current_max_error_over_reference_peak=float(np.max(abs(dj)) / np.max(abs(jref))),
            density_relative_error=float(np.linalg.norm(drho) / np.linalg.norm(ref['rho'])),
            density_error_over_induced_density=float(np.linalg.norm(drho) / response_scale),
            gauge_aligned_wavefunction_relative_error=gauge_invariant_wavefunction_error(run['u'], ref['u']),
            electron_number=float(2 * np.trace(run['gram'], axis1=1, axis2=2).real.mean()),
            gram_max_error=float(np.max(abs(run['gram'] - np.eye(model.no)))),
            post_impulse_energy_range_Ha=float(np.ptp(e)),
            final_energy_Ha=float(e[-1]))
    return dict(reference=str(reference), results=results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('directories', nargs='+')
    args = parser.parse_args()
    result = compare(args.export, args.reference, args.directories)
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
