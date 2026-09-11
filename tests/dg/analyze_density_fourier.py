"""Read-only Fourier diagnostic for a uniform cubic density export.

Bands describe density differences, not the orbital PW cutoff. Density Fourier
components cannot be identified one-to-one with missing orbital basis modes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from density_error_decomposition import read_snapshot_set


def spectrum(delta, length):
    delta = np.asarray(delta, dtype=float)
    if (delta.ndim != 3 or len(set(delta.shape)) != 1 or delta.shape[0] < 2 or
            not np.isfinite(delta).all() or not np.isfinite(length) or length <= 0):
        raise ValueError('finite cubic grid and positive cell length required')
    n = delta.shape[0]
    indices = np.rint(np.fft.fftfreq(n)*n).astype(int)
    modes = np.meshgrid(indices, indices, indices, indexing='ij')
    shell = sum(component**2 for component in modes)
    amplitudes = np.fft.fftn(delta)/delta.size
    power = abs(amplitudes)**2
    total = float(power.sum())
    real_total = float(np.mean(delta**2))
    def fraction(mask):
        return float(power[mask].sum()/total) if total else 0.0
    shell_power = np.bincount(shell.ravel(), weights=power.ravel())
    order = np.argsort(shell_power)[::-1][:12]
    return {
        'mean_square': real_total,
        'parseval_relative_defect': abs(total-real_total)/max(total, np.finfo(float).tiny),
        'zero_fraction': fraction(shell == 0),
        'long_fraction': fraction((shell > 0) & (shell <= 4)),
        'middle_fraction': fraction((shell > 4) & (shell <= 16)),
        'short_fraction': fraction(shell > 16),
        'strongest_shells': [
            {'n_squared': int(s), 'wavelength_bohr': float(length/np.sqrt(s)) if s else None,
             'power_fraction': float(shell_power[s]/total) if total else 0.0}
            for s in order if shell_power[s] > 0],
    }


def analyze(prefix, n, length):
    evidence = read_snapshot_set(prefix)
    rows = np.empty((n**3, 6))
    occupied = np.zeros(n**3, dtype=bool)
    fragment_counts = []
    for name, digest in evidence['sha256'].items():
        path = Path(prefix).parent / name
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('snapshot changed after validation')
        lines = raw.decode('ascii').splitlines()
        fragment = int(lines[1].split()[2])
        data = np.array([[float(v) for v in line.split()] for line in lines[4:-1]])
        ids = data[:, 0].astype(int)-1
        if np.any(ids < 0) or np.any(ids >= n**3):
            raise ValueError('snapshot grid differs from supplied cubic geometry')
        rows[ids] = data[:, 1:]
        occupied[ids] = True
        fragment_counts.append({'fragment': fragment, **{
            label: float(np.dot(data[:, 1], data[:, column]))
            for label, column in (('dc', 2), ('conventional', 3), ('dg_frozen', 4), ('dg_final', 5))}})
    if not occupied.all() or not np.allclose(rows[:, 0], (length/n)**3, rtol=1e-12, atol=0):
        raise ValueError('grid coverage or volume weights disagree with cubic geometry')
    dc, conventional, frozen, final = (rows[:, i].reshape((n, n, n), order='F') for i in range(1, 5))
    differences = {'dc_lcfo_correction': conventional-dc,
                   'dg_increment': frozen-conventional,
                   'frozen_total_displacement': frozen-dc,
                   'relaxation_displacement': final-frozen}
    return {
        'grid': [n]*3, 'cell_length_bohr': length,
        'bands': {'long': 'L/2 <= wavelength <= L',
                  'middle': 'L/4 <= wavelength < L/2', 'short': 'wavelength < L/4'},
        'spectra': {name: spectrum(delta, length) for name, delta in differences.items()},
        'fragment_electron_counts': fragment_counts,
        'source_sha256': evidence['sha256'],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--grid', required=True, type=int, help='points per cubic cell axis')
    parser.add_argument('--length', required=True, type=float, help='cubic cell length in Bohr')
    args = parser.parse_args()
    print(json.dumps(analyze(args.prefix, args.grid, args.length), indent=2, allow_nan=False))
