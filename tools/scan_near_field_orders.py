"""Compare radial-complete saved E/H expansion orders on held-out spatial samples; no solver launch."""
import argparse
import ctypes as c
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
from salmonpy import _material_handles
from salmonpy.field_sampling import SampledField
from salmonpy.near_field import NearFieldSeries


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('samples', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--kernels', type=Path, required=True)
    p.add_argument('--center', type=float, nargs=3, required=True)
    p.add_argument('--radius', type=float, required=True)
    p.add_argument('--orders', type=int, nargs='+', default=[0, 1, 2, 3, 4, 6, 8])
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    with np.load(a.samples, allow_pickle=False) as data:
        times, positions, E, H = (data[k] for k in ('times', 'positions', 'E', 'H'))
        # Deterministic split; held-out points never enter the coefficient fit.
        held = np.arange(len(positions)) % 3 == 1
        rows = []
        with patch.object(_material_handles, '_library_path', return_value=a.kernels):
            api = _material_handles.default_sbe_api()
            compare = api.library.sp_convergence_waveforms
            compare.argtypes = [c.c_int, c.c_void_p, c.c_void_p, c.c_double, c.c_double, c.c_void_p, c.c_void_p]
            compare.restype = c.c_int
            for order in a.orders:
                row = {'lmax': order}
                try:
                    series = NearFieldSeries.from_samples(
                        SampledField('E', times, positions[~held], E[:, ~held], 'E:au'),
                        SampledField('H', times, positions[~held], H[:, ~held], 'H:au'),
                        center=a.center, radius=a.radius, lmax=order)
                    predicted = [series.fields(positions[held], step_index=n) for n in range(len(times))]
                    for index, (name, reference) in enumerate((('E', E), ('H', H))):
                        fit = np.ascontiguousarray([v[index] for v in predicted])
                        ref = np.ascontiguousarray(reference[:, held])
                        stats = np.empty(4)
                        accepted = c.c_int()
                        code = compare(ref.size, fit.ctypes.data, ref.ctypes.data, 0., .05, stats.ctypes.data, c.byref(accepted))
                        if code:
                            raise RuntimeError(f'Fortran residual comparison failed: {code}')
                        row[name] = dict(zip(('absolute_rms', 'relative_rms', 'reference_rms', 'maximum_absolute'), map(float, stats)))
                except (ValueError, RuntimeError) as exc:
                    row['error'] = str(exc)
                rows.append(row)
    report = dict(schema='salmonpy.near-field-order-scan/2', basis='r^(l+2k)Ylm; lmax is total degree', source=str(a.samples),
                  center=a.center, radius=a.radius, training_points=int((~held).sum()),
                  held_out_points=int(held.sum()), rows=rows, physical_convergence_certified=False)
    a.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
