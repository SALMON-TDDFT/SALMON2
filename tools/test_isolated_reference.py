"""Generate a small isolated C2H2 midpoint reference, not a convergence test.

Explicit executable and new external output directory required. Uses one process
and one thread; never downloads, builds, or enables the interactive worker.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def section(text, name, body):
    pattern = rf"(?ms)^&{name}\s*\n.*?^/\s*$"
    result, count = re.subn(pattern, f"&{name}\n{body}\n/", text)
    if count != 1:
        raise ValueError(f"expected one {name} section, found {count}")
    return result


def run(executable, directory):
    executable = Path(executable).resolve(strict=True)
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    samples = Path(__file__).resolve().parents[1] / 'runtime/sources/solver/samples'
    source = samples / 'exercise_01_C2H2_gs'
    gs = root / 'gs'
    gs.mkdir()
    for name in ('C_rps.dat', 'H_rps.dat'):
        shutil.copy2(source / name, gs / name)
    text = (source / 'C2H2_gs.inp').read_text()
    text = section(text, 'rgrid', '  dl(1:3) = 0.35d0, 0.35d0, 0.35d0\n  num_rgrid(1:3) = 24,24,32')
    text = section(text, 'analysis', "  yn_out_psi='n'\n  yn_out_dns='n'")
    text = section(text, 'scf', '  nscf=300\n  threshold=1.d-9')
    env = dict(os.environ)
    for key in tuple(env):
        if key.startswith('SALMONPY_'):
            env.pop(key)
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    started = time.monotonic()

    def execute(folder, contents):
        (folder / 'input.inp').write_text(contents)
        with (folder / 'input.inp').open('rb') as inp, (folder / 'run.log').open('wb') as out:
            subprocess.run([str(executable)], cwd=folder, stdin=inp, stdout=out,
                           stderr=subprocess.STDOUT, env=env, timeout=180, check=True)

    execute(gs, text)
    if '#GS converged at' not in (gs / 'run.log').read_text():
        raise RuntimeError('SCF did not converge; inspect gs/run.log')
    restart = gs / 'data_for_restart'
    if not restart.is_dir():
        raise RuntimeError('ground-state restart missing; inspect gs/run.log')
    rt_text = (samples / 'exercise_03_C2H2_rt/C2H2_rt_pulse.inp').read_text()
    rt_text = section(rt_text, 'rgrid', '  dl(1:3)=0.35d0,0.35d0,0.35d0\n  num_rgrid(1:3)=24,24,32')
    rt_text = section(rt_text, 'tgrid', '  dt=0.001d0\n  nt=20')
    rt_text = section(rt_text, 'analysis', '  de=0.1d0\n  nenergy=10')
    rt_text += "\n&propagation\n propagator='middlepoint'\n/\n"
    runs = []
    for label in ('first', 'repeat'):
        folder = root / label
        folder.mkdir()
        for name in ('C_rps.dat', 'H_rps.dat'):
            shutil.copy2(source / name, folder / name)
        shutil.copytree(restart, folder / 'restart')
        execute(folder, rt_text)
        runs.append(folder)
    outputs = ['C2H2_rt.data', 'C2H2_rt_energy.data']
    if not outputs:
        raise RuntimeError('no rt.data output; inspect logs')
    errors = {}
    for name in outputs:
        first = np.loadtxt(runs[0] / name)
        second = np.loadtxt(runs[1] / name)
        if not np.isfinite(first).all():
            raise AssertionError(f'nonfinite output: {name}')
        np.testing.assert_array_equal(first, second)
        errors[name] = 0.0
    data = np.loadtxt(runs[0] / 'C2H2_rt.data')
    if data.shape != (20, 19):
        raise AssertionError(f'unexpected isolated rt layout: {data.shape}')
    np.testing.assert_allclose(data[:, 0], np.arange(1, 21)*0.001, rtol=0, atol=1.e-12)
    # Isolated writer's E_ext columns are zero in this runtime; E_tot retains
    # the field used by calcVbox. Do not use the header alone to select a port.
    forward_E = -np.diff(data[:, 1:4], axis=0)/0.001
    np.testing.assert_allclose(data[:-1, 10:13], forward_E, rtol=1.e-6, atol=1.e-12)
    if not np.any(data[:, 10:13]) or not np.any(data[:, 13:16]):
        raise AssertionError('reference must have nonzero E and dipole response')
    report = dict(scope='isolated midpoint batch reference; no E worker or physical convergence claim',
                  executable_sha256=digest(executable),
                  pseudopotentials={name: digest(source / name) for name in ('C_rps.dat', 'H_rps.dat')},
                  input_sha256={label: digest(root / label / 'input.inp') for label in ('gs', 'first', 'repeat')},
                  repeat_max_abs_error=errors, scf_converged=True,
                  forward_difference_E_max_abs_error=float(np.max(np.abs(data[:-1, 10:13]-forward_E))),
                  native_E_ext_columns_zero=bool(np.all(data[:, 4:7] == 0)),
                  elapsed_seconds=time.monotonic()-started,
                  threads=1, grid=[24,24,32], steps=20)
    (root / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('directory')
    args = parser.parse_args()
    run(args.executable, args.directory)
