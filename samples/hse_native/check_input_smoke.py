#!/usr/bin/env python3
"""Bounded native HSE input checks using an existing converged Si8 4^3 GS.

Example: python3 samples/hse_native/check_input_smoke.py --exe /path/to/salmon \
  --gs /path/to/data_for_restart --work /tmp/hse-input-smoke

The GS must match gs.inp (Si8, 12^3 real grid, shifted full 4^3 k mesh,
16 occupied orbitals, omega=.11 bohr^-1). No SCF is run. Custom-omega
trajectories test wiring only and are not physical production calculations.
Optional --off-exe checks a serial HSE=OFF PZ odd-step restart. It reuses the
HSE GS solely for a numerical restart comparison, not as a physical PZ GS.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import time

LENGTH = 0.52917721067
ENERGY = 27.21138505
TIME = 0.02418884326505


def quoted(path):
    return "'" + str(path).replace("'", "''") + "'"


def make_input(gs, pseudo, ranks, *, afs=False, omega=None, propagation='', restart=False,
               functional='hse06', steps=None):
    length, energy, time_unit = (LENGTH, ENERGY, TIME) if afs else (1., 1., 1.)
    atoms = [(0, 0, 0), (1, 1, 1), (2, 0, 2), (0, 2, 2),
             (2, 2, 0), (3, 1, 3), (1, 3, 3), (3, 3, 1)]
    coordinates = '\n'.join("'Si' " + ' '.join(format(x * 2.565 * length, '.17g') for x in row) + ' 1'
                            for row in atoms)
    omega_line = '' if omega is None else f'hse_omega={omega}'
    return f"""&calculation
 theory='tddft_response'
/
&control
 sysname='Si'
 directory_read_data={quoted(str(gs) + '/')}
 yn_restart='{'y' if restart else 'n'}'
 checkpoint_interval=1
/
&units
 unit_system='{'A_eV_fs' if afs else 'a.u.'}'
/
&system
 yn_periodic='y'
 al={10.26 * length:.17g},{10.26 * length:.17g},{10.26 * length:.17g}
 nstate=16
 nelec=32
 nelem=1
 natom=8
/
&pseudo
 izatom(1)=14
 file_pseudo(1)={quoted(pseudo)}
 lmax_ps(1)=2
 lloc_ps(1)=2
/
&functional
 xc='{functional}'
 {omega_line}
/
&rgrid
 dl={.855 * length:.17g},{.855 * length:.17g},{.855 * length:.17g}
/
&kgrid
 num_kgrid=4,4,4
/
&tgrid
 nt={steps if steps is not None else (2 if restart else 1)}
 dt={.16 * time_unit:.17g}
/
&scf
 ncg=4
 nscf=1
 threshold=1d-9
/
&emfield
 trans_longi='tr'
 e_impulse={1e-4 * energy / length * time_unit:.17g}
 ae_shape1='impulse'
 epdir_re1=0,0,1
/
&analysis
 out_rt_energy_step=1
 nenergy=8
 de={.001 * energy:.17g}
/
&atomic_coor
{coordinates}
/
&parallel
 nproc_k={ranks}
 nproc_ob=1
 nproc_rgrid=1,1,1
/
{propagation}
"""


def numeric_rows(path):
    rows = [[float(x.replace('D', 'E')) for x in line.split()]
            for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
    if not rows or not all(math.isfinite(x) for row in rows for x in row):
        raise AssertionError(f'Missing/nonfinite data: {path}')
    return rows


def energy_difference(left, right, right_scale=1.):
    a, b = numeric_rows(left / 'Si_rt_energy.data'), numeric_rows(right / 'Si_rt_energy.data')
    if len(a) != len(b):
        raise AssertionError('Energy output row counts differ')
    return max(abs(x[1] - y[1] / right_scale) for x, y in zip(a, b))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', type=Path, required=True)
    p.add_argument('--off-exe', type=Path, help='Optional serial USE_HSE=OFF binary for PZ restart checks')
    p.add_argument('--gs', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True, help='New, nonexisting output directory')
    p.add_argument('--ranks', type=int, default=8)
    p.add_argument('--mpiexec', default='mpiexec')
    p.add_argument('--timeout', type=float, default=180, help='Seconds per case, including MPI startup')
    args = p.parse_args()
    args.exe, args.gs, args.work = (x.resolve() for x in (args.exe, args.gs, args.work))
    pseudo = Path(__file__).resolve().parents[2] / 'testsuites/pseudo/Si_rps.dat'
    if not args.exe.is_file() or not args.gs.is_dir() or not pseudo.is_file():
        p.error('Executable, GS directory and repository Si pseudopotential must exist')
    if args.off_exe is not None:
        args.off_exe = args.off_exe.resolve()
        if not args.off_exe.is_file():
            p.error('--off-exe must exist')
    if args.ranks < 1 or args.ranks > 64 or 64 % args.ranks or args.timeout <= 0:
        p.error('ranks must divide 64; timeout must be positive')
    args.work.mkdir(parents=True, exist_ok=False)
    metrics = {'exe': str(args.exe), 'gs': str(args.gs), 'ranks': args.ranks, 'cases': {}, 'passed': False}
    if args.off_exe:
        metrics['off_exe'] = str(args.off_exe)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')

    def run(name, expected=None, off=False, stale_metadata=False, **options):
        folder = args.work / name
        folder.mkdir()
        ranks = 1 if off else args.ranks
        if stale_metadata:
            checkpoint = folder / 'checkpoint_rt_000001'
            checkpoint.mkdir()
            (checkpoint / 'hse_restart.bin').write_bytes(b'stale HSE metadata\n')
        text = make_input(options.pop('gs', args.gs), pseudo, ranks, **options)
        (folder / 'inputfile').write_text(text)
        started = time.monotonic()
        with (folder / 'inputfile').open() as inp, (folder / 'run.log').open('w') as log:
            command = [str(args.off_exe)] if off else [args.mpiexec, '-n', str(ranks), str(args.exe)]
            process = subprocess.Popen(command,
                                       stdin=inp, stdout=log, stderr=subprocess.STDOUT,
                                       cwd=folder, env=env, start_new_session=True)
            try:
                code = process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                raise AssertionError(f'{name}: per-case timeout reached')
        output = (folder / 'run.log').read_text()
        okay = (code != 0 and expected in output) if expected else code == 0
        metrics['cases'][name] = {'exit_code': code, 'seconds': time.monotonic() - started, 'passed': okay}
        if not okay:
            raise AssertionError(f'{name}: unexpected outcome; inspect {folder / "run.log"}')
        return folder

    try:
        base = run('omitted')
        explicit = run('explicit', omega='.11d0', propagation="&propagation\n propagator='hse_taylor4'\n n_hamil=4\n yn_predictor_corrector='y'\n/")
        for filename in ('Si_rt.data', 'Si_rt_energy.data'):
            if numeric_rows(base / filename) != numeric_rows(explicit / filename):
                raise AssertionError(f'Omitted/explicit defaults differ: {filename}')
        for name, omega in [('afs_default', None), ('afs_explicit', format(.11 / LENGTH, '.17g'))]:
            folder = run(name, afs=True, omega=omega)
            delta = energy_difference(base, folder, ENERGY)
            metrics['cases'][name]['max_energy_difference_ha'] = delta
            if delta > 1e-10:
                raise AssertionError(f'{name}: atomic-unit parity failed ({delta})')
        custom = run('afs_custom', afs=True, omega='.2d0')
        logged = re.search(r'hse_omega\s*\(bohr\^-1\)\s*=\s*([\d.EeDd+\-]+)', (custom / 'variables.log').read_text())
        if logged is None or abs(float(logged[1].replace('D', 'E')) - .2 * LENGTH) > 1e-13:
            raise AssertionError('Custom inverse-Angstrom omega was not converted correctly')
        if energy_difference(base, custom, ENERGY) < 1e-8:
            raise AssertionError('Custom omega did not change the computed energy')
        run('invalid_pc', expected='predictor-corrector', propagation="&propagation\n yn_predictor_corrector='n'\n/")
        for name, omega in [('zero_omega', '0d0'), ('negative_omega', '-.1d0')]:
            run(name, expected='hse_omega must be finite and positive', omega=omega)
        checkpoint = base / 'checkpoint_rt_000001'
        if not (checkpoint / 'hse_restart.bin').is_file():
            raise AssertionError('One-step run did not create HSE RT metadata')
        run('restart_same', gs=checkpoint, restart=True)
        run('restart_changed_omega', gs=checkpoint, restart=True, omega='.2d0', expected='restart physics')
        if args.off_exe:
            fresh = run('pz_fresh_two', off=True, functional='pz', steps=2)
            first = run('pz_first_one', off=True, functional='pz', stale_metadata=True)
            checkpoint = first / 'checkpoint_rt_000001'
            if (checkpoint / 'hse_restart.bin').exists():
                raise AssertionError('PZ checkpoint retained stale HSE metadata')
            resumed = run('pz_odd_restart', off=True, functional='pz', gs=checkpoint, restart=True)
            a, b = (numeric_rows(x / 'Si_rt.data')[-1] for x in (fresh, resumed))
            if a[0] != b[0]:
                raise AssertionError('PZ restart endpoint times differ')
            current_error = max(abs(x - y) for x, y in zip(a[13:16], b[13:16]))
            ea, eb = (numeric_rows(x / 'Si_rt_energy.data')[-1] for x in (fresh, resumed))
            energy_error = abs(ea[1] - eb[1])
            metrics['cases']['pz_odd_restart'].update(max_current_difference_au=current_error,
                                                     energy_difference_ha=energy_error,
                                                     stale_hse_metadata_removed=True)
            if current_error > 1e-10 or energy_error > 1e-10:
                raise AssertionError(f'PZ odd restart differs: J={current_error}, E={energy_error}')
        metrics['passed'] = True
    finally:
        (args.work / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
