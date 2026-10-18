#!/usr/bin/env python3
"""Reproducible local small-cell HSE06 SOC/NLCC integration validation.

This is a developer verification, not converged Si or WSe2 optical data.
The output directory must not exist. The production binary is never rebuilt.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np

HERE = Path(__file__).resolve().parent

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--salmon', type=Path, required=True)
    ap.add_argument('--pseudo', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--threads', type=int, default=1)
    args = ap.parse_args()
    binary, pseudo, output = args.salmon.resolve(), args.pseudo.resolve(), args.output.resolve()
    expected = json.loads((HERE / 'pseudo.json').read_text())
    digest = hashlib.sha256(pseudo.read_bytes()).hexdigest()
    if digest != expected['sha256']:
        raise ValueError('Expected the pinned original Si_PBE19.vps; see pseudo.json')
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(pseudo, output / 'Si_PBE19.vps')
    if args.threads < 1:
        raise ValueError('threads must be positive')
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1')
    report = {'threads': args.threads, 'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
              'pseudo_sha256': digest, 'scope': 'coarse two-atom periodic developer fixture'}

    def run(name, text):
        folder = output / name
        folder.mkdir()
        (folder / 'inputfile').write_text(text)
        with (folder / 'inputfile').open() as inp, (folder / 'run.log').open('w') as log:
            rc = subprocess.run([str(binary)], cwd=folder, stdin=inp, stdout=log,
                                stderr=subprocess.STDOUT, env=env, timeout=1800).returncode
        log = (folder / 'run.log').read_text()
        if rc or 'end SALMON' not in log or 'ERROR STOP' in log:
            raise RuntimeError(f'{name} failed: inspect {folder / "run.log"}')
        return folder, log

    gs = (HERE / 'gs.inp').read_text()
    scalar = gs.replace("spin='noncollinear',yn_spinorbit='y'", "spin='unpolarized',yn_spinorbit='n'")
    scalar = scalar.replace('nstate=8', 'nstate=4').replace("method_mixing='simple_dm'", "method_mixing='simple'")
    for name, text in [('scalar-gs', scalar), ('gs', gs),
                       ('gs-pre', gs.replace('exx_pre_scf_threshold=0d0', 'exx_pre_scf_threshold=1d-4'))]:
        folder, log = run(name, text)
        residuals = re.findall(r'iter and int_x.*?=\s*\d+\s+([0-9.Ee+-]+)', log)
        if not residuals or float(residuals[-1]) >= 1e-7 or '#GS converged at' not in log:
            raise AssertionError(f'{name}: unconverged GS must not be used for RT')
        electrons = re.findall(r'Ne=\s*([0-9.Ee+-]+)', log)
        if not electrons or abs(float(electrons[-1]) - 8.) > 1e-10:
            raise AssertionError(f'{name}: incorrect electron count')
        if name == 'gs-pre' and 'EXX_PRE_SCF switch to target hybrid' not in log:
            raise AssertionError('PBE warmup must switch to HSE before accepting the GS')
        report[name] = {'residual': float(residuals[-1]), 'electrons': float(electrons[-1])}
    base = (HERE / 'rt.inp').read_text()
    for mode in ('ace', 'direct', 'zero', 'halfdt'):
        text = base
        if mode == 'direct':
            text += "\n&propagation propagator='hse_taylor4_full' /\n"
        if mode == 'zero':
            text = text.replace('e_impulse=1d-4', 'e_impulse=0d0')
        if mode == 'halfdt':
            text = text.replace('dt=.02d0,nt=100', 'dt=.01d0,nt=200')
        run('rt-' + mode, text)
    run('rt-pre', base.replace('../gs/data_for_restart/', '../gs-pre/data_for_restart/'))
    scalar_rt = base.replace("spin='noncollinear',yn_spinorbit='y'", "spin='unpolarized',yn_spinorbit='n'")
    scalar_rt = scalar_rt.replace("method_mixing='simple_dm'", "method_mixing='simple'")
    scalar_rt = scalar_rt.replace('nstate=8', 'nstate=4').replace('../gs/data_for_restart/', '../scalar-gs/data_for_restart/')
    run('scalar-rt', scalar_rt)

    def data(name, suffix):
        a = np.loadtxt(output / name / ('Si2_' + suffix + '.data'), ndmin=2)
        if not np.isfinite(a).all():
            raise AssertionError(f'{name}: nonfinite observables')
        return a

    for name in ('rt-ace', 'rt-direct', 'rt-zero', 'rt-halfdt', 'scalar-rt'):
        current = data(name, 'rt')
        count = 200 if name == 'rt-halfdt' else 100
        assert len(current) == count and abs(current[-1, 0] - 2.) < 1e-12
        energy = data(name, 'rt_energy')
        report[name] = {'steps': count, 'final_time_au': float(current[-1, 0]),
                        'energy_span_ha': float(np.ptp(energy[:, 1]))}
    ace, direct, half = data('rt-ace', 'rt'), data('rt-direct', 'rt'), data('rt-halfdt', 'rt')[1::2]
    np.testing.assert_allclose(ace[:, 0], half[:, 0], atol=1e-12, rtol=0)
    report['ace_direct_current_max_au'] = float(np.max(abs(ace[:, 13:16] - direct[:, 13:16])))
    report['halfdt_current_max_au'] = float(np.max(abs(ace[:, 13:16] - half[:, 13:16])))
    assert report['ace_direct_current_max_au'] < 1e-9
    assert report['halfdt_current_max_au'] < 1e-8
    assert report['rt-zero']['energy_span_ha'] < 1e-8
    report['pre_scf_initial_energy_difference_ha'] = float(abs(data('rt-pre', 'rt_energy')[0, 1] -
                                                           data('rt-ace', 'rt_energy')[0, 1]))
    assert report['pre_scf_initial_energy_difference_ha'] < 1e-6
    # A finite-response refinement test catches the endpoint-XC bug that the
    # shared ACE/direct path and stationary zero-field tests cannot detect.
    refined = []
    for name, dt, nt in [('coarse', '.08', 50), ('medium', '.04', 100), ('fine', '.02', 200)]:
        text = base.replace('dt=.02d0,nt=100', f'dt={dt}d0,nt={nt}').replace('e_impulse=1d-4', 'e_impulse=1d-2')
        run('order-' + name, text)
        a = data('order-' + name, 'rt')
        assert len(a) == nt
        refined.append(a[:, 13:16])
    e1 = float(np.max(abs(refined[0] - refined[1][1::2])))
    e2 = float(np.max(abs(refined[1] - refined[2][1::2])))
    report['time_order'] = {'coarse_error': e1, 'fine_error': e2, 'ratio': e1 / e2}
    assert e1 / e2 > 3., 'SOC XC must be evaluated at the corrector midpoint'
    report['result'] = 'PASS'
    (output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
