#!/usr/bin/env python3
"""Opt-in physical-Si regression for independent LDA W in both SALMON RT modes.

Use existing converged TB-mBJ GS and native TB/LDA W files made with the same
crystal, pseudopotential and meshes. This Python test harness does not generate W.
The executable must already have passed source/application-object audits; its
recorded SHA256 is verified before execution. All generated files stay in scratch.
"""
import argparse
import hashlib
import os
from pathlib import Path
import subprocess

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('executable', 'tb-gs', 'tb-w', 'lda-w', 'pseudo', 'work'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--grid', type=int, default=12)
    parser.add_argument('--mesh', type=int, default=2)
    parser.add_argument('--ranks', type=int, default=4)
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    exe = args.executable.resolve()
    if hashlib.sha256(exe.read_bytes()).hexdigest() != args.sha256:
        raise AssertionError('Executable differs from its audited hash')
    if args.ranks < 1 or args.mesh**3 % args.ranks:
        raise ValueError('MPI ranks must divide the k mesh')
    args.work.mkdir(exist_ok=False)
    source = (Path(__file__).with_name('internal_si')/'rt.inp').read_text()
    source = source.replace('@RESTART@', str(args.tb_gs.resolve())+'/')
    source = source.replace('@UPDATE@', "yn_hsex_update_w='n'")
    source = source.replace('nproc_k=4', f'nproc_k={args.ranks}')
    source = source.replace('nstate = 8', 'nstate = 4')
    source = source.replace('8, 8, 8', ', '.join([str(args.grid)]*3))
    source = source.replace('num_kgrid = 2, 2, 2', 'num_kgrid = '+', '.join([str(args.mesh)]*3))
    source = source.replace('hsex_w_qgrid=2,2,2', 'hsex_w_qgrid='+','.join([str(args.mesh)]*3))
    source = source.replace('checkpoint_interval=1', 'checkpoint_interval=0')
    source = source.replace('dt = 0.01', 'dt = 0.05').replace('e_impulse=0.00001', 'e_impulse=0.0001')
    source = source.replace("projection_option='td', out_projection_step=1, yn_out_intraband_current='y'", '')
    start = source.index(' ! Delayed verification probe:')
    end = source.index('\n/', start)
    source = source[:start]+source[end:]
    source = source.replace("kli_mode='hsex'", "kli_mode='hsex', hsex_screening_w='"+str(args.lda_w.resolve())+"'")
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1')
    currents = []
    for theory in ('tddft_response', 'tddft_pulse'):
        directory = args.work/theory
        directory.mkdir()
        (directory/'input').write_text(source.replace("theory='tddft_pulse'", "theory='"+theory+"'"))
        (directory/'Si_rps.dat').symlink_to(args.pseudo.resolve())
        (directory/'hsex-w.bin').symlink_to(args.tb_w.resolve())
        with (directory/'input').open('rb') as data, (directory/'run.log').open('wb') as log:
            subprocess.run(['mpiexec', '-n', str(args.ranks), str(exe)], cwd=directory,
                           env=env, stdin=data, stdout=log, stderr=subprocess.STDOUT, check=True)
        text = (directory/'run.log').read_text()
        if 'end SALMON' not in text or 'HSEX_INDEPENDENT_LDA_W' not in text:
            raise AssertionError('Independent-W RT did not complete: '+theory)
        a = np.loadtxt(directory/'Si_rt.data')
        assert a.shape == (3, 16) and np.isfinite(a).all()
        np.testing.assert_allclose(a[:, 0], [.05, .10, .15], rtol=0, atol=1e-10)
        currents.append(a[:, 13:16])
        if theory == 'tddft_response':
            spectrum = np.loadtxt(directory/'Si_response.data')
            assert spectrum.shape == (1000, 13) and np.isfinite(spectrum).all()
    np.testing.assert_allclose(currents[0], currents[1], rtol=0, atol=1e-12)
    print('INDEPENDENT_LDA_RESPONSE_PULSE_CURRENT_MATCH_PASS max_error=',
          float(np.max(abs(currents[0]-currents[1]))))


if __name__ == '__main__':
    main()
