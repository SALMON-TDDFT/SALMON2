"""Synthetic symmetry-preserving TDDFT fixture, not a converged material.

No solver execution/download/build. Used by the explicit internal fast-PR probe.
Usage: python tools/prepare_fast_hse_symmetry_fixture.py PSEUDO NEW_DIR --group c2z
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pseudo', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--group', choices=('c2z', 'mirror-x', 'full'), default='c2z')
    parser.add_argument('--propagator', choices=('hse_taylor4', 'hse_taylor4_full'), default='hse_taylor4')
    args = parser.parse_args()
    pseudo = args.pseudo.resolve(strict=True)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    restart = out / 'data_for_restart'
    restart.mkdir()
    nk = 8 if args.group == 'full' else 4

    def record(stream, payload):
        stream.write(struct.pack('=i', len(payload)))
        stream.write(payload)
        stream.write(struct.pack('=i', len(payload)))

    with (restart / 'info.bin').open('wb') as stream:
        for n in (nk, 2, 0, 1, 0):
            record(stream, struct.pack('=i', n))
    with (restart / 'occupation.bin').open('wb') as stream:
        record(stream, np.full((2, nk), 2., dtype='=f8').tobytes(order='F'))
    z = np.arange(8)[None, None, :]
    psi = np.empty((8, 8, 8, 2), dtype='=c16', order='F')
    psi[..., 0] = 1 / np.sqrt(512.)
    psi[..., 1] = np.sqrt(2 / 512.) * np.cos(2 * np.pi * z / 8)
    flat = psi.reshape(512, 2, order='F')
    np.testing.assert_allclose(flat.conj().T @ flat, np.eye(2), atol=1e-14)
    (restart / 'wfn.bin').write_bytes(np.repeat(psi[..., None], nk, axis=4).tobytes(order='F'))
    shutil.copy2(pseudo, out / 'pseudo1_rps.dat')
    operations = np.zeros((2, 3, 4))
    operations[:, :, :3] = np.eye(3)
    operations[1, 0, 0] = -1
    if args.group == 'c2z':
        operations[1, 1, 1] = -1
    np.savetxt(out / 'sym.dat', operations.reshape(6, 4), fmt='%.1f')
    text = f"""&calculation
 theory='tddft_response'
/
&control
 sysname='rt',directory_read_data='./data_for_restart/',yn_restart='n'
/
&units
 unit_system='au'
/
&system
 yn_periodic='y',al=8,8,8,nelem=1,natom=1,nelec=4,nstate=2
 yn_symmetry='{'n' if args.group == 'full' else 'yyn'}'
/
&pseudo
 file_pseudo(1)='pseudo1_rps.dat',izatom(1)=14,lloc_ps(1)=2
/
&functional
 xc='hse06',exx_factor_mode='off'
/
&rgrid
 num_rgrid=8,8,8
/
&kgrid
 num_kgrid=2,2,2
/
&tgrid
 dt=.001,nt=12
/
&propagation
 propagator='{args.propagator}',yn_predictor_corrector='y',n_hamil=4
/
&emfield
 ae_shape1='impulse',e_impulse=0,epdir_re1=0,0,1,trans_longi='tr'
/
&analysis
 out_rt_energy_step=1
/
&parallel
 nproc_k=1,nproc_ob=1,nproc_rgrid=1,1,1
/
&atomic_red_coor
 'Si' 0 0 0 1
/
"""
    (out / 'rt.inp').write_text(text)
    (out / 'fixture.json').write_text(json.dumps({
        'scope': 'synthetic symmetry-preserving orbitals; not converged GS',
        'group': args.group, 'grid': [8, 8, 8], 'full_kpoints': 8, 'representatives': nk,
        'files': {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in out.rglob('*') if p.is_file()},
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
