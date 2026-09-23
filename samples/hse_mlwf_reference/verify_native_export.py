"""Validate a completed native export, its file lengths, and independent H action.

Usage: python verify_native_export.py /path/to/export
The rho field is SCF mixed density used by exported potentials, not necessarily
exactly the density of exported orbitals. This checker does not equate them.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from model import NativeModel


def verify(directory):
    path = Path(directory)
    marker = path / 'complete.txt'
    if not marker.exists() or marker.read_text().strip() != 'SALMON_HSE_REFERENCE_V1_COMPLETE':
        raise ValueError('Native export is incomplete: missing or invalid complete.txt')
    lines = (path / 'metadata.txt').read_text().splitlines()
    metadata = {c[0]: c[1:] for line in lines[1:] if (c := line.split()) and not c[0].startswith('#')}
    if metadata.get('endian_little') != ['T' if sys.byteorder == 'little' else 'F']:
        raise ValueError('Native export endianness does not match this reader')
    if metadata.get('rho_semantics') != ['scf_mixed_potential_density']:
        raise ValueError('Unspecified native density semantics')
    shape = tuple(map(int, metadata['grid']))
    ng = int(np.prod(shape))
    no, nk, nlma, nps, ni = (int(metadata[k][0]) for k in ('no', 'nk', 'nlma', 'nps', 'nion'))
    # Each size is in bytes, checked before any numerical reader is invoked.
    sizes = {'psi': 16*ng*no*nk, 'hpsi': 16*ng*no*nk,
             'projectors': 16*ng*nlma*nk, 'rinv_uvu': 8*nlma,
             'k': 8*3*nk, 'weights': 8*nk, 'occupations': 8*no*nk,
             'eigenvalues': 8*no*nk, 'geometry': 8*(21+3*ni),
             'coordinates': 8*sum(shape), 'stencil': 8*25, 'energies': 8*7,
             'raw_projectors': 8*nps*nlma, 'projector_positions': 8*3*nps*nlma,
             'projector_indices': 4*3*nps*nlma, 'projector_counts': 4*nlma}
    sizes.update({name: 8*ng for name in ('rho', 'vlocal', 'vh', 'vxc', 'vpsl')})
    if metadata.get('nlcc_available') == ['T']:
        sizes['rho_nlcc'] = 8*ng
    for name, expected in sizes.items():
        actual = (path / (name+'.bin')).stat().st_size
        if actual != expected:
            raise ValueError(f'{name}.bin: {actual} bytes, expected {expected}')
    model = NativeModel(path)
    uv = np.fromfile(path/'raw_projectors.bin').reshape((nps, nlma), order='F')
    xyz = np.fromfile(path/'projector_positions.bin').reshape((3, nps, nlma), order='F')
    idx = np.fromfile(path/'projector_indices.bin', dtype=np.int32).reshape((3, nps, nlma), order='F')-1
    counts = np.fromfile(path/'projector_counts.bin', dtype=np.int32)
    if not np.isfinite(uv).all() or not np.isfinite(xyz).all():
        raise ValueError('Nonfinite raw projector data')
    if np.any(counts < 0) or np.any(counts > nps):
        raise ValueError('Invalid raw projector counts')
    exported = model.read('projectors', shape+(nlma,nk), complex).transpose(4,3,0,1,2).reshape(nk,nlma,-1)
    error = 0.
    for ik, k in enumerate(model.k):
        dense = np.zeros((nlma, *shape), complex)
        for channel in range(nlma):
            n = counts[channel]
            indices = idx[:, :n, channel]
            if np.any(indices < 0) or np.any(indices >= np.array(shape)[:, None]):
                raise ValueError('Raw projector indices outside physical grid')
            values = uv[:n, channel]*np.exp(-1j*k@xyz[:, :n, channel])
            np.add.at(dense[channel], tuple(indices), values)
        error = max(error, float(np.max(abs(dense.reshape(nlma, -1)-exported[ik]))))
    if not np.isfinite(error) or error > 1e-12:
        raise AssertionError(f'Raw projector reconstruction error {error}')
    result = model.native_parity()
    for key in ('hpsi_relative_error', 'hartree_potential_max_error', 'hartree_energy_error_Ha'):
        if not np.isfinite(result[key]) or abs(result[key]) > 1e-10:
            raise AssertionError(f'Native parity failed: {key}={result[key]}')
    result.update(raw_projector_max_error=error, checked_binary_files=len(sizes),
                  rho_semantics=metadata['rho_semantics'][0])
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    print(json.dumps(verify(parser.parse_args().directory), indent=2))
