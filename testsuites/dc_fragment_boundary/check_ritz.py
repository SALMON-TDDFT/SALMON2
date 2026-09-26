"""Require returned fixed-H orbitals to be a sorted orthonormal Ritz basis."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'samples/dc_hse'))
from locality import read_eigen_pair
for frag in range(1,9):
    p=Path(sys.argv[1])/f'data_dcdft/fragments/{frag:06d}/hse_eigen_after_orthogonalization.bin'
    psi,hp,dv,occ=read_eigen_pair(p)
    h=psi.conj().T@hp*dv
    off=h-np.diag(np.diag(h))
    assert np.max(abs(off)) < 1e-10, (frag,'non-diagonal projected Hamiltonian',np.max(abs(off)))
    assert np.min(np.diff(np.diag(h).real)) > -1e-10, 'unsorted Ritz states'
    assert np.max(abs(psi.conj().T@psi*dv-np.eye(len(occ)))) < 1e-10
print('PASS: returned fixed-H states are orthonormal, diagonal and sorted')
