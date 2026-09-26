"""Independently recover the global DC chemical potential from current Ritz data.
Fixture: Si64 8x1x1, core 16^3, fragment 32x16x16, T=300 K, 256 electrons.
"""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'samples/dc_hse'))
from locality import read_eigen_pair
root=Path(sys.argv[1]);energies=[];weights=[];actual=[]
for i in range(1,9):
    p=root/f'data_dcdft/fragments/{i:06d}'
    psi,hp,dv,_=read_eigen_pair(p/'hse_eigen_after_orthogonalization.bin')
    energies.append(np.sum(psi.conj()*hp,axis=0).real*dv)
    weights.append(np.sum(abs(psi.reshape(32,16,16,64,order='F')[:16])**2,axis=(0,1,2))*dv)
    actual.append(read_eigen_pair(p/'hse_eigen_pair.bin')[3])
e=np.array(energies);w=np.array(weights);actual=np.array(actual)
kt=300*3.166811563e-6
lo=e.min()-1;hi=e.max()+1
for _ in range(100):
    mu=(lo+hi)/2
    occ=2/(1+np.exp(np.clip((e-mu)/kt,-700,700)))
    if np.sum(occ*w)>256:hi=mu
    else:lo=mu
error=float(np.max(abs(actual-occ)))
print('max occupation error against current Ritz spectrum:',error)
assert error < 1e-5, error
assert abs(np.sum(actual*w)-256) < 1e-7
print('PASS: occupations and global electron count agree with the current fixed-H Ritz spectrum')
