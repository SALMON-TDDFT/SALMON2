"""Dense dyadic Hadamard Hamiltonian with an exact spectral oracle (2048 grid points)."""
import os
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import hadamard
from scipy.linalg.lapack import zhetrd,zungqr
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from q_openmp import QKernel
n=2048
a=hadamard(n).astype(float)
levels=np.random.default_rng(27).permutation(np.arange(1,n+1))/16
# All entries in H and the RHS are exactly representable dyadic rationals.
h=np.asarray((a*levels)@a.T/n,complex)
coeff=np.array([[1.,.25+.5j,0.],[.25,0.,1.],[0.,.125, .5j]])
b=(a[:,[0,2,5]]/32)@coeff
f,d,e,tau,info=zhetrd(h,lower=1,lwork=64*n)
assert info==0
rot=np.eye(n,dtype=complex)
rot[1:,1:],_,info=zungqr(f[1:,:-1],tau,lwork=64*n)
assert info==0
kernel=QKernel(os.environ['HSEX_Q_KERNEL'])
for eta in (5e-4,1e-6,1e-9,1e-13):
    delta=levels[[0,2,5]]-levels[0]
    weight=delta/(delta**2+eta**2)
    expected=(n/1024)*(coeff.conj().T@(weight[:,None]*coeff))
    try:
        got,res,team=kernel.solve(rot.conj().T,d,e,np.empty((n,0)),np.ones((n,1,1)),b[:,:,None],
            [3],np.array([[levels[0]]]),np.array([[2.]]),[eta],12,True,original_h=h)
    except ValueError as exc:
        if eta!=1e-13:
            raise
        print('LARGE_ORIGINAL_H_EXPLICIT_REJECTION',n,eta,str(exc),flush=True)
        continue
    ratio=np.max(abs(got[0][0]-expected)/(2e-12+2e-12*abs(expected)))
    print('LARGE_ORIGINAL_H',n,eta,'max_error',np.max(abs(got[0][0]-expected)),'ratio',ratio,'residual',res,flush=True)
    np.testing.assert_allclose(got[0][0],expected,atol=2e-12,rtol=2e-12)
print('LARGE_ORIGINAL_H_PASS')
