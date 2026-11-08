"""Independent 70-digit original-H oracle, including tridiagonal backward error."""
import os
import sys
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.linalg.lapack import zhetrd, zungqr
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from q_openmp import QKernel

rng = np.random.default_rng(12)
n = 24
u = np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]
levels = np.linspace(.1,5,n)
h = (u*levels)@u.conj().T
h = (h+h.conj().T)/2
b = np.column_stack((u[:,0]+.2*u[:,2], u[:,0]*(.3+.7j)+.1*u[:,4], u[:,1]+.3j*u[:,3]))
mp.mp.dps = 90
cv = lambda z: mp.mpc(float(z.real),float(z.imag))
kernel = QKernel(os.environ['HSEX_Q_KERNEL'])
from precise_response import original_gram
references=[]
etas=(5e-4,1e-6,1e-9,1e-13)
for eta in etas:
    a = mp.matrix([[cv(h[i,j])-(mp.mpf(float(levels[0]))+mp.j*mp.mpf(eta) if i==j else 0)
                    for j in range(n)] for i in range(n)])
    bm = mp.matrix([[cv(z) for z in row] for row in b])
    xm = mp.matrix(n,b.shape[1])
    for j in range(b.shape[1]): xm[:,j] = mp.lu_solve(a,bm[:,j])
    gram = bm.H*xm
    gram = (gram+gram.H)/2
    expected = np.array([[complex(z) for z in row] for row in gram.tolist()])
    references.append(expected)
    f,d,e,tau,info = zhetrd(h+2e-15*np.eye(n),lower=1,lwork=64*n)
    assert info==0
    rot = np.eye(n,dtype=complex)
    rot[1:,1:],_,info = zungqr(f[1:,:-1],tau,lwork=64*n)
    assert info==0
    args = (rot.conj().T,d,e,np.empty((n,0)),np.ones((n,1,1)),b[:,:,None],
            [3],np.array([[levels[0]]]),np.array([[2.]]),[eta],12,True)
    got,res,team = kernel.solve(*args,original_h=h)
    print('ORIGINAL_H_ORACLE eta',eta,'max_error',np.max(abs(got[0][0]-expected)), 'residual',res,flush=True)
    np.testing.assert_allclose(got[0][0],expected,rtol=2e-12,atol=2e-12)
    assert team==12
    ref,err=original_gram(h,b,levels[0],eta)
    print('PYTHON_ORIGINAL_H_ORACLE eta',eta,'max_error',np.max(abs(ref-expected)),flush=True)
    np.testing.assert_allclose(ref,expected,rtol=2e-12,atol=2e-12)
for value in (np.nan,np.inf):
    bad=h.copy();bad[-1,-1]=value
    try:
        kernel.solve(*args,original_h=bad)
    except ValueError:
        pass
    else:
        raise AssertionError('nonfinite original H accepted')
# An inaccurate preconditioner must not silently return an unconverged response.
try:
    kernel.solve(*args,original_h=h+.1*np.eye(n))
except ValueError:
    pass
else:
    raise AssertionError('unconverged original-H refinement accepted')
print('ORIGINAL_H_REFINEMENT_PASS')

if os.environ.get('HSEX_ORACLE_FIXTURE'):
    with Path(os.environ['HSEX_ORACLE_FIXTURE']).open('xb') as output:
        np.array([n,b.shape[1],len(etas)],dtype='i4').tofile(output)
        for a in (h,b,rot.conj().T): np.asarray(a.T,dtype='c16').tofile(output)
        for a in (d,e,[levels[0]],etas): np.asarray(a,dtype='f8').tofile(output)
        np.asarray(references,dtype='c16').transpose(0,2,1).tofile(output)
