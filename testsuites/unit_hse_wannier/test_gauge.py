import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]

class GaugeGradient(unittest.TestCase):
 def test_complex_gradient(self):
  with tempfile.TemporaryDirectory() as directory:
   directory=Path(directory);exe=directory/'probe';path=directory/'fixture'
   subprocess.run([os.environ.get('FC','gfortran'),str(ROOT/'src/xc/hse_wannier_gauge.f90'),str(Path(__file__).with_name('gauge_probe.f90')),
                   '-L'+str(Path(os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas'))/'lib'),'-lopenblas','-o',str(exe)],check=True,cwd=directory,capture_output=True)
   rng=np.random.default_rng(721);n=4;nk=2;nb=6;ng=35
   psi=np.array([np.linalg.qr(rng.normal(size=(ng,n))+1j*rng.normal(size=(ng,n)))[0] for _ in range(nk)])
   u=np.array([np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0] for _ in range(nk)])
   b=np.vstack((np.eye(3),-np.eye(3)));wt=np.ones(6)/2;neighbors=np.array([[0,1,0,0,1,0],[1,0,1,1,0,1]])
   raw=np.empty((nk,nb,n,n),complex);position=rng.uniform(0,2*np.pi,(ng,3))
   for k in range(nk):
    for axis in range(3):
     j=neighbors[k,axis]
     raw[k,axis]=psi[k].conj().T@(psi[j]*np.exp(-1j*position[:,axis])[:,None])
     raw[j,axis+3]=raw[k,axis].conj().T
   with path.open('wb') as f:
    f.write(np.array([n,nk,nb],np.int32).tobytes())
    for a in (u.transpose(1,2,0),raw.transpose(2,3,1,0),(neighbors.T+1).astype(np.int32),b.T,wt):
     f.write(a.tobytes(order='F'))
   subprocess.run([str(exe),str(path)],check=True)
   with Path(str(path)+'.out').open('rb') as f:
    values=np.fromfile(f,np.float64,2)
    d=np.fromfile(f,np.complex128).reshape((n,n,nk),order='F').transpose(2,0,1)
   def functional(gauge):
    m=gauge.conj().transpose(0,2,1)[:,None]@raw@gauge[neighbors]
    diag=np.diagonal(m,axis1=2,axis2=3);theta=np.angle(diag)
    centers=-np.einsum('b,ba,kbn->na',wt,b,theta)/nk
    residual=theta+np.einsum('ba,na->bn',b,centers)[None]
    return np.sum(wt[None,:,None]*(1-abs(diag)**2+residual**2))/nk
   self.assertAlmostEqual(values[0],functional(u),places=12)
   z=rng.normal(size=u.shape)+1j*rng.normal(size=u.shape);direction=(z-z.conj().transpose(0,2,1))/2
   def rotate(eps):
    e,v=np.linalg.eigh(-1j*direction)
    return u@((v*np.exp(1j*eps*e)[:,None,:])@v.conj().transpose(0,2,1))
   eps=1e-6;derivative=(functional(rotate(eps))-functional(rotate(-eps)))/(2*eps)
   self.assertAlmostEqual(derivative,-np.vdot(d,direction).real,places=7)
if __name__=='__main__':unittest.main()
