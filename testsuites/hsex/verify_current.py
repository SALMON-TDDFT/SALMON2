"""Check the nonlocal current against a fixed-orbital energy derivative at t=0+.

Use the same GS and impulse for hsex and hsex_hartree, dt <= 1e-6 au, nt=1.
This checks the native RT current including its factor of two and sign.
"""
import argparse
from pathlib import Path
import numpy as np


def verify(gs, screened, control):
    with (screened/'hsex-w.bin').open('rb') as f:
        assert f.read(4)==b'HXR1'
        v=np.fromfile(f,'=i4',9)
        n=v[:3]; mesh=v[3:6]; b,nv,ni=v[6:];nk=int(mesh.prod());ng=int(n.prod());na=b**3
        cell=np.fromfile(f,'=f8',9).reshape(3,3);k=np.fromfile(f,'=f8',3*nk).reshape(nk,3)
        np.fromfile(f,'=f8',3*ni+ng+1)
        R=np.fromfile(f,'=f8',3*nk).reshape(nk,3)
        W=np.fromfile(f,'=f8').reshape(nk,na,na).transpose(0,2,1)
    pos=np.indices(tuple(n[::-1])).reshape(3,-1)[::-1].T
    label=np.ravel_multi_index(tuple((pos//(n//b)).T),(b,)*3,order='F')
    r=pos/n@cell;vol=abs(np.linalg.det(cell))
    psi=np.fromfile(gs/'data_for_restart/wfn.bin',dtype='=c16').reshape(nk,-1,ng)[:,:nv]
    def density(A):
        amp=psi*np.exp(1j*(k+A)@r.T)[:,None,:]
        p=np.stack([amp[:,:,label==j].sum(-1) for j in range(na)],axis=1)*np.sqrt(vol*na)/ng
        pk=p@p.conj().transpose(0,2,1)
        return (np.exp(-1j*R@(k+A).T)@pk.reshape(nk,-1)).reshape(nk,na,na)/nk
    ref=density(np.zeros(3))
    def energy(A):
        return -np.sum(W*abs(density(A)-ref)**2)
    on,off=[np.loadtxt(d/'Si_rt.data',ndmin=2)[0] for d in (screened,control)]
    np.testing.assert_allclose(on[:13],off[:13])
    if on[0]>1e-6: raise ValueError('test requires negligible orbital evolution')
    field=on[1:4];h=1e-7
    dE=np.array([(energy(field+h*e)-energy(field-h*e))/(2*h*vol) for e in np.eye(3)])
    diff=on[13:16]-off[13:16]
    np.testing.assert_allclose(diff,dE,rtol=3e-5,atol=1e-13)
    print('finite-difference and native current',dE,diff)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for a in ('gs','screened','control'):p.add_argument(a,type=Path)
    a=p.parse_args();verify(a.gs,a.screened,a.control)
