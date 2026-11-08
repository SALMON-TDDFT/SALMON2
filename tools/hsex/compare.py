"""Matched Si transverse impulse spectra from the two HSEX modes and their zero-field controls."""
import argparse
from pathlib import Path
import numpy as np


def read_case(path):
    raw=(path/'Si_rt.data').read_bytes()
    return np.loadtxt(raw[:raw.rfind(b'\n')+1].splitlines(),usecols=(0,3,15),ndmin=2)


def compare(root,out,nsteps=None):
    cases=[read_case(root/d) for d in ('hartree','hartree0','sex','sex0')]
    n=min(map(len,cases)) if nsteps is None else nsteps
    if n<2 or min(map(len,cases))<n: raise ValueError('insufficient common trajectory')
    t=cases[0][:n,0];dt=t[0]
    np.testing.assert_allclose(t,dt*np.arange(1,n+1),atol=1e-8,rtol=0)
    kick=cases[0][0,1]
    if kick==0: raise ValueError('missing impulse')
    for a,k in zip(cases,(kick,0,kick,0)):
        np.testing.assert_array_equal(a[:n,0],t);np.testing.assert_allclose(a[:n,1],k)
    j=np.stack([a[:n,2]-b[:n,2] for a,b in zip(cases[::2],cases[1::2])],axis=1)
    w=np.arange(1,251)*.001;e=w*27.211386245988
    u=t/t[-1];z=1+4j*np.pi*dt/(kick*w[:,None])*(np.exp(1j*w[:,None]*t)@((1-3*u*u+2*u**3)[:,None]*j))
    if not np.isfinite(z).all(): raise ValueError('nonfinite dielectric')
    m=(e>=1)&(e<=6)
    np.savetxt(out.with_suffix('.tsv'),np.c_[e,z[:,0].real,z[:,0].imag,z[:,1].real,z[:,1].imag][m],
               delimiter='\t',header='E_eV Re_H Im_H Re_HSEX Im_HSEX')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,sharex=True,figsize=(8.5,6.6),layout='constrained')
    is_tb='tbmbj' in (root/'sex/input').read_text().lower()
    labels=('TB-mBJ + Hartree','TB-mBJ + Hartree + screened exchange') if is_tb else (
        'Induced Hartree','Induced Hartree + screened exchange')
    for i,label in enumerate(labels):
        for ax,y in zip(axes,(z.real,z.imag)): ax.plot(e[m],y[m,i],['-','--'][i],label=label)
        y=z.imag[:,i];ii=np.flatnonzero((y[1:-1]>y[:-2])&(y[1:-1]>y[2:]))+1
        print(label,[(round(e[k],4),round(y[k],4)) for k in ii if 2.5<e[k]<5])
    for ax in axes:
        ax.axvspan(3.15,3.65,color='0.88');ax.grid(alpha=.25);ax.set_xlim(1,6)
    ref='time-dependent TB-mBJ' if is_tb else 'static KLI-SIC reference'
    axes[0].set_title(f'Si · {ref} · {n:,} steps · {t[-1]*.0241888433:.2f} fs')
    axes[0].set_ylabel('Re εzz');axes[0].legend()
    axes[1].set(ylabel='Im εzz',xlabel='Energy (eV)',ylim=(0,None));fig.savefig(out.with_suffix('.png'),dpi=180)
    print('steps',n,'h/T (eV)',2*np.pi*27.211386245988/t[-1])
    return z


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--steps',type=int);a=p.parse_args();compare(a.root,a.output,a.steps)
