#!/usr/bin/env python3
"""Independent angular-shell identities and generated V2 primitive closure."""
import argparse
import json
import math
from pathlib import Path
import unittest
import numpy as np


def ylm(l,m,v):
    r=np.linalg.norm(v);x=v[2]/r;phi=math.atan2(v[1],v[0]);am=abs(m)
    # Associated Legendre with Condon-Shortley phase.
    p=(-1)**am*math.prod(range(1,2*am,2))*max(0.,1-x*x)**(am/2)
    if l>am:
        pm=p;p=x*(2*am+1)*p
        for ll in range(am+2,l+1):pm,p=p,((2*ll-1)*x*p-(ll+am-1)*pm)/(ll-am)
    value=math.sqrt((2*l+1)/(4*math.pi)*math.factorial(l-am)/math.factorial(l+am))*p*np.exp(1j*am*phi)
    return value if m>=0 else (-1)**am*value.conjugate()


def spinor(l,plus,m,v):
    a=m-.5;b=m+.5
    ya=ylm(l,round(a),v) if abs(a)<=l else 0j;yb=ylm(l,round(b),v) if abs(b)<=l else 0j
    if plus:return np.array([math.sqrt(max(0.,(l+m+.5)/(2*l+1)))*ya,math.sqrt(max(0.,(l-m+.5)/(2*l+1)))*yb])
    return np.array([-math.sqrt(max(0.,(l-m+.5)/(2*l+1)))*ya,math.sqrt(max(0.,(l+m+.5)/(2*l+1)))*yb])


class Angular(unittest.TestCase):
    def test_explicit_spin_spherical_harmonics(self):
        sig=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
        v=np.array([.4,-.3,.7]);r=np.linalg.norm(v);step=2e-6
        for l,plus in [(0,True),(1,False),(1,True),(2,False),(2,True)]:
            jj=l+(.5 if plus else -.5);q=2*jj+1;lam=l if plus else -(l+1)
            def wave(m,x):
                rr=np.linalg.norm(x);return rr**l*np.exp(-.7*rr*rr)*spinor(l,plus,m,x)
            density=0.;tau=0.;current=np.zeros((3,3));charge=np.zeros(3);mag=np.zeros(3);magtau=np.zeros(3)
            for m in np.arange(-jj,jj+1):
                w=wave(m,v);dw=[]
                for i in range(3):
                    shift=np.eye(3)[i]*step;dw.append((wave(m,v+shift)-wave(m,v-shift))/(2*step))
                density+=np.vdot(w,w).real;mag+=np.array([np.vdot(w,s@w).real for s in sig])
                for i in range(3):
                    tau+=np.vdot(dw[i],dw[i]).real;charge[i]+=np.vdot(w,dw[i]).imag
                    magtau+=np.array([np.vdot(dw[i],ss@dw[i]).real for ss in sig])
                    for a in range(3):current[i,a]+=np.vdot(w,sig[a]@dw[i]).imag
            amp=r**l*np.exp(-.7*r*r);dramp=(l/r-1.4*r)*amp
            self.assertAlmostEqual(density,q*amp*amp/(4*np.pi),delta=2e-10)
            self.assertAlmostEqual(tau,q*(dramp**2+l*(l+1)*(amp/r)**2)/(4*np.pi),delta=2e-8)
            expected=np.zeros((3,3));c=q*lam*amp*amp/(8*np.pi*r)
            # eps_iab rhat_b
            for i in range(3):
                for a in range(3):expected[i,a]=c*np.cross(np.eye(3)[i],np.eye(3)[a])@v/r
            np.testing.assert_allclose(current,expected,atol=2e-9)
            np.testing.assert_allclose(charge,0,atol=2e-9);np.testing.assert_allclose(mag,0,atol=2e-9);np.testing.assert_allclose(magtau,0,atol=2e-9)


def validate_v2(path,v1):
    lines=path.read_text().splitlines();n=int(lines[0].split()[0]);nc=int(lines[3].split()[2]);start=4+nc;end=start+n+1
    assert lines[1]=='KY_SOC_NLCC_V2' and lines[-1]=='END_KY_SOC_NLCC_V2'
    old=v1.read_text().splitlines();assert lines[0]==old[0] and lines[2:end]==old[2:-1],'legacy numeric rows changed'
    assert lines[end]=='CORE_AUX_PAULI_V1' and lines[end+1]=='PAULI_CLOSED_SHELL_V1 ATOMIC_UNITS TAU_NO_HALF'
    nr,nb,core=map(float,lines[end+2].split());assert int(nr)==n+1
    nb=int(nb);meta=np.array([list(map(float,x.split())) for x in lines[end+3:end+3+nb]])
    aux=np.array([list(map(float,x.split())) for x in lines[end+3+nb:-1]])
    assert aux.shape==(n+1,4+2*nb) and np.isfinite(aux).all()
    r=aux[:,0];q=meta[:,3];l=meta[:,1];lam=meta[:,2];amp=aux[:,4::2];der=aux[:,5::2]
    den=(amp*amp@q)/(4*np.pi);ar=np.zeros_like(amp);ar[1:]=amp[1:]/r[1:,None];ar[0,l==1]=der[0,l==1]
    tau=(der*der+ar*ar*l*(l+1))@q/(4*np.pi)
    cur=np.zeros(len(r));cur[1:]=(amp[1:]**2@(q*lam))/(8*np.pi*r[1:])
    for calc,stored in [(den,aux[:,1]),(tau,aux[:,2]),(cur,aux[:,3])]:np.testing.assert_allclose(calc,stored,rtol=2e-12,atol=1e-100)
    dn=2*(amp*der@q)/(4*np.pi);vw=np.zeros(len(r));nz=den>0;vw[nz]=dn[nz]**2/(4*den[nz])
    assert np.min(tau-vw)>-1e-11
    # Stored derivative identity against existing V1 finite-difference density.
    radial=np.array([list(map(float,x.split())) for x in lines[start:end]])
    np.testing.assert_allclose(dn,radial[:,3]*.529177210903**4,rtol=1e-9,atol=1e-65)
    assert np.max(abs(aux[-1,1:]))<1e-10
    derivative_errors=[]
    for stride in [4,2,1]:
        z=amp[::stride];dd=der[::stride];hh=(r[1]-r[0])*stride
        calc=(z[:-4]-8*z[1:-3]+8*z[3:-1]-z[4:])/(12*hh)
        derivative_errors.append(float(np.max(abs(calc-dd[2:-2]))/np.max(abs(der))))
    assert derivative_errors[2]<derivative_errors[1]<derivative_errors[0]
    assert derivative_errors[2]<1e-3
    return dict(nrows=n+1,blocks=nb,atomic_core_occupancy=core,integrated_pseudo_core_charge=float(np.trapezoid(4*np.pi*r*r*den,r)),minimum_tau_minus_vw=float(np.min(tau-vw)),max_core_tau=float(max(tau)),max_abs_core_current=float(max(abs(cur))),legacy_numeric_rows_identical=True,primitive_derivative_relative_errors_stride4_2_1=derivative_errors)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--v2',type=Path);p.add_argument('--v1',type=Path);a=p.parse_args()
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(Angular);result=unittest.TextTestRunner().run(suite)
    if not result.wasSuccessful():raise SystemExit(1)
    if a.v2:print(json.dumps(validate_v2(a.v2,a.v1),indent=2))
