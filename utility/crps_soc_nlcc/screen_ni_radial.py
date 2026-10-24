#!/usr/bin/env python3
"""Independent radial semilocal/reader-equivalent KB screen, not AE certification.

Uses NumPy dense Hermitian diagonalization of second-order radial finite differences.
All lengths Angstrom and energies eV. Reader mapping: beta=(V_j-Vlocal)*u_j,
N=integral u_j beta dr; angular j channels decouple in spherical symmetry.
Includes the spherical atomic screening reconstructed from the generator's
occupations/waves/NLCC; separately checks unscreened ionic operators.
"""
import argparse
import json
from pathlib import Path
import numpy as np

KINETIC=1.97329e3**2/(2*.511003e6)


def interp(x, r, y):
    # Cubic interpolant through four uniform-grid points; linear at endpoints.
    out=np.interp(x,r,y); h=r[1]; k=np.floor(x/h).astype(int); mask=(k>=1)&(k+2<len(r))
    t=x[mask]/h-k[mask]; j=k[mask]
    out[mask]=(-t*(t-1)*(t-2)*y[j-1]+3*(t+1)*(t-1)*(t-2)*y[j]-3*(t+1)*t*(t-2)*y[j+1]+(t+1)*t*(t-1)*y[j+2])/6
    return out


def screening(rows, channels):
    r=rows[:,0]; h=r[1]; u=rows[:,6::2]
    rho=np.zeros(len(r));rho[1:]=(u[1:]**2@channels[:,3])/(4*np.pi*r[1:]**2);rho[0]=rho[1]
    a=r*r*rho;b=r*rho
    vh=np.zeros(len(r));vh[1:]=4*np.pi*14.39976*h*(np.cumsum(a)[1:]/r[1:]+(np.sum(b)-np.cumsum(b))[1:]);vh[0]=vh[1]
    rs=(3/(4*np.pi*np.maximum(rho+rows[:,2],1e-80)))**(1/3)/.529177
    beta=1/(71.4036*rs);eta=np.sqrt(1+beta*beta)
    # asinh(beta) avoids cancellation in log(beta+sqrt(1+beta^2)).
    psix=-.5+1.5*np.arcsinh(beta)/beta/eta
    xc=-4/3*.4582/rs*psix
    large=rs>1;ss=rs[large]
    xc[large]+=-.1423*(1+7/6*1.0529*np.sqrt(ss)+4/3*.3334*ss)/(1+1.0529*np.sqrt(ss)+.3334*ss)**2
    ss=rs[~large];xc[~large]+=.0311*np.log(ss)+(-.048-.0311/3)+2/3*.002*ss*np.log(ss)+(2*(-.0116)-.002)/3*ss
    return vh+2*13.6058*xc


def run(path, grids):
    lines=path.read_text().splitlines();nc=int(lines[3].split()[2]);channels=np.array([list(map(float,x.split())) for x in lines[4:4+nc]])
    rows=np.array([list(map(float,x.split())) for x in lines[4+nc:-1]]);r=rows[:,0];h=r[1];common=screening(rows,channels)
    results=[]
    for box,step in grids:
        nr=round(box/step);x=np.arange(1,nr)*step;size=len(x)
        kinetic=np.diag(np.full(size,2*KINETIC/step**2))+np.diag(np.full(size-1,-KINETIC/step**2),1)+np.diag(np.full(size-1,-KINETIC/step**2),-1)
        for k,ch in enumerate(channels):
            ll=int(ch[0]);v=rows[:,5+2*k];u=rows[:,6+2*k];vl=rows[:,1];beta=(v-vl)*u
            norm=h*(np.sum(u*beta)-.5*(u[0]*beta[0]+u[-1]*beta[-1]))
            bx=interp(x,r,beta);ux=interp(x,r,u);base=kinetic+np.diag(KINETIC*ll*(ll+1)/x**2)
            for screened in [False,True]:
                extra=common if screened else np.zeros_like(r)
                semi=base+np.diag(interp(x,r,v+extra))
                kb=base+np.diag(interp(x,r,vl+extra))
                if np.max(np.abs(v-vl))>1e-12:
                    assert abs(norm)>1e-12
                    kb+=step*np.outer(bx,bx)/norm
                es=np.linalg.eigvalsh(semi);ek=np.linalg.eigvalsh(kb)
                d=dict(box_angstrom=box,step_angstrom=step,l=ll,two_j=int(ch[1]),screened=screened,kb_norm=norm,
                    semilocal_low_eV=es[:8].tolist(),kb_low_eV=ek[:8].tolist(),semilocal_bound_count=int(np.sum(es<-.05)),kb_bound_count=int(np.sum(ek<-.05)),
                    reference_eV=ch[4],reference_wave_norm=float(step*np.sum(ux**2)))
                if screened:
                    d['semilocal_nearest_reference_eV']=float(es[np.argmin(abs(es-ch[4]))]);d['kb_nearest_reference_eV']=float(ek[np.argmin(abs(ek-ch[4]))])
                    d['semilocal_reference_residual_eV']=float(np.linalg.norm((semi@ux-ch[4]*ux)[3:-3])/np.linalg.norm(ux[3:-3]))
                    d['kb_reference_residual_eV']=float(np.linalg.norm((kb@ux-ch[4]*ux)[3:-3])/np.linalg.norm(ux[3:-3]))
                if screened and box==8 and step==.01:
                    # Driven boundary u(R)=1. The tiny endpoint projector term is
                    # included with trapezoidal half-weight; no AE comparison.
                    values=[]
                    for energy in [0.,2.,5.,10.,20.,50.]:
                        rhs=np.zeros(size);rhs[-1]=KINETIC/step**2
                        us=np.linalg.solve(semi-energy*np.eye(size),rhs)
                        rhs_k=rhs.copy()
                        if abs(norm)>1e-12:
                            rhs_k-=.5*step*bx*interp(np.array([box]),r,beta)[0]/norm
                        uk=np.linalg.solve(kb-energy*np.eye(size),rhs_k)
                        values.append(dict(energy_eV=energy,semilocal_log_derivative_inv_angstrom=float((3-4*us[-1]+us[-2])/(2*step)),kb_log_derivative_inv_angstrom=float((3-4*uk[-1]+uk[-2])/(2*step))))
                    d['finite_radius_scattering_log_derivatives']=values
                d['additional_bound_state_flag']=bool(np.sum(ek<-.05)>np.sum(es<-.05))
                d['deep_extra_state_flag']=bool(ek[0]<es[0]-.5)
                results.append(d)
    return dict(method=__doc__,pseudo=str(path),checks=results,absence_of_ghosts_certified=False,all_electron_scattering_validated=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('pseudo',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--grids',default='8:.04,12:.04,8:.02');a=p.parse_args()
    report=run(a.pseudo,[tuple(map(float,x.split(':'))) for x in a.grids.split(',')]);a.output.write_text(json.dumps(report,indent=2)+'\n')
    print('records',len(report['checks']),'deep-extra-state flags',sum(x['deep_extra_state_flag'] for x in report['checks']))
