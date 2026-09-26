"""Frozen Gamma-density-factor screening diagnostic; never a production Fock operator.

NPZ: q[nfactor,nx,ny,nz] complex density factors (Psi sqrt(f/2) U),
spacing[3] in bohr, omega scalar in bohr^-1. No unit-occupation assumption.
Energy is fragment-periodic HSE exchange (alpha=.25, both spins), not DC core energy.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def analyze(q,spacing,omega,budget):
    return analyze_budgets(q,spacing,omega,[budget])[0]


def analyze_budgets(q,spacing,omega,budgets):
    q=np.asarray(q,dtype=complex);spacing=np.asarray(spacing,dtype=float)
    if q.ndim!=4 or min(q.shape)<1 or not np.isfinite(q).all():
        raise ValueError('q must be finite nonempty (nfactor,nx,ny,nz) density factors')
    if spacing.shape!=(3,) or not np.isfinite(spacing).all() or np.any(spacing<=0):
        raise ValueError('spacing must contain three positive finite bohr lengths')
    budgets=np.asarray(budgets,dtype=float)
    if budgets.ndim!=1 or not len(budgets) or not np.isfinite(budgets).all() or np.any(budgets<0):
        raise ValueError('budgets must be finite nonnegative values')
    if not np.isfinite(omega) or omega<=0:
        raise ValueError('omega must be positive and budget nonnegative, both finite')
    dv=float(np.prod(spacing));n=len(q);alpha=.25
    axes=[2*np.pi*np.fft.fftfreq(s,d=h) for s,h in zip(q.shape[1:],spacing)]
    g2=sum(a*a for a in np.meshgrid(*axes,indexing='ij'))
    kernel=np.full(q.shape[1:],np.pi/omega**2)
    np.divide(4*np.pi*(-np.expm1(-g2/(4*omega**2))),g2,out=kernel,where=g2>0)
    density=abs(q.reshape(n,-1))**2
    overlap=dv*(density@density.T)
    # Discrete Parseval: <rho,V rho> <= max(V_G) ||rho||_2^2.
    # Each off-diagonal pair occurs twice in the spin-summed exchange energy.
    bounds=2*alpha*(np.pi/omega**2)*overlap
    keep=np.ones((len(budgets),n,n),bool);discarded=np.zeros(len(budgets))
    candidates=sorted((float(bounds[i,j]),i,j) for i in range(n) for j in range(i+1,n))
    for bound,i,j in candidates:
        accepted=discarded+bound<=budgets
        keep[accepted,i,j]=keep[accepted,j,i]=False
        discarded[accepted]+=bound
    full=np.zeros_like(q);screened=np.zeros((len(budgets),)+q.shape,dtype=complex)
    for i in range(n):
        for j in range(i,n):
            rho=q[i].conj()*q[j]
            potential=np.fft.ifftn(np.fft.fftn(rho)*kernel)
            full[j]-=q[i]*potential
            for index in np.flatnonzero(keep[:,i,j]):
                screened[index,j]-=q[i]*potential
            if i!=j:
                full[i]-=q[j]*potential.conj()
                for index in np.flatnonzero(keep[:,i,j]):
                    screened[index,i]-=q[j]*potential.conj()
    energy=alpha*dv*np.vdot(q,full).real
    reports=[]
    for index,budget in enumerate(budgets):
        screened_one=screened[index]
        pruned_energy=alpha*dv*np.vdot(q,screened_one).real
        metric=-dv*q.reshape(n,-1).conj()@screened_one.reshape(n,-1).T
        anti=float(np.linalg.norm(metric-metric.conj().T)/max(np.linalg.norm(metric),1e-300))
        hermitian_part=(metric+metric.conj().T)/2
        ev=np.linalg.eigvalsh(hermitian_part)
        minimum=float(ev[0]);maximum=float(ev[-1])
        zero=bool(np.all(screened_one==0))
        # A necessary ACE gate only: passing does not prove variational SCF consistency.
        gate=anti<=1e-10 and (zero or (maximum>0 and minimum>1e-12*maximum))
        kept=int(keep[index].sum())
        reports.append(dict(budget_Ha=float(budget),total_ordered_pairs=n*n,kept_ordered_pairs=kept,
                    kept_unique_pairs=(kept+n)//2,removed_pair_fraction=1-kept/(n*n),
                    full_exchange_Ha=float(energy),screened_exchange_Ha=float(pruned_energy),
                    exchange_error_Ha=float(pruned_energy-energy),discarded_energy_bound_Ha=float(discarded[index]),
                    relative_action_error=float(np.linalg.norm(screened_one-full)/max(np.linalg.norm(full),1e-300)),
                    screened_metric_antihermitian_relative=anti,screened_metric_min_eigenvalue=minimum,
                    necessary_ace_metric_gate_passed=bool(gate),certified_for_scf=False))
    return reports


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture',type=Path)
    parser.add_argument('--budgets',nargs='+',type=float,default=[0,1e-6,1e-5,1e-4])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with np.load(args.fixture,allow_pickle=False) as data:
        q=data['q'];spacing=data['spacing'];omega=float(data['omega'])
        results=analyze_budgets(q,spacing,omega,args.budgets)
    payload=dict(fixture=str(args.fixture.resolve()),scope='frozen Gamma fragment-periodic exchange, not DC core energy',
                 shape=list(q.shape),spacing_bohr=spacing.tolist(),omega_bohr_inverse=omega,results=results)
    args.output.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print(json.dumps(payload,indent=2,allow_nan=False))

if __name__=='__main__': main()
