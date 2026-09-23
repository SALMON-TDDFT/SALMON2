"""Self-consistent occupied-orbital HSE reference on an exported SALMON grid.

Full support is used for the energy and gradient. Pair screening is a controlled
intermediate approximation, disabled for final residual/energy verification.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,time,json
import numpy as np
from numpy.fft import fftn,ifftn
from model import NativeModel,hartree,semilocal_potential,orthonormalize
from semilocal import Semilocal
from exchange import ScreenedKernel,cell_shifts,bloch_to_wannier,wannier_to_bloch,symmetric_exchange

class HSEFunctional:
 def __init__(self,model,xc,fftw=False):
  self.model=model;self.xc=xc;mesh=round(model.nk**(1/3));n=model.shape[0]
  self.kernel=ScreenedKernel((n*mesh,)*3,model.h,.11);self.shifts=cell_shifts(mesh,n)
  self.backend=None
  if fftw:
   from fftw_backend import FFTWConvolution
   self.backend=FFTWConvolution(self.kernel.multiplier);self.kernel.convolve=self.backend.convolve
 def close(self):
  if self.backend is not None:self.backend.close()
 def evaluate(self,u,gauge,pair_tolerance=0.):
  m=self.model;start=time.perf_counter();rho=m.density(u);vh,eh=hartree(rho,m.h)
  if np.max(abs(m.nlcc))>0:raise ValueError('This initial reference requires no NLCC')
  vxc,esl=semilocal_potential(rho,m.nab,self.xc,m.dv)
  core=m.core(u);hu=core+(vh+vxc)*u;local_seconds=time.perf_counter()-start
  ku,stats=self.exchange(u,gauge,pair_tolerance);ex=.5*m.expectation(u,ku)
  hu+=.25*ku
  energy=dict(core=m.expectation(u,core),hartree=eh,semilocal=esl,screened_fock=ex,hse_fock=.25*ex,ion_ion=float(m.native_energies[4]))
  energy['total']=sum(energy[n] for n in ['core','hartree','semilocal','hse_fock','ion_ion'])
  stats.update(local_seconds=local_seconds)
  return hu,energy,stats

 def exchange(self,u,gauge,pair_tolerance=0.):
  m=self.model
  start=time.perf_counter();localized=np.einsum('knxyz,knm->kmxyz',u,gauge,optimize=True)
  w,twist=bloch_to_wannier(localized,m.k,m.h);transform_seconds=time.perf_counter()-start
  action,stats=symmetric_exchange(w,self.shifts,self.kernel,pair_tolerance)
  start=time.perf_counter();vaction=wannier_to_bloch(action,m.k,m.h,twist)
  ku=np.einsum('kmxyz,knm->knxyz',vaction,gauge.conj(),optimize=True)
  transform_seconds+=time.perf_counter()-start
  stats.update(transform_seconds=transform_seconds)
  return ku,stats


def tangent(u,v,dv):
 shape=u.shape;a=u.reshape(shape[0],shape[1],-1);b=v.reshape(a.shape)
 return (b-(b@a.conj().transpose(0,2,1)*dv)@a).reshape(shape)


def to_matrix(u):return np.array([a.reshape(a.shape[0],-1,order='F').T for a in u])

class Localizer:
 def __init__(self,model):
  root=Path(__file__).resolve().parents[2];old=root/'docs/results/si-time-wannier';sys.path.insert(0,str(old))
  from wannier import read_u,geometry
  from mlwf import transport_gauge,overlap_mesh,minimize
  self.transport=transport_gauge;self.overlap=overlap_mesh;self.minimize=minimize;self.model=model
  self.r,_,k,_,self.length,_=geometry(root)
  np.testing.assert_allclose(k,model.k,atol=1e-12)
  old_u,_=read_u(root/'calculations/si_tdcdft_k4/gs/data_for_restart');old_g=np.load(old/'mlwf_initial.npz')['gauge']
  self.previous=old_u@old_g
 def update(self,u,minimize=True):
  start=time.perf_counter();m=self.model;matrix=to_matrix(u)
  g,s=self.transport(matrix,self.previous,m.dv);log=dict(converged=True,iterations=0)
  if minimize:
   raw,nb,b,w=self.overlap(matrix,m.k,self.r,self.length,m.dv)
   g,log=self.minimize(g,raw,nb,b,w,maxiter=600,tolerance=1e-6)
   if not log['converged']:raise RuntimeError(f'MLWF did not converge: {log}')
  self.previous=matrix@g
  return g,dict(seconds=time.perf_counter()-start,iterations=log['iterations'],converged=log['converged'])


def accept_trial(energy,residual,orbitals,previous_energy,tolerance):
 if not np.isfinite(energy) or not np.isfinite(residual) or not np.isfinite(orbitals).all():
  raise FloatingPointError('Nonfinite SCF state')
 return previous_energy is None or energy<=previous_energy+tolerance


def run(directory,output,max_iterations=180):
 m=NativeModel(directory);parity=m.native_parity()
 if any(abs(value)>1e-10 for value in parity.values()):raise RuntimeError(f'Native parity failed: {parity}')
 if not np.allclose(m.weights,1/m.nk):raise ValueError('Equal weights required')
 if not np.allclose(m.occupations[:,:16],2) or np.any(m.occupations[:,16:]!=0):raise ValueError('Si16 occupied orbitals required')
 cache=Path(output);cache.mkdir(parents=True,exist_ok=True);u=orthonormalize(m.psi[:,:16],m.dv);loc=Localizer(m)
 phases=[(1e-4,3e-4),(1e-6,2e-5),(0.,1e-6)];phase=0;history=[];step=1.;previous=None;start=time.perf_counter()
 with Semilocal('hse06') as xc:
  functional=HSEFunctional(m,xc,fftw=True)
  for iteration in range(max_iterations):
   begin=time.perf_counter();g,localization=loc.update(u,minimize=(iteration%10==0));cut,tolerance=phases[phase]
   hu,energy,timing=functional.evaluate(u,g,cut);r=tangent(u,hu,m.dv)
   residual=float(np.sqrt(.5*m.expectation(r,r)));accepted=True
   if not accept_trial(energy['total'],residual,u,previous['energy'] if previous else None,max(1e-10,cut*.5)):
    accepted=False;step*=.5
    if step<1e-7:raise RuntimeError('SCF line search failed to find descent')
    u=orthonormalize(previous['u']+step*previous['direction'],m.dv)
   row=dict(iteration=iteration,phase=phase,pair_tolerance=cut,energy=energy,residual_Ha=residual,step=step,accepted=accepted,
     localization=localization,timing=timing,iteration_seconds=time.perf_counter()-begin,elapsed_seconds=time.perf_counter()-start)
   history.append(row);print(json.dumps(row),flush=True);(cache/'history.json').write_text(json.dumps(history,indent=2)+'\n')
   if not accepted:continue
   np.savez_compressed(cache/'state.npz',u=u,gauge=g,phase=phase,iteration=iteration)
   if residual<tolerance:
    if phase==len(phases)-1:
     result=dict(converged=True,parity=parity,energy=energy,residual_Ha=residual,iterations=iteration+1,wall_seconds=time.perf_counter()-start,
       scope='Self-consistent HSE06 ground state on exported SALMON Si grid; no RT spectrum yet',final_pair_tolerance=cut)
     (cache/'result.json').write_text(json.dumps(result,indent=2)+'\n');functional.close();return result
    phase+=1;previous=None;step=1.;continue
   z=ifftn(fftn(r,axes=(-3,-2,-1))/(.5+m.tsymbol[:,None]),axes=(-3,-2,-1));z=tangent(u,z,m.dv)
   direction=-z
   if previous is not None:
    numerator=m.expectation(r,z-previous['z']);denominator=previous['rz']
    beta=max(0.,min(.8,numerator/max(denominator,1e-30)))
    direction+=beta*tangent(u,previous['direction'],m.dv)
    if m.expectation(r,direction)>=0:direction=-z
    step=min(step*1.15,3.)
   if not np.isfinite(direction).all():raise FloatingPointError('Nonfinite SCF direction')
   previous=dict(u=u.copy(),direction=direction,z=z,rz=m.expectation(r,z),energy=energy['total'])
   u=orthonormalize(u+step*direction,m.dv)
 result=dict(converged=False,iterations=max_iterations,wall_seconds=time.perf_counter()-start,last=history[-1])
 (cache/'result.json').write_text(json.dumps(result,indent=2)+'\n');functional.close();return result

if __name__=='__main__':
 print(json.dumps(run(sys.argv[1],sys.argv[2],int(sys.argv[3]) if len(sys.argv)>3 else 180),indent=2))
