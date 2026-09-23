"""Short Si PT-CN-ACE validation with exact endpoint exchange caches."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,json,time,platform
import numpy as np
from numpy.fft import fftn,ifftn
from model import NativeModel,hartree,semilocal_potential
from scf import HSEFunctional,Localizer
from semilocal import Semilocal
from ace import ACE
from ptcn import ptcn_step,pt_residual
from rt import validate_request

def run(export,state,output,steps=1,dt=.32,amplitude=1e-4,exchange_backend=None):
 validate_request(steps,dt,amplitude,0.)
 state=Path(state);provenance=json.loads((state.parent/'result.json').read_text())
 if not provenance.get('converged') or provenance.get('final_pair_tolerance')!=0:raise ValueError('Converged full-pair SCF state required')
 m=NativeModel(export);data=np.load(state);u=data['u'].copy();g=data['gauge'].copy()
 if u.shape!=(64,16,12,12,12) or np.max(abs(m.nlcc))>0:raise ValueError('Si reference with zero NLCC required')
 m.set_field(np.array([0.,0.,amplitude]));out=Path(output);out.mkdir(parents=True,exist_ok=True)
 loc=None
 if exchange_backend is None:
  loc=Localizer(m);loc.previous=np.array([a.reshape(16,-1,order='F').T for a in u])@g
 cache={};rows=[];build_log=[];start=time.perf_counter()
 with Semilocal('hse06') as xc:
  functional=HSEFunctional(m,xc,fftw=True,exchange_backend=exchange_backend)
  try:
   def local(x,energy=False):
    rho=m.density(x);vh,eh=hartree(rho,m.h);vsl,esl=semilocal_potential(rho,m.nab,xc,m.dv);core=m.core(x)
    if energy:return m.expectation(x,core)+eh+esl+.5*m.expectation(x,cache['full'])+float(m.native_energies[4])
    return core+(vh+vsl)*x
   def build(x):
    tick=time.perf_counter();w,_=functional.exchange(x,g,0.);full_seconds=time.perf_counter()-tick
    tick=time.perf_counter();operator=ACE(x,w,m.dv);compression_seconds=time.perf_counter()-tick
    err=float(np.linalg.norm(operator.apply(x)-w)/np.linalg.norm(w))
    if err>1e-10:raise RuntimeError('ACE interpolation failure')
    build_log.append(dict(full_exchange_seconds=full_seconds,compression_seconds=compression_seconds,interpolation_error=err))
    cache['apply']=lambda y:.25*operator.apply(y);cache['full']=.25*w
    return cache['apply'],cache['full']
   tick=time.perf_counter();build(u);initial_action=local(u)+cache['full'];initial_energy=local(u,energy=True)
   bootstrap_seconds=time.perf_counter()-tick
   initial_residual=pt_residual(u,initial_action,m.dv)
   pt_norm=float(np.linalg.norm(initial_residual));sch_norm=float(np.linalg.norm(initial_action))
   initial_residual_Ha=float(np.sqrt(.5*m.expectation(initial_residual,initial_residual)))
   for step in range(steps):
    tick=time.perf_counter()
    if loc is not None and step and step%10==0:g,_=loc.update(u,minimize=True)
    def precondition(r):return ifftn(fftn(r,axes=(-3,-2,-1))/(1+.5j*dt*m.tsymbol[:,None]),axes=(-3,-2,-1))
    u,info=ptcn_step(u,dt,local,build,m.dv,precondition=precondition,
      initial_action=initial_action,initial_exchange=cache['apply'])
    initial_action=local(u)+cache['full']
    flat=u.reshape(64,16,-1);gram=flat.conj()@flat.transpose(0,2,1)*m.dv
    ne=float(m.density(u).sum()*m.dv);gram_error=float(np.max(abs(gram-np.eye(16))))
    if abs(ne-32)>1e-7 or gram_error>1e-8:raise RuntimeError('PT-CN norm/orthogonality gate exceeded')
    row=dict(step=step+1,time_au=(step+1)*dt,current=m.current(u).tolist(),electron_number=ne,
      gram_error=gram_error,step_wall_seconds=time.perf_counter()-tick,**info)
    rows.append(row);print(json.dumps(row),flush=True)
    (out/'trajectory.json').write_text(json.dumps(rows,indent=2)+'\n')
   final_energy=local(u,energy=True)
  finally:functional.close()
 result=dict(completed=True,steps=steps,dt_au=dt,amplitude=amplitude,method='PT-CN-ACE, Jia/Lin Eq.12',
  bootstrap_seconds=bootstrap_seconds,builds=build_log,total_full_builds=len(build_log),
  initial_pt_to_schrodinger_derivative_norm=pt_norm/sch_norm,initial_pt_residual_Ha=initial_residual_Ha,
  exchange_backend=type(exchange_backend).__name__ if exchange_backend is not None else 'full_mlwf',
  radius_bohr=getattr(exchange_backend,'radius',None),energy_change_Ha=final_energy-initial_energy,
  propagation_seconds=sum(r['step_wall_seconds'] for r in rows),mean_step_seconds=float(np.mean([r['step_wall_seconds'] for r in rows])),
  wall_seconds=time.perf_counter()-start,blas_threads=os.environ.get('OPENBLAS_NUM_THREADS'),
  fftw_threads=1 if exchange_backend is None else 0,
  platform=platform.platform(),scope='Short constant-A pilot only; no optical spectrum or long-time stability claim')
 np.savez_compressed(out/'state.npz',u=u,gauge=g)
 (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
 print(json.dumps(run(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]) if len(sys.argv)>4 else 1,
  float(sys.argv[5]) if len(sys.argv)>5 else .32,float(sys.argv[6]) if len(sys.argv)>6 else 1e-4),indent=2))
