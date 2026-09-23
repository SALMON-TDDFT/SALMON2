"""Si ACE interpolation/application timings and short full-residual RT pilots."""
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
from ace_rt import midpoint_step
from rt import validate_request

def run(export,state,output,steps=1,dt=.08,amplitude=1e-4):
 validate_request(steps,dt,amplitude,0.)
 state=Path(state);provenance=json.loads((state.parent/'result.json').read_text())
 if not provenance.get('converged') or provenance.get('final_pair_tolerance')!=0:raise ValueError('Full-pair converged SCF state required')
 m=NativeModel(export);d=np.load(state);u=d['u'].copy();g=d['gauge'].copy()
 if u.shape!=(64,16,12,12,12):raise ValueError('Si reference shape required')
 if np.max(abs(m.nlcc))>0:raise ValueError('No nonzero NLCC supported')
 m.set_field(np.array([0.,0.,amplitude]));out=Path(output);out.mkdir(parents=True,exist_ok=True)
 loc=Localizer(m);loc.previous=np.array([a.reshape(16,-1,order='F').T for a in u])@g
 rows=[];build_log=[];cache={};start=time.perf_counter()
 with Semilocal('hse06') as xc:
  functional=HSEFunctional(m,xc,fftw=True)
  try:
   def local(x):
    rho=m.density(x);vh,_=hartree(rho,m.h);vsl,_=semilocal_potential(rho,m.nab,xc,m.dv)
    return m.core(x)+(vh+vsl)*x
   def build(x):
    tick=time.perf_counter();w,stats=functional.exchange(x,g,0.);full_seconds=time.perf_counter()-tick
    tick=time.perf_counter();operator=ACE(x,w,m.dv);compression_seconds=time.perf_counter()-tick
    exact_error=float(np.linalg.norm(operator.apply(x)-w)/np.linalg.norm(w))
    if exact_error>1e-10:raise RuntimeError('ACE construction interpolation failure')
    build_log.append(dict(full_exchange_seconds=full_seconds,compression_seconds=compression_seconds,
      interpolation_relative_error=exact_error,metric_condition=operator.condition_max))
    cache['apply']=lambda y:.25*operator.apply(y)
    return cache['apply'],.25*w
   # Standalone repeated application benchmark, also supplies initial energy.
   tick=time.perf_counter();w,stats=functional.exchange(u,g,0.);full_seconds=time.perf_counter()-tick
   tick=time.perf_counter();operator=ACE(u,w,m.dv);compress_seconds=time.perf_counter()-tick
   operator.apply(u);timings=[]
   for repeat in range(7):
    tick=time.perf_counter();wu=operator.apply(u);timings.append(time.perf_counter()-tick)
   interpolation=float(np.linalg.norm(wu-w)/np.linalg.norm(w))
   energy_initial=functional.evaluate(u,g,0.)[1]['total']
   # The first step separately rebuilds ACE: timings include actual build cost.
   for step in range(steps):
    tick=time.perf_counter()
    if step and step%10==0:g,_=loc.update(u,minimize=True)
    def precondition(r):return ifftn(fftn(r,axes=(-3,-2,-1))/(1+.5j*dt*m.tsymbol[:,None]),axes=(-3,-2,-1))
    u,info=midpoint_step(u,dt,local,build,precondition=precondition,initial_exchange=cache.get('apply'))
    flat=u.reshape(64,16,-1);gram=flat.conj()@flat.transpose(0,2,1)*m.dv
    ne=float(m.density(u).sum()*m.dv);gram_error=float(np.max(abs(gram-np.eye(16))))
    if abs(ne-32)>1e-7 or gram_error>1e-8:raise RuntimeError('ACE RT norm/orthogonality failure')
    row=dict(step=step+1,time_au=(step+1)*dt,current=m.current(u).tolist(),electron_number=ne,
      gram_error=gram_error,step_wall_seconds=time.perf_counter()-tick,**info)
    rows.append(row);print(json.dumps(row),flush=True)
    (out/'trajectory.json').write_text(json.dumps(rows,indent=2)+'\n')
   energy_final=functional.evaluate(u,g,0.)[1]['total']
  finally:functional.close()
 result=dict(steps=steps,dt_au=dt,amplitude=amplitude,method='ACE self-consistent implicit midpoint in original Bloch gauge (not PT-CN)',
  full_exchange_seconds=full_seconds,compression_seconds=compress_seconds,apply_seconds=timings,
  median_apply_seconds=float(np.median(timings)),interpolation_relative_error=interpolation,builds=build_log,
  energy_change_Ha=energy_final-energy_initial,wall_seconds=time.perf_counter()-start,
  mean_step_seconds=float(np.mean([r['step_wall_seconds'] for r in rows])),blas_threads=os.environ.get('OPENBLAS_NUM_THREADS'),fftw_threads=1,platform=platform.platform(),
  scope='Short pilot only; no optical spectrum or long-time stability claim')
 np.savez_compressed(out/'state.npz',u=u,gauge=g)
 (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
 print(json.dumps(run(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]) if len(sys.argv)>4 else 1,
  float(sys.argv[5]) if len(sys.argv)>5 else .08,float(sys.argv[6]) if len(sys.argv)>6 else 1e-4),indent=2))
