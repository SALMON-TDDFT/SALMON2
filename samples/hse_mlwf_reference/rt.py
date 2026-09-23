"""Direct self-consistent RT-HSE pilot; four current-orbital H evaluations/step.

The RK4 stages rebuild Hartree, HSE semilocal and screened exchange from the
stage orbitals. No frozen exchange and no occupied-only action used as a linear
Hamiltonian on Taylor powers. This independent reference favors explicitness;
long-time stability and the time-step error must be checked separately.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,json,time
import numpy as np
from model import NativeModel
from semilocal import Semilocal
from scf import HSEFunctional,Localizer


def rk4_step(u,t,dt,rhs):
 k1=rhs(t,u);k2=rhs(t+dt/2,u+dt*k1/2);k3=rhs(t+dt/2,u+dt*k2/2);k4=rhs(t+dt,u+dt*k3)
 return u+dt*(k1+2*k2+2*k3+k4)/6


def validate_request(steps,dt,amplitude,pair_tolerance):
 if not isinstance(steps,int) or steps<1 or not np.isfinite([dt,amplitude,pair_tolerance]).all() or dt<=0 or pair_tolerance<0:
  raise ValueError('Positive steps/dt and finite amplitude/nonnegative cutoff required')


def run(export,state,output,steps=1,dt=.08,amplitude=1e-4,pair_tolerance=0.):
 validate_request(steps,dt,amplitude,pair_tolerance)
 provenance=json.loads((Path(state).parent/'result.json').read_text())
 if not provenance.get('converged') or provenance.get('final_pair_tolerance')!=0.:raise ValueError('A converged full-pair HSE initial state is required')
 m=NativeModel(export);data=np.load(state);u=data['u'].copy();g=data['gauge'].copy();out=Path(output);out.mkdir(parents=True,exist_ok=True)
 if u.shape!=(64,16,12,12,12):raise ValueError('Si reference shape required')
 localization=Localizer(m);localization.previous=np.array([a.reshape(16,-1,order='F').T for a in u])@g
 m.set_field(np.array([0.,0.,amplitude]));records=[];rhs_times=[];initial_energy=None;start=time.perf_counter()
 with Semilocal('hse06') as xc:
  functional=HSEFunctional(m,xc,fftw=True)
  def rhs(t,x):
   nonlocal initial_energy
   tick=time.perf_counter();hx,e,stats=functional.evaluate(x,g,pair_tolerance)
   rhs_times.append(dict(seconds=time.perf_counter()-tick,**stats))
   if initial_energy is None:initial_energy=e['total']
   return -1j*hx
  initial_current=m.current(u).tolist()
  for it in range(steps):
   tick=time.perf_counter()
   if it and it%10==0:g,_=localization.update(u,minimize=True)
   u=rk4_step(u,it*dt,dt,rhs)
   flat=u.reshape(64,16,-1);gram=flat.conj()@flat.transpose(0,2,1)*m.dv
   norm=m.density(u).sum()*m.dv;orthogonality=float(np.max(abs(gram-np.eye(16))))
   if not np.isfinite(u).all() or abs(norm-32)>1e-4 or orthogonality>1e-5:raise RuntimeError('RT norm/orthogonality failure')
   row=dict(step=it+1,time_au=(it+1)*dt,electron_number=float(norm),max_gram_error=orthogonality,current=m.current(u).tolist(),wall_seconds=time.perf_counter()-tick)
   records.append(row);print(json.dumps(row),flush=True)
   (out/'trajectory.json').write_text(json.dumps(records,indent=2)+'\n')
  _,final_energy,_=functional.evaluate(u,g,pair_tolerance)
 functional.close()
 result=dict(completed=True,steps=steps,dt_au=dt,amplitude=amplitude,pair_tolerance=pair_tolerance,initial_current=initial_current,
  energy_change_Ha=final_energy['total']-initial_energy,wall_seconds=time.perf_counter()-start,rhs_timings=rhs_times,
  mean_step_seconds=float(np.mean([r['wall_seconds'] for r in records])),scope='Short RT-HSE pilot; not an optical spectrum')
 np.savez_compressed(out/'state.npz',u=u,gauge=g)
 (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
 print(json.dumps(run(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]) if len(sys.argv)>4 else 1,
   float(sys.argv[5]) if len(sys.argv)>5 else .08,float(sys.argv[6]) if len(sys.argv)>6 else 1e-4),indent=2))
