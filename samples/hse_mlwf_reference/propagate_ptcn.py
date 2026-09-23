"""Resumable constant-A Si PT-CN-ACE propagation, exact exchange acceptance."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse,json,time,fcntl
import numpy as np
from numpy.fft import fftn,ifftn
from model import NativeModel,hartree,semilocal_potential
from scf import HSEFunctional,Localizer,to_matrix
from semilocal import Semilocal
from ace import ACE
from ptcn import ptcn_step
from rt import validate_request
from checkpoint import FORMAT,fingerprint,save_checkpoint,load_checkpoint,write_json

def run(export,state,output,target_steps,dt=.32,amplitude=1e-4,resume=False,checkpoint_every=5):
 validate_request(target_steps,dt,amplitude,0.)
 if not isinstance(checkpoint_every,int) or checkpoint_every<1:raise ValueError('Positive checkpoint interval required')
 export=Path(export);state=Path(state);out=Path(output)
 if out.resolve()==state.parent.resolve():raise ValueError('Output must differ from ground-state input')
 out.mkdir(parents=True,exist_ok=True)
 with (out/'.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  return _run(export,state,out,target_steps,dt,amplitude,resume,checkpoint_every)

def _run(export,state,out,target,dt,amplitude,resume,interval):
 provenance=json.loads((state.parent/'result.json').read_text())
 if not provenance.get('converged') or provenance.get('final_pair_tolerance')!=0:raise ValueError('Full-pair converged SCF required')
 expected=dict(dt=dt,amplitude=amplitude,export_hash=fingerprint([*export.glob('*.bin'),export/'metadata.txt',export/'complete.txt']),initial_hash=fingerprint([state,state.parent/'result.json']))
 checkpoint=out/'restart.npz';m=NativeModel(export)
 if resume:
  u,g,metadata,rows,previous=load_checkpoint(checkpoint,expected)
 else:
  if checkpoint.exists() or (out/'status.json').exists():raise ValueError('Existing trajectory: use --resume or a new output directory')
  with np.load(state) as data:u=data['u'].copy();g=data['gauge'].copy()
  metadata=dict(format=FORMAT,step=0,initial_energy=None,**expected);rows=[];previous=to_matrix(u)@g
 if target<=metadata['step']:raise ValueError('Target must exceed the saved step')
 if u.shape!=(64,16,12,12,12) or np.max(abs(m.nlcc))>0:raise ValueError('Si reference with zero NLCC required')
 def norm_errors(x):
  flat=x.reshape(64,16,-1);gram=flat.conj()@flat.transpose(0,2,1)*m.dv
  ne=float(m.density(x).sum()*m.dv);ge=float(np.max(abs(gram-np.eye(16))))
  if not np.isfinite([ne,ge]).all() or abs(ne-32)>1e-7 or ge>1e-8:raise RuntimeError('PT-CN norm/orthogonality gate exceeded')
  return ne,ge
 norm_errors(u)
 m.set_field(np.array([0.,0.,amplitude]));loc=Localizer(m);loc.previous=previous
 cache={};start=time.perf_counter();start_step=metadata['step'];bootstrap_ready=False;builds=0;accepted=None
 def status(phase,error=None):
  saved_step=accepted[2]['step'] if accepted is not None else metadata['step']
  value=dict(status=phase,pid=os.getpid(),accepted_step=saved_step,target_steps=target,dt_au=dt,
   time_fs=saved_step*dt*.024188843265857,amplitude=amplitude,checkpoint=str(checkpoint.resolve()),
   elapsed_this_run_seconds=time.perf_counter()-start,error=error)
  write_json(out/'status.json',value);return value
 def save():save_checkpoint(checkpoint,*accepted)
 status('running')
 functional=None;xc=None
 try:
  xc=Semilocal('hse06');functional=HSEFunctional(m,xc,fftw=True)
  def local(x,energy=False):
   rho=m.density(x);vh,eh=hartree(rho,m.h);vsl,esl=semilocal_potential(rho,m.nab,xc,m.dv);core=m.core(x)
   if energy:return m.expectation(x,core)+eh+esl+.5*m.expectation(x,cache['full'])+float(m.native_energies[4])
   return core+(vh+vsl)*x
  def build(x):
   nonlocal builds
   w,_=functional.exchange(x,g,0.);op=ACE(x,w,m.dv);builds+=1
   if np.linalg.norm(op.apply(x)-w)>1e-10*np.linalg.norm(w):raise RuntimeError('ACE interpolation failure')
   cache['apply']=lambda y:.25*op.apply(y);cache['full']=.25*w
   return cache['apply'],cache['full']
  build(u);initial_action=local(u)+cache['full'];start_energy=local(u,energy=True)
  if metadata['initial_energy'] is None:metadata['initial_energy']=start_energy
  if not np.isfinite(start_energy):raise RuntimeError('Nonfinite energy')
  accepted=(u,g,metadata,rows,loc.previous);bootstrap_ready=True;save()
  for index in range(start_step,target):
   tick=time.perf_counter();localization=None
   if index and index%10==0:g,localization=loc.update(u,minimize=True)
   def precondition(r):return ifftn(fftn(r,axes=(-3,-2,-1))/(1+.5j*dt*m.tsymbol[:,None]),axes=(-3,-2,-1))
   candidate,info=ptcn_step(u,dt,local,build,m.dv,precondition=precondition,
     initial_action=initial_action,initial_exchange=cache['apply'])
   ne,ge=norm_errors(candidate);energy=local(candidate,energy=True);current=m.current(candidate)
   if not np.isfinite(energy) or not np.isfinite(current).all():raise RuntimeError('Nonfinite endpoint observable')
   next_action=local(candidate)+cache['full']
   row=dict(step=index+1,time_au=(index+1)*dt,current=current.tolist(),electron_number=ne,gram_error=ge,
     energy_Ha=energy,energy_change_Ha=energy-metadata['initial_energy'],localization=localization,
     step_wall_seconds=time.perf_counter()-tick,**info)
   # One accepted snapshot owns state, history and the localization seed.
   accepted=(candidate,g,dict(metadata,step=index+1),rows+[row],loc.previous)
   u,g,metadata,rows,_=accepted;initial_action=next_action
   write_json(out/'trajectory.json',rows)
   if (index+1)%interval==0:save()
   status('running');print(json.dumps(row),flush=True)
  save();result=status('completed')
  result.update(total_full_builds_this_run=builds,start_step=start_step,energy_change_Ha=rows[-1]['energy_change_Ha'],
    max_gram_error=max(r['gram_error'] for r in rows),max_electron_error=max(abs(r['electron_number']-32) for r in rows),
    scope='Constant-A trajectory; no spectral claim until adequate duration/analysis')
  write_json(out/'result.json',result);return result
 except BaseException as error:
  message=f'{type(error).__name__}: {error}'
  if bootstrap_ready:
   try:save()
   except Exception as checkpoint_error:message+=f'; checkpoint write also failed: {checkpoint_error}'
  try:status('failed',message)
  except Exception:pass
  raise
 finally:
  if functional is not None:functional.close()
  if xc is not None:xc.close()

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('export');p.add_argument('state');p.add_argument('output');p.add_argument('target_steps',type=int)
 p.add_argument('--dt',type=float,default=.32);p.add_argument('--amplitude',type=float,default=1e-4);p.add_argument('--resume',action='store_true');p.add_argument('--checkpoint-every',type=int,default=5)
 a=p.parse_args();print(json.dumps(run(**vars(a)),indent=2))
