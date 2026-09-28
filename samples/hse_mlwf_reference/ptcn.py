"""Parallel-transport Crank–Nicolson with ACE inner/full-exchange outer loops.

Jia & Lin, arXiv:1809.09609 Eq.12; row orbital storage, quadrature dv.
Unlike implicit midpoint, nonlinear trapezoidal PT is not asserted to preserve
endpoint orthogonality exactly. No reorthonormalization conceals its error.
"""
import time
import numpy as np

def pt_residual(u,hu,dv):
 shape=u.shape;a=u.reshape(*shape[:2],-1);b=hu.reshape(a.shape)
 return (b-(b@a.conj().transpose(0,2,1)*dv)@a).reshape(shape)

def ptcn_step(u,dt,local_action,build_exchange,dv,precondition=None,
              tolerance=1e-10,inner_tolerance=1e-12,max_inner=150,max_outer=12,
              initial_action=None,initial_exchange=None):
 """Return accepted endpoint and timings; exact initial H[U]U may be cached.

Any supplied initial_action must be full H[U]U at this initial time/field.
An initial ACE is merely a seed. Fresh full exchange certifies the endpoint.
 """
 if not np.isfinite([dt,dv,tolerance,inner_tolerance]).all() or dt<=0 or dv<=0 or not 0<inner_tolerance<tolerance:
  raise ValueError('Positive dt/dv and 0 < inner tolerance < outer tolerance required')
 if any(not isinstance(v,int) or isinstance(v,bool) or v<1 for v in (max_inner,max_outer)):
  raise ValueError('Positive integer iteration limits required')
 u=np.asarray(u,dtype=complex);scale=np.linalg.norm(u)
 if u.ndim<3 or not np.isfinite(u).all() or scale==0:raise ValueError('Finite nonzero k/orbital/grid state required')
 start=time.perf_counter();build_seconds=0.;inner_seconds=0.;builds=0;applications=0
 compressed=initial_exchange
 if initial_action is None:
  tick=time.perf_counter();compressed,full=build_exchange(u);build_seconds+=time.perf_counter()-tick;builds+=1
  initial_action=local_action(u)+full
 initial_action=np.asarray(initial_action)
 if initial_action.shape!=u.shape or not np.isfinite(initial_action).all():raise ValueError('Invalid initial full action')
 if compressed is None:
  tick=time.perf_counter();compressed,_=build_exchange(u);build_seconds+=time.perf_counter()-tick;builds+=1
 rhs=u-.5j*dt*pt_residual(u,initial_action,dv);x=u.copy();history=[];inexact_events=[]
 for outer in range(max_outer):
  tick=time.perf_counter();inner_errors=[];event=None
  for inner in range(max_inner):
   r=x+.5j*dt*pt_residual(x,local_action(x)+compressed(x),dv)-rhs;applications+=1
   error=float(np.linalg.norm(r)/scale)
   if not np.isfinite(error):raise RuntimeError('Nonfinite PT-CN inner residual')
   inner_errors.append(error)
   if error<inner_tolerance:break
   # An inexact inner solve may hand control to the strict full-exchange gate.
   if inner>=8 and error<.5*tolerance and min(inner_errors[-4:])>.5*min(inner_errors[-8:-4]):
    event=dict(outer_iteration=outer,inner_iteration=inner,inner_residual=error)
    inexact_events.append(event);break
   x-=r if precondition is None else precondition(r)
  else:raise RuntimeError(f'PT-CN inner iteration did not converge: residual={error:.6e}')
  inner_seconds+=time.perf_counter()-tick
  tick=time.perf_counter();compressed,full=build_exchange(x);build_seconds+=time.perf_counter()-tick;builds+=1
  r=x+.5j*dt*pt_residual(x,local_action(x)+full,dv)-rhs
  error=float(np.linalg.norm(r)/scale);history.append(error)
  if event is not None:event['full_residual']=error
  if not np.isfinite(error):raise RuntimeError('Nonfinite PT-CN full residual')
  if error<tolerance:
   return x,dict(builds=builds,inner_applications=applications,full_residual=error,residual_history=history,
     build_seconds=build_seconds,inner_seconds=inner_seconds,inexact_inner_exits=len(inexact_events),inexact_events=inexact_events,wall_seconds=time.perf_counter()-start)
 raise RuntimeError('PT-CN outer iteration did not converge')
