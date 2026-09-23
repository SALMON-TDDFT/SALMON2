"""Self-consistent implicit midpoint with ACE inner / full-exchange outer loops.

Original Bloch gauge, not PT-CN. No exchange update interval approximation:
acceptance requires a fresh full nonlinear midpoint residual every step.
"""
import numpy as np
import time

def midpoint_step(u,dt,local_action,build_exchange,precondition=None,
                  tolerance=1e-10,inner_tolerance=1e-12,max_inner=100,max_outer=12,initial_exchange=None):
 """build_exchange(x) returns (linear compressed action, full action on x).

Both include the chosen mixing coefficient. local_action rebuilds all local
nonlinear terms from x. Optional precondition approximates (I+i dt H/2)^-1.
 """
 if not np.isfinite([dt,tolerance,inner_tolerance]).all() or dt<=0 or not 0<inner_tolerance<tolerance:
  raise ValueError('Positive dt and 0 < inner tolerance < outer tolerance required')
 if any(not isinstance(v,int) or isinstance(v,bool) or v<1 for v in (max_inner,max_outer)):
  raise ValueError('Positive integer iteration limits required')
 u=np.asarray(u,dtype=complex);scale=np.linalg.norm(u)
 if not np.isfinite(u).all() or scale==0:raise ValueError('Finite nonzero state required')
 x=u.copy();start=time.perf_counter();build_seconds=0.;inner_seconds=0.;applications=0
 compressed=initial_exchange;builds=0
 if compressed is None:
  tick=time.perf_counter();compressed,_=build_exchange(x);build_seconds+=time.perf_counter()-tick;builds=1
 history=[]
 for outer in range(max_outer):
  tick=time.perf_counter()
  for inner in range(max_inner):
   residual=x-u+.5j*dt*(local_action(x)+compressed(x));applications+=1
   error=float(np.linalg.norm(residual)/scale)
   if not np.isfinite(error):raise RuntimeError('Nonfinite midpoint residual')
   if error<inner_tolerance:break
   x-=residual if precondition is None else precondition(residual)
  else:raise RuntimeError('ACE midpoint inner iteration did not converge')
  inner_seconds+=time.perf_counter()-tick
  tick=time.perf_counter();compressed,full=build_exchange(x);build_seconds+=time.perf_counter()-tick;builds+=1
  residual=x-u+.5j*dt*(local_action(x)+full)
  error=float(np.linalg.norm(residual)/scale);history.append(error)
  if not np.isfinite(error):raise RuntimeError('Nonfinite full midpoint residual')
  if error<tolerance:
   return 2*x-u,dict(builds=builds,inner_applications=applications,full_residual=error,residual_history=history,
     build_seconds=build_seconds,inner_seconds=inner_seconds,wall_seconds=time.perf_counter()-start)
 raise RuntimeError('ACE midpoint outer iteration did not converge')
