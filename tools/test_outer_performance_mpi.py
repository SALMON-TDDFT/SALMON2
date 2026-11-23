"""Synthetic outer schedule comparison, no GS or native solver build.

Run: mpiexec -n 3 (or 5) python tools/test_outer_performance_mpi.py
Field rank plus two collective custom-material groups. Physical kernels remain
in the installed compiled field module. The custom response is a test model.
"""
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelContext,ParallelPlan,SolverGroup
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.coupled import CoupledStepper
from salmonpy.maxwell import Multiscale1D
from salmonpy.stepping import StepTarget,StepRecoveryError
world=MPI.COMM_WORLD
assert world.size in (3,5)
n=(world.size-1)//2
plan=ParallelPlan((SolverGroup('field',(0,),role='field'),
 SolverGroup('left',tuple(range(1,1+n)),point_ids=(10,)),
 SolverGroup('right',tuple(range(1+n,world.size)),point_ids=(20,))))
class Provider:
 identity='outer-performance-test/1'
 def prepare(self,context):return None
class Backend:
 identity='outer-performance-linear/1';execution='collective';reduction='leader'
 def create_state(self,initial,point_id,dt,context):
  class State:
   current=np.zeros(3)
   fail=False
   def step(self,A):
    self.current=np.array(A)*.01
    if self.fail and context.solver_comm.rank==0:raise ValueError('injected callback failure')
    return self.current
   def close(self):pass
  return State()
traces=[]
for mode in ('checked','performance'):
 with ParallelContext.create(world,plan) as context:
  bindings={p:MaterialBinding(str(p),Provider(),Backend()) for p in plan.point_ids}
  def incident(t):return np.outer(np.maximum(t,0),[0.,.0007,.0003])
  field=(Multiscale1D(2,1.,vacuum=(2,2),incident=incident).new_state(dt=.001,steps=24,execution_mode=mode)
         if context.group.role=='field' else None)
  coupled=CoupledStepper.from_multiscale(context,bindings,field,execution_mode=mode)
  response=coupled.response;local=[]
  try:
   assert response.execution_mode==mode
   try:response.execution_mode='checked'
   except AttributeError:pass
   else:raise AssertionError('mode must be read-only')
   for i in range(24):
    if i==8:
     targets=[]
     if coupled.field is not None:
      f=coupled.field
      targets.append(StepTarget('field',f,f.state.snapshot,f.state.restore))
     for p,state in response.states.items():
      state.fail=True
      targets.append(StepTarget(str(p),state,lambda s=state:s.current.copy(),lambda token,s=state:setattr(s,'current',token)))
     try:coupled.step_protected(targets)
     except StepRecoveryError as error:assert error.original_preserved and not error.cleanup_errors
     else:raise AssertionError('injected failure was not returned')
     assert coupled.step_index==8 and response.step_index==8 and not response.failed
     for state in response.states.values():state.fail=False
    result=coupled.step()
    if result is not None:local.append(np.stack((result['A_over_c'],result['J_charge'])))
   if world.rank==0:traces.append(np.array(local))
  finally:
   response.close()
   if coupled.field is not None:coupled.field.close()
if world.rank==0:
 np.testing.assert_array_equal(*traces)
 assert np.max(np.abs(traces[0]))>0
 print('outer performance MPI PASS:',world.size,'ranks; 24 nonzero updates; checked/performance exact; protected rollback/retry both modes')
