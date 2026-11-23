"""E/H+LD and custom A/E/H callbacks: checked/performance and protected retry.

mpiexec -n 3|6 python tools/test_packet_performance_mpi.py BUILD OUT MODE FAULT
No GS, new solver builds, physical accuracy or speed claim. Custom arithmetic is
only a synthetic test callback; FDTD/LD updates use the compiled solver.
"""
import json,sys,hashlib
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelContext,ParallelPlan,SolverGroup
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.coupled import NativeEHField,PacketCoupledStepper
from salmonpy.field_input import FieldContract
from salmonpy.timing import TimeSpec,InputPort,OutputPort
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.ld import LorentzDrude
from salmonpy.native_steps import native_step_targets
from salmonpy.stepping import StepTarget,StepRecoveryError
record,output,mode,fault=sys.argv[1:];out=Path(output);world=MPI.COMM_WORLD
assert world.size in (3,6) and mode in ('checked','performance') and fault in ('0','1')
r=json.loads(Path(record).read_text());assert hashlib.sha256(Path(r['library']).read_bytes()).hexdigest()==r['sha256']
n=world.size//3;dt=.003;steps=96;grid=Grid((40,40,40),(2.,)*3);points={1:(19,20,20),2:(20,20,20)}
contracts={p:FieldContract(TimeSpec(dt,(InputPort('A',offset='.5'),InputPort('E'),InputPort('H',offset='.5')),(OutputPort('j',offset='.5'),)),
 {'A':'A/c:au','E':'E:au','H':'H:au'},[grid.coordinates[d][i[d]] for d in range(3)],'custom') for p,i in points.items()}
class Provider:
 identity='packet-test/1'
 def prepare(self,context):return None
class Backend:
 identity='packet-test-response/1';execution='collective';reduction='leader'
 def create_state(self,initial,point_id,dt,context):
  class State:
   current=np.zeros(3);count=0;fail=False;timing=contracts[point_id].timing
   def step_fields(self,packet):
    self.current=.001*(packet.A.value+packet.E.value+packet.H.value);self.count+=1
    if self.fail and context.solver_comm.rank==0:raise ValueError('injected multiport callback failure')
    return self.current
   def close(self):pass
  return State()
try:
 if world.rank==0:out.mkdir(parents=True,exist_ok=False)
 world.Barrier()
 plan=ParallelPlan((SolverGroup('field',tuple(range(n)),role='field'),SolverGroup('left',tuple(range(n,2*n)),point_ids=(1,)),SolverGroup('right',tuple(range(2*n,3*n)),point_ids=(2,))))
 with ParallelContext.create(world,plan) as context:
  field=response=None;history=[];retried=False
  try:
   if context.group.role=='field':
    calc=FDTD(Geometry(grid,(Shape('box',(4.,4.,4.),(9.,1.,1.)),)),(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),),dt,steps,source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.))
    field=NativeEHField(context,calc,record,out/'field',points=points,grid=(n,1,1),execution_mode=mode)
   b=MaterialBinding('custom',Provider(),Backend());response=ResponseSession(context,{1:b,2:b},dt=dt,field_contracts=contracts,execution_mode=mode)
   driver=PacketCoupledStepper(response,field)
   custom=[]
   for p,s in response.states.items():
    def capture(s=s):return s.current.copy(),s.count
    def restore(token,s=s):s.current,s.count=token
    custom.append(StepTarget(str(p),s,capture,restore))
   targets=native_step_targets(driver,custom_targets=custom)
   for i in range(steps):
    if fault=='1' and i==64:
     for s in response.states.values():s.fail=True
     try:driver.step_protected(targets)
     except StepRecoveryError as e:assert e.original_preserved and not e.cleanup_errors;retried=True
     else:raise AssertionError('failed update accepted')
     for s in response.states.values():s.fail=False;assert s.count==64
     assert response.step_index==driver.step_index==64
    frame=driver.step()
    if world.rank==0:history.append(np.stack([frame['samples'][q].values[0] for q in ('A','E','H')]+[frame['J_charge']]))
   if world.rank==0:
    data=np.asarray(history)
    for q in range(4):assert np.max(np.abs(data[:,q]))>0
    np.savez(out/'traces.npz',values=data)
    (out/'summary.json').write_text(json.dumps(dict(mode=mode,ranks=world.size,steps=steps,retried=retried,ports=['A','E','H']),indent=2)+'\n')
  finally:
   if response is not None:response.close()
   if field is not None:field.close()
 if world.rank==0:print('A/E/H packet performance PASS',mode,world.size,fault,flush=True)
except BaseException:
 import traceback
 traceback.print_exc();world.Abort(1)
