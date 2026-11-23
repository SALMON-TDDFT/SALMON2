"""Actual C2H2 worker + total/reference E/H fields: outer performance regression.

Arguments: FIELD_BUILD WORKER_RECEIPT EXISTING_GS NEW_OUT MODE FAULT MEDIUM
Run with4 or6 ranks, one OMP/BLAS thread and a launcher timeout. This internal
worker validation deliberately does not issue a public runtime-validation receipt.
Only test storage uses GS read-only links; physical kernels/protocols are unchanged.
MEDIUM vacuum/ld, MODE checked/performance, FAULT0/1 (returned frame rejection).
No new GS, GS copies, or electronic disk checkpoints are collected.
"""
import hashlib,json,shutil,sys
from pathlib import Path
from unittest.mock import patch
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelContext,ParallelPlan,SolverGroup
from salmonpy.coupled import NativeEHField
from salmonpy.dipole_source import DipoleMap,MolecularFeedback
from salmonpy.fdtd import FDTD,PlaneWave,SalmonExecutable
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.ld import LorentzDrude
from salmonpy.tddft import TDDFTInitialState
from salmonpy.isolated_tddft import _EEngine,IsolatedTDDFTState
from salmonpy.stepping import StepRecoveryError
record,worker,fixture,out,mode,fault,medium=sys.argv[1:];out=Path(out);fixture=Path(fixture).resolve();world=MPI.COMM_WORLD
assert world.size in (4,6) and mode in ('checked','performance') and fault in ('0','1') and medium in ('vacuum','ld')
n=world.size//2-1;group_size=n+1;pair=world.rank//group_size;local=world.Split(pair,world.rank);dt=.003;steps=64
w=json.loads(Path(worker).read_text());exe=Path(w['executable'])
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert w['status']=='built-and-object-audited' and digest(exe)==w['sha256']
r=json.loads(Path(record).read_text());assert digest(Path(r['library']))==r['sha256']
gs=fixture/'data_for_restart';before={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
def link_copy(source,dest):
 dest.mkdir();records=[]
 for p in source.iterdir():
  if p.is_file() and p.name.endswith(('.inp','.dat')):shutil.copy2(p,dest/p.name);records.append((p.name,digest(dest/p.name)))
 target=(source/'data_for_restart').resolve();(dest/'data_for_restart').symlink_to(target,target_is_directory=True)
 records.extend(('data_for_restart/'+str(p.relative_to(target)),digest(p)) for p in target.rglob('*') if p.is_file())
 return tuple(sorted(records))
try:
 if world.rank==0:out.mkdir(parents=True,exist_ok=False)
 world.Barrier()
 plan=ParallelPlan((SolverGroup('field',tuple(range(n)),role='field'),SolverGroup('labels',(n,),point_ids=(10,))))
 with patch('salmonpy.tddft._copy_bundle',side_effect=link_copy),ParallelContext.create(local,plan) as context:
  model=state=field=None;history=[];rejected=False
  try:
   if context.group.role=='field':
    shapes=();materials=()
    if medium=='ld' and pair==0:
     shapes=(Shape('box',(6.,12.,12.),(6.,0.,0.)),);materials=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),)
    calc=FDTD(Geometry(Grid((40,40,40),(2.,)*3),shapes),materials,dt,steps,source=PlaneWave(AcosPulse(.1 if pair==0 else 0.,20.,.15),0,-16.),observations=((-1.,1.,1.),))
    field=NativeEHField(context,calc,record,out/('field-'+str(pair)),points={10:(19,20,20)},grid=(n,1,1),execution_mode=mode)
   if world.rank==0:
    model=TDDFTInitialState.from_salmon(fixture,'input.inp',executable=SalmonExecutable(exe,w['sha256']),max_steps=steps,threads=1,timeout=60)
    state=IsolatedTDDFTState(_EEngine(model,dt,(-1.,1.,1.),execution_mode=mode),(-1.,1.,1.))
   session=MolecularFeedback(world,field,state=state,mapping=DipoleMap(('molecule',),(10,),[[1.]],[8.]),total_root=0,reference_root=group_size,policy='environment_scattered_only',execution_mode=mode)
   assert session.execution_mode==mode
   def collect(frame):
    if frame is not None:
     history.append(np.r_[frame.field_step,frame.field_time,frame.total_E,frame.own_E,frame.selected_E,frame.drive_E,frame.current.ravel(),state._engine.current])
     assert not frame.total_E.flags.writeable and not frame.current.flags.writeable
   collect(session.start())
   if world.rank==0:initial=state.worker_statistics()
   while session.step_index<steps:
    inject=fault=='1' and session.step_index==32 and not rejected
    def validate(frame):
     if inject:raise ValueError('injected molecular accepted-frame rejection')
    try:frame=session.step(validate=validate)
    except StepRecoveryError as e:
     assert inject and e.original_preserved and not e.cleanup_errors
     assert session.step_index==32 and not session.failed
     if world.rank==0:
      assert state._engine._steps==31 and state.usable and state.worker_statistics()['pid']==initial['pid']
     rejected=True;continue
    collect(frame)
   if world.rank==0:
    stats=state.worker_statistics();assert stats['pid']==initial['pid'] and stats['step']==63 and stats['updates']==63+int(fault)
    np.savez(out/'traces.npz',values=np.asarray(history))
    (out/'summary.json').write_text(json.dumps(dict(mode=mode,ranks=world.size,medium=medium,steps=64,rejected=rejected,worker=stats,GS_collected=False,scope='internal worker/protocol validation, no public runtime receipt'),indent=2)+'\n')
  finally:
   if state is not None:state.close()
   if model is not None:model.close()
   if field is not None:field.close()
 local.Free()
 assert before=={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
 if world.rank==0:print('molecular outer performance PASS',world.size,mode,fault,medium,flush=True)
except BaseException:
 import traceback
 traceback.print_exc();world.Abort(1)
