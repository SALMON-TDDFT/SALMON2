"""Native E/H + LD + two TDDFT/fast-PR HSE points: outer performance regression.

Arguments: BUILD_JSON EXISTING_SYNTHETIC_GS NEW_OUTPUT MODE FAULT
MODE checked/performance; FAULT 0/1. Run with 3 or 6 ranks, a launcher timeout,
and one OMP/BLAS thread. Only the test copier and matching inventory substitute read-only GS links;
production copying and numerical updates remain unchanged. No GS is collected.
"""
import hashlib,json,shutil,sys
from pathlib import Path
from unittest.mock import patch
import numpy as np
from mpi4py import MPI
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.ld import LorentzDrude
from salmonpy.distributed import ParallelContext,ParallelPlan,SolverGroup
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider,NativeTDDFTBackend
from salmonpy.coupled import NativeEHField,CoupledStepper
from salmonpy.native_steps import native_step_targets
from salmonpy.stepping import StepRecoveryError
world=MPI.COMM_WORLD
record,fixture,output,mode,fault=sys.argv[1:];fixture=Path(fixture).resolve();out=Path(output)
assert mode in ('checked','performance') and fault in ('0','1') and world.size in (3,6)
n=world.size//3;steps=96;dt=.003
record_data=json.loads(Path(record).read_text())
hse=record_data.get('schema')=='salmonpy.fast-hse-mpi/1'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert digest(Path(record_data['library']))==record_data['sha256']
gs=fixture/'data_for_restart';before={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
def link_copy(source,dest):
 dest.mkdir();records=[]
 for p in sorted(source.rglob('*')):
  rel=p.relative_to(source)
  if p.is_file():records.append((str(rel),digest(p)))
  if rel.parts[0]=='data_for_restart':continue
  if p.is_dir():(dest/rel).mkdir(exist_ok=True)
  elif p.is_file():shutil.copy2(p,dest/rel)
 (dest/'data_for_restart').symlink_to((source/'data_for_restart').resolve(),target_is_directory=True)
 return test_manifest(source)
from salmonpy._hse_assets import _path as original_asset_path
def test_asset_path(root,name):
 rel=Path(name)
 if not rel.is_absolute() and '..' not in rel.parts and rel.parts and rel.parts[0]=='data_for_restart':
  target=(root/rel).resolve()
  assert target==gs.resolve() or target.is_relative_to(gs.resolve())
  return target
 return original_asset_path(root,name)
def test_manifest(root):
 records=[]
 for p in root.rglob('*'):
  if p.is_symlink():
   assert p==root/'data_for_restart' and p.resolve()==gs.resolve()
  elif p.is_file():records.append((str(p.relative_to(root)),digest(p)))
 if (root/'data_for_restart').is_symlink():
  records.extend(('data_for_restart/'+name,sha) for name,sha in before.items())
 return tuple(sorted(records))
try:
 if world.rank==0:out.mkdir(parents=True,exist_ok=False)
 world.Barrier()
 plan=ParallelPlan((SolverGroup('field',tuple(range(n)),role='field'),
  SolverGroup('left',tuple(range(n,2*n)),point_ids=(1,)),
  SolverGroup('right',tuple(range(2*n,3*n)),point_ids=(2,))))
 grid=Grid((40,40,40),(2.,2.,2.));points={1:(19,20,20),2:(20,20,20)}
 geometry=Geometry(grid,(Shape('box',(4.,4.,4.),(9.,1.,1.)),))
 ld=LorentzDrude(plasma_frequency=20.,strength=[.7,.3],resonance=[0.,25.],damping=[2.,3.])
 calc=FDTD(geometry,(ld,),dt,steps,source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.),observations=((9.,1.,1.),))
 with patch('salmonpy.embedded_tddft._copy_bundle',side_effect=link_copy),patch('salmonpy.fast_hse_mpi._copy_bundle',side_effect=link_copy),patch('salmonpy.embedded_tddft._manifest',side_effect=test_manifest),patch('salmonpy.fast_hse._manifest',side_effect=test_manifest),patch('salmonpy._hse_assets._path',side_effect=test_asset_path),ParallelContext.create(world,plan) as context:
  field=session=None;history=[];restored=False
  try:
   if context.group.role=='field':field=NativeEHField(context,calc,record,out/'field',points=points,grid=(n,1,1),execution_mode=mode)
   if hse:
    from salmonpy.fast_hse import FastHSEMaterialProvider
    from salmonpy.fast_hse_mpi import NativeFastHSEBackend
    provider=FastHSEMaterialProvider(fixture,'rt.inp',max_steps=steps)
    backend=NativeFastHSEBackend(record,out,allow_symmetry=True,execution_mode=mode)
   else:
    provider=NativeTDDFTMaterialProvider(fixture,'rt.inp',max_steps=steps)
    backend=NativeTDDFTBackend(record,out,decomposition=(n,1,1,1,1),execution_mode=mode)
   binding=MaterialBinding('same-initial',provider,backend)
   session=ResponseSession(context,{1:binding,2:binding},dt=dt,execution_mode=mode)
   coupled=CoupledStepper(session,field)
   targets=native_step_targets(coupled)
   for i in range(steps):
    if fault=='1' and i==64:
     # Mutate the actual native electronic state before returning a known error.
     originals={}
     for p,adapter in session._adapters.items():
      original=adapter.step;originals[p]=original
      def fail(A,original=original):
       result=original(A)
       if context.solver_comm.rank==0:raise ValueError('test failure after native propagation')
       return result
      adapter.step=fail
     try:coupled.step_protected(targets)
     except StepRecoveryError as error:
      assert error.original_preserved and not error.cleanup_errors;restored=True
     else:raise AssertionError('fault did not return collectively')
     finally:
      for p,original in originals.items():session._adapters[p].step=original
     assert coupled.step_index==64 and session.step_index==64
    frame=coupled.step()
    if world.rank==0:history.append(np.stack([frame[k] for k in ('A_over_c','E_at_points','J_charge')]))
   if world.rank==0:
    values=np.array(history);assert np.max(np.abs(values[:,0]))>0 and np.max(np.abs(values[:,2]))>0
    probe=np.loadtxt(out/'field/obs1_at_point_rt.data')
    assert np.max(np.abs(probe[:,1:4]))>0
    np.savez(out/'traces.npz',values=values,ld_probe=probe)
    (out/'summary.json').write_text(json.dumps(dict(mode=mode,fault=fault=='1',ranks=world.size,field_ranks=n,material_ranks=n,steps=steps,protected_restore=restored,library_sha256=record_data['sha256'],GS_collected=False,electronic_backend='fast-pr-hse' if hse else 'standard-tddft'),indent=2)+'\n')
  finally:
   if session is not None:session.close()
   if field is not None:field.close()
 assert before=={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
 if world.rank==0:print('native outer performance PASS',mode,'fault',fault,'ranks',world.size,flush=True)
except BaseException:
 import traceback
 traceback.print_exc();world.Abort(1)
