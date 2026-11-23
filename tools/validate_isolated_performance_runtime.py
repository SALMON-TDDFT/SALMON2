"""Issue measured revision5 isolated-worker evidence, including performance mode.

Arguments: AUDITED_WORKER_JSON EXISTING_REFERENCE NEW_OUTPUT
REFERENCE contains first/input.inp, existing first/data_for_restart and dipoles.npz.
Uses read-only GS links only in test storage. Numerical updates, protocols and
public entry checks use the actual solver. Temporary excited-state checkpoints
are deleted after successful validation. No SCF or GS copying is performed.
"""
import hashlib,json,shutil,sys,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
from salmonpy.fdtd import SalmonExecutable
from salmonpy.tddft import TDDFTInitialState
from salmonpy.isolated_tddft import (_EEngine,IsolatedTDDFTState,_CAPABILITIES,
                                    load_isolated_runtime,IsolatedTDDFTInitialState)
from salmonpy.errors import SnapshotError
worker,reference,output=map(Path,sys.argv[1:]);reference=reference.resolve();out=output.resolve();out.mkdir(parents=True,exist_ok=False)
w=json.loads(worker.read_text());exe=Path(w['executable'])
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert w['status']=='built-and-object-audited' and digest(exe)==w['sha256']
with np.load(reference/'dipoles.npz') as data:fields=data['fields'].copy();dt=float(data['dt']);expected=data['values'].copy()
assert fields.shape==(20,3) and expected.shape==(21,9)
gs=reference/'first/data_for_restart';before={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
def links(source,dest):
 dest.mkdir();records=[]
 for p in source.iterdir():
  if p.is_file() and p.name.endswith(('.inp','.dat')):shutil.copy2(p,dest/p.name);records.append((p.name,digest(dest/p.name)))
 target=(source/'data_for_restart').resolve();(dest/'data_for_restart').symlink_to(target,target_is_directory=True)
 records.extend(('data_for_restart/'+str(p.relative_to(target)),digest(p)) for p in target.rglob('*') if p.is_file())
 return tuple(sorted(records))
def equal(a,b):np.testing.assert_array_equal(a,b)
with patch('salmonpy.tddft._copy_bundle',side_effect=links),tempfile.TemporaryDirectory(prefix='isolated-performance-proof-') as temp:
 model=TDDFTInitialState.from_salmon(reference/'first','input.inp',executable=SalmonExecutable(exe,w['sha256']),max_steps=20,threads=1,timeout=60)
 states=[]
 def state(mode):
  s=IsolatedTDDFTState(_EEngine(model,dt,(0.,0.,0.),execution_mode=mode),(0.,0.,0.));states.append(s);return s
 try:
  checked,fast,other=state('checked'),state('performance'),state('performance');traces=[[],[],[]]
  for i,E in enumerate(fields,1):
   for index,s in enumerate((checked,fast,other)):
    previous=[(x.time,x._engine.current.copy()) for x in (checked,fast,other)]
    s.step(E if index<2 else -E);traces[index].append(s._engine.current.copy())
    for k,x in enumerate((checked,fast,other)):
     if k!=index:assert x.time==previous[k][0];equal(x._engine.current,previous[k][1])
   if i==8:
    replay=checked.snapshot();fast.save_direct_checkpoint(Path(temp)/'direct');token=fast.savepoint();pid=fast.worker_statistics()['pid']
  equal(traces[0],traces[1]);control=state('checked')
  for i,E in enumerate(fields):control.step(-E);equal(control._engine.current,traces[2][i])
  values=np.array(traces[0]);reference_error=float(np.max(np.abs(values[:,:6]-expected[1:,:6])))
  mask=expected[1:,8]==1.;energy_error=float(np.max(np.abs(values[mask,7]-expected[1:][mask,7])))
  assert reference_error<=1e-10 and energy_error<=1e-10
  previous=fast.worker_statistics();fast.restore_savepoint(token);stats=fast.worker_statistics()
  assert stats['pid']==pid and stats['step']==8 and stats['updates']==previous['updates']
  for i in range(8,20):fast.step(fields[i]);equal(fast._engine.current,values[i])
  fast.release_savepoint(token)
  try:fast.restore_savepoint(token)
  except SnapshotError:pass
  else:raise AssertionError('released token accepted')
  fast.close()
  direct_stats=[]
  for mode in ('checked','performance'):
   s=state(mode);s.restore_direct_checkpoint(Path(temp)/'direct');stats=s.worker_statistics()
   assert stats['step']==8 and stats['updates']==0 and stats['pid']!=pid and s.execution_mode==mode
   direct_stats.append(stats)
   for i in range(8,20):s.step(fields[i]);equal(s._engine.current,values[i])
  s=state('performance');s.restore(replay)
  for i in range(8,20):s.step(fields[i]);equal(s._engine.current,values[i])
 finally:
  for s in states:s.close()
  model.close()
 report=dict(schema='salmonpy.isolated-worker-validation/5',protocol='SALMONPY_TDDFT_E_V3',executable_sha256=w['sha256'],steps=20,capabilities=sorted(_CAPABILITIES),interleaved_independence=True,direct_disk_process_exit_exact=True,direct_disk_extra_evolutions=0,savepoint_same_process=True,savepoint_extra_evolutions=0,savepoint_stale_rejected=True,repeat_max_abs_error=0.,replay_continuation_max_abs_error=0.,reference_max_abs_error=reference_error,reference_energy_max_abs_error=energy_error,savepoint_continuation_max_abs_error=0.,execution_modes=['checked','performance'],performance_protocol='SALMONPY_TDDFT_E_PERFORMANCE_V3',performance_direct_disk_exact=True,mode_retained_after_restore=True)
 candidate=out/'candidate.json';candidate.write_text(json.dumps(report,indent=2)+'\n')
 # Validate the public loader and model entry before publishing the final record.
 runtime=load_isolated_runtime(exe,candidate)
 with IsolatedTDDFTInitialState.from_salmon(reference/'first','input.inp',runtime=runtime,max_steps=20,threads=1) as initial:
  for mode in ('checked','performance'):
   with initial.new_state(dt=dt,execution_mode=mode) as s:
    for i,E in enumerate(fields):s.step(E);equal(s._engine.current,values[i])
 candidate.rename(out/'report.json')
 np.savez(out/'observables.npz',values=values,fields=fields,dt=dt)
 (out/'summary.json').write_text(json.dumps(dict(public_entry_exact=True,direct_restore_statistics=direct_stats,GS_collected=False,scope='short fixed-ion PZ C2H2 protocol regression; not physical convergence',report=report),indent=2)+'\n')
assert before=={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
print('public isolated performance runtime evidence PASS',w['sha256'],reference_error,energy_error,flush=True)
