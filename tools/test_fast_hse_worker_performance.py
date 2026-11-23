"""Actual serial fast-PR workers: checked/performance vs pre-change worker.

Arguments: NEW_BUILD OLD_BUILD EXISTING_FIXTURE NEW_OUT
FIXTURE contains rt.inp, pseudo1_rps.dat, sym.dat and existing data_for_restart.
Only test copy/inventory functions use GS read-only links. No SCF/GS copies.
Temporary time-evolved checkpoint files are removed on completion. Four profiles
(learned Taylor/PT-CN and full/reduced k) should be supplied as separate runs.
"""
import hashlib,json,shutil,sys,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
from salmonpy.fast_hse import FastHSEMaterialProvider,FastHSEResponseBackend
from salmonpy.errors import RuntimeFailure
new,old,fixture,output=map(Path,sys.argv[1:]);fixture=fixture.resolve();out=output.resolve();out.mkdir(parents=True,exist_ok=False)
gs=(fixture/'data_for_restart').resolve()
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
before={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
def manifest(root):
 records=[]
 for p in root.rglob('*'):
  if p.is_symlink():assert p==root/'data_for_restart' and p.resolve()==gs
  elif p.is_file():records.append((str(p.relative_to(root)),digest(p)))
 if (root/'data_for_restart').is_symlink():records.extend(('data_for_restart/'+k,v) for k,v in before.items())
 return tuple(sorted(records))
def links(source,dest):
 dest.mkdir()
 for p in source.iterdir():
  if p.is_file():shutil.copy2(p,dest/p.name)
 (dest/'data_for_restart').symlink_to(gs,target_is_directory=True)
 return manifest(source)
symmetric="yn_symmetry='y" in (fixture/'rt.inp').read_text().replace(' ','')
def drive(i):return (np.array([0.,0.,.001*np.sin(.17*i)]) if symmetric else np.array([.001*np.sin(.2*i),.0007*(np.cos(.13*i)-1),.0009*np.sin(.17*i)]))
with patch('salmonpy.fast_hse._manifest',side_effect=manifest),patch('salmonpy.fast_hse._copy_bundle',side_effect=links),tempfile.TemporaryDirectory(prefix='hse-worker-checkpoint-') as temp:
 initial=FastHSEMaterialProvider(fixture,max_steps=32).prepare(None);states=[]
 def state(label,record,mode='checked'):
  work=out/label;work.mkdir();s=FastHSEResponseBackend(record,work,execution_mode=mode,threads=1,timeout=60).create_state(initial,1,.001,None);states.append(s);return s
 try:
  if symmetric:
   for label,record,mode in (('old-reduced',old,'checked'),('new-reduced',new,'checked'),('performance-reduced',new,'performance')):
    try:state(label,record,mode)
    except RuntimeFailure as error:assert 'OPEN 1' in str(error)
    else:raise AssertionError('unsupported reduced serial initialization accepted')
   assert before=={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
   (out/'summary.json').write_text(json.dumps(dict(symmetry=True,status='existing unsupported profile unchanged',all_three_initializations_rejected=True,GS_collected=False),indent=2)+'\n')
   print('reduced serial profile unchanged: use caller-MPI symmetry route',flush=True)
   sys.exit(0)
  base=state('baseline',old);checked=state('checked',new);fast=state('performance',new,'performance');other=state('sibling',new,'performance');active=(base,checked,fast,other);traces=[[],[],[],[]]
  for i in range(1,33):
   for k,s in enumerate(active):
    previous=[(x.A.copy(),x.current.copy(),x.steps,x.energy) for x in active]
    s.step(-.7*drive(i) if k==3 else drive(i));traces[k].append(np.r_[s.current,s.energy])
    for j,x in enumerate(active):
     if j!=k:
      np.testing.assert_array_equal(x.A,previous[j][0]);np.testing.assert_array_equal(x.current,previous[j][1]);assert (x.steps,x.energy)==previous[j][2:]
   if i==8:
    fast.checkpoint(Path(temp)/'state.bin',False);sha=digest(Path(temp)/'state.bin');A=fast.A.copy()
  for values in traces[1:3]:np.testing.assert_array_equal(values,traces[0])
  control=state('negative-control',new)
  for i in range(1,33):control.step(-.7*drive(i));np.testing.assert_array_equal(np.r_[control.current,control.energy],traces[3][i-1])
  fast.close()
  for mode in ('checked','performance'):
   restored=state('restored-'+mode,new,mode);assert digest(Path(temp)/'state.bin')==sha
   restored.checkpoint(Path(temp)/'state.bin',True);restored.A=A.copy();assert restored.steps==8 and restored.execution_mode==mode
   np.testing.assert_array_equal(restored.current,np.asarray(traces[0])[7,:3])
   for i in range(9,33):restored.step(drive(i));np.testing.assert_array_equal(np.r_[restored.current,restored.energy],traces[0][i-1])
  # Mode must be negotiated before exposing a state, never silently downgraded.
  try:state('old-reject',old,'performance')
  except RuntimeFailure:pass
  else:raise AssertionError('old worker accepted performance mode')
  failed=state('failed-worker',new,'performance');failed._process.kill();failed._process.wait()
  with patch('subprocess.Popen',side_effect=AssertionError('automatic respawn forbidden')):
   try:failed.step(drive(1))
   except (RuntimeFailure,OSError):pass
   else:raise AssertionError('dead worker accepted step')
  assert failed.closed and failed.steps==0
  np.savez(out/'traces.npz',values=np.array(traces[0]),sibling=np.array(traces[3]))
 finally:
  for s in states:s.close()
assert before=={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
(out/'summary.json').write_text(json.dumps(dict(steps=32,symmetry=symmetric,current_energy_vs_checked_and_previous='array_equal',interleaved_independence=True,process_exit_direct_restore_exact=True,receiving_mode_preserved=True,old_worker_rejected=True,failed_worker_no_respawn=True,GS_bytes_unchanged=True,GS_collected=False,scope='synthetic initial orbitals; no physical convergence/speed claim'),indent=2)+'\n')
print('fast PR HSE process performance PASS',fixture.name,flush=True)
