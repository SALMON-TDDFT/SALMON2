"""Run the public molecular_feedback sample with test-only read-only GS storage.

Arguments/flags are identical to python -m salmonpy.examples.molecular_feedback.
The public runtime loader and initial-state factory are used unchanged. During
restore, any field or electronic time-advance call fails the test. Source GS
files are hashed before/after; no GS copying/SCF is performed. Saved checkpoint
payloads are the requested time-evolved molecular and field state, not GS copies.
"""
import hashlib,json,sys,shutil
from pathlib import Path
from unittest.mock import patch
from mpi4py import MPI
from salmonpy.examples.molecular_feedback import main
from salmonpy.dipole_source import MolecularFeedback
from salmonpy.coupled import NativeEHField
from salmonpy.isolated_tddft import IsolatedTDDFTState,_EEngine
bundle=Path(sys.argv[4]).resolve();output=Path(sys.argv[5]);gs=bundle/'data_for_restart'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
before={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
def links(source,dest):
 dest.mkdir();records=[]
 for p in source.iterdir():
  if p.is_file() and p.name.endswith(('.inp','.dat')):shutil.copy2(p,dest/p.name);records.append((p.name,digest(dest/p.name)))
 target=(source/'data_for_restart').resolve();(dest/'data_for_restart').symlink_to(target,target_is_directory=True)
 records.extend(('data_for_restart/'+str(p.relative_to(target)),digest(p)) for p in target.rglob('*') if p.is_file())
 return tuple(sorted(records))
original=MolecularFeedback.restore_checkpoint;guarded=False;calls=0

def forbidden(*a,**kw):
 global calls
 calls+=1
 raise AssertionError('time evolution forbidden during direct restore')
def restore(self,*a,**kw):
 global guarded
 mode=None if self.state is None else self.state.execution_mode
 with patch.object(NativeEHField,'advance',forbidden),patch.object(NativeEHField,'advance_owned',forbidden),patch.object(IsolatedTDDFTState,'step',forbidden),patch.object(_EEngine,'step',forbidden):
  result=original(self,*a,**kw)
 guarded=True
 if self.state is not None:assert self.state.execution_mode==mode
 return result
try:
 with patch('salmonpy.tddft._copy_bundle',side_effect=links),patch.object(MolecularFeedback,'restore_checkpoint',restore):main()
 assert before=={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
 rows=MPI.COMM_WORLD.allgather(dict(guarded=guarded,forbidden_calls=calls))
 if '--resume' in sys.argv:assert all(r['guarded'] and r['forbidden_calls']==0 for r in rows)
 if MPI.COMM_WORLD.rank==0:(output/'restore-call-guard.json').write_text(json.dumps(dict(ranks=rows,GS_bytes_unchanged=True),indent=2)+'\n')
except BaseException:
 import traceback
 traceback.print_exc();MPI.COMM_WORLD.Abort(1)
