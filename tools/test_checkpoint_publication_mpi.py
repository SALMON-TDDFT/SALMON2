"""Three-rank coordinator publication failure; deterministic rename I/O error."""
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import sys
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelContext, ParallelPlan, SolverGroup
from salmonpy.examples._direct_checkpoint import checkpoint
from salmonpy.errors import RuntimeFailure

comm=MPI.COMM_WORLD
assert comm.size==3
root=Path(sys.argv[1]).resolve()
plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('material',(1,2),point_ids=(1,))))
class State:
    checkpoint_protocol='salmonpy.direct-state/1'
    checkpoint_format='publication-fixture/1'
    checkpoint_modes=('save','restore')
    checkpoint_decomposition=(2,)
    steps=1
    A=np.zeros(3)
    current=np.zeros(3)
    def checkpoint(self,path,restore):
        assert not restore
        Path(path).write_bytes(b'user-owned-payload')

def field_save(path,restore):
    assert not restore
    Path(path.decode()).write_bytes(b'field-payload')
    return 0
with ParallelContext.create(comm,plan) as context:
    field=(SimpleNamespace(lib=SimpleNamespace(sp_eh_checkpoint=field_save),
        _call=lambda fn,*args:fn(*args),checkpoint_decomposition=(1,),step_index=1,
        A=np.zeros((1,3)),E=np.zeros((1,3))) if comm.rank==0 else None)
    response=SimpleNamespace(context=context,states={} if field else {1:State()},failed=False,
        step_index=1,_A=np.zeros((1,3)))
    coupled=SimpleNamespace(response=response,field=field,step_index=1,J=np.zeros((1,3)))
    publication_failure=(patch.object(Path,'replace',side_effect=OSError('publication disk error'))
                         if comm.rank==0 else nullcontext())
    rejected=False
    with publication_failure:
        try: checkpoint(comm,coupled,root,{'steps':4},frames=[])
        except RuntimeFailure as exc:
            assert 'publication' in str(exc) and 'disk error' in str(exc)
            rejected=True
    assert all(comm.allgather(rejected)) and response.failed
    assert not (root/'complete.json').exists()
assert not MPI.Is_finalized()
if comm.rank==0: print('checkpoint publication MPI PASS: rename failure agreed on all ranks; no completed manifest')
