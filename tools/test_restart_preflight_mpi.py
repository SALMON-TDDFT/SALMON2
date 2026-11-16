"""Three-rank preflight: incompatible archives must never call a LOAD callback.

mpiexec -n 3 python tools/test_restart_preflight_mpi.py NEW_OUTPUT
No native runtime required; callback spies stand in for field/material LOAD.
"""
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelPlan,SolverGroup,ParallelContext
from salmonpy.examples._direct_checkpoint import _compatibility,checkpoint,digest

world=MPI.COMM_WORLD
assert world.size==3
root=Path(sys.argv[1]).resolve()
if world.rank==0: root.mkdir(parents=True,exist_ok=False)
world.Barrier()
plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('matter',(1,2),point_ids=(1,))))
calls=[]
def forbidden(*args,**kwargs):
    calls.append('LOAD');raise AssertionError('LOAD must not run during preflight')
with ParallelContext.create(world,plan) as context:
    state=SimpleNamespace(checkpoint_protocol='salmonpy.direct-state/1',checkpoint_format='custom-numeric/1',
        checkpoint_modes=('save','restore'),checkpoint_decomposition=(2,1,1,1,1),checkpoint=forbidden)
    field=(SimpleNamespace(lib=SimpleNamespace(sp_eh_checkpoint=forbidden),_call=forbidden,
        checkpoint_decomposition=(1,1,1),step_index=0) if context.group.role=='field' else None)
    response=SimpleNamespace(context=context,states={} if field else {1:state},failed=False)
    coupled=SimpleNamespace(response=response,field=field,step_index=0)
    record=_compatibility(world,coupled,restore=False)
    assert record==_compatibility(world,coupled,restore=True,saved=record)
    for case in ('inner-layout','ownership','format','missing','unsupported','missing-payload'):
        path=root/case
        if world.rank==0:
            path.mkdir()
            for rank in range(world.size):
                data={'step':1,'J':[[0.,0.,0.]],'response_A':[[0.,0.,0.]],'states':{}}
                if rank==0: data.update(field_A=[[0.,0.,0.]],field_E=[[0.,0.,0.]])
                else: data['states']={'1':{'A':[0.,0.,0.],'current':[0.,0.,0.]}}
                (path/f'rank-{rank}.json').write_text(json.dumps(data))
            (path/'frames.json').write_text('[]')
            saved=json.loads(json.dumps(record))
            if case=='inner-layout': saved['participants'][1]['states']['1']['decomposition']=[1,2,1,1,1]
            if case=='ownership': saved['groups'][1]['ranks']=[2,1]
            if case=='format': saved['participants'][2]['states']['1']['format']='other-numeric/1'
            manifest={'schema':'salmonpy.native-direct/2','ranks':world.size,'signature':{'steps':4},'step':1,
                      'compatibility':saved,'files':{p.name:digest(p) for p in path.iterdir()}}
            if case=='missing':
                manifest['schema']='salmonpy.native-direct/1';manifest.pop('compatibility')
            (path/'complete.json').write_text(json.dumps(manifest))
        world.Barrier()
        state.checkpoint_modes=('save',) if case=='unsupported' and world.rank==2 else ('save','restore')
        rejected=False
        try: checkpoint(world,coupled,path,{'steps':4},restore=True)
        except ValueError as exc:
            assert 'compatibility' in str(exc) or 'capability' in str(exc),str(exc)
            rejected=True
        assert all(world.allgather(rejected))
        assert not calls and not response.failed and coupled.step_index==0
        if field is not None: assert field.step_index==0
    state.checkpoint_modes=('save','restore')
    if world.rank==2: state.checkpoint_format=None
    rejected=False
    try: checkpoint(world,coupled,root/'unsupported-save',{'steps':4},frames=[])
    except ValueError: rejected=True
    assert all(world.allgather(rejected)) and not calls
    assert not (root/'unsupported-save').exists()
    state.checkpoint_format='other-numeric/1' if world.rank==2 else 'custom-numeric/1'
    rejected=False
    try: checkpoint(world,coupled,root/'mixed-format-save',{'steps':4},frames=[])
    except ValueError: rejected=True
    assert all(world.allgather(rejected)) and not calls
    assert not (root/'mixed-format-save').exists()
if world.rank==0: print('restart preflight MPI PASS: exact layout admitted; 8 collective rejections; zero LOAD calls')
