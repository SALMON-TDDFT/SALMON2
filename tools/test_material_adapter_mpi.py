"""Actual2/3rank Python adapter regression; no native solver or GS data.

mpiexec -n 3 python tools/test_material_adapter_mpi.py NEW_JSON
"""
from dataclasses import dataclass
import json
from pathlib import Path
import sys
import numpy as np
import salmonpy
from mpi4py import MPI
from salmonpy.crystal import respond
from salmonpy.materials import bind_material
from salmonpy.material_response import ResponseSession
from salmonpy.distributed import ParallelContext,ParallelPlan,SolverGroup
from salmonpy.examples.custom_material import ExampleMaterial,waveform,bindings,plan_for,run_mpi

world=MPI.COMM_WORLD
result=run_mpi(world)
if world.rank==0:
    expected=np.stack([respond(ExampleMaterial(gain=1. if p in (10,30) else 2.),waveform(p)).current[1:]
                       for p in result['point_ids']],axis=1)
    actual=np.asarray(result['history'])
    assert actual.tobytes()==expected.tobytes()
    result.update(ranks=world.size,serial_maximum_difference=float(np.max(np.abs(actual-expected))),bitwise_serial_parity=True)

plan=plan_for(world.size)
with ParallelContext.create(world,plan) as context:
    values=bindings()
    if world.rank==0:
        wrong=bind_material(ExampleMaterial(gain=1.),key='one',identity='wrong-rank-config')
        values[10]=values[30]=wrong
    try:ResponseSession(context,values,dt=.02)
    except ValueError as exc:assert 'binding descriptors differ across ranks' in str(exc)
    else:raise AssertionError('mismatched identity accepted')

@dataclass(frozen=True)
class Failing(ExampleMaterial):
    backend: str = 'test-failing-material/1'
    def new_state(self,**kwargs):
        state=super().new_state(**kwargs)
        def fail(A):raise RuntimeError('injected model failure')
        state.step=fail
        return state

with ParallelContext.create(world,plan) as context:
    values=bindings();values[10]=bind_material(Failing(),key='failure',identity='intentional-test/1')
    session=ResponseSession(context,values,dt=.02)
    for attempt in range(2):
        try:session.advance(np.zeros((4,3)) if world.rank==0 else None)
        except ValueError:pass
        else:raise AssertionError('failed material session accepted')
        assert session.failed and session.step_index==0
    session.close()

if world.size==3:
    multi=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('material',(1,2),point_ids=(10,20,30,40))))
    with ParallelContext.create(world,multi) as context:
        try:ResponseSession(context,bindings(),dt=.02)
        except ValueError:pass
        else:raise AssertionError('serial material accepted internally parallel group')

world.Barrier()
if world.rank==0:
    result.update(identity_mismatch_rejected=True,material_failure_collective=True,
                  multirank_serial_group_rejected=world.size==3,
                  python=sys.executable,package=salmonpy.__file__)
    with Path(sys.argv[1]).open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print('serial/MPI parity, independent points and collective rejection PASS',world.size)
