"""Persistent SBE through existing point-MPI adapter, with used-model broadcast.

mpiexec -n 3 python tools/test_sbe_handles_mpi.py LIBRARY NEW_JSON
No GS collection or native build. Matching library checksum sidecar required.
"""
import json,sys
from pathlib import Path
import numpy as np
from mpi4py import MPI
import salmonpy
import salmonpy.sbe as sbe
from salmonpy._sbe_handles import API
from salmonpy.materials import bind_material
from salmonpy.distributed import ParallelContext
from salmonpy.material_response import ResponseSession
from salmonpy.examples.custom_material import plan_for,waveform
from salmonpy.crystal import respond
api=API(sys.argv[1]);sbe._handles_api=lambda:api;sbe._native=lambda:api.legacy
world=MPI.COMM_WORLD
if world.rank==0:
    p=np.zeros((2,2,3,1),complex);p[0,1,2,0]=p[1,0,2,0]=.3
    initial=sbe.SBEInitialState([[-.5],[.4]],p,p*.1,[1.],np.eye(3)*4,2,1,True)
    with initial.new_state(dt=.02) as warm:warm.step([0,0,.01])
else:initial=None
initial=world.bcast(initial,root=0)
assert not initial.momentum.flags.writeable
binding=bind_material(initial,key='sbe',identity='synthetic-2band-vnl-order1/1');plan=plan_for(world.size)
history=[]
with ParallelContext.create(world,plan) as context:
    session=ResponseSession(context,{p:binding for p in plan.point_ids},dt=.02)
    try:
        if context.group.point_ids:
            assert len({id(state._raw._native_state.model) for state in session.states.values()})==1
        for n in range(1,9):
            A=np.array([waveform(p).values[n] for p in plan.point_ids]) if world.rank==0 else None
            result=session.advance(A)
            if world.rank==0:history.append(result.copy())
    finally:session.close()
if world.rank==0:
    expected=np.stack([respond(initial,waveform(p)).current[1:] for p in plan.point_ids],axis=1)
    actual=np.asarray(history);assert actual.tobytes()==expected.tobytes()
    summary={'ranks':world.size,'maximum_current_difference':float(np.max(abs(actual-expected))),
             'bitwise_serial_parity':True,'used_initial_model_broadcast':True,'native_model_shared_per_initial':True,
             'history':actual.tolist(),'library_sha256':api.sha256,'package':salmonpy.__file__}
    with Path(sys.argv[2]).open('x') as f:json.dump(summary,f,indent=2,allow_nan=False)
    print('SBE handles, used-model broadcast and point-MPI parity PASS',world.size)
