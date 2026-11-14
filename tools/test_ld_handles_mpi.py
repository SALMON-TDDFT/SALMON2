"""Persistent LD point-MPI parity and explicit session cleanup (2/3 ranks).

mpiexec -n 3 python tools/test_ld_handles_mpi.py LIBRARY NEW_JSON
Synthetic fields only; no native build or GS collection.
"""
import json,sys
from pathlib import Path
import numpy as np
from mpi4py import MPI
import salmonpy
import salmonpy.ld as ld
from salmonpy._ld_handles import API
from salmonpy.distributed import ParallelContext
from salmonpy.material_response import MaterialBinding,ResponseSession,LDMaterialProvider,LDResponseBackend
from salmonpy.examples.custom_material import plan_for,waveform
api=API(sys.argv[1]);ld._handles_api=lambda:api;ld._native=lambda:api.legacy
world=MPI.COMM_WORLD
if world.rank==0:
    initial=ld.LorentzDrude(plasma_frequency=1.,strength=[.7,.3],resonance=[0.,1.4],damping=[.1,.2])
    with initial.new_state(dt=.02) as warm:warm.step([.01,.02,.03])
else:initial=None
initial=world.bcast(initial,root=0)
plan=plan_for(world.size);binding=MaterialBinding('ld',LDMaterialProvider(initial),LDResponseBackend())
finalizers=[];original=ld.LorentzDrude.new_state

def tracked(self,**kwargs):
    state=original(self,**kwargs);finalizers.append(state._native_state.finalizer);return state

ld.LorentzDrude.new_state=tracked
history=[]
with ParallelContext.create(world,plan) as context:
    session=ResponseSession(context,{p:binding for p in plan.point_ids},dt=.02)
    try:
        for n in range(1,9):
            A=np.array([waveform(p).values[n] for p in plan.point_ids]) if world.rank==0 else None
            current=session.advance(A)
            if world.rank==0:history.append(current.copy())
    finally:session.close()
assert all(not f.alive for f in finalizers)
closed=world.gather(len(finalizers),root=0)
ld.LorentzDrude.new_state=original
if world.rank==0:
    reference=[]
    for p in plan.point_ids:
        A=waveform(p).values;E=np.zeros_like(A);E[1:]=-(A[1:]-A[:-1])/.02
        reference.append(initial.respond(E,dt=.02).current[1:])
    expected=np.stack(reference,axis=1);actual=np.asarray(history)
    assert actual.tobytes()==expected.tobytes()
    result={'ranks':world.size,'bitwise_serial_parity':True,'maximum_current_difference':float(np.max(abs(actual-expected))),
            'used_material_broadcast':True,'closed_native_states_per_rank':closed,
            'history':actual.tolist(),'library_sha256':api.sha256,'package':salmonpy.__file__}
    with Path(sys.argv[2]).open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print('LD handles, session cleanup and point-MPI parity PASS',world.size)
