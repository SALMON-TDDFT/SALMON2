"""Three-rank collective admission regression; no native runtime is required.

Run with a launcher timeout: mpiexec -n 3 python tools/test_response_admission_mpi.py
"""
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelContext, ParallelPlan, SolverGroup
from salmonpy.material_response import MaterialBinding, ResponseSession

world=MPI.COMM_WORLD
assert world.size==3
plan=ParallelPlan((SolverGroup('field',(0,),role='field'),
                   SolverGroup('custom',(1,2),point_ids=(10,))))
calls=[]
class Provider:
    identity='admission-initial/1'
    def prepare(self,context): calls.append('prepare');return None
class Backend:
    identity='admission-backend/1';execution='collective';reduction='leader'
    def create_state(self,initial,point_id,dt,context):
        calls.append('create')
        class State:
            current=np.zeros(3)
            def step(self,A): self.current=np.asarray(A)*2;return self.current
            def close(self): pass
        return State()

with ParallelContext.create(world,plan) as context:
    # Valid locally on every rank, but contradictory declarations across groups
    # or within one material group must also fail before callbacks.
    cases=[(1,),(),(True,),(1,1),'1,2',
           (1,2) if world.rank==0 else (2,3),
           (1,2) if world.rank!=2 else (2,3),
           (1,) if world.rank==1 else (1,2),
           (1,2) if world.rank==0 else None]
    for sizes in cases:
        backend=Backend();backend.supported_group_sizes=sizes
        binding=MaterialBinding('custom',Provider(),backend)
        rejected=False
        try: ResponseSession(context,{10:binding},dt=.01)
        except ValueError: rejected=True
        assert all(world.allgather(rejected)), 'admission error not collective'
        assert not calls, 'callbacks ran before admission agreement'
        assert context._dependents==0
    # Equivalent orderings agree. Absence retains the legacy custom contract.
    for declared in (True,False):
        backend=Backend()
        if declared: backend.supported_group_sizes=(1,2) if world.rank%2 else (2,1)
        session=ResponseSession(context,{10:MaterialBinding('custom',Provider(),backend)},dt=.01)
        for step in range(1,5):
            A=np.array([[.001*step,0,0]])
            j=session.advance(A if context.group.role=='field' else None)
            if j is not None: np.testing.assert_array_equal(j,A*2)
        session.close()
        assert context._dependents==0
    assert world.allgather(len(calls))==[0,4,4], 'prepare/create count changed'
if world.rank==0: print('response admission MPI PASS: 9 collective rejections, 2 exact response trajectories')
