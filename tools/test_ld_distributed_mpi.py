"""mpiexec -n 4 python tools/test_ld_distributed_mpi.py
Two LD materials, two independent points each; field role has two ranks.
This is a prescribed-A service test, not implicit E/H FDTD coupling.
"""
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import SolverGroup,ParallelPlan,ParallelContext
from salmonpy.material_response import MaterialBinding,ResponseSession,LDMaterialProvider,LDResponseBackend
from salmonpy.ld import LorentzDrude
assert MPI.COMM_WORLD.size==4
plan=ParallelPlan((SolverGroup('left',(0,),point_ids=(9,3)),SolverGroup('field',(1,2),role='field'),SolverGroup('right',(3,),point_ids=(20,15))))
a=LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[0.],damping=[.1])
b=LorentzDrude(plasma_frequency=2.,strength=[.3,.7],resonance=[.5,1.5],damping=[.1,.2])
x=MaterialBinding('left',LDMaterialProvider(a),LDResponseBackend())
y=MaterialBinding('right',LDMaterialProvider(b),LDResponseBackend())
bindings={3:x,9:x,15:y,20:y}
with ParallelContext.create(MPI.COMM_WORLD,plan) as context:
    session=ResponseSession(context,bindings,dt=.01)
    reference=[(a if p<10 else b).new_state(dt=.01) for p in plan.point_ids]
    old=np.zeros((4,3))
    for step in range(1,6):
        new=np.array([[.001*p*step**2,0.,0.] for p in plan.point_ids])
        actual=session.advance(new if context.group.role=='field' and context.is_leader else None)
        if actual is not None:
            expected=np.array([s.step(-(v-u)/.01)[1] for s,u,v in zip(reference,old,new)])
            np.testing.assert_array_equal(actual,expected)
        old=new
    session.close()
if MPI.COMM_WORLD.rank==0: print('distributed two-material/four-point LD matches serial exactly')
