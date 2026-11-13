"""User-defined A/c-to-J response on two internally parallel material groups.
mpiexec -n 5 python -m salmonpy.examples.mpi_response NEW_JSON
Optional: --bridge PATH --sha256 HASH to exercise the host-MPI Fortran bridge.
One material/one point per group, two ranks inside each material group.
Synthetic memory law only: field role supplies A; this is not an FDTD solver.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import SolverGroup,ParallelPlan,ParallelContext,NativeParallelContext
from salmonpy.material_response import MaterialBinding,ResponseSession

class MyMaterial:
    def __init__(self, scale):
        self.scale = float(scale)
        self.identity = f'example-memory-material/1/scale={self.scale}'
    def prepare(self,context): return self.scale  # Replace with your loader/preparation.

class MyResponse:
    identity='example-memory-response/1'
    execution='collective'
    reduction='sum'  # Every rank returns a local contribution; summed once.
    def create_state(self,initial,point_id,dt,context):
        class State:
            def __init__(self): self.current=np.zeros(3)
            def advance(self,A_old,A_new):
                self.current+=initial*(A_new-A_old)+dt*point_id
                return self.current/context.solver_comm.size
            def close(self): pass
        return State()  # A fresh, independent state for this point.

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('output');p.add_argument('--bridge');p.add_argument('--sha256')
a=p.parse_args()
world=MPI.COMM_WORLD
if world.size!=5: raise ValueError('launch exactly 5 MPI ranks')
if bool(a.bridge)!=bool(a.sha256): raise ValueError('bridge and sha256 must be supplied together')
plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('left',(1,2),point_ids=(10,)),SolverGroup('right',(3,4),point_ids=(20,))))
with ParallelContext.create(world,plan) as context:
    if a.bridge:
        with NativeParallelContext(context,a.bridge,a.sha256) as native:
            assert native.info()==(context.solver_comm.rank,context.solver_comm.size)
            native.barrier()
    bindings={10:MaterialBinding('left-material',MyMaterial(2.),MyResponse()),
              20:MaterialBinding('right-material',MyMaterial(3.),MyResponse())}
    session=ResponseSession(context,bindings,dt=.1)
    history=[]
    for step in range(1,5):
        A=np.array([[point*.01*step,0.,0.] for point in plan.point_ids]) if world.rank==0 else None
        J=session.advance(A)
        if world.rank==0: history.append({'time':step*.1,'A_over_c':A.tolist(),'J_charge':J.tolist()})
    session.close()
# Finish MPI work before root-only file I/O; an I/O error cannot strand peers.
if world.rank==0:
    with Path(a.output).open('x') as f:
        json.dump({'model':'synthetic demonstration, no FDTD','units':'atomic','point_ids':plan.point_ids,'history':history},f,indent=2)
    print(a.output)
