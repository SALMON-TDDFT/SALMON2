"""Prescribed-A TDDFT at two independent macro points (three Python MPI ranks).
mpiexec -n 3 python -m salmonpy.examples.mpi_tddft RUNTIME GS_BUNDLE INPUT NEW_JSON
Requires an audited serial runtime and supported prepared GS/RT bundle.
Two one-rank material groups each retain one serial worker process; no inner
TDDFT MPI or FDTD propagation. Size memory for all two resident electronic states.
"""
import argparse,json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import ParallelPlan,SolverGroup,ParallelContext
from salmonpy.material_response import MaterialBinding,ResponseSession,TDDFTMaterialProvider,TDDFTResponseBackend

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('runtime');p.add_argument('bundle');p.add_argument('input');p.add_argument('output')
p.add_argument('--dt',type=float,default=.02);p.add_argument('--steps',type=int,default=3)
a=p.parse_args()
if MPI.COMM_WORLD.size!=3: raise ValueError('launch exactly three MPI ranks')
if a.steps<1: raise ValueError('positive steps required')
plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('left',(1,),point_ids=(1,)),SolverGroup('right',(2,),point_ids=(2,))))
# MPI failure/partial startup requires job-wide cancellation. This example aborts
# the whole demonstration on errors rather than stranding the other ranks.
try:
    with TDDFTMaterialProvider(a.bundle,a.input,runtime=a.runtime,max_steps=a.steps,threads=1) as provider:
        binding=MaterialBinding('crystal',provider,TDDFTResponseBackend())
        with ParallelContext.create(MPI.COMM_WORLD,plan) as context:
            session=ResponseSession(context,{point:binding for point in plan.point_ids},dt=a.dt)
            history=[]
            try:
                for step in range(1,a.steps+1):
                    A=np.array([[point*step*1e-4,0.,0.] for point in plan.point_ids]) if context.group.role=='field' else None
                    J=session.advance(A)
                    if context.group.role=='field': history.append({'time':step*a.dt,'A_over_c':A.tolist(),'J_charge':J.tolist()})
            finally: session.close()
    if MPI.COMM_WORLD.rank==0:
        with Path(a.output).open('x') as out:
            json.dump({'units':'atomic','point_ids':plan.point_ids,'history':history,
                       'mode':'prescribed A, process-isolated TDDFT; no FDTD'},out,indent=2)
except BaseException:
    import traceback
    traceback.print_exc()
    MPI.COMM_WORLD.Abort(1)
