"""Native E/H FDTD coupled to two internally MPI-parallel TDDFT points.
mpiexec -n 6 python -m salmonpy.examples.mpi_coupled_tddft BUILD_JSON GS_BUNDLE RT_INPUT NEW_OUTPUT
One native lifetime/process; host-MPI library and shared writable output required.
This small finite two-cell geometry is a connection example, not converged optics.
For production TDDFT prefer 1D macro propagation; recommend 3D electronic feedback
with SBE/custom current solvers. This sample is not a 1D reduction.
Set OpenMP/BLAS threads in the launch environment; example budget assumes one.
"""
import argparse,json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.geometry import Grid,Geometry
from salmonpy.fields import AcosPulse
from salmonpy.distributed import ParallelPlan,SolverGroup,ParallelContext
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider,NativeTDDFTBackend
from salmonpy.coupled import NativeEHField,CoupledStepper

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('build_record');p.add_argument('bundle');p.add_argument('input');p.add_argument('output')
p.add_argument('--field-ranks',type=int,choices=(1,2),default=2)
p.add_argument('--material-ranks',type=int,choices=(1,2),default=2)
p.add_argument('--steps',type=int,default=64)
a=p.parse_args();world=MPI.COMM_WORLD
try:
    nf=a.field_ranks;nm=a.material_ranks
    if world.size!=nf+2*nm:raise ValueError('launch field-ranks + 2*material-ranks MPI ranks')
    if a.steps<1:raise ValueError('positive steps required')
    output=Path(a.output).resolve()
    if world.rank==0:
        output.mkdir(parents=True,exist_ok=False)
        print(f'Field ranks={nf}; two material groups of {nm}; two electronic states.')
        print('Memory unmeasured. Use a small supported GS first; no automatic memory limit.')
    world.Barrier()
    plan=ParallelPlan((SolverGroup('field',tuple(range(nf)),role='field'),
        SolverGroup('left',tuple(range(nf,nf+nm)),point_ids=(1,)),
        SolverGroup('right',tuple(range(nf+nm,nf+2*nm)),point_ids=(2,))))
    context=ParallelContext.create(world,plan)
    dt=.003;grid=Grid((40,40,40),(2.,2.,2.))
    points={1:(19,20,20),2:(20,20,20)}
    calculation=FDTD(Geometry(grid),(),dt,a.steps,
        source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.),observations=((-1.,1.,1.),(1.,1.,1.)))
    field=None
    if context.group.role=='field':
        field=NativeEHField(context,calculation,a.build_record,output/'field',points=points,grid=(nf,1,1))
    # Each material rank enters the same native electronic state on its group.
    # Replace these provider/backend objects with a user implementation to change physics.
    provider=NativeTDDFTMaterialProvider(a.bundle,a.input,max_steps=a.steps)
    backend=NativeTDDFTBackend(a.build_record,output,decomposition=(1,nm,1,1,1))
    binding=MaterialBinding('crystal',provider,backend)
    response=ResponseSession(context,{1:binding,2:binding},dt=dt)
    coupled=CoupledStepper(response,field);history=[]
    for _ in range(a.steps):
        frame=coupled.step()
        if world.rank==0:
            history.append({k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in frame.items()})
    response.close()
    if field is not None:field.close()
    context.close()
    if world.rank==0:
        (output/'response.json').write_text(json.dumps({'scheme':NativeEHField.scheme,'units':'atomic',
            'point_ids':plan.point_ids,'history':history},indent=2)+'\n')
        print(output/'response.json')
except BaseException:
    # Native failures cannot be rolled back collectively; do not enter context
    # cleanup collectives on only the failing rank before stopping this sample job.
    import traceback
    traceback.print_exc();world.Abort(1)
