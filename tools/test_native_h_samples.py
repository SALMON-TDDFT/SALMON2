"""2/3 ranks: raw H sampling vs native saved centered-H observations."""
import sys,json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.coupled import NativeEHField
from salmonpy.distributed import ParallelContext,ParallelPlan,SolverGroup
from salmonpy.geometry import Grid,Geometry
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.errors import InputError

world=MPI.COMM_WORLD
build,output=sys.argv[1:3];output=Path(output)
stride=int(sys.argv[3]) if len(sys.argv)>3 else 1
if stride not in (1,2,64):raise ValueError('stride must be 1, 2 or 64')
nf=world.size-1
if nf not in (1,2):raise ValueError('2/3 ranks required')
try:
    dt=.003;steps=64;grid=Grid((40,40,40),(2.,)*3)
    points={10:(19,20,20),20:(22,20,20)}
    sites=list(points.values())
    observations=tuple(tuple(grid.coordinates[d][q[d]] for d in range(3)) for q in sites)
    calculation=FDTD(Geometry(grid),(),dt,steps,source=PlaneWave(AcosPulse(0.,20.,dt*steps),0,-16.),observations=observations,observation_stride=stride)
    plan=ParallelPlan((SolverGroup('field',tuple(range(nf)),role='field'),SolverGroup('labels',(nf,),point_ids=(10,20))))
    if world.rank==0:output.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    with ParallelContext.create(world,plan) as context:
        if context.group.role=='field':
            with NativeEHField(context,calculation,build,output/'field',points=points,grid=(nf,1,1)) as field:
                H=[np.zeros((2,3))];E=[field.E.copy()];A=[field.A.copy()]
                for n in range(1,steps+1):
                    j=np.zeros((2,3))
                    if n==2:j[0,2]=1.e-4
                    field.advance(j)
                    before=(field.A.copy(),field.E.copy(),field.step_index)
                    if n==1 and nf==2:
                        original_live=field.native._live
                        if context.solver_comm.rank==1:
                            def fail_live():raise RuntimeError('injected native liveness rejection')
                            field.native._live=fail_live
                        try:
                            try:field.sample_magnetic_field()
                            except InputError as exc:assert 'injected native liveness rejection' in str(exc)
                            else:raise AssertionError('rank-local liveness failure accepted')
                        finally:field.native._live=original_live
                    sample=field.sample_magnetic_field()
                    repeated=field.sample_magnetic_field()
                    np.testing.assert_array_equal(sample.values,repeated.values)
                    assert sample.times[0]==(n+.5)*dt
                    np.testing.assert_array_equal(field.A,before[0]);np.testing.assert_array_equal(field.E,before[1])
                    assert field.step_index==before[2]
                    H.append(sample.values[0]);E.append(field.E.copy());A.append(field.A.copy())
            if world.rank==0:
                H=np.array(H);error=None
                if stride:
                    expected=np.stack([np.atleast_2d(np.loadtxt(output/'field'/f'obs{i+1}_at_point_rt.data'))[:,4:7]
                                       for i in range(len(sites))],axis=1)
                    centered=(.5*(H[:-1]+H[1:]))[stride-1::stride]
                    np.testing.assert_allclose(centered,expected,rtol=1.e-12,atol=1.e-18)
                    error=float(np.max(np.abs(centered-expected)))
                assert np.max(np.abs(H))>0
                np.savez(output/'samples.npz',H=H,E=E,A=A)
                (output/'report.json').write_text(json.dumps({'field_ranks':nf,'centered_H_max_error':error,'observation_stride':stride,
                    'nonmutating_repeated_read':True,'time':'raw H at (n+1/2)dt'},indent=2))
    world.Barrier()
    if world.rank==0:print('PASS native raw H vs saved centered H',nf)
except BaseException:
    import traceback
    traceback.print_exc();world.Abort(1)
