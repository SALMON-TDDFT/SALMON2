"""One-way native dipole radiation: 2/3 total ranks, 1/2 field ranks.

Arguments: BUILD_JSON NEW_OUTPUT [ISOLATED_WORKER_DIPOLES_NPZ]
The optional NPZ is produced by test_isolated_worker.py. Source trace is delayed
one dt by explicitly repeating p0; no response is fed back to the molecule.
"""
import json
from pathlib import Path
import sys
import numpy as np
from mpi4py import MPI
from salmonpy.dipole_source import DipoleMap, drive_native_dipoles
from salmonpy.fdtd import FDTD, PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid, Geometry
from salmonpy.distributed import ParallelPlan, SolverGroup, ParallelContext
from salmonpy.coupled import NativeEHField
from salmonpy.errors import InputError

world = MPI.COMM_WORLD
record, output = sys.argv[1:3]
output = Path(output)
field_size = world.size-1
if field_size not in (1,2):
    raise ValueError('use 2 or 3 total ranks')
try:
    dt = .003
    spacing = 2.
    if len(sys.argv) == 4:
        with np.load(sys.argv[3], allow_pickle=False) as a:
            dt = float(a['dt'])
            # Columns 3:6 are total charge dipole, not electronic matter current.
            original = a['values'][:,3:6]
        p = np.concatenate((original[:1], original), axis=0)[:,None,:]
        spacing = 20.
    else:
        p = np.zeros((65,1,3))
        # Static displacement after one nonzero current interval.
        p[2:,0,2] = 1.e-4*dt*spacing**3
    steps = len(p)-1
    grid = Grid((40,40,40),(spacing,)*3)
    points = {10:(19,20,20),20:(22,20,20)}
    observations = tuple(tuple(grid.coordinates[d][q[d]] for d in range(3)) for q in points.values())
    calculation = FDTD(Geometry(grid),(),dt,steps,
                       source=PlaneWave(AcosPulse(0.,20.,dt*steps),0,-8*spacing),observations=observations)
    mapping = DipoleMap(('molecule',),(10,20),[[1.],[0.]],[spacing**3]*2)
    plan = ParallelPlan((SolverGroup('field',tuple(range(field_size)),role='field'),
                         SolverGroup('labels-only',(field_size,),point_ids=(10,20))))
    if world.rank == 0:
        output.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    with ParallelContext.create(world,plan) as context:
        if context.group.role == 'field':
            with NativeEHField(context,calculation,record,output/'field',points=points,grid=(field_size,1,1)) as field:
                bad = p.copy()
                bad[1,0,0] += 1
                try:
                    drive_native_dipoles(field,mapping,bad,dt=dt)
                except InputError:
                    pass
                else:
                    raise AssertionError('nonzero first interval accepted')
                assert field.step_index == 0
                if field_size == 2:
                    different = p.copy()
                    if context.solver_comm.rank == 1:
                        different[-1,0,0] += 1
                    try:
                        drive_native_dipoles(field,mapping,different,dt=dt)
                    except InputError:
                        pass
                    else:
                        raise AssertionError('rank-local trajectory mismatch accepted')
                    assert field.step_index == 0
                result = drive_native_dipoles(field,mapping,p,dt=dt)
                samples=field.sample_fields()
                assert set(samples)=={'A','E'}
                assert samples['A'].times[0]==(steps+.5)*dt
                assert samples['E'].times[0]==steps*dt
                np.testing.assert_array_equal(samples['E'].values[0],result['E'][-1])
                np.testing.assert_array_equal(samples['A'].values[0],field.A)
                np.testing.assert_array_equal(samples['E'].positions,observations)
                sample_records=context.solver_comm.allgather(samples['E'].values.tolist())
                assert all(r==sample_records[0] for r in sample_records)
                if len(sys.argv) == 3:
                    expected = np.array([0,0,-4*np.pi*dt*.5e-4])
                    np.testing.assert_allclose(result['E'][2,0],expected,rtol=2.e-13,atol=1.e-18)
                remote = float(np.max(np.abs(result['E'][:,1])))
                if remote <= 1.e-25:
                    raise AssertionError('no radiated field at remote point')
                if world.rank == 0:
                    np.savez(output/'radiation.npz',**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
                    (output/'report.json').write_text(json.dumps(dict(field_ranks=field_size,
                        source='recorded molecular charge dipole' if len(sys.argv)==4 else 'synthetic impulse',
                        steps=steps, dt=dt, explicit_trace_delay=dt, remote_E_max=remote,
                        initial_current_rejected=True, one_way=True),indent=2)+'\n')
    world.Barrier()
    if world.rank == 0:
        print('dipole radiation PASS',field_size,'field ranks')
except BaseException:
    import traceback
    traceback.print_exc()
    world.Abort(1)
