"""Fixed-geometry dipole field and source-work refinement (2 MPI ranks).

Usage: mpiexec -n 2 python tools/test_dipole_refinement_mpi.py BUILD_JSON OUT DX DT
One field rank; one rank supplies point labels. Explicit external output, no build.
"""
import json
from pathlib import Path
import sys
import numpy as np

C = 137.035999084
AMPLITUDE = 1.e-3
START = .03
DURATION = .12
DISTANCE = 6.
END = .30


def pulse(t):
    t = np.asarray(t, dtype=float)
    x = np.pi*(t-START)/DURATION
    inside = (t > START) & (t < START+DURATION)
    x = np.where(inside, x, 0.)
    s, co = np.sin(x), np.cos(x)
    a = np.pi/DURATION
    return (AMPLITUDE*s**4, 4*AMPLITUDE*a*s**3*co,
            4*AMPLITUDE*a*a*(3*s*s*co*co-s**4))


def reference(t):
    p, v, a = pulse(np.asarray(t)-DISTANCE/C)
    return -(p/DISTANCE**3+v/(C*DISTANCE**2)+a/(C*C*DISTANCE))



def check_reference():
    t = np.linspace(START,START+DURATION,10001)
    p,v,a = pulse(t)
    h = t[1]-t[0]
    np.testing.assert_allclose(np.gradient(p,h)[2:-2],v[2:-2],atol=1.e-8,rtol=1.e-6)
    np.testing.assert_allclose(np.gradient(v,h)[2:-2],a[2:-2],atol=1.e-6,rtol=1.e-5)
    numerical = 2/(3*C**3)*np.trapezoid(a*a,t)
    exact = 8*np.pi**4*AMPLITUDE**2/(3*C**3*DURATION**3)
    np.testing.assert_allclose(numerical,exact,rtol=1.e-13)
    return float(exact)


def main():
    from mpi4py import MPI
    from salmonpy.dipole_source import DipoleMap, drive_native_dipoles
    from salmonpy.fdtd import FDTD, PlaneWave
    from salmonpy.fields import AcosPulse
    from salmonpy.geometry import Grid, Geometry
    from salmonpy.distributed import ParallelPlan, SolverGroup, ParallelContext
    from salmonpy.coupled import NativeEHField
    check_reference()
    world = MPI.COMM_WORLD
    if world.size != 2:
        raise ValueError('exactly two ranks required')
    build, out, dx, dt = sys.argv[1:]
    out = Path(out)
    dx, dt = float(dx), float(dt)
    if not np.isfinite(dx+dt) or dx < .75 or dx > 2 or dt < .00075 or dt > .003:
        raise ValueError('bounded study requires .75 <= dx <= 2 and .00075 <= dt <= .003')
    n = 2*int(np.ceil(40/dx))+1
    center = n//2
    offset = int(round(DISTANCE/dx))
    if not np.isclose(offset*dx,DISTANCE,atol=1.e-12,rtol=0):
        raise ValueError('observation distance must be integer number of cells')
    steps = int(round(END/dt))
    if not np.isclose(steps*dt,END,atol=1.e-12,rtol=0):
        raise ValueError('END must be integer number of steps')
    t = np.arange(steps+1)*dt
    dipoles = np.zeros((steps+1,1,3))
    dipoles[:,0,2] = pulse(t)[0]
    grid = Grid((n,)*3,(dx,)*3)
    points = {10:(center,)*3,20:(center+offset,center,center)}
    obs = tuple(tuple(grid.coordinates[d][q[d]] for d in range(3)) for q in points.values())
    calc = FDTD(Geometry(grid),(),dt,steps,source=PlaneWave(AcosPulse(0.,20.,END),0,-20.),observations=obs)
    mapping = DipoleMap(('analytic-charge-dipole',),(10,20),[[1.],[0.]],[dx**3]*2)
    plan = ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('labels',(1,),point_ids=(10,20))))
    if world.rank == 0:
        out.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    with ParallelContext.create(world,plan) as context:
        if context.group.role == 'field':
            with NativeEHField(context,calc,build,out/'field',points=points,grid=(1,1,1)) as field:
                result = drive_native_dipoles(field,mapping,dipoles,dt=dt)
            expected = reference(t)
            actual = result['E'][:,1,2]
            error = float(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
            work = float(-dt*np.sum(result['j']*(result['E'][:-1]+result['E'][1:])*.5)*dx**3)
            analytic_energy = 8*np.pi**4*AMPLITUDE**2/(3*C**3*DURATION**3)
            if not np.isfinite(error+work) or work <= 0:
                raise AssertionError('nonfinite field error or nonpositive net source work')
            report = dict(dx=dx,dt=dt,counts=[n]*3,steps=steps,source_position=list(obs[0]),
                          observation_position=list(obs[1]),relative_L2_E=error,source_work_Ha=work,
                          continuum_radiation_energy_Ha=analytic_energy,
                          source_work_relative_error=abs(work/analytic_energy-1),
                          pulse_start=START,pulse_duration=DURATION,observation_end=END,
                          source='compact sin^4 physical charge dipole',self_feedback=False)
            np.savez(out/'comparison.npz',time=t,expected_Ez=expected,actual_Ez=actual,
                     E=result['E'],j=result['j'],dipoles=dipoles)
            (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps(report,indent=2),flush=True)
    world.Barrier()


if __name__ == '__main__':
    if sys.argv[1:] == ['--check-reference']:
        print('analytic derivatives and radiation integral PASS',check_reference())
        raise SystemExit(0)
    from mpi4py import MPI
    try:
        main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
