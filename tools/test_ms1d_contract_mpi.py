"""MPI contract regressions without a native runtime.

PYTHONPATH=. mpiexec -n 3 python tools/test_ms1d_contract_mpi.py
Exercises mixed step/legacy advance callbacks, ordering, preflight rejection,
invalidated failed sessions, and resource-measurement failures independently.
"""
import numpy as np
from mpi4py import MPI
from salmonpy.coupled import CoupledStepper
from salmonpy.distributed import ParallelContext, ParallelPlan, SolverGroup
from salmonpy.material_response import MaterialBinding, _response_timing
from salmonpy.maxwell import Multiscale1D
import salmonpy.material_response as response_module


class Provider:
    identity = 'manufactured-provider/1'
    def __init__(self): self.calls = 0
    def prepare(self, context):
        self.calls += 1
        return None


class CanonicalState:
    def __init__(self, point, dt, fail=False):
        self.point, self.n, self.fail = point, 0, fail
        self.closed = False
        self.timing = _response_timing(dt, 'manufactured')
    def step(self, A):
        if self.fail: raise ValueError('injected backend failure')
        self.n += 1
        return 0.25*A + np.array([0., self.point*1e-4, self.n*1e-4])
    def close(self): self.closed = True


class LegacyState:
    def __init__(self, point, dt):
        self.inner = CanonicalState(point, dt)
        self.timing = self.inner.timing
        self.previous = np.zeros(3)
    def advance(self, old, new):
        np.testing.assert_array_equal(old, self.previous)
        value = self.inner.step(new)
        self.previous = new.copy()
        return value
    def close(self): self.inner.close()


class Backend:
    execution, reduction = 'collective', 'leader'
    def __init__(self, legacy=False, fail=False, wrong_time=False):
        self.legacy, self.fail, self.wrong_time = legacy, fail, wrong_time
        self.identity = f'manufactured/{legacy}/{fail}/{wrong_time}'
    def create_state(self, prepared, point, dt, context):
        state = LegacyState(point, dt) if self.legacy else CanonicalState(point, dt, self.fail)
        if self.wrong_time: state.timing = state.timing.shifted('1/2')
        return state


def expect_collective_error(world, callback, text):
    try: callback()
    except ValueError as exc:
        assert text in str(exc), str(exc)
    else: raise AssertionError('expected collective error')
    assert world.allreduce(1) == 3


def run(mode):
    world = MPI.COMM_WORLD
    # Group order differs from rank order and point ID order.
    plan = ParallelPlan((SolverGroup('left', (1,), point_ids=(19,)),
                        SolverGroup('field', (0,), role='field'),
                        SolverGroup('right', (2,), point_ids=(4,))))
    provider = Provider()
    bindings = {19: MaterialBinding('canonical', provider, Backend(fail=mode == 'backend_failure', wrong_time=mode == 'timing_failure')),
                4: MaterialBinding('legacy', provider, Backend(legacy=True))}
    with ParallelContext.create(world, plan) as context:
        field = None
        if world.rank == 0:
            incident = (lambda t: np.column_stack((np.zeros_like(t), np.zeros_like(t), np.ones_like(t)))) if mode == 'initial_A' else None
            field = Multiscale1D(2, spacing=.1, vacuum=(3, 2), incident=incident).new_state(dt=.0005, steps=4)
        if mode == 'initial_A':
            expect_collective_error(world, lambda: CoupledStepper.from_multiscale(context, bindings, field), 'vacuum')
            assert world.allreduce(provider.calls) == 0, 'preflight must precede material callbacks'
            return None
        if mode == 'timing_failure':
            expect_collective_error(world, lambda: CoupledStepper.from_multiscale(context, bindings, field), 'point 19')
            assert context._dependents == 0, 'rejected timing must release the response'
            return None
        coupled = CoupledStepper.from_multiscale(context, bindings, field)
        frames = []
        old_clock = response_module.time.perf_counter
        try:
            if mode == 'clock_failure':
                def broken_clock(): raise OSError('injected resource clock failure')
                response_module.time.perf_counter = broken_clock
            if mode == 'backend_failure':
                expect_collective_error(world, coupled.step, 'injected backend failure')
                assert coupled.response.failed
                expect_collective_error(world, coupled.step, 'not usable')
                return None
            for n in range(1, 5):
                frame = coupled.step()
                if frame is not None:
                    expected = .25*frame['A_over_c'] + np.array([[0., 4e-4, n*1e-4], [0., 19e-4, n*1e-4]])
                    np.testing.assert_allclose(frame['J_charge'], expected, rtol=0, atol=1e-18)
                    np.testing.assert_array_equal(frame['E_at_points'], field.grid_E[3:5])
                    frames.append(np.concatenate((field.grid_A.ravel(), frame['J_charge'].ravel())))
            report = coupled.response.communication_report()
            assert report['steps'] == 4
            for rank in report['ranks']:
                assert all(v == 4 for v in rank['calls'].values())
                if mode == 'clock_failure':
                    assert all(v is None for v in rank['seconds'].values())
                    assert rank['measurement_errors']
                else: assert all(v >= 0 for v in rank['seconds'].values())
            expect_collective_error(world, coupled.step, 'steps')
            assert coupled.response.failed
            expect_collective_error(world, coupled.step, 'not usable')
            return np.asarray(frames) if world.rank == 0 else None
        finally:
            response_module.time.perf_counter = old_clock
            coupled.response.close()
            if coupled.field is not None: coupled.field.close()


def main():
    world = MPI.COMM_WORLD
    if world.size != 3: raise ValueError('exactly 3 ranks required')
    run('initial_A')
    run('timing_failure')
    normal = run('normal')
    broken = run('clock_failure')
    if world.rank == 0: np.testing.assert_array_equal(normal, broken)
    run('backend_failure')
    assert world.allreduce(1) == 3
    if world.rank == 0: print('PASS: mixed callbacks/order, initial A/timing preflight, capacity/backend failure, measurement isolation', flush=True)


if __name__ == '__main__':
    try: main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
