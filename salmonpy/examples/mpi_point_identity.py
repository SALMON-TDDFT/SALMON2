"""Synthetic A/c response with persistent user IDs (3 or 4 MPI ranks).

mpiexec -n 3 python -m salmonpy.examples.mpi_point_identity NEW_OUTPUT
One field rank and two material groups; four ranks exercise subgroup reduction.
"""
import argparse
import hashlib
from pathlib import Path
import numpy as np


def run(directory, *, scenario='success'):
    from mpi4py import MPI
    from salmonpy.distributed import ParallelContext, ParallelPlan, SolverGroup
    from salmonpy.material_response import MaterialBinding, ResponseSession
    from salmonpy.errors import InputError
    from salmonpy.point_runs import open_point_arrays
    world = MPI.COMM_WORLD
    if world.size not in (3, 4):
        raise InputError('launch 3 or 4 MPI ranks')
    plan = ParallelPlan((SolverGroup('field', (0,), role='field'),
                         SolverGroup('group-A', tuple(range(1, world.size-1)), point_ids=(10,)),
                         SolverGroup('group-B', (world.size-1,), point_ids=(90,))))
    mapping = {10: 'right/α', 90: 'left/β'}
    entries = [{'point_id': mapping[p], 'material_id': f'film-{p}',
                'definition_digest': hashlib.sha256(f'synthetic-scale={p}'.encode()).hexdigest(),
                'initial_state_id': f'zero-{p}', 'initial_state_digest': hashlib.sha256(b'zero').hexdigest()}
               for p in plan.point_ids]
    record = {'schema': 'salmonpy.point-run/1', 'run_id': 'mpi-example', 'points': entries}
    if scenario == 'mismatch' and world.rank == world.size-1:
        record['run_id'] = 'different-run'
    if scenario == 'mismatch_internal' and world.rank == 2:
        record['run_id'] = 'different-within-group'
    if scenario == 'malformed' and world.rank == world.size-1:
        mapping[90] = mapping[10]
    calls = []

    class Provider:
        def __init__(self, point):
            self.point = point
            self.identity = f'synthetic-provider/{point}'
        def prepare(self, context):
            calls.append('prepare')
            if context.point_records is not None:
                entry = context.point_records[self.point]
                assert entry['point_id'] == mapping[self.point]
                assert entry['initial_state_id'] == f'zero-{self.point}'
                try:
                    entry['material_id'] = 'mutation'
                except TypeError:
                    pass
                else:
                    raise AssertionError('mutable callback metadata')
            return self.point

    class Backend:
        identity = 'synthetic-linear-response/1'
        execution = 'collective'
        reduction = 'sum'
        def create_state(self, initial, point_id, dt, context):
            class State:
                current = np.zeros(3)
                def step(self, A):
                    self.current = initial * A / context.solver_comm.size
                    return self.current.copy()
                def close(self):
                    pass
            return State()

    with ParallelContext.create(world, plan) as context:
        bindings = {p: MaterialBinding(f'cache-{p}', Provider(p), Backend()) for p in plan.point_ids}
        if scenario in ('mismatch', 'mismatch_internal', 'malformed'):
            rejected = False
            try:
                ResponseSession(context, bindings, dt=.1, point_record=record, point_map=mapping)
            except InputError:
                rejected = True
            assert world.allreduce(int(rejected), op=MPI.SUM) == world.size
            assert not any(world.allgather(bool(calls))), 'callbacks must not start on metadata disagreement'
            return
        with_identity = ResponseSession(context, bindings, dt=.1, point_record=record, point_map=mapping)
        currents = []
        drives = []
        try:
            for n in range(1, 5):
                A = np.array([[.01*n, 0, 0], [0, .02*n, 0]]) if world.rank == 0 else None
                current = with_identity.step(A)
                if world.rank == 0:
                    np.testing.assert_allclose(current, np.array([10, 90])[:, None]*A, atol=1e-14, rtol=0)
                    currents.append(current.copy())
                    drives.append(A.copy())
            if scenario in ('writefailure', 'ioerror'):
                # Preflight failure or genuine FileExistsError, collectively propagated.
                if scenario == 'ioerror' and world.rank == 0:
                    Path(directory).mkdir()
                try:
                    with_identity.save_arrays(directory if world.rank == 0 else None,
                        {'current': currents} if world.rank == 0 else None,
                        point_axes={'current': 0 if scenario == 'writefailure' else 1} if world.rank == 0 else None)
                except InputError:
                    rejected = True
                else:
                    rejected = False
                assert world.allreduce(int(rejected), op=MPI.SUM) == world.size
            else:
                archive = with_identity.save_arrays(directory if world.rank == 0 else None,
                    {'current': currents, 'A_over_c': drives, 'time': np.arange(1, 5)*.1} if world.rank == 0 else None,
                    point_axes={'current': 1, 'A_over_c': 1} if world.rank == 0 else None,
                    metadata={'units': 'atomic', 'current': 'charge j', 'scope': 'synthetic regression'} if world.rank == 0 else None)
                if world.rank == 0:
                    checked = open_point_arrays(archive.directory)
                    assert checked.metadata['point_run'] == record
                    np.testing.assert_array_equal(checked.array('current'), currents)
        finally:
            with_identity.close()
        legacy = ResponseSession(context, bindings, dt=.1)
        try:
            for index in range(4):
                current = legacy.step(drives[index] if world.rank == 0 else None)
                if world.rank == 0:
                    np.testing.assert_array_equal(current, currents[index])
        finally:
            legacy.close()
    if world.rank == 0:
        print('PASS identity, callback metadata, direct/legacy current and output:', scenario, directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    run(args.directory)


if __name__ == '__main__':
    main()
