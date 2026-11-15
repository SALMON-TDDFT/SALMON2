"""MPI fault worker: empty-message callback errors must stop all three ranks."""
import argparse
import numpy as np
from mpi4py import MPI
from salmonpy.examples.mpi_recorded_material import RecordedExample
from salmonpy.materials import bind_recorded_material, response_identity
from salmonpy.material_response import ResponseSession
from salmonpy.distributed import ParallelContext, ParallelPlan, SolverGroup

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('failure', choices=('create', 'step'))
args = parser.parse_args()
world = MPI.COMM_WORLD
assert world.size == 3

class Broken(RecordedExample):
    def new_state(self, **kwargs):
        if args.failure == 'create' and world.rank == 1:
            raise RuntimeError()
        state = super().new_state(**kwargs)
        if args.failure == 'step' and world.rank == 1:
            def fail(A): raise RuntimeError()
            state.step = fail
        return state

plan = ParallelPlan((SolverGroup('field', (0,), role='field'),
                     SolverGroup('one', (1,), point_ids=(10,)),
                     SolverGroup('two', (2,), point_ids=(20,))))
with ParallelContext.create(world, plan) as context:
    binding = bind_recorded_material(Broken(), key='custom')
    bindings = {10:binding, 20:binding}
    mapping = {10:'one', 20:'two'}
    record = response_identity(bindings, mapping, run_id='fault')
    session = None
    rejected = False
    try:
        session = ResponseSession(context, bindings, dt=.02, point_record=record, point_map=mapping)
        session.step(np.zeros((2,3)) if world.rank == 0 else None)
    except Exception as exc:
        assert 'RuntimeError' in str(exc), repr(exc)
        rejected = True
    finally:
        if session is not None: session.close()
    assert world.allreduce(int(rejected)) == 3
if world.rank == 0: print('PASS collective empty-message', args.failure)
