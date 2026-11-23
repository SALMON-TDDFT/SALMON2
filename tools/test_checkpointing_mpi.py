"""Host-owned MPI checkpoint admission using only the public common API.

mpiexec -n 3 python tools/test_checkpointing_mpi.py
No calculation runtime, SALMON ParallelContext or file I/O; save callbacks are memory-only.
"""
import copy
from mpi4py import MPI
from salmonpy.checkpointing import CheckpointCapability, check_checkpoint_compatibility

comm = MPI.COMM_WORLD
assert comm.size == 3
GROUPS = [{'name': 'field', 'ranks': [0], 'role': 'field', 'point_ids': []},
          {'name': 'user', 'ranks': [1, 2], 'role': 'material', 'point_ids': ['film']}]
field = CheckpointCapability('user-field/1', ('save', 'restore'), (1,), identity='grid-identity')
material = CheckpointCapability('user-current/1', ('save', 'restore'), (2, 1), identity='model-identity')


def arguments():
    return {'field': field if comm.rank == 0 else None,
            'states': {} if comm.rank == 0 else {'film': material}}


saved = check_checkpoint_compatibility(comm, GROUPS, **arguments())
assert check_checkpoint_compatibility(comm, GROUPS, operation='restore', saved=saved, **arguments()) == saved
cases = ('unsupported', 'partition', 'operation', 'saved-disagreement', 'inventory',
         'local-error', 'plan', 'identity', 'missing-declaration', 'empty-local-exception')
for case in cases:
    kw = arguments()
    kw.update(operation='restore', saved=copy.deepcopy(saved))
    groups = copy.deepcopy(GROUPS)
    if case == 'unsupported' and comm.rank == 2:
        kw['states']['film'] = {'format': 'user-current/1', 'modes': ['save'], 'decomposition': [2, 1]}
    if case == 'partition' and comm.rank != 0:
        kw['states']['film'] = CheckpointCapability('user-current/1', ('save', 'restore'), (1, 2), identity='model-identity')
    if case == 'operation' and comm.rank == 0:
        kw.update(operation='save', saved=None)
    if case == 'saved-disagreement' and comm.rank == 2:
        kw['saved']['participants'][1]['states']['film']['identity'] = 'different-record'
    if case == 'inventory':
        kw['required_files'] = (f'state/rank-{comm.rank}.bin',)
        kw['available_files'] = () if comm.rank == 2 else kw['required_files']
    if case == 'local-error' and comm.rank == 1:
        kw['local_error'] = 'user configuration could not be inspected'
    if case == 'plan' and comm.rank == 1:
        groups[1]['ranks'] = [2, 1]
    if case == 'identity' and comm.rank != 0:
        kw['states']['film'] = CheckpointCapability('user-current/1', ('save', 'restore'), (2, 1), identity='changed-model')
    if case == 'missing-declaration' and comm.rank == 2:
        kw['states']['film'] = None
    if case == 'empty-local-exception' and comm.rank == 2:
        class BrokenGroup(dict):
            def __iter__(self): raise ValueError()
        groups[1] = BrokenGroup(groups[1])
    rejected = False
    try:
        check_checkpoint_compatibility(comm, groups, **kw)
    except ValueError:
        rejected = True
    assert all(comm.allgather(rejected)), case
assert not MPI.Is_finalized()
if comm.rank == 0:
    print('common checkpoint API MPI PASS: positive roundtrip; 10 collective rejections; host MPI preserved')


# A single root clock decides the accepted step for every rank. No GS or I/O.
from unittest.mock import patch
from salmonpy.checkpointing import CheckpointSchedule
from salmonpy.errors import RuntimeFailure
with patch('salmonpy.checkpointing.time.monotonic', return_value=0.):
    timer=CheckpointSchedule(interval_steps=2,stop_after_seconds=10.,comm=comm)
saved=[]
with patch('salmonpy.checkpointing.time.monotonic', return_value=1. if comm.rank==0 else 999.):
    assert not timer.after_step(1,saved.append)
    assert not timer.after_step(2,saved.append)
assert saved==[2]
with patch('salmonpy.checkpointing.time.monotonic', return_value=11. if comm.rank==0 else 0.):
    def broken_save(step):
        if comm.rank==2:raise OSError('injected save failure')
    rejected=False
    try:timer.after_step(3,broken_save)
    except RuntimeFailure:rejected=True
    assert all(comm.allgather(rejected))
    assert not timer.stopped and timer.last_saved_step==2
    assert timer.after_step(3,saved.append)
assert all(comm.allgather(saved==[2,3] and timer.stopped))
if comm.rank==0:print('checkpoint schedule MPI PASS: root clock; periodic save; failed save; retry; agreed stop')
