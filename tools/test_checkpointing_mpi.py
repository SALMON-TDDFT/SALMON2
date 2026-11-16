"""Host-owned MPI checkpoint admission using only the public common API.

mpiexec -n 3 python tools/test_checkpointing_mpi.py
No native runtime, SALMON ParallelContext, persistence callbacks or file I/O.
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
