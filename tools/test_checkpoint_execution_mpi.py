"""Four-rank execution: unequal point counts, opaque results and failure stages."""
from mpi4py import MPI
from salmonpy.checkpointing import CheckpointCapability, CheckpointTarget, execute_checkpoint
from salmonpy.errors import InputError, RuntimeFailure

comm=MPI.COMM_WORLD
assert comm.size==4
GROUPS=[{'name':'field','role':'field','ranks':[0],'point_ids':[]},
        {'name':'one','ranks':[1],'point_ids':['a']},
        {'name':'two','ranks':[2,3],'point_ids':['b','c']}]
owner=next(g for g in GROUPS if comm.rank in g['ranks'])
group_comm=comm.Split(0 if comm.rank==0 else 1 if comm.rank==1 else 2,comm.rank)
calls=[];restored={}

def callback(point,operation,*,fail=False):
    def run():
        assert len(set(group_comm.allgather((point,operation))))==1
        calls.append((point,operation))
        if fail and comm.rank==3: raise ValueError('backend read failed')
        if operation=='save': return lambda:(comm.rank,point)  # Must never be MPI-pickled.
        restored[point]=payloads[point]()
    return run

def bindings(operation,*,fail=False,missing=False):
    cap=CheckpointCapability('user-owned/1',('save','restore'),(len(owner['ranks']),),identity='model1')
    if owner.get('role','material')=='field':
        target=CheckpointTarget(cap,save=callback('field','save'),restore=callback('field','restore'))
        return {'field':target}
    targets={}
    for p in owner['point_ids']:
        targets[p]=CheckpointTarget(cap,save=None if missing and comm.rank==3 and p=='c' else callback(p,'save'),
            restore=callback(p,'restore',fail=fail and p=='b'))
    return {'states':targets}

saved=execute_checkpoint(comm,GROUPS,operation='save',**bindings('save'))
payloads=saved['states'] if owner.get('role','material')=='material' else {'field':saved['field']}
assert all(callable(x) for x in payloads.values())
execute_checkpoint(comm,GROUPS,operation='restore',saved=saved['compatibility'],**bindings('restore'))
assert all(restored[p]==(comm.rank,p) for p in payloads)
assert len(calls)==2*(1 if owner.get('role','material')=='field' else len(owner['point_ids']))

calls.clear()
rejected=False
try: execute_checkpoint(comm,GROUPS,operation='save',**bindings('save',missing=True))
except InputError: rejected=True
assert all(comm.allgather(rejected)) and not calls

calls.clear()
failed=False
try: execute_checkpoint(comm,GROUPS,operation='restore',saved=saved['compatibility'],**bindings('restore',fail=True))
except RuntimeFailure as exc:
    assert 'backend read failed' in str(exc)
    failed=True
assert all(comm.allgather(failed))
assert not any(point=='c' for point,operation in calls), 'next point ran after another rank failed'
assert len(calls)==1

group_comm.Free()
assert not MPI.Is_finalized()
if comm.rank==0: print('checkpoint execution MPI PASS: unequal groups, opaque results, exact restore, zero-call admission failure, later slots stopped')
