"""Four-rank protected restore with field/custom points and an injected late failure."""
from pathlib import Path
import sys
import hashlib
import numpy as np
from mpi4py import MPI
from salmonpy.checkpointing import (CheckpointCapability, CheckpointTarget, CheckpointState,
    execute_checkpoint, restore_checkpoint, AtomicRestoreError)
from salmonpy.examples.custom_checkpoint import UserCurrent

comm=MPI.COMM_WORLD
assert comm.size==4
root=Path(sys.argv[1]).resolve()
if comm.rank==0: root.mkdir(parents=True,exist_ok=False)
comm.Barrier()
groups=[{'name':'field','role':'field','ranks':[0],'point_ids':[]},
        {'name':'one','ranks':[1],'point_ids':['a']},
        {'name':'two','ranks':[2,3],'point_ids':['b','c']}]
owner=next(g for g in groups if comm.rank in g['ranks'])
field=owner.get('role')=='field'
points=['field'] if field else owner['point_ids']
sub=comm.Split(0 if comm.rank==0 else 1 if comm.rank==1 else 2,comm.rank)
wave=np.array([[.001*np.sin(.2*n),0.,0.] for n in range(1,33)])
decays={'field':.99,'a':.97,'b':.95,'c':.92}
checkpoint_states={};originals={};caps={};hashes={};cleanup_calls=[]
for p in points:
    checkpoint_states[p]=UserCurrent(decays[p])
    originals[p]=UserCurrent(decays[p])
    for A in wave[:17]: checkpoint_states[p].step(A)
    for A in wave[:5]: originals[p].step(A)
    caps[p]=CheckpointCapability('user-json-mpi/1',('save','restore'),(len(owner['ranks']),),
        identity=checkpoint_states[p].checkpoint_capability.identity)
live=CheckpointState(field=originals['field'] if field else None,states={} if field else originals)

def path(p): return root/f'rank-{comm.rank}-{p}.json'
def arguments(targets): return {'field':targets['field']} if field else {'states':targets}
def save(p):
    assert len(set(sub.allgather(('save',p))))==1
    checkpoint_states[p].save(path(p))
saved=execute_checkpoint(comm,groups,operation='save',**arguments({p:CheckpointTarget(caps[p],save=lambda p=p:save(p)) for p in points}))
for p in points: hashes[p]=hashlib.sha256(path(p).read_bytes()).hexdigest()

before={p:(s.A.tobytes(),s.current.tobytes(),s.steps,s.load_calls) for p,s in originals.items()}
def stage(p,fail):
    assert len(set(sub.allgather(('stage',p))))==1
    assert hashlib.sha256(path(p).read_bytes()).hexdigest()==hashes[p]
    candidate=UserCurrent(decays[p]);candidate.load(path(p))
    if fail and comm.rank==3 and p=='c':
        candidate.current[:]=12345.  # Partial candidate mutation, never the live state.
        raise OSError('late candidate load failure')
    return candidate

def discard(p,candidate):
    peers=sub.allgather(('discard',p,candidate is None))
    assert all(item[:2]==('discard',p) for item in peers)
    cleanup_calls.append((p,candidate is None))

def targets(fail):
    return {p:CheckpointTarget(caps[p],stage_restore=lambda p=p:stage(p,fail),
        discard_restore=lambda c,p=p:discard(p,c)) for p in points}

failed=False
try:
    restore_checkpoint(comm,groups,live=live,saved=saved['compatibility'],**arguments(targets(True)))
except AtomicRestoreError as exc:
    assert exc.original_preserved and not exc.cleanup_errors
    failed=True
assert all(comm.allgather(failed))
for p,old in originals.items():
    actual=live.field if field else live.states[p]
    assert actual is old
    assert (old.A.tobytes(),old.current.tobytes(),old.steps,old.load_calls)==before[p]
    reference=UserCurrent(decays[p])
    for A in wave[:6]: expected=reference.step(A)
    np.testing.assert_array_equal(actual.step(wave[5]),expected)
assert len(cleanup_calls)==len(points)
if comm.rank==3: assert ('c',True) in cleanup_calls

# Direct-only backends must be refused, not used as a destructive fallback.
calls=[]
unsupported={p:CheckpointTarget(caps[p],restore=lambda:calls.append('LOAD')) for p in points}
rejected=False
try: restore_checkpoint(comm,groups,live=live,saved=saved['compatibility'],**arguments(unsupported))
except ValueError: rejected=True
assert all(comm.allgather(rejected)) and not calls

# Policy or cleanup-hook disagreement must be rejected before candidates run.
from salmonpy.checkpointing import execute_checkpoint as execute
mixed={p:CheckpointTarget(caps[p],restore=lambda:calls.append('direct'),
    stage_restore=lambda p=p:stage(p,False),discard_restore=lambda c,p=p:discard(p,c)) for p in points}
rejected=False
try:
    execute(comm,groups,operation='restore',atomic=comm.rank!=0,live=live if comm.rank!=0 else None,
        saved=saved['compatibility'],**arguments(mixed))
except ValueError:rejected=True
assert all(comm.allgather(rejected)) and not calls
mismatched=targets(False)
if comm.rank==3:
    mismatched['c']=CheckpointTarget(caps['c'],stage_restore=lambda:stage('c',False))
rejected=False
try:restore_checkpoint(comm,groups,live=live,saved=saved['compatibility'],**arguments(mismatched))
except ValueError:rejected=True
assert all(comm.allgather(rejected))

run=restore_checkpoint(comm,groups,live=live,saved=saved['compatibility'],**arguments(targets(False)))
for p,old in originals.items():
    current=live.field if field else live.states[p]
    assert current is not old and current.steps==17 and old.steps==6
    assert (run['previous']['field'] if field else run['previous']['states'][p]) is old
    reference=UserCurrent(decays[p]);expected=np.array([reference.step(A) for A in wave])
    actual=np.array([current.step(A) for A in wave[17:]])
    np.testing.assert_array_equal(actual,expected[17:])
    assert not np.shares_memory(current.current,old.current)
sub.Free()
assert not MPI.Is_finalized()
if comm.rank==0: print('atomic restore MPI PASS: all originals preserved after late rank failure; old evolution continues; successful retry exact; direct fallback refused')
