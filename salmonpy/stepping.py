"""Explicit backend-owned step backups; no MPI initialization or native probing.

Callbacks own numerical state AND side effects. An exception can be recovered
only after every participating rank has returned to the orchestration layer.
"""
from dataclasses import dataclass
from .errors import InputError,RuntimeFailure

__all__=['StepTarget','StepRecoveryError','guarded_step']


@dataclass(frozen=True,eq=False)
class StepTarget:
    """Local owner with capture(), rollback(token), optional release(token).

    capture is read-only; release(None) must clean/idempotently accept a
    failed capture on that rank. Without release, capture must clean its own
    failed preparation and tokens must need no explicit retirement. Tokens
    may be None and are never copied/transferred. Rollback covers accepted clocks,
    latches, caches and external effects; release only retires backup resources.
    Callbacks may use their solver group, never outer orchestration collectives.
    This declaration does not certify hidden state or unsupported native owners.
    """
    key: str
    owner: object
    capture: object
    rollback: object
    release: object = None

    def __post_init__(self):
        if not isinstance(self.key,str) or not self.key.strip():
            raise InputError('nonempty step target key required')
        if self.owner is None or not callable(self.capture) or not callable(self.rollback):
            raise InputError('step target requires owner and capture/rollback callbacks')
        if self.release is not None and not callable(self.release):
            raise InputError('step target release must be callable')


class StepRecoveryError(RuntimeFailure):
    """Failed protected operation; preservation is established by collective vote.

    original_preserved=False forbids retry. Even when True, cleanup_errors need
    host attention before another transaction can use the affected resources.
    """
    def __init__(self,phase,errors,*,original_preserved,rollback_errors=(),cleanup_errors=()):
        self.phase=phase
        self.errors=tuple(errors)
        self.original_preserved=bool(original_preserved)
        self.rollback_errors=tuple(rollback_errors)
        self.cleanup_errors=tuple(cleanup_errors)
        super().__init__(f'protected step failed during {phase}; original_preserved={self.original_preserved}: '+str(self.errors))


def _error(exc):
    return type(exc).__name__+': '+(str(exc) or '<no message>')


def _vote(context,error):
    if context is None:return (error,) if error else ()
    group=context.solver_comm.allgather(error)
    outer=context.leaders.allgather(group) if context.is_leader else None
    outer=context.solver_comm.bcast(outer,root=0)
    return tuple(e for g in outer for e in g if e)


def guarded_step(operation,targets,*,context=None):
    """Run one operation, or restore all explicit targets on a handled exception.

    With ParallelContext every rank calls with the same local key order inside
    its solver group; different groups may have different counts. Operation and
    callbacks must return quiescently on all ranks. No fallback for missing
    capability, process/MPI/allocation recovery or external-effect inference.
    Result and opaque backups never travel through orchestration communicators.
    """
    error=None;items=()
    try:
        if not callable(operation):raise InputError('callable step operation required')
        items=tuple(targets)
        if any(not isinstance(t,StepTarget) for t in items):raise InputError('explicit StepTarget objects required')
        if len({t.key for t in items})!=len(items) or len({id(t.owner) for t in items})!=len(items):
            raise InputError('step targets must have unique keys and owners')
    except Exception as exc:error=_error(exc)
    errors=_vote(context,error)
    if errors:raise StepRecoveryError('admission',errors,original_preserved=True)
    if context is not None:
        orders=context.solver_comm.allgather(tuple((t.key,t.release is not None) for t in items))
        error='target order/release policy differs inside solver group' if any(x!=orders[0] for x in orders) else None
        errors=_vote(context,error)
        if errors:raise StepRecoveryError('admission',errors,original_preserved=True)
        counts=context.leaders.allgather(len(items)) if context.is_leader else None
        counts=context.solver_comm.bcast(counts,root=0)
        slots=max(counts,default=0)
    else:slots=len(items)
    if slots==0:
        raise StepRecoveryError('admission',('at least one explicit step target required',),original_preserved=True)
    captured=[]

    def retire():
        errors=[]
        for slot in reversed(range(slots)):
            failure=None
            if slot<len(captured):
                t,token=captured[slot]
                try:
                    if t.release is not None:t.release(token)
                except Exception as exc:failure=_error(exc)
            errors.extend(_vote(context,failure))
        return tuple(errors)

    for slot in range(slots):
        error=None
        if slot<len(items):
            t=items[slot]
            try:
                captured.append((t,None))
                captured[-1]=(t,t.capture())
            except Exception as exc:error=_error(exc)
        errors=_vote(context,error)
        if errors:
            cleanup=retire()
            raise StepRecoveryError('capture',errors,original_preserved=True,cleanup_errors=cleanup)
    result=None;error=None;cause=None
    try:result=operation()
    except Exception as exc:error=_error(exc);cause=exc
    errors=_vote(context,error)
    if errors:
        rollback=[]
        for slot in reversed(range(slots)):
            failure=None
            if slot<len(captured):
                t,token=captured[slot]
                try:t.rollback(token)
                except Exception as exc:failure=_error(exc)
            rollback.extend(_vote(context,failure))
        cleanup=retire()
        raise StepRecoveryError('operation',errors,original_preserved=not rollback,
                                rollback_errors=rollback,cleanup_errors=cleanup) from cause
    cleanup=retire()
    if cleanup:raise StepRecoveryError('release',cleanup,original_preserved=False,cleanup_errors=cleanup)
    return result
