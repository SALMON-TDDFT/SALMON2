"""Explicit two-level MPI ownership. Importing this module does not initialize MPI.

All create/response calls are collective at their documented scope. Native MPI
process failure is not recoverable here; use a launcher timeout/job cancellation.
"""
from dataclasses import dataclass
import threading


def _ids(values, name):
    values=tuple(values)
    if any(type(i) is not int or i<0 for i in values) or len(set(values))!=len(values):
        raise ValueError(name+' must contain distinct nonnegative integers')
    return values


@dataclass(frozen=True)
class SolverGroup:
    name: str
    ranks: tuple
    point_ids: tuple=()
    role: str='material'

    def __post_init__(self):
        if not isinstance(self.name,str) or not self.name: raise ValueError('group name required')
        object.__setattr__(self,'ranks',_ids(self.ranks,'ranks'))
        object.__setattr__(self,'point_ids',tuple(sorted(_ids(self.point_ids,'point_ids'))))
        if not self.ranks or self.role not in ('field','material'): raise ValueError('invalid group')
        if self.role=='field' and self.point_ids: raise ValueError('field group cannot own material points')
        if self.role=='material' and not self.point_ids: raise ValueError('material points required')


@dataclass(frozen=True)
class ParallelPlan:
    groups: tuple

    def __post_init__(self):
        object.__setattr__(self,'groups',tuple(self.groups))

    @property
    def point_ids(self):
        return tuple(sorted(p for g in self.groups for p in g.point_ids))

    def validate(self, world_size):
        if not self.groups or any(not isinstance(g,SolverGroup) for g in self.groups):
            raise ValueError('SolverGroup entries required')
        if len({g.name for g in self.groups})!=len(self.groups): raise ValueError('duplicate group name')
        if sum(g.role=='field' for g in self.groups)!=1: raise ValueError('exactly one field group required')
        ranks=[r for g in self.groups for r in g.ranks]
        if sorted(ranks)!=list(range(world_size)): raise ValueError('each world rank must belong to exactly one group')
        if len(set(self.point_ids))!=len(self.point_ids): raise ValueError('each point must have exactly one owner')


class ParallelContext:
    @classmethod
    def create(cls, world, plan):
        """Collective on world. Validate agreement before creating communicators."""
        from mpi4py import MPI
        if threading.current_thread() is not threading.main_thread(): raise RuntimeError('MPI calls require main thread')
        if MPI.Is_finalized() or world==MPI.COMM_NULL or world.Is_inter(): raise ValueError('live intracommunicator required')
        error=None
        try: plan.validate(world.size)
        except Exception as exc: error=str(exc)
        entries=world.allgather((error,plan if error is None else None))
        if any(e for e,p in entries): raise ValueError('invalid parallel plan: '+str([e for e,p in entries if e]))
        if any(p!=plan for e,p in entries): raise ValueError('parallel plans differ across ranks')
        self=cls();self.plan=plan
        self.group_id=next(i for i,g in enumerate(plan.groups) if world.rank in g.ranks)
        self.group=plan.groups[self.group_id]
        self.solver_comm=world.Split(self.group_id,self.group.ranks.index(world.rank))
        self.is_leader=self.solver_comm.rank==0
        self.leaders=world.Split(0 if self.is_leader else MPI.UNDEFINED,self.group_id)
        self.field_group=next(i for i,g in enumerate(plan.groups) if g.role=='field')
        self.closed=False;self._dependents=0
        return self

    def close(self):
        """Call after native detach, collectively in all groups; no MPI_Finalize."""
        if self.closed: return
        if threading.current_thread() is not threading.main_thread(): raise RuntimeError('MPI calls require main thread')
        busy=any(self.solver_comm.allgather(self._dependents))
        if self.is_leader: busy=any(self.leaders.allgather(busy))
        busy=self.solver_comm.bcast(busy,root=0)
        if busy:
            raise RuntimeError('close response sessions and detach native contexts first')
        if self.is_leader: self.leaders.Free()
        self.solver_comm.Free();self.closed=True

    def __enter__(self): return self
    def __exit__(self,*args): self.close()


class NativeParallelContext:
    """Attach the external SALMON communication layer to this solver group.

    A host-MPI bridge must be built explicitly. This attaches communication only;
    it does not make TDDFT globals reentrant or expose a full solver library.
    One attached native context per process/library is supported.
    """
    def __init__(self, context, library, sha256):
        import ctypes as c
        import hashlib
        from pathlib import Path
        from mpi4py import MPI
        self.context=context;self.closed=True
        error=None
        try:
            if context.closed: raise ValueError('closed Python context')
            if threading.current_thread() is not threading.main_thread(): raise RuntimeError('main thread required')
            path=Path(library).resolve()
            if hashlib.sha256(path.read_bytes()).hexdigest()!=sha256: raise ValueError('bridge checksum mismatch')
            lib=c.CDLL(str(path))
            lib.sp_attach.argtypes=[c.c_int64];lib.sp_attach.restype=c.c_int
            lib.sp_detach.argtypes=[];lib.sp_detach.restype=c.c_int
            lib.sp_info.argtypes=[c.POINTER(c.c_int),c.POINTER(c.c_int)];lib.sp_info.restype=None
            lib.sp_barrier.argtypes=[];lib.sp_barrier.restype=None
            lib.sp_sum.argtypes=[c.c_double];lib.sp_sum.restype=c.c_double
            lib.sp_library.argtypes=[c.POINTER(c.c_char),c.c_int,c.POINTER(c.c_int)];lib.sp_library.restype=None
            text=c.create_string_buffer(8192);active=c.c_int()
            lib.sp_library(text,len(text),c.byref(active))
            if not active.value or text.value.decode().strip()!=MPI.Get_library_version().rstrip('\x00').strip():
                raise ValueError('bridge and mpi4py must use the same initialized MPI runtime')
            self._lib=lib
        except Exception as exc: error=str(exc)
        errors=context.solver_comm.allgather(error)
        if any(errors): raise ValueError('native attachment rejected: '+str(errors))
        code=self._lib.sp_attach(context.solver_comm.py2f())
        if code: raise RuntimeError('native attachment failed: '+str(code))
        self.closed=False
        context._dependents+=1

    def _live(self):
        if self.closed or self.context.closed: raise RuntimeError('native context is closed')
        if threading.current_thread() is not threading.main_thread(): raise RuntimeError('main thread required')

    def info(self):
        import ctypes as c
        self._live();rank=c.c_int();size=c.c_int()
        self._lib.sp_info(c.byref(rank),c.byref(size));return rank.value,size.value

    def sum(self,value):
        """Collective sum within the native solver group (attachment diagnostic)."""
        self._live();return self._lib.sp_sum(float(value))

    def barrier(self):
        self._live();self._lib.sp_barrier()

    def close(self):
        if self.closed: return
        self._live()
        code=self._lib.sp_detach()
        if code: raise RuntimeError('native detach failed: '+str(code))
        self.closed=True
        self.context._dependents-=1

    def __enter__(self): return self
    def __exit__(self,*args): self.close()
