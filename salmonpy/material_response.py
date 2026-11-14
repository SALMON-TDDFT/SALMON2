"""Experimental distributed A/c -> physical charge-current response exchange.

This is a response service, not an E/H FDTD integrator. dt is fixed. A failed
backend step invalidates the session; no cross-group rollback is promised.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SolverContext:
    """Callback view: group communicator only, no orchestration collectives."""
    solver_comm: object
    group_name: str
    point_ids: tuple


@dataclass(frozen=True)
class MaterialBinding:
    key: str
    provider: object
    backend: object


def _response_timing(dt,integrator):
    """Response-service contract: endpoint A at n*dt in, charge current j at n*dt out."""
    from .timing import TimeSpec,InputPort,OutputPort
    return TimeSpec(dt,(InputPort('A'),),(OutputPort('j'),),integrator=integrator)


class ResponseSession:
    """All context ranks participate, including the field group.

    Providers implement identity, prepare(context). Backends implement identity,
    execution ('collective' or 'leader_only'), reduction ('sum' or 'leader'),
    create_state(initial, point_id, dt, context). States implement
    advance(A_old,A_new) -> 3-vector physical charge-current density and close().
    Every result is tagged implicitly by the fixed plan's point IDs and accepted
    step index. Initial A is zero; current is at the new endpoint. All units a.u.
    """
    def __init__(self,context,bindings,*,dt):
        self.context=context;self.dt=dt;self.step_index=0;self.failed=False;self.closed=False;self._registered=False
        self.states={};self.bindings=dict(bindings);self._A=np.zeros((len(context.plan.point_ids),3))
        error=None;descriptors=None
        try:
            if not np.isfinite(dt) or dt<=0: raise ValueError('positive finite dt required')
            if set(bindings)!=set(context.plan.point_ids): raise ValueError('binding required for every point')
            descriptors=[]
            keys={}
            for p in context.plan.point_ids:
                b=bindings[p]
                if b.key in keys and keys[b.key] is not b:
                    raise ValueError('reuse the same MaterialBinding instance for one material key')
                keys[b.key]=b
                mode=b.backend.execution;reduction=b.backend.reduction
                if not b.key or not b.provider.identity or not b.backend.identity: raise ValueError('stable binding identities required')
                if mode not in ('collective','leader_only') or reduction not in ('sum','leader'): raise ValueError('unsupported execution/reduction')
                owner=next(g for g in context.plan.groups if p in g.point_ids)
                if mode=='leader_only' and len(owner.ranks)!=1: raise ValueError('leader-only backends require one rank per group')
                descriptors.append((p,b.key,b.provider.identity,b.backend.identity,mode,reduction,float(dt)))
        except Exception as exc: error=str(exc)
        self._agree(error)
        # Metadata agreement catches mismatched backend selection before callbacks.
        local=context.solver_comm.allgather(descriptors)
        mismatch=any(v!=descriptors for v in local)
        if context.is_leader:
            all_descriptors=context.leaders.allgather(descriptors)
            mismatch=mismatch or any(v!=descriptors for v in all_descriptors)
        self._agree('binding descriptors differ across ranks' if mismatch else None)
        cache={}
        callback_context=SolverContext(context.solver_comm,context.group.name,context.group.point_ids)
        for point in context.group.point_ids:
            b=bindings[point];error=None
            try:
                if b.key not in cache: cache[b.key]=b.provider.prepare(callback_context)
            except Exception as exc: error=str(exc)
            # Group-only agreement: groups can own different numbers of points.
            errors=context.solver_comm.allgather(error)
            if any(errors): error=str(errors);break
            try:
                state=b.backend.create_state(cache[b.key],point,dt,callback_context)
                if any(state is old for old in self.states.values()): raise ValueError('backend reused mutable state across points')
            except Exception as exc: error=str(exc)
            errors=context.solver_comm.allgather(error)
            if any(errors): error=str(errors);break
            self.states[point]=state
        try: self._agree(error)
        except Exception:
            self.failed=True
            # Only collectively completed states are owned by this session.
            # A failed create callback must release its own partial resources.
            self.close()
            raise
        context._dependents+=1;self._registered=True

    def _agree(self,error):
        c=self.context
        errors=c.solver_comm.allgather(error)
        group_error=next((e for e in errors if e),None)
        if c.is_leader:
            errors=c.leaders.allgather(group_error)
            failure=next((e for e in errors if e),None)
        else: failure=None
        failure=c.solver_comm.bcast(failure,root=0)
        if failure: raise ValueError(failure)

    def advance(self,A):
        """Collective everywhere. Only field leader supplies (points,3) A.

        Returns ordered physical current on every field rank, None on material
        ranks. User callbacks may use solver_comm but not any outer communicator.
        """
        from mpi4py import MPI
        c=self.context;error=None
        try:
            if self.closed or self.failed or c.closed: raise ValueError('session is not usable')
            if c.group.role=='field' and c.is_leader:
                values=np.array(A,dtype=np.float64,order='C',copy=True)
                if values.shape!=self._A.shape or not np.isfinite(values).all(): raise ValueError('finite A array of shape (points,3) required')
            else:
                if A is not None: raise ValueError('only field leader supplies A')
                values=np.empty_like(self._A)
        except Exception as exc: error=str(exc)
        self._agree(error)
        if c.is_leader: c.leaders.Bcast(values,root=c.field_group)
        c.solver_comm.Bcast(values,root=0)
        currents=np.zeros((len(c.group.point_ids),3));error=None
        for i,point in enumerate(c.group.point_ids):
            k=c.plan.point_ids.index(point);b=self.bindings[point]
            try:
                j=np.array(self.states[point].advance(self._A[k].copy(),values[k].copy()),dtype=np.float64,order="C",copy=True)
                if j.shape!=(3,) or not np.isfinite(j).all(): raise ValueError('backend must return finite physical current density (3,)')
            except Exception as exc: error=str(exc)
            errors=c.solver_comm.allgather(error)
            if any(errors): error=str(errors);break
            if b.backend.reduction=='sum': c.solver_comm.Allreduce(j,currents[i],op=MPI.SUM)
            elif c.is_leader: currents[i]=j
            if not np.isfinite(currents[i]).all(): error='nonfinite reduced current'
            errors=c.solver_comm.allgather(error)
            if any(errors): error=str(errors);break
        try: self._agree(error)
        except ValueError:
            self.failed=True
            raise
        result=None
        if c.is_leader:
            counts=[len(g.point_ids)*3 for g in c.plan.groups]
            packed=np.empty(sum(counts),dtype=np.float64) if c.group.role=='field' else None
            c.leaders.Gatherv(currents.ravel(),[packed,counts,MPI.DOUBLE] if packed is not None else None,root=c.field_group)
            if packed is not None:
                order=[p for g in c.plan.groups for p in g.point_ids]
                result=np.array([packed.reshape(-1,3)[order.index(p)] for p in c.plan.point_ids])
        if c.group.role=='field': result=c.solver_comm.bcast(result,root=0)
        self._A=values;self.step_index+=1
        return result

    def close(self):
        if self.closed: return
        error=None
        for point in self.context.group.point_ids:
            if point in self.states:
                try: self.states[point].close()
                except Exception as exc: error=str(exc)
        self.closed=True
        if self._registered:
            self.context._dependents-=1;self._registered=False
        self._agree(error)


class LDMaterialProvider:
    """Immutable LD parameters; no material state is shared across macro points."""
    def __init__(self,model):
        import hashlib,json
        from dataclasses import asdict
        from .ld import LorentzDrude
        if not isinstance(model,LorentzDrude): raise TypeError('LorentzDrude model required')
        self._config=model.config
        self.identity='salmonpy.ld/'+hashlib.sha256(json.dumps(asdict(self._config),sort_keys=True).encode()).hexdigest()

    def prepare(self,context):
        from .ld import LorentzDrude
        c=self._config
        return LorentzDrude(plasma_frequency=c.plasma_frequency,strength=c.strength,
                            resonance=c.resonance,damping=c.damping)


class LDResponseBackend:
    """Standalone endpoint-A driver; NOT legacy semi-implicit E/H coupling.

    E_new=-(A_new-A_old)/dt drives the existing LDState.step. J is its physical
    current density. Initial E/P/J are zero. No FDTD Q or coefficient injection.
    """
    identity='salmonpy.ld.a-endpoints/1'
    execution='leader_only'
    reduction='leader'

    def create_state(self,initial,point_id,dt,context):
        state=initial.new_state(dt=dt)
        class State:
            timing=_response_timing(dt,'LD driven by E_new = -(A_new-A_old)/dt; physical polarization current')
            def advance(self,old,new): return state.step(-(new-old)/dt)[1]
            def close(self): pass
        return State()


@dataclass(frozen=True)
class TDDFTPreparedState:
    """Provider-owned initial snapshot paired with its managed serial runtime.

    Custom providers can return this descriptor. The response backend independently
    verifies the runtime and matches the snapshot's executable checksum before use.
    """
    initial: object
    runtime: object


class TDDFTMaterialProvider:
    """Lazy immutable GS snapshot for process-isolated, serial native workers.

    A verified managed single-rank runtime is required. Material preparation does
    not invoke DFT: supply an existing GS bundle. Close the provider after sessions.
    Each macro point's worker copies the snapshot and owns its evolution history.
    """
    def __init__(self,directory,input_file,*,runtime,max_steps=1000,threads=1,timeout=300.):
        from pathlib import Path
        import hashlib,json
        from .runtime import load_runtime
        from .sbe import _integer
        from .ld import _positive
        from .tddft import _digest
        self._runtime=Path(runtime).expanduser().resolve()
        self._executable=load_runtime(self._runtime,required=['tddft','single-rank'])
        self._directory=Path(directory).expanduser().resolve()
        self._input=Path(input_file)
        if self._input.is_absolute() or '..' in self._input.parts: raise ValueError('relative bundle input required')
        if not (self._directory/self._input).is_file(): raise ValueError('missing TDDFT input file')
        self._max_steps=_integer(max_steps,'max_steps',1,2**31-2)
        self._threads=_integer(threads,'threads',1,2**31-1)
        self._timeout=_positive(timeout,'timeout')
        manifest=[]
        for p in sorted(self._directory.rglob('*')):
            if p.is_symlink(): raise ValueError('bundle symlinks are unsupported')
            if p.is_file(): manifest.append((str(p.relative_to(self._directory)),_digest(p)))
            elif not p.is_dir(): raise ValueError('nonregular bundle entry')
        self._manifest=tuple(manifest)
        record={'files':self._manifest,'input':str(self._input),'binary':self._executable.sha256,
                'max_steps':self._max_steps,'threads':self._threads,'timeout':self._timeout}
        self._identity='salmonpy.tddft.process/'+hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()
        self._initial=None;self._closed=False

    @property
    def identity(self): return self._identity

    def prepare(self,context):
        from .tddft import TDDFTInitialState
        if self._closed: raise RuntimeError('TDDFT material provider is closed')
        if self._initial is None:
            initial=TDDFTInitialState.from_salmon(self._directory,self._input,executable=self._executable,
                      max_steps=self._max_steps,threads=self._threads,timeout=self._timeout)
            if initial.manifest!=self._manifest:
                initial.close()
                raise ValueError('initial bundle changed after material identity was established')
            self._initial=initial
        return TDDFTPreparedState(self._initial,self._runtime)

    def close(self):
        if self._initial is not None: self._initial.close()
        self._closed=True

    def __enter__(self): return self
    def __exit__(self,*args): self.close()


class TDDFTResponseBackend:
    """Existing one-process-per-point TDDFT, not embedded subgroup MPI TDDFT.

    Returns the physical charge-current density that TDDFTState already reports
    (j = -SALMON electronic J); no further sign change here.
    Legacy middlepoint order/worker recovery remain owned by TDDFTState.
    """
    identity='salmonpy.tddft.middlepoint.process.charge/1'
    execution='leader_only'
    reduction='leader'

    def create_state(self,initial,point_id,dt,context):
        from .tddft import TDDFTInitialState
        from .runtime import load_runtime
        if not isinstance(initial,TDDFTPreparedState) or not isinstance(initial.initial,TDDFTInitialState):
            raise TypeError('TDDFTPreparedState with managed serial runtime required')
        verified=load_runtime(initial.runtime,required=['tddft','single-rank'])
        if verified.sha256!=initial.initial.executable.sha256:
            raise ValueError('initial-state executable differs from verified serial runtime')
        from dataclasses import replace
        worker=replace(initial.initial,executable=verified).new_state(dt=dt)
        class State:
            timing=_response_timing(dt,'SALMON middlepoint TDDFT worker; charge current j')
            def __init__(self): self.A=np.zeros(3)
            def advance(self,old,new):
                if not np.array_equal(old,self.A): raise ValueError('old A differs from accepted TDDFT state')
                current=worker.step(new)  # TDDFTState returns charge current j
                self.A=np.array(new,copy=True)
                return current
            def close(self): worker.close()
        return State()
