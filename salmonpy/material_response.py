"""Experimental distributed A/c -> physical charge-current response exchange.

This is a response service, not an E/H FDTD integrator. dt is fixed. A failed
backend step invalidates the session; no cross-group rollback is promised.
"""
from .errors import InputError, InputTypeError, RuntimeFailure
from dataclasses import dataclass
import numpy as np
import time
from collections.abc import Mapping
from .config import _copy_json, _freeze, _identifier


@dataclass(frozen=True)
class SolverContext:
    """Callback view: group communicator only, no orchestration collectives."""
    solver_comm: object
    group_name: str
    point_ids: tuple
    point_records: object = None


@dataclass(frozen=True)
class MaterialBinding:
    key: str
    provider: object
    backend: object
    material_record: object = None

    def __post_init__(self):
        if self.material_record is not None:
            from .point_runs import _identity
            if not isinstance(self.material_record, Mapping):
                raise InputError('material_record must be an identity mapping')
            object.__setattr__(self, 'material_record', _freeze(_identity(dict(self.material_record))))


def _point_metadata(plan, record, mapping):
    if record is None and mapping is None:
        return None
    from .point_runs import _record
    if record is None or not isinstance(mapping, Mapping):
        raise InputError('point_record and point_map must be supplied together')
    checked = _record(record)
    if any(type(p) is not int for p in mapping) or set(mapping) != set(plan.point_ids):
        raise InputError('point_map must cover exactly the integer transport points')
    for value in mapping.values():
        _identifier(value, 'point_map value')
    names = [entry['point_id'] for entry in checked['points']]
    if [mapping[p] for p in plan.point_ids] != names:
        raise InputError('point record order must match explicit point_map in transport order')
    owners = []
    for p in plan.point_ids:
        group = next(g for g in plan.groups if p in g.point_ids)
        owners.append({'transport_point_id': p, 'point_id': mapping[p],
                       'group_name': group.name, 'ranks': list(group.ranks)})
    return {'point_run': checked, 'mpi_point_map': owners}


def _response_timing(dt,integrator):
    """Response-service contract: endpoint A at n*dt in, charge current j at n*dt out."""
    from .timing import TimeSpec,InputPort,OutputPort
    return TimeSpec(dt,(InputPort('A'),),(OutputPort('j'),),integrator=integrator)


class _StepAdapter:
    """Retain accepted A for legacy two-endpoint callbacks."""
    def __init__(self, state):
        self.state = state
        self.A = np.zeros(3)

    def step(self, A):
        callback = getattr(self.state, 'step', None)
        result = callback(A.copy()) if callback is not None else self.state.advance(self.A.copy(), A.copy())
        self.A = A.copy()
        return result


class _Measurements:
    """Optional wall times; broken clocks never change the collective schedule."""
    def __init__(self):
        self.seconds = dict.fromkeys(('A_distribution', 'material', 'j_collection'), 0.)
        self.calls = dict.fromkeys(self.seconds, 0)
        self.errors = []

    def start(self):
        try:
            return time.perf_counter()
        except Exception as exc:
            self.errors.append(str(exc))
            return None

    def stop(self, name, start):
        self.calls[name] += 1
        end = self.start()
        if start is None or end is None or not np.isfinite(end-start) or end < start:
            self.seconds[name] = None
        elif self.seconds[name] is not None:
            self.seconds[name] += end-start


class ResponseSession:
    """All context ranks participate, including the field group.

    Providers implement identity, prepare(context). Backends implement identity,
    execution ('collective' or 'leader_only'), reduction ('sum' or 'leader'),
    create_state(initial, point_id, dt, context). States implement
    step(A_new) -> 3-vector physical charge-current density and close().
    Optional field_contracts selects explicit step_fields(FieldPacket) -> j.
    Legacy advance(A_old,A_new) callbacks are adapted with owned accepted A.
    Every result is tagged implicitly by the fixed plan's point IDs and accepted
    step index. Initial A is zero; current is at the new endpoint. All units a.u.
    """
    def __init__(self,context,bindings,*,dt,point_record=None,point_map=None,field_contracts=None):
        self.context=context;self.dt=dt;self.step_index=0;self.failed=False;self.closed=False;self._registered=False
        self._measurements = _Measurements()
        self._adapters = {}
        self._point_metadata = None
        self._field_contracts = None
        self.states={};self.bindings=dict(bindings);self._A=np.zeros((len(context.plan.point_ids),3))
        error=None;descriptors=None
        try:
            if not np.isfinite(dt) or dt<=0: raise InputError('positive finite dt required')
            if set(bindings)!=set(context.plan.point_ids): raise InputError('binding required for every point')
            self._point_metadata = _point_metadata(context.plan, point_record, point_map)
            contract_record = None
            if field_contracts is not None:
                from .field_input import FieldContract
                if not isinstance(field_contracts, Mapping) or set(field_contracts)!=set(context.plan.point_ids):
                    raise InputError('field contracts required for every point')
                self._field_contracts = dict(field_contracts)
                contract_record = []
                for point in context.plan.point_ids:
                    contract = self._field_contracts[point]
                    if not isinstance(contract,FieldContract) or contract.timing.dt != dt:
                        raise InputError('FieldContract with session dt required')
                    outputs=contract.timing.outputs
                    if (len(outputs)!=1 or outputs[0].quantity!='j' or outputs[0].offset not in (0,.5)
                            or outputs[0].sampling!='point' or contract.timing.clock_offset!=0
                            or contract.timing.t0!=0):
                        raise InputError('packet response requires point charge j at n*dt or (n+1/2)*dt from zero')
                    contract_record.append((point,contract.timing.to_dict(),dict(contract.units),
                                            contract.position.tolist(),contract.position_unit,contract.gauge))
            descriptors=[]
            keys={}
            identities={}
            for p in context.plan.point_ids:
                b=bindings[p]
                if b.key in keys and keys[b.key] is not b:
                    raise InputError('reuse the same MaterialBinding instance for one material key')
                keys[b.key]=b
                if self._point_metadata is not None:
                    entry = self._point_metadata['point_run']['points'][len(descriptors)]
                    identity = {k: v for k, v in entry.items() if k != 'point_id'}
                    if b.material_record is not None and dict(b.material_record) != identity:
                        raise InputError('point identity differs from recorded binding')
                    if identities.setdefault(b.key, identity) != identity:
                        raise InputError('one cached binding key cannot identify different initial/material data')
                mode=b.backend.execution;reduction=b.backend.reduction
                if not b.key or not b.provider.identity or not b.backend.identity: raise InputError('stable binding identities required')
                if mode not in ('collective','leader_only') or reduction not in ('sum','leader'): raise InputError('unsupported execution/reduction')
                owner=next(g for g in context.plan.groups if p in g.point_ids)
                if mode=='leader_only' and len(owner.ranks)!=1: raise InputError('leader-only backends require one rank per group')
                descriptors.append((p,b.key,b.provider.identity,b.backend.identity,mode,reduction,float(dt)))
        except Exception as exc: error=f'{type(exc).__name__}: {exc}'
        self._agree(error)
        # Metadata agreement catches mismatched backend selection before callbacks.
        descriptors = (descriptors, self._point_metadata, contract_record)
        local=context.solver_comm.allgather(descriptors)
        mismatch=any(v!=descriptors for v in local)
        if context.is_leader:
            all_descriptors=context.leaders.allgather(descriptors)
            mismatch=mismatch or any(v!=descriptors for v in all_descriptors)
        self._agree('binding descriptors differ across ranks' if mismatch else None)
        cache={}
        callback_records = None
        if self._point_metadata is not None:
            record = self._point_metadata['point_run']
            callback_records = _freeze({p: {'run_id': record['run_id'], **entry}
                                        for p, entry in zip(context.plan.point_ids, record['points'])
                                        if p in context.group.point_ids})
        callback_context=SolverContext(context.solver_comm,context.group.name,context.group.point_ids,callback_records)
        for point in context.group.point_ids:
            b=bindings[point];error=None
            try:
                if b.key not in cache: cache[b.key]=b.provider.prepare(callback_context)
            except Exception as exc: error=f'{type(exc).__name__}: {exc}'
            # Group-only agreement: groups can own different numbers of points.
            errors=context.solver_comm.allgather(error)
            if any(errors): error=str(errors);break
            try:
                state=b.backend.create_state(cache[b.key],point,dt,callback_context)
                if any(state is old for old in self.states.values()): raise InputError('backend reused mutable state across points')
            except Exception as exc: error=f'{type(exc).__name__}: {exc}'
            errors=context.solver_comm.allgather(error)
            if any(errors): error=str(errors);break
            self.states[point]=state
            self._adapters[point] = _StepAdapter(state)
            if self._field_contracts is not None:
                try:
                    from .field_input import FieldAdapter
                    self._adapters[point]=FieldAdapter(state,self._field_contracts[point])
                    declared=getattr(state,'timing',None)
                    if declared is not None and declared.outputs!=self._field_contracts[point].timing.outputs:
                        raise InputError('backend output timing differs from charge-current contract')
                except Exception as exc: error=f'{type(exc).__name__}: {exc}'
            errors=context.solver_comm.allgather(error)
            if any(errors): error=str(errors);break
        try: self._agree(error)
        except Exception:
            self.failed=True
            # Only collectively completed states are owned by this session.
            # A failed create callback must release its own partial resources.
            self.close()
            raise
        context._dependents+=1;self._registered=True
        try:
            self.initial_current = self._initial_current_records()
        except BaseException:
            self.failed = True
            self.close()
            raise

    def agree(self,error):
        """Collective validation agreement on all context ranks."""
        c=self.context
        errors=c.solver_comm.allgather(error)
        group_error=next((e for e in errors if e),None)
        if c.is_leader:
            errors=c.leaders.allgather(group_error)
            failure=next((e for e in errors if e),None)
        else: failure=None
        failure=c.solver_comm.bcast(failure,root=0)
        if failure: raise InputError(failure)

    _agree = agree  # Compatibility for existing callers.

    def identity_record(self):
        """Detached agreed point/ownership metadata, or None on the legacy path."""
        return _copy_json(self._point_metadata)

    def save_arrays(self, directory=None, arrays=None, *, point_axes=None, metadata=None):
        """Collective output; only field leader provides arguments and receives archive.

        No additional stepping or units conversion. Errors propagate to every rank.
        Output must follow plan.point_ids order. Does not create a restart.
        """
        from .point_runs import _axes
        from .outputs import save_arrays
        c = self.context
        writer = c.group.role == 'field' and c.is_leader
        error = None
        checked = info = None
        try:
            if self.closed or self.failed or self._point_metadata is None:
                raise InputError('saving point arrays requires an active identified session')
            if writer:
                if directory is None or not isinstance(arrays, Mapping):
                    raise InputError('field leader must supply directory and arrays')
                checked = {name: np.asarray(value) for name, value in arrays.items()}
                axes = _axes(point_axes, {name: value.shape for name, value in checked.items()}, len(c.plan.point_ids))
                info = _copy_json({} if metadata is None else metadata)
                if type(info) is not dict or set(info) & {'point_run', 'point_axes', 'mpi_point_map'}:
                    raise InputError('reserved point identity metadata keys')
                info.update(self.identity_record(), point_axes=axes)
            elif any(v is not None for v in (directory, arrays, point_axes, metadata)):
                raise InputError('only field leader supplies output arguments')
        except Exception as exc:
            error = f'{type(exc).__name__}: {exc}'
        self.agree(error)
        archive = None
        error = None
        if writer:
            try:
                archive = save_arrays(directory, checked, metadata=info)
            except Exception as exc:
                error = f'{type(exc).__name__}: {exc}'
        self.agree(error)
        return archive

    def timing_records(self):
        """Collective point timing declarations, returned on every rank."""
        from .timing import matter_timing
        records = {}; error = None
        try:
            records = {p: (self._field_contracts[p].timing if self._field_contracts is not None
                           else matter_timing(s, self.dt)).to_dict() for p, s in self.states.items()}
        except Exception as exc:
            error = f'{type(exc).__name__}: {exc}'
        self.agree(error)
        c = self.context
        groups = c.leaders.allgather(records) if c.is_leader else None
        groups = c.solver_comm.bcast(groups, root=0)
        return {p: spec for group in groups for p, spec in group.items()}

    def _initial_current_records(self):
        """Read physical t=0 current without stepping; None means unavailable."""
        from mpi4py import MPI
        c = self.context
        records = {}; error = None
        for p, state in self.states.items():
            value = None
            try:
                value = getattr(state, 'initial_current', None)
                if value is None:
                    value = getattr(state, 'current', None)
                if value is not None:
                    value = np.array(value, dtype=float, copy=True)
                    if value.shape != (3,) or not np.isfinite(value).all():
                        raise InputError('initial current must be a finite 3-vector')
            except Exception as exc:
                error = f'{type(exc).__name__}: {exc}'
            statuses = c.solver_comm.allgather((error, value is not None))
            if any(e for e, _ in statuses):
                error = next(e for e, _ in statuses if e)
                break
            if len({present for _, present in statuses}) != 1:
                error = 'initial current availability differs within material group'
                break
            if value is not None and self.bindings[p].backend.reduction == 'sum':
                value = c.solver_comm.allreduce(value, op=MPI.SUM)
            records[p] = None if value is None else value.tolist()
        self.agree(error)
        groups = c.leaders.allgather(records) if c.is_leader else None
        groups = c.solver_comm.bcast(groups, root=0)
        return {p: value for group in groups for p, value in group.items()}

    def communication_report(self):
        """Collective per-rank wall times (including collective wait, no barriers).

        Payload bytes are logical arrays per step, not network traffic. Material
        time includes backend work, its subgroup reductions and error agreement.
        Clock failure yields null times and does not fail a calculation.
        """
        c = self.context
        m = self._measurements
        local = {'group': c.group.name, 'solver_rank': c.solver_comm.rank,
                 'seconds': dict(m.seconds), 'calls': dict(m.calls),
                 'measurement_errors': list(m.errors)}
        rows = c.solver_comm.gather(local, root=0)
        groups = c.leaders.allgather(rows) if c.is_leader else None
        groups = c.solver_comm.bcast(groups, root=0)
        return {'schema': 'salmonpy.communication/1', 'steps': self.step_index,
                'logical_A_bytes_per_step': len(c.plan.point_ids)*3*8 if self._field_contracts is None else None,
                'logical_field_value_bytes_per_step': (len(c.plan.point_ids)*3*8 if self._field_contracts is None
                    else sum(len(v.units) for v in self._field_contracts.values())*3*8),
                'field_metadata_bytes_counted': False,
                'logical_j_bytes_per_step': len(c.plan.point_ids)*3*8,
                'ranks': [row for group in groups for row in group]}

    def _step(self,A,*,field_packets=False):
        """Collective everywhere. Only field leader supplies (points,3) A.

        Returns ordered physical current on every field rank, None on material
        ranks. User callbacks may use solver_comm but not any outer communicator.
        """
        from mpi4py import MPI
        c=self.context;error=None
        try:
            if self.closed or self.failed or c.closed: raise InputError('session is not usable')
            if field_packets != (self._field_contracts is not None):
                raise InputError('use step_fields with declared field contracts, otherwise step(A)')
            if field_packets:
                values=self._pack_fields(A)
            elif c.group.role=='field' and c.is_leader:
                values=np.array(A,dtype=np.float64,order='C',copy=True)
                if values.shape!=self._A.shape or not np.isfinite(values).all(): raise InputError('finite A array of shape (points,3) required')
            else:
                if A is not None: raise InputError('only field leader supplies A')
                values=np.empty_like(self._A)
        except Exception as exc: error=f'{type(exc).__name__}: {exc}'
        self._agree(error)
        started = self._measurements.start()
        if field_packets:
            if c.is_leader: values=c.leaders.bcast(values,root=c.field_group)
            values=c.solver_comm.bcast(values,root=0)
            from .field_input import FieldPacket, FieldValue
            error=None
            try:
                values={p:FieldPacket(**{q:FieldValue(**v) for q,v in packet.items()})
                        for p,packet in values.items()}
                for p,contract in self._field_contracts.items():
                    contract.validate(values[p],step_index=self.step_index+1)
            except Exception as exc: error=f'{type(exc).__name__}: {exc}'
            self._agree(error)
        else:
            if c.is_leader: c.leaders.Bcast(values,root=c.field_group)
            c.solver_comm.Bcast(values,root=0)
        self._measurements.stop('A_distribution', started)
        started = self._measurements.start()
        currents=np.zeros((len(c.group.point_ids),3));error=None
        for i,point in enumerate(c.group.point_ids):
            k=c.plan.point_ids.index(point);b=self.bindings[point]
            try:
                raw = self._adapters[point].step_fields(values[point]) if field_packets else self._adapters[point].step(values[k].copy())
                j=np.array(raw,dtype=np.float64,order="C",copy=True)
                if j.shape!=(3,) or not np.isfinite(j).all(): raise InputError('backend must return finite physical current density (3,)')
            except Exception as exc: error=f'{type(exc).__name__}: {exc}'
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
        self._measurements.stop('material', started)
        started = self._measurements.start()
        result=None
        if c.is_leader:
            counts=[len(g.point_ids)*3 for g in c.plan.groups]
            packed=np.empty(sum(counts),dtype=np.float64) if c.group.role=='field' else None
            c.leaders.Gatherv(currents.ravel(),[packed,counts,MPI.DOUBLE] if packed is not None else None,root=c.field_group)
            if packed is not None:
                order=[p for g in c.plan.groups for p in g.point_ids]
                result=np.array([packed.reshape(-1,3)[order.index(p)] for p in c.plan.point_ids])
        if c.group.role=='field': result=c.solver_comm.bcast(result,root=0)
        self._measurements.stop('j_collection', started)
        if not field_packets: self._A=values
        self.step_index+=1
        return result

    def _pack_fields(self,fields):
        from .field_input import FieldPacket
        c=self.context
        if c.group.role!='field' or not c.is_leader:
            if fields is not None: raise InputError('only field leader supplies fields')
            return None
        if not isinstance(fields,Mapping) or set(fields)!=set(c.plan.point_ids):
            raise InputError('FieldPacket mapping must cover exactly transport point IDs')
        packed={}
        for point,contract in self._field_contracts.items():
            packet=fields[point]
            contract.validate(packet,step_index=self.step_index+1)
            packed[point]={}
            for q in contract.units:
                v=getattr(packet,q)
                packed[point][q]=dict(value=v.value.tolist(),time=v.time,position=v.position.tolist(),
                    unit=v.unit,position_unit=v.position_unit,sampling=v.sampling)
        return packed

    def step(self,A):
        """Collective legacy A/c endpoint input, supplied only by field leader."""
        return self._step(A)

    def step_fields(self,fields):
        """Collective explicit per-point A/E/H packets; returns charge current j.

        Contracts are required at construction. No interpolation, conversion or
        dipole-to-current normalization is performed. Only field leader supplies
        a mapping keyed by transport point IDs; all other ranks supply None.
        """
        return self._step(fields,field_packets=True)

    advance = step  # Existing prescribed-A entry point.

    def close(self):
        if self.closed: return
        error=None
        for point in self.context.group.point_ids:
            if point in self.states:
                try: self.states[point].close()
                except Exception as exc: error=f'{type(exc).__name__}: {exc}'
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
        if not isinstance(model,LorentzDrude): raise InputTypeError('LorentzDrude model required')
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
            def close(self): state.close()
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
        if self._input.is_absolute() or '..' in self._input.parts: raise InputError('relative bundle input required')
        if not (self._directory/self._input).is_file(): raise InputError('missing TDDFT input file')
        self._max_steps=_integer(max_steps,'max_steps',1,2**31-2)
        self._threads=_integer(threads,'threads',1,2**31-1)
        self._timeout=_positive(timeout,'timeout')
        manifest=[]
        for p in sorted(self._directory.rglob('*')):
            if p.is_symlink(): raise InputError('bundle symlinks are unsupported')
            if p.is_file(): manifest.append((str(p.relative_to(self._directory)),_digest(p)))
            elif not p.is_dir(): raise InputError('nonregular bundle entry')
        self._manifest=tuple(manifest)
        record={'files':self._manifest,'input':str(self._input),'binary':self._executable.sha256,
                'max_steps':self._max_steps,'threads':self._threads,'timeout':self._timeout}
        self._identity='salmonpy.tddft.process/'+hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()
        self._initial=None;self._closed=False

    @property
    def identity(self): return self._identity

    def prepare(self,context):
        from .tddft import TDDFTInitialState
        if self._closed: raise RuntimeFailure('TDDFT material provider is closed')
        if self._initial is None:
            initial=TDDFTInitialState.from_salmon(self._directory,self._input,executable=self._executable,
                      max_steps=self._max_steps,threads=self._threads,timeout=self._timeout)
            if initial.manifest!=self._manifest:
                initial.close()
                raise InputError('initial bundle changed after material identity was established')
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
            raise InputTypeError('TDDFTPreparedState with managed serial runtime required')
        verified=load_runtime(initial.runtime,required=['tddft','single-rank'])
        if verified.sha256!=initial.initial.executable.sha256:
            raise InputError('initial-state executable differs from verified serial runtime')
        from dataclasses import replace
        worker=replace(initial.initial,executable=verified).new_state(dt=dt)
        class State:
            timing=_response_timing(dt,'SALMON middlepoint TDDFT worker; charge current j')
            def __init__(self): self.A=np.zeros(3)
            def advance(self,old,new):
                if not np.array_equal(old,self.A): raise InputError('old A differs from accepted TDDFT state')
                current=worker.step(new)  # TDDFTState returns charge current j
                self.A=np.array(new,copy=True)
                return current
            def close(self): worker.close()
        return State()
