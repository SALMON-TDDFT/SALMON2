"""Experimental native E/H-to-A/J coupling on separate field/material MPI groups.

E is at native time n*dt. Electronic local time n*dt corresponds to field clock
(n+1/2)*dt; native saved H is centered at E time, not the raw H endpoint.
"""
from .errors import InputError, InputTypeError, NumericalError
from pathlib import Path
import ctypes as c
import json
import os
import threading
from fractions import Fraction
import numpy as np
from . import embedded_tddft as embedded
from .distributed import NativeParallelContext
from .fdtd import FDTD,FDTDMPI,PlaneWave
from .tddft import _digest
from .timing import TimeSpec,InputPort,OutputPort,check_loop

#: Electronic local time n*dt is field clock (n+1/2)*dt in the native E/H scheme.
ELECTRONIC_CLOCK_OFFSET='1/2'


def _agree(comm,error):
    errors=comm.allgather(error)
    if any(errors):raise InputError(str(errors))


def _require_field_record(record):
    """Admit the legacy field record or an explicitly declared fast field ABI."""
    if record.get('schema')=='salmonpy.embedded-tddft/1':return
    field=record.get('field_abi')
    if (record.get('schema')=='salmonpy.fast-hse-mpi/1' and isinstance(field,dict)
            and type(field.get('protocol')) is int and field['protocol']==1
            and field.get('scheme')=='salmonpy.eh-centered-current/1'):
        return
    raise InputError('embedded host library record or explicit fast-runtime E/H ABI required')


class NativeEHField:
    """Collective field-group context; points map IDs to zero-based grid indices.

    Each macro point represents one uniform field-cell volume. Coupled sites and
    their constitutive neighbors must be vacuum in the background geometry; LD
    may occupy other regions. Only zero-initialized plane-wave runs are supported.
    """
    scheme='salmonpy.eh-centered-current/1'
    electronic_clock_offset = ELECTRONIC_CLOCK_OFFSET

    def __init__(self,context,calculation,build_record,directory,*,points,grid):
        comm=context.solver_comm;error=None
        try:
            if context.group.role!='field':raise InputError('field group required')
            if embedded._used:raise InputError('one native lifetime per process')
            if threading.current_thread() is not threading.main_thread():raise InputError('main thread required')
            if not isinstance(calculation,FDTD) or not isinstance(calculation.source,PlaneWave):
                raise InputError('zero-initialized plane-wave FDTD required')
            parallel=FDTDMPI(tuple(grid))
            if parallel.ranks!=comm.size:raise InputError('field decomposition must match group size')
            text=parallel.input_text(calculation)
            self.point_ids=context.plan.point_ids
            if set(points)!=set(self.point_ids):raise InputError('one field index per material point required')
            indices=np.asarray([points[p] for p in self.point_ids])
            if indices.shape!=(len(self.point_ids),3) or indices.dtype.kind not in 'iu':raise InputError('integer grid triples required')
            counts=np.array(calculation.geometry.grid.counts)
            if np.any(indices<10) or np.any(indices>=counts-10):raise InputError('coupled cells must be at least ten cells from boundaries')
            if len({tuple(p) for p in indices})!=len(indices):raise InputError('one macro state per field cell required')
            axis=calculation.source.axis
            source_index=int(np.argmin(np.abs(calculation.geometry.grid.coordinates[axis]-calculation.source.position)))
            if np.any(np.abs(indices[:,axis]-source_index)<=2):
                raise InputError('coupled cells must be clear of the incident source plane and its support')
            ids=calculation.geometry.material_ids()
            for point in indices:
                # Conservative guard covers all constitutive cells touching mapped edges.
                i,j,k=point
                if np.any(ids[i-1:i+2,j-1:j+2,k-1:k+2]):raise InputError('coupled cells require a vacuum constitutive neighborhood')
            record=json.loads(Path(build_record).read_text())
            _require_field_record(record)
            library=Path(record['library']);sha=record['sha256']
            if _digest(library)!=sha:raise InputError('native library checksum mismatch')
            lib=c.CDLL(str(library))
            ints=np.ctypeslib.ndpointer(dtype=np.int32,ndim=2,flags='C_CONTIGUOUS')
            floats=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
            lib.sp_eh_init.argtypes=[c.c_char_p,c.c_int,ints];lib.sp_eh_init.restype=c.c_int
            lib.sp_eh_step.argtypes=[c.c_int,floats,floats];lib.sp_eh_step.restype=c.c_int
            lib.sp_eh_close.argtypes=[];lib.sp_eh_close.restype=c.c_int
            self.directory=Path(directory).expanduser().resolve()
            descriptor=(text,indices.tolist(),sha,str(self.directory))
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        descriptors=comm.allgather(descriptor)
        _agree(comm,'field configuration differs across ranks' if any(x!=descriptor for x in descriptors) else None)
        error=None
        if comm.rank==0:
            try:
                self.directory.mkdir(parents=True,exist_ok=False)
                (self.directory/'inputfile').write_text(text)
            except Exception as exc:error=str(exc)
        _agree(comm,error)
        error=None
        try:
            if (self.directory/'inputfile').read_text()!=text:raise InputError('shared field directory required')
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        self.checkpoint_decomposition=tuple(grid)
        self.context=context;self.lib=lib;self.dt=calculation.dt;self.capacity=calculation.steps
        self.step_index=0;self.A=np.zeros((len(self.point_ids),3));self.E=np.zeros_like(self.A)
        self.closed=False;self.initialized=False;self.failed=False
        from .field_input import _real
        xyz=np.column_stack([calculation.geometry.grid.coordinates[d][indices[:,d]] for d in range(3)])
        self.sample_positions=_real(xyz,(len(self.point_ids),3),'macro sample positions')
        self.manifest={'schema':self.scheme,'status':'started','point_ids':self.point_ids,
            'grid_indices':indices.tolist(),'sample_positions_bohr':self.sample_positions.tolist(),'cell_volume':float(np.prod(calculation.geometry.grid.spacing)),
            'dt':self.dt,'steps':self.capacity,'field_ranks':comm.size,'field_grid':list(grid),
            'library_sha256':sha,'sample_map':'half of same-index and previous-axis E edges',
            'deposit_map':'transpose of sample_map, equal cell volumes',
            'native_E_time':'n*dt','native_saved_H_time':'n*dt',
            'electronic_local_time':'n*dt','electronic_field_clock':'(n+0.5)*dt',
            'initial_feedback_current':'zero; not substituted with initial GS current',
            'timing':self.timing.to_dict(),'electronic_clock_offset':self.electronic_clock_offset}
        self._write_manifest()
        self.native=NativeParallelContext(context,library,sha)
        embedded._used=True
        try:
            code=self._call(lib.sp_eh_init,b'inputfile',len(self.point_ids),np.array(indices,dtype=np.int32,order='C'))
            self.initialized=(code==0)
            reason = {1001: 'invalid unit_system; expected au, a.u. or A_eV_fs',
                      1002: 'explicit atomic coordinates conflict with an MD restart',
                      1003: 'atomic coordinate source is missing or specified more than once',
                      1004: 'reduced atomic coordinates require a periodic system',
                      1005: 'cannot open atomic coordinate file; check the input path',
                      1006: 'cannot read atomic coordinates; check row count and column format',
                      1010: 'unsupported propagator',
                      1011: 'source-support ACE is deferred',
                      1012: 'hse_sr_tolerance must be zero in the supported profile',
                      1013: "a yn_* option only accepts y or n (case-insensitive)"}.get(
                          code, 'native field initialization rejected')
            _agree(comm,f'{reason} (native status {code})' if code else None)
        except BaseException:
            self.failed=True
            self.close();raise

    @property
    def timing(self):
        """salmonpy.time/1: j at (n+1/2)dt (lag 1, t=0 unused) -> A at (n+1/2)dt from E at n*dt."""
        return TimeSpec(self.dt,(InputPort('j',offset='1/2',lag=1,uses_t0=False),),(OutputPort('A',offset='1/2'),),
            history='zero-initialized plane-wave field; initial feedback current zero',
            integrator='native Yee: E at n*dt, saved H centered at E time; A_new = A_old - dt*S(E)')

    def sample_fields(self):
        """Read detached A/E point samples after a step; no collectives or update.

        A is at (n+1/2)dt, E at n*dt. Read H via sample_magnetic_field collectively.
        Positions follow plan.point_ids and the mapped macro-cell centers.
        """
        if self.closed or self.failed or self.step_index<1:
            raise InputError('field samples require a usable field after its first step')
        from .field_sampling import SampledField
        return {'A':SampledField('A',[(self.step_index+.5)*self.dt],self.sample_positions,
                                 self.A[None,:,:],'A/c:au'),
                'E':SampledField('E',[self.step_index*self.dt],self.sample_positions,
                                 self.E[None,:,:],'E:au')}

    def sample_magnetic_field(self):
        """Collective field-group read of macro-centered raw H at (n+1/2)dt.

        Requires the optional sp_eh_sample_h ABI; old libraries reject explicitly.
        Every field rank must call. No time advance or saved-H averaging occurs.
        """
        comm=self.context.solver_comm;error=None;fn=None
        try:
            if self.closed or self.failed or self.step_index<1:
                raise InputError('H samples require a usable field after its first step')
            fn=getattr(self.lib,'sp_eh_sample_h',None)
            if fn is None:raise InputError('runtime lacks sp_eh_sample_h; explicitly build an H-enabled runtime')
            self.native._live()
            fn.argtypes=[c.c_int,np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')]
            fn.restype=c.c_int
            H=np.empty_like(self.A)
        except Exception as exc:error=f'{type(exc).__name__}: {exc}'
        _agree(comm,error)
        error=None
        try:
            code=self._call(fn,len(self.point_ids),H)
            if code or not np.isfinite(H).all():error='native H sampling failed'
        except Exception as exc:error=f'{type(exc).__name__}: {exc}'
        _agree(comm,error)
        from .field_sampling import SampledField
        return SampledField('H',[(self.step_index+.5)*self.dt],self.sample_positions,H[None,:,:],'H:au')

    def _write_manifest(self):
        error=None
        if self.context.is_leader:
            try:(self.directory/'salmonpy-coupled-field.json').write_text(json.dumps(self.manifest,indent=2)+'\n')
            except Exception as exc:error=str(exc)
        _agree(self.context.solver_comm,error)

    def _call(self,fn,*args):
        from ._native_group import _directory_call
        self.native._live()
        return _directory_call(self.lib,self.context.solver_comm,self.directory,fn,*args,minimum_protocol=3)

    def advance(self,J_previous):
        comm=self.context.solver_comm;error=None
        try:
            if self.closed or self.failed or self.step_index>=self.capacity:raise InputError('closed, failed or exhausted field')
            J=np.array(J_previous,dtype=np.float64,order='C')
            if J.shape!=self.A.shape or not np.isfinite(J).all():raise InputError('finite point currents required')
            if self.step_index==0 and np.any(J):raise InputError('initial feedback current must be zero')
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        values=comm.allgather(J.tolist())
        _agree(comm,'field currents differ across ranks' if any(v!=values[0] for v in values) else None)
        E=np.empty_like(self.A)
        try:
            checked=getattr(self.lib,'sp_eh_step_checked',None)
            if checked is not None:
                floats=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
                checked.argtypes=[c.c_int,floats,floats,floats]
                checked.restype=c.c_int
                code=self._call(checked,len(self.point_ids),J,self.A,E)
            else:
                code=self._call(self.lib.sp_eh_step,len(self.point_ids),J,E)
            outcomes=comm.allgather((int(code),bool(np.isfinite(E).all())))
            if embedded._admission_rejected(self.lib,outcomes,allowed=(20,21,22,23),query_name='sp_eh_recovery_protocol'):
                raise NumericalError('native field numerical state restored without replay; retry permitted: '+str(outcomes))
            _agree(comm,'native field advance failed' if code or not np.isfinite(E).all() else None)
            with np.errstate(over="ignore",invalid="ignore"):
                updated_A=self.A-self.dt*E
            _agree(comm,'nonfinite integrated A' if not np.isfinite(updated_A).all() else None)
        except NumericalError:
            raise
        except BaseException as exc:
            self.failed=True
            self.manifest['failure']=str(exc)
            raise
        self.E=E;self.A=updated_A;self.step_index+=1
        return self.A.copy()

    def close(self):
        if self.closed:return
        code=0
        if self.initialized:
            code=self._call(self.lib.sp_eh_close)
        _agree(self.context.solver_comm,'native field close failed' if code else None)
        self.native.close();self.closed=True
        self.manifest['completed_steps']=self.step_index
        self.manifest['status']='failed' if self.failed else ('completed' if self.step_index==self.capacity else 'stopped')
        self._write_manifest()

    def __enter__(self):return self
    def __exit__(self,*args):self.close()


class MultiscaleField:
    """Single-rank execution of an existing Multiscale1D state.

    Plan point order maps to material cells 1..N. The original field owns its
    arrays, boundary update and initial zero-feedback rule.
    """
    electronic_clock_offset = '0'

    def __init__(self, context, state):
        from .maxwell import _Multiscale1DState
        if context.group.role != 'field' or context.solver_comm.size != 1:
            raise InputError('Multiscale1D requires exactly one field rank')
        if not isinstance(state, _Multiscale1DState):
            raise InputTypeError('Multiscale1D.new_state result required')
        if state.step_index != 0 or state.points != len(context.plan.point_ids):
            raise InputError('fresh field with one material cell per plan point required')
        self.context = context
        self.state = state
        self.point_ids = context.plan.point_ids
        self.closed = self.failed = False
        for point, A in zip(self.point_ids, state.A):
            if np.any(A != 0):
                raise InputError(f'point {point} has nonzero A(0)={A.tolist()}; '
                                 'increase left vacuum cells or delay the incident pulse')

    @property
    def dt(self): return self.state.dt
    @property
    def step_index(self): return self.state.step_index
    @property
    def timing(self): return self.state.timing
    @property
    def A(self): return self.state.A
    @property
    def capacity(self): return self.state.config.steps
    @property
    def E(self):
        start = dict(self.state.config.parameters)['vacuum'][0]
        return self.state.grid_E[start:start+len(self.point_ids)].copy()

    def sample_fields(self):
        """Read detached A(n dt)/E((n-1)dt) samples in point order; H absent."""
        if self.closed or self.failed or self.step_index<1:
            raise InputError('field samples require a usable field after its first step')
        from .field_sampling import SampledField
        start=dict(self.state.config.parameters)['vacuum'][0]
        positions=np.zeros((len(self.point_ids),3))
        positions[:,0]=self.state.x[start:start+len(self.point_ids)]
        return {'A':SampledField('A',[self.step_index*self.dt],positions,self.A[None,:,:],'A/c:au'),
                'E':SampledField('E',[(self.step_index-1)*self.dt],positions,self.E[None,:,:],'E:au')}

    def advance(self, j):
        if self.closed or self.failed:
            raise InputError('field is not usable')
        return self.state.step(j)

    def close(self):
        self.closed = True


class CoupledStepper:
    """Collective outer schedule; material implementations remain replaceable.

    Construct ResponseSession on every rank; pass NativeEHField on field ranks
    and None on material ranks, or use from_multiscale. Call step everywhere. Close the response, field,
    and ParallelContext in that order. A step failure invalidates the response.
    """
    @classmethod
    def from_multiscale(cls, context, bindings, field=None, *, point_record=None, point_map=None):
        """Collective factory. Only the single field rank supplies a fresh state.

        Validate A(0) before creating any material state. All ranks call step,
        response.close and field/context cleanup as for native E/H coupling.
        """
        from .material_response import ResponseSession
        adapter = None; error = None; dt = None
        try:
            if context.group.role == 'field':
                adapter = MultiscaleField(context, field)
                dt = adapter.dt
            elif field is not None:
                raise InputError('material ranks must not provide a field')
        except Exception as exc:
            error = str(exc)
        errors = context.solver_comm.allgather(error)
        local = next((e for e in errors if e), None)
        failures = context.leaders.allgather(local) if context.is_leader else None
        failures = context.solver_comm.bcast(failures, root=0)
        if any(failures):
            raise InputError(next(e for e in failures if e))
        if context.is_leader:
            dt = context.leaders.bcast(dt, root=context.field_group)
        dt = context.solver_comm.bcast(dt, root=0)
        response = ResponseSession(context, bindings, dt=dt, point_record=point_record, point_map=point_map)
        try:
            return cls(response, adapter)
        except BaseException:
            response.close()
            raise

    def __init__(self,response,field=None):
        self.response=response;self.field=field;self.step_index=0
        context=response.context;error=None
        if response.closed or response.failed or response.step_index!=0:
            error='coupling requires a fresh response session'
        if context.group.role=='field':
            if not isinstance(field,(NativeEHField,MultiscaleField)) or field.context is not context:
                error='matching field adapter required on field ranks'
            elif field.closed or field.failed or field.step_index!=0:
                error='coupling requires a fresh field'
            elif field.dt!=response.dt or field.point_ids!=context.plan.point_ids:
                error='field/response dt and point order must match'
        elif field is not None:error='material ranks must not provide a field'
        if getattr(response,'_field_contracts',None) is not None:
            error='field packet sessions require explicit sampling; CoupledStepper currently requires legacy A response'
        response.agree(error)
        self.timing=self._check_timing(response,field)
        self.J=np.zeros((len(context.plan.point_ids),3))

    @staticmethod
    def _check_timing(response,field):
        """Collect every point's matter declaration and check it against the field's.

        The field adapter declares the electronic clock offset (0 or 1/2).
        Returns {'field', 'matter' (by point), 'electronic_clock_offset'} on all ranks.
        """
        context=response.context
        matter = response.timing_records()
        result=None;error=None
        try:
            if context.group.role=='field':
                if set(matter)!=set(context.plan.point_ids):raise InputError('timing missing for some points')
                for point,record in matter.items():
                    try:
                        spec = TimeSpec.from_dict(record)
                        if spec.output('j') is None:
                            raise InputError('response service requires physical charge current j')
                        check_loop(field.timing,spec.shifted(field.electronic_clock_offset))
                    except ValueError as exc:raise InputError(f'point {point}: {exc}') from None
                result={'field':field.timing.to_dict(),'matter':matter,'electronic_clock_offset':field.electronic_clock_offset}
        except Exception as exc:error=str(exc)
        response.agree(error)
        if context.is_leader:result=context.leaders.bcast(result,root=context.field_group)
        return context.solver_comm.bcast(result,root=0)

    def step_protected(self, targets):
        """Opt-in mixed step with explicit full-state/side-effect rollback targets.

        Provide exactly one StepTarget per local field/material owner. No native
        capability is inferred. Known callback failures must return quiescently
        everywhere. Retry only when StepRecoveryError.original_preserved is true
        and cleanup_errors is empty. Ordinary step retains its old behavior.
        """
        from .stepping import StepTarget, StepRecoveryError, guarded_step
        r=self.response;error=None;items=()
        try:
            if r.closed or r.failed or r.context.closed:
                raise InputError('usable response required for protected step')
            items=tuple(targets)
            expected=([self.field] if self.field is not None else [])+list(r.states.values())
            if any(not isinstance(t,StepTarget) for t in items):
                raise InputError('explicit StepTarget required for each local owner')
            if (len(items)!=len(expected) or len({id(t.owner) for t in items})!=len(items)
                    or {id(t.owner) for t in items}!={id(s) for s in expected}):
                raise InputError('protected step must cover exactly every local field/material owner')
        except Exception as exc:error=f'{type(exc).__name__}: {exc}'
        r.agree(error)

        def capture():
            m=r._measurements
            arrays={'J':(self.J,self.J.copy()),'A':(r._A,r._A.copy())}
            adapters={p:(a.A,a.A.copy()) for p,a in r._adapters.items()}
            return {'step':self.step_index,'response_step':r.step_index,'failed':r.failed,
                    'arrays':arrays,'adapters':adapters,
                    'seconds':(m.seconds,dict(m.seconds)),
                    'calls':(m.calls,dict(m.calls)),'errors':(m.errors,list(m.errors))}

        def rollback(token):
            # Preserve original array/dictionary aliases as well as values.
            self.step_index=token['step'];r.step_index=token['response_step']
            r.failed=token['failed']
            for original,values in token['arrays'].values():
                if not np.array_equal(original,values):original[:]=values
            self.J=token['arrays']['J'][0];r._A=token['arrays']['A'][0]
            for p,(original,values) in token['adapters'].items():
                if not np.array_equal(original,values):original[:]=values
                r._adapters[p].A=original
            for name in ('seconds','calls'):
                original,values=token[name]
                original.clear();original.update(values)
                setattr(r._measurements,name,original)
            original,values=token['errors'];original[:]=values
            r._measurements.errors=original

        metadata=StepTarget('__coupling_metadata__',self,capture,rollback)
        try:
            return guarded_step(self.step,(metadata,)+items,context=r.context)
        except StepRecoveryError as exc:
            if not exc.original_preserved or exc.cleanup_errors:
                r.failed=True
                if self.field is not None:self.field.failed=True
            raise

    def step(self):
        r=self.response;context=r.context;A=None;error=None
        if r.closed or r.failed:error='response is not usable'
        if r.step_index!=self.step_index or (self.field is not None and self.field.step_index!=self.step_index):
            error='field and response histories are not synchronized'
        try:r.agree(error)
        except Exception:
            r.failed=True
            raise
        if self.field is not None:
            try:A=self.field.advance(self.J)
            except Exception as exc:error=str(exc)
        try:r.agree(error)
        except Exception:
            r.failed=True
            raise
        try:J=r.step(A if context.group.role=='field' and context.is_leader else None)
        except BaseException:
            r.failed=True
            raise
        self.step_index+=1
        if self.field is not None:
            previous=self.J.copy();self.J=J
            n=self.field.step_index;dt=self.field.dt
            return {'field_time':n*dt,'electronic_local_time':n*dt,
                    'electronic_field_clock':(n+float(Fraction(self.timing['electronic_clock_offset'])))*dt,'A_over_c':A.copy(),
                    'E_at_points':self.field.E.copy(),
                    'E_time':(n-1 if isinstance(self.field,MultiscaleField) else n)*dt,
                    'J_previous':previous,'J_charge':J.copy()}
        return None


class PacketCoupledStepper:
    """Exact native samples -> explicit packet current -> next field update.

    All context ranks participate. No interpolation, implicit dipole conversion
    or whole-session rollback. Initial feedback is zero, as in CoupledStepper.
    """
    def __init__(self,response,field=None):
        from .field_input import FieldContract
        self.response=response;self.field=field;self.step_index=0
        cxt=response.context;error=None;record=None
        try:
            contracts=response._field_contracts
            if response.closed or response.failed or response.step_index!=0 or contracts is None:
                raise InputError('fresh packet response session required')
            if cxt.group.role=='field':
                if not isinstance(field,(NativeEHField,MultiscaleField)) or field.context is not cxt:
                    raise InputError('matching native/1D field required')
                if field.closed or field.failed or field.step_index!=0 or field.dt!=response.dt or field.point_ids!=cxt.plan.point_ids:
                    raise InputError('fresh field with matching dt/point order required')
                native=isinstance(field,NativeEHField)
                if native:
                    positions=field.sample_positions
                    phase=.5
                    available={'A':(.5,0),'E':(0,0),'H':(.5,0)}
                else:
                    start=dict(field.state.config.parameters)['vacuum'][0]
                    positions=np.zeros((len(field.point_ids),3))
                    positions[:,0]=field.state.x[start:start+len(field.point_ids)]
                    phase=0.
                    available={'A':(0,0),'E':(0,1)}
                needs_H=False
                for i,point in enumerate(field.point_ids):
                    contract=contracts[point]
                    if (contract.position_unit!='bohr' or not np.array_equal(contract.position,positions[i])
                            or contract.timing.t0!=0 or contract.timing.clock_offset!=0):
                        raise InputError('packet coupling requires exact mapped position and zero local clock origin')
                    output=contract.timing.output('j')
                    if output is None or float(output.offset)!=phase or output.sampling!='point':
                        raise InputError('current output phase differs from field feedback schedule')
                    for port in contract.timing.inputs:
                        if (port.quantity not in available or (float(port.offset),port.lag)!=available[port.quantity]
                                or port.sampling!='point' or contract.units[port.quantity]!={'A':'A/c:au','E':'E:au','H':'H:au'}[port.quantity]):
                            raise InputError('port clock/unit is unavailable without explicit interpolation')
                        needs_H |= port.quantity=='H'
                if needs_H:
                    field.native._live()
                    if getattr(field.lib,'sp_eh_sample_h',None) is None:
                        raise InputError('H-enabled runtime required before packet coupling')
                record={'schema':'salmonpy.packet-coupling/1','current_offset':phase,'needs_H':needs_H,
                        'field_kind':'native-EH' if native else 'Multiscale1D',
                        'initial_feedback':'zero; not a measured ground-state current',
                        'positions_bohr':positions.tolist(),
                        'point_contracts':{p:contracts[p].timing.to_dict() for p in field.point_ids}}
            elif field is not None:raise InputError('material ranks must not supply field')
        except Exception as exc:error=f'{type(exc).__name__}: {exc}'
        response.agree(error)
        if cxt.is_leader:record=cxt.leaders.bcast(record,root=cxt.field_group)
        self.timing=cxt.solver_comm.bcast(record,root=0)
        self.J=np.zeros((len(cxt.plan.point_ids),3))

    def step(self):
        from .field_input import FieldPacket,FieldValue
        r=self.response;cxt=r.context;error=None;packets=None;samples=None
        try:
            if r.closed or r.failed or r.step_index!=self.step_index:
                raise InputError('packet response is unusable or out of sync')
            if self.field is not None:
                if self.field.closed or self.field.failed or self.field.step_index!=self.step_index or self.step_index>=self.field.capacity:
                    raise InputError('packet field is unusable, exhausted or out of sync')
                if self.timing['needs_H']:self.field.native._live()
        except Exception as exc:error=f'{type(exc).__name__}: {exc}'
        r.agree(error)
        try:
            error=None
            if self.field is not None:
                try:
                    self.field.advance(self.J)
                    samples=self.field.sample_fields()
                except Exception as exc:error=f'{type(exc).__name__}: {exc}'
            r.agree(error)
            # Do not enter the H collective if any field rank failed local sampling.
            error=None
            if self.field is not None and self.timing['needs_H']:
                try:samples['H']=self.field.sample_magnetic_field()
                except Exception as exc:error=f'{type(exc).__name__}: {exc}'
            r.agree(error)
            error=None
            if self.field is not None and cxt.is_leader:
                try:
                    packets={}
                    for i,point in enumerate(self.field.point_ids):
                        components={}
                        for q in r._field_contracts[point].units:
                            sample=samples[q]
                            components[q]=FieldValue(sample.values[0,i],time=sample.times[0],
                                position=sample.positions[i],unit=sample.unit)
                        packets[point]=FieldPacket(**components)
                except Exception as exc:error=f'{type(exc).__name__}: {exc}'
            r.agree(error)
            current=r.step_fields(packets)
        except BaseException:
            r.failed=True
            raise
        self.step_index+=1
        if self.field is not None:
            previous=self.J.copy();self.J=current
            return {'field_time':self.step_index*r.dt,
                    'current_time':(self.step_index+self.timing['current_offset'])*r.dt,
                    'samples':samples,'J_previous':previous,'J_charge':current.copy()}
        return None
