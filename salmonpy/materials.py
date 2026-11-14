"""Author one endpoint-A material model for serial and point-distributed use.

No MPI import/initialization or implicit conformance run. The checker exercises
only the explicitly supplied drive and cannot establish a material law's physics.
"""
from dataclasses import dataclass, field
import sys
import numpy as np
from .errors import InputError, InputTypeError, RuntimeFailure
from .fields import FieldSamples
from .ld import _array
from .timing import TimeSpec, InputPort, OutputPort, charge_sign

__all__ = ['MaterialCheck', 'bind_material', 'check_material']


@dataclass(frozen=True, slots=True)
class MaterialCheck:
    backend: str
    steps: int
    dt: float
    maximum_repeat_difference: float
    checks: tuple = ('endpoint timing and clocks', 'finite three-vector currents',
                     'input ownership', 'independent states', 'repeatability')
    not_tested: tuple = ('physical accuracy or units', 'all possible inputs',
                         'failure atomicity', 'model immutability', 'MPI collectives')


def _model(initial):
    if not isinstance(getattr(initial,'backend',None),str) or not initial.backend.strip():
        raise InputError('material requires a nonempty backend identifier')
    if not callable(getattr(initial,'new_state',None)):
        raise InputTypeError('material requires new_state(dt, time0, initial_A)')


def _spec(state, dt, time0):
    spec=getattr(state,'timing',None)
    if not isinstance(spec,TimeSpec):raise InputTypeError('explicit material TimeSpec required')
    if not np.isclose(spec.dt,dt,rtol=1e-12,atol=0) or not np.isclose(spec.t0,time0,rtol=0,atol=abs(dt)*1e-9):
        raise InputError('material timing dt or time origin differs')
    inp=spec.input('A')
    outputs=[p for p in spec.outputs if p.quantity in ('j','J_matter')]
    if (len(spec.inputs)!=1 or inp is None or len(spec.outputs)!=1 or len(outputs)!=1
            or spec.clock_offset+inp.offset!=0 or inp.sampling!='point'
            or spec.clock_offset+outputs[0].offset!=0 or outputs[0].sampling!='point'
            or outputs[0].at_t0!='physical'):
        raise InputError('material must exchange endpoint A and physical endpoint current')
    return spec


def _clock(state, expected, dt, step, time0):
    value=float(_array(state.time,'material time',()))
    scale=max(abs(expected),abs(time0),abs(dt))
    tolerance=min(abs(dt)*1e-3,max(abs(dt)*1e-8,8*np.finfo(float).eps*(step+1)*scale))
    if abs(value-expected)>tolerance:raise InputError('material time differs from endpoint clock')


def _close(state):
    close=getattr(state,'close',None)
    if close is not None:close()


@dataclass
class _Prepared:
    initial: object
    live: list = field(default_factory=list)


class _EndpointState:
    def __init__(self,raw,prepared,dt,time0):
        self._raw=raw;self._prepared=prepared;self._dt=dt;self._time0=time0
        self._steps=0;self._closed=False;self._failed=False
        self._spec=_spec(raw,dt,time0);self._sign=charge_sign(self._spec)
        _clock(raw,time0,dt,0,time0)
        if not callable(getattr(raw,'step',None)):raise InputTypeError('material state requires step(A)')
        self._initial_current=self.current.copy()
        self.timing=TimeSpec(dt,self._spec.inputs,(OutputPort('j'),),t0=time0,
            history=self._spec.history,integrator='serial material adapter: '+self._spec.integrator)

    @property
    def time(self):return self._raw.time
    @property
    def current(self):return self._sign*_array(self._raw.current,'material current',(3,))
    @property
    def initial_current(self):return self._initial_current.copy()

    def step(self,A):
        if self._closed or self._failed:raise InputError('material state is closed or failed')
        values=_array(A,'A',(3,));owned=values.copy()
        try:
            result=self._sign*_array(self._raw.step(owned),'returned current',(3,))
            if not np.array_equal(owned,values):raise InputError('material mutated its input A')
            if getattr(self._raw,'timing',None)!=self._spec:raise InputError('material timing changed during evolution')
            _clock(self._raw,self._time0+(self._steps+1)*self._dt,self._dt,self._steps+1,self._time0)
            if not np.allclose(result,self.current,rtol=1e-12,atol=1e-14):
                raise InputError('returned material current differs from state.current')
        except BaseException:
            self._failed=True
            raise
        self._steps+=1
        return result.copy()

    def close(self):
        if self._closed:return
        self._closed=True
        try:_close(self._raw)
        finally:self._prepared.live[:]=[s for s in self._prepared.live if s is not self._raw]


def _create(prepared,dt,time0,initial_A):
    values=_array(initial_A,'initial_A',(3,));owned=values.copy()
    raw=prepared.initial.new_state(dt=dt,time0=time0,initial_A=owned)
    if any(raw is old for old in prepared.live):
        # The existing owner will close it; never destroy an earlier live point.
        raise InputError('material reused one mutable state for multiple points')
    try:
        if not np.array_equal(owned,values):raise InputError('material mutated initial_A')
        state=_EndpointState(raw,prepared,dt,time0)
    except BaseException as exc:
        try:_close(raw)
        except Exception as cleanup:
            if hasattr(exc,'add_note'):exc.add_note('material cleanup also failed: '+str(cleanup))
            else:exc.material_cleanup_error=str(cleanup)
        raise
    prepared.live.append(raw)
    return state


@dataclass(frozen=True)
class _Provider:
    initial: object
    identity: str
    def prepare(self,context):return _Prepared(self.initial)


class _Backend:
    identity='salmonpy.serial-material.endpoint-charge/1'
    execution='leader_only'
    reduction='leader'
    def create_state(self,prepared,point_id,dt,context):
        if context is not None and context.solver_comm.size!=1:
            raise InputError('serial material adapter requires one rank per material group')
        return _create(prepared,dt,0.,np.zeros(3))


def bind_material(initial, *, key, identity):
    """Adapt a serial crystal model to existing MaterialBinding point exchange.

    identity is the author's stable algorithm/parameter/data identity; it is not
    automatically inferred from Python object addresses. Reuse one binding for
    all points with its key. One MPI rank per material group; no internal MPI.
    """
    _model(initial)
    for label,value in (('key',key),('identity',identity)):
        if not isinstance(value,str) or not value.strip():raise InputError('nonempty material '+label+' required')
    from .tddft import TDDFTInitialState
    if isinstance(initial,TDDFTInitialState):
        raise InputError('use TDDFTMaterialProvider and TDDFTResponseBackend for managed native workers')
    import hashlib,json
    declared=hashlib.sha256(json.dumps([initial.backend,identity],ensure_ascii=True).encode()).hexdigest()
    from .material_response import MaterialBinding
    return MaterialBinding(key,_Provider(initial,'salmonpy.serial-material/'+declared),_Backend())


def check_material(initial, drive, *, rtol=1e-12, atol=1e-14):
    """Exercise two independent states with the caller's representative waveform.

    Returns a MaterialCheck on success; closes only states created by this call.
    It does not prove physical accuracy, failure atomicity or MPI collectives.
    Never invoked implicitly by bind_material or a production session.
    """
    _model(initial)
    if not isinstance(drive,FieldSamples) or drive.kind!='A':raise InputError('FieldSamples(kind="A") required')
    deltas=np.diff(drive.time);dt=float(deltas[0]);time0=float(drive.time[0])
    if not np.allclose(deltas,dt,rtol=1e-8,atol=1e-12*abs(dt)):raise InputError('uniform A samples required')
    tolerances={}
    for name,value in (('rtol',rtol),('atol',atol)):
        raw=np.asarray(value)
        if raw.shape!=() or raw.dtype.kind not in 'iuf':raise InputError('real scalar '+name+' required')
        parsed=float(raw)
        if not np.isfinite(parsed) or parsed<0:raise InputError('finite nonnegative '+name+' required')
        tolerances[name]=parsed
    rtol,atol=tolerances['rtol'],tolerances['atol']
    prepared=_Prepared(initial);states=[];maximum=0.
    try:
        a=_create(prepared,dt,time0,drive.values[0]);states.append(a)
        before=a.current.copy();clock=a.time
        b=_create(prepared,dt,time0,drive.values[0]);states.append(b)
        if a.time!=clock or not np.array_equal(a.current,before):
            raise InputError('creating a material state changed another state')
        if a.timing!=b.timing:raise InputError('material declarations differ between independent states')
        if not np.allclose(a.current,b.current,rtol=rtol,atol=atol):raise InputError('initial material currents are not reproducible')
        maximum=float(np.max(np.abs(a.current-b.current)))
        for index,A in enumerate(drive.values[1:],1):
            before=b.current.copy();clock=b.time
            one=a.step(A)
            _clock(a,float(drive.time[index]),dt,index,time0)
            if b.time!=clock or not np.array_equal(b.current,before):raise InputError('material states share mutable history')
            before=a.current.copy();clock=a.time
            two=b.step(A)
            _clock(b,float(drive.time[index]),dt,index,time0)
            if a.time!=clock or not np.array_equal(a.current,before):raise InputError('material states share mutable history')
            if not np.allclose(one,two,rtol=rtol,atol=atol):raise InputError('material evolution is not reproducible')
            maximum=max(maximum,float(np.max(np.abs(one-two))))
        return MaterialCheck(initial.backend,len(drive.time)-1,dt,maximum)
    finally:
        active_error=sys.exc_info()[1];errors=[]
        for state in states:
            try:state.close()
            except Exception as exc:errors.append(str(exc))
        if errors:
            message='material checker cleanup failed: '+'; '.join(errors)
            if active_error is None:raise RuntimeFailure(message)
            if hasattr(active_error,'add_note'):active_error.add_note(message)
            else:active_error.material_cleanup_error=message
