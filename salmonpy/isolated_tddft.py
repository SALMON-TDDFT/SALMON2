"""Explicit isolated PZ length-gauge E -> charge-dipole process interface."""
from dataclasses import dataclass, field
import json
from pathlib import Path
import numpy as np
from .errors import InputError, NumericalError
from .fdtd import SalmonExecutable
from .field_input import FieldContract, FieldPacket, FieldValue, _real
from .timing import TimeSpec, InputPort
from .snapshots import _fingerprint
from .tddft import TDDFTInitialState, TDDFTState, _digest

__all__ = ['IsolatedRuntime','load_isolated_runtime','DipoleObservables',
           'IsolatedTDDFTInitialState','IsolatedTDDFTState']
_CAPABILITIES = {'isolated','pz','fixed-ions','middlepoint','single-rank','E-to-dipole','measured-electron-number'}


def _receipt(payload, executable):
    if type(payload) is not bytes or len(payload) > 65536:
        raise InputError('bounded isolated-worker validation bytes required')
    def unique(pairs):
        result = {}
        for k,v in pairs:
            if k in result: raise InputError('duplicate receipt key')
            result[k] = v
        return result
    try:
        data = json.loads(payload, object_pairs_hook=unique)
        valid = (data['schema']=='salmonpy.isolated-worker-validation/2'
                 and data['protocol']=='SALMONPY_TDDFT_E_V2'
                 and data['executable_sha256']==executable.sha256
                 and type(data['steps']) is int and data['steps']>=20
                 and data['interleaved_independence'] is True
                 and set(data['capabilities'])==_CAPABILITIES)
        for key in ('repeat_max_abs_error','replay_continuation_max_abs_error','reference_max_abs_error'):
            value=data[key]
            valid = valid and type(value) in (int,float) and np.isfinite(value) and 0 <= value <= 1.e-10
    except (KeyError, TypeError, ValueError, RecursionError) as exc:
        raise InputError('invalid isolated-worker validation receipt') from exc
    if not valid or _digest(executable.path)!=executable.sha256:
        raise InputError('isolated E-worker evidence or binary hash does not match')


@dataclass(frozen=True, slots=True)
class IsolatedRuntime:
    """Locally validated E2 binary; a receipt is evidence, not a signature or audit."""
    executable: SalmonExecutable
    validation: bytes = field(repr=False)

    def __post_init__(self):
        if not isinstance(self.executable,SalmonExecutable):
            raise InputError('SalmonExecutable required')
        _receipt(self.validation,self.executable)


def load_isolated_runtime(executable, validation):
    """Verify local E2 validation JSON and exact executable hash; never build/download."""
    binary = Path(executable).expanduser().resolve(strict=True)
    path = Path(validation).expanduser()
    with path.open('rb') as stream:
        payload = stream.read(65537)
    return IsolatedRuntime(SalmonExecutable(binary,_digest(binary)),payload)


@dataclass(frozen=True, slots=True, eq=False)
class DipoleObservables:
    time: float
    electronic_change: object
    dipole: object
    electron_number: float

    def __post_init__(self):
        object.__setattr__(self,'time',float(_real(self.time,(),'time')))
        object.__setattr__(self,'electronic_change',_real(self.electronic_change,(3,),'electronic dipole change'))
        object.__setattr__(self,'dipole',_real(self.dipole,(3,),'total charge dipole'))
        object.__setattr__(self,'electron_number',float(_real(self.electron_number,(),'electron number')))


class _EEngine(TDDFTState):
    _worker_mode = 'E2'
    _reply_prefix = 'SALMONPY_TDDFT_E_V2 '
    _init_command = 'INIT_E'
    _step_command = 'STEP_E'
    _reply_size = 7
    _snapshot_kind = 'isolated-TDDFT-E-dipole-replay'

    def __init__(self,model,dt,position):
        self._position = position
        super().__init__(model,dt,np.zeros(3))

    @property
    def current(self):
        # Private transport hook: seven raw observables, never exposed as current.
        return _real(self._current,(7,),'dipole/electron observables')

    def _receive(self,index):
        try:
            return super()._receive(index)
        except InputError as exc:
            raise NumericalError('isolated worker returned nonfinite observables') from exc

    def _snapshot_settings(self):
        return _fingerprint((super()._snapshot_settings(),tuple(self._position)))


@dataclass(frozen=True, slots=True, eq=False)
class IsolatedTDDFTInitialState:
    """Owned native GS snapshot; states lease it and own independent processes."""
    _model: TDDFTInitialState = field(repr=False)
    runtime: IsolatedRuntime

    @classmethod
    def from_salmon(cls,directory,input_file,*,runtime,max_steps=1000,threads=1,timeout=300.):
        if not isinstance(runtime,IsolatedRuntime):
            raise InputError('validated IsolatedRuntime required; periodic A runtime is not sufficient')
        _receipt(runtime.validation,runtime.executable)
        model=TDDFTInitialState.from_salmon(directory,input_file,executable=runtime.executable,
                                          max_steps=max_steps,threads=threads,timeout=timeout)
        return cls(model,runtime)

    def new_state(self,*,dt,position=(0.,0.,0.)):
        dt=float(_real(dt,(),'dt'))
        if dt<=0: raise InputError('positive dt required')
        position=_real(position,(3,),'position')
        return IsolatedTDDFTState(_EEngine(self._model,dt,position),position)

    def close(self): self._model.close()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()


class IsolatedTDDFTState:
    """E-driven dipole state. In-memory replay restore, not a current-density model."""
    def __init__(self,engine,position):
        self._engine=engine
        self._contract=FieldContract(TimeSpec(engine.dt,(InputPort('E'),),(),
                            integrator='SALMON isolated middlepoint; E at native step index'),
                            units={'E':'E:au'},position=position,gauge='length')

    @property
    def dt(self): return self._engine.dt
    @property
    def time(self): return self._engine.time
    @property
    def timing(self): return self._contract.timing
    @property
    def observables(self):
        values=self._engine.current
        return DipoleObservables(self.time,values[:3],values[3:6],values[6])

    @property
    def electron_number(self): return self.observables.electron_number

    def step(self,E):
        E=_real(E,(3,),'E')
        self._engine.step(E)
        return self.observables

    def step_fields(self,fields):
        self._contract.validate(fields,step_index=self._engine._steps+1)
        return self.step(fields.E.value)

    def snapshot(self): return self._engine.snapshot()
    def restore(self,snapshot): self._engine.restore(snapshot)
    def save_checkpoint(self,directory):
        """Save accepted E/dipole replay history; GS files remain separate."""
        self._engine.save_checkpoint(directory)

    def restore_checkpoint(self,directory,*,max_bytes=None):
        """Validate matching initial data/settings and atomically replay disk history."""
        self._engine.restore_checkpoint(directory,max_bytes=max_bytes)

    def close(self): self._engine.close()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()
