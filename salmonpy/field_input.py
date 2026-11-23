"""Explicit local A/c, E and H inputs; no conversions or MPI initialization."""
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np
from .errors import InputError
from .timing import TimeSpec, SAMPLINGS

__all__ = ['FieldValue', 'FieldPacket', 'FieldContract', 'FieldAdapter']
_QUANTITIES = ('A','E','H')


def _label(value, name):
    if type(value) is not str or not value or value != value.strip():
        raise InputError(f'{name} must be an explicit nonempty label')
    return value


def _real(value, shape, name):
    try:
        raw = np.asarray(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise InputError(f'invalid {name} numeric data') from exc
    if raw.shape != shape or raw.dtype.kind not in 'iuf' or not np.isfinite(raw).all():
        raise InputError(f'{name} must be finite real numeric data of shape {shape}')
    a = np.asarray(raw, dtype=float)
    if not np.isfinite(a).all():
        raise InputError(f'{name} overflows float64')
    return np.frombuffer(a.tobytes(), dtype=float).reshape(shape)


def _index(value):
    if type(value) is not int or value < 0:
        raise InputError('step_index must be a nonnegative integer')
    return value


@dataclass(frozen=True, slots=True, eq=False)
class FieldValue:
    """One real 3-vector at a declared time (atomic units) and spatial point.

    unit/position_unit are exact labels, not requests for automatic conversion.
    H denotes magnetic field intensity, never B. Packet.A denotes A/c.
    """
    value: object
    time: float
    position: object
    unit: str
    position_unit: str = 'bohr'
    sampling: str = 'point'

    def __post_init__(self):
        object.__setattr__(self, 'value', _real(self.value, (3,), 'field'))
        object.__setattr__(self, 'time', float(_real(self.time, (), 'time')))
        object.__setattr__(self, 'position', _real(self.position, (3,), 'position'))
        _label(self.unit, 'unit')
        _label(self.position_unit, 'position_unit')
        if self.sampling not in SAMPLINGS:
            raise InputError('unsupported field sampling')

    @classmethod
    def _from_transport(cls,value,*,time,position,unit,position_unit='bohr',sampling='point',checked=True):
        """Owned transport value; performance skips vector scans, not metadata."""
        if checked:return cls(value,time,position,unit,position_unit,sampling)
        from .ld import _array_unchecked
        from .sbe import _frozen
        result=object.__new__(cls)
        object.__setattr__(result,'value',_frozen(_array_unchecked(value,'field',(3,)),check=False))
        object.__setattr__(result,'time',float(_real(time,(),'time')))
        object.__setattr__(result,'position',_real(position,(3,),'position'))
        object.__setattr__(result,'unit',_label(unit,'unit'))
        object.__setattr__(result,'position_unit',_label(position_unit,'position_unit'))
        if sampling not in SAMPLINGS:raise InputError('unsupported field sampling')
        object.__setattr__(result,'sampling',sampling)
        return result


@dataclass(frozen=True, slots=True)
class FieldPacket:
    """Optional components: None means absent, distinct from an explicit zero."""
    A: FieldValue | None = None
    E: FieldValue | None = None
    H: FieldValue | None = None

    def __post_init__(self):
        if any(getattr(self, q) is not None and not isinstance(getattr(self, q), FieldValue) for q in _QUANTITIES):
            raise InputError('A/E/H must be FieldValue or None')


@dataclass(frozen=True, slots=True, eq=False)
class FieldContract:
    """Required ports, declared units/gauge and one model sample position.

    TimeSpec supplies absolute times: t0+(n+clock_offset+offset-lag)*dt.
    Caller explicitly interpolates spatially staggered field data to this point.
    """
    timing: TimeSpec
    units: object
    position: object
    gauge: str
    position_unit: str = 'bohr'

    def __post_init__(self):
        if not isinstance(self.timing, TimeSpec):
            raise InputError('explicit TimeSpec required')
        names = tuple(p.quantity for p in self.timing.inputs)
        if not names or any(q not in _QUANTITIES for q in names):
            raise InputError('field inputs must be nonempty A/E/H ports')
        if not isinstance(self.units, Mapping) or set(self.units) != set(names):
            raise InputError('unit declarations must match required field ports exactly')
        units = {q:_label(self.units[q], q+' unit') for q in names}
        object.__setattr__(self, 'units', MappingProxyType(units))
        object.__setattr__(self, 'position', _real(self.position, (3,), 'position'))
        _label(self.position_unit, 'position_unit')
        if self.gauge not in ('velocity','length','phenomenological','custom'):
            raise InputError('unsupported gauge declaration')
        if self.gauge == 'velocity' and 'A' not in names:
            raise InputError('velocity gauge requires A/c')
        if self.gauge == 'length' and 'E' not in names:
            raise InputError('length gauge requires E')

    def validate(self, fields, *, step_index):
        n = _index(step_index)
        if not isinstance(fields, FieldPacket):
            raise InputError('FieldPacket required')
        for port in self.timing.inputs:
            if n == 0 and not port.uses_t0:
                continue
            component = getattr(fields, port.quantity)
            if component is None:
                raise InputError(f'missing required {port.quantity} field')
            expected = self.timing.t0 + (n + float(self.timing.clock_offset + port.offset) - port.lag)*self.timing.dt
            tolerance = min(self.timing.dt*1e-3, max(self.timing.dt*1e-9, 8*np.finfo(float).eps*max(abs(expected),self.timing.dt)))
            if not np.isfinite(expected) or abs(component.time-expected) > tolerance:
                raise InputError(f'{port.quantity} time differs from declared port clock')
            if component.unit != self.units[port.quantity] or component.sampling != port.sampling:
                raise InputError(f'{port.quantity} units or sampling differ')
            if component.position_unit != self.position_unit or not np.array_equal(component.position, self.position):
                raise InputError(f'{port.quantity} sample position differs; interpolate explicitly')


class FieldAdapter:
    """Borrow a state, validating input before explicit legacy or multiport calls.

    legacy_quantity A/E delegates step(vector) and preserves the return value.
    Otherwise the backend must implement step_fields(packet). This is a local
    adapter, not a native isolated-system implementation or MPI field transport.
    """
    def __init__(self, state, contract, *, legacy_quantity=None, step_index=0):
        if not isinstance(contract, FieldContract):
            raise InputError('FieldContract required')
        if legacy_quantity is not None:
            if legacy_quantity not in ('A','E') or tuple(contract.units) != (legacy_quantity,):
                raise InputError('legacy adapter requires exactly its explicit A or E port')
            port = contract.timing.inputs[0]
            canonical = {'A':'A/c:au', 'E':'E:au'}[legacy_quantity]
            if (contract.units[legacy_quantity] != canonical or port.offset != 0 or
                    port.lag != 0 or port.sampling != 'point' or contract.timing.clock_offset != 0):
                raise InputError('legacy step requires canonical atomic units at the new endpoint')
            callback = getattr(state, 'step', None)
        else:
            callback = getattr(state, 'step_fields', None)
        if not callable(callback):
            raise InputError('backend lacks the selected stepping method')
        self._state = state
        self._contract = contract
        self._quantity = legacy_quantity
        self._steps = _index(step_index)
        self._failed = False
        self._check_backend_state()

    def _check_backend_state(self, step_index=None):
        n = self._steps if step_index is None else step_index
        timing = self._contract.timing
        raw_timing = getattr(self._state, 'timing', None)
        if raw_timing is not None:
            if (not isinstance(raw_timing, TimeSpec) or raw_timing.inputs != timing.inputs or
                    raw_timing.dt != timing.dt or raw_timing.t0 != timing.t0 or
                    raw_timing.clock_offset != timing.clock_offset):
                raise InputError('backend time/port declaration differs from field contract')
        dt = getattr(self._state, 'dt', timing.dt)
        if float(_real(dt, (), 'backend dt')) != timing.dt:
            raise InputError('backend dt differs')
        expected = timing.t0+(n+float(timing.clock_offset))*timing.dt
        time = getattr(self._state, 'time', expected)
        tolerance = min(timing.dt*1e-3, max(timing.dt*1e-9, 8*np.finfo(float).eps*max(abs(expected),timing.dt)))
        if not np.isfinite(expected) or abs(float(_real(time, (), 'backend time'))-expected) > tolerance:
            raise InputError('backend clock differs; rebuild adapter after external advance/restore')

    @property
    def contract(self): return self._contract

    @property
    def step_index(self): return self._steps

    def step_fields(self, fields):
        if self._failed:
            raise InputError('adapter invalidated by backend failure')
        self._check_backend_state()
        self._contract.validate(fields, step_index=self._steps+1)
        try:
            if self._quantity is not None:
                result = self._state.step(getattr(fields, self._quantity).value.copy())
            else:
                selected = FieldPacket(**{q:getattr(fields,q) for q in self._contract.units})
                result = self._state.step_fields(selected)
            self._check_backend_state(self._steps+1)
        except BaseException:
            self._failed = True
            raise
        self._steps += 1
        return result
