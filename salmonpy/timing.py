"""Time-layout declarations (schema ``salmonpy.time/1``) and coupling checks.

Every component declares where its exchanged quantities live in time. A port's
time is t0 + (n + clock_offset + offset)*dt; offsets are 0 or 1/2 so that
comparisons are exact. ``check_loop`` rejects couplings whose declarations do
not line up, instead of silently shifting by half a step.
"""
from dataclasses import dataclass, field, asdict
from fractions import Fraction
import warnings
import numpy as np

SCHEMA = 'salmonpy.time/1'
QUANTITIES = ('A', 'E', 'H', 'j', 'J_matter')
SAMPLINGS = ('point', 'interval_average', 'centered_difference')
AT_T0 = ('physical', 'zero')

__all__ = ['SCHEMA', 'InputPort', 'OutputPort', 'TimeSpec', 'crystal_default', 'matter_timing',
           'charge_sign', 'check_loop']


def _half(value, name):
    if isinstance(value, (bool, np.bool_)): raise ValueError(f'{name} must be 0 or 1/2')
    try:
        f = Fraction(value)  # exact for int, float, Fraction and strings such as '1/2'
    except (TypeError, ValueError, ZeroDivisionError):
        raise ValueError(f'{name} must be 0 or 1/2') from None
    if f not in (0, Fraction(1, 2)): raise ValueError(f'{name} must be 0 or 1/2')
    return f


def _choice(value, allowed, name):
    if value not in allowed: raise ValueError(f'{name} must be one of {allowed}')
    return value


@dataclass(frozen=True, slots=True)
class InputPort:
    quantity: str
    offset: Fraction = Fraction(0)
    sampling: str = 'point'
    lag: int = 0
    uses_t0: bool = True

    def __post_init__(self):
        _choice(self.quantity, QUANTITIES, 'quantity')
        object.__setattr__(self, 'offset', _half(self.offset, 'offset'))
        _choice(self.sampling, SAMPLINGS, 'sampling')
        if isinstance(self.lag, (bool, np.bool_)) or not isinstance(self.lag, (int, np.integer)) or self.lag < 0:
            raise ValueError('lag must be a nonnegative integer')
        object.__setattr__(self, 'lag', int(self.lag))
        if not isinstance(self.uses_t0, (bool, np.bool_)): raise ValueError('uses_t0 must be boolean')
        object.__setattr__(self, 'uses_t0', bool(self.uses_t0))


@dataclass(frozen=True, slots=True)
class OutputPort:
    quantity: str
    offset: Fraction = Fraction(0)
    sampling: str = 'point'
    at_t0: str = 'physical'

    def __post_init__(self):
        _choice(self.quantity, QUANTITIES, 'quantity')
        object.__setattr__(self, 'offset', _half(self.offset, 'offset'))
        _choice(self.sampling, SAMPLINGS, 'sampling')
        _choice(self.at_t0, AT_T0, 'at_t0')


@dataclass(frozen=True, slots=True)
class TimeSpec:
    """One component's declaration. ``inputs``/``outputs`` map quantity -> port."""
    dt: float
    inputs: tuple
    outputs: tuple
    t0: float = 0.0
    clock_offset: Fraction = Fraction(0)
    history: str = ''
    integrator: str = ''
    schema: str = field(default=SCHEMA, init=False)

    def __post_init__(self):
        dt = float(self.dt)
        if not np.isfinite(dt) or dt <= 0: raise ValueError('dt must be positive and finite')
        object.__setattr__(self, 'dt', dt)
        t0 = float(self.t0)
        if not np.isfinite(t0): raise ValueError('t0 must be finite')
        object.__setattr__(self, 't0', t0)
        object.__setattr__(self, 'clock_offset', _half(self.clock_offset, 'clock_offset'))
        for name, kind in (('inputs', InputPort), ('outputs', OutputPort)):
            ports = tuple(getattr(self, name))
            if not all(isinstance(p, kind) for p in ports): raise TypeError(f'{name} must be {kind.__name__} objects')
            if len({p.quantity for p in ports}) != len(ports): raise ValueError(f'duplicate {name} quantity')
            object.__setattr__(self, name, ports)
        for name in ('history', 'integrator'):
            if not isinstance(getattr(self, name), str): raise ValueError(f'{name} must be text')

    def input(self, quantity):
        return next((p for p in self.inputs if p.quantity == quantity), None)

    def output(self, quantity):
        return next((p for p in self.outputs if p.quantity == quantity), None)

    @classmethod
    def from_dict(cls, record):
        """Inverse of to_dict (used to exchange declarations between MPI ranks)."""
        if record.get('schema') != SCHEMA: raise ValueError('unsupported timing schema')
        return cls(record['dt'], tuple(InputPort(**{**p, 'offset': p['offset']}) for p in record['inputs'].values()),
                   tuple(OutputPort(**p) for p in record['outputs'].values()), t0=record['t0'],
                   clock_offset=record['clock_offset'], history=record['history'], integrator=record['integrator'])

    def shifted(self, clock_offset):
        """Copy placed on another clock: clock_offset is added (result must be 0 or 1/2)."""
        from dataclasses import replace
        return replace(self, clock_offset=self.clock_offset + _half(clock_offset, 'clock_offset'))

    def to_dict(self):
        """JSON-ready record; offsets are written as strings such as '0' and '1/2'."""
        def port(p):
            d = asdict(p); d['offset'] = str(p.offset); return d
        return {'schema': self.schema, 'dt': self.dt, 't0': self.t0, 'clock_offset': str(self.clock_offset),
                'inputs': {p.quantity: port(p) for p in self.inputs},
                'outputs': {p.quantity: port(p) for p in self.outputs},
                'history': self.history, 'integrator': self.integrator}


def crystal_default(dt, t0=0.0):
    """The documented ``crystal`` contract, for matter backends without a declaration."""
    return TimeSpec(dt, (InputPort('A'),), (OutputPort('j'),), t0=t0,
                    integrator='crystal protocol default: step(A at n+1) -> charge current j at n+1')


def matter_timing(state, dt, t0=0.0):
    """The state's declaration, or the crystal default with a warning when it has none."""
    spec = getattr(state, 'timing', None)
    if spec is None:
        warnings.warn(f'{type(state).__name__} has no timing declaration; assuming the crystal '
                      'contract (A at n*dt -> charge current j). Until 2026-11-14 the contract '
                      'was SALMON matter current J = -j: declare timing to silence this.',
                      stacklevel=3)
        return crystal_default(dt, t0)
    if not isinstance(spec, TimeSpec): raise TypeError('timing must be a salmonpy.timing.TimeSpec')
    return spec


def charge_sign(spec):
    """+1 if the matter output is the charge current j, -1 for SALMON matter current."""
    if spec.output('j') is not None: return 1.0
    if spec.output('J_matter') is not None: return -1.0
    raise ValueError('matter component must output j or J_matter')


#: Explicit adapters allowed between an output and an input quantity.
ADAPTERS = {('J_matter', 'j'): 'j = -J_matter (electronic charge -1)'}


def _connect(producer, out_q, consumer, in_q, label):
    out = producer.output(out_q); inp = consumer.input(in_q)
    if out is None: raise ValueError(f'{label}: producer has no output {out_q!r}')
    if inp is None: raise ValueError(f'{label}: consumer has no input {in_q!r}')
    if out_q != in_q and (out_q, in_q) not in ADAPTERS:
        raise ValueError(f'{label}: no adapter from {out_q} to {in_q}')
    t_out = producer.clock_offset + out.offset; t_in = consumer.clock_offset + inp.offset
    if t_out != t_in:
        raise ValueError(f'{label}: {out_q} is produced at offset {t_out} but consumed at offset {t_in}')
    if out.sampling != inp.sampling:
        raise ValueError(f'{label}: sampling {out.sampling} cannot feed {inp.sampling}')
    if inp.uses_t0 and out.at_t0 != 'physical':
        raise ValueError(f'{label}: consumer uses t=0 but the producer gives a placeholder')
    return inp.lag


def check_loop(field_spec, matter_spec, *, field_in='j', field_out='A', matter_in='A', matter_out=None):
    """Validate field <-> matter coupling (rules 1-4 of salmonpy.time/1). Returns adapters used.

    ``matter_out`` defaults to the matter's declared j, else J_matter.
    """
    if matter_out is None:
        matter_out = 'j' if matter_spec.output('j') is not None else 'J_matter'
    if not np.isclose(field_spec.dt, matter_spec.dt, rtol=1e-12, atol=0):
        raise ValueError('field and matter dt differ')
    if not np.isclose(field_spec.t0, matter_spec.t0, rtol=0, atol=1e-9*field_spec.dt):
        raise ValueError('field and matter t0 differ')
    lag = _connect(field_spec, field_out, matter_spec, matter_in, 'field->matter')
    lag += _connect(matter_spec, matter_out, field_spec, field_in, 'matter->field')
    if lag < 1:
        raise ValueError('coupling loop has total lag 0: an implicit solve would be required')
    return tuple(ADAPTERS[(matter_out, field_in)] for _ in [0] if matter_out != field_in)
