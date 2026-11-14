"""Macroscopic field components with the J -> A protocol.

A field component closes the light-matter loop from the field side. Its state
receives the electric (charge) current density j at time n*dt and returns
Ac=A/c at time (n+1)*dt at each macroscopic point. Matter components (TDDFT,
SBE, user models) keep the A -> J protocol of ``salmonpy.crystal``. ``couple``
alternates the two:

    j_0 = matter(A_0);  A_1 = field(j_0);  j_1 = matter(A_1);  A_2 = field(j_1); ...

All quantities are SALMON atomic units in the laboratory frame. A is SALMON
Ac=A/c and E=-dAc/dt. Interfaces exchange the electric current density j.
Matter modules keep SALMON's electronic matter current J internally and return
j = -J (crystal contract). A backend that declares J_matter output is converted
by ``couple`` at exactly one place.

The update formulas, operation order and initial conditions follow SALMON:
``rt/em_field.f90`` (trans_longi = 'tr', 'lo', '2d') and the 1D Weyl-gauge
multiscale solver (``maxwell/fdtd_weyl_gauge.f90``, ``ms/main_ms.f90``). As in
SALMON, the current at t=0 is not fed back: the first step returns A at dt from
the external/incident field alone.
"""
from dataclasses import dataclass, field, asdict
import numpy as np
from .fields import AcosPulse, FieldSamples
from .ld import _array, _positive
from .timing import TimeSpec, InputPort, OutputPort, matter_timing, charge_sign, check_loop

#: SALMON ``cspeed_au`` (inverse fine-structure constant), atomic units.
SPEED_OF_LIGHT = 137.035999084
_PI = 3.141592653589793


__all__ = ['SPEED_OF_LIGHT', 'Transverse', 'Longitudinal', 'ThinFilm', 'Multiscale1D',
           'FieldStateConfig', 'CoupledResponse', 'couple']


def _steps(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not 1 <= value <= 2**31-3:
        raise ValueError('steps must be an integer in [1, 2**31-3]')
    return int(value)


def _drive(value, name):
    """Accept None, an AcosPulse, a sequence of AcosPulse, FieldSamples(kind='A') or a callable."""
    if value is None or isinstance(value, (AcosPulse, FieldSamples)) or callable(value):
        if isinstance(value, FieldSamples) and value.kind != 'A':
            raise ValueError(f'{name} samples must be A, not E')
        return value
    if isinstance(value, (tuple, list)) and value and all(isinstance(p, AcosPulse) for p in value):
        return tuple(value)
    raise TypeError(f'{name} must be None, AcosPulse, a sequence of AcosPulse, FieldSamples(kind="A") or a callable')


def _describe(value):
    if value is None: return 'none'
    if isinstance(value, AcosPulse): return {'AcosPulse': asdict(value)}
    if isinstance(value, tuple): return {'AcosPulse': [asdict(p) for p in value]}
    if isinstance(value, FieldSamples): return 'FieldSamples(kind="A")'
    return 'callable'


def _pulse_values(pulse, t):
    # AcosPulse.sample needs a strictly increasing grid with two or more samples.
    order = np.argsort(t, kind='stable')
    ts = t[order]
    if np.any(np.diff(ts) <= 0):
        raise ValueError('pulse evaluation times must be distinct')
    padded = ts if len(ts) > 1 else np.array([ts[0], ts[0]+1.])
    v = pulse.sample(padded).values[:len(ts)]
    out = np.empty((len(t), 3)); out[order] = v
    return out


def _evaluate(drive, t, *, dt=None):
    """A(t) for the times t (1D). FieldSamples must lie exactly on the n*dt grid."""
    t = np.asarray(t, dtype=float)
    if drive is None:
        return np.zeros((len(t), 3))
    if isinstance(drive, AcosPulse):
        return _pulse_values(drive, t)
    if isinstance(drive, tuple):
        return sum(_pulse_values(p, t) for p in drive)
    if isinstance(drive, FieldSamples):
        if dt is None:
            raise ValueError('FieldSamples are accepted only for uniform (single-point) fields')
        index = np.rint(t/dt).astype(int)
        if np.any(index < 0) or np.any(index >= len(drive.time)):
            raise ValueError(f'external A samples must cover t=0..{t[-1]:.17g} on the dt grid')
        if not np.allclose(drive.time[index], t, rtol=0, atol=1e-9*dt):
            raise ValueError('external A samples must lie on the n*dt grid starting at t=0')
        return np.array(drive.values[index])
    values = np.array(drive(t.copy()), dtype=float)
    if values.shape != (len(t), 3) or not np.isfinite(values).all():
        raise ValueError('callable drive must return finite A with shape (len(t), 3)')
    return values


def _readonly(a):
    a = np.array(a, dtype=float)
    a.setflags(write=False)
    return a


@dataclass(frozen=True, slots=True)
class FieldStateConfig:
    """Resolved settings of one field state (immutable record, not a checkpoint)."""
    model: str
    dt: float
    steps: int
    points: int
    parameters: tuple
    units: str = field(default='atomic units; Ac=A/c; laboratory frame', init=False)


class _FieldState:
    """Common step bookkeeping. Subclasses implement _advance(J) -> (A_new, commit)."""
    __slots__ = ('_config', '_n', '_A')
    _HISTORY = ''
    _INTEGRATOR = ''

    @property
    def timing(self):
        """salmonpy.time/1 declaration: j at n*dt (lag 1, t=0 unused) -> A at (n+1)*dt."""
        return TimeSpec(self._config.dt, (InputPort('j', lag=1, uses_t0=False),), (OutputPort('A'),),
                        history=self._HISTORY, integrator=self._INTEGRATOR)

    @property
    def config(self): return self._config
    @property
    def dt(self): return self._config.dt
    @property
    def points(self): return self._config.points
    @property
    def time(self): return self._n*self._config.dt
    @property
    def step_index(self): return self._n
    @property
    def A(self):
        """Ac at the current time, shape (points, 3)."""
        return _readonly(self._A)

    def step(self, J):
        """Consume charge current j at the current time; return Ac at the next time, (points, 3).

        Validation or numerical failure leaves the state unchanged.
        """
        J = _array(J, 'J', (self.points, 3))
        if self._n >= self._config.steps:
            raise ValueError('field state capacity (steps) is exhausted')
        A_new, commit = self._advance(J)
        if not np.isfinite(A_new).all():
            raise FloatingPointError('field update became nonfinite; state unchanged')
        commit()
        self._A = A_new
        self._n += 1
        return _readonly(A_new)


# ---------------------------------------------------------------------------
# Uniform (single macroscopic point) fields: SALMON trans_longi = tr / lo / 2d
# ---------------------------------------------------------------------------

def _uniform_external(external, dt, steps):
    return _evaluate(external, np.arange(steps+2)*dt, dt=dt)  # Ac_ext(0:nt+1)


@dataclass(frozen=True, slots=True)
class Transverse:
    """Prescribed total field, SALMON trans_longi='tr': A_tot(n) = A_ext(n)."""
    external: object = None

    def __post_init__(self):
        object.__setattr__(self, 'external', _drive(self.external, 'external'))

    def new_state(self, *, dt, steps):
        dt = _positive(dt, 'dt'); steps = _steps(steps)
        return _TransverseState(self, dt, steps)


class _TransverseState(_FieldState):
    __slots__ = ('_ext',)
    _INTEGRATOR = "SALMON trans_longi='tr': A(n) = A_ext(n)"
    def __init__(self, model, dt, steps):
        self._config = FieldStateConfig('transverse', dt, steps, 1,
                                        (('external', _describe(model.external)),))
        self._ext = _uniform_external(model.external, dt, steps)
        self._n = 0; self._A = self._ext[0][None].copy()

    def _advance(self, J):
        return self._ext[self._n+1][None].copy(), (lambda: None)


@dataclass(frozen=True, slots=True)
class Longitudinal:
    """Bulk longitudinal geometry, SALMON trans_longi='lo'.

    Ac_ind(n+1) = 2 Ac_ind(n) - Ac_ind(n-1) + 4 pi j(n) dt**2  (SALMON: -4 pi curr dt**2),
    A_tot(n) = A_ext(n) + Ac_ind(n), with Ac_ind(0) = Ac_ind(1) = 0.
    Explicit leapfrog: the diamagnetic part of J limits stability to roughly
    omega_p * dt < 2 (omega_p**2 = 4 pi n_e).
    """
    external: object = None

    def __post_init__(self):
        object.__setattr__(self, 'external', _drive(self.external, 'external'))

    def new_state(self, *, dt, steps):
        dt = _positive(dt, 'dt'); steps = _steps(steps)
        return _LongitudinalState(self, dt, steps)


class _LongitudinalState(_FieldState):
    __slots__ = ('_ext', '_ind')
    _HISTORY = 'Ac_ind(0) = Ac_ind(1) = 0'
    _INTEGRATOR = "SALMON trans_longi='lo' leapfrog"
    def __init__(self, model, dt, steps):
        self._config = FieldStateConfig('longitudinal', dt, steps, 1,
                                        (('external', _describe(model.external)),))
        self._ext = _uniform_external(model.external, dt, steps)
        self._ind = np.zeros((2, 3))  # Ac_ind(n-1), Ac_ind(n)
        self._n = 0; self._A = self._ext[0][None].copy()

    @property
    def induced(self):
        """Ac_ind at the current time."""
        return _readonly(self._ind[1])

    def _advance(self, J):
        dt = self.dt
        if self._n == 0:
            new = np.zeros(3)
        else:
            new = 2.0*self._ind[1] - self._ind[0] + 4.0*_PI*J[0]*dt**2
        ind = np.stack((self._ind[1], new))
        def commit(): self._ind = ind
        return (self._ext[self._n+1] + new)[None], commit


@dataclass(frozen=True, slots=True)
class ThinFilm:
    """Thin film between two dielectrics at normal incidence, SALMON trans_longi='2d'.

    epsilon = (incident side, substrate side); thickness in bohr. The returned
    A is SALMON's Ac_tot = A_ext + A_ref, the field transmitted into the film.
    dA_ref/dt = r dA_ext/dt + 4 pi d j / (c (n1+n2)),  r = (n1-n2)/(n1+n2),
    integrated with SALMON's 100 equal substeps per dt and A_ref(0)=A_ref(1)=0.
    """
    thickness: float
    external: object = None
    epsilon: tuple = (1.0, 1.0)

    def __post_init__(self):
        object.__setattr__(self, 'thickness', _positive(self.thickness, 'thickness'))
        object.__setattr__(self, 'external', _drive(self.external, 'external'))
        eps = _array(self.epsilon, 'epsilon', (2,))
        if np.any(eps <= 0): raise ValueError('epsilon must be positive')
        object.__setattr__(self, 'epsilon', tuple(map(float, eps)))

    def new_state(self, *, dt, steps):
        dt = _positive(dt, 'dt'); steps = _steps(steps)
        return _ThinFilmState(self, dt, steps)


class _ThinFilmState(_FieldState):
    __slots__ = ('_ext', '_ref', '_d', '_n1', '_n2')
    _SUBSTEPS = 100  # SALMON em_field.f90: integer, parameter :: m=100
    _HISTORY = 'A_ref(0) = A_ref(1) = 0'
    _INTEGRATOR = "SALMON trans_longi='2d', 100 equal substeps"

    def __init__(self, model, dt, steps):
        self._config = FieldStateConfig('thin_film', dt, steps, 1, (
            ('external', _describe(model.external)), ('thickness', model.thickness),
            ('epsilon', model.epsilon), ('substeps', self._SUBSTEPS)))
        self._ext = _uniform_external(model.external, dt, steps)
        self._d = model.thickness
        self._n1 = np.sqrt(model.epsilon[0]); self._n2 = np.sqrt(model.epsilon[1])
        self._ref = np.zeros(3)
        self._n = 0; self._A = self._ext[0][None].copy()

    @property
    def reflected(self):
        """Reflected Ac (SALMON Ac_ind) at the current time."""
        return _readonly(self._ref)

    @property
    def transmitted(self):
        """Transmitted Ac (SALMON Ac_tot) at the current time."""
        return _readonly(self._A[0])

    def _advance(self, J):
        n = self._n; dt = self.dt
        if n == 0:
            ref = np.zeros(3)
        else:
            n1, n2, m = self._n1, self._n2, self._SUBSTEPS
            rate = ((n1-n2)/(n1+n2))*(self._ext[n+1]-self._ext[n])/dt \
                + (4.0*_PI/(SPEED_OF_LIGHT*(n1+n2)))*self._d*J[0]
            ref = self._ref.copy()
            for _ in range(m):  # same accumulation order as SALMON
                ref = ref + (dt/m)*rate
        def commit(): self._ref = ref
        return (self._ext[n+1] + ref)[None], commit


# ---------------------------------------------------------------------------
# 1D Weyl-gauge propagation with macroscopic points (SALMON multiscale, fdtddim='1d')
# ---------------------------------------------------------------------------

_BOUNDARIES = ('abc', 'pec', 'periodic')


@dataclass(frozen=True, slots=True)
class Multiscale1D:
    """1D Maxwell propagation along x with one matter point per cell 1..cells.

    Grid cells i = 1-vacuum[0] .. cells+vacuum[1] at x = i*spacing (bohr), one
    ghost cell each side. Matter currents enter at cells 1..cells. The incident
    wave A_inc(t - x/c) enters through the left boundary, as SALMON's incident():
    cells i <= 0 are prefilled at times -3dt, -2dt, -dt, one update produces the
    t=0 field, and the boundary uses A_inc at x=(1-vacuum[0]-0.5)*spacing.
    Boundaries are 'abc' (1st-order Mur), 'pec' or 'periodic', per side.
    The incident wave must be transverse (no x component). Stability requires
    c*dt <= spacing; c*dt == spacing is dispersion-free in vacuum.
    """
    cells: int
    spacing: float
    vacuum: tuple = (0, 0)
    incident: object = None
    boundaries: tuple = ('abc', 'abc')

    def __post_init__(self):
        if isinstance(self.cells, (bool, np.bool_)) or not isinstance(self.cells, (int, np.integer)) or self.cells < 1:
            raise ValueError('cells must be a positive integer')
        object.__setattr__(self, 'cells', int(self.cells))
        object.__setattr__(self, 'spacing', _positive(self.spacing, 'spacing'))
        vac = tuple(self.vacuum)
        if len(vac) != 2 or any(isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer)) or v < 0 for v in vac):
            raise ValueError('vacuum must be two nonnegative integers (left, right)')
        object.__setattr__(self, 'vacuum', tuple(map(int, vac)))
        inc = _drive(self.incident, 'incident')
        if isinstance(inc, FieldSamples):
            raise TypeError('incident must be evaluable at arbitrary times: AcosPulse(s) or a callable')
        object.__setattr__(self, 'incident', inc)
        bc = tuple(self.boundaries)
        if len(bc) != 2 or any(b not in _BOUNDARIES for b in bc):
            raise ValueError(f'boundaries must be two of {_BOUNDARIES}')
        if (bc[0] == 'periodic') != (bc[1] == 'periodic'):
            raise ValueError('periodic boundaries must be used on both sides')
        object.__setattr__(self, 'boundaries', bc)

    @property
    def x(self):
        """Cell coordinates (bohr) of the interior grid, cells 1-vacuum[0]..cells+vacuum[1]."""
        i = np.arange(1-self.vacuum[0], self.cells+self.vacuum[1]+1)
        return i*self.spacing

    def new_state(self, *, dt, steps):
        dt = _positive(dt, 'dt'); steps = _steps(steps)
        if SPEED_OF_LIGHT*dt > self.spacing*(1+1e-12):
            raise ValueError('unstable: c*dt must not exceed spacing')
        return _Multiscale1DState(self, dt, steps)


class _Multiscale1DState(_FieldState):
    __slots__ = ('_model', '_is', '_ie', '_old', '_cur', '_new', '_j', '_inc')
    _HISTORY = 'SALMON incident(): cells <= 0 prefilled at -3dt..-dt, one update to t=0; j_em = 0'
    _INTEGRATOR = "SALMON 1D Weyl gauge (fdtd_weyl_gauge.f90), centered E at (n-1)*dt"

    def __init__(self, model, dt, steps):
        self._config = FieldStateConfig('multiscale_1d', dt, steps, model.cells, (
            ('cells', model.cells), ('spacing', model.spacing), ('vacuum', model.vacuum),
            ('incident', _describe(model.incident)), ('boundaries', model.boundaries)))
        self._model = model
        h = model.spacing; c = SPEED_OF_LIGHT
        self._is = 1-model.vacuum[0]; self._ie = model.cells+model.vacuum[1]
        size = self._ie-self._is+3  # interior plus one ghost cell each side
        # Incident samples at the left boundary, n = -1 .. steps+2 (SALMON Ac_inc(:, -1:nt+2)).
        n = np.arange(-1, steps+3)
        self._inc = _evaluate(model.incident, n*dt - (self._is-0.5)*h/c)
        if np.any(self._inc[:, 0] != 0):
            raise ValueError('incident wave must be transverse: x is the propagation axis')
        # Prefill i = is-1 .. 0 with the incoming wave at t = -dt, -2dt, -3dt.
        self._new = np.zeros((size, 3)); self._cur = np.zeros((size, 3)); self._old = np.zeros((size, 3))
        left = np.arange(self._is-1, 1)
        k = left - (self._is-1)
        for array, t0 in ((self._new, -dt), (self._cur, -2*dt), (self._old, -3*dt)):
            array[k] = _evaluate(model.incident, t0 - left*h/c)
        self._j = np.zeros((size, 3))
        # One update to the t=0 profile (SALMON: itt=0 weyl_calc with j_em=0).
        self._old, self._cur, self._new = self._update(self._cur, self._new, self._j, 0)
        self._n = 0
        self._A = self._macro(self._new)

    def _macro(self, a):
        start = 1-(self._is-1)
        return a[start:start+self._model.cells].copy()

    def _inc_at(self, itt):
        return self._inc[itt+1]  # array starts at n=-1

    def _update(self, cur_prev, new_prev, j, itt):
        """One SALMON dt_evolve_Ac_1d call; returns (old, cur, new) after the update."""
        old = cur_prev; cur = new_prev
        dt = self.dt; h = self._model.spacing; c = SPEED_OF_LIGHT
        r = 1.0/h
        tmp = np.empty_like(cur)
        rot2 = np.zeros((len(cur)-2, 3))
        rot2[:, 1:3] = -((cur[2:, 1:3] - 2.0*cur[1:-1, 1:3]) + cur[:-2, 1:3]) * r**2
        tmp[1:-1] = (2.0*cur[1:-1] - old[1:-1]) + j[1:-1]*4.0*_PI*(dt**2)/1.0 \
            - rot2*(c*dt)**2/1.0
        k = (c*dt-h)/(c*dt+h)
        left, right = self._model.boundaries
        if left == 'periodic':
            tmp[0] = tmp[-2]; tmp[-1] = tmp[1]
        else:
            if left == 'pec':
                tmp[0] = cur[0]
            else:
                tmp[0] = cur[1] + k*tmp[1] - k*cur[0] \
                    + (4.0*h)/(c*dt+h)*(self._inc_at(itt)-self._inc_at(itt-1))
            if right == 'pec':
                tmp[-1] = cur[-1]
            else:
                tmp[-1] = cur[-2] + k*tmp[-2] - k*cur[-1]
        return old, cur, tmp

    @property
    def x(self):
        return self._model.x

    @property
    def grid_A(self):
        """Interior Ac on the whole 1D grid at the current time."""
        return _readonly(self._new[1:-1])

    @property
    def grid_E(self):
        """SALMON's centered E = -(A(n) - A(n-2))/(2 dt), located at time (n-1)*dt."""
        return _readonly(-(self._new[1:-1]-self._old[1:-1])*0.5/self.dt)

    def _advance(self, J):
        itt = self._n+1
        j = np.zeros_like(self._j)
        if self._n > 0:
            start = 1-(self._is-1)
            j[start:start+self._model.cells] = J  # SALMON j_em = -curr: charge current
        old, cur, new = self._update(self._cur, self._new, j, itt)
        def commit():
            self._old, self._cur, self._new, self._j = old, cur, new, j
        return self._macro(new), commit


# ---------------------------------------------------------------------------
# Coupling driver
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class CoupledResponse:
    """Coupled history. time (steps+1,), A and current (steps+1, points, 3)."""
    time: np.ndarray
    A: np.ndarray
    current: np.ndarray
    field_config: FieldStateConfig
    backends: tuple
    timing: dict
    convention: str = field(default='A and current at n*dt; current is the electric (charge) current '
                                     'density j (= -SALMON matter current); j(0) not fed back', init=False)
    output_directory: str | None = None


def couple(matter, field_model, *, dt, steps, output_directory=None):
    """Alternate matter A->J and field J->A for ``steps`` steps.

    ``matter`` is one AInitialState used for every point (each point gets its own
    independent state) or a sequence with one AInitialState per field point.
    Every matter state starts at the field's A(0). Before stepping, each matter
    state's ``timing`` declaration (or the crystal default, with a warning) is
    checked against the field's with ``salmonpy.timing.check_loop``. Built-in SBE
    and TDDFT states output the charge current j; a declared J_matter output is
    converted j = -J_matter, the only adapter.
    Matter states are closed on return, also on failure.
    """
    fstate = field_model.new_state(dt=dt, steps=steps)
    dt = fstate.dt; steps = fstate.config.steps; P = fstate.points
    if isinstance(matter, (list, tuple)):
        if len(matter) != P: raise ValueError(f'expected {P} matter initial states, one per field point')
        initials = tuple(matter)
    else:
        initials = (matter,)*P
    backends = tuple(str(getattr(m, 'backend', type(m).__name__)) for m in initials)
    output_path = None
    if output_directory is not None:
        from pathlib import Path
        output_path = str(Path(output_directory).expanduser().resolve())
        if Path(output_path).exists(): raise FileExistsError(output_path)
    states = []; signs = []; matter_records = []
    field_timing = fstate.timing
    try:
        A = np.empty((steps+1, P, 3)); J = np.empty((steps+1, P, 3))
        A[0] = fstate.A
        for p, initial in enumerate(initials):
            state = initial.new_state(dt=dt, time0=0., initial_A=A[0, p].copy())
            if any(state is s for s in states): raise ValueError('matter backend reused one state for two points')
            states.append(state)
            spec = matter_timing(state, dt)
            sign = charge_sign(spec)
            check_loop(field_timing, spec)
            # The only conversion: a backend declaring SALMON matter current J -> electric j.
            signs.append(sign)
            matter_records.append(spec.to_dict())
            J[0, p] = signs[p]*_array(state.current, 'initial current', (3,))
        for n in range(steps):
            A[n+1] = fstate.step(J[n])
            t = (n+1)*dt
            for p, state in enumerate(states):
                J[n+1, p] = signs[p]*_array(state.step(A[n+1, p].copy()), 'current', (3,))
                if not np.isfinite(state.time) or not np.isclose(state.time, t, rtol=0, atol=max(1e-8*dt, 16*np.finfo(float).eps*(n+1)*max(t, dt))):
                    raise ValueError('matter state time does not match the field clock')
        time = np.arange(steps+1)*dt
        timing = {'field': field_timing.to_dict(), 'matter': matter_records}
        if output_path is not None:
            from .outputs import save_arrays
            save_arrays(output_path, {'time': time, 'A': A, 'current': J}, metadata={
                'field': asdict(fstate.config), 'backends': list(backends), 'dt': dt,
                'units': {'time': 'a.u.', 'A': 'Ac=A/c, a.u.', 'current': 'electric (charge) current density, a.u.'},
                'time_location': 'A and current at n*dt', 'matter_current': '-current',
                'timing': timing, 'restart_checkpoint': False})
    finally:
        for state in states:
            close = getattr(state, 'close', None)
            if close is not None: close()
    return CoupledResponse(_readonly(time), _readonly(A), _readonly(J), fstate.config, backends,
                           timing, output_directory=output_path)
