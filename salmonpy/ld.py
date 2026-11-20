"""Isotropic passive Lorentz–Drude model; angular frequencies and atomic units."""
from .errors import InputError, NumericalError, RuntimeFailure, SnapshotError
import ctypes
from .snapshots import _SnapshotMixin, _fingerprint
from dataclasses import dataclass, field
import hashlib
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np


def _array(value, name, shape=None):
    if np.iscomplexobj(value):
        raise InputError(f'{name} must be real-valued')
    a = np.array(value, dtype=np.float64, order='C', copy=True)
    if not np.isfinite(a).all() or (shape is not None and a.shape != shape):
        raise InputError(f'{name}: expected finite values with shape {shape}')
    return a


def _positive(value, name, zero=False):
    a = np.asarray(value)
    if a.ndim != 0 or a.dtype.kind not in 'iuf':
        raise InputError(f'{name} must be a real numeric scalar')
    value = float(value)
    if not np.isfinite(value) or (value < 0 if zero else value <= 0):
        raise InputError(f'{name} must be finite and {"nonnegative" if zero else "positive"}')
    return value


@dataclass(frozen=True, slots=True)
class LDModelConfig:
    """Resolved immutable material definition; obtained through material.config."""
    plasma_frequency: float
    strength: tuple[float, ...]
    resonance: tuple[float, ...]
    damping: tuple[float, ...]
    units: str = field(default='atomic_gaussian', init=False)
    epsilon_infinity: float = 1.0

    @property
    def background_permittivity(self): return self.epsilon_infinity


@dataclass(frozen=True, slots=True)
class LDStateConfig:
    """Resolved initialization settings, not the evolving P/J or current field."""
    model: LDModelConfig
    dt: float
    initial_field: tuple[float, ...]
    initial_time: float = field(default=0.0, init=False)
    initial_polarization: str = field(default='zero_per_pole', init=False)
    initial_current: str = field(default='zero_per_pole', init=False)
    method: str = field(default='salmon_ld_trapezoidal', init=False)


@dataclass(frozen=True)
class Response:
    time: np.ndarray
    polarization: np.ndarray
    current: np.ndarray


class LorentzDrude:
    """Passive isotropic poles. Frequencies are angular frequencies in atomic units.

    epsilon(omega) = epsilon_infinity + sum(f * wp**2 / (w0**2 - omega**2 - i*gamma*omega)).
    Time convention: exp(-i*omega*t). P/J describe dispersive poles; D=epsilon_infinity*E+4*pi*P.
    """
    __slots__ = ('_config', '_wp', '_f', '_w', '_g')

    def __init__(self, *, plasma_frequency, strength, resonance, damping, epsilon_infinity=1.0):
        self._wp = _positive(plasma_frequency, 'plasma_frequency', zero=True)
        self._f = _array(strength, 'strength')
        if self._f.ndim != 1 or not self._f.size:
            raise InputError('strength must be a nonempty one-dimensional array')
        self._w = _array(resonance, 'resonance', self._f.shape)
        self._g = _array(damping, 'damping', self._f.shape)
        if any(np.any(a < 0) for a in (self._f, self._w, self._g)):
            raise InputError('strength, resonance and damping must be nonnegative')
        self._config = LDModelConfig(self._wp, tuple(map(float, self._f)),
                                     tuple(map(float, self._w)), tuple(map(float, self._g)),
                                     epsilon_infinity=_positive(epsilon_infinity,"epsilon_infinity"))
        # Immutable byte-backed arrays cannot be made writable by NumPy callers.
        self._f, self._w, self._g = (np.frombuffer(a.tobytes(), dtype=np.float64)
                                    for a in (self._f, self._w, self._g))

    @property
    def config(self):
        return self._config

    def permittivity(self, omega):
        """Continuous-model dielectric function. Raise at a real-axis pole."""
        omega = _array(omega, 'omega')
        import ctypes as c
        from .errors import RuntimeFailure
        result=np.empty(omega.shape,dtype=np.complex128)
        try: fn=_handles_api().library.sp_ld_permittivity
        except AttributeError as exc:
            raise RuntimeFailure('native LD observables ABI missing; explicitly rebuild kernels') from exc
        fn.argtypes=[c.c_int,c.c_int,c.c_double,c.c_double]+[c.c_void_p]*5
        fn.restype=c.c_int
        status=fn(len(self._f),omega.size,self._wp,self.config.epsilon_infinity,
                  *(a.ctypes.data for a in (self._f,self._w,self._g,omega,result)))
        if status: raise InputError('permittivity is singular or nonfinite at a requested frequency')
        return result

    def _background(self,field,polarization,magnetic_field):
        import ctypes as c
        from .errors import RuntimeFailure
        e=_array(field,'field');p=_array(polarization,'polarization',e.shape)
        h=_array(magnetic_field,'magnetic_field',e.shape)
        if e.ndim<1 or e.shape[-1]!=3: raise InputError('field arrays must have shape (...,3)')
        d=np.empty_like(e);u=np.empty(e.shape[:-1])
        try: fn=_handles_api().library.sp_ld_background
        except AttributeError as exc:
            raise RuntimeFailure('native LD background ABI missing; explicitly rebuild kernels') from exc
        fn.argtypes=[c.c_int,c.c_double]+[c.c_void_p]*5;fn.restype=c.c_int
        status=fn(e.size//3,self.config.epsilon_infinity,*(a.ctypes.data for a in (e,p,h,d,u)))
        if status: raise NumericalError('nonfinite LD background observable')
        return d,u

    def displacement(self,field,polarization):
        """D=epsilon_infinity*E+4*pi*P_LD, with simultaneous (...,3) arrays."""
        return self._background(field,polarization,np.zeros_like(field))[0]

    def background_energy(self,electric_field,magnetic_field):
        """(epsilon_infinity*E^2+H^2)/(8*pi), mu=1; excludes LD pole energy."""
        return self._background(electric_field,np.zeros_like(electric_field),magnetic_field)[1]

    def new_state(self, *, dt, initial_field=(0., 0., 0.)):
        """Independent initially unpolarized state at t=0, with E(0) specified."""
        return LDState(self, dt, initial_field)

    def respond(self, field, *, dt):
        """E samples of shape (N,3) at t=0,dt,...; output includes initial zeros."""
        field = _array(field, 'field')
        if field.ndim != 2 or field.shape[1] != 3 or len(field) == 0:
            raise InputError('field must have shape (N,3), N >= 1')
        state = self.new_state(dt=dt, initial_field=field[0])
        try:
            p, j = state._advance(field)
            return Response(np.arange(len(field))*state.dt, p, j)
        finally:
            state.close()


class LDState(_SnapshotMixin):
    """Independent native per-pole history. A state is not thread-shareable."""
    __slots__ = ('_config', '_material', '_field', '_polarization', '_current',
                 '_steps', '_native_state', '_closed')
    _snapshot_kind = 'ld'
    _snapshot_counter = '_steps'

    def _snapshot_owner(self): return self._material
    def _snapshot_settings(self): return _fingerprint(self.config)
    def _snapshot_arrays(self):
        p,j=self._native_state.export()
        return {'_field':self._field,'_p':p,'_j':j}

    def __init__(self, material, dt, initial_field):
        dt = _positive(dt, 'dt')
        self._material = material
        self._field = _array(initial_field, 'initial_field', (3,))
        self._config = LDStateConfig(material.config, dt, tuple(map(float, self._field)))
        self._polarization=np.zeros(3);self._current=np.zeros(3)
        self._steps = 0;self._closed=False
        from ._ld_handles import NativeState
        self._native_state=NativeState(_handles_api(),material,dt,self._field)

    @property
    def dt(self): return self._config.dt
    @property
    def config(self): return self._config
    @property
    def time(self): return self._steps*self.dt
    @property
    def current(self): return self._current.copy()
    @property
    def polarization(self): return self._polarization.copy()

    def _require_open(self):
        if self._closed:raise RuntimeFailure('LD state is closed')

    def step(self, field):
        """Advance dt using E at the new time; return (P, J), not FDTD Q."""
        self._require_open()
        field = _array(field, 'field', (3,))
        p=np.empty((2,3));j=np.empty_like(p)
        # Allocate accepted metadata and return values before committing native history.
        accepted_p=np.empty(3);accepted_j=np.empty(3);steps=self._steps+1
        result=(p[-1],j[-1]);incoming=field.reshape(1,3)
        self._native_state.advance(incoming,p,j)
        np.copyto(accepted_p,result[0]);np.copyto(accepted_j,result[1])
        self._polarization=accepted_p;self._current=accepted_j
        self._field=field;self._steps=steps
        return result

    def _advance(self, field):
        self._require_open()
        if np.iscomplexobj(field):raise InputError('field must be real-valued')
        field=np.ascontiguousarray(field,dtype=np.float64)
        if field.ndim!=2 or field.shape[1]!=3 or not len(field):
            raise InputError('field must have shape (N,3), N >= 1')
        if not np.array_equal(field[0],self._field):
            raise InputError('batch must begin at the accepted field')
        if len(field)>2**31-1:raise InputError('too many LD samples')
        p=np.empty_like(field);j=np.empty_like(field)
        accepted_p=np.empty(3);accepted_j=np.empty(3);last=field[-1].copy();steps=self._steps+len(field)-1
        last_p=p[-1];last_j=j[-1];result=(p,j)
        self._native_state.advance(field[1:],p,j)
        np.copyto(accepted_p,last_p);np.copyto(accepted_j,last_j)
        # Separate totals prevent returned response arrays from changing accepted history.
        self._polarization=accepted_p;self._current=accepted_j
        self._field=last;self._steps=steps
        return result

    def snapshot(self):
        if self._closed:raise SnapshotError('cannot snapshot a closed LD state')
        return super().snapshot()

    def restore(self,snapshot):
        if self._closed:raise SnapshotError('cannot restore a closed LD state')
        from types import SimpleNamespace
        proxy=SimpleNamespace(dt=self.dt,_steps=self._steps,
            _snapshot_kind=self._snapshot_kind,_snapshot_counter=self._snapshot_counter,
            _snapshot_owner=self._snapshot_owner,_snapshot_settings=self._snapshot_settings,
            _snapshot_arrays=lambda:{'_field':self._field,'_p':np.empty((self._native_state.n,3)),
                                     '_j':np.empty((self._native_state.n,3))})
        _SnapshotMixin.restore(proxy,snapshot)
        try:p,j=self._native_state.replace(proxy._field,proxy._p,proxy._j)
        except (InputError,NumericalError) as exc:raise SnapshotError('invalid LD snapshot history') from exc
        self._field=proxy._field;self._steps=proxy._steps
        self._polarization=p;self._current=j

    def close(self):
        if not self._closed:
            self._native_state.close();self._closed=True

    def __enter__(self):self._require_open();return self
    def __exit__(self,*args):self.close()


def _handles_api():
    from ._ld_handles import default_api
    return default_api()


@lru_cache(maxsize=1)
def _native():
    suffix = {'darwin': '.dylib', 'win32': '.dll'}.get(sys.platform, '.so')
    path = Path(__file__).parent / '_native' / ('libsalmonpy_ld'+suffix)
    if not path.exists():
        raise RuntimeFailure('Build the LD library first: python3 -m salmonpy.build_native')
    record = path.with_suffix(path.suffix+'.sha256')
    if not record.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != record.read_text().strip():
        raise RuntimeFailure('LD library checksum differs; rebuild with salmonpy.build_native')
    lib = ctypes.CDLL(str(path))
    fn = lib.salmonpy_ld
    fn.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double] + [ctypes.POINTER(ctypes.c_double)]*8
    fn.restype = None
    return fn
