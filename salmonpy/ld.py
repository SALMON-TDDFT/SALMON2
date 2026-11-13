"""Isotropic passive Lorentz–Drude model; angular frequencies and atomic units."""
import ctypes
from dataclasses import dataclass, field
import hashlib
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np


def _array(value, name, shape=None):
    if np.iscomplexobj(value):
        raise ValueError(f'{name} must be real-valued')
    a = np.array(value, dtype=np.float64, order='C', copy=True)
    if not np.isfinite(a).all() or (shape is not None and a.shape != shape):
        raise ValueError(f'{name}: expected finite values with shape {shape}')
    return a


def _positive(value, name, zero=False):
    a = np.asarray(value)
    if a.ndim != 0 or a.dtype.kind not in 'iuf':
        raise ValueError(f'{name} must be a real numeric scalar')
    value = float(value)
    if not np.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f'{name} must be finite and {"nonnegative" if zero else "positive"}')
    return value


@dataclass(frozen=True, slots=True)
class LDModelConfig:
    """Resolved immutable material definition; obtained through material.config."""
    plasma_frequency: float
    strength: tuple[float, ...]
    resonance: tuple[float, ...]
    damping: tuple[float, ...]
    units: str = field(default='atomic_gaussian', init=False)
    background_permittivity: float = field(default=1.0, init=False)


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

    epsilon(omega) = 1 + sum(f * wp**2 / (w0**2 - omega**2 - i*gamma*omega)).
    Time convention: exp(-i*omega*t). Instantaneous background epsilon is 1.
    """
    __slots__ = ('_config', '_wp', '_f', '_w', '_g')

    def __init__(self, *, plasma_frequency, strength, resonance, damping):
        self._wp = _positive(plasma_frequency, 'plasma_frequency', zero=True)
        self._f = _array(strength, 'strength')
        if self._f.ndim != 1 or not self._f.size:
            raise ValueError('strength must be a nonempty one-dimensional array')
        self._w = _array(resonance, 'resonance', self._f.shape)
        self._g = _array(damping, 'damping', self._f.shape)
        if any(np.any(a < 0) for a in (self._f, self._w, self._g)):
            raise ValueError('strength, resonance and damping must be nonnegative')
        self._config = LDModelConfig(self._wp, tuple(map(float, self._f)),
                                     tuple(map(float, self._w)), tuple(map(float, self._g)))
        # Immutable byte-backed arrays cannot be made writable by NumPy callers.
        self._f, self._w, self._g = (np.frombuffer(a.tobytes(), dtype=np.float64)
                                    for a in (self._f, self._w, self._g))

    @property
    def config(self):
        return self._config

    def permittivity(self, omega):
        """Continuous-model dielectric function. Raise at a real-axis pole."""
        omega = _array(omega, 'omega')
        result = np.ones(omega.shape, dtype=complex)
        for f, w, g in zip(self._f, self._w, self._g):
            if f == 0 or self._wp == 0:
                continue
            denominator = w*w - omega*omega - 1j*g*omega
            if np.any(denominator == 0):
                raise ValueError('permittivity is singular at a requested frequency')
            result += f*self._wp**2/denominator
        return result

    def new_state(self, *, dt, initial_field=(0., 0., 0.)):
        """Independent initially unpolarized state at t=0, with E(0) specified."""
        return LDState(self, dt, initial_field)

    def respond(self, field, *, dt):
        """E samples of shape (N,3) at t=0,dt,...; output includes initial zeros."""
        field = _array(field, 'field')
        if field.ndim != 2 or field.shape[1] != 3 or len(field) == 0:
            raise ValueError('field must have shape (N,3), N >= 1')
        state = self.new_state(dt=dt, initial_field=field[0])
        p, j = state._advance(field)
        return Response(np.arange(len(field))*state.dt, p, j)


class LDState:
    """Independent explicit per-pole history. A state is not thread-shareable."""
    __slots__ = ('_config', '_material', '_field', '_p', '_j', '_steps')

    def __init__(self, material, dt, initial_field):
        dt = _positive(dt, 'dt')
        self._material = material
        self._field = _array(initial_field, 'initial_field', (3,))
        self._config = LDStateConfig(material.config, dt, tuple(map(float, self._field)))
        self._p = np.zeros((len(material._f), 3))
        self._j = np.zeros_like(self._p)
        self._steps = 0

    @property
    def dt(self): return self._config.dt

    @property
    def config(self): return self._config

    @property
    def time(self): return self._steps*self.dt

    @property
    def current(self): return self._j.sum(axis=0)

    @property
    def polarization(self): return self._p.sum(axis=0)

    def step(self, field):
        """Advance dt using E at the new time; return (P, J), not FDTD Q."""
        field = _array(field, 'field', (3,))
        p, j = self._advance(np.stack((self._field, field)))
        return p[-1].copy(), j[-1].copy()

    def _advance(self, field):
        m = self._material
        p = np.empty_like(field); j = np.empty_like(field)
        # Work on copies so a rejected nonfinite output does not corrupt the state.
        pp = self._p.copy(); jj = self._j.copy()
        native = _native()
        ptr = ctypes.POINTER(ctypes.c_double)
        args = [a.ctypes.data_as(ptr) for a in (m._f,m._w,m._g,field,pp,jj,p,j)]
        native(len(m._f), len(field), self.dt, m._wp, *args)
        if not all(np.isfinite(a).all() for a in (pp,jj,p,j)):
            raise FloatingPointError('LD response became nonfinite; state unchanged')
        self._p, self._j = pp, jj
        self._field = field[-1].copy()
        self._steps += len(field)-1
        return p, j


@lru_cache(maxsize=1)
def _native():
    suffix = {'darwin': '.dylib', 'win32': '.dll'}.get(sys.platform, '.so')
    path = Path(__file__).parent / '_native' / ('libsalmonpy_ld'+suffix)
    if not path.exists():
        raise RuntimeError('Build the LD library first: python3 -m salmonpy.build_native')
    record = path.with_suffix(path.suffix+'.sha256')
    if not record.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != record.read_text().strip():
        raise RuntimeError('LD library checksum differs; rebuild with salmonpy.build_native')
    lib = ctypes.CDLL(str(path))
    fn = lib.salmonpy_ld
    fn.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double] + [ctypes.POINTER(ctypes.c_double)]*8
    fn.restype = None
    return fn
