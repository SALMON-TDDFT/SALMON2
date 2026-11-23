"""Isotropic passive Lorentz–Drude model; angular frequencies and atomic units."""
from .errors import InputError, NumericalError, RuntimeFailure, SnapshotError
import ctypes
from .snapshots import _SnapshotMixin, _fingerprint
from dataclasses import dataclass, field, asdict
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


def _execution_mode_name(value):
    if not isinstance(value,str) or value not in ('checked','performance'):
        raise InputError('execution_mode must be checked or performance')
    return value


def _array_unchecked(value,name,shape=None):
    """Owned real ABI buffer; intentionally no value/finite scan."""
    if np.iscomplexobj(value):raise InputError(f'{name} must be real-valued')
    a=np.array(value,dtype=np.float64,order='C',copy=True)
    if shape is not None and a.shape!=shape:raise InputError(f'{name}: expected shape {shape}')
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
    valid_frequency_range: tuple | None = None

    @property
    def background_permittivity(self): return self.epsilon_infinity


    @property
    def tensor_config(self):
        return LDTensorModelConfig(self)


@dataclass(frozen=True,slots=True)
class LDTensorModelConfig:
    """Compressed isotropic component or reciprocal tensor-pole definition.

    General strengths have shape(poles,3,3). Axes are lab xyz; all numerical
    matrix/positivity checks are native. Independent axes are never averaged.
    """
    component: LDModelConfig | None = None
    plasma_frequency: float = 1.
    strength: tuple = ()
    resonance: tuple = ()
    damping: tuple = ()
    epsilon_infinity: tuple = ((1.,0.,0.),(0.,1.,0.),(0.,0.,1.))
    valid_frequency_range: tuple | None = None
    units: str = field(default='atomic_gaussian',init=False)

    def __post_init__(self):
        if self.component is not None:
            if not isinstance(self.component,LDModelConfig):raise InputError('LDModelConfig component required')
            if len(self.strength) or len(self.resonance) or len(self.damping) or self.plasma_frequency!=1.:
                raise InputError('do not mix isotropic component and tensor parameters')
            if not np.array_equal(self.epsilon_infinity,((1.,0.,0.),(0.,1.,0.),(0.,0.,1.))):
                raise InputError('background belongs to the isotropic component')
            if self.valid_frequency_range is not None:raise InputError('band belongs to the isotropic component')
            component=LorentzDrude.from_config(self.component).config
            object.__setattr__(self,'component',component)
            return
        wp=_positive(self.plasma_frequency,'plasma_frequency',zero=True)
        f=_array(self.strength,'strength');w=_array(self.resonance,'resonance');g=_array(self.damping,'damping')
        if f.ndim!=3 or f.shape[1:]!=(3,3) or not 1<=len(f)<=(2**31-1)//9:
            raise InputError('tensor strength must have shape(poles,3,3)')
        if w.shape!=(len(f),) or g.shape!=w.shape:raise InputError('resonance/damping must match tensor pole count')
        eps=_array(self.epsilon_infinity,'epsilon_infinity',(3,3))
        import ctypes as c
        try:fn=_handles_api().library.sp_ld_tensor_validate
        except AttributeError as exc:raise RuntimeFailure('native anisotropic LD ABI missing; install matching wheel or explicitly rebuild') from exc
        fn.argtypes=[c.c_int,c.c_double]+[c.c_void_p]*4;fn.restype=c.c_int
        if fn(len(f),wp,*(a.ctypes.data for a in (f,w,g,eps))):
            raise InputError('tensor LD requires reciprocal symmetric PSD strengths, positive-definite background and nonnegative resonance/damping')
        band=self.valid_frequency_range
        if band is not None:
            if len(band)!=2:raise InputError('valid_frequency_range needs two endpoints')
            band=tuple(_positive(v,'valid_frequency_range',zero=True) for v in band)
            if band[0]>=band[1]:raise InputError('valid_frequency_range must increase')
        object.__setattr__(self,'plasma_frequency',wp)
        object.__setattr__(self,'strength',tuple(tuple(tuple(float(x) for x in row) for row in pole) for pole in f))
        object.__setattr__(self,'resonance',tuple(map(float,w)));object.__setattr__(self,'damping',tuple(map(float,g)))
        object.__setattr__(self,'epsilon_infinity',tuple(tuple(float(x) for x in row) for row in eps))
        object.__setattr__(self,'valid_frequency_range',band)

    @property
    def representation(self):
        if self.component is not None:return 'isotropic_diagonal'
        matrices=(*self.strength,self.epsilon_infinity)
        return 'reciprocal_tensor' if any(a[i][j]!=0. for a in matrices for i in range(3) for j in range(3) if i!=j) else 'diagonal'
    @property
    def independent_components(self):return 1 if self.component is not None else 6 if self.representation=='reciprocal_tensor' else 3
    @property
    def axes(self):return ('x','y','z')
    @property
    def tensor_config(self):return self

    def permittivity_tensor(self,omega):
        return LorentzDrude.from_config(self).permittivity_tensor(omega)


@dataclass(frozen=True, slots=True)
class LDStateConfig:
    """Resolved initialization settings, not the evolving P/J or current field."""
    model: LDModelConfig | LDTensorModelConfig
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

    def __init__(self, *, plasma_frequency, strength, resonance, damping, epsilon_infinity=1.0, valid_frequency_range=None):
        self._wp = _positive(plasma_frequency, 'plasma_frequency', zero=True)
        self._f = _array(strength, 'strength')
        if self._f.ndim != 1 or not self._f.size:
            raise InputError('strength must be a nonempty one-dimensional array')
        self._w = _array(resonance, 'resonance', self._f.shape)
        self._g = _array(damping, 'damping', self._f.shape)
        if any(np.any(a < 0) for a in (self._f, self._w, self._g)):
            raise InputError('strength, resonance and damping must be nonnegative')
        band=None
        if valid_frequency_range is not None:
            bounds=_array(valid_frequency_range,'valid_frequency_range',(2,))
            if bounds[0]<0 or bounds[1]<=bounds[0]:raise InputError('valid_frequency_range must be increasing and nonnegative')
            band=tuple(map(float,bounds))
        self._config = LDModelConfig(self._wp, tuple(map(float, self._f)),
                                     tuple(map(float, self._w)), tuple(map(float, self._g)),
                                     epsilon_infinity=_positive(epsilon_infinity,"epsilon_infinity"),valid_frequency_range=band)
        # Immutable byte-backed arrays cannot be made writable by NumPy callers.
        self._f, self._w, self._g = (np.frombuffer(a.tobytes(), dtype=np.float64)
                                    for a in (self._f, self._w, self._g))

    @classmethod
    def from_config(cls,config):
        """Construct an owned immutable model from a resolved LDModelConfig."""
        if isinstance(config,LDTensorModelConfig):
            if config.component is None:
                return TensorLorentzDrude(plasma_frequency=config.plasma_frequency,strength=config.strength,
                    resonance=config.resonance,damping=config.damping,epsilon_infinity=config.epsilon_infinity,
                    valid_frequency_range=config.valid_frequency_range)
            config=config.component
        if not isinstance(config,LDModelConfig):raise InputError('LDModelConfig or LDTensorModelConfig required')
        return cls(plasma_frequency=config.plasma_frequency,strength=config.strength,
                   resonance=config.resonance,damping=config.damping,
                   epsilon_infinity=config.epsilon_infinity,valid_frequency_range=config.valid_frequency_range)

    @property
    def config(self):
        return self._config

    @property
    def tensor_config(self):
        return self.config.tensor_config

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

    def permittivity_tensor(self, omega):
        """Relative epsilon[..., i, j] in lab xyz; this model is isotropic.

        Shape is omega.shape+(3,3). No anisotropic propagation is implied.
        The scalar permittivity API remains available for compatibility.
        """
        return _isotropic_tensor(self.permittivity(omega))

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

    def check_excitation(self,spectrum,*,rtol=.01):
        """Warn on supplied drive-spectrum weight outside the calibrated band.

        This uses raw equally weighted FFT-bin squared amplitudes, not absorbed
        energy. Supply an E spectrum in atomic time; no new waveform/run is made.
        """
        import ctypes as c
        import warnings
        from .analysis import Spectrum
        band=self.config.valid_frequency_range
        if band is None:return dict(known=False,reason='no calibrated frequency range')
        if not isinstance(spectrum,Spectrum) or spectrum.config.get('time_unit')!='a.u.':
            raise InputError('supply analysis.spectrum with time_unit="a.u."')
        tol=_positive(rtol,'rtol',zero=True)
        omega=np.ascontiguousarray(spectrum.omega,dtype=np.float64)
        values=np.ascontiguousarray(spectrum.values,dtype=np.complex128)
        if values.ndim not in (1,2) or values.shape[0]!=len(omega):raise InputError('invalid spectrum shape')
        components=1 if values.ndim==1 else values.shape[1]
        stats=np.zeros(3)
        try: fn=_handles_api().library.sp_ld_band_fraction
        except AttributeError as exc:raise RuntimeFailure('native LD band ABI missing; rebuild kernels explicitly') from exc
        fn.argtypes=[c.c_int,c.c_int,c.c_void_p,c.c_void_p,c.c_double,c.c_double,c.c_void_p];fn.restype=c.c_int
        if fn(len(omega),components,omega.ctypes.data,values.ctypes.data,*band,stats.ctypes.data):
            raise InputError('invalid drive spectrum or calibrated band')
        within=bool(stats[0]<=tol) if stats[2]>0 else None
        if within is False:warnings.warn('excitation spectrum extends beyond the calibrated LD band; extrapolation is unvalidated',RuntimeWarning,stacklevel=2)
        return dict(known=True,valid_frequency_range=list(band),outside_fraction=float(stats[0]),rtol=tol,
                    within_band=within,scope='raw FFT-bin squared amplitudes, not absorbed energy; finite band only')

    def new_state(self, *, dt, initial_field=(0., 0., 0.), execution_mode="checked"):
        """Independent initially unpolarized state at t=0, with E(0) specified."""
        return LDState(self, dt, initial_field,execution_mode=execution_mode)

    def respond(self, field, *, dt, execution_mode="checked"):
        """E samples of shape (N,3) at t=0,dt,...; output includes initial zeros."""
        execution_mode=_execution_mode_name(execution_mode)
        parser=_array_unchecked if execution_mode=='performance' else _array
        field = parser(field, 'field')
        if field.ndim != 2 or field.shape[1] != 3 or len(field) == 0:
            raise InputError('field must have shape (N,3), N >= 1')
        state = self.new_state(dt=dt, initial_field=field[0],execution_mode=execution_mode)
        try:
            p, j = state._advance(field)
            return Response(np.arange(len(field))*state.dt, p, j)
        finally:
            state.close()


class LDState(_SnapshotMixin):
    """Independent native per-pole history. A state is not thread-shareable."""
    __slots__ = ('_config', '_material', '_field', '_polarization', '_current',
                 '_steps', '_native_state', '_closed', '_execution_mode', '_parse_field')
    _snapshot_kind = 'ld'
    _snapshot_counter = '_steps'

    def _snapshot_owner(self): return self._material
    def _snapshot_settings(self): return _fingerprint(self.config)
    def _snapshot_arrays(self):
        p,j=self._native_state.export()
        return {'_field':self._field,'_p':p,'_j':j}

    def __init__(self, material, dt, initial_field, *, execution_mode="checked"):
        self._execution_mode=_execution_mode_name(execution_mode)
        self._parse_field=_array_unchecked if execution_mode=='performance' else _array
        dt = _positive(dt, 'dt')
        self._material = material
        self._field = _array(initial_field, 'initial_field', (3,))
        self._config = LDStateConfig(material.config, dt, tuple(map(float, self._field)))
        self._polarization=np.zeros(3);self._current=np.zeros(3)
        self._steps = 0;self._closed=False
        from ._ld_handles import NativeState
        self._native_state=NativeState(_handles_api(),material,dt,self._field,execution_mode=execution_mode)

    @property
    def execution_mode(self):return self._execution_mode
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
        field = self._parse_field(field, 'field', (3,))
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
        if self.execution_mode=='checked' and not np.array_equal(field[0],self._field):
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


@dataclass(frozen=True,slots=True)
class FitReport:
    accepted: bool
    optimizer_converged: bool
    status: str
    relative_rms_error: float | None
    weighted_rms_error: float | None
    rtol: float
    max_iterations: int
    kernel_seconds: float
    valid_frequency_range: tuple
    poles: int
    drude: bool
    fitted_background: bool
    iterations: int
    starts: int
    recommended_dt: float | None
    data_kind: str
    data_sha256: str
    source_sha256: str | None
    source_name: str | None
    reference: str | None
    warnings: tuple
    _metadata_json: str = field(repr=False)

    def as_dict(self):
        """Independent JSON-ready report; no guarantee beyond the supplied band."""
        import json
        from dataclasses import asdict
        result=asdict(self)
        result['provenance']=json.loads(result.pop('_metadata_json'))
        from ._version import __version__
        result.update(schema='salmonpy.ld-fit/1',software_version=__version__,units='atomic_gaussian',
                      algorithm=json.loads(self._metadata_json).get('optimizer',{}).get('algorithm','projected LM; forward finite differences; deterministic multistart'),
                      fitted_background_floor=1e-9,iterations_scope='best start only; SciPy counts function evaluations, not iterations',
                      time_convention='exp(-i omega t)',objective='sum |delta epsilon / max(1,abs(epsilon))|^2',
                      timestep_scope='0.1/fastest active pole/plasma/damping scale: accuracy heuristic only',
                      timing_scope='kernel_seconds measures total fitting wall time including optimizer dispatch',
                      stability_margin=None,stability_scope='requires field grid/CFL; not certified by fit',
                      constraints='strength,damping,resonance >= 0; epsilon_infinity > 0',
                      normalization='nonmagnetic mu=1; plasma_frequency=1; fitted oscillator weight is strength',
                      source_scope='source_file is a caller-supplied reference; arrays are identified separately by data_sha256')
        return json.loads(json.dumps(result,allow_nan=False))


def fit_lorentz_drude(omega,epsilon=None,*,poles,units,nk=None,drude=True,
                      epsilon_infinity=True,initial=None,max_iterations=300,starts=4,
                      rtol=.05,source_file=None,reference=None,metadata=None,optimizer="native"):
    """Passive isotropic LD fitting; native objective with native or optional SciPy optimizer.

    poles includes the optional zero-resonance Drude pole. epsilon_infinity=True
    fits a positive background; False fixes 1, a positive real fixes that value.
    Input is epsilon(N,), nk(N,2), or an explicitly masked dielectric OpticalSpectrum.
    Failed or inaccurate fits warn, report diagnostics and never run a solver.
    A returned candidate is passive but is not automatically certified accurate.
    """
    import ctypes as c
    import json,warnings
    if optimizer not in ("native","scipy"):raise InputError("optimizer must be native or scipy")
    if units!='atomic_gaussian':raise InputError('convert angular frequencies to atomic_gaussian units explicitly')
    if type(drude) is not bool:raise InputError('drude must be boolean')
    integers={'poles':poles,'max_iterations':max_iterations,'starts':starts}
    for name,value in integers.items():
        if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or not 1<=value<=(2**31-2)//3:
            raise InputError(name+' must be a positive supported integer')
    poles=int(poles);max_iterations=int(max_iterations);starts=int(starts)
    tol=_positive(rtol,'rtol',zero=True)
    fit_eps=type(epsilon_infinity) is bool and epsilon_infinity
    fixed_eps=1. if type(epsilon_infinity) is bool else _positive(epsilon_infinity,'epsilon_infinity')
    provenance={} if metadata is None else json.loads(json.dumps(metadata,allow_nan=False))
    if not isinstance(provenance,dict):raise InputError('metadata must be a JSON object')
    from .analysis import OpticalSpectrum
    data_kind='epsilon'
    if isinstance(omega,OpticalSpectrum):
        if epsilon is not None or nk is not None or omega.config.get('observable')!='dielectric_function':
            raise InputError('supply a dielectric OpticalSpectrum alone')
        if omega.config.get('units')!=units or omega.config.get('sign')!='exp(+i omega t)':
            raise InputError('dielectric spectrum must use the declared atomic units and Fourier convention')
        provenance['optical_config']=dict(omega.config)
        provenance['excluded_masked_samples']=int(len(omega.valid)-np.count_nonzero(omega.valid))
        epsilon=omega.values[omega.valid];omega=omega.omega[omega.valid];data_kind='optical_spectrum'
    w=_array(omega,'omega')
    if w.ndim!=1 or len(w)<2 or np.any(w<0) or np.any(np.diff(w)<=0) or (drude and np.any(w==0)):
        raise InputError('omega must increase and be nonnegative; Drude fitting requires positive frequencies')
    if len(w)>((2**31-1)//2):raise InputError('too many fit samples')
    width=1+3*poles
    entries=np.iinfo(np.intp).max//np.dtype(np.float64).itemsize
    if width>entries//width or 2*len(w)>entries//width:
        raise InputError('fit workspace dimensions exceed addressable capacity')
    lib=_handles_api().library
    if nk is not None:
        if epsilon is not None:raise InputError('supply epsilon or nk, not both')
        data=_array(nk,'nk',(len(w),2));nr=np.ascontiguousarray(data[:,0]);ki=np.ascontiguousarray(data[:,1])
        target=np.empty(len(w),dtype=np.complex128)
        try: fn=lib.sp_ld_nk
        except AttributeError as exc:raise RuntimeFailure('native LD fit ABI missing; explicitly rebuild kernels') from exc
        fn.argtypes=[c.c_int]+[c.c_void_p]*3;fn.restype=c.c_int
        if fn(len(w),nr.ctypes.data,ki.ctypes.data,target.ctypes.data):
            raise InputError('nk requires finite nonnegative n,k with exp(-i omega t) absorption convention')
        data_kind='nk'
        provenance['nk_assumption_mu']=1.0
    else:
        target=np.array(epsilon,dtype=np.complex128,copy=True,order='C')
        if target.shape!=w.shape or not np.isfinite(target).all():raise InputError('epsilon must be finite with shape (N,)')
    seed=np.zeros(2+3*poles)
    if initial is not None:
        initial=initial.config if isinstance(initial,LorentzDrude) else initial
        model=LorentzDrude.from_config(initial);cfg=model.config
        if isinstance(cfg,LDTensorModelConfig):raise InputError('scalar LD fitting requires one explicitly selected component, not a tensor initial model')
        if len(cfg.strength)!=poles or (drude and cfg.resonance[0]!=0):
            raise InputError('initial model must match pole count and first Drude resonance')
        seed[0]=cfg.epsilon_infinity;seed[1]=cfg.plasma_frequency
        seed[2::3]=cfg.strength;seed[3::3]=cfg.resonance;seed[4::3]=cfg.damping
        from dataclasses import asdict
        provenance['initial_model']=asdict(cfg)
    source_sha=None;source_name=None
    if source_file is not None:
        source=Path(source_file);before=source.stat();digest=hashlib.sha256()
        with source.open('rb') as stream:
            for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
        after=source.stat()
        if (before.st_mtime_ns,before.st_size)!=(after.st_mtime_ns,after.st_size):raise InputError('source file changed while hashing')
        source_sha=digest.hexdigest();source_name=source.name
    if reference is not None and not isinstance(reference,str):raise InputError('reference must be text')
    digest=hashlib.sha256(b'salmonpy.ld-fit-data/1\0')
    digest.update(w.astype('<f8').tobytes());digest.update(target.astype('<c16').tobytes())
    params=np.empty(1+3*poles);stats=np.zeros(4);details=np.zeros(3,dtype=np.int32)
    try: fn=lib.sp_ld_fit
    except AttributeError as exc:raise RuntimeFailure('native LD fit ABI missing; explicitly rebuild kernels') from exc
    fn.argtypes=[c.c_int]*7+[c.c_double]+[c.c_void_p]*6;fn.restype=c.c_int
    from time import perf_counter
    started=perf_counter()
    if optimizer == "scipy":
        status,optimizer_record=_scipy_fit(lib,w,target,poles,fit_eps,drude,fixed_eps,
            seed if initial is not None else None,max_iterations,starts,params,stats,details)
        provenance['optimizer']=optimizer_record
    else:
        status=fn(len(w),poles,int(fit_eps),int(drude),max_iterations,starts,int(initial is not None),fixed_eps,
                  *(a.ctypes.data for a in (w,target,seed,params,stats,details)))
    elapsed=perf_counter()-started
    messages=[];cfg=None;band=(float(w[0]),float(w[-1]))
    if np.any(target.imag<0):messages.append('negative Im epsilon is present; a passive model cannot reproduce gain exactly')
    if 1+3*poles-int(not fit_eps)-int(drude)>2*len(w):messages.append('underdetermined fit; more independent samples are needed')
    if not status:
        cfg=LorentzDrude(plasma_frequency=1.,strength=params[1::3],resonance=params[2::3],damping=params[3::3],
                         epsilon_infinity=params[0],valid_frequency_range=band).config
        accepted=bool(stats[0]<=tol)
        if not accepted:messages.append('relative RMS exceeds the requested tolerance; candidate is passive but inaccurate')
        if not details[1]:messages.append('optimizer did not reach its convergence criterion; inspect residual and pole choices')
    else:
        accepted=False;messages.append('native fit failed with status '+str(status)+'; no model emitted')
    report=FitReport(accepted=accepted,optimizer_converged=bool(details[1]),status='failed' if status else ('accepted' if accepted else 'inaccurate'),
        relative_rms_error=None if status else float(stats[0]),weighted_rms_error=None if status else float(stats[1]),
        rtol=tol,max_iterations=max_iterations,kernel_seconds=elapsed,
        valid_frequency_range=band,poles=poles,drude=drude,fitted_background=fit_eps,iterations=int(details[0]),starts=int(details[2]),
        recommended_dt=None if status or stats[3]==0 else float(stats[3]),data_kind=data_kind,data_sha256=digest.hexdigest(),
        source_sha256=source_sha,source_name=source_name,reference=reference,warnings=tuple(messages),
        _metadata_json=json.dumps(provenance,allow_nan=False))
    if messages:warnings.warn('; '.join(messages),RuntimeWarning,stacklevel=2)
    return cfg,report


def _scipy_fit(lib,w,target,poles,fit_eps,drude,fixed_eps,initial,max_iterations,starts,params,stats,details):
    """SciPy dispatch only: Fortran owns the objective, seeds and diagnostics."""
    import ctypes as c
    try:
        import scipy
        from scipy.optimize import least_squares
    except ImportError as exc:
        raise RuntimeFailure('SciPy fitting requires the optional salmonpy[fit] dependency') from exc
    try:
        seed_fn=lib.sp_ld_fit_seed
        evaluate=lib.sp_ld_fit_evaluate
    except AttributeError as exc:
        raise RuntimeFailure('SciPy native objective ABI missing; explicitly rebuild kernels') from exc
    seed_fn.argtypes=[c.c_int]*5+[c.c_double]+[c.c_void_p]*3
    evaluate.argtypes=[c.c_int]*2+[c.c_void_p]*6
    seed_fn.restype=evaluate.restype=c.c_int
    full=np.empty(1+3*poles);residual=np.empty(2*len(w))
    physical=np.empty_like(full);scores=np.empty(4)
    free=np.ones(len(full),dtype=bool)
    free[0]=fit_eps
    if drude:free[2]=False
    lower=np.zeros(int(np.count_nonzero(free)))
    if fit_eps:lower[0]=1e-9
    attempts=[];best=None
    for start in range(1,starts+1):
        if seed_fn(len(w),poles,int(fit_eps),int(drude),start,fixed_eps,w.ctypes.data,target.ctypes.data,full.ctypes.data):
            continue
        if initial is not None and start==1:
            scale=float(w[-1]);full[0]=initial[0]
            full[1::3]=np.multiply(initial[2::3],np.square(np.divide(initial[1],scale)))
            full[2::3]=np.divide(initial[3::3],scale)
            full[3::3]=np.divide(initial[4::3],scale)
        if not fit_eps:full[0]=fixed_eps
        if drude:full[2]=0.
        def objective(values):
            full[free]=values
            status=evaluate(len(w),poles,*(a.ctypes.data for a in (w,target,full,residual,physical,scores)))
            if status:raise FloatingPointError('nonfinite or singular native LD residual')
            return residual.copy()
        try:
            result=least_squares(objective,np.maximum(full[free],lower),bounds=(lower,np.inf),
                method='trf',jac='2-point',x_scale='jac',max_nfev=max_iterations,
                ftol=1e-10,xtol=1e-10,gtol=1e-10)
            objective(result.x)
        except (FloatingPointError,ValueError) as exc:
            attempts.append(dict(start=start,status='failed',message=str(exc)))
            continue
        attempt=dict(start=start,status=int(result.status),success=bool(result.success),
            message=str(result.message),nfev=int(result.nfev),njev=None if result.njev is None else int(result.njev),
            optimality=float(result.optimality),weighted_cost=float(scores[2]))
        attempts.append(attempt)
        if best is None or scores[2]<best[0]:
            best=(float(scores[2]),physical.copy(),scores.copy(),attempt)
    record=dict(name='scipy',version=scipy.__version__,algorithm='scipy.optimize.least_squares TRF; native objective; deterministic multistart',
        method='trf',jacobian='2-point',x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-10,
        budget_scope='max_iterations is max_nfev per start; numerical Jacobian calls are additional',attempts=attempts)
    if best is None:return 4,record
    _,best_params,best_stats,chosen=best
    params[:]=best_params
    stats[:]=best_stats
    details[:]=[chosen['nfev'],int(chosen['success']),len(attempts)]
    record['selected_start']=chosen['start']
    return 0,record


def _isotropic_tensor(epsilon):
    """Native isotropic embedding; never discard complex absorption components."""
    import ctypes as c
    values=np.asarray(epsilon,dtype=np.complex128,order='C')
    output=np.empty(values.shape+(3,3),dtype=np.complex128)
    if not 0<values.size<=2**31-1:raise InputError('nonempty supported tensor sample count required')
    try:fn=_handles_api().library.sp_ld_tensor
    except AttributeError as exc:
        raise RuntimeFailure('native dielectric tensor ABI missing; install matching wheel or explicitly rebuild kernels') from exc
    fn.argtypes=[c.c_int,c.c_void_p,c.c_void_p];fn.restype=c.c_int
    if fn(values.size,values.ctypes.data,output.ctypes.data):
        raise InputError('dielectric tensor values must be finite')
    return np.frombuffer(output.tobytes(),dtype=np.complex128).reshape(output.shape)


class TensorLorentzDrude(LorentzDrude):
    """Reciprocal passive LD tensor with independent owned vector histories.

    epsilon = epsilon_infinity + sum(wp**2 * strength_q /
              (resonance_q**2 - omega**2 - i*damping_q*omega)).
    Real symmetric pole strengths are PSD; background is positive definite.
    """
    __slots__=('_strength_matrix','_epsilon_matrix')

    def __init__(self,*,plasma_frequency,strength,resonance,damping,
                 epsilon_infinity=((1.,0.,0.),(0.,1.,0.),(0.,0.,1.)),valid_frequency_range=None):
        config=LDTensorModelConfig(plasma_frequency=plasma_frequency,strength=strength,
            resonance=resonance,damping=damping,epsilon_infinity=epsilon_infinity,
            valid_frequency_range=valid_frequency_range)
        self._config=config;self._wp=config.plasma_frequency
        self._w=np.array(config.resonance);self._g=np.array(config.damping)
        self._strength_matrix=np.array(config.strength);self._epsilon_matrix=np.array(config.epsilon_infinity)
        self._f=np.zeros(len(self._w))  # legacy n/layout slot; tensor creator uses matrices
        for name in ('_w','_g','_f','_strength_matrix','_epsilon_matrix'):
            a=getattr(self,name);setattr(self,name,np.frombuffer(a.tobytes(),dtype=float).reshape(a.shape))

    def permittivity(self,omega):
        raise InputError('anisotropic response requires permittivity_tensor; select components explicitly')

    def permittivity_tensor(self,omega):
        import ctypes as c
        values=_array(omega,'omega');output=np.empty(values.shape+(3,3),dtype=np.complex128)
        try:fn=_handles_api().library.sp_ld_tensor_permittivity
        except AttributeError as exc:raise RuntimeFailure('native tensor observable ABI missing; explicitly rebuild') from exc
        fn.argtypes=[c.c_int,c.c_int,c.c_double]+[c.c_void_p]*6;fn.restype=c.c_int
        if values.size>2**31-1 or fn(len(self._w),values.size,self._wp,
                *(a.ctypes.data for a in (self._epsilon_matrix,self._strength_matrix,self._w,self._g,values,output))):
            raise InputError('tensor permittivity is singular or nonfinite')
        return np.frombuffer(output.tobytes(),dtype=np.complex128).reshape(output.shape)

    def _background(self,field,polarization,magnetic_field):
        import ctypes as c
        e=_array(field,'field');p=_array(polarization,'polarization',e.shape)
        h=_array(magnetic_field,'magnetic_field',e.shape)
        if e.ndim<1 or e.shape[-1]!=3:raise InputError('field arrays require(...,3)')
        d=np.empty_like(e);u=np.empty(e.shape[:-1])
        fn=_handles_api().library.sp_ld_tensor_background
        fn.argtypes=[c.c_int]+[c.c_void_p]*6;fn.restype=c.c_int
        if e.size//3>2**31-1 or fn(e.size//3,*(a.ctypes.data for a in (self._epsilon_matrix,e,p,h,d,u))):
            raise NumericalError('nonfinite tensor background observable')
        return d,u


def fit_lorentz_drude_tensor(omega, epsilon=None, *, initial, units,
                             max_evaluations=1000, rtol=.05, metadata=None):
    """Fit all nine components with shared poles and PSD oscillator strengths.

    Requires an explicit TensorLorentzDrude starting model and optional SciPy.
    The initial SPD background is fixed; pole frequencies/damping and strength
    factors are fitted. Native Fortran computes matrices and weighted residuals.
    Antisymmetric data are retained in the residual, never symmetrized away.
    """
    import json
    import warnings
    try:
        from scipy.optimize import least_squares
    except ImportError as exc:
        raise RuntimeFailure("tensor fitting requires the optional scipy dependency") from exc
    from .analysis import OpticalSpectrum
    if units != 'atomic_gaussian' or not isinstance(initial, TensorLorentzDrude):
        raise InputError('declare atomic_gaussian units and provide a TensorLorentzDrude initial model')
    spectrum_provenance = None
    if isinstance(omega, OpticalSpectrum):
        spectrum_provenance = dict(omega.config)
        if epsilon is not None or omega.config.get('observable') != 'dielectric_tensor':
            raise InputError('supply a dielectric_tensor spectrum alone')
        if omega.config.get('units') != units:
            raise InputError('spectrum units disagree')
        keep = np.all(omega.valid, axis=(1, 2))
        epsilon = omega.values[keep]
        omega = omega.omega[keep]
    w = _array(omega, 'omega')
    data = np.array(epsilon, dtype=np.complex128, copy=True)
    if w.ndim != 1 or len(w) < 2 or data.shape != (len(w), 3, 3):
        raise InputError('positive frequencies and all nine epsilon components are required')
    if not np.isfinite(data).all() or np.any(w <= 0) or np.any(np.diff(w) <= 0):
        raise InputError('frequencies must increase and tensor data must be finite')
    tol = _positive(rtol, 'rtol', zero=True)
    if isinstance(max_evaluations, bool) or not isinstance(max_evaluations, (int, np.integer)) or max_evaluations < 1:
        raise InputError('max_evaluations must be a positive integer')
    provenance = json.loads(json.dumps({} if metadata is None else metadata, allow_nan=False))
    if not isinstance(provenance, dict):
        raise InputError('metadata must be a JSON object')
    # Eigh decomposition is compiled LAPACK; PSD roundoff is clipped only within
    # the already validated initial model. QR produces a lower strength factor.
    parameters = []
    lower = []
    upper = []
    for strength, resonance, damping in zip(initial.config.strength, initial.config.resonance, initial.config.damping):
        values, vectors = np.linalg.eigh(np.asarray(strength) * initial.config.plasma_frequency**2)
        factor = vectors * np.sqrt(np.maximum(values, 0))
        _, triangular = np.linalg.qr(factor.T)
        factor = triangular.T
        parameters.extend([factor[0, 0], factor[1, 0], factor[1, 1], factor[2, 0],
                           factor[2, 1], factor[2, 2], resonance, damping])
        lower.extend([-np.inf]*6 + [0., 0.])
        upper.extend([np.inf]*8)
    p = len(initial.config.resonance)
    target = np.asfortranarray(data.transpose(1, 2, 0))
    background = np.asfortranarray(initial.config.epsilon_infinity, dtype=np.float64)
    residual = np.empty(18*len(w))
    strength = np.empty((3, 3, p), order='F')
    resonance = np.empty(p)
    damping = np.empty(p)
    try:
        fn = _handles_api().library.sp_ld_tensor_fit
    except AttributeError as exc:
        raise RuntimeFailure('native tensor fit ABI unavailable; explicitly rebuild kernels') from exc
    fn.argtypes = [ctypes.c_int, ctypes.c_int] + [ctypes.c_void_p]*8
    fn.restype = ctypes.c_int
    def objective(parameters):
        parameters = np.ascontiguousarray(parameters)
        if fn(len(w), p, *(a.ctypes.data for a in
                           (w, target, background, parameters, residual, strength, resonance, damping))):
            raise NumericalError('nonfinite tensor fit evaluation')
        return residual.copy()
    result = least_squares(objective, parameters, bounds=(lower, upper),
                           max_nfev=int(max_evaluations), x_scale='jac')
    final_residual = objective(result.x)
    candidate = TensorLorentzDrude(plasma_frequency=1., strength=strength.transpose(2, 0, 1),
                    resonance=resonance, damping=damping, epsilon_infinity=background,
                    valid_frequency_range=(float(w[0]), float(w[-1])))
    error = float(np.sqrt(np.mean(final_residual**2)*18))
    accepted = error <= tol
    if not accepted:
        warnings.warn('tensor LD fit exceeds requested RMS tolerance; inspect report before use', RuntimeWarning, stacklevel=2)
    report = dict(schema='salmonpy.tensor-ld-fit/1', accepted=accepted,
                  optimizer_converged=bool(result.success), evaluations=int(result.nfev),
                  weighted_rms_error=error, rtol=tol, fixed_background=True,
                  constraints='shared nonnegative resonance/damping; F=L L^T; reciprocal passive model',
                  objective='RMS Frobenius error / max(1, Frobenius norm of epsilon)',
                  valid_frequency_range=[float(w[0]), float(w[-1])],
                  data_sha256=hashlib.sha256(w.tobytes()+data.tobytes()).hexdigest(),
                  provenance=provenance, input_spectrum=spectrum_provenance,
                  initial_model=asdict(initial.config))
    if not result.success:
        warnings.warn('tensor optimizer did not converge; inspect residual and pole choices', RuntimeWarning, stacklevel=2)
    return candidate, report
