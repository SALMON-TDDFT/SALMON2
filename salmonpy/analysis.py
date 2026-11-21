"""SALMON table reading and explicitly normalized Fourier transforms."""
from .errors import InputError
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
import re
import numpy as np
from .ld import _array
from .fields import _readonly, _times


@dataclass(frozen=True,slots=True)
class Table:
    columns: tuple
    units: tuple
    values: np.ndarray
    headers: tuple
    source: str

    def column(self,name):
        return self.values[:,self.columns.index(name)]


def read_table(path):
    """Read numbered SALMON headers, Fortran D exponents and comment lines.

    Units are retained, never inferred or converted. Unlabeled columns are named
    column_1, column_2,... with unknown units represented by None.
    """
    path=Path(path);headers=[];rows=[];labels={}
    for line in path.read_text().splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith(('#','!')):
            headers.append(line)
            for n,label,unit in re.findall(r'(\d+):([^\s\[]+)\[([^\]]+)\]',line):
                labels[int(n)-1]=(label,unit)
            continue
        clean=re.split(r'[#!]',line,maxsplit=1)[0]
        rows.append([float(x.replace('D','E').replace('d','e')) for x in clean.split()])
    if not rows: raise InputError('table contains no samples')
    data=_array(rows,'table')
    if data.ndim!=2 or data.shape[1]==0: raise InputError('table must be rectangular')
    names=tuple(labels.get(i,(f'column_{i+1}',None))[0] for i in range(data.shape[1]))
    if len(set(names))!=len(names): raise InputError('column names must be unique')
    units=tuple(labels.get(i,('',None))[1] for i in range(data.shape[1]))
    return Table(names,units,_readonly(data),tuple(headers),str(path.resolve()))


@dataclass(frozen=True,slots=True)
class Spectrum:
    omega: np.ndarray
    values: np.ndarray
    config: object


def spectrum(time,values,*,window='hann',time_unit='a.u.',value_unit='a.u.'):
    """dt * sum(w[n]*x[n]*exp(+i*omega*t[n])), positive frequencies only.

    No detrending, zero-padding, one-sided doubling or window-gain correction.
    This is a signal spectrum, not an automatically normalized optical observable.
    """
    t=_times(time);x=_array(values,'values');dt=np.diff(t)
    if x.ndim not in (1,2) or len(x)!=len(t): raise InputError('values must have shape (N,) or (N,M)')
    if not np.allclose(dt,dt[0],rtol=1e-8,atol=1e-12*abs(dt[0])):
        raise InputError('uniform time samples are required')
    if window=='none': w=np.ones(len(t))
    elif window=='hann':
        if len(t)<3: raise InputError('hann requires at least three samples')
        w=np.hanning(len(t))
    else: raise InputError('window must be none or hann')
    xw=x*w if x.ndim==1 else x*w[:,None]
    omega=2*np.pi*np.fft.rfftfreq(len(t),dt[0])
    values=np.conj(np.fft.rfft(xw,axis=0))*dt[0]
    phase=np.exp(1j*omega*t[0]);values*=phase if x.ndim==1 else phase[:,None]
    cfg=MappingProxyType(dict(window=window,sign='exp(+i omega t)',normalization='dt * sum',
                              dt=float(dt[0]),start=float(t[0]),samples=len(t),
                              time_unit=time_unit,value_unit=value_unit,
                              frequency_unit=f'radian/({time_unit})',spectrum_unit=f'({value_unit})*({time_unit})'))
    omega=_readonly(omega)
    values=np.frombuffer(values.tobytes(),dtype=values.dtype).reshape(values.shape)
    return Spectrum(omega,values,cfg)


@dataclass(frozen=True, slots=True)
class OpticalSpectrum:
    """Scalar or tensor optical response; invalid entries carry NaN and valid=False."""
    omega: np.ndarray
    values: np.ndarray
    valid: np.ndarray
    config: object
    _dielectric_values: object = field(default=None, repr=False)
    _dielectric_valid: object = field(default=None, repr=False)


def optical_conductivity(time, current, electric_field, *, current_convention,
                         field_kind, units, window='hann', drive_rtol=1e-8,
                         metadata=None):
    """Return sigma=j/E for one explicitly selected polarization, in Gaussian AU.

    Inputs are simultaneous scalar time series. j must be charge-current density
    and E the total macroscopic electric field, not external E or vector potential.
    This is a finite-time estimate, not a certified linear-response tensor.
    Existing spectrum() handles FFT in NumPy's compiled backend; native Fortran
    performs thresholding and observable conversion. No solver is launched.
    """
    import ctypes as c
    import json
    from .errors import RuntimeFailure
    from .ld import _positive
    from ._sbe_handles import default_api
    if current_convention != 'charge':
        raise InputError('use charge-current density j; convert SALMON matter current J=-j explicitly')
    if field_kind != 'total':
        raise InputError('optical conductivity requires total E; external E or A is not interchangeable')
    if units != 'atomic_gaussian':
        raise InputError('optical observables require atomic_gaussian units; convert explicitly')
    threshold=_positive(drive_rtol,'drive_rtol',zero=True)
    if threshold>=1: raise InputError('drive_rtol must be less than one')
    j=_array(current,'current');e=_array(electric_field,'electric_field')
    if j.ndim!=1 or e.shape!=j.shape:
        raise InputError('select one current component and one driving E polarization, both shape (N,)')
    provenance=json.loads(json.dumps({} if metadata is None else metadata,allow_nan=False))
    if not isinstance(provenance,dict): raise InputError('metadata must be a JSON object')
    js=spectrum(time,j,window=window);es=spectrum(time,e,window=window)
    n=len(js.omega)
    sigma=np.empty(n,dtype=np.complex128);eps=np.empty_like(sigma)
    mask=np.empty(n,dtype=np.int32)
    try: fn=default_api().library.sp_optical_ratio
    except AttributeError as exc:
        raise RuntimeFailure('native optical ABI unavailable; install a matching wheel or explicitly rebuild kernels') from exc
    fn.argtypes=[c.c_int,c.c_void_p,c.c_void_p,c.c_void_p,c.c_double,c.c_void_p,c.c_void_p,c.c_void_p]
    fn.restype=c.c_int
    status=fn(n,js.omega.ctypes.data,js.values.ctypes.data,es.values.ctypes.data,
              threshold,sigma.ctypes.data,eps.ctypes.data,mask.ctypes.data)
    if status: raise InputError('native optical conversion rejected input')
    cfg=dict(js.config,observable='optical_conductivity',units=units,
             current_convention=current_convention,field_kind=field_kind,
             gauge='electric_field',drive_rtol=threshold,formula='sigma=j/E',
             value_unit='1/atomic_time',duration=float(time[-1]-time[0]),
             peak_abs_E=float(np.max(np.abs(e))),
             spectrum_unit='1/atomic_time',
             response_scope='single polarization; linearity and convergence not established',
             metadata=provenance)
    # Byte-backed arrays preserve the same immutability as spectrum().
    def frozen(a): return np.frombuffer(a.tobytes(),dtype=a.dtype).reshape(a.shape)
    return OpticalSpectrum(js.omega,frozen(sigma),frozen(mask!=0),MappingProxyType(cfg),
                           frozen(eps),frozen(mask==3))


def dielectric_function(conductivity):
    """epsilon=1+4*pi*i*sigma/omega, exp(-i*omega*t); DC is undefined.

    Accepts optical_conductivity() output, preserving its drive and provenance.
    A separate incident/external-field ratio is not a bulk dielectric function.
    """
    if not isinstance(conductivity,OpticalSpectrum) or conductivity.config.get('observable') not in ('optical_conductivity','optical_conductivity_tensor'):
        raise InputError('pass the result of optical_conductivity')
    if conductivity._dielectric_values is None or conductivity._dielectric_valid is None:
        raise InputError('conductivity lacks native dielectric conversion data')
    cfg=dict(conductivity.config)
    values=conductivity._dielectric_values;valid=conductivity._dielectric_valid
    cfg.update(observable='dielectric_tensor' if conductivity.config.get('observable')=='optical_conductivity_tensor' else 'dielectric_function',formula='epsilon=1+4*pi*i*sigma/omega',value_unit='dimensionless',spectrum_unit='dimensionless')
    return OpticalSpectrum(conductivity.omega,values,valid,MappingProxyType(cfg))


def _optical_diagnostic(reference, candidate):
    import ctypes as c
    from ._sbe_handles import default_api
    from .errors import RuntimeFailure
    for result in (reference,candidate):
        if not isinstance(result,OpticalSpectrum) or result.config.get('observable')!='optical_conductivity':
            raise InputError('diagnostics require optical_conductivity results')
    if not np.array_equal(reference.omega,candidate.omega):
        raise InputError('diagnostics require identical frequency grids; no implicit interpolation')
    a=np.ascontiguousarray(reference.values);b=np.ascontiguousarray(candidate.values)
    ma=np.ascontiguousarray(reference.valid,dtype=np.int32)
    mb=np.ascontiguousarray(candidate.valid,dtype=np.int32);stats=np.zeros(4)
    try: fn=default_api().library.sp_optical_diagnostic
    except AttributeError as exc: raise RuntimeFailure('rebuild the native optical kernel explicitly') from exc
    fn.argtypes=[c.c_int]+[c.c_void_p]*6;fn.restype=c.c_int
    status=fn(len(a),reference.omega.ctypes.data,a.ctypes.data,b.ctypes.data,
              ma.ctypes.data,mb.ctypes.data,stats.ctypes.data)
    if status: raise InputError('invalid optical diagnostic input')
    return stats


def optical_linearity(reference, candidate, *, rtol):
    """Compare supplied amplitude-normalized conductivities; never reruns a solver.

    Same dt/window/frequency grid and polarization/physical configuration are
    required. A passing result is evidence for these two drives, not a universal
    linear-response certificate. The caller owns physical metadata consistency.
    """
    from .ld import _positive
    tol=_positive(rtol,'rtol',zero=True)
    if not all(isinstance(x,OpticalSpectrum) for x in (reference,candidate)):
        raise InputError('linearity requires optical_conductivity results')
    for key in ('dt','start','samples','window','units','field_kind','current_convention'):
        if reference.config.get(key)!=candidate.config.get(key):
            raise InputError('linearity requires matching '+key)
    stats=_optical_diagnostic(reference,candidate)
    error=float(stats[1]/stats[2]) if stats[2]>0 else None
    return dict(relative_rms_error=error,rtol=tol,
                accepted=None if error is None else bool(error<=tol),
                compared_frequencies=int(stats[0]),
                reference_peak_abs_E=reference.config.get('peak_abs_E'),
                candidate_peak_abs_E=candidate.config.get('peak_abs_E'),
                scope='supplied drives and common valid frequencies only; check matching material/polarization')


def optical_sum_rule(conductivity, *, expected_weight, normalization):
    """Finite-band trapezoidal integral of Re sigma over adjacent valid bins.

    expected_weight and its normalization are explicitly supplied by the caller.
    For all-electron Gaussian AU the positive-frequency f-sum is pi*n/2 (m=1).
    Pseudopotentials, nonlocal current operators and missing bands require an
    appropriate reference. No universal pass/fail threshold is imposed.
    """
    from .ld import _positive
    expected=_positive(expected_weight,'expected_weight')
    if not isinstance(normalization,str) or not normalization.strip():
        raise InputError('describe the applicable density/current-operator normalization')
    stats=_optical_diagnostic(conductivity,conductivity)
    return dict(integrated_weight=float(stats[3]),expected_weight=expected,
                captured_fraction=float(stats[3]/expected),normalization=normalization,
                frequency_min=float(conductivity.omega[0]),frequency_max=float(conductivity.omega[-1]),
                valid_frequencies=int(stats[0]),
                scope='finite sampled band; masked gaps omitted, no endpoint extrapolation or missing-band correction',
                certified=False)


def waveplate_diagnostic(time,electric_field,*,omega,indices,thickness,units,angle=0.):
    """Native carrier phasors/Stokes vs a lossless normal-incidence vacuum slab.

    angle is the principal-axis rotation about lab z, in radians. At angle=0,
    E[:,0:2] are the principal transverse components; incident amplitudes
    are equal and in phase. Full pulse/ringdown must lie in the supplied time
    window. No windowing or automatic convergence claim is made. Fresnel
    reflections and the slab's multiple reflections are included analytically.
    """
    import ctypes as c
    from .ld import _handles_api,_positive
    from .errors import RuntimeFailure
    if units!='atomic_gaussian':raise InputError('waveplate diagnostics require atomic_gaussian')
    t=_times(time);e=_array(electric_field,'electric_field',(len(t),3))
    if len(t)>2**31-1:raise InputError('too many diagnostic samples')
    frequency=_positive(omega,'omega');length=_positive(thickness,'thickness')
    if len(indices)!=2:raise InputError('two positive real principal refractive indices required')
    nx,ny=(_positive(v,'indices') for v in indices)
    z=np.empty(2,dtype=np.complex128);s=np.empty(3);reference=np.empty(3);stats=np.empty(4)
    rotation=_array(angle,'angle',())
    try:
        fn=(_handles_api().library.sp_waveplate_diagnostic if float(rotation)==0. else
            _handles_api().library.sp_waveplate_rotated)
    except AttributeError as exc:raise RuntimeFailure('native waveplate diagnostic ABI missing; explicitly rebuild kernels') from exc
    scalar_args=[frequency,nx,ny,length]
    if float(rotation)!=0.:scalar_args.append(float(rotation))
    fn.argtypes=[c.c_int,c.c_void_p,c.c_void_p]+[c.c_double]*len(scalar_args)+[c.c_void_p]*4;fn.restype=c.c_int
    if fn(len(t),t.ctypes.data,e.ctypes.data,*scalar_args,*(a.ctypes.data for a in (z,s,reference,stats))):
        raise InputError('nonfinite, zero-field or invalid waveplate diagnostic')
    return dict(phasor=z,phasor_unit='E:au * time:au',phasor_normalization='trapezoidal Fourier integral, not peak field amplitude',stokes_normalized=s,analytic_stokes_normalized=reference,
        geometric_retardance=float(stats[0]),measured_phase_y_minus_x=float(stats[1]),
        analytic_phase_y_minus_x=float(stats[2]),stokes_error_norm=float(stats[3]),
        stokes_convention='S3=2 Im(Ex conj(Ey))/S0; phasor uses exp(+i omega t)',
        angle_radians=float(rotation),phase_axes='principal',
        scope='lossless rotated slab in vacuum; equal in-phase incident principal amplitudes; Stokes/phasor in lab xy; not a convergence certificate')


def dielectric_tensor(time, current, electric_field, *, current_convention,
                      field_kind, units, window='hann', drive_rtol=1e-8,
                      metadata=None):
    """Three independent runs -> epsilon[frequency, current axis, field axis].

    Histories have shape (3, N, 3): experiment, simultaneous sample, Cartesian
    component. Native inversion solves sigma = J E^-1, using total E and charge
    current. DC, weak or rank-deficient excitation is masked; missing components
    are never inferred. Linear response and physical convergence remain unchecked.
    """
    import ctypes as c
    import json
    from ._sbe_handles import default_api
    from .ld import _positive
    from .errors import RuntimeFailure
    if current_convention != 'charge' or field_kind != 'total' or units != 'atomic_gaussian':
        raise InputError('declare charge-current density, total E and atomic_gaussian units explicitly')
    threshold = _positive(drive_rtol, 'drive_rtol', zero=True)
    if threshold >= 1:
        raise InputError('drive_rtol must be less than one')
    j = _array(current, 'current')
    e = _array(electric_field, 'electric_field')
    if j.shape != (3, len(time), 3) or e.shape != j.shape:
        raise InputError('provide all three independent runs with shape (3, N, 3)')
    provenance = json.loads(json.dumps({} if metadata is None else metadata, allow_nan=False))
    if not isinstance(provenance, dict):
        raise InputError('metadata must be a JSON object')
    js = [spectrum(time, run, window=window) for run in j]
    es = [spectrum(time, run, window=window) for run in e]
    omega = js[0].omega
    # Fortran axes: current/field component, experiment, frequency.
    jf = np.asfortranarray(np.stack([run.values for run in js], axis=2).transpose(1, 2, 0))
    ef = np.asfortranarray(np.stack([run.values for run in es], axis=2).transpose(1, 2, 0))
    sigma = np.empty(jf.shape, dtype=np.complex128, order='F')
    epsilon = np.empty_like(sigma, order='F')
    mask = np.empty(len(omega), dtype=np.int32)
    try:
        fn = default_api().library.sp_optical_tensor
    except AttributeError as exc:
        raise RuntimeFailure('native tensor optical ABI unavailable; explicitly rebuild kernels') from exc
    fn.argtypes = [c.c_int, c.c_void_p, c.c_void_p, c.c_void_p, c.c_double,
                   c.c_void_p, c.c_void_p, c.c_void_p]
    fn.restype = c.c_int
    if fn(len(omega), omega.ctypes.data, jf.ctypes.data, ef.ctypes.data, threshold,
          sigma.ctypes.data, epsilon.ctypes.data, mask.ctypes.data):
        raise InputError('native tensor optical conversion rejected input')
    values = epsilon.transpose(2, 0, 1).copy()
    valid = np.broadcast_to((mask != 0)[:, None, None], values.shape).copy()
    values[~valid] = complex(float('nan'), float('nan'))
    def frozen(a):
        return np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape)
    config = dict(js[0].config, observable='dielectric_tensor', units=units,
                  current_convention=current_convention, field_kind=field_kind,
                  axes=('x', 'y', 'z'), independent_runs=3, drive_rtol=threshold,
                  formula='epsilon=I+4*pi*i*(J E^-1)/omega',
                  value_unit='dimensionless', metadata=provenance,
                  response_scope='finite-time three-drive estimate; linearity and convergence unchecked')
    return OpticalSpectrum(omega, frozen(values), frozen(valid), MappingProxyType(config))


@dataclass(frozen=True, slots=True)
class LongWavelengthCorrection:
    """Independent observation results; neither solver nor input is modified."""
    raw_conductivity: OpticalSpectrum
    raw_dielectric: OpticalSpectrum
    conductivity: OpticalSpectrum
    dielectric: OpticalSpectrum
    applied: np.ndarray
    report: object


def long_wavelength_correction(observation, *, fit_range, system, units,
                               cutoff=None, absorption='preserve'):
    """Observation-only PRB85,045134 Appendix A4 insulating correction.

    Fit Re epsilon=c0+c2*omega**2 in a caller-selected trusted frequency band.
    Below cutoff, correct Im sigma consistently. Default preserves absorption
    as in the paper. absorption='analytic' additionally fits Re sigma starting
    at omega**2: an explicit extra insulating assumption, not the paper's recipe.
    Only declared ordinary zero-DC gapped response is supported; never call for
    Drude/Hall conductivity or photoexcited free carriers. No automatic band-gap
    inference, experimental static preset or global causality claim is made.
    All numerical conversion/fitting/correction is Fortran. Inputs are spectra,
    never TDDFT states or feedback currents. Frequencies are angular Gaussian AU.
    """
    import ctypes as c
    import json
    import warnings
    from .ld import _positive
    from ._sbe_handles import default_api
    from .errors import RuntimeFailure
    if system != 'gapped_insulator':
        raise InputError('A4 correction requires a declared zero-DC gapped_insulator; retain metal/free-carrier/Hall data')
    if units != 'atomic_gaussian':
        raise InputError('convert frequency bounds to atomic_gaussian angular frequencies explicitly')
    if absorption not in ('preserve','analytic'):
        raise InputError('absorption must be preserve or analytic')
    if not isinstance(observation,OpticalSpectrum):
        raise InputError('pass an OpticalSpectrum conductivity or dielectric result')
    name=observation.config.get('observable')
    kinds={'optical_conductivity':0,'optical_conductivity_tensor':0,'dielectric_function':1,'dielectric_tensor':1}
    if name not in kinds:raise InputError('conductivity or dielectric observation required')
    if observation.config.get('units')!=units or observation.config.get('sign')!='exp(+i omega t)':
        raise InputError('observation units/Fourier convention disagree')
    kind=kinds[name]
    if kind==0 and (observation.config.get('current_convention')!='charge' or observation.config.get('field_kind')!='total'):
        raise InputError('use charge-current conductivity referred to total E')
    bounds=_array(fit_range,'fit_range',(2,))
    if bounds[0]<=0 or bounds[1]<=bounds[0]:raise InputError('require 0 < fit_range[0] < fit_range[1]')
    edge=float(bounds[0]) if cutoff is None else _positive(cutoff,'cutoff')
    omega=_array(observation.omega,'omega')
    data=np.array(observation.values,dtype=np.complex128,copy=True)
    mask=np.array(observation.valid,dtype=bool,copy=True)
    if omega.ndim!=1 or len(omega)==0 or data.shape not in ((len(omega),),(len(omega),3,3)) or mask.shape!=data.shape:
        raise InputError('spectrum shape must be (N,) or (N,3,3) with matching masks and omega')
    shape=data.shape;channels=1 if data.ndim==1 else 9
    if len(omega)>2**31-1:raise InputError('too many frequency bins')
    values=np.asfortranarray(data.reshape(len(omega),channels).T)
    valid=np.asfortranarray(mask.reshape(len(omega),channels).T,dtype=np.int32)
    diagonal=np.ones(1,dtype=np.int32) if channels==1 else np.eye(3,dtype=np.int32).ravel()
    raw_sigma=np.empty_like(values,order='F');raw_eps=np.empty_like(values,order='F')
    sigma=np.empty_like(values,order='F');eps=np.empty_like(values,order='F')
    raw_bits=np.empty(values.shape,dtype=np.int32,order='F');bits=np.empty_like(raw_bits,order='F')
    applied=np.empty_like(raw_bits,order='F');coeff=np.empty((4,channels),order='F')
    stats=np.empty((4,channels),order='F');counts=np.empty((3,channels),dtype=np.int32,order='F')
    try:fn=default_api().library.sp_optical_a4
    except AttributeError as exc:raise RuntimeFailure('native observation A4 ABI unavailable; explicitly rebuild/install matching kernels') from exc
    fn.argtypes=[c.c_int]*4+[c.c_void_p]*4+[c.c_double]*3+[c.c_void_p]*10
    fn.restype=c.c_int
    status=fn(len(omega),channels,kind,int(absorption=='analytic'),
        *(x.ctypes.data for x in (omega,values,valid,diagonal)),float(bounds[0]),float(bounds[1]),edge,
        *(x.ctypes.data for x in (raw_sigma,raw_eps,sigma,eps,raw_bits,bits,applied,coeff,stats,counts)))
    if status:raise InputError('native A4 observation correction rejected inconsistent frequencies, masks or valid data')
    notices=[]
    if edge>bounds[1]:notices.append('correction extends beyond the upper fitted frequency; inspect extrapolation')
    if np.any(counts[1]!=0):notices.append('some components lack a usable fit or correction; their unavailable bins retain raw values/masks')
    for notice in notices:warnings.warn(notice,RuntimeWarning,stacklevel=2)
    def number(x):return float(x) if np.isfinite(x) else None
    components=[]
    for j in range(channels):
        fitted=counts[0,j]>=2 and counts[1,j]!=1
        components.append(dict(index=[j//3,j%3] if channels==9 else [],fit_points=int(counts[0,j]),
            status=('fitted','insufficient_fit_data','numerical_correction_unavailable')[int(counts[1,j])],
            corrected_bins=int(counts[2,j]),
            static_dielectric_real=number(coeff[0,j]) if fitted else None,
            epsilon_real_omega2=number(coeff[1,j]) if fitted else None,
            sigma_imaginary_linear_slope=number(stats[3,j]) if fitted else None,
            sigma_real_omega2=number(coeff[2,j]) if fitted and absorption=='analytic' else None,
            sigma_real_omega4=number(coeff[3,j]) if fitted and absorption=='analytic' else None,
            fit_epsilon_real_rms=number(stats[0,j]) if fitted else None,
            max_epsilon_real_change=number(stats[1,j]),max_sigma_real_change=number(stats[2,j])))
    report=dict(schema='salmonpy.optical-a4/1',source_observable=name,reference='doi:10.1103/PhysRevB.85.045134 Appendix A4 and following text',
        system_assumption='ordinary zero-DC gapped optical response; no Drude/free-carrier/Hall contribution',
        assumption_verified=False,units=units,fit_range=bounds.tolist(),cutoff=edge,
        absorption=absorption,components=components,warnings=notices,
        method='native even quadratic Re epsilon; Im sigma=-omega*(Re epsilon-I)/(4*pi)',
        dc_scope='explicit analytic model extrapolation' if absorption=='analytic' else 'original DC validity retained; static real fit reported separately',
        absorption_scope='extra analytic insulating approximation' if absorption=='analytic' else 'raw Re sigma and Im epsilon preserved',
        numerical_scope='observation-only; no solver feedback, global KK reconstruction or scissor shift')
    base=json.loads(json.dumps(dict(observation.config),allow_nan=False))
    base.update(units=units,current_convention='charge',field_kind='total')
    correction=json.loads(json.dumps(report,allow_nan=False))
    def frozen(x):return np.frombuffer(x.tobytes(),dtype=x.dtype).reshape(x.shape)
    def unpack(x):return frozen(np.asarray(x.T).reshape(shape))
    w=frozen(omega)
    def pair(s,e,b,corrected):
        sv=unpack((b&1)!=0);ev=unpack((b&2)!=0)
        sc=dict(base,formula='sigma=omega*(epsilon-I)/(4*pi*i)',observable='optical_conductivity' if channels==1 else 'optical_conductivity_tensor',value_unit='1/atomic_time',spectrum_unit='1/atomic_time')
        ec=dict(base,formula='epsilon=I+4*pi*i*sigma/omega',observable='dielectric_function' if channels==1 else 'dielectric_tensor',value_unit='dimensionless',spectrum_unit='dimensionless')
        if corrected:
            sc['long_wavelength_correction']=json.loads(json.dumps(correction))
            ec['long_wavelength_correction']=json.loads(json.dumps(correction))
        return (OpticalSpectrum(w,unpack(s),sv,MappingProxyType(sc),unpack(e),ev),
                OpticalSpectrum(w,unpack(e),ev,MappingProxyType(ec)))
    raw_s,raw_e=pair(raw_sigma,raw_eps,raw_bits,False)
    corrected_s,corrected_e=pair(sigma,eps,bits,True)
    return LongWavelengthCorrection(raw_s,raw_e,corrected_s,corrected_e,unpack(applied!=0),MappingProxyType(correction))


@dataclass(frozen=True, slots=True)
class OpticalSensitivity:
    spectra: tuple
    reference_index: int
    report: object


@dataclass(frozen=True, slots=True)
class KCurrentContributions:
    time: np.ndarray
    k_ids: np.ndarray
    current: np.ndarray
    total: np.ndarray
    spectrum: object
    report: object


def _observed_transform(time, data, frequencies, window):
    """Native requested-frequency transform; Python only arranges storage."""
    import ctypes as c
    import warnings
    from ._sbe_handles import default_api
    from .errors import RuntimeFailure
    windows={'none':0,'hann':1,'cubic':2}
    if window not in windows:raise InputError('window must be none, hann or cubic')
    t=_times(time);x=_array(data,'data');w=_array(frequencies,'frequencies')
    if x.ndim<1 or len(x)!=len(t) or w.ndim!=1 or not len(w):
        raise InputError('data must start with the time axis and frequencies must be a nonempty vector')
    shape=x.shape[1:];channels=int(np.prod(shape)) if shape else 1
    if max(len(t),channels,len(w))>2**31-1:raise InputError('observation dimensions exceed native ABI')
    if len(t)*channels*len(w)>100_000_000:
        warnings.warn('explicit-frequency observation transform has large N*channels*frequencies cost; select fewer frequencies if needed',RuntimeWarning,stacklevel=2)
    values=np.asfortranarray(x.reshape(len(t),channels).T)
    out=np.empty((channels,len(w)),dtype=np.complex128,order='F')
    valid=np.empty(len(w),dtype=np.int32);stats=np.empty(4)
    try:fn=default_api().library.sp_observation_spectrum
    except AttributeError as exc:raise RuntimeFailure('native trajectory observation ABI missing; explicitly rebuild or install matching kernels') from exc
    fn.argtypes=[c.c_int]*3+[c.c_void_p]*3+[c.c_int]+[c.c_void_p]*3;fn.restype=c.c_int
    if fn(len(t),channels,len(w),t.ctypes.data,values.ctypes.data,w.ctypes.data,windows[window],out.ctypes.data,valid.ctypes.data,stats.ctypes.data):
        raise InputError('native observation transform requires finite uniform times, ordered nonnegative frequencies and finite data')
    return out.T.reshape((len(w),)+shape).copy(),valid!=0,stats


def optical_sensitivity(time,current,electric_field,*,frequencies,sample_counts,
                        current_convention,field_kind,units,windows=('none','cubic'),drive_rtol=1e-8):
    """Compare prefixes/windows on explicit COMMON angular frequencies, no solver.

    Longest prefix and first requested window are the descriptive reference.
    Prefixes keep the same initial time; no interpolation, zero padding, k scan
    or physical convergence certificate is implied. Numerical transforms and
    optical diagnostics are native Fortran; selected-frequency cost is O(N*F).
    """
    import ctypes as c
    import operator
    from time import perf_counter
    from .ld import _positive
    from ._sbe_handles import default_api
    if current_convention!='charge' or field_kind!='total' or units!='atomic_gaussian':
        raise InputError('declare charge current, total E and atomic_gaussian units')
    t=_times(time);j=_array(current,'current');e=_array(electric_field,'electric_field')
    if j.shape!=(len(t),) or e.shape!=j.shape:raise InputError('select one current and total-E component with shape (N,)')
    try:
        counts=tuple(operator.index(v) for v in sample_counts)
    except (TypeError,ValueError) as exc:raise InputError('sample_counts must contain integers') from exc
    if not counts or any(v<2 or v>len(t) for v in counts):raise InputError('sample counts must lie in [2, available samples]')
    if isinstance(windows,str):raise InputError('windows must be a sequence of window names')
    windows=tuple(windows)
    if not windows:raise InputError('at least one window required')
    threshold=_positive(drive_rtol,'drive_rtol',zero=True)
    if threshold>=1:raise InputError('drive_rtol must be below one')
    def freeze(x):return np.frombuffer(x.tobytes(),dtype=x.dtype).reshape(x.shape)
    w=freeze(_array(frequencies,'frequencies'))
    results=[];cases=[]
    for n in counts:
        for window in windows:
            started=perf_counter()
            transformed,usable,stats=_observed_transform(t[:n],np.column_stack((j[:n],e[:n])),w,window)
            transform_seconds=perf_counter()-started
            transformed[~usable]=0j
            js=np.ascontiguousarray(transformed[:,0]);es=np.ascontiguousarray(transformed[:,1])
            sigma=np.empty(len(w),complex);eps=np.empty_like(sigma);mask=np.empty(len(w),np.int32)
            fn=default_api().library.sp_optical_ratio
            fn.argtypes=[c.c_int]+[c.c_void_p]*3+[c.c_double]+[c.c_void_p]*3;fn.restype=c.c_int
            if fn(len(w),w.ctypes.data,js.ctypes.data,es.ctypes.data,threshold,sigma.ctypes.data,eps.ctypes.data,mask.ctypes.data):
                raise InputError('native optical conversion rejected transformed trajectory')
            cfg=dict(observable='optical_conductivity',units=units,current_convention='charge',field_kind='total',
                     sign='exp(+i omega t)',formula='sigma=j/E',value_unit='1/atomic_time',samples=n,window=window,
                     dt=float(stats[0]),duration=float(stats[1]),start=float(t[0]),drive_rtol=threshold,
                     nominal_frequency_resolution=float(stats[2]),nyquist=float(stats[3]),
                     transform='native dt*sum; explicit shared frequencies; no interpolation or window-gain correction')
            results.append(OpticalSpectrum(w,freeze(sigma),freeze(mask!=0),MappingProxyType(cfg),freeze(eps),freeze(mask==3)))
            cases.append(dict(samples=n,window=window,transform_seconds=transform_seconds,transform_work_items=n*2*len(w),duration=float(stats[1]),nominal_frequency_resolution=float(stats[2]),
                              nyquist=float(stats[3]),masked_frequencies=int(np.count_nonzero(mask==0))))
    ref=counts.index(max(counts))*len(windows)
    for case,result in zip(cases,results):
        a=np.ascontiguousarray(results[ref].values);b=np.ascontiguousarray(result.values)
        ma=np.ascontiguousarray(results[ref].valid,dtype=np.int32);mb=np.ascontiguousarray(result.valid,dtype=np.int32)
        stats=np.empty(5);fn=default_api().library.sp_observation_difference
        fn.argtypes=[c.c_int]+[c.c_void_p]*6;fn.restype=c.c_int
        if fn(len(w),w.ctypes.data,a.ctypes.data,b.ctypes.data,ma.ctypes.data,mb.ctypes.data,stats.ctypes.data):
            raise InputError('native stored-spectrum comparison rejected data')
        case.update(compared_frequencies=int(stats[0]),relative_rms_difference=float(stats[4]) if np.isfinite(stats[4]) else None)
    import hashlib
    digest=hashlib.sha256()
    for array in (t,j,e,w):digest.update(np.asarray(array,dtype='<f8').tobytes())
    return OpticalSensitivity(tuple(results),ref,MappingProxyType(dict(cases=cases,reference_index=ref,input_sha256=digest.hexdigest(),
        certified=False,scope='same stored trajectory only; window/time sensitivity is not grid/k/linearity convergence',
        frequency_grid='explicit common angular frequencies; bins above Nyquist or weak drive are masked')))


def k_current_contributions(decomposed,*,dt,current_convention,units,time0=0.,
                            frequencies=None,window='cubic',total_current=None):
    """Sum already weighted band/spin entries by k ID, with NO new k/occ weight.

    Input convention applies to decomposed entries only. Optional total_current
    must be simultaneous PHYSICAL CHARGE current, shape (N,3), in the same units.
    Partial entries remain partial; no missing band/k weight is reconstructed.
    Spectrum is charge current, not per-k dielectric functions. Relative RMS may
    exceed one when signed contributions cancel; it is not a probability.
    """
    import ctypes as c
    from .outputs import DecomposedCurrent
    from .ld import _positive
    from ._sbe_handles import default_api
    from .errors import RuntimeFailure
    if not isinstance(decomposed,DecomposedCurrent):raise InputError('pass read_decomposed_current output')
    if current_convention not in ('charge','matter') or units!='atomic_gaussian':
        raise InputError('declare source charge/matter current and atomic_gaussian units explicitly')
    if any(u not in (None,'a.u.','au','a.u','atomic units') for u in decomposed.units):
        raise InputError('convert decomposed current units to atomic_gaussian explicitly')
    delta=_positive(dt,'dt');start=float(_array(time0,'time0',()))
    raw=_array(decomposed.current,'decomposed current');ids=np.asarray(decomposed.indices)
    steps=np.asarray(decomposed.step)
    if steps.ndim!=1 or steps.dtype.kind not in 'iu' or ids.ndim!=2 or ids.shape[1]!=3 or ids.dtype.kind not in 'iu':
        raise InputError('integer step and (k,band,spin) index arrays required')
    if raw.shape!=(len(steps),len(ids),3) or len(steps)<2 or not len(ids):raise InputError('inconsistent decomposed array shapes')
    labels=[tuple(int(v) for v in row) for row in ids]
    if any(v<1 or v>2**63-1 for row in labels for v in row) or len(set(labels))!=len(labels):raise InputError('positive unique k/band/spin entries required')
    if any(int(s)<0 or int(s)>2**63-1 for s in steps):raise InputError('nonnegative signed-int64 steps required')
    kids=sorted({row[0] for row in labels});lookup={k:i+1 for i,k in enumerate(kids)}
    groups=np.array([lookup[row[0]] for row in labels],np.int32)
    if max(len(steps),len(labels),len(kids))>2**31-1:raise InputError('decomposition dimensions exceed native ABI')
    steps=np.ascontiguousarray(steps,dtype=np.int64);values=np.asfortranarray(raw.transpose(2,1,0))
    time=np.empty(len(steps));current=np.empty((3,len(kids),len(steps)),order='F')
    total=np.empty((3,len(steps)),order='F');stats=np.empty((2,len(kids)),order='F')
    try:fn=default_api().library.sp_k_current
    except AttributeError as exc:raise RuntimeFailure('native k-current observation ABI missing; explicitly rebuild matching kernels') from exc
    fn.argtypes=[c.c_int]*3+[c.c_void_p]*3+[c.c_double]*2+[c.c_int]+[c.c_void_p]*4;fn.restype=c.c_int
    if fn(len(steps),len(ids),len(kids),steps.ctypes.data,groups.ctypes.data,values.ctypes.data,delta,start,
          1 if current_convention=='charge' else -1,time.ctypes.data,current.ctypes.data,total.ctypes.data,stats.ctypes.data):
        raise InputError('native k-current reduction rejected nonfinite data, invalid times or numerical overflow')
    def freeze(x):return np.frombuffer(x.tobytes(),dtype=x.dtype).reshape(x.shape)
    current=freeze(current.transpose(2,1,0).copy());total=freeze(total.T.copy());time=freeze(time)
    comparison=None
    if total_current is not None:
        ref=np.ascontiguousarray(_array(total_current,'total_current',total.shape));candidate=np.ascontiguousarray(total)
        if ref.size>2**31-1:raise InputError('total-current comparison exceeds native ABI')
        error=np.empty(4);accepted=np.empty(1,np.int32);cmp=default_api().library.sp_convergence_waveforms
        cmp.argtypes=[c.c_int]+[c.c_void_p]*2+[c.c_double]*2+[c.c_void_p]*2;cmp.restype=c.c_int
        if cmp(ref.size,candidate.ctypes.data,ref.ctypes.data,0.,0.,error.ctypes.data,accepted.ctypes.data):raise InputError('invalid total current comparison')
        comparison=dict(rms_difference=float(error[0]),relative_rms_error=float(error[1]) if error[2]>0 else None,
                        maximum_difference=float(error[3]),scope='same sampled charge-current convention; no missing entries inferred')
    spectrum_result=None
    if frequencies is not None:
        transformed,usable,transform_stats=_observed_transform(time,current,frequencies,window)
        spectrum_result=Spectrum(freeze(_array(frequencies,'frequencies')),freeze(transformed),MappingProxyType(dict(
            observable='weighted_k_charge_current_spectrum',window=window,sign='exp(+i omega t)',units=units,
            frequency_valid=usable.tolist(),dt=float(transform_stats[0]),normalization='dt*sum',
            axes=('frequency','k','xyz'),value_unit='charge-current density * atomic time')))
    components=[dict(k_id=k,rms=float(stats[0,i]),relative_rms_to_total=float(stats[1,i]) if np.isfinite(stats[1,i]) else None) for i,k in enumerate(kids)]
    import hashlib
    digest=hashlib.sha256()
    for array in (steps,np.asarray(ids,dtype='<i8'),np.asarray(raw,dtype='<f8')):digest.update(array.tobytes())
    return KCurrentContributions(time,freeze(np.asarray(kids,dtype=np.int64)),current,total,spectrum_result,MappingProxyType(dict(input_sha256=digest.hexdigest(),
        components=components,total_current_comparison=comparison,weights='already applied; never reweighted',
        source_convention=current_convention,output_convention='charge',units=units,dt=delta,time0=start,
        entries=len(labels),certified=False,scope='only supplied entries; no proof of k/band completeness or k-sampling convergence')))


@dataclass(frozen=True, slots=True)
class PortSpectrum:
    """Directional flux ratios and unwindowed time-domain fluences."""
    omega: np.ndarray
    R: np.ndarray
    T: np.ndarray
    valid: np.ndarray
    incident_fluence: float
    reflected_fluence: float
    transmitted_fluence: float
    config: object


def _energy_native(name, signature):
    import ctypes as c
    from ._sbe_handles import default_api
    from .errors import RuntimeFailure
    try:
        fn = getattr(default_api().library, name)
    except AttributeError as exc:
        raise RuntimeFailure('Energy diagnostic ABI missing; install a matching wheel or explicitly rebuild kernels') from exc
    fn.argtypes = signature
    fn.restype = c.c_int
    return fn


def absorbed_energy_density(time, electric_field, current, *, units, current_convention,
                            field_kind):
    """Integral of collocated total E dot charge j, independently per macro point.

    Arrays are (time, point, xyz); returns Hartree/bohr**3 in Gaussian AU.
    This is net work including reversible material storage, NOT heat alone.
    Negative work is retained. No solver state or observations are modified.
    """
    import ctypes as c
    if units != 'atomic_gaussian' or current_convention != 'charge' or field_kind != 'total':
        raise InputError('require atomic_gaussian, charge current j=-J, and total E')
    t = np.ascontiguousarray(_times(time), dtype=np.float64)
    e = np.ascontiguousarray(_array(electric_field,'electric_field'), dtype=np.float64)
    j = np.ascontiguousarray(_array(current,'current'), dtype=np.float64)
    if e.ndim != 3 or e.shape != j.shape or e.shape[0] != len(t) or e.shape[1] < 1 or e.shape[2] != 3:
        raise InputError('E and j must both have shape (time, macro_point, 3)')
    if max(len(t),e.shape[1]) > 2**31-1: raise InputError('observation dimensions exceed native ABI')
    work = np.empty(e.shape[1], dtype=np.float64)
    fn = _energy_native('sp_absorbed_work', [c.c_int,c.c_int]+[c.c_void_p]*4)
    if fn(len(t),e.shape[1],t.ctypes.data,e.ctypes.data,j.ctypes.data,work.ctypes.data):
        raise InputError('native material work rejected input or nonfinite arithmetic')
    return _readonly(work)


def reflection_transmission(time, incident_E, incident_H, reflected_E, reflected_H,
                            transmitted_E, transmitted_H, *, normal, units,
                            window='none', drive_rtol=1e-8):
    """R=-Re(Er x Hr*)·n / Pi and T=Re(Et x Ht*)·n / Pi.

    Explicit direction-separated plane-port fields, shape (N,3), equal port area.
    n points from incident to transmitted side. Match E/H time and position first;
    never pass staggered Yee samples, total fields as reflected fields, or A/J.
    Does not infer incident-wave subtraction or integrate a nonuniform surface.
    """
    import ctypes as c
    import warnings
    from .ld import _positive
    if units != 'atomic_gaussian': raise InputError('require atomic_gaussian E/H units')
    threshold = _positive(drive_rtol,'drive_rtol',zero=True)
    if threshold >= 1: raise InputError('drive_rtol must be less than one')
    t = np.ascontiguousarray(_times(time), dtype=np.float64)
    n = np.ascontiguousarray(_array(normal,'normal'), dtype=np.float64)
    if n.shape != (3,): raise InputError('normal must have shape (3,)')
    arrays = [np.ascontiguousarray(_array(v,'port field'), dtype=np.float64)
              for v in (incident_E,incident_H,reflected_E,reflected_H,transmitted_E,transmitted_H)]
    if any(a.shape != (len(t),3) for a in arrays): raise InputError('port fields must have shape (N,3)')
    if len(t) > 2**31-1: raise InputError('observation dimensions exceed native ABI')
    transforms = [spectrum(t,a,window=window) for a in arrays]
    values = [np.ascontiguousarray(a.values, dtype=np.complex128) for a in transforms]
    size = len(transforms[0].omega)
    R = np.empty(size); T = np.empty(size); valid = np.empty(size,dtype=np.int32)
    fn = _energy_native('sp_port_spectra',[c.c_int]+[c.c_void_p]*7+[c.c_double]+[c.c_void_p]*3)
    if fn(size,*(a.ctypes.data for a in values),n.ctypes.data,threshold,R.ctypes.data,T.ctypes.data,valid.ctypes.data):
        raise InputError('native spectral flux rejected input or nonfinite arithmetic')
    fn = _energy_native('sp_port_fluence',[c.c_int]+[c.c_void_p]*5)
    fluences = []
    for E,H in zip(arrays[::2],arrays[1::2]):
        value = c.c_double()
        if fn(len(t),t.ctypes.data,E.ctypes.data,H.ctypes.data,n.ctypes.data,c.byref(value)):
            raise InputError('native time-domain fluence rejected input or nonfinite arithmetic')
        fluences.append(value.value)
    if not np.any(valid): warnings.warn('Incident spectral flux is unresolved; R/T are invalid',UserWarning,stacklevel=2)
    if np.any(R[valid!=0] < 0) or np.any(T[valid!=0] < 0):
        warnings.warn('Negative directional R/T: inspect port decomposition and normal',UserWarning,stacklevel=2)
    config = dict(transforms[0].config, observable='reflection_transmission',units=units,
                  drive_rtol=threshold,normal=tuple(n),fluence_unit='hartree/bohr^2',
                  fluence_window='none',ratio_normalization='directional E cross conjugate(H)',
                  scope='equal-area uniform plane ports; complete time trace and separated incident/reflected waves required')
    return PortSpectrum(transforms[0].omega,_readonly(R),_readonly(T),np.frombuffer((valid!=0).tobytes(),dtype=np.bool_),
                        fluences[0],-fluences[1],fluences[2],MappingProxyType(config))
