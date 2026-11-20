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
    """Single-polarization response; invalid frequencies carry NaN and valid=False."""
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
    if not isinstance(conductivity,OpticalSpectrum) or conductivity.config.get('observable')!='optical_conductivity':
        raise InputError('pass the result of optical_conductivity')
    if conductivity._dielectric_values is None or conductivity._dielectric_valid is None:
        raise InputError('conductivity lacks native dielectric conversion data')
    cfg=dict(conductivity.config)
    values=conductivity._dielectric_values;valid=conductivity._dielectric_valid
    cfg.update(observable='dielectric_function',formula='epsilon=1+4*pi*i*sigma/omega',value_unit='dimensionless',spectrum_unit='dimensionless')
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
                reference_peak_abs_E=reference.config['peak_abs_E'],
                candidate_peak_abs_E=candidate.config['peak_abs_E'],
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
