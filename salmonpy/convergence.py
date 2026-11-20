"""Explicit SALMON grid/dt campaigns; no native execution on import or planning.

Norm conservation is the primary dt diagnostic. It does not certify phase or
response accuracy; optional adjacent-current comparisons are reported separately.
"""
import ctypes as c
from functools import lru_cache
from copy import deepcopy
import math
import warnings
import numpy as np
from .errors import InputError, NumericalError, RuntimeFailure
from .structure import grid_from_spacing

__all__ = ['plan_grid_scan','run_grid_scan','plan_dt_scan','run_dt_scan','run_current_grid_scan']
_AU_TIME_FS = .02418884326505  # SALMON input convention.


@lru_cache(maxsize=1)
def _kernel():
    # Existing checksum-verified bundled library; no compilation or fallback.
    from ._sbe_handles import default_api
    lib=default_api().library
    try:
        lib.sp_convergence_check.argtypes=[c.c_int,c.c_void_p,c.c_int]
        lib.sp_convergence_norm.argtypes=[c.c_int,c.c_void_p,c.c_double,c.c_double,c.c_double,c.c_void_p,c.POINTER(c.c_int)]
        lib.sp_convergence_compare.argtypes=[c.c_int,c.c_void_p,c.c_void_p,c.c_double,c.c_double,c.c_void_p,c.POINTER(c.c_int)]
        lib.sp_convergence_waveforms.argtypes=lib.sp_convergence_compare.argtypes
        lib.sp_convergence_waveforms.restype=c.c_int
        for name in ('sp_convergence_check','sp_convergence_norm','sp_convergence_compare'):
            getattr(lib,name).restype=c.c_int
    except AttributeError as exc:
        raise RuntimeFailure('native convergence ABI unavailable; install a matching wheel or explicitly rebuild kernels') from exc
    return lib


def _checked(values,positive=False):
    array=np.ascontiguousarray(values,dtype=np.float64)
    if _kernel().sp_convergence_check(array.size,array.ctypes.data,int(positive)):
        raise NumericalError('nonfinite or invalid numerical sample')
    return array


def _comparison(coarse,fine,atol,rtol):
    a=np.ascontiguousarray(coarse,dtype=np.float64);b=np.ascontiguousarray(fine,dtype=np.float64)
    if a.shape!=b.shape:raise InputError('comparison sampling shapes differ')
    stats=np.empty(2);accepted=c.c_int()
    if _kernel().sp_convergence_compare(a.size,a.ctypes.data,b.ctypes.data,atol,rtol,stats.ctypes.data,c.byref(accepted)):
        raise NumericalError('native comparison rejected nonfinite arithmetic')
    return {'absolute_error':float(stats[0]),'relative_error':float(stats[1]),'within_tolerance':bool(accepted.value)}


def _norm_report(values,reference,atol,rtol):
    array=np.ascontiguousarray(values,dtype=np.float64);stats=np.empty(6);accepted=c.c_int()
    if _kernel().sp_convergence_norm(array.size,array.ctypes.data,reference,atol,rtol,stats.ctypes.data,c.byref(accepted)):
        raise NumericalError('native norm diagnostics rejected nonfinite arithmetic')
    return dict(zip(('norm_max_loss','norm_max_gain','norm_max_loss_fraction','norm_max_gain_fraction',
                     'norm_max_absolute_error','norm_max_relative_error'),map(float,stats)),
                norm_within_tolerance=bool(accepted.value))


def _number(value, name, *, positive=False):
    if isinstance(value,(bool,np.bool_)):
        raise InputError(name+' must be a finite number')
    try:
        result=float(value)
    except (TypeError,ValueError) as exc:
        raise InputError(name+' must be a finite number') from exc
    if not math.isfinite(result) or (result<=0 if positive else result<0):
        raise InputError(name+' must be finite and '+('positive' if positive else 'nonnegative'))
    return result


def _identity(material_id,metadata):
    if not isinstance(material_id,str) or not material_id.strip():
        raise InputError('nonempty material_id required')
    return {'material_id':material_id,'metadata':deepcopy(metadata or {})}


def plan_grid_scan(cell_angstrom, *, radii, material_id, lower_spacing=.1,
                   factor=.8, minimum=8, multiples=(1,1,1), metadata=None):
    """PP radius/3 to lower_spacing in angstrom, with explicit radius provenance.

    radii maps species to radius_angstrom, definition and pp_sha256. The caller
    selects the physical cutoff channel/radius; radial-table extents are not cutoffs.
    Finer than0.1 angstrom is allowed. Fixed low-k settings belong in metadata.
    """
    identity=_identity(material_id,metadata)
    if not isinstance(radii,dict) or not radii:
        raise InputError('explicit species radii required')
    records=deepcopy(radii)
    for species,record in records.items():
        if not isinstance(species,str) or not isinstance(record,dict):
            raise InputError('species radius records required')
        record['radius_angstrom']=_number(record.get('radius_angstrom'),'radius',positive=True)
        digest=record.get('pp_sha256','')
        if not isinstance(record.get('definition'),str) or not record['definition'].strip():
            raise InputError('explicit radius definition required; do not infer a table extent')
        if not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
            raise InputError('pseudopotential SHA256 required for each species')
    lower=_number(lower_spacing,'lower_spacing',positive=True)
    factor=_number(factor,'factor',positive=True)
    if factor>=1:raise InputError('refinement factor must be below1')
    initial=min(v['radius_angstrom'] for v in records.values())/3
    widths=[initial]
    while widths[-1]>lower:
        next_width=max(lower,widths[-1]*factor)
        if next_width==widths[-1]:raise InputError('refinement factor cannot reduce spacing')
        widths.append(next_width)
    trials=[];seen=set()
    for width in widths:
        grid=grid_from_spacing(cell_angstrom,spacing=width,unit='angstrom',minimum=minimum,multiples=multiples)
        if grid.counts not in seen:
            trials.append({'index':len(trials),'grid':grid.as_dict(),**deepcopy(identity)})
            seen.add(grid.counts)
    return {'schema':'salmonpy.grid-scan/1',**identity,'radii':records,
            'cell_angstrom':np.asarray(cell_angstrom,dtype=float).tolist(),
            'initial_spacing_angstrom':initial,'lower_spacing_angstrom':lower,'trials':trials}


def _tolerances(atol,rtol):
    return _number(atol,'atol'),_number(rtol,'rtol')


def run_grid_scan(plan, evaluate, *, atol, rtol, stop_when_converged=True):
    """Evaluate a fixed scalar observable (e.g. GS Hartree/atom) on each grid.

    evaluate(trial) returns value and scf_converged; nonconvergence warns and
    retains provisional evidence. Units, k points and SCF conditions are caller-owned.
    """
    atol,rtol=_tolerances(atol,rtol)
    if plan.get('schema')!='salmonpy.grid-scan/1':raise InputError('grid scan plan required')
    if type(stop_when_converged) is not bool:raise InputError('stop_when_converged must be boolean')
    _kernel()
    rows=[];selected=None;provisional=False;selected_provisional=False
    for trial in deepcopy(plan['trials']):
        row={**trial,'status':'complete','comparison':None}
        try:
            result=evaluate(deepcopy(trial));value=float(result['value'])
            _checked([value])
            if type(result.get('scf_converged')) is not bool:raise InputError('explicit scf_converged boolean required')
            row.update(value=value,scf_converged=result['scf_converged'],result_metadata=deepcopy(result.get('metadata',{})))
            if not row['scf_converged']:
                provisional=True;warnings.warn('SCF not converged; grid evidence is provisional',UserWarning,stacklevel=2)
            if rows and rows[-1]['status']=='complete':
                coarse=rows[-1]
                row['comparison']={**_comparison([coarse['value']],[value],atol,rtol),
                    'coarse_index':coarse['index'],'provisional':not (coarse['scf_converged'] and row['scf_converged'])}
                if row['comparison']['within_tolerance'] and (selected is None or
                        (selected_provisional and not row['comparison']['provisional'])):
                    selected=max(coarse['grid']['actual_spacing_angstrom'])
                    selected_provisional=row['comparison']['provisional']
        except (NumericalError,FloatingPointError) as exc:
            row.update(status='numerical_failure',value=None,error=str(exc))
        except (RuntimeFailure,OSError) as exc:
            row.update(status='runtime_failure',value=None,error=str(exc))
        rows.append(row)
        if (stop_when_converged and row.get('comparison')
                and row['comparison']['within_tolerance'] and not row['comparison']['provisional']):
            break
    return {'schema':'salmonpy.grid-scan-report/1','plan':deepcopy(plan),'material_id':plan['material_id'],
            'atol':atol,'rtol':rtol,'trials':rows,'selected_spacing_angstrom':selected,
            'provisional':provisional,'kpoint_convergence_established':False,
            'stopped_early':len(rows)<len(plan['trials']),'selected_provisional':selected_provisional}


def plan_dt_scan(*, duration_fs=3., levels=7, initial_dt=.2, material_id,
                 impulse=(0.,0.,1e-4), expected_electrons=None, metadata=None):
    """Halved dt, same physical duration and same A/c impulse on every trial.

    Initial A/c is zero. Every endpoint step receives the same post-impulse A/c;
    the integrated E is -impulse, independent of dt. End time is rounded up to a
    coarse step and then identical across trials. All quantities except fs are au.
    """
    identity=_identity(material_id,metadata)
    duration=_number(duration_fs,'duration_fs',positive=True)
    dt=_number(initial_dt,'initial_dt',positive=True)
    if type(levels) is not int or levels<1:raise InputError('positive integer levels required')
    A=np.asarray(impulse,dtype=float)
    if A.shape!=(3,) or not np.isfinite(A).all():raise InputError('finite three-component impulse required')
    if expected_electrons is not None:expected_electrons=_number(expected_electrons,'expected_electrons',positive=True)
    steps=math.ceil(duration/_AU_TIME_FS/dt)
    trials=[]
    for i in range(levels):
        value=dt/(2**i)
        if value==0:raise InputError('dt underflow in refinement plan')
        trials.append({'index':i,'dt_au':value,'steps':steps*(2**i)})
    return {'schema':'salmonpy.dt-scan/1',**identity,'trials':trials,
            'requested_duration_fs':duration,'duration_au':steps*dt,
            'duration_fs':steps*dt*_AU_TIME_FS,'impulse_A_over_c':A.tolist(),
            'expected_electrons':expected_electrons,'time_unit_fs':_AU_TIME_FS}


def run_dt_scan(plan, new_state, *, electron_number=None, atol=0., rtol=1e-5,
                current_atol=None,current_rtol=None):
    """Stream independent states and stop a nonfinite trial at its first bad step.

    new_state(dt) starts at t=0,A=0. electron_number(state) must measure the actual
    occupied-wavefunction/density norm, not return the configured nelec. It is
    sampled at t=0 and after every step. Norm and optional current tolerances are
    independent; neither is a material-independent physical accuracy certificate.
    Default relative norm tolerance is1e-5. Exceeding it records a result and
    continues the trial; only numerical breakdown ends a trial immediately.
    """
    atol,rtol=_tolerances(atol,rtol)
    compare=current_atol is not None or current_rtol is not None
    if compare:
        if current_atol is None or current_rtol is None:raise InputError('both current tolerances required')
        current_atol,current_rtol=_tolerances(current_atol,current_rtol)
    if plan.get('schema')!='salmonpy.dt-scan/1':raise InputError('dt scan plan required')
    _kernel()
    A=np.asarray(plan['impulse_A_over_c'],dtype=float);A.setflags(write=False)
    rows=[];nonbreaking=None;acceptable=None;current_choice=None
    def measured(state):
        value=float(state.electron_number if electron_number is None else electron_number(state))
        _checked([value],positive=True)
        return value
    for trial in deepcopy(plan['trials']):
        state=None;row={**trial,'status':'complete','current':None,'electron_number':None,'cleanup_error':None,'completed_steps':0,'attempted_step':0,'current_comparison':None}
        currents=[];norms=[]
        try:
            state=new_state(trial['dt_au']);norms=[measured(state)]
            reference=norms[0] if plan['expected_electrons'] is None else plan['expected_electrons']
            initial=np.asarray(state.current,dtype=float)
            if initial.shape!=(3,):raise RuntimeFailure('three-component current required')
            _checked(initial)
            currents=[initial.copy()]
            for step in range(1,trial['steps']+1):
                row['attempted_step']=step
                value=np.asarray(state.step(A),dtype=float)
                if value.shape!=(3,):raise RuntimeFailure('three-component current required')
                _checked(value)
                norms.append(measured(state));currents.append(value.copy());row['completed_steps']=step
            row.update(current=np.asarray(currents).tolist(),electron_number=norms,reference_electrons=reference,
                       **_norm_report(norms,reference,atol,rtol))
            nonbreaking=trial['dt_au'] if nonbreaking is None else max(nonbreaking,trial['dt_au'])
            if row['norm_within_tolerance']:acceptable=trial['dt_au'] if acceptable is None else max(acceptable,trial['dt_au'])
            if compare and rows and rows[-1]['status']=='complete':
                coarse=rows[-1];fine=np.asarray(currents)[::2];previous=np.asarray(coarse['current'])
                comparison=_comparison(previous,fine,current_atol,current_rtol)
                within=comparison['within_tolerance']
                row['current_comparison']={'coarse_index':coarse['index'],
                    'max_absolute_error':comparison['absolute_error'],
                    'relative_error':comparison['relative_error'],'within_tolerance':within}
                if within and current_choice is None:current_choice=coarse['dt_au']
        except (NumericalError,FloatingPointError) as exc:
            row.update(status='numerical_failure',error=str(exc),failed_step=row['attempted_step'],current=None,electron_number=None)
        except (RuntimeFailure,OSError) as exc:
            row.update(status='runtime_failure',error=str(exc),current=None,electron_number=None)
        finally:
            if state is not None:
                try:state.close()
                except Exception as exc:row['cleanup_error']=str(exc)
        rows.append(row)
    return {'schema':'salmonpy.dt-scan-report/1','plan':deepcopy(plan),'material_id':plan['material_id'],
            'metadata':deepcopy(plan['metadata']),'duration_au':plan['duration_au'],
            'impulse_A_over_c':plan['impulse_A_over_c'].copy(),'atol':atol,'rtol':rtol,
            'current_atol':current_atol,'current_rtol':current_rtol,'trials':rows,
            'nonbreaking_dt_au':nonbreaking,'norm_acceptable_dt_au':acceptable,
            'current_converged_dt_au':current_choice,'physical_accuracy_certified':False}


def _waveform_comparison(coarse, fine, atol, rtol):
    a=np.ascontiguousarray(coarse,dtype=np.float64)
    b=np.ascontiguousarray(fine,dtype=np.float64)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or not len(a):
        raise InputError('current traces must share shape (time,3)')
    stats=np.empty(4);accepted=c.c_int()
    code=_kernel().sp_convergence_waveforms(a.size,a.ctypes.data,b.ctypes.data,
                                           atol,rtol,stats.ctypes.data,c.byref(accepted))
    if code:raise NumericalError('native current comparison rejected nonfinite arithmetic')
    return {'absolute_rms_error':float(stats[0]),
            'relative_rms_error':float(stats[1]) if stats[2]>0 else None,
            'reference_rms':float(stats[2]),'maximum_absolute_error':float(stats[3]),
            'within_tolerance':bool(accepted.value)}


def run_current_grid_scan(plan, evaluate, *, rtol=.05, atol=0., stop_when_converged=True):
    """Screen grid currents at norm-qualified dt; default relative RMS is5%.

    This screening default does not certify dielectric or high-harmonic spectra.
    Select tighter observable-specific tolerances and verify spectral convergence
    when the scientific target requires them.

    evaluate(trial) returns scf_converged and dt_report from run_dt_scan. Hold
    material, excitation, k sampling and physical duration fixed. Nested endpoint
    grids are aligned by subsampling; no interpolation changes the observables.
    """
    atol,rtol=_tolerances(atol,rtol)
    if plan.get('schema')!='salmonpy.grid-scan/1':raise InputError('grid scan plan required')
    if type(stop_when_converged) is not bool:raise InputError('stop_when_converged must be boolean')
    _kernel()
    rows=[];selected=None;selected_provisional=False;provisional=False
    for trial in deepcopy(plan['trials']):
        row={**trial,'status':'complete','comparison':None}
        try:
            value=evaluate(deepcopy(trial))
            if type(value.get('scf_converged')) is not bool:raise InputError('explicit scf_converged boolean required')
            report=deepcopy(value['dt_report'])
            if report['material_id']!=plan['material_id']:raise InputError('grid and dt material IDs differ')
            row.update(scf_converged=value['scf_converged'],dt_report=report)
            if not row['scf_converged']:
                provisional=True
                warnings.warn('SCF not converged; current-grid evidence is provisional',UserWarning,stacklevel=2)
            dt=report.get('norm_acceptable_dt_au')
            eligible=[t for t in report['trials'] if t['dt_au']==dt and t['status']=='complete' and t.get('norm_within_tolerance')]
            if not eligible:
                row['status']='dt_unresolved'
            else:
                row['selected_dt_au']=dt
                trace=np.asarray(eligible[0]['current'],dtype=float)
                if trace.ndim!=2 or trace.shape[1]!=3 or len(trace)<2:raise InputError('current trace must include initial and final three-component samples')
                if not math.isclose((len(trace)-1)*dt,report['duration_au'],rel_tol=1e-12):raise InputError('current samples do not span the declared duration')
                row['current']=_checked(trace).tolist()
                if rows and rows[-1]['status']=='complete':
                    coarse=rows[-1];previous=coarse['dt_report']
                    if (report['material_id']!=previous['material_id'] or
                        report['impulse_A_over_c']!=previous['impulse_A_over_c'] or
                        report['atol']!=previous['atol'] or report['rtol']!=previous['rtol'] or
                        report['duration_au']!=previous['duration_au']):
                        raise InputError('current grid comparison requires matching material, impulse, norm target and duration')
                    a=np.asarray(coarse['current'],dtype=float);b=np.asarray(row['current'],dtype=float)
                    if coarse['selected_dt_au']>=dt:
                        ratio=coarse['selected_dt_au']/dt;stride=round(ratio)
                        if stride<1 or not math.isclose(ratio,stride,rel_tol=1e-12):raise InputError('nested current sample times required')
                        b=b[::stride]
                    else:
                        ratio=dt/coarse['selected_dt_au'];stride=round(ratio)
                        if stride<1 or not math.isclose(ratio,stride,rel_tol=1e-12):raise InputError('nested current sample times required')
                        a=a[::stride]
                    row['comparison']={**_waveform_comparison(a,b,atol,rtol),
                        'coarse_index':coarse['index'],'sample_count':len(a),
                        'provisional':not (coarse['scf_converged'] and row['scf_converged'])}
                    comparison=row['comparison']
                    if comparison['within_tolerance'] and (selected is None or
                        (selected_provisional and not comparison['provisional'])):
                        selected=max(coarse['grid']['actual_spacing_angstrom'])
                        selected_provisional=comparison['provisional']
        except (NumericalError,FloatingPointError) as exc:
            row.update(status='numerical_failure',error=str(exc))
        except (RuntimeFailure,OSError) as exc:
            row.update(status='runtime_failure',error=str(exc))
        rows.append(row)
        if (stop_when_converged and row.get('comparison') and
            row['comparison']['within_tolerance'] and not row['comparison']['provisional']):break
    return {'schema':'salmonpy.current-grid-report/1','plan':deepcopy(plan),'material_id':plan['material_id'],
            'metric':'current relative RMS at common physical times','rtol':rtol,'atol':atol,
            'trials':rows,'selected_spacing_angstrom':selected,'selected_provisional':selected_provisional,
            'provisional':provisional,'stopped_early':len(rows)<len(plan['trials']),
            'kpoint_convergence_established':False}
