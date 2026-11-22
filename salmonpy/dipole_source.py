"""One-way charge-dipole sources, in atomic units, for native E/H fields."""
from dataclasses import dataclass
import numpy as np
from .errors import InputError, NumericalError, RuntimeFailure
from .field_input import _real

__all__ = ['DipoleMap', 'drive_native_dipoles', 'MolecularDrive', 'select_molecular_field']


@dataclass(frozen=True, slots=True, eq=False)
class DipoleMap:
    """Map molecular dipoles (e*bohr) to macro-cell polarization (e/bohr**2).

    weights[point,molecule] sums to one for each molecule. Volumes are physical
    FDTD cell volumes in bohr**3, not molecular electronic-box volumes. Dipoles
    must already be rotated into laboratory axes and include the intended count
    of molecules. The native field separately deposits macro currents on Yee edges.
    """
    molecule_ids: tuple
    point_ids: tuple
    weights: object
    cell_volumes: object

    def __post_init__(self):
        molecules = tuple(self.molecule_ids)
        points = tuple(self.point_ids)
        if (not molecules or any(type(x) is not str or not x.strip() for x in molecules)
                or len(set(molecules)) != len(molecules)):
            raise InputError('distinct nonempty molecule IDs required')
        if (not points or any(type(x) is not int or x < 0 for x in points)
                or len(set(points)) != len(points)):
            raise InputError('distinct nonnegative integer point IDs required')
        weights = _real(self.weights, (len(points), len(molecules)), 'weights')
        volumes = _real(self.cell_volumes, (len(points),), 'cell volumes')
        if np.any(weights < 0) or not np.allclose(weights.sum(axis=0), 1, rtol=0, atol=1.e-12):
            raise InputError('nonnegative weights must sum to one per molecule')
        if np.any(volumes <= 0):
            raise InputError('positive physical cell volumes required')
        for name, value in (('molecule_ids', molecules), ('point_ids', points),
                            ('weights', weights), ('cell_volumes', volumes)):
            object.__setattr__(self, name, value)

    def polarization(self, dipoles):
        p = _real(dipoles, (len(self.molecule_ids), 3), 'charge dipoles')
        return _mapped_dipoles(self,p,p,dt=1.,mode=0,name='polarization')

    def current(self, old_dipoles, new_dipoles, *, dt):
        """Interval-average charge current; timestamp is midpoint of the interval."""
        dt = float(_real(dt, (), 'dt'))
        if dt <= 0:
            raise InputError('positive dt required')
        old = _real(old_dipoles, (len(self.molecule_ids), 3), 'old dipoles')
        new = _real(new_dipoles, old.shape, 'new dipoles')
        return _mapped_dipoles(self,new,old,dt=dt,mode=1,name='dipole current')


def drive_native_dipoles(field, mapping, dipoles, *, dt):
    """Collectively drive a fresh NativeEHField with precomputed charge dipoles.

    dipoles has shape (intervals+1,molecules,3), at times n*dt starting at zero.
    The first interval must produce zero current. To play a nonstationary trace,
    explicitly repeat its initial sample, introducing a declared one-dt delay.
    Call on every field-group rank with identical inputs. Field lifetime is owned
    by caller; backend failure is not a whole-trajectory rollback. One-way only.
    """
    from .coupled import NativeEHField
    if not isinstance(field, NativeEHField):
        raise InputError('NativeEHField required')
    comm = field.context.solver_comm
    error = None
    try:
        if not isinstance(mapping, DipoleMap):
            raise InputError('DipoleMap required')
        dt = float(_real(dt, (), 'dt'))
        if dt != field.dt or dt <= 0:
            raise InputError('dipole dt must equal field dt; resample explicitly')
        if field.closed or field.failed or field.step_index != 0:
            raise InputError('fresh usable field required')
        if mapping.point_ids != tuple(field.point_ids):
            raise InputError('mapping point order differs from native field')
        if not np.allclose(mapping.cell_volumes, field.manifest['cell_volume'], rtol=1.e-14, atol=0):
            raise InputError('mapping volumes differ from native field cells')
        raw = np.asarray(dipoles)
        if raw.ndim != 3 or raw.shape[0] < 2:
            raise InputError('at least two dipole time samples required')
        p = _real(raw, (len(raw), len(mapping.molecule_ids), 3), 'dipole trajectory')
        if len(p)-1 > field.capacity:
            raise InputError('dipole trajectory exceeds field capacity')
        currents = np.array([mapping.current(a,b,dt=dt) for a,b in zip(p[:-1],p[1:])])
        if np.any(currents[0]):
            raise InputError('first interval current must be zero; explicitly hold initial dipoles for one step')
        descriptor = (dt, mapping.molecule_ids, mapping.point_ids, mapping.weights.tobytes(),
                      mapping.cell_volumes.tobytes(), p.tobytes())
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    errors = comm.allgather(error)
    if any(errors):
        raise InputError(next(e for e in errors if e))
    # Agree the entire trajectory before advancing even the first field step.
    import hashlib
    key = (descriptor[:3], hashlib.sha256(b''.join(descriptor[3:])).hexdigest(), p.shape)
    keys = comm.allgather(key)
    if any(k != keys[0] for k in keys):
        raise InputError('dipole source differs across field ranks')
    electric = [field.E.copy()]
    for current in currents:
        field.advance(current)
        electric.append(field.E.copy())
    shape = (len(p), len(mapping.point_ids), 3)
    return {'E': _real(electric, shape, 'sampled E'),
            'E_time': _real(np.arange(len(p))*dt, (len(p),), 'E time'),
            'j': _real(currents, (len(p)-1, len(mapping.point_ids), 3), 'charge current'),
            'j_time': _real((np.arange(len(p)-1)+.5)*dt, (len(p)-1,), 'current time'),
            'point_ids': mapping.point_ids, 'molecule_ids': mapping.molecule_ids}


@dataclass(frozen=True, slots=True)
class MolecularDrive:
    """Immutable sampled drive and caller-declared selection provenance."""
    field: object
    metadata: object

    def __post_init__(self):
        from .field_input import FieldValue
        from .config import _copy_json, _freeze
        if not isinstance(self.field, FieldValue):
            raise InputError('FieldValue required')
        object.__setattr__(self, 'metadata', _freeze(_copy_json(self.metadata)))

    def record(self):
        """Return a detached JSON-compatible selection record."""
        from .config import _thaw
        return _thaw(self.metadata)


def select_molecular_field(total, *, molecule_id, policy, provenance,
                           reference=None, reference_provenance=None):
    """Select E explicitly; never start a solver or infer a local-field model.

    provided_drive preserves an already selected field with a required meaning.
    environment_scattered_only subtracts this molecule's declared own-source
    vacuum reference. It removes vacuum radiation reaction too. Provenance is
    caller-declared, not independent verification of Maxwell simulations.
    """
    from .field_input import FieldValue
    from .config import _copy_json, _identifier
    _identifier(molecule_id, 'molecule_id')

    def electric(sample):
        if (not isinstance(sample, FieldValue) or sample.unit != 'E:au'
                or sample.position_unit != 'bohr' or sample.sampling != 'point'):
            raise InputError('point-sampled E:au FieldValue at bohr positions required')

    def info(raw):
        result = _copy_json(raw)
        if type(result) is not dict:
            raise InputError('provenance must be a JSON object')
        return result

    electric(total)
    declaration = info(provenance)
    ref_info = None
    if policy == 'provided_drive':
        if reference is not None or reference_provenance is not None:
            raise InputError('provided_drive does not subtract a reference')
        meaning = declaration.get('meaning')
        if type(meaning) is not str or not meaning.strip():
            raise InputError('provided_drive requires a provenance meaning')
        result = total.value
    elif policy == 'environment_scattered_only':
        electric(reference)
        ref_info = info(reference_provenance)
        if (reference.time != total.time or not np.array_equal(reference.position, total.position)):
            raise InputError('reference time and position must match exactly; interpolate explicitly')
        if (ref_info.get('source_ids') != [molecule_id] or ref_info.get('medium') != 'vacuum'
                or ref_info.get('incident') is not False):
            raise InputError('reference must contain only this molecule in vacuum, with no incident field')
        for key in ('footprint_id', 'sample_map_id', 'own_source_history_id'):
            _identifier(declaration.get(key), key)
            if ref_info.get(key) != declaration[key]:
                raise InputError('reference differs in '+key)
        result = _selected_molecular_field(total.value,reference.value)
    else:
        raise InputError('explicit provided_drive or environment_scattered_only policy required')
    selected = FieldValue(result, total.time, total.position, 'E:au')
    metadata = {'schema':'salmonpy.molecular-drive/1', 'molecule_id':molecule_id,
                'policy':policy, 'time':total.time, 'position':total.position.tolist(),
                'unit':'E:au', 'position_unit':'bohr', 'sampling':'point',
                'input_E':total.value.tolist(), 'selected_E':selected.value.tolist(),
                'reference_E':None if reference is None else reference.value.tolist(),
                'provenance':declaration, 'reference_provenance':ref_info,
                'provenance_status':'caller-declared; not solver-certified'}
    return MolecularDrive(selected, metadata)


def _molecular_function(name,arguments):
    import ctypes as c
    from ._sbe_handles import default_api
    from .errors import RuntimeFailure
    try:
        library=default_api().library;probe=library.sp_molecular_protocol
        probe.argtypes=[];probe.restype=c.c_int
        if probe()!=1:raise RuntimeFailure('unsupported molecular development-kernel protocol')
        fn=getattr(library,name)
    except AttributeError as exc:raise RuntimeFailure('matching molecular test kernel required; explicitly build/install it') from exc
    fn.argtypes=arguments;fn.restype=c.c_int
    return fn


def _midpoint_oscillator(p,v,total_left,total_right,*,dt,omega,volume=1.,own_left=None,own_right=None):
    """Development-only unit-mass/charge test response; all numerics in Fortran.

    Returns new p,v, interval-average j and statistics [H0,H1,Wtotal,Wown,
    Wdrive,H1-H0-Wdrive,Wtotal-Wown-Wdrive,Wheld-Wdrive]. No input mutation.
    This is not an isolated-TDDFT propagator or a public two-way driver.
    """
    import ctypes as c
    vectors=[np.ascontiguousarray(_real(value,(3,),name)) for value,name in (
        (p,'p'),(v,'v'),(total_left,'total left'),(total_right,'total right'),
        (np.zeros(3) if own_left is None else own_left,'own left'),
        (np.zeros(3) if own_right is None else own_right,'own right'))]
    scalars=[float(_real(value,(),name)) for value,name in ((dt,'dt'),(omega,'omega'),(volume,'volume'))]
    outputs=[np.empty(3),np.empty(3),np.empty(3),np.empty(8)]
    fn=_molecular_function('sp_molecular_midpoint',[c.c_double]*3+[c.c_void_p]*10)
    status=fn(*scalars,*(a.ctypes.data for a in vectors+outputs))
    if status==1:raise InputError('invalid molecular midpoint test input')
    if status:raise NumericalError('molecular midpoint numerical failure; inputs unchanged')
    return tuple(outputs)


def _molecular_work(p,v,q,w,total_left,total_right,*,omega,own_left,own_right):
    """Actual accepted-state work diagnostics for the development fixture."""
    import ctypes as c
    arrays=[np.ascontiguousarray(_real(a,(3,),'molecular work vector'))
            for a in (p,v,q,w,total_left,total_right,own_left,own_right)]
    output=np.empty(8);fn=_molecular_function('sp_molecular_work',[c.c_double]+[c.c_void_p]*9)
    code=fn(float(omega),*(a.ctypes.data for a in arrays),output.ctypes.data)
    if code:raise NumericalError('molecular work diagnostic rejected nonfinite values')
    return output


def _molecular_summary(statistics):
    import ctypes as c
    a=np.ascontiguousarray(statistics,dtype=np.float64)
    if a.ndim!=2 or a.shape[1]!=8 or not len(a):raise InputError('nonempty molecular statistics[steps,8] required')
    output=np.empty(11);fn=_molecular_function('sp_molecular_summary',[c.c_int,c.c_void_p,c.c_void_p])
    if fn(len(a),a.ctypes.data,output.ctypes.data):raise NumericalError('invalid molecular summary')
    return output


def _selected_molecular_field(total,own):
    import ctypes as c
    total=np.ascontiguousarray(_real(total,(3,),'total E'))
    own=np.ascontiguousarray(_real(own,(3,),'own-source E'));out=np.empty(3)
    fn=_molecular_function('sp_molecular_select',[c.c_void_p]*3)
    if fn(total.ctypes.data,own.ctypes.data,out.ctypes.data):raise NumericalError('invalid molecular selected field')
    return out


def _mapped_dipoles(mapping,p,old,*,dt,mode,name):
    """Owned map arrays are immutable; result is published only after native success."""
    import ctypes as c
    p=np.ascontiguousarray(p);old=np.ascontiguousarray(old)
    result=np.empty((len(mapping.point_ids),3))
    fn=_molecular_function('sp_molecular_map',[c.c_int]*3+[c.c_double]+[c.c_void_p]*5)
    code=fn(mode,len(mapping.point_ids),len(mapping.molecule_ids),dt,
            mapping.weights.ctypes.data,mapping.cell_volumes.ctypes.data,
            p.ctypes.data,old.ctypes.data,result.ctypes.data)
    if code==1:raise InputError('invalid dipole map input')
    if code:raise NumericalError(name+' overflow')
    return _real(result,result.shape,name)


class _IsolatedDipoleTrial:
    """Private single-molecule interval-E transaction for the FDTD fixture.

    The isolated worker owns its length-gauge integration and applied-E history.
    Cancellation uses its supported replay restore, not direct electronic restore.
    Field backups and collective acceptance belong to the caller. No mutation of
    the molecule is allowed between propose and commit/cancel.
    """
    def __init__(self, state, mapping, *, sampling='midpoint', initial_field=(0.,0.,0.)):
        from .isolated_tddft import IsolatedTDDFTState
        if not isinstance(state, IsolatedTDDFTState):
            raise InputError('isolated E-driven TDDFT state required')
        if not isinstance(mapping, DipoleMap) or len(mapping.molecule_ids) != 1:
            raise InputError('single-molecule dipole map required')
        if sampling not in ('midpoint','endpoint'):
            raise InputError('explicit midpoint or endpoint molecular sampling required')
        self.sampling = sampling
        self._left = _real(initial_field,(3,),'initial selected E')
        self._accepted_right = None
        self.state = state
        self.mapping = mapping
        self._saved = None
        self._drive = None
        self._validated = False

    def propose(self, total, own):
        if self._saved is not None:
            raise InputError('finish the pending molecular trial first')
        right = _selected_molecular_field(total, own)
        drive = _interval_field(self._left,right) if self.sampling=='midpoint' else right
        self._accepted_right = None
        old = self.state.observables.dipole
        self._saved = self.state.snapshot()
        self._drive = drive
        try:
            new = self.state.step(drive).dipole
            return self.mapping.current([old], [new], dt=self.state.dt)
        except Exception:
            self.cancel()
            raise

    def validate(self, total, own):
        if self._saved is None:
            raise InputError('no pending molecular trial')
        self._validated = False
        right = _selected_molecular_field(total, own)
        actual = _interval_field(self._left,right) if self.sampling=='midpoint' else right
        if not np.allclose(actual, self._drive, rtol=1.e-10, atol=1.e-14):
            raise NumericalError('accepted molecular interval E differs from prediction; cancel and retry')
        self._accepted_right = right
        self._validated = True

    def commit(self):
        if self._saved is None or not self._validated:
            raise InputError('validate the accepted molecular endpoint before commit')
        self._left = self._accepted_right
        self._accepted_right = None
        self._saved = None
        self._drive = None
        self._validated = False

    def cancel(self):
        if self._saved is None:
            raise InputError('no pending molecular trial')
        # Keep the token if replay itself fails; the caller can retry restoration.
        self.state.restore(self._saved)
        self._accepted_right = None
        self._saved = None
        self._drive = None
        self._validated = False


def _energy_work(fields,dipoles,energies,fresh):
    """Fortran cumulative held/midpoint work and fresh internal-energy residuals."""
    import ctypes as c
    fields=np.asarray(fields)
    if fields.ndim!=2 or fields.shape[1]!=3 or not len(fields):
        raise InputError('nonempty endpoint E history required')
    n=len(fields)
    arrays=[np.ascontiguousarray(_real(a,shape,name)) for a,shape,name in (
        (fields,(n,3),'E history'),(dipoles,(n,3),'dipole history'),
        (energies,(n,),'energy history'),(fresh,(n,),'energy freshness'))]
    out=np.empty((n,6))
    fn=_molecular_function('sp_molecular_energy_work',[c.c_int]+[c.c_void_p]*5)
    status=fn(n,*[a.ctypes.data for a in arrays],out.ctypes.data)
    if status==1:raise InputError('invalid energy/work input or initial energy not fresh')
    if status==2:raise RuntimeFailure('energy/work buffer allocation failed')
    if status:raise NumericalError('energy/work arithmetic overflow')
    return _real(out,out.shape,'energy/work history')


def _interval_field(left,right):
    """Fortran interval-average selected E; the solver still receives constant E."""
    import ctypes as c
    left=np.ascontiguousarray(_real(left,(3,),'left selected E'))
    right=np.ascontiguousarray(_real(right,(3,),'right selected E'))
    out=np.empty(3)
    fn=_molecular_function('sp_molecular_interval_field',[c.c_void_p]*3)
    status=fn(left.ctypes.data,right.ctypes.data,out.ctypes.data)
    if status:raise NumericalError('invalid interval field')
    return out
