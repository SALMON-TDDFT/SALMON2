"""Charge-dipole sources and explicit molecular feedback, in atomic units, for caller-owned FDTD fields."""
from dataclasses import dataclass
import numpy as np
from .errors import InputError, NumericalError, RuntimeFailure
from .field_input import _real

__all__ = ['DipoleMap', 'drive_native_dipoles', 'MolecularDrive', 'select_molecular_field', 'MolecularFrame', 'MolecularFeedback']


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
        if probe()!=1:raise RuntimeFailure('unsupported molecular kernel protocol')
        fn=getattr(library,name)
    except AttributeError as exc:raise RuntimeFailure('matching molecular kernel required; explicitly build/install it') from exc
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
    Cancellation uses its direct live savepoint and never starts a replay worker.
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
        self._saved = self.state.savepoint()
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
        self.state.release_savepoint(self._saved)
        self._left = self._accepted_right
        self._accepted_right = None
        self._saved = None
        self._drive = None
        self._validated = False

    def cancel(self):
        if self._saved is None:
            raise InputError('no pending molecular trial')
        # Keep the live token if restoration/release fails; never replay.
        self.state.restore_savepoint(self._saved)
        self.state.release_savepoint(self._saved)
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

@dataclass(frozen=True, slots=True, eq=False)
class MolecularFrame:
    """Accepted molecular/field samples on the reporting rank, in atomic units."""
    field_step: int
    field_time: float
    molecule: object
    total_E: object
    own_E: object
    selected_E: object
    drive_E: object
    current: object
    field_energy: object = None

    def __post_init__(self):
        from .config import _copy_json, _freeze
        for name in ('total_E','own_E','selected_E','drive_E'):
            object.__setattr__(self,name,_real(getattr(self,name),(3,),name))
        object.__setattr__(self,'current',_real(self.current,(1,3),'charge current'))
        object.__setattr__(self,'field_energy',_freeze(_copy_json(self.field_energy)))

    def record(self):
        """Detached JSON-compatible accepted frame; no file I/O."""
        from .config import _thaw
        return dict(field_step=self.field_step,field_time=self.field_time,
            electronic_time=self.molecule.time,dipole=self.molecule.dipole.tolist(),
            electron_number=self.molecule.electron_number,internal_energy=self.molecule.internal_energy,
            total_E=self.total_E.tolist(),own_E=self.own_E.tolist(),selected_E=self.selected_E.tolist(),
            drive_E=self.drive_E.tolist(),current=self.current.tolist(),field_energy=_thaw(self.field_energy))


class MolecularFeedback:
    """Caller-owned single-molecule MPI feedback with own-vacuum subtraction.

    Call collectively on comm. Supply a local field on every field-group rank,
    None on bookkeeping ranks, and state only on total_root. Both roots must be
    field-group leaders. The reference must be vacuum driven only by this dipole,
    with no illumination; this physical declaration is the caller's responsibility.
    All objects/MPI lifetimes remain caller-owned. Numerical operations are Fortran.
    A recoverable rejected interval raises StepRecoveryError(original_preserved=True).
    Electron restoration uses direct live savepoints. No automatic retry or I/O.
    """
    def __init__(self,comm,field,*,state=None,mapping,total_root,reference_root,policy,observer=None):
        from .native_steps import _NativeBackup
        error=None;descriptor=None;backup=None;footprint=None
        try:
            if policy!='environment_scattered_only':raise InputError('explicit environment_scattered_only policy required')
            if (type(total_root) is not int or type(reference_root) is not int
                    or total_root==reference_root or not 0<=total_root<comm.size
                    or not 0<=reference_root<comm.size):raise InputError('distinct field leader ranks required')
            if not isinstance(mapping,DipoleMap) or len(mapping.molecule_ids)!=1 or len(mapping.point_ids)!=1:
                raise InputError('single molecule and single macro point required')
            if (state is not None)!=(comm.rank==total_root):raise InputError('molecular state belongs only on total_root')
            if state is not None and state.time!=0.:raise InputError('fresh molecular state required')
            if field is not None:
                if field.closed or field.failed or field.step_index!=0:raise InputError('fresh usable field required')
                grid=field.calculation.geometry.grid
                footprint=(tuple(grid.counts),tuple(grid.spacing),tuple(field.calculation.boundaries.faces),
                    tuple(map(tuple,field.manifest['grid_indices'])),
                    tuple(map(tuple,field.manifest['sample_positions_bohr'])),
                    field.manifest['sample_map'],field.manifest['deposit_map'])
                if tuple(field.point_ids)!=mapping.point_ids:raise InputError('field and dipole point IDs differ')
                if not np.allclose(mapping.cell_volumes,field.manifest['cell_volume'],rtol=1e-14,atol=0.):
                    raise InputError('mapping volumes differ from field cells')
            elif observer is not None:raise InputError('energy observer requires its local field')
            if observer is not None and observer.field is not field:raise InputError('observer belongs to another field')
            if field is not None:
                backup=_NativeBackup(field,field.context.solver_comm,'sp_eh_step_transaction',
                    field.directory,('A','E','step_index','failed','closed'))
            descriptor=dict(error=None,key=(total_root,reference_root,policy,mapping.molecule_ids,
                mapping.point_ids,mapping.weights.tobytes(),mapping.cell_volumes.tobytes()),
                field=field is not None,leader=field is not None and field.context.is_leader,
                dt=None if field is None else field.dt,capacity=None if field is None else field.capacity,
                state_dt=None if state is None else state.dt,observer=observer is not None,footprint=footprint)
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        records=comm.allgather(descriptor if error is None else dict(error=error))
        if any(r['error'] for r in records):raise InputError('molecular feedback setup: '+str(records))
        if any(r['key']!=records[0]['key'] for r in records):raise InputError('molecular feedback settings differ across ranks')
        leaders=[i for i,r in enumerate(records) if r['leader']]
        if sorted(leaders)!=sorted((total_root,reference_root)):raise InputError('roots must identify the two field leaders')
        if any(r['footprint']!=records[total_root]['footprint'] for r in records if r['field']):
            raise InputError('own-vacuum reference must use matching grid, boundaries and sampling/deposition footprint')
        dt=records[total_root]['dt']
        if len({r['observer'] for r in records if r['field']})!=1:
            raise InputError('enable energy observation on both complete field groups or neither')
        if any(r['dt']!=dt for r in records if r['field']) or records[total_root]['state_dt']!=dt:
            raise InputError('molecular and field time steps must match')
        self.comm=comm;self.field=field;self.state=state;self.mapping=mapping
        self.total_root=total_root;self.reference_root=reference_root;self.dt=dt
        self.capacity=min(r['capacity'] for r in records if r['field'])
        self.observer=observer;self.started=False;self.failed=False;self.step_index=0
        self._bridge=None;self._needs_rollback=False
        self._backup=backup

    def _field_step_indices(self):
        from .api import _bind_field_operations
        _bind_field_operations(self.field.lib,('sp_eh_step_indices',))
        values=np.empty(2,dtype=np.int32)
        if self.field._call(self.field.lib.sp_eh_step_indices,values):
            raise RuntimeFailure('cannot query actual Fortran FDTD step indices')
        return tuple(int(v) for v in values)

    def _bundle_rows(self, *, saving):
        import json,hashlib
        from .fdtd import FDTDMPI
        row=None;error=None
        try:
            if self.failed or (saving and not self.started):
                raise InputError('save requires an accepted usable feedback state')
            if self._backup is not None and self._backup.active:
                raise InputError('finish the pending field update before checkpointing')
            if self._bridge is not None and self._bridge._saved is not None:
                raise InputError('finish the pending molecular update before checkpointing')
            field=self.field
            contract=dict(rank=self.comm.rank,roots=[self.total_root,self.reference_root],
                dt=self.dt,capacity=self.capacity,molecules=list(self.mapping.molecule_ids),
                points=list(self.mapping.point_ids),weights=self.mapping.weights.tolist(),
                volumes=self.mapping.cell_volumes.tolist(),field=None,molecule=None)
            row=dict(contract=contract,session_step=self.step_index,started=self.started,field_step=None,molecular_step=None)
            if field is not None:
                if field.closed or field.failed or field.step_index!=self.step_index:
                    raise InputError('field and feedback accepted step indices differ')
                from .api import _bind_field_operations
                _bind_field_operations(field.lib,('sp_eh_step_indices','sp_eh_checkpoint'))
                if getattr(field,'_restore_candidate',None) is not None:
                    raise InputError('field restore candidate already pending')
                assets=field.calculation.input_assets()
                specification=dict(input=FDTDMPI(field.checkpoint_decomposition).input_text(field.calculation),
                    assets=assets,library=field.manifest['library_sha256'],
                    indices=field.manifest['grid_indices'],positions=field.manifest['sample_positions_bohr'],
                    sample=field.manifest['sample_map'],deposit=field.manifest['deposit_map'],
                    decomposition=list(field.checkpoint_decomposition),
                    observer=None if self.observer is None else self.observer.box.tolist())
                contract['field']=hashlib.sha256(json.dumps(specification,sort_keys=True).encode()).hexdigest()
                row.update(field_step=field.step_index,A=field.A.tolist(),E=field.E.tolist(),leader=field.context.is_leader)
            if self.state is not None:
                if not self.state.usable or self.state._engine._savepoint is not None:
                    raise InputError('usable molecule without a pending savepoint required')
                statistics=self.state.worker_statistics()
                contract['molecule']=self.state._engine._snapshot_settings()
                row['molecular_step']=statistics['step']
                if self.started and statistics['step']!=self.step_index-1:
                    raise InputError('molecule must be one step behind the accepted FDTD step')
                row['left']=None if self._bridge is None else self._bridge._left.tolist()
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        records=self.comm.allgather((error,row))
        if any(e for e,_ in records):raise InputError('molecular checkpoint admission: '+str([e for e,_ in records]))
        # All ranks passed local admission before any owning-field-group query.
        error=None
        try:
            if self.field is not None:
                actual,pending=self._field_step_indices()
                if actual!=self.field.step_index or pending!=-1:
                    raise InputError('Fortran field step differs or a restore is pending')
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        errors=self.comm.allgather(error)
        if any(errors):raise InputError('molecular checkpoint field indices: '+str(errors))
        rows=[v for _,v in records]
        if len({v['started'] for v in rows})!=1:raise InputError('feedback acceptance phases differ across ranks')
        if len({v['session_step'] for v in rows})!=1:raise InputError('feedback step numbers differ across ranks')
        if self.started:
            _check_molecular_steps(rows[self.total_root]['field_step'],rows[self.reference_root]['field_step'],
                                   rows[self.total_root]['molecular_step'])
        return rows

    def save_checkpoint(self, directory, *, history=None):
        """Collectively save both FDTD groups, their LD state and the molecule.

        Call after accepting a step. Optional JSON frame history belongs only on
        total_root. Writes a new directory; no update, replay or hidden SCF.
        """
        import json,os,tempfile,shutil
        from pathlib import Path
        from .tddft import _digest
        from .api import _bind_field_operations
        rows=self._bundle_rows(saving=True)
        destination=Path(directory).expanduser().absolute()
        paths=self.comm.allgather(str(destination))
        if any(v!=paths[0] for v in paths):raise InputError('shared checkpoint destination required')
        prepared=None
        if self.comm.rank==self.total_root:
            try:
                if destination.exists() or destination.is_symlink():raise InputError('checkpoint destination exists')
                if history is not None and (not isinstance(history,list) or len(history)!=self.step_index):
                    raise InputError('optional history needs one frame per accepted FDTD step')
                json.dumps(history,allow_nan=False)
                destination.parent.mkdir(parents=True,exist_ok=True)
                prepared=(None,tempfile.mkdtemp(prefix='.molecular-bundle-',dir=destination.parent))
            except Exception as exc:prepared=(str(exc),None)
        error,temporary=self.comm.bcast(prepared,root=self.total_root)
        if error:raise InputError(error)
        root=Path(temporary);files={};error=None
        try:
            if self.field is not None:
                _bind_field_operations(self.field.lib,('sp_eh_checkpoint',))
                path=root/('field-'+str(self.comm.rank)+'.bin')
                if self.field._call(self.field.lib.sp_eh_checkpoint,os.fsencode(path),0):
                    raise RuntimeFailure('FDTD direct save failed')
                files[path.name]=_digest(path)
                if self.field.context.is_leader:
                    stream=root/('streams-'+str(self.comm.rank));stream.mkdir()
                    for path in self.field.directory.glob('*.data'):
                        if path.is_symlink() or not path.is_file():raise InputError('regular private field stream required')
                        target=stream/path.name;shutil.copyfile(path,target)
                        files[target.relative_to(root).as_posix()]=_digest(target)
            if self.state is not None:
                self.state.save_direct_checkpoint(root/'molecule')
                for path in (root/'molecule').iterdir():
                    files[path.relative_to(root).as_posix()]=_digest(path)
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        reports=self.comm.allgather((error,files))
        error=None
        if self.comm.rank==self.total_root:
            try:
                if any(e for e,_ in reports):raise RuntimeFailure(str([e for e,_ in reports]))
                inventory={k:v for _,part in reports for k,v in part.items()}
                record=dict(schema='salmonpy.molecular-feedback-direct/1',field_step=self.step_index,
                    molecular_step=self.step_index-1,rows=rows,files=inventory,history=history)
                (root/'complete.json').write_text(json.dumps(record,allow_nan=False,indent=2)+'\n')
                root.rename(destination)
            except Exception as exc:
                error=type(exc).__name__+': '+str(exc)
                shutil.rmtree(root,ignore_errors=True)
        error=self.comm.bcast(error,root=self.total_root)
        if error:raise RuntimeFailure('molecular checkpoint save: '+error)

    def restore_checkpoint(self, directory):
        """Collectively stage all three calculations before publishing any state.

        Returns saved history/indices on total_root, None elsewhere. Requires the
        same inputs/runtime/MPI layout. Never calls start() or elapsed step().
        """
        import json
        from pathlib import Path, PurePosixPath
        from .tddft import _digest
        from ._native_checkpoint import NativeRestoreCandidate
        from .isolated_tddft import _EEngine
        rows=self._bundle_rows(saving=False)
        root=Path(directory).expanduser().absolute()
        paths=self.comm.allgather(str(root))
        if any(v!=paths[0] for v in paths):raise InputError('shared checkpoint source required')
        payload=None
        if self.comm.rank==self.total_root:
            try:
                record=json.loads((root/'complete.json').read_text())
                if record['schema']!='salmonpy.molecular-feedback-direct/1':raise InputError('unsupported molecular checkpoint')
                _check_molecular_steps(record['field_step'],record['field_step'],record['molecular_step'])
                if len(record['rows'])!=self.comm.size:raise InputError('saved MPI layout differs')
                for current,saved in zip(rows,record['rows']):
                    if current['contract']!=saved['contract']:raise InputError('saved model/runtime/MPI configuration differs')
                    if saved['session_step']!=record['field_step']:raise InputError('saved session step differs')
                    if saved['field_step'] is not None and saved['field_step']!=record['field_step']:
                        raise InputError('saved FDTD step differs')
                _check_molecular_steps(record['rows'][self.total_root]['field_step'],
                    record['rows'][self.reference_root]['field_step'],record['rows'][self.total_root]['molecular_step'])
                history=record['history']
                if history is not None and (not isinstance(history,list) or len(history)!=record['field_step']):
                    raise InputError('saved history count differs from the accepted field step')
                if not isinstance(record['files'],dict) or not record['files']:raise InputError('missing payload inventory')
                required={'molecule/manifest.json'}
                required.update('field-'+str(i)+'.bin' for i,row in enumerate(rows) if row['field_step'] is not None)
                if not required.issubset(record['files']):raise InputError('required numerical payload missing')
                for name in record['files']:
                    path=PurePosixPath(name)
                    if path.is_absolute() or '..' in path.parts or str(path)!=name or '\\' in name:
                        raise InputError('unsafe checkpoint path')
                    if any((root/Path(*path.parts[:i])).is_symlink() for i in range(1,len(path.parts)+1)):
                        raise InputError('checkpoint symlinks unsupported')
                    if not (root/name).is_file():raise InputError('regular checkpoint payload required')
                if set(record['files'])!={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p!=root/'complete.json'}:
                    raise InputError('checkpoint payload inventory differs')
                payload=(None,record)
            except Exception as exc:payload=(type(exc).__name__+': '+str(exc),None)
        error,record=self.comm.bcast(payload,root=self.total_root)
        if error:raise InputError('molecular checkpoint metadata: '+error)
        error=None
        try:
            for index,(name,sha) in enumerate(sorted(record['files'].items())):
                if index%self.comm.size==self.comm.rank and _digest(root/name)!=sha:
                    raise InputError('checkpoint checksum mismatch: '+name)
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        errors=self.comm.allgather(error)
        if any(errors):raise InputError('molecular checkpoint integrity: '+str(errors))
        field_candidate=molecule_candidate=streams=bridge=None
        error=None
        saved=record['rows'][self.comm.rank]
        try:
            if self.field is not None:
                field_candidate=NativeRestoreCandidate(self.field,'sp_eh_checkpoint',root/('field-'+str(self.comm.rank)+'.bin'))
                field_candidate.stage()
                _,pending=self._field_step_indices()
                if pending!=record['field_step']:raise InputError('staged Fortran FDTD step differs from manifest')
                A=_real(saved['A'],self.field.A.shape,'saved A')
                E=_real(saved['E'],self.field.E.shape,'saved E')
                if self.field.context.is_leader:
                    streams=_MolecularStreams(self.field,root/('streams-'+str(self.comm.rank)))
            if self.state is not None:
                engine=self.state._engine
                molecule_candidate=_EEngine(engine._model,engine.dt,engine._position)
                molecule_candidate.restore_direct_checkpoint(root/'molecule')
                stats=molecule_candidate.worker_statistics()
                if stats['step']!=record['molecular_step'] or stats['updates']!=0:
                    raise InputError('molecular restore step differs or performed propagation')
                left=_real(saved['left'],(3,),'saved selected E')
                selected=_selected_molecular_field(record['rows'][self.total_root]['E'][0],
                                                   record['rows'][self.reference_root]['E'][0])
                if not np.array_equal(left,selected):raise InputError('saved selected E differs from the two field samples')
                bridge=_IsolatedDipoleTrial(self.state,self.mapping,initial_field=left)
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        errors=self.comm.allgather(error)
        if not any(errors):
            try:
                if streams is not None:streams.publish()
            except Exception as exc:error=type(exc).__name__+': '+str(exc)
            errors=self.comm.allgather(error)
        if any(errors):
            cleanup=[]
            for candidate,method in ((field_candidate,'discard'),(molecule_candidate,'close'),(streams,'rollback')):
                if candidate is not None:
                    try:getattr(candidate,method)()
                    except Exception as exc:cleanup.append(str(exc))
            if streams is not None:
                try:streams.close()
                except Exception as exc:cleanup.append(str(exc))
            failures=self.comm.allgather(cleanup)
            if any(failures):self.failed=True
            from .checkpointing import AtomicRestoreError
            if any(failures):
                raise RuntimeFailure('molecular restore cleanup failed; preservation not guaranteed: '+str(failures))
            raise AtomicRestoreError('molecular candidate restore failed: '+str(errors))
        # All I/O and native payload reads finished before this publication phase.
        error=None
        try:
            if field_candidate is not None:field_candidate.publish()
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        errors=self.comm.allgather(error)
        if any(errors):
            self.failed=True
            if molecule_candidate is not None:molecule_candidate.close()
            if streams is not None:streams.close()
            raise RuntimeFailure('native field publication failed after preparation; preservation not guaranteed: '+str(errors))
        if self.field is not None:
            self.field.A=np.array(A,copy=True);self.field.E=np.array(E,copy=True)
            self.field.step_index=record['field_step']
        retired=None
        if self.state is not None:
            retired=self.state._engine;self.state._engine=molecule_candidate
            self._bridge=bridge
        self.step_index=record['field_step'];self.started=True
        cleanup=[]
        for resource in (streams,retired):
            if resource is not None:
                try:resource.close()
                except Exception as exc:cleanup.append(str(exc))
        error=None;verified=dict(rank=self.comm.rank,field_step=None,molecular_step=None)
        try:
            if self.field is not None:
                actual,pending=self._field_step_indices()
                if actual!=self.step_index or pending!=-1:raise RuntimeFailure('published Fortran field step differs')
                verified['field_step']=actual
            if self.state is not None:
                statistics=self.state.worker_statistics()
                if statistics['step']!=record['molecular_step'] or statistics['updates']!=0:
                    raise RuntimeFailure('published molecule step differs or advanced during restore')
                verified.update(molecular_step=statistics['step'],molecular_updates=statistics['updates'])
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        checked=self.comm.allgather((error,verified,cleanup))
        if any(e for e,_,_ in checked):
            self.failed=True
            raise RuntimeFailure('published step verification failed: '+str(checked))
        cleanup_errors=[e for _,_,errors in checked for e in errors]
        if cleanup_errors and self.comm.rank==self.total_root:
            import warnings
            warnings.warn('checkpoint restored; retired resource cleanup: '+str(cleanup_errors),RuntimeWarning)
        return dict(field_step=self.step_index,molecular_step=record['molecular_step'],
                    verified_steps=[v for _,v,_ in checked],history=record['history']) if self.comm.rank==self.total_root else None

    def _capture(self):
        self._needs_rollback=False
        token=None;error=None
        if self._backup is None:return token,error
        try:
            if self.field.step_index!=self.step_index:raise InputError('caller changed the field clock during feedback')
            token=self._backup.capture()
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        # A leader-only append-journal failure must not let peers advance.
        errors=self.field.context.solver_comm.allgather(error)
        return token,next((e for e in errors if e),None)

    def _packet(self,error=None,*,observe=False):
        try:
            energy=None if self.observer is None or not observe else self.observer.sample_step()
            E=None if self.field is None else self.field.E[0].copy()
        except Exception as exc:
            error=type(exc).__name__+': '+str(exc);E=None;energy=None
        return dict(error=error,E=E,energy=energy)

    @staticmethod
    def _errors(packets):
        return tuple(p['error'] for p in packets if p['error'])

    def _recover(self,token,errors):
        from .stepping import StepRecoveryError
        recovery=[]
        if self._backup is not None and self._backup.active and self._needs_rollback:
            try:self._backup.rollback(token)
            except Exception as exc:recovery.append(str(exc))
        if self._bridge is not None and self._bridge._saved is not None:
            try:self._bridge.cancel()
            except Exception as exc:recovery.append(str(exc))
        if self._backup is not None:
            try:self._backup.release(token)
            except Exception as exc:recovery.append(str(exc))
        if self.state is not None and not self.state.usable:
            recovery.append('isolated worker is not usable; no replay recovery attempted')
        # Recovery-only vote: successful steps retain the existing four exchanges.
        failures=tuple(e for group in self.comm.allgather(recovery) for e in group)
        self.failed=bool(failures)
        raise StepRecoveryError('molecular feedback',errors,original_preserved=not failures,rollback_errors=failures)

    def _frame(self,packets,drive,current):
        total=packets[self.total_root];own=packets[self.reference_root]
        return MolecularFrame(self.step_index+1,(self.step_index+1)*self.dt,self.state.observables,
            total['E'],own['E'],_selected_molecular_field(total['E'],own['E']),drive,current,
            [total['energy'],own['energy']])

    def start(self):
        """Collective zero-current hold; report the immutable first frame on total_root."""
        if self.started or self.failed:raise InputError('start requires a fresh usable feedback session')
        token,error=self._capture()
        try:
            if error is None and self.field is not None:
                self._needs_rollback=True
                self.field.advance(np.zeros((1,3)))
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        packets=self.comm.gather(self._packet(error,observe=error is None),root=self.total_root)
        error=None;frame=None
        if self.comm.rank==self.total_root:
            try:
                failures=self._errors(packets)
                if failures:raise RuntimeFailure(str(failures))
                selected=_selected_molecular_field(packets[self.total_root]['E'],packets[self.reference_root]['E'])
                if np.any(selected):raise InputError('initial molecular E must be zero; increase source separation or delay')
                self._bridge=_IsolatedDipoleTrial(self.state,self.mapping,initial_field=selected)
                frame=self._frame(packets,np.zeros(3),np.zeros((1,3)))
            except Exception as exc:error=type(exc).__name__+': '+str(exc)
        error=self.comm.bcast(error,root=self.total_root)
        if error:self._recover(token,(error,))
        if self._backup is not None:self._backup.release(token)
        self.started=True;self.step_index=1
        return frame

    def step(self,*,validate=None):
        """Collective accepted interval; optional root-only read-only frame validator.

        A validator exception rejects the interval. Catch StepRecoveryError on
        every rank and retry only if original_preserved is true. The callback
        must not mutate states or perform collectives. No automatic retries.
        """
        if not self.started or self.failed or self.step_index>=self.capacity:
            raise InputError('started usable feedback with remaining field capacity required')
        token,error=self._capture()
        try:
            if error is None and self.field is not None:
                self._needs_rollback=True
                self.field.advance(np.zeros((1,3)))
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        packets=self.comm.gather(self._packet(error),root=self.total_root)
        payload=None
        if self.comm.rank==self.total_root:
            try:
                failures=self._errors(packets)
                if failures:raise RuntimeFailure(str(failures))
                if self.state.time!=(self.step_index-1)*self.dt:
                    raise InputError('caller changed the molecular clock during feedback')
                current=self._bridge.propose(packets[self.total_root]['E'],packets[self.reference_root]['E'])
                payload=(None,current)
            except Exception as exc:payload=(type(exc).__name__+': '+str(exc),None)
        error,current=self.comm.bcast(payload,root=self.total_root)
        if error:self._recover(token,(error,))
        error=None
        try:
            if self._backup is not None:
                self._backup.rollback(token);self._backup.release(token);token=None
                token,error=self._capture()
            if error is None and self.field is not None:
                self._needs_rollback=True
                self.field.advance(current)
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        packets=self.comm.gather(self._packet(error,observe=error is None),root=self.total_root)
        error=None;frame=None
        if self.comm.rank==self.total_root:
            try:
                failures=self._errors(packets)
                if failures:raise RuntimeFailure(str(failures))
                self._bridge.validate(packets[self.total_root]['E'],packets[self.reference_root]['E'])
                frame=self._frame(packets,self._bridge._drive,current)
                if validate is not None:validate(frame)
                self._bridge.commit()
            except Exception as exc:error=type(exc).__name__+': '+str(exc)
        error=self.comm.bcast(error,root=self.total_root)
        if error:self._recover(token,(error,))
        if self._backup is not None:self._backup.release(token)
        self.step_index+=1
        return frame


def _check_molecular_steps(total_step, reference_step, molecular_step):
    """Accepted FDTD n and molecular n-1; never infer or shift a saved index."""
    if (any(type(v) is not int for v in (total_step,reference_step,molecular_step))
            or total_step < 1 or reference_step != total_step or molecular_step != total_step-1):
        raise InputError('accepted steps must be total FDTD n, reference FDTD n, molecule n-1')
    return total_step,molecular_step


class _MolecularStreams:
    """Stage closed-per-update FDTD point streams; retain originals until commit."""
    def __init__(self, field, source):
        import tempfile
        import shutil
        from pathlib import Path
        self.root=field.directory
        self.work=Path(tempfile.mkdtemp(prefix='.molecular-streams-',dir=self.root))
        self.changed=False
        try:
            self.originals={p.name:p.read_bytes() for p in self.root.glob('*.data')}
            self.incoming={p.name:p.read_bytes() for p in source.glob('*.data')} if source.exists() else {}
            for name,raw in self.incoming.items():(self.work/name).write_bytes(raw)
        except BaseException:
            shutil.rmtree(self.work,ignore_errors=True)
            raise

    def publish(self):
        import os
        self.changed=True
        for name in self.incoming:os.replace(self.work/name,self.root/name)
        for name in self.originals.keys()-self.incoming.keys():(self.root/name).unlink()

    def rollback(self):
        import os
        if not self.changed:return
        for name,raw in self.originals.items():
            path=self.work/name;path.write_bytes(raw);os.replace(path,self.root/name)
        for name in self.incoming.keys()-self.originals.keys():(self.root/name).unlink(missing_ok=True)
        self.changed=False

    def close(self):
        import shutil
        shutil.rmtree(self.work)
