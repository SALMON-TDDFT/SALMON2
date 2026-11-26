"""Explicit isolated PZ length-gauge E -> charge-dipole process interface."""
from dataclasses import dataclass, field
import json
from pathlib import Path
import numpy as np
from .errors import InputError, NumericalError, RuntimeFailure, SnapshotError
from .fdtd import SalmonExecutable
from .field_input import FieldContract, FieldPacket, FieldValue, _real
from .timing import TimeSpec, InputPort
from .snapshots import _fingerprint
from .tddft import TDDFTInitialState, TDDFTState, _digest
from .ld import _array_unchecked
from .sbe import _frozen,_integer

__all__ = ['IsolatedRuntime','load_isolated_runtime','DipoleObservables','ChargeDensitySamples','OrbitalMagneticResponse',
           'IsolatedTDDFTInitialState','IsolatedTDDFTState','IsolatedSavepoint','run_molecular_mpi','MolecularMPIController']
_CAPABILITIES = {'isolated','pz','fixed-ions','middlepoint','single-rank','E-to-dipole','measured-electron-number','internal-energy-freshness','live-savepoint','direct-disk-state'}


def _receipt(payload, executable):
    if type(payload) is not bytes or len(payload) > 65536:
        raise InputError('bounded isolated-worker validation bytes required')
    def unique(pairs):
        result = {}
        for k,v in pairs:
            if k in result: raise InputError('duplicate receipt key')
            result[k] = v
        return result
    try:
        data = json.loads(payload, object_pairs_hook=unique)
        valid = (data['schema']=='salmonpy.isolated-worker-validation/5'
                 and data['protocol']=='SALMONPY_TDDFT_E_V3'
                 and data['executable_sha256']==executable.sha256
                 and type(data['steps']) is int and data['steps']>=20
                 and data['interleaved_independence'] is True
                 and set(data['capabilities'])==_CAPABILITIES
                 and data['direct_disk_process_exit_exact'] is True
                 and type(data['direct_disk_extra_evolutions']) is int and data['direct_disk_extra_evolutions']==0
                 and data['savepoint_same_process'] is True
                 and type(data['savepoint_extra_evolutions']) is int and data['savepoint_extra_evolutions']==0
                 and data['savepoint_stale_rejected'] is True)
        for key in ('repeat_max_abs_error','replay_continuation_max_abs_error','reference_max_abs_error','reference_energy_max_abs_error','savepoint_continuation_max_abs_error'):
            value=data[key]
            valid = valid and type(value) in (int,float) and np.isfinite(value) and 0 <= value <= 1.e-10
    except (KeyError, TypeError, ValueError, RecursionError) as exc:
        raise InputError('invalid isolated-worker revision5 validation receipt') from exc
    if 'native_near_field' in data:
        evidence = data['native_near_field']
        try:
            near_valid = (evidence['schema']=='salmonpy.isolated-near-field-validation/1'
                          and type(evidence['steps']) is int and evidence['steps']>=20
                          and evidence['interleaved_independence'] is True
                          and evidence['savepoint_same_process'] is True
                          and evidence['direct_disk_process_exit_exact'] is True)
            for key in ('reference_max_abs_error','zero_source_max_abs_error'):
                value=evidence[key]
                near_valid = near_valid and type(value) in (int,float) and np.isfinite(value) and 0<=value<=1.e-10
            near_valid = near_valid and evidence['zero_source_max_abs_error']==0
        except (KeyError,TypeError,ValueError):
            near_valid=False
        if not near_valid:
            raise InputError('invalid near-field validation evidence')
    if 'solid_field' in data:
        try:
            solid=data['solid_field']
            solid_valid=(solid['schema']=='salmonpy.solid-field-validation/1'
                         and solid['direct_disk_process_exit_exact'] is True
                         and solid['savepoint_same_process'] is True
                         and solid['uniform_max_abs_error']>=0
                         and solid['uniform_max_abs_error']<=1.e-10
                         and solid['quadrupole_response_nonzero'] is True)
        except (KeyError,TypeError,ValueError):
            solid_valid=False
        if not solid_valid:raise InputError('invalid spherical-harmonic field validation evidence')
    if 'charge_density' in data:
        evidence=data['charge_density']
        if (not isinstance(evidence,dict) or evidence.get('schema')!='salmonpy.charge-density-validation/1'
            or not all(evidence.get(k) is True for k in ('signed_charge','read_only','direct_restart_exact'))):
            raise InputError('invalid charge-density validation evidence')
    if 'orbital_magnetic' in data:
        evidence=data['orbital_magnetic']
        if (not isinstance(evidence,dict) or evidence.get('schema')!='salmonpy.orbital-magnetic-validation/1'
            or not all(evidence.get(k) is True for k in ('hamiltonian_derivative','nonlocal_projectors','direct_restart_exact'))):
            raise InputError('invalid orbital magnetic validation evidence')
    if 'grid_stream' in data:
        evidence=data['grid_stream']
        if (not isinstance(evidence,dict) or evidence.get('schema')!='salmonpy.grid-stream-validation/1'
            or not all(evidence.get(k) is True for k in ('prescribed_matches','direct_restart_exact','rejection_preserves_state'))):
            raise InputError('invalid live grid validation evidence')
    if 'near_eh' in data:
        try:
            evidence=data['near_eh']
            near_valid=(evidence['schema']=='salmonpy.near-eh-validation/3'
                        and evidence['radial_terms'] is True
                        and evidence['exponential_envelope'] is True
                        and evidence['uniform_matches'] is True
                        and evidence['magnetic_effect'] is True
                        and evidence['electric_quadrupole'] is True
                        and evidence['checkpoint_process_exit_exact'] is True
                        and evidence['center_preserved'] is True
                        and type(evidence['restore_extra_evolutions']) is int
                        and evidence['restore_extra_evolutions']==0)
        except (KeyError,TypeError,ValueError):
            near_valid=False
        if not near_valid:raise InputError('invalid spatial E/H validation evidence')
    if 'grid_eh' in data:
        try:
            evidence=data['grid_eh']
            grid_valid=(evidence['schema']=='salmonpy.grid-field-validation/1'
                and all(evidence[k] is True for k in ('uniform_matches','off_grid_linear_exact',
                    'staggered_components','cell_crossing_integral','checkpoint_process_exit_exact'))
                and type(evidence['restore_extra_evolutions']) is int and evidence['restore_extra_evolutions']==0)
        except (KeyError,TypeError):grid_valid=False
        if not grid_valid:raise InputError('invalid direct-grid field validation evidence')
    if not valid or _digest(executable.path)!=executable.sha256:
        raise InputError('isolated E-worker revision5 evidence or binary hash does not match')


@dataclass(frozen=True, slots=True)
class IsolatedRuntime:
    """Locally validated E3 binary; a receipt is evidence, not a signature or audit."""
    executable: SalmonExecutable
    validation: bytes = field(repr=False)

    def __post_init__(self):
        if not isinstance(self.executable,SalmonExecutable):
            raise InputError('SalmonExecutable required')
        _receipt(self.validation,self.executable)


def load_isolated_runtime(executable, validation):
    """Verify current isolated-worker evidence and executable hash; never build/download."""
    binary = Path(executable).expanduser().resolve(strict=True)
    path = Path(validation).expanduser()
    with path.open('rb') as stream:
        payload = stream.read(65537)
    return IsolatedRuntime(SalmonExecutable(binary,_digest(binary)),payload)


@dataclass(frozen=True, slots=True, eq=False)
class DipoleObservables:
    time: float
    electronic_change: object
    dipole: object
    electron_number: float
    internal_energy: object

    def __post_init__(self):
        object.__setattr__(self,'time',float(_real(self.time,(),'time')))
        object.__setattr__(self,'electronic_change',_real(self.electronic_change,(3,),'electronic dipole change'))
        object.__setattr__(self,'dipole',_real(self.dipole,(3,),'total charge dipole'))
        object.__setattr__(self,'electron_number',float(_real(self.electron_number,(),'electron number')))
        if self.internal_energy is not None:
            object.__setattr__(self,'internal_energy',float(_real(self.internal_energy,(),'internal energy')))


    @classmethod
    def _from_worker(cls,time,values,*,checked):
        if checked:return cls(time,values[:3],values[3:6],values[6],values[7] if values[8]==1. else None)
        # Transport already fixed the nine-element ABI; immutable owned views.
        values=_frozen(values,check=False)
        result=object.__new__(cls)
        for key,value in (('time',float(time)),('electronic_change',values[:3]),
            ('dipole',values[3:6]),('electron_number',float(values[6])),
            ('internal_energy',float(values[7]) if values[8]==1. else None)):
            object.__setattr__(result,key,value)
        return result


@dataclass(frozen=True, slots=True, eq=False)
class IsolatedSavepoint:
    """Opaque live-worker token; explicitly release it, never serialize it."""
    step_index: int
    time: float
    worker_pid: int
    updates: int
    _owner: object = field(repr=False)
    _observables: object = field(repr=False)
    _journal: object = field(repr=False)

    def __post_init__(self):
        object.__setattr__(self, '_observables', _real(self._observables, (9,), 'savepoint observations'))


def _near_feature(path):
    with Path(path).open() as stream:magic=stream.readline().strip()
    if magic=='SALMONPY_NEAR_GRID_V1':return 'grid_eh'
    if magic=='SALMONPY_NEAR_EH_V3':return 'near_eh'
    raise InputError('unsupported prepared spatial E/H format; regenerate the field input')


class _EEngine(TDDFTState):
    _worker_mode = 'E3'
    _reply_prefix = ('SALMONPY_TDDFT_E_V3 ', 'SALMONPY_E_STATE_V1 ')
    _init_command = 'INIT_E'
    _performance_init_command = 'INIT_E_PERFORMANCE'
    _performance_reply_prefix = ('SALMONPY_TDDFT_E_PERFORMANCE_V3 ', 'SALMONPY_E_STATE_V1 ')
    _step_command = 'STEP_E'
    _reply_size = 9
    _snapshot_kind = 'isolated-TDDFT-E-dipole-replay'

    def __init__(self,model,dt,position,*,execution_mode="checked",native_near_field=False,electric_potential=False,near_field=False):
        self._worker_mode = "E3M" if near_field else ("E3Y" if electric_potential else ("E3N" if native_near_field else "E3"))
        self._near_spec = None
        if near_field:
            with (model._snapshot/'salmonpy-near-eh.dat').open() as stream:
                magic=stream.readline().strip()
                if magic=='SALMONPY_NEAR_GRID_V1':
                    self._near_spec=('grid',)
                else:
                    lmax=int(stream.readline().split()[0])
                    geometry=tuple(float(v) for v in stream.readline().split()[:7])
                    self._near_spec=(lmax,geometry)
        self._position = position
        super().__init__(model,dt,np.zeros(3),execution_mode=execution_mode)

    @property
    def current(self):
        # Private transport hook: nine raw observables, never exposed as current.
        if self.execution_mode=='performance':return _frozen(self._current,check=False)
        return _real(self._current,(9,),'dipole/electron/energy observables')

    def _receive(self,index):
        try:
            values=super()._receive(index)
            if values[8] not in (0.,1.):raise InputError('invalid energy freshness flag')
            return values
        except InputError as exc:
            raise NumericalError('isolated worker returned nonfinite observables') from exc

    def charge_density(self):
        if not self.usable:raise SnapshotError('density read requires an accepted usable worker')
        status,step,_=self._command('DENSITY_E',8)
        if status:raise RuntimeFailure('density transfer rejected without propagation: '+str(status))
        if step!=self._steps:raise RuntimeFailure('density transfer step differs from state')
        path=Path(self._work.name)/'run/salmonpy-charge-samples.bin'
        return _decode_charge_density(path.read_bytes(),step=step,dt=self.dt)

    def orbital_magnetic_response(self):
        if not self.usable:raise SnapshotError('magnetic read requires an accepted usable worker')
        status,step,_=self._command('ORBITAL_E',10)
        if status:raise RuntimeFailure('orbital magnetic query rejected without propagation: '+str(status))
        if step!=self._steps:raise RuntimeFailure('orbital response step differs from state')
        raw=(Path(self._work.name)/'run/salmonpy-orbital-response.bin').read_bytes()
        if len(raw)<192:raise RuntimeFailure('truncated orbital response')
        magic,version,actual,nx,ny,nz=map(int,np.frombuffer(raw[:24],dtype=np.int32))
        if magic!=16909060 or version!=1 or actual!=step or min(nx,ny,nz)<2 or len(raw)!=192+nx*ny*nz*24:
            raise RuntimeFailure('orbital response geometry/step/size differs')
        metadata=np.frombuffer(raw,offset=24,count=21,dtype=np.float64)
        M=np.frombuffer(raw,offset=192,dtype=np.float64).reshape(nz,ny,nx,3)
        return OrbitalMagneticResponse(step,step*self.dt,metadata[:18].reshape(2,3,3),metadata[18:],M)

    def _stage_grid_field(self,window,n):
        if not self.usable:raise SnapshotError('live grid staging requires an accepted usable worker')
        _write_grid_packet(Path(self._work.name)/'run/salmonpy-grid-frames.bin',window,n,self.dt)
        status,step,_=self._command('GRID_E',9)
        if status:raise InputError('live grid transfer rejected; check fixed geometry, target step and complete finite E/H: '+str(status))
        if step!=self._steps:raise RuntimeFailure('live grid staging advanced or changed the worker step')

    def _snapshot_settings(self):
        return _fingerprint((super()._snapshot_settings(),tuple(self._position),self._worker_mode))

    def _spawn(self):
        self._savepoint = None
        self._restore_pending = None
        self._restore_updates = None
        self._direct_failed = False
        self._direct_restores = 0
        self._direct_origin = False
        super()._spawn()
        self._command('PROBE_E', 5)

    def _command(self, command, operation):
        """Receive a state-control reply, with no spawn, replay or propagation."""
        import queue
        if self._closed or self._direct_failed or self._proc is None or self._proc.poll() is not None:
            raise RuntimeFailure('isolated worker is not live; direct savepoint cannot replay it')
        try:
            self._send(command)
            line = self._events.get(timeout=self._timeout)
            if line is None:
                raise RuntimeFailure('isolated worker stopped during savepoint operation')
            values = line.split()
            if len(values) != 5 or values[0] != 'SALMONPY_E_STATE_V1':
                raise RuntimeFailure('isolated savepoint protocol mismatch')
            op, status, step, updates = map(int, values[1:])
            if op != operation or not 0 <= step <= self._capacity or updates < 0:
                raise RuntimeFailure('isolated savepoint sequence mismatch')
            return status, step, updates
        except Exception as exc:
            self._direct_failed = True
            raise RuntimeFailure('isolated savepoint communication failed; no replay attempted') from exc

    @property
    def usable(self):
        return (not self._closed and not self._direct_failed and self._restore_pending is None
                and self._proc is not None and self._proc.poll() is None)

    def worker_statistics(self):
        status, step, updates = self._command('PROBE_E', 5)
        if status:
            raise RuntimeFailure('isolated worker status rejected')
        return dict(pid=self._proc.pid, step=step, updates=updates,
                    direct_restores=self._direct_restores, suspended=self._restore_pending is not None)

    def _token(self, token):
        if (type(token) is not IsolatedSavepoint or token is not self._savepoint
                or token._owner is not self or self._proc is None or token.worker_pid != self._proc.pid):
            raise SnapshotError('savepoint is released, stale or belongs to another worker')

    def savepoint(self):
        from .native_steps import _AppendJournal
        if self._savepoint is not None:
            raise SnapshotError('one live savepoint per isolated worker; release it explicitly')
        status, step, updates = self._command('CAPTURE_E', 1)
        if status:
            raise SnapshotError('worker declined savepoint creation: '+str(status))
        if step != self._steps:
            self._direct_failed = True
            raise RuntimeFailure('savepoint clock differs from the accepted state')
        try:
            journal = _AppendJournal.capture(Path(self._work.name)/'run')
            token = IsolatedSavepoint(self._steps, self.time, self._proc.pid, updates,
                                      self, self.current, journal)
        except Exception:
            status, _, _ = self._command('RELEASE_E', 3)
            if status:
                self._direct_failed = True
                raise RuntimeFailure('failed savepoint preparation could not release its native backup')
            raise
        self._savepoint = token
        return token

    def restore_savepoint(self, token):
        self._token(token)
        if self._restore_pending is None:
            status, _, before = self._command('PROBE_E', 5)
            if status:
                raise RuntimeFailure('cannot inspect native restore counter')
            status, step, updates = self._command('RESTORE_E', 2)
            if status or step != token.step_index or updates != before:
                self._direct_failed = True
                raise RuntimeFailure('native restore failed; no replay attempted')
            self._restore_pending = token
            self._restore_updates = before
        # If truncation fails, preserve the token/phase and forbid stepping until
        # explicit restoration completes. Do not repeat the native copy command.
        token._journal.restore()
        status, step, updates = self._command('RESUME_E', 4)
        if status or step != token.step_index or updates != self._restore_updates:
            self._direct_failed = True
            raise RuntimeFailure('native output reconnection failed; state is not usable')
        self._current = np.array(token._observables, copy=True)
        self._steps = token.step_index
        self._history = self._history[:self._steps]
        self._curr_history = self._curr_history[:self._steps]
        self._restore_pending = None
        self._restore_updates = None
        self._direct_restores += 1

    def release_savepoint(self, token):
        self._token(token)
        if self._restore_pending is not None:
            raise SnapshotError('complete the pending output restore before releasing the savepoint')
        status, step, _ = self._command('RELEASE_E', 3)
        if status or step != self._steps:
            self._direct_failed = True
            raise RuntimeFailure('native savepoint release failed')
        self._savepoint = None

    def step(self, E, *, spatial_command=None):
        if self._direct_failed or self._restore_pending is not None:
            raise RuntimeFailure('isolated worker has an incomplete or failed direct restore')
        if self.execution_mode=='performance' and self._worker_mode!='E3M':return super().step(E)
        if self._savepoint is None and not self._direct_origin and self._worker_mode!='E3M':
            return super().step(E)
        # The live-savepoint route never enters TDDFTState's process/replay recovery.
        E = (_array_unchecked(E,'E',(3,)) if self.execution_mode=='performance' else _real(E, (3,), 'E'))
        if self._steps >= self._capacity:
            raise InputError('isolated worker capacity exhausted')
        if not np.isfinite((self._steps+1)*self.dt):raise InputError('isolated time overflow')
        try:
            self._send(spatial_command if spatial_command is not None else self._step_command+' '+' '.join(f'{v:.17e}' for v in E))
            current = self._receive(self._steps+1)
        except Exception as exc:
            self._direct_failed = True
            if isinstance(exc,NumericalError):
                raise NumericalError('isolated live-savepoint numerical failure; no replay attempted') from exc
            raise RuntimeFailure('isolated live-savepoint step failed; no replay attempted') from exc
        self._current = current
        self._history.append(E.copy())
        self._curr_history.append(current.copy())
        self._steps += 1
        return self.current

    def save_direct_checkpoint(self, directory):
        """Persist Fortran numerical state and output streams without propagation."""
        import shutil
        import tempfile
        if self._savepoint is not None or not self.usable:
            raise SnapshotError('release the live savepoint and complete restoration before disk save')
        destination = Path(directory).expanduser().absolute()
        if destination.exists() or destination.is_symlink():
            raise SnapshotError('checkpoint destination already exists')
        destination.parent.mkdir(parents=True, exist_ok=True)
        run = Path(self._work.name) / 'run'
        binary = run / 'isolated-state.bin'
        if binary.exists():
            binary.unlink()
        status, step, _ = self._command('SAVE_E', 6)
        if status or step != self._steps:
            raise SnapshotError('worker declined direct disk save: '+str(status))
        temporary = Path(tempfile.mkdtemp(prefix='.isolated-checkpoint-', dir=destination.parent))
        try:
            from .native_steps import _AppendJournal
            streams = sorted(_AppendJournal.capture(run).files)
            files = {}
            for name in ['isolated-state.bin', *streams]:
                shutil.copyfile(run/name, temporary/name)
                files[name] = dict(sha256=_digest(temporary/name), bytes=(temporary/name).stat().st_size)
            metadata = dict(schema='salmonpy.isolated-direct/1', settings=self._snapshot_settings(),
                            executable_sha256=self._executable.sha256, step=self._steps,
                            observables=self.current.tolist(), history=[v.tolist() for v in self._history],
                            observable_history=[v.tolist() for v in self._curr_history], files=files)
            (temporary/'manifest.json').write_text(json.dumps(metadata, allow_nan=False, indent=2)+'\n')
            temporary.rename(destination)
        finally:
            if temporary.exists():
                shutil.rmtree(temporary)
            binary.unlink(missing_ok=True)

    def restore_direct_checkpoint(self, directory, *, max_bytes=None):
        """Validate, load in an independent fresh worker, then replace this worker."""
        import shutil
        if self._savepoint is not None or not self.usable:
            raise SnapshotError('release the live savepoint before disk restoration')
        try:
            root = Path(directory).expanduser().resolve(strict=True)
            manifest = root/'manifest.json'
            if manifest.is_symlink() or not manifest.is_file():
                raise ValueError('checkpoint manifest is not a regular file')
            if max_bytes is not None and manifest.stat().st_size > max_bytes:
                raise ValueError('checkpoint exceeds explicit read limit')
            data = json.loads(manifest.read_text())
            if (data['schema'] != 'salmonpy.isolated-direct/1'
                    or data['settings'] != self._snapshot_settings()
                    or data['executable_sha256'] != self._executable.sha256):
                raise ValueError('checkpoint model, runtime or settings differ')
            step = data['step']
            if type(step) is not int or not 0 <= step <= self._capacity:
                raise ValueError('invalid saved step')
            observations = _real(data['observables'], (9,), 'saved observations')
            history = _real(data['history'], (step, 3), 'saved E history') if step else np.empty((0,3))
            currents = _real(data['observable_history'], (step, 9), 'saved observations history') if step else np.empty((0,9))
            files = data['files']
            if not isinstance(files, dict) or 'isolated-state.bin' not in files:
                raise ValueError('missing numerical state')
            if {v.name for v in root.iterdir()} != {'manifest.json', *files}:
                raise ValueError('checkpoint file inventory differs')
            total = manifest.stat().st_size
            for name, record in files.items():
                if (not isinstance(name,str) or Path(name).name != name or name in ('.','..')
                        or not (name == 'isolated-state.bin' or name.endswith('.data'))):
                    raise ValueError('invalid checkpoint filename')
                path = root/name
                if path.is_symlink() or not path.is_file() or path.stat().st_size != record['bytes']:
                    raise ValueError('checkpoint file type or size differs')
                total += path.stat().st_size
            if max_bytes is not None and total > max_bytes:
                raise ValueError('checkpoint exceeds explicit read limit')
            for name, record in files.items():
                if _digest(root/name) != record['sha256']:
                    raise ValueError('checkpoint checksum mismatch')
        except (KeyError, ValueError, TypeError, OSError) as exc:
            raise SnapshotError('invalid isolated direct checkpoint: '+str(exc)) from exc
        candidate = None
        try:
            candidate = _EEngine(self._model, self.dt, self._position,execution_mode=self.execution_mode,
                                 native_near_field=self._worker_mode=="E3N",electric_potential=self._worker_mode=="E3Y",near_field=self._worker_mode=="E3M")
            run = Path(candidate._work.name)/'run'
            shutil.copyfile(root/'isolated-state.bin', run/'isolated-state.bin')
            status, actual_step, updates = candidate._command('LOAD_E', 7)
            if status or actual_step != step or updates != 0:
                raise SnapshotError('worker rejected direct checkpoint; original worker retained')
            # LOAD_E leaves native output writers closed until imported bytes exist.
            expected_streams = {v.name for v in run.glob('*.data')}
            if expected_streams != {name for name in files if name.endswith('.data')}:
                raise SnapshotError('saved output stream inventory differs')
            for name in files:
                if name.endswith('.data'):
                    if not (run/name).is_file() or (run/name).is_symlink():
                        raise SnapshotError('saved output stream differs from initialized streams')
                    shutil.copyfile(root/name, run/name)
            status, actual_step, updates = candidate._command('RELEASE_E', 3)
            if status or actual_step != step or updates != 0:
                raise SnapshotError('restored worker could not reconnect outputs')
            candidate._steps = step
            candidate._current = np.array(observations, copy=True)
            candidate._history = [v.copy() for v in history]
            candidate._curr_history = [v.copy() for v in currents]
            candidate._direct_origin = True
            self.close()
            self.__dict__.update(candidate.__dict__)
            # Transfer ownership; candidate must never close the installed process.
            candidate._closed = True
            candidate._proc = None
            candidate._work = None
            candidate._owner = None
            candidate = None
        finally:
            if candidate is not None:
                candidate.close()

    def snapshot(self):
        if self._worker_mode=='E3M':raise SnapshotError('near-field states use live savepoints or direct checkpoints')
        if self._savepoint is not None:
            raise SnapshotError('release the live savepoint before creating a replay snapshot')
        return super().snapshot()

    def restore_checkpoint(self, directory, *, max_bytes=None):
        if self._savepoint is not None:
            raise SnapshotError('release the live savepoint before replay checkpoint restoration')
        return super().restore_checkpoint(directory, max_bytes=max_bytes)

    def restore(self, snapshot):
        if self._worker_mode=='E3M':raise SnapshotError('near-field states use live savepoints or direct checkpoints')
        if self._savepoint is not None:
            raise SnapshotError('release the live savepoint before replay restoration')
        return super().restore(snapshot)

    def close(self):
        self._savepoint = None
        self._restore_pending = None
        super().close()



@dataclass(frozen=True, slots=True, eq=False)
class IsolatedTDDFTInitialState:
    """Owned native GS snapshot; states lease it and own independent processes."""
    _model: TDDFTInitialState = field(repr=False)
    runtime: IsolatedRuntime
    native_near_field: bool = False
    electric_potential: bool = False
    near_field: bool = False

    def __post_init__(self):
        if type(self.near_field) is not bool:
            raise InputError('near_field must be bool')
        if self.near_field and (self.electric_potential or not ({'near_eh','grid_eh'} & set(json.loads(self.runtime.validation)))):
            raise InputError('select validated near_field separately from electric_potential')
        if self.near_field and _near_feature(self._model._snapshot/'salmonpy-near-eh.dat') not in json.loads(self.runtime.validation):
            raise InputError('validation evidence for this field representation required')
        if type(self.electric_potential) is not bool:
            raise InputError('electric_potential must be bool')
        if self.electric_potential and 'solid_field' not in json.loads(self.runtime.validation):
            raise InputError('spherical-harmonic field validation evidence required')
        if type(self.native_near_field) is not bool:
            raise InputError('native_near_field must be bool')
        if self.native_near_field and 'native_near_field' not in json.loads(self.runtime.validation):
            raise InputError('near-field validation evidence required for this executable')

    @classmethod
    def from_salmon(cls,directory,input_file,*,runtime,max_steps=1000,threads=1,timeout=300.,native_near_field=False,electric_potential=False,near_field=False):
        if not isinstance(runtime,IsolatedRuntime):
            raise InputError('validated IsolatedRuntime required; periodic A runtime is not sufficient')
        _receipt(runtime.validation,runtime.executable)
        if type(native_near_field) is not bool:
            raise InputError('native_near_field must be bool')
        if native_near_field and 'native_near_field' not in json.loads(runtime.validation):
            raise InputError('near-field validation evidence required for this executable')
        if type(electric_potential) is not bool:
            raise InputError('electric_potential must be bool')
        if electric_potential:
            if 'solid_field' not in json.loads(runtime.validation):
                raise InputError('spherical-harmonic field validation evidence required')
            path=Path(directory)/'salmonpy-ylm.dat'
            if not path.is_file() or path.is_symlink():
                raise InputError('prepared salmonpy-ylm.dat file required in the input bundle')
        if type(near_field) is not bool:
            raise InputError('near_field must be bool')
        if near_field:
            if electric_potential or not ({'near_eh','grid_eh'} & set(json.loads(runtime.validation))):
                raise InputError('select validated near_field separately from electric_potential')
            path=Path(directory)/'salmonpy-near-eh.dat'
            if not path.is_file() or path.is_symlink():
                raise InputError('prepared salmonpy-near-eh.dat required in the input bundle')
            feature=_near_feature(path)
            if feature not in json.loads(runtime.validation):raise InputError('validation evidence for this field representation required')

        model=TDDFTInitialState.from_salmon(directory,input_file,executable=runtime.executable,
                                          max_steps=max_steps,threads=threads,timeout=timeout)
        return cls(model,runtime,native_near_field,electric_potential,near_field)

    def new_state(self,*,dt,position=(0.,0.,0.),execution_mode="checked"):
        dt=float(_real(dt,(),'dt'))
        if dt<=0: raise InputError('positive dt required')
        position=_real(position,(3,),'position')
        result=IsolatedTDDFTState(_EEngine(self._model,dt,position,execution_mode=execution_mode,
                                                native_near_field=self.native_near_field,electric_potential=self.electric_potential,near_field=self.near_field),position)
        result._charge_density_enabled='charge_density' in json.loads(self.runtime.validation)
        result._orbital_magnetic_enabled='orbital_magnetic' in json.loads(self.runtime.validation)
        result._grid_stream_enabled='grid_stream' in json.loads(self.runtime.validation)
        return result

    def close(self): self._model.close()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()


@dataclass(frozen=True,slots=True,eq=False)
class ChargeDensitySamples:
    """Read-only electronic charge density (negative electron number density).

    Positions are TDDFT coordinates in bohr, densities e/bohr^3, volumes bohr^3.
    Ionic charges are not included; fixed ions cancel from density differences.
    """
    step_index: int
    time: float
    positions: object
    charge_density: object
    volume_elements: object

    def __post_init__(self):
        if type(self.step_index) is not int or self.step_index<0:raise InputError('nonnegative integer density step required')
        time=float(_real(self.time,(),'density time'))
        if time<0:raise InputError('nonnegative density time required')
        object.__setattr__(self,'time',time)
        points=np.asarray(self.positions)
        if points.ndim!=2 or points.shape[1]!=3 or not len(points):raise InputError('density positions[N,3] required')
        object.__setattr__(self,'positions',_real(points,points.shape,'density positions'))
        object.__setattr__(self,'charge_density',_real(self.charge_density,(len(points),),'signed charge density'))
        volumes=_real(self.volume_elements,(len(points),),'density volume elements')
        if np.any(volumes<=0):raise InputError('positive density volume elements required')
        object.__setattr__(self,'volume_elements',volumes)


@dataclass(frozen=True,slots=True,eq=False)
class OrbitalMagneticResponse:
    """Spinless orbital M=-d<E>/dH on the supplied component H grids.

    Gaussian/atomic units, H=B in vacuum. Does not contain electric polarization
    current or spin magnetization. The Hamiltonian's applied endpoint is used.
    """
    step_index: int
    time: float
    origin: object
    spacing: object
    magnetization: object

    def __post_init__(self):
        _integer(self.step_index,'orbital response step',0,2**31-1)
        time=float(_real(self.time,(),'orbital response time'))
        if time<0:raise InputError('nonnegative response time required')
        object.__setattr__(self,'time',time)
        object.__setattr__(self,'origin',_real(self.origin,(2,3,3),'component origins'))
        spacing=_real(self.spacing,(3,),'spacing')
        if np.any(spacing<=0):raise InputError('positive spacing required')
        object.__setattr__(self,'spacing',spacing)
        M=np.asarray(self.magnetization)
        if M.ndim!=4 or M.shape[-1]!=3 or min(M.shape[:3])<2:raise InputError('magnetization[z,y,x,3] required')
        object.__setattr__(self,'magnetization',_real(M,M.shape,'orbital magnetization'))

    def current_density(self,*,previous=None):
        """Only c*curl(M), on E components. Add once to dP/dt, not a full orbital current.

        With previous, average adjacent accepted magnetizations in Fortran for
        the interval current. Otherwise return the instantaneous endpoint current.
        Requires Yee component origins and a zero positive-face padding layer.
        The patch must cover the molecular support and its surface current.
        """
        from .near_field import _kernel,_ptr
        expected=self.origin[0]+.5*self.spacing-np.diag(self.spacing)
        tolerance=64*np.finfo(float).eps*max(1.,float(np.max(np.abs(expected))))
        if not np.allclose(self.origin[1],expected,rtol=0,atol=tolerance):
            raise InputError('magnetic return requires Yee-aligned E/H component origins')
        M=self.magnetization
        if previous is not None:
            if (not isinstance(previous,OrbitalMagneticResponse) or previous.step_index+1!=self.step_index
                or previous.time>=self.time or previous.magnetization.shape!=M.shape
                or not np.array_equal(previous.origin,self.origin) or not np.array_equal(previous.spacing,self.spacing)):
                raise InputError('adjacent accepted magnetic responses with identical geometry required')
            n=_integer(M.size,'magnetic sample count',1,2**31-1);centered=np.empty_like(M)
            code=_kernel('sp_grid_average')(n,_ptr(previous.magnetization),_ptr(M),_ptr(centered))
            if code:raise NumericalError('magnetic interval averaging rejected: '+str(code))
            M=centered
        counts=np.array(M.shape[2::-1],dtype=np.int32);result=np.empty_like(M)
        status=_kernel('sp_grid_magnetization_current')(_ptr(counts),_ptr(self.spacing),_ptr(M),_ptr(result))
        if status==1:raise InputError('magnetic current requires finite samples and a zero positive-face padding layer')
        if status==4:raise NumericalError('nonfinite orbital magnetic curl')
        if status:raise RuntimeFailure('orbital magnetic curl allocation failed')
        result.setflags(write=False);return result

    def total_current(self,electric_current,*,previous=None):
        """Add c*curl(M) once to the electric polarization current, in Fortran.

        Do not pass a current that already contains orbital magnetic response.
        """
        from .near_field import _kernel,_ptr
        magnetic=self.current_density(previous=previous)
        electric=_real(electric_current,magnetic.shape,'electric polarization current')
        n=_integer(electric.size,'current element count',1,2**31-1);result=np.empty_like(electric)
        status=_kernel('sp_grid_current_sum')(n,_ptr(electric),_ptr(magnetic),_ptr(result))
        if status==4:raise NumericalError('nonfinite combined current')
        if status:raise RuntimeFailure('current combination rejected: '+str(status))
        result.setflags(write=False);return result


def _decode_charge_density(raw,*,step,dt):
    if len(raw)<16:raise RuntimeFailure('truncated density transfer header')
    magic,version,actual,n=np.frombuffer(raw[:16],dtype=np.int32)
    if magic!=16909060 or version!=1 or actual!=step or n<1 or len(raw)!=16+int(n)*40:
        raise RuntimeFailure('density transfer format/step differs')
    data=np.frombuffer(raw,offset=16,dtype=np.float64).reshape(int(n),5)
    return ChargeDensitySamples(step,step*dt,data[:,:3],data[:,3],data[:,4])


class IsolatedTDDFTState:
    """E-driven charge-dipole state with live savepoints and separate replay archives."""
    def __init__(self,engine,position):
        self._engine=engine
        self._contract=FieldContract(TimeSpec(engine.dt,(InputPort('E'),),(),
                            integrator='SALMON isolated middlepoint; E at native step index'),
                            units={'E':'E:au'},position=position,gauge='length')

    @property
    def execution_mode(self):return self._engine.execution_mode
    @property
    def dt(self): return self._engine.dt
    @property
    def time(self): return self._engine.time
    @property
    def timing(self): return self._contract.timing
    @property
    def observables(self):
        values=self._engine.current
        return DipoleObservables._from_worker(self.time,values,checked=self.execution_mode=="checked")

    @property
    def electron_number(self): return self.observables.electron_number

    @property
    def internal_energy(self): return self.observables.internal_energy

    def orbital_magnetic_response(self):
        """Read only -d<E>/dH; no propagation and no electric/spin current."""
        if not getattr(self,'_orbital_magnetic_enabled',False):
            raise InputError('runtime validation does not declare orbital magnetic return')
        return self._engine.orbital_magnetic_response()

    def charge_density(self):
        """Explicit current-step density read; no time evolution or replay."""
        if not getattr(self,'_charge_density_enabled',False):
            raise InputError('runtime validation does not declare charge-density readout')
        return self._engine.charge_density()

    def step(self,E):
        E=(_array_unchecked(E,'E',(3,)) if self.execution_mode=='performance' else _real(E,(3,),'E'))
        self._engine.step(E)
        return self.observables

    def step_fields(self,fields):
        self._contract.validate(fields,step_index=self._engine._steps+1)
        return self.step(fields.E.value)

    @property
    def usable(self): return self._engine.usable

    def worker_statistics(self):
        """Explicit local worker query: PID, actual evolution count and restore count."""
        return self._engine.worker_statistics()

    def savepoint(self):
        """Retain one direct live-worker state; no propagation or GS copying."""
        return self._engine.savepoint()

    def restore_savepoint(self, token):
        """Restore this worker directly; retain token for repeated restore or release."""
        self._engine.restore_savepoint(token)

    def release_savepoint(self, token):
        """Explicitly retire a live savepoint; leave the current state unchanged."""
        self._engine.release_savepoint(token)

    def step_grid_field(self,window,*,step_index,E_extra=(0.,0.,0.)):
        """Advance target step n with grid frames labeled n,n+1.

        window has two relative times 0,dt and time_origin=n*dt. Geometry is
        fixed by the prepared initial grid. Direct saves retain streamed fields.
        """
        from .near_field import GridFieldSeries
        if not getattr(self,'_grid_stream_enabled',False):raise InputError('runtime validation does not declare live grid input')
        if self._engine._near_spec!=('grid',):raise InputError('prepared direct-grid near_field=True state required')
        n=self._engine._steps+1
        step_index=_integer(step_index,'step_index',1,self._engine._capacity)
        if step_index!=n or n>self._engine._capacity:
            raise InputError('live grid target must equal the next molecular step')
        if not isinstance(window,GridFieldSeries) or len(window.times)!=2:
            raise InputError('two-frame GridFieldSeries required')
        if not np.isclose(window.times[1],self.dt,rtol=1e-12,atol=0.) or not np.isclose(window.time_origin,n*self.dt,rtol=1e-12,atol=1e-15):
            raise InputError('live grid frames must carry field-step labels n,n+1 at fixed dt')
        E=_real(E_extra,(3,),'additional uniform E')
        self._engine._stage_grid_field(window,n)
        self._engine.step(E)
        return self.observables

    def step_near_field(self, left, right, *, E_extra=(0.,0.,0.)):
        from .near_field import NearFieldFrame
        spec=self._engine._near_spec
        if spec is None:raise InputError('near_field=True state required')
        if spec[0]=='grid':raise InputError('use step(E_extra) for prepared grids or step_grid_field for live grids, not coefficient frames')
        n=self._engine._steps+1
        for frame,index in ((left,n),(right,n+1)):
            if not isinstance(frame,NearFieldFrame) or frame.step_index!=index:
                raise InputError('spatial fields must carry step indices n and n+1')
            if (frame.lmax!=spec[0] or tuple((*frame.center,frame.radius,*frame.decay_vector))!=spec[1]
                    or not np.isclose(frame.time,index*self.dt,rtol=1.e-12,atol=1.e-15)):
                raise InputError('spatial field center/radius/decay/order/time differs from initialized model')
        E=_real(E_extra,(3,),'additional uniform E')
        values=(*E,n,*left.E_coefficients.ravel(),*right.E_coefficients.ravel(),
                *left.H_coefficients.ravel(),*right.H_coefficients.ravel())
        command='STEP_EH '+' '.join(format(float(v),'.17e') for v in values[:3])+' '+str(n)+' '+' '.join(format(float(v),'.17e') for v in values[4:])
        self._engine.step(E,spatial_command=command)
        return self.observables

    def snapshot(self):
        if self._engine._worker_mode=='E3M':raise SnapshotError('near-field states use live savepoints or direct checkpoints')
        return self._engine.snapshot()
    def restore(self,snapshot):
        if self._engine._worker_mode=='E3M':raise SnapshotError('near-field states use live savepoints or direct checkpoints')
        self._engine.restore(snapshot)
    def save_direct_checkpoint(self, directory):
        """Save numerical state directly; no E-history replay is needed to resume."""
        self._engine.save_direct_checkpoint(directory)

    def restore_direct_checkpoint(self, directory, *, max_bytes=None):
        """Resume in a new process, retaining the original state if import fails."""
        self._engine.restore_direct_checkpoint(directory, max_bytes=max_bytes)

    def save_checkpoint(self,directory):
        """Save accepted E/dipole replay history; GS files remain separate."""
        self._engine.save_checkpoint(directory)

    def restore_checkpoint(self,directory,*,max_bytes=None):
        """Validate matching initial data/settings and atomically replay disk history."""
        self._engine.restore_checkpoint(directory,max_bytes=max_bytes)

    def close(self): self._engine.close()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()


def _caller_mpi():
    import sys
    mpi=sys.modules.get('mpi4py.MPI')
    if mpi is None or not mpi.Is_initialized() or mpi.Is_finalized():
        raise InputError('initialize MPI in the caller before using a molecular MPI group')
    return mpi


def run_molecular_mpi(build_record, solver_comm, control_comm, *, controller_rank,
                      work_directory, input_file='input.inp', decomposition):
    """Run one molecular lifetime collectively on caller-owned MPI ranks.

    Blocking call on dedicated main-thread solver ranks. The external controller
    sends INIT_E, STEP_E, GRID_E, PROBE_E and CLOSE on control_comm. No MPI launcher,
    MPI initialization/finalization or process creation occurs here. Current host
    scope is fixed-ion scalar PZ/direct E/H. Direct disk restore requires an
    unchanged group layout and binary. Deep-error recovery remains separate.
    """
    import ctypes as c
    import os
    import threading
    MPI=_caller_mpi()
    error=None;descriptor=None
    try:
        if not MPI.Is_initialized() or MPI.Is_finalized():raise InputError('caller-owned live MPI required')
        if threading.current_thread() is not threading.main_thread():raise InputError('main thread required')
        rank=_integer(controller_rank,'controller_rank',0,control_comm.size-1)
        if control_comm.rank==rank:raise InputError('controller must be outside the molecular group')
        root=Path(work_directory).resolve(strict=True)
        inp=Path(input_file)
        if inp.is_absolute() or '..' in inp.parts:raise InputError('relative molecular input path required')
        if not (root/inp).is_file():raise InputError('molecular input file missing')
        layout=tuple(_integer(v,'decomposition',1,2**31-1) for v in decomposition)
        if len(layout)!=5 or layout[0]!=1 or np.prod(np.asarray(layout,dtype=object))!=solver_comm.size:
            raise InputError('isolated (k=1,orbital,x,y,z) decomposition must match solver group')
        record=json.loads(Path(build_record).read_text())
        library=Path(record['library']).resolve(strict=True)
        if _digest(library)!=record['sha256']:raise InputError('molecular library checksum mismatch')
        lib=c.CDLL(str(library))
        from .api import query_host_abi,_bind_host_operations
        query_host_abi(lib);_bind_host_operations(lib)
        text=c.create_string_buffer(8192);active=c.c_int()
        lib.sp_library(text,len(text),c.byref(active))
        if not active.value or text.value.decode().strip()!=MPI.Get_library_version().rstrip('\x00').strip():
            raise InputError('calculation library and mpi4py must use the same initialized MPI runtime')
        try:probe=lib.sp_molecular_mpi_protocol
        except AttributeError as exc:
            raise InputError('calculation library lacks the molecular MPI entry; explicitly prepare a matching build') from exc
        probe.argtypes=[];probe.restype=c.c_int
        if probe()!=2:raise InputError('CPU molecular host-MPI protocol 2 required')
        fn=lib.sp_molecular_mpi_run
        fn.restype=c.c_int;fn.argtypes=[c.c_int,c.c_int,c.c_int,c.c_char_p,c.c_void_p]
        descriptor=(str(root),str(inp),layout,record['sha256'],rank)
    except Exception as exc:error=f'{type(exc).__name__}: {exc}'
    admissions=solver_comm.allgather((error,descriptor))
    if any(item[0] for item in admissions):raise InputError(next(item[0] for item in admissions if item[0]))
    if any(item[1]!=descriptor for item in admissions):raise InputError('molecular group configuration differs')
    initial_error=None
    if solver_comm.rank==0:
        try:
            marker=root/'salmonpy-molecular-initial.json'
            if marker.exists():raise InputError('fresh molecular working directory required')
            assets=[]
            for path in sorted(root.rglob('*')):
                if path.is_symlink():raise InputError('molecular input symlinks unsupported')
                if path.is_file():assets.append((path.relative_to(root).as_posix(),_digest(path)))
            identity=_fingerprint((tuple(assets),record['sha256']))
            with marker.open('x') as stream:
                stream.write(json.dumps({'identity':identity,'library_sha256':record['sha256']})+'\n')
        except Exception as exc:initial_error=f'{type(exc).__name__}: {exc}'
    errors=solver_comm.allgather(initial_error)
    if any(errors):raise InputError(next(error for error in errors if error))
    raw=np.ascontiguousarray(layout,dtype=np.int32)
    previous=Path.cwd()
    try:
        os.chdir(root)
        status=fn(solver_comm.py2f(),control_comm.py2f(),rank,os.fsencode(str(inp)),raw.ctypes.data)
    finally:os.chdir(previous)
    statuses=solver_comm.allgather(status)
    if any(statuses):raise RuntimeFailure('molecular host-MPI exit status '+str(statuses))


class MolecularMPIController:
    """Root-only E/H controller for a dedicated run_molecular_mpi group.

    Uses caller-owned MPI, shared working files and blocking point-to-point
    commands. No process launch and no application-wide step collectives.
    The caller owns job termination if a solver process or MPI fails. This is
    not the process-isolated savepoint/recovery interface.
    """
    def __init__(self,control_comm,*,solver_root,work_directory,dt,max_steps,decomposition):
        _caller_mpi()
        self.comm=control_comm
        self.root=_integer(solver_root,'solver_root',0,control_comm.size-1)
        if self.root==control_comm.rank:raise InputError('separate molecular root required')
        self.directory=Path(work_directory).resolve(strict=True)
        self.dt=float(_real(dt,(),'dt'))
        if self.dt<=0:raise InputError('positive dt required')
        self.capacity=_integer(max_steps,'max_steps',1,2**31-2)
        self.decomposition=tuple(_integer(v,'decomposition',1,2**31-1) for v in decomposition)
        if len(self.decomposition)!=5 or self.decomposition[0]!=1:raise InputError('isolated decomposition required')
        self.step_index=0;self.closed=False
        self._send(f'INIT_E {self.dt:.17e} {self.capacity} 0 0 0')
        self.observables=self._receive(0)
        initial=json.loads((self.directory/'salmonpy-molecular-initial.json').read_text())
        self._identity=_fingerprint((initial,self.dt,self.capacity,self.decomposition))

    def _send(self,line):
        MPI=_caller_mpi()
        if self.closed:raise RuntimeFailure('molecular MPI controller is closed')
        self.comm.Send([line.encode('ascii'),MPI.CHAR],dest=self.root,tag=4101)

    def _receive(self,n):
        value=np.empty(10)
        self.comm.Recv(value,source=self.root,tag=4102)
        value=_real(value,(10,),'molecular MPI response')
        if value[0]!=n or value[9] not in (0.,1.):raise RuntimeFailure('molecular MPI response step/flag differs')
        return DipoleObservables._from_worker(n*self.dt,value[1:],checked=True)

    def step_grid_field(self,window,*,step_index,E_extra=(0.,0.,0.)):
        from .near_field import GridFieldSeries
        n=_integer(step_index,'step_index',1,self.capacity)
        if n!=self.step_index+1:raise InputError('next molecular step required')
        if not isinstance(window,GridFieldSeries) or len(window.times)!=2:
            raise InputError('two-frame GridFieldSeries required')
        if not np.isclose(window.times[1],self.dt,rtol=1e-12,atol=0.) or not np.isclose(window.time_origin,n*self.dt,rtol=1e-12,atol=1e-15):
            raise InputError('field frames n,n+1 at the molecular dt required')
        extra=_real(E_extra,(3,),'additional uniform E')
        _write_grid_packet(self.directory/'salmonpy-grid-frames.bin',window,n,self.dt)
        self._send('GRID_E');status=self._status(9)
        if status['status'] or status['step']!=self.step_index:
            raise RuntimeFailure('molecular MPI grid staging failed: '+str(status))
        self._send('STEP_E '+' '.join(f'{v:.17e}' for v in extra))
        result=self._receive(n);self.step_index=n;self.observables=result
        return result

    def _status(self,operation):
        value=np.empty(4,dtype=np.int32)
        self.comm.Recv(value,source=self.root,tag=4103)
        if value[0]!=operation:raise RuntimeFailure('molecular MPI operation reply differs')
        return dict(status=int(value[1]),step=int(value[2]),updates=int(value[3]))

    def worker_statistics(self):
        self._send('PROBE_E')
        result=self._status(5)
        if result.pop('status'):raise RuntimeFailure('molecular MPI counter query failed')
        return result

    def save_direct_checkpoint(self,directory):
        """Save rank-owned electronic payloads plus the root output streams."""
        import shutil
        root=Path(directory).resolve()
        root.mkdir(parents=True,exist_ok=False)
        self._send('SAVE_E');result=self._status(6)
        if result['status'] or result['step']!=self.step_index:
            raise RuntimeFailure('molecular MPI save failed: '+str(result))
        files={}
        ranks=int(np.prod(np.asarray(self.decomposition,dtype=object)))
        # Retire only the transient files owned by this completed SAVE_E.
        for rank in range(ranks):
            source=self.directory/f'isolated-state-{rank:08d}.bin'
            target=root/source.name
            shutil.copy2(source,target);files[target.name]=_digest(target)
        streams={}
        for source in sorted(self.directory.glob('*.data')):
            target=root/source.name
            shutil.copy2(source,target);streams[target.name]=_digest(target)
        metadata={'schema':'salmonpy.molecular-mpi-checkpoint/1','identity':self._identity,
                  'step':self.step_index,'ranks':ranks,'files':files,'streams':streams}
        (root/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
        # The archive remains intact; remove working SAVE files so a later save
        # can use Fortran status=new without overwriting a prior archive.
        for name in files:(self.directory/name).unlink()
        return metadata

    def restore_direct_checkpoint(self,directory):
        """Restore into a fresh group without replay; errors are not group rollback."""
        import shutil
        root=Path(directory).resolve()
        if self.step_index!=0:raise InputError('fresh molecular MPI group required for restore')
        metadata=json.loads((root/'manifest.json').read_text())
        ranks=int(np.prod(np.asarray(self.decomposition,dtype=object)))
        expected={f'isolated-state-{rank:08d}.bin' for rank in range(ranks)}
        if (metadata.get('schema')!='salmonpy.molecular-mpi-checkpoint/1' or
                metadata.get('identity')!=self._identity or metadata.get('ranks')!=ranks or
                set(metadata.get('files',{}))!=expected):
            raise InputError('molecular MPI checkpoint input/binary/layout differs')
        step=_integer(metadata.get('step'),'checkpoint step',0,self.capacity)
        streams=metadata.get('streams',{})
        if not isinstance(streams,dict):raise InputError('molecular stream inventory required')
        for name,sha in {**metadata['files'],**streams}.items():
            if (not isinstance(name,str) or Path(name).name!=name or name in ('.','..') or
                    (name not in expected and not name.endswith('.data'))):
                raise InputError('invalid molecular checkpoint filename')
            path=root/name
            if path.is_symlink() or not path.is_file() or _digest(path)!=sha:
                raise InputError('molecular MPI checkpoint checksum differs')
        for name in expected:shutil.copy2(root/name,self.directory/name)
        self._send('LOAD_E');result=self._status(7)
        if result['status'] or result['step']!=step or result['updates']!=0:
            raise RuntimeFailure('molecular MPI restore failed: '+str(result))
        # LOAD_E closed the fresh output streams; replace the initial prefix
        # before reconnecting the existing owning units for append.
        for name in streams:shutil.copy2(root/name,self.directory/name)
        self._send('RELEASE_E');result=self._status(3)
        if result['status']:raise RuntimeFailure('molecular MPI stream reconnection failed: '+str(result))
        self._send('OBSERVE_E')
        self.observables=self._receive(step);self.step_index=step
        for name in expected:(self.directory/name).unlink()
        return self.observables

    def close(self):
        if not self.closed:self._send('CLOSE');self.closed=True


def _write_grid_packet(path,window,n,dt):
    header=np.array([16909060,1,n,*window.counts],dtype=np.int32)
    geometry=np.array([*window.center,*window.molecule_position,dt,*window.origin.ravel(),*window.spacing])
    frames=np.stack((np.stack((window.E[0],window.H[0])),np.stack((window.E[1],window.H[1]))))
    with Path(path).open('wb') as stream:
        stream.write(header.tobytes());stream.write(geometry.tobytes());stream.write(frames.tobytes())
