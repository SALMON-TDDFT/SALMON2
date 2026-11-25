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
from .sbe import _frozen

__all__ = ['IsolatedRuntime','load_isolated_runtime','DipoleObservables',
           'IsolatedTDDFTInitialState','IsolatedTDDFTState','IsolatedSavepoint']
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
        return IsolatedTDDFTState(_EEngine(self._model,dt,position,execution_mode=execution_mode,
                                                native_near_field=self.native_near_field,electric_potential=self.electric_potential,near_field=self.near_field),position)

    def close(self): self._model.close()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()


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

    def step_near_field(self, left, right, *, E_extra=(0.,0.,0.)):
        from .near_field import NearFieldFrame
        spec=self._engine._near_spec
        if spec is None:raise InputError('near_field=True state required')
        if spec[0]=='grid':raise InputError('prepared grid histories advance with step(E_extra), not coefficient frames')
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
