"""Process-isolated native SALMON TDDFT, with the common Ac=A/c → J protocol.

A separately built, audited SALMON executable with worker support is required.
The initial bundle contains a valid RT input, data_for_restart and relative assets.
States return the electric (charge) current density j; the worker's SALMON
electronic matter current J is kept internally (``matter_current``), j = -J.
"""
from .errors import InputError, InputTypeError, RuntimeFailure, SnapshotError
from .snapshots import _SnapshotMixin, _fingerprint
from dataclasses import dataclass, field
import hashlib
import json
import os
from pathlib import Path
import queue
import shutil
import subprocess
import tempfile
import threading
import numpy as np
from .fdtd import SalmonExecutable
from .ld import _array, _positive, _array_unchecked, _execution_mode_name
from .sbe import _integer, _frozen


class OutputSaveError(RuntimeFailure):
    """Output export failed; raw worker files remain at recovery_directory."""
    def __init__(self,directory):
        self.recovery_directory=Path(directory)
        super().__init__(f'Output save failed; worker files retained at {directory}')


class _WorkingDirectory:
    # Explicit lifetime lets failed exports survive both close() and GC.
    def __init__(self): self.name=tempfile.mkdtemp(prefix='salmonpy-tddft-worker-')
    def cleanup(self): shutil.rmtree(self.name,ignore_errors=True)


class _BundleStore:
    """One private baseline, retained until loader and every lease release it."""
    def __init__(self):
        self._lock=threading.RLock()
        self._directory=tempfile.TemporaryDirectory(prefix='salmonpy-tddft-initial-')
        self.name=self._directory.name
        self._closed=False;self._users=0

    @property
    def closed(self):
        with self._lock:return self._closed

    def _acquire_locked(self):
        count=self._users+1
        lease=_BundleLease(self)
        self._users=count;lease._released=False
        return lease

    def acquire(self):
        with self._lock:
            if self._closed:raise RuntimeFailure('initial snapshot is closed')
            return self._acquire_locked()

    def _cleanup_unused(self):
        if self._closed and self._users==0:self._directory.cleanup()

    def cleanup(self):
        with self._lock:
            self._closed=True
            self._cleanup_unused()


class _BundleLease:
    """Idempotent lifetime reference; an existing state can fork after loader close."""
    def __init__(self,store):
        self._released=True;self._store=store;self.name=store.name

    def fork(self):
        with self._store._lock:
            if self._released:raise RuntimeFailure('TDDFT bundle lease is closed')
            return self._store._acquire_locked()

    def cleanup(self):
        with self._store._lock:
            if not self._released:
                self._released=True;self._store._users-=1
            self._store._cleanup_unused()

    def __del__(self):
        try:
            if not getattr(self,'_released',True):self.cleanup()
        except Exception:pass


def _digest(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()


def _copy_bundle(source,dest):
    if dest.resolve().is_relative_to(source.resolve()):
        raise InputError('snapshot destination must be outside the input bundle')
    dest.mkdir()
    records=[]
    for path in sorted(source.rglob('*')):
        if path.is_symlink(): raise InputError('bundle symlinks are unsupported; supply actual files')
        target=dest/path.relative_to(source)
        if path.is_dir(): target.mkdir(exist_ok=True)
        elif path.is_file():
            shutil.copy2(path,target)
            records.append((str(path.relative_to(source)),_digest(target)))
        else: raise InputError('bundle contains a non-regular file')
    return tuple(records)


@dataclass(frozen=True,slots=True,eq=False)
class TDDFTInitialState:
    """Owned immutable snapshot; loading native binary records is done by SALMON.

    Use from_salmon, and close the snapshot when no further new states are needed.
    States lease this baseline and own independent writable worker copies.
    Existing states survive loader close; the final release removes the baseline.
    """
    executable: SalmonExecutable
    max_steps: int
    threads: int
    timeout: float
    manifest: tuple
    _input: str=field(repr=False)
    _snapshot: Path=field(repr=False)
    _owner: object=field(repr=False)
    backend: str=field(default='SALMON-TDDFT-middlepoint',init=False)
    _source_root: Path | None = field(default=None,repr=False)
    _record_environment: tuple | None = field(default=None,repr=False)

    @classmethod
    def from_salmon(cls,directory,input_file,*,executable,max_steps=1000,threads=1,timeout=300.):
        if not isinstance(executable,SalmonExecutable): raise InputTypeError('expected SalmonExecutable')
        max_steps=_integer(max_steps,'max_steps',1,2**31-2)
        threads=_integer(threads,'threads',1,2**31-1);timeout=_positive(timeout,'timeout')
        root=Path(directory).expanduser().resolve(); inp=Path(input_file)
        if inp.is_absolute() or '..' in inp.parts: raise InputError('input_file must be relative to the bundle')
        for name in ('info.bin','occupation.bin','wfn.bin'):
            if not (root/'data_for_restart'/name).is_file(): raise InputError(f'missing GS record: {name}')
        text=(root/inp).read_text()
        if '!SALMONPY_END_INPUT' in text or any(len(line)>256 for line in text.splitlines()):
            raise InputError('input contains worker marker or a line longer than SALMON accepts')
        owner=_BundleStore()
        snapshot=Path(owner.name)/'bundle'
        try:
            manifest=_copy_bundle(root,snapshot)
            # Read the copied input so its bytes and the manifest describe one snapshot.
            text=(snapshot/inp).read_text()
            from .gs_native import _receipt_environment
            record_environment=_receipt_environment(snapshot,executable)
            return cls(executable,max_steps,threads,timeout,manifest,text,snapshot,owner,
                       _source_root=root,_record_environment=record_environment)
        except BaseException:
            owner.cleanup();raise

    def new_state(self,*,dt,time0=0.,initial_A=(0.,0.,0.),execution_mode="checked"):
        dt=_positive(dt,'dt');time0=float(_array(time0,'time0',()))
        A=_array(initial_A,'initial_A',(3,))
        if time0!=0.: raise InputError('fresh TDDFT GS states start at time0=0')
        if self._owner.closed or not self._snapshot.exists(): raise RuntimeFailure('initial snapshot is closed')
        return TDDFTState(self,dt,A,execution_mode=execution_mode)

    def close(self): self._owner.cleanup()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()


def _read_output(stream,events,log,prefix='SALMONPY_TDDFT_V1 '):
    try:
        with log.open('w') as output:
            for line in stream:
                output.write(line);output.flush()
                if line.startswith(prefix): events.put(line)
    finally:
        stream.close();events.put(None)


class TDDFTState(_SnapshotMixin):
    """One live SALMON process. Use a context manager or close explicitly.

    Only step(Anew) mutates the calculation. In checked mode, on worker failure the accepted
    A history is replayed into a new process before the error is returned. If
    replay itself fails the state is closed and cannot silently continue.
    """
    _worker_mode = '1'
    _reply_prefix = 'SALMONPY_TDDFT_V1 '
    _init_command = 'INIT'
    _performance_init_command = 'INIT_PERFORMANCE'
    _performance_reply_prefix = 'SALMONPY_TDDFT_PERFORMANCE_V1 '
    _step_command = 'STEP'
    _reply_size = 3

    def __init__(self,initial,dt,A,*,execution_mode="checked"):
        self._execution_mode=_execution_mode_name(execution_mode)
        if execution_mode=='performance':
            self._init_command=self._performance_init_command
            self._reply_prefix=self._performance_reply_prefix
        self._model=initial
        self._dt=dt;self._steps=0;self._history=[];self._A0=A.copy()
        self._executable=initial.executable;self._capacity=initial.max_steps
        self._threads=initial.threads;self._timeout=initial.timeout;self._input=initial._input
        self._closed=False;self._proc=None;self._work=None;self._reader=None
        self._retain_outputs=False
        self._initial_hashes=dict(initial.manifest);self._curr_history=[]
        self._owner=None;self._bundle=initial._snapshot
        try:
            self._owner=initial._owner.acquire()
            self._spawn()
            self._J0=self._current.copy()
        except BaseException:
            self.close();raise

    @property
    def execution_mode(self):return self._execution_mode
    @property
    def dt(self): return self._dt
    @property
    def timing(self):
        """salmonpy.time/1 declaration (crystal contract; SALMON middlepoint)."""
        from .timing import TimeSpec,InputPort,OutputPort
        return TimeSpec(self._dt,(InputPort('A'),),(OutputPort('j'),),
            integrator='SALMON middlepoint: H at (A_n+A_{n+1})/2; J at A_{n+1}')
    @property
    def time(self): return self._steps*self.dt
    @property
    def current(self):
        """Electric (charge) current density j = -J_matter at the current time."""
        return _frozen(-self._current,check=self.execution_mode=="checked")
    @property
    def electron_number(self):
        """Actual Fortran-integrated occupied electron density at this endpoint."""
        if self._closed: raise RuntimeFailure("TDDFT state is closed")
        value=getattr(self,"_electron_number",None)
        if value is None: raise RuntimeFailure("runtime does not expose measured electron number; explicitly build/use a matching runtime")
        return value

    @property
    def matter_current(self):
        """SALMON electronic matter current J (internal representation)."""
        return _frozen(self._current,check=self.execution_mode=="checked")
    @property
    def max_steps(self): return self._capacity

    def _spawn(self):
        # A process worker is an independent single-rank MPI application, not
        # another member of the Python host's communicator. Do not pass parent
        # runtime connection descriptors/rank identities into its MPI_Init.
        environment={k:v for k,v in os.environ.items()
            if not k.startswith(('OMPI_','PMI_','PMIX_','PRTE_'))
            and k not in ('MPI_LOCALRANKID','MPI_LOCALNRANKS',
                          'SLURM_PROCID','SLURM_LOCALID')}
        expected=self._model._record_environment
        if expected is not None and any(environment.get(k,'')!=v for k,v in expected):
            raise InputError('native GS record I/O environment differs from conversion')
        if _digest(self._executable.path)!=self._executable.sha256:
            raise InputError('executable SHA256 differs from the audited binary')
        self._work=_WorkingDirectory()
        work=Path(self._work.name)/'run'
        if _copy_bundle(self._bundle,work)!=tuple(sorted(self._initial_hashes.items())):
            raise InputError('initial snapshot checksum mismatch')
        self._events=queue.Queue()
        self._proc=subprocess.Popen([str(self._executable.path)],cwd=work,
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1,
            env={**environment,'SALMONPY_TDDFT_WORKER':self._worker_mode,'SALMONPY_TDDFT_NORM':'1','SALMONPY_TDDFT_EXECUTION_MODE':self.execution_mode,'OMP_NUM_THREADS':str(self._threads)})
        self._reader=threading.Thread(target=_read_output,args=(self._proc.stdout,self._events,work/'worker.log',self._reply_prefix),daemon=True)
        self._reader.start()
        self._send(self._input.rstrip()+'\n!SALMONPY_END_INPUT\n'+
                   f'{self._init_command} {self.dt:.17e} {self._capacity} '+' '.join(f'{v:.17e}' for v in self._A0))
        self._current=self._receive(0)

    def _send(self,text):
        # A monitored writer works with pipes on Windows as well as POSIX.
        # Capture this process: recovery must never redirect an old writer.
        proc = self._proc
        done = queue.Queue()
        def write():
            try:
                proc.stdin.write(text+'\n')
                proc.stdin.flush()
                done.put(None)
            except BaseException as error:
                done.put(error)
        writer = threading.Thread(target=write, daemon=True)
        writer.start()
        try:
            error = done.get(timeout=self._timeout)
        except queue.Empty as exc:
            # Do not send CLOSE to a worker whose input pipe is already blocked.
            if proc.poll() is None:
                proc.kill()
            proc.wait()
            writer.join()
            raise TimeoutError('SALMON TDDFT worker input timed out') from exc
        writer.join()
        if error is not None:
            raise error

    def _receive(self,index):
        try: line=self._events.get(timeout=self._timeout)
        except queue.Empty as exc: raise TimeoutError('SALMON TDDFT worker response timed out') from exc
        if line is None:
            log=Path(self._work.name)/'run/worker.log'
            message=log.read_text()[-6000:]
            from .errors import NumericalError
            if 'SALMONPY_NUMERICAL_FAILURE' in message: raise NumericalError(message)
            raise RuntimeFailure('SALMON TDDFT worker exited before reply:\n'+message)
        fields=line.split()
        self._electron_number=None
        if self._reply_size==3:
            if len(fields)!=6:raise RuntimeFailure('TDDFT worker protocol requires measured electron number')
            self._electron_number=float(fields[-1]);fields=fields[:-1]
        if len(fields)!=2+self._reply_size or int(fields[1])!=index: raise RuntimeFailure('TDDFT worker protocol sequence mismatch')
        parser=_array if self.execution_mode=='checked' else _array_unchecked
        return parser([float(x) for x in fields[2:]],'worker observables',(self._reply_size,))

    def _stop(self,*,cleanup=True):
        returncode=None
        if self._proc is not None:
            proc=self._proc
            if proc.poll() is None:
                try:
                    self._send('CLOSE');proc.wait(timeout=min(self._timeout,5.))
                except (OSError,TimeoutError,subprocess.TimeoutExpired):
                    proc.kill();proc.wait()
            if proc.stdin:
                try: proc.stdin.close()
                except OSError: pass
            if self._reader: self._reader.join(timeout=5.)
            returncode=proc.returncode
            self._proc=None
        if cleanup and self._work is not None and not self._retain_outputs: self._work.cleanup();self._work=None
        return returncode

    def _step_performance(self,Anew):
        if self._closed:raise RuntimeFailure('TDDFT state is closed')
        A=_array_unchecked(Anew,'Anew',(3,))
        if self._steps>=self._capacity:raise InputError('TDDFT max_steps reached')
        try:
            self._send(self._step_command+' '+' '.join(f'{v:.17e}' for v in A))
            current=self._receive(self._steps+1)
        except Exception as exc:
            self.close()
            raise RuntimeFailure('performance TDDFT worker failed; state closed without automatic replay') from exc
        self._current=current;self._history.append(A);self._curr_history.append(current.copy());self._steps+=1
        return self.current

    def step(self,Anew):
        if self.execution_mode=="performance":return self._step_performance(Anew)
        if self._closed: raise RuntimeFailure('TDDFT state is closed')
        A=_array(Anew,'Anew',(3,))
        if self._steps>=self._capacity: raise InputError('TDDFT max_steps reached; create a new state with larger capacity')
        if not np.isfinite((self._steps+1)*self.dt): raise InputError('time overflow')
        previous=self._current.copy()
        try:
            self._send(self._step_command+' '+' '.join(f'{v:.17e}' for v in A))
            current=self._receive(self._steps+1)
        except Exception as exc:
            self._stop()
            try:
                self._spawn()
                for index,oldA in enumerate(self._history,1):
                    self._send(self._step_command+' '+' '.join(f'{v:.17e}' for v in oldA))
                    self._current=self._receive(index)
                if not np.allclose(self._current,previous,rtol=1e-12,atol=1e-14):
                    raise RuntimeFailure('TDDFT replay differs from the accepted state')
                self._current=previous
            except Exception as recovery:
                self.close()
                raise RuntimeFailure('TDDFT worker failed and recovery failed; state closed') from recovery
            from .errors import NumericalError
            if isinstance(exc,NumericalError): raise NumericalError('TDDFT numerical step failed; previous accepted state restored') from exc
            raise RuntimeFailure('TDDFT step failed; previous accepted state restored') from exc
        self._current=current;self._history.append(A.copy());self._curr_history.append(current.copy());self._steps+=1
        return self.current

    _snapshot_kind = 'TDDFT-process-replay'
    _snapshot_counter = '_steps'

    def _snapshot_owner(self): return self._model

    def _snapshot_settings(self):
        return _fingerprint((self.dt, self._capacity, self._threads,
                             self._executable.sha256, self._input,
                             tuple(sorted(self._initial_hashes.items()))), (self._A0,))

    def _snapshot_arrays(self):
        return {'_A0': self._A0, '_J0': self._J0, '_current': self._current,
                '_history': np.asarray(self._history, dtype=float).reshape(-1, 3),
                '_curr_history': np.asarray(self._curr_history, dtype=float).reshape(-1, self._reply_size)}

    def snapshot(self):
        """Save immutable accepted A/current history for in-memory replay.

        No worker is advanced. This is not a native wavefunction checkpoint.
        """
        if self._closed: raise SnapshotError('TDDFT state is closed')
        return super().snapshot()

    def save_checkpoint(self, directory):
        """Write accepted replay history to a new Zarr directory, without GS data."""
        from ._tddft_checkpoint import save
        save(self, directory)

    def restore_checkpoint(self, directory, *, max_bytes=None):
        """Validate a disk checkpoint and replay atomically in this supplied model.

        Identical model contents, executable, dt, initial A, threads and capacity
        are required; model object identity may differ across Python processes.
        Initial GS files must be supplied separately. Cost grows with history.
        The default file bound follows configured capacity and archive layout.
        max_bytes adds a caller budget for files/decoded payload, not total RSS.
        """
        if self._closed: raise SnapshotError('TDDFT state is closed')
        from ._tddft_checkpoint import load
        self.restore(load(self, directory, max_bytes))

    @staticmethod
    def _retire_replay_worker(state):
        """Best-effort retirement cannot turn a committed restore into failure.

        Keep diagnostics on the surviving state if ordinary close failed. The
        captured handles allow cleanup even after close partially changed state.
        """
        proc, work, owner = state._proc, state._work, state._owner
        errors = []
        try:
            state.close()
            return ()
        except Exception as exc:
            errors.append(repr(exc))
        actions = []
        if proc is not None:
            def terminate():
                if proc.poll() is None: proc.kill()
                proc.wait(timeout=5.)
            actions.append(terminate)
            if proc.stdin is not None: actions.append(proc.stdin.close)
        if state._reader is not None:
            actions.append(lambda: state._reader.join(timeout=5.))
        if proc is not None and proc.stdout is not None: actions.append(proc.stdout.close)
        if work is not None: actions.append(work.cleanup)
        if owner is not None: actions.append(owner.cleanup)
        for action in actions:
            try: action()
            except Exception as exc: errors.append(repr(exc))
        state._proc = state._work = state._owner = state._reader = None
        state._closed = True
        return tuple(errors)

    def restore(self, snapshot):
        """Replay a token in a new worker, then replace this accepted branch.

        Requires the same initial model object, dt, initial A, capacity, threads
        and executable. Cost grows with saved step count; both workers coexist
        during validation. Replay currents must match within rtol=1e-12,
        atol=1e-14. A failed replay leaves this worker and history untouched.
        """
        from types import SimpleNamespace
        from .snapshots import StateSnapshot
        if self._closed: raise SnapshotError('TDDFT state is closed')
        if type(snapshot) is not StateSnapshot:
            raise SnapshotError('expected StateSnapshot from state.snapshot()')
        n = snapshot.step_index
        if type(n) is not int or not 0 <= n <= self._capacity:
            raise SnapshotError('invalid snapshot step index')
        # Reuse the common schema/integrity checks, allowing variable history size.
        expected = {'_A0': np.empty(3), '_J0': np.empty(self._reply_size),
                    '_current': np.empty(self._reply_size), '_history': np.empty((n,3)),
                    '_curr_history': np.empty((n,self._reply_size))}
        validated = SimpleNamespace(dt=self.dt, _snapshot_kind=self._snapshot_kind,
            _snapshot_counter='_steps', _snapshot_capacity=self._capacity,
            _snapshot_owner=self._snapshot_owner, _snapshot_settings=self._snapshot_settings,
            _snapshot_arrays=lambda: expected)
        _SnapshotMixin.restore(validated, snapshot)
        if (not np.array_equal(validated._A0, self._A0)
                or not np.allclose(validated._J0, self._J0, rtol=1e-12, atol=1e-14)
                or not np.array_equal(validated._current,
                                      validated._curr_history[-1] if n else validated._J0)):
            raise SnapshotError('inconsistent TDDFT snapshot history')
        candidate = object.__new__(type(self))
        candidate.__dict__ = self.__dict__.copy()
        candidate._proc = candidate._work = candidate._reader = candidate._owner = None
        candidate._retain_outputs = False
        try:
            candidate._owner = self._owner.fork()
            candidate._bundle = self._bundle
            candidate._spawn()
            if not np.allclose(candidate._current, validated._J0, rtol=1e-12, atol=1e-14):
                raise SnapshotError('initial replay current differs')
            candidate._J0 = candidate._current.copy()
            replayed_currents = []
            for index, (A, J) in enumerate(zip(validated._history, validated._curr_history), 1):
                candidate._send(candidate._step_command+' '+' '.join(format(v, '.17e') for v in A))
                candidate._current = candidate._receive(index)
                if not np.allclose(candidate._current, J, rtol=1e-12, atol=1e-14):
                    raise SnapshotError(f'replay current differs at step {index}')
                replayed_currents.append(candidate._current.copy())
            candidate._steps = n
            candidate._history = [a.copy() for a in validated._history]
            candidate._curr_history = replayed_currents
            # Keep the replayed native endpoint current, rather than hiding drift.
        except BaseException as exc:
            self._restore_cleanup_errors = self._retire_replay_worker(candidate)
            if not isinstance(exc, Exception): raise
            raise SnapshotError('TDDFT replay failed; original state preserved') from exc
        # Reader threads capture their state's queue/log. Transfer them without
        # rebinding their state object: each reader was started with explicit queue/log arguments.
        self.__dict__, candidate.__dict__ = candidate.__dict__, self.__dict__
        self._restore_cleanup_errors = self._retire_replay_worker(candidate)

    def finish(self,directory):
        """Close cleanly and archive accepted A/J history and generated native files.

        Spatial outputs/strides are selected in the initial input's &analysis.
        Unchanged GS restart files are indexed by hash, not copied. Other files
        are retained, including outputs already present in the initial bundle.
        This archive is for postprocessing, not a restart checkpoint.
        """
        if self._closed: raise RuntimeFailure('TDDFT state is closed')
        dest=Path(directory).expanduser().resolve()
        if dest.exists(): raise FileExistsError(dest)
        for parent in (Path(self._work.name),Path(self._owner.name)):
            if dest.is_relative_to(parent): raise InputError('archive must be outside temporary state directories')
        from .outputs import _write_archive
        try:
            if self._stop(cleanup=False)!=0: raise RuntimeFailure('worker did not exit cleanly; output not finalized')
            metadata={'backend':'SALMON-TDDFT-middlepoint','execution_mode':self.execution_mode,'dt':self.dt,'steps':self._steps,
                'units':{'time':'a.u.','A':'Ac=A/c, a.u.','current':'electric (charge) current density, a.u.'},
                'time_location':'A and current are accepted endpoint samples; first sample is the initial state',
                'matter_current':'-current','executable_sha256':self._executable.sha256,'threads':self._threads,
                'initial_files':self._initial_hashes,'native_output_semantics':'original SALMON headers/input; no inferred units or time shifts',
                'restart_checkpoint':False}
            return _write_archive(dest,{'time':np.arange(self._steps+1)*self.dt,
                'A':np.array([self._A0]+self._history),'current':-np.array([self._J0]+self._curr_history)},
                metadata,native_source=Path(self._work.name)/'run',initial_hashes=self._initial_hashes,input_text=self._input)
        except BaseException as exc:
            self._retain_outputs=True
            recovery=Path(self._work.name)
            # Raw outputs must survive even if the filesystem also rejects these
            # small recovery records. Never delete the only copy after save failure.
            try:
                np.savez_compressed(recovery/'accepted-response.npz',time=np.arange(self._steps+1)*self.dt,
                    A=np.array([self._A0]+self._history),J_matter=np.array([self._J0]+self._curr_history))
                (recovery/'recovery.json').write_text(json.dumps({'error':str(exc),
                    'executable_sha256':self._executable.sha256,'steps':self._steps,
                    'dt':self.dt,'initial_files':self._initial_hashes,'input':self._input},indent=2))
            except OSError:
                pass
            if not isinstance(exc,Exception):
                exc.recovery_directory=recovery
                if hasattr(exc,'add_note'): exc.add_note(f'Worker output files retained at {recovery}')
                raise
            raise OutputSaveError(recovery) from exc
        finally:
            self.close()

    def close(self):
        self._closed=True
        try: self._stop()
        finally:
            if getattr(self,'_owner',None) is not None: self._owner.cleanup()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()
    def __del__(self):
        try: self.close()
        except Exception: pass
