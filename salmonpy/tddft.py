"""Process-isolated native SALMON TDDFT, with the common Ac=A/c → J protocol.

A separately built, audited SALMON executable with worker support is required.
The initial bundle contains a valid RT input, data_for_restart and relative assets.
States return the electric (charge) current density j; the worker's SALMON
electronic matter current J is kept internally (``matter_current``), j = -J.
"""
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
from .ld import _array, _positive
from .sbe import _integer, _frozen


class OutputSaveError(RuntimeError):
    """Output export failed; raw worker files remain at recovery_directory."""
    def __init__(self,directory):
        self.recovery_directory=Path(directory)
        super().__init__(f'Output save failed; worker files retained at {directory}')


class _WorkingDirectory:
    # Explicit lifetime lets failed exports survive both close() and GC.
    def __init__(self): self.name=tempfile.mkdtemp(prefix='salmonpy-tddft-worker-')
    def cleanup(self): shutil.rmtree(self.name,ignore_errors=True)


def _digest(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()


def _copy_bundle(source,dest):
    if dest.resolve().is_relative_to(source.resolve()):
        raise ValueError('snapshot destination must be outside the input bundle')
    dest.mkdir()
    records=[]
    for path in sorted(source.rglob('*')):
        if path.is_symlink(): raise ValueError('bundle symlinks are unsupported; supply actual files')
        target=dest/path.relative_to(source)
        if path.is_dir(): target.mkdir(exist_ok=True)
        elif path.is_file():
            shutil.copy2(path,target)
            records.append((str(path.relative_to(source)),_digest(target)))
        else: raise ValueError('bundle contains a non-regular file')
    return tuple(records)


@dataclass(frozen=True,slots=True,eq=False)
class TDDFTInitialState:
    """Owned immutable snapshot; loading native binary records is done by SALMON.

    Use from_salmon, and close the snapshot when no further new states are needed.
    Each existing state owns another independent copy and survives snapshot close.
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

    @classmethod
    def from_salmon(cls,directory,input_file,*,executable,max_steps=1000,threads=1,timeout=300.):
        if not isinstance(executable,SalmonExecutable): raise TypeError('expected SalmonExecutable')
        max_steps=_integer(max_steps,'max_steps',1,2**31-2)
        threads=_integer(threads,'threads',1,2**31-1);timeout=_positive(timeout,'timeout')
        root=Path(directory).expanduser().resolve(); inp=Path(input_file)
        if inp.is_absolute() or '..' in inp.parts: raise ValueError('input_file must be relative to the bundle')
        for name in ('info.bin','occupation.bin','wfn.bin'):
            if not (root/'data_for_restart'/name).is_file(): raise ValueError(f'missing GS record: {name}')
        text=(root/inp).read_text()
        if '!SALMONPY_END_INPUT' in text or any(len(line)>256 for line in text.splitlines()):
            raise ValueError('input contains worker marker or a line longer than SALMON accepts')
        owner=tempfile.TemporaryDirectory(prefix='salmonpy-tddft-initial-')
        snapshot=Path(owner.name)/'bundle'
        try:
            manifest=_copy_bundle(root,snapshot)
            # Read the copied input so its bytes and the manifest describe one snapshot.
            text=(snapshot/inp).read_text()
            return cls(executable,max_steps,threads,timeout,manifest,text,snapshot,owner)
        except BaseException:
            owner.cleanup();raise

    def new_state(self,*,dt,time0=0.,initial_A=(0.,0.,0.)):
        dt=_positive(dt,'dt');time0=float(_array(time0,'time0',()))
        A=_array(initial_A,'initial_A',(3,))
        if time0!=0.: raise ValueError('fresh TDDFT GS states start at time0=0')
        if not self._snapshot.exists(): raise RuntimeError('initial snapshot is closed')
        return TDDFTState(self,dt,A)

    def close(self): self._owner.cleanup()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()


def _read_output(stream,events,log):
    try:
        with log.open('w') as output:
            for line in stream:
                output.write(line);output.flush()
                if line.startswith('SALMONPY_TDDFT_V1 '): events.put(line)
    finally:
        stream.close();events.put(None)


class TDDFTState:
    """One live SALMON process. Use a context manager or close explicitly.

    Only step(Anew) mutates the calculation. On a worker failure the accepted
    A history is replayed into a new process before the error is returned. If
    replay itself fails the state is closed and cannot silently continue.
    """
    def __init__(self,initial,dt,A):
        self._dt=dt;self._steps=0;self._history=[];self._A0=A.copy()
        self._executable=initial.executable;self._capacity=initial.max_steps
        self._threads=initial.threads;self._timeout=initial.timeout;self._input=initial._input
        self._closed=False;self._proc=None;self._work=None;self._reader=None
        self._retain_outputs=False
        self._initial_hashes=dict(initial.manifest);self._curr_history=[]
        self._owner=tempfile.TemporaryDirectory(prefix='salmonpy-tddft-state-')
        self._bundle=Path(self._owner.name)/'bundle'
        try:
            if _copy_bundle(initial._snapshot,self._bundle)!=initial.manifest:
                raise ValueError('initial snapshot checksum mismatch')
            self._spawn()
            self._J0=self._current.copy()
        except BaseException:
            self.close();raise

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
        return _frozen(-self._current)
    @property
    def matter_current(self):
        """SALMON electronic matter current J (internal representation)."""
        return _frozen(self._current)
    @property
    def max_steps(self): return self._capacity

    def _spawn(self):
        if _digest(self._executable.path)!=self._executable.sha256:
            raise ValueError('executable SHA256 differs from the audited binary')
        self._work=_WorkingDirectory()
        work=Path(self._work.name)/'run';_copy_bundle(self._bundle,work)
        self._events=queue.Queue()
        self._proc=subprocess.Popen([str(self._executable.path)],cwd=work,
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1,
            env={**os.environ,'SALMONPY_TDDFT_WORKER':'1','OMP_NUM_THREADS':str(self._threads)})
        self._reader=threading.Thread(target=_read_output,args=(self._proc.stdout,self._events,work/'worker.log'),daemon=True)
        self._reader.start()
        self._send(self._input.rstrip()+'\n!SALMONPY_END_INPUT\n'+
                   f'INIT {self.dt:.17e} {self._capacity} '+' '.join(f'{v:.17e}' for v in self._A0))
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
            raise RuntimeError('SALMON TDDFT worker exited before reply:\n'+log.read_text()[-6000:])
        fields=line.split()
        if len(fields)!=5 or int(fields[1])!=index: raise RuntimeError('TDDFT worker protocol sequence mismatch')
        return _array([float(x) for x in fields[2:]],'worker current',(3,))

    def _stop(self,*,cleanup=True):
        returncode=None
        if self._proc is not None:
            proc=self._proc
            if proc.poll() is None:
                try:
                    self._send('CLOSE');proc.wait(timeout=min(self._timeout,5.))
                except (OSError,subprocess.TimeoutExpired):
                    proc.kill();proc.wait()
            if proc.stdin:
                try: proc.stdin.close()
                except OSError: pass
            if self._reader: self._reader.join(timeout=5.)
            returncode=proc.returncode
            self._proc=None
        if cleanup and self._work is not None and not self._retain_outputs: self._work.cleanup();self._work=None
        return returncode

    def step(self,Anew):
        if self._closed: raise RuntimeError('TDDFT state is closed')
        A=_array(Anew,'Anew',(3,))
        if self._steps>=self._capacity: raise ValueError('TDDFT max_steps reached; create a new state with larger capacity')
        if not np.isfinite((self._steps+1)*self.dt): raise ValueError('time overflow')
        previous=self._current.copy()
        try:
            self._send('STEP '+' '.join(f'{v:.17e}' for v in A))
            current=self._receive(self._steps+1)
        except Exception as exc:
            self._stop()
            try:
                self._spawn()
                for index,oldA in enumerate(self._history,1):
                    self._send('STEP '+' '.join(f'{v:.17e}' for v in oldA))
                    self._current=self._receive(index)
                if not np.allclose(self._current,previous,rtol=1e-12,atol=1e-14):
                    raise RuntimeError('TDDFT replay differs from the accepted state')
                self._current=previous
            except Exception as recovery:
                self.close()
                raise RuntimeError('TDDFT worker failed and recovery failed; state closed') from recovery
            raise RuntimeError('TDDFT step failed; previous accepted state restored') from exc
        self._current=current;self._history.append(A.copy());self._curr_history.append(current.copy());self._steps+=1
        return self.current

    def finish(self,directory):
        """Close cleanly and archive accepted A/J history and generated native files.

        Spatial outputs/strides are selected in the initial input's &analysis.
        Unchanged GS restart files are indexed by hash, not copied. Other files
        are retained, including outputs already present in the initial bundle.
        This archive is for postprocessing, not a restart checkpoint.
        """
        if self._closed: raise RuntimeError('TDDFT state is closed')
        dest=Path(directory).expanduser().resolve()
        if dest.exists(): raise FileExistsError(dest)
        for parent in (Path(self._work.name),Path(self._owner.name)):
            if dest.is_relative_to(parent): raise ValueError('archive must be outside temporary state directories')
        from .outputs import _write_archive
        try:
            if self._stop(cleanup=False)!=0: raise RuntimeError('worker did not exit cleanly; output not finalized')
            metadata={'backend':'SALMON-TDDFT-middlepoint','dt':self.dt,'steps':self._steps,
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
            if hasattr(self,'_owner'): self._owner.cleanup()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()
    def __del__(self):
        try: self.close()
        except Exception: pass
