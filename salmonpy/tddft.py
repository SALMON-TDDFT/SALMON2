"""Process-isolated native SALMON TDDFT, with the common Ac=A/c → J protocol.

A separately built, audited SALMON executable with worker support is required.
The initial bundle contains a valid RT input, data_for_restart and relative assets.
J is SALMON's electronic matter current (charge current is its negative).
"""
from dataclasses import dataclass, field
import hashlib
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


def _digest(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()


def _copy_bundle(source,dest):
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
        self._owner=tempfile.TemporaryDirectory(prefix='salmonpy-tddft-state-')
        self._bundle=Path(self._owner.name)/'bundle'
        try:
            if _copy_bundle(initial._snapshot,self._bundle)!=initial.manifest:
                raise ValueError('initial snapshot checksum mismatch')
            self._spawn()
        except BaseException:
            self.close();raise

    @property
    def dt(self): return self._dt
    @property
    def time(self): return self._steps*self.dt
    @property
    def current(self): return _frozen(self._current)
    @property
    def max_steps(self): return self._capacity

    def _spawn(self):
        if _digest(self._executable.path)!=self._executable.sha256:
            raise ValueError('executable SHA256 differs from the audited binary')
        self._work=tempfile.TemporaryDirectory(prefix='salmonpy-tddft-worker-')
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
        self._proc.stdin.write(text+'\n');self._proc.stdin.flush()

    def _receive(self,index):
        try: line=self._events.get(timeout=self._timeout)
        except queue.Empty as exc: raise TimeoutError('SALMON TDDFT worker response timed out') from exc
        if line is None:
            log=Path(self._work.name)/'run/worker.log'
            raise RuntimeError('SALMON TDDFT worker exited before reply:\n'+log.read_text()[-6000:])
        fields=line.split()
        if len(fields)!=5 or int(fields[1])!=index: raise RuntimeError('TDDFT worker protocol sequence mismatch')
        return _array([float(x) for x in fields[2:]],'worker current',(3,))

    def _stop(self):
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
            self._proc=None
        if self._work is not None: self._work.cleanup();self._work=None

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
        self._current=current;self._history.append(A.copy());self._steps+=1
        return self.current

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
