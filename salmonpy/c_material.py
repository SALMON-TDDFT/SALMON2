"""Explicitly load a checksum-verified serial C material plugin (preview ABI1).

No native build, MPI initialization or library load on import. Native plugins
are trusted executable code; ABI checks do not isolate faults in native code.
"""
import ctypes as c
from pathlib import Path
import hashlib
import json
import threading
import numpy as np
from .errors import InputError, InputTypeError, RuntimeFailure
from .ld import _array, _positive
from .timing import TimeSpec,InputPort,OutputPort

__all__=['CMaterial','CABIError']
H=c.c_int64;P=c.POINTER(c.c_double);E=c.POINTER(c.c_char)
CreateModel=c.CFUNCTYPE(c.c_int32,c.c_void_p,c.c_size_t,c.POINTER(H),E,c.c_size_t)
Destroy=c.CFUNCTYPE(c.c_int32,H,E,c.c_size_t)
CreateState=c.CFUNCTYPE(c.c_int32,H,c.c_double,c.c_double,P,c.POINTER(H),E,c.c_size_t)
Observe=c.CFUNCTYPE(c.c_int32,H,P,P,E,c.c_size_t)
Step=c.CFUNCTYPE(c.c_int32,H,P,P,P,E,c.c_size_t)


class _API(c.Structure):
    _fields_=[('abi_version',c.c_int32),('struct_size',c.c_int32),('backend',c.c_char_p),
              ('input_lag',c.c_int32),('uses_initial_A',c.c_int32),
              ('model_create',CreateModel),('model_destroy',Destroy),('state_create',CreateState),
              ('state_observe',Observe),('state_step',Step),('state_destroy',Destroy)]


class CABIError(RuntimeFailure):
    """A native callback returned nonzero status; status retains the ABI code."""
    def __init__(self,status,message):
        self.status=int(status)
        super().__init__(f'C material status {self.status}: {message}')


_abi_lock=threading.RLock()


class _Library:
    def __init__(self,path,sha256,loaded=None):
        self.path=path;self.sha256=sha256;self.lock=_abi_lock
        self.library=c.CDLL(str(path)) if loaded is None else loaded
        try:query=self.library.salmonpy_material_query
        except AttributeError as exc:raise InputError('missing salmonpy_material_query') from exc
        query.argtypes=[c.c_int32,c.c_size_t,c.POINTER(_API)];query.restype=c.c_int32
        self.api=_API();status=query(1,c.sizeof(_API),c.byref(self.api))
        if status or self.api.abi_version!=1 or self.api.struct_size!=c.sizeof(_API):
            raise InputError('unsupported C material ABI version or table size')
        if (not self.api.backend or self.api.input_lag<0 or self.api.uses_initial_A not in (0,1)
                or any(not getattr(self.api,name) for name in ('model_create','model_destroy','state_create','state_observe','state_step','state_destroy'))):
            raise InputError('incomplete C material ABI metadata/callbacks')
        try:self.backend=self.api.backend.decode('utf-8')
        except UnicodeError as exc:raise InputError('C material backend must be UTF-8') from exc
        if not self.backend.strip():raise InputError('empty C material backend')

    def call(self,name,*args):
        with self.lock:
            error=c.create_string_buffer(1024)
            status=getattr(self.api,name)(*args,error,len(error))
            if status:raise CABIError(status,error.value.decode('utf-8',errors='replace'))


_libraries={}
_pinned_paths={}
_pinned_files={}
_pinned_images={}


def _load(path,sha256):
    path=Path(path).expanduser().resolve()
    if hashlib.sha256(path.read_bytes()).hexdigest()!=sha256:raise InputError('C material library checksum mismatch')
    info=path.stat();file_id=(info.st_dev,info.st_ino) if info.st_ino else None
    with _abi_lock:
        if (path in _pinned_paths and _pinned_paths[path]!=sha256
                or file_id in _pinned_files and _pinned_files[file_id]!=sha256):
            raise InputError('loaded library is pinned; use a new path/file for a different build')
        if path in _libraries:
            library=_libraries[path]
            if library.sha256!=sha256:raise InputError('loaded library path is pinned; use a new path for a different build')
            return library
        raw=c.CDLL(str(path))
        # dlopen succeeded: pin even if subsequent ABI/query validation fails.
        mapped=_pinned_images.get(raw._handle)
        actual=sha256 if mapped is None else mapped[0]
        _pinned_paths[path]=actual
        if file_id is not None:_pinned_files[file_id]=actual
        if actual!=sha256:raise InputError('OS loader reused another image; use a new library filename')
        _pinned_images[raw._handle]=(actual,raw)
        library=_Library(path,sha256,loaded=raw);_libraries[path]=library
        return library


def _finish(owner,primary):
    try:owner.close()
    except Exception as cleanup:
        if primary is None:raise
        if hasattr(primary,'add_note'):primary.add_note('C material cleanup also failed: '+str(cleanup))
        else:primary.material_cleanup_error=str(cleanup)


class CMaterial:
    """Library-owned read-only model; configuration bytes are plugin-defined.

    Calls must be serialized by a host; Python serializes native calls process-wide;
    Python model/state lifecycle operations are not thread-shareable.
    close prevents new states and defers native model destruction until all
    existing states close. Use context managers for deterministic ownership.
    """
    def __init__(self,library,sha256,configuration):
        if type(configuration) is not bytes:raise InputTypeError('configuration must be immutable bytes')
        self._library=_load(library,sha256);self._handle=H(0);self._states=0;self._closed=False
        self._configuration=configuration
        buffer=c.create_string_buffer(configuration)
        self._library.call('model_create',c.cast(buffer,c.c_void_p),len(configuration),c.byref(self._handle))
        if self._handle.value<=0:raise RuntimeFailure('C material returned an invalid model handle')
        self._backend=self._library.backend
        self._identity='salmonpy.c-material/1/'+hashlib.sha256(json.dumps(
            [sha256,hashlib.sha256(configuration).hexdigest(),self._backend]).encode()).hexdigest()

    @property
    def backend(self):return self._backend
    @property
    def identity(self):return self._identity

    def new_state(self,*,dt,time0=0.,initial_A=(0.,0.,0.)):
        if self._closed:raise RuntimeFailure('C material model is closed')
        return _State(self,dt,time0,initial_A)

    def _release(self):
        if self._closed and not self._states and self._handle.value:
            self._library.call('model_destroy',self._handle)
            self._handle=H(0)

    def close(self):
        self._closed=True;self._release()
    def __enter__(self):return self
    def __exit__(self,kind,error,traceback):_finish(self,error)


class _State:
    def __init__(self,model,dt,time0,A):
        self._model=model;self._dt=_positive(dt,'dt');self._time0=float(_array(time0,'time0',()))
        values=_array(A,'initial_A',(3,));self._steps=0;self._closed=False;self._failed=False
        self._handle=H(0);self._library=model._library
        self._library.call('state_create',model._handle,self._dt,self._time0,values.ctypes.data_as(P),c.byref(self._handle))
        if self._handle.value<=0:raise RuntimeFailure('C material returned an invalid state handle')
        model._states+=1
        try:
            clock=c.c_double();current=np.empty(3)
            self._library.call('state_observe',self._handle,c.byref(clock),current.ctypes.data_as(P))
            self._accept(clock.value,current,0)
        except BaseException as exc:
            _finish(self,exc);raise
        self.timing=TimeSpec(self._dt,(InputPort('A',lag=self._library.api.input_lag,
            uses_t0=bool(self._library.api.uses_initial_A)),),(OutputPort('j'),),t0=self._time0,
            integrator='C material ABI1: '+model.backend)

    @property
    def dt(self):return self._dt
    @property
    def time(self):return self._time
    @property
    def current(self):return self._current.copy()

    def _accept(self,clock,current,index):
        expected=self._time0+index*self._dt
        tolerance=min(abs(self._dt)*1e-3,max(abs(self._dt)*1e-8,
            8*np.finfo(float).eps*(index+1)*max(abs(expected),abs(self._time0),abs(self._dt))))
        if not np.isfinite(clock) or not np.isfinite(current).all() or abs(clock-expected)>tolerance:
            self._failed=True
            raise RuntimeFailure('C material returned an invalid clock/current; state invalidated')
        self._time=float(clock);self._current=current.copy();self._steps=index

    def step(self,A):
        if self._closed or self._failed:raise RuntimeFailure('C material state is closed or invalidated')
        values=_array(A,'A',(3,));clock=c.c_double();current=np.empty(3)
        self._library.call('state_step',self._handle,values.ctypes.data_as(P),c.byref(clock),current.ctypes.data_as(P))
        self._accept(clock.value,current,self._steps+1)
        return self.current

    def close(self):
        if self._closed:return
        self._library.call('state_destroy',self._handle)
        self._closed=True;self._handle=H(0);self._model._states-=1;self._model._release()
    def __enter__(self):return self
    def __exit__(self,kind,error,traceback):_finish(self,error)
