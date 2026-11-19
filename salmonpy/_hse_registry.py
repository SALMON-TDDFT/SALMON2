"""Explicit different-input fast-PR HSE models on one host-MPI attachment."""
import ctypes as c
import os
from pathlib import Path
import threading
from types import SimpleNamespace
import numpy as np
from .errors import InputError,RuntimeFailure
from ._shared_hse import NativeFastHSEModel


def _configure(lib):
    vector=np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,shape=(3,),flags='C_CONTIGUOUS')
    specs={'sp_fast_hse_registry_protocol':[],
           'sp_fast_hse_model_register':[c.c_char_p,c.c_double,c.c_int,c.c_int,c.POINTER(c.c_int64),c.POINTER(c.c_int64),vector],
           'sp_fast_hse_model_spawn':[c.c_int64,c.c_double,c.c_int,c.POINTER(c.c_int64),vector],
           'sp_fast_hse_model_close':[c.c_int64]}
    for name,args in specs.items():
        fn=getattr(lib,name,None)
        if fn is None:raise InputError('heterogeneous HSE registry requires an explicitly built matching runtime')
        fn.argtypes=args;fn.restype=c.c_int
    if lib.sp_fast_hse_registry_protocol()!=1:raise InputError('unsupported HSE registry protocol')


def open_registry(backend,context):
    from . import fast_hse_mpi as hse
    from .distributed import NativeParallelContext
    comm=context.solver_comm;error=None
    try:
        if hse._used or hse.embedded._used:raise InputError('one native electronic attachment lifetime per process')
        if threading.current_thread() is not threading.main_thread():raise InputError('main thread required')
        if comm.size not in backend.supported_group_sizes:raise InputError('supported HSE subgroup size required')
        if hse._digest(backend.library)!=backend.sha256:raise InputError('native HSE library changed')
        lib=c.CDLL(str(backend.library));hse._configure(lib);_configure(lib)
        from ._shared_hse import _configure_handles
        _configure_handles(lib)
    except Exception as exc:error=str(exc)
    hse._agree(comm,error)
    descriptor=(backend.sha256,str(backend.work_root),backend.allow_symmetry)
    hse._agree(comm,'HSE registry settings differ' if any(v!=descriptor for v in comm.allgather(descriptor)) else None)
    owner=context if hasattr(context,'closed') and hasattr(context,'_dependents') else SimpleNamespace(
        solver_comm=comm,closed=False,_dependents=0)
    native=NativeParallelContext(owner,backend.library,backend.sha256)
    hse._used=hse.embedded._used=True
    return NativeFastHSERegistry(backend,native,lib)


class NativeFastHSERegistry:
    def __init__(self,backend,native,lib):
        self.backend=backend;self._native=native;self.lib=lib
        self.closed=False;self._models=set();self._serial=0
    def _live(self):self._native._live()
    def open_model(self,initial,dt):
        from .fast_hse_mpi import _agree
        self._live()
        _agree(self._native.context.solver_comm,'HSE registry is closed' if self.closed else None)
        model=self.backend._create(initial,None,dt,self._native.context,shared=True,registry=self)
        model._registry_id=self._serial
        self._serial+=1
        return model
    def _detach_if_unused(self):
        if self.closed and not self._models and not self._native.closed:self._native.close()
    def close(self):
        if self.closed:return
        self._live();errors=[]
        for model in sorted(tuple(self._models),key=lambda m:m._registry_id):
            try:
                if not model.closed:model.close()
            except Exception as exc:errors.append(str(exc))
        if errors:raise RuntimeFailure('HSE registry model cleanup failed: '+str(errors))
        self.closed=True;self._detach_if_unused()
    def __enter__(self):return self
    def __exit__(self,*args):self.close()


class _ModelLease:
    def __init__(self,registry):
        self.registry=registry;self.context=registry._native.context;self.closed=False;self.model=None
    def _live(self):
        self.registry._live()
        if self.closed:raise InputError('HSE model attachment lease is closed')
    def close(self):
        if self.closed:return
        self.closed=True
        self.registry._models.discard(self.model)
        self.registry._detach_if_unused()


class _ModelLibrary:
    def __init__(self,registry):self.registry=registry;self.token=0
    def __getattr__(self,name):return getattr(self.registry.lib,name)
    def sp_fast_hse_handle_init(self,path,dt,capacity,symmetry,state,current):
        token=c.c_int64()
        code=self.registry.lib.sp_fast_hse_model_register(path,dt,capacity,symmetry,c.byref(token),state,current)
        if code==0:self.token=token.value
        return code
    def sp_fast_hse_model_new(self,dt,capacity,state,current):
        return self.registry.lib.sp_fast_hse_model_spawn(self.token,dt,capacity,state,current)


class _RegisteredModel(NativeFastHSEModel):
    def _detach_if_unused(self):
        # An idle child never releases the registry's process-wide reservation.
        if self.closed and not self._states and not self._native.closed:self._native.close()
    def _selection(self):
        from .fast_hse_mpi import _agree
        comm=self._native.context.solver_comm
        descriptor=(self._registry_id,self.lib.token,self.closed)
        _agree(comm,'HSE model selections differ' if any(v!=descriptor for v in comm.allgather(descriptor)) else None)
    def _call(self,fn,*args):
        from .fast_hse_mpi import _agree
        self._native._live();previous=Path.cwd();error=None
        try:
            try:os.chdir(self.work)
            except Exception as exc:error=str(exc)
            _agree(self._native.context.solver_comm,error)
            return fn(*args)
        finally:os.chdir(previous)
    def new_state(self,*args,**kwargs):
        self._selection()
        return super().new_state(*args,**kwargs)
    def close(self):
        self._selection()
        if self.closed:return
        self._native._live()
        if self.lib.token:
            from .fast_hse_mpi import _agree
            code=self._call(self.lib.sp_fast_hse_model_close,self.lib.token)
            _agree(self._native.context.solver_comm,'HSE native model close failed' if code else None,runtime=True)
        super().close()
