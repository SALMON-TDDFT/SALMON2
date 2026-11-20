"""Explicit heterogeneous standard-model registry. No MPI/loading on import."""
import ctypes as c
import math
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import threading
import numpy as np
from .errors import InputError,RuntimeFailure
from .fast_hse import _positive
from .sbe import _integer
from ._shared_tddft import NativeTDDFTModel,_configure,_profile,_agree


def _configure_registry(lib):
    _configure(lib)
    vec=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
    dec=np.ctypeslib.ndpointer(dtype=np.int32,shape=(5,),flags='C_CONTIGUOUS')
    specs={'registry_protocol':[],
           'model_register':[c.c_char_p,c.c_double,c.c_int,vec,dec,c.POINTER(c.c_int64),c.POINTER(c.c_int64),vec],
           'model_spawn':[c.c_int64,c.POINTER(c.c_int64),vec],
           'model_close':[c.c_int64]}
    for name,args in specs.items():
        fn=getattr(lib,'sp_tddft_'+name,None)
        if fn is None:raise InputError('heterogeneous standard TDDFT registry ABI missing; use an explicitly matching runtime')
        fn.argtypes=args;fn.restype=c.c_int
    if lib.sp_tddft_registry_protocol()!=1:raise InputError('unsupported TDDFT registry protocol')


def open_registry(backend,context):
    from . import embedded_tddft as embedded
    from .distributed import NativeParallelContext
    comm=context.solver_comm;error=None
    try:
        if embedded._used:raise InputError('one native electronic attachment lifetime per process')
        if threading.current_thread() is not threading.main_thread():raise InputError('main thread required')
        if embedded._digest(backend.library)!=backend.sha256:raise InputError('native library changed')
        lib=c.CDLL(str(backend.library));_configure_registry(lib)
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    values=comm.allgather((backend.identity,str(backend.work_root)))
    _agree(comm,'registry parameters differ within subgroup' if any(x!=values[0] for x in values) else None)
    owner=context if hasattr(context,'closed') and hasattr(context,'_dependents') else SimpleNamespace(
        solver_comm=comm,closed=False,_dependents=0)
    native=NativeParallelContext(owner,backend.library,backend.sha256)
    embedded._used=True
    return NativeTDDFTRegistry(backend,native,lib)


class _RegistryLease:
    def __init__(self,registry):self.registry=registry;self.context=registry._native.context
    def _live(self):self.registry._native._live()
    def close(self):pass  # The model controls attachment retirement, not each state.


class _RegisteredLibrary:
    def __init__(self,registry,model_token,state_token,current):
        self.registry=registry;self.model_token=model_token
        self.pending=(state_token,current.copy())
    def __getattr__(self,name):return getattr(self.registry.lib,name)
    def sp_tddft_model_new(self,token,current):
        if self.pending is not None:
            handle,initial=self.pending;self.pending=None
            token._obj.value=handle;current[:]=initial
            return 0
        return self.registry.lib.sp_tddft_model_spawn(self.model_token,token,current)


class _RegisteredModel(NativeTDDFTModel):
    def new_state(self,point_id=None):
        if self.registry._native.closed:raise InputError('TDDFT registry is closed')
        self._native._live();comm=self._native.context.solver_comm
        _agree(comm,'TDDFT registry is closed' if self.registry.closed else None)
        tokens=comm.allgather(self.model_token)
        _agree(comm,'model selections differ within subgroup' if any(t!=self.model_token for t in tokens) else None)
        return super().new_state(point_id)
    def _detach_if_unused(self):
        if self.closed and not self._states:
            self.registry._models.discard(self)
            self.registry._detach_if_unused()
    def close(self):
        if self.registry._native.closed:return
        self._native._live();comm=self._native.context.solver_comm
        descriptor=(self.model_token,self.closed,self.lib.pending is not None)
        descriptors=comm.allgather(descriptor)
        _agree(comm,'model close selections differ within subgroup' if any(v!=descriptor for v in descriptors) else None)
        if self.closed:return
        if self.lib.pending is not None:
            handle,_=self.lib.pending
            code=self._call(self.registry.lib.sp_tddft_handle_close,handle)
            _agree(comm,'queued state close failed: '+str(code) if code else None,runtime=True)
            self.lib.pending=None
        code=self._call(self.registry.lib.sp_tddft_model_close,self.model_token)
        _agree(comm,'native model close failed: '+str(code) if code else None,runtime=True)
        self.closed=True;self._detach_if_unused()


class NativeTDDFTRegistry:
    """One host-MPI attachment; independent immutable full-k standard models.

    Each open_model initializes a first state, held until its first new_state call.
    Model/registry close prevents creation and lets existing states finish.
    """
    def __init__(self,backend,native,lib):
        self.backend=backend;self._native=native;self.lib=lib
        self.closed=False;self._attempted=False;self._models=set()
    def _call(self,work,fn,*args):
        from ._native_group import _directory_call
        self._native._live()
        return _directory_call(self.lib,self._native.context.solver_comm,work,fn,*args,minimum_protocol=2)

    def open_model(self,initial,dt,*,decomposition=None):
        from .embedded_tddft import NativeTDDFTInitialState,_manifest,_copy_bundle,_digest
        comm=self._native.context.solver_comm;error=None
        try:
            self._native._live()
            if self.closed:raise InputError('TDDFT registry is closed')
            if not isinstance(initial,NativeTDDFTInitialState):raise InputError('NativeTDDFTInitialState required')
            filename=Path(initial.input_file)
            if filename.is_absolute() or '..' in filename.parts or any(x in initial.input_file for x in ('\0','\n','\r')):
                raise InputError('relative input filename required')
            dt=_positive(dt,'dt')
            count=_integer(initial.max_steps,'max_steps',1,2**31-2)
            dec=tuple(_integer(v,'decomposition',1,2**31-1) for v in
                      (self.backend.decomposition if decomposition is None else decomposition))
            if len(dec)!=5 or math.prod(dec)!=comm.size:raise InputError('decomposition must match registry subgroup size')
            if self.backend.work_root.is_relative_to(initial.directory):raise InputError('work root must be outside the initial bundle')
            if _manifest(initial.directory)!=initial.manifest:raise InputError('initial bundle changed')
            if _digest(self.backend.library)!=self.backend.sha256:raise InputError('native library changed')
            base=_profile(initial)
            descriptor=(initial.manifest,initial.input_file,dt,count,dec)
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        values=comm.allgather(descriptor)
        _agree(comm,'model parameters differ within subgroup' if any(x!=values[0] for x in values) else None)
        work=None;error=None
        if comm.rank==0:
            try:
                work=Path(tempfile.mkdtemp(prefix='model-',dir=self.backend.work_root))
                if _copy_bundle(initial.directory,work/'bundle')!=initial.manifest:raise InputError('initial bundle changed during copy')
            except Exception as exc:
                error=str(exc)
                if work is not None:shutil.rmtree(work)
        error,name=comm.bcast((error,str(work) if work is not None else None),root=0)
        _agree(comm,error)
        work=Path(name)/'bundle';error=None
        try:
            if _manifest(work)!=initial.manifest:raise InputError('work root must be shared within subgroup')
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        model_token=c.c_int64();state_token=c.c_int64();current=np.zeros(3)
        self._attempted=True
        code=self._call(work,self.lib.sp_tddft_model_register,initial.input_file.encode(),dt,count,
                        np.zeros(3),np.asarray(dec,dtype=np.int32),c.byref(model_token),c.byref(state_token),current)
        outcomes=comm.allgather((int(code),model_token.value,state_token.value,bool(np.isfinite(current).all())))
        if any(code or mid>=0 or sid<=0 or not finite for code,mid,sid,finite in outcomes) or len({(m,s) for _,m,s,_ in outcomes})!=1:
            if all(code==0 and mid<0 and sid>0 for code,mid,sid,_ in outcomes) and len({(m,s) for _,m,s,_ in outcomes})==1:
                self._call(work,self.lib.sp_tddft_handle_close,state_token.value)
                self._call(work,self.lib.sp_tddft_model_close,model_token.value)
            raise RuntimeFailure('TDDFT model registration failed: '+str(outcomes))
        proxy=_RegisteredLibrary(self,model_token.value,state_token.value,current)
        identity='registry-v1/'+hashlib.sha256(json.dumps([self.backend.sha256,initial.input_file,dec]).encode()).hexdigest()
        model=_RegisteredModel(_RegistryLease(self),proxy,work,initial,dt,dec,identity,base)
        model.registry=self;model.model_token=model_token.value
        model._initialized=True;model._attempted=True
        self._models.add(model)
        return model
    def _detach_if_unused(self):
        if self.closed and not self._models and not self._native.closed:
            self._native.close()
            if not self._attempted:
                from . import embedded_tddft
                embedded_tddft._used=False
    def close(self):
        if self._native.closed:return
        self._native._live();self.closed=True
        errors=[]
        for model in sorted(self._models,key=lambda value:value.model_token):
            try:model.close()
            except Exception as exc:errors.append(str(exc))
        try:self._detach_if_unused()
        except Exception as exc:errors.append(str(exc))
        if errors:raise RuntimeFailure('TDDFT registry cleanup failed: '+str(errors))
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
