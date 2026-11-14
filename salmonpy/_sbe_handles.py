"""Owned typed SBE handles; no library load/build on import.

Only small vectors cross the step call. Density export/import is explicit.
"""
import ctypes as c
import hashlib
from functools import lru_cache
from pathlib import Path
import sys
import threading
import weakref
import numpy as np
from .errors import InputError,NumericalError,RuntimeFailure

H=c.c_int64
_lock=threading.RLock()


class API:
    def __init__(self,path,sha256=None):
        path=Path(path).resolve()
        if sha256 is None:
            record=path.with_suffix(path.suffix+'.sha256')
            if not record.is_file():raise RuntimeFailure('missing native SBE checksum record')
            sha256=record.read_text().strip()
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=sha256:
            raise RuntimeFailure('native SBE library checksum mismatch')
        self.library=c.CDLL(str(path));self.lock=_lock;self.sha256=sha256;self.path=path
        signatures={
            'sp_sbe_model_create':[c.c_int32]*5+[c.c_double]+[c.c_void_p]*4+[c.POINTER(H)],
            'sp_sbe_model_destroy':[H],
            'sp_sbe_state_create':[H,c.c_double,c.c_void_p,c.POINTER(H),c.c_void_p],
            'sp_sbe_step':[H,c.c_void_p,c.c_void_p],
            'sp_sbe_density':[H,c.c_void_p],
            'sp_sbe_set':[H,c.c_void_p,c.c_void_p,c.c_void_p],
            'sp_sbe_state_destroy':[H]}
        try:
            for name,args in signatures.items():
                fn=getattr(self.library,name);fn.argtypes=args;fn.restype=c.c_int32
        except AttributeError as exc:
            raise RuntimeFailure('native SBE handle ABI missing; explicitly rebuild kernels or install a matching wheel') from exc
        self.legacy=self.library.salmonpy_sbe
        self.legacy.argtypes=[c.c_int]*5+[c.c_double]*2+[c.c_void_p]*7+[c.c_int,c.c_void_p]
        self.legacy.restype=None

    def call(self,name,*args):
        with self.lock:
            status=getattr(self.library,name)(*args)
        if status==4:raise NumericalError('native SBE rejected nonfinite arithmetic; state unchanged')
        if status==1:raise InputError('invalid native SBE parameters')
        if status:raise RuntimeFailure(f'native SBE {name} failed with status {status}')


@lru_cache(maxsize=1)
def default_api():
    suffix={'darwin':'.dylib','win32':'.dll'}.get(sys.platform,'.so')
    return API(Path(__file__).parent/'_native'/('libsalmonpy_ld'+suffix))


def _cleanup_failure(api,name,handle,primary):
    try:api.call(name,handle)
    except Exception as cleanup:
        if hasattr(primary,'add_note'):primary.add_note('SBE cleanup also failed: '+str(cleanup))
        else:primary.sbe_cleanup_error=str(cleanup)


def _release_state(api,handle,model):
    # Capturing model keeps it alive until the native state has been released.
    api.call('sp_sbe_state_destroy',handle)


_models=weakref.WeakKeyDictionary()


def model_for(initial,api):
    with _lock:
        model=_models.get(initial)
        if model is None:
            model=NativeModel(initial,api);_models[initial]=model
        return model


class NativeModel:
    def __init__(self,initial,api):
        self.api=api;self.shape=initial.eigen.shape;self.handle=H(0)
        if len(self.shape)!=2 or min(self.shape)<1 or max(self.shape)>2**31-1:
            raise InputError('unsupported SBE model dimensions')
        nb,nk=self.shape
        for a,shape,dtype in ((initial.eigen,(nb,nk),np.dtype(float)),
                (initial.weights,(nk,),np.dtype(float)),
                (initial.momentum,(nb,nb,3,nk),np.dtype(complex)),
                (initial.nonlocal_velocity,(nb,nb,3,nk),np.dtype(complex))):
            if a.shape!=shape or a.dtype!=dtype or not a.flags.f_contiguous:
                raise InputError('SBE initial array layout changed after construction')
        api.call('sp_sbe_model_create',nb,nk,initial.electrons,initial.correction_order,
                 int(initial.nonlocal_correction),initial.volume,
                 *(a.ctypes.data for a in (initial.eigen,initial.weights,initial.momentum,initial.nonlocal_velocity)),
                 c.byref(self.handle))
        try:self.finalizer=weakref.finalize(self,api.call,'sp_sbe_model_destroy',int(self.handle.value))
        except BaseException as exc:
            _cleanup_failure(api,'sp_sbe_model_destroy',int(self.handle.value),exc);raise

    def new_state(self,dt,A):return NativeState(self,dt,A)


class NativeState:
    def __init__(self,model,dt,A):
        self.model=model;self.handle=H(0);self.current=np.empty(3)
        model.api.call('sp_sbe_state_create',model.handle,dt,A.ctypes.data,c.byref(self.handle),self.current.ctypes.data)
        try:self.finalizer=weakref.finalize(self,_release_state,model.api,int(self.handle.value),model)
        except BaseException as exc:
            _cleanup_failure(model.api,'sp_sbe_state_destroy',int(self.handle.value),exc);raise

    @property
    def alive(self):return self.finalizer.alive

    def step(self,A):
        if not self.alive:raise RuntimeFailure('SBE state is closed')
        current=np.empty(3)
        self.model.api.call('sp_sbe_step',self.handle,A.ctypes.data,current.ctypes.data)
        self.current=current
        return current

    def density(self):
        if not self.alive:raise RuntimeFailure('SBE state is closed')
        nb,nk=self.model.shape;rho=np.empty((nb,nb,nk),dtype=complex,order='F')
        self.model.api.call('sp_sbe_density',self.handle,rho.ctypes.data)
        return rho

    def replace(self,rho,A,J):
        if not self.alive:raise RuntimeFailure('SBE state is closed')
        nb,nk=self.model.shape
        rho=np.asfortranarray(rho,dtype=complex)
        A=np.ascontiguousarray(A,dtype=float);J=np.ascontiguousarray(J,dtype=float)
        if rho.shape!=(nb,nb,nk) or A.shape!=(3,) or J.shape!=(3,):
            raise InputError('SBE restore array layout mismatch')
        accepted=J.copy()
        self.model.api.call('sp_sbe_set',self.handle,rho.ctypes.data,A.ctypes.data,J.ctypes.data)
        self.current=accepted

    def close(self):
        if self.finalizer.alive:self.finalizer()
