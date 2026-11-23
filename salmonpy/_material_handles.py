"""LD/SBE library bindings and persistent state ownership.

Importing does not load a library. Each model family retains its own lock/cache.
Numerical updates and state storage remain in the compiled library.
"""
import ctypes as c
from functools import lru_cache
import hashlib
from pathlib import Path
import sys
import threading
import weakref
import numpy as np
from .errors import InputError, NumericalError, RuntimeFailure

H = c.c_int64
_ld_lock = threading.RLock()
_sbe_lock = threading.RLock()


class _Library:
    """One checksum/load/signature/status contract for both calculation families."""
    def __init__(self, path, sha256=None):
        self.path = Path(path).resolve()
        if sha256 is None:
            record = self.path.with_suffix(self.path.suffix + '.sha256')
            if not record.is_file():
                raise RuntimeFailure(f'missing {self.label} library checksum record')
            sha256 = record.read_text().strip()
        if not self.path.is_file() or hashlib.sha256(self.path.read_bytes()).hexdigest() != sha256:
            raise RuntimeFailure(f'{self.label} library checksum mismatch')
        self.sha256 = sha256
        self.library = c.CDLL(str(self.path))
        self.lock = self.family_lock

    def call(self, name, *args):
        with self.lock:
            status = getattr(self.library, name)(*args)
        if status == 4:
            raise NumericalError(f'{self.label} rejected nonfinite arithmetic; state unchanged')
        if status == 1:
            raise InputError(f'invalid {self.label} parameters')
        if status:
            raise RuntimeFailure(f'{self.label} {name} failed with status {status}')

    def bind(self, signatures, result_type):
        try:
            for name, args in signatures.items():
                fn = getattr(self.library, name)
                fn.argtypes = args
                fn.restype = result_type
        except AttributeError as exc:
            raise RuntimeFailure(f'{self.label} handle ABI missing; explicitly rebuild kernels or install matching wheel') from exc


def _library_path():
    suffix = {'darwin': '.dylib', 'win32': '.dll'}.get(sys.platform, '.so')
    return Path(__file__).parent / '_native' / ('libsalmonpy_ld' + suffix)


class LDAPI(_Library):
    label = "native LD"
    family_lock = _ld_lock

    def __init__(self, path, sha256=None):
        super().__init__(path, sha256)
        signatures={
            'sp_ld_create':[c.c_int,c.c_double,c.c_double]+[c.c_void_p]*4+[c.POINTER(H)],
            'sp_ld_advance':[H,c.c_int]+[c.c_void_p]*3,
            'sp_ld_advance_performance':[H,c.c_int]+[c.c_void_p]*3,
            'sp_ld_export':[H]+[c.c_void_p]*2,
            'sp_ld_set':[H]+[c.c_void_p]*5,
            'sp_ld_destroy':[H]}
        self.bind(signatures, c.c_int)
        self.legacy=self.library.salmonpy_ld
        self.legacy.argtypes=[c.c_int,c.c_int,c.c_double,c.c_double]+[c.POINTER(c.c_double)]*8
        self.legacy.restype=None


class SBEAPI(_Library):
    label = "native SBE"
    family_lock = _sbe_lock

    def __init__(self, path, sha256=None):
        super().__init__(path, sha256)
        signatures={
            'sp_sbe_model_create':[c.c_int32]*5+[c.c_double]+[c.c_void_p]*4+[c.POINTER(H)],
            'sp_sbe_model_destroy':[H],
            'sp_sbe_state_create':[H,c.c_double,c.c_void_p,c.POINTER(H),c.c_void_p],
            'sp_sbe_state_create_performance':[H,c.c_double,c.c_void_p,c.POINTER(H),c.c_void_p],
            'sp_sbe_step':[H,c.c_void_p,c.c_void_p],
            'sp_sbe_step_performance':[H,c.c_void_p,c.c_void_p],
            'sp_sbe_density':[H,c.c_void_p],
            'sp_sbe_set':[H,c.c_void_p,c.c_void_p,c.c_void_p],
            'sp_sbe_state_destroy':[H]}
        self.bind(signatures, c.c_int32)
        self.legacy=self.library.salmonpy_sbe
        self.legacy.argtypes=[c.c_int]*5+[c.c_double]*2+[c.c_void_p]*7+[c.c_int,c.c_void_p]
        self.legacy.restype=None


@lru_cache(maxsize=1)
def default_ld_api():
    return LDAPI(_library_path())


@lru_cache(maxsize=1)
def default_sbe_api():
    return SBEAPI(_library_path())


class LDState:
    def __init__(self,api,material,dt,field,*,execution_mode="checked"):
        self.api=api;self.n=len(material._f);self.handle=H(0)
        self._advance_name="sp_ld_advance_performance" if execution_mode=="performance" else "sp_ld_advance"
        if not 1<=self.n<=(2**31-1)//3:raise InputError('unsupported LD pole count')
        from .ld import TensorLorentzDrude
        if isinstance(material,TensorLorentzDrude):
            try:fn=api.library.sp_ld_create_tensor
            except AttributeError as exc:raise RuntimeFailure('native tensor state ABI missing; explicitly rebuild') from exc
            fn.argtypes=[c.c_int,c.c_double,c.c_double]+[c.c_void_p]*5+[c.POINTER(H)];fn.restype=c.c_int
            api.call('sp_ld_create_tensor',self.n,dt,material._wp,
                *(a.ctypes.data for a in (material._strength_matrix,material._w,material._g,material._epsilon_matrix,field)),c.byref(self.handle))
        else:
            for a in (material._f,material._w,material._g):
                if a.shape!=(self.n,) or a.dtype!=np.dtype(float) or not a.flags.c_contiguous:
                    raise InputError('LD material array layout changed after construction')
            api.call('sp_ld_create',self.n,dt,material._wp,
                     *(a.ctypes.data for a in (material._f,material._w,material._g,field)),c.byref(self.handle))
        try:self.finalizer=weakref.finalize(self,api.call,'sp_ld_destroy',int(self.handle.value))
        except BaseException as exc:
            try:api.call('sp_ld_destroy',int(self.handle.value))
            except Exception as cleanup:
                if hasattr(exc,'add_note'):exc.add_note('LD cleanup also failed: '+str(cleanup))
            raise

    @property
    def alive(self):return self.finalizer.alive

    def close(self):self.finalizer()

    def advance(self,field,p,j):
        self.api.call(self._advance_name,self.handle,len(field),field.ctypes.data,p.ctypes.data,j.ctypes.data)

    def export(self):
        p=np.empty((self.n,3));j=np.empty_like(p)
        self.api.call('sp_ld_export',self.handle,p.ctypes.data,j.ctypes.data)
        return p,j

    def replace(self,field,p,j):
        field=np.ascontiguousarray(field,dtype=np.float64)
        p=np.ascontiguousarray(p,dtype=np.float64);j=np.ascontiguousarray(j,dtype=np.float64)
        if field.shape!=(3,) or p.shape!=(self.n,3) or j.shape!=(self.n,3):
            raise InputError('invalid LD history shape')
        pt=np.empty(3);jt=np.empty(3);result=(pt,jt)
        self.api.call('sp_ld_set',self.handle,*(a.ctypes.data for a in (field,p,j,pt,jt)))
        return result


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
    with _sbe_lock:
        model=_models.get(initial)
        if model is None:
            model=SBEModel(initial,api);_models[initial]=model
        return model


class SBEModel:
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

    def new_state(self,dt,A,*,execution_mode="checked"):return SBEState(self,dt,A,execution_mode=execution_mode)


class SBEState:
    def __init__(self,model,dt,A,*,execution_mode="checked"):
        self.model=model;self.handle=H(0);self.current=np.empty(3)
        self._step_name='sp_sbe_step_performance' if execution_mode=='performance' else 'sp_sbe_step'
        create='sp_sbe_state_create_performance' if execution_mode=='performance' else 'sp_sbe_state_create'
        model.api.call(create,model.handle,dt,A.ctypes.data,c.byref(self.handle),self.current.ctypes.data)
        try:self.finalizer=weakref.finalize(self,_release_state,model.api,int(self.handle.value),model)
        except BaseException as exc:
            _cleanup_failure(model.api,'sp_sbe_state_destroy',int(self.handle.value),exc);raise

    @property
    def alive(self):return self.finalizer.alive

    def step(self,A):
        if not self.alive:raise RuntimeFailure('SBE state is closed')
        current=np.empty(3)
        self.model.api.call(self._step_name,self.handle,A.ctypes.data,current.ctypes.data)
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
