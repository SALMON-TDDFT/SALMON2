"""Private persistent LD ownership; explicit verified library loading only."""
import ctypes as c
from functools import lru_cache
import hashlib
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
        self.path=Path(path).resolve()
        if sha256 is None:
            record=self.path.with_suffix(self.path.suffix+'.sha256')
            if not record.is_file():raise RuntimeFailure('missing native LD checksum record')
            sha256=record.read_text().strip()
        if not self.path.is_file() or hashlib.sha256(self.path.read_bytes()).hexdigest()!=sha256:
            raise RuntimeFailure('native LD library checksum mismatch')
        self.sha256=sha256;self.library=c.CDLL(str(self.path));self.lock=_lock
        signatures={
            'sp_ld_create':[c.c_int,c.c_double,c.c_double]+[c.c_void_p]*4+[c.POINTER(H)],
            'sp_ld_advance':[H,c.c_int]+[c.c_void_p]*3,
            'sp_ld_export':[H]+[c.c_void_p]*2,
            'sp_ld_set':[H]+[c.c_void_p]*5,
            'sp_ld_destroy':[H]}
        try:
            for name,args in signatures.items():
                fn=getattr(self.library,name);fn.argtypes=args;fn.restype=c.c_int
        except AttributeError as exc:
            raise RuntimeFailure('native LD handle ABI missing; explicitly rebuild kernels or install matching wheel') from exc
        self.legacy=self.library.salmonpy_ld
        self.legacy.argtypes=[c.c_int,c.c_int,c.c_double,c.c_double]+[c.POINTER(c.c_double)]*8
        self.legacy.restype=None

    def call(self,name,*args):
        with self.lock:status=getattr(self.library,name)(*args)
        if status==4:raise NumericalError('native LD rejected nonfinite arithmetic; state unchanged')
        if status==1:raise InputError('invalid native LD parameters')
        if status:raise RuntimeFailure(f'native LD {name} failed with status {status}')


@lru_cache(maxsize=1)
def default_api():
    suffix={'darwin':'.dylib','win32':'.dll'}.get(sys.platform,'.so')
    return API(Path(__file__).parent/'_native'/('libsalmonpy_ld'+suffix))


class NativeState:
    def __init__(self,api,material,dt,field):
        self.api=api;self.n=len(material._f);self.handle=H(0)
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
        self.api.call('sp_ld_advance',self.handle,len(field),field.ctypes.data,p.ctypes.data,j.ctypes.data)

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
