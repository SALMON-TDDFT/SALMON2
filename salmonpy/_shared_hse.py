"""Shared-model leases and token-bound compatibility for native fast-PR HSE.

No import-time MPI attachment, native loading or compilation.
"""
import ctypes as c
import hashlib
import json
import numpy as np

from .errors import InputError, InputTypeError, RuntimeFailure
from .fast_hse import _positive
from .sbe import _integer


def _configure_handles(lib, *, symmetry=False):
    vector=np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,shape=(3,),flags='C_CONTIGUOUS')
    specs={'sp_fast_hse_handle_protocol':[], 'sp_fast_hse_model_protocol':[],
           'sp_fast_hse_handle_init':[c.c_char_p,c.c_double,c.c_int,c.c_int,c.POINTER(c.c_int64),vector],
           'sp_fast_hse_model_new':[c.c_double,c.c_int,c.POINTER(c.c_int64),vector],
           'sp_fast_hse_handle_step':[c.c_int64,vector,vector,vector,c.POINTER(c.c_double)],
           'sp_fast_hse_handle_info':[c.c_int64,c.POINTER(c.c_int),vector,vector],
           'sp_fast_hse_handle_history':[c.c_int64,np.ctypeslib.ndpointer(dtype=np.int32,ndim=1,shape=(8,),flags='C_CONTIGUOUS')],
           'sp_fast_hse_handle_transaction':[c.c_int64,c.c_int],
           'sp_fast_hse_handle_checkpoint':[c.c_int64,c.c_char_p,c.c_int],
           'sp_fast_hse_handle_close':[c.c_int64],
           'sp_fast_hse_handle_label':[c.c_int64,c.c_void_p,c.c_int]}
    for name,args in specs.items():
        fn=getattr(lib,name,None)
        if fn is None:raise InputError('shared HSE model ABI missing; explicitly rebuild current runtime')
        fn.argtypes=args;fn.restype=c.c_int
    if lib.sp_fast_hse_handle_protocol()!=1 or lib.sp_fast_hse_model_protocol()!=1:
        raise InputError('unsupported shared HSE model protocol')

    if symmetry:
        query=getattr(lib,'sp_fast_hse_shared_symmetry_protocol',None)
        if query is None:raise InputError('shared symmetry HSE requires an explicit current runtime build')
        query.argtypes=[];query.restype=c.c_int
        if query()!=1:raise InputError('unsupported shared symmetry HSE protocol')


class _Bound:
    def __init__(self,fn,token):self.fn,self.token=fn,token
    def __call__(self,*args):return self.fn(self.token,*args)


class _HandleLibrary:
    def __init__(self,lib,token):
        self._lib=lib
        for old,new in [('step','step'),('checkpoint','checkpoint'),('info','info'),
                        ('close','close'),('history_info','history'),('step_transaction','transaction')]:
            setattr(self,'sp_fast_hse_'+old,_Bound(getattr(lib,'sp_fast_hse_handle_'+new),token))
    def __getattr__(self,name):return getattr(self._lib,name)


class _Lease:
    def __init__(self,model):
        self.model=model;self.context=model._native.context;self.closed=False;self.state=None
    def _live(self):
        if self.closed:raise RuntimeFailure('native HSE state lease is closed')
        self.model._native._live()
    def close(self):
        if self.closed:return
        self.closed=True
        self.model._states.discard(self.state)
        self.model._detach_if_unused()


class NativeFastHSEModel:
    """One immutable material/symmetry model and shared host-MPI lease; states evolve independently.

    Construct with NativeFastHSEBackend.open_model. close prevents new states;
    existing states keep their connection until individually closed. All calls
    are collective/main-thread. No state or file is finalized by Python GC.
    """
    def __init__(self,native,lib,work,dt,capacity,propagator,backend_identity,manifest,*,symmetry=False):
        self._native=native;self.lib=lib;self.work=work
        self.dt=dt;self.capacity=capacity;self.propagator=propagator
        self.symmetry_enabled=symmetry
        self.closed=False;self.failed=False;self._attempted=False;self._initialized=False
        self._states=set();self._serial=0
        self.identity='salmonpy.shared-fast-hse/1/'+hashlib.sha256(
            json.dumps([backend_identity,manifest],sort_keys=True).encode()).hexdigest()

    def _call(self,fn,*args):
        from ._native_group import _directory_call
        self._native._live()
        return _directory_call(self.lib,self._native.context.solver_comm,self.work,fn,*args)

    def new_state(self,point_id=None,*,dt=None,capacity=None):
        """Fresh original-GS state, never a copy of another state's excitation."""
        from .fast_hse_mpi import _State, _agree
        comm=self._native.context.solver_comm
        error=None
        try:
            self._native._live()
            if self.closed or self.failed:raise InputError('shared HSE model is closed or failed')
            timestep=self.dt if dt is None else _positive(dt,'dt')
            count=self.capacity if capacity is None else _integer(capacity,'capacity',1,2**31-2)
            point=self._serial if point_id is None else point_id
            if isinstance(point,(bool,np.bool_)) or not isinstance(point,(str,int,np.integer)):
                raise InputTypeError('point_id must be a string or integer')
            point=int(point) if isinstance(point,np.integer) else point
            descriptor=(point,timestep,count,self._serial,self._initialized)
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        descriptors=comm.allgather(descriptor)
        _agree(comm,'shared HSE state parameters differ within subgroup' if any(x!=descriptor for x in descriptors) else None)
        token=c.c_int64();current=np.zeros(3)
        first=not self._initialized
        self._attempted=True
        if first:
            code=self._call(self.lib.sp_fast_hse_handle_init,b'inputfile',timestep,count,int(self.symmetry_enabled),c.byref(token),current)
        else:
            code=self._call(self.lib.sp_fast_hse_model_new,timestep,count,c.byref(token),current)
        statuses=comm.allgather((int(code),token.value,bool(np.isfinite(current).all())))
        if any(code or handle<=0 or not finite for code,handle,finite in statuses) or len({h for _,h,_ in statuses})!=1:
            if first:self.failed=True
            raise RuntimeFailure('shared native HSE state initialization failed: '+str(statuses))
        self._initialized=True
        bound=_HandleLibrary(self.lib,token.value)
        lease=_Lease(self)
        state=_State(lease,bound,self.work,timestep,count,self.propagator)
        lease.state=state
        state.handle=token.value;state.point_id=point;state.model_identity=self.identity
        state.initialized=True;state.symmetry_enabled=self.symmetry_enabled
        self._states.add(state)
        try:
            label=c.create_string_buffer(81)
            code=self._call(self.lib.sp_fast_hse_handle_label,token.value,label,len(label))
            _agree(comm,'native HSE output label query failed' if code else None,runtime=True)
            prefix=label.value.decode('ascii')
            if not prefix or any(ch not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_' for ch in prefix):
                raise RuntimeFailure('invalid shared native output basename')
            state.output_files=(prefix+'_rt.data',prefix+'_rt_energy.data')
            step=c.c_int();A=np.zeros(3);j=np.zeros(3)
            code=state._call(bound.sp_fast_hse_info,c.byref(step),A,j)
            _agree(comm,'new native HSE endpoint invalid' if code or step.value!=0 or not np.isfinite(A).all() or not np.isfinite(j).all() else None,runtime=True)
            state.A=A.copy();state.current=j.copy();state.initial_current=current.copy()
            self._serial+=1
            return state
        except BaseException:
            state.close()
            raise

    def _detach_if_unused(self):
        if self.closed and not self._states and not self._native.closed:
            self._native.close()
            if not self._attempted:
                from . import fast_hse_mpi
                fast_hse_mpi._used=fast_hse_mpi.embedded._used=False

    def close(self):
        if self.closed:return
        self._native._live()
        self.closed=True
        self._detach_if_unused()

    def __enter__(self):return self
    def __exit__(self,*args):self.close()


from ._native_group import _SessionGroup
