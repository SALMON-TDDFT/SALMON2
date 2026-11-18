"""Borrow a retained-GS session attachment; never write/read a GS restart bundle."""
import ctypes as c
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
from .errors import InputError,RuntimeFailure
from .fast_hse import _positive
from .sbe import _integer
from ._shared_tddft import _agree


class _BorrowedAttachment:
    def __init__(self,native):self.native=native;self.context=native.context;self._closed=False
    @property
    def closed(self):return self._closed or self.native.closed
    def _live(self):
        self.native._live()
        if self._closed:raise InputError('borrowed GS/TDDFT attachment is closed')
    def close(self):self._closed=True


def open_model(gs,dt,max_steps):
    from ._tddft_registry import NativeTDDFTRegistry,_RegisteredModel,_RegisteredLibrary,_RegistryLease,_configure_registry
    from .embedded_tddft import NativeTDDFTInitialState
    gs._live();session=gs.session;comm=session.comm;error=None
    try:
        dt=_positive(dt,'dt');count=_integer(max_steps,'max_steps',1,2**31-2)
        if session.backend.decomposition[1]>gs._electrons//2:
            raise InputError('RT occupied orbital count is smaller than orbital ranks; use spatial/k splitting')
        _configure_registry(session.lib)
        query=getattr(session.lib,'sp_tddft_memory_gs_protocol',None)
        if query is None or query()!=1:raise InputError('in-memory GS→TDDFT protocol missing; explicitly use a matching runtime')
        fn=session.lib.sp_tddft_model_from_gs
        fn.argtypes=[c.c_int64,c.c_double,c.c_int,c.c_void_p,c.POINTER(c.c_int64),c.POINTER(c.c_int64),c.c_void_p]
        fn.restype=c.c_int
        descriptor=(gs._token,dt,count,gs.material_id,session.backend.decomposition)
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    _agree(comm,'GS/TDDFT model choices differ' if any(v!=descriptor for v in comm.allgather(descriptor)) else None)
    work=None;error=None
    if comm.rank==0:
        try:
            work=Path(tempfile.mkdtemp(prefix='gs-rt-',dir=session.backend.work_root))
            for name,data in gs._assets:
                (work/name).parent.mkdir(parents=True,exist_ok=True)
                (work/name).write_bytes(data)
        except Exception as exc:error=str(exc)
    error=comm.bcast(error,root=0);work=comm.bcast(work,root=0)
    _agree(comm,error)
    try:
        for name,data in gs._assets:
            if (work/name).read_bytes()!=data:
                raise InputError('RT pseudopotential is not shared or differs')
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    dec=np.asarray(session.backend.decomposition,dtype=np.int32)
    model_token=c.c_int64();state_token=c.c_int64();current=np.zeros(3)
    code=session._call(work,fn,gs._token,dt,count,dec.ctypes.data,c.byref(model_token),c.byref(state_token),current.ctypes.data)
    outcomes=comm.allgather((int(code),model_token.value,state_token.value,bool(np.isfinite(current).all())))
    if any(code or mid>=0 or sid<=0 or not finite for code,mid,sid,finite in outcomes) or len({(m,s) for _,m,s,_ in outcomes})!=1:
        if all(code==0 and mid<0 and sid>0 for code,mid,sid,_ in outcomes):
            session._call(work,session.lib.sp_tddft_handle_close,state_token.value)
            session._call(work,session.lib.sp_tddft_model_close,model_token.value)
        raise RuntimeFailure('memory GS→TDDFT registration failed: '+str(outcomes))
    registry=session._td_registry
    if registry is None:
        registry=NativeTDDFTRegistry(session.backend,_BorrowedAttachment(session._native),session.lib)
        session._td_registry=registry
    registry._attempted=True
    source=hashlib.sha256(gs.identity.encode()).hexdigest()
    initial=NativeTDDFTInitialState(work,'memory-gs',(('in-memory-gs',source),),count)
    proxy=_RegisteredLibrary(registry,model_token.value,state_token.value,current)
    model=_RegisteredModel(_RegistryLease(registry),proxy,work,initial,dt,tuple(dec.tolist()),
        'memory-gs/'+source,Path('.'))
    model.registry=registry;model.model_token=model_token.value;model.material_id=gs.material_id
    model._initialized=True;model._attempted=True
    registry._models.add(model)
    return model
