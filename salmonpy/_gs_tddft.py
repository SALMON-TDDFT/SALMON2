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
    if getattr(gs.session.backend,'_fast',False) and not getattr(gs.session.backend,'_ordinary_rt',False):
        return _open_fast_model(gs,dt,max_steps)
    from ._tddft_registry import NativeTDDFTRegistry,_RegisteredModel,_RegisteredLibrary,_RegistryLease,_configure_registry
    from .embedded_tddft import NativeTDDFTInitialState
    gs._live();session=gs.session;comm=session.comm;error=None
    try:
        dt=_positive(dt,'dt');count=_integer(max_steps,'max_steps',1,2**31-2)
        if getattr(session.backend,'_ordinary_rt',False):
            from .bands import _input_groups,_scalar
            groups=_input_groups(gs._input_text)
            if _scalar(groups.get('functional',''),'xc','pz').lower() not in ('pz','tbmbj','libxc_pz','libxc_pbe'):
                raise InputError('ordinary TDDFT supports PZ/TB-mBJ/Libxc-PZ/PBE; GS arrays remain available')
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


def _open_fast_model(gs,dt,max_steps):
    """Value-owned GS seed on the existing fast-PR attachment; no GS file copy."""
    from . import fast_hse_mpi as hse
    from .fast_hse import _input_text
    from .bands import _input_groups,_scalar
    from ._hse_registry import NativeFastHSERegistry,_ModelLease,_ModelLibrary,_RegisteredModel,_configure
    from ._shared_hse import _configure_handles
    from ._hse_geometry import preflight
    gs._live();session=gs.session;comm=session.comm;lib=session.lib;error=None
    try:
        dt=_positive(dt,'dt');count=_integer(max_steps,'max_steps',1,2**31-2)
        if session.backend.decomposition!=(comm.size,1,1,1,1):
            raise InputError('fast memory GS→TDDFT requires k-only decomposition')
        query=getattr(lib,'sp_fast_hse_memory_gs_protocol',None)
        if query is None:raise InputError('memory GS→fast HSE requires an explicitly matching runtime')
        query.argtypes=[];query.restype=c.c_int
        if query()!=1:raise InputError('unsupported memory GS→fast HSE protocol')
        groups=_input_groups(gs._input_text)
        groups['calculation']="\n theory='tddft_response'\n"
        groups['control']=("\n sysname='memory_rt', yn_restart='n', base_directory='./',"
                           "\n directory_read_data='./in-memory-no-gs/'\n")
        groups['tgrid']=f'\n dt={dt:.17g},nt={count}\n'
        groups['parallel']=f'\n nproc_k={comm.size},nproc_ob=1,nproc_rgrid=1,1,1\n'
        groups['emfield']="\n ae_shape1='impulse',e_impulse=0,epdir_re1=0,0,1,trans_longi='tr'\n"
        groups['analysis']=groups.get('analysis','')+'\n out_rt_energy_step=1\n'
        groups.setdefault('propagation',"\n propagator='hse_taylor4',yn_predictor_corrector='n',n_hamil=4\n")
        text=''.join('&'+key+body+'/\n' for key,body in groups.items())
        groups,propagator=_input_text(text)
        symmetry='y' in _scalar(groups.get('system',''),'yn_symmetry','n')
        hse._configure(lib);_configure(lib);_configure_handles(lib,symmetry=symmetry)
        fn=lib.sp_fast_hse_model_from_gs
        vector=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
        fn.argtypes=[c.c_int64,c.c_char_p,c.c_double,c.c_int,c.c_int,c.POINTER(c.c_int64),c.POINTER(c.c_int64),vector]
        fn.restype=c.c_int
        boundary=lib.sp_fast_hse_gs_boundary;boundary.argtypes=[c.c_int];boundary.restype=c.c_int
        descriptor=(gs._token,dt,count,gs.material_id,text,gs.identity)
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    _agree(comm,'GS/HSE model choices differ' if any(v!=descriptor for v in comm.allgather(descriptor)) else None)
    work=None;error=None
    if comm.rank==0:
        try:
            work=Path(tempfile.mkdtemp(prefix='gs-hse-',dir=session.backend.work_root))
            (work/'inputfile').write_text(text)
            for name,data in gs._assets:
                (work/name).parent.mkdir(parents=True,exist_ok=True)
                (work/name).write_bytes(data)
        except Exception as exc:error=str(exc)
    error,work=comm.bcast((error,work),root=0);_agree(comm,error)
    try:
        if (work/'inputfile').read_text()!=text:raise InputError('RT inputs differ across ranks')
        for name,data in gs._assets:
            if (work/name).read_bytes()!=data:raise InputError('RT assets differ across ranks')
        operations=hse._symmetry_preflight(work/'sym.dat',_scalar(groups['system'],'yn_symmetry')) if symmetry else None
        preflight(work,groups,operations)
        hse._require_cell_protocol(lib,groups,symmetry,operations)
        if symmetry:hse._require_history_protocol(lib,hse._symmetry_factor(groups,propagator),propagator)
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    mid=c.c_int64();sid=c.c_int64();current=np.zeros(3)
    code=session._call(work,fn,gs._token,b'inputfile',dt,count,int(symmetry),c.byref(mid),c.byref(sid),current)
    outcomes=comm.allgather((int(code),mid.value,sid.value,bool(np.isfinite(current).all())))
    if any(code or m>=0 or h<=0 or not finite for code,m,h,finite in outcomes) or len({(m,h) for _,m,h,_ in outcomes})!=1:
        if all(code==0 and m<0 and h>0 for code,m,h,_ in outcomes):
            session._call(work,lib.sp_fast_hse_handle_close,sid.value)
            session._call(work,lib.sp_fast_hse_model_close,mid.value)
        raise RuntimeFailure('memory GS→fast HSE registration failed: '+str(outcomes))
    registry=session._td_registry
    if registry is None:
        registry=NativeFastHSERegistry(session.backend,_BorrowedAttachment(session._native),lib)
        session._td_registry=registry

    class SeedLibrary(_ModelLibrary):
        def __init__(self):
            super().__init__(registry);self.token=mid.value;self.bootstrap=sid.value
        def sp_fast_hse_handle_init(self,path,timestep,capacity,symmetry,state,j):
            if self.bootstrap and (timestep,capacity)==(dt,count):
                c.cast(state,c.POINTER(c.c_int64))[0]=self.bootstrap;j[:]=current
                self.bootstrap=0;return 0
            if self.bootstrap:
                code=lib.sp_fast_hse_handle_close(self.bootstrap)
                if code:return code
                self.bootstrap=0
            return self.sp_fast_hse_model_new(timestep,capacity,state,j)

    class SeedModel(_RegisteredModel):
        def close(self):
            if self.closed:return
            self._selection()
            if self.lib.bootstrap:
                code=self._call(lib.sp_fast_hse_handle_close,self.lib.bootstrap)
                _agree(comm,'memory GS bootstrap cleanup failed' if code else None,runtime=True)
                self.lib.bootstrap=0
            super().close()

    lease=_ModelLease(registry)
    model=SeedModel(lease,SeedLibrary(),work,dt,count,propagator,session.backend.sha256,
                    ('memory-gs',gs.identity),symmetry=symmetry)
    lease.model=model;model.material_id=gs.material_id;model.model_token=mid.value
    model._registry_id=registry._serial;registry._serial+=1
    registry._models.add(model)
    return model
