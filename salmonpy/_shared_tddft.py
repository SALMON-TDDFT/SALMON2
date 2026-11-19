"""Same-model standard TDDFT leases. No loading/MPI/build on import."""
import ctypes as c
import hashlib
import json
from pathlib import Path
import numpy as np
from .errors import InputError,InputTypeError,RuntimeFailure
from .sbe import _integer
from .fast_hse import _positive,_relative
from .bands import _input_groups
from ._hse_geometry import _string


def _agree(comm,error,runtime=False):
    errors=comm.allgather(error)
    if any(errors):raise (RuntimeFailure if runtime else InputError)(str(errors))


def _profile(initial):
    groups=_input_groups((initial.directory/initial.input_file).read_text())
    body='\n'.join(text for name,text in groups.items() if name not in ('atomic_coor','atomic_red_coor','reentrance'))
    required={'yn_periodic':('y','n'),'spin':('unpolarized','unpolarized'),
              'yn_symmetry':('n','n'),'yn_md':('n','n'),'yn_spinorbit':('n','n'),
              'yn_jm':('n','n'),'propagator':('middlepoint','middlepoint'),
              'trans_longi':('tr','tr'),'projection_option':('no','no'),'dm_unfold_option':('no','no')}
    for name,(wanted,default) in required.items():
        value=_string(body,name,default).lower()
        if name=='yn_symmetry' and value and set(value)=={'n'}:continue
        if value!=wanted:raise InputError('shared standard TDDFT requires '+name+'='+wanted)
    if _string(body,'xc','pz').lower() not in ('pz','tbmbj','libxc_pz','libxc_pbe'):
        raise InputError('shared standard TDDFT requires PZ/TB-mBJ/Libxc-PZ/PBE')
    for name in ('yn_out_dns_rt','yn_out_dns_ac_je','yn_out_micro_je','yn_out_elf_rt',
                 'yn_out_estatic_rt','yn_out_current_decomposed','yn_out_intraband_current',
                 'yn_out_mag_decomposed_rt','yn_out_mag_micro_rt','yn_out_spin_current_decomposed',
                 'yn_out_spin_current_micro','yn_out_mom_distr_rt','yn_out_rt_energy_components','yn_out_rvf_rt'):
        if _string(body,name,'n').lower()!='n':raise InputError('shared TDDFT stream ownership does not support '+name)
    base=_relative(_string(body,'base_directory','./'))
    _relative(_string(body,'sysname','default'))
    return base


def _configure(lib):
    vec=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
    dec=np.ctypeslib.ndpointer(dtype=np.int32,shape=(5,),flags='C_CONTIGUOUS')
    specs={'handle_protocol':[],'model_protocol':[],
           'handle_init':[c.c_char_p,c.c_double,c.c_int,vec,dec,c.POINTER(c.c_int64),vec],
           'model_new':[c.POINTER(c.c_int64),vec],'handle_step':[c.c_int64,vec,vec,vec],
           'handle_info':[c.c_int64,c.POINTER(c.c_int),vec,vec],
           'handle_checkpoint':[c.c_int64,c.c_char_p,c.c_int],
           'handle_transaction':[c.c_int64,c.c_int],'handle_close':[c.c_int64],
           'handle_label':[c.c_int64,c.c_void_p,c.c_int]}
    for name,args in specs.items():
        fn=getattr(lib,'sp_tddft_'+name,None)
        if fn is None:raise InputError('standard TDDFT handle ABI missing; explicitly use/build a matching runtime')
        fn.argtypes=args;fn.restype=c.c_int
    if lib.sp_tddft_handle_protocol()!=1 or lib.sp_tddft_model_protocol()!=1:
        raise InputError('unsupported standard TDDFT handle/model protocol')


class _Bound:
    def __init__(self,fn,token):self.fn,self.token=fn,token
    def __call__(self,*args):return self.fn(self.token,*args)


class _Library:
    def __init__(self,lib,token):
        self._lib=lib
        for old,new in (('step','step'),('close','close'),('checkpoint','checkpoint'),('step_transaction','transaction')):
            setattr(self,'sp_tddft_'+old,_Bound(getattr(lib,'sp_tddft_handle_'+new),token))
    def __getattr__(self,name):return getattr(self._lib,name)


class _Lease:
    def __init__(self,model):
        self.model=model;self.context=model._native.context;self.closed=False;self.state=None
    def _live(self):
        if self.closed:raise RuntimeFailure('native TDDFT lease is closed')
        self.model._native._live()
    def close(self):
        if self.closed:return
        self.closed=True;self.model._states.discard(self.state);self.model._detach_if_unused()


class NativeTDDFTModel:
    """One immutable full-k model, common dt/capacity/decomposition, independent states."""
    def __init__(self,native,lib,work,initial,dt,decomposition,identity,base):
        self._native=native;self.lib=lib;self.work=work;self._initial=initial
        self._dt=float(dt);self._capacity=initial.max_steps;self._decomposition=decomposition;self.base=base
        self.closed=False;self.failed=False;self._initialized=False;self._attempted=False;self._serial=0;self._states=set()
        self.identity='salmonpy.shared-tddft/1/'+hashlib.sha256(json.dumps(
            [identity,initial.manifest,self.dt,self.capacity],sort_keys=True).encode()).hexdigest()
    @property
    def dt(self):return self._dt
    @property
    def capacity(self):return self._capacity
    @property
    def decomposition(self):return self._decomposition
    @property
    def initial(self):return self._initial
    def _call(self,fn,*args):
        self._native._live()
        import os
        old=Path.cwd()
        try:os.chdir(self.work);return fn(*args)
        finally:os.chdir(old)
    def new_state(self,point_id=None):
        from .embedded_tddft import _NativeState
        from .material_response import _response_timing
        comm=self._native.context.solver_comm;error=None
        try:
            self._native._live()
            if self.closed or self.failed:raise InputError('shared TDDFT model is closed or failed')
            point=self._serial if point_id is None else point_id
            if isinstance(point,(bool,np.bool_)) or not isinstance(point,(str,int,np.integer)):
                raise InputTypeError('point_id must be a string or integer')
            point=int(point) if isinstance(point,np.integer) else point
            descriptor=(point,self._serial,self._initialized)
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        descriptors=comm.allgather(descriptor)
        _agree(comm,'TDDFT state parameters differ within subgroup' if any(x!=descriptor for x in descriptors) else None)
        token=c.c_int64();j=np.zeros(3);self._attempted=True
        if not self._initialized:
            code=self._call(self.lib.sp_tddft_handle_init,self.initial.input_file.encode(),self.dt,self.capacity,
                            np.zeros(3),np.array(self.decomposition,dtype=np.int32),c.byref(token),j)
        else:code=self._call(self.lib.sp_tddft_model_new,c.byref(token),j)
        outcomes=comm.allgather((int(code),token.value,bool(np.isfinite(j).all())))
        if any(code or handle<=0 or not finite for code,handle,finite in outcomes) or len({h for _,h,_ in outcomes})!=1:
            if all(code==0 and h>0 for code,h,_ in outcomes) and len({h for _,h,_ in outcomes})==1:
                self._call(self.lib.sp_tddft_handle_close,token.value)
            # Status1 from model_new is the pre-allocation output-name guard.
            # Its failed reservation advances generation without touching siblings.
            recoverable=(self._initialized and all(code==1 and h==0 for code,h,_ in outcomes))
            if not recoverable:self.failed=True
            raise RuntimeFailure('shared native TDDFT creation failed: '+str(outcomes))
        self._initialized=True
        lease=_Lease(self)
        state=_SharedState(lease,_Library(self.lib,token.value),self.work,self.capacity)
        lease.state=state;state.handle=token.value;state.point_id=point;state.model_identity=self.identity
        state._current_sign=1;state.initialized=True;state.checkpoint_decomposition=self.decomposition
        state.timing=_response_timing(self.dt,'shared native standard TDDFT; charge current j')
        self._states.add(state)
        try:
            label=c.create_string_buffer(257)
            code=self._call(self.lib.sp_tddft_handle_label,token.value,label,len(label))
            _agree(comm,'native TDDFT label failed' if code else None,runtime=True)
            prefix=label.value.decode()
            files=tuple(str(self.base/(prefix+suffix)) for suffix in ('_rt.data','_rt_energy.data'))
            for name in files:_relative(name)
            state.output_files=files
            state._refresh();state.initial_current=j.copy();self._serial+=1
            return state
        except BaseException:state.close();raise
    def _detach_if_unused(self):
        if self.closed and not self._states and not self._native.closed:
            self._native.close()
            if not self._attempted:
                from . import embedded_tddft
                embedded_tddft._used=False
    def close(self):
        if self.closed:return
        self._native._live();self.closed=True;self._detach_if_unused()
    def __enter__(self):return self
    def __exit__(self,*args):self.close()


# Import after model definitions: embedded_tddft exposes NativeTDDFTModel lazily.
from .embedded_tddft import _NativeState
class _SharedState(_NativeState):
    @property
    def electron_number(self):
        """Collective read-only norm of this specific state handle."""
        from .errors import NumericalError
        if self.closed or self.failed:raise InputError('closed or failed native state')
        fn=getattr(self.native.model.lib,'sp_tddft_handle_electron_number',None)
        if fn is None:raise RuntimeFailure('runtime lacks measured electron-number ABI; explicitly rebuild/use a matching runtime')
        fn.argtypes=[c.c_int64,c.POINTER(c.c_double)];fn.restype=c.c_int
        value=c.c_double();code=self._call(fn,self.handle,c.byref(value))
        codes=self.native.context.solver_comm.allgather(int(code))
        if any(v==4 for v in codes):raise NumericalError('nonfinite/nonpositive measured electron number')
        if any(codes):raise RuntimeFailure('electron-number query failed: '+str(codes))
        return value.value

    def _refresh(self):
        step=c.c_int();A=np.zeros(3);j=np.zeros(3)
        model=self.native.model
        code=self._call(model.lib.sp_tddft_handle_info,self.handle,c.byref(step),A,j)
        bad=code or not 0<=step.value<=self.capacity or not np.isfinite(A).all() or not np.isfinite(j).all()
        _agree(self.native.context.solver_comm,'native endpoint query failed' if bad else None,runtime=True)
        endpoint=(step.value,tuple(A),tuple(j))
        endpoints=self.native.context.solver_comm.allgather(endpoint)
        if any(value!=endpoint for value in endpoints):
            self.failed=True
            raise RuntimeFailure('native restored endpoints differ within subgroup')
        self.steps=step.value;self.A=A.copy();self.current=j.copy()
    @property
    def dt(self):return self.native.model.dt
    def advance(self,old,new):
        _agree(self.native.context.solver_comm,'real A vectors required' if np.iscomplexobj(old) or np.iscomplexobj(new) else None)
        result=super().advance(old,new);self.current=result.copy();return result
    def step(self,A):return self.advance(self.A.copy(),A)
    @property
    def time(self):return self.steps*self.dt
    def checkpoint(self,path,restore=False):
        import os
        comm=self.native.context.solver_comm;error=None
        try:
            self.native._live();p=Path(path).expanduser().resolve()
            if self.closed or self.failed:raise InputError('closed or failed state')
            if type(restore) is not bool:raise InputTypeError('restore must be boolean')
            if restore and (self.steps!=0 or not p.is_file()):raise InputError('restore requires a fresh state and existing file')
            if not restore and (p.exists() or not p.parent.is_dir()):raise InputError('new checkpoint path required')
            if any(x in str(p) for x in ('\0','\n','\r')) or len(os.fsencode(p))>4000:raise InputError('invalid checkpoint path')
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        modes=comm.allgather((restore,self.steps));_agree(comm,'checkpoint epoch differs' if any(x!=modes[0] for x in modes) else None)
        fn=self.lib.sp_tddft_checkpoint;fn.argtypes=[c.c_char_p,c.c_int];fn.restype=c.c_int
        code=self._call(fn,os.fsencode(p),int(restore))
        _agree(comm,'native checkpoint rejected: '+str(code) if code else None,runtime=True)
        if restore:self._refresh()
