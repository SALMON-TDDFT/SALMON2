"""Experimental one-point TDDFT on a supplied host-MPI communicator.

One native electronic lifetime per Python process, on its main thread. SALMON
module configuration is not reentrant: start a fresh MPI job for another run.
A shared writable work root is required across members of each material group.
"""
from .errors import InputError, InputTypeError, RuntimeFailure, NumericalError
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import os
import threading
import numpy as np
from .tddft import _digest, _copy_bundle
from .sbe import _integer

_used = False


def _manifest(root):
    records=[]
    for path in sorted(root.rglob('*')):
        if path.is_symlink(): raise InputError('bundle symlinks are unsupported')
        if path.is_file(): records.append((str(path.relative_to(root)),_digest(path)))
        elif not path.is_dir(): raise InputError('non-regular bundle asset')
    return tuple(records)


@dataclass(frozen=True)
class NativeTDDFTInitialState:
    directory: Path
    input_file: str
    manifest: tuple
    max_steps: int


class NativeTDDFTMaterialProvider:
    """Validate a saved GS/RT bundle; each solver group copies its own state.

    The caller keeps the source bundle unchanged until state creation. This does
    not run DFT or silently convert restart formats. All members need file access.
    """
    def __init__(self, directory, input_file, *, max_steps=1000):
        root=Path(directory).expanduser().resolve();inp=Path(input_file)
        if inp.is_absolute() or '..' in inp.parts: raise InputError('relative input filename required')
        for name in ('info.bin','occupation.bin','wfn.bin'):
            if not (root/'data_for_restart'/name).is_file(): raise InputError('missing GS record: '+name)
        text=(root/inp).read_text()
        if '!SALMONPY_END_INPUT' in text or any(len(line)>256 for line in text.splitlines()):
            raise InputError('unsupported input marker or line length')
        steps=_integer(max_steps,'max_steps',1,2**31-2)
        self.initial=NativeTDDFTInitialState(root,str(inp),_manifest(root),steps)
        self.identity='salmonpy.embedded-initial/1/'+hashlib.sha256(json.dumps(
            [str(inp),self.initial.manifest,steps],sort_keys=True).encode()).hexdigest()

    def prepare(self, context):
        if _manifest(self.initial.directory)!=self.initial.manifest:
            raise InputError('initial bundle changed after preparation')
        return self.initial


def _admission_rejected(lib, outcomes, allowed=(10, 11, 12, 13, 14, 15), query_name='sp_tddft_error_protocol'):
    """Only a declared protocol and unanimous status establish no mutation."""
    codes = [status for status, finite in outcomes]
    if not codes or codes[0] not in allowed or len(set(codes)) != 1 or not all(finite for _, finite in outcomes):
        return False
    query = getattr(lib, query_name, None)
    if query is None:
        return False
    import ctypes as c
    query.argtypes = []
    query.restype = c.c_int
    return query() == 1


class NativeTDDFTBackend:
    """Collective native A/c-to-charge-current adapter, one point per group.

    decomposition is (k, orbital, x, y, z); its product must equal group size.
    Declared deep numerical errors restore accepted state before returning.
    Unknown native errors may terminate the process; retain an MPI timeout policy.
    """
    execution='collective'
    reduction='leader'

    def __init__(self, build_record, work_root, *, decomposition, shared_model=False):
        if type(shared_model) is not bool:raise InputTypeError("shared_model must be boolean")
        self.shared_model=shared_model
        record=json.loads(Path(build_record).read_text())
        if record.get('schema')!='salmonpy.embedded-tddft/1': raise InputError('embedded TDDFT build record required')
        self.library=Path(record['library']).resolve();self.sha256=record['sha256']
        if _digest(self.library)!=self.sha256: raise InputError('native library checksum mismatch')
        self.work_root=Path(work_root).expanduser().resolve()
        if not self.work_root.is_dir(): raise InputError('existing shared work root required')
        self.decomposition=tuple(_integer(v,'decomposition',1,2**31-1) for v in decomposition)
        if len(self.decomposition)!=5: raise InputError('decomposition requires k, orbital, x, y, z')
        identity=[self.sha256,self.decomposition,str(self.work_root)]
        if shared_model:identity.append('shared-model-v1')
        self.identity='salmonpy.embedded-tddft/1/'+hashlib.sha256(json.dumps(identity).encode()).hexdigest()

    def open_group(self,initial,dt,context):
        if not self.shared_model:return None
        from ._native_group import _SessionGroup
        return _SessionGroup(self.open_model(initial,dt,context))

    def open_model(self,initial,dt,context):
        return self._create(initial,None,dt,context,shared=True)

    def create_state(self,initial,point_id,dt,context):
        if self.shared_model:raise InputError('shared_model uses ResponseSession/open_group or open_model')
        return self._create(initial,point_id,dt,context,shared=False)

    def _create(self, initial, point_id, dt, context, *, shared):
        global _used
        from types import SimpleNamespace
        from .distributed import NativeParallelContext
        import ctypes as c
        import math, shutil, tempfile
        comm=context.solver_comm
        error=None
        try:
            if _used: raise InputError('one native electronic lifetime per process; start a fresh MPI job')
            if threading.current_thread() is not threading.main_thread(): raise InputError('main thread required')
            if not shared and len(context.point_ids)!=1: raise InputError('native TDDFT requires one point per group')
            if math.prod(self.decomposition)!=comm.size: raise InputError('inner decomposition must match group size')
            if not isinstance(initial,NativeTDDFTInitialState): raise InputTypeError('NativeTDDFTInitialState required')
            _integer(initial.max_steps,'max_steps',1,2**31-2)
            inp=Path(initial.input_file)
            if inp.is_absolute() or '..' in inp.parts or '\x00' in initial.input_file:
                raise InputError('relative input filename required')
            if not (initial.directory/inp).is_file():raise InputError('missing native input')
            if not np.isfinite(dt) or dt<=0: raise InputError('positive finite dt required')
            if _manifest(initial.directory)!=initial.manifest: raise InputError('initial bundle changed')
            if _digest(self.library)!=self.sha256: raise InputError('native library changed')
            if shared:
                from ._shared_tddft import _profile,_configure
                from .fast_hse import _positive
                if self.work_root.is_relative_to(initial.directory):raise InputError('work root must be outside the initial bundle')
                dt=_positive(dt,'dt');base=_profile(initial)
            lib=c.CDLL(str(self.library))
            if shared:_configure(lib)
            vec=np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,shape=(3,),flags='C_CONTIGUOUS')
            dec=np.ctypeslib.ndpointer(dtype=np.int32,ndim=1,shape=(5,),flags='C_CONTIGUOUS')
            lib.sp_tddft_init.argtypes=[c.c_char_p,c.c_double,c.c_int,vec,dec,vec];lib.sp_tddft_init.restype=c.c_int
            lib.sp_tddft_step.argtypes=[vec,vec,vec];lib.sp_tddft_step.restype=c.c_int
            lib.sp_tddft_close.argtypes=[];lib.sp_tddft_close.restype=c.c_int
        except Exception as exc: error=str(exc)
        errors=comm.allgather(error)
        if any(errors): raise InputError(str(errors))
        descriptor=(self.identity,initial.input_file,initial.manifest,dt,initial.max_steps,self.decomposition,shared)
        descriptors=comm.allgather(descriptor)
        if any(x!=descriptor for x in descriptors):raise InputError('TDDFT configuration differs within subgroup')
        work=None;error=None
        if comm.rank==0:
            try:
                work=Path(tempfile.mkdtemp(prefix='point-'+str(point_id)+'-',dir=self.work_root))
                if _copy_bundle(initial.directory,work/'bundle')!=initial.manifest:
                    raise InputError('initial bundle changed during copy')
            except Exception as exc:
                error=str(exc)
                if work is not None: shutil.rmtree(work)
        error,work=comm.bcast((error,str(work) if work else None),root=0)
        if error: raise InputError(error)
        work=Path(work)/'bundle'
        error=None
        try:
            if _manifest(work)!=initial.manifest: raise InputError('work root must be shared by all group ranks')
        except Exception as exc:error=str(exc)
        errors=comm.allgather(error)
        if any(errors):
            if comm.rank==0:shutil.rmtree(work.parent)
            raise InputError(str(errors))
        owner=(context if shared and hasattr(context,'closed') and hasattr(context,'_dependents')
               else SimpleNamespace(solver_comm=comm,closed=False,_dependents=0))
        native=NativeParallelContext(owner,self.library,self.sha256)
        _used=True
        if shared:
            from ._shared_tddft import NativeTDDFTModel
            return NativeTDDFTModel(native,lib,work,initial,dt,self.decomposition,self.identity,base)
        state=_NativeState(native,lib,work,initial.max_steps)
        state.checkpoint_decomposition=self.decomposition
        retryable_initialization=False
        try:
            J=np.empty(3)
            code=state._call(lib.sp_tddft_init,os.fsencode(initial.input_file),float(dt),initial.max_steps,
                             np.zeros(3),np.array(self.decomposition,dtype=np.int32),J)
            state.initialized=(code==0)
            from .material_response import _response_timing
            state.timing=_response_timing(float(dt),'embedded SALMON middlepoint TDDFT; charge current -J_native')
            outcomes=comm.allgather((int(code),bool(np.isfinite(J).all())))
            if any(code or not finite for code,finite in outcomes):
                retryable_initialization=_admission_rejected(lib,outcomes,allowed=(10,11))
                if retryable_initialization:
                    raise InputError('native initialization rejected before setup; retry permitted: '+str(outcomes))
                raise RuntimeFailure('native initialization failed: '+str(outcomes))
            state.initial_current=-J
            return state
        except BaseException:
            state.close()
            if retryable_initialization:
                _used=False
            raise


class _NativeState:
    def __init__(self,native,lib,work,capacity):
        self._current_sign=-1
        self.timing=None  # declared once dt is known (create_state)
        self.lib=lib
        self.native=native;self.work=work;self.capacity=capacity
        self.closed=False;self.failed=False;self.initialized=False;self.steps=0;self.A=np.zeros(3)

    def _call(self,fn,*args):
        self.native._live()
        # Native relative I/O is process-global; calls are main-thread only.
        # Token shims select shared owners. Restore the caller's directory.
        previous=Path.cwd()
        try:
            os.chdir(self.work)
            return fn(*args)
        finally:os.chdir(previous)

    def advance(self,old,new):
        error=None
        try:
            old=np.array(old,dtype=np.float64,order='C');new=np.array(new,dtype=np.float64,order='C')
            if self.closed or self.failed or self.steps>=self.capacity:raise InputError('closed, failed or exhausted native state')
            if old.shape!=(3,) or new.shape!=(3,) or not np.isfinite(new).all():raise InputError('finite A vectors required')
            if not np.array_equal(old,self.A):raise InputError('stale A endpoint')
        except Exception as exc:error=str(exc)
        comm=self.native.context.solver_comm
        errors=comm.allgather(error)
        if any(errors):raise InputError(str(errors))
        values=comm.allgather((old.tolist(),new.tolist()))
        if any(v!=values[0] for v in values):raise InputError('A endpoints differ within group')
        J=np.empty(3)
        code=self._call(self.lib.sp_tddft_step,old,new,J)
        outcomes=comm.allgather((int(code),bool(np.isfinite(J).all())))
        if any(code or not finite for code,finite in outcomes):
            if _admission_rejected(self.lib,outcomes):
                raise InputError('native step rejected before update; accepted state unchanged: '+str(outcomes))
            if _admission_rejected(self.lib,outcomes,allowed=(20,21,22,23,24,25),
                                   query_name='sp_tddft_recovery_protocol'):
                raise NumericalError('native TDDFT step failed; accepted state restored without replay; retry permitted: '+str(outcomes))
            self.failed=True
            raise RuntimeFailure('native step failed; numerical state is not guaranteed unchanged: '+str(outcomes))
        self.A=new.copy();self.steps+=1
        return self._current_sign*J  # Legacy matter J or handle charge j, converted once.

    def close(self):
        if self.closed:return
        code=0
        if self.initialized:
            code=self._call(self.lib.sp_tddft_close)
        codes=self.native.context.solver_comm.allgather(int(code))
        if any(codes):raise RuntimeFailure('native close failed: '+str(codes))
        self.native.close();self.closed=True
        # Keep native files for diagnosis/postprocessing; caller owns work_root.


def __getattr__(name):
    if name=='NativeTDDFTModel':
        from ._shared_tddft import NativeTDDFTModel
        return NativeTDDFTModel
    raise AttributeError(name)
