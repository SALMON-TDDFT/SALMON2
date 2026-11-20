"""Experimental fast-PR HSE on an explicitly supplied k-only MPI subgroup.

Build explicitly against the host MPI. One native lifetime on the main thread
per process on the legacy route; open_model provides shared same-model full/reduced-k state leases.
Input and shared work directories are caller-owned. No MPI import,
initialization, finalization, compilation or download occurs on module import.
"""
import ctypes as c
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import threading
from types import SimpleNamespace

import numpy as np

from . import embedded_tddft as embedded
from .distributed import NativeParallelContext
from .errors import InputError, InputTypeError, RuntimeFailure, NumericalError
from .bands import _scalar
from .fast_hse import _Initial, _input_text, _positive, _vector
from .material_response import _response_timing
from .tddft import _copy_bundle, _digest

__all__ = ['NativeFastHSEBackend', 'NativeFastHSEModel', 'NativeFastHSERegistry']
_used = False


def _agree(comm, error, *, runtime=False):
    errors = comm.allgather(error)
    if any(errors):
        raise (RuntimeFailure if runtime else InputError)(str(errors))



def _symmetry_preflight(path, setting):
    """Read-only checks of the retained fast-PR group, before native lifetime.

    This is NOT a cell/atom/grid/k-star/restart validator. Native checks remain.
    Match symmetry.f90 record/filter semantics and exx_symmetry.f90 group rules.
    """
    def fail(reason):
        raise InputError('sym.dat: '+reason)

    if not isinstance(setting,str) or re.fullmatch(r'[yn]{1,3}',setting) is None:
        fail('yn_symmetry must contain one to three lowercase y/n characters')
    try:
        lines = Path(path).read_text().splitlines()
    except (OSError, UnicodeError) as exc:
        raise InputError('sym.dat: cannot read symmetry operations') from exc
    if not lines or len(lines) % 3:
        fail('require three records per operation')
    rows = []
    for line in lines:
        row = []
        # List-directed real literals, including Fortran D exponents/repeats.
        content = line.split('!',1)[0]
        if re.search(r',\s*,',content) or content.strip().startswith(',') or content.rstrip().endswith(','):
            fail('null values are not allowed')
        try:
            for token in content.replace(',',' ').split():
                count = 1
                if '*' in token:
                    repeat, token = token.split('*',1)
                    if re.fullmatch(r'[0-9]+',repeat) is None:
                        fail('invalid Fortran repeat count')
                    count = int(repeat)
                if count < 1 or count > 4-len(row):
                    fail('each record must contain exactly four real values')
                if re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?',token) is None:
                    fail('invalid Fortran real literal')
                value = float(token.replace('D','e').replace('d','e'))
                row.extend([value]*count)
        except (ValueError, OverflowError) as exc:
            raise InputError('sym.dat: invalid real literal') from exc
        if len(row) != 4 or not np.isfinite(row).all():
            fail('each record must contain four finite real values')
        rows.append(row)
    operations = np.asarray(rows,dtype=float).reshape(-1,3,4)
    identity = np.eye(3)
    # read_sw_symmetry initializes axes on when any y is present; each n
    # disables its own column using exact comparisons, before validation.
    disabled = [i for i,c in enumerate(setting) if c == 'n']
    retained = [op for op in operations if all(np.array_equal(op[:,i],identity[:,i]) for i in disabled)]
    if not retained:
        fail('no retained operations')
    ops = np.asarray(retained)
    rotation, translation = ops[:,:,:3],ops[:,:,3]
    tolerance = 1e-10
    for r in rotation:
        if np.max(abs(r-np.rint(r))) > tolerance or np.max(abs(r.T@r-identity)) > tolerance:
            fail('fast-PR rotations must be fractional signed permutations')
    if not any(np.max(abs(r-identity)) < tolerance and np.max(abs(t)) < tolerance
               for r,t in zip(rotation,translation)):
        fail('retained group has no identity')
    for a in range(len(ops)):
        delta = translation[a]-translation
        same_rotation = np.max(abs(rotation[a]-rotation),axis=(1,2)) < tolerance
        same_translation = np.max(abs(delta-np.rint(delta)),axis=1) < tolerance
        if np.count_nonzero(same_rotation & same_translation) != 1:
            fail('duplicate operation modulo lattice translations')
        for b in range(len(ops)):
            product = rotation[a]@rotation[b]
            shift = rotation[a]@translation[b]+translation[a]
            delta = shift-translation
            matched = (np.max(abs(rotation-product),axis=(1,2)) < tolerance) & (
                np.max(abs(delta-np.rint(delta)),axis=1) < tolerance)
            if np.count_nonzero(matched) != 1:
                fail('retained operations do not form a closed group')
    ops.setflags(write=False)
    return ops


def _configure(lib):
    vec = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, shape=(3,), flags='C_CONTIGUOUS')
    lib.sp_fast_hse_init.argtypes = [c.c_char_p, c.c_double, c.c_int, vec]
    lib.sp_fast_hse_step.argtypes = [vec, vec, vec, c.POINTER(c.c_double)]
    lib.sp_fast_hse_checkpoint.argtypes = [c.c_char_p, c.c_int]
    lib.sp_fast_hse_info.argtypes = [c.POINTER(c.c_int), vec, vec]
    lib.sp_fast_hse_close.argtypes = []
    for name in ('init','step','checkpoint','info','close'):
        getattr(lib, 'sp_fast_hse_'+name).restype = c.c_int


def _require_cell_protocol(lib,groups,symmetry,operations=None):
    if not symmetry:return
    from ._hse_geometry import _numeric
    vectors=np.column_stack([_numeric(groups.get('system',''),'al_vec'+str(i),3) for i in range(1,4)])
    general=bool(np.any(vectors-np.diag(np.diag(vectors))))
    def require(name,message):
        query=getattr(lib,name,None)
        if query is None:raise InputError(message+'; rebuild explicitly with current sources')
        query.argtypes=[];query.restype=c.c_int
        if query()!=1:raise InputError('unsupported '+message+' protocol')
    if general:
        require('sp_fast_hse_general_cell_protocol','general-cell symmetry requires a corrected fast-PR core')
    if operations is not None:
        for operation in operations:
            rotation=operation[:,:3]
            if general:rotation=vectors@rotation@np.linalg.inv(vectors)
            if np.max(abs(rotation[:,2]-[0.,0.,1.]))>1e-10:
                require('sp_fast_hse_symmetry_field_protocol','arbitrary-field symmetry requires a corrected fast-PR core')
                break


def _symmetry_factor(groups,propagator):
    body=groups.get('functional','')
    factor=_scalar(body,'exx_factor_mode',default='off')
    if factor not in ('off','strict','learned') or (propagator=='hse_taylor4_full' and factor!='off'):
        raise InputError('symmetry requires a supported fast-PR factor mode; full-action remains factor-off test-only')
    if propagator=='hse_ptcn' and factor!='off':
        from ._hse_geometry import _numeric
        if factor!='learned' or _numeric(body,'exx_factor_exact_interval',integer=True,default=8)[0]!=1:
            raise InputError('learned PT-CN requires exx_factor_mode=learned and exx_factor_exact_interval=1, matching the native update contract')
    return factor


def _history_info(lib):
    query=getattr(lib,'sp_fast_hse_history_info',None)
    if query is None:raise InputError('native history diagnostics unavailable')
    values=np.zeros(8,dtype=np.int32)
    query.argtypes=[np.ctypeslib.ndpointer(dtype=np.int32,shape=(8,),flags='C_CONTIGUOUS')]
    query.restype=c.c_int
    if query(values)!=0:raise RuntimeFailure('native history diagnostics failed')
    return values


def _require_history_protocol(lib,factor,propagator):
    if factor=='off' and propagator in ('hse_taylor4','hse_taylor4_full'):return
    query=getattr(lib,'sp_fast_hse_symmetry_history_protocol',None)
    if query is None:raise InputError('native core has no symmetry history protocol; rebuild explicitly with current sources')
    query.argtypes=[];query.restype=c.c_int
    if query()!=1:raise InputError('unsupported native symmetry history protocol')


def _initializer(lib, allow_symmetry, group_size=1):
    if not allow_symmetry:
        return lib.sp_fast_hse_init
    query = getattr(lib, 'sp_fast_hse_symmetry_protocol', None)
    entry = getattr(lib, 'sp_fast_hse_init_symmetry', None)
    if query is None or entry is None:
        raise InputError('native build has no symmetry opt-in protocol; build explicitly with current sources')
    query.argtypes = []
    query.restype = c.c_int
    if query() != 1:
        raise InputError('unsupported native symmetry opt-in protocol')
    if group_size > 1:
        kquery = getattr(lib, 'sp_fast_hse_symmetry_kmpi_protocol', None)
        if kquery is None:
            raise InputError('native build has no symmetry k-MPI protocol')
        kquery.argtypes = []
        kquery.restype = c.c_int
        if kquery() != 1:
            raise InputError('unsupported native symmetry k-MPI protocol')
    entry.argtypes = lib.sp_fast_hse_init.argtypes
    entry.restype = c.c_int
    return entry


class NativeFastHSEBackend:
    """Collective fast-PR HSE response; k-only groups of one or two ranks.

    The explicit build receipt binds a host-MPI native library. No conventional
    HSE engine or fallback is supported. Material points require separate groups;
    the legacy route has one point per group; open_model shares one verified copy
    among independent full-k states. Errors can invalidate the
    whole MPI job; native process-terminating failures are not recoverable.
    allow_symmetry is explicit, identity-bound and restricted to one/two internal
    ranks; Taylor strict/learned and PT-CN need the matching history capability.
    See FAST_HSE_SYMMETRY for numerical limitations.
    """
    supported_group_sizes = (1,2)
    execution = 'collective'
    reduction = 'leader'

    def __init__(self, build_record, work_root, *, allow_symmetry=False, shared_model=False, registry_model=False):
        if type(allow_symmetry) is not bool:
            raise InputTypeError('allow_symmetry must be boolean')
        if type(shared_model) is not bool:
            raise InputTypeError('shared_model must be boolean')
        if type(registry_model) is not bool:raise InputTypeError("registry_model must be boolean")
        if shared_model and registry_model:raise InputError("choose shared_model or registry_model")
        self.registry_model=registry_model
        self.shared_model = shared_model
        self.allow_symmetry = allow_symmetry
        self.supported_group_sizes = (1,2)
        record = json.loads(Path(build_record).expanduser().resolve().read_text())
        if record.get('schema') != 'salmonpy.fast-hse-mpi/1' or record.get('status') != 'built-and-audited':
            raise InputError('audited fast-PR HSE host-MPI build record required')
        self.library = Path(record['library']).expanduser().resolve()
        self.sha256 = record['sha256']
        self.origin_sha256 = record.get('origin_sha256')
        for value in (self.sha256, self.origin_sha256):
            if not isinstance(value,str) or re.fullmatch('[0-9a-f]{64}',value) is None:
                raise InputError('native library and fast-PR origin hashes required')
        if not self.library.is_file() or _digest(self.library) != self.sha256:
            raise InputError('fast HSE MPI library checksum mismatch')
        self.work_root = Path(work_root).expanduser().resolve()
        if not self.work_root.is_dir():
            raise InputError('existing shared work root required')
        identity = [self.sha256,self.origin_sha256]
        if allow_symmetry:
            identity.append('symmetry-v1')
        if shared_model:
            identity.append('shared-model-v1')
        if registry_model:identity.append('registry-model-v1')
        self.identity = 'salmonpy.fast-hse-mpi/1/' + hashlib.sha256(json.dumps(identity).encode()).hexdigest()

    def open_group(self, initial, dt, context):
        """Optional ResponseSession-owned factory; default route stays unchanged."""
        if getattr(self,'registry_model',False):
            resources=getattr(context,'group_resources',None)
            if resources is None:raise InputError('registry_model requires session-owned group resources')
            key='salmonpy.fast-hse.registry/1'
            registry=resources.get(key)
            if registry is None:
                registry=self.open_registry(context)
                resources[key]=registry
            elif (registry.backend.sha256,registry.backend.work_root,registry.backend.allow_symmetry)!=(self.sha256,self.work_root,self.allow_symmetry):
                raise InputError('one HSE registry group requires matching library, work root and symmetry permission')
            from ._shared_hse import _SessionGroup
            return _SessionGroup(registry.open_model(initial,dt))
        if not self.shared_model:
            return None
        from ._shared_hse import _SessionGroup
        return _SessionGroup(self.open_model(initial,dt,context))

    def open_registry(self,context):
        """Explicit different-input models sharing one host-MPI attachment."""
        from ._hse_registry import open_registry
        return open_registry(self,context)

    def open_model(self, initial, dt, context):
        """Collectively share one immutable full-k GS and MPI lease among states."""
        return self._create(initial, None, dt, context, shared=True)

    def create_state(self, initial, point_id, dt, context):
        if self.shared_model or getattr(self,"registry_model",False):
            raise InputError('shared/registry model uses ResponseSession/open_group or explicit model API')
        return self._create(initial, point_id, dt, context, shared=False)

    def _create(self, initial, point_id, dt, context, *, shared, registry=None):
        global _used
        comm = context.solver_comm
        error = None
        try:
            if registry is None and (_used or embedded._used):
                raise InputError('one native lifetime per process; use a fresh MPI job')
            if threading.current_thread() is not threading.main_thread():
                raise InputError('native fast HSE requires the main thread')
            if not shared and (len(context.point_ids) != 1 or context.point_ids[0] != point_id):
                raise InputError('one fast HSE point per solver group required')
            if comm.size not in self.supported_group_sizes:
                raise InputError('validated fast HSE k-only group requires one or two ranks')
            if not isinstance(initial, _Initial):
                raise InputTypeError('prepared FastHSEMaterialProvider state required')
            if self.work_root.is_relative_to(initial.directory):
                raise InputError('work root must be outside the initial bundle')
            dt = _positive(dt,'dt')
            if embedded._manifest(initial.directory) != initial.manifest:
                raise InputError('initial fast HSE bundle changed')
            if _digest(self.library) != self.sha256:
                raise InputError('fast HSE MPI library changed')
            groups, propagator = _input_text((initial.directory/initial.input_file).read_text())
            from ._hse_geometry import preflight, _string
            symmetry_setting = _string(groups.get('system',''), 'yn_symmetry', 'n')
            if re.fullmatch(r'[yn]{1,3}',symmetry_setting) is None:
                raise InputError('yn_symmetry requires one to three lowercase y/n characters')
            symmetry = 'y' in symmetry_setting
            if symmetry and not self.allow_symmetry:
                raise InputError('symmetry input requires explicit allow_symmetry=True')
            operations = None
            if symmetry:
                factor = _symmetry_factor(groups,propagator)
                operations = _symmetry_preflight(initial.directory/'sym.dat',
                    symmetry_setting)
            grid=preflight(initial.directory,groups,operations)
            from ._hse_assets import preflight as assets_preflight
            assets_preflight(initial.directory,groups,operations,grid,comm.size)
            groups['tgrid'] = groups.get('tgrid','') + '\n dt='+format(dt,'.17g')+', nt='+str(initial.max_steps)+'\n'
            groups['parallel'] = groups.get('parallel','') + '\n nproc_k='+str(comm.size)+', nproc_ob=1, nproc_rgrid=1,1,1\n'
            text = '\n'.join('&'+k+v+'/' for k,v in groups.items())+'\n'
            descriptor = (self.sha256,initial.manifest,point_id,dt,initial.max_steps,text,self.allow_symmetry,shared)
        except Exception as exc:
            error = str(exc)
        _agree(comm,error)
        descriptors = comm.allgather(descriptor)
        _agree(comm,'fast HSE configuration differs within subgroup' if any(x!=descriptor for x in descriptors) else None)
        work = None; error = None; expected = None
        if comm.rank == 0:
            try:
                tag = hashlib.sha256(str(point_id).encode()).hexdigest()[:12]
                parent = Path(tempfile.mkdtemp(prefix='native-fast-hse-'+tag+'-',dir=self.work_root))
                work = parent/'bundle'
                if _copy_bundle(initial.directory,work) != initial.manifest:
                    raise InputError('initial bundle changed during copy')
                (work/'inputfile').write_text(text)
                expected = embedded._manifest(work)
            except Exception as exc:
                error = str(exc)
        error, work, expected = comm.bcast((error,str(work) if work else None,expected),root=0)
        if error:
            raise InputError(error)
        work = Path(work); error = None
        try:
            if embedded._manifest(work) != expected:
                raise InputError('fast HSE work directory is not shared identically')
            lib = c.CDLL(str(self.library))
            _configure(lib)
            _require_cell_protocol(lib,groups,symmetry,operations)
            if symmetry:_require_history_protocol(lib,factor,propagator)
            if shared:
                from ._shared_hse import _configure_handles
                _configure_handles(lib,symmetry=symmetry)
            else:
                initializer = _initializer(lib,self.allow_symmetry,comm.size)
        except Exception as exc:
            error = str(exc)
        _agree(comm,error)
        owner = (context if shared and hasattr(context,'closed') and hasattr(context,'_dependents')
                 else SimpleNamespace(solver_comm=comm,closed=False,_dependents=0))
        if registry is not None:
            from ._hse_registry import _ModelLease,_ModelLibrary,_RegisteredModel
            native=_ModelLease(registry)
            model=_RegisteredModel(native,_ModelLibrary(registry),work,dt,initial.max_steps,
                                  propagator,self.identity,initial.manifest,symmetry=symmetry)
            native.model=model
            registry._models.add(model)
            return model
        native = NativeParallelContext(owner,self.library,self.sha256)
        _used = embedded._used = True
        if shared:
            return NativeFastHSEModel(native,lib,work,dt,initial.max_steps,propagator,self.identity,initial.manifest,symmetry=symmetry)
        state = _State(native,lib,work,dt,initial.max_steps,propagator)
        state.symmetry_enabled = symmetry
        try:
            current = np.empty(3)
            code = state._call(initializer,b'inputfile',dt,initial.max_steps,current)
            state.initialized = code == 0
            statuses = comm.allgather((int(code),bool(np.isfinite(current).all())))
            if any(code or not finite for code,finite in statuses):
                raise RuntimeFailure('native fast HSE initialization failed: '+str(statuses))
            state.current = current.copy()
            state.initial_current = current.copy()
            return state
        except BaseException:
            state.failed = True
            state.close()
            raise


class _State:
    checkpoint_protocol = 'salmonpy.direct-state/1'
    checkpoint_format = 'salmonpy.fast-hse-kmpi-native/2'
    checkpoint_modes = ('save','restore')

    def __init__(self,native,lib,work,dt,capacity,propagator):
        self.checkpoint_decomposition=(native.context.solver_comm.size,1,1,1,1)
        self._native = native  # Avoid the legacy sp_tddft checkpoint dispatch.
        self.lib = lib
        self.work = work
        self.dt = dt
        self.capacity = capacity
        self.steps = 0
        self.A = np.zeros(3)
        self.current = np.zeros(3)
        self.initial_current = np.zeros(3)
        self.energy = None
        self.closed = False
        self.failed = False
        self.initialized = False
        self.timing = _response_timing(dt,'native fast-PR HSE '+propagator+'; k-only MPI; charge current j')

    @property
    def time(self):
        return self.steps*self.dt

    def _call(self,fn,*args):
        from ._native_group import _directory_call
        self._native._live()
        return _directory_call(self.lib,self._native.context.solver_comm,self.work,fn,*args,minimum_protocol=2)

    def _validate(self,callback):
        error = None; result = None
        try:
            if self.closed or self.failed:
                raise InputError('fast HSE state is closed or failed')
            self._native._live()
            result = callback()
        except Exception as exc:
            error = str(exc)
        _agree(self._native.context.solver_comm,error)
        return result

    def advance(self,old,new):
        def validate():
            a = _vector(old, 'old A')
            b = _vector(new, 'new A')
            if not np.array_equal(a,self.A):
                raise InputError('stale fast HSE A endpoint')
            if self.steps>=self.capacity:
                raise InputError('fast HSE propagation capacity exhausted')
            return a,b
        a,b = self._validate(validate)
        comm = self._native.context.solver_comm
        endpoints = comm.allgather((a.tolist(),b.tolist(),self.steps))
        _agree(comm,'A endpoints differ within fast HSE subgroup' if any(x!=endpoints[0] for x in endpoints) else None)
        current = np.empty(3); energy = c.c_double()
        code = self._call(self.lib.sp_fast_hse_step,a,b,current,c.byref(energy))
        statuses = comm.allgather((int(code),bool(np.isfinite(current).all() and np.isfinite(energy.value))))
        if any(code or not finite for code,finite in statuses):
            from .embedded_tddft import _admission_rejected
            if _admission_rejected(self.lib,statuses,allowed=(34,),
                                   query_name='sp_fast_hse_symmetry_protocol'):
                raise InputError('A breaks the retained symmetry little group; accepted state unchanged; correct A and retry')
            if _admission_rejected(self.lib,statuses,allowed=(20,21,22,23,24,25,30,32,33),
                                   query_name='sp_fast_hse_recovery_protocol'):
                raise NumericalError('fast-PR HSE step failed; accepted electronic/exchange state restored without replay; retry permitted: '+str(statuses))
            self.failed = True
            raise RuntimeFailure('fast HSE native step failed: '+str(statuses))
        self.steps += 1; self.A = b.copy(); self.current = current.copy(); self.energy = energy.value
        return self.current.copy()

    def step(self,A):
        return self.advance(self.A.copy(),A)

    def checkpoint(self,path,restore):
        def validate():
            if type(restore) is not bool:
                raise InputTypeError('restore must be boolean')
            p = Path(path).expanduser().resolve()
            if any(x in str(p) for x in ('\x00','\n','\r')) or len(os.fsencode(p))>4000:
                raise InputError('invalid native checkpoint filename')
            if restore and (self.steps!=0 or not p.is_file()):
                raise InputError('restore requires a fresh state and an existing checkpoint')
            if not restore and (p.exists() or not p.parent.is_dir()):
                raise InputError('new checkpoint filename in an existing directory required')
            return p
        p = self._validate(validate)
        comm = self._native.context.solver_comm
        modes = comm.allgather((restore,self.steps))
        _agree(comm,'checkpoint epoch/mode differs within subgroup' if any(x!=modes[0] for x in modes) else None)
        code = self._call(self.lib.sp_fast_hse_checkpoint,os.fsencode(p),int(restore))
        codes = comm.allgather(int(code))
        if any(codes):
            self.failed = True
            raise RuntimeFailure('fast HSE native checkpoint failed: '+str(codes))
        if restore:
            step = c.c_int(); A = np.empty(3); current = np.empty(3)
            code = self._call(self.lib.sp_fast_hse_info,c.byref(step),A,current)
            statuses = comm.allgather((int(code),int(step.value),bool(np.isfinite(A).all() and np.isfinite(current).all())))
            if (any(code or not finite or not 0<=n<=self.capacity for code,n,finite in statuses)
                    or len({n for _,n,_ in statuses})!=1):
                self.failed = True
                raise RuntimeFailure('restored native fast HSE metadata invalid: '+str(statuses))
            self.steps = step.value; self.A = A.copy(); self.current = current.copy(); self.energy = None

    def close(self):
        if self.closed:
            return
        code = 0
        if self.initialized:
            code = self._call(self.lib.sp_fast_hse_close)
        codes = self._native.context.solver_comm.allgather(int(code))
        if any(codes):
            self.failed = True
            raise RuntimeFailure('fast HSE native close failed: '+str(codes))
        self._native.close()
        self.closed = True


from ._shared_hse import NativeFastHSEModel

from ._hse_registry import NativeFastHSERegistry
