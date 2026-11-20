"""Retained GS on explicit host MPI. Import never loads/builds a runtime."""
import ctypes as c
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import threading
from types import SimpleNamespace
import warnings
import numpy as np
from .errors import InputError, InputTypeError, RuntimeFailure, NumericalError
from .sbe import _integer

__all__=['NativeGSBackend','NativeGSSession','NativeGSState','NativeGSArray','NativeGSInfo','NativeSCFSnapshot']


def _agree(comm,error):
    errors=comm.allgather(error)
    if any(errors):raise InputError('retained GS admission: '+str(errors))


def _configure(lib):
    specs={'protocol':[], 'create':[c.c_char_p,c.POINTER(c.c_int64)],
           'info':[c.c_int64,c.c_void_p,c.c_void_p],
           'array_info':[c.c_int64,c.c_int,c.POINTER(c.c_int),c.c_void_p,c.c_void_p],
           'read':[c.c_int64,c.c_int,c.c_void_p,c.c_void_p,c.c_int64],
           'close':[c.c_int64], 'resources':[c.c_void_p]}
    for name,args in specs.items():
        fn=getattr(lib,'sp_gs_'+name,None)
        if fn is None:raise InputError('retained GS protocol missing; explicitly build a matching runtime')
        fn.argtypes=args;fn.restype=c.c_int
    if lib.sp_gs_protocol()!=1:raise InputError('unsupported retained GS protocol')


def _configure_scf(lib):
    specs={'scf_protocol':[],
           'start':[c.c_char_p,c.c_int,c.POINTER(c.c_int64)],
           'continue':[c.c_int64,c.c_int],
           'scf_info':[c.c_int64,c.c_void_p],
           'scf_snapshot':[c.c_int64,c.POINTER(c.c_int64)],
           'scf_restore':[c.c_int64,c.c_int64],
           'scf_snapshot_close':[c.c_int64,c.c_int64]}
    for name,args in specs.items():
        fn=getattr(lib,'sp_gs_'+name,None)
        if fn is None:raise InputError('active SCF protocol unavailable; explicitly use/build a matching runtime')
        fn.argtypes=args;fn.restype=c.c_int
    if lib.sp_gs_scf_protocol()!=1:raise InputError('unsupported active SCF protocol')


def _scf_disk_api(lib):
    probe=getattr(lib,'sp_gs_scf_disk_protocol',None)
    if probe is None:raise InputError('SCF disk restart requires a matching explicitly built runtime')
    probe.argtypes=[];probe.restype=c.c_int
    if probe()!=1:raise InputError('unsupported SCF disk protocol')
    if _scf_candidate_function(lib) is None:raise InputError('SCF disk restart requires protected publication support')
    lib.sp_gs_scf_bootstrap.argtypes=[c.c_char_p,c.POINTER(c.c_int64)]
    lib.sp_gs_scf_bootstrap.restype=c.c_int
    lib.sp_gs_scf_disk.argtypes=[c.c_int64,c.c_char_p,c.c_int]
    lib.sp_gs_scf_disk.restype=c.c_int
    return lib.sp_gs_scf_disk


def _scf_k_repartition_api(lib,*,hse=False):
    probe=getattr(lib,'sp_gs_scf_k_repartition_protocol',None)
    if probe is None:raise InputError('k-only SCF repartition requires a matching explicitly built runtime')
    probe.argtypes=[];probe.restype=c.c_int
    version=probe()
    if version not in (1,2) or (hse and version<2):
        raise InputError('fast HSE k repartition requires protocol2' if hse else 'unsupported SCF k repartition protocol')
    fn=lib.sp_gs_scf_k_repartition
    fn.argtypes=[c.c_int64,c.c_char_p,c.c_int,c.c_void_p];fn.restype=c.c_int
    return fn


def _scf_o_repartition_api(lib,*,hse=False):
    probe=getattr(lib,'sp_gs_scf_o_repartition_protocol',None)
    if probe is None:raise InputError('orbital-only SCF repartition requires a matching explicitly built runtime')
    probe.argtypes=[];probe.restype=c.c_int
    version=probe()
    if version not in (1,2) or (hse and version<2):
        raise InputError('fast HSE orbital repartition requires protocol2' if hse else 'unsupported SCF orbital repartition protocol')
    fn=lib.sp_gs_scf_o_repartition
    fn.argtypes=[c.c_int64,c.c_char_p,c.c_int,c.c_void_p];fn.restype=c.c_int
    return fn


def _scf_r_repartition_api(lib,*,hse=False):
    probe=getattr(lib,'sp_gs_scf_r_repartition_protocol',None)
    if probe is None:raise InputError('spatial SCF repartition requires a matching explicitly built runtime')
    probe.argtypes=[];probe.restype=c.c_int
    version=probe()
    if version not in (1,2) or (hse and version<2):
        raise InputError('fast HSE spatial repartition requires protocol2' if hse else 'unsupported SCF spatial repartition protocol')
    fn=lib.sp_gs_scf_r_repartition
    fn.argtypes=[c.c_int64,c.c_char_p,c.c_int,c.c_void_p];fn.restype=c.c_int
    return fn


def _scf_m_repartition_api(lib,*,hse=False,kspace=False):
    probe=getattr(lib,'sp_gs_scf_m_repartition_protocol',None)
    if probe is None:raise InputError('mixed SCF repartition requires a matching explicitly built runtime')
    probe.argtypes=[];probe.restype=c.c_int
    version=probe()
    if version not in (1,2,3) or (hse and version<2) or (kspace and version<3):raise InputError('unsupported SCF mixed repartition protocol')
    fn=lib.sp_gs_scf_m_repartition
    fn.argtypes=[c.c_int64,c.c_char_p,c.c_int,c.c_void_p];fn.restype=c.c_int
    return fn


def _scf_repartition_kind(old,new):
    if old[1:]==[1,1,1,1] and tuple(new[1:])==(1,1,1,1):return 'k'
    if old[0]==1 and old[2:]==[1,1,1] and new[0]==1 and tuple(new[2:])==(1,1,1):return 'o'
    if old[:2]==[1,1] and tuple(new[:2])==(1,1):return 'r'
    return 'm'


def _scf_saved_boxes(entry,mixed=False):
    box=entry.get('grid_bounds')
    if mixed and isinstance(box,list):box=box+[entry.get('orbital_bounds'),entry.get('k_bounds')]
    if (not isinstance(box,list) or len(box)!=(5 if mixed else 3) or any(not isinstance(pair,list) or len(pair)!=2
            or any(type(v) is not int or not 1<=v<2**31 for v in pair) or pair[1]<pair[0] for pair in box)):
        raise InputError('SCF repartition requires valid saved grid_bounds and orbital/k bounds')
    return box


def _scf_k_input(text,decomposition):
    # Only the generated parallel block may differ. Every physical/input byte in
    # other groups is retained; do not ignore arbitrary parallel algorithm flags.
    import re
    from .bands import _input_groups,_scalar
    groups=_input_groups(text)
    body=groups.pop('parallel','')
    match=re.fullmatch(r'\s*nproc_k\s*=\s*(\d+)\s*,\s*nproc_ob\s*=\s*(\d+)\s*,\s*nproc_rgrid\s*=\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*',body,re.I)
    if match is None or tuple(map(int,match.groups()))!=tuple(decomposition):
        raise InputError('SCF repartition requires an explicit matching parallel block')
    if _scalar(groups.get('functional',''),'xc','').lower() not in ('pz','libxc_pz','libxc_pbe','hse06'):
        raise InputError('SCF k repartition requires a supported ordinary/fast-PR HSE profile')
    return groups


def _scf_disk_manifest(directory,session,text,assets,material_id,*,repartition=False):
    directory=Path(directory).expanduser().resolve()
    record=json.loads((directory/'manifest.json').read_text())
    expected={'schema':'salmonpy.native-scf-disk/1','runtime_sha256':session.backend.sha256,
              'decomposition':list(session.backend.decomposition),'ranks':session.comm.size,
              'input_text':text,'material_id':material_id,
              'assets':{name:hashlib.sha256(data).hexdigest() for name,data in assets}}
    migrating=repartition and (record.get('ranks')!=session.comm.size or record.get('decomposition')!=list(session.backend.decomposition))
    for key,value in expected.items():
        if migrating and key in ('decomposition','ranks','input_text'):continue
        if record.get(key)!=value:raise InputError('SCF disk checkpoint incompatible: '+key)
    old_size=record.get('ranks')
    if type(old_size) is not int or not 1<=old_size<2**31:raise InputError('invalid SCF saved rank count')
    if migrating:
        old_decomposition=record.get('decomposition')
        orbital=(old_decomposition==[1,old_size,1,1,1] and tuple(session.backend.decomposition)==(1,session.comm.size,1,1,1))
        k_only=(old_decomposition==[old_size,1,1,1,1] and tuple(session.backend.decomposition)==(session.comm.size,1,1,1,1))
        spatial=(isinstance(old_decomposition,list) and len(old_decomposition)==5
                 and old_decomposition[:2]==[1,1] and math.prod(old_decomposition[2:])==old_size
                 and session.backend.decomposition[:2]==(1,1)
                 and math.prod(session.backend.decomposition[2:])==session.comm.size)
        kspace=getattr(session.backend,"hse_kspace",False)
        mixed=kspace or not (orbital or k_only or spatial)
        if kspace:orbital=k_only=spatial=False
        if (not isinstance(old_decomposition,list) or len(old_decomposition)!=5
                or any(type(v) is not int or v<1 for v in old_decomposition)
                or math.prod(old_decomposition)!=old_size):
            raise InputError('invalid saved MPI decomposition')
        new_groups=_scf_k_input(text,session.backend.decomposition)
        if _scf_k_input(record['input_text'],old_decomposition)!=new_groups:
            raise InputError('SCF repartition physical inputs differ')
        from .bands import _scalar
        hse=_scalar(new_groups['functional'],'xc','').lower()=='hse06'
        if mixed:
            if hse and not kspace and (old_decomposition[0]!=1 or session.backend.decomposition[0]!=1):
                raise InputError('mixed HSE SCF repartition requires Gamma k groups')
            _scf_m_repartition_api(session.lib,hse=hse,kspace=kspace)
        elif spatial:
            _scf_r_repartition_api(session.lib,hse=hse)
        elif orbital:
            _scf_o_repartition_api(session.lib,hse=hse)
        else:_scf_k_repartition_api(session.lib,hse=hse)
    files=record.get('files')
    if not isinstance(files,list) or len(files)!=old_size:raise InputError('SCF checkpoint rank files missing')
    for rank,entry in enumerate(files):
        if entry.get('name')!=f'rank-{rank:08d}.bin':raise InputError('invalid SCF rank filename')
    if migrating:
        if spatial or mixed:
            boxes=[_scf_saved_boxes(entry,mixed) for entry in files]
            for i,box in enumerate(boxes):
                if any(all(max(x[0],y[0])<=min(x[1],y[1]) for x,y in zip(box,other)) for other in boxes[:i]):
                    raise InputError('SCF saved spatial domains overlap')
            extent=[max(box[i][1] for box in boxes) for i in range(5 if mixed else 3)]
            if sum(math.prod(hi-lo+1 for lo,hi in box) for box in boxes)!=math.prod(extent):
                raise InputError('SCF saved spatial domains have gaps')
        else:
            bounds=[]
            for entry in files:
                pair=entry.get('orbital_bounds' if orbital else 'k_bounds')
                if not isinstance(pair,list) or len(pair)!=2 or any(type(v) is not int or not 1<=v<2**31 for v in pair) or pair[1]<pair[0]:
                    raise InputError('SCF repartition requires valid saved global orbital/k bounds')
                bounds.append(pair)
            ordered=sorted(bounds)
            if ordered[0][0]!=1 or any(b[0]!=a[1]+1 for a,b in zip(ordered,ordered[1:])):
                raise InputError('SCF saved orbital/k ranges overlap or have gaps')
            if session.comm.size>ordered[-1][1]:
                raise InputError('more orbital ranks than saved bands; choose fewer ranks' if orbital else 'more k ranks than saved k points; choose fewer k ranks')
        selected=[entry for index,entry in enumerate(files) if index%session.comm.size==session.comm.rank]
    else:selected=[files[session.comm.rank]]
    for entry in selected:
        path=directory/entry['name']
        with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
        if path.stat().st_size!=entry['bytes'] or digest!=entry['sha256']:
            raise InputError('SCF checkpoint rank file checksum mismatch')
    return directory,record


class NativeSCFSnapshot:
    """Collective native in-memory numerical snapshot, tied to its live GS owner."""
    def __init__(self,owner,token):
        self._owner=owner;self._token=token;self.closed=False
    def close(self):
        if not self._owner.closed:
            self._owner._live()
            candidate=getattr(self._owner,'_restore_candidate',None)
            owns_candidate=candidate is not None and candidate.snapshot is self
            selection=(self._owner._token,self._token,self.closed,owns_candidate)
            _agree(self._owner.session.comm,'SCF snapshot close selections differ' if any(
                v!=selection for v in self._owner.session.comm.allgather(selection)) else None)
            if self.closed:return
            if owns_candidate:candidate.close()
            self._owner.session._check(self._owner.session.lib.sp_gs_scf_snapshot_close(
                self._owner._token,self._token))
        self.closed=True
    def __reduce__(self):
        raise TypeError('native SCF snapshots are live in-memory tokens, not portable files')
    def __enter__(self):return self
    def __exit__(self,*args):self.close()


def _scf_candidate_function(lib):
    probe=getattr(lib,'sp_gs_scf_candidate_protocol',None)
    if probe is None:return None  # old runtimes retain direct execution
    probe.argtypes=[];probe.restype=c.c_int
    if probe()!=1:raise InputError('unsupported native SCF candidate protocol')
    fn=lib.sp_gs_scf_candidate
    fn.argtypes=[c.c_int64,c.c_int64,c.c_int];fn.restype=c.c_int
    return fn


class _SCFRestoreCandidate:
    """Detached native numerical payload; publish only after the host ready vote.

    The original owner/resources stay live and reserved. Memory or authenticated
    disk input supplies the numerical payload; this is not a second runnable solver. Allocation/process failures during
    native publication are outside handled staging-error recovery.
    """
    def __init__(self,state,snapshot,fn,*,directory=None):
        self.state=state;self.snapshot=snapshot;self.fn=fn;self.directory=directory
        self.ready=False;self.reserved=False;self.finished=False

    def stage(self):
        self.state._scf_ready()
        error=None
        if self.ready or self.finished:error='SCF candidate already staged or retired'
        if getattr(self.state,'_restore_candidate',None) is not None:error='SCF owner already reserved'
        if self.snapshot is not None and self.snapshot.closed:error='SCF candidate snapshot closed before staging'
        _agree(self.state.session.comm,error)
        path=None;error=None
        if self.directory is not None:
            try:
                root,_=_scf_disk_manifest(self.directory,self.state.session,self.state._input_text,
                                          self.state._assets,self.state.material_id)
                path=root/f'rank-{self.state.session.comm.rank:08d}.bin'
            except Exception as exc:error=str(exc)
            _agree(self.state.session.comm,error)
        self.state._restore_candidate=self;self.reserved=True
        if path is None:code=self.fn(self.state._token,self.snapshot._token,2)
        else:code=self.state.session.lib.sp_gs_scf_disk(self.state._token,os.fsencode(path),2)
        if code:
            # Candidate1 errors are admitted before allocation. Never discard a
            # pre-existing raw-C reservation that this failed call did not own.
            self.reserved=False;self.state._restore_candidate=None
        self.state.session._check(code)
        self.ready=True
        return self

    def discard(self):
        if self.finished or not self.reserved:return
        self.state._scf_ready()
        self.state.session._check(self.fn(self.state._token,0,4))
        self.finished=True;self.state._restore_candidate=None

    def publish(self):
        self.state._scf_ready()
        _agree(self.state.session.comm,None if self.ready and not self.finished
            else 'prepared live SCF candidate required')
        self.state.session._check(self.fn(self.state._token,0,3))
        self.finished=True;self.state._restore_candidate=None
        self.state._metadata_converged=self.state.info.converged
        return self.state

    def close(self):
        self.discard()

    def __enter__(self):return self
    def __exit__(self,*args):self.close()

    def __reduce__(self):
        raise TypeError('native SCF candidates are rank-local live resources')


@dataclass(frozen=True)
class NativeGSInfo:
    energy_hartree: float | None
    residual: float | None
    tolerance: float | None
    iterations: int
    converged: bool | None
    atoms: int
    bands: int


@dataclass(frozen=True)
class NativeGSArray:
    values: np.ndarray
    lower_indices: tuple
    quantity: str
    unit: str


class NativeGSBackend:
    def __init__(self,build_record,work_root,*,decomposition=(1,1,1,1,1),hse_orbital=False,hse_spatial=False,hse_kspace=False):
        record=json.loads(Path(build_record).read_text())
        self._ordinary_rt=type(record.get('ordinary_rt_abi')) is int and record['ordinary_rt_abi']==1
        self._fast=record.get('schema')=='salmonpy.fast-hse-mpi/1' or record.get('source_profile')=='fast-pr'
        if self._fast:
            abi=record.get('gs_abi')
            if (not isinstance(abi,dict) or type(abi.get('protocol')) is not int or abi['protocol']!=1
                    or not isinstance(abi.get('modes'),list) or 'import' not in abi['modes']
                    or not isinstance(abi.get('arrays'),list)
                    or not all(q in abi['arrays'] for q in ('density','occupations','orbitals'))):
                raise InputError('explicit fast GS import ABI declaration required')
        elif record.get('schema')!='salmonpy.embedded-tddft/1':
            raise InputError('verified embedded or declared fast GS build record required')
        self._fast_hse=record.get('schema')=='salmonpy.fast-hse-mpi/1'
        self._supports_compute=not self._fast or 'run' in record.get('gs_abi',{}).get('modes',[])
        if type(hse_orbital) is not bool or (hse_orbital and not self._fast_hse):
            raise InputError('hse_orbital must be a boolean for a fast-PR HSE runtime')
        self._hse_orbital=hse_orbital
        if type(hse_spatial) is not bool or (hse_spatial and not self._fast_hse) or (hse_spatial and hse_orbital):
            raise InputError('hse_spatial requires a fast-PR HSE runtime and excludes hse_orbital')
        self._hse_spatial=hse_spatial
        if type(hse_kspace) is not bool or (hse_kspace and (not self._fast_hse or hse_spatial or hse_orbital)):
            raise InputError('hse_kspace requires a fast-PR HSE runtime and a separate explicit profile')
        self._hse_kspace=hse_kspace
        from ._gs_archive import record_profile
        self._io_profile=record_profile(record)
        self._library=Path(record['library']).resolve();self._sha256=record['sha256']
        if hashlib.sha256(self.library.read_bytes()).hexdigest()!=self.sha256:raise InputError('native library checksum mismatch')
        self._work_root=Path(work_root).expanduser().resolve()
        if not self.work_root.is_dir():raise InputError('existing shared work root required')
        self._decomposition=tuple(_integer(v,'decomposition',1,2**31-1) for v in decomposition)
        if len(self.decomposition)!=5:raise InputError('five decomposition values required')

    @property
    def hse_kspace(self):return self._hse_kspace
    @property
    def hse_spatial(self):return self._hse_spatial
    @property
    def hse_orbital(self):return self._hse_orbital
    @property
    def library(self):return self._library
    @property
    def sha256(self):return self._sha256
    @property
    def work_root(self):return self._work_root
    @property
    def decomposition(self):return self._decomposition

    def open(self,context):
        from . import embedded_tddft
        from .distributed import NativeParallelContext
        comm=context.solver_comm;error=None
        try:
            if embedded_tddft._used:raise InputError('one electronic native attachment lifetime per process')
            if threading.current_thread() is not threading.main_thread():raise InputError('main thread required')
            if math.prod(self.decomposition)!=comm.size:raise InputError('decomposition must match subgroup size')
            if hashlib.sha256(self.library.read_bytes()).hexdigest()!=self.sha256:raise InputError('native library changed')
            lib=c.CDLL(str(self.library));_configure(lib)
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        descriptor=(self.sha256,str(self.work_root),self.decomposition)
        _agree(comm,'backend parameters differ' if any(v!=descriptor for v in comm.allgather(descriptor)) else None)
        owner=context if hasattr(context,'closed') and hasattr(context,'_dependents') else SimpleNamespace(
            solver_comm=comm,closed=False,_dependents=0)
        native=NativeParallelContext(owner,self.library,self.sha256)
        embedded_tddft._used=True
        return NativeGSSession(self,native,lib)


class NativeGSSession:
    def __init__(self,backend,native,lib):
        self.backend=backend;self._native=native;self.lib=lib;self.closed=False;self._states={};self._td_registry=None
    @property
    def comm(self):return self._native.context.solver_comm
    def _live(self):
        self._native._live()
        if self.closed:raise InputError('GS session is closed')
    def _call(self,work,fn,*args):
        self._live();previous=Path.cwd();error=None
        try:
            try:os.chdir(work)
            except Exception as exc:error=str(exc)
            _agree(self.comm,error)
            return fn(*args)
        finally:os.chdir(previous)
    def _check(self,code):
        codes=self.comm.allgather(code)
        if any(codes):raise RuntimeFailure('retained GS native call failed: '+str(codes))

    def _gs_call(self,work,fn,*args):
        boundary=None
        if getattr(self.backend,'_fast',False) and self._td_registry is not None:
            if self.backend._ordinary_rt:
                pause=self.lib.sp_tddft_park_for_gs;pause.argtypes=[];pause.restype=c.c_int
                resume=self.lib.sp_tddft_resume_after_gs;resume.argtypes=[];resume.restype=None
                self._check(pause())
                boundary=lambda mode: (resume() or 0)
            else:
                boundary=self.lib.sp_fast_hse_gs_boundary
                self._check(boundary(0))
        try:
            return self._call(work,fn,*args)
        finally:
            if boundary is not None:self._check(boundary(1))

    def run(self,calculation,*,scf=None,material_id='material',steps=None,restart=None,repartition=False):
        from .dft import DFT,SCFSettings
        from .config import _identifier
        self._live();error=None
        try:
            if not getattr(self.backend,'_supports_compute',not getattr(self.backend,'_fast',False)):
                raise InputError('this fast GS profile imports existing states; a declared compute-capable runtime is required')
            if not isinstance(calculation,DFT):raise InputTypeError('immutable DFT calculation required')
            if calculation.xc=='hse06':
                if not getattr(self.backend,'_fast_hse',False):
                    raise InputError('HSE06 GS requires a declared fast-PR HSE runtime')
                spatial_profile=getattr(self.backend,'hse_spatial',False)
                if getattr(self.backend,'hse_kspace',False):
                    probe=getattr(self.lib,'sp_gs_hse_kspace_protocol',None)
                    if probe is None:raise InputError('multi-k spatial HSE needs a matching explicit runtime')
                    probe.argtypes=[];probe.restype=c.c_int
                    version=probe()
                    if (version not in (1,2,3) or (self.backend.decomposition[1]!=1 and version<3)
                            or self.backend.decomposition[1]>calculation.states or math.prod(calculation.kmesh)<=1
                            or (calculation.symmetry is not None and version<2) or calculation.electrons!=2*calculation.states):
                        raise InputError('multi-k spatial HSE requires full mesh, occupied spin pairs and no orbital split')
                elif spatial_profile or getattr(self.backend,'hse_orbital',False) or self.backend.decomposition!=(self.comm.size,1,1,1,1):
                    name='sp_gs_hse_spatial_protocol' if spatial_profile else 'sp_gs_hse_orbital_protocol'
                    probe=getattr(self.lib,name,None)
                    if probe is None:raise InputError('fast HSE GS requires k-only MPI or a matching explicitly selected Gamma runtime')
                    probe.argtypes=[];probe.restype=c.c_int
                    version=probe()
                    valid_layout=(self.backend.decomposition[0]==1 and (version>=2 or self.backend.decomposition[1]==1)) if spatial_profile else self.backend.decomposition==(1,self.comm.size,1,1,1)
                    if (version not in (1,2) or not valid_layout
                            or (getattr(self.backend,'hse_orbital',False) and version<2)
                            or tuple(calculation.kmesh)!=(1,1,1) or calculation.symmetry is not None
                            or calculation.electrons!=2*calculation.states):
                        raise InputError('fast HSE GS requires k-only MPI or a Gamma orbital/spatial layout and a matching runtime')
            if type(repartition) is not bool or (repartition and restart is None):
                raise InputError('repartition is an explicit boolean option requiring restart')
            if restart is not None and steps is not None:
                raise InputError('restart loads the saved boundary; call continue_scf separately')
            if restart is not None:
                _configure_scf(self.lib);_scf_disk_api(self.lib)
            if steps is not None:
                steps=_integer(steps,'steps',1,2**31-2)
                if calculation.xc not in ('pz','libxc_pz','libxc_pbe','hse06'):
                    raise InputError('active SCF currently supports periodic scalar PZ/Libxc-PZ/PBE/fast-PR HSE06')
                _configure_scf(self.lib)
            scf=SCFSettings() if scf is None else scf
            _identifier(material_id,'material_id')
            text=calculation.input_text(scf)
            if calculation.xc=='hse06' and (getattr(self.backend,'hse_orbital',False) or getattr(self.backend,'hse_spatial',False)):
                text=text.replace('&functional',"&functional\n yn_hse_wannier='y',",1).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
            if calculation.xc=='hse06' and getattr(self.backend,'hse_kspace',False):
                text=text.replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
            k,ob,x,y,z=self.backend.decomposition
            text=text.replace('&system',f'&parallel\n nproc_k={k},nproc_ob={ob},nproc_rgrid={x},{y},{z}\n/\n&system',1)
            assets=calculation.input_assets
            descriptor=(material_id,text,tuple(p.sha256 for p in calculation.pseudopotentials))
            if calculation.symmetry is not None:descriptor+=(calculation.symmetry.sha256,)
            if restart is not None:
                restart,saved=_scf_disk_manifest(restart,self,text,assets,material_id,repartition=repartition)
                descriptor+=(str(restart),saved)
        except Exception as exc:error=str(exc)
        _agree(self.comm,error)
        _agree(self.comm,'DFT inputs/material IDs differ' if any(v!=(descriptor,steps) for v in self.comm.allgather((descriptor,steps))) else None)
        work=None;error=None
        if self.comm.rank==0:
            try:
                work=Path(tempfile.mkdtemp(prefix='gs-',dir=self.backend.work_root))
                (work/'gs.inp').write_text(text)
                for name,data in assets:(work/name).write_bytes(data)
            except Exception as exc:error=str(exc)
        error=self.comm.bcast(error,root=0);work=self.comm.bcast(work,root=0)
        _agree(self.comm,error)
        error=None
        try:
            if not work.is_dir() or (work/'gs.inp').read_text()!=text:
                raise InputError('GS work/input is not visible or differs on this rank')
            for name,data in assets:
                if (work/name).read_bytes()!=data:raise InputError('GS input asset differs on this rank: '+name)
            if calculation.symmetry is not None and calculation.xc=='hse06':
                from .fast_hse_mpi import _symmetry_preflight
                from ._hse_geometry import preflight
                from .bands import _input_groups
                operations=_symmetry_preflight(work/'sym.dat',calculation.symmetry.directions)
                preflight(work,_input_groups(text),operations)
        except Exception as exc:error=str(exc)
        _agree(self.comm,error)
        token=c.c_int64()
        if restart is not None:
            code=self._gs_call(work,self.lib.sp_gs_scf_bootstrap,b'gs.inp',c.byref(token))
        elif steps is None:
            code=self._gs_call(work,self.lib.sp_gs_create,b'gs.inp',c.byref(token))
        else:
            code=self._gs_call(work,self.lib.sp_gs_start,b'gs.inp',steps,c.byref(token))
        self._check(code)
        state=NativeGSState(self,token.value,work,material_id,calculation)
        state._scf_controlled=steps is not None or restart is not None
        state._input_text=text;state._input_origin='computed'
        state._electrons=calculation.electrons
        state._assets=assets
        self._states[token.value]=state
        try:
            if restart is not None:
                fn=_scf_disk_api(self.lib)
                if repartition and saved['decomposition']!=list(self.backend.decomposition):
                    kind='m' if getattr(self.backend,'hse_kspace',False) else _scf_repartition_kind(saved['decomposition'],self.backend.decomposition)
                    if kind=='m':
                        bounds=np.asarray([_scf_saved_boxes(entry,True) for entry in saved['files']],dtype=np.int32)
                        migrate=_scf_m_repartition_api(self.lib,hse=calculation.xc=='hse06',kspace=getattr(self.backend,'hse_kspace',False))
                    else:
                        key={'r':'grid_bounds','o':'orbital_bounds','k':'k_bounds'}[kind]
                        bounds=np.asarray([entry[key] for entry in saved['files']],dtype=np.int32)
                        migrate=(_scf_r_repartition_api(self.lib,hse=calculation.xc=='hse06') if kind=='r' else
                                 _scf_o_repartition_api(self.lib,hse=calculation.xc=='hse06') if kind=='o' else
                                 _scf_k_repartition_api(self.lib,hse=calculation.xc=='hse06'))
                    self._check(migrate(state._token,os.fsencode(restart),saved['ranks'],bounds.ctypes.data))
                else:
                    path=restart/f'rank-{self.comm.rank:08d}.bin'
                    self._check(fn(state._token,os.fsencode(path),2))
                candidate=_scf_candidate_function(self.lib)
                if candidate is None:raise InputError('SCF disk publication protocol missing')
                self._check(candidate(state._token,0,3))
                state._metadata_converged=state.info.converged
            info=state.info
        except BaseException:
            state.close();raise
        state._metadata_converged=info.converged
        state._identity='salmonpy.retained-gs/1/'+hashlib.sha256(json.dumps(
            [descriptor,self.backend.sha256,self.backend.decomposition,vars(info)],sort_keys=True).encode()).hexdigest()
        if not info.converged and ((steps is None and restart is None) or not state.scf_status['paused']):
            warnings.warn(f'SCF not converged after {info.iterations} iterations '
                          f'(residual={info.residual:g}, tolerance={info.tolerance:g}); retaining the last iterate',
                          RuntimeWarning,stacklevel=2)
        return state

    def load(self,directory,input_file='rt.inp',*,material_id='material'):
        """Read an existing GS checkpoint without SCF or a new restart bundle.

        The input must describe that checkpoint. Integrity checks cannot establish
        the physical provenance of unlabelled legacy binary data.
        """
        import re
        from .embedded_tddft import NativeTDDFTMaterialProvider
        from .bands import _input_groups,_scalar
        from .config import _identifier
        self._live();error=None
        try:
            _identifier(material_id,'material_id')
            provider=NativeTDDFTMaterialProvider(directory,input_file)
            initial=provider.prepare(None)
            source_input=(initial.directory/initial.input_file).read_text()
            groups=_input_groups(source_input)
            if (getattr(self.backend,'_ordinary_rt',False)
                    and _scalar(groups.get('functional',''),'xc','pz').lower()=='hse06'):
                raise InputError('HSE06 requires the fast-HSE runtime, not the HSE-disabled ordinary profile')
            electrons=int(_scalar(groups['system'],'nelec'))
            if electrons<2:raise InputError('positive occupied electronic state required')
            allowed={'calculation','control','units','system','pseudo','functional','rgrid','kgrid',
                     'atomic_coor','atomic_red_coor','scf','analysis','parallel','tgrid','propagation','emfield'}
            if set(groups)-allowed:raise InputError('unsupported checkpoint input groups: '+str(set(groups)-allowed))
            groups={key:value for key,value in groups.items() if key not in ('tgrid','propagation','emfield','analysis')}
            groups['calculation']="\n theory='dft'\n"
            if getattr(self.backend,'_fast',False):
                groups['functional']=groups.get('functional','')+"\n exx_factor_mode='off',exx_pre_scf_threshold=0\n"
            source=str(initial.directory/'data_for_restart').replace("'","''")+'/'
            groups['control']=("\n sysname='gs', yn_restart='y', read_gs_restart_data='wfn', "
                               "write_gs_restart_data='no', yn_self_checkpoint='n', base_directory='./', \n"
                               f"directory_read_data='{source}'\n")
            k,o,x,y,z=self.backend.decomposition
            groups['parallel']=f'\n nproc_k={k},nproc_ob={o},nproc_rgrid={x},{y},{z}\n'
            assets=[]
            for body in groups.values():
                clean=re.sub(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*",
                             lambda m:'' if m[0].startswith('!') else m[0],body)
                for match in re.finditer(r"\b(file_pseudo|file_atom_coor|file_atom_red_coor|file_kw|file_kgrid)\s*(?:\([^)]*\))?\s*=\s*(['\"])(.*?)\2",clean,re.I):
                    name=match[3]
                    if name=='none':continue
                    path=Path(name)
                    if path.is_absolute() or '..' in path.parts:raise InputError('bundle-relative input asset required')
                    assets.append((name,(initial.directory/path).read_bytes()))
            if getattr(self.backend,'_fast',False) and 'y' in _scalar(groups.get('system',''),'yn_symmetry','n'):
                assets.append(('sym.dat',(initial.directory/'sym.dat').read_bytes()))
            if not assets:raise InputError('explicit pseudopotential assets required')
            text=''.join(f'&{key}\n{body}\n/\n' for key,body in groups.items())
            if any(len(line)>256 for line in text.splitlines()):raise InputError('import input line exceeds native 256-character limit')
            fn=getattr(self.lib,'sp_gs_import',None)
            if fn is None:raise InputError('existing-GS import requires a matching runtime')
            fn.argtypes=[c.c_char_p,c.POINTER(c.c_int64)];fn.restype=c.c_int
            descriptor=(material_id,text,initial.manifest,tuple((n,hashlib.sha256(d).hexdigest()) for n,d in assets))
        except Exception as exc:error=str(exc)
        _agree(self.comm,error)
        _agree(self.comm,'checkpoint inputs differ' if any(v!=descriptor for v in self.comm.allgather(descriptor)) else None)
        work=None;error=None
        if self.comm.rank==0:
            try:
                work=Path(tempfile.mkdtemp(prefix='gs-import-',dir=self.backend.work_root))
                (work/'gs.inp').write_text(text)
                for name,data in assets:
                    (work/name).parent.mkdir(parents=True,exist_ok=True)
                    (work/name).write_bytes(data)
            except Exception as exc:error=str(exc)
        error=self.comm.bcast(error,root=0);work=self.comm.bcast(work,root=0)
        _agree(self.comm,error)
        try:
            provider.prepare(None)
            if (work/'gs.inp').read_text()!=text:raise InputError('import input differs on this rank')
            for name,data in assets:
                if (work/name).read_bytes()!=data:raise InputError('import asset differs on this rank')
        except Exception as exc:error=str(exc)
        _agree(self.comm,error)
        token=c.c_int64()
        self._check(self._gs_call(work,fn,b'gs.inp',c.byref(token)))
        state=NativeGSState(self,token.value,work,material_id,None)
        state._input_text=source_input;state._input_origin='imported'
        state._electrons=electrons;state._assets=tuple(assets)
        state._identity='salmonpy.imported-gs/1/'+hashlib.sha256(json.dumps(
            [descriptor,self.backend.sha256,self.backend.decomposition],sort_keys=True).encode()).hexdigest()
        self._states[token.value]=state
        return state

    def load_archive(self,archive,directory,*,converter,executable,pseudopotentials,
                     material_id=None,timeout=300.):
        """Collectively convert GS/2 once on root, then import without SCF.

        Explicit directory retains the converted bundle. This initializes arrays,
        not a complete SCF/native-state checkpoint. Conversion failures precede
        native import; the existing native loader recovery limits still apply.
        """
        from ._gs_archive import load_archive
        return load_archive(self,archive,directory,converter=converter,executable=executable,
                            pseudopotentials=pseudopotentials,material_id=material_id,timeout=timeout)

    def resources(self):
        self._live();counts=np.zeros(3,dtype=np.int32)
        self._check(self.lib.sp_gs_resources(counts.ctypes.data))
        return {'states':int(counts[0]),'live_child_communicators':int(counts[1]),'communicator_slots':int(counts[2])}

    def close(self):
        if self.closed:return
        errors=[]
        if self._td_registry is not None:
            for model in sorted(self._td_registry._models,key=lambda value:value.model_token):
                for state in sorted(model._states,key=lambda value:value.handle):
                    try:state.close()
                    except Exception as exc:errors.append(str(exc))
            try:self._td_registry.close()
            except Exception as exc:errors.append(str(exc))
        for token in sorted(self._states):
            try:self._states[token].close()
            except Exception as exc:errors.append(str(exc))
        if errors:raise RuntimeFailure('GS session cleanup failed: '+str(errors))
        self._native.close();self.closed=True
    def __enter__(self):return self
    def __exit__(self,*args):self.close()


_QUANTITIES={'density':(1,'bohr^-3'),'occupations':(2,'1'),'eigenvalues':(3,'hartree'),
             'forces':(4,'hartree/bohr'),'orbitals':(5,'native orbital amplitude'),'scf_history':(6,'mixed; see columns')}


class NativeGSState:
    def __init__(self,session,token,work,material_id,calculation):
        self.session=session;self._token=token;self.directory=work
        self._material_id=material_id;self._calculation=calculation;self._identity=None;self.closed=False
        self._electrons=None;self._assets=()
        self._input_text=None;self._input_origin='unknown';self._metadata_converged=None
        self._archive_origin=None
        self._archive_export_context=None;self._archive_export_potentials=()
    @property
    def identity(self):return self._identity
    @property
    def material_id(self):return self._material_id
    @property
    def calculation(self):return self._calculation
    @property
    def archive_origin(self):return self._archive_origin

    def input_metadata(self):
        """Return detached declared-input metadata, even after native state close.

        Optional f90nml is imported only here; no native/MPI operation is made.
        """
        from .gs_input import describe_gs_input
        if self._input_text is None:raise InputError('no captured GS input text')
        return describe_gs_input(self._input_text,assets=self._assets,material_id=self.material_id,
                                 origin=self._input_origin,converged=self._metadata_converged)

    def _live(self):
        self.session._live()
        _agree(self.session.comm,'GS state is closed' if self.closed else None)
    @property
    def info(self):
        self._live();stats=np.zeros(4,dtype=np.int32);values=np.zeros(3,dtype=np.float64)
        self.session._check(self.session.lib.sp_gs_info(self._token,stats.ctypes.data,values.ctypes.data))
        if any(self.session.comm.allgather(not bool(np.isfinite(values).all()))):
            raise NumericalError('nonfinite retained GS metadata')
        if getattr(getattr(self.session,'backend',None),'_fast',False) and stats[1]<0:
            return NativeGSInfo(None,None,None,int(stats[0]),None,int(stats[2]),int(stats[3]))
        return NativeGSInfo(*map(float,values),int(stats[0]),None if stats[1]<0 else bool(stats[1]),int(stats[2]),int(stats[3]))
    @property
    def energy_hartree(self):return self.info.energy_hartree
    @property
    def converged(self):return self.info.converged
    @property
    def iterations(self):return self.info.iterations
    def _scf_ready(self):
        self._live();error=None
        try:
            if self.calculation is None or self.calculation.xc not in ('pz','libxc_pz','libxc_pbe','hse06'):
                raise InputError('active SCF requires a computed periodic scalar PZ/Libxc-PZ/PBE/fast-PR HSE06 owner')
            _configure_scf(self.session.lib)
        except Exception as exc:error=str(exc)
        _agree(self.session.comm,error)
    @property
    def scf_status(self):
        self._scf_ready();values=np.zeros(6,dtype=np.int32)
        self.session._check(self.session.lib.sp_gs_scf_info(self._token,values.ctypes.data))
        result=dict(zip(('paused','iterations','max_iterations','initializations','executed_iterations','snapshots'),map(int,values)))
        result['paused']=bool(result['paused'])
        return result
    @property
    def checkpoint_capability(self):
        """Collective declaration for memory snapshots of this live SCF owner."""
        from .checkpointing import CheckpointCapability
        # Query the actual live profile: a computed/imported or symmetry owner
        # must not advertise capabilities that its numerical entry rejects.
        _ = self.scf_status
        identity=hashlib.sha256(json.dumps(
            [self.identity,self.material_id,self._token],sort_keys=True).encode()).hexdigest()
        return CheckpointCapability('salmonpy.native-scf-memory/1',('save','restore'),
            tuple(self.session.backend.decomposition),identity=identity)

    def checkpoint_target(self,snapshot=None,*,directory=None):
        """Bind common execution; use snapshot for memory or directory for disk.

        No snapshot is taken here. With a snapshot and candidate protocol1,
        protected staging retains independent numerical data. After the common
        ready vote, publish each returned candidate collectively. Old runtimes
        expose direct-only callbacks, never an in-place protected fallback.
        A disk binding performs no file I/O; stage authenticates files later.
        """
        from .checkpointing import CheckpointTarget
        if directory is not None:
            error=None
            if snapshot is not None:error='select a memory snapshot or a disk directory, not both'
            try:directory=Path(directory).expanduser().resolve()
            except Exception as exc:error=str(exc)
            _agree(self.session.comm,error)
            paths=self.session.comm.allgather(str(directory))
            _agree(self.session.comm,'SCF disk target paths differ' if any(p!=paths[0] for p in paths) else None)
            capability=self.scf_disk_capability
            fn=_scf_candidate_function(self.session.lib)
            candidate=_SCFRestoreCandidate(self,None,fn,directory=directory)
            return CheckpointTarget(capability,save=lambda:self.save_scf(directory),
                restore=lambda:self.restore_scf(directory),stage_restore=candidate.stage,
                discard_restore=lambda unused:candidate.discard())
        capability=self.checkpoint_capability
        error=None
        if snapshot is not None and (not isinstance(snapshot,NativeSCFSnapshot)
                or snapshot._owner is not self or snapshot.closed):
            error='SCF checkpoint target requires an open snapshot of the same live owner'
        _agree(self.session.comm,error)
        fn=None;error=None
        try:fn=_scf_candidate_function(self.session.lib)
        except Exception as exc:error=str(exc)
        _agree(self.session.comm,error)
        candidate=None if snapshot is None or fn is None else _SCFRestoreCandidate(self,snapshot,fn)
        return CheckpointTarget(capability,save=self.snapshot,
            restore=None if snapshot is None else lambda:self.restore(snapshot),
            stage_restore=None if candidate is None else candidate.stage,
            discard_restore=None if candidate is None else lambda unused:candidate.discard())

    @property
    def scf_disk_capability(self):
        """Metadata-only same-build/model/topology disk declaration; no file reads."""
        from .checkpointing import CheckpointCapability
        self._scf_ready();error=None
        try:_scf_disk_api(self.session.lib)
        except Exception as exc:error=str(exc)
        _agree(self.session.comm,error)
        identity=hashlib.sha256(json.dumps([self.material_id,self._input_text,
            {name:hashlib.sha256(data).hexdigest() for name,data in self._assets},
            self.session.backend.sha256],sort_keys=True).encode()).hexdigest()
        return CheckpointCapability('salmonpy.native-scf-disk/1',('save','restore'),
            tuple(self.session.backend.decomposition),identity=identity)

    def restore_scf(self,directory):
        """Collectively restore this live owner from authenticated disk data.

        A single owner is staged before publication. Use common protected restore
        for multiple owners; direct restore across several callbacks is non-atomic.
        """
        target=self.checkpoint_target(directory=directory)
        candidate=target.stage_restore()
        return candidate.publish()

    def continue_scf(self,*,steps):
        self._live();error=None
        try:steps=_integer(steps,'steps',1,2**31-2)
        except Exception as exc:error=str(exc)
        _agree(self.session.comm,error)
        _agree(self.session.comm,'SCF step budgets differ' if any(x!=steps for x in self.session.comm.allgather(steps)) else None)
        self._scf_ready()
        self.session._check(self.session._gs_call(self.directory,self.session.lib.sp_gs_continue,self._token,steps))
        self._metadata_converged=self.info.converged
        if not self.scf_status['paused'] and not self._metadata_converged:
            warnings.warn('SCF iteration cap reached without convergence; retaining the last iterate',RuntimeWarning,stacklevel=2)
        return self.info
    def snapshot(self):
        self._scf_ready();token=c.c_int64()
        self.session._check(self.session.lib.sp_gs_scf_snapshot(self._token,c.byref(token)))
        return NativeSCFSnapshot(self,token.value)
    def restore(self,snapshot):
        self._live();error=None
        if not isinstance(snapshot,NativeSCFSnapshot) or snapshot._owner is not self or snapshot.closed:
            error='SCF snapshot requires the same live GS owner'
        _agree(self.session.comm,error)
        self._scf_ready()
        self.session._check(self.session._gs_call(self.directory,self.session.lib.sp_gs_scf_restore,self._token,snapshot._token))
        self._metadata_converged=self.info.converged

    def write_scf_report(self,directory):
        """Save residual history as CSV and an optional logarithmic plot.

        Plotting failure warns and leaves the numerical report available.
        New output directory only; collective, with root writing the files.
        """
        self._scf_ready();comm=self.session.comm;error=None
        try:directory=Path(directory).expanduser().resolve()
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        history=self.array('scf_history').values
        info=self.info
        _agree(comm,'SCF report paths differ' if any(p!=str(directory) for p in comm.allgather(str(directory))) else None)
        report=None
        if comm.rank==0:
            try:
                directory.mkdir(parents=True,exist_ok=False)
                import csv
                with (directory/'residuals.csv').open('w',newline='') as stream:
                    writer=csv.writer(stream);writer.writerow(['iteration','residual','threshold','PBE_pre_scf'])
                    writer.writerows(history.T)
                report={'info':vars(info),'columns':['iteration','residual','threshold','PBE_pre_scf'],
                        'metric':self._input_text.split("convergence='")[-1].split("'")[0]
                            if "convergence='" in self._input_text else 'native SCF convergence criterion',
                        'plot_status':'not attempted'}
            except Exception as exc:error=str(exc)
        _agree(comm,error)
        if comm.rank==0:
            figure=None
            try:
                from matplotlib.figure import Figure
                from matplotlib.backends.backend_agg import FigureCanvasAgg
                from matplotlib.ticker import MaxNLocator
                figure=Figure(figsize=(7,6));FigureCanvasAgg(figure)
                detail,axis=figure.subplots(2,1,sharex=True)
                detail.semilogy(history[0],history[1],marker='o')
                detail.set(ylabel=report['metric'],title='Residual detail')
                axis.semilogy(history[0],history[1],marker='o',label='Residual')
                if any(value>0 for value in history[2]):axis.semilogy(history[0],history[2],linestyle='--',label='Threshold')
                axis.set(xlabel='SCF iteration',ylabel=report['metric'])
                figure.suptitle('SCF converged' if info.converged else 'SCF not converged — last iterate retained')
                for panel in (detail,axis):
                    panel.grid(True,which='both',alpha=.25)
                    panel.xaxis.set_major_locator(MaxNLocator(integer=True))
                axis.legend();figure.tight_layout()
                figure.savefig(directory/'residuals.png',dpi=150)
                report['plot_status']='saved'
            except Exception as exc:
                report['plot_status']='unavailable';report['plot_error']=str(exc)
                warnings.warn('SCF data saved; residual plot unavailable: '+str(exc),UserWarning,stacklevel=2)
            finally:
                if figure is not None:figure.clear()
            try:(directory/'report.json').write_text(json.dumps(report,indent=2)+'\n')
            except Exception as exc:
                warnings.warn('Residual CSV saved; report metadata unavailable: '+str(exc),UserWarning,stacklevel=2)
        return directory

    def save_scf(self,directory):
        """Collectively save a complete same-build/topology numerical SCF boundary.

        The manifest is committed last. Failed output stays incomplete; never
        overwrite an existing checkpoint. No SCF iterations or exchange calls.
        """
        self._scf_ready();comm=self.session.comm;error=None
        try:
            fn=_scf_disk_api(self.session.lib)
            directory=Path(directory).expanduser().resolve()
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        _agree(comm,'SCF output paths differ' if any(p!=str(directory) for p in comm.allgather(str(directory))) else None)
        if comm.rank==0:
            try:directory.mkdir(parents=True,exist_ok=False)
            except Exception as exc:error=str(exc)
        _agree(comm,error)
        path=directory/f'rank-{comm.rank:08d}.bin'
        self.session._check(fn(self._token,os.fsencode(path),0))
        shape=np.ones(7,dtype=np.int64);lower=np.ones(7,dtype=np.int64);ndim=c.c_int()
        self.session._check(self.session.lib.sp_gs_array_info(self._token,5,c.byref(ndim),shape.ctypes.data,lower.ctypes.data))
        grid_shape=np.ones(7,dtype=np.int64);grid_lower=np.ones(7,dtype=np.int64)
        self.session._check(self.session.lib.sp_gs_array_info(self._token,1,c.byref(ndim),grid_shape.ctypes.data,grid_lower.ctypes.data))
        entry=None;error=None
        try:
            with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
            entry={'name':path.name,'bytes':path.stat().st_size,'sha256':digest,
                   'k_bounds':[int(lower[5]),int(lower[5]+shape[5]-1)],
                   'orbital_bounds':[int(lower[4]),int(lower[4]+shape[4]-1)]}
            entry['grid_bounds']=[[int(grid_lower[i]),int(grid_lower[i]+grid_shape[i]-1)] for i in range(3)]
        except Exception as exc:error=str(exc)
        _agree(comm,error);files=comm.allgather(entry)
        info=self.info;progress=self.scf_status
        record={'schema':'salmonpy.native-scf-disk/1','runtime_sha256':self.session.backend.sha256,
                'decomposition':list(self.session.backend.decomposition),'ranks':comm.size,
                'material_id':self.material_id,'input_text':self._input_text,
                'assets':{name:hashlib.sha256(data).hexdigest() for name,data in self._assets},
                'files':files,'info':vars(info),'progress':progress}
        if comm.rank==0:
            try:(directory/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
            except Exception as exc:error=str(exc)
        _agree(comm,error)
        try:self.write_scf_report(directory/'diagnostics')
        except Exception as exc:
            warnings.warn('SCF checkpoint saved; residual report unavailable: '+str(exc),UserWarning,stacklevel=2)
        return directory

    def tddft_model(self,*,dt,max_steps=1000):
        from ._gs_tddft import open_model
        return open_model(self,dt,max_steps)
    def export_archive(self,directory,*,orbital_block=16):
        """Collectively export supported GS/2 arrays; root returns the archive.

        Workers return None. This is not a complete SCF/native-handle checkpoint.
        """
        if getattr(getattr(self.session,'backend',None),'_fast',False) and self.calculation is None:
            raise InputError('fast imported GS archive export requires explicit model metadata')
        from ._gs_archive import export_retained
        return export_retained(self,directory,orbital_block)

    def array(self,quantity):
        self._live();error=None
        if not isinstance(quantity,str) or quantity not in _QUANTITIES:
            error='unknown GS array quantity; supported quantity string required'
        elif getattr(getattr(self.session,'backend',None),'_fast',False) and self.calculation is None and quantity=='eigenvalues':
            error='eigenvalues were not computed by fast GS import'
        elif quantity=='forces' and self._calculation is None:
            error='forces were not computed by existing-GS import'
        _agree(self.session.comm,error)
        quantities=self.session.comm.allgather(quantity)
        _agree(self.session.comm,'GS array quantities differ across ranks'
               if any(value!=quantity for value in quantities) else None)
        if quantity=='forces' and getattr(self,'_scf_controlled',False) and self.scf_status['paused']:
            raise InputError('forces are not evaluated at a paused SCF boundary')
        code,unit=_QUANTITIES[quantity];ndim=c.c_int()
        symmetry=getattr(self._calculation,'symmetry',None)
        if quantity=='forces' and symmetry is not None and 'y' in symmetry.directions:
            query=getattr(self.session.lib,'sp_gs_cartesian_force_protocol',None)
            if query is not None:query.argtypes=[];query.restype=c.c_int
            if query is None or query()!=1:
                unit='native force components (Cartesian basis unverified)'
                warnings.warn('Runtime lacks the Cartesian symmetric-force basis declaration; values remain '
                    'available with an unverified basis. Explicitly use a corrected runtime for hartree/bohr.',
                    UserWarning,stacklevel=2)
        shape=np.ones(7,dtype=np.int64);lower=np.ones(7,dtype=np.int64)
        self.session._check(self.session.lib.sp_gs_array_info(self._token,code,c.byref(ndim),shape.ctypes.data,lower.ctypes.data))
        dims=tuple(map(int,shape[:ndim.value]));size=math.prod(dims)
        real=np.empty(size,dtype=np.float64);imag=np.empty(size,dtype=np.float64) if code==5 else None
        self.session._check(self.session.lib.sp_gs_read(self._token,code,real.ctypes.data,
                            None if imag is None else imag.ctypes.data,size))
        if imag is None:values=real
        else:
            values=np.empty(size,dtype=np.complex128)
            values.real=real;values.imag=imag
        if any(self.session.comm.allgather(not bool(np.isfinite(values).all()))):
            raise NumericalError('nonfinite retained GS array: '+quantity)
        frozen=np.frombuffer(values.tobytes(),dtype=values.dtype).reshape(dims,order='F')
        return NativeGSArray(frozen,tuple(map(int,lower[:ndim.value])),quantity,unit)
    def close(self):
        if self.session.closed and self.closed:return
        self.session._native._live()
        candidate=getattr(self,'_restore_candidate',None)
        descriptor=(self._token,self.closed,candidate is not None)
        _agree(self.session.comm,'state close selections differ' if any(
            v!=descriptor for v in self.session.comm.allgather(descriptor)) else None)
        if self.closed:return
        if candidate is not None:candidate.close()
        self.session._check(self.session.lib.sp_gs_close(self._token))
        self.closed=True;self.session._states.pop(self._token,None)
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
