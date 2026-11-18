"""Retained standard GS on explicit host MPI. Import never loads/builds a runtime."""
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

__all__=['NativeGSBackend','NativeGSSession','NativeGSState','NativeGSArray','NativeGSInfo']


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


@dataclass(frozen=True)
class NativeGSInfo:
    energy_hartree: float
    residual: float
    tolerance: float
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
    def __init__(self,build_record,work_root,*,decomposition=(1,1,1,1,1)):
        record=json.loads(Path(build_record).read_text())
        if record.get('schema')!='salmonpy.embedded-tddft/1':raise InputError('verified standard embedded build record required')
        self._library=Path(record['library']).resolve();self._sha256=record['sha256']
        if hashlib.sha256(self.library.read_bytes()).hexdigest()!=self.sha256:raise InputError('native library checksum mismatch')
        self._work_root=Path(work_root).expanduser().resolve()
        if not self.work_root.is_dir():raise InputError('existing shared work root required')
        self._decomposition=tuple(_integer(v,'decomposition',1,2**31-1) for v in decomposition)
        if len(self.decomposition)!=5:raise InputError('five decomposition values required')

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

    def run(self,calculation,*,scf=None,material_id='material'):
        from .dft import DFT,SCFSettings
        from .config import _identifier
        self._live();error=None
        try:
            if not isinstance(calculation,DFT):raise InputTypeError('immutable DFT calculation required')
            scf=SCFSettings() if scf is None else scf
            _identifier(material_id,'material_id')
            text=calculation.input_text(scf).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
            k,ob,x,y,z=self.backend.decomposition
            text=text.replace('&system',f'&parallel\n nproc_k={k},nproc_ob={ob},nproc_rgrid={x},{y},{z}\n/\n&system',1)
            descriptor=(material_id,text,tuple(p.sha256 for p in calculation.pseudopotentials))
        except Exception as exc:error=str(exc)
        _agree(self.comm,error)
        _agree(self.comm,'DFT inputs/material IDs differ' if any(v!=descriptor for v in self.comm.allgather(descriptor)) else None)
        work=None;error=None
        if self.comm.rank==0:
            try:
                work=Path(tempfile.mkdtemp(prefix='gs-',dir=self.backend.work_root))
                (work/'gs.inp').write_text(text)
                for i,pseudo in enumerate(calculation.pseudopotentials):
                    (work/calculation._pseudo_name(i)).write_bytes(pseudo.data)
            except Exception as exc:error=str(exc)
        error=self.comm.bcast(error,root=0);work=self.comm.bcast(work,root=0)
        _agree(self.comm,error)
        error=None
        try:
            if not work.is_dir() or (work/'gs.inp').read_text()!=text:
                raise InputError('GS work/input is not visible or differs on this rank')
            for i,pseudo in enumerate(calculation.pseudopotentials):
                if hashlib.sha256((work/calculation._pseudo_name(i)).read_bytes()).hexdigest()!=pseudo.sha256:
                    raise InputError('GS pseudopotential differs on this rank')
        except Exception as exc:error=str(exc)
        _agree(self.comm,error)
        token=c.c_int64()
        code=self._call(work,self.lib.sp_gs_create,b'gs.inp',c.byref(token))
        self._check(code)
        state=NativeGSState(self,token.value,work,material_id,calculation)
        state._electrons=calculation.electrons
        state._assets=tuple((calculation._pseudo_name(i),pp.data) for i,pp in enumerate(calculation.pseudopotentials))
        self._states[token.value]=state
        try:info=state.info
        except BaseException:
            state.close();raise
        state._identity='salmonpy.retained-gs/1/'+hashlib.sha256(json.dumps(
            [descriptor,self.backend.sha256,self.backend.decomposition,vars(info)],sort_keys=True).encode()).hexdigest()
        if not info.converged:
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
            groups=_input_groups((initial.directory/initial.input_file).read_text())
            electrons=int(_scalar(groups['system'],'nelec'))
            if electrons<2:raise InputError('positive occupied electronic state required')
            allowed={'calculation','control','units','system','pseudo','functional','rgrid','kgrid',
                     'atomic_coor','atomic_red_coor','scf','analysis','parallel','tgrid','propagation','emfield'}
            if set(groups)-allowed:raise InputError('unsupported checkpoint input groups: '+str(set(groups)-allowed))
            groups={key:value for key,value in groups.items() if key not in ('tgrid','propagation','emfield','analysis')}
            groups['calculation']="\n theory='dft'\n"
            source=str(initial.directory/'data_for_restart').replace("'","''")+'/'
            groups['control']=("\n sysname='gs', yn_restart='y', read_gs_restart_data='wfn', "
                               "write_gs_restart_data='no', yn_self_checkpoint='n', base_directory='./', \n"
                               f"directory_read_data='{source}'\n")
            k,o,x,y,z=self.backend.decomposition
            groups['parallel']=f'\n nproc_k={k},nproc_ob={o},nproc_rgrid={x},{y},{z}\n'
            assets=[]
            for body in groups.values():
                clean=re.sub(r"![^\n]*",'',body)
                for match in re.finditer(r"\b(file_pseudo|file_atom_coor|file_atom_red_coor|file_kgrid)\s*(?:\([^)]*\))?\s*=\s*(['\"])(.*?)\2",clean,re.I):
                    name=match[3]
                    if name=='none':continue
                    path=Path(name)
                    if path.is_absolute() or '..' in path.parts:raise InputError('bundle-relative input asset required')
                    assets.append((name,(initial.directory/path).read_bytes()))
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
        self._check(self._call(work,fn,b'gs.inp',c.byref(token)))
        state=NativeGSState(self,token.value,work,material_id,None)
        state._electrons=electrons;state._assets=tuple(assets)
        state._identity='salmonpy.imported-gs/1/'+hashlib.sha256(json.dumps(
            [descriptor,self.backend.sha256,self.backend.decomposition],sort_keys=True).encode()).hexdigest()
        self._states[token.value]=state
        return state

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
             'forces':(4,'hartree/bohr'),'orbitals':(5,'native orbital amplitude')}


class NativeGSState:
    def __init__(self,session,token,work,material_id,calculation):
        self.session=session;self._token=token;self.directory=work
        self._material_id=material_id;self._calculation=calculation;self._identity=None;self.closed=False
        self._electrons=None;self._assets=()
    @property
    def identity(self):return self._identity
    @property
    def material_id(self):return self._material_id
    @property
    def calculation(self):return self._calculation
    def _live(self):
        self.session._live()
        _agree(self.session.comm,'GS state is closed' if self.closed else None)
    @property
    def info(self):
        self._live();stats=np.zeros(4,dtype=np.int32);values=np.zeros(3,dtype=np.float64)
        self.session._check(self.session.lib.sp_gs_info(self._token,stats.ctypes.data,values.ctypes.data))
        if any(self.session.comm.allgather(not bool(np.isfinite(values).all()))):
            raise NumericalError('nonfinite retained GS metadata')
        return NativeGSInfo(*map(float,values),int(stats[0]),None if stats[1]<0 else bool(stats[1]),int(stats[2]),int(stats[3]))
    @property
    def energy_hartree(self):return self.info.energy_hartree
    @property
    def converged(self):return self.info.converged
    @property
    def iterations(self):return self.info.iterations
    def tddft_model(self,*,dt,max_steps=1000):
        from ._gs_tddft import open_model
        return open_model(self,dt,max_steps)
    def export_archive(self,directory,*,orbital_block=16):
        """Export a serial Gamma DFT-owned GS to portable GS/2 arrays.

        Reconstruct density from orbitals. This is not a complete SCF checkpoint;
        mixing history, live native resources and forces are not serialized.
        """
        from .dft import DFT
        from .gs_archive import save_ground_state
        self._live()
        model=self.calculation
        if not isinstance(model,DFT):
            raise InputError('export requires the immutable DFT model of a retained computed GS')
        if self.session.comm.size!=1 or self.session.backend.decomposition!=(1,1,1,1,1):
            raise InputError('retained GS export currently requires a serial subgroup')
        if model.kmesh!=(1,1,1):
            raise InputError('retained GS export currently supports Gamma only')
        destination=Path(directory).expanduser()
        if destination.exists():raise FileExistsError(destination)
        info=self.info
        if info.atoms!=len(model.numbers) or info.bands!=model.states:
            raise InputError('retained GS metadata differs from the DFT model')
        density=self.array('density')
        orbitals=self.array('orbitals')
        occupations=self.array('occupations')
        eigenvalues=self.array('eigenvalues')
        if density.values.shape!=model.grid or density.lower_indices!=(1,1,1):
            raise InputError('unsupported retained periodic density grid layout')
        expected=(model.states,1,1)
        for array in (occupations,eigenvalues):
            if array.values.shape!=expected or array.lower_indices!=(1,1,1):
                raise InputError('unsupported retained orbital/k/spin layout')
        values=orbitals.values
        lower=orbitals.lower_indices
        if (values.ndim!=7 or len(lower)!=7 or values.shape[3:]!=(1,model.states,1,1)
                or lower[3:]!=(1,1,1,1)):
            raise InputError('unsupported retained spin/orbital/k/image layout')
        # Density bounds define the physical region; orbitals additionally carry halos.
        starts=tuple(density.lower_indices[i]-lower[i] for i in range(3))
        if any(start<0 or start+n>values.shape[i] for i,(start,n) in enumerate(zip(starts,model.grid))):
            raise InputError('retained orbital halos do not cover the physical grid')
        region=tuple(slice(start,start+n) for start,n in zip(starts,model.grid))
        waves=values[region+(0,slice(None),0,0)].transpose(3,2,1,0)[None]
        occ=occupations.values[:,:,0].T
        # Match the existing GS/1 exporter: preserve waves/occupations and derive rho.
        rho=np.zeros(tuple(reversed(model.grid)))
        try:
            with np.errstate(over='raise',invalid='raise'):
                for band in range(model.states):
                    rho+=occ[0,band]*(waves[0,band].real**2+waves[0,band].imag**2)
        except FloatingPointError as exc:
            raise NumericalError('nonfinite reconstructed retained GS density') from exc
        if not info.converged:
            warnings.warn('SCF not converged; exporting the retained last iterate without '
                          'certifying physical convergence',RuntimeWarning,stacklevel=2)
        provenance=dict(schema='salmonpy.retained-gs-export/1',material_id=self.material_id,
            state_identity=self.identity,target_sha256=self.session.backend.sha256,
            converged=info.converged,iterations=info.iterations,
            density_source='reconstructed_from_native_orbitals',
            local_channels=[dict(atomic_number=p.atomic_number,local_channel=p.local_channel)
                            for p in model.pseudopotentials])
        references=[dict(atomic_number=p.atomic_number,filename=model._pseudo_name(i),
                         sha256=p.sha256,source='retained DFT model')
                    for i,p in enumerate(model.pseudopotentials)]
        return save_ground_state(destination,wavefunctions=waves,density=rho,
            eigenvalues=eigenvalues.values[:,:,0].T,occupations=occ,kpoints=[[0.,0.,0.]],weights=[1.],
            cell_bohr=model.cell,fractional_positions=model.fractional_positions,
            atomic_numbers=list(model.numbers),electrons=model.electrons,energy_hartree=info.energy_hartree,
            xc=model.xc,pseudopotentials=references,orbital_block=orbital_block,provenance=provenance)

    def array(self,quantity):
        self._live();error=None
        if quantity not in _QUANTITIES:error='unknown GS array quantity'
        if quantity=='forces' and self._calculation is None:error='forces were not computed by existing-GS import'
        _agree(self.session.comm,error)
        code,unit=_QUANTITIES[quantity];ndim=c.c_int()
        shape=np.ones(7,dtype=np.int64);lower=np.ones(7,dtype=np.int64)
        self.session._check(self.session.lib.sp_gs_array_info(self._token,code,c.byref(ndim),shape.ctypes.data,lower.ctypes.data))
        dims=tuple(map(int,shape[:ndim.value]));size=math.prod(dims)
        real=np.empty(size,dtype=np.float64);imag=np.empty(size,dtype=np.float64) if code==5 else None
        self.session._check(self.session.lib.sp_gs_read(self._token,code,real.ctypes.data,
                            None if imag is None else imag.ctypes.data,size))
        values=real if imag is None else real+1j*imag
        if any(self.session.comm.allgather(not bool(np.isfinite(values).all()))):
            raise NumericalError('nonfinite retained GS array: '+quantity)
        frozen=np.frombuffer(values.tobytes(),dtype=values.dtype).reshape(dims,order='F')
        return NativeGSArray(frozen,tuple(map(int,lower[:ndim.value])),quantity,unit)
    def close(self):
        if self.session.closed and self.closed:return
        self.session._native._live()
        descriptor=(self._token,self.closed)
        _agree(self.session.comm,'state close selections differ' if any(
            v!=descriptor for v in self.session.comm.allgather(descriptor)) else None)
        if self.closed:return
        self.session._check(self.session.lib.sp_gs_close(self._token))
        self.closed=True;self.session._states.pop(self._token,None)
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
