"""Immutable isolated FDTD inputs and batch execution of an existing SALMON binary.

No FDTD equations are reimplemented here. All physical inputs are atomic units.
"""
from .errors import InputError, InputTypeError, RuntimeFailure
from dataclasses import dataclass, asdict, field
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import numpy as np
from .ld import LorentzDrude, LDModelConfig, LDTensorModelConfig, _positive
from .geometry import Geometry, vector
from .analysis import read_table
from .fields import AcosPulse
from .subpixel import SubpixelSmoothing, _material_data, _profile_asset, _is_ld_material


def _integer(v,name):
    if isinstance(v,bool) or not isinstance(v,(int,np.integer)) or v<1:
        raise InputError(f'{name} must be a positive integer')
    return int(v)


def _number(v): return format(float(v),'.17g')
def _vector(v): return ', '.join(_number(x) for x in v)


@dataclass(frozen=True,slots=True)
class Boundaries:
    """Faces in x-,x+,y-,y+,z-,z+ order; periodic faces occur in opposite pairs.

    Coordinates remain centered, including mixed periodic/PML domains.
    PML thickness is fixed by SALMON (8 cells)."""
    faces: tuple=('pml',)*6

    def __post_init__(self):
        faces=tuple(self.faces)
        if len(faces)!=6 or any(x not in ('pml','pec','periodic') for x in faces):
            raise InputError('six pml/pec/periodic faces required')
        for axis in range(3):
            if (faces[2*axis]=='periodic')!=(faces[2*axis+1]=='periodic'):
                raise InputError('periodic faces must be specified in opposite pairs')
        object.__setattr__(self,'faces',faces)


@dataclass(frozen=True,slots=True)
class Dielectric:
    epsilon: float=1.
    mu: float=1.
    conductivity: float=0.

    def __post_init__(self):
        for key in ('epsilon','mu','conductivity'):
            object.__setattr__(self,key,_positive(getattr(self,key),key,zero=(key=='conductivity')))
        # Restrict initial API to media no faster than vacuum; use vacuum CFL.
        if self.epsilon*self.mu<1: raise InputError('initial API requires epsilon*mu >= 1')


@dataclass(frozen=True,slots=True)
class Impulse:
    amplitude: float=1e-3
    polarization: tuple=(0.,0.,1.)

    def __post_init__(self):
        object.__setattr__(self,'amplitude',_positive(self.amplitude,'amplitude'))
        p=vector(self.polarization,'polarization')
        if not np.isclose(np.linalg.norm(p),1.,atol=1e-12,rtol=0):
            raise InputError('polarization must be a unit vector')
        object.__setattr__(self,'polarization',p)


@dataclass(frozen=True,slots=True)
class PlaneWave:
    """Acos2 pulse source on a plane normal to the positive propagation axis."""
    pulse: AcosPulse
    axis: int
    position: float

    def __post_init__(self):
        if not isinstance(self.pulse,AcosPulse): raise InputTypeError('pulse must be AcosPulse')
        if isinstance(self.axis,bool) or not isinstance(self.axis,(int,np.integer)) or self.axis not in (0,1,2): raise InputError('axis must be 0,1,2')
        object.__setattr__(self,'axis',int(self.axis))
        if self.pulse.polarization[self.axis]!=0: raise InputError('polarization must be transverse')
        if np.asarray(self.position).shape!=() or np.asarray(self.position).dtype.kind not in 'iuf' or not np.isfinite(self.position):
            raise InputError('position must be a finite real numeric scalar')
        object.__setattr__(self,'position',float(self.position))


@dataclass(frozen=True,slots=True)
class PlaneMonitor:
    """Complete transverse-periodic plane; nearest grid location is recorded."""
    axis: int
    position: float

    def __post_init__(self):
        if isinstance(self.axis,(bool,np.bool_)) or not isinstance(self.axis,(int,np.integer)) or self.axis not in (0,1,2):
            raise InputError('plane axis must be 0, 1 or 2')
        from .ld import _array
        if isinstance(self.position,(bool,np.bool_)):raise InputError('plane position must be a real scalar')
        position=_array(self.position,'plane position')
        if position.shape!=():raise InputError('plane position must be scalar')
        object.__setattr__(self,'axis',int(self.axis))
        object.__setattr__(self,'position',float(position))


@dataclass(frozen=True,slots=True)
class FDTD:
    geometry: Geometry
    materials: tuple
    dt: float
    steps: int
    boundaries: Boundaries=Boundaries()
    source: Impulse | PlaneWave=Impulse()
    observations: tuple=((0.,0.,0.),)
    observation_stride: int=1
    name: str='fdtd'
    tensor_rtol: float=1e-10
    tensor_maxiter: int=200
    smoothing: SubpixelSmoothing | None=None
    planes: tuple=()

    def __post_init__(self):
        if not isinstance(self.geometry,Geometry): raise InputTypeError('geometry must be Geometry')
        if not isinstance(self.boundaries,Boundaries): raise InputTypeError('boundaries must be Boundaries')
        if not isinstance(self.source,(Impulse,PlaneWave)): raise InputTypeError('source must be Impulse or PlaneWave')
        materials=tuple(m.config if isinstance(m,LorentzDrude) else m for m in self.materials)
        materials=tuple(m.component if isinstance(m,LDTensorModelConfig) and m.component is not None else m for m in materials)
        if len(materials)>200 or any(not isinstance(m,(Dielectric,LDModelConfig,LDTensorModelConfig)) for m in materials):
            raise InputError('up to200 Dielectric, LD or reciprocal tensor materials required')
        for m in materials:
            if isinstance(m,LDTensorModelConfig):
                if len(m.strength)>100:raise InputError('at most100 LD poles supported by SALMON')
        # Revalidate configs constructed directly rather than obtained from material.config.
        for m in materials:
            if isinstance(m,LDModelConfig):
                LorentzDrude(plasma_frequency=m.plasma_frequency,strength=m.strength,
                             resonance=m.resonance,damping=m.damping,epsilon_infinity=m.epsilon_infinity,valid_frequency_range=m.valid_frequency_range)
                if len(m.strength)>100: raise InputError('at most 100 LD poles supported by SALMON')
        materials=tuple(LorentzDrude(plasma_frequency=m.plasma_frequency,strength=m.strength,
                                    resonance=m.resonance,damping=m.damping,epsilon_infinity=m.epsilon_infinity,valid_frequency_range=m.valid_frequency_range).config
                        if isinstance(m,LDModelConfig) else m for m in materials)
        object.__setattr__(self,'materials',materials)
        if self.smoothing is not None:
            if not isinstance(self.smoothing,SubpixelSmoothing):raise InputTypeError('smoothing must be SubpixelSmoothing or None')
            carrier=self.source.pulse.omega if isinstance(self.source,PlaneWave) else None
            import warnings
            try:
                if not any(_is_ld_material(material) for material in materials):
                    raise InputError('this mode is enabled only for explicitly specified LD materials')
                for material in materials:_material_data(material)
            except InputError as exc:
                warnings.warn('subpixel smoothing not applied; retaining original material coefficients: '+str(exc),RuntimeWarning,stacklevel=2)
            else:
                self.smoothing.resolve(carrier)
                if any(len(_material_data(material)[2]) for material in materials):
                    warnings.warn('dispersive subpixel coefficients match only reference_omega; inspect response away from that frequency',RuntimeWarning,stacklevel=2)
        for m in materials:
            band=m.valid_frequency_range if isinstance(m,(LDModelConfig,LDTensorModelConfig)) else None
            if band is not None:
                import warnings
                if isinstance(self.source,Impulse):
                    warnings.warn('broadband impulse extends beyond the calibrated LD band; extrapolation is unvalidated',RuntimeWarning,stacklevel=2)
                elif not band[0]<=self.source.pulse.omega<=band[1]:
                    warnings.warn('source carrier lies outside the calibrated LD band; extrapolation is unvalidated',RuntimeWarning,stacklevel=2)
                else:
                    warnings.warn('source carrier is in the calibrated LD band; finite-pulse spectral tails are not certified. Supply an E spectrum to check_excitation.',RuntimeWarning,stacklevel=2)

        if isinstance(self.source,Impulse) and any(isinstance(m,Dielectric) for m in materials):
            raise InputError('SALMON impulse excitation requires LD materials; use PlaneWave for dielectrics')
        if isinstance(self.source,Impulse) and any(((isinstance(m,LDModelConfig) and m.epsilon_infinity!=1.0) or isinstance(m,LDTensorModelConfig)) for m in materials):
            import warnings
            warnings.warn('Non-unit LD background: impulse remains a dispersive-pole current kick. Time evolution and raw data continue, but legacy normalized optical spectra are omitted; use PlaneWave for optical response.',RuntimeWarning,stacklevel=2)
        object.__setattr__(self,'tensor_rtol',_positive(self.tensor_rtol,'tensor_rtol'))
        object.__setattr__(self,'tensor_maxiter',_integer(self.tensor_maxiter,'tensor_maxiter'))
        if self.tensor_rtol>=1:raise InputError('tensor_rtol must be less than one')
        dt=_positive(self.dt,'dt');object.__setattr__(self,'dt',dt)
        for key in ('steps','observation_stride'):
            object.__setattr__(self,key,_integer(getattr(self,key),key))
        if self.observation_stride>self.steps: raise InputError('observation_stride exceeds steps')
        cfl=1/(137.035999084*np.sqrt(sum(1/d**2 for d in self.geometry.grid.spacing)))
        fastest=min([1.0]+[m.epsilon_infinity if isinstance(m,LDModelConfig) else (1/np.max(np.sum(np.abs(np.linalg.inv(m.epsilon_infinity)),axis=1)) if m.representation=='reciprocal_tensor' else min(m.epsilon_infinity[i][i] for i in range(3))) if isinstance(m,LDTensorModelConfig) else m.epsilon*m.mu for m in materials])
        cfl*=np.sqrt(fastest)
        if dt>=cfl: raise InputError(f'dt must be below whole-domain CFL limit {cfl:.17g} a.u.')
        if any(s.material>len(materials) for s in self.geometry.shapes):
            raise InputError('shape refers to an undefined material')
        contacts=self.geometry.boundary_contacts(self.boundaries.faces)
        if contacts: raise InputError(f'matter intersects vacuum PML on {contacts}; enlarge the grid')
        observations=tuple(vector(x,'observation') for x in self.observations)
        if not observations or len(observations)>200: raise InputError('1..200 observation points required')
        coords=self.geometry.grid.coordinates
        if isinstance(self.source,PlaneWave):
            axis=self.source.axis;axis_coords=coords[axis]
            i=int(np.argmin(np.abs(axis_coords-self.source.position)))
            low=9 if self.boundaries.faces[2*axis]=='pml' else 1
            high=len(axis_coords)-10 if self.boundaries.faces[2*axis+1]=='pml' else len(axis_coords)-2
            if not axis_coords[low]<=self.source.position<=axis_coords[high]:
                raise InputError('source plane must lie inside the domain, clear of boundaries')
            if np.any(np.take(self.geometry.material_ids(),[i-1,i,i+1],axis=axis)):
                raise InputError('source plane and adjacent grid planes must be vacuum')
        if any(any(x<a[0] or x>a[-1] for x,a in zip(p,coords)) for p in observations):
            raise InputError('observation is outside the material grid')
        object.__setattr__(self,'observations',observations)
        planes=tuple(self.planes)
        for plane in planes:
            if not isinstance(plane,PlaneMonitor):raise InputTypeError('planes require PlaneMonitor values')
            if any(self.boundaries.faces[2*d:2*d+2]!=('periodic','periodic') for d in range(3) if d!=plane.axis):
                raise InputError('automatic plane colocation requires transverse periodic faces')
            axis_coords=self.geometry.grid.coordinates[plane.axis]
            if plane.position<axis_coords[0] or plane.position>axis_coords[-1]:
                raise InputError('plane position is outside the grid')
            index=int(np.argmin(np.abs(axis_coords-plane.position)))
            low=9 if self.boundaries.faces[2*plane.axis]=='pml' else 1
            high=len(axis_coords)-10 if self.boundaries.faces[2*plane.axis+1]=='pml' else len(axis_coords)-2
            if index<low or index>high:raise InputError('plane and its preceding neighbor must clear normal boundaries/PML')
        object.__setattr__(self,'planes',planes)
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]{0,39}',self.name): raise InputError('invalid calculation name')

    @property
    def config(self): return self

    @property
    def tensor_materials(self):
        """LD tensor records; isotropic one-component compression is exact.

        Dielectric compatibility records are retained. Native E/H consumes
        diagonal and reciprocal off-diagonal records without dropping entries.
        """
        return tuple(m.tensor_config if isinstance(m,LDModelConfig) else m for m in self.materials)

    @property
    def smoothing_reference(self):
        if self.smoothing is None or not any(_is_ld_material(m) for m in self.materials):return None
        try:
            for material in self.materials:_material_data(material)
        except InputError:return None
        return self.smoothing.resolve(self.source.pulse.omega if isinstance(self.source,PlaneWave) else None)

    def _plane_assets(self,assets):
        if self.planes:
            rows=['SALMONPY_NATIVE_PLANES_V1',str(len(self.planes))]
            for plane in self.planes:
                coords=self.geometry.grid.coordinates[plane.axis]
                index=int(np.argmin(np.abs(coords-plane.position)))
                rows.append(f'{plane.axis+1} {index}')
            assets['salmonpy-planes.dat']='\n'.join(rows)+'\n'
        return assets

    def input_assets(self):
        """Explicit native tensor records; no compilation or downloads."""
        reference=self.smoothing_reference
        if reference is not None:
            return self._plane_assets({'salmonpy-subpixel.dat':_profile_asset(self.materials,self.geometry,reference,self.tensor_rtol,self.tensor_maxiter)})
        assets={}
        for full in (False,True):
            entries=[(i,m) for i,m in enumerate(self.materials,1)
                     if isinstance(m,LDTensorModelConfig) and (m.representation=='reciprocal_tensor')==full]
            if not entries:continue
            rows=[('SALMONPY_LD_TENSOR 1' if full else 'SALMONPY_LD_DIAGONAL 1'),str(len(entries))]
            if full:rows.append(_number(self.tensor_rtol)+' '+str(self.tensor_maxiter))
            for index,m in entries:
                rows.append(f'{index} {len(m.strength)}')
                if full:
                    rows.append(_vector(v for row in m.epsilon_infinity for v in row))
                    rows.extend(_vector(v for row in pole for v in row) for pole in m.strength)
                else:
                    rows.append(_vector(m.epsilon_infinity[i][i] for i in range(3)))
                    rows.extend(_vector(pole[i][i] for i in range(3)) for pole in m.strength)
            assets['salmonpy-ld-tensor.dat' if full else 'salmonpy-ld-diagonal.dat']='\n'.join(rows)+'\n'
        return self._plane_assets(assets)

    def input_text(self):
        """Generate an inspectable namelist; tensor materials also need input_assets()."""
        g=self.geometry.grid
        if isinstance(self.source,Impulse):
            emfield="&emfield\n ae_shape1 = 'impulse'\n e_impulse = "+_number(self.source.amplitude)+ '\n epdir_re1 = '+_vector(self.source.polarization)+'\n/'
        else:
            p=self.source.pulse
            # FDTD uses cos-carrier Ac; TDDFT/common AcosPulse uses -sin.
            emfield="&emfield\n ae_shape1 = 'Acos2'\n I_wcm2_1 = -1\n E_amplitude1 = "+_number(p.amplitude)+ '\n omega1 = '+_number(p.omega)+ '\n tw1 = '+_number(p.duration)+ '\n t1_start = '+_number(p.start)+ '\n phi_cep1 = '+_number(p.cep+0.25)+ '\n epdir_re1 = '+_vector(p.polarization)+'\n/'
        parts=["&calculation\n theory = 'maxwell'\n/",f"&control\n sysname = '{self.name}'\n/",
               "&units\n unit_system = 'au'\n/","&system\n yn_periodic = 'n'\n/",
               emfield,
               '&analysis\n nenergy = 100\n de = 0.001\n/']
        lines=['&maxwell',' al_em = '+_vector(g.size),
               ' dl_em = '+_vector(g.spacing),' dt_em = '+_number(self.dt),f' nt_em = {self.steps}',
               " yn_make_shape = 'y'"," yn_output_shape = 'y'",f' n_s = {len(self.geometry.shapes)}',
               f' media_num = {len(self.materials)}',f' obs_num_em = {len(self.observations)}',
               f' obs_samp_em = {self.observation_stride}']
        if isinstance(self.source,PlaneWave):
            direction=[0,0,0];direction[self.source.axis]=1
            position=[0,0,0];position[self.source.axis]=self.source.position
            lines.extend([" wave_input = 'source'",' ek_dir1 = '+_vector(direction),
                          ' source_loc1 = '+_vector(position)])
        for i,face in enumerate(self.boundaries.faces):
            lines.append(f" boundary_em({i//2+1},{i%2+1}) = '{'abc' if face=='pml' else face}'")
        for i,s in enumerate(self.geometry.shapes,1):
            kind='rectangular-cylinder' if s.kind=='box' else 'ellipsoid'
            lines.extend([f" typ_s({i}) = '{kind}'",f' id_s({i}) = {s.material}',
                          f' inf_s({i},1:3) = '+_vector(s.size),f' ori_s({i},:) = '+_vector(s.center)])
        for i,tensor in enumerate(self.tensor_materials,1):
            m=tensor.component if isinstance(tensor,LDTensorModelConfig) and tensor.component is not None else tensor
            if isinstance(m,Dielectric) and self.smoothing_reference is not None:
                m=LDModelConfig(0.,(0.,),(0.,),(0.,),epsilon_infinity=m.epsilon)
            if isinstance(m,Dielectric):
                lines.extend([f" media_type({i}) = 'vacuum'",f' epsilon_em({i}) = '+_number(m.epsilon),
                              f' mu_em({i}) = '+_number(m.mu),f' sigma_em({i}) = '+_number(m.conductivity)])
            else:
                n=len(m.strength)
                kind=('ld-tensor' if m.representation=='reciprocal_tensor' else 'ld-diagonal') if isinstance(m,LDTensorModelConfig) else 'lorentz-drude'
                if self.smoothing_reference is not None:kind='ld-subpixel'
                strength=[pole[0][0] for pole in m.strength] if isinstance(m,LDTensorModelConfig) else m.strength
                background=m.epsilon_infinity[0][0] if isinstance(m,LDTensorModelConfig) else m.epsilon_infinity
                lines.extend([f" media_type({i}) = '{kind}'",f' pole_num_ld({i}) = {n}',
                              f' omega_p_ld({i}) = '+_number(m.plasma_frequency),
                              f' f_ld({i},1:{n}) = '+_vector(strength),
                              f' omega_ld({i},1:{n}) = '+_vector(m.resonance),
                              f' gamma_ld({i},1:{n}) = '+_vector(m.damping)])
                if background!=1.0:
                    lines.append(f' epsilon_em({i}) = '+_number(background))
        for i,p in enumerate(self.observations,1): lines.append(f' obs_loc_em({i},:) = '+_vector(p))
        result='\n\n'.join(parts+['\n'.join(lines+['/'])])+'\n'
        for name,text in self.input_assets().items():
            result+='! '+name+' SHA256='+hashlib.sha256(text.encode()).hexdigest()+'\n'
        return result


@dataclass(frozen=True,slots=True)
class FDTDResult:
    directory: Path
    observations: tuple
    planes: tuple=()

    def reflection_transmission(self,reference,*,front=0,back=1,window='none',drive_rtol=1e-8):
        """Use matched native planes; subtract reference E/H before flux ratios."""
        from .analysis import reflection_transmission
        if not isinstance(reference,FDTDResult):raise InputTypeError('matched vacuum FDTDResult required')
        try:incoming=self.planes[front];outgoing=self.planes[back];ref=reference.planes[front]
        except IndexError as exc:raise InputError('complete front/back native plane data required') from exc
        for plane in (outgoing,ref):
            if (plane.axis!=incoming.axis or not np.array_equal(plane.time,incoming.time)
                    or not np.array_equal(plane.area_weights,incoming.area_weights)
                    or any(not np.array_equal(a,b) for a,b in zip(plane.coordinates,incoming.coordinates))):
                raise InputError('plane geometry/sampling/clocks differ')
        if ref.position!=incoming.position:raise InputError('reference front plane differs')
        if outgoing.position==incoming.position:raise InputError('front and back planes must be distinct')
        normal=[0.,0.,0.];normal[incoming.axis]=1. if outgoing.position>incoming.position else -1.
        return reflection_transmission(incoming.time,ref.E,ref.H,incoming.E-ref.E,incoming.H-ref.H,
            outgoing.E,outgoing.H,normal=normal,units='atomic_gaussian',window=window,drive_rtol=drive_rtol,
            area_weights=incoming.area_weights)

    def energy_balance(self,reference,absorbed,volumes,*,material_ids,field_energy_change,
                       other_outgoing_energy=0.,front=0,back=1,window='none',drive_rtol=1e-8,rtol=.01,atol=0.):
        """Connect collected/collocated planes to a budget; storage/work stay explicit."""
        ports=self.reflection_transmission(reference,front=front,back=back,window=window,drive_rtol=drive_rtol)
        return energy_balance(absorbed,volumes,ports,material_ids=material_ids,
            transverse_area=ports.config['surface_area'],field_energy_change=field_energy_change,
            other_outgoing_energy=other_outgoing_energy,units='atomic_gaussian',rtol=rtol,atol=atol)


@dataclass(frozen=True,slots=True)
class FDTDMPI:
    """Spatial decomposition of native FDTD; LD states stay on owning grid ranks."""
    grid: tuple
    launcher: tuple=('mpiexec',)
    rank_arguments: tuple=('-n','{ranks}')

    def __post_init__(self):
        grid=tuple(self.grid)
        if len(grid)!=3: raise InputError('three MPI grid dimensions required')
        object.__setattr__(self,'grid',tuple(_integer(v,'MPI grid') for v in grid))
        launcher=tuple(self.launcher)
        if isinstance(self.launcher,str) or not launcher or any(not isinstance(x,str) or not x or '\x00' in x for x in launcher):
            raise InputError('launcher must be a nonempty argument tuple')
        object.__setattr__(self,'launcher',launcher)
        arguments=tuple(self.rank_arguments)
        if isinstance(self.rank_arguments,str) or any(not isinstance(x,str) or not x or '\x00' in x for x in arguments):
            raise InputError('rank_arguments must be an argument sequence (empty is allowed)')
        object.__setattr__(self,'rank_arguments',arguments)
        if self.ranks>2**31-1: raise InputError('too many MPI ranks')

    @property
    def ranks(self):
        import math
        return math.prod(self.grid)

    def command(self,executable):
        """Return argv only; never launch, interpret a shell or inspect a binary."""
        path=os.fspath(executable)
        if not isinstance(path,str) or not path or '\x00' in path:
            raise InputError('executable must be a nonempty text path')
        arguments=tuple(token.replace('{ranks}',str(self.ranks)) for token in self.rank_arguments)
        return (*self.launcher,*arguments,path)

    def input_text(self,calculation):
        if any(n//p<2 for n,p in zip(calculation.geometry.grid.counts,self.grid)):
            raise InputError('each rank needs at least two grid cells per direction')
        return calculation.input_text()+'\n&parallel\n nproc_k = 1\n nproc_ob = 1\n nproc_rgrid = '+', '.join(map(str,self.grid))+'\n/\n'


def _run_fdtd_mpi(command,inp,*,cwd,log,env,timeout):
    import signal
    with subprocess.Popen(command,stdin=subprocess.PIPE,stdout=log,stderr=subprocess.STDOUT,
                          text=True,cwd=cwd,env=env,start_new_session=True) as process:
        try: process.communicate(inp,timeout=timeout)
        except BaseException:
            try: os.killpg(process.pid,signal.SIGTERM)
            except ProcessLookupError: pass
            try: process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try: os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError: pass
                process.wait()
            raise
        return process.returncode


@dataclass(frozen=True,slots=True)
class SalmonExecutable:
    """An externally built/validated SALMON executable identified by SHA256.

    If a runtime bundle is available, use load_runtime or from_runtime to verify all files.
    A hash identifies the binary; it does not replace source/object build audits.
    """
    path: Path
    sha256: str
    mpi_enabled: bool=False
    _nlcc_force_xc: tuple = field(default=(), init=False, repr=False, compare=False)
    _runtime_validation_stats: tuple = field(default=(), init=False, repr=False, compare=False)

    @property
    def runtime_validation_report(self):
        """Local I/O verification report; no job token or scientific claims."""
        return dict(self._runtime_validation_stats)

    @property
    def nlcc_force_xc(self):
        """XC names with matching completed NLCC implementation validation."""
        return self._nlcc_force_xc

    def __post_init__(self):
        object.__setattr__(self,'path',Path(self.path).expanduser().resolve())
        if type(self.mpi_enabled) is not bool: raise InputTypeError('mpi_enabled must be boolean')
        if not re.fullmatch(r'[0-9a-f]{64}',self.sha256): raise InputError('expected lowercase SHA256')

    @classmethod
    def from_runtime(cls,directory,*,required=(),host=None,validation=None):
        """Delegate bundle integrity/capability checks to runtime.load_runtime."""
        from .runtime import load_runtime
        return load_runtime(directory,required=required,host=host,validation=validation)

    @classmethod
    def from_mpi_build(cls,build_record):
        """Load a trusted explicit runtime build record; reject serial builds."""
        record_path=Path(build_record).expanduser().resolve()
        record=json.loads(record_path.read_text())
        executable=Path(record.get('executable','')).resolve()
        if record.get('status')!='built-and-audited' or executable!=record_path.parent/'build/salmon':
            raise InputError('expected an audited runtime build.json')
        if record.get('build_configuration',{}).get('mpi_enabled') is not True:
            raise InputError('runtime was not built with MPI or lacks completed MPI capability record; rebuild explicitly')
        if hashlib.sha256(executable.read_bytes()).hexdigest()!=record.get('sha256'):
            raise InputError('MPI runtime differs from audited build')
        return cls(executable,record['sha256'],mpi_enabled=True)

    def run(self,calculation,directory,*,timeout=300.,threads=1,mpi=None):
        if not isinstance(calculation,FDTD): raise InputTypeError('calculation must be FDTD')
        timeout=_positive(timeout,'timeout');threads=_integer(threads,'threads')
        executable_bytes=self.path.read_bytes()
        if hashlib.sha256(executable_bytes).hexdigest()!=self.sha256:
            raise InputError('executable SHA256 does not match the validated binary')
        assets=calculation.input_assets()
        plane_supported=b'SALMONPY_NATIVE_PLANES_V1' in executable_bytes
        if calculation.planes and not plane_supported:
            import warnings
            warnings.warn('Runtime lacks native plane collection; FDTD calculation continues, plane diagnostics unavailable. '
                          'Explicitly obtain/build a matching runtime for plane diagnostics.',UserWarning,stacklevel=2)
        if 'salmonpy-ld-diagonal.dat' in assets and b'SALMONPY_FDTD_LD_DIAGONAL_V1' not in executable_bytes:
            raise InputError('runtime lacks diagonal tensor LD support; explicitly build or obtain a matching runtime')
        if 'salmonpy-ld-tensor.dat' in assets and b'SALMONPY_FDTD_LD_TENSOR_V1' not in executable_bytes:
            raise InputError('native runtime lacks off-diagonal tensor FDTD capability; explicitly rebuild it')
        if 'salmonpy-subpixel.dat' in assets and b'SALMONPY_SUBPIXEL_REFERENCE_V2' not in executable_bytes:
            raise InputError('native runtime lacks reference-frequency subpixel module; explicitly rebuild it')
        if mpi is not None:
            if not isinstance(mpi,FDTDMPI): raise InputTypeError('mpi must be FDTDMPI')
            if not self.mpi_enabled: raise InputError('MPI requires an audited MPI-enabled executable; use from_mpi_build')
            if os.name!='posix': raise InputError('MPI launcher currently supports macOS/Linux')
        if self.mpi_enabled:
            import sys
            imported_mpi=sys.modules.get('mpi4py.MPI')
            mpi_active=imported_mpi is not None and imported_mpi.Is_initialized() and not imported_mpi.Is_finalized()
            if mpi_active or any(k in os.environ for k in ('OMPI_COMM_WORLD_RANK','PMI_RANK','PMIX_RANK','MPI_LOCALRANKID','SLURM_PROCID')):
                raise InputError('launch FDTD MPI from a fresh serial Python process, not inside an MPI rank or initialized MPI')
        inp=calculation.input_text() if mpi is None else mpi.input_text(calculation)
        command=[str(self.path)] if mpi is None else list(mpi.command(self.path))
        directory=Path(directory).expanduser().resolve()
        try:directory.mkdir(parents=True,exist_ok=False)
        except FileExistsError as exc:
            raise FileExistsError(exc.errno,'Output directory already exists. Choose a new directory.',str(directory)) from exc
        (directory/'inputfile').write_text(inp)
        for name,text in assets.items():(directory/name).write_bytes(text.encode('utf-8'))
        record=dict(config=asdict(calculation),units='atomic_gaussian',executable=str(self.path),
                    executable_sha256=self.sha256,input_sha256=hashlib.sha256(inp.encode()).hexdigest(),
                    threads=threads,timeout=timeout,status='started',command=command,
                    mpi=None if mpi is None else asdict(mpi),mpi_ranks=1 if mpi is None else mpi.ranks,
                    ld_tensor_materials=[None if not isinstance(t,LDTensorModelConfig) else dict(representation=t.representation,independent_components=t.independent_components,axes=t.axes,parameters=asdict(t)) for t in calculation.tensor_materials],
                    input_asset_sha256={name:hashlib.sha256(text.encode()).hexdigest() for name,text in assets.items()})
        manifest=directory/'salmonpy-run.json'
        manifest.write_text(json.dumps(record,indent=2))
        try:
            with (directory/'stdout.log').open('w') as log:
                env={**os.environ,'OMP_NUM_THREADS':str(threads)}
                if mpi is None:
                    completed=subprocess.run(command,input=inp,text=True,cwd=directory,stdout=log,
                                             stderr=subprocess.STDOUT,timeout=timeout,env=env,check=False)
                    returncode=completed.returncode
                else:
                    returncode=_run_fdtd_mpi(command,inp,cwd=directory,log=log,env=env,timeout=timeout)
            output=(directory/'stdout.log').read_text()
            if returncode or 'FDTD end' not in output or 'end SALMON' not in output:
                raise RuntimeFailure(f'SALMON did not finish; see {directory / "stdout.log"}')
            if len(calculation.materials)>1 and 'SALMONPY_FDTD_INTERFACE_V2' not in output:
                raise RuntimeFailure('multiple materials require the material-ID-independent v2 runtime; rebuild or obtain a new bundle')
            record['interface_assignment']='geometric-staircase-v2' if 'SALMONPY_FDTD_INTERFACE_V2' in output else 'legacy-single-material'
            tables=tuple(read_table(directory/f'obs{i}_at_point_rt.data')
                         for i in range(1,len(calculation.observations)+1))
            expected=calculation.steps//calculation.observation_stride
            if any(len(t.values)!=expected for t in tables): raise RuntimeFailure('incomplete observation output')
            record['status']='complete'
            planes=()
            if calculation.planes and not plane_supported:
                record['plane_observations_status']='unavailable'
                record['plane_observations_error']='runtime lacks native plane collection'
            if calculation.planes and plane_supported:
                try:
                    planes=read_native_planes(directory,calculation,ranks=1 if mpi is None else mpi.ranks)
                    record['plane_observations_status']='complete'
                except Exception as exc:
                    import warnings
                    record['plane_observations_status']='unavailable'
                    record['plane_observations_error']=str(exc)
                    warnings.warn('FDTD completed; plane diagnostics unavailable: '+str(exc),UserWarning,stacklevel=2)
            return FDTDResult(directory,tables,planes)
        except BaseException as exc:
            record['status']='failed' if isinstance(exc,Exception) else 'cancelled'
            record['error']=str(exc) or type(exc).__name__
            raise
        finally:
            manifest.write_text(json.dumps(record,indent=2))


@dataclass(frozen=True, slots=True)
class EnergyBalance:
    material_ids: tuple
    absorbed_energy_density: np.ndarray
    incident_energy: float
    reflected_energy: float
    transmitted_energy: float
    material_work: float
    field_energy_change: float
    other_outgoing_energy: float
    residual: float
    relative_residual: float
    within_tolerance: bool
    spectrum: object
    config: object


def energy_balance(absorbed, volumes, ports, *, material_ids, transverse_area,
                   field_energy_change, other_outgoing_energy=0., units,
                   rtol=.01, atol=0.):
    """FDTD/multiscale energy budget, all energies Hartree, volumes bohr**3.

    absorbed is absorbed_energy_density() output, point order matches material_ids
    and volumes (one representative volume per macro point). Repeated material IDs
    are valid; each entry is an independent macro point. No MPI or solver calls.
    Work includes dispersive/electronic storage. field_energy_change must include
    instantaneous background field storage, NOT that work again. Supply lateral
    flux/PML losses not represented by the ports via other_outgoing_energy.
    A nonconserved runnable calculation warns; it is never stopped or altered.
    """
    import ctypes as c
    import warnings
    from types import MappingProxyType
    from .analysis import PortSpectrum, _energy_native
    from .fields import _readonly
    from .ld import _array
    if units != 'atomic_gaussian' or not isinstance(ports,PortSpectrum):
        raise InputError('require atomic_gaussian units and reflection_transmission() output')
    if ports.config.get('units') != units: raise InputError('port units differ')
    w = np.ascontiguousarray(_array(absorbed,'absorbed'),dtype=np.float64)
    v = np.ascontiguousarray(_array(volumes,'volumes'),dtype=np.float64)
    if w.ndim != 1 or w.size < 1 or v.shape != w.shape:
        raise InputError('absorbed and volumes require one entry per macro point')
    if isinstance(material_ids,(str,bytes)): raise InputError('material_ids must be a sequence')
    ids = tuple(material_ids)
    if len(ids) != len(w) or any(not isinstance(i,str) or not i for i in ids):
        raise InputError('one nonempty material ID per macro point required')
    if len(w) > 2**31-1: raise InputError('observation dimensions exceed native ABI')
    area = _positive(transverse_area,'transverse_area')
    observed_area = ports.config.get('surface_area')
    if observed_area is not None and not np.isclose(area,observed_area,rtol=1e-12,atol=0.):
        raise InputError('transverse_area differs from the integrated port surface area')
    tolerance = _positive(rtol,'rtol',zero=True); absolute = _positive(atol,'atol',zero=True)
    # Scalar inputs are checked before ctypes sees them; numerical balance is native.
    delta = _array(field_energy_change,'field_energy_change')
    other = _array(other_outgoing_energy,'other_outgoing_energy')
    if delta.shape != () or other.shape != (): raise InputError('field/other energy must be scalar')
    terms = np.ascontiguousarray([ports.incident_fluence,ports.reflected_fluence,
                                 ports.transmitted_fluence,float(delta),float(other)],dtype=np.float64)
    stats = np.empty(6); accepted = c.c_int()
    fn = _energy_native('sp_energy_balance',[c.c_int]+[c.c_void_p]*3+[c.c_double]*3+[c.c_void_p]*2)
    if fn(len(w),w.ctypes.data,v.ctypes.data,terms.ctypes.data,area,tolerance,absolute,stats.ctypes.data,c.byref(accepted)):
        raise InputError('native energy budget requires positive incident energy/volumes and finite arithmetic')
    if not accepted.value:
        warnings.warn(f'FDTD energy imbalance {stats[4]:.6g} hartree ({stats[5]:.3%} of incident); '
                      'inspect time coverage, collocation, volume weights, ports and stored field energy',
                      UserWarning,stacklevel=2)
    config = MappingProxyType(dict(units=units,energy_unit='hartree',point_volume_unit='bohr^3',
                    transverse_area=area,rtol=tolerance,atol=absolute,
                    formula='incident-reflected-transmitted-material_work-field_change-other_outgoing',
                    diagnostic_only=True,scope='finite-time net material work; not irreversible heat or a convergence certificate'))
    return EnergyBalance(ids,_readonly(w),*map(float,stats[:3]),float(stats[3]),float(delta),float(other),
                         float(stats[4]),float(stats[5]),bool(accepted.value),ports,config)


@dataclass(frozen=True,slots=True)
class PlaneData:
    axis: int
    position: float
    time: np.ndarray
    coordinates: tuple
    E: np.ndarray
    H: np.ndarray
    area_weights: np.ndarray
    config: object

    def flux(self,normal=None):
        from .analysis import surface_flux
        if normal is None:
            normal=[0.,0.,0.];normal[self.axis]=1.
        return surface_flux(self.time,self.E,self.H,normal=normal,area_weights=self.area_weights,units='atomic_gaussian')


def read_native_planes(directory,calculation,*,ranks=1):
    """Read owner-local stream tiles and collocate periodic Yee planes in Fortran.

    Batch fresh-run profile only. H is already at E time; no temporal averaging.
    No solver launch, compilation or MPI call; incomplete observations are rejected.
    """
    import ctypes as c
    from types import MappingProxyType
    from .analysis import _energy_native
    from .fields import _readonly
    if not isinstance(calculation,FDTD):raise InputTypeError('FDTD configuration required')
    ranks=_integer(ranks,'ranks');directory=Path(directory).expanduser().resolve()
    times=np.arange(0,calculation.steps+1,calculation.observation_stride,dtype=np.int64)
    result=[]
    for number,monitor in enumerate(calculation.planes,1):
        tangent=tuple(d for d in range(3) if d!=monitor.axis)
        shape=tuple(calculation.geometry.grid.counts[d] for d in tangent)
        raw=np.empty((2,len(times),*shape,6),dtype=np.float64)
        covered=np.zeros((2,*shape),dtype=bool);seen=set();global_box=None
        files=sorted(directory.glob(f'salmonpy_plane_{number:04d}_rank_*.bin'))
        index=int(np.argmin(abs(calculation.geometry.grid.coordinates[monitor.axis]-monitor.position)))
        for path in files:
            with path.open('rb') as stream:
                if stream.read(32).rstrip(b' \0')!=b'SALMONPY_NATIVE_PLANES_V1':raise InputError('invalid plane stream magic')
                header=stream.read(84)
                if len(header)!=84:raise InputError('truncated plane header')
                endian='<' if int.from_bytes(header[:4],'little')==16909060 else '>'
                h=np.frombuffer(header,dtype=endian+'i4').astype(np.int64)
                dtbytes=stream.read(8)
                if len(dtbytes)!=8:raise InputError('truncated plane timestep')
                dt=float(np.frombuffer(dtbytes,dtype=endian+'f8')[0])
                if (h[0]!=16909060 or h[1]!=1 or h[2]!=monitor.axis+1 or h[20]!=0
                    or h[17]!=calculation.observation_stride or h[18]!=calculation.steps or dt!=calculation.dt):
                    raise InputError('plane stream version/configuration/initial clock differs')
                if global_box is None:global_box=tuple(h[5:11])
                if tuple(h[5:11])!=global_box:raise InputError('native global plane bounds differ')
                if tuple(h[8:11]-h[5:8]+1)!=calculation.geometry.grid.counts:
                    raise InputError('native plane grid differs')
                if h[3]!=h[5+monitor.axis]+index or h[4]!=h[3]-1:raise InputError('native plane location differs')
                rank=int(h[19])
                if rank in seen or rank<0 or rank>=ranks:raise InputError('invalid/duplicate native plane owner')
                seen.add(rank)
                lo=h[11:14];hi=h[14:17]
                if np.any(lo<h[5:8]) or np.any(hi>h[8:11]) or np.any(hi<lo):raise InputError('invalid plane tile bounds')
                spans=tuple(slice(int(lo[d]-h[5+d]),int(hi[d]-h[5+d]+1)) for d in tangent)
                sizes=tuple(int(hi[d]-lo[d]+1) for d in tangent)
                owns=[q for q,location in enumerate(h[3:5]) if lo[monitor.axis]<=location<=hi[monitor.axis]]
                if not owns:raise InputError('plane stream has no owned plane')
                for q in owns:
                    if np.any(covered[(q,*spans)]):raise InputError('overlapping native plane tiles')
                    covered[(q,*spans)]=True
                size=6*sizes[0]*sizes[1]*8
                for ti,step in enumerate(times):
                    tag=stream.read(4)
                    if len(tag)!=4 or int.from_bytes(tag,'little' if endian=='<' else 'big',signed=True)!=step:
                        raise InputError('incomplete or inconsistent plane sample clock')
                    for q in owns:
                        data=stream.read(size)
                        if len(data)!=size:raise InputError('truncated plane sample')
                        tile=np.frombuffer(data,dtype=endian+'f8').reshape(sizes[1],sizes[0],6).transpose(1,0,2)
                        raw[(q,ti,*spans,slice(None))]=tile
                if stream.read(1):raise InputError('unexpected trailing plane data')
        if not np.all(covered):raise InputError('native plane coverage is incomplete')
        fields=[np.ascontiguousarray(raw[q,:,:,:,i:i+3]) for q,i in ((0,0),(0,3),(1,0),(1,3))]
        E=np.empty_like(fields[0]);H=np.empty_like(E)
        fn=_energy_native('sp_collocate_yee_plane',[c.c_int]*4+[c.c_void_p]*6)
        if fn(len(times),*shape,monitor.axis+1,*(a.ctypes.data for a in fields),E.ctypes.data,H.ctypes.data):
            raise InputError('native plane colocation rejected nonfinite fields or geometry')
        coordinates=tuple(_readonly(calculation.geometry.grid.coordinates[d]) for d in tangent)
        weights=np.full(shape,np.prod([calculation.geometry.grid.spacing[d] for d in tangent]),dtype=np.float64)
        position=float(calculation.geometry.grid.coordinates[monitor.axis][index])
        result.append(PlaneData(monitor.axis,position,_readonly(times*calculation.dt),coordinates,
            _readonly(E),_readonly(H),_readonly(weights),MappingProxyType(dict(
                units='atomic_gaussian',native_H_time='already centered at E time',
                spatial_interpolation='native Yee: two E edges, four H faces; periodic transverse wrap',
                raw_tiles=len(files),requested_position=monitor.position,source='native owner-local plane streams'))))
    return tuple(result)
