"""Immutable isolated FDTD inputs and batch execution of an existing SALMON binary.

No FDTD equations are reimplemented here. All physical inputs are atomic units.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import numpy as np
from .ld import LorentzDrude, LDModelConfig, _positive
from .geometry import Geometry, vector
from .analysis import read_table
from .fields import AcosPulse


def _integer(v,name):
    if isinstance(v,bool) or not isinstance(v,(int,np.integer)) or v<1:
        raise ValueError(f'{name} must be a positive integer')
    return int(v)


def _number(v): return format(float(v),'.17g')
def _vector(v): return ', '.join(_number(x) for x in v)


@dataclass(frozen=True,slots=True)
class Boundaries:
    """Faces in x-,x+,y-,y+,z-,z+ order. PML thickness is fixed by SALMON (8 cells)."""
    faces: tuple=('pml',)*6

    def __post_init__(self):
        faces=tuple(self.faces)
        if len(faces)!=6 or any(x not in ('pml','pec') for x in faces):
            raise ValueError('six pml/pec faces required; periodic support is not yet exposed')
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
        if self.epsilon*self.mu<1: raise ValueError('initial API requires epsilon*mu >= 1')


@dataclass(frozen=True,slots=True)
class Impulse:
    amplitude: float=1e-3
    polarization: tuple=(0.,0.,1.)

    def __post_init__(self):
        object.__setattr__(self,'amplitude',_positive(self.amplitude,'amplitude'))
        p=vector(self.polarization,'polarization')
        if not np.isclose(np.linalg.norm(p),1.,atol=1e-12,rtol=0):
            raise ValueError('polarization must be a unit vector')
        object.__setattr__(self,'polarization',p)


@dataclass(frozen=True,slots=True)
class PlaneWave:
    """Acos2 pulse source on a plane normal to the positive propagation axis."""
    pulse: AcosPulse
    axis: int
    position: float

    def __post_init__(self):
        if not isinstance(self.pulse,AcosPulse): raise TypeError('pulse must be AcosPulse')
        if isinstance(self.axis,bool) or not isinstance(self.axis,(int,np.integer)) or self.axis not in (0,1,2): raise ValueError('axis must be 0,1,2')
        object.__setattr__(self,'axis',int(self.axis))
        if self.pulse.polarization[self.axis]!=0: raise ValueError('polarization must be transverse')
        if np.asarray(self.position).shape!=() or np.asarray(self.position).dtype.kind not in 'iuf' or not np.isfinite(self.position):
            raise ValueError('position must be a finite real numeric scalar')
        object.__setattr__(self,'position',float(self.position))


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

    def __post_init__(self):
        if not isinstance(self.geometry,Geometry): raise TypeError('geometry must be Geometry')
        if not isinstance(self.boundaries,Boundaries): raise TypeError('boundaries must be Boundaries')
        if not isinstance(self.source,(Impulse,PlaneWave)): raise TypeError('source must be Impulse or PlaneWave')
        materials=tuple(m.config if isinstance(m,LorentzDrude) else m for m in self.materials)
        if len(materials)>200 or any(not isinstance(m,(Dielectric,LDModelConfig)) for m in materials):
            raise ValueError('up to 200 Dielectric or LD materials required')
        # Revalidate configs constructed directly rather than obtained from material.config.
        for m in materials:
            if isinstance(m,LDModelConfig):
                LorentzDrude(plasma_frequency=m.plasma_frequency,strength=m.strength,
                             resonance=m.resonance,damping=m.damping)
                if len(m.strength)>100: raise ValueError('at most 100 LD poles supported by SALMON')
        materials=tuple(LorentzDrude(plasma_frequency=m.plasma_frequency,strength=m.strength,
                                    resonance=m.resonance,damping=m.damping).config
                        if isinstance(m,LDModelConfig) else m for m in materials)
        object.__setattr__(self,'materials',materials)
        if isinstance(self.source,Impulse) and any(isinstance(m,Dielectric) for m in materials):
            raise ValueError('SALMON impulse excitation requires LD materials; use PlaneWave for dielectrics')
        dt=_positive(self.dt,'dt');object.__setattr__(self,'dt',dt)
        for key in ('steps','observation_stride'):
            object.__setattr__(self,key,_integer(getattr(self,key),key))
        if self.observation_stride>self.steps: raise ValueError('observation_stride exceeds steps')
        cfl=1/(137.035999084*np.sqrt(sum(1/d**2 for d in self.geometry.grid.spacing)))
        if dt>=cfl: raise ValueError(f'dt must be below vacuum CFL limit {cfl:.17g} a.u.')
        if any(s.material>len(materials) for s in self.geometry.shapes):
            raise ValueError('shape refers to an undefined material')
        contacts=self.geometry.boundary_contacts(self.boundaries.faces)
        if contacts: raise ValueError(f'matter intersects vacuum PML on {contacts}; enlarge the grid')
        observations=tuple(vector(x,'observation') for x in self.observations)
        if not observations or len(observations)>200: raise ValueError('1..200 observation points required')
        coords=self.geometry.grid.coordinates
        if isinstance(self.source,PlaneWave):
            axis=self.source.axis;axis_coords=coords[axis]
            i=int(np.argmin(np.abs(axis_coords-self.source.position)))
            low=9 if self.boundaries.faces[2*axis]=='pml' else 1
            high=len(axis_coords)-10 if self.boundaries.faces[2*axis+1]=='pml' else len(axis_coords)-2
            if not axis_coords[low]<=self.source.position<=axis_coords[high]:
                raise ValueError('source plane must lie inside the domain, clear of boundaries')
            if np.any(np.take(self.geometry.material_ids(),[i-1,i,i+1],axis=axis)):
                raise ValueError('source plane and adjacent grid planes must be vacuum')
        if any(any(x<a[0] or x>a[-1] for x,a in zip(p,coords)) for p in observations):
            raise ValueError('observation is outside the material grid')
        object.__setattr__(self,'observations',observations)
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]{0,39}',self.name): raise ValueError('invalid calculation name')

    @property
    def config(self): return self

    def input_text(self):
        """Generate a complete, inspectable SALMON namelist in atomic units."""
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
        for i,m in enumerate(self.materials,1):
            if isinstance(m,Dielectric):
                lines.extend([f" media_type({i}) = 'vacuum'",f' epsilon_em({i}) = '+_number(m.epsilon),
                              f' mu_em({i}) = '+_number(m.mu),f' sigma_em({i}) = '+_number(m.conductivity)])
            else:
                n=len(m.strength)
                lines.extend([f" media_type({i}) = 'lorentz-drude'",f' pole_num_ld({i}) = {n}',
                              f' omega_p_ld({i}) = '+_number(m.plasma_frequency),
                              f' f_ld({i},1:{n}) = '+_vector(m.strength),
                              f' omega_ld({i},1:{n}) = '+_vector(m.resonance),
                              f' gamma_ld({i},1:{n}) = '+_vector(m.damping)])
        for i,p in enumerate(self.observations,1): lines.append(f' obs_loc_em({i},:) = '+_vector(p))
        return '\n\n'.join(parts+['\n'.join(lines+['/'])])+'\n'


@dataclass(frozen=True,slots=True)
class FDTDResult:
    directory: Path
    observations: tuple


@dataclass(frozen=True,slots=True)
class FDTDMPI:
    """Spatial decomposition of native FDTD; LD states stay on owning grid ranks."""
    grid: tuple
    launcher: tuple=('mpiexec',)

    def __post_init__(self):
        grid=tuple(self.grid)
        if len(grid)!=3: raise ValueError('three MPI grid dimensions required')
        object.__setattr__(self,'grid',tuple(_integer(v,'MPI grid') for v in grid))
        launcher=tuple(self.launcher)
        if isinstance(self.launcher,str) or not launcher or any(not isinstance(x,str) or not x or '\x00' in x for x in launcher):
            raise ValueError('launcher must be a nonempty argument tuple')
        object.__setattr__(self,'launcher',launcher)
        if self.ranks>2**31-1: raise ValueError('too many MPI ranks')

    @property
    def ranks(self):
        import math
        return math.prod(self.grid)

    def input_text(self,calculation):
        if any(n//p<2 for n,p in zip(calculation.geometry.grid.counts,self.grid)):
            raise ValueError('each rank needs at least two grid cells per direction')
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

    A hash identifies the binary; it does not replace source/object build audits.
    """
    path: Path
    sha256: str
    mpi_enabled: bool=False

    def __post_init__(self):
        object.__setattr__(self,'path',Path(self.path).expanduser().resolve())
        if type(self.mpi_enabled) is not bool: raise TypeError('mpi_enabled must be boolean')
        if not re.fullmatch(r'[0-9a-f]{64}',self.sha256): raise ValueError('expected lowercase SHA256')

    @classmethod
    def from_mpi_build(cls,build_record):
        """Load a trusted explicit runtime build record; reject serial builds."""
        record_path=Path(build_record).expanduser().resolve()
        record=json.loads(record_path.read_text())
        executable=Path(record.get('executable','')).resolve()
        if record.get('status')!='built-and-audited' or executable!=record_path.parent/'build/salmon':
            raise ValueError('expected an audited runtime build.json')
        if record.get('build_configuration',{}).get('mpi_enabled') is not True:
            raise ValueError('runtime was not built with MPI or lacks completed MPI capability record; rebuild explicitly')
        if hashlib.sha256(executable.read_bytes()).hexdigest()!=record.get('sha256'):
            raise ValueError('MPI runtime differs from audited build')
        return cls(executable,record['sha256'],mpi_enabled=True)

    def run(self,calculation,directory,*,timeout=300.,threads=1,mpi=None):
        if not isinstance(calculation,FDTD): raise TypeError('calculation must be FDTD')
        timeout=_positive(timeout,'timeout');threads=_integer(threads,'threads')
        if hashlib.sha256(self.path.read_bytes()).hexdigest()!=self.sha256:
            raise ValueError('executable SHA256 does not match the validated binary')
        if mpi is not None:
            if not isinstance(mpi,FDTDMPI): raise TypeError('mpi must be FDTDMPI')
            if not self.mpi_enabled: raise ValueError('MPI requires an audited MPI-enabled executable; use from_mpi_build')
            if os.name!='posix': raise ValueError('MPI launcher currently supports macOS/Linux')
        if self.mpi_enabled:
            import sys
            imported_mpi=sys.modules.get('mpi4py.MPI')
            mpi_active=imported_mpi is not None and imported_mpi.Is_initialized() and not imported_mpi.Is_finalized()
            if mpi_active or any(k in os.environ for k in ('OMPI_COMM_WORLD_RANK','PMI_RANK','PMIX_RANK','MPI_LOCALRANKID','SLURM_PROCID')):
                raise ValueError('launch FDTD MPI from a fresh serial Python process, not inside an MPI rank or initialized MPI')
        inp=calculation.input_text() if mpi is None else mpi.input_text(calculation)
        command=[str(self.path)] if mpi is None else [*mpi.launcher,'-n',str(mpi.ranks),str(self.path)]
        directory=Path(directory).expanduser().resolve();directory.mkdir(parents=True,exist_ok=False)
        (directory/'inputfile').write_text(inp)
        record=dict(config=asdict(calculation),units='atomic_gaussian',executable=str(self.path),
                    executable_sha256=self.sha256,input_sha256=hashlib.sha256(inp.encode()).hexdigest(),
                    threads=threads,timeout=timeout,status='started',command=command,
                    mpi=None if mpi is None else asdict(mpi),mpi_ranks=1 if mpi is None else mpi.ranks)
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
                raise RuntimeError(f'SALMON did not finish; see {directory / "stdout.log"}')
            if len(calculation.materials)>1 and 'SALMONPY_FDTD_INTERFACE_V2' not in output:
                raise RuntimeError('multiple materials require the material-ID-independent v2 runtime; rebuild or obtain a new bundle')
            record['interface_assignment']='geometric-staircase-v2' if 'SALMONPY_FDTD_INTERFACE_V2' in output else 'legacy-single-material'
            tables=tuple(read_table(directory/f'obs{i}_at_point_rt.data')
                         for i in range(1,len(calculation.observations)+1))
            expected=calculation.steps//calculation.observation_stride
            if any(len(t.values)!=expected for t in tables): raise RuntimeError('incomplete observation output')
            record['status']='complete'
            return FDTDResult(directory,tables)
        except BaseException as exc:
            record['status']='failed' if isinstance(exc,Exception) else 'cancelled'
            record['error']=str(exc) or type(exc).__name__
            raise
        finally:
            manifest.write_text(json.dumps(record,indent=2))
