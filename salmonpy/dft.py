"""ASE-independent, fixed-ion periodic ground-state workflow (one MPI rank)."""
from .errors import InputError, InputTypeError, RuntimeFailure
from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import json
import warnings
import numpy as np
from .fdtd import SalmonExecutable
from .sbe import _integer
from .ld import _positive
from .tddft import _digest, _copy_bundle


def _frozen(value,shape):
    a=np.array(value,dtype=float,copy=True)
    if a.shape!=shape or not np.isfinite(a).all(): raise InputError(f'expected finite array of shape {shape}')
    return np.frombuffer(a.tobytes(),dtype=a.dtype).reshape(shape)


@dataclass(frozen=True)
class Pseudopotential:
    atomic_number: int
    local_channel: int
    suffix: str
    data: bytes=field(repr=False)

    def __post_init__(self):
        object.__setattr__(self,'atomic_number',_integer(self.atomic_number,'atomic_number',1,118))
        object.__setattr__(self,'local_channel',_integer(self.local_channel,'local_channel',0,4))
        if self.suffix not in ('_rps.dat','.upf'): raise InputError('supported formats: KY _rps.dat and NC UPF')
        if not isinstance(self.data,bytes) or not self.data: raise InputError('nonempty immutable pseudopotential bytes required')
        if self.suffix=='.upf':
            import xml.etree.ElementTree as ET
            try:
                header=ET.fromstring(self.data).find('PP_HEADER')
                if header is None or header.get('pseudo_type','').strip().upper()!='NC': raise InputError('NC UPF required')
                for name in ('is_paw','is_ultrasoft','has_so'):
                    if header.get(name,'F').strip().lower() not in ('f','false','.false.'): raise InputError('scalar NC UPF required')
            except ET.ParseError as exc: raise InputError('UPF XML could not be parsed') from exc

    @classmethod
    def from_file(cls,filename,*,atomic_number,local_channel):
        path=Path(filename)
        suffix='_rps.dat' if path.name.endswith('_rps.dat') else path.suffix.lower()
        return cls(atomic_number,local_channel,suffix,path.read_bytes())

    @property
    def sha256(self): return hashlib.sha256(self.data).hexdigest()

    @property
    def valence_electrons(self):
        """Integer ionic charge from KY or NC-UPF metadata, never atomic Z.

        This checks the charge header only, not the entire potential or its XC
        provenance. Fractional charges are unsupported by the neutral ASE path.
        """
        try:
            if self.suffix == '_rps.dat':
                # SALMON Read_PS_KY: mr, step, mlps, rzps.
                header = self.data.decode('utf-8').splitlines()[0]
                value = header.split('!', 1)[0].replace(',', ' ').split()[3]
            else:
                import xml.etree.ElementTree as ET
                value = ET.fromstring(self.data).find('PP_HEADER').get('z_valence')
            charge = float(value.replace('D', 'e').replace('d', 'e'))
            if not np.isfinite(charge) or not charge.is_integer() or not 1 <= charge <= self.atomic_number:
                raise InputError('unsupported ionic charge')
        except (ValueError, IndexError, AttributeError, UnicodeError) as exc:
            raise InputError('pseudopotential requires a valid integer valence charge') from exc
        return int(charge)


@dataclass(frozen=True)
class SCFSettings:
    max_iterations: int=300
    tolerance: float=1e-9

    def __post_init__(self):
        object.__setattr__(self,'max_iterations',_integer(self.max_iterations,'max_iterations',1,2**31-1))
        object.__setattr__(self,'tolerance',_positive(self.tolerance,'tolerance'))


def _read_status(path):
    try:
        rows=Path(path).read_text().splitlines()
        if len(rows)!=4 or rows[0]!='SALMONPY_GS_V1': raise InputError('unsupported status report')
        flag,iterations=map(int,rows[1].split());residual,tolerance,energy=map(float,rows[2].split())
        if flag not in (0,1) or iterations<1 or rows[3]!='rho_dne' or not np.isfinite([residual,tolerance,energy]).all():
            raise InputError('invalid status')
        if residual<0 or tolerance<=0: raise InputError('invalid residual or tolerance')
        converged=bool(flag==1 and residual<tolerance)
        if not converged:
            warnings.warn(f'SCF not converged after {iterations} iterations (residual={residual:.6g}, tolerance={tolerance:.6g}); continuing with the saved state', RuntimeWarning, stacklevel=2)
        return dict(iterations=iterations,residual=residual,tolerance=tolerance,energy_hartree=energy,converged=converged)
    except (OSError,ValueError) as exc:
        raise RuntimeFailure('Invalid or missing GS status; inspect gs_status.dat and stdout.log') from exc


def _read_forces(path,natoms):
    rows=Path(path).read_text().splitlines()
    if len(rows)!=natoms+2 or rows[0]!='SALMONPY_FORCES_V1' or int(rows[1])!=natoms:
        raise InputError('invalid native force report')
    return _frozen([[float(x) for x in row.split()] for row in rows[2:]],(natoms,3))


@dataclass(frozen=True,eq=False)
class DFT:
    """Geometry (cell rows in bohr, fractional positions) and fixed model settings.

    Electron count and scalar pseudopotential/XC consistency are explicit caller
    responsibilities. Native defaults govern the solver beyond SCF stopping criteria.
    """
    numbers: tuple
    cell: np.ndarray
    fractional_positions: np.ndarray
    pseudopotentials: tuple
    electrons: int
    states: int
    grid: tuple
    kmesh: tuple
    xc: str='pz'

    def __post_init__(self):
        numbers=tuple(_integer(n,'atomic number',1,118) for n in self.numbers)
        if not numbers: raise InputError('at least one atom required')
        pseudos=tuple(self.pseudopotentials)
        if any(not isinstance(p,Pseudopotential) for p in pseudos): raise InputTypeError('expected Pseudopotential objects')
        if len({p.atomic_number for p in pseudos})!=len(pseudos) or {p.atomic_number for p in pseudos}!=set(numbers):
            raise InputError('exactly one pseudopotential per atomic species required')
        electrons=_integer(self.electrons,'electrons',2,2**31-1)
        states=_integer(self.states,'states',1,2**31-1)
        if electrons%2 or states<electrons//2: raise InputError('closed-shell even electrons and enough states required')
        cell=_frozen(self.cell,(3,3))
        if np.linalg.det(cell)<=1e-12: raise InputError('positive-volume full periodic cell required')
        positions=_frozen(self.fractional_positions,(len(numbers),3))
        grid=tuple(_integer(n,'grid',8,2**20) for n in self.grid)
        kmesh=tuple(_integer(n,'kmesh',1,2**20) for n in self.kmesh)
        if len(grid)!=3 or len(kmesh)!=3: raise InputError('three grid and kmesh dimensions required')
        if self.xc.lower() not in ('pz','libxc_pz','libxc_pbe'): raise InputError('supported xc: pz, libxc_pz, libxc_pbe')
        for name,value in [('numbers',numbers),('pseudopotentials',pseudos),('cell',cell),('fractional_positions',positions),
                           ('electrons',electrons),('states',states),('grid',grid),('kmesh',kmesh),('xc',self.xc.lower())]:
            object.__setattr__(self,name,value)

    def _pseudo_name(self,i): return f'pseudo{i+1}'+self.pseudopotentials[i].suffix

    def input_text(self,scf=SCFSettings()):
        if not isinstance(scf,SCFSettings): raise InputTypeError('expected SCFSettings')
        vector=lambda x: ','.join(f'{v:.17e}' for v in x)
        lattice='\n'.join(f' al_vec{i+1}='+vector(row)+',' for i,row in enumerate(self.cell))
        pseudo='\n'.join(f" file_pseudo({i+1})='{self._pseudo_name(i)}', izatom({i+1})={p.atomic_number}, lloc_ps({i+1})={p.local_channel},"
                         for i,p in enumerate(self.pseudopotentials))
        index={p.atomic_number:i+1 for i,p in enumerate(self.pseudopotentials)}
        atoms='\n'.join(f"'Z{z}' "+' '.join(f'{v:.17e}' for v in row)+f' {index[z]}'
                        for z,row in zip(self.numbers,self.fractional_positions))
        return f"""&calculation
 theory='dft'
/
&control
 sysname='gs', base_directory='./', yn_restart='n', write_gs_restart_data='all', yn_self_checkpoint='n'
/
&units
 unit_system='au'
/
&system
 yn_periodic='y', spin='unpolarized', nelem={len(index)}, natom={len(self.numbers)}, nelec={self.electrons}, nstate={self.states},
{lattice}
/
&pseudo
{pseudo}
/
&functional
 xc='{self.xc}'
/
&rgrid
 num_rgrid={','.join(map(str,self.grid))}
/
&kgrid
 num_kgrid={','.join(map(str,self.kmesh))}
/
&scf
 nscf={scf.max_iterations}, threshold={scf.tolerance:.17e}, convergence='rho_dne'
/
&atomic_red_coor
{atoms}
/
"""

    def _restart_signature(self):
        """Physical/discretization identity; only ionic positions may change."""
        return {'xc': self.xc, 'electrons': self.electrons, 'states': self.states,
                'grid': list(self.grid), 'kmesh': list(self.kmesh),
                'cell_bohr': self.cell.tolist(), 'numbers': list(self.numbers),
                'pseudopotentials': [{'atomic_number': p.atomic_number,
                                     'sha256': p.sha256, 'local_channel': p.local_channel}
                                    for p in self.pseudopotentials]}

    def run(self,executable,directory,*,scf=SCFSettings(),threads=1,timeout=600.,initial_guess=None):
        import os,subprocess
        import shutil
        if not isinstance(executable,SalmonExecutable): raise InputTypeError('expected SalmonExecutable')
        text=self.input_text(scf)
        threads=_integer(threads,'threads',1,2**31-1);timeout=_positive(timeout,'timeout')
        if _digest(executable.path)!=executable.sha256: raise InputError('executable checksum mismatch')
        dest=Path(directory).expanduser().resolve()
        source = None
        if initial_guess is not None:
            if not isinstance(initial_guess, GroundStateResult):
                raise InputTypeError('initial_guess must be a completed GroundStateResult')
            source = GroundStateResult.open(initial_guess.directory)
            source_record = json.loads((source.directory/'ground-state.json').read_text())
            signature = dict(self._restart_signature(), executable_sha256=executable.sha256)
            if any(source_record.get(key) != value for key, value in signature.items()):
                raise InputError('initial guess requires matching executable, cell, grid, kmesh, species, potentials, electrons, states and XC')
            if dest.is_relative_to(source.directory):
                raise InputError('initial guess output must be outside its source directory')
            text = text.replace("yn_restart='n'", "yn_restart='y', read_gs_restart_data='wfn', "
                                "yn_reset_step_restart='y', directory_read_data='./initial_guess/'", 1)
        dest.mkdir(parents=True,exist_ok=False)
        work=dest/'native';work.mkdir()
        record={'schema':'salmonpy.gs/1','status':'started','executable_sha256':executable.sha256,
                'xc':self.xc,'electrons':self.electrons,'states':self.states,'grid':self.grid,'kmesh':self.kmesh,
                'cell_bohr':self.cell.tolist(),'fractional_positions':self.fractional_positions.tolist(),'numbers':self.numbers,
                'pseudopotentials':[{'atomic_number':p.atomic_number,'sha256':p.sha256,'local_channel':p.local_channel} for p in self.pseudopotentials],
                'max_iterations':scf.max_iterations,'tolerance':scf.tolerance,'threads':threads}
        manifest=dest/'ground-state.json'
        try:
            if source is not None:
                guess = work/'initial_guess'
                guess.mkdir()
                for name in ('info.bin', 'occupation.bin', 'wfn.bin'):
                    relative = 'native/data_for_restart/'+name
                    shutil.copy2(source.directory/relative, guess/name)
                    if _digest(guess/name) != source_record['files'][relative]:
                        raise InputError('initial guess changed while copying')
                record['initial_guess'] = {'mode': 'wavefunction',
                    'source_directory': str(source.directory),
                    'manifest_sha256': _digest(source.directory/'ground-state.json')}
            (work/'gs.inp').write_text(text)
            for i,p in enumerate(self.pseudopotentials): (work/self._pseudo_name(i)).write_bytes(p.data)
            record['input_sha256']=_digest(work/'gs.inp');manifest.write_text(json.dumps(record,indent=2))
            env={k:v for k,v in os.environ.items() if k!='SALMONPY_TDDFT_WORKER'}
            env.update(OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1')
            with (work/'stdout.log').open('w') as log:
                proc=subprocess.run([str(executable.path)],input=text,text=True,cwd=work,stdout=log,
                                    stderr=subprocess.STDOUT,env=env,timeout=timeout,check=False)
            if proc.returncode or 'end SALMON' not in (work/'stdout.log').read_text(): raise RuntimeFailure('SALMON did not finish; inspect native/stdout.log')
            result=_read_status(work/'gs_status.dat')
            _read_forces(work/'gs_forces.dat',len(self.numbers))
            if result['tolerance']!=scf.tolerance: raise RuntimeFailure('native convergence threshold differs from request')
            for name in ('info.bin','occupation.bin','wfn.bin'):
                if not (work/'data_for_restart'/name).is_file(): raise RuntimeFailure('missing restart record: '+name)
            record.update(result,status='complete',files={p.relative_to(dest).as_posix():_digest(p) for p in work.rglob('*') if p.is_file()})
        except BaseException as exc:
            record.update(status='failed',error=str(exc));raise
        finally:
            manifest.write_text(json.dumps(record,indent=2))
        return GroundStateResult.open(dest)


@dataclass(frozen=True)
class GroundStateResult:
    directory: Path
    energy_hartree: float
    iterations: int
    converged: bool=True

    @classmethod
    def open(cls,directory):
        from pathlib import PureWindowsPath
        root=Path(directory).expanduser().resolve()
        record=json.loads((root/'ground-state.json').read_text())
        if record.get('schema')!='salmonpy.gs/1' or record.get('status')!='complete': raise InputError('not a completed GS result')
        files=record.get('files',{})
        required={'native/gs.inp','native/gs_status.dat',*(f'native/data_for_restart/{n}' for n in ('info.bin','occupation.bin','wfn.bin'))}
        if not required<=files.keys(): raise InputError('incomplete GS manifest')
        paths=list((root/'native').rglob('*'))
        if any(p.is_symlink() for p in paths) or {p.relative_to(root).as_posix() for p in paths if p.is_file()}!=set(files):
            raise InputError('GS native files differ from the saved manifest')
        for name,digest in files.items():
            path=Path(name)
            if path.is_absolute() or PureWindowsPath(name).drive or '..' in path.parts or '\\' in name or not (root/path).resolve().is_relative_to(root):
                raise InputError('invalid manifest path')
            if _digest(root/path)!=digest: raise InputError('GS checksum mismatch: '+name)
        status=_read_status(root/'native/gs_status.dat')
        if status['energy_hartree']!=record['energy_hartree'] or status['iterations']!=record['iterations']:
            raise InputError('inconsistent GS metadata')
        if 'converged' in record and record['converged'] != status['converged']:
            raise InputError('inconsistent GS convergence metadata')
        return cls(root,status['energy_hartree'],status['iterations'],status['converged'])

    @property
    def forces_hartree_per_bohr(self):
        verified=self.open(self.directory)
        record=json.loads((verified.directory/'ground-state.json').read_text())
        return _read_forces(verified.directory/'native/gs_forces.dat',len(record['numbers']))

    def export_archive(self,directory,*,converter,orbital_block=16,timeout=300.):
        """Export verified serial native GS in converter scope; density is reconstructed.

        Requires a matching explicit compiler-paired exporter and AU eigenvalue
        output. Does not run SCF or preserve density-mixing/RT histories.
        """
        from .gs_native import GammaRestartConverter
        if not isinstance(converter,GammaRestartConverter):raise InputTypeError('GammaRestartConverter required')
        return converter.export_archive(self,directory,orbital_block=orbital_block,timeout=timeout)

    def band_initial_state(self,*,executable):
        from .bands import BandInitialState
        verified=self.open(self.directory)
        return BandInitialState.from_salmon(verified.directory/'native','gs.inp',executable=executable)

    def tddft_initial_state(self,*,executable,max_steps=1000,threads=1,timeout=300.):
        import tempfile
        from .bands import _input_groups
        from .tddft import TDDFTInitialState
        verified=self.open(self.directory)
        with tempfile.TemporaryDirectory(prefix='salmonpy-gs-to-rt-') as tmp:
            bundle=Path(tmp)/'bundle';_copy_bundle(verified.directory/'native',bundle)
            groups=_input_groups((bundle/'gs.inp').read_text())
            groups['calculation']="\n theory='tddft_pulse'\n"
            groups['control']="\n sysname='rt', directory_read_data='./data_for_restart/', yn_restart='n'\n"
            groups.pop('scf',None)
            groups['tgrid']='\n dt=0.02, nt=1\n'
            groups['emfield']="\n ae_shape1='none', trans_longi='tr'\n"
            groups['propagation']="\n propagator='middlepoint', yn_predictor_corrector='n'\n"
            (bundle/'rt.inp').write_text('\n'.join('&'+key+body+'/' for key,body in groups.items())+'\n')
            return TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=executable,max_steps=max_steps,threads=threads,timeout=timeout)
