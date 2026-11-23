"""ASE-independent, fixed-ion periodic ground-state workflow (one MPI rank)."""
from .errors import InputError, InputTypeError, RuntimeFailure, SCFReviewRequired
from dataclasses import dataclass, field, asdict
from pathlib import Path
import hashlib
import json
import warnings
import numpy as np
from .fdtd import SalmonExecutable
from .sbe import _integer
from .ld import _positive
from .tddft import _digest, _copy_bundle



def _review_option(value):
    if type(value) is not bool:
        raise InputTypeError('allow_unconverged must be boolean')


def _review_scf(result, allow_unconverged=False):
    """Retain completed data/owners and pause downstream work for review."""
    _review_option(allow_unconverged)
    if result.converged:
        return
    if allow_unconverged:
        warnings.warn('SCF not converged; explicitly reviewed continuation uses the last iterate',
                      RuntimeWarning, stacklevel=2)
        return
    # Deterministic sibling name also agrees across ranks of a native owner.
    directory = result.directory.parent / (result.directory.name + '-scf-review')
    suffix = 0
    while directory.exists():
        suffix += 1
        directory = result.directory.parent / (result.directory.name + f'-scf-review-{suffix}')
    report = None
    try:
        report = result.write_scf_report(directory)
    except Exception as exc:
        warnings.warn('SCF report unavailable; inspect retained native output: ' + str(exc),
                      RuntimeWarning, stacklevel=2)
    pause = SCFReviewRequired(result, report)
    warnings.warn(str(pause), RuntimeWarning, stacklevel=2)
    raise pause


def _restart_output_text(text,save_restart):
    """Keep fresh-run and preparation output controls byte-identical."""
    if type(save_restart) is not bool:raise InputTypeError('save_restart must be boolean')
    return text if save_restart else text.replace("write_gs_restart_data='all'","write_gs_restart_data='no'")


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
    reference: object=field(default=None,repr=False)

    def __post_init__(self):
        object.__setattr__(self,'atomic_number',_integer(self.atomic_number,'atomic_number',1,118))
        object.__setattr__(self,'local_channel',_integer(self.local_channel,'local_channel',0,4))
        if self.suffix not in ('_rps.dat','.upf'): raise InputError('supported formats: KY _rps.dat and NC UPF')
        if not isinstance(self.data,bytes) or not self.data: raise InputError('nonempty immutable pseudopotential bytes required')
        if self.reference is not None:
            from .pseudopotentials import PseudopotentialReference, _check_upf
            if not isinstance(self.reference, PseudopotentialReference):
                raise InputTypeError('reference must be a PseudopotentialReference')
            if self.reference.atomic_number != self.atomic_number or self.reference.sha256 != self.sha256:
                raise InputError('pseudopotential identity differs from its reference')
            if self.reference.file_format == 'UPF2':
                _check_upf(self.reference, self.data)
        if self.suffix=='.upf':
            import xml.etree.ElementTree as ET
            try:
                header=ET.fromstring(self.data).find('PP_HEADER')
                if header is None or header.get('pseudo_type','').strip().upper()!='NC': raise InputError('NC UPF required')
                for name in ('is_paw','is_ultrasoft'):
                    if header.get(name,'F').strip().lower() not in ('f','false','.false.'): raise InputError('norm-conserving UPF required')
            except ET.ParseError as exc: raise InputError('UPF XML could not be parsed') from exc

    @classmethod
    def from_file(cls,filename,*,atomic_number,local_channel,reference=None):
        path=Path(filename)
        suffix='_rps.dat' if path.name.endswith('_rps.dat') else path.suffix.lower()
        return cls(atomic_number,local_channel,suffix,path.read_bytes(),reference)

    @property
    def has_spin_orbit(self):
        if self.suffix != '.upf':
            return False
        import xml.etree.ElementTree as ET
        header = ET.fromstring(self.data).find('PP_HEADER')
        return header.get('has_so', 'F').strip().lower() in ('t', 'true', '.true.')

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
    mixing_method: str | None=None
    mixing_rate: float | None=0.3
    adaptive_mixing: bool=False
    initial_wavefunction: str | None=None
    random_seed: int | None=None

    @classmethod
    def for_magnetic(cls, *, tolerance=1e-5, **settings):
        """Explicit magnetic-material starting point, not a convergence guarantee.

        Inspect energy and magnetic-moment changes as well as the residual.
        SOC alone does not select this preset. Explicit tolerances take priority.
        """
        return cls(tolerance=tolerance, **settings)

    def __post_init__(self):
        if type(self.adaptive_mixing) is not bool:raise InputError('adaptive_mixing must be bool')
        object.__setattr__(self,'max_iterations',_integer(self.max_iterations,'max_iterations',1,2**31-1))
        object.__setattr__(self,'tolerance',_positive(self.tolerance,'tolerance'))
        if self.mixing_method is not None and self.mixing_method not in ('simple','simple_dm','broyden','pulay','simple_potential'):
            raise InputError('unsupported mixing_method')
        if self.initial_wavefunction is not None and self.initial_wavefunction not in ('random','gauss','gauss2','gauss3','gauss4','gauss5','gauss10','pp'):
            raise InputError('unsupported initial_wavefunction')
        if self.mixing_rate is not None:
            rate=np.asarray(self.mixing_rate)
            if rate.shape!=() or rate.dtype.kind not in 'iuf' or not np.isfinite(rate):
                raise InputError('mixing_rate must be a finite real scalar')
            object.__setattr__(self,'mixing_rate',float(rate))
            if not 0<=float(rate)<=1:
                warnings.warn('mixing_rate outside [0,1]; solver acceptance and stability must be checked',RuntimeWarning,stacklevel=2)
        if self.random_seed is not None:
            object.__setattr__(self,'random_seed',_integer(self.random_seed,'random_seed',0,2**31-1))


_GUESS_FILES=('info.bin','occupation.bin','wfn.bin')


def _initial_guess_source(calculation,initial_guess,destination,executable_sha256=None):
    """Authenticate the existing wavefunction seed; never execute or modify it."""
    if not isinstance(initial_guess,GroundStateResult):
        raise InputTypeError('initial_guess must be a completed GroundStateResult')
    source=initial_guess._with_restart()
    record=json.loads((source.directory/'ground-state.json').read_text())
    signature=calculation._restart_signature()
    if executable_sha256 is not None:signature['executable_sha256']=executable_sha256
    if any(record.get(key)!=value for key,value in signature.items()):
        raise InputError('initial guess requires matching executable, cell, grid, kmesh, species, potentials, electrons, states and XC')
    if Path(destination).expanduser().resolve().is_relative_to(source.directory):
        raise InputError('initial guess output must be outside its source directory')
    return source,record


def _initial_guess_text(text):
    return text.replace("yn_restart='n'", "yn_restart='y', read_gs_restart_data='wfn', "
                        "yn_reset_step_restart='y', directory_read_data='./initial_guess/'",1)


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
            warnings.warn(f'SCF not converged after {iterations} iterations (residual={residual:.6g}, tolerance={tolerance:.6g}); last iterate retained; inspect residual history before continuing', RuntimeWarning, stacklevel=2)
        return dict(iterations=iterations,residual=residual,tolerance=tolerance,energy_hartree=energy,converged=converged)
    except (OSError,ValueError) as exc:
        raise RuntimeFailure('Invalid or missing GS status; inspect gs_status.dat and stdout.log') from exc


def _read_forces(path,natoms):
    rows=Path(path).read_text().splitlines()
    if len(rows)!=natoms+2 or rows[0]!='SALMONPY_FORCES_V1' or int(rows[1])!=natoms:
        raise InputError('invalid native force report')
    return _frozen([[float(x) for x in row.split()] for row in rows[2:]],(natoms,3))


@dataclass(frozen=True,eq=False)
class SymmetryOperations:
    """Explicit fractional rotations/translations, shape (operations,3,4).

    This records user-supplied symmetry, not symmetry discovery or a certificate
    that atoms/grid/k stars preserve it. Native/geometry admission still applies.
    """
    operations: np.ndarray
    directions: str='yyy'

    def __post_init__(self):
        values=np.asarray(self.operations)
        if np.iscomplexobj(values):raise InputError('symmetry operations must be real')
        if values.ndim!=3 or values.shape[0]<1 or values.shape[1:]!=(3,4):
            raise InputError('symmetry operations require shape (n,3,4), n>=1')
        if not isinstance(self.directions,str) or not 1<=len(self.directions)<=3 or any(c not in 'yn' for c in self.directions):
            raise InputError('symmetry directions require one to three y/n characters')
        object.__setattr__(self,'operations',_frozen(values,values.shape))

    @property
    def data(self):
        return (''.join(' '.join(format(float(x),'.17g') for x in row)+'\n'
                        for row in self.operations.reshape(-1,4))).encode('ascii')

    @property
    def sha256(self):return hashlib.sha256(self.data).hexdigest()


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
    symmetry: SymmetryOperations | None=None
    cval: float | None=None
    spin_orbit: bool=False
    electronic_temperature: float | None=None
    electronic_temperature_k: float | None=None

    def __post_init__(self):
        if self.symmetry is not None and not isinstance(self.symmetry,SymmetryOperations):
            raise InputTypeError('symmetry requires explicit SymmetryOperations or None')
        numbers=tuple(_integer(n,'atomic number',1,118) for n in self.numbers)
        if not numbers: raise InputError('at least one atom required')
        pseudos=tuple(self.pseudopotentials)
        if any(not isinstance(p,Pseudopotential) for p in pseudos): raise InputTypeError('expected Pseudopotential objects')
        if len({p.atomic_number for p in pseudos})!=len(pseudos) or {p.atomic_number for p in pseudos}!=set(numbers):
            raise InputError('exactly one pseudopotential per atomic species required')
        electrons=_integer(self.electrons,'electrons',2,2**31-1)
        states=_integer(2*electrons if self.spin_orbit and self.states is None else self.states,'states',1,2**31-1)
        if type(self.spin_orbit) is not bool:
            raise InputTypeError('spin_orbit must be boolean')
        if any(p.has_spin_orbit != self.spin_orbit for p in pseudos):
            raise InputError('spin_orbit must match every pseudopotential; no scalar/SOC substitution')
        if self.spin_orbit:
            if states < electrons:
                raise InputError('SOC requires at least electrons spinor states')
            if states<2*electrons:
                warnings.warn('SOC default is twice the electron count; retaining explicitly supplied states',RuntimeWarning,stacklevel=2)
        elif electrons%2 or states<electrons//2:
            raise InputError('closed-shell even electrons and enough states required')
        cell=_frozen(self.cell,(3,3))
        if np.linalg.det(cell)<=1e-12: raise InputError('positive-volume full periodic cell required')
        positions=_frozen(self.fractional_positions,(len(numbers),3))
        grid=tuple(_integer(n,'grid',8,2**20) for n in self.grid)
        kmesh=tuple(_integer(n,'kmesh',1,2**20) for n in self.kmesh)
        if len(grid)!=3 or len(kmesh)!=3: raise InputError('three grid and kmesh dimensions required')
        if self.xc.lower() not in ('pz','tbmbj','libxc_pz','libxc_pbe','hse06'): raise InputError('supported xc: pz, tbmbj, libxc_pz, libxc_pbe, hse06 (fast native GS)')
        if self.electronic_temperature is not None and self.electronic_temperature_k is not None:
            raise InputError('specify electronic temperature in Hartree or kelvin, not both')
        if self.electronic_temperature_k is not None:
            object.__setattr__(self,'electronic_temperature_k',_positive(self.electronic_temperature_k,'electronic_temperature_k'))
        if self.electronic_temperature is not None:
            object.__setattr__(self,'electronic_temperature',_positive(self.electronic_temperature,'electronic_temperature'))
        cval=self.cval
        if cval is not None:
            value=np.asarray(cval)
            if value.shape!=() or value.dtype.kind not in 'iuf' or not np.isfinite(value):
                raise InputError('cval must be a finite real scalar or None')
            cval=float(value)
            if self.xc.lower()!='tbmbj':
                warnings.warn('cval is only used by TB-mBJ; ignoring it for '+self.xc,RuntimeWarning,stacklevel=2)
                cval=None
            elif cval<=1e-10:
                cval=None  # Native automatic-c sentinel, normalized for restart identity.
        object.__setattr__(self,'cval',cval)
        for name,value in [('numbers',numbers),('pseudopotentials',pseudos),('cell',cell),('fractional_positions',positions),
                           ('electrons',electrons),('states',states),('grid',grid),('kmesh',kmesh),('xc',self.xc.lower())]:
            object.__setattr__(self,name,value)

    def _pseudo_name(self,i): return f'pseudo{i+1}'+self.pseudopotentials[i].suffix

    @property
    def input_assets(self):
        values=tuple((self._pseudo_name(i),p.data) for i,p in enumerate(self.pseudopotentials))
        return values if self.symmetry is None else values+(('sym.dat',self.symmetry.data),)

    def input_text(self,scf=SCFSettings()):
        if not isinstance(scf,SCFSettings): raise InputTypeError('expected SCFSettings')
        vector=lambda x: ','.join(f'{v:.17e}' for v in x)
        lattice='\n'.join(f' al_vec{i+1}='+vector(row)+',' for i,row in enumerate(self.cell))
        if self.symmetry is not None:lattice+=f"\n yn_symmetry='{self.symmetry.directions}',"
        pseudo='\n'.join(f" file_pseudo({i+1})='{self._pseudo_name(i)}', izatom({i+1})={p.atomic_number}, lloc_ps({i+1})={p.local_channel},"
                         for i,p in enumerate(self.pseudopotentials))
        index={p.atomic_number:i+1 for i,p in enumerate(self.pseudopotentials)}
        atoms='\n'.join(f"'Z{z}' "+' '.join(f'{v:.17e}' for v in row)+f' {index[z]}'
                        for z,row in zip(self.numbers,self.fractional_positions))
        hybrid_options=', exx_pre_scf_threshold=0' if self.xc=='hse06' else ''
        if self.cval is not None:hybrid_options+=f', cval={self.cval:.17e}'
        system_options='' if self.electronic_temperature is None else f', temperature={self.electronic_temperature:.17e}'
        if self.electronic_temperature_k is not None:system_options=f', temperature_k={self.electronic_temperature_k:.17e}'
        options=[]
        if scf.adaptive_mixing:options.append("yn_auto_mixing='y', update_mixing_ratio=1.2")
        if scf.mixing_method is not None:options.append(f"method_mixing='{scf.mixing_method}'")
        if scf.mixing_rate is not None:
            coefficient={'broyden':'alpha_mb','pulay':'beta_p'}.get(scf.mixing_method or 'broyden','mixrate')
            options.append(f'{coefficient}={scf.mixing_rate:.17e}')
        if scf.initial_wavefunction is not None:options.append(f"method_init_wf='{scf.initial_wavefunction}'")
        if scf.random_seed is not None:options.append(f'iseed_number_change={scf.random_seed}')
        scf_extra='' if not options else '\n '+', '.join(options)
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
 yn_periodic='y', {"spin='noncollinear', yn_spinorbit='y'" if self.spin_orbit else "spin='unpolarized'"}, nelem={len(index)}, natom={len(self.numbers)}, nelec={self.electrons}, nstate={self.states}{system_options},
{lattice}
/
&pseudo
{pseudo}
/
&functional
 xc='{self.xc}'{hybrid_options}
/
&rgrid
 num_rgrid={','.join(map(str,self.grid))}
/
&kgrid
 num_kgrid={','.join(map(str,self.kmesh))}
/
&scf
 nscf={scf.max_iterations}, threshold={scf.tolerance:.17e}, convergence='rho_dne'{scf_extra}
/
&atomic_red_coor
{atoms}
/
"""

    def _restart_signature(self):
        """Physical/discretization identity; only ionic positions may change."""
        value={'xc': self.xc, 'electrons': self.electrons, 'states': self.states,
                'grid': list(self.grid), 'kmesh': list(self.kmesh),
                'cell_bohr': self.cell.tolist(), 'numbers': list(self.numbers),
                'pseudopotentials': [{'atomic_number': p.atomic_number,
                                     'sha256': p.sha256, 'local_channel': p.local_channel}
                                    for p in self.pseudopotentials]}
        if self.spin_orbit:value['spin_orbit']=True
        if self.electronic_temperature is not None:value['electronic_temperature']=self.electronic_temperature
        if self.electronic_temperature_k is not None:value['electronic_temperature_k']=self.electronic_temperature_k
        if self.xc=='tbmbj':value['cval']=self.cval
        if self.symmetry is not None:value['symmetry']={'directions':self.symmetry.directions,'sha256':self.symmetry.sha256}
        return value

    def run(self,executable,directory,*,scf=SCFSettings(),threads=1,timeout=600.,initial_guess=None,save_restart=True,allow_unconverged=False):
        """Run SCF; terminal nonconvergence retains results and pauses for review.

        allow_unconverged=True requires explicit review, never automatic retry.
        """
        import os,subprocess
        import shutil
        _review_option(allow_unconverged)
        if not isinstance(executable,SalmonExecutable): raise InputTypeError('expected SalmonExecutable')
        if self.xc=='hse06':
            raise InputError('HSE06 GS uses NativeGSBackend with a declared fast-PR compute runtime; batch HSE selection requires fast-PR provenance')
        if type(save_restart) is not bool:raise InputTypeError('save_restart must be boolean')
        if self.xc=='tbmbj':
            warnings.warn('TB-mBJ is a model potential for states/bands/response; total energy and forces '
                          'are unvalidated diagnostics, not an ASE optimization capability.',RuntimeWarning,stacklevel=2)
        text=_restart_output_text(self.input_text(scf),save_restart)
        threads=_integer(threads,'threads',1,2**31-1);timeout=_positive(timeout,'timeout')
        if _digest(executable.path)!=executable.sha256: raise InputError('executable checksum mismatch')
        dest=Path(directory).expanduser().resolve()
        source = None
        if initial_guess is not None:
            source,source_record=_initial_guess_source(self,initial_guess,dest,executable.sha256)
            text=_initial_guess_text(text)
        dest.mkdir(parents=True,exist_ok=False)
        work=dest/'native';work.mkdir()
        record={'schema':'salmonpy.gs/1','status':'started','executable_sha256':executable.sha256,
                'xc':self.xc,'electrons':self.electrons,'states':self.states,'grid':self.grid,'kmesh':self.kmesh,
                'cell_bohr':self.cell.tolist(),'fractional_positions':self.fractional_positions.tolist(),'numbers':self.numbers,
                'pseudopotentials':[{'atomic_number':p.atomic_number,'sha256':p.sha256,'local_channel':p.local_channel} for p in self.pseudopotentials],
                'max_iterations':scf.max_iterations,'tolerance':scf.tolerance,'threads':threads,'restart_saved':save_restart}
        if self.spin_orbit:record['spin_orbit']=True
        if self.electronic_temperature is not None:record['electronic_temperature']=self.electronic_temperature
        if self.electronic_temperature_k is not None:record['electronic_temperature_k']=self.electronic_temperature_k
        if scf.adaptive_mixing or any(getattr(scf,k) is not None for k in ('mixing_method','mixing_rate','initial_wavefunction','random_seed')):
            record['scf_settings']=asdict(scf)
        if any(p.reference is not None for p in self.pseudopotentials):
            record['pseudopotential_sources']=[None if p.reference is None else asdict(p.reference)
                                               for p in self.pseudopotentials]
        if self.xc=='tbmbj':
            record.update(cval=self.cval,cval_mode='automatic' if self.cval is None else 'fixed',
                          observable_scope='states/bands/response; energy/forces are unvalidated diagnostics')
        if self.symmetry is not None:record['symmetry']=self._restart_signature()['symmetry']
        manifest=dest/'ground-state.json'
        try:
            if source is not None:
                guess = work/'initial_guess'
                guess.mkdir()
                for name in _GUESS_FILES:
                    relative = 'native/data_for_restart/'+name
                    shutil.copy2(source.directory/relative, guess/name)
                    if _digest(guess/name) != source_record['files'][relative]:
                        raise InputError('initial guess changed while copying')
                record['initial_guess'] = {'mode': 'wavefunction',
                    'source_directory': str(source.directory),
                    'manifest_sha256': _digest(source.directory/'ground-state.json')}
            (work/'gs.inp').write_text(text)
            for name,data in self.input_assets:(work/name).write_bytes(data)
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
                exists=(work/'data_for_restart'/name).is_file()
                if save_restart and not exists:raise RuntimeFailure('missing restart record: '+name)
            present=[(work/'data_for_restart'/name).is_file() for name in ('info.bin','occupation.bin','wfn.bin')]
            if not save_restart and any(present):
                warnings.warn('runtime wrote restart output despite save_restart=False; calculation results are retained',RuntimeWarning,stacklevel=2)
            record['restart_saved']=all(present)
            record.update(result,status='complete',files={p.relative_to(dest).as_posix():_digest(p) for p in work.rglob('*') if p.is_file()})
        except BaseException as exc:
            record.update(status='failed',error=str(exc));raise
        finally:
            manifest.write_text(json.dumps(record,indent=2))
        result = GroundStateResult.open(dest)
        _review_scf(result, allow_unconverged)
        return result


@dataclass(frozen=True)
class GroundStateResult:
    directory: Path
    energy_hartree: float
    iterations: int
    converged: bool=True
    restart_saved: bool=True

    @classmethod
    def open(cls,directory):
        from pathlib import PureWindowsPath
        root=Path(directory).expanduser().resolve()
        record=json.loads((root/'ground-state.json').read_text())
        if record.get('schema')!='salmonpy.gs/1' or record.get('status')!='complete': raise InputError('not a completed GS result')
        files=record.get('files',{})
        saved=record.get('restart_saved')
        if type(saved) is not bool:raise InputError('invalid restart_saved flag')
        restart={f'native/data_for_restart/{n}' for n in ('info.bin','occupation.bin','wfn.bin')}
        required={'native/gs.inp','native/gs_status.dat'}
        if saved:required.update(restart)
        if not required<=files.keys(): raise InputError('incomplete GS manifest')
        paths=list((root/'native').rglob('*'))
        if any(p.is_symlink() for p in paths) or {p.relative_to(root).as_posix() for p in paths if p.is_file()}!=set(files):
            raise InputError('GS native files differ from the saved manifest')
        for name,digest in files.items():
            path=Path(name)
            if path.is_absolute() or PureWindowsPath(name).drive or '..' in path.parts or '\\' in name or not (root/path).resolve().is_relative_to(root):
                raise InputError('invalid manifest path')
            if _digest(root/path)!=digest: raise InputError('GS checksum mismatch: '+name)
        if record.get('xc')=='tbmbj':
            from .bands import _input_groups,_scalar
            functional=_input_groups((root/'native/gs.inp').read_text()).get('functional','')
            try:value=float(_scalar(functional,'cval','-1').replace('d','e'))
            except ValueError as exc:raise InputError('invalid TB-mBJ cval input') from exc
            expected=None if value<=1e-10 else value
            if (_scalar(functional,'xc')!='tbmbj' or not np.isfinite(value)
                    or record.get('cval')!=expected):
                raise InputError('TB-mBJ cval/XC metadata differs from GS input')
        status=_read_status(root/'native/gs_status.dat')
        if status['energy_hartree']!=record['energy_hartree'] or status['iterations']!=record['iterations']:
            raise InputError('inconsistent GS metadata')
        if 'converged' in record and record['converged'] != status['converged']:
            raise InputError('inconsistent GS convergence metadata')
        return cls(root,status['energy_hartree'],status['iterations'],status['converged'],saved)

    def write_scf_report(self, directory):
        """Write actual native residual history and an optional plot, without SCF."""
        import csv
        verified = self.open(self.directory)
        status = _read_status(verified.directory/'native/gs_status.dat')
        rows = []
        log = verified.directory/'native/stdout.log'
        if log.is_file():
            for line in log.read_text().splitlines():
                if line.strip().startswith('iter and int_x|rho_i(x)-rho_i-1(x)|dx/nelec'):
                    fields = line.rsplit('=', 1)[-1].split()
                    if len(fields) == 2:
                        iteration, residual = int(fields[0]), float(fields[1].replace('D','e').replace('d','e'))
                        if iteration > 0 and np.isfinite(residual) and residual >= 0:
                            rows.append((iteration, residual, status['tolerance']))
        source = 'native-stdout-rho_dne' if rows else 'final-status-only'
        if not rows:
            rows = [(status['iterations'], status['residual'], status['tolerance'])]
        directory = Path(directory).expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=False)
        with (directory/'residuals.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['iteration','residual','threshold'])
            writer.writerows(rows)
        report = dict(status, metric='rho_dne', history_source=source, plot_status='unavailable')
        figure = None
        try:
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            figure = Figure(figsize=(7,4)); FigureCanvasAgg(figure)
            axis = figure.subplots()
            axis.semilogy([r[0] for r in rows], [r[1] for r in rows], marker='o', label='Residual')
            axis.axhline(status['tolerance'], linestyle='--', label='Threshold')
            axis.set(xlabel='SCF iteration', ylabel='rho_dne',
                     title='SCF converged' if status['converged'] else 'SCF not converged — review required')
            axis.grid(True, which='both', alpha=.25); axis.legend(); figure.tight_layout()
            figure.savefig(directory/'residuals.png', dpi=150)
            report['plot_status'] = 'saved'
        except Exception as exc:
            report['plot_error'] = str(exc)
            warnings.warn('SCF residual CSV saved; plot unavailable: ' + str(exc), UserWarning, stacklevel=2)
        finally:
            if figure is not None: figure.clear()
        (directory/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        return directory

    def _with_restart(self):
        verified=self.open(self.directory)
        if not verified.restart_saved:
            raise InputError('restart arrays were not saved; run a new GS with save_restart=True for this operation')
        return verified

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
        return converter.export_archive(self._with_restart(),directory,orbital_block=orbital_block,timeout=timeout)

    def band_initial_state(self,*,executable,allow_unconverged=False):
        from .bands import BandInitialState
        verified=self._with_restart()
        _review_scf(verified, allow_unconverged)
        return BandInitialState.from_salmon(verified.directory/'native','gs.inp',executable=executable)

    def tddft_initial_state(self,*,executable,max_steps=1000,threads=1,timeout=300.,allow_unconverged=False):
        import tempfile
        from .bands import _input_groups
        from .tddft import TDDFTInitialState
        verified=self._with_restart()
        _review_scf(verified, allow_unconverged)
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
