"""Explicit compiler-bound serial GS/2→native TDDFT bundle bridge.

No compilation, downloads or solver launch during conversion. An independently
built, audited and probed converter must match the caller's runtime hash.
"""
from dataclasses import dataclass,replace
import hashlib,json,os,shutil,subprocess,tempfile
from pathlib import Path
import numpy as np
from .errors import InputError,InputTypeError,RuntimeFailure
from .fdtd import SalmonExecutable

__all__=['GammaRestartConverter','SerialRestartConverter']
RECEIPT='salmonpy-native-gs.json'
ENVIRONMENT=('GFORTRAN_CONVERT_UNIT',)


def _digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def _environment(environ=None):
    value=(os.environ if environ is None else environ).get('GFORTRAN_CONVERT_UNIT','')
    if value.lower() not in ('','native','swap','big_endian','little_endian'):
        raise InputError('unit-specific native record conversion settings are unsupported')
    return (('GFORTRAN_CONVERT_UNIT',value),)


def _receipt_environment(root,engine):
    """Called on the copied immutable initial bundle, never loads executable paths."""
    p=root/RECEIPT
    if not p.exists():return None
    try:
        r=json.loads(p.read_text())
        if r['schema'] not in ('salmonpy.gs-gamma-bundle/1','salmonpy.gs-serial-bundle/1') or r['target_sha256']!=engine.sha256:raise InputError('native GS receipt requires its paired runtime')
        env=r['record_environment']
        if type(env) is not dict or set(env)!=set(ENVIRONMENT) or any(type(v) is not str for v in env.values()):raise InputError('invalid native GS record environment')
        expected=r['native_files']
        if type(expected) is not dict or set(expected)!={'info.bin','occupation.bin','wfn.bin'}:
            raise InputError('invalid native GS file receipt')
        if any(_digest(root/'data_for_restart'/name)!=value for name,value in expected.items()):
            raise InputError('native GS files differ from conversion receipt')
        if env['GFORTRAN_CONVERT_UNIT'].lower() not in ('','native','swap','big_endian','little_endian'):
            raise InputError('unit-specific record settings are unsupported')
        if r['schema']=='salmonpy.gs-serial-bundle/1':
            if r.get('kpoints_file')!='kpoints.dat' or _digest(root/'kpoints.dat')!=r.get('kpoints_sha256'):
                raise InputError('native k-point file differs from conversion receipt')
        return tuple(sorted(env.items()))
    except (KeyError,TypeError,ValueError) as exc:
        if isinstance(exc,InputError):raise
        raise InputError('invalid native GS conversion receipt') from exc


@dataclass(frozen=True)
class GammaRestartConverter:
    """Loaded from a trusted explicit build record; hashes rechecked before use."""
    path:Path
    sha256:str
    target_sha256:str
    build_record:Path
    exporter:Path | None=None
    exporter_sha256:str | None=None
    _multi_k=False

    @classmethod
    def from_build(cls,record):
        p=Path(record).expanduser().resolve();r=json.loads(p.read_text())
        if (r.get('schema')!='salmonpy.gs-native-build/1' or r.get('status')!='built-audited-and-probed'
                or r.get('native_byte_parity') is not True):raise InputError('audited native-reader/writer-probed converter build required')
        converter=Path(r['converter']).resolve()
        if not converter.is_file() or _digest(converter)!=r['sha256']:raise InputError('converter differs from audited executable')
        engine=Path(r['target_executable']).resolve()
        if not engine.is_file() or _digest(engine)!=r['target_sha256']:raise InputError('paired runtime differs from build record')
        exporter=None;exporter_sha=None
        if 'exporter' in r:
            exporter=Path(r['exporter']).resolve();exporter_sha=r.get('exporter_sha256')
            if r.get('native_export_parity') is not True or not exporter.is_file() or _digest(exporter)!=exporter_sha:
                raise InputError('native exporter differs from audited/probed build')
        return cls(converter,r['sha256'],r['target_sha256'],p,exporter,exporter_sha)

    def write_bundle(self,archive,directory,*,executable,pseudopotentials,timeout=300.):
        """Create a new self-contained RT bundle; does not start TDDFT.

        Γ only, one serial CPU rank, fixed ions/unpolarized; grid>=8, closed-shell
        even electrons and native DFT XC scope. Density/eigenvalue arrays are not
        native restart records: SALMON reconstructs density/potentials from u/occ.
        SCF convergence/norm/orthogonality remain caller responsibilities.
        """
        from .gs_archive import GroundStateArchive
        from .dft import DFT,Pseudopotential
        from .ld import _positive
        if not isinstance(archive,GroundStateArchive):raise InputTypeError('GroundStateArchive required')
        if not isinstance(executable,SalmonExecutable):raise InputTypeError('SalmonExecutable required')
        timeout=_positive(timeout,'timeout')
        if executable.mpi_enabled or executable.sha256!=self.target_sha256 or _digest(executable.path)!=self.target_sha256:raise InputError('converter requires its paired serial runtime')
        if _digest(self.path)!=self.sha256:raise InputError('converter hash changed')
        # Reverify all source bytes, including non-wave metadata, before producing a bundle.
        archive=GroundStateArchive.open(archive.directory);m=archive.metadata
        nk=archive.wavefunction_shape[0];kpoints=archive.array('kpoints');weights=archive.array('weights')
        if not self._multi_k:
            if nk!=1 or not np.array_equal(kpoints,np.zeros((1,3))) or not np.array_equal(weights,np.ones(1)):
                raise InputError('native bridge supports exact Gamma with weight one only')
        else:
            cell=archive.array('cell_bohr')
            if not np.array_equal(cell,np.diag(np.diag(cell))) or np.any(np.diag(cell)<=0):
                raise InputError('serial multi-k bridge requires aligned positive orthogonal cell')
        n=m['electrons']
        if not float(n).is_integer():raise InputError('native bridge requires integer closed-shell electrons')
        potentials=tuple(pseudopotentials)
        if any(not isinstance(p,Pseudopotential) for p in potentials):raise InputTypeError('explicit Pseudopotential objects required')
        refs={p['atomic_number']:p['sha256'] for p in m['pseudopotentials']}
        if len(potentials)!=len(refs) or {p.atomic_number:p.sha256 for p in potentials}!=refs:raise InputError('supplied pseudopotentials do not match archive species/hashes')
        if 'provenance' in m:
            channels={c['atomic_number']:c['local_channel'] for c in m['provenance']['local_channels']}
            if {p.atomic_number:p.local_channel for p in potentials}!=channels:raise InputError('local channels differ from source GS')
        model=DFT(tuple(map(int,archive.array('atomic_numbers'))),archive.array('cell_bohr'),archive.array('fractional_positions'),
                  potentials,int(n),m['n_bands'],tuple(reversed(m['grid_zyx'])),(1,1,1),m['xc'])
        from .bands import _input_groups
        groups=_input_groups(model.input_text());groups['calculation']="\n theory='tddft_pulse'\n"
        groups['control']="\n sysname='rt', base_directory='./', directory_read_data='./data_for_restart/', yn_restart='n', yn_self_checkpoint='n'\n"
        if self._multi_k:groups['kgrid']="\n file_kw='kpoints.dat'\n"
        groups.pop('scf',None);groups['tgrid']='\n dt=0.02, nt=1\n'
        groups['emfield']="\n ae_shape1='none', trans_longi='tr'\n"
        groups['propagation']="\n propagator='middlepoint', yn_predictor_corrector='n'\n"
        text='\n'.join('&'+key+body+'/' for key,body in groups.items())+'\n'
        root=Path(directory).expanduser()
        if root.resolve().is_relative_to(archive.directory):raise InputError('native bundle must be outside portable archive')
        process_environment=dict(os.environ);environment=_environment(process_environment)
        root.mkdir(parents=True,exist_ok=False)
        try:
            (root/'data_for_restart').mkdir()
            nx,ny,nz=model.grid;nb=model.states
            (root/'request.txt').write_text(f'{nx} {ny} {nz} {nb}\n{int(n)}\n'+(f'{nk}\n' if self._multi_k else ''))
            (root/'portable_occupations.bin').write_bytes(np.asarray(archive.array('occupations'),dtype='<f8').tobytes())
            if self._multi_k:
                (root/'portable_weights.bin').write_bytes(np.asarray(weights,dtype='<f8').tobytes())
                lines=[str(nk)]+[str(k+1)+' '+' '.join(format(float(x),'.17e') for x in (*kpoints[k],weights[k])) for k in range(nk)]
                (root/'kpoints.dat').write_text('\n'.join(lines)+'\n')
            block=m['orbital_block']
            with (root/'portable_real.bin').open('wb') as real,(root/'portable_imag.bin').open('wb') as imag:
                for k in range(nk):
                    for lo in range(0,nb,block):
                        waves=archive.wavefunctions(k,lo,min(lo+block,nb))
                        real.write(np.asarray(waves.real,dtype='<f8').tobytes(order='C'))
                        imag.write(np.asarray(waves.imag,dtype='<f8').tobytes(order='C'))
            del waves
            with (root/'converter.log').open('w') as log:
                result=subprocess.run([str(self.path)],cwd=root,stdout=log,stderr=subprocess.STDOUT,timeout=timeout,check=False,env=process_environment)
            if result.returncode:raise RuntimeFailure('native GS conversion failed: '+(root/'converter.log').read_text()[-4000:])
            files={name:_digest(root/'data_for_restart'/name) for name in ('info.bin','occupation.bin','wfn.bin')}
            if (root/'data_for_restart/wfn.bin').stat().st_size!=16*nx*ny*nz*nb*nk:raise RuntimeFailure('converted complex wavefunction size mismatch')
            if self._multi_k:(root/'portable_weights.bin').unlink()
            for name in ('request.txt','portable_occupations.bin','portable_real.bin','portable_imag.bin'):(root/name).unlink()
            (root/'rt.inp').write_text(text)
            for i,p in enumerate(potentials):(root/model._pseudo_name(i)).write_bytes(p.data)
            receipt={'schema':'salmonpy.gs-gamma-bundle/1','target_sha256':self.target_sha256,'converter_sha256':self.sha256,
                     'portable_manifest_sha256':_digest(archive.directory/'manifest.json'),'native_files':files,
                     'record_environment':dict(environment),'scope':'Gamma, serial CPU, periodic unpolarized fixed ions; not SCF certification'}
            if self._multi_k:receipt.update(schema='salmonpy.gs-serial-bundle/1',kpoints_file='kpoints.dat',kpoints_sha256=_digest(root/'kpoints.dat'),scope='aligned orthogonal cell, serial CPU, periodic unpolarized fixed ions; not SCF certification')
            (root/RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n')
            return root.resolve()
        except BaseException:
            shutil.rmtree(root);raise

    def export_archive(self,result,directory,*,orbital_block=16,timeout=300.):
        """Export verified native GS/1 in this converter scope; reconstruct density."""
        from ._gs_archive import export_gamma
        return export_gamma(self,result,directory,orbital_block=orbital_block,timeout=timeout)

    def tddft_initial_state(self,archive,*,executable,pseudopotentials,max_steps=1000,threads=1,timeout=300.):
        """Convert into owned immutable TDDFT baseline; no evolving solver starts yet."""
        from .tddft import TDDFTInitialState
        with tempfile.TemporaryDirectory(prefix='salmonpy-gamma-to-rt-') as tmp:
            bundle=self.write_bundle(archive,Path(tmp)/'bundle',executable=executable,pseudopotentials=pseudopotentials,timeout=timeout)
            initial=TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=executable,max_steps=max_steps,threads=threads,timeout=timeout)
            return replace(initial,_source_root=archive.directory)


class SerialRestartConverter(GammaRestartConverter):
    """Opt-in aligned orthogonal multi-k bridge; Gamma API remains restricted."""
    _multi_k=True

    @classmethod
    def from_build(cls,record):
        converter=super().from_build(record)
        if json.loads(converter.build_record.read_text()).get('multi_k_byte_parity') is not True:
            raise InputError('explicitly build a converter with probed multi-k capability')
        return converter
