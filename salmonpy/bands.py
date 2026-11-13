"""Read native dft_band output and plot explicitly referenced band energies."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np

HARTREE_EV = 27.211386245988


@dataclass(frozen=True)
class BandStructure:
    """Cartesian k in bohr^-1; energies shaped (k, band, spin) in Hartree.

    No automatic Fermi level or symmetry labels are inferred. Native output may
    contain repeated endpoint padding. Convergence must be checked in the run log.
    """
    k_cartesian: np.ndarray
    energies: np.ndarray

    def __post_init__(self):
        k = np.array(self.k_cartesian, dtype=float, copy=True)
        e = np.array(self.energies, dtype=float, copy=True)
        if k.ndim != 2 or k.shape[1] != 3 or len(k) < 2:
            raise ValueError('k_cartesian must have shape (at least 2, 3)')
        if e.ndim != 3 or e.shape[0] != len(k) or min(e.shape[1:]) < 1:
            raise ValueError('energies must have shape (k, band, spin)')
        if not np.isfinite(k).all() or not np.isfinite(e).all():
            raise ValueError('band data must be finite')
        for name, arr in [('k_cartesian', k), ('energies', e)]:
            frozen = np.frombuffer(arr.tobytes(), dtype=arr.dtype).reshape(arr.shape)
            object.__setattr__(self, name, frozen)

    @classmethod
    def from_salmon(cls, filename):
        """Read block-format band.dat, including spin channels and actual k.

        Supports legacy and corrected writers. Actual Cartesian columns are used,
        never the legacy first coordinate columns (which can be transformed twice).
        """
        lines = iter(Path(filename).read_text().splitlines())
        try:
            dims = []
            for key in ('Number_of_Bands:', 'Number_of_kpt_in_each_block:', 'Number_of_blocks:'):
                row = next(lines).split()
                if len(row) != 2 or row[0] != key:
                    raise ValueError('invalid band.dat header')
                dims.append(int(row[1]))
            nb, nk, blocks = dims
            if min(dims) < 1:
                raise ValueError('band dimensions must be positive')
            coords, energies = [], []
            for block in range(blocks):
                for ik in range(nk):
                    row = next(lines).split()
                    if len(row) != 7 or int(row[0]) != block*nk+ik+1:
                        raise ValueError('invalid k-point record')
                    coords.append([float(x.replace('D', 'E')) for x in row[4:]])
                for ik in range(nk):
                    bands = []
                    for ib in range(nb):
                        row = next(lines).split()
                        if len(row) not in (3, 4) or [int(x) for x in row[:2]] != [ik+1, ib+1]:
                            raise ValueError('invalid eigenvalue record')
                        bands.append([float(x.replace('D', 'E')) for x in row[2:]])
                    energies.append(bands)
            if any(line.strip() for line in lines):
                raise ValueError('unexpected trailing band data')
            return cls(coords, energies)
        except StopIteration as exc:
            raise ValueError('truncated band.dat') from exc

    @property
    def distance(self):
        return np.r_[0., np.cumsum(np.linalg.norm(np.diff(self.k_cartesian, axis=0), axis=1))]

    def plot(self, *, reference_ev, ticks=(), ax=None, title=None):
        """Plot using caller-supplied energy zero and (k-index, label) ticks.

        reference_ev is an absolute energy in eV, not a band index. For an
        insulator one may explicitly use the highest occupied energy on this path;
        this is not necessarily the global valence-band maximum.
        """
        if not np.isfinite(reference_ev):
            raise ValueError('reference_ev must be finite')
        ticks = list(ticks)
        for index, label in ticks:
            if not isinstance(index, (int, np.integer)) or not 0 <= index < len(self.energies):
                raise ValueError('tick index outside band path')
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        x = self.distance
        for spin in range(self.energies.shape[2]):
            ax.plot(x, self.energies[:, :, spin]*HARTREE_EV-reference_ev,
                    color=('C0', 'C1')[spin % 2], linewidth=.9)
        if ticks:
            ax.set_xticks([x[i] for i, _ in ticks], [label for _, label in ticks])
            for i, _ in ticks:
                ax.axvline(x[i], color='0.85', linewidth=.6)
        ax.axhline(0, color='0.5', linewidth=.6, linestyle='--')
        ax.set(xlabel='k-path distance (bohr⁻¹)', ylabel=f'Energy − ({reference_ev:.6f} eV) [eV]', title=title)
        return ax


@dataclass(frozen=True, eq=False)
class BandPath:
    """Continuous path in fractional reciprocal coordinates; both endpoints included."""
    vertices: np.ndarray
    intervals: tuple
    labels: tuple

    def __post_init__(self):
        from .sbe import _integer
        v=np.asarray(self.vertices,dtype=float)
        if v.ndim!=2 or v.shape[1]!=3 or not 2<=len(v)<=11 or not np.isfinite(v).all():
            raise ValueError('vertices must be finite (2..11, 3)')
        n=tuple(_integer(x,'intervals',1,1000000) for x in self.intervals)
        labels=tuple(self.labels)
        if len(n)!=len(v)-1 or len(labels)!=len(v) or any(not isinstance(x,str) or not x for x in labels):
            raise ValueError('one interval count per segment and one label per vertex required')
        if sum(n)>1000000: raise ValueError('path exceeds one million intervals')
        object.__setattr__(self,'vertices',np.frombuffer(v.tobytes(),dtype=float).reshape(v.shape))
        object.__setattr__(self,'intervals',n)
        object.__setattr__(self,'labels',labels)

    @property
    def fractional(self):
        return np.concatenate([a+(b-a)*np.arange(n)[:,None]/n for a,b,n in
                               zip(self.vertices[:-1],self.vertices[1:],self.intervals)]+[self.vertices[-1:]])

    def cartesian(self, reciprocal):
        """reciprocal[:, i] is reciprocal basis vector i, in bohr^-1 (includes 2π)."""
        b=np.asarray(reciprocal,dtype=float)
        if b.shape!=(3,3) or not np.isfinite(b).all() or abs(np.linalg.det(b))<1e-15:
            raise ValueError('reciprocal basis must be finite and nonsingular')
        return self.fractional @ b.T

    @property
    def ticks(self):
        return tuple(zip([0,*np.cumsum(self.intervals).tolist()],self.labels))


@dataclass(frozen=True, eq=False)
class ExplicitKPoints:
    """Exact fractional k-point sequence, independent of ASE or path conventions."""
    kpoints: np.ndarray

    def __post_init__(self):
        v=np.asarray(self.kpoints,dtype=float)
        if v.ndim!=2 or v.shape[1]!=3 or not 2<=len(v)<=1000000 or not np.isfinite(v).all():
            raise ValueError('kpoints must be finite (2..1000000, 3)')
        object.__setattr__(self,'kpoints',np.frombuffer(v.tobytes(),dtype=float).reshape(v.shape))

    @property
    def fractional(self): return self.kpoints

    def cartesian(self,reciprocal):
        b=np.asarray(reciprocal,dtype=float)
        if b.shape!=(3,3) or not np.isfinite(b).all() or abs(np.linalg.det(b))<1e-15:
            raise ValueError('reciprocal basis must be finite and nonsingular')
        return self.fractional @ b.T


def _check_report(filename, stored_points, bands):
    """Require versioned native block convergence, not exit code or log heuristics."""
    try:
        rows=Path(filename).read_text().splitlines()
        if rows[0]!='SALMONPY_BANDS_V1': raise ValueError('unsupported report')
        reciprocal=np.array([[float(x) for x in line.split()] for line in rows[1:4]]).T
        nk,nref=map(int,rows[4].split())
        if nk<1 or nref!=bands or stored_points%nk or len(rows)!=5+stored_points//nk:
            raise ValueError('incomplete report')
        for first,line in zip(range(1,stored_points+1,nk),rows[5:]):
            if list(map(int,line.split()))!=[first,nk,1]:
                raise ValueError('unconverged or inconsistent block')
        if reciprocal.shape!=(3,3) or not np.isfinite(reciprocal).all(): raise ValueError('invalid basis')
        return reciprocal
    except (OSError,ValueError,IndexError) as exc:
        raise RuntimeError('Band convergence report missing, incomplete or unconverged; inspect native run') from exc


def _input_groups(text):
    """Split SALMON groups without treating slashes/comments inside strings as syntax.

    Keep group bodies verbatim, including SALMON's non-namelist atom records.
    This is intentionally not a general Fortran namelist evaluator.
    """
    import re
    token=re.compile(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*|&[A-Za-z_][A-Za-z_0-9]*|/")
    groups={};name=None;start=0
    for m in token.finditer(text):
        v=m.group()
        if v.startswith(('!',"'",'"')): continue
        if v.startswith('&'):
            if name is not None: raise ValueError('nested or unterminated input group')
            name=v[1:].lower();start=m.end()
            if name in groups: raise ValueError('duplicate input group')
        elif name is not None:
            groups[name]=text[start:m.start()];name=None
    if name is not None: raise ValueError('unterminated input group')
    return groups


def _scalar(body,name,default=None):
    import re
    # Mask comments while preserving quoted values.
    clean=re.sub(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*",
                 lambda m: '' if m[0].startswith('!') else m[0],body)
    matches=re.findall(r'\b'+re.escape(name)+r'''\s*=\s*('[^']*'|"[^"]*"|[^,\s/]+)''',clean,re.I)
    return matches[-1].strip("'\"").lower() if matches else default


def _validate_gs_input(text,root):
    import re
    from pathlib import PureWindowsPath
    groups=_input_groups(text)
    allowed={'calculation','control','units','system','pseudo','functional','rgrid','kgrid',
             'atomic_coor','atomic_red_coor','scf','analysis','parallel','band'}
    if set(groups)-allowed: raise ValueError('unsupported GS input groups: '+str(set(groups)-allowed))
    system=groups.get('system','');calc=groups.get('calculation','')
    if _scalar(calc,'theory')!='dft': raise ValueError('supply the original DFT ground-state input')
    if any(_scalar(calc,k,'n')!='n' for k in ['yn_md','yn_opt','yn_dc','yn_conventional_from_dcdft']):
        raise ValueError('fixed-ion conventional DFT required')
    if _scalar(system,'yn_periodic')!='y' or _scalar(system,'spin','unpolarized')!='unpolarized':
        raise ValueError('periodic unpolarized ground state required')
    if _scalar(system,'yn_spinorbit','n')!='n': raise ValueError('spin-orbit unsupported')
    if _scalar(groups.get('functional',''),'xc','pz') not in ('pz','libxc_pz','libxc_pbe'):
        raise ValueError('band runner supports PZ, Libxc-PZ and Libxc-PBE')
    # Validate every supplied asset path in retained groups; no implicit external reads.
    for key,body in groups.items():
        if key in ('control','analysis','band','parallel'): continue
        for m in re.finditer(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*",body):
            if m[0].startswith('!'): continue
            value=m[0][1:-1]
            if any(x in value for x in ('/','\\')) or value.startswith('~'):
                path=Path(value)
                if path.is_absolute() or PureWindowsPath(value).drive or '\\' in value or '..' in path.parts or value.startswith('~'):
                    raise ValueError('input assets must be bundle-relative paths')
                if not (root/path).is_file(): raise ValueError('missing relative asset: '+value)
    # Bare filenames are paths too (file_pseudo array or scalar coordinate files).
    for body in groups.values():
        body=re.sub(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*",
                    lambda m: '' if m[0].startswith('!') else m[0],body)
        for m in re.finditer(r'''\b(file_pseudo|file_atom_coor|file_atom_red_coor|file_kgrid)\s*(?:\([^)]*\))?\s*=\s*(['"])(.*?)\2''',body,re.I):
            if m[1].lower() in ('file_atom_coor','file_atom_red_coor') and m[3]=='none': continue
            if not (root/m[3]).is_file(): raise ValueError('missing input asset: '+m[3])
    return groups


@dataclass(frozen=True)
class BandRun:
    directory: Path
    bands: BandStructure
    path: BandPath | ExplicitKPoints

    def plot(self, *, reference_ev, **kwargs):
        if isinstance(self.path,ExplicitKPoints):
            raise ValueError('explicit k points need external path metadata; use the ASE result plot')
        return self.bands.plot(reference_ev=reference_ev,ticks=self.path.ticks,**kwargs)


@dataclass(frozen=True, eq=False)
class BandInitialState:
    """Owned GS snapshot. Use from_salmon and a context manager.

    Binary restart provenance must match the supplied original DFT input; this
    adapter verifies file integrity, not the physical provenance of native binaries.
    """
    executable: object
    manifest: tuple
    _snapshot: Path
    _owner: object
    _origin: Path
    _input: str

    @classmethod
    def from_salmon(cls,directory,input_file,*,executable):
        import tempfile
        from .fdtd import SalmonExecutable
        from .tddft import _copy_bundle
        if not isinstance(executable,SalmonExecutable): raise TypeError('expected SalmonExecutable')
        root=Path(directory).expanduser().resolve();inp=Path(input_file)
        if inp.is_absolute() or '..' in inp.parts: raise ValueError('input_file must be bundle-relative')
        for name in ('info.bin','occupation.bin','wfn.bin'):
            if not (root/'data_for_restart'/name).is_file(): raise ValueError('missing GS record: '+name)
        owner=tempfile.TemporaryDirectory(prefix='salmonpy-band-initial-');snapshot=Path(owner.name)/'bundle'
        try:
            manifest=_copy_bundle(root,snapshot);text=(snapshot/inp).read_text()
            _validate_gs_input(text,snapshot)
            return cls(executable,manifest,snapshot,owner,root,text)
        except BaseException:
            owner.cleanup();raise

    def close(self): self._owner.cleanup()
    def __enter__(self): return self
    def __exit__(self,*args): self.close()

    def run(self,path,directory,*,bands,tolerance=1e-5,max_iterations=300,threads=1,timeout=300.):
        import json,os,subprocess
        from .tddft import _copy_bundle,_digest
        from .sbe import _integer
        from .ld import _positive
        from .outputs import save_arrays
        if not isinstance(path,(BandPath,ExplicitKPoints)): raise TypeError('expected BandPath or ExplicitKPoints')
        bands=_integer(bands,'bands',1,2**31-1)
        iterations=_integer(max_iterations,'max_iterations',1,2**31-1)
        threads=_integer(threads,'threads',1,2**31-1)
        tolerance=_positive(tolerance,'tolerance');timeout=_positive(timeout,'timeout')
        if not self._snapshot.exists(): raise RuntimeError('initial state is closed')
        if _digest(self.executable.path)!=self.executable.sha256: raise ValueError('executable checksum mismatch')
        groups=_validate_gs_input(self._input,self._snapshot)
        nstate=_scalar(groups.get('system',''),'nstate')
        if nstate is None or bands>int(nstate): raise ValueError('bands must not exceed explicit nstate')
        dest=Path(directory).expanduser().resolve()
        if dest.is_relative_to(self._origin) or dest.is_relative_to(self._snapshot):
            raise ValueError('output must be outside the initial bundle')
        dest.mkdir(parents=True,exist_ok=False)
        record={'schema':'salmonpy.band-run/1','status':'started','executable_sha256':self.executable.sha256,
                'initial_files':dict(self.manifest),'fractional_kpoints':path.fractional.tolist(),
                'bands':bands,'tolerance_Ha':tolerance,'max_iterations':iterations,
                'threads':threads,'timeout':timeout}
        if isinstance(path,BandPath):
            record.update(vertices=path.vertices.tolist(),intervals=path.intervals,labels=path.labels)
        manifest=dest/'band-run.json'
        try:
            work=dest/'native'
            if _copy_bundle(self._snapshot,work)!=self.manifest: raise ValueError('initial snapshot changed')
            for name in ('band.dat','band_status.dat'):
                (work/name).unlink(missing_ok=True)
            groups['calculation']="\n theory='dft_band'\n"
            groups['control']="\n sysname='bands', base_directory='./', yn_restart='y', directory_read_data='./data_for_restart/',\n read_gs_restart_data='all', yn_reset_step_restart='y'\n"
            groups['parallel']='\n nproc_k=1, nproc_ob=1, nproc_rgrid=1,1,1\n'
            groups['analysis']='\n'
            groups['scf']=groups.get('scf','').rstrip().rstrip(',')+f'\n nscf={iterations}\n'
            if isinstance(path,ExplicitKPoints):
                point_file=work/'salmonpy-kpoints.data'
                np.savetxt(point_file,path.fractional,fmt='%.17e',header=str(len(path.fractional)),comments='')
                record['kpoints_sha256']=_digest(point_file)
                body=f"\n file_band_kpoints='salmonpy-kpoints.data', nref_band={bands}, tol_esp_diff={tolerance:.17e}\n"
            else:
                body=f"\n lattice='non', nref_band={bands}, num_of_segments={len(path.intervals)}, tol_esp_diff={tolerance:.17e},\n"
                body+='\n'.join(f' ndiv_segment({i})={n},' for i,n in enumerate(path.intervals,1))+'\n'
                body+='\n'.join(f' kpt(:,{i})='+','.join(f'{x:.17e}' for x in v)+',' for i,v in enumerate(path.vertices,1))+'\n'
            groups['band']=body
            text='\n'.join('&'+key+body+'/' for key,body in groups.items())+'\n'
            (work/'band-input.inp').write_text(text)
            record['input_sha256']=_digest(work/'band-input.inp')
            manifest.write_text(json.dumps(record,indent=2))
            env={k:v for k,v in os.environ.items() if k not in ('SALMONPY_TDDFT_WORKER',)}
            env.update(OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1')
            with (work/'stdout.log').open('w') as log:
                completed=subprocess.run([str(self.executable.path)],input=text,text=True,cwd=work,stdout=log,
                                         stderr=subprocess.STDOUT,timeout=timeout,env=env,check=False)
            if completed.returncode or 'end SALMON' not in (work/'stdout.log').read_text():
                raise RuntimeError('SALMON did not finish; inspect '+str(work/'stdout.log'))
            raw=BandStructure.from_salmon(work/'band.dat')
            if raw.energies.shape[1]<bands or raw.energies.shape[2]!=1:
                raise RuntimeError('unexpected native band/spin dimensions')
            reciprocal=_check_report(work/'band_status.dat',len(raw.energies),bands)
            expected=path.cartesian(reciprocal);count=len(expected)
            nk=int((work/'band_status.dat').read_text().splitlines()[4].split()[0])
            if len(raw.energies)!=((count+nk-1)//nk)*nk: raise RuntimeError('unexpected path/padding size')
            if not np.allclose(raw.k_cartesian[:count],expected,rtol=1e-12,atol=1e-12):
                raise RuntimeError('native k path differs from requested reciprocal path')
            if not np.allclose(raw.k_cartesian[count:],expected[-1],rtol=1e-12,atol=1e-12):
                raise RuntimeError('invalid endpoint padding')
            result=BandStructure(raw.k_cartesian[:count],raw.energies[:count,:bands])
            record.update(status='complete',physical_points=count,stored_points=len(raw.energies),
                          reciprocal_columns_bohr_inverse=reciprocal.tolist())
            save_arrays(dest/'result',{'k_cartesian':result.k_cartesian,'energies':result.energies,
                                      'distance':result.distance},metadata={**record,'units':{'energies':'Hartree','k_cartesian':'bohr^-1','distance':'bohr^-1'}})
            return BandRun(dest,result,path)
        except BaseException as exc:
            record.update(status='failed',error=str(exc));raise
        finally:
            manifest.write_text(json.dumps(record,indent=2))
