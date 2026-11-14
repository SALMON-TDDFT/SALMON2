"""Explicit portable GS array interchange (Zarr v3); no native solver invocation.

This does not decode SALMON GS/1 records or reconstruct an SCF/RT restart.
Hashes detect corruption, not publisher authenticity. Do not concurrently edit an
archive while reading, or mutate input arrays while writing.
"""
from dataclasses import dataclass
import hashlib,itertools,json,math
from pathlib import Path,PureWindowsPath
import shutil
from types import MappingProxyType
import numpy as np
from .errors import InputError

__all__=['GroundStateArchive','save_ground_state']
SCHEMA='salmonpy.gs/2'
DTYPES={'float64':np.dtype('<f8'),'int64':np.dtype('<i8')}
UNITS={'wavefunction_real':'bohr^-3/2','wavefunction_imag':'bohr^-3/2','density':'bohr^-3',
       'eigenvalues':'hartree','occupations':'electron','weights':'dimensionless',
       'kpoints':'reduced_reciprocal_2pi','cell_bohr':'bohr','fractional_positions':'fractional','atomic_numbers':'atomic_number'}
AXES={'wavefunction_real':['k','band','z','y','x'],'wavefunction_imag':['k','band','z','y','x'],
      'density':['z','y','x'],'eigenvalues':['k','band'],'occupations':['k','band'],
      'weights':['k'],'kpoints':['k','reciprocal_lattice_vector'],'cell_bohr':['lattice_vector','cartesian'],
      'fractional_positions':['atom','lattice_vector'],'atomic_numbers':['atom']}


def _json(value):return (json.dumps(value,sort_keys=True,allow_nan=False,indent=2)+'\n').encode()

def _write(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)


def _budget(value):
    if value is not None and (type(value) is not int or value<=0):raise InputError('max_bytes must be None or a positive integer')
    return value


def _integer(x,name):
    if type(x) is not int or x<1:raise InputError(name+' must be a positive integer')
    return x


def _real(value,name,shape=None):
    if np.iscomplexobj(value):raise InputError(name+' must be real')
    a=np.asarray(value,dtype=np.float64)
    if shape is not None and a.shape!=shape:raise InputError('invalid '+name+' shape')
    if not np.isfinite(a).all():raise InputError('nonfinite '+name)
    return a


def _scalar(x,name):
    if isinstance(x,(bool,np.bool_)):raise InputError('invalid '+name)
    a=_real(x,name,())
    return float(a)


def _sha(value):return type(value) is str and len(value)==64 and all(c in '0123456789abcdef' for c in value)


def _references(values,numbers):
    if type(values) is not list:raise InputError('pseudopotentials must be a list of references')
    seen=set();result=[]
    for p in values:
        if type(p) is not dict or set(p)!={'atomic_number','filename','sha256','source'}:raise InputError('invalid pseudopotential reference')
        z=p['atomic_number'];name=p['filename'];source=p['source']
        if type(z) is not int or z not in numbers or z in seen:raise InputError('invalid/duplicate pseudopotential species')
        if (type(name) is not str or not name or name in ('.','..') or '/' in name or '\\' in name
                or PureWindowsPath(name).drive or not _sha(p['sha256']) or type(source) is not str or not source.strip()):
            raise InputError('invalid pseudopotential filename/hash/source')
        seen.add(z);result.append(dict(p))
    if seen!=set(numbers):raise InputError('pseudopotential references must cover every species')
    return result


def _layout(name,shape,chunk):
    return dict(zarr_format=3,node_type='array',shape=list(shape),data_type='int64' if name=='atomic_numbers' else 'float64',
                chunk_grid={'name':'regular','configuration':{'chunk_shape':list(chunk)}},
                chunk_key_encoding={'name':'default','configuration':{'separator':'/'}},fill_value=0,
                codecs=[{'name':'bytes','configuration':{'endian':'little'}}],
                attributes={'units':UNITS[name],'axes':AXES[name]})


def _specs(m):
    nk=m['n_kpoints'];nb=m['n_bands'];na=m['n_atoms'];grid=tuple(m['grid_zyx']);block=m['orbital_block']
    shapes={'wavefunction_real':(nk,nb,*grid),'wavefunction_imag':(nk,nb,*grid),'density':grid,
            'eigenvalues':(nk,nb),'occupations':(nk,nb),'kpoints':(nk,3),'weights':(nk,),
            'cell_bohr':(3,3),'fractional_positions':(na,3),'atomic_numbers':(na,)}
    return {name:_layout(name,shape,(1,block,*grid) if name.startswith('wavefunction')
                        else (1,block) if name in ('eigenvalues','occupations') else shape)
            for name,shape in shapes.items()}


def _chunks(spec):
    shape=tuple(spec['shape']);chunk=tuple(spec['chunk_grid']['configuration']['chunk_shape'])
    for index in itertools.product(*(range((n+c-1)//c) for n,c in zip(shape,chunk))):
        rows=tuple(slice(i*c,min((i+1)*c,n)) for i,c,n in zip(index,chunk,shape))
        yield 'c/'+'/'.join(map(str,index)),rows,chunk


def _group():return {'zarr_format':3,'node_type':'group','attributes':{'schema':SCHEMA}}


def _freeze(x):
    if isinstance(x,dict):return MappingProxyType({k:_freeze(v) for k,v in x.items()})
    if isinstance(x,list):return tuple(_freeze(v) for v in x)
    return x


def _readonly(a):return np.frombuffer(a.tobytes(order='C'),dtype=a.dtype).reshape(a.shape)


def _metadata(m):
    keys={'n_kpoints','n_bands','n_atoms','grid_zyx','orbital_block','electrons','energy_hartree','xc',
          'pseudopotentials','boundary','spin','wavefunction_convention','cell_convention',
          'grid_convention','bloch_phase','reciprocal_convention'}
    if type(m) is not dict or not keys<=set(m) or set(m)-keys-{'provenance'}:raise InputError('invalid GS metadata')
    if 'provenance' in m:_export_provenance(m['provenance'])
    for name in ('n_kpoints','n_bands','n_atoms','orbital_block'):_integer(m[name],name)
    if type(m['grid_zyx']) is not list or len(m['grid_zyx'])!=3:raise InputError('invalid GS grid')
    for n in m['grid_zyx']:_integer(n,'grid dimension')
    if m['orbital_block']>m['n_bands']:raise InputError('orbital block exceeds band count')
    if (m['boundary']!='periodic' or m['spin']!='unpolarized' or m['wavefunction_convention']!='cell_periodic_u'
            or m['cell_convention']!='row_lattice_vectors'
            or m['grid_convention']!='fractional_origin_zero_endpoint_excluded'
            or m['bloch_phase']!='exp(+i*k_cartesian_dot_r)'
            or m['reciprocal_convention']!='2pi_inverse_transpose_row_lattice'):raise InputError('unsupported GS conventions')
    if type(m['xc']) is not str or not m['xc'].strip():raise InputError('nonempty XC identifier required')
    if _scalar(m['electrons'],'electrons')<=0:raise InputError('electrons must be positive')
    _scalar(m['energy_hartree'],'energy')


def _export_provenance(p):
    keys={'schema','source_manifest_sha256','target_sha256','exporter_sha256','density_source',
          'eigenvalue_source','native_real_orbital','local_channels','record_environment'}
    if type(p) is not dict or set(p)!=keys or p['schema']!='salmonpy.gs-native-export/1':raise InputError('invalid GS export provenance')
    if any(not _sha(p[k]) for k in ('source_manifest_sha256','target_sha256','exporter_sha256')):raise InputError('invalid GS export provenance hashes')
    if (p['density_source']!='reconstructed_from_native_orbitals' or p['eigenvalue_source']!='native/gs_eigen.data'
            or type(p['native_real_orbital']) is not bool):raise InputError('unsupported GS export semantics')
    channels=p['local_channels'];seen=set()
    if type(channels) is not list or not channels:raise InputError('invalid local channels')
    for c in channels:
        if (type(c) is not dict or set(c)!={'atomic_number','local_channel'} or type(c['atomic_number']) is not int
                or not 1<=c['atomic_number']<=118 or c['atomic_number'] in seen or type(c['local_channel']) is not int
                or not 0<=c['local_channel']<=4):raise InputError('invalid local channels')
        seen.add(c['atomic_number'])
    env=p['record_environment']
    if type(env) is not dict or set(env)!={'GFORTRAN_CONVERT_UNIT'} or type(env['GFORTRAN_CONVERT_UNIT']) is not str:raise InputError('invalid GS record environment')
    if env['GFORTRAN_CONVERT_UNIT'].lower() not in ('','native','swap','big_endian','little_endian'):raise InputError('unsupported GS record environment')


def _physical(arrays,m):
    cell=arrays['cell_bohr'];weights=arrays['weights'];occ=arrays['occupations'];numbers=arrays['atomic_numbers']
    if not np.isfinite(np.linalg.det(cell)) or np.linalg.det(cell)<=0:raise InputError('cell must have positive finite volume')
    if np.any(weights<0) or not np.isclose(weights.sum(),1.,rtol=0,atol=1e-12):raise InputError('k weights must be nonnegative and sum to one')
    if np.any(occ<0) or np.any(occ>2):raise InputError('unpolarized occupations must lie in [0,2]')
    count=float(np.sum(weights[:,None]*occ))
    if not np.isclose(count,m['electrons'],rtol=1e-12,atol=1e-12):raise InputError('occupations/weights disagree with electrons')
    if np.any(numbers<1) or np.any(numbers>118):raise InputError('invalid atomic numbers')
    if np.any(arrays['density']<0):raise InputError('density must be nonnegative')
    _references(m['pseudopotentials'],set(map(int,numbers)))
    if 'provenance' in m and {c['atomic_number'] for c in m['provenance']['local_channels']}!=set(map(int,numbers)):
        raise InputError('local channels do not cover archive species')


def save_ground_state(directory,*,wavefunctions,density,eigenvalues,occupations,kpoints,weights,
                      cell_bohr,fractional_positions,atomic_numbers,electrons,energy_hartree,xc,pseudopotentials,
                      orbital_block=16,provenance=None):
    """Write caller-supplied GS arrays, manifest last; return verified lazy archive.

    Waves are complex128 cell-periodic u[k,band,z,y,x]. No norm/orthogonality,
    SCF convergence or physical density consistency is certified. No overwrite,
    downloads, native-record parsing or solver execution. Inputs must not change.
    """
    waves=np.asarray(wavefunctions,dtype=np.complex128)
    if waves.ndim!=5 or min(waves.shape)<1:raise InputError('wavefunctions must have shape (k,band,z,y,x)')
    nk,nb,*grid=waves.shape
    block=min(_integer(orbital_block,'orbital_block'),nb)
    positions=_real(fractional_positions,'fractional_positions')
    if positions.ndim!=2 or positions.shape[1]!=3 or len(positions)<1:raise InputError('positions must be (atom,3)')
    raw_numbers=np.asarray(atomic_numbers)
    if raw_numbers.dtype.kind not in 'iu' or any(type(n) in (bool,np.bool_) for n in atomic_numbers):raise InputError('atomic_numbers must be integers')
    numbers=np.asarray(atomic_numbers,dtype=np.int64)
    if numbers.shape!=(len(positions),):raise InputError('invalid atomic_numbers shape')
    m=dict(n_kpoints=nk,n_bands=nb,n_atoms=len(positions),grid_zyx=grid,orbital_block=block,
           electrons=_scalar(electrons,'electrons'),energy_hartree=_scalar(energy_hartree,'energy'),xc=xc,
           pseudopotentials=_references(pseudopotentials,set(map(int,numbers))),boundary='periodic',spin='unpolarized',
           wavefunction_convention='cell_periodic_u',cell_convention='row_lattice_vectors',
           grid_convention='fractional_origin_zero_endpoint_excluded',
           bloch_phase='exp(+i*k_cartesian_dot_r)',reciprocal_convention='2pi_inverse_transpose_row_lattice')
    if provenance is not None:m['provenance']=json.loads(_json(provenance))
    _metadata(m);specs=_specs(m)
    arrays={'wavefunction_real':waves.real,'wavefunction_imag':waves.imag,'atomic_numbers':numbers}
    supplied=dict(density=density,eigenvalues=eigenvalues,occupations=occupations,kpoints=kpoints,weights=weights,
                  cell_bohr=cell_bohr,fractional_positions=positions)
    for name,value in supplied.items():arrays[name]=_real(value,name,tuple(specs[name]['shape']))
    _physical(arrays,m)
    root=Path(directory).expanduser();root.mkdir(parents=True,exist_ok=False)
    try:
        manifest={'schema':SCHEMA,'metadata':m,'arrays':specs,'files':{}}
        def put(name,raw):
            _write(root/name,raw);manifest['files'][name]=hashlib.sha256(raw).hexdigest()
        put('zarr.json',_json(_group()))
        for name,spec in specs.items():
            put(name+'/zarr.json',_json(spec))
            for key,rows,extent in _chunks(spec):
                part=arrays[name][rows]
                if not np.isfinite(part).all():raise InputError('nonfinite '+name)
                chunk=np.zeros(extent,dtype=DTYPES[spec['data_type']]);chunk[tuple(slice(0,n) for n in part.shape)]=part
                put(name+'/'+key,chunk.tobytes(order='C'))
        _write(root/'manifest.json',_json(manifest))
        return GroundStateArchive.open(root)
    except BaseException:
        shutil.rmtree(root);raise


@dataclass(frozen=True)
class GroundStateArchive:
    """Verified immutable metadata + lazy block reads; hashes rechecked on access."""
    directory:Path
    metadata:object
    _arrays:object
    _hashes:object
    schema:str=SCHEMA

    @property
    def wavefunction_shape(self):return tuple(self._arrays['wavefunction_real']['shape'])

    @classmethod
    def open(cls,directory,*,max_bytes=None):
        """Verify inventory/size/hash/layout and physical metadata, without loading waves.

        Optional max_bytes limits total archive file bytes, not process RSS. There
        is no machine-independent default ceiling; payload sizes derive from the
        declared supported layout and are checked before array allocation.
        """
        limit=_budget(max_bytes);root=Path(directory).expanduser()
        try:
            if root.is_symlink() or not root.is_dir():raise InputError('GS archive must be a real directory')
            mp=root/'manifest.json'
            if mp.is_symlink() or not mp.is_file():raise InputError('missing GS manifest')
            if limit is not None and mp.stat().st_size>limit:raise InputError('archive exceeds max_bytes')
            manifest=json.loads(mp.read_bytes())
            if type(manifest) is not dict or set(manifest)!={'schema','metadata','arrays','files'} or manifest['schema']!=SCHEMA:raise InputError('unsupported GS manifest')
            m=manifest['metadata'];_metadata(m);specs=_specs(m)
            if manifest['arrays']!=specs:raise InputError('unsupported GS array layout')
            hashes=manifest['files']
            count=m['n_kpoints']*((m['n_bands']+m['orbital_block']-1)//m['orbital_block'])
            if type(hashes) is not dict or len(hashes)!=17+4*count:
                raise InputError('GS archive file count disagrees with dimensions')
            expected={'zarr.json':None}
            for name,spec in specs.items():
                expected[name+'/zarr.json']=None
                for key,_,extent in _chunks(spec):expected[name+'/'+key]=math.prod(extent)*8
            if set(hashes)!=set(expected):raise InputError('GS archive file inventory differs')
            required=set(expected)|{'manifest.json'};allowed=required.copy()
            for name in required:allowed.update(p.as_posix() for p in Path(name).parents if p!=Path('.'))
            paths={};total=0
            for p in root.rglob('*'):
                name=p.relative_to(root).as_posix()
                # Reject an unknown directory before the iterator descends into it.
                if name not in allowed or p.is_symlink():raise InputError('unexpected GS entry or symlink')
                if p.is_dir():
                    if name in required:raise InputError('GS file is a directory')
                    continue
                if not p.is_file() or name not in required:raise InputError('non-regular GS entry')
                total+=p.stat().st_size
                if limit is not None and total>limit:raise InputError('archive exceeds max_bytes')
                paths[name]=p
            if set(paths)!=required:raise InputError('GS archive file inventory differs')
            for name,size in expected.items():
                if not _sha(hashes[name]):raise InputError('invalid GS checksum')
                if size is not None and paths[name].stat().st_size!=size:raise InputError('invalid GS chunk size')
                data=paths[name].read_bytes()
                if hashlib.sha256(data).hexdigest()!=hashes[name]:raise InputError('GS checksum mismatch: '+name)
                if name=='zarr.json' and json.loads(data)!=_group():raise InputError('invalid GS root metadata')
                if name.endswith('/zarr.json') and json.loads(data)!=specs[name.split('/')[0]]:raise InputError('GS array metadata differs')
                if size is not None and not np.isfinite(np.frombuffer(data,dtype=DTYPES[specs[name.split('/')[0]]['data_type']])).all():raise InputError('nonfinite GS chunk')
            result=cls(root.resolve(),_freeze(m),_freeze(specs),_freeze(hashes))
            # No full wavefunction allocation. Small metadata and one density grid are checked.
            _physical({name:result.array(name) for name in ('cell_bohr','weights','occupations','atomic_numbers','density')},m)
            return result
        except (OSError,ValueError,TypeError,KeyError,OverflowError,RecursionError) as exc:
            if isinstance(exc,InputError):raise
            raise InputError('invalid or unreadable GS archive') from exc

    def tddft_initial_state(self,*,converter,executable,pseudopotentials,max_steps=1000,threads=1,timeout=300.):
        """Explicit Γ-only native conversion; requires paired audited converter/runtime."""
        from .gs_native import GammaRestartConverter
        if not isinstance(converter,GammaRestartConverter):raise InputError('GammaRestartConverter required')
        return converter.tddft_initial_state(self,executable=executable,pseudopotentials=pseudopotentials,
                                            max_steps=max_steps,threads=threads,timeout=timeout)

    def _read(self,name,key):
        p=self.directory/name/key
        if p.is_symlink() or any(parent.is_symlink() for parent in p.parents):raise InputError('archive path changed to symlink')
        spec=self._arrays[name];extent=spec['chunk_grid']['configuration']['chunk_shape']
        if not p.is_file() or p.stat().st_size!=math.prod(extent)*8:raise InputError('GS chunk size changed after open')
        data=p.read_bytes()
        if hashlib.sha256(data).hexdigest()!=self._hashes[name+'/'+key]:raise InputError('GS chunk changed after open')
        return data

    def array(self,name,*,max_bytes=None):
        """Read one non-wavefunction array; immutable result. Use wavefunctions for waves."""
        limit=_budget(max_bytes)
        if name not in self._arrays or name.startswith('wavefunction'):raise InputError('unknown or wavefunction array name')
        spec=self._arrays[name];shape=tuple(spec['shape']);size=math.prod(shape)*8
        if limit is not None and size>limit:raise InputError('decoded array exceeds max_bytes')
        a=np.empty(shape,dtype=DTYPES[spec['data_type']])
        for key,rows,extent in _chunks(spec):
            part=np.frombuffer(self._read(name,key),dtype=a.dtype).reshape(extent)
            a[rows]=part[tuple(slice(0,s.stop-s.start) for s in rows)]
        return _readonly(a)

    def wavefunctions(self,k,start=0,stop=None,*,max_bytes=None):
        """Read bands [start,stop) at one k; immutable complex128 (band,z,y,x)."""
        limit=_budget(max_bytes);nk,nb,*grid=self.wavefunction_shape
        if stop is None:stop=nb
        if any(type(n) is not int for n in (k,start,stop)) or not (0<=k<nk and 0<=start<stop<=nb):raise InputError('invalid k/band range')
        spec=self._arrays['wavefunction_real'];block=spec['chunk_grid']['configuration']['chunk_shape'][1]
        size=(stop-start)*math.prod(grid)*16
        if limit is not None and size>limit:raise InputError('decoded wavefunctions exceed max_bytes')
        result=np.empty((stop-start,*grid),dtype=np.complex128)
        for name,target in (('wavefunction_real',result.real),('wavefunction_imag',result.imag)):
            for b in range(start//block,(stop-1)//block+1):
                key=f'c/{k}/{b}/0/0/0';extent=(1,block,*grid)
                data=np.frombuffer(self._read(name,key),dtype='<f8').reshape(extent)[0]
                lo=max(start,b*block);hi=min(stop,(b+1)*block)
                target[lo-start:hi-start]=data[lo-b*block:hi-b*block]
        return _readonly(result)
