"""Opt-in, versioned postprocessing archives; no native runtime is needed to read."""
from .errors import InputError
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import shutil
import numpy as np
from .analysis import read_table
from .fields import _readonly


@dataclass(frozen=True,slots=True)
class SpinDynamics:
    """Native SOC time series; no unit, sign or spin-normalization conversion.

    pauli_current axes: (sample, Pauli index 0/x/y/z, flow direction x/y/z).
    Magnetization retains native mag_x/y/z and its component-sign convention,
    not a converted magnetic moment in Bohr magnetons. Spin-current columns retain the
    native Pauli-matrix convention, without an extra hbar/2 factor.
    """
    time: np.ndarray
    magnetization: np.ndarray
    pauli_current: np.ndarray
    time_unit: str
    magnetization_unit: str
    current_unit: str
    headers: tuple
    source: str

    @property
    def spin_current(self):
        """(sample, spin polarization x/y/z, flow direction x/y/z)."""
        return self.pauli_current[:,1:,:]


def read_spin_dynamics(path):
    """Read SALMON *_rt_spin.data as owned, immutable arrays.

    Accepts one path or an ordered list/tuple of restart output paths.
    Matching overlap rows are removed; times are never shifted.
    Sampling times belong to these files. Do not align with current/energy by
    row number: output strides and inclusion of step0 may differ.
    Channel0 is retained, not silently relabeled as charge current j=-Jm.
    All physical observables are computed by the existing Fortran writer.
    """
    table=read_table(path)
    names=('Time','mag_x','mag_y','mag_z')+tuple(
        'spin_curr_'+spin+'_'+direction for spin in ('0','x','y','z') for direction in ('x','y','z'))
    if set(table.columns)!=set(names):raise InputError('expected native spin dynamics columns')
    indices=[table.columns.index(n) for n in names]
    units=tuple(table.units[i] for i in indices)
    if any(u is None for u in units) or len(set(units[1:4]))!=1 or len(set(units[4:]))!=1:
        raise InputError('consistent labeled spin dynamics units required')
    values=table.values[:,indices]
    return SpinDynamics(_readonly(values[:,0]),_readonly(values[:,1:4]),
                        _readonly(values[:,4:].reshape(-1,4,3)),units[0],units[1],units[4],
                        table.headers,table.source)


def _indices(value):
    a=np.asarray(value,dtype=np.int64)
    return np.frombuffer(a.tobytes(),dtype=a.dtype).reshape(a.shape)


@dataclass(frozen=True)
class CubeGrid:
    """Scalar grid; axes are row vectors: position = origin + index @ axes."""
    values: np.ndarray
    origin: np.ndarray
    axes: np.ndarray
    atoms: np.ndarray
    comments: tuple
    unit_note: str


def read_cube(path):
    """Read SALMON scalar CUBE output; preserve its unit comment without conversion.

    Multi-orbital and negative-count cube variants are deliberately rejected.
    """
    with Path(path).open() as stream:
        comments=(stream.readline().rstrip(),stream.readline().rstrip())
        header=stream.readline().split()
        if len(header)!=4: raise InputError('scalar cube header required')
        nat=int(header[0]);origin=[float(x) for x in header[1:]]
        grid=[stream.readline().split() for _ in range(3)]
        if nat<0 or any(len(row)!=4 or int(row[0])<=0 for row in grid):
            raise InputError('only positive-count scalar cubes supported')
        shape=tuple(int(row[0]) for row in grid);axes=[[float(x) for x in row[1:]] for row in grid]
        atoms=[[float(x) for x in stream.readline().split()] for _ in range(nat)]
        if any(len(row)!=5 for row in atoms): raise InputError('invalid cube atom records')
        values=np.array([float(x.replace('D','E').replace('d','e')) for x in stream.read().split()])
    if values.size!=int(np.prod(shape)): raise InputError('cube scalar count does not match grid')
    arrays=[values.reshape(shape),np.asarray(origin),np.asarray(axes),np.asarray(atoms).reshape(nat,5)]
    if not all(np.isfinite(a).all() for a in arrays): raise InputError('cube contains nonfinite values')
    return CubeGrid(*[_readonly(a) for a in arrays],comments,comments[1])


@dataclass(frozen=True)
class DecomposedCurrent:
    step: np.ndarray
    indices: np.ndarray
    current: np.ndarray
    units: tuple
    headers: tuple


def read_decomposed_current(path):
    """Read weighted k/band/spin current blocks, retaining native step indices.

    current shape is (sample, entry, xyz); indices columns are (k, band, spin).
    Multiply step by the run's dt to obtain time. No additional occupation or
    k-point weight should be applied: SALMON's writer has already applied both.
    """
    steps=[];blocks=[];headers=[]
    for line in Path(path).read_text().splitlines():
        text=line.strip()
        if not text: continue
        if text.startswith(('#','!')): headers.append(text);continue
        values=text.split()
        if len(values)==1:
            steps.append(int(values[0]));blocks.append([])
        elif len(values)==6 and blocks:
            blocks[-1].append([float(x.replace('D','E').replace('d','e')) for x in values])
        else: raise InputError('invalid decomposed-current block')
    if not blocks or any(not b or len(b)!=len(blocks[0]) for b in blocks):
        raise InputError('incomplete decomposed-current blocks')
    a=np.asarray(blocks)
    if not np.isfinite(a).all() or np.any(np.diff(steps)<=0): raise InputError('invalid block values or steps')
    ids=a[0,:,:3]
    if np.any(ids!=np.floor(ids)) or np.any(ids<1) or not np.all(a[:,:,:3]==ids):
        raise InputError('inconsistent k/band/spin indices')
    units={int(n):u for n,_,u in re.findall(r'(\d+):([^\s\[]+)\[([^\]]+)\]','\n'.join(headers))}
    return DecomposedCurrent(_indices(steps),_indices(ids),
        _readonly(a[:,:,3:]),tuple(units.get(i) for i in (4,5,6)),tuple(headers))


def _digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def _member(root,name):
    p=Path(name)
    if p.is_absolute() or '..' in p.parts: raise InputError('archive-relative path required')
    target=root/p
    if not target.resolve().is_relative_to(root.resolve()): raise InputError('path leaves archive')
    return target


def _validate_array(value,spec,name):
    if list(value.shape)!=spec['shape'] or str(value.dtype)!=spec['dtype']:
        raise InputError('array shape or dtype differs from manifest: '+name)
    if value.dtype.kind not in 'biufc' or not np.isfinite(value).all():
        raise InputError('finite numerical array required: '+name)


@dataclass(frozen=True)
class OutputArchive:
    directory: Path
    metadata: dict
    files: tuple
    arrays: dict

    @classmethod
    def open(cls,directory,*,verify=True):
        root=Path(directory).expanduser().resolve()
        manifest=json.loads((root/'manifest.json').read_text())
        if manifest.get('schema')!='salmonpy.outputs/1': raise InputError('unsupported output schema')
        entries=manifest.get('files')
        specs=manifest.get('arrays')
        if not isinstance(entries,list) or not isinstance(specs,dict):
            raise InputError('archive file and array indexes are required')
        names=[]
        for entry in entries:
            if not isinstance(entry,dict) or not isinstance(entry.get('path'),str):
                raise InputError('invalid archive file entry')
            name=entry['path']
            if name in names: raise InputError('duplicate archive file entry: '+name)
            names.append(name)
            p=_member(root,name)
            if not p.is_file(): raise InputError('missing archive file: '+name)
            if not isinstance(entry.get('sha256'),str) or not re.fullmatch('[0-9a-f]{64}',entry['sha256']):
                raise InputError('invalid archive checksum: '+name)
            if verify and _digest(p)!=entry['sha256']: raise InputError('output checksum mismatch: '+name)
        if 'arrays.npz' not in names:
            raise InputError('arrays.npz must be indexed exactly once')
        for name,spec in specs.items():
            if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',name) or not isinstance(spec,dict):
                raise InputError('invalid array index')
            shape=spec.get('shape')
            if not isinstance(shape,list) or any(type(n) is not int or n<0 for n in shape):
                raise InputError('invalid array shape: '+name)
            if not isinstance(spec.get('dtype'),str): raise InputError('invalid array dtype: '+name)
        with np.load(root/'arrays.npz',allow_pickle=False) as data:
            if len(data.files)!=len(specs) or set(data.files)!=set(specs):
                raise InputError('array names differ from manifest')
            for name in data.files:
                _validate_array(data[name],specs[name],name)
        return cls(root,manifest['metadata'],tuple(entries),specs)

    def array(self,name):
        if name not in self.arrays: raise InputError('array is not indexed: '+str(name))
        path=_member(self.directory,'arrays.npz')
        if 'arrays.npz' not in {entry['path'] for entry in self.files}:
            raise InputError('arrays.npz is not indexed')
        with np.load(path,allow_pickle=False) as data:
            value=data[name].copy()
        _validate_array(value,self.arrays[name],name)
        value.setflags(write=False)
        return value

    def table(self,name):
        path=_member(self.directory,name)
        if name not in {f['path'] for f in self.files}: raise InputError('file is not indexed in this archive')
        return read_table(path)

    def cube(self,name):
        path=_member(self.directory,name)
        if name not in {f['path'] for f in self.files}: raise InputError('file is not indexed in this archive')
        return read_cube(path)

    def decomposed_current(self,name):
        path=_member(self.directory,name)
        if name not in {f['path'] for f in self.files}: raise InputError('file is not indexed in this archive')
        return read_decomposed_current(path)


def save_arrays(directory,arrays,*,metadata):
    """Save selected numerical arrays with caller-specified units and conventions.

    Suitable for LD P/J, FDTD fields, SBE density matrices and analysis results.
    No units or physical interpretation are inferred. Existing runs are never overwritten.
    """
    return _write_archive(directory,arrays,metadata)


def _write_archive(directory,arrays,metadata,*,native_source=None,initial_hashes=None,input_text=None):
    checked={}
    for name,value in arrays.items():
        if not isinstance(name,str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',name):
            raise InputError('array names must be simple identifiers')
        a=np.asarray(value)
        if a.dtype.kind not in 'biufc' or not np.isfinite(a).all(): raise InputError('finite numerical arrays required')
        checked[name]=a
    json.dumps(metadata,allow_nan=False)  # validate before creating any output
    root=Path(directory).expanduser().resolve()
    if native_source is not None and root.is_relative_to(Path(native_source).resolve()):
        raise InputError('output directory must be outside native working directory')
    root.mkdir(parents=True,exist_ok=False)
    try:
        np.savez_compressed(root/'arrays.npz',**checked)
        if input_text is not None: (root/'input.inp').write_text(input_text)
        if native_source is not None:
            source=Path(native_source);initial_hashes={Path(k).as_posix():v for k,v in (initial_hashes or {}).items()}
            for p in sorted(source.rglob('*')):
                if p.is_symlink(): raise InputError('native output symlinks are unsupported')
                if not p.is_file(): continue
                rel=p.relative_to(source).as_posix()
                if rel.startswith('data_for_restart/') and initial_hashes.get(rel)==_digest(p): continue
                target=root/'native'/rel;target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(p,target)
        files=[]
        for p in sorted(root.rglob('*')):
            if not p.is_file(): continue
            record={'path':p.relative_to(root).as_posix(),'bytes':p.stat().st_size,'sha256':_digest(p)}
            if record['path'].startswith('native/'):
                rel=record['path'][7:]
                record['content_matches_initial_bundle']=(initial_hashes or {}).get(rel)==record['sha256']
            if p.suffix in ('.data','.cube'):
                headers=[]
                with p.open(errors='replace') as stream:
                    for _ in range(32):
                        line=stream.readline()
                        if not line: break
                        if p.suffix=='.data' and not line.lstrip().startswith(('#','!')): break
                        if p.suffix=='.cube' and len(headers)==2: break
                        headers.append(line.rstrip())
                record['headers']=headers
                record['columns']=[{'index':int(n),'name':label,'unit':unit}
                    for n,label,unit in re.findall(r'(\d+):([^\s\[]+)\[([^\]]+)\]','\n'.join(headers))]
            files.append(record)
        manifest={'schema':'salmonpy.outputs/1','metadata':metadata,
                  'arrays':{k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in checked.items()},'files':files}
        (root/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False))
        return OutputArchive.open(root)
    except BaseException:
        shutil.rmtree(root)
        raise
