"""Strict Zarr v3 subset for TDDFT replay history, without GS payloads.

A caller supplies the model and executable; this module never loads executable
paths from disk. Hashes detect corruption, not an untrusted publisher's identity.
"""
import hashlib
import json
import math
from pathlib import Path
import shutil
import numpy as np
from .errors import SnapshotError
from .snapshots import StateSnapshot, _SCHEMA, _integrity
from ._version import __version__

SCHEMA = 'salmonpy.tddft-checkpoint/1'
CHUNK_ROWS = 1024
NAMES = {'A0':'_A0', 'J0':'_J0', 'current':'_current',
         'history':'_history', 'current_history':'_curr_history'}



def _format(state):
    if state._snapshot_kind == 'isolated-TDDFT-E-dipole-replay':
        return ('salmonpy.isolated-tddft-checkpoint/2',
                {'E0':'_A0', 'observables0':'_J0', 'observables':'_current',
                 'E_history':'_history', 'observable_history':'_curr_history'}, 7)
    return SCHEMA, NAMES, 3


def _shape(private, rows, width):
    columns = 3 if private in ('_A0', '_history') else width
    return (rows, columns) if private in ('_history', '_curr_history') else (columns,)


def _json(data):
    return (json.dumps(data, sort_keys=True, allow_nan=False, indent=2)+'\n').encode('utf-8')


def _write(path, raw):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)


def _layout(shape, chunk):
    return dict(zarr_format=3, node_type='array', shape=list(shape), data_type='float64',
                chunk_grid={'name':'regular','configuration':{'chunk_shape':list(chunk)}},
                chunk_key_encoding={'name':'default','configuration':{'separator':'/'}},
                fill_value=0., codecs=[{'name':'bytes','configuration':{'endian':'little'}}],
                attributes={})


def _binding(state):
    return dict(settings_sha256=state._snapshot_settings(), dt=state.dt,
                capacity=state._capacity, threads=state._threads,
                executable_sha256=state._executable.sha256,
                initial_files=dict(sorted(state._initial_hashes.items())))


def _limits(state):
    """Format-derived byte ceilings for this explicitly configured capacity.

    Input and observation float64 values, three fixed vectors, full edge chunks,
    metadata and every hash entry are counted. No machine/RSS estimate is used.
    """
    schema,names,width=_format(state)
    capacity=state._capacity
    count=(capacity+CHUNK_ROWS-1)//CHUNK_ROWS
    file_count=9+2*count
    arrays={}
    for name,private in names.items():
        shape=_shape(private,capacity,width)
        arrays[name]=_layout(shape,(CHUNK_ROWS,shape[-1]) if len(shape)==2 else shape)
    template=dict(schema=schema,version=__version__,kind=state._snapshot_kind,
                  step_index=capacity,time=0.,binding=_binding(state),arrays=arrays,files={})
    longest=max([*[name+'/c/'+str(max(count-1,0))+'/0' for name in names],
                 *[name+'/zarr.json' for name in names]],key=len)
    # Each JSON hash entry is shorter than this standalone one-entry object.
    entry_size=len(_json({longest:'f'*64}))
    # A finite IEEE binary64 JSON clock fits within 32 characters.
    manifest_size=len(_json(template))+32+file_count*entry_size
    group_size=len(_json({'zarr_format':3,'node_type':'group','attributes':{'schema':schema}}))
    metadata_size=sum(len(_json(spec)) for spec in arrays.values())
    padded=capacity if capacity<=CHUNK_ROWS else count*CHUNK_ROWS
    return manifest_size+group_size+metadata_size+8*(3+2*width)+8*(3+width)*padded,manifest_size


def _chunks(shape, chunk):
    if len(shape)==1:
        yield 'c/0', slice(None), shape
    else:
        for i, start in enumerate(range(0, shape[0], chunk[0])):
            yield f'c/{i}/0', slice(start, min(start+chunk[0], shape[0])), chunk


def save(state, directory):
    schema,names,width=_format(state)
    token=state.snapshot()
    root=Path(directory).expanduser()
    for source in (state._bundle,state._model._snapshot,state._model._source_root):
        if source is not None and root.resolve().is_relative_to(source.resolve()):
            raise SnapshotError('checkpoint destination must be outside initial bundles')
    # mkdir claims only a new destination; never clean somebody else's directory.
    root.mkdir(parents=True, exist_ok=False)
    try:
        manifest=dict(schema=schema, version=__version__, kind=token.kind,
                      step_index=token.step_index, time=token.time, binding=_binding(state),
                      arrays={}, files={})
        def put(name,raw):
            _write(root/name,raw)
            manifest['files'][name]=hashlib.sha256(raw).hexdigest()
        put('zarr.json',_json({'zarr_format':3,'node_type':'group','attributes':{'schema':schema}}))
        # Use the immutable token, rather than rereading live arrays during I/O.
        frozen={name:np.frombuffer(raw,dtype=dtype).reshape(shape,order=order)
                for name,shape,dtype,order,raw in token._data}
        for public, private in names.items():
            a=frozen[private]
            chunk=a.shape if a.ndim==1 else (min(max(a.shape[0],1),CHUNK_ROWS),a.shape[1])
            spec=_layout(a.shape,chunk);manifest['arrays'][public]=spec
            put(f'{public}/zarr.json',_json(spec))
            for key, rows, extent in _chunks(a.shape,chunk):
                block=np.zeros(extent,dtype='<f8')
                part=a[rows]
                if a.ndim==1:block[:]=part
                else:block[:len(part)]=part
                put(f'{public}/{key}',block.tobytes(order='C'))
        # Readers require this final commit marker. Interrupted writes cannot load.
        _write(root/'manifest.json',_json(manifest))
    except BaseException:
        shutil.rmtree(root)
        raise


def load(state, directory, max_bytes):
    if max_bytes is not None and (type(max_bytes) is not int or max_bytes<=0):
        raise SnapshotError('max_bytes must be None or a positive integer')
    schema,names,width=_format(state)
    root=Path(directory).expanduser()
    try:
        ceiling,manifest_ceiling=_limits(state)
        budget=ceiling if max_bytes is None else min(ceiling,max_bytes)
        if root.is_symlink() or not root.is_dir():raise SnapshotError('checkpoint must be a real directory')
        manifest_path=root/'manifest.json'
        if (manifest_path.is_symlink() or not manifest_path.is_file()
                or manifest_path.stat().st_size>min(budget,manifest_ceiling)):
            raise SnapshotError('missing or oversized checkpoint manifest')
        m=json.loads(manifest_path.read_bytes())
        if (type(m) is not dict or m.get('schema')!=schema or m.get('version')!=__version__
                or m.get('kind')!=state._snapshot_kind or m.get('binding')!=_binding(state)):
            raise SnapshotError('checkpoint requires identical model, dt, capacity, threads and executable')
        n=m['step_index'];clock=m['time']
        if (type(n) is not int or not 0<=n<=state._capacity or type(clock) not in (int,float)
                or not math.isfinite(clock) or clock!=n*state.dt):
            raise SnapshotError('invalid checkpoint clock')
        shapes={name:_shape(private,n,width) for name,private in names.items()}
        if type(m['arrays']) is not dict or set(m['arrays'])!=set(names):raise SnapshotError('invalid checkpoint array names')
        hashes=m['files']
        count=(n+CHUNK_ROWS-1)//CHUNK_ROWS
        if type(hashes) is not dict or len(hashes)!=9+2*count:
            raise SnapshotError('checkpoint declares an invalid file count')
        if 8*(3+2*width)+8*(3+width)*n>budget:raise SnapshotError('checkpoint arrays exceed max_bytes')
        expected={'zarr.json':None};specs={}
        for name,shape in shapes.items():
            chunk=shape if len(shape)==1 else (min(max(n,1),CHUNK_ROWS),shape[1])
            if m['arrays'][name]!=_layout(shape,chunk):raise SnapshotError('unsupported checkpoint array layout')
            expected[f'{name}/zarr.json']=None;specs[name]=(shape,chunk)
            for key,_,extent in _chunks(shape,chunk):expected[f'{name}/{key}']=math.prod(extent)*8
        if set(hashes)!=set(expected):raise SnapshotError('checkpoint file inventory differs')
        # Only canonical paths and their parents are allowed. An unexpected
        # directory is rejected before descending into it, including empty ones.
        required=set(expected)|{'manifest.json'}
        allowed=required.copy()
        for name in expected:
            allowed.update(p.as_posix() for p in Path(name).parents if p!=Path('.'))
        files={};total=0
        for path in root.rglob('*'):
            relative=path.relative_to(root).as_posix()
            if relative not in allowed or path.is_symlink():raise SnapshotError('unexpected checkpoint entry or symlink')
            if path.is_dir():
                if relative in expected or relative=='manifest.json':raise SnapshotError('checkpoint file is a directory')
                continue
            if not path.is_file() or relative not in required:
                raise SnapshotError('checkpoint contains a non-regular entry')
            total+=path.stat().st_size
            if total>budget:raise SnapshotError('checkpoint exceeds max_bytes')
            files[relative]=path
        if set(files)!=required:raise SnapshotError('checkpoint file inventory differs')
        # Verify all sizes and hashes before allocating replay arrays or spawning.
        for name,size in expected.items():
            path=files[name]
            if size is not None and path.stat().st_size!=size:raise SnapshotError('invalid checkpoint chunk length')
            raw=path.read_bytes()
            if hashlib.sha256(raw).hexdigest()!=hashes[name]:raise SnapshotError('checkpoint checksum mismatch')
            if name=='zarr.json':
                if json.loads(raw)!={'zarr_format':3,'node_type':'group','attributes':{'schema':schema}}:
                    raise SnapshotError('invalid checkpoint group metadata')
            elif name.endswith('/zarr.json') and json.loads(raw)!=m['arrays'][name.split('/')[0]]:
                raise SnapshotError('checkpoint array metadata differs')
        data=[]
        for name,private in names.items():
            shape,chunk=specs[name];a=np.empty(shape,dtype=float)
            for key,rows,extent in _chunks(shape,chunk):
                raw=files[f'{name}/{key}'].read_bytes()
                if hashlib.sha256(raw).hexdigest()!=hashes[f'{name}/{key}']:raise SnapshotError('checkpoint changed during load')
                part=np.frombuffer(raw,dtype='<f8').reshape(extent)
                if len(shape)==1:a[:]=part
                else:a[rows]=part[:len(a[rows])]
            if not np.isfinite(a).all():raise SnapshotError('nonfinite checkpoint history')
            data.append((private,shape,a.dtype.str,'C',a.tobytes(order='C')))
        data=tuple(data);settings=state._snapshot_settings()
        return StateSnapshot(kind=state._snapshot_kind,step_index=n,time=clock,schema=_SCHEMA,
            version=__version__,_owner=state._snapshot_owner(),_settings=settings,_data=data,
            _integrity=_integrity(state._snapshot_kind,n,clock,settings,data))
    except (OSError,ValueError,TypeError,KeyError,OverflowError,RecursionError) as exc:
        if isinstance(exc,SnapshotError):raise
        raise SnapshotError('invalid or unreadable TDDFT checkpoint') from exc
