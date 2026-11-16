"""Same-topology direct checkpoint coordinator for the native mixed fixture.

No pickle, no propagation during load. A completed manifest binds every rank file.
This internal example helper is not a general user-backend archive API.
"""
import ctypes as c
import hashlib
import json
from pathlib import Path
from dataclasses import replace
import numpy as np


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''): h.update(block)
    return h.hexdigest()


def checkpoint(world, coupled, directory, signature, *, restore=False, frames=None):
    r = coupled.response
    field = coupled.field
    root = Path(directory).resolve()
    error = None
    try:
        for state in r.states.values():
            if hasattr(state, 'native'):
                getattr(state.lib, 'sp_tddft_checkpoint')
            elif not callable(getattr(getattr(state, '_raw', None), 'snapshot', None)):
                raise ValueError('backend has no direct checkpoint capability')
        if field is not None: getattr(field.lib, 'sp_eh_checkpoint')
        if restore:
            manifest = json.loads((root/'complete.json').read_text())
            if manifest['schema'] != 'salmonpy.native-direct/1' or manifest['signature'] != signature or manifest['ranks'] != world.size:
                raise ValueError('checkpoint input/library/topology mismatch')
            if not 0 < manifest['step'] < signature['steps']:
                raise ValueError('invalid checkpoint step')
            if not ({f'rank-{i}.json' for i in range(world.size)} | {'frames.json'}) <= set(manifest['files']):
                raise ValueError('missing rank checkpoint')
            for name, sha in manifest['files'].items():
                p = root/name
                if p.parent != root or digest(p) != sha: raise ValueError('checkpoint checksum mismatch')
            if set(manifest['files']) != {p.name for p in root.iterdir() if p.name != 'complete.json'}:
                raise ValueError('checkpoint file inventory mismatch')
            data = json.loads((root/f'rank-{world.rank}.json').read_text())
            if data['step'] != manifest['step']:
                raise ValueError('checkpoint epochs differ')
            count = len(r.context.plan.point_ids)
            for name in ('J', 'response_A'):
                a = np.asarray(data[name], dtype=float)
                if a.shape != (count,3) or not np.isfinite(a).all():
                    raise ValueError('invalid checkpoint coupling arrays')
            if field is not None:
                for name in ('field_A','field_E'):
                    a = np.asarray(data[name], dtype=float)
                    if a.shape != (count,3) or not np.isfinite(a).all():
                        raise ValueError('invalid checkpoint field arrays')
            if set(data['states']) != {str(p) for p in r.states}:
                raise ValueError('checkpoint point IDs differ')
        elif world.rank == 0:
            root.mkdir(parents=True, exist_ok=False)
    except Exception as exc: error = str(exc)
    errors = world.allgather(error)
    if any(errors): raise ValueError(str(errors))
    world.Barrier()
    path = root/f'rank-{world.rank}.json'
    data = data if restore else {
        'step': coupled.step_index, 'J': coupled.J.tolist(), 'response_A': r._A.tolist(), 'states': {}}

    def native_call(owner, symbol, stem):
        fn = getattr(owner.lib, symbol)
        fn.argtypes = [c.c_char_p, c.c_int]; fn.restype = c.c_int
        # Per-outer-rank filenames avoid collisions between independent solver groups.
        code = owner._call(fn, str(root/stem).encode(), int(restore))
        if code: raise RuntimeError(f'{symbol} failed with code {code}')

    try:
        if field is not None:
            native_call(field, 'sp_eh_checkpoint', f'field-{world.rank}.bin')
            if restore:
                field.A = np.asarray(data['field_A']); field.E = np.asarray(data['field_E'])
                field.step_index = data['step']
            else:
                data.update(field_A=field.A.tolist(), field_E=field.E.tolist())
        for point, state in r.states.items():
            key = str(point)
            if hasattr(state, 'native'):
                native_call(state, 'sp_tddft_checkpoint', f'electronic-{world.rank}.bin')
                if restore:
                    state.steps = data['step']; state.A = np.asarray(data['states'][key]['A'])
                else: data['states'][key] = {'A': state.A.tolist()}
            else:
                raw = state._raw
                token = raw.snapshot()
                if restore:
                    item = data['states'][key]
                    payload = tuple((n,tuple(shape),dtype,order,bytes.fromhex(hexdata))
                                    for n,shape,dtype,order,hexdata in item['data'])
                    raw.restore(replace(token, step_index=item['step'], time=item['time'],
                                        _settings=item['settings'], _data=payload, _integrity=item['integrity']))
                    state._steps = data['step']
                else:
                    data['states'][key] = {'step': token.step_index, 'time': token.time,
                        'settings': token._settings, 'integrity': token._integrity,
                        'data': [(n,shape,dtype,order,payload.hex()) for n,shape,dtype,order,payload in token._data]}
        if restore:
            coupled.step_index = r.step_index = data['step']
            coupled.J = np.asarray(data['J']); r._A = np.asarray(data['response_A'])
            for point, adapter in r._adapters.items():
                adapter.A = r._A[list(r.context.plan.point_ids).index(point)].copy()
        else: path.write_text(json.dumps(data, allow_nan=False)+'\n')
    except Exception as exc: error = str(exc)
    errors = world.allgather(error)
    if any(errors):
        r.failed = True
        raise RuntimeError(str(errors))
    world.Barrier()
    if not restore and world.rank == 0:
        (root/'frames.json').write_text(json.dumps([
            {k: v.tolist() if isinstance(v,np.ndarray) else v for k,v in h.items()}
            for h in frames], allow_nan=False)+'\n')
        files = {p.name: digest(p) for p in sorted(root.iterdir())}
        (root/'complete.tmp').write_text(json.dumps({'schema':'salmonpy.native-direct/1',
            'ranks':world.size,'signature':signature,'step':coupled.step_index,'files':files}, indent=2)+'\n')
        (root/'complete.tmp').replace(root/'complete.json')
    world.Barrier()
    return data['step']
