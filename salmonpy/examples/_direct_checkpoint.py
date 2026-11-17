"""Same-topology direct checkpoint coordinator for the native mixed fixture.

No pickle, no propagation during load. A completed manifest binds every rank file.
This internal example helper is not a general user-backend archive API.
"""
import ctypes as c
import hashlib
import json
from pathlib import Path
from dataclasses import replace
from functools import partial
from salmonpy.checkpointing import CheckpointTarget, execute_checkpoint
from salmonpy.errors import RuntimeFailure
import numpy as np


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''): h.update(block)
    return h.hexdigest()


def _verify_files(root, files, rank, ranks):
    """Hash each payload on one rank; caller agrees errors before any load."""
    if not isinstance(files, dict) or not files:
        raise ValueError('invalid checkpoint inventory')
    checked = 0
    size = 0
    for index, (name, sha) in enumerate(sorted(files.items())):
        if not isinstance(name, str) or Path(name).name != name or name in ('.','..','complete.json'):
            raise ValueError('invalid checkpoint filename')
        path = root/name
        if path.is_symlink() or not path.is_file():
            raise ValueError('invalid checkpoint payload')
        if not isinstance(sha, str) or len(sha) != 64 or any(c not in '0123456789abcdef' for c in sha):
            raise ValueError('checkpoint checksum mismatch')
        if index % ranks == rank:
            if digest(path) != sha:
                raise ValueError('checkpoint checksum mismatch')
            checked += 1
            size += path.stat().st_size
    if set(files) != {p.name for p in root.iterdir() if p.name != 'complete.json'}:
        raise ValueError('checkpoint file inventory mismatch')
    return {'rank': rank, 'files_checked': checked, 'bytes_checked': size}


def _compatibility(world, coupled, *, restore, saved=None, available_files=None, protected=False):
    """Solver-specific admission adapted onto the public metadata validator."""
    from dataclasses import asdict
    from salmonpy.checkpointing import CheckpointCapability, check_checkpoint_compatibility
    from salmonpy.sbe import SBEState
    from salmonpy.ld import LDState

    def component(state, *, symbol=None, native_format=None, ranks=1):
        if getattr(state,'closed',False) or getattr(state,'failed',False):
            raise ValueError('checkpoint capability requires a live state')
        if restore and not protected and getattr(state,'steps',getattr(state,'step_index',0))!=0:
            raise ValueError('checkpoint compatibility requires a fresh unadvanced state')
        if symbol is not None:
            if not callable(getattr(state.lib,symbol,None)):
                raise ValueError('missing native checkpoint capability '+symbol)
            return CheckpointCapability(native_format,('save','restore'),
                getattr(state,'checkpoint_decomposition',None))
        if getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
            if not callable(getattr(state,'checkpoint',None)):
                raise ValueError('missing direct checkpoint capability')
            declared=getattr(state,'checkpoint_capability',None)
            if declared is not None:
                return declared if isinstance(declared,CheckpointCapability) else CheckpointCapability.from_dict(declared)
            return CheckpointCapability(getattr(state,'checkpoint_format',None),
                getattr(state,'checkpoint_modes',None),getattr(state,'checkpoint_decomposition',None))
        raw=getattr(state,'_raw',None)
        if not isinstance(raw,(SBEState,LDState)):
            raise ValueError('backend has no declared direct checkpoint capability; snapshot may require replay')
        if not callable(getattr(raw,'snapshot',None)) or not callable(getattr(raw,'restore',None)):
            raise ValueError('snapshot/restore capability required')
        if restore and not protected and raw._steps!=0:
            raise ValueError('checkpoint compatibility requires a fresh snapshot state')
        if ranks!=1:
            raise ValueError('snapshot checkpoint capability requires one rank per group')
        return CheckpointCapability('salmonpy.'+raw._snapshot_kind+'-snapshot/1',('save','restore'),(1,))

    error=None;groups=None;group=None;field=None;states={};required=[]
    try:
        r=coupled.response;c=r.context
        if getattr(r,'failed',False):
            raise ValueError('checkpoint capability unavailable after session failure')
        if restore and not protected and (getattr(coupled,'step_index',0)!=0 or getattr(r,'step_index',0)!=0):
            raise ValueError('checkpoint compatibility requires a fresh unadvanced session')
        c.plan.validate(world.size)
        groups=[asdict(g) for g in c.plan.groups]
        group=c.group.name
        if world.rank not in c.group.ranks or set(r.states)!=set(c.group.point_ids):
            raise ValueError('checkpoint compatibility ownership differs from plan')
        ranks=len(c.group.ranks)
        if (coupled.field is not None)!=(c.group.role=='field'):
            raise ValueError('checkpoint compatibility field role differs')
        if coupled.field is not None:
            field=component(coupled.field,symbol='sp_eh_checkpoint',native_format='salmonpy.eh-native/1',ranks=ranks)
            required.append(f'field-{world.rank}.bin')
        for point,state in r.states.items():
            if hasattr(state,'native'):
                states[point]=component(state,symbol='sp_tddft_checkpoint',native_format='salmonpy.tddft-native/2',ranks=ranks)
                required.append(f'electronic-{world.rank}.bin')
            else:
                states[point]=component(state,ranks=ranks)
                if getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
                    required.append(f'worker-{world.rank}.bin')
    except Exception as exc:
        error=str(exc)
    return check_checkpoint_compatibility(world,groups,group=group,field=field,states=states,
        operation='restore' if restore else 'save',saved=saved,local_error=error,
        required_files=required if available_files is not None else (),available_files=available_files)


class _SnapshotCandidate:
    """Detached endpoint-model state; no original handle is restored in place."""
    def __init__(self, state, item, step, A, adapter, owner_plan):
        self.original = state
        self.item = item
        self.step = step
        self.A = A
        self.adapter = adapter
        self.owner_plan = owner_plan
        self.raw = None
        self.state = None

    def stage(self):
        from copy import copy
        from salmonpy.sbe import SBEState
        raw = self.original._raw
        if isinstance(raw, SBEState):
            self.raw = raw._initial.new_state(dt=raw.dt,time0=raw._time0)
        else:
            self.raw = raw._material.new_state(dt=raw.dt,initial_field=raw.config.initial_field)
        token = raw.snapshot()
        item = self.item
        payload = tuple((n,tuple(shape),dtype,order,bytes.fromhex(hexdata))
                        for n,shape,dtype,order,hexdata in item['data'])
        self.raw.restore(replace(token,step_index=item['step'],time=item['time'],
            _settings=item['settings'],_data=payload,_integrity=item['integrity']))
        if self.raw._steps != self.step:
            raise ValueError('snapshot checkpoint clock differs')
        if isinstance(raw,SBEState) and not np.array_equal(self.raw._A,self.A):
            raise ValueError('snapshot A differs from response history')
        self.state = copy(self.original)
        self.state._raw = self.raw
        self.state._steps = self.step
        self.new_adapter = copy(self.adapter)
        self.new_adapter.state = self.state
        self.new_adapter.A = self.A
        for index, old in enumerate(self.owner_plan):
            if old is raw:
                self.owner_plan[index] = self.raw
                break
        else:
            raise ValueError('snapshot state has no matching live owner')
        return self

    def discard(self):
        if self.raw is not None:
            self.raw.close()


def _restore_protected(world,coupled,root,data,compatibility,verification,history=None):
    """Stage all numerical and Python state before publishing any live owner.

    Native publication may allocate; fatal native/process/MPI failures are outside
    handled-exception recovery. There is no archive I/O after the ready vote.
    """
    from copy import copy
    from salmonpy._native_checkpoint import NativeRestoreCandidate
    from salmonpy.checkpointing import CheckpointState, restore_checkpoint
    from salmonpy.materials import _EndpointState
    from salmonpy.material_response import _StepAdapter
    from salmonpy.sbe import SBEState
    from salmonpy.ld import LDState
    r=coupled.response;field=coupled.field
    participant=compatibility['participants'][world.rank]
    error=None;field_target=None;targets={};native=[];snapshots=[];updates=[];owners={}
    try:
        step=data['step']
        J=np.array(data['J'],dtype=float,copy=True)
        response_A=np.array(data['response_A'],dtype=float,copy=True)
        states=dict(r.states);adapters=dict(r._adapters)
        if field is not None:
            candidate=NativeRestoreCandidate(field,'sp_eh_checkpoint',root/f'field-{world.rank}.bin')
            native.append(candidate)
            field_target=CheckpointTarget(participant['field'],stage_restore=candidate.stage,
                discard_restore=lambda unused,candidate=candidate:candidate.discard())
            updates.extend(((field,'A',np.array(data['field_A'],dtype=float)),
                            (field,'E',np.array(data['field_E'],dtype=float)),(field,'step_index',step)))
        for point,state in r.states.items():
            A=response_A[list(r.context.plan.point_ids).index(point)].copy()
            adapter=r._adapters[point]
            if type(adapter) is not _StepAdapter:
                raise ValueError('protected native sample requires endpoint response adapters')
            if hasattr(state,'native'):
                candidate=NativeRestoreCandidate(state,'sp_tddft_checkpoint',root/f'electronic-{world.rank}.bin')
                updates.extend(((state,'steps',step),(state,'A',A)))
                native.append(candidate)
            elif getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
                candidate=NativeRestoreCandidate(state,'sp_fast_hse_checkpoint',root/f'worker-{world.rank}.bin')
                updates.extend(((state,'steps',step),(state,'A',A),
                    (state,'current',np.array(data['states'][str(point)]['current'],dtype=float)),(state,'energy',None)))
                native.append(candidate)
            else:
                if not isinstance(state,_EndpointState) or not isinstance(state._raw,(SBEState,LDState)):
                    raise ValueError('protected snapshot requires an owned built-in endpoint state')
                prepared=state._prepared
                if id(prepared) not in owners:
                    owners[id(prepared)]=(prepared,list(prepared.live))
                candidate=_SnapshotCandidate(state,data['states'][str(point)],step,A,adapter,owners[id(prepared)][1])
                snapshots.append((point,candidate))
            targets[point]=CheckpointTarget(participant['states'][str(point)],stage_restore=candidate.stage,
                discard_restore=lambda unused,candidate=candidate:candidate.discard())
            replacement=copy(adapter);replacement.A=A
            adapters[point]=replacement
        holder=CheckpointState(field=field,states=r.states)
        validation=verification
    except Exception as exc:
        error=type(exc).__name__+': '+str(exc)
    # Admission also agrees local preparation/protocol errors; no old ABI fallback.
    restore_checkpoint(world,compatibility['groups'],group=r.context.group.name,live=holder if error is None else None,
        field=field_target,states=targets,saved=compatibility,local_error=error)
    # Snapshot-dependent bindings are constructed by the stage callbacks, before
    # their global ready vote. These assignments invoke no user/backend callbacks.
    for point,candidate in snapshots:
        states[point]=candidate.state
        adapters[point]=candidate.new_adapter
    publication_error=None
    try:
        for candidate in native:
            candidate.publish()
    except Exception as exc:
        publication_error=type(exc).__name__+': '+str(exc)
    failures=world.allgather(publication_error)
    if any(failures):
        r.failed=True
        raise RuntimeFailure('native publication failed after ready vote; original preservation is not guaranteed: '+str(failures))
    for owner,name,value in updates:
        setattr(owner,name,value)
    for prepared,planned in owners.values():
        prepared.live=planned
    r.states=states;r._adapters=adapters
    coupled.step_index=r.step_index=step
    coupled.J=J;r._A=response_A
    coupled.checkpoint_validation=validation
    coupled.checkpoint_frames=history
    cleanup=[]
    for point,candidate in snapshots:
        try:candidate.original.close()
        except Exception as exc:cleanup.append((str(point),type(exc).__name__+': '+str(exc)))
    coupled.checkpoint_cleanup_errors=world.allgather(cleanup)
    return step


def checkpoint(world, coupled, directory, signature, *, restore=False, frames=None, protected=True):
    r = coupled.response
    field = coupled.field
    root = Path(directory).resolve()
    error = None
    manifest_id = None
    try:
        if type(restore) is not bool or type(protected) is not bool:
            raise ValueError("restore and protected must be booleans")
        if sum(getattr(s,'checkpoint_protocol',None)=='salmonpy.direct-state/1' for s in r.states.values())>1:
            raise ValueError('one worker checkpoint per rank required')
        if sum(hasattr(s,'native') for s in r.states.values())>1:
            raise ValueError('one native electronic checkpoint per rank required')
        for state in r.states.values():
            if hasattr(state, 'native'):
                getattr(state.lib, 'sp_tddft_checkpoint')
            elif getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1' and callable(getattr(state,'checkpoint',None)):
                pass
            elif not callable(getattr(getattr(state, '_raw', None), 'snapshot', None)):
                raise ValueError('backend has no direct checkpoint capability')
        if field is not None: getattr(field.lib, 'sp_eh_checkpoint')
        if restore:
            manifest = json.loads((root/'complete.json').read_text())
            manifest_id = hashlib.sha256(json.dumps(manifest,sort_keys=True,allow_nan=False).encode()).hexdigest()
            if manifest['schema'] not in ('salmonpy.native-direct/1','salmonpy.native-direct/2') or manifest['signature'] != signature or manifest['ranks'] != world.size:
                raise ValueError('checkpoint input/library/topology mismatch')
            if type(manifest['step']) is not int or not 0 < manifest['step'] < signature['steps']:
                raise ValueError('invalid checkpoint step')
            if not ({f'rank-{i}.json' for i in range(world.size)} | {'frames.json'}) <= set(manifest['files']):
                raise ValueError('missing rank checkpoint')
            verification = _verify_files(root, manifest['files'], world.rank, world.size)
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
            for point,state in r.states.items():
                if hasattr(state,'native') or getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
                    A = np.asarray(data['states'][str(point)]['A'],dtype=float)
                    index = list(r.context.plan.point_ids).index(point)
                    if A.shape!=(3,) or not np.isfinite(A).all() or not np.array_equal(A,np.asarray(data['response_A'])[index]):
                        raise ValueError('checkpoint material A differs from response history')
                    if getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
                        current = np.asarray(data['states'][str(point)]['current'],dtype=float)
                        if current.shape!=(3,) or not np.isfinite(current).all():
                            raise ValueError('invalid checkpoint material current')
    except Exception as exc: error = str(exc)
    policy=(restore if type(restore) is bool else None,protected if type(protected) is bool else None)
    agreement = world.allgather((error, manifest_id, policy))
    errors = [item[0] for item in agreement]
    if any(errors): raise ValueError(str(errors))
    if any(item[2]!=policy for item in agreement):
        raise ValueError('checkpoint policies differ across ranks')
    if restore and len({item[1] for item in agreement}) != 1:
        raise ValueError('checkpoint manifests differ across ranks')
    compatibility=_compatibility(world,coupled,restore=restore,
        saved=manifest.get('compatibility') if restore else None,
        available_files=manifest['files'] if restore else None,protected=protected and restore)
    if not restore:
        error=None
        if world.rank==0:
            try: root.mkdir(parents=True,exist_ok=False)
            except Exception as exc: error=str(exc)
        errors=world.allgather(error)
        if any(errors): raise ValueError(str(errors))
    if restore and protected:
        from salmonpy._memory_advice import restore_advice
        restore_advice(world,coupled,root,data)
        error=None;history=None
        try:
            history=json.loads((root/'frames.json').read_text())
            if not isinstance(history,list) or len(history)!=data['step']:
                raise ValueError('checkpoint frame count mismatch')
            for frame in history:
                for name in ('A_over_c','E_at_points','J_charge','J_previous'):
                    values=np.asarray(frame[name],dtype=float)
                    if values.shape!=(count,3) or not np.isfinite(values).all():
                        raise ValueError('invalid checkpoint frame arrays')
                for name in ('field_time','electronic_local_time','electronic_field_clock','E_time'):
                    value=np.asarray(frame[name],dtype=float)
                    if value.shape!=() or not np.isfinite(value):
                        raise ValueError('invalid checkpoint frame clock')
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        errors=world.allgather(error)
        if any(errors):raise ValueError('checkpoint history admission failed: '+str(errors))
        verification=world.allgather(verification)
        return _restore_protected(world,coupled,root,data,compatibility,verification,history)
    if restore:
        coupled.checkpoint_validation = world.allgather(verification)
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

    def field_io():
        native_call(field, 'sp_eh_checkpoint', f'field-{world.rank}.bin')
        if restore:
            field.A = np.asarray(data['field_A']); field.E = np.asarray(data['field_E'])
            field.step_index = data['step']
        else:
            data.update(field_A=field.A.tolist(), field_E=field.E.tolist())

    def state_io(point,state):
        key = str(point)
        if hasattr(state, 'native'):
            native_call(state, 'sp_tddft_checkpoint', f'electronic-{world.rank}.bin')
            if restore:
                state.steps = data['step']; state.A = np.asarray(data['states'][key]['A'])
            else: data['states'][key] = {'A': state.A.tolist()}
        elif getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
            state.checkpoint(root/f'worker-{world.rank}.bin', restore=restore)
            if restore:
                if state.steps != data['step']: raise ValueError('worker checkpoint clock differs')
                state.A = np.asarray(data['states'][key]['A'])
                np.testing.assert_array_equal(state.current,np.asarray(data['states'][key]['current']))
            else: data['states'][key] = {'A':state.A.tolist(),'current':state.current.tolist()}
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

    participant=compatibility['participants'][world.rank]
    field_target=(CheckpointTarget(participant['field'],save=field_io,restore=field_io)
                  if field is not None else None)
    targets={point:CheckpointTarget(participant['states'][str(point)],
        save=partial(state_io,point,state),restore=partial(state_io,point,state))
        for point,state in r.states.items()}
    try:
        run=execute_checkpoint(world,compatibility['groups'],group=r.context.group.name,
            field=field_target,states=targets,operation='restore' if restore else 'save',atomic=False,
            saved=manifest.get('compatibility') if restore else None)
        compatibility=run['compatibility']
    except RuntimeFailure:
        r.failed=True
        raise
    try:
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
    error=None
    try:
        if not restore and world.rank == 0:
            (root/'frames.json').write_text(json.dumps([
                {k: v.tolist() if isinstance(v,np.ndarray) else v for k,v in h.items()}
                for h in frames], allow_nan=False)+'\n')
            files = {p.name: digest(p) for p in sorted(root.iterdir())}
            (root/'complete.tmp').write_text(json.dumps({'schema':'salmonpy.native-direct/2','compatibility':compatibility,
                'ranks':world.size,'signature':signature,'step':coupled.step_index,'files':files}, indent=2)+'\n')
            (root/'complete.tmp').replace(root/'complete.json')
    except Exception as exc:
        error=type(exc).__name__+': '+str(exc)
    errors=world.allgather(error)
    if any(e is not None for e in errors):
        r.failed=True
        raise RuntimeFailure('checkpoint completion publication failed: '+str(errors))
    world.Barrier()
    return data['step']
