"""Checkpoint admission and explicit backend-owned execution.

Metadata checks perform no I/O or callbacks. Execution calls only explicitly
bound backend functions. Hosts own MPI, storage, authentication and lifetimes.
"""
from collections.abc import Mapping
from dataclasses import dataclass
import math
import json
from pathlib import PurePosixPath
from types import MappingProxyType
import threading

from .config import _copy_json
from .errors import InputError, RuntimeFailure

__all__ = ['CheckpointCapability', 'check_checkpoint_compatibility',
           'CheckpointTarget', 'execute_checkpoint', 'CheckpointState',
           'restore_checkpoint', 'AtomicRestoreError']


def _text(value, name, *, empty=False):
    if not isinstance(value, str) or (not empty and not value.strip()):
        raise InputError(name + ': string required')
    try:
        value.encode('utf-8')
    except UnicodeError as exc:
        raise InputError(name + ': valid Unicode required') from exc
    return value


@dataclass(frozen=True)
class CheckpointCapability:
    """Immutable declaration, not a certification of a backend implementation.

    format is a versioned caller-defined identifier. decomposition describes
    the actual subgroup layout; its product must match that group's rank count.
    Optional identity can pin a model/input/runtime digest chosen by the caller.
    """
    format: str
    modes: tuple
    decomposition: tuple
    identity: str = ''

    def __post_init__(self):
        _text(self.format, 'checkpoint format')
        _text(self.identity, 'checkpoint identity', empty=True)
        if (not isinstance(self.modes, tuple) or not self.modes
                or any(type(m) is not str or m not in ('save', 'restore') for m in self.modes)
                or len(set(self.modes)) != len(self.modes)):
            raise InputError('checkpoint capability requires unique save/restore modes')
        if (not isinstance(self.decomposition, tuple) or not self.decomposition
                or any(type(n) is not int or n < 1 for n in self.decomposition)):
            raise InputError('checkpoint capability requires a positive integer decomposition tuple')
        object.__setattr__(self, 'modes', tuple(sorted(self.modes)))

    def to_dict(self):
        """Detached JSON values; empty identity is omitted for record compatibility."""
        result = {'format': self.format, 'modes': list(self.modes),
                  'decomposition': list(self.decomposition)}
        if self.identity:
            result['identity'] = self.identity
        return result

    @classmethod
    def from_dict(cls, value):
        """Validate explicit JSON/mapping declarations; never infer capabilities."""
        required = {'format', 'modes', 'decomposition'}
        if (not isinstance(value, Mapping) or not required <= set(value)
                or set(value) - required - {'identity'}):
            raise InputError('explicit checkpoint capability format/modes/decomposition required')
        if (not isinstance(value['modes'], (list, tuple))
                or not isinstance(value['decomposition'], (list, tuple))):
            raise InputError('checkpoint capability modes/decomposition must be sequences')
        return cls(value['format'], tuple(value['modes']), tuple(value['decomposition']),
                   value.get('identity', ''))


def _sequence(value, name):
    if not isinstance(value, (list, tuple)):
        raise InputError(name + ': list or tuple required')
    return value


def _point_key(value):
    if type(value) is int and value >= 0:
        return str(value)
    if type(value) is str:
        return _text(value, 'point ID')
    raise InputError('point IDs must be nonnegative integers or nonempty strings')


def _groups(groups, size):
    result = []
    names = set()
    ranks_seen = []
    points_seen = set()
    for item in _sequence(groups, 'checkpoint groups'):
        if (not isinstance(item, Mapping) or not {'name', 'ranks'} <= set(item)
                or set(item) - {'name', 'ranks', 'role', 'point_ids'}):
            raise InputError('checkpoint group requires name/ranks and optional role/point_ids')
        name = _text(item['name'], 'group name')
        ranks = list(_sequence(item['ranks'], 'group ranks'))
        points = list(_sequence(item.get('point_ids', ()), 'group point IDs'))
        role = item.get('role', 'material')
        if (name in names or role not in ('field', 'material') or not ranks
                or any(type(n) is not int or n < 0 for n in ranks)):
            raise InputError('invalid checkpoint group name/role/ranks')
        keys = [_point_key(p) for p in points]
        if len(set(keys)) != len(keys) or points_seen.intersection(keys):
            raise InputError('checkpoint points require unique owners and string representations')
        if role == 'field' and points:
            raise InputError('field group cannot own material point IDs')
        if role == 'material' and not points:
            raise InputError('material group requires point IDs')
        names.add(name)
        ranks_seen.extend(ranks)
        points_seen.update(keys)
        result.append({'name': name, 'ranks': ranks, 'point_ids': points, 'role': role})
    if not result or sorted(ranks_seen) != list(range(size)):
        raise InputError('checkpoint groups must partition communicator ranks exactly once')
    return result


def _files(value):
    if isinstance(value, Mapping):
        value = list(value)
    if not isinstance(value, (list, tuple, set, frozenset)):
        raise InputError('checkpoint file inventory must be a collection of relative names')
    names = []
    for name in value:
        _text(name, 'checkpoint filename')
        path = PurePosixPath(name)
        if (path.is_absolute() or str(path) != name or name in ('.', '..')
                or '..' in path.parts or any(c in name for c in ('\\', ':', '\x00', '\n', '\r'))):
            raise InputError('checkpoint filename must be a safe relative POSIX path')
        names.append(name)
    if len(set(names)) != len(names):
        raise InputError('duplicate required checkpoint files')
    return set(names)


def check_checkpoint_compatibility(comm, groups, *, group=None, field=None, states=None,
                                   operation='save', saved=None, required_files=(),
                                   available_files=None, local_error=None):
    """Collectively admit metadata, returning a detached direct-compatibility/1 record.

    comm is a live host-owned communicator with rank/size/allgather, or None for
    serial use. groups are ordinary mappings; no SALMON context is required.
    states maps this group's point IDs to capabilities or explicit mappings.
    Mapping validation occurs inside the collective; local_error propagates an
    earlier caller-side error. All ranks must request the same operation/saved
    record. Required-file checks inspect names only, never contents or hashes.
    Host code may perform SAVE/LOAD only after this call returns successfully.
    """
    rank, size = (0, 1) if comm is None else (comm.rank, comm.size)
    if type(rank) is not int or type(size) is not int or not 0 <= rank < size:
        raise InputError('live checkpoint communicator rank/size required')
    gather = (lambda value: [value]) if comm is None else comm.allgather
    error = topology = participant = request = saved_record = None
    try:
        if local_error is not None:
            raise InputError(_text(local_error, 'local checkpoint error'))
        if type(operation) is not str or operation not in ('save', 'restore'):
            raise InputError('checkpoint operation must be save or restore')
        request = operation
        if operation == 'save' and saved is not None:
            raise InputError('saved compatibility record is only accepted for restore')
        if saved is not None:
            saved_record = _copy_json(saved)
        topology = _groups(groups, size)
        owner = next(g for g in topology if rank in g['ranks'])
        if group is not None and group != owner['name']:
            raise InputError('checkpoint compatibility group ownership differs')
        if (field is not None) != (owner['role'] == 'field'):
            raise InputError('checkpoint capability required exactly on field ranks')
        material_states = {} if states is None else states
        if not isinstance(material_states, Mapping):
            raise InputError('checkpoint states must map point IDs to capabilities')
        keys = [_point_key(p) for p in material_states]
        if len(set(keys)) != len(keys) or set(keys) != {_point_key(p) for p in owner['point_ids']}:
            raise InputError('checkpoint capability required for every owned point')

        def declaration(value):
            cap = value if isinstance(value, CheckpointCapability) else CheckpointCapability.from_dict(value)
            if operation not in cap.modes:
                raise InputError('checkpoint capability does not support ' + operation)
            if math.prod(cap.decomposition) != len(owner['ranks']):
                raise InputError('checkpoint compatibility decomposition does not match group ranks')
            return cap.to_dict()

        participant = {'rank': rank, 'group': owner['name'],
                       'field': declaration(field) if field is not None else None,
                       'states': {_point_key(p): declaration(v) for p, v in material_states.items()}}
        required = _files(required_files)
        inventory = None if available_files is None else _files(available_files)
        if operation == 'restore' and required and (inventory is None or not required <= inventory):
            missing = required if inventory is None else required - inventory
            raise InputError('checkpoint compatibility missing required payload inventory: ' + str(sorted(missing)))
    except Exception as exc:
        error = type(exc).__name__ + ': ' + str(exc)
    gathered = gather((error, topology, participant, request, saved_record))
    if len(gathered) != size:
        raise InputError('checkpoint compatibility communicator inventory differs')
    if any(row[0] for row in gathered):
        raise InputError('checkpoint capability/compatibility: ' + str([row[0] for row in gathered]))
    if any(row[3] != request for row in gathered):
        raise InputError('checkpoint compatibility operation differs across ranks')
    if any(row[1] != topology for row in gathered):
        raise InputError('checkpoint compatibility plans differ across ranks')
    canonical = lambda value: json.dumps(value, sort_keys=True, allow_nan=False)
    if any(canonical(row[4]) != canonical(saved_record) for row in gathered):
        raise InputError('saved checkpoint compatibility records differ across ranks')
    participants = [row[2] for row in gathered]
    if [p['rank'] for p in participants] != list(range(size)):
        raise InputError('checkpoint compatibility rank inventory differs')
    for owner in topology:
        members = [participants[n] for n in owner['ranks']]
        declarations = [(p['field'], p['states']) for p in members]
        if any(p['group'] != owner['name'] for p in members) or any(d != declarations[0] for d in declarations):
            raise InputError('checkpoint capability/compatibility declarations differ within MPI group ' + owner['name'])
    current = {'schema': 'salmonpy.direct-compatibility/1',
               'groups': topology, 'participants': participants}
    if operation == 'restore' and (saved_record is None or canonical(saved_record) != canonical(current)):
        raise InputError('checkpoint compatibility metadata missing or MPI layout/format/capability/identity differs; legacy archives require their original reader')
    return current


class AtomicRestoreError(RuntimeFailure):
    """A staged restore failed before publication; live holders are unchanged."""
    def __init__(self, message, *, cleanup_errors=()):
        super().__init__(message)
        self.original_preserved = True
        self.cleanup_errors = tuple(cleanup_errors)


class CheckpointState:
    """Process-local state holder, published by one internal bundle assignment.

    Access field/states through this holder after restore. Existing object aliases
    remain references to the old states. Pause computation during restoration;
    stage callbacks must never mutate these original objects or share mutable
    numerical resources with their candidates.
    """
    __slots__ = ('_current', '_lock')

    def __init__(self, *, field=None, states=None):
        values = {} if states is None else states
        if not isinstance(values, Mapping):
            raise InputError('CheckpointState states must be a point mapping')
        keys = [_point_key(p) for p in values]
        if len(set(keys)) != len(keys):
            raise InputError('CheckpointState point representations must be distinct')
        self._current = (field, MappingProxyType({_point_key(p): s for p, s in values.items()}))
        self._lock = threading.Lock()

    @property
    def field(self):
        return self._current[0]

    @property
    def states(self):
        return self._current[1]


@dataclass(frozen=True)
class CheckpointTarget:
    """Bind an explicit capability to backend-owned zero-argument I/O functions.

    Validation occurs inside execute_checkpoint's collective. Callback results
    remain local and may be opaque tokens, arrays, receipts or None. Bind paths
    using closures or functools.partial; the common API prescribes no payload.
    """
    capability: object
    save: object = None
    restore: object = None
    stage_restore: object = None
    discard_restore: object = None


def execute_checkpoint(comm, groups, *, group=None, field=None, states=None,
                       operation='save', saved=None, required_files=(),
                       available_files=None, local_error=None, atomic=False, live=None):
    """Admit then execute backend callbacks; return compatibility/field/states.

    Arguments match check_checkpoint_compatibility, using CheckpointTarget
    bindings for field and states. All ranks agree callback availability before
    any I/O. One group-local component runs per slot, followed by a global error
    agreement; unequal point counts are allowed. Collective callbacks
    must complete their group protocol before raising; this function cannot
    unblock peers stranded inside a backend collective. RuntimeFailure
    stops later slots but does not roll back completed callbacks. Opaque return
    values are never included in MPI messages. Preflight raises InputError.
    """
    declarations = {}
    field_declaration = None
    callbacks = {}
    field_callback = None
    discarders = {}
    field_discard = None
    error = local_error
    try:
        if error is None:
            if type(operation) is not str or operation not in ('save', 'restore'):
                raise InputError('checkpoint operation must be save or restore')

            if type(atomic) is not bool:
                raise InputError('atomic must be boolean')
            if atomic and operation != 'restore':
                raise InputError('atomic execution requires operation=restore')
            if atomic and type(live) is not CheckpointState:
                raise InputError('atomic restore requires CheckpointState')
            if not atomic and live is not None:
                raise InputError('use restore_checkpoint or atomic=True with a live holder')

            def binding(value):
                if not isinstance(value, CheckpointTarget):
                    raise InputError('CheckpointTarget required for checkpoint execution')
                callback = value.stage_restore if atomic else getattr(value, operation)
                if not callable(callback):
                    raise InputError('checkpoint ' + ('stage_restore' if atomic else operation) + ' callback required; no destructive fallback')
                if atomic and value.discard_restore is not None and not callable(value.discard_restore):
                    raise InputError('discard_restore must be callable or None')
                return value.capability, callback

            if field is not None:
                field_declaration, field_callback = binding(field)
                field_discard = field.discard_restore
            targets = {} if states is None else states
            if not isinstance(targets, Mapping):
                raise InputError('checkpoint states must map point IDs to targets')
            for point, target in targets.items():
                declarations[point], callbacks[point] = binding(target)
                discarders[_point_key(point)] = target.discard_restore
            if atomic:
                if set(live.states) != {_point_key(p) for p in targets}:
                    raise InputError('live holder must cover exactly the owned points')
                if field is None and live.field is not None:
                    raise InputError('live holder has an unbound field state')
    except Exception as exc:
        error = type(exc).__name__ + ': ' + str(exc)
    record = check_checkpoint_compatibility(comm, groups, group=group,
        field=field_declaration, states=declarations, operation=operation,
        saved=saved, required_files=required_files, available_files=available_files,
        local_error=error)
    rank = 0 if comm is None else comm.rank
    gather = (lambda value: [value]) if comm is None else comm.allgather
    owner = next(g for g in record['groups'] if rank in g['ranks'])
    protocol = {'field': field_discard is not None,
                'states': {k: callback is not None for k, callback in discarders.items()}} if atomic else None
    policies = gather((atomic, protocol))
    if any(policy[0] != atomic for policy in policies):
        raise InputError('atomic restore policy differs across ranks')
    if atomic:
        for item in record['groups']:
            if any(policies[r][1] != policies[item['ranks'][0]][1] for r in item['ranks']):
                raise InputError('discard_restore availability differs within an MPI group')
    canonical_callbacks = {_point_key(p): callback for p, callback in callbacks.items()}
    if owner['role'] == 'field':
        schedule = [('field', field_callback)]
    else:
        schedule = [(str(p), canonical_callbacks[str(p)]) for p in owner['point_ids']]
    slots = max(1 if g['role'] == 'field' else len(g['point_ids']) for g in record['groups'])
    if atomic:
        cleanup = [field_discard] if owner['role'] == 'field' else [discarders[str(p)] for p in owner['point_ids']]
        return _restore_candidates(gather, rank, record, owner, schedule, cleanup, live, slots)
    result = {'compatibility': record, 'field': None, 'states': {}}
    for slot in range(slots):
        error = None
        cause = None
        component = None
        if slot < len(schedule):
            component, callback = schedule[slot]
            try:
                value = callback()
                if owner['role'] == 'field':
                    result['field'] = value
                else:
                    result['states'][component] = value
            except Exception as exc:
                cause = exc
                error = type(exc).__name__ + ': ' + str(exc)
        statuses = gather((rank, component, error))
        failures = [status for status in statuses if status[2] is not None]
        if failures:
            raise RuntimeFailure('checkpoint ' + operation + ' callback failure at slot '
                + str(slot) + ': ' + str(failures)
                + '; later callbacks stopped; completed state changes are not rolled back') from cause
    return result


def _restore_candidates(gather, rank, record, owner, schedule, cleanup, live, slots):
    reserved = live._lock.acquire(blocking=False)
    errors = gather(None if reserved else 'CheckpointState is already in use')
    if any(error is not None for error in errors):
        if reserved:
            live._lock.release()
        raise InputError('atomic restore reservation failed: ' + str(errors))
    initialization_error = None
    try:
        original = live._current
        candidates = [None] * len(schedule)
        candidate_ids = set()
        live_ids = {id(s) for s in (original[0], *original[1].values()) if s is not None}
        result = {'compatibility': record, 'field': None, 'states': {},
                  'previous': {'field': original[0], 'states': original[1]}, 'atomic': True}
    except Exception as exc:
        initialization_error = type(exc).__name__ + ': ' + str(exc)
    errors = gather(initialization_error)
    if any(error is not None for error in errors):
        live._lock.release()
        raise AtomicRestoreError('atomic restore preparation failed: ' + str(errors)
            + '; original states preserved')

    def discard(attempted):
        failures = []
        for slot in range(attempted, -1, -1):
            error = None
            component = None
            if slot < len(schedule):
                component = schedule[slot][0]
                callback = cleanup[slot]
                if callback is not None:
                    try:
                        callback(candidates[slot])
                    except Exception as exc:
                        error = type(exc).__name__ + ': ' + str(exc)
            failures.extend(row for row in gather((rank, component, error)) if row[2] is not None)
        return failures

    published = False
    try:
        for slot in range(slots):
            error = None
            cause = None
            component = None
            if slot < len(schedule):
                component, callback = schedule[slot]
                try:
                    candidate = callback()
                    if candidate is None or id(candidate) in live_ids or id(candidate) in candidate_ids:
                        raise InputError('stage_restore must return an independent distinct candidate')
                    candidates[slot] = candidate
                    candidate_ids.add(id(candidate))
                    if owner['role'] == 'field':
                        result['field'] = candidate
                    else:
                        result['states'][component] = candidate
                except Exception as exc:
                    cause = exc
                    error = type(exc).__name__ + ': ' + str(exc)
            failures = [row for row in gather((rank, component, error)) if row[2] is not None]
            if failures:
                cleanup_errors = discard(slot)
                raise AtomicRestoreError('atomic restore staging failed: ' + str(failures)
                    + '; original states preserved', cleanup_errors=cleanup_errors) from cause
        error = None
        next_bundle = None
        try:
            if live._current is not original:
                raise InputError('live holder changed during candidate staging')
            next_bundle = (result['field'], MappingProxyType(dict(result['states'])))
        except Exception as exc:
            error = type(exc).__name__ + ': ' + str(exc)
        errors = gather(error)
        if any(item is not None for item in errors):
            cleanup_errors = discard(slots-1)
            raise AtomicRestoreError('atomic restore publication vote failed: ' + str(errors)
                + '; original states preserved', cleanup_errors=cleanup_errors)
        # Exact CheckpointState type, one existing slot, no user callback here.
        # result['previous'] keeps old resources alive throughout publication.
        object.__setattr__(live, '_current', next_bundle)
        published = True
        gather(None)
        return result
    finally:
        if not published:
            object.__setattr__(live, '_current', original)
        live._lock.release()


def restore_checkpoint(comm, groups, *, live, group=None, field=None, states=None,
                       saved=None, required_files=(), available_files=None, local_error=None):
    """Protected restore: stage independent candidates and publish after all succeed.

    Direct-only targets are rejected before callbacks. On handled staging errors,
    AtomicRestoreError reports preserved original holders and candidate-cleanup
    errors. Fatal process/communicator failure and an isolation-violating backend
    are outside this contract. Retire previous states explicitly after success.
    """
    return execute_checkpoint(comm, groups, live=live, atomic=True, operation='restore',
        group=group, field=field, states=states, saved=saved, required_files=required_files,
        available_files=available_files, local_error=local_error)
