"""Declarative checkpoint admission for serial and host-owned MPI programs.

No solver imports, MPI initialization, filesystem access or persistence callbacks.
The caller owns actual I/O, file authentication and state lifetime validation.
"""
from collections.abc import Mapping
from dataclasses import dataclass
import math
import json
from pathlib import PurePosixPath

from .config import _copy_json
from .errors import InputError

__all__ = ['CheckpointCapability', 'check_checkpoint_compatibility']


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
