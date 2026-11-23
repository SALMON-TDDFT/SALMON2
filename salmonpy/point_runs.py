"""Explicit serial point bindings and identity-bearing numerical archives."""
from collections.abc import Mapping
from dataclasses import dataclass
import re
from types import MappingProxyType
from uuid import uuid4
import numpy as np

from .config import _identifier, _copy_json, _freeze
from .errors import InputError, SnapshotError
from ._version import __version__
from .outputs import OutputArchive, save_arrays

__all__ = ['PointRun', 'PointState', 'open_point_arrays']


def _identity(value):
    result = _copy_json(value)
    required = {'material_id', 'definition_digest'}
    optional = {'initial_state_id', 'initial_state_digest'}
    if type(result) is not dict or set(result) not in (required, required | optional):
        raise InputError('identity requires material_id/definition_digest and optionally both initial-state fields')
    _identifier(result['material_id'], 'material_id')
    if 'initial_state_id' in result:
        _identifier(result['initial_state_id'], 'initial_state_id')
    for key in ('definition_digest', 'initial_state_digest'):
        if key in result and (type(result[key]) is not str or not re.fullmatch('[0-9a-f]{64}', result[key])):
            raise InputError(f'{key} must be lowercase SHA256')
    return result


def _record(value):
    record = _copy_json(value)
    if type(record) is not dict or set(record) != {'schema', 'run_id', 'points'} or record['schema'] != 'salmonpy.point-run/1':
        raise InputError('invalid point-run record')
    _identifier(record['run_id'], 'run_id')
    if type(record['points']) is not list or not record['points']:
        raise InputError('point-run requires a nonempty ordered point list')
    points, materials, initials = set(), {}, {}
    for entry in record['points']:
        if type(entry) is not dict or 'point_id' not in entry:
            raise InputError('invalid point entry')
        point = _identifier(entry['point_id'], 'point_id')
        if point in points:
            raise InputError('duplicate point_id')
        points.add(point)
        identity = _identity({key: val for key, val in entry.items() if key != 'point_id'})
        material, digest = identity['material_id'], identity['definition_digest']
        if materials.setdefault(material, digest) != digest:
            raise InputError('conflicting definition digests for material_id')
        if 'initial_state_id' in identity:
            initial = identity['initial_state_id']
            signature = (material, identity['initial_state_digest'])
            if initials.setdefault(initial, signature) != signature:
                raise InputError('conflicting initial_state_id')
    return record


def _axes(axes, shapes, count):
    if not isinstance(axes, Mapping) or not axes:
        raise InputError('explicit nonempty point_axes mapping required')
    result = dict(axes)
    for name, axis in result.items():
        if name not in shapes or type(axis) is not int or axis < 0 or axis >= len(shapes[name]):
            raise InputError('invalid point axis or unknown array')
        if shapes[name][axis] != count:
            raise InputError('point axis length differs from ordered point mapping')
    return result


@dataclass(frozen=True, slots=True, eq=False)
class _PointSnapshot:
    """Process-local identity envelope; the backend owns the opaque payload."""
    _identity: tuple
    _payload: object
    schema: str = 'salmonpy.point-snapshot/1'
    version: str = __version__

    def record(self):
        return dict(self._identity)

    @property
    def step_index(self):
        return self._payload.step_index

    @property
    def time(self):
        return self._payload.time


@dataclass(frozen=True, init=False)
class PointState:
    """One serial state; its prepared model is borrowed, not owned."""
    _state: object
    _identity: object
    _closed: bool

    def __init__(self):
        raise InputError('use PointRun.new_state')

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        if name in self._identity:
            return self._identity[name]
        if self._closed:
            raise InputError('point state is closed')
        return getattr(self._state, name)

    def record(self):
        return dict(self._identity)

    def snapshot(self):
        """Save a process-local backend token with this exact point/run identity."""
        if self._closed:
            raise SnapshotError('point state is closed')
        callback = getattr(self._state, 'snapshot', None)
        if not callable(callback):
            raise SnapshotError('point backend does not support snapshots')
        return _PointSnapshot(tuple(sorted(self._identity.items())), callback())

    def restore(self, snapshot):
        """Reject mismatched identity before the backend's own restore checks.

        Backend snapshot semantics, settings checks and failure atomicity remain
        its responsibility. This is not a branch, clone or MPI transaction API.
        """
        if self._closed:
            raise SnapshotError('point state is closed')
        if (type(snapshot) is not _PointSnapshot
                or snapshot.schema != 'salmonpy.point-snapshot/1'
                or snapshot.version != __version__
                or snapshot._identity != tuple(sorted(self._identity.items()))):
            raise SnapshotError('point snapshot identity/version mismatch')
        callback = getattr(self._state, 'restore', None)
        if not callable(callback):
            raise SnapshotError('point backend does not support restore')
        return callback(snapshot._payload)

    def save_checkpoint(self, directory):
        """Save known TDDFT replay state with exact point/run identity; new path only."""
        from ._point_checkpoint import save
        save(self, directory)

    def restore_checkpoint(self, directory, *, max_bytes=None):
        """Verify exact point/run/material/initial identity before TDDFT replay.

        Supply the same run ID for an explicit resume. This is not a branch API
        or MPI whole-session checkpoint. max_bytes applies to the inner archive.
        """
        from ._point_checkpoint import restore
        restore(self, directory, max_bytes=max_bytes)

    def close(self):
        if not self._closed:
            object.__setattr__(self, '_closed', True)
            close = getattr(self._state, 'close', None)
            if close is not None:
                close()

    def __enter__(self):
        if self._closed:
            raise InputError('point state is closed')
        return self

    def __exit__(self, *args):
        self.close()


class PointRun:
    """Ordered point_id -> model bindings. Models expose record/new_state.

    A run owns each created state once, but never closes borrowed models. Custom
    backends must return independent states; this wrapper cannot inspect hidden
    shared native storage. No MPI, step scheduling, checkpoint or unit conversion.
    """
    def __init__(self, points, *, run_id=None):
        if not isinstance(points, Mapping) or not points:
            raise InputError('points must be a nonempty ordered mapping')
        models = dict(points)
        entries = []
        for point, model in models.items():
            _identifier(point, 'point_id')
            if not callable(getattr(model, 'record', None)) or not callable(getattr(model, 'new_state', None)):
                raise InputError('point models require record() and new_state()')
            entries.append({'point_id': point, **_identity(model.record())})
        self._record = _freeze(_record({'schema': 'salmonpy.point-run/1',
                                       'run_id': str(uuid4()) if run_id is None else run_id,
                                       'points': entries}))
        self._models = MappingProxyType(models)
        self._states = {}
        self._closed = False

    @property
    def run_id(self):
        return self._record['run_id']

    @property
    def point_ids(self):
        return tuple(self._models)

    def record(self):
        return {'schema': self._record['schema'], 'run_id': self.run_id,
                'points': [dict(entry) for entry in self._record['points']]}

    def new_state(self, point_id, **kwargs):
        _identifier(point_id, 'point_id')
        if self._closed or point_id not in self._models or point_id in self._states:
            raise InputError('closed run, unknown point or point state already created')
        model = self._models[point_id]
        entry = next(entry for entry in self._record['points'] if entry['point_id'] == point_id)
        identity = {key: val for key, val in entry.items() if key != 'point_id'}
        if _identity(model.record()) != identity:
            raise InputError('model identity changed after point binding')
        raw = model.new_state(**kwargs)
        if any(raw is state._state for state in self._states.values()):
            raise InputError('backend returned the same state object for multiple points')
        state = object.__new__(PointState)
        object.__setattr__(state, '_state', raw)
        object.__setattr__(state, '_identity', MappingProxyType({'run_id': self.run_id, **entry}))
        object.__setattr__(state, '_closed', False)
        self._states[point_id] = state
        return state

    def save_arrays(self, directory, arrays, *, point_axes, metadata=None):
        if not isinstance(arrays, Mapping):
            raise InputError('arrays must be a mapping')
        checked = {name: np.asarray(value) for name, value in arrays.items()}
        axes = _axes(point_axes, {name: value.shape for name, value in checked.items()}, len(self.point_ids))
        info = _copy_json({} if metadata is None else metadata)
        if type(info) is not dict or 'point_run' in info or 'point_axes' in info or 'mpi_point_map' in info:
            raise InputError('metadata must be an object without reserved point_run/point_axes')
        info.update(point_run=self.record(), point_axes=axes)
        return save_arrays(directory, checked, metadata=info)

    def close(self):
        if self._closed:
            return
        self._closed = True
        failures = []
        for state in self._states.values():
            try:
                state.close()
            except BaseException as exc:
                failures.append(exc)
        if failures:
            raise failures[0]

    def __enter__(self):
        if self._closed:
            raise InputError('point run is closed')
        return self

    def __exit__(self, *args):
        self.close()


def _mpi_mapping(mapping, record):
    if type(mapping) is not list or len(mapping) != len(record['points']):
        raise InputError('invalid MPI point mapping')
    transport, groups, owners = [], {}, {}
    for row, entry in zip(mapping, record['points']):
        if type(row) is not dict or set(row) != {'transport_point_id', 'point_id', 'group_name', 'ranks'}:
            raise InputError('invalid MPI point mapping entry')
        number = row['transport_point_id']
        if type(number) is not int or number < 0 or row['point_id'] != entry['point_id']:
            raise InputError('MPI transport/user point identity mismatch')
        transport.append(number)
        group, ranks = row['group_name'], row['ranks']
        if type(group) is not str or not group or type(ranks) is not list or not ranks:
            raise InputError('invalid MPI group ownership')
        if any(type(rank) is not int or rank < 0 for rank in ranks) or len(set(ranks)) != len(ranks):
            raise InputError('invalid MPI owner ranks')
        if groups.setdefault(group, ranks) != ranks:
            raise InputError('inconsistent MPI group ranks')
        for rank in ranks:
            if owners.setdefault(rank, group) != group:
                raise InputError('MPI rank assigned to different groups')
    if transport != sorted(set(transport)):
        raise InputError('MPI transport points must be unique in sorted transport order')


def open_point_arrays(directory):
    """Verify outputs/1 arrays and point metadata; does not restore solver state."""
    archive = OutputArchive.open(directory)
    metadata = archive.metadata
    if type(metadata) is not dict or 'point_run' not in metadata or 'point_axes' not in metadata:
        raise InputError('point-run metadata missing')
    record = _record(metadata['point_run'])
    if 'mpi_point_map' in metadata:
        _mpi_mapping(metadata['mpi_point_map'], record)
    _axes(metadata['point_axes'], {name: spec['shape'] for name, spec in archive.arrays.items()}, len(record['points']))
    return archive
