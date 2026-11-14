"""Versioned, immutable, process-local snapshots of accepted states.

Tokens are created by state.snapshot(), not constructed or edited by callers.
They retain the originating model. They are not a disk checkpoint format.
"""
from dataclasses import dataclass, field
import hashlib
import json
import math
import numpy as np
from ._version import __version__
from .errors import SnapshotError

__all__ = ['StateSnapshot']
_SCHEMA = 'salmonpy.snapshot/1'


@dataclass(frozen=True, slots=True, eq=False)
class StateSnapshot:
    """Opaque immutable token; public metadata describes its accepted time."""
    kind: str
    step_index: int
    time: float
    schema: str
    version: str
    _owner: object = field(repr=False)
    _settings: str = field(repr=False)
    _data: tuple = field(repr=False)
    _integrity: str = field(repr=False)


def _fingerprint(settings, arrays=()):
    digest = hashlib.sha256(repr(settings).encode('utf-8'))
    for value in arrays:
        a = np.asarray(value)
        digest.update(repr((a.shape, a.dtype.str)).encode('ascii'))
        digest.update(a.tobytes(order='C'))
    return digest.hexdigest()


def _integrity(kind, step, time, settings, data):
    digest = hashlib.sha256(json.dumps([kind, step, time, settings], allow_nan=False).encode('utf-8'))
    for name, shape, dtype, order, payload in data:
        digest.update(repr((name, shape, dtype, order)).encode('utf-8'))
        digest.update(payload)
    return digest.hexdigest()


class _SnapshotMixin:
    __slots__ = ()

    def snapshot(self):
        """Copy accepted history without advancing; token retains this model."""
        data = []
        for name, a in self._snapshot_arrays().items():
            if not np.isfinite(a).all():
                raise SnapshotError('cannot snapshot a nonfinite state')
            order = 'F' if a.flags.f_contiguous and not a.flags.c_contiguous else 'C'
            data.append((name, a.shape, a.dtype.str, order, a.tobytes(order=order)))
        data = tuple(data)
        step = int(getattr(self, self._snapshot_counter))
        time = self.time
        if not math.isfinite(time):
            raise SnapshotError('cannot snapshot a nonfinite time')
        settings = self._snapshot_settings()
        return StateSnapshot(kind=self._snapshot_kind, step_index=step, time=time,
                             schema=_SCHEMA, version=__version__,
                             _owner=self._snapshot_owner(), _settings=settings, _data=data,
                             _integrity=_integrity(self._snapshot_kind, step, time, settings, data))

    def restore(self, snapshot):
        """Restore atomically into this state; incompatible tokens are rejected.

        A compatible state uses the same model object, dt, time origin and field
        capacity/drive. No native solver is invoked and no model is modified.
        """
        if type(snapshot) is not StateSnapshot:
            raise SnapshotError('expected StateSnapshot from state.snapshot()')
        if any(type(getattr(snapshot, name)) is not str
               for name in ('schema', 'version', 'kind', '_settings', '_integrity')):
            raise SnapshotError('invalid snapshot metadata')
        if (snapshot.schema != _SCHEMA or snapshot.version != __version__
                or snapshot.kind != self._snapshot_kind
                or snapshot._owner is not self._snapshot_owner()
                or snapshot._settings != self._snapshot_settings()):
            raise SnapshotError('snapshot requires the same version, model, dt, time origin and field capacity/drive')
        n = snapshot.step_index
        if type(n) is not int or not 0 <= n <= getattr(self, '_snapshot_capacity', 2**63-1):
            raise SnapshotError('invalid snapshot step index')
        origin = getattr(self, '_time0', 0.)
        if type(snapshot.time) not in (int, float) or not math.isfinite(snapshot.time) or snapshot.time != origin+n*self.dt:
            raise SnapshotError('invalid snapshot clock')
        expected = self._snapshot_arrays()
        restored = {}
        try:
            if type(snapshot._data) is not tuple or len(snapshot._data) != len(expected):
                raise SnapshotError('invalid snapshot history')
            if snapshot._integrity != _integrity(snapshot.kind, n, snapshot.time, snapshot._settings, snapshot._data):
                raise SnapshotError('snapshot integrity mismatch')
            for name, shape, dtype, order, payload in snapshot._data:
                if name not in expected or name in restored:
                    raise SnapshotError('invalid snapshot history names')
                target = expected[name]
                if shape != target.shape or dtype != target.dtype.str or order not in ('C', 'F') or type(payload) is not bytes:
                    raise SnapshotError('snapshot array layout mismatch')
                # Every allocation and validity check precedes state mutation.
                array = np.frombuffer(payload, dtype=target.dtype).reshape(shape, order=order).copy(order=order)
                if not np.isfinite(array).all():
                    raise SnapshotError('nonfinite snapshot history')
                restored[name] = array
        except (ValueError, TypeError, OverflowError) as exc:
            raise SnapshotError('invalid snapshot payload') from exc
        for name, array in restored.items():
            setattr(self, name, array)
        setattr(self, self._snapshot_counter, n)
