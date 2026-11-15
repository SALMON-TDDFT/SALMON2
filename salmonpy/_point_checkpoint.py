"""Identity envelope around the existing managed TDDFT replay checkpoint."""
import hashlib
import json
from pathlib import Path
import re
import shutil
from .errors import SnapshotError
from ._version import __version__

_SCHEMA = 'salmonpy.point-checkpoint/1'
_LIMIT = 65536
_MANIFEST_LIMIT = 16 * 1024 * 1024


def _native_state(point):
    from .builtin_materials import ConfiguredState
    from .tddft import TDDFTState
    if point._closed:
        raise SnapshotError('point state is closed')
    state = point._state
    if isinstance(state, ConfiguredState):
        state = state._state
    if not isinstance(state, TDDFTState):
        raise SnapshotError('point disk checkpoints currently require managed TDDFTState')
    return state


def _hash(path, limit=_MANIFEST_LIMIT):
    count = 0
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            count += len(block)
            if count > limit:
                raise SnapshotError('point checkpoint inner manifest is oversized')
            digest.update(block)
    return digest.hexdigest()


def save(point, directory):
    state = _native_state(point)
    root = Path(directory).expanduser()
    root.mkdir(parents=True, exist_ok=False)
    try:
        state.save_checkpoint(root/'state')
        record = {'schema': _SCHEMA, 'version': __version__, 'identity': point.record(),
                  'state_manifest_sha256': _hash(root/'state/manifest.json')}
        content = (json.dumps(record, sort_keys=True, allow_nan=False, ensure_ascii=True)+'\n').encode()
        if len(content) > _LIMIT:
            raise SnapshotError('point identity envelope exceeds 64 KiB')
        (root/'point.json').write_bytes(content)
    except BaseException:
        shutil.rmtree(root)
        raise


def restore(point, directory, *, max_bytes=None):
    state = _native_state(point)
    if max_bytes is not None and (type(max_bytes) is not int or max_bytes <= 0):
        raise SnapshotError('max_bytes must be None or a positive integer')
    root = Path(directory).expanduser()
    manifest = root/'point.json'
    payload = root/'state'
    if (root.is_symlink() or not root.is_dir() or
            {p.name for p in root.iterdir()} != {'point.json','state'} or
            manifest.is_symlink() or not manifest.is_file() or
            payload.is_symlink() or not payload.is_dir()):
        raise SnapshotError('invalid point checkpoint layout')
    with manifest.open('rb') as stream:
        content = stream.read(_LIMIT+1)
    if len(content) > _LIMIT:
        raise SnapshotError('point identity envelope exceeds 64 KiB')
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise SnapshotError('duplicate point checkpoint key')
            result[key] = value
        return result
    try:
        record = json.loads(content, object_pairs_hook=pairs)
    except (ValueError, RecursionError) as exc:
        raise SnapshotError('invalid point checkpoint JSON') from exc
    if (type(record) is not dict or set(record) != {'schema','version','identity','state_manifest_sha256'} or
            record['schema'] != _SCHEMA or record['version'] != __version__):
        raise SnapshotError('unsupported point checkpoint schema/version')
    if record['identity'] != point.record():
        raise SnapshotError('checkpoint run/point/material/initial identity differs')
    digest = record['state_manifest_sha256']
    inner = payload/'manifest.json'
    if (type(digest) is not str or not re.fullmatch('[0-9a-f]{64}', digest) or
            inner.is_symlink() or not inner.is_file() or
            inner.stat().st_size > min(_MANIFEST_LIMIT, max_bytes if max_bytes is not None else _MANIFEST_LIMIT) or
            _hash(inner, min(_MANIFEST_LIMIT, max_bytes if max_bytes is not None else _MANIFEST_LIMIT)) != digest):
        raise SnapshotError('point checkpoint inner manifest checksum mismatch')
    # Existing loader validates native binding/history and provides atomic replay.
    state.restore_checkpoint(payload, max_bytes=max_bytes)
