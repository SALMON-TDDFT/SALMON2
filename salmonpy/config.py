"""Configuration-only material catalog; no solver creation or plugin discovery.

Backend resolvers are explicitly registered trusted Python code. They validate
model-specific settings and return JSON-compatible resolved parameter objects.
"""
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from types import MappingProxyType
import unicodedata

from .errors import InputError

__all__ = ['MaterialRegistry', 'MaterialCatalog', 'resolve_materials', 'load_materials']
_SCHEMA = 'salmonpy.materials/1'


def _copy_json(value, path='$', active=None):
    """Detach JSON values, rejecting lossy conversions and recursive containers."""
    if active is None:
        active = set()
    if value is None or type(value) is bool:
        return value
    if type(value) is int:
        try:
            str(value)  # Honor the host interpreter's integer serialization limit.
        except ValueError as exc:
            raise InputError(f'{path}: integer exceeds JSON serialization limit') from exc
        return value
    if type(value) is float and math.isfinite(value):
        return value
    if type(value) is str:
        try:
            value.encode('utf-8')
        except UnicodeEncodeError as exc:
            raise InputError(f'{path}: invalid Unicode string') from exc
        return value
    if type(value) not in (dict, list):
        raise InputError(f'{path}: expected JSON values with finite numbers')
    if id(value) in active:
        raise InputError(f'{path}: cyclic container')
    active.add(id(value))
    try:
        if isinstance(value, list):
            return [_copy_json(v, f'{path}[{i}]', active) for i, v in enumerate(value)]
        result = {}
        for key, item in value.items():
            if type(key) is not str:
                raise InputError(f'{path}: object keys must be strings')
            _copy_json(key, path, active)
            result[key] = _copy_json(item, f'{path}[{key!r}]', active)
        return result
    finally:
        active.remove(id(value))


def _object(value, path, required=(), optional=()):
    if type(value) is not dict:
        raise InputError(f'{path}: expected an object')
    missing = set(required) - value.keys()
    unknown = value.keys() - set(required) - set(optional)
    if missing or unknown:
        raise InputError(f'{path}: missing keys {sorted(missing)}; unknown keys {sorted(unknown)}')


def _identifier(value, path):
    if (type(value) is not str or not value or value != value.strip()
            or any(unicodedata.category(c) in ('Cc', 'Cs') for c in value)):
        raise InputError(f'{path}: expected a nonempty ID without edge whitespace or control characters')
    return value


def _freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _thaw(value):
    if isinstance(value, MappingProxyType):
        return {key: _thaw(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw(item) for item in value]
    return value


@dataclass(frozen=True)
class _Backend:
    version: str
    models: frozenset
    resolver: object


class MaterialRegistry:
    """Explicit local registry; registration never calls the resolver.

    register(name, *, version, models, resolve_parameters) binds a trusted
    callable ``resolve_parameters(model, parameters, *, base_directory)``.
    It must validate unknown parameters, units, ranges and combinations and
    return a plain JSON object. Do not allocate solver states in this callback.
    Registration does not certify execution capabilities or physical accuracy.
    """

    def __init__(self):
        self._backends = {}

    def register(self, name, *, version, models, resolve_parameters):
        _identifier(name, 'backend')
        _identifier(version, 'backend.version')
        if name in self._backends:
            raise InputError(f'backend {name!r}: already registered')
        if not isinstance(models, (list, tuple)) or not models:
            raise InputError('backend.models: expected a nonempty list or tuple')
        names = tuple(_identifier(model, 'backend.models') for model in models)
        if len(set(names)) != len(names):
            raise InputError('backend.models: duplicate model')
        if not callable(resolve_parameters):
            raise InputError('backend.resolve_parameters: expected a callable')
        self._backends[name] = _Backend(version, frozenset(names), resolve_parameters)


@dataclass(frozen=True, init=False)
class MaterialCatalog:
    """Read-only resolved definitions, reference recipes and definition hashes.

    Obtain via resolve_materials/load_materials. Initial-state recipes are only
    structurally checked, not opened/prepared or certified compatible. Hashes
    cover resolved definitions, not referenced file contents or callback code.
    """
    materials: object
    initial_states: object
    definition_digests: object

    def __init__(self):
        raise InputError('use resolve_materials or load_materials to create a catalog')

    def as_dict(self):
        """Independent JSON-safe record, not a loadable recipe or checkpoint."""
        return {'schema': 'salmonpy.materials-resolved/1',
                'materials': _thaw(self.materials),
                'initial_states': _thaw(self.initial_states),
                'definition_digests': _thaw(self.definition_digests)}


def resolve_materials(data, *, registry, base_directory):
    """Validate a plain Python JSON object through explicitly registered resolvers.

    IDs/references and schema are checked before any callback. ValueError or
    TypeError from a resolver is wrapped in InputError with the material path;
    other callback and OS exceptions propagate. No partial catalog is returned.
    base_directory is explicit; relative file access belongs to the resolver.
    """
    if not isinstance(registry, MaterialRegistry):
        raise InputError('registry: expected MaterialRegistry')
    try:
        data = _copy_json(data)
    except RecursionError as exc:
        raise InputError('configuration nesting is too deep') from exc
    _object(data, '$', ('schema', 'materials'), ('initial_states',))
    if data['schema'] != _SCHEMA:
        raise InputError(f'schema: expected {_SCHEMA}')
    materials = data['materials']
    initial = data.get('initial_states', {})
    if type(materials) is not dict or not materials:
        raise InputError('materials: expected a nonempty object')
    if type(initial) is not dict:
        raise InputError('initial_states: expected an object')
    selected = {}
    for key, definition in materials.items():
        path = f'materials[{key!r}]'
        _identifier(key, path)
        _object(definition, path, ('model', 'backend', 'parameters'), ('label',))
        _identifier(definition['model'], path + '.model')
        _identifier(definition['backend'], path + '.backend')
        if 'label' in definition and type(definition['label']) is not str:
            raise InputError(f'{path}.label: expected string')
        if type(definition['parameters']) is not dict:
            raise InputError(f'{path}.parameters: expected an object')
        backend = registry._backends.get(definition['backend'])
        if backend is None or definition['model'] not in backend.models:
            raise InputError(f'{path}: unregistered backend/model pair')
        selected[key] = backend
    for key, entry in initial.items():
        path = f'initial_states[{key!r}]'
        _identifier(key, path)
        _object(entry, path, ('material_id',), ('source', 'prepare'))
        _identifier(entry['material_id'], path + '.material_id')
        if entry['material_id'] not in materials:
            raise InputError(f'{path}.material_id: unknown material')
        if ('source' in entry) == ('prepare' in entry):
            raise InputError(f'{path}: specify exactly one of source and prepare')
        if type(entry.get('source', entry.get('prepare'))) is not dict:
            raise InputError(f'{path}: source/prepare must be an object')
    try:
        base = Path(base_directory).resolve()
    except (TypeError, ValueError) as exc:
        raise InputError('base_directory: expected a path') from exc
    digests = {}
    for key, definition in materials.items():
        backend = selected[key]
        path = f'materials[{key!r}].parameters'
        try:
            parameters = backend.resolver(definition['model'], definition['parameters'],
                                          base_directory=base)
            parameters = _copy_json(parameters, path)
        except RecursionError as exc:
            raise InputError(f'{path}: configuration nesting is too deep') from exc
        except (ValueError, TypeError) as exc:
            raise InputError(f'{path}: {exc}') from exc
        if type(parameters) is not dict:
            raise InputError(f'{path}: resolver must return an object')
        definition['parameters'] = parameters
        definition['backend_version'] = backend.version
        hashed = {k: v for k, v in definition.items() if k != 'label'}
        encoded = json.dumps({'schema': 'salmonpy.material-definition/1', **hashed},
                             sort_keys=True, separators=(',', ':'), ensure_ascii=False,
                             allow_nan=False).encode('utf-8')
        digests[key] = hashlib.sha256(encoded).hexdigest()
    result = object.__new__(MaterialCatalog)
    for name, value in (('materials', materials), ('initial_states', initial),
                        ('definition_digests', digests)):
        object.__setattr__(result, name, _freeze(value))
    return result


def load_materials(path, *, registry):
    """Read UTF-8 JSON, reject duplicate keys, resolve relative to its directory.

    OS errors propagate. JSON content never imports modules or launches solvers.
    """
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise InputError(f'JSON: duplicate key {key!r}')
            result[key] = value
        return result

    def constant(value):
        raise InputError(f'JSON: non-finite number {value}')

    path = Path(path).absolute()
    try:
        data = json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=pairs,
                          parse_constant=constant)
    except (ValueError, RecursionError) as exc:
        raise InputError(f'{path.name}: invalid JSON: {exc}') from exc
    return resolve_materials(data, registry=registry, base_directory=path.parent)
