"""Two user material IDs with separate validated settings; no solver execution."""
from pathlib import Path
import json
import math

from salmonpy.config import MaterialRegistry, resolve_materials
from salmonpy.errors import InputError


def resolve_response(model, parameters, *, base_directory):
    # This illustrative model uses a dimensionless coefficient, not a real material.
    if set(parameters) != {'coefficient'}:
        raise InputError('expected coefficient only')
    value = parameters['coefficient']
    if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
        raise InputError('coefficient must be finite and positive')
    return {'coefficient': float(value)}


def main():
    registry = MaterialRegistry()
    registry.register('local', version='example-1', models=['lab:response'],
                      resolve_parameters=resolve_response)
    catalog = resolve_materials({
        'schema': 'salmonpy.materials/1',
        'materials': {
            name: {'model': 'lab:response', 'backend': 'local',
                   'parameters': {'coefficient': value}}
            for name, value in [('film-left', 2), ('film-right', 3)]
        },
        'initial_states': {
            'left-ground': {'material_id': 'film-left', 'prepare': {'kind': 'equilibrium'}}
        }
    }, registry=registry, base_directory=Path.cwd())
    print(json.dumps(catalog.as_dict(), indent=2))


if __name__ == '__main__':
    main()
