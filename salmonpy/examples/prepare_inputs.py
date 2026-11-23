"""Stage a user-provided synthetic response input for offline transfer.

Usage: python -m salmonpy.examples.prepare_inputs COEFFICIENT_JSON NEW_DIRECTORY
The input JSON contains only {"coefficient": 2.0}; this is not a real material.
Replace the explicit resolver with your model's validator for actual workflows.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

from salmonpy.config import MaterialRegistry
from salmonpy.errors import InputError
from salmonpy.preparation import prepare_inputs, verify_inputs


def resolve_response(model, parameters, *, base_directory):
    if set(parameters) != {'coefficient_file'} or parameters['coefficient_file'] != 'data/response.json':
        raise InputError('expected staged coefficient_file=data/response.json')
    path = base_directory/parameters['coefficient_file']
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise InputError('duplicate coefficient key')
            result[key] = value
        return result
    content = json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=pairs)
    if type(content) is not dict or set(content) != {'coefficient'}:
        raise InputError('input file must contain coefficient only')
    coefficient = content['coefficient']
    if type(coefficient) not in (int, float) or not math.isfinite(coefficient) or coefficient <= 0:
        raise InputError('coefficient must be finite and positive')
    return {'coefficient_file': 'data/response.json', 'coefficient': float(coefficient),
            'input_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file', type=Path)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    registry = MaterialRegistry()
    registry.register('local', version='example-1', models=['lab:response'],
                      resolve_parameters=resolve_response)
    bundle = prepare_inputs({
        'schema': 'salmonpy.materials/1',
        'materials': {'film': {'model': 'lab:response', 'backend': 'local',
                               'parameters': {'coefficient_file': 'data/response.json'}}}
    }, args.directory, registry=registry, base_directory=Path.cwd(),
        files={'data/response.json': args.input_file},
        provenance={'data/response.json': {'source': 'user-supplied synthetic example',
                                          'license': 'user-managed'}})
    manifest_hash = hashlib.sha256((bundle/'manifest.json').read_bytes()).hexdigest()
    verify_inputs(bundle, manifest_sha256=manifest_hash)
    print(json.dumps({'directory': str(bundle), 'manifest_sha256': manifest_hash,
                      'scope': 'synthetic input preparation, no solver execution'}, indent=2))


if __name__ == '__main__':
    main()
