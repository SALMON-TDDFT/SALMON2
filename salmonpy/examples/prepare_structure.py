"""Prepare a local structure and grid report for offline transfer; no DFT run."""
import argparse
import hashlib
import json
from pathlib import Path

from salmonpy.config import MaterialRegistry
from salmonpy.errors import InputError
from salmonpy.preparation import prepare_inputs
from salmonpy.structure import StructureInput, grid_from_spacing


def resolve_structure(model, parameters, *, base_directory):
    expected = {'structure_file', 'format', 'index', 'grid_spacing', 'grid_spacing_unit', 'kmesh'}
    if set(parameters) != expected or parameters['structure_file'] != 'structures/input':
        raise InputError('expected explicit structure/frame/spacing/unit/kmesh parameters')
    kmesh = parameters['kmesh']
    if (type(kmesh) is not list or len(kmesh) != 3
            or any(type(n) is not int or not 1 <= n <= 2**20 for n in kmesh)):
        raise InputError('kmesh must contain three positive integers <= 2**20')
    structure = StructureInput.from_file(base_directory/'structures/input',
                                        format=parameters['format'], index=parameters['index'])
    atoms = structure.atoms
    if not all(atoms.pbc):
        raise InputError('this example requires three periodic axes')
    grid = grid_from_spacing(atoms.cell.array, spacing=parameters['grid_spacing'],
                             unit=parameters['grid_spacing_unit'])
    return {**parameters, 'structure': structure.record(), 'grid': grid.as_dict()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('structure_file', type=Path)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--format', required=True, choices=['cif', 'extxyz'])
    parser.add_argument('--index', type=int, default=0)
    parser.add_argument('--material-id', default='material')
    parser.add_argument('--spacing', type=float, required=True)
    parser.add_argument('--unit', required=True, choices=['angstrom', 'bohr'])
    parser.add_argument('--kmesh', type=int, nargs=3, required=True)
    args = parser.parse_args()
    registry = MaterialRegistry()
    registry.register('structure-preview', version='example-1', models=['structure-grid'],
                      resolve_parameters=resolve_structure)
    data = {'schema': 'salmonpy.materials/1', 'materials': {
        args.material_id: {'backend': 'structure-preview', 'model': 'structure-grid',
                          'parameters': {'structure_file': 'structures/input',
                                         'format': args.format, 'index': args.index,
                                         'grid_spacing': args.spacing,
                                         'grid_spacing_unit': args.unit, 'kmesh': args.kmesh}}}}
    directory = prepare_inputs(data, args.directory, registry=registry,
        base_directory=Path.cwd(), files={'structures/input': args.structure_file},
        provenance={'structures/input': {'source': 'user-provided structure', 'license': 'user-managed'}})
    record = json.loads((directory/'resolved.json').read_text())
    print(json.dumps({'directory': str(directory), 'materials': record['materials'],
                     'manifest_sha256': hashlib.sha256((directory/'manifest.json').read_bytes()).hexdigest(),
                     'scope': 'structure and grid preparation only, no physical convergence or solver run'}, indent=2))


if __name__ == '__main__':
    main()
