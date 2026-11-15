"""Run two explicitly configured synthetic LD materials; outputs outside inputs."""
import argparse
import json
from pathlib import Path
import numpy as np

from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import load_materials
from salmonpy.preparation import prepare_inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    registry = builtin_material_registry()
    recipe = {'schema': 'salmonpy.materials/1', 'materials': {
        key: {'backend': 'salmonpy.ld', 'model': 'ld', 'parameters': {
            'units': 'atomic_gaussian', 'plasma_frequency': frequency,
            'strength': [1.0], 'resonance': [.3], 'damping': [.1]}}
        for key, frequency in [('film-left', .5), ('film-right', .7)]}}
    args.directory.mkdir(parents=True, exist_ok=False)
    bundle = prepare_inputs(recipe, args.directory/'inputs', registry=registry,
                            base_directory=Path.cwd(), files={})
    catalog = load_materials(bundle/'materials.json', registry=registry)
    records = {}
    for key in catalog.materials:
        material = material_from_catalog(catalog, key, base_directory=bundle)
        with material.new_state(dt=.02) as a, material.new_state(dt=.02) as b:
            currents = []
            for i in range(1, 21):
                field = [.01*np.sin(i*.1), 0, 0]
                before = b.current.copy()
                a.step(field)
                np.testing.assert_array_equal(b.current, before)
                b.step(field)
                np.testing.assert_array_equal(a.current, b.current)
                currents.append(a.current.tolist())
            records[key] = {**material.record(), 'current': currents, 'dt': .02,
                            'drive': 'E', 'units': 'atomic_gaussian', 'independent_states': True}
    result = {'scope': 'synthetic LD connection test, not validated material data', 'materials': records}
    (args.directory/'response.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'directory': str(args.directory), 'materials': list(records), 'independence': 'passed'}))


if __name__ == '__main__':
    main()
