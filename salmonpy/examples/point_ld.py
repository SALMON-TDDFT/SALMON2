"""Two independent points of one synthetic LD material, with saved ID mapping."""
import argparse
from pathlib import Path
import numpy as np
from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import resolve_materials
from salmonpy.point_runs import PointRun, open_point_arrays


def run(directory):
    recipe = {'schema': 'salmonpy.materials/1', 'materials': {'film/α': {
        'backend': 'salmonpy.ld', 'model': 'ld', 'parameters': {
            'units': 'atomic_gaussian', 'plasma_frequency': .5,
            'strength': [1.], 'resonance': [.3], 'damping': [.1]}}}}
    catalog = resolve_materials(recipe, registry=builtin_material_registry(), base_directory=Path.cwd())
    model = material_from_catalog(catalog, 'film/α', base_directory=Path.cwd())
    dt = .02
    with PointRun({'right': model, 'left': model}) as calculation:
        right = calculation.new_state('right', dt=dt)
        left = calculation.new_state('left', dt=dt)
        current = [np.stack([right.current, left.current])]
        electric = [np.zeros((2, 3))]
        with model.new_state(dt=dt) as reference:
            for index in range(1, 21):
                field = np.array([.01*np.sin(index*.1), 0, 0])
                before = left.current.copy()
                right.step(field)
                np.testing.assert_array_equal(before, left.current)
                reference.step(field)
                np.testing.assert_array_equal(right.current, reference.current)
                left.step(-field)
                current.append(np.stack([right.current, left.current]))
                electric.append(np.stack([field, -field]))
        if np.array_equal(right.current, left.current):
            raise AssertionError('opposite drives produced identical currents')
        calculation.save_arrays(directory, {'current': current, 'E': electric, 'time': np.arange(21)*dt},
                                point_axes={'current': 1, 'E': 1},
                                metadata={'units': 'atomic_gaussian', 'drive': 'E',
                                          'scope': 'synthetic LD connection check',
                                          'direct_reference_exact': True, 'independent_states': True})
        record = calculation.record()
    loaded = open_point_arrays(directory)
    np.testing.assert_array_equal(loaded.array('current'), current)
    if loaded.metadata['point_run'] != record:
        raise AssertionError('identity roundtrip failed')
    return loaded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    args = parser.parse_args()
    result = run(args.directory)
    print('Passed direct-reference, independence and ID roundtrip:', result.directory)


if __name__ == '__main__':
    main()
