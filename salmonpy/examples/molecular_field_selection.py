"""Pure electric-drive selection; no solver, MPI or physical convergence claim."""
import json
import numpy as np
from salmonpy.field_input import FieldValue
from salmonpy.dipole_source import select_molecular_field


def main():
    # Artificial additive fields: own vacuum contribution + external + other molecule.
    own = np.array([4., 5., 6.])
    external = np.array([1., 0., 0.])
    other = np.array([0., 2., 0.])
    def sample(values):
        return FieldValue(values, time=.1, position=[1., 2., 3.], unit='E:au')
    provenance = dict(footprint_id='example-footprint/1', sample_map_id='example-map/1',
                      own_source_history_id='molecule-1-history/1')
    selected = select_molecular_field(sample(own+external+other), molecule_id='molecule-1',
        policy='environment_scattered_only', provenance=provenance,
        reference=sample(own), reference_provenance=dict(provenance,
            source_ids=['molecule-1'], medium='vacuum', incident=False))
    np.testing.assert_array_equal(selected.field.value, external+other)
    print(json.dumps(selected.record(), indent=2))


if __name__ == '__main__':
    main()
