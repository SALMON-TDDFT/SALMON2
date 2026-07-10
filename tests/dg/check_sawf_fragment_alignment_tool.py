#!/usr/bin/env python3
"""Regression tests for symmetry-compatible periodic fragment maps."""

from fractions import Fraction
import importlib.util
import itertools
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "align_periodic_structure_to_fragments.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("sawf_fragment_alignment", TOOL)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {TOOL}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FragmentAlignmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.alignment = load_tool()

    def inversion(self, tau):
        return self.alignment.SymOp(
            rotation=((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
            translation=(tau, tau, tau),
        )

    def test_half_grid_inversion_maps_each_fragment_to_itself(self):
        grid_map = self.alignment.integer_grid_map(
            self.inversion(Fraction(15, 32)), (32, 32, 32)
        )

        self.assertEqual(grid_map.shift, (15, 15, 15))
        self.assertEqual(
            self.alignment.inversion_center_diagnostic(grid_map),
            (Fraction(15, 2), Fraction(15, 2), Fraction(15, 2)),
        )
        targets = self.alignment.fragment_target_enumeration(
            grid_map, (16, 16, 16)
        )
        self.assertEqual(len(targets), 8)
        for source, target_set in targets.items():
            self.assertEqual(target_set, frozenset((source,)))

    def test_half_grid_center_need_not_be_an_integer_grid_index(self):
        grid_map = self.alignment.integer_grid_map(
            self.inversion(Fraction(15, 32)), (32, 32, 32)
        )
        centers = self.alignment.inversion_center_diagnostic(grid_map)

        self.assertTrue(all(center.denominator == 2 for center in centers))
        self.assertEqual(grid_map.map_index((0, 0, 0)), (15, 15, 15))
        self.assertEqual(grid_map.map_index((15, 15, 15)), (0, 0, 0))

    def test_original_inversion_splits_each_fragment_across_eight_targets(self):
        grid_map = self.alignment.integer_grid_map(
            self.inversion(Fraction(1, 8)), (32, 32, 32)
        )
        targets = self.alignment.fragment_target_enumeration(
            grid_map, (16, 16, 16)
        )

        self.assertEqual(grid_map.shift, (4, 4, 4))
        expected_targets = frozenset(itertools.product(range(2), repeat=3))
        for target_set in targets.values():
            self.assertEqual(target_set, expected_targets)

    def test_identity_is_accepted(self):
        identity = self.alignment.SymOp(
            rotation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            translation=(0, 0, 0),
        )
        grid_map = self.alignment.integer_grid_map(identity, (12, 16, 20))

        self.assertEqual(grid_map.map_index((11, 15, 19)), (11, 15, 19))

    def test_nonintegral_grid_translation_is_rejected(self):
        operation = self.alignment.SymOp(
            rotation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            translation=(Fraction(1, 7), 0, 0),
        )

        with self.assertRaisesRegex(ValueError, "non-integral.*shift"):
            self.alignment.integer_grid_map(operation, (32, 32, 32))

    def test_fragment_shape_must_divide_mesh(self):
        identity = self.alignment.SymOp(
            rotation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            translation=(0, 0, 0),
        )
        grid_map = self.alignment.integer_grid_map(identity, (32, 32, 32))

        with self.assertRaisesRegex(ValueError, "divide"):
            self.alignment.fragment_target_enumeration(grid_map, (12, 16, 16))

    def test_mesh_dimensions_must_be_positive_integers(self):
        identity = self.alignment.SymOp(
            rotation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            translation=(0, 0, 0),
        )

        with self.assertRaisesRegex(ValueError, "positive integers"):
            self.alignment.integer_grid_map(identity, (32, 0, 32))

    def test_unequal_mesh_axis_swap_is_rejected(self):
        swap_xy = self.alignment.SymOp(
            rotation=((0, 1, 0), (1, 0, 0), (0, 0, 1)),
            translation=(0, 0, 0),
        )

        with self.assertRaisesRegex(ValueError, "non-integral.*coefficient"):
            self.alignment.integer_grid_map(swap_xy, (4, 8, 6))

    def test_unequal_mesh_compatible_finite_order_operation_is_a_bijection(self):
        sign_change = self.alignment.SymOp(
            rotation=((-1, 0, 0), (0, 1, 0), (0, 0, -1)),
            translation=(0, 0, 0),
        )
        grid_map = self.alignment.integer_grid_map(sign_change, (4, 8, 6))

        self.assertEqual(
            grid_map.coefficients,
            ((-1, 0, 0), (0, 1, 0), (0, 0, -1)),
        )
        images = {
            grid_map.map_index(index)
            for index in itertools.product(range(4), range(8), range(6))
        }
        self.assertEqual(len(images), 4 * 8 * 6)

    def test_non_unimodular_point_operation_is_rejected(self):
        with self.assertRaisesRegex(ValueError, r"determinant.*\+/-1"):
            self.alignment.SymOp(
                rotation=((3, 0, 0), (0, 1, 0), (0, 0, 1)),
                translation=(0, 0, 0),
            )

    def test_infinite_order_unipotent_shear_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "finite order"):
            self.alignment.SymOp(
                rotation=((1, 0, 0), (1, 1, 0), (0, 0, 1)),
                translation=(0, 0, 0),
            )

    def test_integer_grid_map_direct_construction_normalizes_immutable_values(self):
        grid_map = self.alignment.IntegerGridMap(
            mesh=[4, 8, 6],
            coefficients=[[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
            shift=[5, -1, 12],
        )

        self.assertEqual(grid_map.mesh, (4, 8, 6))
        self.assertEqual(
            grid_map.coefficients,
            ((-1, 0, 0), (0, 1, 0), (0, 0, -1)),
        )
        self.assertEqual(grid_map.shift, (1, 7, 0))
        with self.assertRaises(TypeError):
            grid_map.shift[0] = 0

    def test_integer_grid_map_direct_construction_rejects_malformed_values(self):
        valid_coefficients = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        invalid_cases = (
            ((4, 0, 6), valid_coefficients, (0, 0, 0)),
            ((4, 8, 6), ((1, 0), (0, 1), (0, 0)), (0, 0, 0)),
            ((4, 8, 6), ((1, 0, 0), (0, 1.5, 0), (0, 0, 1)), (0, 0, 0)),
            ((4, 8, 6), ((0, 0, 0), (0, 0, 0), (0, 0, 0)), (0, 0, 0)),
            ((4, 8, 6), valid_coefficients, (0, 0)),
            ((4, 8, 6), valid_coefficients, (0, "bad", 0)),
        )

        for mesh, coefficients, shift in invalid_cases:
            with self.subTest(mesh=mesh, coefficients=coefficients, shift=shift):
                with self.assertRaises(ValueError):
                    self.alignment.IntegerGridMap(mesh, coefficients, shift)

    def test_periodic_same_species_atom_bijection(self):
        inversion = self.alignment.SymOp(
            rotation=((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
            translation=(0, 0, 0),
        )
        atoms = (
            self.alignment.PeriodicAtom("Si", (Fraction(1, 8), 0, 0)),
            self.alignment.PeriodicAtom("Si", (Fraction(7, 8), 0, 0)),
        )

        permutation = self.alignment.periodic_same_species_atom_bijection(
            inversion, atoms
        )
        self.assertEqual(permutation, (1, 0))

    def test_periodic_atom_bijection_rejects_wrong_species(self):
        inversion = self.alignment.SymOp(
            rotation=((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
            translation=(0, 0, 0),
        )
        atoms = (
            self.alignment.PeriodicAtom("Si", (Fraction(1, 8), 0, 0)),
            self.alignment.PeriodicAtom("C", (Fraction(7, 8), 0, 0)),
        )

        with self.assertRaisesRegex(ValueError, "same-species.*bijection"):
            self.alignment.periodic_same_species_atom_bijection(inversion, atoms)

    def test_periodic_atom_bijection_respects_duplicate_multiplicity(self):
        inversion = self.alignment.SymOp(
            rotation=((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
            translation=(0, 0, 0),
        )
        atoms = (
            self.alignment.PeriodicAtom("Si", (Fraction(1, 8), 0, 0)),
            self.alignment.PeriodicAtom("Si", (Fraction(1, 8), 0, 0)),
            self.alignment.PeriodicAtom("Si", (Fraction(7, 8), 0, 0)),
        )

        with self.assertRaisesRegex(ValueError, "same-species.*bijection"):
            self.alignment.periodic_same_species_atom_bijection(inversion, atoms)


if __name__ == "__main__":
    unittest.main(verbosity=2)
