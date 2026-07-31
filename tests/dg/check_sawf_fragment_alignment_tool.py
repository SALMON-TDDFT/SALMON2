#!/usr/bin/env python3
"""Regression tests for symmetry-compatible periodic fragment maps."""

from fractions import Fraction
import importlib.util
import itertools
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "align_periodic_structure_to_fragments.py"
SAMPLE = ROOT / "tests/dg/data/si64_sawf_alignment"
SOURCE_ATOM = SAMPLE / "atom.dat"
ALIGNED_ATOM = SAMPLE / "atom_sawf_aligned.dat"
ALIGNED_SYMMETRY = SAMPLE / "sym_sawf_aligned.dat"
ALIGNED_INPUT = SAMPLE / "inputfile_gs_w90_pseudo_sawf_aligned_nw576_nb664"
C64_SYMMETRY_TEXT = (
    "# Identity and inversion about fractional position (1/16, 1/16, 1/16).\n"
    " 1  0  0  0.000\n"
    " 0  1  0  0.000\n"
    " 0  0  1  0.000\n"
    "-1  0  0  0.125\n"
    " 0 -1  0  0.125\n"
    " 0  0 -1  0.125\n"
)


def write_c64_symmetry(directory, name="sym.dat"):
    path = Path(directory) / name
    path.write_text(C64_SYMMETRY_TEXT, encoding="ascii")
    return path


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

    def test_metric_incompatible_finite_order_operation_is_rejected(self):
        swap_xy = self.alignment.SymOp(
            rotation=((0, 1, 0), (1, 0, 0), (0, 0, 1)),
            translation=(0, 0, 0),
        )

        with self.assertRaisesRegex(ValueError, "lattice metric"):
            self.alignment.validate_lattice_metric(swap_xy, (10, 12, 14))

    def test_buffer_must_follow_signed_axis_permutation(self):
        swap_xy = self.alignment.SymOp(
            rotation=((0, 1, 0), (1, 0, 0), (0, 0, 1)),
            translation=(0, 0, 0),
        )

        with self.assertRaisesRegex(ValueError, "buffer.*axis"):
            self.alignment.validate_buffer_geometry(swap_xy, (4, 6, 8))
        self.alignment.validate_buffer_geometry(swap_xy, (6, 6, 8))

    def test_symmetry_parser_preserves_order_and_validates_groups(self):
        with tempfile.TemporaryDirectory() as temporary:
            sym_path = Path(temporary) / "sym.dat"
            sym_path.write_text(
                "# identity then inversion\n"
                " 1 0 0 0\n 0 1 0 0\n 0 0 1 0\n"
                "-1 0 0 1/8\n 0 -1 0 1/8\n 0 0 -1 1/8\n",
                encoding="ascii",
            )
            operations = self.alignment.parse_symmetry_file(sym_path)

        self.assertEqual(len(operations), 2)
        self.assertEqual(operations[0].rotation, self.alignment._IDENTITY)
        self.assertEqual(operations[1].translation, (Fraction(1, 8),) * 3)
        self.alignment.validate_symmetry_group(operations)

    def test_symmetry_group_rejects_missing_products_and_duplicate_identity(self):
        identity = self.alignment.SymOp(self.alignment._IDENTITY, (0, 0, 0))
        quarter_turn = self.alignment.SymOp(
            ((0, -1, 0), (1, 0, 0), (0, 0, 1)), (0, 0, 0)
        )
        with self.assertRaisesRegex(ValueError, "closure"):
            self.alignment.validate_symmetry_group((identity, quarter_turn))
        with self.assertRaisesRegex(ValueError, "identity.*exactly once"):
            self.alignment.validate_symmetry_group((identity, identity))

    def test_nonintegral_base_shift_can_be_aligned_by_a_half_integer_candidate(self):
        identity = self.alignment.SymOp(self.alignment._IDENTITY, (0, 0, 0))
        inversion = self.inversion(Fraction(1, 8))
        atom = self.alignment.PeriodicAtom("C", (Fraction(1, 16),) * 3)

        result = self.alignment.find_fragment_compatible_translation(
            (atom,), (identity, inversion), (1, 1, 1),
            (4, 4, 4), (2, 2, 2), (0, 0, 0)
        )

        self.assertEqual(result.translation, (Fraction(1, 16),) * 3)
        self.assertEqual(result.grid_maps[1].shift, (1, 1, 1))
        for source, target in result.fragment_maps[1]:
            self.assertEqual(source, target)

    def test_tolerant_group_validation_accepts_decimal_thirds(self):
        translations = tuple(
            self.alignment.SymOp(
                self.alignment._IDENTITY, (Fraction(index, 3), 0, 0)
            )
            for index in range(3)
        )
        rendered = self.alignment._render_symmetry(
            translations, (0, 0, 0), (0, 0, 0)
        )
        with tempfile.TemporaryDirectory() as temporary:
            sym_path = Path(temporary) / "thirds.dat"
            sym_path.write_text(rendered, encoding="ascii")
            reread = self.alignment.parse_symmetry_file(sym_path)

        with self.assertRaisesRegex(ValueError, "closure"):
            self.alignment.validate_symmetry_group(reread)
        self.alignment.validate_symmetry_group_tolerant(
            reread, Fraction(1, 10**12)
        )

        inconsistent = (
            reread[0],
            reread[1],
            self.alignment.SymOp(self.alignment._IDENTITY, (Fraction(67, 100), 0, 0)),
        )
        with self.assertRaisesRegex(ValueError, "closure"):
            self.alignment.validate_symmetry_group_tolerant(
                inconsistent, Fraction(1, 10**12)
            )

    def test_validated_publication_accepts_decimal_thirds_roundtrip(self):
        cell = (Fraction(1), Fraction(1), Fraction(1))
        atoms = tuple(
            self.alignment.PeriodicAtom(
                "C", (Fraction(index, 3), 0, 0), "'C'", "1"
            )
            for index in range(3)
        )
        operations = tuple(
            self.alignment.SymOp(
                self.alignment._IDENTITY, (Fraction(index, 3), 0, 0)
            )
            for index in range(3)
        )
        atom_text = self.alignment._render_atoms(atoms, cell)
        self.assertIn("0.333333333333333333333333", atom_text)
        self.assertIn("0.666666666666666666666667", atom_text)
        with tempfile.TemporaryDirectory() as temporary:
            output_atom = Path(temporary) / "atom.dat"
            output_sym = Path(temporary) / "sym.dat"
            self.alignment._write_validated_outputs(
                output_atom, output_sym,
                atom_text,
                self.alignment._render_symmetry(operations, (0, 0, 0), (0, 0, 0)),
                cell, atoms, operations, Fraction(1, 10**12), False
            )

            self.assertTrue(output_atom.exists())
            self.assertNotIn("/", "\n".join(
                line for line in output_sym.read_text(encoding="ascii").splitlines()
                if not line.startswith("#")
            ))

    def test_denominator_three_search_uses_residue_filtering(self):
        identity = self.alignment.SymOp(self.alignment._IDENTITY, (0, 0, 0))
        inversion = self.inversion(Fraction(1, 96))
        atom = self.alignment.PeriodicAtom("C", (Fraction(1, 192),) * 3)

        started = time.perf_counter()
        result = self.alignment.find_fragment_compatible_translation(
            (atom,), (identity, inversion), (1, 1, 1),
            (32, 32, 32), (2, 2, 2), (0, 0, 0)
        )
        elapsed = time.perf_counter() - started

        self.assertEqual(result.translation, (Fraction(11, 48),) * 3)
        self.assertEqual(result.grid_maps[1].shift, (15, 15, 15))
        self.assertLess(elapsed, 10.0, f"denominator-three search took {elapsed:.3f} s")

    def test_already_compatible_identity_uses_zero_translation_fast_path(self):
        identity = self.alignment.SymOp(self.alignment._IDENTITY, (0, 0, 0))
        atom = self.alignment.PeriodicAtom("C", (0, 0, 0))

        with mock.patch.object(
            self.alignment,
            "_candidate_grid_shifts",
            wraps=self.alignment._candidate_grid_shifts,
        ) as grid_shifts:
            result = self.alignment.find_fragment_compatible_translation(
                (atom,), (identity,), (1, 1, 1),
                (32, 32, 32), (2, 2, 2), (0, 0, 0)
            )

        self.assertEqual(result.translation, (0, 0, 0))
        self.assertEqual(grid_shifts.call_count, 1)

    def test_c64_search_selects_half_grid_center_and_is_deterministic(self):
        atom_path = (
            ROOT
            / "tests/dg/data/si64_sawf_alignment/atom.dat"
        )
        atoms = self.alignment.parse_atom_file(atom_path, (Fraction("13.44"),) * 3)
        with tempfile.TemporaryDirectory() as temporary:
            operations = self.alignment.parse_symmetry_file(
                write_c64_symmetry(temporary)
            )

        started = time.perf_counter()
        first = self.alignment.find_fragment_compatible_translation(
            atoms, operations, (Fraction("13.44"),) * 3,
            (32, 32, 32), (2, 2, 2), (6, 6, 6)
        )
        second = self.alignment.find_fragment_compatible_translation(
            atoms, operations, (Fraction("13.44"),) * 3,
            (32, 32, 32), (2, 2, 2), (6, 6, 6)
        )
        elapsed = time.perf_counter() - started

        self.assertEqual(len(atoms), 64)
        self.assertEqual(first.translation, (Fraction(11, 64),) * 3)
        self.assertEqual(first.operations[1].translation, (Fraction(15, 32),) * 3)
        self.assertEqual(first.grid_maps[1].shift, (15, 15, 15))
        self.assertEqual(
            self.alignment.inversion_center_diagnostic(first.grid_maps[1]),
            (Fraction(15, 2),) * 3,
        )
        self.assertEqual(first, second)
        self.assertLess(elapsed, 10.0, f"two C64 searches took {elapsed:.3f} s")

    def test_committed_c64_aligned_files_are_exact_tool_outputs(self):
        source_atom_before = SOURCE_ATOM.read_bytes()
        legacy_sym = SAMPLE / "sym.dat"
        legacy_sym_before = legacy_sym.read_bytes() if legacy_sym.exists() else None

        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            source_sym = write_c64_symmetry(temporary, "source_sym.dat")
            regenerated_atom = temporary / ALIGNED_ATOM.name
            regenerated_symmetry = temporary / ALIGNED_SYMMETRY.name
            completed = subprocess.run(
                [
                    sys.executable, str(TOOL),
                    "--input-atom", str(SOURCE_ATOM),
                    "--input-sym", str(source_sym),
                    "--output-atom", str(regenerated_atom),
                    "--output-sym", str(regenerated_symmetry),
                    "--cell", "13.44", "13.44", "13.44",
                    "--mesh", "32", "32", "32",
                    "--fragments", "2", "2", "2",
                    "--buffer", "6", "6", "6",
                ],
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(regenerated_atom.read_bytes(), ALIGNED_ATOM.read_bytes())
            self.assertEqual(
                regenerated_symmetry.read_bytes(), ALIGNED_SYMMETRY.read_bytes()
            )

        aligned_atoms = self.alignment.parse_atom_file(
            ALIGNED_ATOM, (Fraction("13.44"),) * 3
        )
        source_atoms = self.alignment.parse_atom_file(
            SOURCE_ATOM, (Fraction("13.44"),) * 3
        )
        operations = self.alignment.parse_symmetry_file(ALIGNED_SYMMETRY)
        self.assertEqual(len(aligned_atoms), 64)
        self.assertTrue(ALIGNED_ATOM.read_text(encoding="ascii").startswith("  'C'"))
        self.assertFalse(ALIGNED_ATOM.read_text(encoding="ascii").startswith("#"))
        self.assertEqual(
            [(atom.label, atom.species_index, atom.trailing) for atom in aligned_atoms],
            [(atom.label, atom.species_index, atom.trailing) for atom in source_atoms],
        )
        expected_translation = (Fraction(11, 64),) * 3
        for source, aligned in zip(source_atoms, aligned_atoms):
            self.assertEqual(
                aligned.position,
                tuple((value + shift) % 1 for value, shift in zip(
                    source.position, expected_translation
                )),
            )
        self.assertEqual(operations[1].translation, (Fraction(15, 32),) * 3)
        self.assertIn("# translation=11/64 11/64 11/64", ALIGNED_SYMMETRY.read_text())
        self.assertIn("# buffer=6 6 6", ALIGNED_SYMMETRY.read_text())
        for operation in operations:
            permutation = self.alignment.periodic_same_species_atom_bijection(
                operation, aligned_atoms
            )
            self.assertEqual(sorted(permutation), list(range(64)))
            grid_map = self.alignment.integer_grid_map(operation, (32, 32, 32))
            targets = self.alignment.fragment_target_enumeration(
                grid_map, (16, 16, 16)
            )
            self.assertEqual(len(targets), 8)
            self.assertTrue(all(targets[source] == frozenset((source,)) for source in targets))
        inversion_map = self.alignment.integer_grid_map(operations[1], (32, 32, 32))
        self.assertEqual(inversion_map.shift, (15, 15, 15))
        self.assertEqual(
            self.alignment.inversion_center_diagnostic(inversion_map),
            (Fraction(15, 2),) * 3,
        )
        self.assertEqual(Fraction("13.44") * Fraction(11, 64), Fraction("2.31"))
        self.assertEqual(SOURCE_ATOM.read_bytes(), source_atom_before)
        if legacy_sym_before is not None:
            self.assertEqual(legacy_sym.read_bytes(), legacy_sym_before)

    def test_aligned_input_uses_only_dedicated_aligned_geometry(self):
        text = ALIGNED_INPUT.read_text(encoding="ascii")
        active = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
        required = (
            "sysname = 'c64_dc_pseudo_sawf_aligned_nw576_nb664'",
            "file_atom_coor = 'atom_sawf_aligned.dat'",
            "wannier_symmetry_file = 'sym_sawf_aligned.dat'",
            "wannier_site_symmetry = 'file'",
            "wannier_symmetry_tolerance = 1.0d-6",
            "wannier_num_wann = 576",
            "wannier_num_bands = 664",
            "nstate = 664",
            "num_rgrid(1:3) = 32,32,32",
            "num_fragment(1:3) = 2,2,2",
            "num_rgrid_buffer(1:3) = 6,6,6",
            "yn_eigenexa = 'y'",
        )
        for setting in required:
            self.assertEqual(active.count(setting), 1, setting)
        self.assertNotIn("file_atom_coor = 'atom.dat'", active)
        self.assertNotIn("wannier_symmetry_file = 'sym.dat'", active)
        self.assertNotRegex(active, r"(?m)^\s*(?:nstate|wannier_num_(?:wann|bands))\s*=\s*256\b")
        self.assertNotIn("yn_conventional_from_dcdft", active.lower())

    def test_cli_publishes_atomically_and_preserves_inputs(self):
        source_atom = (
            ROOT
            / "tests/dg/data/si64_sawf_alignment/atom.dat"
        )
        atom_before = source_atom.read_bytes()
        with tempfile.TemporaryDirectory() as temporary:
            source_sym = write_c64_symmetry(temporary, "input_sym.dat")
            sym_before = source_sym.read_bytes()
            output_atom = Path(temporary) / "atom_sawf_aligned.dat"
            output_sym = Path(temporary) / "sym_sawf_aligned.dat"
            command = [
                sys.executable, str(TOOL),
                "--input-atom", str(source_atom),
                "--input-sym", str(source_sym),
                "--output-atom", str(output_atom),
                "--output-sym", str(output_sym),
                "--cell", "13.44", "13.44", "13.44",
                "--mesh", "32", "32", "32",
                "--fragments", "2", "2", "2",
                "--buffer", "6", "6", "6",
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("translation=11/64 11/64 11/64", completed.stdout)
            self.assertIn("buffer=6 6 6", completed.stdout)
            self.assertTrue(output_atom.read_text(encoding="ascii").startswith("  'C'"))
            self.assertIn("buffer=6 6 6", output_sym.read_text(encoding="ascii"))
            self.assertEqual(len(self.alignment.parse_atom_file(
                output_atom, (Fraction("13.44"),) * 3
            )), 64)

            refused = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("already exists", refused.stderr)
            self.assertEqual(source_sym.read_bytes(), sym_before)

        self.assertEqual(source_atom.read_bytes(), atom_before)

    def test_cli_failure_creates_no_partial_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            atom_path = temporary / "atom.dat"
            sym_path = temporary / "sym.dat"
            output_atom = temporary / "aligned_atom.dat"
            output_sym = temporary / "aligned_sym.dat"
            atom_path.write_text("  'Si' 0.5 0 0 1\n  'C' 1.5 0 0 2\n", encoding="ascii")
            sym_path.write_text(
                "1 0 0 0\n0 1 0 0\n0 0 1 0\n"
                "-1 0 0 0\n0 -1 0 0\n0 0 -1 0\n",
                encoding="ascii",
            )
            command = [
                sys.executable, str(TOOL),
                "--input-atom", str(atom_path), "--input-sym", str(sym_path),
                "--output-atom", str(output_atom), "--output-sym", str(output_sym),
                "--cell", "2", "2", "2", "--mesh", "4", "4", "4",
                "--fragments", "2", "2", "2", "--buffer", "1", "1", "1",
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)

            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("same-species", completed.stderr)
            self.assertFalse(output_atom.exists())
            self.assertFalse(output_sym.exists())

    def test_cli_force_still_rejects_an_output_that_aliases_an_input(self):
        source_atom = (
            ROOT
            / "tests/dg/data/si64_sawf_alignment/atom.dat"
        )
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            atom_path = temporary / "atom.dat"
            sym_path = temporary / "sym.dat"
            output_sym = temporary / "aligned_sym.dat"
            atom_path.write_bytes(source_atom.read_bytes())
            sym_path.write_text(C64_SYMMETRY_TEXT, encoding="ascii")
            atom_before = atom_path.read_bytes()
            command = [
                sys.executable, str(TOOL),
                "--input-atom", str(atom_path), "--input-sym", str(sym_path),
                "--output-atom", str(atom_path), "--output-sym", str(output_sym),
                "--cell", "13.44", "13.44", "13.44",
                "--mesh", "32", "32", "32", "--fragments", "2", "2", "2",
                "--buffer", "6", "6", "6", "--force",
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)

            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("dedicated", completed.stderr)
            self.assertEqual(atom_path.read_bytes(), atom_before)
            self.assertFalse(output_sym.exists())

    def test_publication_failure_rolls_back_the_first_output(self):
        atom = self.alignment.PeriodicAtom("C", (0, 0, 0), "'C'", "1")
        identity = self.alignment.SymOp(self.alignment._IDENTITY, (0, 0, 0))
        cell = (Fraction(2),) * 3
        with tempfile.TemporaryDirectory() as temporary:
            output_atom = Path(temporary) / "atom.dat"
            output_sym = Path(temporary) / "sym.dat"
            real_replace = self.alignment.os.replace
            calls = 0

            def fail_second_replace(source, destination):
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise OSError("injected publication failure")
                return real_replace(source, destination)

            with mock.patch.object(self.alignment.os, "replace", fail_second_replace):
                with self.assertRaisesRegex(OSError, "injected"):
                    self.alignment._write_validated_outputs(
                        output_atom,
                        output_sym,
                        self.alignment._render_atoms((atom,), cell),
                        self.alignment._render_symmetry((identity,), (0, 0, 0), (1, 1, 1)),
                        cell,
                        (atom,),
                        (identity,),
                        Fraction(1, 10**12),
                        False,
                    )

            self.assertFalse(output_atom.exists())
            self.assertFalse(output_sym.exists())

    def test_keyboard_interrupt_during_publication_rolls_back_first_output(self):
        atom = self.alignment.PeriodicAtom("C", (0, 0, 0), "'C'", "1")
        identity = self.alignment.SymOp(self.alignment._IDENTITY, (0, 0, 0))
        cell = (Fraction(2),) * 3
        with tempfile.TemporaryDirectory() as temporary:
            output_atom = Path(temporary) / "atom.dat"
            output_sym = Path(temporary) / "sym.dat"
            real_replace = self.alignment.os.replace
            calls = 0

            def interrupt_second_replace(source, destination):
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise KeyboardInterrupt()
                return real_replace(source, destination)

            with mock.patch.object(self.alignment.os, "replace", interrupt_second_replace):
                with self.assertRaises(KeyboardInterrupt):
                    self.alignment._write_validated_outputs(
                        output_atom, output_sym,
                        self.alignment._render_atoms((atom,), cell),
                        self.alignment._render_symmetry(
                            (identity,), (0, 0, 0), (1, 1, 1)
                        ),
                        cell, (atom,), (identity,), Fraction(1, 10**12), False
                    )

            self.assertFalse(output_atom.exists())
            self.assertFalse(output_sym.exists())

    def test_cli_rejects_zero_fragment_count_without_traceback_or_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            atom_path = temporary / "atom.dat"
            sym_path = temporary / "sym.dat"
            output_atom = temporary / "aligned_atom.dat"
            output_sym = temporary / "aligned_sym.dat"
            atom_path.write_text("  'C' 0 0 0 1\n", encoding="ascii")
            sym_path.write_text(
                "1 0 0 0\n0 1 0 0\n0 0 1 0\n", encoding="ascii"
            )
            completed = subprocess.run(
                [
                    sys.executable, str(TOOL),
                    "--input-atom", str(atom_path), "--input-sym", str(sym_path),
                    "--output-atom", str(output_atom), "--output-sym", str(output_sym),
                    "--cell", "1", "1", "1", "--mesh", "4", "4", "4",
                    "--fragments", "0", "2", "2", "--buffer", "0", "0", "0",
                ],
                text=True, capture_output=True, check=False,
            )

            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("positive integers", completed.stderr)
            self.assertNotIn("Traceback", completed.stderr)
            self.assertFalse(output_atom.exists())
            self.assertFalse(output_sym.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
