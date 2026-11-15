#!/usr/bin/env python3

import gzip
import math
from pathlib import Path
import sys
import tempfile
import unittest


CIF2SALMON_DIRECTORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CIF2SALMON_DIRECTORY))

import cif2salmon


def cif_text(atom_rows, extra_header="_space_group_name_H-M_alt 'P 1'", extra_tags=""):
    return """data_test
_cell_length_a 4.0(1)
_cell_length_b 5.0
_cell_length_c 6.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
{}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
{}{}
""".format(extra_header, extra_tags, atom_rows)


class Cif2SalmonTest(unittest.TestCase):
    def read_text(self, text, block=None, tolerance=1.0e-5):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.cif"
            path.write_text(text, encoding="utf-8")
            return cif2salmon.read_cif_structure(path, block, tolerance)

    def test_p1_fractional_structure_and_rendering(self):
        structure = self.read_text(cif_text("Si1 Si 0.125(2) 0.25 0.5\n"))

        self.assertEqual(structure.block_name, "test")
        self.assertEqual(len(structure.atoms), 1)
        self.assertEqual(structure.atoms[0].symbol, "Si")
        self.assertEqual(structure.atoms[0].fractional, (0.125, 0.25, 0.5))
        output = cif2salmon.render_salmon(structure, "/some/path/input.cif")
        self.assertIn("unit_system = 'A_eV_fs'", output)
        self.assertIn("al = 4.000000000000, 5.000000000000, 6.000000000000", output)
        self.assertNotIn("al_vec", output)
        self.assertIn("natom = 1", output)
        self.assertIn("! species index: 1 = Si (Z=14), count = 1", output)
        self.assertIn("! total atom count: natom = 1", output)
        self.assertIn("'Si'  0.125000000000", output)
        self.assertNotIn("/some/path", output)

    def test_species_comments_include_each_count_and_total_natom(self):
        structure = self.read_text(
            cif_text(
                "O1 O 0 0 0\n"
                "Si1 Si 0.25 0.25 0.25\n"
                "O2 O 0.5 0.5 0.5\n"
            )
        )
        output = cif2salmon.render_salmon(structure, "input.cif")

        self.assertIn("! species index: 1 = O (Z=8), count = 2", output)
        self.assertIn("! species index: 2 = Si (Z=14), count = 1", output)
        self.assertIn("! total atom count: natom = 3", output)

    def test_spacegroup_number_expands_and_deduplicates_special_position(self):
        structure = self.read_text(
            cif_text(
                "C1 C 0.1 0.2 0.3\nC2 C 0 0 0\n",
                extra_header="_space_group_IT_number 2",
            )
        )

        self.assertEqual(len(structure.atoms), 3)
        positions = {atom.fractional for atom in structure.atoms}
        self.assertIn((0.1, 0.2, 0.3), positions)
        self.assertIn((0.9, 0.8, 0.7), positions)
        self.assertIn((0.0, 0.0, 0.0), positions)

    def test_cartesian_coordinates_are_converted_to_fractional(self):
        text = """data_cartesian
_cell_length_a 4
_cell_length_b 5
_cell_length_c 6
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_Cartn_x
_atom_site_Cartn_y
_atom_site_Cartn_z
C1 C 2 2.5 3
"""
        structure = self.read_text(text)
        self.assertEqual(structure.atoms[0].fractional, (0.5, 0.5, 0.5))

    def test_triclinic_cell_vectors_use_standard_geometry(self):
        text = cif_text("C1 C 0 0 0\n").replace(
            "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90",
            "_cell_angle_alpha 70\n_cell_angle_beta 80\n_cell_angle_gamma 75",
        )
        structure = self.read_text(text)
        a, b, c = structure.vectors

        self.assertAlmostEqual(math.sqrt(sum(x * x for x in a)), 4.0)
        self.assertAlmostEqual(math.sqrt(sum(x * x for x in b)), 5.0)
        self.assertAlmostEqual(math.sqrt(sum(x * x for x in c)), 6.0)
        self.assertAlmostEqual(
            sum(a[i] * b[i] for i in range(3)) / 20.0,
            math.cos(math.radians(75)),
        )
        self.assertAlmostEqual(
            sum(a[i] * c[i] for i in range(3)) / 24.0,
            math.cos(math.radians(80)),
        )
        self.assertAlmostEqual(
            sum(b[i] * c[i] for i in range(3)) / 30.0,
            math.cos(math.radians(70)),
        )
        output = cif2salmon.render_salmon(structure, "triclinic.cif")
        self.assertIn("al_vec1", output)
        self.assertIn("al_vec2", output)
        self.assertIn("al_vec3", output)

    def test_unknown_and_partial_occupancies_are_rejected(self):
        for occupancy in ("?", ".", "0.5"):
            with self.subTest(occupancy=occupancy):
                text = cif_text(
                    "C1 C 0 0 0 {}\n".format(occupancy),
                    extra_tags="_atom_site_occupancy\n",
                )
                with self.assertRaisesRegex(cif2salmon.ConversionError, "occupancy"):
                    self.read_text(text)

    def test_missing_spacegroup_is_not_silently_treated_as_p1(self):
        with self.assertRaisesRegex(cif2salmon.ConversionError, "silently assume P1"):
            self.read_text(cif_text("C1 C 0 0 0\n", extra_header=""))

    def test_multiple_structure_blocks_require_selection(self):
        first = cif_text("C1 C 0 0 0\n").replace("data_test", "data_first")
        second = cif_text("O1 O 0 0 0\n").replace("data_test", "data_second")
        with self.assertRaisesRegex(cif2salmon.ConversionError, "multiple structure blocks"):
            self.read_text(first + second)

        structure = self.read_text(first + second, block="data_second")
        self.assertEqual(structure.block_name, "second")
        self.assertEqual(structure.atoms[0].symbol, "O")

    def test_atoms_format_omits_companion_namelists(self):
        structure = self.read_text(cif_text("C1 C 0 0 0\n"))
        output = cif2salmon.render_salmon(structure, "input.cif", atoms_only=True)

        self.assertNotIn("&units", output)
        self.assertNotIn("&system", output)
        self.assertIn("&atomic_red_coor", output)

    def test_cartesian_output_uses_atomic_coor_in_angstrom(self):
        structure = self.read_text(cif_text("Si1 Si 0.125 0.25 0.5\n"))
        output = cif2salmon.render_salmon(
            structure,
            "input.cif",
            coordinate_system="cartesian",
        )

        self.assertIn("&atomic_coor", output)
        self.assertNotIn("&atomic_red_coor", output)
        self.assertIn("Cartesian atomic coordinates are in angstrom", output)
        self.assertIn(
            "'Si'  0.500000000000  1.250000000000  3.000000000000  1",
            output,
        )

    def test_unknown_output_coordinate_system_is_rejected(self):
        structure = self.read_text(cif_text("C1 C 0 0 0\n"))
        with self.assertRaisesRegex(cif2salmon.ConversionError, "coordinate system"):
            cif2salmon.render_salmon(
                structure,
                "input.cif",
                coordinate_system="spherical",
            )

    def test_legacy_repository_sample_expands_to_eight_atoms(self):
        sample = CIF2SALMON_DIRECTORY.parent / "sym" / "Si_sample.cif"
        structure = cif2salmon.read_cif_structure(sample)

        self.assertEqual(len(structure.atoms), 8)
        self.assertEqual(structure.species, (("Si", 14),))
        self.assertTrue(
            all(0.0 <= value < 1.0 for atom in structure.atoms for value in atom.fractional)
        )

    def test_gzip_compressed_cif_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.cif.gz"
            with gzip.open(path, "wt", encoding="utf-8") as stream:
                stream.write(cif_text("C1 C 0 0 0\n"))

            structure = cif2salmon.read_cif_structure(path)

        self.assertEqual(len(structure.atoms), 1)
        self.assertEqual(structure.atoms[0].symbol, "C")


if __name__ == "__main__":
    unittest.main()
