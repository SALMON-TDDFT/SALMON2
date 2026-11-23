import tempfile
from pathlib import Path
import unittest

from salmonpy.errors import InputError


IDENTITY = '1 0 0 0\n0 1 0 0\n0 0 1 0\n'
C2Z = '-1 0 0 0\n0 -1 0 0\n0 0 1 0\n'


class SymmetryInputAdmissionTests(unittest.TestCase):
    def check(self, text, setting='yyn'):
        from salmonpy.fast_hse_mpi import _symmetry_preflight
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'sym.dat'
            path.write_text(text)
            before=path.read_bytes()
            try:
                _symmetry_preflight(path,setting)
            finally:
                self.assertEqual(path.read_bytes(),before)

    def test_valid_group_and_fortran_numbers(self):
        self.check(IDENTITY+C2Z)
        self.check(IDENTITY+IDENTITY.replace('1 0 0 0','1 0 0 .5'))
        self.check('1d0, 2*0, 0 ! first row\n0 1D0 0 0\n2*0 1 0\n'+C2Z)

    def test_malformed_records_and_nonfinite(self):
        for text in ['',IDENTITY+'1 0 0 0\n',IDENTITY.replace('1 0 0 0','1 0 0'),
                     IDENTITY.replace('1 0 0 0','NaN 0 0 0'),
                     IDENTITY.replace('1 0 0 0','1e999 0 0 0'),
                     IDENTITY.replace('1 0 0 0','1_0 0 0 0'),
                     IDENTITY.replace('1 0 0 0','0_1*1 0 0 0'),
                     IDENTITY.replace('1 0 0 0','+1*1 0 0 0')]:
            with self.subTest(text=text),self.assertRaisesRegex(InputError,'sym.dat'):
                self.check(text)

    def test_group_guards(self):
        quarter='0 -1 0 0\n1 0 0 0\n0 0 1 0\n'
        for text in [C2Z,IDENTITY+IDENTITY,IDENTITY+quarter,
                     IDENTITY.replace('1 0 0 0','2 0 0 0')]:
            with self.subTest(text=text),self.assertRaisesRegex(InputError,'sym.dat'):
                self.check(text)

    def test_filter_matches_disabled_axes_before_group_validation(self):
        self.check(IDENTITY+C2Z,setting='nnn')
        with self.assertRaisesRegex(InputError,'retained'):
            self.check(C2Z,setting='nnn')

    def test_duplicate_translations_modulo_cell(self):
        with self.assertRaisesRegex(InputError,'duplicate'):
            self.check(IDENTITY+IDENTITY.replace('1 0 0 0','1 0 0 1'))

    def test_invalid_axis_flags_are_not_silently_normalized(self):
        for flags in ['Yyn','yyyy','y?','']:
            with self.subTest(flags=flags),self.assertRaisesRegex(InputError,'yn_symmetry'):
                self.check(IDENTITY,setting=flags)
