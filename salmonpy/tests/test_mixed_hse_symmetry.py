import contextlib
import io
import unittest
from unittest.mock import patch

from salmonpy.examples.mpi_mixed_native import main


class MixedSymmetryAdmissionTests(unittest.TestCase):
    def rejected(self, extra, message):
        args = ['sample', 'build', 'runtime', 'bundle', 'output', '--layers', 'layout',
                '--fast-hse-build', 'worker', '--fast-hse-mpi-build', 'native', *extra]
        errors = io.StringIO()
        with patch('sys.argv', args), contextlib.redirect_stderr(errors):
            with self.assertRaises(SystemExit) as raised:
                main()
        self.assertEqual(raised.exception.code, 2)
        self.assertIn(message, errors.getvalue())

    def test_requires_explicit_full_k_reference(self):
        self.rejected(['--hse-symmetry'], '--hse-reference-bundle')

    def test_rejects_unvalidated_symmetry_inner_k2(self):
        self.rejected(['--hse-symmetry', '--hse-reference-bundle', 'full', '--hse-ranks', '2'], 'group1')

    def test_reference_does_not_implicitly_enable_symmetry(self):
        self.rejected(['--hse-reference-bundle', 'full'], 'requires --hse-symmetry')
