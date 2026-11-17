import contextlib
import io
import unittest
from types import SimpleNamespace
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

    def test_accepts_inner_k2_without_starting_real_mpi(self):
        class MPIBoundary:
            @property
            def COMM_WORLD(self):
                raise RuntimeError('MPI boundary reached')
        args = ['sample', 'build', 'runtime', 'bundle', 'output', '--layers', 'layout',
                '--fast-hse-build', 'worker', '--fast-hse-mpi-build', 'native',
                '--hse-symmetry', '--hse-reference-bundle', 'full', '--hse-ranks', '2']
        with patch('sys.argv', args), patch.dict('sys.modules', {'mpi4py': SimpleNamespace(MPI=MPIBoundary())}):
            with self.assertRaisesRegex(RuntimeError, 'MPI boundary reached'):
                main()

    def test_reference_does_not_implicitly_enable_symmetry(self):
        self.rejected(['--hse-reference-bundle', 'full'], 'requires --hse-symmetry')
