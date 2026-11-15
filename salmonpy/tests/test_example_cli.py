"""MPI help must be usable before installing/initializing MPI."""
from pathlib import Path
import subprocess
import sys
import unittest


class ExampleCLITests(unittest.TestCase):
    def test_mpi_help_without_mpi(self):
        for name in ('mpi_response', 'mpi_tddft', 'mpi_coupled_tddft',
                     'mpi_multiscale_tddft', 'mpi_recorded_material',
                     'mpi_configured_native'):
            with self.subTest(name=name):
                # Block MPI even on developer machines where it is installed.
                script = '''import runpy, sys
sys.modules['mpi4py'] = None
sys.argv = [sys.argv[1], '--help']
runpy.run_module(sys.argv[0], run_name='__main__')
'''
                result = subprocess.run(
                    [sys.executable, '-c', script, 'salmonpy.examples.' + name],
                    capture_output=True, text=True, timeout=20)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)
                self.assertNotIn('Traceback', result.stderr)
