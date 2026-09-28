"""Actual communication wrappers, MPI decomposition parity and retained-norm probes."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class AdaptiveSupportTest(unittest.TestCase):
    def test_periodic_norm_support(self):
        module = ROOT / 'src/xc/exx_adaptive_support.f90'
        self.assertTrue(module.exists(), 'adaptive source support is not implemented')
        with tempfile.TemporaryDirectory() as tmp:
            exe = Path(tmp) / 'probe'
            sources = [ROOT / 'src/misc/nvtx_wrapper.f90',
                       ROOT / 'src/parallel/communication.f90', module,
                       Path(__file__).with_name('adaptive_support_probe.f90')]
            build = subprocess.run([os.environ.get('MPIFC', 'mpifort'), '-cpp', '-O0', '-g',
                                    '-fcheck=all', '-ffree-line-length-none', '-fallow-argument-mismatch',
                                    *map(str, sources), '-o', str(exe)],
                                   cwd=tmp, text=True, capture_output=True)
            self.assertEqual(build.returncode, 0, build.stderr)
            reference = None
            for ranks in (1, 2, 4):
                run = subprocess.run(['mpiexec', '-n', str(ranks), str(exe)], cwd=tmp,
                                     text=True, capture_output=True,
                                     env=dict(os.environ, OMP_NUM_THREADS='1'))
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                values = list(map(float, run.stdout.split()))
                self.assertEqual(len(values), 8)
                if reference is not None:
                    for value, expected in zip(values, reference):
                        self.assertAlmostEqual(value, expected, delta=1e-10)
                reference = values


if __name__ == '__main__':
    unittest.main()
