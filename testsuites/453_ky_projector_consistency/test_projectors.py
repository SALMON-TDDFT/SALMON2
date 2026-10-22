"""Production-object regression for scalar and equal-j KY projector consistency.

GNU/Makefiles CPU builds only; no model implementation replaces input_pp.
"""
import argparse
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
P = argparse.ArgumentParser()
P.add_argument('--build-dir', required=True, type=Path)
P.add_argument('--mpi')
P.add_argument('--np-flag', default='-n')
ARGS = P.parse_args()
BUILD = ARGS.build_dir.resolve()


def checked(command, **kwargs):
    result = subprocess.run(command, capture_output=True, text=True, **kwargs)
    if result.returncode:
        raise RuntimeError('command failed: '+shlex.join(map(str, command))+'\n'+result.stdout+result.stderr)
    return result


def fixture(path):
    # Equal plus/minus p,d channels must remain equal after all processing.
    h = .03
    lines = ['100 .03 2 4', 'KY_SOC_NLCC_V1',
             'ANGSTROM EV U_R_INV_SQRT_ANGSTROM RHO_INV_ANGSTROM3', '0 0 5',
             '0 1 .6 2 -1', '1 1 .6 0 -1', '1 3 .6 2 -1',
             '2 3 .6 0 -1', '2 5 .6 0 -1']
    for i in range(101):
        r = i*h
        row = [r, -2., 0., 0., 0.]
        for l in [0, 1, 1, 2, 2]:
            row += [-2.+math.exp(-r*r), r**(l+1)*math.exp(-r*r)]
        lines.append(' '.join('%.17g' % x for x in row))
    lines.append('END_KY_SOC_NLCC_V1')
    path.write_text('\n'.join(lines)+'\n')


def read_tables(path):
    lines = iter(path.read_text().splitlines())
    n, lm, radius = next(lines).split()
    channels = []
    for _ in range(int(lm)+1):
        meta = list(map(float, next(lines).split()))
        channels.append((meta, [list(map(float, next(lines).split())) for _ in range(int(n))]))
    return float(radius), channels


class Projectors(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        checked([sys.executable, str(ROOT/'tools/check_forbidden_eigensolvers.py'),
                 '--source-root', str(ROOT/'src')])
        cls.tmp = tempfile.TemporaryDirectory(prefix='ky-projector-')
        cls.directory = Path(cls.tmp.name)
        # Use the exact production link, replacing only the main program object.
        command = shlex.split((BUILD/'src/CMakeFiles/salmon.dir/link.txt').read_text())
        compiler = command[0]
        obj = cls.directory/'driver.o'
        checked([compiler, '-fopenmp', '-fcheck=all', '-I'+str(BUILD), '-c',
                 str(HERE/'driver.f90'), '-o', str(obj)], cwd=cls.directory)
        main = next(i for i, s in enumerate(command) if s.endswith('/main.f90.o'))
        command[main] = str(obj)
        cls.binary = cls.directory/'probe'
        command[command.index('-o')+1] = str(cls.binary)
        application_tree = BUILD/'src/CMakeFiles/salmon.dir'
        trees = {application_tree, cls.directory}
        # Also audit application objects supplied from outside the target tree.
        for argument in command[1:]:
            if argument.endswith(('.o', '.obj')):
                linked = (BUILD/'src'/argument).resolve()
                if application_tree not in linked.parents:
                    trees.add(linked.parent)
        for tree in sorted(trees):
            checked([sys.executable, str(ROOT/'tools/check_forbidden_eigensolvers.py'),
                     '--object-tree', str(tree)])
        checked(command, cwd=BUILD/'src')
        cls.equal = cls.directory/'equal_rps.dat'
        fixture(cls.equal)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_case(self, pseudo, mask, soc='y'):
        directory = self.directory/('run-'+str(len(list(self.directory.glob('run-*')))))
        directory.mkdir()
        command = [str(self.binary), str(pseudo), mask, soc]
        if ARGS.mpi:
            command = [ARGS.mpi, ARGS.np_flag, '2']+command
        env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1')
        checked(command, cwd=directory, env=env, timeout=90)
        tables = [read_tables(p) for p in sorted(directory.glob('tables.rank*'))]
        self.assertEqual(len(tables), 2 if ARGS.mpi else 1)
        self.assertTrue(all(x == tables[0] for x in tables), 'MPI broadcast mismatch')
        return tables[0]

    def assert_equal_j(self, channels):
        self.assertTrue(any(meta[1] != 0 for meta, _ in channels[1:]),
                        'fixture must have at least one nonlocal l>0 channel')
        for meta, rows in channels[1:]:
            self.assertEqual(meta[1], meta[2], 'minus-j projector sign missing')
            self.assertAlmostEqual(meta[3], meta[4], delta=1e-12*max(1, abs(meta[3])))
            for r, plus, minus, dplus, dminus in rows:
                self.assertAlmostEqual(plus, minus, delta=1e-11*max(1, abs(plus)), msg=str((meta[0], r)))
                self.assertAlmostEqual(dplus, dminus, delta=1e-10*max(1, abs(dplus)), msg=str((meta[0], r)))

    def test_legacy_scalar_is_spin_independent(self):
        pseudo = ROOT/'testsuites/pseudo/Si_rps.dat'
        for mask in ['n', 'y']:
            with self.subTest(mask=mask):
                scalar = self.run_case(pseudo, mask, 'n')
                spinor = self.run_case(pseudo, mask)
                self.assertEqual(scalar[0], spinor[0])
                self.assert_equal_j(spinor[1])
                for (m, a), (n, b) in zip(scalar[1], spinor[1]):
                    self.assertEqual(m[1], n[1])
                    self.assertEqual(m[3], n[3])
                    self.assertEqual([r[1::2] for r in a], [r[1::2] for r in b])

    def test_equal_j_channels_stay_equal_when_masked(self):
        self.assert_equal_j(self.run_case(self.equal, 'n')[1])
        self.assert_equal_j(self.run_case(self.equal, 'y')[1])

    def test_mask_radius_is_expanded_once(self):
        raw, _ = self.run_case(self.equal, 'n')
        radius, channels = self.run_case(self.equal, 'y')
        h = channels[0][1][1][0]-channels[0][1][0][0]
        # First radial point above gamma*rc, allowing rounding on an exact
        # boundary. A second expansion is many radial steps outside this band.
        self.assertGreaterEqual(radius, 1.8*raw-1e-12)
        self.assertLessEqual(radius, 1.8*raw+h+1e-12)


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]])
