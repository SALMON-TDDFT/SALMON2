"""Real Fortran vector/density regressions against independent Cartesian/SU(2) references.

Run with Python 3 + NumPy and a GNU-compatible Fortran compiler (FC, default gfortran).
Only the tiny fixtures execute; the full application source first passes the ban checker.
"""
from pathlib import Path
import itertools
import os
import shlex
import subprocess
import sys
import tempfile
import unittest
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def models():
    """Direct primitive bases, including improper cubic and Se screw operations."""
    W = np.array([[-1., 1., 1.], [1., -1., 1.], [1., 1., -1.]])
    cubic = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([-1, 1], repeat=3):
            R = np.zeros((3, 3)); R[range(3), perm] = signs
            if np.array_equal(R @ [0, 0, 1], [0, 0, 1]):
                cubic.append(np.column_stack((np.linalg.solve(W, R @ W), np.zeros(3))))
    Se = np.array([[.5, .5, 0], [-np.sqrt(3)/2, np.sqrt(3)/2, 0], [0, 0, 1.2]])
    C = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
    screw = np.array([np.column_stack((np.linalg.matrix_power(C, j), [0, 0, j/3])) for j in range(3)])
    C2 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    transverse = np.array([np.column_stack((M, np.zeros(3))) for M in [np.eye(3), C2]])
    return [('bcc-z', W, np.array(cubic), (4, 4, 4)),
            ('Se-screw-z', Se, screw, (4, 4, 6)),
            ('Se-C2-x', Se, transverse, (4, 4, 6)),
            ('orthogonal-z', np.eye(3), np.array([np.column_stack((R, np.zeros(3))) for R in
              [np.eye(3), np.diag([-1., -1., 1.])]]), (4, 4, 4))]


def data(A, ops):
    inverse = np.linalg.inv(A)
    B = 2*np.pi*inverse.T
    cart = np.array([A @ M[:, :3] @ inverse for M in ops])
    reciprocal = ops.copy()
    for i, R in enumerate(cart):
        reciprocal[i, :, :3] = np.linalg.solve(B, R @ B)
    return inverse, cart, reciprocal


def encode(n, arrays):
    return '\n'.join([str(n)] + [str(x) for a in arrays for x in np.asarray(a).flatten(order='F')]) + '\n'


def su2(Q):
    """Independent spinor rotation; compare densities, which ignore the double-cover sign."""
    angle = np.arccos(np.clip((np.trace(Q)-1)/2, -1, 1))
    if abs(angle) < 1e-10:
        return np.eye(2, dtype=complex)
    values, vectors = np.linalg.eig(Q)
    axis = np.real(vectors[:, np.argmin(abs(values-1))]); axis /= np.linalg.norm(axis)
    if abs(np.sin(angle)) > 1e-8:
        skew = np.array([Q[2, 1]-Q[1, 2], Q[0, 2]-Q[2, 0], Q[1, 0]-Q[0, 1]])
        if axis @ skew < 0:
            axis = -axis
    pauli = np.array([[axis[2], axis[0]-1j*axis[1]], [axis[0]+1j*axis[1], -axis[2]]])
    return np.cos(angle/2)*np.eye(2) - 1j*np.sin(angle/2)*pauli


class Covariance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT/'tools/check_forbidden_eigensolvers.py'),
                        '--source-root', str(ROOT/'src')], check=True)
        cls.temp = tempfile.TemporaryDirectory()
        cls.build = Path(cls.temp.name)
        compiler = shlex.split(os.environ.get('FC', 'gfortran'))
        for target, source, driver in [('vector', 'symmetry_vector.f90', 'driver.f90'),
                                       ('density', 'symmetry_spin_density.f90', 'density_driver.f90')]:
            subprocess.run(compiler + ['-fcheck=all', str(HERE/'sym_stub.f90'),
              str(ROOT/'src/symmetry'/source), str(HERE/driver), '-o', str(cls.build/target)],
              cwd=cls.build, check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_cartesian_current_in_skew_and_rotated_bases(self):
        for name, base, ops, _ in models():
            for Q in [np.eye(3), np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])]:
                A = Q @ base
                inverse, cart, reciprocal = data(A, ops)
                for v in [np.array([0., 0., 1.]), np.array([1., 2., 3.])]:
                    with self.subTest(model=name, rotated=bool(Q[0, 2]), vector=v):
                        text = encode(len(ops), [A, inverse, ops.transpose(1, 2, 0),
                                                reciprocal.transpose(1, 2, 0), v])
                        got = np.fromstring(subprocess.check_output([str(self.build/'vector')], input=text, text=True), sep=' ')
                        np.testing.assert_allclose(got, np.mean([R @ v for R in cart], axis=0), atol=1e-12, rtol=0)

    def run_density(self, A, ops, g, D, check=True):
        inverse, _, reciprocal = data(A, ops)
        pairs = np.array([[z.real, z.imag] for z in D.flatten(order='F')]).T
        text = encode(len(ops), [g, A, inverse, ops.transpose(1, 2, 0),
                                reciprocal.transpose(1, 2, 0), pairs])
        return subprocess.run([str(self.build/'density')], input=text, text=True, capture_output=True, check=check)

    def test_spin_density_matches_explicit_rotated_spinors(self):
        rng = np.random.default_rng(942)
        for name, A, ops, dims in models():
            with self.subTest(model=name):
                g = np.array(dims); shape = dims+(2, 2)
                spinors = rng.normal(size=dims+(2,)) + 1j*rng.normal(size=dims+(2,))
                D = np.einsum('...i,...j->...ij', spinors, spinors.conj())
                expected = np.zeros_like(D)
                _, cart, _ = data(A, ops)
                for M, R in zip(ops, cart):
                    U = su2(round(np.linalg.det(R))*R)
                    for ijk in np.ndindex(*g):
                        point = np.linalg.solve(M[:, :3], np.array(ijk)/g-M[:, 3])*g
                        np.testing.assert_allclose(point, np.rint(point), atol=1e-12)
                        source = tuple(np.rint(point).astype(int) % g)
                        v = U @ spinors[source]
                        expected[ijk] += np.outer(v, v.conj())/len(ops)
                values = np.fromstring(self.run_density(A, ops, g, D).stdout, sep=' ')
                got = (values[::2]+1j*values[1::2]).reshape(shape, order='F')
                np.testing.assert_allclose(got, expected, atol=1e-12, rtol=0)
                self.assertLess(abs(np.trace(got, axis1=3, axis2=4).sum()-np.trace(D, axis1=3, axis2=4).sum()), 1e-10)

    def test_noncommensurate_screw_grid_rejected(self):
        _, A, ops, _ = models()[1]; g = np.array([4, 4, 5])
        D = np.zeros(tuple(g)+(2, 2), dtype=complex)
        result = self.run_density(A, ops, g, D, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('noncommensurate grid', result.stderr)

    def test_nonfinite_translation_rejected(self):
        _, A, ops, dims = models()[2]; ops[1, 0, 3] = np.nan
        result = self.run_density(A, ops, np.array(dims), np.zeros(dims+(2, 2), complex), check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('nonfinite operation', result.stderr)

    def test_nonhermitian_density_rejected(self):
        _, A, ops, dims = models()[0]; D = np.zeros(dims+(2, 2), complex); D[..., 0, 1] = 1
        result = self.run_density(A, ops, np.array(dims), D, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('non-Hermitian', result.stderr)

if __name__ == '__main__':
    unittest.main()
