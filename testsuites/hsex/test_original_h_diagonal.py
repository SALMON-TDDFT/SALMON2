"""Exact and adjacent diagonal poles: no false response from cancelling E and Hii."""
import os
import sys
from pathlib import Path
import numpy as np
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from q_openmp import QKernel
from precise_response import original_gram, original_residual

n = 64
band = 21
shift = 1e-7*np.random.default_rng(544).uniform(-1, 1, 216)[0]
d = np.linspace(-1e-4, 1e-4, n) + shift*np.r_[-np.ones(n//2), np.ones(n//2)]
h = np.diag(d).astype(complex)
psi = np.eye(n, dtype=complex)[:, band:band+1, None]
phase = np.ones((n, 3, 1), complex)
phase[:, 1, 0] = np.exp(.37j)
phase[:, 2, 0] = np.exp(.91j)
rhs = psi[:, :, 0]*phase[:, :, 0]
kernel = QKernel(os.environ['HSEX_Q_KERNEL'])
mp.mp.dps = 90
for eta in (1e-13, 5e-14):
    for energy in (d[band], np.nextafter(d[band], np.inf)):
        delta = mp.mpf(float(d[band]))-mp.mpf(float(energy))
        weight = float(delta/(delta*delta+mp.mpf(eta)**2))
        expected = weight*(rhs.conj().T@rhs)
        got, residual, team = kernel.solve(
            np.eye(n), d, np.zeros(n-1), np.empty((n, 0)), psi, phase, [3],
            np.array([[energy]]), np.array([[2.]]), [eta], 12, True, original_h=h)
        np.testing.assert_allclose(got[0][0], expected, rtol=2e-12, atol=2e-12)
        assert team == 12
        reference, _ = original_gram(h, rhs, energy, eta)
        np.testing.assert_allclose(reference, expected, rtol=2e-12, atol=2e-12)
        print('DIAGONAL_POLE eta', eta, 'delta', float(delta), 'native_max_error',
              np.max(abs(got[0][0]-expected)), 'python_max_error',
              np.max(abs(reference-expected)), 'residual', residual, flush=True)
# A rounded E-Hii would lose E entirely here; its low part must survive.
# Keep an imaginary diagonal as well, testing the residual helper's full H.
wide_h = np.diag([1e16+.25j, -1e16-.5j])
wide_x = np.ones((2, 1), complex)
wide_b = wide_h.diagonal().real[:, None].astype(complex)
wide_residual = original_residual(wide_h, 1., .01, wide_x, wide_b)
np.testing.assert_allclose(wide_residual[:, 0], [1.-.24j, 1.+.51j], atol=2e-16, rtol=0)
print('ORIGINAL_H_DIAGONAL_PASS')
