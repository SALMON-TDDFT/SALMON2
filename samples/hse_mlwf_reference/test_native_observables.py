import unittest
import numpy as np
from compare_native_propagators import gauge_invariant_wavefunction_error


class GaugeComparisonTest(unittest.TestCase):
    def test_occupied_rotation_is_removed_but_excitation_is_not(self):
        rng = np.random.default_rng(199)
        u = np.zeros((4, 2, 2), complex)
        u[:2, :, :] = np.eye(2)[:, :, None]
        v = u.copy()
        for k in range(2):
            q, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))
            v[:, :, k] = u[:, :, k] @ q
        self.assertLess(gauge_invariant_wavefunction_error(u, v), 1e-14)
        v[2, 0, 0] += .01
        self.assertGreater(gauge_invariant_wavefunction_error(u, v), .001)
