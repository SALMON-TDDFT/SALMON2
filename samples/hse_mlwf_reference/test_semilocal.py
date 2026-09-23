import importlib.util
import unittest
import numpy as np

class SemilocalTests(unittest.TestCase):
    def wrapper(self, name):
        self.assertIsNotNone(importlib.util.find_spec('semilocal'), 'Libxc semilocal wrapper is missing')
        from semilocal import Semilocal
        return Semilocal(name)

    def test_hse_coefficients(self):
        with self.wrapper('hse06') as xc:
            self.assertAlmostEqual(xc.coefficients['omega'], .11)
            self.assertAlmostEqual(xc.coefficients['short_range_exact'], .25)
            self.assertAlmostEqual(xc.coefficients['long_range_exact'], 0.)

    def test_density_and_gradient_derivatives(self):
        rho = np.array([.001, .02, .3])
        sigma = np.array([1e-7, .0003, .04])
        for name in ('hse06', 'pbe'):
            with self.wrapper(name) as xc:
                eps, vrho, vsigma = xc.evaluate(rho, sigma)
                dr = rho * 1e-5
                ds = sigma * 1e-4
                ep = xc.evaluate(rho + dr, sigma)[0] * (rho + dr)
                em = xc.evaluate(rho - dr, sigma)[0] * (rho - dr)
                np.testing.assert_allclose((ep-em)/(2*dr), vrho, rtol=2e-6, atol=1e-8)
                ep = xc.evaluate(rho, sigma + ds)[0] * rho
                em = xc.evaluate(rho, sigma - ds)[0] * rho
                np.testing.assert_allclose((ep-em)/(2*ds), vsigma, rtol=2e-5, atol=1e-7)

    def test_hse_independent_component_weighting(self):
        # Libxc7 HSE uses WPBEH(omega=0) as its full-range exchange,
        # not GGA_X_PBE(101): the exchange-hole model differs at finite gradient.
        import ctypes as ct
        rho = np.array([.001, .01, .1, 1.])
        sigma = np.array([1e-7, .0003, .04, .2])
        with self.wrapper('hse06') as xc:
            lib = xc.lib
            lib.xc_func_set_ext_params_name.argtypes = [ct.c_void_p, ct.c_char_p, ct.c_double]
            lib.xc_func_set_ext_params_name.restype = None
            expected = np.zeros((3, 4))
            ptr = ct.POINTER(ct.c_double)
            for identifier, omega, weight in ((524, 0., 1.), (524, .11, -.25), (130, None, 1.)):
                f = lib.xc_func_alloc()
                self.assertEqual(lib.xc_func_init(f, identifier, 1), 0)
                try:
                    if omega is not None:
                        lib.xc_func_set_ext_params_name(f, b'_omega', omega)
                    values = np.empty((3, 4))
                    lib.xc_gga_exc_vxc(f, 4, *[a.ctypes.data_as(ptr) for a in (rho, sigma, *values)])
                    expected += weight * values
                finally:
                    lib.xc_func_end(f)
                    lib.xc_func_free(f)
            np.testing.assert_allclose(xc.evaluate(rho, sigma), expected, rtol=2e-13, atol=1e-14)

    def test_vacuum_and_low_density(self):
        for name in ('hse06', 'pbe'):
            with self.wrapper(name) as xc:
                result = xc.evaluate(np.array([0., 1e-30, 1e-12, .1]), np.zeros(4))
                for component in result:
                    self.assertTrue(np.isfinite(component).all())
                    self.assertEqual(component[0], 0.)
                for rho, sigma in (([-1.], [0.]), ([1.], [-1.]), ([np.nan], [0.]), ([1., 2.], [0.])):
                    with self.assertRaises(ValueError):
                        xc.evaluate(rho, sigma)

    def test_pbe_uniform_exchange_correlation_and_hse_weighting(self):
        # At homogeneous density, HSE remainder must be less negative than PBE,
        # since its short-range PBE exchange portion was already subtracted.
        rho = np.array([.01, .1, 1.])
        with self.wrapper('pbe') as pbe, self.wrapper('hse06') as hse:
            ep = pbe.evaluate(rho, np.zeros(3))[0]
            eh = hse.evaluate(rho, np.zeros(3))[0]
            self.assertTrue(np.all(ep < eh))
            self.assertTrue(np.all(eh < 0.))

if __name__ == '__main__':
    unittest.main()
