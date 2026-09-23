"""Numerical and lifetime checks for cached single-thread FFTW convolution."""
import importlib.util
import unittest
import numpy as np
from exchange import screened_multiplier

class FFTWTests(unittest.TestCase):
    def backend(self, multiplier):
        self.assertIsNotNone(importlib.util.find_spec('fftw_backend'), 'FFTW backend is missing')
        from fftw_backend import FFTWConvolution
        return FFTWConvolution(multiplier)

    def test_complex_noncubic_convolution(self):
        rng = np.random.default_rng(109)
        shape = (9, 8, 7)
        multiplier = rng.normal(size=shape) + 1j*rng.normal(size=shape)
        rho = rng.normal(size=shape) + 1j*rng.normal(size=shape)
        original = rho.copy()
        with self.backend(multiplier) as kernel:
            actual = kernel.convolve(rho)
            expected = np.fft.ifftn(np.fft.fftn(rho)*multiplier)
            np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=3e-13)
            np.testing.assert_array_equal(rho, original)
            self.assertEqual(kernel.shape, shape)
            self.assertGreaterEqual(kernel.plan_seconds, 0.)

    def test_screened_kernel_including_qzero(self):
        shape = (12, 12, 12)
        multiplier = screened_multiplier(shape, .855, .11)
        with self.backend(multiplier) as kernel:
            np.testing.assert_allclose(kernel.convolve(np.ones(shape)), np.pi/.11**2, rtol=1e-13)
            rng = np.random.default_rng(14)
            rho = rng.normal(size=shape) + 1j*rng.normal(size=shape)
            np.testing.assert_allclose(kernel.convolve(rho),
                np.fft.ifftn(np.fft.fftn(rho)*multiplier), rtol=2e-13, atol=2e-13)

    def test_repeated_calls_are_independent_and_close_rejects(self):
        shape = (6, 5, 4)
        with self.backend(np.ones(shape)) as kernel:
            first = kernel.convolve(np.ones(shape))
            saved = first.copy()
            second = kernel.convolve(np.full(shape, 3j))
            np.testing.assert_array_equal(first, saved)
            self.assertFalse(np.shares_memory(first, second))
            np.testing.assert_allclose(second, 3j, atol=1e-14)
            with self.assertRaises(ValueError):
                kernel.convolve(np.zeros((2, 3, 4)))
        kernel.close()
        with self.assertRaises(RuntimeError):
            kernel.convolve(np.zeros(shape))

    def test_invalid_multiplier(self):
        for bad in (np.ones((3, 3)), np.ones((0, 3, 3)), np.full((3, 3, 3), np.nan)):
            with self.assertRaises(ValueError):
                with self.backend(bad):
                    pass

if __name__ == '__main__':
    unittest.main()
