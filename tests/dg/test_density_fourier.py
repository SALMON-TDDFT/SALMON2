import unittest
import importlib.util
from pathlib import Path
import numpy as np


class FourierTest(unittest.TestCase):
    def test_known_modes(self):
        path = Path(__file__).with_name('analyze_density_fourier.py')
        self.assertTrue(path.exists(), 'Fourier analyzer missing')
        spec = importlib.util.spec_from_file_location('fourier', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        x = np.arange(16)[:, None, None]
        delta = np.broadcast_to(np.cos(2*np.pi*x/16) + 2*np.cos(2*np.pi*5*x/16), (16, 16, 16))
        result = module.spectrum(delta, 10.0)
        self.assertAlmostEqual(result['mean_square'], 2.5)
        self.assertAlmostEqual(result['parseval_relative_defect'], 0, places=14)
        self.assertAlmostEqual(result['long_fraction'], 0.2)
        self.assertAlmostEqual(result['middle_fraction'], 0)
        self.assertAlmostEqual(result['short_fraction'], 0.8)
        self.assertAlmostEqual(result['zero_fraction'], 0)
        zero = module.spectrum(np.zeros((16, 16, 16)), 10.0)
        self.assertEqual(zero['long_fraction'], 0)
        self.assertAlmostEqual(module.spectrum(np.ones((16, 16, 16)), 10.0)['zero_fraction'], 1)


if __name__ == '__main__':
    unittest.main()
