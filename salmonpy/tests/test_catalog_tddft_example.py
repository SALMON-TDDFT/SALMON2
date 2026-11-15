import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from salmonpy.errors import InputError
from salmonpy.examples.catalog_tddft import run


class CatalogExampleTests(unittest.TestCase):
    def test_nonreal_waveforms_rejected_before_runtime_loading(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for A, dt in ((np.zeros((2, 3), complex), 0.1),
                          (np.zeros((2, 3), bool), 0.1),
                          (np.zeros((2, 3)), True),
                          (np.zeros((2, 3)), '0.1')):
                np.savez(root/'wave.npz', A=A, dt=dt)
                with self.subTest(dtype=A.dtype, dt=dt), patch(
                        'salmonpy.examples.catalog_tddft.load_runtime',
                        side_effect=AssertionError('must not load runtime')):
                    with self.assertRaises(InputError):
                        run(root/'materials.json', 'ground', root/'runtime',
                            root/'wave.npz', root/'output')
            self.assertFalse((root/'output').exists())
