"""Velocity-gauge field/current integration tests against a native Si export."""
import os
from pathlib import Path
import tempfile
import unittest
import numpy as np
from model import NativeModel

EXPORT = Path(os.environ.get('SALMON_NATIVE_EXPORT_TEST', '/private/tmp/salmon-hse-export-check/export'))

class ExportValidationTests(unittest.TestCase):
    def test_missing_completion_marker_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory)/'metadata.txt').write_text('SALMON_HSE_REFERENCE_V1\n')
            with self.assertRaisesRegex(ValueError, 'incomplete|completion'):
                NativeModel(directory)

    def test_endianness_and_density_semantics_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path/'complete.txt').write_text('SALMON_HSE_REFERENCE_V1_COMPLETE\n')
            for line, error in [('endian_little invalid\n', 'endian'),
                                ('endian_little T\nrho_semantics incorrect\n', 'density|rho')]:
                (path/'metadata.txt').write_text('SALMON_HSE_REFERENCE_V1\n'+line)
                with self.assertRaisesRegex(ValueError, error):
                    NativeModel(directory)

@unittest.skipUnless((EXPORT/'complete.txt').exists(), 'Native export integration fixture not available')
class FieldCurrentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = NativeModel(EXPORT)

    def setUp(self):
        self.assertTrue(hasattr(self.model, 'set_field'), 'NativeModel.set_field is missing')
        self.assertTrue(hasattr(self.model, 'current'), 'NativeModel.current is missing')
        self.model.set_field(np.zeros(3))

    def test_zero_field_recovers_native_projectors_and_hamiltonian(self):
        model = self.model
        reference = model.read('projectors',model.shape+(model.nlma,model.nk),complex).transpose(4,3,0,1,2).reshape(model.nk,model.nlma,-1)
        model.set_field([.014,-.02,.008])
        self.assertGreater(np.linalg.norm(model.projectors-reference), 1e-5)
        model.set_field(np.zeros(3))
        np.testing.assert_allclose(model.projectors,reference,rtol=3e-14,atol=1e-15)
        self.assertLess(model.native_parity()['hpsi_relative_error'],1e-12)

    def test_current_matches_core_energy_derivative(self):
        model = self.model
        # Mix one occupied orbital with a conduction orbital with complex phase,
        # creating a coherent nonequilibrium state without changing its norm.
        u = model.psi[:,:16].copy()
        theta = .19
        u[:,0] = np.cos(theta)*model.psi[:,0]+1j*np.sin(theta)*model.psi[:,16]
        field = np.array([.017,-.013,.021])
        model.set_field(field)
        current = model.current(u)
        fd = []
        for axis in range(3):
            step=np.zeros(3);step[axis]=2e-5
            model.set_field(field+step);ep=model.expectation(u,model.core(u))
            model.set_field(field-step);em=model.expectation(u,model.core(u))
            fd.append((ep-em)/(2*step[axis])/(model.dv*np.prod(model.shape)))
        model.set_field(field)
        self.assertGreater(np.linalg.norm(current), 1e-5)
        np.testing.assert_allclose(current,fd,rtol=1e-7,atol=1e-10)
        # A repeat must not modify the state or accumulate gauge phases.
        model.set_field(field)
        np.testing.assert_allclose(model.current(u),current,rtol=2e-14,atol=1e-15)

    def test_anisotropic_or_nonorthogonal_geometry_rejected(self):
        geometry = np.fromfile(EXPORT/'geometry.bin')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            for source in EXPORT.iterdir():
                if source.name != 'geometry.bin':
                    (path/source.name).symlink_to(source)
            for index in (0, 4):
                changed = geometry.copy()
                changed[index] += .1
                changed.tofile(path/'geometry.bin')
                with self.assertRaisesRegex(ValueError, 'geometry|isotropic|orthogonal|cubic'):
                    NativeModel(path)

    def test_invalid_field(self):
        for value in ([0.,0.], [0.,0.,np.nan]):
            with self.assertRaises(ValueError):
                self.model.set_field(value)

if __name__ == '__main__':
    unittest.main()
