import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('compare_spin',ROOT/'tools/compare_spin.py')
mod=importlib.util.module_from_spec(spec)

class CompareSpin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not spec.origin or not Path(spec.origin).exists():
            raise AssertionError('magnetization comparison tool is missing')
        spec.loader.exec_module(mod)

    def load(self, text):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'spin.data'
            path.write_text('# 1:Time[fs] 2:mag_x[none] 3:mag_y[none] 4:mag_z[none]\n'+text)
            return mod.read_spin(path)

    def test_rotation_does_not_look_like_demagnetization(self):
        run=self.load('0 0 0 1\n1 1 0 0\n')
        result=mod.compare({'LSDA':run,'MGGA':run})
        self.assertEqual(result['runs']['LSDA']['magnitude'],[1.,1.])
        self.assertEqual(result['runs']['LSDA']['delta_magnitude'],[0.,0.])

    def test_zero_initial_magnetization_is_not_divided(self):
        run=self.load('0 0 0 0\n1 0 0 0.1\n')
        self.assertIsNone(mod.compare({'AFM':run})['runs']['AFM']['normalized_magnitude'])

    def test_bad_data_and_time_mismatch_rejected(self):
        for rows in ('0 nan 0 1\n','0 0 0 1\n0 0 0 1\n','1 0 0 1\n0 0 0 1\n'):
            with self.assertRaises(ValueError): self.load(rows)
        a=self.load('0 0 0 1\n1 0 0 .9\n')
        b=self.load('0 0 0 1\n2 0 0 .9\n')
        with self.assertRaises(ValueError): mod.compare({'A':a,'B':b})

if __name__=='__main__': unittest.main()
