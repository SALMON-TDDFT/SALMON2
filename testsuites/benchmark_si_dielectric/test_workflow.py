import unittest
from pathlib import Path
from run import make_input


class Inputs(unittest.TestCase):
    def test_k4_common_conditions(self):
        for xc in ('pbe', 'hse06', 'pbe0', 'pbeh40'):
            gs = make_input(xc, 'gs', Path('/tmp/producer'))
            rt = make_input(xc, 'rt', Path('/tmp/producer'))
            for text in (gs, rt):
                self.assertIn('num_kgrid=4,4,4', text)
                self.assertIn('num_rgrid=16,16,16', text)
                self.assertIn('natom=8', text)
                self.assertIn(f"xc='{"libxc_pbe" if xc == "pbe" else xc}'", text)
                self.assertNotIn('temperature', text)
                self.assertNotIn("yn_dc='y'", text)
            self.assertIn('e_impulse=1d-4', rt)
            self.assertIn('nt=4375', rt)
            self.assertNotIn('exx_pre_scf_threshold', rt)
            self.assertIn("directory_read_data='/tmp/producer/data_for_restart/'", rt)

    def test_equal_time_probes(self):
        a = make_input('pbe', 'rt', Path('/tmp/gs'), .08, 16)
        b = make_input('pbe', 'rt', Path('/tmp/gs'), .04, 32)
        self.assertIn('nt=16', a)
        self.assertIn('nt=32', b)
        self.assertIn('dt=0.08', a)
        self.assertIn('dt=0.04', b)


if __name__ == '__main__':
    unittest.main()
