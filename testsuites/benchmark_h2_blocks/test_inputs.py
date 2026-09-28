"""Benchmark controls must not silently change the full reference or DC seeds."""
import importlib.util
from pathlib import Path
import unittest
spec=importlib.util.spec_from_file_location('block_benchmark',Path(__file__).with_name('run.py'))
bench=importlib.util.module_from_spec(spec);spec.loader.exec_module(bench)

class Inputs(unittest.TestCase):
    def test_default_unchanged(self):
        self.assertNotIn('exx_pair_screening',bench.rt_block((2,1,1),2,.999))
        self.assertNotIn('exx_pair_screening',bench.gs_input((2,1,1)))

    def test_screened_adaptive_only(self):
        full=bench.rt_block((2,1,1),2,1.,pair_tolerance=1e-6)
        adaptive=bench.rt_block((2,1,1),2,.999,pair_tolerance=1e-6)
        self.assertNotIn('exx_pair_screening',full)
        self.assertIn("exx_pair_screening='on'",adaptive)
        self.assertIn('exx_pair_tolerance=9.99999999999999955d-07',adaptive)
        self.assertIn('nt=16',adaptive)

    def test_invalid_budget(self):
        for value in (-1,float('nan'),float('inf')):
            with self.assertRaises(ValueError):bench.rt_block((2,1,1),2,.999,pair_tolerance=value)

if __name__=='__main__':unittest.main()
