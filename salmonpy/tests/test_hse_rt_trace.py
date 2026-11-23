"""Successful comparison must prove that the intended fast path was executed."""
import importlib.util
from pathlib import Path
import unittest
import salmonpy

package=Path(salmonpy.__file__).resolve().parent
root=package/'_runtime_source' if (package/'_runtime_source').is_dir() else package.parent
spec=importlib.util.spec_from_file_location('hse_rt_trace',root/'tools/test_hse_rt_state_runtime.py')
trace=importlib.util.module_from_spec(spec);spec.loader.exec_module(trace)


class HSETraceTests(unittest.TestCase):
    def test_learned_taylor_requires_installed_prediction(self):
        with self.assertRaisesRegex(ValueError,'prediction'):
            trace._check_scientific_trace([],[],'hse_taylor4','learned')
        exact=['FACTOR_HISTORY step/predicted: 7 F']
        with self.assertRaisesRegex(ValueError,'prediction'):
            trace._check_scientific_trace(exact,exact,'hse_taylor4','learned')
        predicted=['FACTOR_HISTORY step/predicted: 7 T']
        self.assertEqual(trace._check_scientific_trace(predicted,predicted,'hse_taylor4','learned'),1)

    def test_ptcn_requires_solver_residual_trace(self):
        with self.assertRaisesRegex(ValueError,'PT-CN'):
            trace._check_scientific_trace([],[],'hse_ptcn','learned')
        residual=['HSE_PT_CN residual= 1e-12 builds=3 applications=10 Gram=1e-15 electrons=4']
        self.assertEqual(trace._check_scientific_trace(residual,residual,'hse_ptcn','learned'),0)

    def test_different_diagnostics_rejected_and_exact_path_allowed(self):
        with self.assertRaisesRegex(ValueError,'differ'):
            trace._check_scientific_trace(['one'],['two'],'hse_taylor4','off')
        self.assertEqual(trace._check_scientific_trace([],[],'hse_taylor4','off'),0)
