import importlib
import json
import subprocess
import sys
import unittest
from pathlib import Path

class PublicAPI(unittest.TestCase):
    def test_version_and_inventory(self):
        import salmonpy
        from salmonpy.api import PUBLIC_API
        self.assertEqual(salmonpy.__version__, '0.1.0a5')
        for module, names in PUBLIC_API.items():
            if module == 'salmonpy.ase_calculator' and importlib.util.find_spec('ase') is None:
                continue
            imported = importlib.import_module(module)
            for name in names:
                self.assertTrue(hasattr(imported, name), (module, name))
        self.assertIn('__version__', salmonpy.__all__)

    def test_errors_keep_builtin_catches(self):
        from salmonpy.errors import SalmonpyError, InputError, InputTypeError, NumericalError, RuntimeFailure
        from salmonpy.ld import LorentzDrude
        from salmonpy.runtime import RuntimeUnavailable
        for cls, old in ((InputError, ValueError), (InputTypeError, TypeError),
                         (NumericalError, FloatingPointError), (RuntimeFailure, RuntimeError),
                         (RuntimeUnavailable, RuntimeError)):
            self.assertTrue(issubclass(cls, SalmonpyError))
            self.assertTrue(issubclass(cls, old))
        with self.assertRaises(InputError):
            LorentzDrude(plasma_frequency=-1, strength=[1], resonance=[0], damping=[0])

    def test_import_is_quiet_and_optional_dependencies_stay_unloaded(self):
        code = '''import sys
import salmonpy
assert not any(x in sys.modules for x in ('mpi4py', 'matplotlib.pyplot', 'ase'))
assert salmonpy.__version__ == '0.1.0a5'
'''
        import salmonpy
        args = [sys.executable]
        if Path(salmonpy.__file__).resolve().parent != Path(__file__).resolve().parents[1]:
            args.append('-I')
        code += '\nassert salmonpy.__file__ == ' + repr(salmonpy.__file__)
        result = subprocess.run(args + ['-c', code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, '')
        self.assertEqual(result.stderr, '')

    def test_version_matches_distribution(self):
        import re
        import salmonpy
        root=Path(__file__).resolve().parents[2]
        self.assertEqual(re.search(r'^version = "([^"]+)"',
                                   (root/'pyproject.toml').read_text(),re.M)[1],salmonpy.__version__)
        if Path(salmonpy.__file__).resolve().parent != root/'salmonpy':
            from importlib.metadata import version
            self.assertEqual(version('salmonpy'),salmonpy.__version__)

    def test_numerical_failure_is_classified_and_atomic(self):
        import numpy as np
        from salmonpy.maxwell import Longitudinal
        from salmonpy.errors import NumericalError
        state=Longitudinal().new_state(dt=1e100,steps=2)
        state.step(np.zeros((1,3)))
        before=state.A.copy(); time=state.time
        with np.errstate(over='ignore',invalid='ignore'):
            with self.assertRaises(NumericalError):state.step(np.full((1,3),1e200))
        np.testing.assert_array_equal(state.A,before)
        self.assertEqual(state.time,time)
