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
        from salmonpy.api import native_contract
        legacy = native_contract('fast_hse', route='legacy')
        shared = native_contract('fast_hse', route='shared')
        self.assertEqual(legacy['model_settings_owner'], 'process compatibility settings')
        self.assertEqual(shared['model_settings_owner'], 'explicit model')
        self.assertFalse(shared['runtime_verified'])
        self.assertIn('matter J', native_contract('tddft', route='legacy')['native_current'])
        self.assertEqual(shared['native_current'], 'physical charge-current density j')
        shared['limitations'].append('host edit')
        self.assertNotIn('host edit', native_contract('fast_hse', route='shared')['limitations'])
        self.assertEqual(native_contract('lab.custom', route='anything')['status'], 'unverified')
        from unittest.mock import patch
        from salmonpy.doctor import probe_native_protocols
        with patch('subprocess.run', side_effect=AssertionError('child must not run')):
            custom = probe_native_protocols('/not-read', component='lab.custom', route='anything')
            self.assertEqual(custom['status'], 'unverified')
            self.assertFalse(custom['library_loaded'])
        import tempfile
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lib = root/'not-a-library'; lib.write_bytes(b'not native code')
            record = root/'record.json'
            record.write_text(json.dumps({'library': str(lib), 'sha256': '0'*64}))
            with patch('subprocess.run', side_effect=AssertionError('mismatched code must not load')):
                result = probe_native_protocols(record, component='tddft', route='shared')
            self.assertFalse(result['library_loaded'])
            self.assertEqual(result['errors'][0]['phase'], 'identity')
            import hashlib
            record.write_text(json.dumps({'library': str(lib), 'sha256': hashlib.sha256(lib.read_bytes()).hexdigest()}))
            child_failure = subprocess.CompletedProcess([], -11, '', 'injected child failure')
            with patch('subprocess.run', return_value=child_failure):
                failed = probe_native_protocols(record, component='tddft', route='shared')
            self.assertIsNone(failed['library_loaded'])
            self.assertEqual(failed['errors'][0]['phase'], 'child-probe')




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
from salmonpy.api import native_contract
assert native_contract('fast_hse', route='legacy')['known']
assert not any(x in sys.modules for x in ('mpi4py', 'matplotlib.pyplot', 'ase'))
assert salmonpy.__version__ == '0.1.0a5'
# Import contracts formerly repeated in six module test files. Keep per-module
# diagnostics, including ASE/f90nml laziness and the missing-MPI checkpoint path.
import importlib
for module in ('distributed', 'embedded_tddft', 'materials', 'c_material', 'ase', 'checkpointing'):
    importlib.import_module('salmonpy.' + module)
    assert 'mpi4py.MPI' not in sys.modules, module
    if module == 'ase':
        assert 'ase' not in sys.modules and 'f90nml' not in sys.modules, module
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
        # Preserve the separate unavailable-optional-dependency environment.
        missing = "import sys; sys.modules['mpi4py']=None; import salmonpy.checkpointing; assert 'mpi4py.MPI' not in sys.modules"
        result = subprocess.run(args + ['-c', missing], capture_output=True, text=True, timeout=20)
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
