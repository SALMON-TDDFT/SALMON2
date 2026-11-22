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
        self.assertEqual(salmonpy.__version__, '0.1.0a7')
        for module, names in PUBLIC_API.items():
            if module == 'salmonpy.ase_calculator' and importlib.util.find_spec('ase') is None:
                continue
            imported = importlib.import_module(module)
            for name in names:
                self.assertTrue(hasattr(imported, name), (module, name))
        self.assertIn('__version__', salmonpy.__all__)
        from salmonpy.api import native_contract
        single = native_contract('fast_hse', route='single')
        shared = native_contract('fast_hse', route='shared')
        self.assertEqual(single['model_settings_owner'], 'explicit single-state settings')
        self.assertEqual(shared['model_settings_owner'], 'explicit model')
        self.assertFalse(shared['runtime_verified'])
        self.assertIn('matter J', native_contract('tddft', route='single')['native_current'])
        self.assertEqual(shared['native_current'], 'physical charge-current density j')
        shared['limitations'].append('host edit')
        self.assertNotIn('host edit', native_contract('fast_hse', route='shared')['limitations'])
        self.assertEqual(native_contract('lab.custom', route='anything')['status'], 'unverified')
        from salmonpy.api import query_host_abi
        from salmonpy.errors import RuntimeFailure
        from types import SimpleNamespace
        class Probe:
            values=(1,32,64,64,32,64,1,1,1,3)
            status=0
            def __call__(self,capacity,buffer):
                self.capacity=capacity
                for i,v in enumerate(self.values):buffer[i]=v
                return self.status
        probe=Probe();report=query_host_abi(SimpleNamespace(sp_host_abi=probe))
        self.assertEqual(probe.capacity,10)
        self.assertEqual(report['version'],1)
        self.assertFalse(report['initializes_MPI'])
        self.assertEqual(report['directory_protocol_required'],3)
        for index,bad in ((0,2),(1,64),(2,32),(3,32),(4,64),(5,32),(6,2),(7,0),(8,2),(9,2)):
            probe.values=tuple(bad if i==index else v for i,v in enumerate(Probe.values))
            with self.subTest(word=index),self.assertRaises(RuntimeFailure):query_host_abi(SimpleNamespace(sp_host_abi=probe))
        probe.values=Probe.values;probe.status=13
        with self.assertRaises(RuntimeFailure):query_host_abi(SimpleNamespace(sp_host_abi=probe))
        with self.assertRaisesRegex(RuntimeFailure,'host ABI'):query_host_abi(SimpleNamespace())
        self.assertFalse(native_contract('fast_hse',route='legacy')['known'])
        from salmonpy.api import _bind_field_operations
        field_ops=native_contract('fdtd',route='single')['field_operations']
        self.assertEqual(field_ops['sp_eh_plane_read']['collective'],'local')
        self.assertEqual(field_ops['sp_eh_current_mask']['collective'],'local')
        self.assertEqual(field_ops['sp_eh_sample_h']['collective'],'solver group')
        self.assertEqual(field_ops['sp_eh_checkpoint']['modes'][2],'stage')
        class FieldFunction:
            def __call__(self,*args):raise AssertionError('binding must not advance the field')
        fns={name:FieldFunction() for name in field_ops}
        for missing in fns:
            with self.subTest(field_operation=missing),self.assertRaisesRegex(RuntimeFailure,missing):
                _bind_field_operations(SimpleNamespace(**{k:v for k,v in fns.items() if k!=missing}),tuple(fns))
            self.assertFalse(any(hasattr(f,'argtypes') for f in fns.values()))
        _bind_field_operations(SimpleNamespace(**fns),tuple(fns))
        self.assertEqual(len(fns['sp_eh_step_owned'].argtypes),5)
        self.assertEqual(len(fns['sp_eh_plane_read'].argtypes),5)

        from salmonpy.api import _bind_host_operations
        operations=native_contract('fdtd',route='single')['host_operations']
        self.assertEqual(set(operations),{'sp_attach','sp_detach','sp_info','sp_sum','sp_barrier','sp_library'})
        self.assertEqual(operations['sp_sum']['collective'],'solver group')
        self.assertEqual(operations['sp_info']['collective'],'local')
        self.assertEqual(operations['sp_detach']['requires'],'attached; all owned states closed')
        operations['sp_sum']['arguments'].append('host edit')
        self.assertNotIn('host edit',native_contract('fdtd',route='single')['host_operations']['sp_sum']['arguments'])
        class Function:
            def __call__(self,*args):raise AssertionError('binding must not execute operations')
        functions={name:Function() for name in operations}
        for missing in functions:
            partial=SimpleNamespace(**{k:v for k,v in functions.items() if k!=missing})
            with self.subTest(missing=missing),self.assertRaisesRegex(RuntimeFailure,missing):
                _bind_host_operations(partial)
            self.assertFalse(any(hasattr(f,'argtypes') for f in functions.values()))
        complete=SimpleNamespace(**functions);_bind_host_operations(complete)
        import ctypes as c
        self.assertEqual(complete.sp_attach.argtypes,[c.c_int64])
        self.assertEqual(complete.sp_info.argtypes,[c.POINTER(c.c_int),c.POINTER(c.c_int)])
        self.assertIsNone(complete.sp_library.restype)
        self.assertEqual(complete.sp_sum.restype,c.c_double)

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
assert native_contract('fast_hse', route='single')['known']
assert not any(x in sys.modules for x in ('mpi4py', 'matplotlib.pyplot', 'ase'))
assert salmonpy.__version__ == '0.1.0a7'
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
