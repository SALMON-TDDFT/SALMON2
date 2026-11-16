"""Protocol/lifecycle tests with an executable mock, not a physics validation."""
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np

from salmonpy.errors import InputError, RuntimeFailure
from salmonpy.fast_hse import FastHSEMaterialProvider, FastHSEResponseBackend


@unittest.skipUnless(os.name == 'posix', 'worker requires POSIX')
class FastHSEProtocolTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.bundle = self.root / 'initial'
        (self.bundle / 'data_for_restart').mkdir(parents=True)
        for name in ('info.bin', 'occupation.bin', 'wfn.bin'):
            (self.bundle / 'data_for_restart' / name).write_bytes(b'fixture')
        self.text = ("&calculation theory='tddft_response'/\n&functional xc='hse06'/\n"
                     "&propagation propagator='hse_taylor4', n_hamil=4/\n&tgrid dt=.1,nt=99/\n")
        (self.bundle / 'rt.inp').write_text(self.text)
        self.executable = self.root / 'worker'
        self.executable.write_text('#!' + sys.executable + '\n' + '''
import os, sys, time
from pathlib import Path
Path('pid').write_text(str(os.getpid()))
step = 0
print('ordinary native output ignored by protocol', flush=True)
print('SALMONPY_HSE_READY 0 1 2 3', flush=True)
for line in sys.stdin:
    command, _, data = line.strip().partition(' ')
    if command == 'STEP':
        if Path('hang').exists():
            time.sleep(60)
        if Path('bad-response').exists():
            print('SALMONPY_HSE_STEP 999 4 5 6 7', flush=True)
            continue
        step += 1
        print('SALMONPY_HSE_STEP', step, '4 5 6 7', flush=True)
    elif command == 'SAVE':
        Path(data).write_text(str(step))
        print('SALMONPY_HSE_SAVED', step, flush=True)
    elif command == 'LOAD':
        step = int(Path(data).read_text())
        print('SALMONPY_HSE_LOADED', step, '4 5 6', flush=True)
    elif command == 'CLOSE':
        print('SALMONPY_HSE_CLOSED', step, flush=True)
        break
''')
        self.executable.chmod(0o700)
        self.record = self.root / 'build.json'
        self.metadata = dict(schema='salmonpy.fast-hse-worker/1', status='built-and-audited',
                             executable=str(self.executable), origin_sha256='a' * 64,
                             sha256=hashlib.sha256(self.executable.read_bytes()).hexdigest())
        self.record.write_text(json.dumps(self.metadata))
        self.provider = FastHSEMaterialProvider(self.bundle, max_steps=3)

    def state(self, *, timeout=2., context=None):
        backend = FastHSEResponseBackend(self.record, self.root, timeout=timeout)
        state = backend.create_state(self.provider.prepare(context), 7, .01, context)
        self.addCleanup(state.close)
        return state

    def test_current_endpoint_and_capacity(self):
        state = self.state()
        self.assertEqual(state.time, 0.)
        np.testing.assert_array_equal(state.initial_current, [1, 2, 3])
        old = np.zeros(3); new = np.array([.1, .2, .3])
        np.testing.assert_array_equal(state.advance(old, new), [4, 5, 6])
        self.assertEqual(state.energy, 7.)
        np.testing.assert_array_equal(old, [0, 0, 0])
        new[:] = 9
        np.testing.assert_array_equal(state.A, [.1, .2, .3])
        with self.assertRaises(InputError):
            state.advance([0, 0, 0], [1, 2, 3])
        self.assertEqual(state.steps, 1)
        self.assertEqual(state.time, .01)
        np.testing.assert_array_equal(state.current, [4, 5, 6])
        state.step([0, 0, 0]); state.step([0, 0, 0])
        with self.assertRaises(InputError):
            state.step([0, 0, 0])
        self.assertEqual(state.steps, 3)
        self.assertIn("propagator='hse_taylor4', n_hamil=4", (state.work / 'inputfile').read_text())
        self.assertIn('nt=3', (state.work / 'inputfile').read_text())

    def test_direct_checkpoint_and_cleanup(self):
        first = self.state(); first.step([.1, .2, .3])
        checkpoint = self.root / 'state with spaces.bin'
        first.checkpoint(checkpoint, False)
        with self.assertRaises(InputError):
            first.checkpoint(checkpoint, True)
        first.close(); first.close()
        self.assertIsNotNone(first._process.poll())
        self.assertFalse(first._reader.is_alive())
        fresh = self.state(); fresh.checkpoint(checkpoint, True)
        self.assertEqual(fresh.checkpoint_protocol, 'salmonpy.direct-state/1')
        self.assertEqual(fresh.steps, 1)
        self.assertEqual(fresh.time, .01)
        np.testing.assert_array_equal(fresh.current, [4, 5, 6])
        np.testing.assert_array_equal(fresh.initial_current, [1, 2, 3])
        fresh.A[:] = [.1, .2, .3]  # Authenticated caller metadata restores endpoint A.
        fresh.step([0, 0, 0]); self.assertEqual(fresh.steps, 2)

    def test_timeout_and_protocol_failure_invalidate_and_reap(self):
        for marker, timeout in [('hang', .3), ('bad-response', 2.)]:
            with self.subTest(marker=marker):
                state = self.state(timeout=timeout)
                (state.work / marker).touch()
                with self.assertRaises(RuntimeFailure):
                    state.step([0, 0, 0])
                self.assertTrue(state.closed)
                self.assertIsNotNone(state._process.poll())
                self.assertFalse(state._reader.is_alive())
                self.assertEqual(state.steps, 0)
                with self.assertRaises(InputError):
                    state.step([0, 0, 0])

    def test_hashes_bundle_and_paths(self):
        changed = dict(self.metadata, sha256='0' * 64)
        self.record.write_text(json.dumps(changed))
        with self.assertRaises(InputError):
            FastHSEResponseBackend(self.record, self.root)
        self.record.write_text(json.dumps(self.metadata))
        backend = FastHSEResponseBackend(self.record, self.root)
        self.executable.write_text(self.executable.read_text() + '\n# changed\n')
        with self.assertRaises(InputError):
            backend.create_state(self.provider.prepare(None), 0, .01, None)
        with self.assertRaises(InputError):
            FastHSEMaterialProvider(self.bundle, '../rt.inp')
        (self.bundle / 'link').symlink_to(self.bundle / 'rt.inp')
        with self.assertRaises(InputError):
            self.provider.prepare(None)
        (self.bundle / 'link').unlink()
        (self.bundle / 'data_for_restart' / 'wfn.bin').write_bytes(b'changed')
        with self.assertRaises(InputError):
            self.provider.prepare(None)

    def test_nested_work_root_preserves_bundle(self):
        from salmonpy.embedded_tddft import _manifest
        nested = self.bundle / 'work'
        nested.mkdir()
        provider = FastHSEMaterialProvider(self.bundle, max_steps=3)
        before = _manifest(self.bundle)
        directories = tuple(sorted(p.relative_to(self.bundle) for p in self.bundle.rglob('*')))
        for work_root in (self.bundle, nested):
            backend = FastHSEResponseBackend(self.record, work_root)
            with self.assertRaisesRegex(InputError, 'outside the initial bundle'):
                backend.create_state(provider.prepare(None), 0, .01, None)
            self.assertEqual(_manifest(self.bundle), before)
            self.assertEqual(tuple(sorted(p.relative_to(self.bundle) for p in self.bundle.rglob('*'))), directories)

    def test_profiles_and_rank_contract(self):
        for text in [self.text.replace("xc='hse06'", "xc='pbe0'"),
                     self.text.replace("theory='tddft_response'", "theory='dft'"),
                     self.text.replace("propagator='hse_taylor4'", "propagator='middlepoint'"),
                     self.text + "&pseudo file_pseudo(1)='../external.upf'/\n"]:
            (self.bundle / 'rt.inp').write_text(text)
            with self.assertRaises(InputError):
                FastHSEMaterialProvider(self.bundle)
        (self.bundle / 'rt.inp').write_text(self.text)
        context = SimpleNamespace(solver_comm=SimpleNamespace(size=2))
        with self.assertRaises(InputError):
            self.state(context=context)

    def test_atomic_units_and_external_A_only(self):
        for extra in ("&units unit_system='A_eV_fs'/", "&emfield ae_shape1='Acos2'/",
                      "&emfield ae_shape1='impulse',e_impulse=.01/", "&emfield trans_longi='lo'/"):
            (self.bundle/'rt.inp').write_text(self.text+extra+'\n')
            with self.assertRaises(InputError):
                FastHSEMaterialProvider(self.bundle,max_steps=3)
        (self.bundle/'rt.inp').write_text(self.text+"&emfield ae_shape1='impulse',e_impulse=0d0/\n")
        FastHSEMaterialProvider(self.bundle,max_steps=3)


if __name__ == '__main__':
    unittest.main()
