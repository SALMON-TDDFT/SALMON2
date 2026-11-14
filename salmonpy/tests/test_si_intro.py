"""Starter rejects changed inputs and misleading numerical comparisons."""
from contextlib import nullcontext
from types import SimpleNamespace
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import numpy as np
from salmonpy.examples import si_intro


class SiIntroTests(unittest.TestCase):
    def reference(self):
        return {'energy_hartree': -1., 'energy_atol': 1e-8,
                'time': [0., .02, .04], 'current': [[0., 0., 0.], [0., 0., .001], [0., 0., .002]],
                'current_atol': 1e-10, 'current_rtol': 1e-6}

    def test_reference_detects_sign_and_time_errors(self):
        ref = self.reference()
        self.assertTrue(si_intro.compare_result(ref, -1., ref['time'], ref['current'])['passed'])
        self.assertFalse(si_intro.compare_result(ref, -1., ref['time'], -np.array(ref['current']))['passed'])
        self.assertFalse(si_intro.compare_result(ref, -1., np.array(ref['time'])+.02, ref['current'])['passed'])
        self.assertFalse(si_intro.compare_result(ref, float('nan'), ref['time'], ref['current'])['passed'])

    def test_invalid_reference_tolerances_are_rejected(self):
        for key, value in [('current_atol', float('inf')), ('energy_atol', -1.), ('current_rtol', float('nan'))]:
            ref = self.reference(); ref[key] = value
            with self.assertRaises(ValueError):
                si_intro.compare_result(ref, -1., ref['time'], ref['current'])

    def test_unknown_artifact_is_rejected_before_install(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root/'outer.zip'; archive.write_bytes(b'wrong artifact')
            with patch.object(si_intro, '_assets', return_value={'provenance.json': json.dumps({'runtime_artifacts': []}).encode()}), \
                 patch.object(si_intro, 'install_runtime') as install:
                with self.assertRaisesRegex(ValueError, 'artifact'):
                    si_intro.prepare(root/'prepared', archive)
                install.assert_not_called()
                self.assertFalse((root/'prepared').exists())

    def test_artifact_members_are_checked_before_extraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); archive = root/'outer.zip'
            with zipfile.ZipFile(archive, 'w') as z:
                z.writestr('catalog.json', '{}')
                z.writestr('../escape', 'bad')
            provenance = {'runtime_artifacts': [{'outer_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
                                                'archive_member': 'runtime.zip'}]}
            with patch.object(si_intro, '_assets', return_value={'provenance.json': json.dumps(provenance).encode()}), \
                 patch.object(si_intro, 'install_runtime') as install:
                with self.assertRaisesRegex(ValueError, 'members'):
                    si_intro.prepare(root/'prepared', archive)
                install.assert_not_called()
                self.assertFalse((root/'escape').exists())

    def test_changed_prepared_inputs_do_not_launch_calculation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            assets = {'input.json': b'{}', 'expected.json': b'{}', 'provenance.json': b'{}', 'Si_rps.dat': b'pseudo'}
            for name, content in assets.items(): (root/name).write_bytes(content)
            (root/'input.json').write_text('{"changed":true}')
            with patch.object(si_intro, '_assets', return_value=assets), patch.object(si_intro, 'load_runtime') as load:
                with self.assertRaisesRegex(ValueError, 'input.json'):
                    si_intro.run(root, threads=1, memory_budget_gib=4)
                load.assert_not_called()
                self.assertFalse((root/'result').exists())

    def prepared_fixture(self, root):
        assets = si_intro._assets()
        for name, data in assets.items(): (root/name).write_bytes(data)
        configuration = json.loads(assets['input.json'])
        model, scf = si_intro._model(configuration, root)
        (root/'gs.inp').write_text(model.input_text(scf))
        entry = json.loads(assets['provenance.json'])['runtime_artifacts'][0]
        (root/'preparation.json').write_text(json.dumps({
            'schema': 'salmonpy.si-intro-preparation/1', 'status': 'complete',
            'asset_sha256': {name: hashlib.sha256(data).hexdigest() for name, data in assets.items()},
            'runtime_artifact': entry, 'executable_sha256': 'fake-engine'}))
        return json.loads(assets['expected.json'])

    def test_native_failure_is_not_measurement_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); self.prepared_fixture(root)
            with patch.object(si_intro, 'load_runtime', return_value=SimpleNamespace(sha256='fake-engine')), \
                 patch.object(si_intro.DFT, 'run', side_effect=RuntimeError('native failed')), \
                 patch.object(si_intro.time, 'perf_counter', side_effect=OSError('no clock')):
                with self.assertRaisesRegex(RuntimeError, 'native failed'):
                    si_intro.run(root, threads=1, memory_budget_gib=4)
            report = json.loads((root/'result/report.json').read_text())
            self.assertEqual(report['calculation']['status'], 'failed')
            self.assertEqual(report['validation']['status'], 'not_performed')
            self.assertEqual(report['resources']['wall_time']['status'], 'unavailable')

    def test_validation_failure_retains_successful_calculation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); expected = self.prepared_fixture(root)
            gs = SimpleNamespace(energy_hartree=expected['energy_hartree'], iterations=54,
                tddft_initial_state=lambda **kwargs: nullcontext(object()))
            response = SimpleNamespace(time=np.array(expected['time']),
                current=-np.array(expected['current']), field_convention='physical charge current j')
            with patch.object(si_intro, 'load_runtime', return_value=SimpleNamespace(sha256='fake-engine')), \
                 patch.object(si_intro.DFT, 'run', return_value=gs), \
                 patch.object(si_intro, 'respond', return_value=response):
                report = si_intro.run(root, threads=1, memory_budget_gib=4)
            self.assertEqual(report['status'], 'validation_failed')
            self.assertEqual(report['calculation']['status'], 'passed')
            self.assertEqual(report['validation']['status'], 'failed')
            self.assertTrue((root/'result/current.npz').is_file())
            with self.assertRaises(FileExistsError):
                with patch.object(si_intro, 'load_runtime', return_value=SimpleNamespace(sha256='fake-engine')):
                    si_intro.run(root, threads=1, memory_budget_gib=4)

    def test_broken_clock_does_not_raise(self):
        with patch.object(si_intro.time, 'perf_counter', side_effect=OSError('clock unavailable')):
            self.assertIsNone(si_intro._clock())
        self.assertEqual(si_intro._elapsed(None), {'status': 'unavailable', 'seconds': None})
