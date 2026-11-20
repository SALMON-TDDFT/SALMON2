"""Offline input collection and transferred-file integrity contracts."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from salmonpy.config import MaterialRegistry
from salmonpy.errors import InputError


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.preparation'),
                             'offline input preparation is missing')
        from salmonpy import preparation
        self.api = preparation
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root/'source.dat').write_text('synthetic input\n')
        self.registry = MaterialRegistry()
        def resolve(model, parameters, *, base_directory):
            if set(parameters) != {'input'}:
                raise InputError('expected input')
            if not (base_directory/parameters['input']).is_file():
                raise InputError('missing material input')
            return parameters
        self.registry.register('local', version='1', models=['test'], resolve_parameters=resolve)
        self.recipe = {'schema': 'salmonpy.materials/1', 'materials': {
            'layer-A': {'backend': 'local', 'model': 'test',
                        'parameters': {'input': 'assets/input.dat'}}}}
        self.assets = {'assets/input.dat': 'source.dat'}

    def prepare(self, **kwargs):
        return self.api.prepare_inputs(self.recipe, self.root/'bundle',
            registry=self.registry, base_directory=self.root,
            files=kwargs.pop('files', self.assets), **kwargs)

    def test_relocation_and_hash_anchor(self):
        dest = self.prepare(provenance={'assets/input.dat': {'source': 'synthetic', 'license': 'CC0'}})
        anchor = hashlib.sha256((dest/'manifest.json').read_bytes()).hexdigest()
        moved = self.root/'transferred'
        shutil.move(dest, moved)
        (self.root/'source.dat').unlink()
        result = self.api.verify_inputs(moved, manifest_sha256=anchor)
        self.assertEqual(result['schema'], 'salmonpy.input-bundle/1')
        self.assertEqual(result['files']['assets/input.dat']['provenance']['source'], 'synthetic')
        resolved = json.loads((moved/'resolved.json').read_text())
        self.assertIn('layer-A', resolved['materials'])
        self.assertEqual(resolved['materials']['layer-A']['parameters']['input'], 'assets/input.dat')
        with self.assertRaises(InputError): self.api.verify_inputs(moved, manifest_sha256='0'*64)

    def test_missing_files_reported_without_output(self):
        with self.assertRaisesRegex(InputError, 'missing.one.*missing.two'):
            self.prepare(files={'a': 'missing.one', 'b': 'missing.two'})
        self.assertFalse((self.root/'bundle').exists())

    def test_no_overwrite(self):
        dest = self.prepare()
        before = (dest/'manifest.json').read_bytes()
        with self.assertRaises(FileExistsError): self.prepare()
        self.assertEqual((dest/'manifest.json').read_bytes(), before)

    def test_resolver_failure_cleans_staging(self):
        self.recipe['materials']['layer-A']['parameters']['input'] = 'not-staged.dat'
        with self.assertRaises(InputError): self.prepare()
        self.assertEqual(sorted(p.name for p in self.root.iterdir()), ['source.dat'])

    def test_unsafe_and_colliding_names(self):
        for name in ('../escape', '/absolute', 'a/../b', 'a\\b', 'a//b', 'manifest.json',
                     'resolved.json', 'materials.json', 'CON', 'a:', 'a.'):
            with self.subTest(name=name), self.assertRaises(InputError):
                self.prepare(files={name: 'source.dat'})
        for files in ({'a': 'source.dat', 'A': 'source.dat'},
                      {'a': 'source.dat', 'a/b': 'source.dat'},
                      {'Dir/a': 'source.dat', 'dir/b': 'source.dat'}):
            with self.assertRaises(InputError): self.prepare(files=files)

    def test_symlinks_and_directories_are_not_assets(self):
        (self.root/'link').symlink_to(self.root/'source.dat')
        for path in ('link', '.'):
            with self.assertRaises(InputError): self.prepare(files={'asset': path})

    def test_corruption_missing_and_extra_files(self):
        dest = self.prepare()
        asset = dest/'assets/input.dat'
        original = asset.read_bytes()
        asset.write_bytes(b'corrupt')
        with self.assertRaises(InputError): self.api.verify_inputs(dest)
        asset.unlink()
        with self.assertRaises(InputError): self.api.verify_inputs(dest)
        asset.write_bytes(original)
        (dest/'unexpected').write_text('extra')
        with self.assertRaises(InputError): self.api.verify_inputs(dest)
        (dest/'unexpected').unlink()
        asset.unlink()
        asset.symlink_to(self.root/'source.dat')
        with self.assertRaises(InputError): self.api.verify_inputs(dest)

    def test_manifest_does_not_enable_path_escape(self):
        dest = self.prepare()
        p = dest/'manifest.json'
        original = json.loads(p.read_text())
        for name in ('../source.dat', '/etc/passwd'):
            data = json.loads(json.dumps(original))
            data['files'][name] = data['files'].pop('assets/input.dat')
            p.write_text(json.dumps(data))
            with self.assertRaises(InputError): self.api.verify_inputs(dest)
        p.write_text('{"schema":"a","schema":"b"}')
        with self.assertRaises(InputError): self.api.verify_inputs(dest)

    def test_unknown_provenance_and_nonfinite_are_rejected(self):
        for provenance in ({'missing': {}}, {'assets/input.dat': {'number': float('nan')}}):
            with self.assertRaises(InputError): self.prepare(provenance=provenance)

    def test_unicode_aliases_rejected_before_resolver(self):
        calls = []
        registry = MaterialRegistry()
        registry.register('local', version='1', models=['test'],
                          resolve_parameters=lambda *a, **k: calls.append(True) or {})
        try:
            self.api.prepare_inputs(self.recipe, self.root/'bundle', registry=registry,
                base_directory=self.root, files={'café/a': 'source.dat', 'cafe\u0301/b': 'source.dat'})
        except InputError:
            pass
        else:
            self.fail('canonically equivalent directory names accepted')
        self.assertEqual(calls, [], 'collision must be rejected before user code runs')

    def test_manifest_case_alias_on_case_sensitive_filesystem(self):
        dest = self.prepare()
        alias = dest/'MANIFEST.JSON'
        if alias.exists():
            self.skipTest('requires case-sensitive filesystem')
        alias.write_text('extra asset')
        p = dest/'manifest.json'
        data = json.loads(p.read_text())
        data['files']['MANIFEST.JSON'] = {
            'sha256': hashlib.sha256(alias.read_bytes()).hexdigest(),
            'size': alias.stat().st_size, 'provenance': {}}
        p.write_text(json.dumps(data))
        with self.assertRaises(InputError): self.api.verify_inputs(dest)


class JobPreflightTests(unittest.TestCase):
    setUp = PreparationTests.setUp
    prepare = PreparationTests.prepare
    def test_report_distinguishes_scope_and_failures(self):
        dest = self.prepare()
        result = self.api.preflight_inputs(dest, registry=self.registry)
        self.assertTrue(result['verified'])
        self.assertFalse(result['solver_started'])
        self.assertEqual(result['errors'], [])
        (dest/'assets/input.dat').write_text('corrupt')
        result = self.api.preflight_inputs(dest, registry=self.registry)
        self.assertFalse(result['verified'])
        self.assertEqual(result['errors'][0]['phase'], 'integrity')

    def test_changed_resolver_rejected_without_modifying_inputs(self):
        dest = self.prepare()
        before = {str(p): p.read_bytes() for p in dest.rglob('*') if p.is_file()}
        registry = MaterialRegistry()
        registry.register('local', version='2', models=['test'],
                          resolve_parameters=lambda model, parameters, **kw: parameters)
        result = self.api.preflight_inputs(dest, registry=registry)
        self.assertFalse(result['verified'])
        self.assertEqual(result['errors'][0]['phase'], 'materials')
        self.assertEqual(before, {str(p): p.read_bytes() for p in dest.rglob('*') if p.is_file()})


class StaticJobCheckTests(unittest.TestCase):
    def test_code_is_never_executed_and_errors_are_independent(self):
        from salmonpy.preparation import check_job
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            code, inputs, sentinel = root/'run.py', root/'input.json', root/'executed'
            code.write_text("from pathlib import Path\nPath("+repr(str(sentinel))+").write_text('BAD')\nraise RuntimeError('must not run')\n")
            inputs.write_text('{"steps": 3}')
            self.assertEqual(check_job(code, inputs)['status'], 'passed')
            self.assertFalse(sentinel.exists())
            inputs.write_text('{"backend":"lab.custom","parameters":{"user_knob":17}}')
            self.assertEqual(check_job(code, inputs)['status'], 'passed')
            code.write_text('return 1\n')
            inputs.write_text('{broken')
            report = check_job(code, inputs)
            self.assertEqual(report['status'], 'failed')
            self.assertEqual({e['phase'] for e in report['errors']}, {'python', 'inputs'})
            self.assertFalse(sentinel.exists())
