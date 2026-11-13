import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from salmonpy import runtime


class RuntimeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.host = {'system': 'linux', 'machine': 'x86_64', 'version': '2.35', 'libc': 'glibc'}
        self.meta = {'schema': 'salmonpy.runtime/1', 'interface': runtime.INTERFACE,
                     'platform': dict(self.host), 'features': ['dft', 'tddft', 'pz'],
                     'executable': 'bin/salmon', 'files': {'bin/salmon': hashlib.sha256(b'binary').hexdigest()},
                     'provenance': {'test': 'synthetic fixture'}}
        self.archive = self.root/'runtime.zip'
        self.make_archive()

    def make_archive(self, extra=None):
        with zipfile.ZipFile(self.archive, 'w') as z:
            z.writestr('runtime.json', json.dumps(self.meta))
            z.writestr('bin/salmon', b'binary')
            if extra: z.writestr(extra, b'bad')
        self.entry = {'id': 'test-runtime', 'interface': runtime.INTERFACE,
                      'platform': self.host, 'features': self.meta['features'],
                      'url': 'https://example.invalid/runtime.zip',
                      'sha256': hashlib.sha256(self.archive.read_bytes()).hexdigest()}
        self.catalog = {'schema': 'salmonpy.runtimes/1', 'artifacts': [self.entry]}

    def install(self):
        def download(url, target): target.write_bytes(self.archive.read_bytes())
        with patch.object(runtime, '_download', side_effect=download):
            return runtime.install_runtime(self.root/'installed', catalog=self.catalog,
                                           required=['tddft'], host=self.host)

    def test_install_and_revalidate_all_files(self):
        engine = self.install()
        self.assertEqual(engine.path, (self.root/'installed/bin/salmon').resolve())
        self.assertEqual(runtime.load_runtime(self.root/'installed', host=self.host).sha256, engine.sha256)
        engine.path.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            runtime.load_runtime(self.root/'installed', host=self.host)

    def test_unsupported_never_downloads(self):
        with patch.object(runtime, '_download') as download:
            with self.assertRaises(runtime.RuntimeUnavailable):
                runtime.install_runtime(self.root/'absent', catalog=self.catalog,
                                        required=['libxc'], host=self.host)
            download.assert_not_called()
        self.assertFalse((self.root/'absent').exists())

    def test_platform_and_interface(self):
        for change in ({'machine': 'arm64'}, {'version': '2.17'}, {'libc': 'musl'}):
            self.assertEqual(runtime.inspect_runtime(self.catalog, host={**self.host, **change})['available'], [])
        self.entry['interface'] = 'future/2'
        self.assertEqual(runtime.inspect_runtime(self.catalog, host=self.host)['available'], [])

    def test_archive_checksum_and_cleanup(self):
        self.entry['sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'checksum'): self.install()
        self.assertFalse((self.root/'installed').exists())

    def test_unsafe_and_unlisted_members(self):
        for name in ('../escape', '/absolute', 'C:/evil', 'extra.txt'):
            self.make_archive(name)
            with self.assertRaises(ValueError): self.install()
            self.assertFalse((self.root/'installed').exists())

    def test_catalog_manifest_mismatch(self):
        self.meta['features'] = ['pz']
        self.make_archive()
        self.entry['features'] = ['tddft']
        with self.assertRaisesRegex(ValueError, 'catalog'): self.install()

    def test_existing_destination_untouched(self):
        self.install()
        with self.assertRaises(FileExistsError): self.install()
        self.assertEqual((self.root/'installed/bin/salmon').read_bytes(), b'binary')

    def test_ci_archive_install_without_network(self):
        self.entry['url'] = None
        with patch.object(runtime, '_download') as download:
            engine = runtime.install_runtime(self.root/'local', catalog=self.catalog,
                                             archive=self.archive, required=['tddft'], host=self.host)
            self.assertTrue(engine.path.exists())
            download.assert_not_called()
        with self.assertRaisesRegex(runtime.RuntimeUnavailable, '--archive'):
            runtime.install_runtime(self.root/'missing', catalog=self.catalog, host=self.host)

    def test_build_invalid_arguments_do_not_create_directory(self):
        for jobs, options in [(0, []), (True, []), (1, ['--bad'])]:
            with self.assertRaises(ValueError):
                runtime.build_runtime(self.root/'upstream', self.root/'build', jobs=jobs, cmake_args=options)
            self.assertFalse((self.root/'build').exists())

    def test_build_failure_keeps_diagnostics(self):
        import subprocess
        with patch.object(runtime.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'fixture')):
            with self.assertRaises(subprocess.CalledProcessError):
                runtime.build_runtime(self.root/'upstream', self.root/'build', jobs=1)
        record = json.loads((self.root/'build/build.json').read_text())
        self.assertEqual(record['status'], 'failed')
        self.assertTrue((self.root/'build/build.log').exists())

    def test_no_implicit_build(self):
        with patch('subprocess.run') as run:
            runtime.inspect_runtime({'schema':'salmonpy.runtimes/1', 'artifacts':[]}, host=self.host)
            run.assert_not_called()

if __name__ == '__main__': unittest.main()
