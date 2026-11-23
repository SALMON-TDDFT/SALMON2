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
        from salmonpy.fdtd import SalmonExecutable
        verified=SalmonExecutable.from_runtime(self.root/'installed',required=['tddft'],host=self.host)
        self.assertEqual(verified.sha256,engine.sha256)
        self.assertEqual(SalmonExecutable.from_runtime(self.root/'installed',host=self.host).sha256,engine.sha256)
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
                runtime.build_runtime(self.root/'build', jobs=jobs, cmake_args=options)
            self.assertFalse((self.root/'build').exists())

    def test_build_failure_keeps_diagnostics(self):
        import subprocess
        with patch.object(runtime.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'fixture')):
            with self.assertRaises(subprocess.CalledProcessError):
                runtime.build_runtime(self.root/'build', jobs=1)
        record = json.loads((self.root/'build/build.json').read_text())
        self.assertEqual(record['status'], 'failed')
        configure=next(command for command in record['commands'] if command[0]=='cmake')
        self.assertIn('-DUSE_HSE=OFF',configure)
        self.assertTrue((self.root/'build/build.log').exists())

    def test_no_implicit_build(self):
        with patch('subprocess.run') as run:
            runtime.inspect_runtime({'schema':'salmonpy.runtimes/1', 'artifacts':[]}, host=self.host)
            run.assert_not_called()

    def test_job_receipt_reuses_auxiliary_hashes_and_hides_token(self):
        import os
        if os.name!='posix':self.skipTest('fast receipt needs POSIX change-time semantics')
        engine=self.install();bundle=self.root/'installed'
        aux=bundle/'lib/data.bin';aux.parent.mkdir();aux.write_bytes(b'library data'*100)
        meta=json.loads((bundle/'runtime.json').read_text());meta['files']['lib/data.bin']=runtime._hash(aux)
        (bundle/'runtime.json').write_text(json.dumps(meta))
        token='job-secret-'*6;receipt=self.root/'job-proof.json'
        validation=runtime.prepare_runtime_job(bundle,receipt,job_id='job-42',job_token=token,host=self.host)
        self.assertNotIn(token,receipt.read_text());self.assertNotIn(token,repr(validation))
        with patch.object(runtime,'_hash',wraps=runtime._hash) as hashed,patch.object(Path,'rglob',side_effect=AssertionError('fast reuse must not scan the tree')):
            verified=runtime.load_runtime(bundle,host=self.host,validation=validation)
        self.assertEqual(verified.sha256,engine.sha256)
        self.assertNotIn(aux,[Path(call.args[0]) for call in hashed.call_args_list])
        self.assertIn(engine.path,[Path(call.args[0]) for call in hashed.call_args_list])
        self.assertEqual(verified.runtime_validation_report['mode'],'job-receipt')
        self.assertNotIn(token,str(verified.runtime_validation_report))
        report=verified.runtime_validation_report;report['mode']='host edit'
        self.assertEqual(verified.runtime_validation_report['mode'],'job-receipt')
        from salmonpy.material_response import TDDFTMaterialProvider,TDDFTResponseBackend
        source=self.root/'initial';source.mkdir();(source/'rt.inp').write_text('synthetic parser-only input')
        with patch.object(runtime,'load_runtime',return_value=verified) as provider_load:
            provider=TDDFTMaterialProvider(source,'rt.inp',runtime=bundle,runtime_validation=validation)
            self.assertIs(provider_load.call_args.kwargs['validation'],validation)
            provider.close()
        self.assertIs(TDDFTResponseBackend(runtime_validation=validation)._runtime_validation,validation)

        from salmonpy.fdtd import SalmonExecutable
        self.assertEqual(SalmonExecutable.from_runtime(bundle,host=self.host,validation=validation).sha256,engine.sha256)
        with self.assertRaises(FileExistsError):runtime.prepare_runtime_job(bundle,receipt,job_id='job-42',job_token=token,host=self.host)
        with self.assertRaises(ValueError):runtime.prepare_runtime_job(bundle,bundle/'bad-proof.json',job_id='job-42',job_token=token,host=self.host)

    def test_job_receipt_falls_back_safely_and_rejects_corrupt_runtime(self):
        import os
        from dataclasses import replace
        self.install();bundle=self.root/'installed';receipt=self.root/'job-proof.json'
        good=runtime.prepare_runtime_job(bundle,receipt,job_id='job-42',job_token='s'*64,host=self.host)
        for bad in (replace(good,job_id='other-job'),replace(good,job_token='x'*64),replace(good,receipt=self.root/'missing.json')):
            with self.assertWarns(RuntimeWarning):loaded=runtime.load_runtime(bundle,host=self.host,validation=bad)
            self.assertEqual(loaded.runtime_validation_report['mode'],'full')
        signed=receipt.read_text();changed=json.loads(signed);changed['payload']['manifest']['features'].append('libxc')
        receipt.write_text(json.dumps(changed))
        with self.assertWarns(RuntimeWarning),self.assertRaises(runtime.RuntimeUnavailable):
            runtime.load_runtime(bundle,host=self.host,required=['libxc'],validation=good)
        receipt.write_text(signed)
        # Metadata-only changes do not block an otherwise valid calculation.
        old=(bundle/'bin/salmon').stat();os.utime(bundle/'bin/salmon',ns=(old.st_atime_ns,old.st_mtime_ns+1_000_000))
        with self.assertWarns(RuntimeWarning):loaded=runtime.load_runtime(bundle,host=self.host,validation=good)
        self.assertEqual(loaded.runtime_validation_report['mode'],'full')
        # Same length and restored mtime still cannot bless changed executable bytes.
        (bundle/'bin/salmon').write_bytes(b'BADBAD');os.utime(bundle/'bin/salmon',ns=(old.st_atime_ns,old.st_mtime_ns))
        with self.assertWarns(RuntimeWarning),self.assertRaisesRegex(ValueError,'checksum'):
            runtime.load_runtime(bundle,host=self.host,validation=good)
        (bundle/'bin/salmon').write_bytes(b'binary')
        aux=bundle/'support.bin';aux.write_bytes(b'original')
        meta=json.loads((bundle/'runtime.json').read_text());meta['files']['support.bin']=runtime._hash(aux)
        (bundle/'runtime.json').write_text(json.dumps(meta))
        fresh=runtime.prepare_runtime_job(bundle,self.root/'fresh-proof.json',job_id='job-42',job_token='s'*64,host=self.host)
        before=aux.stat();aux.write_bytes(b'changed!');os.utime(aux,ns=(before.st_atime_ns,before.st_mtime_ns))
        with self.assertWarns(RuntimeWarning),self.assertRaisesRegex(ValueError,'checksum'):
            runtime.load_runtime(bundle,host=self.host,validation=fresh)
        aux.write_bytes(b'original')
        (bundle/'unlisted').write_text('unexpected')
        with self.assertWarns(RuntimeWarning),self.assertRaisesRegex(ValueError,'unlisted'):
            runtime.load_runtime(bundle,host=self.host,validation=good)


if __name__ == '__main__': unittest.main()
