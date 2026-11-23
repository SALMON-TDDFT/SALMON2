"""CI policy checks; Python3.12+ and PyYAML, no solver or GS execution."""
import fnmatch
import importlib.util
import json
import sys
import tempfile
from unittest.mock import patch
from pathlib import Path
import tomllib
import unittest
import yaml

ROOT=Path(__file__).resolve().parents[1]


def workflow(name):
    # BaseLoader follows YAML string keys, including GitHub's unquoted `on`.
    return yaml.load((ROOT/'.github/workflows'/name).read_text(),Loader=yaml.BaseLoader)


class CIConfiguration(unittest.TestCase):
    def test_all_platforms_and_python_versions_keep_wheel_checks(self):
        wheel=workflow('ld-wheels.yml')['jobs']['wheels']
        self.assertEqual(set(wheel['strategy']['matrix']['os']),
                         {'ubuntu-22.04','windows-2022','macos-15','macos-15-intel'})
        config=tomllib.loads((ROOT/'pyproject.toml').read_text())['tool']['cibuildwheel']
        self.assertIn('--suite core',config['test-command'])
        self.assertEqual(config.get('test-extras',[]),[])
        for version in ('310','311','312','313','314'):
            identifier='cp'+version+'-manylinux_x86_64'
            self.assertTrue(any(fnmatch.fnmatchcase(identifier,p) for p in config['build'].split()))
            selected=dict(config)
            for override in config['overrides']:
                if fnmatch.fnmatchcase(identifier,override['select']):selected.update(override)
            self.assertIn('--suite '+('full' if version=='312' else 'core'),selected['test-command'])
            self.assertEqual(selected.get('test-extras',[]),['ase','plot'] if version=='312' else [])

    def test_matrix_builds_the_shared_checked_sdist(self):
        jobs=workflow('ld-wheels.yml')['jobs']
        self.assertIn('sdist',jobs['wheels']['needs'])
        download=next(s for s in jobs['wheels']['steps'] if s.get('uses','').startswith('actions/download-artifact@'))
        self.assertEqual(download['with']['name'],'ld-sdist')
        build=next(s for s in jobs['wheels']['steps'] if s.get('uses','').startswith('pypa/cibuildwheel@'))
        self.assertIn('needs.sdist.outputs.archive',build['with']['package-dir'])
        serial='\n'.join(s.get('run','') for s in jobs['sdist']['steps'])
        self.assertIn('test_wheel_native_hash.py',serial)
        self.assertIn('test_installed_ld.py',serial)
        self.assertIn('test_ci_configuration.py',serial)

    def test_physics_campaigns_remain_manual_and_worker_is_consolidated(self):
        for name in ('runtime-bundles.yml','nlcc-consistency.yml'):
            self.assertEqual(set(workflow(name)['on']),{'workflow_dispatch'})
        nickel='\n'.join(s.get('run','') for s in workflow('nlcc-consistency.yml')['jobs']['nlcc']['steps'])
        self.assertIn('test_nickel_soc.py',nickel)
        self.assertEqual(nickel.count('python -m salmonpy.runtime build'),2)
        self.assertFalse((ROOT/'.github/workflows/tddft-worker.yml').exists())
        bundled=workflow('runtime-bundles.yml')
        self.assertEqual(bundled['on']['workflow_dispatch']['inputs']['scope']['options'],['bundle','worker'])
        steps=bundled['jobs']['bundle']['steps']
        worker=next(s for s in steps if s.get('name')=='Standalone TDDFT worker regression')
        self.assertEqual(worker['if'],"inputs.scope == 'worker'")
        self.assertEqual(sum(s.get('name')=='Explicit build using the public command' for s in steps),1)
        commands='\n'.join(s.get('run','') for s in steps)
        self.assertIn('test_runtime_preparation.py',commands)
        self.assertIn('relocation-test',commands)
        self.assertIn('SALMONPY_TEST_PC=y',commands)
        logs=next(s for s in steps if s.get('name')=='Build diagnostics')['with']['path']
        self.assertIn('worker-check/**/run.log',logs)
        self.assertNotIn('worker-check/*\n',logs)

    def test_installed_runner_selects_core_and_rejects_missing_tests(self):
        spec=importlib.util.spec_from_file_location('installed_runner',ROOT/'tools/test_installed_ld.py')
        runner=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)  # No installed package or numerical import.
        self.assertTrue({'test_public_api.py','test_ld.py','test_sbe.py','test_manual.py'}<=set(runner.CORE_PATTERNS))
        names=[Path(p).stem for p in runner.CORE_PATTERNS]+['test_extra']
        with tempfile.TemporaryDirectory() as tmp,patch.dict(sys.modules):
            folder=Path(tmp)/'salmonpy/tests';folder.mkdir(parents=True)
            for name in names:
                sys.modules.pop(name,None)
                (folder/(name+'.py')).write_text('import unittest\nclass Example(unittest.TestCase):\n def test_ok(self):pass\n')
            previous=list(sys.path)
            try:
                self.assertEqual(runner.build_suite(Path(tmp),'core').countTestCases(),len(names)-1)
                self.assertEqual(runner.build_suite(Path(tmp),'full').countTestCases(),len(names))
                (folder/runner.CORE_PATTERNS[0]).unlink()
                with self.assertRaisesRegex(RuntimeError,'missing'):
                    runner.build_suite(Path(tmp),'core')
            finally:sys.path[:]=previous

    def test_source_contract_and_manual_runs_are_preserved(self):
        native=workflow('native-contract.yml')
        self.assertIn('push',native['on'])
        self.assertIn('pull_request',native['on'])
        for event in ('push','pull_request'):
            paths=native['on'][event]['paths']
            required_paths=['runtime/source.json','tools/build_embedded_tddft.py',
                            'tools/check_forbidden_eigensolvers.py','tests/native/host_abi.c']
            required_paths+=list(json.loads((ROOT/'runtime/source.json').read_text())['shared_kernels'])
            for path in required_paths:
                self.assertTrue(any(fnmatch.fnmatchcase(path,pattern) for pattern in paths),path)
        commands='\n'.join(s.get('run','') for s in native['jobs']['query-and-parser']['steps'])
        for required in ('git archive HEAD','test_runtime_sources.py','check_forbidden_eigensolvers.py',
                         'build_embedded_tddft.py','--abi-contract','--parser-host --ranks 2'):
            self.assertIn(required,commands)
        for name in ('ld-wheels.yml','native-contract.yml'):
            concurrency=workflow(name)['concurrency']
            self.assertIn('workflow_dispatch',concurrency['cancel-in-progress'])
            self.assertIn('github.run_id',concurrency['group'])
            self.assertIn('github.workflow',concurrency['group'])


if __name__=='__main__':unittest.main()
