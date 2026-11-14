"""Pinned fast-PR preparation refuses overwrites and keeps failed evidence."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import salmonpy

package=Path(salmonpy.__file__).resolve().parent
root=package/'_runtime_source' if (package/'_runtime_source').is_dir() else package.parent
spec=importlib.util.spec_from_file_location('hse_rt_prepare',root/'tools/prepare_hse_rt_state.py')
prepare=importlib.util.module_from_spec(spec);spec.loader.exec_module(prepare)
PROFILE=json.loads((root/'runtime/hse-rt-state.json').read_text())


class HSEStatePreparation(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.source=Path(self.tmp.name)/'source';self.source.mkdir()
        self.output=Path(self.tmp.name)/'output'

    def test_recursive_output_rejected_before_git(self):
        with patch.object(prepare.subprocess,'check_output',side_effect=AssertionError('must not query')):
            with self.assertRaisesRegex(ValueError,'external'):
                prepare.prepare(self.source,self.source/'nested')
        self.assertFalse((self.source/'nested').exists())

    def test_wrong_pin_rejected_before_creating_output(self):
        with patch.object(prepare.subprocess,'check_output',return_value='0'*40):
            with self.assertRaisesRegex(ValueError,'pinned PR'):
                prepare.prepare(self.source,self.output)
        self.assertFalse(self.output.exists())

    def test_existing_evidence_is_preserved(self):
        self.output.mkdir();keep=self.output/'keep';keep.write_text('existing evidence')
        with patch.object(prepare.subprocess,'check_output',return_value=PROFILE['upstream_commit']):
            with self.assertRaises(FileExistsError):prepare.prepare(self.source,self.output)
        self.assertEqual(keep.read_text(),'existing evidence')
        self.assertEqual(list(self.output.iterdir()),[keep])

    def test_failed_archive_preserved_and_not_approved(self):
        with patch.object(prepare.subprocess,'check_output',return_value=PROFILE['upstream_commit']), \
             patch.object(prepare.subprocess,'run',side_effect=OSError('archive failure')):
            with self.assertRaisesRegex(OSError,'archive failure'):prepare.prepare(self.source,self.output)
        marker=json.loads((self.output/'salmonpy-hse-rt-state.json').read_text())
        self.assertEqual(marker['status'],'failed-not-approved')
        self.assertEqual(list(self.source.iterdir()),[])
