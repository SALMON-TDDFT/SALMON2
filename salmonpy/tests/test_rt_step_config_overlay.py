"""Opt-in preparation validation and failed-snapshot preservation (no build)."""
import importlib.util
import json
import shutil
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import salmonpy

package=Path(salmonpy.__file__).resolve().parent
root=package/'_runtime_source' if (package/'_runtime_source').is_dir() else package.parent
spec=importlib.util.spec_from_file_location('rt_step_overlay',root/'tools/prepare_rt_step_config.py')
overlay=importlib.util.module_from_spec(spec);spec.loader.exec_module(overlay)
PROFILE=json.loads((root/'runtime/rt-step-config.json').read_text())


class StepConfigPreparation(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.source=Path(self.tmp.name)/'source'
        base=json.loads((root/'runtime/source.json').read_text())
        shutil.copytree(root/'runtime'/base['source_tree'],self.source)
        (self.source/'salmonpy-runtime-source.json').write_text(json.dumps(base))
        self.original=(self.source/'src/main.f90').read_bytes()
        self.output=Path(self.tmp.name)/'output'

    def test_profile_and_recursive_output_rejected_before_copy(self):
        with self.assertRaises(ValueError):overlay.prepare(self.source,self.source/'nested')
        record=self.source/'salmonpy-runtime-source.json';record.write_text('{}')
        with self.assertRaisesRegex(ValueError,'pinned base'):overlay.prepare(self.source,self.output)
        self.assertFalse(self.output.exists());self.assertEqual(record.read_text(),'{}')

    def test_existing_destination_preserved(self):
        self.output.mkdir();keep=self.output/'keep';keep.write_text('evidence')
        with self.assertRaises(FileExistsError):overlay.prepare(self.source,self.output)
        self.assertEqual(keep.read_text(),'evidence')
        self.assertEqual(list(self.output.iterdir()),[keep])

    def test_failed_copy_keeps_explicit_unapproved_record(self):
        with patch.object(overlay.store,'copy_files',side_effect=OSError('simulated copy failure')):
            with self.assertRaisesRegex(OSError,'simulated copy failure'):
                overlay.prepare(self.source,self.output)
        marker=json.loads((self.output/'salmonpy-rt-step-config.json').read_text())
        self.assertEqual(marker['status'],'failed-not-approved')
        self.assertEqual((self.source/'src/main.f90').read_bytes(),self.original)

    def test_already_modified_source_rejected(self):
        (self.source/'salmonpy-rt-step-config.json').write_text('existing evidence')
        with self.assertRaisesRegex(ValueError,'already has'):overlay.prepare(self.source,self.output)
        self.assertFalse(self.output.exists())

    def test_modified_unreplaced_source_is_rejected_before_copy(self):
        (self.source/'src/main.f90').write_text('changed')
        with self.assertRaisesRegex(ValueError,'source hash changed'):
            overlay.prepare(self.source,self.output)
        self.assertFalse(self.output.exists())

    def test_failed_preparation_cannot_be_used(self):
        p=self.source/'salmonpy-runtime-source.json';record=json.loads(p.read_text())
        record['status']='failed-not-approved';p.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError,'verified current base'):
            overlay.prepare(self.source,self.output)
        self.assertFalse(self.output.exists())
