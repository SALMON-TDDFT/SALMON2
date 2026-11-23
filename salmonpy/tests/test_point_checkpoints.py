import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from salmonpy.errors import SnapshotError
from salmonpy.point_runs import PointRun

class Model:
    def record(self): return {'material_id':'Si','definition_digest':'a'*64,'initial_state_id':'gs','initial_state_digest':'b'*64}
    def new_state(self, **kwargs): return Fake()

class Fake:
    def save_checkpoint(self, directory):
        root=Path(directory);root.mkdir()
        (root/'manifest.json').write_text('{"test":true}')
    def restore_checkpoint(self, directory, **kwargs): self.restored=Path(directory)
    def close(self): pass

class PointCheckpointTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)/'checkpoint'
        self.run=PointRun({'p':Model(),'q':Model()},run_id='run')
        self.addCleanup(self.run.close)
        self.a=self.run.new_state('p');self.b=self.run.new_state('q')
        # Exercise envelope validation without invoking a native worker.
        self.check=patch('salmonpy._point_checkpoint._native_state',side_effect=lambda state:state._state)

    def test_roundtrip_identity_and_mismatch_no_restore(self):
        with self.check:
            self.a.save_checkpoint(self.root)
            self.a.restore_checkpoint(self.root)
            self.assertEqual(self.a._state.restored,self.root/'state')
            with self.assertRaises(SnapshotError):self.b.restore_checkpoint(self.root)
            self.assertFalse(hasattr(self.b._state,'restored'))
            with PointRun({'p':Model()},run_id='other') as other:
                with self.assertRaises(SnapshotError):other.new_state('p').restore_checkpoint(self.root)

    def test_corruption_and_extra_members_rejected(self):
        with self.check:
            self.a.save_checkpoint(self.root)
            (self.root/'state/manifest.json').write_text('changed')
            with self.assertRaises(SnapshotError):self.a.restore_checkpoint(self.root)
            self.assertFalse(hasattr(self.a._state,'restored'))

    def test_failed_save_cleans_only_new_directory(self):
        with self.check:
            with patch.object(self.a._state,'save_checkpoint',side_effect=RuntimeError('failed')):
                with self.assertRaises(RuntimeError):self.a.save_checkpoint(self.root)
            self.assertFalse(self.root.exists())
            self.root.mkdir();(self.root/'keep').write_text('keep')
            with self.assertRaises(FileExistsError):self.a.save_checkpoint(self.root)
            self.assertEqual((self.root/'keep').read_text(),'keep')

    def test_bad_envelope_rejected(self):
        with self.check:
            self.a.save_checkpoint(self.root)
            manifest=self.root/'point.json'
            original=manifest.read_text()
            for text in ('[]', '{"schema":"one","schema":"two"}', original+' '*70000):
                manifest.write_text(text)
                with self.assertRaises(SnapshotError):self.a.restore_checkpoint(self.root)
            manifest.write_text(original)
            (self.root/'extra').write_text('extra')
            with self.assertRaises(SnapshotError):self.a.restore_checkpoint(self.root)

    def test_custom_backend_not_implicitly_serialized(self):
        with self.assertRaises(SnapshotError):self.a.save_checkpoint(self.root)
        self.assertFalse(self.root.exists())

    def test_invalid_budget_rejected_before_native_restore(self):
        with self.check:
            self.a.save_checkpoint(self.root)
            for budget in (0, True, -1, '10', 1):
                with self.subTest(budget=budget), self.assertRaises(SnapshotError):
                    self.a.restore_checkpoint(self.root, max_bytes=budget)
            self.assertFalse(hasattr(self.a._state,'restored'))

    def test_symlink_root_and_inner_manifest_rejected(self):
        with self.check:
            self.a.save_checkpoint(self.root)
            alias=self.root.parent/'alias';alias.symlink_to(self.root,target_is_directory=True)
            with self.assertRaises(SnapshotError):self.a.restore_checkpoint(alias)
            manifest=self.root/'state/manifest.json'
            saved=self.root.parent/'saved.json';manifest.rename(saved);manifest.symlink_to(saved)
            with self.assertRaises(SnapshotError):self.a.restore_checkpoint(self.root)
            self.assertFalse(hasattr(self.a._state,'restored'))

    def test_closed_state_rejected(self):
        self.a.close()
        with self.assertRaisesRegex(SnapshotError,'closed'):self.a.save_checkpoint(self.root)
        with self.assertRaisesRegex(SnapshotError,'closed'):self.a.restore_checkpoint(self.root)
