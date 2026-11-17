from pathlib import Path
import tempfile
import unittest
from salmonpy.native_steps import _AppendJournal
from salmonpy.errors import InputError,RuntimeFailure


class NativeStepTests(unittest.TestCase):
    def test_append_rollback_preserves_inode_and_bytes(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'rt.data';p.write_bytes(b'accepted\n');inode=p.stat().st_ino
            saved=_AppendJournal.capture(root)
            with p.open('ab') as f:f.write(b'failed\n')
            saved.restore()
            self.assertEqual(p.read_bytes(),b'accepted\n');self.assertEqual(p.stat().st_ino,inode)

    def test_replaced_writer_is_rejected_before_truncation(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'rt.data';p.write_bytes(b'accepted\n');saved=_AppendJournal.capture(root)
            p.rename(root/'old');p.write_bytes(b'replaced')
            with self.assertRaises(RuntimeFailure):saved.restore()
            self.assertEqual(p.read_bytes(),b'replaced')

    def test_new_or_missing_stream_is_not_silently_deleted(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'rt.data';p.write_bytes(b'accepted\n');saved=_AppendJournal.capture(root)
            (root/'new.data').write_bytes(b'unexpected')
            with self.assertRaises(RuntimeFailure):saved.restore()
            self.assertTrue((root/'new.data').exists())

    def test_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'real').write_bytes(b'original')
            try:(root/'rt.data').symlink_to(root/'real')
            except (OSError,NotImplementedError):self.skipTest('symlinks unavailable on this platform')
            with self.assertRaises(InputError):_AppendJournal.capture(root)

    def backup(self,root,code=0):
        from types import SimpleNamespace
        import numpy as np
        from salmonpy.native_steps import _NativeBackup
        calls=[]
        def run(mode):calls.append(mode);return code if mode==0 else 0
        lib=SimpleNamespace(tx_protocol=lambda:1,tx=run)
        owner=SimpleNamespace(lib=lib,native=SimpleNamespace(_live=lambda:None),work=root,
            closed=False,failed=False,A=np.zeros(3),steps=0,_call=lambda f,*args:f(*args))
        comm=SimpleNamespace(rank=0,allgather=lambda v:[v])
        return _NativeBackup(owner,comm,'tx',root,('A','steps','failed')),owner,calls

    def test_native_python_aliases_and_output_positions_restore(self):
        import numpy as np
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'rt.data';p.write_bytes(b'accepted\n')
            backup,owner,calls=self.backup(root);A=owner.A
            saved=backup.capture();owner.A=np.ones(3);owner.steps=1
            with p.open('ab') as f:f.write(b'failed\n')
            backup.rollback(saved);backup.release(saved)
            self.assertIs(owner.A,A);np.testing.assert_array_equal(owner.A,0)
            self.assertEqual(owner.steps,0);self.assertFalse(owner.failed)
            self.assertEqual(p.read_bytes(),b'accepted\n')
            self.assertEqual(calls,[0,1,3,2]);self.assertIsNone(owner._step_guard)

    def test_partial_capture_retires_native_backup_with_none(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'rt.data';p.write_bytes(b'old')
            backup,owner,calls=self.backup(root)
            from unittest.mock import patch
            with patch.object(_AppendJournal,'capture',side_effect=OSError('journal unavailable')):
                with self.assertRaises(OSError):backup.capture()
            backup.release(None);backup.release(None)
            self.assertEqual(calls,[0,2]);self.assertIsNone(owner._step_guard)

    def test_rejected_begin_does_not_retire_other_native_reservation(self):
        with tempfile.TemporaryDirectory() as d:
            backup,owner,calls=self.backup(Path(d),14)
            with self.assertRaises(InputError):backup.capture()
            backup.release(None)
            self.assertEqual(calls,[0]);self.assertFalse(owner.failed)

    def test_old_protocol_is_not_inferred(self):
        from salmonpy.native_steps import _NativeBackup
        with tempfile.TemporaryDirectory() as d:
            _,owner,calls=self.backup(Path(d));del owner.lib.tx_protocol
            with self.assertRaises(InputError):_NativeBackup(owner,None,'tx',d,())
            self.assertEqual(calls,[])

    def test_standard_and_adapter_output_helpers_match(self):
        import salmonpy
        package=Path(salmonpy.__file__).resolve().parent
        root=package/'_runtime_source' if (package/'_runtime_source').is_dir() else package.parent
        self.assertEqual((package/'native_adapters/step_output.f90').read_bytes(),
                         (root/'runtime/sources/standard/src/io/step_output.f90').read_bytes())
