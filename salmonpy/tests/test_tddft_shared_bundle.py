"""Shared baseline lifetime; real isolated protocol workers, no physical GS."""
from pathlib import Path
import gc
import unittest
from unittest.mock import patch
import numpy as np
import test_tddft_snapshots as fixture
from salmonpy.errors import RuntimeFailure,SnapshotError
from salmonpy.tddft import TDDFTState


class TDDFTSharedBundle(unittest.TestCase):
    setUp=fixture.TDDFTSnapshots.setUp
    state=fixture.TDDFTSnapshots.state

    def test_one_baseline_independent_workers_and_last_release(self):
        a,b=self.state(),self.state()
        self.assertEqual(a._bundle,self.model._snapshot,'states should lease the one baseline')
        self.assertEqual(a._bundle,b._bundle)
        baseline=a._bundle/'data_for_restart/wfn.bin'
        aw=Path(a._work.name)/'run/data_for_restart/wfn.bin'
        bw=Path(b._work.name)/'run/data_for_restart/wfn.bin'
        self.assertFalse(aw.samefile(baseline));self.assertFalse(bw.samefile(baseline));self.assertFalse(aw.samefile(bw))
        raw=baseline.read_bytes();aw.write_bytes(b'worker-local output')
        self.assertEqual(baseline.read_bytes(),raw);self.assertEqual(bw.read_bytes(),raw)
        self.model.close();self.model.close();self.assertTrue(baseline.exists())
        with self.assertRaises(RuntimeFailure):self.model.new_state(dt=.02)
        a.step([.1,0,0]);token=a.snapshot();a.step([.2,0,0]);a.restore(token)
        expected=a.step([.2,0,0]);b.restore(token)
        self.assertEqual(b.step([.2,0,0]).tobytes(),expected.tobytes())
        a.close();a.close();self.assertTrue(baseline.exists())
        b.close();b.close();self.assertFalse(baseline.exists())

    def test_failed_creation_and_restore_release_only_their_leases(self):
        a=self.state();a.step([.1,0,0]);saved=a.snapshot();pid=a._proc.pid
        baseline=a._bundle
        with patch.object(TDDFTState,'_spawn',side_effect=OSError('injected startup failure')):
            with self.assertRaises(OSError):self.model.new_state(dt=.02)
            with self.assertRaises(SnapshotError):a.restore(saved)
        self.assertEqual(a._proc.pid,pid);self.assertTrue(baseline.exists())
        self.model.close();a.step([.2,0,0]);a.close()
        self.assertFalse(baseline.exists(),'failed acquisitions must not leak baseline references')

    def test_corruption_refused_before_start_and_live_worker_unchanged(self):
        a=self.state();saved=a.snapshot();pid=a._proc.pid;before=a.current.copy()
        (a._bundle/'data_for_restart/wfn.bin').write_bytes(b'corrupt shared baseline')
        with patch('salmonpy.tddft.subprocess.Popen',side_effect=AssertionError('must not launch corrupted bundle')):
            with self.assertRaises(ValueError):self.model.new_state(dt=.02)
            with self.assertRaises(SnapshotError):a.restore(saved)
        self.assertEqual(a._proc.pid,pid);np.testing.assert_array_equal(a.current,before)
        a.step([.1,0,0])

    def test_failed_output_retention_survives_final_baseline_release(self):
        import shutil
        from salmonpy.tddft import OutputSaveError
        a=self.state();a.step([.1,0,0]);baseline=a._bundle
        work=Path(a._work.name);self.addCleanup(shutil.rmtree,work,ignore_errors=True)
        (work/'run/valuable.data').write_text('worker output must survive')
        self.model.close()
        with patch('salmonpy.outputs._write_archive',side_effect=OSError('disk full')):
            with self.assertRaises(OutputSaveError) as caught:a.finish(Path(self.tmp.name)/'archive')
        self.assertEqual(caught.exception.recovery_directory,work)
        self.assertFalse(baseline.exists());self.assertTrue(work.exists())
        self.assertEqual((work/'run/valuable.data').read_text(),'worker output must survive')
        a.close();self.assertTrue(work.exists())

    def test_gc_and_released_lease_cleanup(self):
        a=self.model.new_state(dt=.02);baseline=a._bundle;proc=a._proc
        self.model.close();self.assertTrue(baseline.exists())
        del a;gc.collect()
        self.assertIsNotNone(proc.poll());self.assertFalse(baseline.exists())

    def test_closed_state_cannot_fork_a_released_lease(self):
        a=self.state();lease=a._owner;a.close()
        with self.assertRaises(RuntimeFailure):lease.fork()
        self.assertTrue(self.model._snapshot.exists())

    def test_closed_owner_without_states_is_deleted(self):
        baseline=self.model._snapshot;self.model.close()
        self.assertFalse(baseline.exists())

if __name__=='__main__':unittest.main()
