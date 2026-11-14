"""Disk persistence contract; deterministic workers test interface, not physics."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import test_tddft_snapshots as fixture
from salmonpy.errors import SnapshotError
from salmonpy.tddft import TDDFTInitialState, TDDFTState


class TDDFTCheckpoints(unittest.TestCase):
    setUp = fixture.TDDFTSnapshots.setUp
    state = fixture.TDDFTSnapshots.state

    def path(self, name='checkpoint'):
        return Path(self.tmp.name)/name

    def independent(self, **kwargs):
        from dataclasses import replace
        model=TDDFTInitialState.from_salmon(self.model._snapshot,'rt.inp',
            executable=self.model.executable,max_steps=30)
        self.addCleanup(model.close)
        if kwargs: model=replace(model,**kwargs)
        state=model.new_state(dt=.02);self.addCleanup(state.close)
        return state

    def test_disk_recreates_independent_model_and_future_steps(self):
        a=self.state();values=np.random.default_rng(37).normal(size=(25,3))*.01
        for v in values[:5]:a.step(v)
        a.save_checkpoint(self.path())
        expected=np.array([a.step(v) for v in values[5:]])
        b=self.independent();self.assertIsNot(a._model,b._model)
        a.close();self.model.close()
        self.assertIsNone(b.restore_checkpoint(self.path()))
        self.assertEqual(b.time,.1)
        actual=np.array([b.step(v) for v in values[5:]])
        self.assertEqual(expected.tobytes(),actual.tobytes())
        files=[str(p.relative_to(self.path())) for p in self.path().rglob('*') if p.is_file()]
        self.assertFalse(any('wfn' in p or 'restart' in p or 'pseudo' in p for p in files))

    def test_save_cannot_modify_initial_bundles(self):
        a=self.state()
        for root in (a._bundle,self.model._snapshot,self.model._source_root):
            target=root/'checkpoint'
            with self.assertRaises(SnapshotError):a.save_checkpoint(target)
            self.assertFalse(target.exists())

    def test_zero_step_archive_and_overwrite_refusal(self):
        a=self.state();a.save_checkpoint(self.path());a.step([1,0,0]);a.restore_checkpoint(self.path())
        self.assertEqual(a.time,0.)
        before=(self.path()/'manifest.json').read_bytes()
        with self.assertRaises(FileExistsError):a.save_checkpoint(self.path())
        self.assertEqual(before,(self.path()/'manifest.json').read_bytes())

    def test_corruption_extra_file_and_incomplete_archive_preserve_state(self):
        a=self.state();a.step([.1,0,0]);a.save_checkpoint(self.path());pid=a._proc.pid;current=a.current.copy()
        chunk=self.path()/'history/c/0/0';raw=chunk.read_bytes();chunk.write_bytes(b'x'+raw[1:])
        with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())
        chunk.write_bytes(raw);(self.path()/'unexpected').write_text('extra')
        with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())
        (self.path()/'unexpected').unlink();(self.path()/'manifest.json').unlink()
        with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())
        self.assertEqual(a._proc.pid,pid);np.testing.assert_array_equal(a.current,current)

    def test_incompatible_settings_and_bundle_fail_before_spawn(self):
        a=self.state();a.save_checkpoint(self.path())
        other=self.state(dt=.03)
        incompatible=[self.independent(**c) for c in ({'threads':2},{'max_steps':31})]
        with patch.object(TDDFTState,'_spawn',side_effect=AssertionError('must not spawn')):
            with self.assertRaises(SnapshotError):other.restore_checkpoint(self.path())
            for b in incompatible:
                with self.assertRaises(SnapshotError):b.restore_checkpoint(self.path())
        b=self.independent();(b._bundle/'data_for_restart/wfn.bin').write_bytes(b'changed')
        with self.assertRaises(SnapshotError):b.restore_checkpoint(self.path())

    def test_bounded_reads_symlinks_and_malformed_manifest(self):
        a=self.state();a.step([.1,0,0]);a.save_checkpoint(self.path())
        for limit in (0,True,1):
            with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path(),max_bytes=limit)
        m=self.path()/'manifest.json';raw=m.read_bytes();m.write_text('[]')
        with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())
        m.write_bytes(raw)
        link=self.path('link');link.symlink_to(self.path(),target_is_directory=True)
        with self.assertRaises(SnapshotError):a.restore_checkpoint(link)
        chunk=self.path()/'history/c/0/0';payload=chunk.read_bytes();chunk.unlink()
        external=self.path('outside');external.write_bytes(payload);chunk.symlink_to(external)
        with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())

    def test_save_failure_cleanup_and_replay_failure_are_atomic(self):
        a=self.state();a.step([.1,0,0])
        from salmonpy import _tddft_checkpoint as disk
        original=disk._write
        def fail(path,raw):
            if path.name=='manifest.json':raise OSError('injected full disk')
            return original(path,raw)
        with patch.object(disk,'_write',side_effect=fail):
            with self.assertRaises(OSError):a.save_checkpoint(self.path())
        self.assertFalse(self.path().exists())
        a.save_checkpoint(self.path());a.step([.2,0,0]);pid=a._proc.pid;before=a.current.copy()
        with patch.object(TDDFTState,'_receive',side_effect=RuntimeError('worker failed')):
            with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())
        self.assertEqual(a._proc.pid,pid);np.testing.assert_array_equal(a.current,before)
        a.step([.3,0,0])

    def test_metadata_recursion_and_inventory_limits_are_normalized(self):
        a=self.state();a.save_checkpoint(self.path())
        from salmonpy import _tddft_checkpoint as disk
        with patch.object(disk.json,'loads',side_effect=RecursionError('nested JSON')):
            with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())
        unexpected=self.path()/'unexpected-empty-directory';unexpected.mkdir()
        with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())

    def test_claimed_chunk_inventory_is_bounded_before_expansion(self):
        from salmonpy import _tddft_checkpoint as disk
        a=self.state();a.save_checkpoint(self.path());a._capacity=1000000
        path=self.path()/'manifest.json';m=json.loads(path.read_text())
        m.update(binding=disk._binding(a),step_index=1000000,time=1000000*a.dt)
        for name in ('history','current_history'):
            m['arrays'][name]=disk._layout((1000000,3),(1024,3))
        path.write_text(json.dumps(m));chunks=disk._chunks
        def guarded(shape,chunk):
            if len(shape)==2 and shape[0]>1000:raise AssertionError('must reject before expanding inventory')
            return chunks(shape,chunk)
        with patch.object(disk,'_chunks',side_effect=guarded):
            with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())

    def test_default_budget_follows_model_capacity_and_explicit_budget(self):
        from salmonpy import _tddft_checkpoint as disk
        a=self.state();budget30=disk._limits(a)[0]
        a._capacity=1000000;budget_large=disk._limits(a)[0]
        self.assertGreater(budget_large,budget30)
        self.assertGreater(budget_large,48*1000000)
        a._capacity=30;a.save_checkpoint(self.path())
        size=sum(p.stat().st_size for p in self.path().rglob('*') if p.is_file())
        a.restore_checkpoint(self.path(),max_bytes=size)
        with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path(),max_bytes=size-1)

    def test_manifest_identity_clock_and_large_claims_are_rejected(self):
        import copy
        a=self.state();a.save_checkpoint(self.path());path=self.path()/'manifest.json'
        good=json.loads(path.read_text());pid=a._proc.pid
        bad=[]
        for key,value in (('step_index',True),('step_index',2**62),('time',float('nan')),('version','future')):
            m=copy.deepcopy(good);m[key]=value;bad.append(m)
        for key,value in (('executable_sha256','0'*64),('settings_sha256','0'*64)):
            m=copy.deepcopy(good);m['binding'][key]=value;bad.append(m)
        m=copy.deepcopy(good);m['arrays']['history']['shape']=[2**62,3];bad.append(m)
        for m in bad:
            path.write_text(json.dumps(m))
            with patch.object(TDDFTState,'_spawn',side_effect=AssertionError('must not spawn')):
                with self.assertRaises(SnapshotError):a.restore_checkpoint(self.path())
        self.assertEqual(a._proc.pid,pid)

    def test_chunk_boundary_and_little_endian_zarr_layout(self):
        a=self.state()
        for _ in range(5):a.step([.1,0,0])
        from salmonpy import _tddft_checkpoint as disk
        with patch.object(disk,'CHUNK_ROWS',2):
            a.save_checkpoint(self.path());b=self.independent();b.restore_checkpoint(self.path())
            np.testing.assert_array_equal(a.current,b.current)
            meta=json.loads((self.path()/'history/zarr.json').read_text())
            self.assertEqual(meta['shape'],[5,3]);self.assertEqual(meta['chunk_grid']['configuration']['chunk_shape'],[2,3])
            self.assertEqual(meta['codecs'],[{'name':'bytes','configuration':{'endian':'little'}}])
            self.assertEqual(len((self.path()/'history/c/2/0').read_bytes()),2*3*8)

if __name__=='__main__':unittest.main()
