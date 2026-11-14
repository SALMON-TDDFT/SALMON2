"""Public replay contract using a real isolated deterministic protocol worker.

This tests Python transaction/lifecycle behavior, not SALMON physics.
"""
from dataclasses import replace
import hashlib
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
from unittest.mock import patch
import numpy as np
from salmonpy.errors import SnapshotError
from salmonpy.fdtd import SalmonExecutable
from salmonpy.tddft import TDDFTInitialState, TDDFTState

WORKER = '''#!/usr/bin/env python3
import sys
j=[0.,0.,0.]; n=0
for line in sys.stdin:
    p=line.split()
    if not p: continue
    if p[0]=='INIT':
        dt=float(p[1]); j=list(map(float,p[3:6]))
    elif p[0]=='STEP':
        n+=1; a=list(map(float,p[1:])); j=[.83*x+dt*y+.001*n for x,y in zip(j,a)]
    elif p[0]=='CLOSE': break
    else: continue
    print('SALMONPY_TDDFT_V1',n,*j,flush=True)
'''

class TDDFTSnapshots(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        root=Path(self.tmp.name); bundle=root/'bundle';bundle.mkdir()
        restart=bundle/'data_for_restart';restart.mkdir()
        for name in ('info.bin','occupation.bin','wfn.bin'): (restart/name).write_bytes(b'protocol fixture')
        (bundle/'rt.inp').write_text('&calculation\n/\n')
        worker=root/'worker';worker.write_text(WORKER);worker.chmod(0o755)
        # Windows CreateProcess cannot execute shebang scripts directly.
        popen=subprocess.Popen
        def launch(args,**kwargs):
            if args==[str(worker)]: args=[sys.executable,str(worker)]
            return popen(args,**kwargs)
        launcher=patch('salmonpy.tddft.subprocess.Popen',side_effect=launch)
        launcher.start();self.addCleanup(launcher.stop)
        exe=SalmonExecutable(worker,hashlib.sha256(worker.read_bytes()).hexdigest())
        self.model=TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=exe,max_steps=30)
        self.addCleanup(self.model.close)
    def state(self,**kw):
        s=self.model.new_state(dt=kw.pop('dt',.02),**kw);self.addCleanup(s.close);return s
    def test_replay_and_independent_branches(self):
        a,b=self.state(),self.state();initial=a.snapshot()
        values=np.random.default_rng(41).normal(size=(25,3))*.01
        for v in values[:5]:a.step(v)
        saved=a.snapshot();at_save=a.current.copy();b.restore(saved)
        self.assertEqual(saved.step_index,5);self.assertEqual(b.time,.1)
        np.testing.assert_array_equal(b.current,at_save)
        expected=[a.step(v) for v in values[5:]]
        a.restore(saved)
        for v,j in zip(values[5:],expected):
            np.testing.assert_array_equal(a.step(v),j)
            np.testing.assert_array_equal(b.step(v),j)
        a.restore(initial);self.assertEqual(a.time,0.)
        np.testing.assert_array_equal(a.current,np.zeros(3))
        # Independent owned bundles permit replay after the initial loader closes.
        self.model.close();b.restore(saved);np.testing.assert_array_equal(b.current,at_save)
    def test_incompatible_and_corrupt_tokens_preserve_worker(self):
        a=self.state();a.step([.1,.2,.3]);saved=a.snapshot()
        for target in (self.state(dt=.03),self.state(initial_A=[.1,0,0])):
            with self.assertRaises(SnapshotError): target.restore(saved)
        for bad in (None,replace(saved,_owner=object()),replace(saved,step_index=-1),
                    replace(saved,_data=()),replace(saved,time=float('nan'))):
            pid=a._proc.pid;current=a.current.copy()
            with self.assertRaises(SnapshotError):a.restore(bad)
            self.assertEqual(a._proc.pid,pid);np.testing.assert_array_equal(a.current,current)
    def test_failed_replay_preserves_live_branch(self):
        a=self.state();a.step([.1,0,0]);saved=a.snapshot();a.step([.2,0,0])
        pid=a._proc.pid;before=a.current.copy()
        real=TDDFTState._receive
        def wrong(state,index):return real(state,index)+1
        with patch.object(TDDFTState,'_receive',wrong):
            with self.assertRaises(SnapshotError):a.restore(saved)
        self.assertEqual(a._proc.pid,pid);self.assertEqual(a.time,.04)
        np.testing.assert_array_equal(a.current,before)
        a.step([.3,0,0])
    def test_spawn_failure_and_intermediate_mismatch_are_atomic(self):
        a=self.state()
        for _ in range(3): a.step([.1,0,0])
        saved=a.snapshot();before=a.current.copy();pid=a._proc.pid
        with patch.object(TDDFTState,'_spawn',side_effect=OSError('cannot start')):
            with self.assertRaises(SnapshotError):a.restore(saved)
        real=TDDFTState._receive
        def wrong_middle(state,index):
            j=real(state,index)
            return j+1 if index==2 else j
        with patch.object(TDDFTState,'_receive',wrong_middle):
            with self.assertRaises(SnapshotError):a.restore(saved)
        self.assertEqual(a._proc.pid,pid);np.testing.assert_array_equal(a.current,before)
        a.step([.2,0,0])

    def test_finish_exports_restored_branch(self):
        a=self.state();a.step([.1,0,0]);saved=a.snapshot()
        a.step([999,999,999]);a.restore(saved);expected=a.step([.2,0,0])
        archive=a.finish(Path(self.tmp.name)/'output')
        np.testing.assert_array_equal(archive.array('A'),[[0,0,0],[.1,0,0],[.2,0,0]])
        np.testing.assert_array_equal(archive.array('current')[-1],expected)
        self.assertEqual(archive.metadata['steps'],2)

    def test_nonfinite_and_layout_tampering_rejected(self):
        from salmonpy.snapshots import _integrity
        a=self.state();a.step([.1,0,0]);saved=a.snapshot();pid=a._proc.pid
        for shape,payload in (((1,3),np.full((1,3),np.nan).tobytes()),
                              ((3,1),np.zeros((1,3)).tobytes())):
            records=tuple((name,shape,dtype,order,payload) if name=='_history' else record
                          for record in saved._data for name,_,dtype,order,_ in [record])
            bad=replace(saved,_data=records,
                _integrity=_integrity(saved.kind,saved.step_index,saved.time,saved._settings,records))
            with self.assertRaises(SnapshotError):a.restore(bad)
        self.assertEqual(a._proc.pid,pid)

    def test_retirement_failure_does_not_report_failed_restore(self):
        a=self.state();a.step([.1,0,0]);saved=a.snapshot();expected=a.current.copy()
        a.step([.2,0,0]);old=a._proc;stop=TDDFTState._stop
        def fail_old(state,**kw):
            if state._proc is old: raise OSError('injected old-worker cleanup failure')
            return stop(state,**kw)
        with patch.object(TDDFTState,'_stop',fail_old):a.restore(saved)
        self.assertEqual(a.time,.02);np.testing.assert_array_equal(a.current,expected)
        self.assertIsNotNone(old.poll());self.assertTrue(a._restore_cleanup_errors)
        a.step([.2,0,0])

    def test_failed_candidate_cleanup_does_not_mask_replay_error(self):
        a=self.state();saved=a.snapshot();old=a._proc
        receive=TDDFTState._receive;stop=TDDFTState._stop
        def wrong(state,index):return receive(state,index)+1
        def fail_candidate(state,**kw):
            if state._proc is not old: raise OSError('candidate cleanup failure')
            return stop(state,**kw)
        with patch.object(TDDFTState,'_receive',wrong),patch.object(TDDFTState,'_stop',fail_candidate):
            with self.assertRaisesRegex(SnapshotError,'original state preserved'):a.restore(saved)
        self.assertIs(a._proc,old);self.assertIsNone(old.poll());a.step([.1,0,0])

    def test_tolerated_replay_roundoff_allows_repeated_restore(self):
        a=self.state();saved=a.snapshot();receive=TDDFTState._receive
        def rounded(state,index):return receive(state,index)+1e-15
        with patch.object(TDDFTState,'_receive',rounded):
            a.restore(saved);a.restore(saved)
            # A new token must also be internally consistent with replayed currents.
            a.step([.1,0,0]);second=a.snapshot();a.restore(second)
        self.assertEqual(a.time,.02)

    def test_closed_rejects_snapshot_and_restore(self):
        a=self.state();s=a.snapshot();a.close()
        with self.assertRaises(SnapshotError):a.snapshot()
        with self.assertRaises(SnapshotError):a.restore(s)
