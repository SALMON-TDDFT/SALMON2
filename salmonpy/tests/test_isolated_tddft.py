import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from salmonpy.errors import InputError, SnapshotError, RuntimeFailure

WORKER = '''import sys
v=[0.]*6;n=0
for line in sys.stdin:
 p=line.split()
 if not p:continue
 if p[0]=='INIT_E':dt=float(p[1])
 elif p[0]=='STEP_E':
  n+=1;e=list(map(float,p[1:]));v=[.9*x+dt*e[k%3]+.001*n for k,x in enumerate(v)]
 elif p[0]=='CLOSE':break
 else:continue
 print('SALMONPY_TDDFT_E_V1',n,*v,flush=True)
'''


def receipt(binary):
    return dict(schema='salmonpy.isolated-worker-validation/1',
                executable_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                protocol='SALMONPY_TDDFT_E_V1',steps=20,interleaved_independence=True,
                replay_continuation_max_abs_error=0.,repeat_max_abs_error=0.,reference_max_abs_error=0.,
                capabilities=['isolated','pz','fixed-ions','middlepoint','single-rank','E-to-dipole'])


class IsolatedTDDFTTests(unittest.TestCase):
    def setUp(self):
        from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        root=Path(self.tmp.name);self.binary=root/'worker.py';self.binary.write_text(WORKER)
        self.validation=root/'validation.json';self.validation.write_text(json.dumps(receipt(self.binary)))
        self.bundle=root/'bundle';self.bundle.mkdir();restart=self.bundle/'data_for_restart';restart.mkdir()
        for n in ('info.bin','occupation.bin','wfn.bin'):(restart/n).write_bytes(b'protocol-only fixture')
        (self.bundle/'rt.inp').write_text('&calculation\n/\n')
        original=subprocess.Popen
        def launch(args,**kw):
            if len(args)==1 and Path(args[0]).resolve()==self.binary.resolve():args=[sys.executable,str(self.binary)]
            return original(args,**kw)
        mock=patch('salmonpy.tddft.subprocess.Popen',side_effect=launch);mock.start();self.addCleanup(mock.stop)
        self.runtime=load_isolated_runtime(self.binary,self.validation)
        self.model=IsolatedTDDFTInitialState.from_salmon(self.bundle,'rt.inp',runtime=self.runtime,max_steps=25)
        self.addCleanup(self.model.close)

    def state(self,dt=.02):
        state=self.model.new_state(dt=dt);self.addCleanup(state.close);return state

    def test_independent_restore_and_loader_lifetime(self):
        a,b=self.state(),self.state()
        a.step([0,0,.1]);saved=a.snapshot()
        np.testing.assert_array_equal(b.observables.dipole,[0,0,0])
        expected=a.step([0,0,.2]).dipole
        a.restore(saved);b.restore(saved)
        self.model.close()
        np.testing.assert_array_equal(a.step([0,0,.2]).dipole,expected)
        np.testing.assert_array_equal(b.step([0,0,.2]).dipole,expected)
        self.assertFalse(hasattr(a,'current'))
        self.assertFalse(hasattr(a,'save_checkpoint'))

    def test_bad_packet_snapshot_and_failed_step(self):
        from salmonpy.field_input import FieldPacket,FieldValue
        a=self.state();bad=self.state(.03).snapshot()
        with self.assertRaises(SnapshotError):a.restore(bad)
        with self.assertRaises(InputError):a.step_fields(FieldPacket())
        with self.assertRaises(InputError):a.step_fields(FieldPacket(E=FieldValue([0,0,1],time=.01,position=(0,0,0),unit='E:au')))
        self.assertEqual(a.time,0.)
        a.step([0,0,.1]);old=a.observables.dipole.copy()
        a._engine._proc.kill();a._engine._proc.wait()
        with self.assertRaises(RuntimeFailure):a.step([0,0,.2])
        np.testing.assert_array_equal(a.observables.dipole,old)
        self.assertEqual(a.time,.02)
        a.step([0,0,.2])

    def test_runtime_capability_and_drift(self):
        from salmonpy.isolated_tddft import load_isolated_runtime
        d=receipt(self.binary);d['protocol']='SALMONPY_TDDFT_V1';self.validation.write_text(json.dumps(d))
        with self.assertRaises(InputError):load_isolated_runtime(self.binary,self.validation)
        self.binary.write_text(WORKER+'\n# changed')
        with self.assertRaises(InputError):self.state()

    def test_failed_restore_preserves_original_worker(self):
        from salmonpy.isolated_tddft import _EEngine
        a=self.state();a.step([0,0,.1]);saved=a.snapshot();a.step([0,0,.2])
        before=a.observables.dipole.copy();proc=a._engine._proc
        receive=_EEngine._receive
        def fail_replay(engine,index):
            if engine is not a._engine and index>0:raise RuntimeFailure('injected replay failure')
            return receive(engine,index)
        with patch.object(_EEngine,'_receive',fail_replay):
            with self.assertRaises(SnapshotError):a.restore(saved)
        self.assertIs(a._engine._proc,proc)
        self.assertIsNone(proc.poll())
        np.testing.assert_array_equal(a.observables.dipole,before)
        self.assertEqual(a.time,.04)
        a.step([0,0,.3])

    def test_receipt_rejects_false_evidence_and_duplicate_keys(self):
        from salmonpy.isolated_tddft import load_isolated_runtime
        for key,value in [('steps',True),('interleaved_independence',False),
                          ('reference_max_abs_error',float('nan')),('capabilities',[])]:
            data=receipt(self.binary);data[key]=value;self.validation.write_text(json.dumps(data))
            with self.subTest(key=key),self.assertRaises(InputError):
                load_isolated_runtime(self.binary,self.validation)
        self.validation.write_text('{"schema":"a","schema":"b"}')
        with self.assertRaises(InputError):load_isolated_runtime(self.binary,self.validation)

    def test_position_and_protocol_snapshot_mismatch_before_spawn(self):
        from dataclasses import replace
        a=self.state()
        with self.model.new_state(dt=.02,position=(1.,0,0)) as b:
            token=b.snapshot()
        proc=a._engine._proc
        with patch('salmonpy.tddft.subprocess.Popen',side_effect=AssertionError('must not spawn')):
            with self.assertRaises(SnapshotError):a.restore(token)
            with self.assertRaises(SnapshotError):a.restore(replace(a.snapshot(),kind='TDDFT-process-replay'))
        self.assertIs(a._engine._proc,proc)
        self.assertEqual(a.time,0.)
