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
v=[0.]*6;n=0;updates=0;saved=None
codes={'CAPTURE_E':1,'RESTORE_E':2,'RELEASE_E':3,'RESUME_E':4,'PROBE_E':5}
for line in sys.stdin:
 p=line.split()
 if not p:continue
 if p[0] in codes:
  op=codes[p[0]];status=0
  if op==1:
   if saved is not None:status=14
   else:saved=(n,list(v))
  elif op==2:
   if saved is None:status=14
   else:n=saved[0];v=list(saved[1])
  elif op==3:
   if saved is None:status=14
   else:saved=None
  elif op==4 and saved is None:status=14
  print('SALMONPY_E_STATE_V1',op,status,n,updates,flush=True)
  continue
 if p[0]=='INIT_E':dt=float(p[1])
 elif p[0]=='STEP_E':
  n+=1;updates+=1;e=list(map(float,p[1:]));v=[.9*x+dt*e[k%3]+.001*n for k,x in enumerate(v)]
 elif p[0]=='CLOSE':break
 else:continue
 print('SALMONPY_TDDFT_E_V3',n,*v,10.-n*1.e-8,-1.+n*.01,1. if n!=2 else 0.,flush=True)
'''


def receipt(binary):
    return dict(schema='salmonpy.isolated-worker-validation/4',
                executable_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                protocol='SALMONPY_TDDFT_E_V3',steps=20,interleaved_independence=True,
                savepoint_same_process=True,savepoint_extra_evolutions=0,savepoint_stale_rejected=True,
                savepoint_continuation_max_abs_error=0.,
                replay_continuation_max_abs_error=0.,repeat_max_abs_error=0.,reference_max_abs_error=0.,reference_energy_max_abs_error=0.,
                capabilities=['isolated','pz','fixed-ions','middlepoint','single-rank','E-to-dipole','measured-electron-number','internal-energy-freshness','live-savepoint'])


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
        self.assertTrue(callable(a.save_checkpoint))

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
                          ('reference_max_abs_error',float('nan')),('capabilities',[]),
                          ('schema','salmonpy.isolated-worker-validation/3'),
                          ('savepoint_extra_evolutions',1),('savepoint_same_process',False),
                          ('savepoint_stale_rejected',False)]:
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

    def test_disk_checkpoint_rebuilt_model_and_rejections(self):
        from salmonpy.isolated_tddft import IsolatedTDDFTInitialState
        checkpoint=Path(self.tmp.name)/'checkpoint'
        with self.model.new_state(dt=.02) as a:
            for i in range(5): a.step([0.,0.,.001*i])
            a.save_checkpoint(checkpoint)
            expected=a.step([0.,0.,.007]).dipole
        manifest=json.loads((checkpoint/'manifest.json').read_text())
        self.assertEqual(manifest['schema'],'salmonpy.isolated-tddft-checkpoint/3')
        self.assertEqual(set(manifest['arrays']),{'E0','observables0','observables','E_history','observable_history'})
        with IsolatedTDDFTInitialState.from_salmon(self.bundle,'rt.inp',runtime=self.runtime,max_steps=25) as rebuilt:
            with rebuilt.new_state(dt=.02) as b:
                b.restore_checkpoint(checkpoint)
                np.testing.assert_array_equal(b.step([0.,0.,.007]).dipole,expected)
                before=b.observables.dipole.copy();clock=b.time
                with patch.object(b._engine,'_spawn',side_effect=AssertionError('must reject before spawn')):
                    with self.assertRaises(SnapshotError): b.restore_checkpoint(checkpoint,max_bytes=1)
                    chunk=checkpoint/'observables/c/0';raw=chunk.read_bytes();chunk.write_bytes(b'x'+raw[1:])
                    with self.assertRaises(SnapshotError): b.restore_checkpoint(checkpoint)
                    chunk.write_bytes(raw)
                np.testing.assert_array_equal(b.observables.dipole,before)
                self.assertEqual(b.time,clock)
            with rebuilt.new_state(dt=.02,position=(1.,0.,0.)) as c:
                with patch.object(c._engine,'_spawn',side_effect=AssertionError('must reject before spawn')):
                    with self.assertRaises(SnapshotError): c.restore_checkpoint(checkpoint)

    def test_disk_empty_and_chunk_boundary(self):
        from salmonpy.isolated_tddft import IsolatedTDDFTInitialState
        with IsolatedTDDFTInitialState.from_salmon(self.bundle,'rt.inp',runtime=self.runtime,max_steps=1025) as model:
            with model.new_state(dt=.02) as state:
                empty=Path(self.tmp.name)/'empty'
                state.save_checkpoint(empty)
                state.step([0.,0.,.001])
                state.restore_checkpoint(empty)
                self.assertEqual(state.time,0.)
                for _ in range(1025): state.step([0.,0.,.001])
                expected=state.observables.dipole.copy()
                full=Path(self.tmp.name)/'full'
                state.save_checkpoint(full)
                state.restore_checkpoint(empty)
                state.restore_checkpoint(full)
                np.testing.assert_array_equal(state.observables.dipole,expected)
                self.assertEqual(state.time,1025*.02)

    def test_fdtd_trial_commit_cancel_and_independence(self):
        from salmonpy.dipole_source import DipoleMap, _IsolatedDipoleTrial
        from salmonpy.errors import NumericalError
        a,b=self.state(),self.state()
        mapping=DipoleMap(('molecule',),(10,),[[1.]],[8.])
        bridge=_IsolatedDipoleTrial(a,mapping,sampling='endpoint')
        drive=np.array([.03,-.01,.02]);own=np.array([.1,.2,.3])
        before=a.observables.dipole.copy()
        current=bridge.propose(drive+own,own)
        from salmonpy.dipole_source import _selected_molecular_field
        expected=b.step(_selected_molecular_field(drive+own,own)).dipole
        np.testing.assert_array_equal(a.observables.dipole,expected)
        np.testing.assert_allclose(current,mapping.current([before],[expected],dt=.02),rtol=0,atol=1e-18)
        with self.assertRaises(InputError):bridge.propose(drive,own)
        with self.assertRaises(NumericalError):bridge.validate(drive+own+[.01,0,0],own)
        bridge.cancel()
        self.assertEqual(a.time,0.)
        np.testing.assert_array_equal(a.observables.dipole,before)
        bridge.propose(drive+own,own)
        bridge.validate(drive+own,own)
        bridge.commit()
        np.testing.assert_array_equal(a.observables.dipole,expected)
        bridge.propose(drive+own,own);bridge.cancel()
        np.testing.assert_array_equal(a.observables.dipole,expected)
        self.assertEqual(a.time,.02)
        with self.assertRaises(InputError):bridge.commit()
        with self.assertRaises(InputError):_IsolatedDipoleTrial(a,DipoleMap(('a','b'),(10,),[[1.,1.]],[8.]))

    def test_worker_does_not_inherit_parent_mpi_connection(self):
        import os
        original=subprocess.Popen
        launched=[]
        def capture(args,**kw):
            launched.append(kw['env'].copy())
            return original(args,**kw)
        parent=dict(OMPI_COMM_WORLD_RANK='3',OMPI_COMM_WORLD_SIZE='4',
                    PMIX_RANK='3',PMIX_SERVER_URI41='parent-only',PMI_FD='123',
                    PRTE_LAUNCHED='1',MPI_LOCALRANKID='3',SLURM_PROCID='3',
                    SLURM_LOCALID='3',SLURM_JOB_ID='retain-allocation',
                    OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        with patch.dict(os.environ,parent),patch('salmonpy.tddft.subprocess.Popen',side_effect=capture):
            a=self.state()
            a.step([0,0,.01])
            self.assertEqual(os.environ['PMIX_RANK'],'3')
        env=launched[0]
        for key in parent:
            if key.startswith(('OMPI_','PMIX_','PMI_','PRTE_')) or key in ('MPI_LOCALRANKID','SLURM_PROCID','SLURM_LOCALID'):
                self.assertTrue(key not in env, key+' leaked into child MPI environment')
        self.assertEqual(env['SLURM_JOB_ID'],'retain-allocation')

    def test_measured_electron_number_survives_replay(self):
        a=self.state()
        self.assertEqual(a.observables.electron_number,10.)
        a.step([0.,0.,.01]);token=a.snapshot()
        expected=a.observables.electron_number
        a.step([0.,0.,.02]);a.restore(token)
        self.assertEqual(a.observables.electron_number,expected)
        self.assertEqual(a.electron_number,expected)

        import queue
        from salmonpy.errors import NumericalError
        events=a._engine._events
        try:
            for bad in ('nan','inf','-inf'):
                a._engine._events=queue.Queue()
                a._engine._events.put('SALMONPY_TDDFT_E_V3 1 0 0 0 0 0 0 '+bad+' -1 1')
                with self.assertRaises(NumericalError):a._engine._receive(1)
        finally:a._engine._events=events
        self.assertEqual(a.electron_number,expected)


    def test_energy_is_explicitly_fresh_and_restores(self):
        a=self.state()
        self.assertEqual(a.observables.internal_energy,-1.)
        a.step([0.,0.,.01]);token=a.snapshot()
        expected=a.internal_energy
        a.step([0.,0.,.02])
        self.assertIsNone(a.internal_energy)
        a.restore(token)
        self.assertEqual(a.internal_energy,expected)


    def test_midpoint_trial_retains_left_field_on_cancel(self):
        from salmonpy.dipole_source import DipoleMap,_IsolatedDipoleTrial
        from salmonpy.errors import NumericalError
        a,b=self.state(),self.state()
        bridge=_IsolatedDipoleTrial(a,DipoleMap(('molecule',),(10,),[[1.]],[8.]),sampling='midpoint')
        bridge.propose([3.,5.,7.],[1.,1.,1.])
        first=b.step([1.,2.,3.]).dipole.copy()
        np.testing.assert_array_equal(a.observables.dipole,first)
        bridge.cancel()
        self.assertEqual(a.time,0.)
        bridge.propose([3.,5.,7.],[1.,1.,1.]);bridge.validate([3.,5.,7.],[1.,1.,1.]);bridge.commit()
        bridge.propose([5.,9.,13.],[1.,1.,1.])
        expected=b.step([3.,6.,9.]).dipole.copy()
        np.testing.assert_array_equal(a.observables.dipole,expected)
        with self.assertRaises(NumericalError):bridge.validate([5.1,9.,13.],[1.,1.,1.])
        bridge.cancel()
        np.testing.assert_array_equal(a.observables.dipole,first)
        bridge.propose([5.,9.,13.],[1.,1.,1.]);bridge.validate([5.,9.,13.],[1.,1.,1.]);bridge.commit()
        np.testing.assert_array_equal(a.observables.dipole,expected)

    def test_public_feedback_start_rejection_retry_and_caller_ownership(self):
        from salmonpy.dipole_source import MolecularFeedback, DipoleMap
        from salmonpy.stepping import StepRecoveryError
        from types import SimpleNamespace
        state=self.state()
        class Comm:
            rank=0;size=2;different_footprint=False
            def allgather(self,value):
                if isinstance(value,dict):
                    other=dict(value,state_dt=None)
                    if self.different_footprint:other['footprint']=('different',)
                    return [value,other]
                return [value,[]]
            def gather(self,value,root=0):
                # A source-free own-vacuum reference; no hidden physical solver.
                return [value,dict(error=None,E=np.zeros(3),energy=None)]
            def bcast(self,value,root=0):return value
        class Backup:
            fail=False
            def __init__(self,*a):self.active=False
            def capture(self):
                self.active=True
                if self.fail:raise OSError('leader journal unavailable')
                return field.step_index
            def rollback(self,token):field.step_index=token
            def release(self,token):self.active=False
        field=SimpleNamespace(dt=state.dt,step_index=0,capacity=20,closed=False,failed=False,
            E=np.zeros((1,3)),point_ids=(10,),manifest={'cell_volume':8.,'grid_indices':[[19,20,20]],
                'sample_positions_bohr':[[-1.,1.,1.]],'sample_map':'same-map','deposit_map':'same-deposit'},
            calculation=SimpleNamespace(geometry=SimpleNamespace(grid=SimpleNamespace(counts=(40,40,40),spacing=(2.,2.,2.))),
                boundaries=SimpleNamespace(faces=('pml',)*6)),directory=self.tmp.name,
            context=SimpleNamespace(is_leader=True,solver_comm=SimpleNamespace(size=1,allgather=lambda value:[value])))
        def advance(current):field.step_index+=1
        field.advance=advance
        mapping=DipoleMap(('molecule',),(10,),[[1.]],[8.])
        with patch('salmonpy.native_steps._NativeBackup',Backup):
            comm=Comm();comm.different_footprint=True
            with self.assertRaises(InputError):
                MolecularFeedback(comm,field,state=state,mapping=mapping,total_root=0,
                    reference_root=1,policy='environment_scattered_only')
            self.assertEqual(field.step_index,0)
            comm.different_footprint=False
            session=MolecularFeedback(comm,field,state=state,mapping=mapping,total_root=0,
                reference_root=1,policy='environment_scattered_only')
            first=session.start();self.assertEqual(first.field_step,1)
            before=state.observables.dipole.copy()
            Backup.fail=True
            with self.assertRaises(StepRecoveryError) as failed_capture:session.step()
            self.assertTrue(failed_capture.exception.original_preserved)
            self.assertEqual(field.step_index,1);self.assertEqual(state.time,0.)
            Backup.fail=False
            def reject(frame):raise ValueError('injected consumer rejection')
            with self.assertRaises(StepRecoveryError) as cm:session.step(validate=reject)
            self.assertTrue(cm.exception.original_preserved)
            self.assertEqual(field.step_index,1);self.assertEqual(state.time,0.)
            np.testing.assert_array_equal(state.observables.dipole,before)
            result=session.step()
            self.assertEqual(result.field_step,2)
            self.assertEqual(result.molecule.time,state.dt)
            with self.assertRaises(ValueError):result.selected_E[0]=1.
            self.assertFalse(field.closed)
            self.assertEqual(state.step([0,0,0]).time,2*state.dt)

    def test_live_savepoint_same_process_no_replay_and_explicit_lifetime(self):
        a,b=self.state(),self.state()
        for E in ([.1,.2,0.],[.3,-.1,.2]):a.step(E)
        token=a.savepoint();pid=a.worker_statistics()['pid']
        saved=a.observables;other=b.observables.dipole.copy()
        with self.assertRaises(SnapshotError):a.savepoint()
        with self.assertRaises(SnapshotError):a.snapshot()
        with self.assertRaises(SnapshotError):b.restore_savepoint(token)
        expected=a.step([.2,.1,-.2])
        updates=a.worker_statistics()['updates']
        with patch.object(a._engine,'_spawn',side_effect=AssertionError('restore must not replay')):
            with patch.object(token._journal,'restore',side_effect=OSError('journal temporarily unavailable')):
                with self.assertRaises(OSError):a.restore_savepoint(token)
            self.assertFalse(a.usable)
            with self.assertRaises(RuntimeFailure):a.step([0.,0.,0.])
            a.restore_savepoint(token)
            self.assertEqual(a.worker_statistics()['updates'],updates)
            self.assertEqual(a.worker_statistics()['pid'],pid)
            self.assertEqual(a.time,saved.time)
            np.testing.assert_array_equal(a.observables.dipole,saved.dipole)
            np.testing.assert_array_equal(a.step([.2,.1,-.2]).dipole,expected.dipole)
            a.restore_savepoint(token)
        np.testing.assert_array_equal(b.observables.dipole,other)
        a.release_savepoint(token)
        with self.assertRaises(SnapshotError):a.restore_savepoint(token)
        newer=a.savepoint()
        with self.assertRaises(SnapshotError):a.restore_savepoint(token)
        a.release_savepoint(newer)
        pending=a.savepoint()
        a._engine._proc.kill();a._engine._proc.wait()
        with patch.object(a._engine,'_spawn',side_effect=AssertionError('dead direct worker must not replay')):
            with self.assertRaises(RuntimeFailure):a.step([0.,0.,0.])
            with self.assertRaises(RuntimeFailure):a.restore_savepoint(pending)
        self.assertFalse(a.usable)
