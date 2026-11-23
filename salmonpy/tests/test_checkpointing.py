"""Common checkpoint metadata validation must not invoke solver code."""
import copy
import unittest
from dataclasses import FrozenInstanceError
from salmonpy.checkpointing import CheckpointCapability, check_checkpoint_compatibility

GROUPS=[{'name':'material','role':'material','ranks':[0],'point_ids':[10]}]

class CheckpointingTests(unittest.TestCase):
    def test_timed_checkpoint_and_periodic_save(self):
        from unittest.mock import patch
        from salmonpy.checkpointing import CheckpointSchedule
        saved=[]
        with patch('salmonpy.checkpointing.time.monotonic',side_effect=[10.,11.,13.,15.]):
            schedule=CheckpointSchedule(interval_steps=2,interval_seconds=3.,stop_after_seconds=5.)
            self.assertFalse(schedule.after_step(1,saved.append))
            self.assertFalse(schedule.after_step(2,saved.append))
            self.assertTrue(schedule.after_step(3,saved.append))
        self.assertEqual(saved,[2,3])
        self.assertEqual(schedule.last_saved_step,3)
        self.assertTrue(schedule.stopped)
        with self.assertRaises(ValueError):schedule.after_step(4,saved.append)
        disabled=CheckpointSchedule()
        self.assertFalse(disabled.after_step(1,saved.append))
        self.assertEqual(saved,[2,3])
        for kwargs in ({'interval_steps':0},{'stop_after_seconds':float('nan')},
                       {'interval_seconds':-1},{'initial_step':True}):
            with self.assertRaises(ValueError):CheckpointSchedule(**kwargs)

    def test_timed_checkpoint_save_failure_can_retry(self):
        from unittest.mock import patch
        from salmonpy.checkpointing import CheckpointSchedule
        from salmonpy.errors import RuntimeFailure
        with patch('salmonpy.checkpointing.time.monotonic',return_value=0.):
            schedule=CheckpointSchedule(stop_after_seconds=1.)
        def fail(step):raise OSError('disk full')
        with patch('salmonpy.checkpointing.time.monotonic',return_value=2.):
            with self.assertRaisesRegex(RuntimeFailure,'disk full'):schedule.after_step(7,fail)
            self.assertFalse(schedule.stopped)
            self.assertIsNone(schedule.last_saved_step)
            saved=[]
            self.assertTrue(schedule.after_step(7,saved.append))
            self.assertEqual(saved,[7])

    def test_serial_roundtrip_is_detached_and_identity_bound(self):
        cap=CheckpointCapability('user-drude/1',('save','restore'),(1,),identity='model-sha256')
        with self.assertRaises(FrozenInstanceError): cap.format='changed'
        item=cap.to_dict();item['decomposition'][0]=2
        self.assertEqual(cap.decomposition,(1,))
        saved=check_checkpoint_compatibility(None,GROUPS,states={10:cap})
        original=copy.deepcopy(saved)
        self.assertEqual(saved,check_checkpoint_compatibility(None,GROUPS,states={10:cap.to_dict()},operation='restore',saved=saved))
        self.assertEqual(saved,original)
        other=CheckpointCapability('user-drude/1',('save','restore'),(1,),identity='different-model')
        with self.assertRaisesRegex(ValueError,'compatibility'):
            check_checkpoint_compatibility(None,GROUPS,states={10:other},operation='restore',saved=saved)

    def test_bad_declaration_and_inventory_rejected_without_io(self):
        cap=CheckpointCapability('user/1',('save',),(1,))
        with self.assertRaisesRegex(ValueError,'capability'):
            check_checkpoint_compatibility(None,GROUPS,states={10:cap},operation='restore',saved={})
        both=CheckpointCapability('user/1',('save','restore'),(1,))
        saved=check_checkpoint_compatibility(None,GROUPS,states={10:both})
        with self.assertRaisesRegex(ValueError,'required'):
            check_checkpoint_compatibility(None,GROUPS,states={10:both},operation='restore',saved=saved,
                required_files=('rank-0.bin',),available_files=())
        for invalid in ('../escape','/absolute','a\\b','./file'):
            with self.assertRaises(ValueError):
                check_checkpoint_compatibility(None,GROUPS,states={10:both},required_files=(invalid,))

    def test_unknown_and_malformed_contracts_are_not_guessed(self):
        for declaration in ({}, {'format':'user/1','modes':['save','restore'],'decomposition':[2]},
                            {'format':'user/1','modes':['save','restore'],'decomposition':[True]},
                            {'format':'user/1','modes':['save','restore'],'decomposition':[1],'typo':1}):
            with self.subTest(declaration=declaration),self.assertRaises(ValueError):
                check_checkpoint_compatibility(None,GROUPS,states={10:declaration})
        with self.assertRaises(ValueError):
            check_checkpoint_compatibility(None,GROUPS,states={})

    def test_custom_backend_example_continues_without_replay(self):
        import tempfile
        from pathlib import Path
        from salmonpy.examples.custom_checkpoint import run
        with tempfile.TemporaryDirectory() as root:
            result=run(Path(root)/'new-run')
        self.assertEqual(result['maximum_error'],0.)
        self.assertEqual(result['new_steps_after_load'],15)
        self.assertTrue(result['independent_states'])

    def test_saved_numeric_type_changes_are_not_silently_coerced(self):
        cap=CheckpointCapability('user/1',('save','restore'),(1,))
        saved=check_checkpoint_compatibility(None,GROUPS,states={10:cap})
        saved['participants'][0]['states']['10']['decomposition']=[True]
        with self.assertRaisesRegex(ValueError,'compatibility'):
            check_checkpoint_compatibility(None,GROUPS,states={10:cap},operation='restore',saved=saved)


    def test_empty_local_exception_is_still_reported_as_input_error(self):
        from salmonpy.errors import InputError
        class BrokenGroup(dict):
            def __iter__(self):
                raise ValueError()
        cap=CheckpointCapability('user/1',('save','restore'),(1,))
        with self.assertRaises(InputError):
            check_checkpoint_compatibility(None,[BrokenGroup(GROUPS[0])],states={10:cap})

    def test_execution_returns_backend_owned_tokens_and_restores(self):
        from salmonpy.checkpointing import CheckpointTarget, execute_checkpoint
        cap=CheckpointCapability('opaque/1',('save','restore'),(1,),identity='model')
        token=object();calls=[]
        target=CheckpointTarget(cap,save=lambda:token,restore=lambda:calls.append('restore'))
        saved=execute_checkpoint(None,GROUPS,states={10:target},operation='save')
        self.assertIs(saved['states']['10'],token)
        execute_checkpoint(None,GROUPS,states={10:target},operation='restore',saved=saved['compatibility'])
        self.assertEqual(calls,['restore'])

    def test_bad_callback_or_compatibility_runs_no_backend(self):
        from salmonpy.checkpointing import CheckpointTarget, execute_checkpoint
        from unittest.mock import Mock
        cap=CheckpointCapability('user/1',('save','restore'),(1,))
        callback=Mock()
        target=CheckpointTarget(cap,save=callback)
        with self.assertRaisesRegex(ValueError,'callback'):
            execute_checkpoint(None,GROUPS,states={10:target},operation='restore',saved={})
        target=CheckpointTarget(cap,save=callback,restore=callback)
        with self.assertRaisesRegex(ValueError,'compatibility'):
            execute_checkpoint(None,GROUPS,states={10:target},operation='restore',saved={})
        callback.assert_not_called()

    def test_callback_failure_stops_later_point_and_preserves_error_family(self):
        from salmonpy.checkpointing import CheckpointTarget, execute_checkpoint
        from salmonpy.errors import RuntimeFailure
        from unittest.mock import Mock
        groups=[{'name':'user','ranks':[0],'point_ids':['a','b']}]
        cap=CheckpointCapability('user/1',('save','restore'),(1,))
        def broken(): raise ValueError('disk full')
        later=Mock()
        with self.assertRaisesRegex(RuntimeFailure,'disk full'):
            execute_checkpoint(None,groups,states={'a':CheckpointTarget(cap,save=broken),
                'b':CheckpointTarget(cap,save=later)},operation='save')
        later.assert_not_called()

    def test_ld_snapshot_payload_remains_backend_owned_and_exact(self):
        import numpy as np
        from salmonpy.ld import LorentzDrude
        from salmonpy.checkpointing import CheckpointTarget, execute_checkpoint
        material=LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[.5],damping=[.1])
        cap=CheckpointCapability('ld-token/1',('save','restore'),(1,))
        with material.new_state(dt=.01) as original, material.new_state(dt=.01) as restored:
            original.step([.01,0,0])
            saved=execute_checkpoint(None,GROUPS,states={10:CheckpointTarget(cap,save=original.snapshot)})
            token=saved['states']['10']
            execute_checkpoint(None,GROUPS,operation='restore',saved=saved['compatibility'],
                states={10:CheckpointTarget(cap,restore=lambda:restored.restore(token))})
            for n in range(2,8):
                expected=original.step([n*.01,0,0])
                actual=restored.step([n*.01,0,0])
                for a,b in zip(actual,expected): np.testing.assert_array_equal(a,b)

    def test_atomic_failure_preserves_all_originals_and_allows_continuation(self):
        from salmonpy.checkpointing import CheckpointState, CheckpointTarget, restore_checkpoint, AtomicRestoreError
        groups=[{'name':'user','ranks':[0],'point_ids':['a','b']}]
        cap=CheckpointCapability('candidate/1',('save','restore'),(1,))
        old_a={'step':7};old_b={'step':9}
        live=CheckpointState(states={'a':old_a,'b':old_b})
        discarded=[]
        def fail():
            candidate={'step':100}
            candidate['partial']=True
            raise ValueError('failed after partial candidate load')
        targets={'a':CheckpointTarget(cap,stage_restore=lambda:{'step':10},discard_restore=lambda c:discarded.append(c)),
                 'b':CheckpointTarget(cap,stage_restore=fail,discard_restore=lambda c:discarded.append(c))}
        saved=check_checkpoint_compatibility(None,groups,states={'a':cap,'b':cap})
        with self.assertRaises(AtomicRestoreError) as caught:
            restore_checkpoint(None,groups,live=live,states=targets,saved=saved)
        self.assertTrue(caught.exception.original_preserved)
        self.assertIs(live.states['a'],old_a);self.assertIs(live.states['b'],old_b)
        self.assertEqual((old_a,old_b),({'step':7},{'step':9}))
        self.assertEqual(discarded,[None,{'step':10}])
        live.states['a']['step']+=1
        self.assertEqual(old_a['step'],8)

    def test_atomic_success_publishes_candidates_together(self):
        from salmonpy.checkpointing import CheckpointState, CheckpointTarget, restore_checkpoint
        cap=CheckpointCapability('candidate/1',('save','restore'),(1,))
        original={'step':3};candidate={'step':17}
        live=CheckpointState(states={10:original})
        saved=check_checkpoint_compatibility(None,GROUPS,states={10:cap})
        run=restore_checkpoint(None,GROUPS,live=live,saved=saved,
            states={10:CheckpointTarget(cap,stage_restore=lambda:candidate)})
        self.assertIs(live.states['10'],candidate)
        self.assertIs(run['previous']['states']['10'],original)
        self.assertEqual(original,{'step':3})

    def test_atomic_rejects_direct_only_and_live_alias(self):
        from salmonpy.checkpointing import CheckpointState, CheckpointTarget, restore_checkpoint, AtomicRestoreError
        from unittest.mock import Mock
        cap=CheckpointCapability('candidate/1',('save','restore'),(1,))
        original={};live=CheckpointState(states={10:original})
        saved=check_checkpoint_compatibility(None,GROUPS,states={10:cap})
        direct=Mock()
        with self.assertRaises(ValueError):
            restore_checkpoint(None,GROUPS,live=live,saved=saved,states={10:CheckpointTarget(cap,restore=direct)})
        direct.assert_not_called()
        with self.assertRaises(AtomicRestoreError):
            restore_checkpoint(None,GROUPS,live=live,saved=saved,
                states={10:CheckpointTarget(cap,stage_restore=lambda:original)})
        self.assertIs(live.states['10'],original)

    def test_atomic_cleanup_failure_reports_preserved_original(self):
        from salmonpy.checkpointing import CheckpointState, CheckpointTarget, restore_checkpoint, AtomicRestoreError
        groups=[{'name':'user','ranks':[0],'point_ids':['a','b']}]
        cap=CheckpointCapability('candidate/1',('save','restore'),(1,))
        live=CheckpointState(states={'a':{'step':1},'b':{'step':2}})
        original=live.states
        saved=check_checkpoint_compatibility(None,groups,states={'a':cap,'b':cap})
        def fail(): raise OSError('candidate load error')
        def cleanup(candidate): raise OSError('candidate cleanup error')
        with self.assertRaises(AtomicRestoreError) as caught:
            restore_checkpoint(None,groups,live=live,saved=saved,states={
                'a':CheckpointTarget(cap,stage_restore=lambda:{'step':17},discard_restore=cleanup),
                'b':CheckpointTarget(cap,stage_restore=fail)})
        self.assertIs(live.states,original)
        self.assertTrue(caught.exception.original_preserved)
        self.assertIn('cleanup error',str(caught.exception.cleanup_errors))

    def test_atomic_ld_state_candidate_continues_exactly(self):
        import numpy as np
        from salmonpy.ld import LorentzDrude
        from salmonpy.checkpointing import CheckpointState, CheckpointTarget, execute_checkpoint, restore_checkpoint
        material=LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[.5],damping=[.1])
        cap=CheckpointCapability('ld-token/1',('save','restore'),(1,))
        old=material.new_state(dt=.01)
        live=CheckpointState(states={10:old})
        old.step([.01,0,0])
        saved=execute_checkpoint(None,GROUPS,states={10:CheckpointTarget(cap,save=old.snapshot)})
        token=saved['states']['10']
        def stage():
            candidate=material.new_state(dt=.01)
            try: candidate.restore(token)
            except BaseException:
                candidate.close();raise
            return candidate
        try:
            restored=restore_checkpoint(None,GROUPS,live=live,saved=saved['compatibility'],
                states={10:CheckpointTarget(cap,stage_restore=stage,discard_restore=lambda c:c.close() if c is not None else None)})
            self.assertIsNot(live.states['10'],old)
            self.assertIs(restored['previous']['states']['10'],old)
            for n in range(2,8):
                for a,b in zip(old.step([n*.01,0,0]),live.states['10'].step([n*.01,0,0])):
                    np.testing.assert_array_equal(a,b)
        finally:
            if live.states['10'] is not old: live.states['10'].close()
            old.close()

    def test_atomic_sbe_abort_and_retry_preserve_future_currents(self):
        import numpy as np
        from salmonpy.examples.mpi_mixed_multiscale import models
        from salmonpy.checkpointing import CheckpointState, CheckpointTarget, restore_checkpoint, AtomicRestoreError
        _,model=models()
        cap=CheckpointCapability('sbe-token/1',('save','restore'),(1,))
        other=CheckpointCapability('counter/1',('save','restore'),(1,))
        groups=[{'name':'user','ranks':[0],'point_ids':[10,20]}]
        old=model.new_state(dt=.01);reference=model.new_state(dt=.01)
        wave=[[0,0,.0001*n] for n in range(1,13)]
        expected=[reference.step(A).copy() for A in wave]
        for A in wave[:5]: old.step(A)
        token=old.snapshot()
        for A in wave[5:8]: old.step(A)
        counter={'step':8};live=CheckpointState(states={10:old,20:counter})
        saved=check_checkpoint_compatibility(None,groups,states={10:cap,20:other})
        def stage():
            c=model.new_state(dt=.01)
            try:c.restore(token)
            except BaseException:c.close();raise
            return c
        def fail():raise OSError('second point failed')
        target=CheckpointTarget(cap,stage_restore=stage,discard_restore=lambda c:c.close() if c is not None else None)
        before=old.snapshot()
        try:
            with self.assertRaises(AtomicRestoreError):
                restore_checkpoint(None,groups,live=live,saved=saved,
                    states={10:target,20:CheckpointTarget(other,stage_restore=fail)})
            self.assertIs(live.states['10'],old)
            self.assertEqual(old.snapshot()._integrity,before._integrity)
            self.assertEqual(counter,{'step':8})
            np.testing.assert_array_equal(old.step(wave[8]),expected[8])
            restore_checkpoint(None,groups,live=live,saved=saved,
                states={10:target,20:CheckpointTarget(other,stage_restore=lambda:{'step':5})})
            for A,j in zip(wave[5:],expected[5:]):
                np.testing.assert_array_equal(live.states['10'].step(A),j)
            self.assertEqual(counter,{'step':8})
        finally:
            if live.states['10'] is not old:live.states['10'].close()
            old.close();reference.close()
