"""Common checkpoint metadata validation must not invoke solver code."""
import copy
import unittest
from dataclasses import FrozenInstanceError
from salmonpy.checkpointing import CheckpointCapability, check_checkpoint_compatibility

GROUPS=[{'name':'material','role':'material','ranks':[0],'point_ids':[10]}]

class CheckpointingTests(unittest.TestCase):
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

    def test_import_without_optional_mpi_is_quiet(self):
        import subprocess
        import sys
        from pathlib import Path
        import salmonpy
        args=[sys.executable]
        if Path(salmonpy.__file__).resolve().parent != Path(__file__).resolve().parents[1]:
            args.append('-I')
        code="import sys; sys.modules['mpi4py']=None; import salmonpy.checkpointing; assert 'mpi4py.MPI' not in sys.modules"
        result=subprocess.run(args+['-c',code],text=True,capture_output=True,timeout=20)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(result.stdout,'')
        self.assertEqual(result.stderr,'')

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
