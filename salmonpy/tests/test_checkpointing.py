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
