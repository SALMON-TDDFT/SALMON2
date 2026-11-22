"""Admission failures must precede native mutation."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from salmonpy.examples._direct_checkpoint import checkpoint, digest


class DirectCheckpointTests(unittest.TestCase):
    def world(self):
        return SimpleNamespace(rank=0, size=1, allgather=lambda x:[x], Barrier=lambda:None)

    def coupled(self):
        raw=SimpleNamespace(snapshot=lambda:None)
        response=SimpleNamespace(states={1:SimpleNamespace(_raw=raw)}, failed=False)
        return SimpleNamespace(response=response, field=None)

    def test_unsupported_backend_before_directory_creation(self):
        c=self.coupled();c.response.states={1:object()}
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'checkpoint'
            with self.assertRaisesRegex(ValueError,'capability'):checkpoint(self.world(),c,path,{})
            self.assertFalse(path.exists())

    def test_incomplete_rejected_before_state_change(self):
        with tempfile.TemporaryDirectory() as d:
            c=self.coupled()
            with self.assertRaises(ValueError):checkpoint(self.world(),c,d,{},restore=True)
            self.assertFalse(c.response.failed)

    def test_changed_signature_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'complete.json').write_text(json.dumps({'schema':'salmonpy.native-direct/4','signature':{'old':1},'ranks':1}))
            with self.assertRaisesRegex(ValueError,'mismatch'):checkpoint(self.world(),self.coupled(),p,{'steps':128},restore=True)

    def test_corrupted_rank_file_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);data=p/'rank-0.json';data.write_text('{}');(p/'frames.json').write_text('[]')
            m={'schema':'salmonpy.native-direct/4','signature':{'steps':128},'step':64,'ranks':1,'files':{data.name:digest(data), 'frames.json':digest(p/'frames.json')}}
            (p/'complete.json').write_text(json.dumps(m));data.write_text('changed')
            with self.assertRaisesRegex(ValueError,'checksum'):checkpoint(self.world(),self.coupled(),p,{'steps':128},restore=True)

    def test_files_hashed_once_across_ranks(self):
        from salmonpy.examples._direct_checkpoint import _verify_files
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            files={}
            for i in range(9):
                p=root/f'payload-{i}.bin';p.write_bytes(bytes([i])*(i+1));files[p.name]=digest(p)
            (root/'complete.json').write_text('{}')
            seen=[]
            def counted(path):
                seen.append(Path(path).name)
                return digest(path)
            with patch('salmonpy.examples._direct_checkpoint.digest',side_effect=counted):
                reports=[_verify_files(root,files,rank,4) for rank in range(4)]
            self.assertCountEqual(seen,files)
            self.assertEqual(sum(x['files_checked'] for x in reports),9)
            self.assertEqual(sum(x['bytes_checked'] for x in reports),45)

    def test_assigned_corruption_detected(self):
        from salmonpy.examples._direct_checkpoint import _verify_files
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'payload';p.write_bytes(b'ok');files={p.name:digest(p)}
            (root/'complete.json').write_text('{}');p.write_bytes(b'bad')
            with self.assertRaisesRegex(ValueError,'checksum'):_verify_files(root,files,0,4)

    def test_reject_path_escape_on_nonassigned_rank(self):
        from salmonpy.examples._direct_checkpoint import _verify_files
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):_verify_files(Path(d),{'../escape':'0'*64},3,4)


    def test_manifest_disagreement_rejected_before_mutation(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            data={'step':64,'J':[[0.,0.,0.]],'response_A':[[0.,0.,0.]],'states':{'1':{}}}
            (p/'rank-0.json').write_text(json.dumps(data));(p/'frames.json').write_text('[]')
            files={n:digest(p/n) for n in ('rank-0.json','frames.json')}
            (p/'complete.json').write_text(json.dumps({'schema':'salmonpy.native-direct/4',
                'signature':{'steps':128},'ranks':1,'step':64,'files':files}))
            world=self.world();world.allgather=lambda value:[value,(None,'different-manifest',value[2])]
            coupled=self.coupled();coupled.response.context=SimpleNamespace(plan=SimpleNamespace(point_ids=(1,)))
            with self.assertRaisesRegex(ValueError,'manifests differ'):
                checkpoint(world,coupled,p,{'steps':128},restore=True)
            self.assertFalse(coupled.response.failed)

    def test_worker_A_mismatch_rejected_before_loading(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            data={'step':64,'J':[[0.,0.,0.]],'response_A':[[0.,0.,0.]],
                  'states':{'1':{'A':[1.,0.,0.],'current':[0.,0.,0.]}}}
            (p/'rank-0.json').write_text(json.dumps(data));(p/'frames.json').write_text('[]')
            files={n:digest(p/n) for n in ('rank-0.json','frames.json')}
            (p/'complete.json').write_text(json.dumps({'schema':'salmonpy.native-direct/4',
                'signature':{'steps':128},'ranks':1,'step':64,'files':files}))
            coupled=self.coupled()
            coupled.response.states={1:SimpleNamespace(checkpoint_protocol='salmonpy.direct-state/1',
                checkpoint=lambda *a,**k:self.fail('must reject before native load'))}
            coupled.response.context=SimpleNamespace(plan=SimpleNamespace(point_ids=(1,)))
            with self.assertRaisesRegex(ValueError,'A differs'):
                checkpoint(self.world(),coupled,p,{'steps':128},restore=True)
            self.assertFalse(coupled.response.failed)

    def test_standard_current_rejected_before_staging(self):
        from unittest.mock import Mock
        import numpy as np
        for bad in (None, [0.,0.], [0.,float('inf'),0.]):
            with self.subTest(current=bad), tempfile.TemporaryDirectory() as d:
                p=Path(d)
                data={'step':64,'J':[[0.,0.,0.]],'response_A':[[0.,0.,0.]],
                      'states':{'1':{'A':[0.,0.,0.],'current':bad}}}
                (p/'rank-0.json').write_text(json.dumps(data));(p/'frames.json').write_text('[]')
                files={n:digest(p/n) for n in ('rank-0.json','frames.json')}
                (p/'complete.json').write_text(json.dumps({'schema':'salmonpy.native-direct/4',
                    'signature':{'steps':128},'ranks':1,'step':64,'files':files}))
                state=SimpleNamespace(native=object(),current=np.zeros(3),checkpoint=Mock(),
                                      lib=SimpleNamespace(sp_tddft_checkpoint=Mock()))
                coupled=self.coupled();coupled.response.states={1:state}
                coupled.response.context=SimpleNamespace(plan=SimpleNamespace(point_ids=(1,)))
                with self.assertRaisesRegex(ValueError,'invalid checkpoint material current'):
                    checkpoint(self.world(),coupled,p,{'steps':128},restore=True,protected=True)
                state.checkpoint.assert_not_called()
                np.testing.assert_array_equal(state.current,np.zeros(3))
                self.assertFalse(coupled.response.failed)

    def test_multiple_worker_files_cannot_collide(self):
        from salmonpy.examples._direct_checkpoint import _state_filename
        state=lambda:SimpleNamespace(checkpoint_protocol='salmonpy.direct-state/1')
        r=SimpleNamespace(states={1:state(),2:state()})
        self.assertEqual({_state_filename(r,p,0) for p in r.states},
                         {'worker-0-point-1.bin','worker-0-point-2.bin'})

    def test_restart_layout_includes_inner_decomposition_and_capability(self):
        from salmonpy.examples._direct_checkpoint import _compatibility
        from salmonpy.distributed import ParallelPlan,SolverGroup
        from unittest.mock import Mock
        group=SolverGroup('electronic',(0,),point_ids=(1,))
        state=SimpleNamespace(checkpoint_protocol='salmonpy.direct-state/1',
            checkpoint_format='custom-numeric/1',checkpoint_modes=('save','restore'),
            checkpoint_decomposition=(1,1,1,1,1),checkpoint=Mock())
        context=SimpleNamespace(group=group,plan=ParallelPlan((group,SolverGroup('field',(1,),role='field'))))
        coupled=SimpleNamespace(response=SimpleNamespace(states={1:state},context=context),field=None)
        record=_compatibility(self.compat_world(),coupled,restore=False)
        self.assertEqual(record,_compatibility(self.compat_world(),coupled,restore=True,saved=record))
        bad=json.loads(json.dumps(record));bad['participants'][0]['states']['1']['decomposition']=[1,2,1,1,1]
        with self.assertRaisesRegex(ValueError,'compatibility'):
            _compatibility(self.compat_world(),coupled,restore=True,saved=bad)
        state.checkpoint_modes=('save',)
        with self.assertRaisesRegex(ValueError,'capability'):
            _compatibility(self.compat_world(),coupled,restore=True,saved=record)
        state.checkpoint_modes=('save','restore');state.steps=1
        with self.assertRaisesRegex(ValueError,'fresh'):
            _compatibility(self.compat_world(),coupled,restore=True,saved=record)
        state.checkpoint.assert_not_called()

    def test_snapshot_alone_does_not_prove_direct_restart_support(self):
        from salmonpy.examples._direct_checkpoint import _compatibility
        from salmonpy.distributed import ParallelPlan,SolverGroup
        from unittest.mock import Mock
        group=SolverGroup('electronic',(0,),point_ids=(1,))
        raw=SimpleNamespace(snapshot=Mock(),restore=Mock())
        coupled=SimpleNamespace(field=None,response=SimpleNamespace(states={1:SimpleNamespace(_raw=raw)},
            context=SimpleNamespace(group=group,plan=ParallelPlan((group,SolverGroup('field',(1,),role='field'))))))
        with self.assertRaisesRegex(ValueError,'capability'):
            _compatibility(self.compat_world(),coupled,restore=False)
        raw.snapshot.assert_not_called();raw.restore.assert_not_called()

    def compat_world(self):
        def gather(value):
            return [value,(None,value[1],{'rank':1,'group':'field','states':{},
                'field':{'format':'field/1','modes':['restore','save'],'decomposition':[1]}},value[3],value[4])]
        return SimpleNamespace(rank=0,size=2,allgather=gather)

    def test_multiple_native_files_cannot_collide(self):
        from salmonpy.examples._direct_checkpoint import _state_filename
        r=SimpleNamespace(states={1:SimpleNamespace(native=object()),2:SimpleNamespace(native=object())})
        self.assertEqual({_state_filename(r,p,0) for p in r.states},
                         {'electronic-0-point-1.bin','electronic-0-point-2.bin'})

    def test_public_capability_object_on_custom_backend(self):
        from salmonpy.checkpointing import CheckpointCapability
        from salmonpy.distributed import ParallelPlan,SolverGroup
        from salmonpy.examples._direct_checkpoint import _compatibility
        from unittest.mock import Mock
        group=SolverGroup('electronic',(0,),point_ids=(1,))
        state=SimpleNamespace(checkpoint_protocol='salmonpy.direct-state/1',
            checkpoint_capability=CheckpointCapability('custom/1',('save','restore'),(1,),identity='model-v1'),
            checkpoint=Mock())
        coupled=SimpleNamespace(field=None,response=SimpleNamespace(states={1:state},
            context=SimpleNamespace(group=group,plan=ParallelPlan((group,SolverGroup('field',(1,),role='field'))))))
        record=_compatibility(self.compat_world(),coupled,restore=False)
        self.assertEqual(record['participants'][0]['states']['1']['identity'],'model-v1')
        state.checkpoint.assert_not_called()

    def test_protected_policy_agreed_before_branch(self):
        from unittest.mock import patch
        world=self.world()
        world.allgather=lambda value:[value,(None,value[1],(False,False))]
        coupled=self.coupled()
        with tempfile.TemporaryDirectory() as root:
            with patch('salmonpy.examples._direct_checkpoint._compatibility') as admission:
                with self.assertRaisesRegex(ValueError,'checkpoint policies differ'):
                    checkpoint(world,coupled,Path(root)/'new',{},protected=True)
                admission.assert_not_called()
        self.assertFalse(coupled.response.failed)

    def test_binary_names_include_rank_and_point(self):
        from salmonpy.examples._direct_checkpoint import _state_filename
        worker=lambda:SimpleNamespace(checkpoint_protocol='salmonpy.direct-state/1')
        a,b=worker(),worker()
        r=SimpleNamespace(states={10:a,20:b})
        self.assertEqual(_state_filename(r,10,3),'worker-3-point-10.bin')
        self.assertNotEqual(_state_filename(r,10,3),_state_filename(r,20,3))
        self.assertEqual(_state_filename(SimpleNamespace(states={10:a}),10,3),'worker-3-point-10.bin')

    def test_version_four_identity_cannot_be_silently_downgraded(self):
        from salmonpy.examples._mixed_identity import check_identity
        from salmonpy.errors import InputError
        identity={'run':'example','points':[10,20]}
        check_identity(identity,identity,'salmonpy.native-direct/4')
        check_identity(None,None,'salmonpy.native-direct/4')
        with self.assertRaises(InputError):check_identity(identity,None,'salmonpy.native-direct/4')
        with self.assertRaises(InputError):check_identity(None,identity,'salmonpy.native-direct/4')
