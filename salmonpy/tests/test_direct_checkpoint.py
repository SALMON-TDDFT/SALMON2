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
            p=Path(d);(p/'complete.json').write_text(json.dumps({'schema':'salmonpy.native-direct/1','signature':{'old':1},'ranks':1}))
            with self.assertRaisesRegex(ValueError,'mismatch'):checkpoint(self.world(),self.coupled(),p,{'steps':128},restore=True)

    def test_corrupted_rank_file_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);data=p/'rank-0.json';data.write_text('{}');(p/'frames.json').write_text('[]')
            m={'schema':'salmonpy.native-direct/1','signature':{'steps':128},'step':64,'ranks':1,'files':{data.name:digest(data), 'frames.json':digest(p/'frames.json')}}
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
            (p/'complete.json').write_text(json.dumps({'schema':'salmonpy.native-direct/1',
                'signature':{'steps':128},'ranks':1,'step':64,'files':files}))
            world=self.world();world.allgather=lambda value:[value,(None,'different-manifest')]
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
            (p/'complete.json').write_text(json.dumps({'schema':'salmonpy.native-direct/1',
                'signature':{'steps':128},'ranks':1,'step':64,'files':files}))
            coupled=self.coupled()
            coupled.response.states={1:SimpleNamespace(checkpoint_protocol='salmonpy.direct-state/1',
                checkpoint=lambda *a,**k:self.fail('must reject before native load'))}
            coupled.response.context=SimpleNamespace(plan=SimpleNamespace(point_ids=(1,)))
            with self.assertRaisesRegex(ValueError,'A differs'):
                checkpoint(self.world(),coupled,p,{'steps':128},restore=True)
            self.assertFalse(coupled.response.failed)

    def test_multiple_worker_files_cannot_collide(self):
        c=self.coupled()
        state=lambda:SimpleNamespace(checkpoint_protocol='salmonpy.direct-state/1',checkpoint=lambda *a,**k:None)
        c.response.states={1:state(),2:state()}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'checkpoint'
            with self.assertRaisesRegex(ValueError,'one worker checkpoint'):
                checkpoint(self.world(),c,p,{})
            self.assertFalse(p.exists())

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
            return [value,(None,value[1],{'rank':1,'group':'field','states':{},'field':{'format':'field/1'}})]
        return SimpleNamespace(rank=0,size=2,allgather=gather)

    def test_multiple_native_files_cannot_collide(self):
        coupled=self.coupled()
        coupled.response.states={1:SimpleNamespace(native=object()),2:SimpleNamespace(native=object())}
        with tempfile.TemporaryDirectory() as d:
            target=Path(d)/'checkpoint'
            with self.assertRaisesRegex(ValueError,'one native electronic checkpoint'):
                checkpoint(self.world(),coupled,target,{})
            self.assertFalse(target.exists())
