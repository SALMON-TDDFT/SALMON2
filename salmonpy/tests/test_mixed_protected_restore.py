"""Protected mixed restore must preserve live scientific states on late failure."""
import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np
from salmonpy.examples.mpi_mixed_multiscale import models
from salmonpy.materials import bind_material
from salmonpy.material_response import _StepAdapter
from salmonpy.checkpointing import CheckpointCapability, check_checkpoint_compatibility, AtomicRestoreError
from salmonpy.examples._direct_checkpoint import _restore_protected


class MixedProtectedRestoreTests(unittest.TestCase):
    def test_late_bad_sbe_candidate_keeps_all_originals_and_retry_exact(self):
        from salmonpy.material_response import SolverContext
        _, model = models()
        binding = bind_material(model,key='test',identity='same-model')
        context = SolverContext(SimpleNamespace(rank=0,size=1),'material',(10,20))
        initial = binding.provider.prepare(context)
        states={p:binding.backend.create_state(initial,p,.003,context) for p in (10,20)}
        groups=[{'name':'material','ranks':[0],'point_ids':[10,20]}]
        cap=CheckpointCapability('salmonpy.sbe-snapshot/1',('save','restore'),(1,))
        compatibility=check_checkpoint_compatibility(None,groups,states={p:cap for p in states})
        world=SimpleNamespace(rank=0,size=1,allgather=lambda x:[x])
        data={'step':3,'states':{},'J':np.zeros((2,3)).tolist(),'response_A':np.zeros((2,3)).tolist()}
        try:
            for n in range(1,4):
                for p,s in states.items():s.step(np.array([n*p*1e-5,0.,0.]))
            for p,s in states.items():
                token=s._raw.snapshot()
                data['response_A'][list(states).index(p)]=[3*p*1e-5,0.,0.]
                data['states'][str(p)]={'step':token.step_index,'time':token.time,'settings':token._settings,
                    'integrity':token._integrity,'data':[(a,b,c,d,e.hex()) for a,b,c,d,e in token._data]}
            for n in range(4,6):
                for p,s in states.items():s.step(np.array([n*p*1e-5,0.,0.]))
            response=SimpleNamespace(states=states,_adapters={p:_StepAdapter(s) for p,s in states.items()},
                _A=np.zeros((2,3)),step_index=5,failed=False,context=SimpleNamespace(group=SimpleNamespace(name='material'),
                plan=SimpleNamespace(point_ids=(10,20))))
            coupled=SimpleNamespace(field=None,response=response,J=np.ones((2,3)),step_index=5)
            original={p:s._raw.snapshot() for p,s in states.items()}
            live=states[10]._prepared.live.copy()
            valid=data['states']['20']['integrity']
            data['states']['20']['integrity']='bad'
            with tempfile.TemporaryDirectory() as out:
                with self.assertRaises(AtomicRestoreError):
                    _restore_protected(world,coupled,Path(out),data,compatibility,[])
                self.assertIs(response.states,states)
                self.assertEqual(response.step_index,5)
                self.assertFalse(response.failed)
                self.assertEqual(states[10]._prepared.live,live)
                for p,s in states.items():self.assertEqual(s._raw.snapshot()._integrity,original[p]._integrity)
                data['states']['20']['integrity']=valid
                _restore_protected(world,coupled,Path(out),data,compatibility,[])
                self.assertEqual(coupled.step_index,3)
                for p,s in response.states.items():
                    reference=model.new_state(dt=.003)
                    try:
                        for n in range(1,7):
                            expected=reference.step([n*p*1e-5,0.,0.])
                            if n>3:np.testing.assert_array_equal(s.step(np.array([n*p*1e-5,0.,0.])),expected)
                    finally:reference.close()
        finally:
            for s in getattr(locals().get('response',None),'states',states).values():s.close()
