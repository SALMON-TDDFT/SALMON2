import importlib.util
import unittest

class DistributedTest(unittest.TestCase):
    def test_plan_validation(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.distributed'), 'MPI group interface is missing')
        from salmonpy.distributed import SolverGroup, ParallelPlan
        p=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('material',(1,2),point_ids=(7,9))))
        p.validate(3)
        self.assertEqual(p.point_ids,(7,9))
        with self.assertRaises(ValueError): p.validate(4)
        with self.assertRaises(ValueError):
            ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('bad',(0,1),point_ids=(1,)))).validate(2)
        with self.assertRaises(ValueError):
            ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('a',(1,),point_ids=(3,)),SolverGroup('b',(2,),point_ids=(3,)))).validate(3)
        with self.assertRaises(ValueError): SolverGroup('bad',(True,))

    def test_import_does_not_initialize_mpi(self):
        import subprocess,sys
        code="import sys; import salmonpy.distributed; assert 'mpi4py.MPI' not in sys.modules"
        self.assertEqual(subprocess.run([sys.executable,'-c',code],capture_output=True).returncode,0)

    def test_ld_independent_a_response(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.material_response'))
        import salmonpy.material_response as m
        self.assertTrue(hasattr(m,'LDMaterialProvider'),'LD provider is missing')
        import numpy as np
        from salmonpy.ld import LorentzDrude
        model=LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[0.],damping=[.1])
        provider=m.LDMaterialProvider(model)
        initial=provider.prepare(None);backend=m.LDResponseBackend()
        one=backend.create_state(initial,1,.01,None);two=backend.create_state(initial,2,.01,None)
        reference=model.new_state(dt=.01)
        for old,new in [(np.zeros(3),np.array([.01,0.,0.])),(np.array([.01,0.,0.]),np.array([.03,0.,0.]))]:
            np.testing.assert_array_equal(one.advance(old,new),reference.step(-(new-old)/.01)[1])
        np.testing.assert_array_equal(two.advance(np.zeros(3),np.zeros(3)),np.zeros(3))
        one.close();two.close()

    def test_tddft_backend_rejects_other_initial_types(self):
        from salmonpy.material_response import TDDFTResponseBackend
        with self.assertRaises(TypeError):
            TDDFTResponseBackend().create_state(object(),1,.02,None)
