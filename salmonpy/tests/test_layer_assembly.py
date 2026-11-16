import copy
import unittest
from types import SimpleNamespace
import numpy as np
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy.fields import AcosPulse
from salmonpy.fdtd import PlaneWave
from salmonpy.ld import LorentzDrude
from salmonpy.material_response import MaterialBinding
from salmonpy.layers import assemble_native_layers


def layout():
    return {'schema': 'salmonpy.native-layers/1', 'length_unit': 'bohr', 'axis': 0,
            'origin_cell': [19,20,20], 'layers': [
                {'id':'s','material_id':'sbe','thickness':2.,'transverse_cells':[1,1]},
                {'id':'t','material_id':'td','thickness':2.,'transverse_cells':[1,1]},
                {'id':'gap','material_id':None,'thickness':4.,'transverse_cells':[1,1]},
                {'id':'a','material_id':'drude','thickness':2.,'transverse_cells':[3,3]},
                {'id':'b','material_id':'lorentz','thickness':2.,'transverse_cells':[3,3]}]}


class LayerAssemblyTests(unittest.TestCase):
    def materials(self):
        provider=SimpleNamespace(identity='fixture',prepare=lambda c:self.fail('must not prepare'))
        backend=SimpleNamespace(identity='fixture',execution='leader_only',reduction='leader')
        return {'sbe':MaterialBinding('sbe',provider,backend),'td':MaterialBinding('td',provider,backend),
                'drude':LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),
                'lorentz':LorentzDrude(plasma_frequency=15.,strength=[.7,.3],resonance=[25.,35.],damping=[2.,3.])}

    def assemble(self, data=None, **kwargs):
        options=dict(field_grid=(1,2,1),ranks_per_point={'sbe':1,'td':1},rank_budget=4)
        options.update(kwargs)
        return assemble_native_layers(layout() if data is None else data,self.materials(),
            grid=Grid((40,40,40),(2.,2.,2.)),dt=.003,steps=128,
            source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.),**options)

    def test_matches_manual_geometry_and_plan(self):
        a=self.assemble();a.parallel_plan.validate(4)
        self.assertEqual(dict(a.points),{1:(19,20,20),2:(20,20,20)})
        expected=Geometry(a.calculation.geometry.grid,(Shape('box',(2.,4.,4.),(7.,1.,1.),material=1),
                                                     Shape('box',(2.,4.,4.),(9.,1.,1.),material=2)))
        np.testing.assert_array_equal(a.calculation.geometry.material_ids(),expected.material_ids())
        self.assertEqual(a.record()['required_ranks'],4)
        with self.assertRaises(TypeError):a.points[1]=(0,0,0)

    def test_rounding_and_identity_are_explicit(self):
        d=layout();d['layers'][0]['thickness']=2.1
        a=self.assemble(d,rank_budget=5)
        self.assertEqual(a.record()['layers'][0]['actual_thickness'],4.)
        self.assertEqual(a.record()['required_ranks'],5)
        self.assertIs(a.bindings[1],a.bindings[2])
        self.assertNotEqual(a.record()['points'][0]['point_id'],a.record()['points'][1]['point_id'])
        d['layers'][0]['id']='changed'
        self.assertEqual(a.record()['layers'][0]['id'],'s')

    def test_budget_never_silently_subsamples(self):
        d=layout();d['layers'][0]['transverse_cells']=[3,3]
        with self.assertRaisesRegex(ValueError,'required.*12|12.*budget'):self.assemble(d)

    def test_forbidden_overlap_and_bad_inputs(self):
        d=layout();d['layers'].pop(2)
        with self.assertRaisesRegex(ValueError,'LD|constitutive'):self.assemble(d)
        for key,value in [('length_unit','angstrom'),('axis',True),('unknown',1)]:
            d=layout();d[key]=value
            with self.assertRaises(ValueError):self.assemble(d)
        d=layout();d['layers'][0]['thickness']=0
        with self.assertRaises(ValueError):self.assemble(d)
        d=layout();d['layers'][1]['id']='s'
        with self.assertRaises(ValueError):self.assemble(d)

    def test_material_and_rank_contracts(self):
        d=layout();d['layers'][0]['material_id']='missing'
        with self.assertRaises(ValueError):self.assemble(d)
        with self.assertRaises(ValueError):self.assemble(ranks_per_point={'sbe':2,'td':1},rank_budget=5)
        with self.assertRaises(ValueError):self.assemble(ranks_per_point={'sbe':1,'td':1,'drude':1})

    def test_declared_collective_group_sizes_are_checked_before_execution(self):
        from unittest.mock import patch
        materials=self.materials()
        backend=SimpleNamespace(identity='custom-k-only', execution='collective', reduction='leader',
                                supported_group_sizes=(1,2))
        materials['td']=MaterialBinding('td',materials['td'].provider,backend)
        with patch.object(self,'materials',return_value=materials):
            with self.assertRaisesRegex(ValueError,'td.*supported_group_sizes'):
                self.assemble(ranks_per_point={'sbe':1,'td':3},rank_budget=6)
            valid=self.assemble(ranks_per_point={'sbe':1,'td':2},rank_budget=5)
            self.assertEqual(valid.record()['required_ranks'],5)
            for invalid in ((),(True,),(0,),(1,1),'1,2'):
                backend.supported_group_sizes=invalid
                with self.assertRaisesRegex(ValueError,'supported_group_sizes'):
                    self.assemble()
