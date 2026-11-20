import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from salmonpy.config import resolve_materials
from salmonpy.errors import InputError


def ld_recipe():
    return {'schema':'salmonpy.materials/1','materials':{
        key:{'backend':'salmonpy.ld','model':'ld','parameters':{
            'units':'atomic_gaussian','plasma_frequency':wp,
            'strength':[1.0],'resonance':[.3],'damping':[.1]}}
        for key,wp in [('metal-A',.5),('metal-B',.7)]}}


def sbe_files(root):
    (root/'si_k.data').write_text('# header\n'*5+'1 0 0 0 1.0\n')
    (root/'si_eigen.data').write_text('# eigen\n# occupations\n# esp[a.u.]\nk= 1, spin= 1\n1 -0.5 2\n2 0.4 0\n')
    text='# header\n'*3
    for part in range(2):
        if part: text+='# nonlocal\n'
        for b in range(2):
            for c in range(2):
                v=.3 if b!=c and part==0 else 0.
                text+=f'1 {b+1} {c+1} {v} 0 0 0 0 0\n'
    (root/'si_tm.data').write_text(text)
    return {'schema':'salmonpy.materials/1','materials':{'silicon':{
        'backend':'salmonpy.sbe','model':'sbe','parameters':{
            'units':'atomic','directory':'.','sysname':'si','bands':2,'kpoints':1,
            'electrons':2,'lattice':(np.eye(3)*4).tolist()}}}}


class BuiltinMaterialsTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.builtin_materials'),
                             'built-in configuration adapters missing')
        from salmonpy import builtin_materials
        self.api=builtin_materials
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.registry=self.api.builtin_material_registry()

    def resolve(self,data):
        return resolve_materials(data,registry=self.registry,base_directory=self.root)

    def build(self,catalog,key):
        return self.api.material_from_catalog(catalog,key,base_directory=self.root)

    def test_ld_background_recipe_and_provider(self):
        from salmonpy.material_response import LDMaterialProvider
        recipe=ld_recipe();recipe['materials']['metal-A']['parameters']['epsilon_infinity']=2.25
        model=self.build(self.resolve(recipe),'metal-A')
        self.assertEqual(model.config.epsilon_infinity,2.25)
        baseline=self.build(self.resolve(ld_recipe()),'metal-A')
        self.assertNotEqual(model.definition_digest,baseline.definition_digest)
        self.assertNotEqual(LDMaterialProvider(model.model).identity,LDMaterialProvider(baseline.model).identity)
        provider=LDMaterialProvider(model.model)
        self.assertEqual(provider.prepare(None).config.epsilon_infinity,2.25)

    def test_resolution_does_not_open_native(self):
        with (patch('salmonpy.ld._handles_api',side_effect=AssertionError('native call')),
              patch('salmonpy.sbe._handles_api',side_effect=AssertionError('native call'))):
            self.build(self.resolve(ld_recipe()),'metal-A')
            self.build(self.resolve(sbe_files(self.root)),'silicon')

    def test_ld_direct_agreement_identity_and_independence(self):
        from salmonpy.ld import LorentzDrude
        catalog=self.resolve(ld_recipe())
        model=self.build(catalog,'metal-A')
        direct=LorentzDrude(plasma_frequency=.5,strength=[1],resonance=[.3],damping=[.1])
        with model.new_state(dt=.02) as a, model.new_state(dt=.02) as b, direct.new_state(dt=.02) as reference:
            for field in ([.1,0,0],[.2,.01,0],[0,0,0]):
                before=b.current.copy()
                actual=a.step(field); expected=reference.step(field)
                np.testing.assert_array_equal(actual,expected)
                np.testing.assert_array_equal(b.current,before)
                np.testing.assert_array_equal(b.step(field),actual)
            self.assertEqual(a.material_id,'metal-A')
            self.assertEqual(a.definition_digest,catalog.definition_digests['metal-A'])
            with self.assertRaises(AttributeError): a.material_id='other'
        self.assertNotEqual(model.definition_digest,self.build(catalog,'metal-B').definition_digest)

    def test_sbe_direct_agreement_independence_and_snapshot(self):
        from salmonpy.sbe import SBEInitialState
        catalog=self.resolve(sbe_files(self.root))
        model=self.build(catalog,'silicon')
        direct=SBEInitialState.from_salmon(self.root,'si',bands=2,kpoints=1,electrons=2,lattice=np.eye(3)*4)
        with model.new_state(dt=.02) as a, model.new_state(dt=.02) as b, direct.new_state(dt=.02) as ref:
            snapshot=a.snapshot()
            for field in ([.01,0,0],[.02,0,0],[0,0,0]):
                before=b.density.copy()
                np.testing.assert_array_equal(a.step(field),ref.step(field))
                np.testing.assert_array_equal(b.density,before)
                np.testing.assert_array_equal(b.step(field),a.current)
            a.restore(snapshot)
            self.assertEqual(a.time,0)
            self.assertEqual(a.material_id,'silicon')
            self.assertEqual(b.time,.06)

    def test_sbe_source_change_rejected(self):
        catalog=self.resolve(sbe_files(self.root))
        p=self.root/'si_eigen.data';p.write_text(p.read_text().replace('-0.5','-0.6'))
        with self.assertRaisesRegex(InputError,'changed|match'):self.build(catalog,'silicon')

    def test_schema_units_and_paths_rejected(self):
        for field,value in [('units','eV'),('typo',1),('strength',[-1])]:
            data=ld_recipe();data['materials']['metal-A']['parameters'][field]=value
            with self.assertRaises(InputError):self.resolve(data)
        for field,value in [('directory','../'),('sysname','../si'),('units','eV'),('nonlocal_correction',1)]:
            data=sbe_files(self.root);data['materials']['silicon']['parameters'][field]=value
            with self.assertRaises(InputError):self.resolve(data)

    def test_unknown_material_and_custom_backend_not_dispatched(self):
        with self.assertRaises(InputError):self.build(self.resolve(ld_recipe()),'absent')
        self.registry.register('mine',version='1',models=['custom'],resolve_parameters=lambda *a,**k:{})
        catalog=self.resolve({'schema':'salmonpy.materials/1','materials':{'x':{'backend':'mine','model':'custom','parameters':{}}}})
        with self.assertRaises(InputError):self.build(catalog,'x')

    def test_relocated_sbe_bundle_keeps_digest_and_options(self):
        import shutil
        from salmonpy.preparation import prepare_inputs, verify_inputs
        from salmonpy.config import load_materials
        data=sbe_files(self.root)
        data['materials']['silicon']['parameters'].update(
            active_bands=1,correction_order=3,nonlocal_correction=True)
        files={name:name for name in ('si_k.data','si_eigen.data','si_tm.data')}
        bundle=prepare_inputs(data,self.root/'bundle',registry=self.registry,
                              base_directory=self.root,files=files)
        catalog=load_materials(bundle/'materials.json',registry=self.registry)
        moved=self.root/'moved';shutil.move(bundle,moved)
        for name in files:(self.root/name).unlink()
        verify_inputs(moved)
        model=self.api.material_from_catalog(catalog,'silicon',base_directory=moved)
        self.assertEqual(model.definition_digest,catalog.definition_digests['silicon'])
        self.assertEqual(model.eigen.shape,(1,1))
        self.assertEqual(model.correction_order,3)
        self.assertTrue(model.nonlocal_correction)
