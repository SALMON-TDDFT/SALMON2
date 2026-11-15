import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import resolve_materials
from salmonpy.errors import InputError


@unittest.skipUnless(importlib.util.find_spec('ase'), 'optional ASE')
class DFTMaterialConfigTests(unittest.TestCase):
    def setUp(self):
        from ase import Atoms
        from ase.io import write
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.atoms=Atoms('Si2',positions=[[0,0,0],[1,2,3]],cell=[5,6,7],pbc=True)
        write(self.root/'cell.xyz',self.atoms,format='extxyz')
        (self.root/'Si_rps.dat').write_text('1 0.1 2 4\n')
        self.recipe={'schema':'salmonpy.materials/1','materials':{'si-film':{
            'backend':'salmonpy.dft','model':'dft','parameters':{
                'structure':{'path':'cell.xyz','format':'extxyz'},
                'pseudopotentials':[{'path':'Si_rps.dat','atomic_number':14,'local_channel':2}],
                'xc':'pz','states':4,'kmesh':[1,1,1],
                'grid':{'spacing':{'value':.4,'unit':'angstrom'}}}}}}

    def resolve(self):
        return resolve_materials(self.recipe,registry=builtin_material_registry(),base_directory=self.root)

    def build(self,catalog):
        return material_from_catalog(catalog,'si-film',base_directory=self.root)

    def test_matches_original_dft_and_keeps_grid_provenance(self):
        from salmonpy.ase import dft_from_atoms
        from salmonpy.dft import Pseudopotential
        with patch('salmonpy.dft.DFT.run',side_effect=AssertionError('unexpected execution')):
            catalog=self.resolve();configured=self.build(catalog)
        direct=dft_from_atoms(self.atoms,pseudopotentials=[Pseudopotential.from_file(
            self.root/'Si_rps.dat',atomic_number=14,local_channel=2)],states=4,kmesh=(1,1,1),grid=(13,15,18))
        self.assertEqual(configured.model.input_text(),direct.input_text())
        p=catalog.materials['si-film']['parameters']
        self.assertEqual(p['electrons'],8)
        self.assertEqual(p['grid_resolution']['counts'],(13,15,18))
        self.assertEqual(set(p['file_sha256']),{'cell.xyz','Si_rps.dat'})
        self.assertEqual(configured.material_id,'si-film')

    def test_count_grid_and_changed_files(self):
        self.recipe['materials']['si-film']['parameters']['grid']={'counts':[16,16,16]}
        catalog=self.resolve();self.assertEqual(self.build(catalog).grid,(16,16,16))
        p=self.root/'Si_rps.dat';p.write_text(p.read_text()+'# changed\n')
        with self.assertRaisesRegex(InputError,'changed|match'):self.build(catalog)

    def test_paths_units_conflicting_grid_and_unknown_fields(self):
        import copy
        for key,value in [('grid',{'counts':[8,8,8],'spacing':{'value':.4,'unit':'angstrom'}}),
                          ('grid',{'spacing':{'value':.4,'unit':'unknown'}}),
                          ('xc','unknown'),('electrons',6),('misspelled',1),
                          ('structure',{'path':'../cell.xyz','format':'extxyz'})]:
            with self.subTest(key=key):
                original=copy.deepcopy(self.recipe)
                self.recipe['materials']['si-film']['parameters'][key]=value
                with self.assertRaises(InputError):self.resolve()
                self.recipe=original

    def test_gs_tddft_identity_delegation_without_claiming_native_execution(self):
        configured=self.build(self.resolve())
        class Initial:
            closed=False
            def close(self): self.closed=True
        initial=Initial()
        class GS:
            def tddft_initial_state(self,**kwargs):
                self.arguments=kwargs
                return initial
        raw=GS()
        with patch('salmonpy.dft.DFT.run',return_value=raw) as run:
            gs=configured.run('engine','new-dir',threads=1)
        run.assert_called_once_with('engine','new-dir',threads=1)
        self.assertEqual(gs.material_id,'si-film')
        self.assertIs(gs.result,raw)
        with gs.tddft_initial_state(executable='engine',max_steps=4) as rt:
            self.assertEqual(rt.definition_digest,configured.definition_digest)
            self.assertEqual(rt.material_id,'si-film')
            self.assertIs(rt.model,initial)
        self.assertTrue(initial.closed)

    def test_wrapped_warm_start_unwraps_before_existing_checks(self):
        configured=self.build(self.resolve())
        raw=object()
        with patch('salmonpy.dft.DFT.run',return_value=raw): previous=configured.run('engine','old')
        with patch('salmonpy.dft.DFT.run',return_value=object()) as run:
            configured.run('engine','new',initial_guess=previous)
        self.assertIs(run.call_args.kwargs['initial_guess'],raw)

    def test_upf_element_mismatch_rejected(self):
        path=self.root/'wrong.upf'
        path.write_text('<UPF><PP_HEADER element="Ge" pseudo_type="NC" z_valence="4"/></UPF>')
        self.recipe['materials']['si-film']['parameters']['pseudopotentials'][0]['path']='wrong.upf'
        with self.assertRaisesRegex(InputError,'element'): self.resolve()
