import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import resolve_materials
from salmonpy.errors import InputError
from salmonpy.fdtd import SalmonExecutable


@unittest.skipUnless(importlib.util.find_spec('ase') and importlib.util.find_spec('f90nml'), 'optional ASE/f90nml')
class InitialStateTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.initial_states'), 'saved GS selection missing')
        from salmonpy.initial_states import tddft_from_catalog
        self.load=tddft_from_catalog
        from ase import Atoms
        from ase.io import write
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        write(self.root/'cell.xyz',Atoms('Si2',positions=[[0,0,0],[1,2,3]],cell=[5,6,7],pbc=True),format='extxyz')
        (self.root/'Si_rps.dat').write_text('1 0.1 2 4\n')
        executable=self.root/'solver';executable.write_bytes(b'not a real native executable')
        self.engine=SalmonExecutable(executable,self.hash(executable))
        self.recipe={'schema':'salmonpy.materials/1','materials':{'film':{
            'backend':'salmonpy.dft','model':'dft','parameters':{
                'structure':{'path':'cell.xyz','format':'extxyz'},
                'pseudopotentials':[{'path':'Si_rps.dat','atomic_number':14,'local_channel':2}],
                'xc':'pz','states':4,'kmesh':[1,1,1],'grid':{'counts':[8,8,8]}}}}}
        model=material_from_catalog(self.catalog(),'film',base_directory=self.root).model
        self.native=self.root/'saved/native';self.native.mkdir(parents=True)
        (self.native/'data_for_restart').mkdir()
        for name in ('info.bin','occupation.bin','wfn.bin'):
            (self.native/'data_for_restart'/name).write_bytes(b'fixture only')
        (self.native/'gs.inp').write_text(model.input_text())
        (self.native/'gs_status.dat').write_text('SALMONPY_GS_V1\n1 2\n1e-10 1e-9 -1.0\nrho_dne\n')
        (self.native/'pseudo1_rps.dat').write_bytes(model.pseudopotentials[0].data)
        self.record={'schema':'salmonpy.gs/1','status':'complete','restart_saved':True,**model._restart_signature(),
                     'fractional_positions':model.fractional_positions.tolist(),
                     'energy_hartree':-1.,'iterations':2,'max_iterations':300,'tolerance':1e-9,
                     'executable_sha256':self.engine.sha256}
        self.recipe['initial_states']={'ground':{'material_id':'film','source':{
            'format':'salmonpy.gs/1','path':'saved','manifest_sha256':''}}}
        self.pin()

    @staticmethod
    def hash(path):return hashlib.sha256(path.read_bytes()).hexdigest()

    def pin(self):
        self.record['files']={p.relative_to(self.native.parent).as_posix():self.hash(p)
                              for p in self.native.rglob('*') if p.is_file()}
        path=self.native.parent/'ground-state.json'
        path.write_text(json.dumps(self.record))
        self.recipe['initial_states']['ground']['source']['manifest_sha256']=self.hash(path)

    def catalog(self):
        return resolve_materials(self.recipe,registry=builtin_material_registry(),base_directory=self.root)

    def initial(self):
        return self.load(self.catalog(),'ground',base_directory=self.root,executable=self.engine,max_steps=3)

    def test_malformed_manifest_metadata_is_input_error(self):
        original = dict(self.record)
        for name, value in [('files', []), ('files', {'native/gs.inp': 4}),
                            ('energy_hartree', None), ('iterations', True),
                            ('schema', 'unknown'), ('status', 'running')]:
            self.record = {**original, name: value}
            path = self.native.parent/'ground-state.json'
            path.write_text(json.dumps(self.record))
            self.recipe['initial_states']['ground']['source']['manifest_sha256'] = self.hash(path)
            with self.subTest(name=name, value=value), self.assertRaises(InputError):
                self.initial()

    def test_selection_records_identity_without_starting_worker(self):
        catalog=self.catalog()
        with patch('salmonpy.tddft.TDDFTState',side_effect=AssertionError('worker')):
            with self.initial() as initial:
                self.assertEqual(initial.material_id,'film')
                self.assertEqual(initial.initial_state_id,'ground')
                self.assertEqual(initial.initial_state_digest,self.recipe['initial_states']['ground']['source']['manifest_sha256'])
                self.assertEqual(initial.definition_digest,catalog.definition_digests['film'])
                self.assertEqual(initial.record()['initial_state_id'],'ground')
                raw=initial.model
        self.assertTrue(raw._owner.closed)

    def test_manifest_anchor_and_native_corruption_rejected(self):
        (self.native/'data_for_restart/wfn.bin').write_bytes(b'changed')
        with self.assertRaises(InputError):self.initial()
        self.pin()
        self.recipe['initial_states']['ground']['source']['manifest_sha256']='0'*64
        with self.assertRaisesRegex(InputError,'manifest'):self.initial()

    def test_scientific_mismatch_rejected(self):
        for name,value in [('grid',[9,8,8]),('kmesh',[2,1,1]),('xc','libxc_pz'),
                           ('fractional_positions',[[0,0,0],[0,0,0]])]:
            before=self.record[name];self.record[name]=value;self.pin()
            with self.subTest(name=name),self.assertRaisesRegex(InputError,'mismatch'):
                self.initial()
            self.record[name]=before
        self.pin()

    def test_actual_input_and_pseudo_checked_not_only_metadata(self):
        path=self.native/'gs.inp';original=path.read_text()
        path.write_text(original.replace('num_rgrid=8,8,8','num_rgrid=9,8,8'));self.pin()
        with self.assertRaisesRegex(InputError,'input'):self.initial()
        path.write_text(original)
        (self.native/'pseudo1_rps.dat').write_bytes(b'other potential');self.pin()
        with self.assertRaisesRegex(InputError,'pseudo'):self.initial()

    def test_source_schema_and_executable_identity(self):
        source=self.recipe['initial_states']['ground']['source']
        source['path']='../outside'
        with self.assertRaises(InputError):self.initial()
        source['path']='saved';source['unknown']=1
        with self.assertRaises(InputError):self.initial()
        source.pop('unknown')
        self.record['executable_sha256']='0'*64;self.pin()
        with self.assertRaisesRegex(InputError,'executable'):self.initial()

    def test_unknown_id_or_prepare_recipe_refused(self):
        with self.assertRaises(InputError):
            self.load(self.catalog(),'missing',base_directory=self.root,executable=self.engine)
        self.recipe['initial_states']['ground']={'material_id':'film','prepare':{}}
        with self.assertRaises(InputError):self.initial()

    def test_failed_snapshot_check_closes_candidate(self):
        from types import SimpleNamespace
        from unittest.mock import Mock
        candidate=SimpleNamespace(manifest=(),close=Mock())
        with patch('salmonpy.dft.GroundStateResult.tddft_initial_state',return_value=candidate):
            with self.assertRaisesRegex(InputError,'snapshot'):self.initial()
        candidate.close.assert_called_once_with()

    def test_initial_identity_propagates_to_state_wrapper(self):
        from types import SimpleNamespace
        from salmonpy.tddft import TDDFTInitialState
        with self.initial() as initial:
            with patch.object(TDDFTInitialState,'new_state',return_value=SimpleNamespace()) as create:
                state=initial.new_state(dt=.02)
            create.assert_called_once_with(dt=.02)
            self.assertEqual(state.initial_state_id,'ground')
            self.assertEqual(state.initial_state_digest,initial.initial_state_digest)
            self.assertEqual(state.material_id,initial.material_id)
