import unittest
from salmonpy.errors import InputError
from salmonpy.examples.custom_material import ExampleMaterial


class Recorded:
    backend = 'recorded-example/1'
    def __init__(self): self.material_id = 'film'
    def record(self): return {'material_id': self.material_id, 'definition_digest': 'a'*64}
    def new_state(self, **kwargs): return ExampleMaterial().new_state(**kwargs)


class RecordedTests(unittest.TestCase):
    def setUp(self):
        from salmonpy.materials import bind_recorded_material, response_identity
        self.bind = bind_recorded_material
        self.record = response_identity

    def test_binding_and_ordered_record(self):
        model = Recorded()
        binding = self.bind(model, key='cache')
        record = self.record({90:binding,10:binding}, {10:'right',90:'left'}, run_id='run')
        self.assertEqual([e['point_id'] for e in record['points']], ['right','left'])
        self.assertEqual(binding.material_record['material_id'], 'film')
        with self.assertRaises(TypeError): binding.material_record['material_id']='other'
        record['points'][0]['material_id']='other'
        self.assertEqual(binding.material_record['material_id'], 'film')
        model.material_id='changed'
        with self.assertRaises(InputError): binding.provider.prepare(None)

    def test_invalid_mapping_and_unrecorded_binding(self):
        from salmonpy.materials import bind_material
        binding = self.bind(Recorded(), key='cache')
        for mapping in ({}, {10:'x',20:'x'}, {10:' x'}):
            with self.assertRaises(InputError): self.record({10:binding}, mapping, run_id='run')
        old = bind_material(ExampleMaterial(), key='old', identity='old')
        with self.assertRaises(InputError): self.record({10:old}, {10:'x'}, run_id='r')

    def test_wrapped_tddft_does_not_bypass_managed_runtime(self):
        from salmonpy.builtin_materials import ConfiguredMaterial, _make
        from salmonpy.tddft import TDDFTInitialState
        from salmonpy.materials import bind_material
        raw = object.__new__(TDDFTInitialState)
        wrapped = _make(ConfiguredMaterial, model=raw, material_id='Si', definition_digest='a'*64)
        with self.assertRaisesRegex(InputError, 'TDDFTMaterialProvider'):
            bind_material(wrapped, key='test', identity='test')

    def test_change_after_prepare_rejected(self):
        model = Recorded()
        binding = self.bind(model, key='cache')
        prepared = binding.provider.prepare(None)
        model.material_id = 'changed'
        with self.assertRaises(InputError): binding.backend.create_state(prepared, 10, .02, None)

    def test_configured_sbe_matches_direct(self):
        import tempfile
        from pathlib import Path
        import numpy as np
        if __package__:
            from .test_builtin_materials import sbe_files
        else:
            from test_builtin_materials import sbe_files
        from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
        from salmonpy.config import resolve_materials
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            catalog = resolve_materials(sbe_files(root), registry=builtin_material_registry(), base_directory=root)
            model = material_from_catalog(catalog, 'silicon', base_directory=root)
            binding = self.bind(model, key='sbe')
            adapted = binding.backend.create_state(binding.provider.prepare(None), 10, .02, None)
            try:
                with model.new_state(dt=.02) as direct:
                    for n in range(1, 5):
                        A = [n*.0001, 0, 0]
                        np.testing.assert_array_equal(adapted.step(A), direct.step(A))
            finally:
                adapted.close()

    def test_direct_binding_copies_record(self):
        from salmonpy.material_response import MaterialBinding
        data = Recorded().record()
        binding = MaterialBinding('key', object(), object(), data)
        data['material_id'] = 'changed'
        self.assertEqual(binding.material_record['material_id'], 'film')
        with self.assertRaises(TypeError): binding.material_record['material_id'] = 'bad'
