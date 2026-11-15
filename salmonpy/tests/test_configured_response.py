import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
from salmonpy.errors import InputError
from salmonpy.builtin_materials import ConfiguredMaterial, _make, builtin_material_registry, material_from_catalog
from salmonpy.config import resolve_materials
if __package__:
    from .test_builtin_materials import ld_recipe
else:
    from test_builtin_materials import ld_recipe


class ConfiguredResponseTests(unittest.TestCase):
    def setUp(self):
        from salmonpy.materials import bind_configured_ld_endpoints, bind_configured_tddft
        self.ld_bind = bind_configured_ld_endpoints
        self.tddft_bind = bind_configured_tddft

    def ld_model(self):
        catalog=resolve_materials(ld_recipe(),registry=builtin_material_registry(),base_directory=Path.cwd())
        return material_from_catalog(catalog,'metal-A',base_directory=Path.cwd())

    def test_ld_binding_no_state_until_create_and_exact_conversion(self):
        model=self.ld_model()
        with patch.object(type(model.model),'new_state',side_effect=AssertionError('state')):
            binding=self.ld_bind(model,key='ld')
            prepared=binding.provider.prepare(None)
        self.assertEqual(dict(binding.material_record),model.record())
        a=binding.backend.create_state(prepared,10,.02,None)
        b=binding.backend.create_state(prepared,20,.02,None)
        try:
            with model.new_state(dt=.02) as reference:
                old=np.zeros(3)
                for n in range(1,5):
                    new=np.array([n*.0001,0,0])
                    expected=reference.step(-(new-old)/.02)[1]
                    np.testing.assert_array_equal(a.advance(old,new),expected)
                    np.testing.assert_array_equal(b.advance(old,new),expected)
                    old=new
        finally:
            a.close();b.close()

    def tddft_model(self):
        from salmonpy.tddft import TDDFTInitialState
        raw=object.__new__(TDDFTInitialState)
        for name,value in {'executable':SimpleNamespace(sha256='a'*64),'manifest':(('rt.inp','b'*64),),
                           '_input':'input','max_steps':20,'threads':1,'timeout':300.}.items():
            object.__setattr__(raw,name,value)
        return _make(ConfiguredMaterial,model=raw,material_id='Si',definition_digest='c'*64,
                     initial_state_id='GS',initial_state_digest='d'*64)

    def test_tddft_runtime_verified_and_snapshot_borrowed(self):
        model=self.tddft_model()
        with patch('salmonpy.runtime.load_runtime',return_value=model.model.executable) as load:
            binding=self.tddft_bind(model,key='tddft',runtime='/tmp/fake-runtime')
            prepared=binding.provider.prepare(None)
            self.assertIs(prepared.initial,model.model)
            self.assertEqual(dict(binding.material_record),model.record())
            self.assertTrue(load.called)
        with patch('salmonpy.runtime.load_runtime',return_value=SimpleNamespace(sha256='f'*64)):
            with self.assertRaises(InputError): self.tddft_bind(model,key='bad',runtime='/tmp/fake-runtime')

    def test_wrong_types_rejected(self):
        with self.assertRaises(InputError): self.ld_bind(self.ld_model().model,key='bad')
        with self.assertRaises(InputError): self.tddft_bind(self.ld_model(),key='bad',runtime='/tmp/unused')
