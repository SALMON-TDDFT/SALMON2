"""Contracts for the configuration-only material catalog."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from salmonpy.errors import InputError


def recipe():
    return {"schema": "salmonpy.materials/1", "materials": {
        "Si": {"model": "lab:test", "backend": "local", "parameters": {"mesh": [2, 3, 4]}},
        "Si-slab": {"model": "lab:test", "backend": "local", "parameters": {"mesh": [4, 4, 4]}}
    }, "initial_states": {"ground": {"material_id": "Si", "source": {"path": "gs/manifest.json"}}}}


class ConfigTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.config'),
                             'material catalog implementation is missing')
        from salmonpy import config
        self.config = config
        self.seen = []
        def resolve(model, parameters, *, base_directory):
            self.seen.append((model, base_directory))
            if set(parameters) != {'mesh'}:
                raise InputError('expected mesh only')
            return {'mesh': parameters['mesh'], 'units': 'dimensionless'}
        self.registry = config.MaterialRegistry()
        self.registry.register('local', version='test-1', models=['lab:test'],
                               resolve_parameters=resolve)

    def resolve(self, data=None):
        return self.config.resolve_materials(recipe() if data is None else data,
                                             registry=self.registry, base_directory=Path.cwd())

    def test_ids_and_parameters_are_independent_and_frozen(self):
        data = recipe()
        catalog = self.resolve(data)
        data['materials']['Si']['parameters']['mesh'][0] = 99
        self.assertEqual(catalog.materials['Si']['parameters']['mesh'], (2, 3, 4))
        self.assertEqual(catalog.materials['Si-slab']['parameters']['mesh'], (4, 4, 4))
        with self.assertRaises(TypeError):
            catalog.materials['Si']['model'] = 'changed'
        with self.assertRaises(AttributeError):
            catalog.materials = {}
        exported = catalog.as_dict()
        exported['materials']['Si']['parameters']['mesh'][0] = 0
        self.assertEqual(catalog.materials['Si']['parameters']['mesh'][0], 2)
        self.assertEqual(catalog.initial_states['ground']['material_id'], 'Si')
        self.assertEqual(catalog.materials['Si']['backend_version'], 'test-1')

    def test_json_and_python_share_resolution_and_file_base(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'materials.json'
            path.write_text(json.dumps(recipe()))
            loaded = self.config.load_materials(path, registry=self.registry)
            direct = self.config.resolve_materials(recipe(), registry=self.registry,
                                                    base_directory=path.parent)
            self.assertEqual(loaded.as_dict(), direct.as_dict())
            self.assertEqual(self.seen[-1][1], path.parent.resolve())

    def test_digest_is_order_independent_and_tracks_parameters(self):
        a = self.resolve()
        data = recipe()
        data['materials'] = dict(reversed(list(data['materials'].items())))
        self.assertEqual(a.definition_digests, self.resolve(data).definition_digests)
        data['materials']['Si']['parameters']['mesh'][0] = 9
        b = self.resolve(data)
        self.assertNotEqual(a.definition_digests['Si'], b.definition_digests['Si'])
        self.assertEqual(a.definition_digests['Si-slab'], b.definition_digests['Si-slab'])

    def test_bad_references_and_unknown_fields_are_rejected_before_callback(self):
        for change in ('reference', 'backend', 'model', 'placement', 'typo'):
            with self.subTest(change=change):
                data = recipe()
                if change == 'reference': data['initial_states']['ground']['material_id'] = 'missing'
                elif change in ('backend', 'model'): data['materials']['Si'][change] = 'missing'
                elif change == 'placement': data['placement'] = {}
                else: data['materials']['Si']['paramters'] = {}
                self.seen.clear()
                with self.assertRaises(InputError): self.resolve(data)
                self.assertEqual(self.seen, [])

    def test_initial_state_has_exactly_one_recipe_or_source(self):
        for entry in ({'material_id': 'Si'}, {'material_id': 'Si', 'source': {}, 'prepare': {}},
                      {'material_id': 'Si', 'source': 'file'}):
            data = recipe()
            data['initial_states']['ground'] = entry
            with self.assertRaises(InputError): self.resolve(data)

    def test_identifier_is_exact_and_not_sanitized(self):
        for identifier in ('', ' Si', 'Si ', 'Si\n', 'Si\x00'):
            data = recipe()
            data['materials'][identifier] = data['materials'].pop('Si')
            with self.assertRaises(InputError): self.resolve(data)
        data = recipe()
        data['materials']['結晶/A'] = data['materials'].pop('Si')
        data['initial_states']['ground']['material_id'] = '結晶/A'
        self.assertIn('結晶/A', self.resolve(data).materials)

    def test_non_json_values_and_cycles_are_rejected(self):
        cyclic = []; cyclic.append(cyclic)
        for value in (float('nan'), float('inf'), object(), {1: 'bad'}, cyclic, '\ud800'):
            data = recipe()
            data['materials']['Si']['parameters']['mesh'] = value
            with self.assertRaises(InputError): self.resolve(data)

    def test_duplicate_keys_in_json_are_not_silently_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'bad.json'
            for content in ('{"schema": "x", "schema": "y"}', '{"x":NaN}', '{invalid'):
                path.write_text(content)
                with self.assertRaises(InputError): self.config.load_materials(path, registry=self.registry)

    def test_registry_requires_explicit_unique_callable_registration(self):
        with self.assertRaises(InputError):
            self.registry.register('local', version='2', models=['other'], resolve_parameters=lambda *a: {})
        with self.assertRaises(InputError):
            self.registry.register('bad', version='1', models=['x'], resolve_parameters=None)

    def test_callback_validation_reports_material_path(self):
        data = recipe()
        data['materials']['Si']['parameters']['typo'] = 1
        with self.assertRaisesRegex(InputError, 'materials.*Si.*parameters'):
            self.resolve(data)

    def test_callback_output_is_validated_and_copied(self):
        returned = {'values': [1]}
        registry = self.config.MaterialRegistry()
        registry.register('local', version='1', models=['lab:test'],
                          resolve_parameters=lambda *a, **kw: returned)
        result = self.config.resolve_materials(recipe(), registry=registry, base_directory=Path.cwd())
        returned['values'][0] = 9
        self.assertEqual(result.materials['Si']['parameters']['values'], (1,))
        returned['values'] = float('nan')
        with self.assertRaises(InputError):
            self.config.resolve_materials(recipe(), registry=registry, base_directory=Path.cwd())

    def test_unserializable_integer_is_an_input_error(self):
        import sys
        if not getattr(sys, 'get_int_max_str_digits', lambda: 0)():
            self.skipTest('interpreter has no integer string conversion limit')
        data = recipe()
        data['materials']['Si']['parameters']['mesh'] = 10 ** (sys.get_int_max_str_digits() + 1)
        with self.assertRaisesRegex(InputError, 'materials.*Si'):
            self.resolve(data)
