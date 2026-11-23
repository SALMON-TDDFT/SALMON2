import unittest
import importlib
from salmonpy.manual import load, search

class ManualTest(unittest.TestCase):
    def test_material_guides_are_packaged(self):
        from importlib.resources import files
        for name in ('FDTD_MATERIALS.md', 'FDTD_MATERIALS.en.md', 'DISTRIBUTED.md', 'DISTRIBUTED.en.md',
                     'FIELDS.md', 'FIELDS.en.md'):
            self.assertTrue(files('salmonpy').joinpath(name).is_file(), name)

    def test_tasks_and_examples(self):
        d=load()
        self.assertEqual(len({t['id'] for t in d['tasks']}),len(d['tasks']))
        for task in d['tasks']:
            self.assertTrue(task['instructions'])
            self.assertIn(task['external_input_source'],d['sources'])
            available=all(importlib.util.find_spec(name) is not None for name in task.get('environment',{}).get('import_requires',[]))
            for api in (task['api'] if available else []):
                parts=api.split('.')
                obj=importlib.import_module('.'.join(parts[:2]))
                for name in parts[2:]: obj=getattr(obj,name)
            for code in task['examples']: compile(code,task['id'],'exec')
    def test_tddft_current_sign_is_charge_current(self):
        task = next(t for t in load()['tasks'] if t['id'] == 'tddft')
        self.assertIn('physical charge-current density j', task['outputs'])
        self.assertIn('J = -j', task['outputs'])
        self.assertIn('do not negate again', task['outputs'])

    def test_search(self):
        for query,expected in [('バンド図を描きたい','bands'),('TDDFT','tddft'),('ドルーデ','ld'),('runtime','runtime')]:
            self.assertEqual(search(query)['tasks'][0]['id'],expected)
        self.assertEqual(search('unrecognizedxyz')['tasks'],[])
        with self.assertRaises(ValueError):search('')

    def test_portable_guidance(self):
        d=load()
        self.assertIn('environment',search('バンド'))
        import json
        for task in d['tasks']:
            text=json.dumps(task,ensure_ascii=False)
            for host_path in ('/Users/', '/tmp/', '/opt/homebrew/', 'canonical calc2026'):
                self.assertNotIn(host_path,text)
            self.assertIn('environment',task)

    def test_beginner_resource_intake_is_always_returned(self):
        for query in ('DFT','LD','初めて','unrecognizedxyz'):
            intake=search(query)['beginner_intake']
            self.assertEqual([q['id'] for q in intake['questions']],['memory_budget','parallel_budget','intent'])
            self.assertTrue(intake['reuse_known_answers'])
            self.assertEqual(intake['implementation_status']['hard_memory_enforcement'],'not_implemented')
            self.assertIn('evidence',intake['capacity_report_fields'])

    def test_first_time_and_unknown_guidance(self):
        for query in ('初めて', '何ができる', 'unrecognizedxyz'):
            result = search(query)
            self.assertTrue(result['getting_started'])
            self.assertIn('dft', {t['id'] for t in result['task_catalog']})
        ids = {t['id'] for t in search('レーザーによる電子の動き')['tasks']}
        self.assertTrue({'sbe', 'tddft'}.issubset(ids))

    def test_english_intents_and_language(self):
        self.assertEqual(search('I want to run TDDFT')['tasks'][0]['id'], 'tddft')
        self.assertEqual(search('What can I do?')['tasks'], [])
        self.assertTrue({'sbe','tddft'} <= {t['id'] for t in search('laser induced current')['tasks']})
        self.assertEqual(search('初めて', lang='en')['language'], 'en')
        self.assertTrue(search('初めて', lang='en')['getting_started']['first_action'].isascii())

    def test_ase_optical_task_and_guides(self):
        from importlib.resources import files
        d = load()
        task = next(t for t in d['tasks'] if t['id'] == 'ase_tddft')
        self.assertIn('salmonpy.ase.tddft_from_atoms', task['api'])
        self.assertEqual(task['guides']['en'], 'ASE.en.md')
        self.assertIn('ase_tddft', {t['id'] for t in search('Atoms current')['tasks']})
        self.assertTrue(files('salmonpy').joinpath('examples/ase_tddft.py').is_file())
        for name in ('ASE.md', 'ASE.en.md', 'DFT.md', 'ASE_CALCULATOR.md', 'AI_USAGE.md'):
            self.assertEqual(d['documents'][name]['sections'][0]['text'],
                             files('salmonpy').joinpath(name).read_text(encoding='utf-8'))

    def test_all_embedded_guides_match_packaged_text(self):
        from importlib.resources import files
        d=load()
        for name,document in d['documents'].items():
            self.assertEqual(document['sections'][0]['text'],
                             files('salmonpy').joinpath(name).read_text(encoding='utf-8'),name)
        for query,id_ in [('snapshot restore','snapshots'),('公開API','public_api')]:
            self.assertIn(id_,{t['id'] for t in search(query)['tasks']})

    def test_common_checkpoint_api_is_searchable_and_packaged(self):
        from importlib.resources import files
        tasks=search('checkpoint compatibility')['tasks']
        self.assertEqual(tasks[0]['id'],'checkpoint_compatibility')
        for name in ('CHECKPOINTING.md','CHECKPOINTING.en.md','examples/custom_checkpoint.py'):
            self.assertTrue(files('salmonpy').joinpath(name).is_file(),name)
