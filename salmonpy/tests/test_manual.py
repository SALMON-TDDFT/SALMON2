import unittest
import importlib
from salmonpy.manual import load, search, task_details

class ManualTest(unittest.TestCase):
    def test_user_facing_hse_names_and_entry_guidance(self):
        from importlib.resources import files
        import json
        d=load()
        for task in d['tasks']:
            fields={k:task.get(k) for k in ('intent','intent_i18n','inputs','outputs','keywords','constraints')}
            text=json.dumps(fields,ensure_ascii=False).casefold()
            for internal_label in ('高速pr','fast-pr','fast pr'):
                self.assertNotIn(internal_label,text,task['id'])
        for name in ('README.md','README.en.md','AI_USAGE.md','FAST_HSE_SYMMETRY.md','FAST_HSE_SYMMETRY.en.md'):
            text=files('salmonpy').joinpath(name).read_text(encoding='utf-8').casefold()
            for internal_label in ('高速pr','fast-pr','fast pr'):
                self.assertNotIn(internal_label,text,name)
        guide=files('salmonpy').joinpath('FAST_HSE_SYMMETRY.md').read_text(encoding='utf-8')
        self.assertIn('allow_symmetry=True',guide)
        self.assertIn('同じMPI分割',guide)
        self.assertIn('円偏光',guide)
        self.assertIn('docs/development/',guide)
        ja=files('salmonpy').joinpath('DFT.md').read_text(encoding='utf-8')
        self.assertNotIn('UPFのPAW、ultrasoft、SOCは拒否',ja)
        self.assertNotIn('擬ポテンシャルをダウンロードする処理はありません',ja)
        self.assertIn('spin_orbit=True',ja)
        self.assertIn('PSEUDOPOTENTIALS.md',ja)
        for name in ('README.md','README.en.md'):
            text=files('salmonpy').joinpath(name).read_text(encoding='utf-8')
            self.assertIn('mpi_mixed_native',text)
        self.assertIn('Configuration of this sample',files('salmonpy').joinpath('README.en.md').read_text(encoding='utf-8'))

        retained=next(t for t in d['tasks'] if t['id']=='retained_native_gs')
        constraints=' '.join(retained['constraints'])
        self.assertNotIn('Full SCF checkpoints, SOC/+U/DC GS',constraints)
        self.assertNotIn('legacy sources/default preparation',constraints)
        self.assertIn('save_scf',constraints)
        from salmonpy.native_gs import NativeGSState
        self.assertTrue(callable(NativeGSState.set_scf_spin_field))
        for id_ in ('checkpoint_compatibility','retained_native_gs','shared_standard_tddft','nickel_soc_smoke'):
            task=next(t for t in d['tasks'] if t['id']==id_)
            self.assertTrue(any('set_scf_spin_field([0., 0., 0.])' in c
                                and 'continue SCF at zero field' in c for c in task['constraints']))

        for name in ('GS_ARCHIVE.md','GS_ARCHIVE.en.md'):
            text=files('salmonpy').joinpath(name).read_text(encoding='utf-8')
            self.assertIn('SerialRestartConverter',text.split('```')[0])
        for name in ('CHECKPOINTING.md','CHECKPOINTING.en.md'):
            text=files('salmonpy').joinpath(name).read_text(encoding='utf-8')
            self.assertIn('execution_mode="checked"',text)
            self.assertIn('SOC',text)


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

    def test_search_routing_and_bounded_details(self):
        import json
        import subprocess
        import sys
        from salmonpy.manual import task_details
        for query in ('LDとSBEとTDDFTを混ぜて多層膜を計算',
                      'combine LD SBE TDDFT in a multilayer'):
            result=search(query)
            self.assertEqual(result['tasks'][0]['id'],'native_coupled')
            self.assertLessEqual(len(result['tasks']),5)
            self.assertLess(len(json.dumps(result,ensure_ascii=False)),100000)
            self.assertTrue(all('instructions' not in t for t in result['tasks']))
        self.assertIn('retained_native_gs',{t['id'] for t in search('retain',lang='en')['tasks']})
        self.assertEqual(search('task_details')['tasks'][0]['id'],'public_api')
        detail=task_details('checkpoint_compatibility',lang='en')
        self.assertEqual([i['heading'] for i in detail['instructions']],['en'])
        self.assertIn('CheckpointSchedule',detail['instructions'][0]['text'])
        with self.assertRaises(ValueError):task_details('unknown-task')
        with self.assertRaises(ValueError):search('test',limit=0)
        out=subprocess.check_output([sys.executable,'-m','salmonpy.manual','--lang','en','--task','checkpoint_compatibility'],text=True)
        self.assertEqual(json.loads(out)['id'],'checkpoint_compatibility')

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
        import json
        task=task_details(tasks[0]['id'])
        text=json.dumps(task,ensure_ascii=False)
        for obsolete in ('salmonpy.native-direct/3','solver-direct/1','solver-direct/2',
                         '--unprotected-direct-restore','`protected=False`'):
            self.assertNotIn(obsolete,text)
        self.assertIn('stop_after_seconds',task['inputs'])
        self.assertIn('True',task['outputs'])
        self.assertTrue(any('after_step' in example for example in task['examples']))
        for example in task['examples']:
            if 'CheckpointSchedule' in example:
                exec(compile(example,'timed-checkpoint-manual','exec'),{})
        self.assertIn('examples/mpi_mixed_native.py',task['sample_files'])
        for name in ('CHECKPOINTING.md','CHECKPOINTING.en.md'):
            guide=files('salmonpy').joinpath(name).read_text(encoding='utf-8')
            self.assertIn('salmonpy.native-direct/4',guide)
            self.assertIn('--stop-after-seconds',guide)
            self.assertIn('--direct-resume',guide)

        for name in ('CHECKPOINTING.md','CHECKPOINTING.en.md','examples/custom_checkpoint.py'):
            self.assertTrue(files('salmonpy').joinpath(name).is_file(),name)
