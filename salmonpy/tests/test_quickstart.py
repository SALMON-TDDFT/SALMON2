"""First calculation orchestration: no solver build or implicit downloads."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


class QuickstartTest(unittest.TestCase):
    def test_examples_success_failure_and_existing_output(self):
        from salmonpy.quickstart import run
        from salmonpy.errors import RuntimeFailure
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def completed(command, **kwargs):
                kwargs['stdout'].write('calculation finished\n')
                return subprocess.CompletedProcess(command, 0)
            with patch('salmonpy.quickstart.subprocess.run', side_effect=completed) as invoke:
                report = run(root/'success')
            self.assertEqual(report['status'], 'passed')
            self.assertEqual(len(invoke.call_args_list), 2)
            for call in invoke.call_args_list:
                command = call.args[0]
                self.assertIn(command[command.index('-m')+1], ('salmonpy.examples.ld_response', 'salmonpy.examples.sbe_response'))
                self.assertEqual(command[command.index('-X')+1], 'utf8')
                self.assertEqual(call.kwargs['env']['OMP_NUM_THREADS'], '1')
                self.assertEqual(call.kwargs['env']['OPENBLAS_NUM_THREADS'], '1')
            with self.assertRaises(FileExistsError):
                run(root/'success')
            self.assertTrue((root/'success/summary.json').is_file())
            with patch('salmonpy.quickstart.subprocess.run', return_value=subprocess.CompletedProcess([], 2)) as invoke:
                with self.assertRaises(RuntimeFailure): run(root/'failure')
                self.assertEqual(invoke.call_count, 1)
            self.assertEqual(json.loads((root/'failure/summary.json').read_text())['status'], 'failed')
            with patch('salmonpy.quickstart.subprocess.run', side_effect=subprocess.TimeoutExpired([], 60)):
                with self.assertRaises(RuntimeFailure): run(root/'timeout')
            self.assertEqual(json.loads((root/'timeout/summary.json').read_text())['status'], 'failed')

    def test_first_instructions_are_short_and_executable(self):
        from salmonpy.manual import load, search
        data = load()
        for lang in ('ja', 'en'):
            first = data['getting_started_i18n'][lang]
            self.assertEqual(first['test_command'], 'python -m salmonpy.quickstart')
            self.assertEqual(first['success_marker'], 'SALMONpy first calculations: PASS')
            self.assertEqual(len(first['capabilities']), 7)
            self.assertTrue(first['explanation_for_user'])
            self.assertTrue(all(c['input'] and c['output'] and c['preparation'] for c in first['capabilities']))
        self.assertEqual(search('first calculation')['tasks'][0]['id'], 'first_calculation')
        import sys
        result = subprocess.run([sys.executable, '-m', 'salmonpy.manual'], capture_output=True, text=True, encoding='utf-8', env=dict(os.environ, PYTHONIOENCODING='cp1252'), timeout=20, check=True)
        self.assertLess(len(result.stdout), 10000)
        self.assertEqual(json.loads(result.stdout)['test_command'], 'python -m salmonpy.quickstart')
