import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from salmonpy.doctor import diagnose


class DoctorTests(unittest.TestCase):
    def test_no_execution_and_no_implied_budget(self):
        with patch('subprocess.run', side_effect=AssertionError('must not launch')):
            result = diagnose()
        self.assertEqual(result['external_salmon']['status'], 'not_configured')
        self.assertEqual(result['resources']['memory_budget'], 'not_specified')
        self.assertEqual(result['resources']['maximum_problem_size'], 'unknown')
        self.assertIn('next_steps', result)

    def test_identity_is_not_capability_certification(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'salmon'
            path.write_bytes(b'not an executable')
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            result = diagnose(executable=path, sha256=digest)
            self.assertEqual(result['external_salmon']['status'], 'hash_matches')
            self.assertEqual(result['external_salmon']['capabilities'], 'unverified')
            result = diagnose(executable=path, sha256='0'*64)
            self.assertEqual(result['external_salmon']['status'], 'hash_mismatch')
            self.assertEqual(diagnose(executable=path)['external_salmon']['status'], 'hash_required')
            self.assertEqual(diagnose(executable=Path(tmp)/'absent')['external_salmon']['status'], 'missing')

    def test_native_failure_is_reported(self):
        with patch('salmonpy.ld._native', side_effect=OSError('missing runtime')):
            result = diagnose()
        self.assertEqual(result['native']['ld']['status'], 'unavailable')
        self.assertIn('missing runtime', result['native']['ld']['detail'])

    def test_english_diagnosis_preserves_status_codes(self):
        report = diagnose(lang='en')
        self.assertEqual(report['external_salmon']['status'], 'not_configured')
        self.assertTrue(all(s.isascii() for s in report['next_steps']))
        with self.assertRaises(ValueError): diagnose(lang='xx')

    def test_english_cli_help(self):
        import subprocess, sys
        result=subprocess.run([sys.executable,'-m','salmonpy.doctor','--lang','en','--help'],
                              capture_output=True,text=True,check=True)
        self.assertIn('External SALMON',result.stdout)
        self.assertTrue(result.stdout.isascii())
