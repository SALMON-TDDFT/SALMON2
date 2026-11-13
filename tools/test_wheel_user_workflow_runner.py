"""Regression: measurement failure must not change the calculation outcome."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from test_wheel_user_workflow import run_case


class MeasurementSeparationTests(unittest.TestCase):
    def test_failed_clock_preserves_success_and_failure(self):
        def unavailable():
            raise OSError('measurement unavailable')
        with tempfile.TemporaryDirectory() as tmp:
            for code in (0, 3):
                r = run_case([sys.executable, '-c', f'raise SystemExit({code})'],
                             Path(tmp), f'case{code}', dict(os.environ), clock=unavailable)
                self.assertEqual(r['exit_code'], code)
                self.assertEqual(r['status'], 'passed' if code == 0 else 'failed')
                self.assertEqual(r['timing']['status'], 'unavailable')
                self.assertIsNone(r['timing']['wall_seconds'])
                self.assertIsNone(r['peak_memory']['bytes'])


if __name__ == '__main__':
    unittest.main()
