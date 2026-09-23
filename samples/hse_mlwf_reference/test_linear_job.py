import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import run_linear_comparison as job


class LinearJobTests(unittest.TestCase):
 def test_analysis_retry_does_not_repeat_completed_propagation(self):
  with tempfile.TemporaryDirectory() as directory:
   hse=Path(directory)
   def propagation(*a,**kw):
    (hse/'status.json').write_text(json.dumps(dict(status='completed',accepted_step=2750,dt_au=.32,amplitude=1e-4)))
   with patch.object(job,'run',side_effect=propagation),patch.object(job,'analyze',side_effect=RuntimeError('analysis failure')):
    with self.assertRaisesRegex(RuntimeError,'analysis failure'):job.execute('export','state',hse,'td','out')
   self.assertEqual(json.loads((hse/'job_status.json').read_text())['status'],'failed')
   with patch.object(job,'fingerprint',create=True,return_value='fixture'),patch.object(job,'load_checkpoint',create=True),patch.object(job,'run',side_effect=AssertionError('must not propagate again')),patch.object(job,'analyze') as analysis:
    job.execute('export','state',hse,'td','out');analysis.assert_called_once()
   self.assertEqual(json.loads((hse/'job_status.json').read_text())['status'],'completed')

 def test_analysis_retry_rejects_wrong_initial_state(self):
  with tempfile.TemporaryDirectory() as directory:
   hse=Path(directory);(hse/'status.json').write_text(json.dumps(dict(status='completed',accepted_step=2750,dt_au=.32,amplitude=1e-4)))
   with patch.object(job,'fingerprint',create=True,return_value='changed'),patch.object(job,'load_checkpoint',create=True,side_effect=ValueError('Restart mismatch: initial_hash')),patch.object(job,'analyze') as analysis:
    with self.assertRaisesRegex(ValueError,'initial_hash'):job.execute('export','state',hse,'td','out')
    analysis.assert_not_called()

if __name__=='__main__':unittest.main()
