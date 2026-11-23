"""Bad campaign inputs must be rejected before spawning any MPI process."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class ScalingCLI(unittest.TestCase):
    def reject(self,extra):
        script=Path(__file__).resolve().parents[2]/'tools/test_mixed_scaling.py'
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'output'
            command=[sys.executable,str(script),'build','runtime','bundle',str(output),
                     '--rank-budget','18','--memory-budget-gib','48','--launcher',str(Path(tmp)/'absent-launcher'),*extra]
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,2,result.stderr)
            self.assertFalse(output.exists())

    def test_bad_memory_budget_before_mpi(self):
        for value in ('nan','inf','0','-1'):
            with self.subTest(value=value):self.reject(['--memory-budget-gib',value])

    def test_incompatible_or_odd_field_shape_before_mpi(self):
        for shape in (('80','80','80'),('95','96','96')):
            with self.subTest(shape=shape):self.reject(['--field-shape',*shape])

    def test_field_observation_campaign_rejects_bad_options_before_mpi(self):
        script=Path(__file__).resolve().parents[2]/'tools/test_diagonal_field_energy_mpi.py'
        for option,value,message in (('--observation-stride','0','observation-stride must be positive'),
                                     ('--grid-count','41','grid-count must be even and at least40'),
                                     ('--stop-after','0','stop-after must be in1..64')):
            with self.subTest(option=option),tempfile.TemporaryDirectory() as tmp:
                output=Path(tmp)/'output'
                result=subprocess.run([sys.executable,str(script),'build',str(output),option,value],
                                      capture_output=True,text=True)
                self.assertEqual(result.returncode,2,result.stderr)
                self.assertIn(message,result.stderr)
                self.assertFalse(output.exists())
