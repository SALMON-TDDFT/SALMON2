"""Frozen-state finite differences of asymmetric core/full projector products."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]

@unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
class DCProjectorTest(unittest.TestCase):
    def test_truncated_complex_projectors(self):
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'probe'
            subprocess.run(['gfortran','-O2','-fcheck=all','-ffree-line-length-none',
                str(ROOT/'src/gs/dc/dc_projector_force.f90'),
                str(Path(__file__).with_name('dc_projector_probe.f90')),'-o',str(exe)],cwd=tmp,check=True,capture_output=True)
            subprocess.run([str(exe)],cwd=tmp,check=True,capture_output=True)
