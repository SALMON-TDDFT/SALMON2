from pathlib import Path
import os
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[3]
class StreamTest(unittest.TestCase):
    def test_payload_streaming(self):
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'probe'
            subprocess.run([os.environ.get('FC','gfortran'),'-O0','-g','-fcheck=all','-ffree-line-length-none',str(ROOT/'src/gs/dc/lcfo_mesh_stream.f90'),str(Path(__file__).with_name('stream_probe.f90')),'-o',str(exe)],cwd=tmp,check=True)
            subprocess.run([str(exe)],check=True,timeout=20)
if __name__=='__main__':unittest.main()
