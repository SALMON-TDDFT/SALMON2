from pathlib import Path
import os
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[3]
class TileTest(unittest.TestCase):
    def test_complex_tiles(self):
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'probe'
            subprocess.run([os.environ.get('FC','gfortran'),'-O0','-g','-fcheck=all','-ffree-line-length-none',str(ROOT/'src/gs/dc/lcfo_mesh_tile.f90'),str(Path(__file__).with_name('tile_probe.f90')),'-o',str(exe)],cwd=tmp,check=True)
            subprocess.run([str(exe)],check=True,timeout=10)
if __name__=='__main__':unittest.main()
