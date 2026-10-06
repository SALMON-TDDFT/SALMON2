import os
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[3]

class LocalFFTTest(unittest.TestCase):
    def test_direct_periodic_convolution(self):
        fftw=Path(os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw'))
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'probe'
            build=subprocess.run([os.environ.get('FC','gfortran'),'-O0','-g','-fcheck=all',
                '-I'+str(fftw/'include'),str(ROOT/'src/xc/exx_local_fft.f90'),
                str(Path(__file__).with_name('local_fft_probe.f90')),
                '-L'+str(fftw/'lib'),'-lfftw3','-o',str(exe)],cwd=tmp,text=True,capture_output=True)
            self.assertEqual(build.returncode,0,build.stderr)
            run=subprocess.run([exe],cwd=tmp,text=True,capture_output=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)

if __name__=='__main__':unittest.main()
