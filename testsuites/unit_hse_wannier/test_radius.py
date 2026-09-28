import os
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]

class RadiusTest(unittest.TestCase):
    def test_periodic_support_hse_and_pbeh(self):
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'probe'
            fftw=Path(os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw'))
            blas=Path(os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas'))
            cmd=[os.environ.get('FC','gfortran'),'-O0','-g','-fcheck=all','-fopenmp','-I'+str(fftw/'include')]
            cmd += [str(ROOT/'src/xc'/f) for f in ['exx_local_fft.f90','lcfo_wf_support.f90','exx_wannier_gauge.f90','exx_wannier.f90']]
            cmd += [str(Path(__file__).with_name('radius_probe.f90')),'-L'+str(fftw/'lib'),'-lfftw3',
                    '-L'+str(blas/'lib'),'-lopenblas','-o',str(exe)]
            build=subprocess.run(cmd,cwd=tmp,text=True,capture_output=True)
            self.assertEqual(build.returncode,0,build.stderr)
            run=subprocess.run([exe],cwd=tmp,text=True,capture_output=True,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)

if __name__=='__main__':unittest.main()
