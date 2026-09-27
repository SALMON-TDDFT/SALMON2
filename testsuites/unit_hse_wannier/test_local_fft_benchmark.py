"""Bounded operator/FFT-volume benchmark; timings are reported, not asserted."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]

class LocalFFTBenchmark(unittest.TestCase):
    def test_volume_and_action(self):
        fftw=Path(os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw'))
        blas=Path(os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas'))
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'probe'
            cmd=[os.environ.get('FC','gfortran'),'-O2','-g','-fcheck=all','-fopenmp',
                 '-fno-tree-loop-vectorize','-I'+str(fftw/'include')]
            cmd += [str(ROOT/'src/xc'/f) for f in ['exx_local_fft.f90','lcfo_wf_support.f90',
                    'hse_wannier_gauge.f90','hse_wannier.f90']]
            cmd += [str(Path(__file__).with_name('local_fft_benchmark.f90')),
                    '-L'+str(fftw/'lib'),'-lfftw3','-L'+str(blas/'lib'),'-lopenblas','-o',str(exe)]
            build=subprocess.run(cmd,cwd=tmp,text=True,capture_output=True)
            self.assertEqual(build.returncode,0,build.stderr)
            run=subprocess.run([exe],cwd=tmp,text=True,capture_output=True,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            print(run.stdout,end='')

if __name__=='__main__':unittest.main()
