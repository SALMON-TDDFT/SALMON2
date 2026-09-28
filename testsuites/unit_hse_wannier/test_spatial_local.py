"""MPI compact convolution parity with the same distributed global multiplier."""
import argparse, os, subprocess, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);p.add_argument('--ranks',type=int,nargs='+',default=[1,2,4]);a=p.parse_args()
b=a.build.resolve();obj=b/'src/CMakeFiles/salmon.dir'
with tempfile.TemporaryDirectory() as tmp:
    exe=Path(tmp)/'probe'
    subprocess.run([os.environ.get('MPIFC','mpifort'),'-O0','-g','-fopenmp','-fcheck=all','-ffree-line-length-none',
        '-I'+str(b),'-I'+os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw')+'/include',
        str(ROOT/'src/xc/exx_local_fft.f90'),str(ROOT/'src/xc/exx_spatial_local.f90'),
        str(Path(__file__).with_name('spatial_local_probe.f90')),
        *[str(obj/s) for s in ['xc/fftw_pencils.f90.o','parallel/communication.f90.o','misc/nvtx_wrapper.f90.o']],
        '-L'+os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw')+'/lib','-lfftw3','-o',str(exe)],cwd=tmp,check=True)
    for rank in a.ranks:
        subprocess.run([os.environ.get('MPIEXEC','mpiexec'),'-n',str(rank),str(exe)],check=True,timeout=60,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
