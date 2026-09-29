"""Compare distributed Gamma exchange to serial MLWF exchange using build objects."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import os
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);args=p.parse_args()
b=args.build.resolve();obj=b/'src/CMakeFiles/salmon.dir'
names=['xc/exx_orbitals.f90.o','xc/exx_ace.f90.o','xc/exx_pair_candidates.f90.o','xc/exx_spatial.f90.o','xc/exx_batch_backend.f90.o','xc/exx_spatial_local.f90.o','xc/exx_wannier_gauge.f90.o','xc/exx_wannier.f90.o',
       'xc/exx_local_fft.f90.o','xc/lcfo_wf_support.f90.o','xc/fftw_pencils.f90.o',
       'parallel/communication.f90.o','misc/nvtx_wrapper.f90.o']
with tempfile.TemporaryDirectory() as tmp:
    exe=Path(tmp)/'probe'
    subprocess.run(['mpifort','-fopenmp','-fcheck=all','-ffree-line-length-none','-I'+str(b),str(Path(__file__).with_name('exchange_driver.f90')),*[str(obj/n) for n in names],'-L/opt/homebrew/opt/fftw/lib','-lfftw3','-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],check=True)
    for ranks, orbitals in ((1,1),(2,1),(4,1),(2,2),(4,2),(8,2),(4,4),(8,4),(8,8)):
        for omega in ('0', '.11', '.3'):
            subprocess.run(['mpiexec','-n',str(ranks),str(exe),omega,str(orbitals)],check=True,timeout=15,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
