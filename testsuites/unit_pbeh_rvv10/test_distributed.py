#!/usr/bin/env python3
"""Direct MPI rVV10 comparisons using the configured native build objects.
Run: python3 test_distributed.py --build /path/to/USE_HSE_MPI_build
"""
import argparse
from pathlib import Path
import os
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark',action='store_true')
    parser.add_argument('--scale',type=int,default=1,choices=range(1,5))
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--mpirun',default='mpiexec')
    parser.add_argument('--mpifort',default='mpifort')
    parser.add_argument('--fftw-prefix',type=Path,default=Path('/opt/homebrew'))
    args=parser.parse_args()
    build=args.build.resolve();obj=build/'src/CMakeFiles/salmon.dir'
    objects=['xc/fftw_pencils.f90.o','xc/rvv10.f90.o','xc/rvv10_distributed.f90.o','xc/rvv10_ffte.f.o',
             'parallel/communication.f90.o','misc/nvtx_wrapper.f90.o']
    objects += ['ext/FFTE/'+name+'.f.o' for name in ['factor','kernel','fft235','pzfft3dv_mod']]
    with tempfile.TemporaryDirectory(prefix='rvv10-distributed-') as tmp:
        exe=Path(tmp)/'probe'
        subprocess.run([args.mpifort,'-fopenmp','-fcheck=all','-ffree-line-length-none','-I'+str(build),
            str(Path(__file__).with_name('fftw_benchmark.f90' if args.benchmark else 'distributed_probe.f90')),
            *[str(obj/name) for name in objects],'-L'+str(args.fftw_prefix/'lib'),'-lfftw3','-o',str(exe)],check=True)
        for ranks in (2,4):
            for threads in ((1,) if args.benchmark else (1,2)):
                env=dict(os.environ,OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1')
                for nq in ((args.scale,) if args.benchmark else (8,16,32)):
                    subprocess.run([args.mpirun,'-n',str(ranks),str(exe),str(nq)],env=env,check=True,timeout=90)
                print(f'Passed {ranks} MPI ranks, {threads} threads',flush=True)

if __name__=='__main__':main()
