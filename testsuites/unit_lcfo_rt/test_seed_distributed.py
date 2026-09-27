"""Distributed seed QR versus independent root QR, including snapshot order."""
from pathlib import Path
import os, subprocess, tempfile
here = Path(__file__).resolve().parent
root = here.parents[1]
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder); (p/'config.h').write_text('')
    subprocess.run(['mpifort', '-cpp', '-DUSE_MPI', '-DUSE_SCALAPACK', '-I'+folder,
                    '-O2', '-fcheck=all', '-fexternal-blas', '-fno-tree-loop-vectorize',
                    str(root/'src/xc/hse_wannier_gauge.f90'), str(root/'src/xc/lcfo_dist_rows.f90'),
                    str(root/'src/xc/lcfo_seed.f90'), str(here/'seed_stream_probe.f90'),
                    '-L/opt/homebrew/opt/scalapack/lib', '-lscalapack',
                    '-L/opt/homebrew/opt/openblas/lib', '-lopenblas', '-o', str(p/'probe')], cwd=p, check=True)
    for ranks in (1, 2, 3, 4):
        result = subprocess.run(['mpirun', '-np', str(ranks), str(p/'probe')], cwd=p,
            check=True, timeout=120, stdout=subprocess.PIPE, universal_newlines=True,
            env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                     SALMON_LCFO_SEED_DISTRIBUTED='1'))
        assert 'LCFO seed: distributed pivoted QR' in result.stdout, 'Distributed backend not used'
        print(ranks, result.stdout.strip())

# A requested backend must never silently fall back in unsupported builds.
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder); (p/'config.h').write_text('')
    (p/'unsupported.f90').write_text('''program unsupported
#ifdef USE_MPI
 use mpi
#endif
 use lcfo_seed, only: lcfo_seed_gamma
 implicit none
 integer :: np,rank,ierr,status
 integer,allocatable :: counts(:)
 complex(8) :: local(2,2),u(2,2)
 np=1;rank=0
#ifdef USE_MPI
 call MPI_Init(ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
#endif
 allocate(counts(np));counts=2;local=1d0
#ifdef USE_MPI
 call lcfo_seed_gamma(local,counts,MPI_COMM_WORLD,u,status)
 call MPI_Finalize(ierr)
#else
 call lcfo_seed_gamma(local,counts,0,u,status)
#endif
end program
''')
    for mpi in (False, True):
        command = ['mpifort' if mpi else 'gfortran', '-cpp', '-I'+folder, '-O0']
        if mpi:command += ['-DUSE_MPI']
        command += [str(root/'src/xc'/name) for name in
                    ('hse_wannier_gauge.f90','lcfo_dist_rows.f90','lcfo_seed.f90')]
        command += [str(p/'unsupported.f90'),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(p/'probe')]
        subprocess.run(command,cwd=p,check=True)
        result = subprocess.run((['mpirun','-np','2'] if mpi else [])+[str(p/'probe')],
            cwd=p,timeout=30,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',SALMON_LCFO_SEED_DISTRIBUTED='1'))
        assert result.returncode != 0 and 'requires MPI and ScaLAPACK' in result.stdout, result.stdout
        print('Unsupported backend rejected:', 'MPI-only' if mpi else 'serial')
