from pathlib import Path
import os, subprocess, tempfile
here = Path(__file__).resolve().parent
root = here.parents[2]
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder);(p/'config.h').write_text('')
    subprocess.run(['mpifort', '-cpp','-ffree-line-length-none','-fallow-argument-mismatch', '-DUSE_MPI', '-I'+folder, '-O2', '-fcheck=all',
                    '-fexternal-blas', '-fno-tree-loop-vectorize',
                    str(root/'src/misc/nvtx_wrapper.f90'), str(root/'src/parallel/communication.f90'),
                    str(root/'src/xc/exx_wannier_gauge.f90'), str(root/'src/xc/lcfo_dist_rows.f90'),
                    str(root/'src/xc/lcfo_seed.f90'), str(here/'seed_stream_probe.f90'),
                    '-L/opt/homebrew/opt/openblas/lib', '-lopenblas', '-o', str(p/'probe')], cwd=p, check=True)
    for ranks in (2, 4):
        subprocess.run(['mpirun', '-np', str(ranks), str(p/'probe')], cwd=p, check=True,
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))

# Compile the same cases without USE_MPI; all communication uses serial branches.
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder);(p/'config.h').write_text('')
    probe = (here/'seed_stream_probe.f90').read_text().replace(' use mpi\n','')
    probe = probe.replace(' call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)', ' rank=0;np=1')
    probe = probe.replace('MPI_COMM_WORLD','0')
    probe = probe.replace('  call MPI_Allreduce(status,lowest,1,MPI_INTEGER,MPI_MIN,0,ierr)', '  lowest=status')
    probe = probe.replace('  call MPI_Allreduce(status,highest,1,MPI_INTEGER,MPI_MAX,0,ierr)', '  highest=status')
    probe = probe.replace(' call MPI_Finalize(ierr)', '')
    assert 'MPI_' not in probe
    (p/'probe.f90').write_text(probe)
    subprocess.run(['gfortran','-cpp','-ffree-line-length-none','-fallow-argument-mismatch','-I'+folder,'-O2','-fcheck=all','-fexternal-blas','-fno-tree-loop-vectorize',
                    str(root/'src/parallel/communication_dummy.f90'),
                    str(root/'src/xc/exx_wannier_gauge.f90'),str(root/'src/xc/lcfo_dist_rows.f90'),
                    str(root/'src/xc/lcfo_seed.f90'),str(p/'probe.f90'),
                    '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(p/'probe')],cwd=p,check=True)
    subprocess.run([str(p/'probe')],cwd=p,check=True,
                   env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
