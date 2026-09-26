"""MPI2/MPI4 algebra and local-storage checks, with optional ScaLAPACK backend."""
from pathlib import Path
import os,subprocess,tempfile,argparse
p=argparse.ArgumentParser();p.add_argument('--compile-only',action='store_true');p.add_argument('--scalapack',action='store_true');p.add_argument('--ranks',type=int,choices=[2,4],default=2);args=p.parse_args()
root=Path(__file__).resolve().parents[2];here=Path(__file__).parent
with tempfile.TemporaryDirectory() as folder:
 exe=Path(folder)/'probe'
 (Path(folder)/'config.h').write_text('')
 command=['mpifort','-I',folder,'-cpp','-DUSE_MPI','-O2','-fexternal-blas','-fno-tree-loop-vectorize','-fcheck=all']
 if args.scalapack:command+=['-DUSE_SCALAPACK']
 command += [str(root/'src/xc'/name) for name in ['hse_ace.f90','lcfo_dist_rows.f90','lcfo_dist_dense.f90']]
 command += [str(here/'distributed_build_probe.f90'),'-L/opt/homebrew/opt/openblas/lib','-lopenblas']
 if args.scalapack:command+=['-L/opt/homebrew/opt/scalapack/lib','-lscalapack']
 subprocess.run(command+['-o',str(exe)],cwd=folder,check=True)
 if not args.compile_only:subprocess.run(['mpirun','-np',str(args.ranks),str(exe)],cwd=folder,check=True,timeout=120,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'))
