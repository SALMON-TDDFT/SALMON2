"""Run unchanged mesh MLWF/QR/localization/polar code; only communication is stubbed in serial."""
from pathlib import Path
import argparse,os,subprocess,tempfile
parser=argparse.ArgumentParser()
parser.add_argument('--mpi',action='store_true')
args=parser.parse_args()
here=Path(__file__).resolve().parent;root=here.parents[1]
with tempfile.TemporaryDirectory() as folder:
 p=Path(folder);(p/'config.h').write_text('')
 stub=root/'testsuites/unit_lcfo_rt'/('local_action_stubs.f90' if args.mpi else 'transport_stubs.f90')
 sources=['hse_wannier_gauge.f90','hse_ace.f90','lcfo_dist_rows.f90','lcfo_dist_dense.f90',
          'lcfo_seed.f90','lcfo_mlwf_links.f90','hse_grid_wannier.f90']
 command=['mpifort' if args.mpi else 'gfortran','-cpp','-I'+folder,'-O2','-fcheck=all',
          '-fexternal-blas','-fno-tree-loop-vectorize']
 if args.mpi:command+=['-DUSE_MPI']
 command += [str(stub)]+[str(root/'src/xc'/name) for name in sources]+[str(here/'probe.f90'),
             '-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(p/'probe')]
 subprocess.run(command,cwd=p,check=True)
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 for ranks in ([2,4] if args.mpi else [1]):
  launch=['mpirun','-np',str(ranks)] if args.mpi else []
  run=subprocess.run(launch+[str(p/'probe')],cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  assert run.returncode==0,run.stdout
  assert 'PASS mesh MLWF' in run.stdout,run.stdout
  assert 'WARNING Grid MLWF radius: sphere norm below 99.9%' in run.stdout,run.stdout
  print('PASS grid MLWF ranks',ranks)
 if not args.mpi:
  run=subprocess.run([str(p/'probe'),'1'],cwd=p,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  assert run.returncode!=0 and 'invalid source' in run.stdout,run.stdout
  print('PASS invalid disabled finite radius')
