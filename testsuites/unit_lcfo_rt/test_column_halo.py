from pathlib import Path
import subprocess,tempfile,os
r=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as d:
 p=Path(d);(p/'config.h').write_text('#define USE_MPI\n')
 subprocess.run(['mpifort','-cpp','-O2','-fcheck=all','-fno-tree-loop-vectorize','-I'+d,str(r/'src/xc/lcfo_dist_rows.f90'),str(Path(__file__).with_name('column_halo_probe.f90')),'-o',str(p/'probe')],cwd=d,check=True)
 for n in (2,4):
  subprocess.run(['mpiexec','-n',str(n),str(p/'probe')],cwd=d,check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
