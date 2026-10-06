"""Sequential MPI2 old/new allocation benchmark; reports process peak RSS."""
from pathlib import Path
import os,subprocess,tempfile
root=Path(__file__).resolve().parents[2];here=root/"testsuites"/"652_dc_lcfo/rt"
with tempfile.TemporaryDirectory() as folder:
 p=Path(folder)
 subprocess.run(['cc','-c',str(here/'peak_rss.c'),'-o',str(p/'rss.o')],check=True)
 subprocess.run(['mpifort','-O2','-fexternal-blas','-fno-tree-loop-vectorize',str(here/'local_action_stubs.f90'),str(root/'src/xc/lcfo_mlwf_links.f90'),str(here/'link_memory_probe.f90'),str(p/'rss.o'),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(p/'probe')],cwd=p,check=True)
 for n in (512,1024):
  for mode in ('dense','tiled'):
   subprocess.run(['mpirun','-np','2',str(p/'probe'),mode,str(n)],cwd=p,check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
