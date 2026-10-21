"""Run the real block-cyclic regression with MPI-only and OpenMP executables.
Example: python3 run_omp_rounds.py /absolute/mpi_round /absolute/omp_round --ranks 4
Elapsed totals include launcher/setup; use MLWF_PARALLEL_SWEEP for localization timing.
"""
import argparse,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('candidate',type=Path)
p.add_argument('--mpiexec',default='mpiexec');p.add_argument('--ranks',type=int,default=4)
a=p.parse_args();assert a.ranks>0
for name,exe,threads in [('MPI-only',a.baseline,1)]+[('OMP',a.candidate,t) for t in (1,2,4,12)]:
 env=dict(os.environ,OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 start=time.monotonic()
 r=subprocess.run([a.mpiexec,'--bind-to','none','--oversubscribe','-n',str(a.ranks),str(exe.resolve())],env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,timeout=60)
 assert r.returncode==0,(name,threads,r.stdout,r.stderr)
 assert r.stdout.count('PARALLEL_ROUND_PASS')==10,r.stdout
 if name=='OMP':
  teams={int(l.split()[-1]) for l in r.stdout.splitlines() if 'MLWF_JACOBI_THREADS' in l}
  assert teams=={threads},(threads,teams)
 print(json.dumps(dict(version=name,omp=threads,mpi=a.ranks,seconds=time.monotonic()-start,pass_=True)),flush=True)
