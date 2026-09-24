"""Uses pre-existing compare_* inputs and reference_restart; see memory report."""
from pathlib import Path
import subprocess,os,time,re,json
import argparse
parser=argparse.ArgumentParser(description='Sequential native HSE old/new MPI8 RSS and timing pilot (macOS ps).')
parser.add_argument('--before',required=True,type=Path)
parser.add_argument('--after',required=True,type=Path)
args=parser.parse_args()
p=Path(__file__).resolve().parents[2]/'calculations/si_hse_native';rows=[]
for mode in ['hse_taylor4_dt0.16','hse_ptcn_dt0.32','hse_taylor4_full_dt0.16']:
 for version in ['old','new']:
  d=p/f'memory_{version}_{mode}';d.mkdir(exist_ok=True)
  s=(p/f'compare_{mode}'/'inputfile').read_text();s=re.sub(r'(?m)^(\s*nt\s*=).*',r'\g<1> 5',s);s=re.sub(r'checkpoint_interval=\d+','checkpoint_interval=5',s);s=s.replace('nproc_k=4','nproc_k=8');(d/'inputfile').write_text(s)
  exe=str((args.before if version=='old' else args.after).resolve())
  with (d/'inputfile').open() as inp,(d/'run.log').open('w') as log:
   proc=subprocess.Popen(['/opt/homebrew/bin/mpiexec','-n','8','/usr/bin/env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1',exe],cwd=d,stdin=inp,stdout=log,stderr=subprocess.STDOUT)
   peak_sum=peak_rank=0;samples=0
   while proc.poll() is None:
    ps=subprocess.run(['/bin/ps','-axo','pid=,ppid=,rss=,comm='],capture_output=True,text=True,check=True)
    records=[]
    for line in ps.stdout.splitlines():
     a=line.split(None,3)
     if len(a)==4:records.append((int(a[0]),int(a[1]),int(a[2]),a[3]))
    descendants={proc.pid}
    for _ in range(5):descendants.update(a[0] for a in records if a[1] in descendants)
    mem=[a[2] for a in records if a[0] in descendants and a[3]==exe]
    if len(mem)==8:peak_sum=max(peak_sum,sum(mem));peak_rank=max(peak_rank,max(mem));samples+=1
    time.sleep(.1)
  txt=(d/'run.log').read_text()
  if proc.returncode or 'end SALMON' not in txt:raise RuntimeError(str(d))
  rt=float(re.search(r'^rt iterations\s+(.+)$',txt,re.M)[1].split()[1]);row=dict(case=mode,version=version,mpi=8,steps=5,rt_seconds=rt,seconds_per_step=rt/5,sampled_peak_sum_rss_KiB=peak_sum,sampled_peak_rank_rss_KiB=peak_rank,samples=samples);rows.append(row);print(row,flush=True);(p/'memory_benchmark.json').write_text(json.dumps(rows,indent=2)+'\n')
