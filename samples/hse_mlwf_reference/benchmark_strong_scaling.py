"""Fixed Si Taylor+ACE strong scaling; macOS RSS sampling, sequential MPI runs."""
import argparse, json, os, re, signal, subprocess, time
from pathlib import Path
import numpy as np

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--executable',type=Path,required=True)
 parser.add_argument('--pause-reference',type=int,help='Verified MPI launcher PID for this repository reference job')
 parser.add_argument('--steps',type=int,default=5)
 args=parser.parse_args();exe=str(args.executable.resolve())
 root=Path(__file__).resolve().parents[2];base=root/'calculations/si_hse_native'
 out=root/'docs/results/si-hse-strong-scaling';out.mkdir(exist_ok=True)
 template=(base/'memory_tuned_hse_taylor4_dt0.16/inputfile').read_text()
 def processes():
  lines=subprocess.check_output(['/bin/ps','-axo','pid=,ppid=,state=,rss=,args='],text=True).splitlines()
  return [dict(pid=int(a[0]),ppid=int(a[1]),state=a[2],rss=int(a[3]),args=a[4]) for a in (s.split(None,4) for s in lines) if len(a)==5]
 def descendants(pid,records):
  group={pid}
  for _ in range(10):group.update(r['pid'] for r in records if r['ppid'] in group)
  return group
 paused=[];runs=[];proc=None
 def interrupted(signum,frame):raise KeyboardInterrupt
 signal.signal(signal.SIGTERM,interrupted)
 try:
  if args.pause_reference:
   records=processes();group=descendants(args.pause_reference,records)
   selected=[r for r in records if r['pid'] in group]
   assert len(selected)==9 and all('samples/hse_mlwf_reference/mpi_run.py' in r['args'] for r in selected)
   assert all('T' not in r['state'] for r in selected),'Reference already suspended; refusing to change its state'
   for r in selected:
    os.kill(r['pid'],signal.SIGSTOP);paused.append(r['pid'])
   (out/'reference_pause.json').write_text(json.dumps({'paused_pids':paused,'resumed':False},indent=2)+'\n')
  for rep in range(2):
   for ranks in ([1,2,4,8,16] if rep==0 else [16,8,4,2,1]):
    name=f'strong_taylor_n{ranks}_r{rep}';d=base/name;d.mkdir(exist_ok=True)
    s=re.sub(r'(?m)^(\s*nt\s*=).*',r'\g<1> '+str(args.steps),template)
    s=re.sub(r'checkpoint_interval=\d+','checkpoint_interval='+str(args.steps),s)
    s=re.sub(r'nproc_k=\d+','nproc_k='+str(ranks),s);(d/'inputfile').write_text(s)
    with (d/'inputfile').open() as inp,(d/'run.log').open('w') as log:
     proc=subprocess.Popen(['mpiexec','-n',str(ranks),'/usr/bin/env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1',exe],cwd=d,stdin=inp,stdout=log,stderr=subprocess.STDOUT)
     peak_sum=peak_rank=0;samples=0;max_rank_vector=[]
     while proc.poll() is None:
      records=processes();group=descendants(proc.pid,records)
      memory=[r['rss'] for r in records if r['pid'] in group and r['args']==exe]
      if len(memory)==ranks:
       if sum(memory)>peak_sum:peak_sum=sum(memory);max_rank_vector=memory
       peak_rank=max(peak_rank,max(memory));samples+=1
      time.sleep(.1)
    txt=(d/'run.log').read_text()
    if proc.returncode or 'end SALMON' not in txt:raise RuntimeError(name)
    assert samples>0,'No RSS samples matched the executable'
    rt=float(re.search(r'^rt iterations\s+(.+)$',txt,re.M)[1].split()[1])
    total=float(re.search(r'total calculation time,([\d.]+)',txt)[1])
    row=dict(case=name,ranks=ranks,repeat=rep,steps=args.steps,rt_seconds=rt,total_seconds=total,seconds_per_step=rt/args.steps,sampled_peak_sum_rss_KiB=peak_sum,sampled_peak_rank_rss_KiB=peak_rank,rss_vector_at_sum_peak_KiB=max_rank_vector,samples=samples)
    ref=base/f'strong_taylor_n1_r0/checkpoint_rt_{args.steps:06d}/wfn.bin'
    u=np.fromfile(ref,np.complex128);v=np.fromfile(d/f'checkpoint_rt_{args.steps:06d}/wfn.bin',np.complex128)
    row['orbital_relative_difference']=float(np.linalg.norm(v-u)/np.linalg.norm(u))
    assert row['orbital_relative_difference']<1e-10
    runs.append(row);(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(row,flush=True)
 finally:
  try:
   if proc is not None and proc.poll() is None:
    proc.terminate();proc.wait(timeout=30)
  finally:
   resumed=[]
   for pid in paused:
    try:os.kill(pid,signal.SIGCONT);resumed.append(pid)
    except ProcessLookupError:pass
   if paused:(out/'reference_pause.json').write_text(json.dumps({'paused_pids':paused,'resumed_pids':resumed,'resumed':len(resumed)==len(paused)},indent=2)+'\n')
if __name__=='__main__':main()
