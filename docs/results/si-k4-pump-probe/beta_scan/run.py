from pathlib import Path
import os,sys,json,time,signal,subprocess,fcntl
out=Path(__file__).resolve().parent
root=out.parents[3]
base=root/'calculations/si_hse_native/k4_pump_probe'
source=base/'td_g40_ground'
def save(p,d):
 tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2));tmp.replace(p)
with (out/'run.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 results=[]
 for beta in [.001,.004,.016]:
  name=f'td_g40_b{round(beta*1000):03d}_ground';d=base/name;d.mkdir(exist_ok=True)
  if (d/'status.json').exists():
   old=json.loads((d/'status.json').read_text())
   if old.get('completed'):results.append(old);continue
   raise RuntimeError(f'Existing incomplete case: {name}')
  (d/'inputfile').write_text((source/'inputfile').read_text().replace('tdcdft_damping=0.0',f'tdcdft_damping={beta}'))
  for link in ['restart','Si_rps.dat']:
   if not (d/link).exists():(d/link).symlink_to((source/link).resolve())
  start=time.monotonic();peak=0;reason=None
  with (d/'inputfile').open() as inp,(d/'run.log').open('w') as log:
   p=subprocess.Popen(['/opt/homebrew/bin/mpiexec','--bind-to','none','-n','4','/private/tmp/salmon-hse-symmetry/salmon'],cwd=d,stdin=inp,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE'))
   save(out/'status.json',dict(status='running',case=name,pid=p.pid,results=results))
   try:
    while p.poll() is None:
     rows=subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True);rows=[list(map(int,r.split())) for r in rows.splitlines()];ids={p.pid}
     for _ in range(5):ids.update(r[0] for r in rows if r[1] in ids)
     rss=sum(r[2] for r in rows if r[0] in ids);peak=max(peak,rss)
     if rss>48*1024**2 or time.monotonic()-start>21600:
      reason='resource_limit';os.killpg(p.pid,signal.SIGTERM);p.wait(timeout=30);break
     time.sleep(5)
   finally:
    if p.poll() is None:os.killpg(p.pid,signal.SIGTERM);p.wait(timeout=30)
  item=dict(case=name,beta=beta,exit_code=p.returncode,completed=p.returncode==0 and 'end SALMON' in (d/'run.log').read_text(),wall_seconds=time.monotonic()-start,peak_RSS_KiB=peak,stop_reason=reason)
  save(d/'status.json',item);results.append(item)
  if not item['completed']:save(out/'status.json',dict(status='failed',results=results));sys.exit(1)
 subprocess.run([sys.executable,str(out/'analyze.py')],check=True)
 save(out/'status.json',dict(status='completed',results=results))
