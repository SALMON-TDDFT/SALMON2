"""A pulse restart must not be accepted by either the pulse or impulse path."""
import argparse,os,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('exe');p.add_argument('work');a=p.parse_args();b=Path(a.work);source=b/'probe';results={}
for name in ['pulse_restart','pulse_as_impulse']:
 d=b/name;d.mkdir(exist_ok=True);s=(source/'inputfile').read_text().replace("directory_read_data='restart/'","directory_read_data='../probe/checkpoint_rt_000016/'").replace('nt=16','nt=17')
 s=s.replace("yn_restart='n'","yn_restart='y'")
 if name=='pulse_as_impulse':s=s.replace("ae_shape1='Acos2'","ae_shape1='impulse'").replace('ae_shape2="impulse"','ae_shape2="none"')
 (d/'inputfile').write_text(s)
 r=subprocess.run(['mpiexec','--bind-to','none','-n','4',str(Path(a.exe).resolve())],input=s,text=True,cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),timeout=90)
 (d/'run.log').write_text(r.stdout)
 expected='laser RT restart' if name=='pulse_restart' else 'restart physics'
 results[name]=dict(exit_code=r.returncode,expected=expected,passed=r.returncode!=0 and expected in r.stdout)
 assert results[name]['passed'],r.stdout[-2000:]
(b/'restart_guard_metrics.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
