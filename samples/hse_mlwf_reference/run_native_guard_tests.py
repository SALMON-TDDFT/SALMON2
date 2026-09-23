"""Run prepared negative native HSE restart fixtures (no accepted steps allowed)."""
import argparse,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('binary');p.add_argument('fixtures');p.add_argument('--mpiexec',default='mpiexec');a=p.parse_args()
expect={'bad_dt':'restart physics','bad_impulse':'restart physics','bad_metadata':'restart physics','bad_missing_metadata':'restart physics','bad_pseudo':'restart physics','bad_frozen':'self-consistent functional'}
result={}
for case,reason in expect.items():
 directory=Path(a.fixtures)/case
 with (directory/'inputfile').open('rb') as inp:
  run=subprocess.run([a.mpiexec,'-n','1',a.binary],stdin=inp,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,cwd=directory,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),timeout=60)
 text=run.stdout.decode(errors='replace');(directory/'run.log').write_text(text)
 assert run.returncode!=0 and reason in text and 'HSE_PT_CN' not in text,(case,run.returncode,text[-1500:])
 result[case]=dict(rejected=True,exit_code=run.returncode,expected_reason=reason)
print(json.dumps(result,indent=2))
