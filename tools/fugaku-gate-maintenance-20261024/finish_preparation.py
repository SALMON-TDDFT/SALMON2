from pathlib import Path
import json,hashlib,subprocess,datetime
old=Path.cwd().resolve();new=old.parent/'parallel-jacobi-omp-plefix-20261024';fix=new/'gate-regression'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def check(gate):
 r=subprocess.run(['python3',str(gate)],cwd=str(fix),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 return r.returncode,r.stdout
red=check(old/'round_gate.py');green=check(new/'round_gate.py');assert red[0]!=0 and green[0]==0
tests={'original_RED':red[0],'sparse_PLE_GREEN':green[0]}
p=fix/'exit.rank-255.txt';v=p.read_text();p.unlink();tests['missing_rank']=check(new/'round_gate.py')[0];assert tests['missing_rank']!=0;p.write_text(v)
p=fix/'launcher-exit.txt';p.write_text('1');tests['failed_launcher']=check(new/'round_gate.py')[0];assert tests['failed_launcher']!=0;p.write_text('0')
p=fix/'exit.rank-255.txt';p.write_text('1');tests['failed_child']=check(new/'round_gate.py')[0];assert tests['failed_child']!=0;p.write_text(v)
root=next(p for p in fix.glob('output.*/**/stdout.*') if p.name.endswith('.0'))
v=root.read_text();root.write_text(v+'FATAL fixture');tests['fatal_stdout']=check(new/'round_gate.py')[0];assert tests['fatal_stdout']!=0;root.write_text(v)
root.write_text(v.replace('MLWF_JACOBI_THREADS 12','MLWF_JACOBI_THREADS 1'));tests['wrong_OMP']=check(new/'round_gate.py')[0];assert tests['wrong_OMP']!=0;root.write_text(v)
root.write_text(v.replace('PARALLEL_ROUND_PASS','MISSING_PASS',1));tests['missing_PASS']=check(new/'round_gate.py')[0];assert tests['missing_PASS']!=0;root.write_text(v)
err=root.parent/'stderr.fixture.0';err.write_text('error fixture');tests['nonempty_stderr']=check(new/'round_gate.py')[0];assert tests['nonempty_stderr']!=0;err.write_text('')
assert check(new/'round_gate.py')[0]==0
(new/'gate-regression.json').write_text(json.dumps({'tests':tests,'old_failure':red[1],'green':green[1]},indent=2))
for p in new.glob('*.py'):compile(p.read_text(),str(p),'exec')
for p in new.glob('*.sh'):subprocess.check_call(['bash','-n',str(p)])
for base in [old,old.parent/'parallel-jacobi-v3']:
 m=json.loads((base/'build-manifest.json').read_text())
 assert sha(base/'salmon')==m['binary_sha256']
 assert sha(base/'round_test')==m['round_binary_sha256']
assert sha(new/'salmon')==sha(old/'salmon') and sha(new/'round_test')==sha(old/'round_test')
files={str(p.relative_to(new)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in new.rglob('*') if p.is_file() and 'data_dcdft' not in p.parts and 'gate-regression' not in p.parts}
record={'time':datetime.datetime.now().isoformat(),'old':str(old),'new':str(new),'files':files,'binary_reused_unchanged':True,'build_or_object_audit_repeated':False,'GS_links_readonly':True,'regression':tests,'resources':'64nodes256MPI12OMP6h','old_jobs_untouched':[52613213,52614422,52615486]}
(new/'plefix-preflight.json').write_text(json.dumps(record,indent=2))
print('PLEFIX_PREPARED',str(new),'REGRESSION',tests,'BINARY',sha(new/'salmon'))

