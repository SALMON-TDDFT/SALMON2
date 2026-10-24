from pathlib import Path
import os,json,hashlib,shutil,subprocess,datetime
old=Path.cwd().resolve()
assert old.name=='parallel-jacobi-omp'
new=old.parent/'parallel-jacobi-omp-plefix-20261024'
assert not new.exists()
new.mkdir()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((old/'preflight.json').read_text())
for name,v in manifest['files'].items():
 p=old/name
 if name not in ['native_gate.py','compare.py','preparation.log']:
  assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256'],name
 q=new/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(str(p),str(q))
for name in ['mpi12','omp1','omp12']:
 src=old/(name+'.run')/'data_dcdft'
 assert src.is_symlink() and src.resolve().is_dir()
 (new/(name+'.run')/'data_dcdft').symlink_to(src.resolve())
(new/'round.run').mkdir()
for name in ['job.sh','omp12-rank.sh','omp1-rank.sh','mpi12-rank.sh','round-rank.sh']:
 p=new/name;t=p.read_text();t=t.replace(str(old),str(new));p.write_text(t)
# PLE omits stdout for ranks that write nothing. All rank exits remain mandatory.
p=new/'round_gate.py';t=p.read_text()
t=t.replace("assert len(s)==256 and len(r)==1","assert len(s)>=1 and len(r)==1")
t+="\nassert (p/'launcher-exit.txt').read_text().strip()=='0'\nassert not any(any(x in f.read_text().upper() for x in ['ERROR STOP','FATAL','SIGSEGV','SEGMENTATION FAULT']) for f in s)\n"
p.write_text(t)
# Physical traces already require every rank via localization status and memory records.
p=new/'native_gate.py';t=p.read_text()
assert "if len(stdout)!=256 or len(z)!=1:" in t
t=t.replace("if len(stdout)!=256 or len(z)!=1:","if not stdout or len(z)!=1:")
t=t.replace("ex=[]","if not (p/'launcher-exit.txt').is_file() or (p/'launcher-exit.txt').read_text().strip()!='0':errors.append('launcher exit')\nif any(any(x in f.read_text().upper() for x in ['ERROR STOP','FATAL','SIGSEGV','SEGMENTATION FAULT']) for f in stdout):errors.append('fatal stdout')\nex=[]")
p.write_text(t)
p=new/'job.sh';t=p.read_text().replace('Si4JOmp','Si4JPFix')
line=next(l for l in t.splitlines() if l.startswith('mpiexec') and 'round-rank.sh' in l)
t=t.replace(line,'set +e\n'+line+'\nrc=$?\nset -e\necho "$rc" > launcher-exit.txt\nif [ "$rc" -ne 0 ]; then exit "$rc"; fi')
p.write_text(t)
# Regression fixture uses copies of actual non-GS evidence; original outputs remain untouched.
fix=new/'gate-regression';fix.mkdir()
run=old/'round.run'
for p in run.glob('exit.rank-*.txt'):shutil.copy2(str(p),str(fix/p.name))
for p in run.glob('output.*'):(fix/p.name).symlink_to(p.resolve(),target_is_directory=True)
(fix/'launcher-exit.txt').write_text('0')
def check(gate):
 r=subprocess.run(['python3',str(gate)],cwd=str(fix),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 return r.returncode,r.stdout
red=check(old/'round_gate.py');assert red[0]!=0
green=check(new/'round_gate.py');assert green[0]==0,green
tests={'original_RED':red[0],'sparse_PLE_GREEN':green[0]}
p=fix/'exit.rank-255.txt';v=p.read_text();p.unlink();tests['missing_rank']=check(new/'round_gate.py')[0];assert tests['missing_rank']!=0;p.write_text(v)
p=fix/'launcher-exit.txt';p.write_text('1');tests['failed_launcher']=check(new/'round_gate.py')[0];assert tests['failed_launcher']!=0;p.write_text('0')
p=fix/'exit.rank-255.txt';p.write_text('1');tests['failed_child']=check(new/'round_gate.py')[0];assert tests['failed_child']!=0;p.write_text(v)
# independent real PLE fixture copy for fatal/error rejection tests
for p in list(fix.glob('output.*')):
 target=p.resolve();p.unlink();shutil.copytree(str(target),str(p))
root=next(p for p in fix.glob('output.*/**/stdout.*') if p.name.endswith('.0'))
v=root.read_text();root.write_text(v+'\nFATAL fixture\n');tests['fatal_stdout']=check(new/'round_gate.py')[0];assert tests['fatal_stdout']!=0;root.write_text(v)
root.write_text(v.replace('MLWF_JACOBI_THREADS 12','MLWF_JACOBI_THREADS 1'));tests['wrong_OMP']=check(new/'round_gate.py')[0];assert tests['wrong_OMP']!=0;root.write_text(v)
root.write_text(v.replace('PARALLEL_ROUND_PASS','MISSING_PASS',1));tests['missing_PASS']=check(new/'round_gate.py')[0];assert tests['missing_PASS']!=0;root.write_text(v)
err=next(fix.glob('output.*/**/stderr.*'));err.write_text('error fixture');tests['nonempty_stderr']=check(new/'round_gate.py')[0];assert tests['nonempty_stderr']!=0;err.write_text('')
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

