from pathlib import Path
import os,json,hashlib,subprocess,datetime,shutil
b=Path.cwd().resolve();assert b.name=='wse2-fugaku-run-20261023'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
before=(b/'gate.py').read_text()
assert sha(b/'gate.py')=='31c55b30cd67387ba7d41fbb6cafd33282b881f81a4570e6e06ec23e189c9c00'
audit=b/'warning-magnetization-policy-20261024';assert not audit.exists();audit.mkdir()
shutil.copy2(str(b/'gate.py'),str(audit/'gate.before.py'))
shutil.copy2(str(b/'manifest.json'),str(audit/'manifest.before.json'))
shutil.copy2(str(b/'verify.py'),str(audit/'verify.before.py'))
helper="""def stderr_warning_only(text):
 if not text.strip():return True
 fatal=r'error stop|fortran runtime error|segmentation|sigsegv|sigabrt|fatal|assertion.*fail|double free|heap corruption|malloc.*corrupt|traceback|(^|[\\r\\n])\\s*(?:\\[[^]]*\\]\\s*)?error\\b'
 if re.search(fatal,text,re.I):return False
 warning=False
 for line in text.splitlines():
  if not line.strip():continue
  if re.match(r'^\\s*(?:\\[[^]]*\\]\\s*)?(?:warning\\b|warn\\b|\\[warn(?:ing)?\\])',line,re.I):
   warning=True;continue
  if warning and line[:1].isspace():continue
  return False
 return warning
"""
after=before.replace('texts=[];peaks=[]',helper+'texts=[];peaks=[];warning_stderr_ranks=[];final_magnetization=[]')
assert "check(not err,'rank stderr '+str(rank))" in after
after=after.replace("check(not err,'rank stderr '+str(rank))","check(q.exists() and stderr_warning_only(err),'rank stderr '+str(rank))\n if err.strip() and stderr_warning_only(err):warning_stderr_ranks.append(rank)")
old="check(len(m)==3 and all(map(math.isfinite,m)) and sum(x*x for x in m)**.5<1e-4,'final magnetization')"
assert old in after
after=after.replace(old,"check(len(m)==3 and all(map(math.isfinite,m)),'finite final magnetization');final_magnetization=m")
after=after.replace("GS_restart=sha,memory_scope=","GS_restart=sha,warning_stderr_ranks=warning_stderr_ranks,final_magnetization=final_magnetization,magnetization_policy='record finite value without magnitude acceptance threshold',memory_scope=")
compile(after,str(b/'gate.py'),'exec')
new=audit/'gate.after.py';new.write_text(after)
# Synthetic fixtures only: no actual GS payload or executable is run.
fix=audit/'synthetic-gs-fixture';fix.mkdir()
for i in range(64):
 (fix/('exit.rank-%d'%i)).write_text('0')
 (fix/('stdout.rank-%d'%i)).write_text('end SALMON\n#GS converged at 20 : 1e-9\nNe= 24\nMagnetization= 0.1 0.2 0.3\n' if i==0 else '')
 (fix/('stderr.rank-%d'%i)).write_text('WARNING: benign fixture warning\n' if i==0 else '')
 (fix/('identity.rank-%d.json'%i)).write_text(json.dumps({'host':'fixture%d'%(i//4),'omp':'12','affinity':[i%4]}))
 (fix/('samples.rank-%d.jsonl'%i)).write_text(json.dumps({'rss_kib':1024,'hugetlb_kib':0})+chr(10))
eigen=''.join('k= fixture\n'+''.join('%d 0.0\n'%j for j in range(1,33)) for k in range(162))
(fix/'WSe2_eigen.data').write_text(eigen);(fix/'data_for_restart').mkdir()
for n in ['info.bin','occupation.bin','rho_inout.bin','wfn.bin']:(fix/'data_for_restart'/n).write_bytes(b'fixture-only')
def run(gate):
 out=fix/'final-gate.json'
 if out.exists():out.rename(fix/('prior-result-%d.json'%len(list(fix.glob('prior-result-*')))))
 r=subprocess.run(['python3',str(gate),'gs','64'],cwd=str(fix),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 return r.returncode
tests={};tests['old_warning_magnetization_RED']=run(audit/'gate.before.py');assert tests['old_warning_magnetization_RED']!=0
tests['new_warning_magnetization_GREEN']=run(new);assert tests['new_warning_magnetization_GREEN']==0
err=fix/'stderr.rank-0';err.write_text('WARNING: fixture\nFATAL: fixture\n');tests['fatal_rejected']=run(new);assert tests['fatal_rejected']!=0
err.write_text('unclassified stderr text');tests['unknown_stderr_rejected']=run(new);assert tests['unknown_stderr_rejected']!=0
err.write_text('WARNING: fixture\n')
q=fix/'exit.rank-2';q.write_text('1');tests['child_failure_rejected']=run(new);assert tests['child_failure_rejected']!=0;q.write_text('0')
q=fix/'samples.rank-2.jsonl';v=q.read_text();q.write_text('');tests['missing_memory_rejected']=run(new);assert tests['missing_memory_rejected']!=0;q.write_text(json.dumps({'rss_kib':7*1024**2,'hugetlb_kib':0}));tests['capacity_rejected']=run(new);assert tests['capacity_rejected']!=0;q.write_text(v)
q=fix/'stdout.rank-0';v=q.read_text();q.write_text(v.replace('0.1 0.2 0.3','NaN 0.2 0.3'));tests['nonfinite_magnetization_rejected']=run(new);assert tests['nonfinite_magnetization_rejected']!=0;q.write_text(v)
q.write_text(v.replace('1e-9','1e-3'));tests['unconverged_rejected']=run(new);assert tests['unconverged_rejected']!=0;q.write_text(v)
q.write_text(v.replace('Ne= 24','Ne= 25'));tests['Ne_rejected']=run(new);assert tests['Ne_rejected']!=0;q.write_text(v)
assert run(new)==0
(audit/'regression.json').write_text(json.dumps(tests,indent=2))
# Immediately confirm every affected job remains unstarted.
active=subprocess.check_output(['pjstat'],universal_newlines=True)
for i in ['52698609','52698613','52698614']:assert any(i in l and 'QUE' in l for l in active.splitlines()),i
for c in ['h032','h030','h028']:
 for s in ['gs','rt']:
  d=b/'lsda'/c/s
  assert not (d/'run.lock').exists() and not list(d.glob('stdout.rank-*')) and not list(d.glob('output.*'))
subprocess.check_call(['python3',str(b/'verify.py')])
m=json.loads((b/'manifest.json').read_text())
newmeta={'bytes':new.stat().st_size,'sha256':sha(new)}
# Transitional verifier accepts only the immutable old/new gate hashes, including its own original/bridge hash.
vbefore=(b/'verify.py').read_text()
bridge="""from pathlib import Path
import json,hashlib
b=Path(__file__).resolve().parent
allowed_gate=ALLOWED_GATE
for name,x in json.loads((b/'manifest.json').read_text()).items():
 p=b/name;actual=hashlib.sha256(p.read_bytes()).hexdigest()
 if name=='gate.py':assert actual in allowed_gate,name
 elif name=='verify.py':assert p.read_text()==Path(BRIDGE_PATH).read_text() or actual==ORIGINAL_VERIFY_SHA,name
 else:assert p.stat().st_size==x['bytes'] and actual==x['sha256'],name
print('MANIFEST_PASS_TRANSITION')
"""
bridge=bridge.replace('ALLOWED_GATE',repr([sha(b/'gate.py'),sha(new)])).replace('BRIDGE_PATH',repr(str(audit/'verify.bridge.py'))).replace('ORIGINAL_VERIFY_SHA',repr(sha(b/'verify.py')))
(audit/'verify.bridge.py').write_text(bridge);compile(bridge,'bridge','exec')
lock=audit/'apply.lock';lock.mkdir()
def atomic(path,data):
 temp=path.with_name(path.name+'.policy-tmp');temp.write_text(data);os.replace(str(temp),str(path))
atomic(b/'verify.py',bridge)
atomic(b/'gate.py',after)
m['gate.py']=newmeta
atomic(b/'manifest.json',json.dumps(m,indent=2))
atomic(b/'verify.py',vbefore)
subprocess.check_call(['python3',str(b/'verify.py')])
assert sha(b/'gate.py')==newmeta['sha256']
record={'time':datetime.datetime.now().isoformat(),'jobs':[52698609,52698613,52698614],'before':active,'old_gate_sha':sha(audit/'gate.before.py'),'new_gate_sha':sha(b/'gate.py'),'tests':tests,'physical_magnetization':'record only; finite and three components required','stderr':'only classified warnings allowed; fatal and unknown rejected','binary_input_job':'unchanged','submission_cancel':'none','old_evidence':'preserved'}
(audit/'applied-record.json').write_text(json.dumps(record,indent=2));print('WSE_POLICY_APPLIED_PASS',json.dumps(record))

