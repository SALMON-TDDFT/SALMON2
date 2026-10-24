from pathlib import Path
import json,math,subprocess,hashlib,shutil,datetime
old=Path.cwd().resolve()
assert old.name=='parallel-jacobi-v3'
w=old.parent/'parallel-jacobi-omp-plefix-20261024'
active=subprocess.check_output(['pjstat'],universal_newlines=True)
assert any('52745397' in l and 'QUE' in l for l in active.splitlines())
for c in ['round','mpi12','omp1','omp12']:
 d=w/(c+'.run');assert not (d/'run.lock').exists() and not list(d.glob('output.*'))
helper="""def memory_samples(rows,child_exit):
 peak=0;count=0;ignored=False;errors=[]
 for j,d in enumerate(rows):
  try:
   rss=d['rss_kib'];huge=d['hugetlb_kib']
   if rss is None or huge is None:
    if rss is None and huge is None and j==len(rows)-1 and count>0 and child_exit=='0':
     ignored=True;continue
    errors.append('missing memory field');continue
   if not all(isinstance(x,(int,float)) and math.isfinite(x) and x>=0 for x in [rss,huge]):
    errors.append('invalid memory sample');continue
   peak=max(peak,rss+huge);count+=1
  except (KeyError,TypeError,ValueError):errors.append('invalid memory sample')
 if not count or peak>5*1024**2:errors.append('rank capacity')
 return peak,count,ignored,errors
"""
scope={'math':math};exec(helper,scope);check=scope['memory_samples']
v={'rss_kib':100,'hugetlb_kib':200};null={'rss_kib':None,'hugetlb_kib':None}
tests={}
cases={'valid':([v],'0',True),'terminal_null':([v,null],'0',True),'middle_null':([v,null,v],'0',False),'only_null':([null],'0',False),'all_empty':([],'0',False),'failed_child':([v,null],'1',False),'partial_null':([v,{'rss_kib':None,'hugetlb_kib':1}],'0',False),'multiple_null':([v,null,null],'0',False),'overcapacity':([{'rss_kib':6*1024**2,'hugetlb_kib':0},null],'0',False),'invalid_nan':([{'rss_kib':float('nan'),'hugetlb_kib':0}],'0',False)}
for name,(rows,exitcode,passed) in cases.items():
 result=check(rows,exitcode);assert (not result[3])==passed,(name,result);tests[name]={'expected_pass':passed,'errors':result[3]}
real=[];rawnull=0
for i in range(256):
 rows=[json.loads(l) for l in (old/'comparison.run'/('samples.rank-'+str(i)+'.jsonl')).read_text().splitlines()]
 result=check(rows,(old/'comparison.run'/('exit.rank-'+str(i)+'.txt')).read_text().strip())
 assert not result[3],(i,result)
 rawnull+=sum(x.get('rss_kib') is None or x.get('hugetlb_kib') is None for x in rows)
 real.append({'rank':i,'peak_KiB':result[0],'valid_count':result[1],'ignored_terminal_null':result[2]})
assert rawnull==82 and sum(x['ignored_terminal_null'] for x in real)==82
p=w/'native_gate.py';before=p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest()=='bef75c4f8ef80f330f5bde5d1f84acea11d5477668c12903cc7b468d1e968996'
start=before.index('peaks=[]');end=before.index("r={'pass':",start)
block="""peaks=[];ignored_terminal_memory_ranks=[]
for i in range(256):
 f=p/('samples.rank-'+str(i)+'.jsonl')
 if not f.is_file():errors.append('missing rank sample');continue
 try:memory_rows=[json.loads(line) for line in f.open()]
 except (ValueError,TypeError):errors.append('invalid memory sample');continue
 peak,count,ignored,memory_errors=memory_samples(memory_rows,ex[i])
 errors.extend(memory_errors);peaks.append(peak)
 if ignored:ignored_terminal_memory_ranks.append(i)
"""
after=before[:start]+helper+block+before[end:]
after=after.replace("'PJM_final':'separate verification required'","'ignored_terminal_memory_ranks':ignored_terminal_memory_ranks,'PJM_final':'separate verification required'")
compile(after,str(p),'exec')
backup=w/'native_gate.py.before-terminal-memory-fix-20261024';assert not backup.exists()
shutil.copy2(str(p),str(backup));p.write_text(after)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((w/'plefix-preflight.json').read_text())
for n,x in m['files'].items():
 if n!='native_gate.py':assert sha(w/n)==x['sha256'],n
record={'time':datetime.datetime.now().isoformat(),'job':52745397,'user_authorized':True,'active_before':active,'before_sha':sha(backup),'after_sha':sha(p),'tests':tests,'real52614422':{'ranks':256,'ignored_final_null_ranks':82,'minimum_valid_samples':min(x['valid_count'] for x in real),'max_sample_KiB':max(x['peak_KiB'] for x in real)},'old_red':'old checker unconditionally flags82 terminal nulls','all_other_static_SHA':'PASS','capacity_threshold_KiB':5*1024**2,'binary_input_job':'unchanged','old52613213':'untouched'}
(w/'terminal-memory-fix-20261024.json').write_text(json.dumps(record,indent=2))
(w/'terminal-memory-real-rank-check-20261024.json').write_text(json.dumps(real,indent=2))
print('TERMINAL_MEMORY_FIX_PASS',json.dumps(record))

