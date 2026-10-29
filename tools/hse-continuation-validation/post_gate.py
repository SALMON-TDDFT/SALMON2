from pathlib import Path
import json,re,math,collections,os
from continuation_checks import check_endpoints,check_pulse
from checkpoint_v3 import validate
p=Path.cwd();cfg=json.loads((p/'preflight.json').read_text());start=cfg['start_step'];target=cfg['nt'];k=cfg['exact_interval']
out=Path(os.environ.get('NR16_GATE_OUTPUT',str(p)))
assert (p/'launcher-exit.txt').read_text().strip()=='0'
streams={}
for f in p.glob('output.*/**/stdout.*'):
 rank=int(f.name.rsplit('.',1)[1]);assert rank not in streams;streams[rank]=f
assert set(streams)==set(range(cfg['mpi'])),('missing or extra PLE stdout',len(streams))
for f in list(p.glob('output.*/**/stderr.*'))+list(p.glob('stderr.log')):
 assert f.read_text().strip()=='',('unexpected stderr',str(f))
root=streams[0].read_text();assert 'end SALMON' in root
stop=[int(x) for x in re.findall(r'shutdown the calculation, iter\s*=\s*(\d+)',root)]
assert len(set(stop))<=1
end=stop[-1] if stop else target;assert start<end<=target
incomplete=end<target;assert incomplete==bool(stop)
if incomplete:
 cps=list(p.glob('checkpoint_rt_*'));assert len(cps)==1 and cps[0].name=='checkpoint_rt_%06d'%end
else:assert not list(p.glob('checkpoint_rt_*'))
hosts=collections.defaultdict(list);svdmax=[0.,0.,0.];svdtotal=0
for rank in range(cfg['mpi']):
 d=json.loads((p/('summary.rank-%d.json'%rank)).read_text());ident=d['identity']
 assert ident['rank']==rank and ident['omp']=='12' and d['exit_code']==0 and d['samples']>0
 assert d['child_rusage_maxrss_units']=='KiB' and d['child_rusage_maxrss']>0
 assert (p/('samples.rank-%d.jsonl.gz'%rank)).stat().st_size>0
 cores=set(ident['affinity']);peak=d['sample_max_RSS_plus_Huge_KiB']
 assert len(cores)==12 and math.isfinite(peak) and 0<peak<6*1024**2
 hosts[ident['host']].append((cores,peak))
 text=streams[rank].read_text(errors='strict')
 assert not re.search(r'ERROR STOP|MPI_ABORT|out of memory|SIGSEGV|severe|fatal|(?<![A-Za-z])nan(?![A-Za-z])|(?<![A-Za-z])[+-]?inf(?:inity)?(?![A-Za-z])|FACTOR_SVD_INFO',text,re.I),rank
 if rank<80:
  teacher=[x.split(':',1)[1].split() for x in text.splitlines() if 'FACTOR_HISTORY k/teacher/step/count/coeff/ready:' in x]
  svds=[x.split(':',1)[1].split() for x in text.splitlines() if 'FACTOR_SVD_QUALITY n/residual/unit_u/unit_v/smin/smax:' in x]
  check_endpoints([int(x[2]) for x in teacher],len(svds),end,k,start=start)
  for x in teacher:assert int(x[0])==rank+1 and int(x[1])==int(x[2])//k-2
  for x in svds:
   n=int(x[0]);v=list(map(float,x[1:]));assert n==16 and len(v)==5 and all(math.isfinite(y) for y in v)
   assert min(v[:3])>=0 and v[0]<=1e-10 and max(v[1:3])<=1e-10*math.sqrt(n) and v[4]>0 and v[3]>v[4]*1e-12
   for j in range(3):svdmax[j]=max(svdmax[j],v[j])
  svdtotal+=len(svds)
 else:assert 'EXX_SERVICE_COMPUTE' in text and 'EXX_SERVICE_STOP' in text,rank
assert len(hosts)==cfg['nodes']
nodepeaks=[]
for group in hosts.values():
 assert len(group)==4
 for i,(cores,peak) in enumerate(group):
  for prev,_ in group[:i]:assert not cores.intersection(prev)
 nodepeaks.append(sum(x[1] for x in group));assert nodepeaks[-1]<24*1024**2
flags=[x.split(':',1)[1].split() for x in root.splitlines() if 'FACTOR_HISTORY step/predicted:' in x]
assert [(int(x[0]),x[1]) for x in flags]==[(i,'F' if i%k==0 else 'T') for i in range(start+1,end+1)]
physical=[]
for line in root.splitlines():
 x=line.split()
 if len(x)==7 and x[0].isdigit() and start<int(x[0])<=end:
  v=list(map(float,x));assert all(math.isfinite(y) for y in v) and abs(v[1]-int(x[0])*.03*.024188843265857)<1e-5 and abs(v[5]-32)<1e-3
  physical.append(int(x[0]))
assert physical==list(range(start+1,end+1)),('Ne/progress coverage',len(physical))
def rows(name,cols=None):
 result=[]
 for line in (p/name).read_text().splitlines():
  if not line.strip() or line.lstrip().startswith(('#','!')):continue
  x=[float(y.replace('D','E').replace('d','e')) for y in line.split()]
  assert all(math.isfinite(y) for y in x) and (cols is None or len(x)==cols),name
  result.append(x)
 assert result,name
 return result
def window(name,cols):
 allrows=rows(name,cols);assert all(x[0]>=start*.03-1e-6 for x in allrows)
 result=[x for x in allrows if x[0]>start*.03+1e-8]
 assert len(result)==end-start,(name,len(result),end-start)
 assert all(abs(x[0]-(start+1+i)*.03)<1e-6 for i,x in enumerate(result))
 return result
rt=window('Si_rt.data',16);energy=window('Si_rt_energy.data',3)
pulse=check_pulse(rt,start=start)
nex=rows('Si_nex.data');ov=rows('Si_ovlp.data');assert len(ov)%81==0
frames=[]
for offset in range(0,len(ov),81):
 h=ov[offset];a=ov[offset+1:offset+81];assert len(h)==1 and h[0].is_integer()
 frames.append(int(h[0]));assert all(len(x)==33 and int(x[0])==i+1 for i,x in enumerate(a))
expected=list(range((start//256+1)*256,end+1,256))
if not incomplete and end%256:expected.append(end)
assert frames==expected and len(nex)==len(frames),(frames,expected)
assert all(abs(x[0]-i*.03)<1e-6 and len(x)>=3 for i,x in zip(frames,nex))
if (p/'variables.log').exists():assert not re.search(r'#namelist: .*status=\s*[1-9]',(p/'variables.log').read_text())
summary={'pass':True,'complete':not incomplete,'checkpoint_incomplete':incomplete,'start_step':start,'last_step':end,'target_step':target,'rows_rt':len(rt),'rows_energy_after_restart':len(energy),'child_exit0':cfg['mpi'],'SVD_records':svdtotal,'SVD_max_residual_unit_u_unit_v':svdmax,'pulse_field_maxabs_error':pulse,'projection_steps':frames,'remaining_steps':target-end}
if frames:
 summary.update(last_projection_step=frames[-1],last_projection_excited=nex[-1][1],last_projection_holes=nex[-1][2])
 (out/'band-k-occupations.json').write_text(json.dumps({'step':frames[-1],'rows':ov[-80:]},indent=2))
if incomplete:
 manifest=validate(p/('checkpoint_rt_%06d'%end),end)
 (out/'checkpoint-manifest.json').write_text(json.dumps(manifest,indent=2))
(out/'rt-gate.json').write_text(json.dumps(summary,indent=2))
(out/'memory-gate.json').write_text(json.dumps({'pass':True,'hosts':len(hosts),'node_max_sum_sampled_rank_peaks_KiB':max(nodepeaks),'scope':'sum of sampled individual rank peaks; not simultaneous or guaranteed instantaneous peak'},indent=2))
print('CONTINUATION_RT_AND_MEMORY_GATE_PASS',json.dumps(summary),flush=True)
if incomplete:raise SystemExit(86)
