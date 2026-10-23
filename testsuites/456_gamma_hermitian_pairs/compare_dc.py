#!/usr/bin/env python3
"""Small production DC comparison; run once per requested orbital MPI layout."""
import argparse,json,os,re,subprocess,time,hashlib,sys,shutil,signal,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',required=True);p.add_argument('--output',required=True);p.add_argument('--ranks',type=int,default=8);p.add_argument('--resume-analysis',action='store_true');a=p.parse_args()
source=Path(__file__).resolve().parents[2];root=Path(a.output).resolve();root.mkdir(exist_ok=a.resume_analysis,parents=True)
assert a.ranks in [2,4,8]
base=(source/'testsuites/422_H_dcdft_hse/inputfile').read_text()
base=base.replace('nproc_k=2','nproc_k=1').replace('nproc_ob=1','nproc_ob=%d'%(a.ranks//2))
base=base.replace('nproc_rgrid_tot=4,1,1','nproc_rgrid_tot=%s'%{2:'2,1,1',4:'4,1,1',8:'4,2,1'}[a.ranks])
base=base.replace('num_kgrid=1,2,1','num_kgrid=1,1,1').replace('nscf=200','nscf=4').replace(" xc='hse06'"," exx_pre_scf_threshold=0d0\n xc='hse06'")
results={}
for mode in ['reference','default']:
 d=root/mode;reuse=d.exists() and a.resume_analysis
 if reuse:
  assert all((d/('rank-%d.json'%i)).exists() for i in range(a.ranks))
  wall=max(json.loads((d/('rank-%d.json'%i)).read_text())['wall'] for i in range(a.ranks));returncode=0
 else:
  d=root/mode;d.mkdir();text=base
  if mode=='default':text=text.replace(" yn_hse_sr_pair_symmetry='N'\n",'')
  (d/'inputfile').write_text(text);(d/'H_rps.dat').write_bytes((source/'testsuites/pseudo/H_rps.dat').read_bytes())
  env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
  for key in list(env):
   if key.startswith('SALMON_'):del env[key]
  with (d/'stdout').open('x') as out,(d/'stderr').open('x') as err:
   t=time.monotonic();q=subprocess.Popen([shutil.which('mpiexec'),'-n',str(a.ranks),sys.executable,str(Path(__file__).with_name('rank_run.py').resolve()),str(Path(a.binary).resolve())],stdin=(d/'inputfile').open(),stdout=out,stderr=err,cwd=d,env=env,start_new_session=True)
   try:q.wait(timeout=240)
   except subprocess.TimeoutExpired:
    os.killpg(q.pid,signal.SIGTERM);q.wait();raise
   wall=time.monotonic()-t
  returncode=q.returncode
 records=[json.loads((d/('rank-%d.json'%i)).read_text()) for i in range(a.ranks)]
 assert all(x['exit']==0 for x in records),records
 t=''.join((d/('stdout.rank-%d'%i)).read_text() for i in range(a.ranks))
 e=(d/'stderr').read_text()+''.join((d/('stderr.rank-%d'%i)).read_text() for i in range(a.ranks))
 assert returncode==0 and 'end SALMON' in t,(mode,returncode,t[-2000:],e[-1000:])
 assert not re.search(r'(?i)error stop|segmentation|sigsegv',t+e)
 vals=[float(x.replace('D','E')) for x in re.findall(r'Total Energy\s*=\s*(\S+)',t)]
 dc=[(int(i),float(v.replace('D','E'))) for i,v in re.findall(r'DC #SCF\s*=\s*(\d+).*?diff\s*=\s*(\S+)',t)]
 pairs=[int(x) for x in re.findall(r'EXX_GAMMA_HERMITIAN local_pairs=(\d+)',t)]
 assert dc and vals
 assert bool(pairs)==(mode=='default'),(mode,pairs,t[-2000:])
 ex=[float(v) for v in re.findall(r'DC_HSE_CORE exchange Ha\s*=\s*(\S+)',t)]
 ne=[float(v) for v in re.findall(r'integral\(rho_tot\)\s*=\s*(\S+)',t)]
 eigfiles=list((d/'data_dcdft/total').glob('*eigen*'));assert len(eigfiles)==1
 eigen=[float(line.split()[-1]) for line in eigfiles[0].read_text().splitlines() if line and not line.startswith('#')]
 assert ex and ne and eigen and all(math.isfinite(v) for v in vals+ex+ne+eigen)
 assert max(abs(v-4) for v in ne)<1e-8
 results[mode]=dict(ranks=records,rank_wall_max_seconds=max(x['wall'] for x in records),launcher_wall_seconds=None if reuse else wall,exchange_Ha=ex,eigen_Ha=eigen,nelec=ne,energy_eV=vals,dc_residual=dc,pair_counts=pairs,stderr=e,stdout_sha256=hashlib.sha256(t.encode()).hexdigest())
x=results['reference'];y=results['default'];assert len(x['energy_eV'])==len(y['energy_eV'])
assert len(x['dc_residual'])==len(y['dc_residual'])
assert all(i[0]==j[0] and abs(i[1]-j[1])<1e-10 for i,j in zip(x['dc_residual'],y['dc_residual']))
diff=max(abs(i-j) for i,j in zip(x['energy_eV'],y['energy_eV']));assert diff<1e-7,diff
assert len(x['exchange_Ha'])==len(y['exchange_Ha']) and len(x['eigen_Ha'])==len(y['eigen_Ha'])
dex=max(abs(i-j) for i,j in zip(x['exchange_Ha'],y['exchange_Ha']));dev=max(abs(i-j) for i,j in zip(x['eigen_Ha'],y['eigen_Ha']))
assert dex<1e-10 and dev<1e-10,(dex,dev)
results.update(exchange_maxabs_Ha=dex,eigen_maxabs_Ha=dev,ranks=a.ranks,orbital_ranks=a.ranks//2,energy_maxabs_eV=diff,PASS=True,scope='4 SCF steps including initialization; not converged GS or robust performance benchmark')
dest=root/('comparison-refined.json' if a.resume_analysis else 'comparison.json')
with dest.open('x') as f:json.dump(results,f,indent=2)
print(json.dumps({k:v for k,v in results.items() if k not in ['reference','default']},indent=2))
