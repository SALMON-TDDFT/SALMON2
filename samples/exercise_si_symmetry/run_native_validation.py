"""Run prepared symmetry fixtures and compare currents; optional local MPI integration test."""
from pathlib import Path
import argparse,json,os,re,subprocess,time
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',required=True,type=Path);args=p.parse_args()
root=Path(__file__).resolve().parents[2];base=root/'calculations/si_hse_native/symmetry_validation'
runs={};exe=str(args.executable.resolve())
for name in ['hse_full','hse_reduced','td_full','td_reduced']:
 d=base/name;n=3 if 'reduced' in name else 4
 with (d/'inputfile').open() as inp,(d/'run.log').open('w') as log:
  r=subprocess.run(['mpiexec','--bind-to','none','-n',str(n),exe],cwd=d,stdin=inp,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),timeout=600)
 text=(d/'run.log').read_text();assert r.returncode==0 and 'end SALMON' in text,name
 runs[name]={'rt_seconds':float(re.search(r'^rt iterations\s+(.+)$',text,re.M)[1].split()[1]),'ranks':n}
 print(name,runs[name],flush=True)
for model in ['hse','td']:
 a=np.loadtxt(base/(model+'_full')/'Si_rt.data');b=np.loadtxt(base/(model+'_reduced')/'Si_rt.data')
 np.testing.assert_allclose(a[:,0],b[:,0],rtol=0,atol=1e-12)
 err=np.linalg.norm(a[:,15]-b[:,15])/np.linalg.norm(a[:,15]);assert err<1e-5,(model,err)
 runs[model+'_relative_current_error']=float(err)
(base/'parity.json').write_text(json.dumps(runs,indent=2)+'\n')
print(json.dumps(runs,indent=2))
