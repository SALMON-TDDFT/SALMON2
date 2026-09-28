#!/usr/bin/env python3
"""Real Gamma LCFO -> complex mesh response regression using an independently generated PBE seed."""
import argparse,json,os,re,shutil,subprocess,tempfile
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--seed',type=Path,required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--mpiexec',default=os.environ.get('SALMON_TEST_MPIEXEC','mpiexec'));a=p.parse_args()
root=Path(tempfile.mkdtemp(prefix='lcfo-real-response-'));seed=a.seed.resolve();inp=a.input.read_text();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
base=re.sub(r'\bnt\s*=\s*\d+','nt=1',inp)
def run(name,ranks=4,mutation=None):
 d=root/name;d.mkdir();shutil.copy2(seed/'H_rps.dat',d/'H_rps.dat');shutil.copytree(seed/'data_dcdft',d/'data_dcdft')
 s=base.replace('nproc_rgrid=1,2,2',f'nproc_rgrid={ranks//4},2,2').replace('nproc_rgrid_tot=1,2,2',f'nproc_rgrid_tot={ranks//4},2,2');(d/'inputfile').write_text(s)
 if mutation:
  q=d/'data_dcdft/fragments/000002/basis_functions.bin'
  if mutation=='mixed':q.write_bytes(b'SLCFO_COMPLEX_V1 '+q.read_bytes()[16:])
  elif mutation=='truncated':q.write_bytes(b'bad')
 with (d/'output').open('wb') as out:
  result=subprocess.run([a.mpiexec,'--bind-to','none','-n',str(ranks),str(a.binary.resolve())],cwd=d,env=env,input=s.encode(),stdout=out,stderr=subprocess.STDOUT,timeout=120)
 log=(d/'output').read_text()
 if mutation:
  assert 'DC-LCFO reconstruction: invalid or mixed file headers' in log,(d,log[-1000:]);assert 'end SALMON' not in log
 else:
  assert result.returncode==0 and 'end SALMON' in log,(d,log[-1000:])
  v=np.loadtxt(d/'h2_rt.data',ndmin=2);assert np.isfinite(v).all();return v
x=run('mpi4');y=run('mpi8',8);np.testing.assert_allclose(x,y,rtol=0,atol=1e-10)
run('mixed',mutation='mixed');run('truncated',mutation='truncated')
print(json.dumps({'passed':4,'work':str(root),'current_difference':float(np.max(abs(x[:,13:]-y[:,13:])))}))
