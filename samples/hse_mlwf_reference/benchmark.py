"""Accuracy-matched Si exchange benchmark; these orbitals are NOT HSE SCF."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,json,time,platform
import numpy as np
from exchange import *
root=Path(__file__).resolve().parents[2];out=root/'docs/results/si-hse-mlwf';out.mkdir(exist_ok=True)
sys.path.insert(0,str(root/'docs/results/si-time-wannier'))
from wannier import read_u,geometry

def prepare(case):
 r,c,k,b,L,h=geometry(root);old=root/'docs/results/si-time-wannier'
 if case=='initial':
  path=root/'calculations/si_tdcdft_k4/gs/data_for_restart';g=np.load(old/'mlwf_initial.npz')['gauge']
 elif case=='strong':
  path=Path('/private/tmp/salmon-si-time-wannier-dense/strong/checkpoint_rt_001600');g=np.load(old/'strong_mlwf.npz')['final_gauge']
 else:raise ValueError(case)
 u,_=read_u(path);v=u@g
 cube=np.array([v[ik].T.reshape(16,12,12,12,order='F') for ik in range(64)])
 start=time.perf_counter();w,twist=bloch_to_wannier(cube,k,h);conversion=time.perf_counter()-start
 return cube,w,k,twist,h,conversion

def main(case):
 cube,w,k,twist,h,conversion=prepare(case);kernel=ScreenedKernel(w.shape[-3:],h,.11);shifts=cell_shifts(4,12)
 print(case,'norm',np.sum(abs(w)**2,axis=(1,2,3))*h**3,flush=True)
 scores=pair_overlaps(w,shifts,h**3)
 print('pairs',[(x,int((scores>=x).sum())) for x in [0,1e-6,1e-5,1e-4,1e-3,.01,.05]],flush=True)
 kernel.convolve(w[0].conj()*w[0])
 metrics=dict(case=case,provenance='PZ GS' if case=='initial' else 'old polarization-TDCDFT pump snapshot at3.096 fs; not HSE',
  numpy=np.__version__,python=sys.version,platform=platform.platform(),threads=os.environ.get('OPENBLAS_NUM_THREADS'),
  conversion_seconds=conversion,mixing=.25,omega=.11,runs=[])
 start=time.perf_counter();full,stats=localized_exchange(w,shifts,kernel,0.);elapsed=time.perf_counter()-start
 e=closed_shell_energy(w,full,h**3);metrics['full']=dict(stats,wall_seconds=elapsed,exchange_energy_Ha=e,hse_exchange_Ha=.25*e)
 cache=Path('/private/tmp/salmon-hse-mlwf');cache.mkdir(exist_ok=True)
 np.savez_compressed(cache/f'{case}_reference.npz',wannier=w,full_action=full,k=k,twist=twist,spacing=h)
 print('full',metrics['full'],flush=True)
 start=time.perf_counter();direct=bloch_exchange(cube,k,h,.11);metrics['bloch_seconds']=time.perf_counter()-start
 back=wannier_to_bloch(full,k,h,twist)
 metrics['bloch_action_relative_error']=float(np.linalg.norm(back-direct)/np.linalg.norm(direct))
 print('bloch',metrics['bloch_seconds'],metrics['bloch_action_relative_error'],flush=True)
 for threshold in [1e-6,1e-5,1e-4,1e-3,.01,.05]:
  row=dict(pair_tolerance=threshold,timings=[])
  for repeat in range(3):
   tick=time.perf_counter();action,s=localized_exchange(w,shifts,kernel,threshold);row['timings'].append(dict(s,wall_seconds=time.perf_counter()-tick))
  ee=closed_shell_energy(w,action,h**3)
  row.update(exchange_energy_Ha=ee,hse_energy_error_meV_per_atom=(ee-e)*.25*27.211386245988*1000/8,
   action_relative_error=float(np.linalg.norm(action-full)/np.linalg.norm(full)),
   median_wall_seconds=float(np.median([x['wall_seconds'] for x in row['timings']])))
  metrics['runs'].append(row);print('threshold',row['pair_tolerance'],row['median_wall_seconds'],row['action_relative_error'],row['hse_energy_error_meV_per_atom'],flush=True)
  (out/f'{case}_benchmark.json').write_text(json.dumps(metrics,indent=2)+'\n')
 print('done',case,flush=True)
if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else 'initial')
