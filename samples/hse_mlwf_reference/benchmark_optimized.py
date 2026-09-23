"""Repeated equal-accuracy timings of reciprocal-pair reuse and MLWF screening."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,json,time,platform,resource
import numpy as np
from exchange import *
from fftw_backend import FFTWConvolution
from benchmark import prepare
root=Path(__file__).resolve().parents[2];out=root/'docs/results/si-hse-mlwf'

def main(case):
 if case=='hse':
  from model import NativeModel
  m=NativeModel('/private/tmp/salmon-hse-export-check/export');d=np.load('/private/tmp/salmon-hse-mlwf/scf/state.npz')
  cube=np.einsum('knxyz,knm->kmxyz',d['u'],d['gauge'],optimize=True);k=m.k;h=m.h
  tick=time.perf_counter();w,twist=bloch_to_wannier(cube,k,h);conversion=time.perf_counter()-tick
 else:cube,w,k,twist,h,conversion=prepare(case)
 kernel=ScreenedKernel(w.shape[-3:],h,.11);shifts=cell_shifts(4,12)
 result=dict(case=case,implementation='Python NumPy + single-thread FFTW, reciprocal pair reuse',conversion_seconds=conversion,
  shape=list(w.shape),orbital_bytes=w.nbytes,platform=platform.platform(),numpy=np.__version__,threads=os.environ['OPENBLAS_NUM_THREADS'],rows=[])
 with FFTWConvolution(kernel.multiplier) as backend:
  kernel.convolve=backend.convolve;result['fftw_plan_seconds']=backend.plan_seconds
  backend.convolve(w[0].conj()*w[0])
  for cut in [0.,1e-7,1e-6,1e-5,1e-4,1e-3]:
   repeats=[]
   for repeat in range(3):
    tick=time.perf_counter();action,stats=symmetric_exchange(w,shifts,kernel,cut);repeats.append(dict(stats,wall_seconds=time.perf_counter()-tick))
   e=closed_shell_energy(w,action,h**3)
   if cut==0:reference=action.copy();e0=e
   row=dict(pair_tolerance=cut,runs=repeats,median_seconds=float(np.median([x['wall_seconds'] for x in repeats])),
    action_relative_error=float(np.linalg.norm(action-reference)/np.linalg.norm(reference)),
    hse_energy_error_meV_per_atom=(e-e0)*.25*27.211386245988*1000/8,hse_exchange_energy_Ha=.25*e)
   result['rows'].append(row);print(case,cut,row['median_seconds'],row['action_relative_error'],row['hse_energy_error_meV_per_atom'],flush=True)
   (out/f'{case}_optimized_benchmark.json').write_text(json.dumps(result,indent=2)+'\n')
 back=wannier_to_bloch(reference,k,h,twist)
 tick=time.perf_counter();direct=bloch_exchange(cube,k,h,.11);result['bloch_reference_seconds']=time.perf_counter()-tick
 result['bloch_action_relative_error']=float(np.linalg.norm(back-direct)/np.linalg.norm(direct))
 result['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
 (out/f'{case}_optimized_benchmark.json').write_text(json.dumps(result,indent=2)+'\n')
 print('done',case,flush=True)
if __name__=='__main__':main(sys.argv[1])
