"""Source-domain truncation and selected-pair timings at fixed grid kernel."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,json,time
import numpy as np
from exchange import *
root=Path(__file__).resolve().parents[2];out=root/'docs/results/si-hse-mlwf'
case=sys.argv[1] if len(sys.argv)>1 else 'initial'
f=np.load(Path('/private/tmp/salmon-hse-mlwf')/f'{case}_reference.npz');w=f['wannier'];ref=f['full_action'];h=float(f['spacing']);dv=h**3
kernel=ScreenedKernel(w.shape[-3:],h,.11);shifts=cell_shifts(4,12);scores=pair_overlaps(w,shifts,dv)
e0=closed_shell_energy(w,ref,dv);rows=[]
for width in [12,16,20,24,32,48]:
 origins=box_origins(w,width)
 discarded=[]
 for i,origin in enumerate(origins):
  ix=np.ix_(*[(np.arange(width)+s)%n for s,n in zip(origin,kernel.shape)])
  discarded.append(float(1-np.sum(abs(w[(i,)+ix])**2)*dv))
 for threshold in [0.,1e-6,1e-5,1e-4,1e-3]:
  timings=[];selected=scores>=threshold
  for repeat in range(3):
   tick=time.perf_counter();a=source_box_exchange(w,w,shifts,origins,width,kernel,selected);timings.append(time.perf_counter()-tick)
  ee=closed_shell_energy(w,a,dv)
  row=dict(width=width,pair_threshold=threshold,retained_pairs=int(selected.sum()),timings=timings,median_seconds=float(np.median(timings)),
    source_tail_norm_max=max(discarded),action_relative_error=float(np.linalg.norm(a-ref)/np.linalg.norm(ref)),
    hse_exchange_trace_error_meV_per_atom=(ee-e0)*.25*27.211386245988*1000/8)
  rows.append(row);print(row,flush=True)
  (out/f'{case}_local_benchmark.json').write_text(json.dumps(rows,indent=2)+'\n')
