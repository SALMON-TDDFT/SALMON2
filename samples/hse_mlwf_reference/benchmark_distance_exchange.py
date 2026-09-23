"""Si fixed-state distance-kernel exchange and ACE amortization measurements."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse
import time
from pathlib import Path
import numpy as np
from distance_exchange import DistanceExchange
from benchmark_local_support import load_snapshot
from model import NativeModel
from exchange import ScreenedKernel,cell_shifts,bloch_to_wannier,wannier_to_bloch,symmetric_exchange
from fftw_backend import FFTWConvolution
from ace import ACE
from checkpoint import write_json


def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('export');p.add_argument('state');p.add_argument('output')
 p.add_argument('--radii',type=float,nargs='+',default=[6.,10.,16.])
 p.add_argument('--repeats',type=int,default=2);p.add_argument('--block-rows',type=int,default=16)
 a=p.parse_args()
 if a.repeats<1:p.error('repeats must be positive')
 m=NativeModel(a.export);u,g,digest=load_snapshot(a.state)
 cube=np.einsum('knxyz,knm->kmxyz',u,g,optimize=True);w,twist=bloch_to_wannier(cube,m.k,m.h)
 kernel=ScreenedKernel(w.shape[-3:],m.h,.11);shifts=cell_shifts(round(m.nk**(1/3)),m.shape[0])
 result=dict(state_sha256=digest,nk=m.nk,primitive_shape=list(m.shape),occupied=u.shape[1],
             timing_scope='One BLAS/FFTW thread, two methods on identical state. Reference RT may run concurrently.',rows=[])
 with FFTWConvolution(kernel.multiplier) as backend:
  kernel.convolve=backend.convolve;backend.convolve(w[0].conj()*w[0]);times=[]
  for _ in range(a.repeats):
   tick=time.perf_counter();action,_=symmetric_exchange(w,shifts,kernel,0.);times.append(time.perf_counter()-tick)
  back=wannier_to_bloch(action,m.k,m.h,twist)
  ref=np.einsum('kmxyz,knm->knxyz',back,g.conj(),optimize=True)
 result['mlwf_reference']=dict(seconds=times,median_seconds=float(np.median(times)),
                              transform_and_localization_included=False)
 e0=.5*m.expectation(u,ref)
 for radius in [None]+a.radii:
  tick=time.perf_counter();op=DistanceExchange(m.shape,m.h,m.k,radius=radius,block_rows=a.block_rows)
  setup=time.perf_counter()-tick;op.apply(u,u);times=[]
  for _ in range(a.repeats):
   tick=time.perf_counter();action,stats=op.apply(u,u);times.append(time.perf_counter()-tick)
  energy=.5*m.expectation(u,action);flat=u.reshape(*u.shape[:2],-1);af=action.reshape(flat.shape)
  metric=flat.conj()@af.transpose(0,2,1)*m.dv
  defect=np.max(np.linalg.norm(metric-metric.conj().transpose(0,2,1),axis=(1,2))/np.linalg.norm(metric,axis=(1,2)))
  row=dict(stats,setup_seconds=setup,seconds=times,median_seconds=float(np.median(times)),
   action_relative_error=float(np.linalg.norm(action-ref)/np.linalg.norm(ref)),
   hse_exchange_error_meV_per_atom=(energy-e0)*.25*27.211386245988*1000/8,
   metric_antihermiticity=float(defect),minimum_kernel_multiplier=float(op.multiplier.min()))
  try:
   tick=time.perf_counter();ace=ACE(u,action,m.dv);compression=time.perf_counter()-tick
   ace.apply(u);applications=[]
   for _ in range(5):
    tick=time.perf_counter();ac=ace.apply(u);applications.append(time.perf_counter()-tick)
   at=float(np.median(applications));cost=row['median_seconds']
   row.update(ace_accepted=True,ace_compression_seconds=compression,ace_apply_seconds=at,
              ace_span_relative_error=float(np.linalg.norm(ac-action)/np.linalg.norm(action)),
              ace_condition_max=ace.condition_max,
              break_even_uses=(cost+compression)/(cost-at) if cost>at else None)
  except ValueError as error:row.update(ace_accepted=False,ace_error=str(error))
  result['rows'].append(row);write_json(Path(a.output),result);print(row,flush=True)


if __name__=='__main__':main()
