"""Fixed-support HSE exchange accuracy/timing and ACE compatibility gate.

Usage: python benchmark_local_support.py EXPORT STATE OUTPUT.json
Works with SCF state.npz or an atomic RT restart.npz. No propagation is changed.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse
import hashlib
import io
import time
from pathlib import Path
import numpy as np
from model import NativeModel
from exchange import (ScreenedKernel,cell_shifts,bloch_to_wannier,wannier_to_bloch,
                      symmetric_exchange,box_origins,closed_shell_energy)
from local_support import LocalSupportExchange
from fftw_backend import FFTWConvolution
from ace import ACE
from checkpoint import write_json


def load_snapshot(path):
 # Atomic replacement of a live checkpoint cannot change this byte snapshot.
 raw=Path(path).read_bytes()
 with np.load(io.BytesIO(raw),allow_pickle=False) as data:
  return data['u'],data['gauge'],hashlib.sha256(raw).hexdigest()


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('export');parser.add_argument('state');parser.add_argument('output')
 parser.add_argument('--widths',type=int,nargs='+',default=[12,16,20,24])
 parser.add_argument('--repeats',type=int,default=2)
 args=parser.parse_args()
 if args.repeats<1:parser.error('repeats must be positive')
 m=NativeModel(args.export)
 u,g,state_digest=load_snapshot(args.state)
 cube=np.einsum('knxyz,knm->kmxyz',u,g,optimize=True)
 w,twist=bloch_to_wannier(cube,m.k,m.h);kernel=ScreenedKernel(w.shape[-3:],m.h,.11)
 shifts=cell_shifts(round(m.nk**(1/3)),m.shape[0])
 result=dict(state=str(Path(args.state).resolve()),state_sha256=state_digest,
             shape=list(w.shape),spacing_bohr=m.h,repeats=args.repeats,rows=[],
             scope='Fixed gauge and supports; variational orbital functional, not an RT implementation',
             timing_caveat='Single-thread FFTW/BLAS; full-support RT job may run concurrently')
 def metric(action):
  back=wannier_to_bloch(action,m.k,m.h,twist)
  ku=np.einsum('kmxyz,knm->knxyz',back,g.conj(),optimize=True)
  a=u.reshape(*u.shape[:2],-1);b=ku.reshape(a.shape)
  matrix=a.conj()@b.transpose(0,2,1)*m.dv
  defect=np.linalg.norm(matrix-matrix.conj().transpose(0,2,1),axis=(1,2))/np.linalg.norm(matrix,axis=(1,2))
  info=dict(max_relative_metric_antihermiticity=float(max(defect)))
  try:
   ace=ACE(u,ku,m.dv)
   info.update(ace_accepted=True,ace_span_relative_error=float(np.linalg.norm(ace.apply(u)-ku)/np.linalg.norm(ku)))
  except ValueError as error:info.update(ace_accepted=False,ace_error=str(error))
  return info
 with FFTWConvolution(kernel.multiplier) as backend:
  kernel.convolve=backend.convolve;backend.convolve(w[0].conj()*w[0]);times=[]
  for _ in range(args.repeats):
   tick=time.perf_counter();ref,stats=symmetric_exchange(w,shifts,kernel,0.);times.append(time.perf_counter()-tick)
  e0=closed_shell_energy(w,ref,m.dv)
  result['reference']=dict(seconds=times,median_seconds=float(np.median(times)),
                           hse_exchange_Ha=.25*e0,**metric(ref))
 for width in args.widths:
  origins=box_origins(w,width);tick=time.perf_counter()
  with LocalSupportExchange(kernel,shifts,origins,width,fftw=True) as op:
   setup=time.perf_counter()-tick
   tails=[float(np.sum(abs(w[i])**2)*m.dv-np.sum(abs(w[(i,)+ix])**2)*m.dv) for i,ix in enumerate(op.indices)]
   op.apply(w);times=[]
   for _ in range(args.repeats):
    tick=time.perf_counter();action,stats=op.apply(w);times.append(time.perf_counter()-tick)
  e=closed_shell_energy(w,action,m.dv)
  row=dict(side_bohr=width*m.h,setup_seconds=setup,seconds=times,
           median_seconds=float(np.median(times)),source_tail_norm_max=max(tails),
           action_relative_error=float(np.linalg.norm(action-ref)/np.linalg.norm(ref)),
           hse_exchange_error_meV_per_atom=(e-e0)*.25*27.211386245988*1000/8,
           **stats,**metric(action))
  result['rows'].append(row);write_json(Path(args.output),result)
  print(row,flush=True)
 write_json(Path(args.output),result)


if __name__=='__main__':main()
