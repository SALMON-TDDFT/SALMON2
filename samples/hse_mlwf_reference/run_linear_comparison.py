"""Sequential resumable HSE propagation followed by matched TDCDFT analysis.

Suitable for a detached local process. No timer/scheduler or notification is
created. The propagator owns atomic physical checkpoints; job_status.json owns
the propagation/analysis stage and records any failure.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse
import fcntl
import json
import math
from pathlib import Path
from checkpoint import write_json,fingerprint,load_checkpoint
from propagate_ptcn import run
from compare_linear import analyze


def execute(export,state,hse_directory,tdcdft_directory,output,target_steps=2750,dt=.32,amplitude=1e-4,exchange_backend=None):
 hse=Path(hse_directory)
 with (hse/'.comparison.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  job=dict(pid=os.getpid(),target_steps=target_steps,end_time_au=target_steps*dt,
           hse_directory=str(hse.resolve()),comparison_directory=str(Path(output).resolve()),
           mpi_ranks=getattr(exchange_backend,'size',1))
  try:
   write_json(hse/'job_status.json',dict(job,status='running',stage='propagation'))
   saved=json.loads((hse/'status.json').read_text()) if (hse/'status.json').exists() else {}
   complete=saved.get('status')=='completed' and saved.get('accepted_step',0)>=target_steps
   if complete:
    if not math.isclose(saved['dt_au'],dt,rel_tol=0,abs_tol=1e-14) or not math.isclose(saved['amplitude'],amplitude,rel_tol=0,abs_tol=1e-14):
     raise ValueError('Completed trajectory physics differs from requested comparison')
    export_path=Path(export);state_path=Path(state)
    expected=dict(dt=dt,amplitude=amplitude,
     export_hash=fingerprint([*export_path.glob('*.bin'),export_path/'metadata.txt',export_path/'complete.txt']),
     initial_hash=fingerprint([state_path,state_path.parent/'result.json']))
    load_checkpoint(hse/'restart.npz',expected)
   else:
    run(export,state,hse_directory,target_steps,dt=dt,amplitude=amplitude,resume=True,exchange_method='blocked',exchange_backend=exchange_backend)
   write_json(hse/'job_status.json',dict(job,status='running',stage='analysis'))
   analyze(hse_directory,tdcdft_directory,output,end_time=target_steps*dt)
   write_json(hse/'job_status.json',dict(job,status='completed',stage='analysis'))
  except BaseException as error:
   try:write_json(hse/'job_status.json',dict(job,status='failed',error=f'{type(error).__name__}: {error}'))
   except Exception:pass
   raise


if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('export','state','hse_directory','tdcdft_directory','output'):p.add_argument(name)
 p.add_argument('--target-steps',type=int,default=2750);p.add_argument('--dt',type=float,default=.32)
 p.add_argument('--amplitude',type=float,default=1e-4)
 execute(**vars(p.parse_args()))
