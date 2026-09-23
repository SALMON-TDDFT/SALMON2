"""MPI entry point: rank0 owns propagation/files, other ranks serve exchange."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse
from mpi4py import MPI
from mpi_exchange import MPIExchange


def configuration_from_export(comm,export):
 packet=None
 if comm.Get_rank()==0:
  try:
   from model import NativeModel
   m=NativeModel(export)
   packet=(dict(shape=m.shape,h=m.h,k=m.k,omega=.11,radius=None,block_rows=16),None)
  except Exception as error:packet=(None,f'{type(error).__name__}: {error}')
 configuration,error=comm.bcast(packet,root=0)
 if error:raise RuntimeError('MPI setup: '+error)
 return configuration


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('export','state','hse_directory','tdcdft_directory','output'):p.add_argument(name)
 p.add_argument('--target-steps',type=int,default=2750);p.add_argument('--dt',type=float,default=.32)
 p.add_argument('--amplitude',type=float,default=1e-4)
 p.add_argument('--propagation-only',action='store_true',help='Validation pilot: skip final comparison')
 a=p.parse_args();comm=MPI.COMM_WORLD
 service=MPIExchange(comm,configuration_from_export(comm,a.export))
 if comm.Get_rank():service.serve();return
 try:
  if a.propagation_only:
   from propagate_ptcn import run
   run(a.export,a.state,a.hse_directory,a.target_steps,dt=a.dt,amplitude=a.amplitude,
       resume=True,exchange_method='blocked',exchange_backend=service)
  else:
   from run_linear_comparison import execute
   values=vars(a);values.pop('propagation_only');execute(**values,exchange_backend=service)
 finally:service.close()


if __name__=='__main__':main()
