"""Run with mpiexec: real-collective parity and worker-error recovery checks."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json
import numpy as np
from mpi4py import MPI
from distance_exchange import DistanceExchange
from mpi_exchange import MPIExchange


def main():
 comm=MPI.COMM_WORLD;rank=comm.Get_rank();size=comm.Get_size();h=.8;n=3;mesh=2
 axis=(np.arange(mesh)+.3-mesh/2)*2*np.pi/(mesh*n*h)
 k=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1).reshape(-1,3)
 config=dict(shape=(n,)*3,h=h,k=k,omega=.3,block_rows=4)
 service=MPIExchange(comm,config);original=service.operator.apply;calls=[0]
 def injected(*args,**kwargs):
  calls[0]+=1
  if rank==min(1,size-1) and calls[0]==3:raise RuntimeError('injected worker computation failure')
  return original(*args,**kwargs)
 service.operator.apply=injected
 if rank:
  service.serve();return
 try:
  rng=np.random.default_rng(91);u=(rng.normal(size=(8,2,3,3,3))+1j*rng.normal(size=(8,2,3,3,3)))*.05
  v=(rng.normal(size=(8,3,3,3,3))+1j*rng.normal(size=(8,3,3,3,3)))*.05
  serial=DistanceExchange(**config);errors=[]
  for target in (u,v):
   ref=serial.apply(u,target)[0];got,stats=service.apply(u,target)
   np.testing.assert_allclose(got,ref,atol=1e-14,rtol=1e-14)
   errors.append(float(np.linalg.norm(got-ref)/np.linalg.norm(ref)))
  try:service.apply(u,u)
  except RuntimeError as error:assert 'injected worker' in str(error)
  else:raise AssertionError('Worker error was not delivered to rank0')
  np.testing.assert_allclose(service.apply(u,u)[0],serial.apply(u,u)[0],atol=1e-14,rtol=1e-14)
  import mpi_exchange
  allocate=mpi_exchange.np.empty_like
  def fail_allocation(*a,**kw):raise MemoryError('injected root result allocation failure')
  mpi_exchange.np.empty_like=fail_allocation
  try:
   try:service.apply(u,u)
   except RuntimeError as error:assert 'injected root result allocation' in str(error)
   else:raise AssertionError('Allocation error was not coordinated')
  finally:mpi_exchange.np.empty_like=allocate
  np.testing.assert_allclose(service.apply(u,u)[0],serial.apply(u,u)[0],atol=1e-14,rtol=1e-14)
  print(json.dumps(dict(passed=True,ranks=size,relative_errors=errors,worker_error_recovered=True,allocation_error_recovered=True)),flush=True)
 finally:service.close()


if __name__=='__main__':main()
