"""Communication-inclusive MPI exchange timing on a fixed Si snapshot."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,time,json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from mpi_exchange import MPIExchange
from mpi_run import configuration_from_export
from distance_exchange import DistanceExchange
from benchmark_local_support import load_snapshot
from checkpoint import write_json


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('export');p.add_argument('state');p.add_argument('output')
 p.add_argument('--repeats',type=int,default=3);a=p.parse_args()
 if a.repeats<1:p.error('Positive repeats required')
 comm=MPI.COMM_WORLD;configuration=configuration_from_export(comm,a.export);service=MPIExchange(comm,configuration)
 if comm.Get_rank():service.serve();return
 try:
  u,_,digest=load_snapshot(a.state);serial=DistanceExchange(**configuration)
  tick=time.perf_counter();reference,_=serial.apply(u,u);serial_seconds=time.perf_counter()-tick
  service.apply(u,u);times=[];runs=[]
  for _ in range(a.repeats):
   tick=time.perf_counter();result,stats=service.apply(u,u);times.append(time.perf_counter()-tick);runs.append(stats)
  error=float(np.linalg.norm(result-reference)/np.linalg.norm(reference))
  if error>1e-12:raise RuntimeError(f'MPI exchange differs from serial: {error}')
  report=dict(mpi_ranks=comm.Get_size(),mpi_library=MPI.Get_library_version(),state_sha256=digest,
    serial_seconds=serial_seconds,mpi_seconds=times,mpi_median_seconds=float(np.median(times)),
    speedup=serial_seconds/float(np.median(times)),relative_error=error,runs=runs,
    scope='Exchange only; communication included; one serial reference timing, three MPI repetitions; serial production job remains active')
  write_json(Path(a.output),report);print(json.dumps(report),flush=True)
 finally:service.close()


if __name__=='__main__':main()
