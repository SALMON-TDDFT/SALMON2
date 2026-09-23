"""Persistent MPI exchange workers with a rank-zero serial solver/controller.

Construct collectively. Non-root ranks call serve(); only root calls apply()
and close(). All ranks replicate source/target arrays but compute disjoint rows.
Python computation failures are reported collectively before reduction, allowing
root to checkpoint and close the service. Broken MPI communication aborts the
communicator rather than leaving orphaned workers in a collective.
"""
import time
import numpy as np
from mpi4py import MPI
from distance_exchange import DistanceExchange


class CollectiveExchangeError(RuntimeError):
 """All ranks observed a recoverable failure before reduction."""


class MPIExchange:
 def __init__(self,comm,configuration):
  self.comm=comm;self.rank=comm.Get_rank();self.size=comm.Get_size();self.closed=False
  self.operator=None;error=None
  try:self.operator=DistanceExchange(**configuration)
  except Exception as exc:error=f'rank {self.rank}: {type(exc).__name__}: {exc}'
  errors=comm.allgather(error)
  if any(errors):raise RuntimeError('MPI exchange setup failed: '+'; '.join(e for e in errors if e))
  self.radius=self.operator.radius

 def validate_for_hse(self):self.operator.validate_for_hse()

 def _perform(self,packet,sources=None,targets=None):
  start=time.perf_counter();_,source_shape,target_shape,same=packet
  if self.rank:sources=np.empty(source_shape,complex)
  self.comm.Bcast(sources,root=0)
  if same:targets=sources
  else:
   if self.rank:targets=np.empty(target_shape,complex)
   self.comm.Bcast(targets,root=0)
  error=None;partial=stats=result=None
  try:
   partial,stats=self.operator.apply(sources,targets,row_rank=self.rank,row_size=self.size)
   result=np.empty_like(partial) if self.rank==0 else None
  except BaseException as exc:error=f'rank {self.rank}: {type(exc).__name__}: {exc}'
  errors=self.comm.allgather(error)
  if any(errors):
   if self.rank==0:raise CollectiveExchangeError('MPI exchange computation failed: '+'; '.join(e for e in errors if e))
   return None
  # Only one rank owns each row. This sum combines disjoint contributions,
  # rather than changing the orbital or k reduction order within a row.
  self.comm.Reduce(partial,result,op=MPI.SUM,root=0)
  rank_stats=self.comm.gather(stats,root=0)
  if self.rank:return None
  entries=sum(s['kernel_entries'] for s in rank_stats)
  info=dict(exchange_seconds=time.perf_counter()-start,mpi_ranks=self.size,
   rank_exchange_seconds=[s['exchange_seconds'] for s in rank_stats],
   workspace_bytes_per_rank=[s['workspace_bytes'] for s in rank_stats],
   primitive_pair_fraction=sum(s['primitive_pair_fraction'] for s in rank_stats),
   nonzero_kernel_fraction=sum(s['nonzero_kernel_entries'] for s in rank_stats)/max(1,entries),
   radius_bohr=self.radius,block_rows=self.operator.block_rows)
  return result,info

 def apply(self,sources,targets):
  if self.rank!=0 or self.closed:raise RuntimeError('Only active rank0 may dispatch exchange')
  same=targets is sources
  sources=np.ascontiguousarray(sources,dtype=np.complex128)
  targets=sources if same else np.ascontiguousarray(targets,dtype=np.complex128)
  packet=('apply',sources.shape,targets.shape,same);start=time.perf_counter()
  try:
   self.comm.bcast(packet,root=0);result,stats=self._perform(packet,sources,targets)
  except CollectiveExchangeError:raise
  except BaseException:
   self.comm.Abort(1);raise
  stats['exchange_seconds']=time.perf_counter()-start
  return result,stats

 def serve(self):
  if self.rank==0:raise RuntimeError('Rank0 must dispatch, not serve')
  try:
   while True:
    packet=self.comm.bcast(None,root=0)
    if packet[0]=='stop':self.closed=True;return
    if packet[0]!='apply':raise RuntimeError('Unknown MPI exchange command')
    self._perform(packet)
  except BaseException:
   self.comm.Abort(1);raise

 def close(self):
  if self.rank!=0:raise RuntimeError('Only rank0 may close the service')
  if not self.closed:
   self.comm.bcast(('stop',),root=0);self.closed=True
