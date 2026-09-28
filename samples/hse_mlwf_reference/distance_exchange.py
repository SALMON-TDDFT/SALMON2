"""Common periodic exchange kernel with an optional relative-distance cutoff.

Stream primitive-cell density-matrix blocks; k FFTs replace explicit k-pair
sums. Sources are not truncated or gauge-localized. The real-space kernel is
the inverse FFT of the same sampled HSE multiplier used by the reference,
including its grid weight and finite-grid aliasing. A finite radius changes
that discrete kernel, not just the algorithm. The radius is minimum-image
electron separation in bohr, NOT distance between Wannier centers.
"""
import time
import numpy as np
from exchange import ScreenedKernel,mesh_indices,cell_shifts


class DistanceExchange:
 def __init__(self,shape,h,k,omega=.11,radius=None,block_rows=16):
  shape=tuple(shape)
  if len(shape)!=3 or len(set(shape))!=1 or any(int(n)!=n or n<1 for n in shape):
   raise ValueError('Positive cubic integer primitive grid required')
  if not np.isfinite(h) or h<=0 or not np.isfinite(omega) or omega<=0:
   raise ValueError('Positive finite spacing and screening required')
  if radius is not None and (not np.isfinite(radius) or radius<=0):
   raise ValueError('Positive finite cutoff radius required')
  if not np.isfinite(block_rows) or int(block_rows)!=block_rows or block_rows<1:
   raise ValueError('Positive integer row block required')
  self.shape=tuple(int(n) for n in shape);self.h=h;self.radius=radius;self.block_rows=int(block_rows)
  self.k=np.array(k,float,copy=True);self.nk=len(self.k);n=self.shape[0]
  indices,self.mesh=mesh_indices(self.k,n*h)
  self.order=np.argsort(np.ravel_multi_index(indices.T,(self.mesh,)*3))
  self.inverse_order=np.argsort(self.order)
  self.points=np.stack(np.meshgrid(*[np.arange(n)]*3,indexing='ij'),axis=-1).reshape(-1,3)
  self.phase=np.exp(1j*((self.k-self.k[0])@(self.points*h).T))
  self.shifts=cell_shifts(self.mesh,n);self.super_n=self.mesh*n
  kernel=ScreenedKernel((self.super_n,)*3,h,omega)
  self.real_kernel=np.fft.ifftn(kernel.multiplier).real
  if radius is not None:
   x=np.arange(self.super_n);x=np.minimum(x,self.super_n-x)*h
   r2=x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2
   self.real_kernel*=r2<=radius**2
  # Real/inversion-even kernel guarantees Hermiticity. ACE additionally checks
  # the occupied metric sign; do not assume arbitrary cutoffs preserve it.
  self.multiplier=np.fft.fftn(self.real_kernel).real

 def validate_for_hse(self):
  # Diagnostic apply() permits indefinite kernels; the HSE integration does not.
  if np.min(self.multiplier)<-1e-12*np.max(abs(self.multiplier)):
   raise ValueError('HSE distance kernel must be positive semidefinite on the active grid')

 def apply(self,sources,targets,row_rank=0,row_size=1):
  if not isinstance(row_rank,(int,np.integer)) or not isinstance(row_size,(int,np.integer)) or row_size<1 or not 0<=row_rank<row_size:
   raise ValueError('Valid integer row rank and size required')
  start=time.perf_counter();sources=np.asarray(sources);targets=np.asarray(targets)
  for a in (sources,targets):
   if a.ndim!=5 or a.shape[0]!=self.nk or a.shape[2:]!=self.shape or a.shape[1]<1 or not np.isfinite(a).all():
    raise ValueError('Finite k/orbital/grid arrays matching the operator required')
  s=sources.reshape(self.nk,sources.shape[1],-1)*self.phase[:,None,:]
  t=targets.reshape(self.nk,targets.shape[1],-1)*self.phase[:,None,:]
  out=np.zeros_like(t,dtype=complex);ng=len(self.points);workspace=0;columns_total=0
  formation=transform=application=0.;nonzero=total_entries=0
  for lo in range(row_rank*self.block_rows,ng,row_size*self.block_rows):
   hi=min(lo+self.block_rows,ng);points=self.points[lo:hi]
   delta=points[:,None,:]-self.points[None,:,:]
   columns=np.arange(ng)
   if self.radius is not None:
    shortest=(delta+self.shape[0]//2)%self.shape[0]-self.shape[0]//2
    columns=np.flatnonzero(np.any(np.sum(shortest**2,axis=-1)*self.h**2<=self.radius**2,axis=0))
   delta=delta[:,columns];columns_total+=(hi-lo)*len(columns)
   tick=time.perf_counter()
   density=s[:,:,lo:hi].transpose(0,2,1)@s[:,:,columns].conj()
   formation+=time.perf_counter()-tick
   tick=time.perf_counter()
   # fft(C)/Nk is the density at translated source cells; inverse transformation
   # to the Bloch operator supplies Nk, so these normalizations cancel.
   grid=density[self.order].reshape((self.mesh,)*3+(hi-lo,len(columns)))
   translated=np.fft.fftn(grid,axes=(0,1,2))
   offsets=(delta[None]-self.shifts[:,None,None,:])%self.super_n
   weights=self.real_kernel[tuple(np.moveaxis(offsets,-1,0))]
   nonzero+=int(np.count_nonzero(weights));total_entries+=weights.size
   translated*=weights.reshape(grid.shape)
   matrix=-np.fft.ifftn(translated,axes=(0,1,2)).reshape(density.shape)[self.inverse_order]
   transform+=time.perf_counter()-tick
   tick=time.perf_counter()
   out[:,:,lo:hi]=(matrix@t[:,:,columns].transpose(0,2,1)).transpose(0,2,1)
   application+=time.perf_counter()-tick
   # Explicit arrays only; NumPy/BLAS internal work buffers are not included.
   workspace=max(workspace,sum(a.nbytes for a in (density,grid,translated,offsets,weights,matrix)))
  out*=self.phase.conj()[:,None,:]
  return out.reshape(targets.shape),dict(exchange_seconds=time.perf_counter()-start,
   density_seconds=formation,transform_seconds=transform,application_seconds=application,
   workspace_bytes=int(workspace+s.nbytes+t.nbytes+out.nbytes),
   primitive_pair_fraction=columns_total/ng**2,nonzero_kernel_fraction=nonzero/max(1,total_entries),
   kernel_entries=total_entries,nonzero_kernel_entries=nonzero,
   radius_bohr=self.radius,block_rows=self.block_rows,row_rank=int(row_rank),row_size=int(row_size))
