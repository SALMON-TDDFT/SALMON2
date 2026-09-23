"""Fixed-support variational MLWF exchange, for accuracy experiments.

For v_i=P_i w_i, returns P_i K[v]v_i, the derivative of E_x[v]/4
in the real directional-derivative convention. P_i and the MLWF gauge must
remain fixed while taking that derivative. Orbital-specific projectors break
occupied-unitary invariance: this is NOT a common Hermitian operator on trial
vectors, and must not be silently passed to ACE or the production propagator.
"""
import time
import numpy as np


class LocalSupportExchange:
 def __init__(self,kernel,shifts,origins,width,fftw=False):
  shape=np.asarray(kernel.shape)
  if len(shape)!=3 or len(set(shape))!=1:
   raise ValueError('Cubic periodic grid required')
  if not np.isfinite(width) or int(width)!=width or not 1<=width<=min(shape):
   raise ValueError('Integer support width within grid required')
  shifts=np.asarray(shifts);origins=np.asarray(origins)
  for a in (shifts,origins):
   if a.ndim!=2 or a.shape[1]!=3 or not len(a) or not np.isfinite(a).all() or not np.equal(a,np.floor(a)).all():
    raise ValueError('Nonempty integer coordinate triples required')
  self.shape=tuple(shape);self.width=int(width);self.kernel=kernel
  self.origins=origins.astype(int)%shape;shifts=shifts.astype(int)%shape
  lookup={tuple(s):i for i,s in enumerate(shifts)}
  if len(lookup)!=len(shifts) or (0,0,0) not in lookup:
   raise ValueError('Unique translation group including zero required')
  for s in shifts:
   if any(tuple((s+t)%shape) not in lookup for t in shifts):
    raise ValueError('Translations must form a periodic group')
  self.inverse=[lookup[tuple((-s)%shape)] for s in shifts]
  self.shifts=shifts;self.backend=None;self.closed=False
  self.indices=[np.ix_(*[(np.arange(self.width)+s)%n for s,n in zip(o,shape)]) for o in self.origins]
  # The source is confined to one box. Only its overlap with the target
  # support contributes to the pair density and both gradient terms.
  self.pairs=[]
  for i,origin in enumerate(self.origins):
   for j in range(i,len(self.origins)):
    for ir,shift in enumerate(shifts):
     if i==j and ir>self.inverse[ir]:continue
     axes=[(np.arange(self.width)+s+d)%n for s,d,n in zip(origin,shift,shape)]
     keep=[((a-o)%n)<self.width for a,o,n in zip(axes,self.origins[j],shape)]
     if not all(np.any(a) for a in keep):continue
     mask=keep[0][:,None,None]*keep[1][None,:,None]*keep[2][None,None,:]
     self.pairs.append((i,j,ir,np.ix_(*axes),mask,i!=j or ir!=self.inverse[ir]))
  # Full-periodic FFT is cheaper once zero padding would exceed the cell.
  self.pad=min(2*self.width,shape[0]);p=self.pad
  offsets=np.arange(p)
  if p<shape[0]:offsets=np.where(offsets<self.width,offsets,offsets-p)
  real=np.fft.ifftn(kernel.multiplier).real
  self.multiplier=np.fft.fftn(real[np.ix_(*[offsets%n for n in shape])])
  self.plan_seconds=0.
  if fftw:
   from fftw_backend import FFTWConvolution
   self.backend=FFTWConvolution(self.multiplier);self.plan_seconds=self.backend.plan_seconds

 def apply(self,w):
  if self.closed:raise RuntimeError('Local support exchange is closed')
  w=np.asarray(w)
  if w.shape!=(len(self.origins),)+self.shape or not np.isfinite(w).all():
   raise ValueError('Finite orbitals matching supports and grid required')
  start=time.perf_counter();out=np.zeros_like(w,dtype=complex)
  local=[w[(i,)+ix] for i,ix in enumerate(self.indices)]
  box=(slice(0,self.width),)*3
  padded=np.zeros((self.pad,)*3,complex)
  for i,j,ir,ix,mask,reciprocal in self.pairs:
   source=local[i];target=w[(j,)+ix]*mask
   padded.fill(0);padded[box]=source.conj()*target
   if self.backend is None:potential=np.fft.ifftn(np.fft.fftn(padded)*self.multiplier)[box]
   else:potential=self.backend.convolve(padded)[box]
   out[(j,)+ix]-=source*potential*mask
   if reciprocal:out[(i,)+self.indices[i]]-=target*potential.conj()
  return out,dict(convolutions=len(self.pairs),possible_directed_pairs=len(w)**2*len(self.shifts),
                 width=self.width,fft_width=int(self.pad),exchange_seconds=time.perf_counter()-start)

 def close(self):
  if self.backend is not None:self.backend.close()
  self.closed=True

 def __enter__(self):return self
 def __exit__(self,*exc):self.close()
