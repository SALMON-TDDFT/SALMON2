"""Periodic short-range Fock reference and translation-aware MLWF pair screening.

One-spin occupied sources have unit occupation; closed-shell Ex=Tr(W^dagger K W).
No HSE mixing factor is included here. Pair-screened occupied actions are NOT an
operator on arbitrary trial vectors. The generic apply_exchange remains linear
in targets for a fixed occupied source set.
"""
import itertools
import time
import numpy as np
from numpy.fft import fftn, ifftn


def screened_multiplier(shape,h,omega,shift=None):
 if omega<=0 or h<=0:raise ValueError('Positive omega and spacing required')
 axes=[2*np.pi*np.fft.fftfreq(n,d=h) for n in shape]
 q=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
 if shift is not None:
  # Fold into the same sampled supercell Nyquist interval as Wannier reference.
  q=(q+np.asarray(shift)+np.pi/h)%(2*np.pi/h)-np.pi/h
 q2=np.sum(q*q,axis=-1);v=np.empty(shape);nonzero=q2>1e-24
 v[nonzero]=4*np.pi*(-np.expm1(-q2[nonzero]/(4*omega**2)))/q2[nonzero]
 v[~nonzero]=np.pi/omega**2
 return v


class ScreenedKernel:
 def __init__(self,shape,h,omega):
  self.shape=tuple(shape);self.h=float(h);self.dv=h**3;self.omega=omega
  self.multiplier=screened_multiplier(shape,h,omega)
 def convolve(self,rho):
  return ifftn(fftn(rho,axes=(-3,-2,-1))*self.multiplier,axes=(-3,-2,-1))


def apply_exchange(sources,targets,kernel):
 sources=np.asarray(sources);targets=np.asarray(targets)
 if sources.shape[1:]!=kernel.shape or targets.shape[1:]!=kernel.shape:raise ValueError('Grid shape differs')
 out=np.zeros_like(targets,dtype=complex)
 for w in sources:
  for j,phi in enumerate(targets):out[j]-=w*kernel.convolve(w.conj()*phi)
 return out


def closed_shell_energy(orbitals,action,dv):
 return float(np.vdot(orbitals,action).real*dv)


def cell_shifts(mesh,primitive_grid):
 return np.array(list(itertools.product(range(mesh),repeat=3)),dtype=int)*primitive_grid


def pair_overlaps(w,shifts,dv):
 """S_ijR=int |w_i(r-R)|² |w_j(r)|² / sqrt(IPR_i IPR_j)."""
 density=abs(w)**2;f=fftn(density,axes=(-3,-2,-1))
 ipr=np.sum(density*density,axis=(1,2,3))*dv
 result=np.empty((len(w),len(w),len(shifts)))
 indices=tuple(shifts.T)
 for i in range(len(w)):
  for j in range(len(w)):
   c=ifftn(f[i].conj()*f[j]).real
   result[i,j]=np.maximum(c[indices]*dv,0)/np.sqrt(ipr[i]*ipr[j])
 return result


def localized_exchange(w,shifts,kernel,pair_tolerance=0.):
 """Occupied action on home-cell Wanniers. Full periodic pair convolution.

Pair screening is symmetric under i,j,R -> j,i,-R when all cell shifts exist.
Spatial convolution optimization is measured separately; this path truncates
pairs only and retains the common reciprocal-grid kernel exactly.
"""
 if pair_tolerance<0:raise ValueError('Negative threshold')
 start=time.perf_counter();scores=pair_overlaps(w,shifts,kernel.dv)
 selected=scores>=pair_tolerance;screen_time=time.perf_counter()-start
 out=np.zeros_like(w,dtype=complex);start=time.perf_counter();convolution_time=0.
 for i in range(len(w)):
  for ir,s in enumerate(shifts):
   js=np.flatnonzero(selected[i,:,ir])
   if not len(js):continue
   source=np.roll(w[i],tuple(s),(0,1,2))
   for j in js:
    rho=source.conj()*w[j];tick=time.perf_counter();v=kernel.convolve(rho)
    convolution_time+=time.perf_counter()-tick;out[j]-=source*v
 return out,dict(total_pairs=int(selected.size),retained_pairs=int(selected.sum()),pair_tolerance=pair_tolerance,
  screening_seconds=screen_time,exchange_seconds=time.perf_counter()-start,convolution_seconds=convolution_time,
  grid_points=int(np.prod(kernel.shape)),wavefunction_bytes=int(w.nbytes))


def mesh_indices(k,primitive_length):
 k=np.asarray(k);mesh=round(len(k)**(1/3))
 if mesh**3!=len(k) or k.shape!=(len(k),3) or not np.isfinite(k).all():raise ValueError('Full cubic k mesh required')
 scaled=(k-k[0])*primitive_length*mesh/(2*np.pi);rounded=np.rint(scaled)
 if np.max(abs(scaled-rounded))>1e-8:raise ValueError('Nonuniform k mesh')
 indices=rounded.astype(int)%mesh
 if len(np.unique(indices,axis=0))!=len(k):raise ValueError('Duplicate/missing k points')
 return indices,mesh


def bloch_to_wannier(u,k,h):
 nk,no,n,ny,nz=u.shape
 if n!=ny or n!=nz:raise ValueError('Cubic primitive grid required')
 indices,mesh=mesh_indices(k,n*h);twist=np.asarray(k[0])
 r=np.stack(np.meshgrid(*(np.arange(n)*h for _ in range(3)),indexing='ij'),axis=-1)
 coefficients=np.empty((mesh,mesh,mesh,no,n,n,n),complex)
 for ik,index in enumerate(indices):coefficients[tuple(index)]=u[ik]*np.exp(1j*np.einsum('xyza,a->xyz',r,k[ik]-twist))
 cells=ifftn(coefficients,axes=(0,1,2))
 return cells.transpose(3,0,4,1,5,2,6).reshape(no,n*mesh,n*mesh,n*mesh),twist


def wannier_to_bloch(w,k,h,twist):
 nk=len(k);mesh=round(nk**(1/3));n=w.shape[-1]//mesh;no=len(w)
 indices,checked=mesh_indices(k,n*h)
 if not np.allclose(twist,k[0],atol=1e-12,rtol=0):raise ValueError('Twist differs from mesh origin')
 r=np.stack(np.meshgrid(*(np.arange(n)*h for _ in range(3)),indexing='ij'),axis=-1)
 cells=w.reshape(no,mesh,n,mesh,n,mesh,n).transpose(1,3,5,0,2,4,6)
 coefficients=fftn(cells,axes=(0,1,2));out=np.empty((nk,no,n,n,n),complex)
 for ik,index in enumerate(indices):out[ik]=coefficients[tuple(index)]*np.exp(-1j*np.einsum('xyza,a->xyz',r,k[ik]-twist))
 return out


def bloch_exchange(u,k,h,omega):
 """Independent primitive-cell occupied exchange, equal k weights."""
 nk,no=u.shape[:2];out=np.zeros_like(u)
 for ik,kv in enumerate(k):
  for iq,qv in enumerate(k):
   multiplier=screened_multiplier(u.shape[-3:],h,omega,kv-qv)
   for source in u[iq]:
    pairs=source.conj()[None,...]*u[ik]
    potential=ifftn(fftn(pairs,axes=(-3,-2,-1))*multiplier,axes=(-3,-2,-1))
    out[ik]-=source[None,...]*potential/nk
 return out

class LocalConvolution:
 """Exact restricted block of the SAME periodic discrete convolution.

The source and requested output lie in one box. Zero padding embeds its
Toeplitz block in a circulant FFT; no isolated/minimum-image kernel is substituted.
"""
 def __init__(self,kernel,width):
  self.width=int(width)
  if width<1 or width>min(kernel.shape):raise ValueError('Invalid source width')
  self.full=width==kernel.shape[0] and len(set(kernel.shape))==1
  self.kernel=kernel
  if self.full:return
  self.pad=2*width
  real_kernel=ifftn(kernel.multiplier).real
  offsets=np.arange(self.pad);offsets=np.where(offsets<width,offsets,offsets-self.pad)
  # The unused offset -width can be assigned consistently; never sampled by
  # output/input indices in [0,width). All used differences are -(width-1)..width-1.
  index=np.ix_(*[offsets%n for n in kernel.shape])
  self.multiplier=fftn(real_kernel[index])
 def convolve(self,rho):
  if self.full:return self.kernel.convolve(rho)
  w=self.width;shape=rho.shape[:-3]+(self.pad,)*3
  padded=np.zeros(shape,complex);padded[...,:w,:w,:w]=rho
  result=ifftn(fftn(padded,axes=(-3,-2,-1))*self.multiplier,axes=(-3,-2,-1))
  return result[...,:w,:w,:w]


def box_origins(w,width):
 """Circular density centers; supports are fixed for one linear application."""
 origins=[]
 for orbital in w:
  density=abs(orbital)**2;center=[]
  for axis,n in enumerate(orbital.shape):
   marginal=density.sum(axis=tuple(a for a in range(3) if a!=axis))
   z=np.sum(marginal*np.exp(2j*np.pi*np.arange(n)/n))
   if abs(z)<1e-12*density.sum():raise ValueError('Delocalized orbital: no unique support center')
   center.append((np.angle(z)%(2*np.pi))*n/(2*np.pi))
  origins.append(np.floor(np.array(center)-width/2).astype(int)%np.array(orbital.shape))
 return np.array(origins)


def source_box_exchange(sources,targets,shifts,origins,width,kernel,selected=None):
 """Linear Hermitian source-support operator when selected is None.

selected[i,j,R] additionally skips occupied-target pairs and is only suitable
for the separately validated occupied-action approximation, not general hpsi.
"""
 local=LocalConvolution(kernel,width);out=np.zeros_like(targets,dtype=complex)
 for i,origin in enumerate(origins):
  index=np.ix_(*[(np.arange(width)+s)%n for s,n in zip(origin,kernel.shape)])
  source=sources[(i,)+index]
  for ir,shift in enumerate(shifts):
   js=np.arange(len(targets)) if selected is None else np.flatnonzero(selected[i,:,ir])
   if not len(js):continue
   index=np.ix_(*[(np.arange(width)+s+d)%n for s,d,n in zip(origin,shift,kernel.shape)])
   for j in js:
    rho=source.conj()*targets[(j,)+index]
    out[(j,)+index]-=source*local.convolve(rho)
 return out


def symmetric_exchange(w,shifts,kernel,pair_tolerance=0.):
 """Reuse i,j,R <-> j,i,-R pair potentials without changing the operator."""
 start=time.perf_counter();scores=pair_overlaps(w,shifts,kernel.dv);n=len(w);nr=len(shifts)
 inverse=[]
 for shift in shifts:
  match=np.where(np.all(shifts==(-shift)%np.array(kernel.shape),axis=1))[0]
  if len(match)!=1:raise ValueError('Translation mesh must contain unique inverses')
  inverse.append(int(match[0]))
 # Explicit symmetrization only removes roundoff in the overlap test.
 scores=.5*(scores+scores.transpose(1,0,2)[:,:,inverse]);selected=scores>=pair_tolerance
 screening=time.perf_counter()-start;out=np.zeros_like(w,dtype=complex);count=0;start=time.perf_counter()
 for i in range(n):
  for j in range(i,n):
   for ir,shift in enumerate(shifts):
    inv=inverse[ir]
    if i==j and ir>inv:continue
    if not selected[i,j,ir]:continue
    source=np.roll(w[i],tuple(shift),(0,1,2));v=kernel.convolve(source.conj()*w[j]);count+=1
    out[j]-=source*v
    if i!=j or ir!=inv:out[i]-=np.roll(w[j]*v.conj(),tuple(-shift),(0,1,2))
 return out,dict(total_pairs=n*n*nr,retained_pairs=int(selected.sum()),convolutions=count,pair_tolerance=pair_tolerance,
  screening_seconds=screening,exchange_seconds=time.perf_counter()-start)
