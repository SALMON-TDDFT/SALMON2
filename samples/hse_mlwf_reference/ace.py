"""Per-k adaptively compressed exchange with grid-weighted inner products.

Input exchange excludes the HSE mixing coefficient. Exact on the construction
span, not generally on new targets or a changed occupied density matrix.
"""
import numpy as np

class ACE:
 def __init__(self,orbitals,exchange_action,dv):
  u=np.asarray(orbitals);w=np.asarray(exchange_action)
  if u.ndim<3 or u.shape!=w.shape or not np.isfinite(dv) or dv<=0:
   raise ValueError('Matching k/orbital/grid arrays and positive volume element required')
  if not np.isfinite(u).all() or not np.isfinite(w).all():raise ValueError('Nonfinite ACE input')
  self.grid_shape=u.shape[2:];self.nk=u.shape[0];self.dv=dv
  a=u.reshape(*u.shape[:2],-1);b=w.reshape(a.shape)
  metric=-(a.conj()@b.transpose(0,2,1))*dv
  scale=np.linalg.norm(metric,axis=(1,2))
  if np.any(scale==0) or np.any(np.linalg.norm(metric-metric.conj().transpose(0,2,1),axis=(1,2))>1e-10*scale):
   raise ValueError('Exchange metric must be nonzero Hermitian')
  metric=(metric+metric.conj().transpose(0,2,1))*.5
  e,v=np.linalg.eigh(metric)
  if np.any(e[:,0]<=1e-12*e[:,-1]) or np.any(e[:,-1]<=0):
   raise ValueError('Exchange metric is not numerically positive definite')
  self.factors=(v.transpose(0,2,1)@b)/np.sqrt(e)[:,:,None]
  self.condition_max=float(np.max(e[:,-1]/e[:,0]))
 def apply(self,targets):
  x=np.asarray(targets)
  if x.ndim<3 or x.shape[0]!=self.nk or x.shape[2:]!=self.grid_shape:
   raise ValueError('Target k/grid shape differs')
  a=x.reshape(*x.shape[:2],-1);f=self.factors
  return (-((a@f.conj().transpose(0,2,1))*self.dv)@f).reshape(x.shape)
