"""SALMON-grid independent reference Hamiltonian, atomic units, fixed ions.

Use exported native pseudopotential projectors, finite-difference coefficients
and k mesh; validate frozen native Hpsi before a self-consistent calculation.
"""
from pathlib import Path
import sys
import numpy as np
from numpy.fft import fftn,ifftn


def kinetic_symbol(shape,h,k,lap0,lap,nab):
 q=np.stack(np.meshgrid(*[2*np.pi*np.fft.fftfreq(n,d=h) for n in shape],indexing='ij'),axis=0)
 result=np.full(shape,lap0+.5*np.dot(k,k))
 for axis in range(3):
  for d in range(1,5):result-=lap[d-1,axis]*np.cos(d*h*q[axis]);result+=2*k[axis]*nab[d-1,axis]*np.sin(d*h*q[axis])
 return result


def gradient(field,nab):
 result=np.zeros((3,)+field.shape,dtype=field.dtype)
 for axis in range(3):
  for d in range(1,5):result[axis]+=nab[d-1,axis]*(np.roll(field,-d,axis=axis)-np.roll(field,d,axis=axis))
 return result


def hartree(rho,h):
 q=np.stack(np.meshgrid(*[2*np.pi*np.fft.fftfreq(n,d=h) for n in rho.shape],indexing='ij'),axis=0)
 q2=np.sum(q*q,axis=0);multiplier=np.zeros_like(q2);np.divide(4*np.pi,q2,out=multiplier,where=q2>0)
 v=ifftn(fftn(rho)*multiplier).real
 return v,float(.5*np.sum(rho*v)*h**3)


def semilocal_potential(rho,nab,xc,dv):
 g=gradient(rho,nab);eps,vrho,vsigma=xc.evaluate(rho,np.sum(g*g,axis=0));v=vrho.copy()
 for axis in range(3):
  field=2*vsigma*g[axis]
  for d in range(1,5):v-=nab[d-1,axis]*(np.roll(field,-d,axis=axis)-np.roll(field,d,axis=axis))
 return v,float(np.sum(rho*eps)*dv)


def orthonormalize(u,dv):
 shape=u.shape;flat=u.reshape(shape[0],shape[1],-1);out=np.empty_like(flat)
 for ik,a in enumerate(flat):
  gram=a@a.conj().T*dv;e,v=np.linalg.eigh(gram)
  if e.min()<1e-12:raise ValueError('Linearly dependent occupied orbitals')
  out[ik]=((v/np.sqrt(e))@v.conj().T)@a
 return out.reshape(shape)


class NativeModel:
 def __init__(self,directory):
  self.path=Path(directory)
  marker=self.path/'complete.txt'
  if not marker.exists() or marker.read_text().strip()!='SALMON_HSE_REFERENCE_V1_COMPLETE':
   raise ValueError('Native export incomplete: missing or invalid completion marker')
  lines=(self.path/'metadata.txt').read_text().splitlines()
  if lines[0].strip()!='SALMON_HSE_REFERENCE_V1':raise ValueError('Unknown export format')
  self.metadata={c[0]:c[1:] for line in lines[1:] if (c:=line.split()) and not c[0].startswith('#')}
  if self.metadata.get('endian_little')!=['T' if sys.byteorder=='little' else 'F']:
   raise ValueError('Native export endianness does not match reader')
  if self.metadata.get('rho_semantics')!=['scf_mixed_potential_density']:
   raise ValueError('Native export rho density semantics unspecified')
  m=self.metadata;self.shape=tuple(map(int,m['grid']));self.nk=int(m['nk'][0]);self.no=int(m['no'][0]);self.nlma=int(m['nlma'][0]);self.dv=float(m['hvol'][0]);self.h=self.dv**(1/3)
  if len(set(self.shape))!=1:raise ValueError('Cubic reference grid required')
  geometry=self.read('geometry',(21+3*int(m['nion'][0]),))
  self.hgs=geometry[:3].copy();self.primitive_a=geometry[3:12].reshape(3,3,order='F').copy()
  expected_cell=np.diag(np.asarray(self.shape)*self.h)
  if not np.isfinite(geometry).all() or not np.allclose(self.hgs,self.h,rtol=1e-12,atol=1e-13):
   raise ValueError('Native reference requires isotropic physical geometry')
  if not np.allclose(self.primitive_a,expected_cell,rtol=1e-12,atol=1e-13):
   raise ValueError('Native reference requires axis-aligned orthogonal cubic geometry consistent with grid')
  self.k=self.read('k',(3,self.nk)).T;self.weights=self.read('weights',(self.nk,));self.occupations=self.read('occupations',(self.no,self.nk)).T
  self.psi=self.read('psi',self.shape+(self.no,self.nk),complex).transpose(4,3,0,1,2)
  self.native_hpsi=self.read('hpsi',self.shape+(self.no,self.nk),complex).transpose(4,3,0,1,2)
  self.projectors=self.read('projectors',self.shape+(self.nlma,self.nk),complex).transpose(4,3,0,1,2).reshape(self.nk,self.nlma,-1)
  self.rinv=self.read('rinv_uvu',(self.nlma,));self.vps=self.read('vpsl',self.shape);self.rho=self.read('rho',self.shape)
  self.native_local=self.read('vlocal',self.shape);self.native_energies=self.read('energies',(7,))
  coeff=self.read('stencil',(25,));self.lap0=coeff[0];self.lap=coeff[1:13].reshape(4,3,order='F');self.nab=coeff[13:].reshape(4,3,order='F')
  self.tsymbol=np.array([kinetic_symbol(self.shape,self.h,k,self.lap0,self.lap,self.nab) for k in self.k])
  self.nlcc=self.read('rho_nlcc',self.shape) if m.get('nlcc_available')==['T'] else np.zeros(self.shape)
  self.volume=self.dv*np.prod(self.shape)
  self._zero_tsymbol=self.tsymbol.copy()
  q=np.stack(np.meshgrid(*[2*np.pi*np.fft.fftfreq(n,d=self.h) for n in self.shape],indexing='ij'),axis=0)
  self._momentum=np.zeros((3,)+self.shape)
  for axis in range(3):
   for d in range(1,5):self._momentum[axis]+=2*self.nab[d-1,axis]*np.sin(d*self.h*q[axis])
  self.nps=int(m['nps'][0])
  self._raw_uv=self.read('raw_projectors',(self.nps,self.nlma)).T.copy()
  self._raw_xyz=self.read('projector_positions',(3,self.nps,self.nlma)).transpose(0,2,1).copy()
  indices=self.read('projector_indices',(3,self.nps,self.nlma),np.int32).transpose(0,2,1)-1
  counts=self.read('projector_counts',(self.nlma,),np.int32)
  if np.any(counts<0) or np.any(counts>self.nps):raise ValueError('Invalid raw projector counts')
  valid=np.arange(self.nps)[None,:]<counts[:,None]
  if np.any(indices[:,valid]<0) or np.any(indices[:,valid]>=np.array(self.shape)[:,None]):
   raise ValueError('Invalid raw projector indices')
  indices[:,~valid]=0
  self._raw_uv[~valid]=0.
  self._raw_xyz[:,~valid]=0.
  self._raw_indices=np.ravel_multi_index(tuple(indices),self.shape)
  self._projector_channels=np.broadcast_to(np.arange(self.nlma)[:,None],self._raw_indices.shape)
  self.set_field(np.zeros(3))
 def read(self,name,shape,dtype=float):
  a=np.fromfile(self.path/(name+'.bin'),dtype=np.complex128 if dtype is complex else np.int32 if dtype is np.int32 else np.float64)
  if a.size!=int(np.prod(shape)):raise ValueError(f'{name}: export size differs')
  return a.reshape(shape,order='F')
 def set_field(self,A):
  """Set uniform A/c in atomic units, rebuilding k+A kinetic/nonlocal terms.

  Each raw projector image retains its own displacement before duplicate grid
  entries are summed; multiplying an already accumulated dense B is incorrect.
  """
  A=np.asarray(A,dtype=float)
  if A.shape!=(3,) or not np.isfinite(A).all():raise ValueError('Field must be a finite three-vector')
  self.A=A.copy()
  self.tsymbol=self._zero_tsymbol+np.einsum('a,axyz->xyz',A,self._momentum)[None]
  self.tsymbol+=(self.k@A+.5*np.dot(A,A))[:,None,None,None]
  self._raw_field=self._raw_uv[None]*np.exp(-1j*np.einsum('ka,alj->klj',self.k+A,self._raw_xyz))
  self.projectors.fill(0.)
  for ik in range(self.nk):
   np.add.at(self.projectors[ik],(self._projector_channels,self._raw_indices),self._raw_field[ik])

 def current(self,u):
  """Electron-number current = (1/V) d E_core / d(A/c), fixed orbitals.

  Includes the finite-difference kinetic velocity and derivative of nonlocal
  projectors. Spin degeneracy is two and supplied orbitals are fully occupied,
  matching density()/expectation(). This is not conventional negative-charge
  current. Density is taken from u, never from the mixed exported SCF density.
  No extra EXX current is added for self-consistent full Fock exchange.
  """
  u=np.asarray(u)
  if u.ndim!=5 or u.shape[0]!=self.nk or u.shape[2:]!=self.shape:
   raise ValueError('Orbital shape differs from native model')
  fourier=fftn(u,axes=(-3,-2,-1))
  momentum_density=np.sum(abs(fourier)**2,axis=1)
  current=2*self.dv/np.prod(self.shape)*np.einsum('k,kxyz,axyz->a',self.weights,momentum_density,self._momentum)
  norms=np.sum(abs(u)**2,axis=(1,2,3,4))
  current+=2*self.dv*np.einsum('k,k,ka->a',self.weights,norms,self.k+self.A)
  flat=u.reshape(self.nk,u.shape[1],-1)
  for ik in range(self.nk):
   selected=flat[ik][:,self._raw_indices]
   projected=np.einsum('nlj,lj->nl',selected,self._raw_field[ik].conj())
   # c=<B|u>; d c/dA = <dB|u> = sum i*r*conj(B_raw)*u.
   dprojected=1j*np.einsum('nlj,lj,alj->anl',selected,self._raw_field[ik].conj(),self._raw_xyz)
   current+=4*self.dv*self.weights[ik]*np.einsum('nl,anl,l->a',projected.conj(),dprojected,self.rinv).real
  return current/self.volume

 def kinetic(self,u):
  return ifftn(fftn(u,axes=(-3,-2,-1))*self.tsymbol[:,None],axes=(-3,-2,-1))
 def nonlocal_potential(self,u):
  flat=u.reshape(self.nk,u.shape[1],-1);out=np.empty_like(flat)
  for ik in range(self.nk):out[ik]=((flat[ik]@self.projectors[ik].conj().T)*self.rinv)@self.projectors[ik]
  return out.reshape(u.shape)
 def density(self,u):
  return 2*np.einsum('k,knxyz->xyz',self.weights,abs(u)**2)
 def expectation(self,u,hu):
  return float(2*np.einsum('k,knxyz,knxyz->',self.weights,u.conj(),hu).real*self.dv)
 def core(self,u):return self.kinetic(u)+self.vps*u+self.nonlocal_potential(u)
 def native_parity(self):
  got=self.kinetic(self.psi)+self.native_local*self.psi+self.nonlocal_potential(self.psi)
  vh,eh=hartree(self.rho,self.h);reference_vh=self.read('vh',self.shape)
  return dict(hpsi_relative_error=float(np.linalg.norm(got-self.native_hpsi)/np.linalg.norm(self.native_hpsi)),
   hartree_potential_max_error=float(abs(vh-reference_vh).max()),hartree_energy_error_Ha=float(eh-self.native_energies[2]),
   nonlocal_energy_error_Ha=float(self.expectation(self.psi[:,:16],self.nonlocal_potential(self.psi[:,:16]))-self.native_energies[6]))
