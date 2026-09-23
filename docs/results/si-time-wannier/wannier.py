"""Projected occupied Wannier diagnostics; primitive-cell normalized SALMON u_nk."""
from pathlib import Path
import itertools,struct
import numpy as np

def projected_gauge(u,trial,dv):
 a=u.conj().T@trial*dv
 left,s,right=np.linalg.svd(a,full_matrices=False)
 if len(s)!=u.shape[1] or s[-1]<1e-8*max(s[0],1):raise ValueError('Singular or incomplete trial projection')
 return left@right,s

def reconstruct(v,k,r,cells):
 """w_R=0 evaluated on cell translations; shifted-k boundary twist retained.

Cell-normalized u requires 1/Nk (not 1/sqrt(Nk)) for normalized supercell w.
"""
 nk=len(k)
 bloch=v*np.exp(1j*k@r.T)[:,:,None]
 return (np.exp(1j*cells@k.T)@bloch.reshape(nk,-1)/nk).reshape(len(cells),len(r),v.shape[-1])

def moments(rho,pos,lengths,dv):
 norm=rho.sum(axis=0)*dv
 p=rho*dv/norm
 z=np.exp(2j*np.pi*pos/lengths).T@p
 center=(np.angle(z).T%(2*np.pi))*lengths/(2*np.pi)
 delta=(pos[:,None,:]-center[None,:,:]+lengths/2)%lengths-lengths/2
 axes=np.einsum('rn,rna->na',p,delta**2)
 edge=np.any(abs(delta)>.4*lengths,axis=2)
 return dict(norm=norm,center=center,spread=axes.sum(axis=1),spread_axes=axes,
             concentration=abs(z.T),boundary_weight=np.sum(p*edge,axis=0))

def read_u(directory,ngrid=1728):
 d=Path(directory);raw=(d/'info.bin').read_bytes()
 if len(raw)!=60:raise ValueError('Only current five-record shared info.bin is supported')
 values=[]
 for i in range(5):
  lead,val,trail=struct.unpack('<iii',raw[i*12:(i+1)*12])
  if lead!=4 or trail!=4:raise ValueError('Unsupported Fortran record format')
  values.append(val)
 nk,no,step,nproc,real=values
 if real:raise ValueError('This diagnostic expects complex Bloch orbitals')
 size=ngrid*no*nk*16
 if (d/'wfn.bin').stat().st_size!=size:raise ValueError('Wavefunction size/layout mismatch')
 u=np.fromfile(d/'wfn.bin',dtype='<c16').reshape((ngrid,no,nk),order='F').transpose(2,0,1)
 occraw=(d/'occupation.bin').read_bytes();nbytes=8*nk*no
 if struct.unpack('<i',occraw[:4])[0]!=nbytes or struct.unpack('<i',occraw[-4:])[0]!=nbytes:raise ValueError('Occupation layout mismatch')
 occ=np.frombuffer(occraw[4:-4],dtype='<f8').reshape(nk,no)
 if not np.allclose(occ[:,:16],2) or not np.allclose(occ[:,16:],0):raise ValueError('Expected sixteen doubly occupied Si bands')
 return np.array(u[:,:,:16]),values

def geometry(root):
 length=10.26;h=.855;grid=12;mesh=4
 # Explicit meshgrid flattening follows x-fast SALMON Fortran storage.
 xyz=np.stack(np.meshgrid(*(np.arange(grid)*h for _ in range(3)),indexing='ij'),axis=-1).reshape(-1,3,order='F')
 cells=np.array(list(itertools.product(range(mesh),repeat=3)))*length
 kdata=np.loadtxt(Path(root)/'calculations/si_tdcdft_k4/gs/Si_k.data',skiprows=5,max_rows=64)
 if len(kdata)!=64 or not np.allclose(kdata[:,4],1/64):raise ValueError('Expected full4³ equal-weight k mesh')
 k=kdata[:,1:4]*2*np.pi/length
 atoms=np.array([[0,0,0],[2.565,2.565,2.565],[5.13,0,5.13],[0,5.13,5.13],[5.13,5.13,0],[7.695,2.565,7.695],[2.565,7.695,7.695],[7.695,7.695,2.565]])
 centers=[]
 for i in range(8):
  for j in range(i+1,8):
   delta=(atoms[j]-atoms[i]+length/2)%length-length/2
   if abs(np.linalg.norm(delta)-length*np.sqrt(3)/4)<1e-8:centers.append((atoms[i]+delta/2)%length)
 assert len(centers)==16
 return xyz,cells,k,np.array(centers),length,h

def initial_gauge(u,k,r,centers,length,dv,sigma=1.2):
 translations=np.array(list(itertools.product((-1,0,1),repeat=3)))*length
 localized=np.exp(-np.sum((r[None,:,None,:]-centers[None,None,:,:]-translations[:,None,None,:])**2,axis=-1)/(2*sigma**2))
 gauges=[];singular=[]
 for ik,kv in enumerate(k):
  trial=np.einsum('t,trn->rn',np.exp(1j*translations@kv),localized)*np.exp(-1j*r@kv)[:,None]
  trial/=np.sqrt(np.sum(abs(trial)**2,axis=0)*dv)
  g,s=projected_gauge(u[ik],trial,dv);gauges.append(g);singular.append(s)
 return np.array(gauges),np.array(singular)
