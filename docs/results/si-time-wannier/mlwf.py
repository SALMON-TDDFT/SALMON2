"""Discrete Marzari-Vanderbilt spread on an orthogonal full k mesh.

Occupied subspace only. Minimize unitary gauge dependence; report convergence,
not a guarantee of the global minimum. Atomic units throughout.
"""
import numpy as np

def functional(u,raw,neighbors,b,weights):
 nk,nb,n,_=raw.shape
 m=u.conj().transpose(0,2,1)[:,None,:,:]@raw@u[neighbors]
 diag=np.diagonal(m,axis1=2,axis2=3)
 if np.min(abs(diag))<1e-12:raise ValueError('Vanishing diagonal overlap: ill-defined phase')
 theta=np.angle(diag)
 centers=-np.einsum('b,ba,kbn->na',weights,b,theta)/nk
 residual=theta+np.einsum('ba,na->bn',b,centers)[None,:,:]
 spread=float(np.sum(weights[None,:,None]*(1-abs(diag)**2+residual**2))/nk)
 omega_i=float(np.sum(weights[None,:]*(n-np.sum(abs(m)**2,axis=(2,3))))/nk)
 # Line search on the nonnegative gauge-dependent part avoids subtracting
 # changes near machine precision from the much larger invariant spread.
 offdiag=m.copy()
 idx=np.arange(n);offdiag[:,:,idx,idx]=0
 variable=float(np.sum(weights[None,:]*(np.sum(abs(offdiag)**2,axis=(2,3))+np.sum(residual**2,axis=2)))/nk)
 q=(-2*diag.conj()-2j*residual/diag)*weights[None,:,None]/nk
 back=-np.sum(m*q[:,:,None,:],axis=1)
 np.add.at(back,neighbors.reshape(-1),(q[:,:,:,None]*m).reshape(-1,n,n))
 descent=(back-back.conj().transpose(0,2,1))/2
 return spread,descent,dict(centers=centers,omega_invariant=omega_i,variable_spread=variable,
                           spreads=np.sum(weights[None,:,None]*(1-abs(diag)**2+residual**2),axis=(0,1))/nk)

def rotate(u,d,step):
 eig,v=np.linalg.eigh(-1j*d)
 exp=(v*np.exp(1j*step*eig)[:,None,:])@v.conj().transpose(0,2,1)
 return u@exp

def minimize(initial,raw,neighbors,b,weights,maxiter=300,tolerance=1e-6):
 u=initial.copy();step=1.;history=[];converged=False
 for iteration in range(maxiter+1):
  f,d,info=functional(u,raw,neighbors,b,weights);norm=float(np.linalg.norm(d));history.append(f)
  if norm<tolerance:converged=True;break
  if iteration==maxiter:break
  for backtrack in range(35):
   candidate=rotate(u,d,step)
   try:trial=functional(candidate,raw,neighbors,b,weights)[2]['variable_spread']
   except ValueError:trial=float('inf')
   if trial<=info['variable_spread']-1e-4*step*norm**2:break
   step*=.5
  else:break
  u=candidate;step=min(step*1.5,1000.)
 return u,dict(converged=converged,iterations=iteration,spread=history[-1],initial_spread=history[0],gradient_norm=norm,history=history)

def overlap_mesh(orbitals,k,r,length,dv):
 nk,nr,n=orbitals.shape;mesh=round(nk**(1/3));delta=2*np.pi/(mesh*length)
 b=np.vstack([np.eye(3)*delta,-np.eye(3)*delta]);neighbors=np.empty((nk,6),int);raw=np.empty((nk,6,n,n),complex)
 for ik,kv in enumerate(k):
  for ib,bv in enumerate(b[:3]):
   target=kv+bv;distance=(k-target+np.pi/length)%(2*np.pi/length)-np.pi/length
   j=int(np.argmin(np.linalg.norm(distance,axis=1)))
   if np.linalg.norm(distance[j])>1e-9:raise ValueError('Nonuniform k mesh')
   gvec=target-k[j]
   if np.linalg.norm(gvec)<1e-10:gvec=np.zeros(3)
   phase=np.exp(-1j*r@gvec)
   mm=orbitals[ik].conj().T@(orbitals[j]*phase[:,None])*dv
   neighbors[ik,ib]=j;raw[ik,ib]=mm
   neighbors[j,ib+3]=ik;raw[j,ib+3]=mm.conj().T
 weights=np.full(6,1/(2*delta**2))
 return raw,neighbors,b,weights


def transport_gauge(current,previous_wannier,dv):
 """Previous U expressed in the closest current occupied frame (polar transport).

Previous_wannier is previous_u @ previous_U. This is an initial guess only;
subsequent spread minimization still uses the instantaneous orbitals.
 """
 overlap=current.conj().transpose(0,2,1)@previous_wannier*dv
 left,singular,right=np.linalg.svd(overlap)
 if singular.min()<1e-8:raise ValueError('Singular temporal occupied-subspace overlap')
 return left@right,singular
