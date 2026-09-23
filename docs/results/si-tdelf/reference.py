"""Finite-mesh cold uniform gas calibration; no propagated-state fitting."""
from pathlib import Path
import sys,json
import numpy as np
out=Path(__file__).resolve().parent;root=out.parents[2]
sys.path.insert(0,str(out.parent/'si-time-wannier'));sys.path.insert(0,str(out.parent/'si-ueg-distance'))
from wannier import geometry
from ueg import global_indices,fermi_reference,mean_momentum
r,c,k,b,L,h=geometry(root);indices,q=global_indices(k,12,L,4)
f=fermi_reference(q,16*64);n=f.sum()/(64*L**3)
tau=np.sum(f*np.sum(q*q,axis=-1))/(64*L**3)
j=n*mean_momentum(f,q);d=tau-j@j/n
ref=3/5*(6*np.pi**2)**(2/3)*n**(5/3)
z=dict(one_spin_density=n,tau=tau,continuum_reference=ref,finite_mesh_elf=1/(1+(d/ref)**2),
       continuum_cold_ueg_elf=.5,note='Finite k-grid quadrature is not exact continuum Fermi sphere; no normalization adjustment applied.')
(out/'reference.json').write_text(json.dumps(z,indent=2)+'\n');print(json.dumps(z,indent=2))
