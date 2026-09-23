"""Central probe derivative and finite-window apparent width."""
import numpy as np

def central_response(time_plus,j_plus,time_minus,j_minus,amplitude):
 t=np.asarray(time_plus);u=np.asarray(time_minus);p=np.asarray(j_plus);m=np.asarray(j_minus)
 if not np.isfinite(amplitude) or amplitude<=0:raise ValueError('Positive finite amplitude required')
 if t.shape!=u.shape or p.shape!=t.shape or m.shape!=u.shape or not np.allclose(t,u,rtol=0,atol=1e-10):raise ValueError('Probe time grids differ')
 if not all(np.isfinite(x).all() for x in (t,u,p,m)):raise ValueError('Nonfinite probe data')
 return (p-m)/(2*amplitude)

def apparent_width(energy,absorption,low=2.,high=4.):
 mask=(energy>=low)&(energy<=high);x=energy[mask];y=absorption[mask]
 if len(x)<3:return None
 i=int(np.argmax(y));half=y[i]/2
 if half<=0:return None
 left=np.where(y[:i]<half)[0];right=np.where(y[i+1:]<half)[0]
 if not len(left) or not len(right):return None
 l=left[-1];r=i+1+right[0]
 xl=x[l]+(half-y[l])/(y[l+1]-y[l])*(x[l+1]-x[l])
 xr=x[r-1]+(half-y[r-1])/(y[r]-y[r-1])*(x[r]-x[r-1])
 return float(xr-xl)
