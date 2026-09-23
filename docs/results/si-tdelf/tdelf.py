"""One-spin TDELF: tau has no factor 1/2; equal k weights, occupied states only."""
import numpy as np

def gradients(u,k,length,grid,zero_nyquist=False):
    gaxis=np.fft.fftfreq(grid)*grid*2*np.pi/length
    if zero_nyquist and grid%2==0:gaxis[grid//2]=0.
    gv=np.stack(np.meshgrid(gaxis,gaxis,gaxis,indexing='ij'),axis=-1)
    result=np.empty((3,)+u.shape,complex)
    for ik in range(len(k)):
        cube=u[ik].reshape((grid,)*3+(u.shape[-1],),order='F')
        coef=np.fft.fftn(cube,axes=(0,1,2))
        for a in range(3):
            result[a,ik]=np.fft.ifftn(1j*(gv[...,a]+k[ik,a])[...,None]*coef,axes=(0,1,2)).reshape(-1,u.shape[-1],order='F')
    return result

def fields(u,grad):
    nk=len(u)
    n=np.sum(abs(u)**2,axis=(0,2))/nk
    tau=np.sum(abs(grad)**2,axis=(0,1,3))/nk
    ug=np.sum(u.conj()[None]*grad,axis=(1,3))/nk
    gradrho=2*ug.real
    jp=ug.imag
    if np.any(n<1e-12):raise ValueError('Low density requires an explicit masking policy')
    weiz=np.sum(gradrho**2,axis=0)/(4*n)
    flow=np.sum(jp**2,axis=0)/n
    d=tau-weiz-flow
    if d.min() < -1e-10*max(1.,tau.max()):raise ValueError('Negative Pauli curvature')
    d=np.maximum(d,0.)
    reference=3/5*(6*np.pi**2)**(2/3)*n**(5/3)
    return dict(n=n,tau=tau,jp=jp,curvature=d,elf=1/(1+(d/reference)**2),
                elf_without_current=1/(1+((tau-weiz)/reference)**2),flow_term=flow)
