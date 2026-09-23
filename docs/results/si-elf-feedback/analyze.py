"""Analyze native ELF feedback runs and independently check checkpoint ELF."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;root=out.parents[2]
sys.path.insert(0,str(out.parent/'si-time-wannier'));sys.path.insert(0,str(out.parent/'si-tdelf'))
from wannier import read_u,geometry
from tdelf import fields,gradients
work=Path('/private/tmp/salmon-si-elf-feedback');au_fs=.02418884326505
_,_,k,_,length,h=geometry(root)
def finite_difference(u):
 result=np.empty((3,)+u.shape,complex)
 for ik in range(len(k)):
  cube=u[ik].reshape((12,12,12,u.shape[-1]),order='F')
  for axis in range(3):
   grad=sum(c*(np.roll(cube,-m,axis)-np.roll(cube,m,axis))/h for m,c in enumerate([.8,-.2,4/105,-1/280],1))
   result[axis,ik]=(grad+1j*k[ik,axis]*cube).reshape(-1,u.shape[-1],order='F')
 return result
def measure(u,grad):
 z=fields(u,grad);return float(z['n']@(z['elf']-.5)**2/z['n'].sum())

def main():
 names=('none','weak','strong','fixed_strong','impulse','fixed_impulse','strong_stride1')
 data={};summary={}
 for name in names:
  p=work/name
  if 'end SALMON' not in (p/'outputfile').read_text():raise RuntimeError(f'{name} is incomplete')
  rt=np.loadtxt(p/'Si_rt.data');xc=np.loadtxt(p/'Si_rt_xc.data');assert len(rt)==len(xc)==1600
  assert np.isfinite(rt).all() and np.isfinite(xc).all()
  data[name]=(rt,xc)
  np.savez_compressed(out/f'{name}.npz',rt=rt,xc=xc)
  alpha=xc[:,7] if xc.shape[1]==13 else np.full(len(xc),.2)
  q=xc[:,8] if xc.shape[1]==13 else np.full(len(xc),np.nan)
  trace=np.column_stack((rt[:,0]*au_fs,alpha,q,rt[:,15],rt[:,3],xc[:,3],xc[:,6]))
  np.savetxt(out/f'{name}.csv',trace[9::10],delimiter=',',header='time_fs,alpha,Q_ELF,Jz,Ac_ext_z,Ac_xc_z,Exc_z',comments='')
  post=rt[:,0]>60
  summary[name]=dict(alpha_initial=.2,alpha_final=float(alpha[-1]),alpha_min=float(alpha.min()),alpha_max=float(alpha.max()),
     max_alpha_deviation=float(np.max(abs(alpha-.2))),Q_initial=float(xc[0,12]) if xc.shape[1]==13 else None,
     Q_final=float(q[-1]) if xc.shape[1]==13 else None,
     current_rms_post=float(np.sqrt(np.mean(rt[post,15]**2))),current_final=float(rt[-1,15]),
     xc_potential_final=float(xc[-1,3]),max_abs_current=float(np.max(abs(rt[:,15]))))
 for name,ref in [('impulse','fixed_impulse'),('strong','fixed_strong'),('strong_stride1','strong')]:
  a,ax=data[name];b,bx=data[ref]
  summary[name]['current_relative_l2_vs_'+ref]=float(np.linalg.norm(a[:,15]-b[:,15])/np.linalg.norm(b[:,15]))
  if ax.shape[1]==bx.shape[1]==13:summary[name]['max_alpha_difference_vs_'+ref]=float(np.max(abs(ax[:,7]-bx[:,7])))
 u0,_=read_u(root/'calculations/si_tdcdft_k4/gs/data_for_restart')
 fd0=measure(u0,finite_difference(u0));fft0=measure(u0,gradients(u0,k,length,12,zero_nyquist=True))
 np.testing.assert_allclose(fd0,data['none'][1][0,12],atol=3e-14,rtol=0)
 checks=[]
 for name in ('none','weak','strong'):
  for step in (200,800,1600):
   u,_=read_u(work/name/f'checkpoint_rt_{step:06d}')
   fd=measure(u,finite_difference(u));fft=measure(u,gradients(u,k,length,12,zero_nyquist=True))
   native=float(data[name][1][step-1,8]);np.testing.assert_allclose(fd,native,atol=3e-14,rtol=0)
   checks.append(dict(case=name,step=step,native_Q=native,offline_fd_Q=fd,offline_fft_Q=fft,
                      alpha_if_fft_normalized=.2*fft/fft0))
 summary['derivative_check']=dict(initial_fd=fd0,initial_fft=fft0,snapshots=checks)
 # Finite-mesh cold Fermi sea: check how closely the discretized reference
 # reaches the continuum ELF=1/2 endpoint, without changing the model.
 sys.path.insert(0,str(out.parent/'si-ueg-distance'))
 from ueg import global_indices,fermi_reference
 indices,q=global_indices(k,12,length,4);occupations=fermi_reference(q,16*64).ravel()[indices]
 gaxis=np.fft.fftfreq(12)*12*2*np.pi/length
 g=np.stack(np.meshgrid(gaxis,gaxis,gaxis,indexing='ij'),axis=-1).reshape(-1,3,order='F')
 gfd=sum(2*c*np.sin(m*g*h)/h for m,c in enumerate([.8,-.2,4/105,-1/280],1))
 p=gfd[None]+k[:,None];n=occupations.sum()/(64*length**3)
 tau=np.sum(occupations*np.sum(p*p,axis=2))/(64*length**3)
 j=np.sum(occupations[:,:,None]*p,axis=(0,1))/(64*length**3)
 reference=.6*(6*np.pi**2)**(2/3)*n**(5/3)
 gas_elf=1/(1+((tau-j@j/n)/reference)**2)
 summary['finite_mesh_cold_gas']=dict(elf=float(gas_elf),alpha_on_Si_normalization=float(.2*(gas_elf-.5)**2/fd0),
     note='Continuum endpoint remains .5; no finite-mesh offset was fitted.')
 summary['conditions']=dict(k_grid=[4,4,4],dt_au=.08,nt=1600,pulse_width_au=60,omega_au=.2,alpha0=.2,elf_stride=10,
                            model='alpha0 * density-weighted (ELF-0.5)^2 / initial value; no upper clipping')
 (out/'metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
 fig,ax=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
 for name,color in [('none','0.5'),('weak','tab:blue'),('strong','tab:red')]:
  rt,xc=data[name];t=rt[:,0]*au_fs
  ax[0].plot(t,xc[:,7],label=name,color=color)
 ax[0].set(ylabel='alpha');ax[0].legend()
 for name,style,label in [('strong','-','ELF, stride 10'),('fixed_strong','--','Fixed alpha = 0.2'),('strong_stride1',':','ELF, stride 1')]:
  rt,xc=data[name];t=rt[:,0]*au_fs
  ax[1].plot(t,rt[:,15],style,label=label);ax[2].plot(t,xc[:,3],style,label=label)
 ax[1].set(ylabel='Jz (a.u.)');ax[1].legend();ax[2].set(ylabel='Axc,z / c (a.u.)',xlabel='Time (fs)')
 for a in ax:a.axvline(60*au_fs,color='0.5',ls=':');a.grid(alpha=.2)
 fig.suptitle('Si 4x4x4 k: self-consistent ELF-dependent alpha\nStrong pulse: 1e13 W/cm2; cold uniform-gas endpoint alpha = 0')
 fig.savefig(out/'feedback.png',dpi=160);fig.savefig(out/'feedback.pdf')
 fig,ax=plt.subplots(2,1,figsize=(9,5.5),sharex=True,layout='constrained')
 rt,xc=data['impulse'];fixed=data['fixed_impulse'][0];t=rt[:,0]*au_fs
 ax[0].plot(t,rt[:,15],label='ELF');ax[0].plot(t,fixed[:,15],'--',label='Fixed alpha');ax[0].legend();ax[0].set(ylabel='Jz (a.u.)')
 ax[1].plot(t,xc[:,7]-.2);ax[1].set(ylabel='alpha - 0.2',xlabel='Time (fs)')
 fig.suptitle('Small impulse = 1e-4 a.u.; short-time response check')
 fig.savefig(out/'impulse.png',dpi=160)
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
