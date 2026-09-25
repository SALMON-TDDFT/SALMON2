from pathlib import Path
import sys,json,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;root=out.parents[3];sys.path.insert(0,str(root/'samples/exercise_si_tdcdft'));from analyze import response
b=root/'calculations/si_hse_native/k4_pump_probe';ev=27.211386245988;energy=np.arange(.05,5.001,.005);metrics={};fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
for label,gamma in [('g10',.001),('g40',.004)]:
 r=np.loadtxt(b/f'td_{label}_ground/Si_rt.data');xc=np.loadtxt(b/f'td_{label}_ground/Si_rt_xc.data');t=r[:,0];sel=t>80;tt=t[sel]-80;j=r[sel,15];a=xc[sel,3];dt=tt[1]-tt[0]
 # Direct recurrence check, away from impulse. Discrete relation removes window boundary terms.
 residual=(xc[2:,3]-2*xc[1:-1,3]+xc[:-2,3])/dt**2+gamma*xc[1:-1,3]-.2*r[1:-1,15]
 m=dict(bare_auxiliary_eV=float(np.sqrt(gamma)*ev),recurrence_relative_L2=float(np.linalg.norm(residual)/np.linalg.norm(.2*r[1:-1,15])),windows={})
 for duration in [440,600,720,880]:
  z=tt<=duration+1e-8;ep=response(tt[z],j[z],1e-4,energy);band=np.where((energy>1.7)&(energy<2.4))[0];i=band[np.argmin(ep.imag[band])];m['windows'][duration]=dict(min_eV=float(energy[i]),min_Im=float(ep.imag[i]))
  if label=='g40':ax[0,0].plot(energy,ep.imag,label=f'{duration*.024188843:.1f} fs')
 ax[1,0].plot(tt*.024188843,a,label=f'gamma={gamma}')
 # Exponential window is a diagnostic only, no claim of physical damping.
 m['exponential_windows']={}
 for eta_eV in [.1,.2,.4]:
  om=energy/ev;weighted=j*np.exp(-eta_eV/ev*tt)*dt/1e-4
  ft=np.array([np.sum(weighted*np.exp(1j*w*tt)) for w in om]);ep=1+4*np.pi*1j*ft/om;band=np.where((energy>1.7)&(energy<2.4))[0];i=band[np.argmin(ep.imag[band])]
  m['exponential_windows'][eta_eV]=dict(min_eV=float(energy[i]),min_Im=float(ep.imag[i]))
  if label=='g40':ax[0,1].plot(energy,ep.imag,label=f'eta={eta_eV} eV')
 if label=='g40':
  # Variable projection fit to three dominant modes plus a constant.
  fitmask=tt<880;tf=tt[fitmask];jf=j[fitmask]
  def design(freq):return np.column_stack([np.ones(len(tf))]+[v for f in freq for v in [np.cos(f/ev*tf),np.sin(f/ev*tf)]])
  def fun(freq):M=design(freq);return (M@np.linalg.lstsq(M,jf,rcond=None)[0]-jf)/np.linalg.norm(jf)
  freq=np.array([2.06,2.76,3.3]);step=.02
  for iteration in range(45):
   improved=False
   for k in range(3):
    candidates=[]
    for shift in [-step,0,step]:
     trial=freq.copy();trial[k]+=shift;candidates.append((np.linalg.norm(fun(trial)),trial))
    best=min(candidates,key=lambda x:x[0])[1];improved=improved or bool(np.any(best!=freq));freq=best
   if not improved:step*=.5
   if step<1e-7:break
  M=design(freq);cj=np.linalg.lstsq(M,jf,rcond=None)[0];ca=np.linalg.lstsq(M,a[fitmask],rcond=None)[0]
  modes=[]
  for k,f in enumerate(freq):
   z=cj[1+2*k]+1j*cj[2+2*k];za=ca[1+2*k]+1j*ca[2+2*k];ratio=za/z
   modes.append(dict(eV=float(f),J_cos=float(z.real),J_sin=float(z.imag),A_over_J_abs=float(abs(ratio)),phase_degrees=float(np.angle(ratio,deg=True)),oscillator_ratio_expected=float(.2/(gamma-(f/ev)**2))))
  m['fit']=dict(relative_residual=float(np.linalg.norm(fun(freq))),modes=modes)
  ax[1,1].plot(tf*.024188843,jf,label='J data');ax[1,1].plot(tf*.024188843,M@cj,'--',label='3-mode fit')
 metrics[label]=m
for x in ax[0]:x.set(xlim=(1.7,2.4),ylim=(-200,40),xlabel='Energy (eV)',ylabel='Im epsilon');x.axhline(0,color='gray',lw=.5)
ax[0,0].set_title('gamma=.004: cubic windows');ax[0,1].set_title('gamma=.004: exponential windows (diagnostic)')
ax[1,0].set(xlabel='Time after probe (fs)',ylabel='Axc/c (a.u.)');ax[1,1].set(xlabel='Time after probe (fs)',ylabel='Jz (a.u.)')
for x in ax.flat:x.legend();x.grid(alpha=.2)
fig.savefig(out/'diagnostic.png',dpi=150);(out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(metrics,indent=2))
