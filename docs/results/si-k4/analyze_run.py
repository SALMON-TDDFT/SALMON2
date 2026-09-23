from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path.cwd();sys.path.insert(0,str(root/'samples/exercise_si_tdcdft'))
from analyze import response,read_current,HARTREE_EV,peak_metrics
out=root/'docs/results/si-k4';runs=root/'calculations/si_tdcdft_k4'
energy=np.arange(.01,8.001,.01)
data={};spectra={}
for mode in ['alda','lrc']:
 p=runs/mode
 if 'end SALMON' not in (p/'outputfile').read_text():raise RuntimeError(f'{mode} incomplete')
 r=read_current(p/'Si_rt.data')
 assert r.shape[0]==12000 and np.isfinite(r).all()
 eps=response(r[:,0],r[:,15],.001,energy)
 data[mode]=r;spectra[mode]=eps
 np.savetxt(out/f'{mode}.csv',np.c_[energy,eps.real,eps.imag],delimiter=',',header='energy_eV,Re_epsilon,Im_epsilon',comments='')
xc=np.loadtxt(runs/'lrc/Si_rt_xc.data');assert np.isfinite(xc).all()
metrics={}
for mode,eps in spectra.items():
 vals=[]
 for i in range(1,len(energy)-1):
  if eps.imag[i]>eps.imag[i-1] and eps.imag[i]>=eps.imag[i+1] and 1.5<=energy[i]<=6:
   vals.append(dict(energy_eV=float(energy[i]),height=float(eps.imag[i])))
 vals.sort(key=lambda p:p['height'],reverse=True)
 metrics[mode]=dict(local_maxima_1p5_to_6_eV=vals[:6],
   peak_2_to_4_eV=peak_metrics(energy,eps.imag,2,4),
   peak_4_to_6_eV=peak_metrics(energy,eps.imag,4,6),
   area_2_to_6_eV=float(np.trapezoid(eps.imag[(energy>=2)&(energy<=6)],energy[(energy>=2)&(energy<=6)])),
   max_abs_current_au=float(np.max(np.abs(data[mode][:,15]))))
metrics['settings']=dict(kgrid=[4,4,4],rgrid=[12,12,12],alpha=.2,dt_au=.08,nt=12000,
 duration_fs=float(data['alda'][-1,0]*.02418884326505),kick_au=.001,polarization='z',
 window='1-3(t/T)^2+2(t/T)^3',energy_shift_eV=0,additional_broadening_eV=0)
metrics['max_abs_Axc_over_c_au']=float(np.max(abs(xc[:,3])))
sensitivity={}
for mode in ['alda','lrc']:
 sensitivity[mode]=[]
 for n in [6000,9000,12000]:
  d=data[mode][:n]
  e=response(d[:,0],d[:,15],.001,energy)
  sensitivity[mode].append(dict(duration_fs=float(d[-1,0]*.02418884326505),
      peak_2_to_4_eV=peak_metrics(energy,e.imag,2,4)))
(out/'window_sensitivity.json').write_text(json.dumps(sensitivity,indent=2)+'\n')
(out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
fig,(ax,dx)=plt.subplots(2,1,figsize=(8,6),sharex=True,gridspec_kw={'height_ratios':[3,1]},layout='constrained')
for mode,label in [('alda','PZ-ALDA'),('lrc',r'PZ-ALDA + LRC ($\alpha=0.2$)')]:
 ax.plot(energy,spectra[mode].imag,label=label,lw=1.5)
ax.set(ylabel=r'Im $\epsilon_{zz}$',title=r'Si impulse response: fixed $4^3$ k points')
visible=(energy>=1.5)&(energy<=6)
yvis=np.concatenate([eps.imag[visible] for eps in spectra.values()])
ax.set_ylim(min(0,float(yvis.min()))-5,float(yvis.max())*1.1)
ax.legend();ax.grid(alpha=.2)
dx.plot(energy,spectra['lrc'].imag-spectra['alda'].imag,color='tab:purple')
dx.axhline(0,color='gray',lw=.6);dx.set(xlabel='Photon energy (eV), no energy shift',ylabel='LRC - ALDA',xlim=(1.5,6));dx.grid(alpha=.2)
delta=(spectra['lrc'].imag-spectra['alda'].imag)[visible]
dx.set_ylim(float(delta.min())-5,float(delta.max())+5)
fig.savefig(out/'absorption.png',dpi=180);fig.savefig(out/'absorption.pdf')
plt.close(fig)
fig,(ax,bx)=plt.subplots(2,1,figsize=(8,5),sharex=True,layout='constrained')
for mode in data:
 r=data[mode];ax.plot(r[:,0]*.02418884326505,r[:,15]/.001,label=mode.upper(),lw=.8)
ax.set(ylabel='j / kick (a.u.)');ax.legend()
bx.plot(xc[:,0]*.02418884326505,xc[:,3]);bx.set(xlabel='Time (fs)',ylabel='Axc / c (a.u.)')
fig.savefig(out/'time_trace.png',dpi=160)
print(json.dumps(metrics,indent=2))
