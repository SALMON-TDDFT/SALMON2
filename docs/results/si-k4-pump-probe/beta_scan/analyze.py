from pathlib import Path
import os,sys,json
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;root=out.parents[3];base=root/'calculations/si_hse_native/k4_pump_probe'
sys.path.insert(0,str(root/'samples/exercise_si_tdcdft'));from analyze import response
energy=np.arange(.1,6.001,.01);metrics={};fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for beta in [0,.001,.004,.016]:
 name='td_g40_ground' if beta==0 else f'td_g40_b{round(beta*1000):03d}_ground';d=base/name
 assert json.loads((d/'status.json').read_text())['completed']
 r=np.loadtxt(d/'Si_rt.data');a=np.loadtxt(d/'Si_rt_xc.data');en=np.loadtxt(d/'Si_rt_energy.data')
 assert np.isfinite(r).all() and np.isfinite(a).all() and np.isfinite(en).all() and r[-1,0]>=959.9
 m={'windows':{},'max_abs_Axc':float(abs(a[:,3]).max()),'max_abs_J':float(abs(r[:,15]).max())}
 for duration in [600,720,880]:
  z=(r[:,0]>80)&(r[:,0]<=80+duration+1e-7);ep=response(r[z,0]-80,r[z,15],1e-4,energy)
  neg=(energy>=1.7)&(energy<=2.4);pos=(energy>=2.4)&(energy<=3.1)
  m['windows'][duration]={'negative_min':float(ep.imag[neg].min()),'negative_area':float(np.trapezoid(np.minimum(ep.imag[neg],0),energy[neg])),'positive_max':float(ep.imag[pos].max()),'positive_area':float(np.trapezoid(np.maximum(ep.imag[pos],0),energy[pos]))}
  np.savetxt(out/f'beta{beta}_window{duration}.csv',np.column_stack([energy,ep.real,ep.imag]),delimiter=',',header='energy_eV,Re_epsilon,Im_epsilon')
  if duration==880:axes[0].plot(energy,ep.imag,label=f'beta={beta}')
 axes[1].plot((a[:,0]-80)*.024188843,a[:,3],label=f'beta={beta}');metrics[beta]=m
axes[0].set(xlim=(1.5,4),xlabel='Energy (eV)',ylabel='Im epsilon',title='Unpumped, alpha=.2 gamma=.004');axes[0].axhline(0,color='gray',lw=.5)
axes[1].set(xlim=(0,21.3),xlabel='Time after impulse (fs)',ylabel='Axc/c (a.u.)')
for ax in axes:ax.legend();ax.grid(alpha=.2)
fig.savefig(out/'comparison.png',dpi=160);(out/'metrics.json').write_text(json.dumps(metrics,indent=2))
