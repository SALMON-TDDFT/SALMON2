"""Matched finite-window central probe comparison; no shifts or extra damping."""
from pathlib import Path
import os,sys,json
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;root=out.parents[2];b=root/'calculations/si_hse_native/k4_pump_probe'
sys.path.insert(0,str(root/'samples/exercise_si_tdcdft'));from analyze import response
energy=np.arange(.01,8.001,.01);metrics={};spectra={}
for method,prefix,dt in [('HSE06','hse',.16),('TDCDFT','td_g40',.08)]:
 data={};m={}
 for kind in ['pump','plus','minus','ground']:
  d=b/f'{prefix}_{kind}';status=json.loads((d/'status.json').read_text());assert status['completed']
  r=np.loadtxt(d/'Si_rt.data');assert len(r)==round(960/dt) and np.isfinite(r).all();data[kind]=r
  norms=[]
  for line in (d/'run.log').read_text().splitlines():
   c=line.split()
   if len(c)!=7:continue
   try:
    step=int(c[0]);v=list(map(float,c[1:]))
    if step>0 and abs(v[0]-step*dt*.02418884326505)<1e-6:norms.append(v[4])
   except ValueError:pass
  assert norms and max(abs(n-32) for n in norms)<1e-4
  m[kind+'_max_norm_error']=max(abs(n-32) for n in norms)
 t=data['pump'][:,0];pre=t<=80;post=~pre
 for kind in ['plus','minus']:
  np.testing.assert_allclose(data[kind][pre],data['pump'][pre],atol=1e-10,rtol=1e-8)
  np.testing.assert_allclose(data[kind][post,3]-data['pump'][post,3],1e-4 if kind=='plus' else -1e-4,atol=1e-13,rtol=0)
 central=(data['plus'][:,15]-data['minus'][:,15])/2e-4
 even=(data['plus'][:,15]+data['minus'][:,15]-2*data['pump'][:,15])/2e-4
 m['even_over_odd_current']=float(np.linalg.norm(even[post])/np.linalg.norm(central[post]));m['windows']={}
 for w in [600,720,880]:
  mask=post&(t<=80+w+1e-8)
  ep=response(t[mask]-80,central[mask],1,energy);eg=response(t[mask]-80,data['ground'][mask,15],1e-4,energy)
  spectra[method,w]=(eg,ep)
  np.savetxt(out/f'{method}_{w}au.csv',np.column_stack([energy,eg.real,eg.imag,ep.real,ep.imag]),delimiter=',',header='eV,ground_re,ground_im,pumped_re,pumped_im')
  band=(energy>=2.5)&(energy<=4.5);indices=np.where(band)[0];i=indices[np.argmax(eg.imag[band])]
  m['windows'][w]=dict(ground_max_eV=float(energy[i]),ground_height=float(eg.imag[i]),pumped_at_ground_max=float(ep.imag[i]),signed_area_ground=float(np.trapezoid(eg.imag[band],energy[band])),signed_area_pumped=float(np.trapezoid(ep.imag[band],energy[band])))
 metrics[method]=m
fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained');mask=(energy>=1.5)&(energy<=6)
for row,method in enumerate(['HSE06','TDCDFT']):
 eg,ep=spectra[method,880];axes[row,0].plot(energy[mask],eg.imag[mask],label='Unpumped');axes[row,0].plot(energy[mask],ep.imag[mask],label='Pumped')
 for w in [600,720,880]:axes[row,1].plot(energy[mask],spectra[method,w][1].imag[mask],label=f'{w*.024188843:.2f} fs')
 for a in axes[row]:a.set(title=method,xlabel='Energy (eV)',ylabel='Im epsilon');a.axhline(0,color='0.4',lw=.6);a.grid(alpha=.2);a.legend()
fig.suptitle('Si 4³: central pump–probe at 1.94 fs; fixed alpha=0.2, gamma=0.004 for TDCDFT')
fig.savefig(out/'comparison_preliminary.png',dpi=150);fig.savefig(out/'comparison_preliminary.pdf')
(out/'comparison_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(metrics,indent=2))
