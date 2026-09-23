"""Differential frozen-screening spectra. No population or binding-energy inference."""
from pathlib import Path
import sys,json,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'samples/exercise_si_tdcdft'))
from analyze import response,peak_metrics
def width_in_band(energy,eps):
 # FWHM relative to zero only if both crossings exist within the peak interval.
 mask=(energy>=2)&(energy<=4);xx=energy[mask];yy=eps.imag[mask];ip=int(np.argmax(yy));half=yy[ip]/2
 left=np.where(yy[:ip]<half)[0];right=np.where(yy[ip+1:]<half)[0]
 width=None
 if len(left) and len(right) and half>0:
  l=left[-1];rr=ip+1+right[0]
  xl=xx[l]+(half-yy[l])/(yy[l+1]-yy[l])*(xx[l+1]-xx[l])
  xr=xx[rr-1]+(half-yy[rr-1])/(yy[rr]-yy[rr-1])*(xx[rr]-xx[rr-1]);width=float(xr-xl)
 return width

fs=.02418884326505;probe=80.;work=root/'calculations/si_tdcdft_k4'
pump=np.loadtxt(work/'bare_k/strong_long/Si_rt.data');px=np.loadtxt(work/'bare_k/strong_long/Si_rt_xc.data')
energy=np.arange(.01,8.001,.01);data={};spectra={};metrics={}

# Reproduce the short full-feedback diagnostic separately from valid frozen spectra.
diag={};normalized=[]
for name,amp in [('diagnostic_full',.001),('diagnostic_half',.0005)]:
 d=work/'bare_k_probe'/name
 rr=np.loadtxt(d/'Si_rt.data');xx=np.loadtxt(d/'Si_rt_xc.data')
 jj=rr[:,15]-pump[:len(rr),15];pre=rr[:,0]<79.8;post=rr[:,0]>80.2
 diag[name]=dict(pre_current_difference=float(abs(jj[pre]).max()),alpha_before_probe=float(xx[np.where(pre)[0][-1],7]),alpha_min_after=float(xx[post,7].min()),alpha_max_after=float(xx[post,7].max()),differential_current_norm=float(np.linalg.norm(jj)/amp))
 normalized.append(jj/amp)
diag['normalized_current_relative_difference']=float(np.linalg.norm(normalized[0]-normalized[1])/np.linalg.norm(normalized[1]))
(out/'diagnostic_metrics.json').write_text(json.dumps(diag,indent=2)+'\n')

for name in ('pumped','pumped_half','pumped_small','ground_one','ground_screened'):
 d=work/'bare_k_probe'/name
 if not (d/'outputfile').exists():continue
 log=(d/'outputfile').read_text()
 if 'end SALMON' not in log:
  status=out/f'{name}_status.json'
  if status.exists():metrics[name]=dict(valid_spectrum=False,status=json.loads(status.read_text()),reason=('Confirmed norm breakdown; no peak extracted' if name=='ground_one' and 'Ne=' in log else 'Incomplete trajectory; no peak extracted'))
  continue
 r=np.loadtxt(d/'Si_rt.data');x=np.loadtxt(d/'Si_rt_xc.data')
 assert len(r)==len(x)==12000 and np.isfinite(r).all() and np.isfinite(x).all()
 np.testing.assert_allclose(r[:,0],pump[:,0],atol=1e-8,rtol=0)
 norm=[]
 for line in log.splitlines():
  cols=line.split()
  if len(cols)!=7:continue
  try:step=int(cols[0]);v=list(map(float,cols[1:]))
  except ValueError:continue
  if step>0 and abs(v[0]-step*.08*fs)<1e-6:norm.append(v[4])
 j=r[:,15]-pump[:,15] if name.startswith('pumped') else r[:,15]
 amp={'pumped_half':.0005,'pumped_small':.0001}.get(name,.001)
 pre=r[:,0]<=probe;post=r[:,0]>probe
 assert abs(j[pre]).max()<1e-9,(name,abs(j[pre]).max())
 assert np.ptp(x[post,7])==0
 if name.startswith('pumped'):np.testing.assert_allclose(x[r[:,0]>60,7],px[r[:,0]>60,7],rtol=1e-9,atol=1e-12)
 eps=response(r[:,0],j,amp,energy,probe)
 spectra[name]=eps;data[name]=(r[:,0],j/amp,x)
 np.savetxt(out/f'{name}.csv',np.column_stack([energy,eps.real,eps.imag]),delimiter=',',header='energy_eV,Re_epsilon,Im_epsilon',comments='')
 m=peak_metrics(energy,eps.imag,2,4)
 width=width_in_band(energy,eps)
 m.update(alpha=float(x[-1,7]),max_norm_error=float(np.max(abs(np.array(norm)-32))),pre_probe_max_difference=float(abs(j[pre]).max()),observation_fs=float((r[-1,0]-probe)*fs),FWHM_eV=width)
 windows={}
 for span in (600.,720.,880.):
  keep=r[:,0]<=probe+span+1e-8
  e=response(r[keep,0],j[keep],amp,energy,probe)
  windows[str(span)]=peak_metrics(energy,e.imag,2,4)
 m['window_checks_au']=windows;metrics[name]=m
if 'pumped_half' in data:
 t,j,x=data['pumped'];_,jh,_=data['pumped_half'];post=t>probe
 metrics['linearity']=dict(normalized_current_relative_l2=float(np.linalg.norm(j[post]-jh[post])/np.linalg.norm(jh[post])),spectrum_relative_l2_2_4eV=float(np.linalg.norm((spectra['pumped']-spectra['pumped_half']).imag[(energy>=2)&(energy<=4)])/np.linalg.norm(spectra['pumped_half'].imag[(energy>=2)&(energy<=4)])))
if 'pumped_small' in data:
 t,j,x=data['pumped_half'];_,js,_=data['pumped_small'];post=t>probe;mask=(energy>=2)&(energy<=4)
 metrics['small_probe_check']=dict(normalized_current_relative_l2=float(np.linalg.norm(j[post]-js[post])/np.linalg.norm(js[post])),spectrum_relative_l2_2_4eV=float(np.linalg.norm((spectra['pumped_half']-spectra['pumped_small']).imag[mask])/np.linalg.norm(spectra['pumped_small'].imag[mask])))
# Earlier fixed-alpha=.2 equilibrium response, truncated to the identical 880-au window.
old=np.loadtxt(work/'lrc/Si_rt.data');keep=old[:,0]<=880+1e-8
old_eps=response(old[keep,0],old[keep,15],.001,energy,0)
spectra['old_alpha02']=old_eps
metrics['old_alpha02']=peak_metrics(energy,old_eps.imag,2,4)
metrics['old_alpha02']['FWHM_eV']=width_in_band(energy,old_eps)
metrics['old_alpha02']['window_checks_au']={}
for span in (600.,720.,880.):
 keep_old=old[:,0]<=span+1e-8
 ee=response(old[keep_old,0],old[keep_old,15],.001,energy,0)
 metrics['old_alpha02']['window_checks_au'][str(span)]=peak_metrics(energy,ee.imag,2,4)
np.savetxt(out/'old_alpha02.csv',np.column_stack([energy,old_eps.real,old_eps.imag]),delimiter=',',header='energy_eV,Re_epsilon,Im_epsilon',comments='')
primary='pumped_small' if 'pumped_small' in spectra else 'pumped_half'
if all(n in spectra for n in ('pumped','pumped_half','ground_screened')):
 metrics['comparisons']={}
 for ref in ('ground_screened','old_alpha02'):
  m=metrics[primary];b=metrics[ref]
  metrics['comparisons'][ref]=dict(primary=primary,peak_shift_eV=m['peak_eV']-b['peak_eV'],height_change_percent=100*(m['height']/b['height']-1),area_change_percent=100*(m['area_eV']/b['area_eV']-1))
 band=(energy>=2)&(energy<=4.5)
 fig,ax=plt.subplots(2,1,figsize=(9,8),layout='constrained')
 labels={'old_alpha02':'Earlier equilibrium, fixed alpha=0.2','ground_screened':'No pump, alpha=0.0002063','pumped':'Pumped, frozen alpha=0.0002063','pumped_half':'Pumped, half probe','pumped_small':'Pumped, probe 0.0001'}
 for n in ('old_alpha02','ground_screened','pumped','pumped_half','pumped_small'):
  if n not in spectra:continue
  ax[0].plot(energy[band],spectra[n].imag[band],label=labels[n],ls='--' if n=='pumped_half' else '-')
 ax[0].set_xlim(2,4.5);ax[0].set_ylabel('Im epsilon (differential response)');ax[0].legend(fontsize=9)
 ax[1].plot(energy[band],(spectra[primary]-spectra['old_alpha02']).imag[band],label='Pump - earlier equilibrium (alpha=0.2)')
 ax[1].plot(energy[band],(spectra[primary]-spectra['ground_screened']).imag[band],label='Pump - no pump (same alpha)')
 ax[1].set_xlim(2,4.5);ax[1].set_ylabel('Delta Im epsilon');ax[1].set_xlabel('Energy (eV)');ax[1].legend(fontsize=9)
 for a in ax:a.axhline(0,color='0.5',lw=.6);a.grid(alpha=.2)
 fig.suptitle('Si 4x4x4 k: frozen-screening pump-probe\nProbe at 1.935 fs, 21.286 fs observation, pump 1e13 W/cm2')
 fig.savefig(out/'spectra.png',dpi=160);fig.savefig(out/'spectra.pdf')
 np.savetxt(out/'differential.csv',np.column_stack([energy,(spectra[primary]-spectra['old_alpha02']).imag,(spectra[primary]-spectra['ground_screened']).imag]),delimiter=',',header='energy_eV,delta_vs_alpha02,delta_vs_same_alpha',comments='')
(out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(metrics,indent=2))
