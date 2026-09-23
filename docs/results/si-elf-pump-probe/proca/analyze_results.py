"""Analyze completed ELF-Proca central pump-probe response for a gamma prefix."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;sys.path.insert(0,str(out.parent))
import analyze_probe as common
from probe_tools import central_response,apparent_width
common.work=Path('/private/tmp/salmon-si-elf-proca-probe');common.out=out
work=common.work;fs=common.fs;probe=80.;eta=1e-4;energy=np.arange(.01,8.001,.01);band=(energy>=2)&(energy<=4)
response=common.response;peak_metrics=common.peak_metrics

def main(prefix):
 data={};metrics={'gamma_prefix':prefix,'trajectories':{}}
 for kind in ('pump','plus','minus','halfplus','halfminus','ground'):
  r,x,s=common.load(prefix+'_'+kind);data[kind]=(r,x);metrics['trajectories'][kind]=s
  if not s['valid_trajectory']:raise RuntimeError(f'{prefix}_{kind}: invalid; no spectrum')
 pump,px=data['pump'];t=pump[:,0];pre=t<=probe;post=t>probe
 for kind,(r,x) in data.items():
  np.testing.assert_allclose(r[:,0],t,atol=1e-8,rtol=0)
  if kind!='ground':
   np.testing.assert_allclose(r[pre],pump[pre],atol=1e-12,rtol=1e-9)
   np.testing.assert_allclose(x[pre],px[pre],atol=1e-12,rtol=1e-9)
 for kind,amp in [('plus',eta),('minus',-eta),('halfplus',eta/2),('halfminus',-eta/2)]:
  np.testing.assert_allclose(data[kind][0][post,3]-pump[post,3],amp,atol=1e-14,rtol=0)
 assert np.max(abs(data['ground'][0][pre,15]))<1e-9
 normalized={'pumped':central_response(t,data['plus'][0][:,15],t,data['minus'][0][:,15],eta),
  'pumped_half':central_response(t,data['halfplus'][0][:,15],t,data['halfminus'][0][:,15],eta/2),
  'one_sided':(data['plus'][0][:,15]-pump[:,15])/eta,
  'ground':data['ground'][0][:,15]/eta}
 spectra={}
 for name,j in normalized.items():
  eps=response(t,j,1.,energy,probe);spectra[name]=eps
  m=peak_metrics(energy,eps.imag,2,4);m['apparent_FWHM_eV']=apparent_width(energy,eps.imag);m['tracked_2p8_3p8eV']=peak_metrics(energy,eps.imag,2.8,3.8);m['windows']={}
  for span in (600.,720.,880.):
   keep=t<=probe+span+1e-8;ee=response(t[keep],j[keep],1.,energy,probe);m['windows'][str(span)]=peak_metrics(energy,ee.imag,2,4);m['windows'][str(span)]['tracked_2p8_3p8eV']=peak_metrics(energy,ee.imag,2.8,3.8)
  metrics[name]=m
  np.savetxt(out/f'{prefix}_{name}_spectrum.csv',np.column_stack((energy,eps.real,eps.imag)),delimiter=',',header='energy_eV,Re_epsilon,Im_epsilon',comments='')
 def compare(a,b):
  return dict(normalized_current_relative_l2=float(np.linalg.norm((normalized[a]-normalized[b])[post])/np.linalg.norm(normalized[b][post])),
   spectrum_relative_l2_2_4eV=float(np.linalg.norm((spectra[a]-spectra[b]).imag[band])/np.linalg.norm(spectra[b].imag[band])),
   peak_shift_eV=metrics[a]['peak_eV']-metrics[b]['peak_eV'],height_change_percent=100*(metrics[a]['height']/metrics[b]['height']-1),
   area_change_percent=100*(metrics[a]['area_eV']/metrics[b]['area_eV']-1))
 metrics['probe_linearity']=compare('pumped','pumped_half');metrics['one_sided_check']=compare('one_sided','pumped_half');metrics['pump_effect']=compare('pumped_half','ground')
 metrics['alpha_response']=dict(max_plus_minus=float(abs(data['plus'][1][:,7]-data['minus'][1][:,7]).max()),max_half_plus_minus=float(abs(data['halfplus'][1][:,7]-data['halfminus'][1][:,7]).max()))
 metrics['observation_fs']=float((t[-1]-probe)*fs)
 np.savez_compressed(out/f'{prefix}_differential_currents.npz',time_au=t,**normalized)
 np.savetxt(out/f'{prefix}_difference_spectrum.csv',np.column_stack((energy,(spectra['pumped_half']-spectra['ground']).imag)),delimiter=',',header='energy_eV,delta_Im_epsilon',comments='')
 (out/f'{prefix}_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
 gamma=metrics['trajectories']['pump']['gamma']
 fig,ax=plt.subplots(2,1,figsize=(9,7),sharex=True,layout='constrained');show=(energy>=2)&(energy<=4.5)
 for name,label,style in [('ground','Unpumped, same beta/gamma','-'),('pumped','Pumped, central probe 1e-4','--'),('pumped_half','Pumped, central probe 5e-5','-')]:
  ax[0].plot(energy[show],spectra[name].imag[show],style,label=label)
 ax[0].set(ylabel='Im epsilon (differential)');ax[0].legend(fontsize=9)
 ax[1].plot(energy[show],(spectra['pumped_half']-spectra['ground']).imag[show],label='Pumped - unpumped');ax[1].set(xlabel='Energy (eV)',ylabel='Delta Im epsilon');ax[1].legend()
 for a in ax:a.set_xlim(2,4.5);a.axhline(0,color='0.5',lw=.7);a.grid(alpha=.2)
 fig.suptitle(f'Si 4x4x4 k: dynamic ELF alpha + Proca, beta=0, gamma={gamma:g}\nProbe at1.935 fs;21.286 fs observation; pump1e13 W/cm2')
 fig.savefig(out/f'{prefix}_spectra.png',dpi=160);fig.savefig(out/f'{prefix}_spectra.pdf')
 fig,ax=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
 for kind in ('pump','plus','minus','ground'):
  r,x=data[kind];ax[0].plot(t*fs,x[:,7],label=kind);ax[1].plot(t*fs,x[:,3])
 for name in ('ground','pumped','pumped_half'):ax[2].plot(t*fs,normalized[name],label=name)
 ax[0].set(ylabel='alpha');ax[0].legend();ax[1].set(ylabel='Axc,z / c (a.u.)');ax[2].set(ylabel='Delta Jz / probe',xlabel='Time (fs)');ax[2].legend()
 for a in ax:a.axvline(probe*fs,color='0.5',ls=':');a.grid(alpha=.2)
 fig.savefig(out/f'{prefix}_dynamics.png',dpi=160)
 fig,ax=plt.subplots(2,1,figsize=(9,7),sharex=True,layout='constrained')
 for span in (600.,720.,880.):
  keep=t<=probe+span+1e-8
  for axis,name in zip(ax,('ground','pumped_half')):
   ee=response(t[keep],normalized[name][keep],1.,energy,probe)
   axis.plot(energy[show],ee.imag[show],label=f'{span*fs:.2f} fs')
 for axis,name in zip(ax,('Unpumped','Pumped, central probe 5e-5')):
  axis.set(ylabel='Im epsilon',title=name,xlim=(2,4.5));axis.legend(title='Observation window');axis.axhline(0,color='0.5',lw=.7);axis.grid(alpha=.2)
 ax[-1].set_xlabel('Energy (eV)');fig.savefig(out/f'{prefix}_windows.png',dpi=160)
 print(json.dumps(metrics,indent=2))
if __name__=='__main__':main(sys.argv[1])
