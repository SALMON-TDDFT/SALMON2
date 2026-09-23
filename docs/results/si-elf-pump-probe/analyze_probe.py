"""Full-feedback, finite-window Si pump-probe response. No scissor shift."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
from pathlib import Path
import sys,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from probe_tools import central_response,apparent_width
out=Path(__file__).resolve().parent;root=out.parents[2]
sys.path.insert(0,str(root/'samples/exercise_si_tdcdft'))
from analyze import response,peak_metrics
work=Path('/private/tmp/salmon-si-elf-pump-probe');fs=.02418884326505;probe=80.;eta=1e-4
energy=np.arange(.01,8.001,.01);band=(energy>=2)&(energy<=4)

def load(name,expected_steps=12000):
 p=work/name;status=json.loads((p/'status.json').read_text())
 r=np.atleast_2d(np.loadtxt(p/'Si_rt.data'));x=np.atleast_2d(np.loadtxt(p/'Si_rt_xc.data'))
 norms=[];norm_times=[];dt=float(r[1,0]-r[0,0])
 for line in (p/'outputfile').read_text().splitlines():
  c=line.split()
  if len(c)!=7:continue
  try:step=int(c[0]);v=list(map(float,c[1:]))
  except ValueError:continue
  if step>0 and abs(v[0]-step*dt*fs)<1e-6:norms.append(v[4]);norm_times.append(v[0])
 status.update(rows=len(r),end_time_fs=float(r[-1,0]*fs),max_norm_error=float(np.max(abs(np.array(norms)-32))) if norms else None,
  alpha_min=float(x[:,7].min()),alpha_max=float(x[:,7].max()),alpha_final=float(x[-1,7]),max_abs_Axc=float(abs(x[:,1:4]).max()),max_abs_J=float(abs(r[:,13:16]).max()))
 status['valid_trajectory']=bool(status['completed'] and len(r)==len(x)==expected_steps and np.isfinite(r).all() and np.isfinite(x).all() and norms and status['max_norm_error']<1e-4)
 bad=np.where(abs(np.array(norms)-32)>1e-4)[0]
 status['first_logged_norm_error_over_1e_4_fs']=float(norm_times[bad[0]]) if len(bad) else None
 log=(p/'outputfile').read_text();failed_steps=re.findall(r'\bt=\s*(\d+)',log)
 status['failure_step']=int(failed_steps[-1]) if failed_steps and not status['completed'] else None
 status['failure_time_fs']=status['failure_step']*dt*fs if status['failure_step'] is not None else None
 status['fatal_electron_counts']=sorted(set(float(v) for v in re.findall(r'Ne=\s*([0-9.Ee+\-]+)',log)))
 np.savez_compressed(out/f'{name}_trace.npz',rt=r,xc=x,norm_time_fs=np.array(norm_times),electron_count=np.array(norms))
 return r,x,status

def main():
 names=[n for n in ('pump','plus','minus','half_plus','half_minus','ground','ground_half','ground_zero') if (work/n/'status.json').exists()]
 data={};metrics={'trajectories':{}}
 for n in names:
  r,x,s=load(n);data[n]=(r,x);metrics['trajectories'][n]=s
 (out/'trajectory_status.json').write_text(json.dumps(metrics,indent=2)+'\n')
 if not all(n in data and metrics['trajectories'][n]['valid_trajectory'] for n in ('pump','plus','minus','half_plus','half_minus','ground')):
  # Explicitly avoid a spectrum from incomplete or failed full-feedback runs.
  fig,ax=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
  for n in names:
   r,x=data[n];t=r[:,0]*fs
   ax[0].plot(t,x[:len(t),7],label=n);ax[1].plot(t,x[:len(t),3]);ax[2].plot(t,r[:,15])
  ax[0].set(ylabel='alpha');ax[0].legend();ax[1].set(ylabel='Axc,z / c (a.u.)');ax[2].set(ylabel='Jz (a.u.)',xlabel='Time (fs)')
  for a in ax:a.axvline(probe*fs,color='0.5',ls=':');a.grid(alpha=.2)
  fig.suptitle('Full ELF-feedback trajectories: spectrum validation incomplete')
  fig.savefig(out/'trajectory_diagnostic.png',dpi=160)
  print(json.dumps(metrics,indent=2));return
 pump,px=data['pump'];t=pump[:,0];pre=t<=probe;post=t>probe
 for n,(r,x) in data.items():
  if not metrics['trajectories'][n]['valid_trajectory']:continue
  np.testing.assert_allclose(r[:,0],t,atol=1e-8,rtol=0)
  if not n.startswith('ground'):
   np.testing.assert_allclose(r[pre],pump[pre],atol=1e-12,rtol=1e-9)
   np.testing.assert_allclose(x[pre],px[pre],atol=1e-12,rtol=1e-9)
 # Check the signed external step; do not infer its sign from a name alone.
 for n,a in [('plus',eta),('minus',-eta),('half_plus',eta/2),('half_minus',-eta/2)]:
  np.testing.assert_allclose(data[n][0][post,3]-pump[post,3],a,atol=1e-14,rtol=0)
 normalized={
  'pumped':central_response(t,data['plus'][0][:,15],t,data['minus'][0][:,15],eta),
  'pumped_half':central_response(t,data['half_plus'][0][:,15],t,data['half_minus'][0][:,15],eta/2),
  'one_sided':(data['plus'][0][:,15]-pump[:,15])/eta,
  'ground':data['ground'][0][:,15]/eta,
 }
 if 'ground_zero' in data and metrics['trajectories']['ground_zero']['valid_trajectory']:
  normalized['ground']=(data['ground'][0][:,15]-data['ground_zero'][0][:,15])/eta
 if 'ground_half' in data and metrics['trajectories']['ground_half']['valid_trajectory']:normalized['ground_half']=data['ground_half'][0][:,15]/(eta/2)
 spectra={}
 for name,j in normalized.items():
  assert np.max(abs(j[pre]))<1e-5,(name,'nonzero response before probe')
  eps=response(t,j,1.,energy,probe);spectra[name]=eps
  m=peak_metrics(energy,eps.imag,2,4);m['apparent_FWHM_eV']=apparent_width(energy,eps.imag)
  m['windows']={}
  for span in (600.,720.,880.):
   keep=t<=probe+span+1e-8;ee=response(t[keep],j[keep],1.,energy,probe)
   m['windows'][str(span)]=peak_metrics(energy,ee.imag,2,4)
  metrics[name]=m
  np.savetxt(out/f'{name}_spectrum.csv',np.column_stack((energy,eps.real,eps.imag)),delimiter=',',header='energy_eV,Re_epsilon,Im_epsilon',comments='')
 def difference_metrics(a,b):
  return dict(normalized_current_relative_l2=float(np.linalg.norm((normalized[a]-normalized[b])[post])/np.linalg.norm(normalized[b][post])),
   spectrum_relative_l2_2_4eV=float(np.linalg.norm((spectra[a]-spectra[b]).imag[band])/np.linalg.norm(spectra[b].imag[band])),
   peak_shift_eV=metrics[a]['peak_eV']-metrics[b]['peak_eV'],height_change_percent=100*(metrics[a]['height']/metrics[b]['height']-1),
   area_change_percent=100*(metrics[a]['area_eV']/metrics[b]['area_eV']-1))
 metrics['probe_linearity']=difference_metrics('pumped','pumped_half')
 metrics['one_sided_check']=difference_metrics('one_sided','pumped_half')
 metrics['pump_effect']=difference_metrics('pumped_half','ground')
 if 'ground_half' in normalized:metrics['ground_linearity']=difference_metrics('ground','ground_half')
 metrics['alpha_probe_response']=dict(max_plus_minus=float(abs(data['plus'][1][:,7]-data['minus'][1][:,7]).max()),
  max_half_plus_minus=float(abs(data['half_plus'][1][:,7]-data['half_minus'][1][:,7]).max()))
 np.savez_compressed(out/'differential_currents.npz',time_au=t,**normalized)
 (out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
 fig,ax=plt.subplots(2,1,figsize=(9,7.5),sharex=True,layout='constrained');show=(energy>=2)&(energy<=4.5)
 for n,label,style in [('ground','No pump, dynamic ELF','-'),('pumped','Pumped, central probe 1e-4','--'),('pumped_half','Pumped, central probe 5e-5','-')]:
  ax[0].plot(energy[show],spectra[n].imag[show],style,label=label)
 ax[0].set(ylabel='Im epsilon (differential)');ax[0].legend()
 ax[1].plot(energy[show],(spectra['pumped_half']-spectra['ground']).imag[show],label='Pump - no pump')
 ax[1].set(xlabel='Energy (eV)',ylabel='Delta Im epsilon');ax[1].legend()
 for a in ax:a.set_xlim(2,4.5);a.axhline(0,color='0.5',lw=.7);a.grid(alpha=.2)
 fig.suptitle('Si 4x4x4 k: full ELF-alpha pump-probe\nProbe at1.935 fs;21.286 fs observation; pump1e13 W/cm2')
 fig.savefig(out/'spectra.png',dpi=160);fig.savefig(out/'spectra.pdf')
 fig,ax=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
 for n in ('pump','plus','minus','ground'):
  r,x=data[n];ax[0].plot(t*fs,x[:,7],label=n);ax[1].plot(t*fs,x[:,3])
 for n in ('ground','pumped','pumped_half'):ax[2].plot(t*fs,normalized[n],label=n)
 ax[0].set(ylabel='alpha');ax[0].legend();ax[1].set(ylabel='Axc,z / c (a.u.)');ax[2].set(ylabel='Delta Jz / probe',xlabel='Time (fs)');ax[2].legend()
 for a in ax:a.axvline(probe*fs,color='0.5',ls=':');a.grid(alpha=.2)
 fig.savefig(out/'dynamics.png',dpi=160);fig.savefig(out/'dynamics.pdf')
 print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
