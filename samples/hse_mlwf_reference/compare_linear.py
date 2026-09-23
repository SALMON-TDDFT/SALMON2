"""Matched finite-window HSE / native TDCDFT impulse response.

No interpolation is used in Fourier transforms. Each trace keeps its native
sampling and both use exactly the same post-impulse duration and cubic window.
Transient-only mode deliberately produces no spectrum or peak measurements.
"""
import argparse
import importlib.util
import json
import os,tempfile
from pathlib import Path
import numpy as np
from checkpoint import write_json,fingerprint

_path=Path(__file__).resolve().parents[1]/'exercise_si_tdcdft/analyze.py'
_spec=importlib.util.spec_from_file_location('salmon_impulse_analysis',_path)
_analysis=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_analysis)
FS_PER_AU=.024188843265857


def uniform(time,current):
 time=np.asarray(time,float);current=np.asarray(current,float)
 if time.ndim!=1 or len(time)<3 or current.shape!=time.shape or not np.isfinite(time+current).all():
  raise ValueError('Matching finite time/current vectors required')
 dt=time[1]-time[0]
 if dt<=0 or not np.allclose(np.diff(time),dt,atol=1e-9,rtol=1e-9):raise ValueError('Uniform increasing trace required')
 return time,current,dt


def post_impulse(time,current,field):
 time,current,dt=uniform(time,current);field=np.asarray(field,float)
 if field.shape!=time.shape or not np.isfinite(field).all() or np.max(abs(field))==0:
  raise ValueError('Finite nonzero impulse field required')
 index=int(np.flatnonzero(abs(field)>.5*np.max(abs(field)))[0]);kick=float(field[index])
 if not np.allclose(field[:index],0.,atol=1e-13,rtol=0) or not np.allclose(field[index:],kick,atol=1e-13,rtol=1e-9):
  raise ValueError('Expected a single external A step, no pump or ramp')
 origin=float(time[index]-dt)
 return time[index:]-origin,current[index:],kick,origin


def same_window(t1,j1,a1,t2,j2,a2,end,energy):
 if not np.isfinite(end) or end<=0:raise ValueError('Positive finite duration required')
 result=[]
 for t,j,a in ((t1,j1,a1),(t2,j2,a2)):
  t,j,dt=uniform(t,j)
  if not np.isclose(t[0],dt,atol=1e-8,rtol=0):raise ValueError('Trace must start one step after impulse')
  keep=t<=end+1e-8
  if np.count_nonzero(keep)<3 or not np.isclose(t[keep][-1],end,atol=1e-8,rtol=0):
   raise ValueError('Both traces must contain the exact requested endpoint')
  result.append(_analysis.response(t[keep],j[keep],a,energy))
 return tuple(result)


def analyze(hse_directory,tdcdft_directory,output,end_time=880.,transient_only=False):
 hse=Path(hse_directory);td=Path(tdcdft_directory);out=Path(output);out.mkdir(parents=True,exist_ok=True)
 hs=json.loads((hse/'status.json').read_text());ts=json.loads((td/'status.json').read_text())
 if hs['status']=='failed' or not ts.get('completed') or ts.get('exit_code')!=0:
  raise ValueError('Failed or incomplete reference trajectory')
 if not transient_only and hs['status']!='completed':raise ValueError('HSE must complete before spectral comparison')
 rows=json.loads((hse/'trajectory.json').read_text())[:hs['accepted_step']]
 if len(rows)!=hs['accepted_step'] or any(r['step']!=i+1 for i,r in enumerate(rows)):
  raise ValueError('HSE history/status mismatch')
 th=np.array([r['time_au'] for r in rows]);jh=np.array([r['current'][2] for r in rows]);ah=float(hs['amplitude'])
 raw=_analysis.read_current(td/'Si_rt.data');tt,jt,at,origin=post_impulse(raw[:,0],raw[:,15],raw[:,3])
 if not np.isclose(ah,at,atol=1e-13,rtol=1e-9):raise ValueError('HSE and TDCDFT impulses must match')
 manifest=json.loads((td/'comparison.json').read_text()) if (td/'comparison.json').exists() else {}
 metadata=dict(hse_status=hs,tdcdft_status=ts,tdcdft_provenance=manifest,impulse_au=ah,
  tdcdft_probe_origin_au=origin,hse_dt_au=th[1]-th[0],tdcdft_dt_au=tt[1]-tt[0],
  energy_shift_eV=0.,additional_broadening_eV=0.,window='1-3(t/T)^2+2(t/T)^3',
  source_hashes=dict(hse=fingerprint([hse/'trajectory.json']),tdcdft=fingerprint([td/'Si_rt.data',td/'inputfile'])))
 os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'salmon-mpl-cache'))
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 duration=min(th[-1],tt[-1],end_time)
 fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
 for t,j,a,label in ((th,jh,ah,'HSE06'),(tt,jt,at,'TDCDFT (ELF-Proca)')):
  keep=t<=duration+1e-8;ax.plot(t[keep]*FS_PER_AU,j[keep]/a,label=label,lw=1)
 ax.set(xlabel='Time after impulse (fs)',ylabel='Electron current / A step (a.u.)',
        title='Si, fixed 4³ k mesh'+(' — preliminary time trace' if transient_only else ''))
 ax.legend();ax.grid(alpha=.2);fig.savefig(out/('transient_preview.png' if transient_only else 'current.png'),dpi=160);plt.close(fig)
 if transient_only:
  metadata.update(status='transient_only',available_time_fs=duration*FS_PER_AU,
                  scope='No spectrum or peak claim from this short preview')
  write_json(out/'preview.json',metadata);return metadata
 energy=np.arange(.01,8.001,.01);windows=[]
 for fraction in (.5,.75,1.):
  end=np.floor(end_time*fraction/float(hs['dt_au'])+1e-9)*float(hs['dt_au'])
  eh,et=same_window(th,jh,ah,tt,jt,at,end,energy)
  values={}
  for label,eps in (('hse',eh),('tdcdft',et)):
   values[label]={f'peak_{lo}_{hi}_eV':_analysis.peak_metrics(energy,eps.imag,lo,hi) for lo,hi in ((2,4),(4,6))}
  windows.append(dict(end_time_au=end,duration_fs=end*FS_PER_AU,metrics=values))
  if fraction==1.:
   np.savetxt(out/'spectra.csv',np.c_[energy,eh.real,eh.imag,et.real,et.imag],delimiter=',',
    header='energy_eV,HSE_Re,HSE_Im,TDCDFT_Re,TDCDFT_Im',comments='')
   fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
   ax.plot(energy,eh.imag,label='HSE06');ax.plot(energy,et.imag,label='TDCDFT (ELF-Proca)')
   ax.set(xlabel='Photon energy (eV), no shift',ylabel='Im epsilon_zz',xlim=(1,7),
    title=f'Si impulse response — common {end*FS_PER_AU:.2f} fs window')
   ax.legend();ax.grid(alpha=.2);fig.savefig(out/'spectra.png',dpi=180);fig.savefig(out/'spectra.pdf');plt.close(fig)
 metadata.update(status='completed',windows=windows,end_time_au=end_time,
  finite_window_scale_eV=2*np.pi*_analysis.HARTREE_EV/end_time,
  scope='Fixed-grid finite-window comparison; peak maxima are not proof of bound excitons or convergence')
 write_json(out/'comparison.json',metadata)
 lines=['# HSE / TDCDFT finite-window comparison','',
  f'Completed common window: {end_time*FS_PER_AU:.5f} fs. No energy shift or additional damping.',
  '','These are optical maxima on this fixed grid, not extracted exciton binding energies.',
  'HSE and TDCDFT use their own ground states; differences include their band structures.',
  '','| Window (fs) | Region (eV) | HSE maximum (eV) | TDCDFT maximum (eV) | HSE minus TDCDFT (eV) | Boundary / nonpositive maximum |',
  '|---:|---|---:|---:|---:|---|']
 for window in windows:
  for lo,hi in ((2,4),(4,6)):
   key=f'peak_{lo}_{hi}_eV';h=window['metrics']['hse'][key];t=window['metrics']['tdcdft'][key]
   flagged=h['boundary_peak'] or t['boundary_peak'] or h['height']<=0 or t['height']<=0
   lines.append(f"| {window['duration_fs']:.3f} | {lo}–{hi} | {h['peak_eV']:.2f} | {t['peak_eV']:.2f} | {h['peak_eV']-t['peak_eV']:+.2f} | {'yes' if flagged else 'no'} |")
 lines+=['','Compare the different windows before assigning a peak shift. A boundary maximum is not a resolved peak.',
  '',f"Finite-window scale 2πℏ/T: {metadata['finite_window_scale_eV']:.3f} eV; energy output spacing is 0.01 eV.",
  'No k-point, long-time or long-time-step convergence claim is made.','',
  '![Absorption comparison](spectra.png)','', '![Current comparison](current.png)']
 (out/'RESULTS.md').write_text('\n'.join(lines)+'\n')
 return metadata


if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('hse_directory');p.add_argument('tdcdft_directory');p.add_argument('output')
 p.add_argument('--end-time',type=float,default=880.);p.add_argument('--transient-only',action='store_true')
 print(json.dumps(analyze(**vars(p.parse_args())),indent=2))
