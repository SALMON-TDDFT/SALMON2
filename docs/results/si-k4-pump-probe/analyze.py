"""Inspect completed trajectories; only form spectra from complete finite runs."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
from pathlib import Path
import json,re,sys
import numpy as np
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'samples/exercise_si_tdcdft'))
from analyze import response
work=root/'calculations/si_hse_native/k4_pump_probe';summary={};arrays={}
for d in sorted(work.glob('td_*')):
 if not (d/'status.json').exists():continue
 s=json.loads((d/'status.json').read_text());rows=[]
 dt=float(re.search(r'(?m)^\s*dt\s*=\s*([0-9.]+)',(d/'inputfile').read_text()).group(1))
 for l in (d/'run.log').read_text().splitlines():
  c=l.split()
  if len(c)!=7:continue
  try:
   step=int(c[0]);v=list(map(float,c[1:]))
   if step>0 and abs(v[0]-step*dt*.02418884326505)<1e-6:rows.append(v)
  except ValueError:pass
 if rows:s.update(last_time_fs=rows[-1][0],max_electron_error=max(abs(v[4]-32) for v in rows))
 if (d/'Si_rt.data').exists():
  r=np.atleast_2d(np.loadtxt(d/'Si_rt.data'));arrays[d.name]=r
  s.update(rows=len(r),finite=bool(np.isfinite(r).all()),max_abs_J=float(np.max(abs(r[:,15]))))
  x=np.atleast_2d(np.loadtxt(d/'Si_rt_xc.data'))
  s['max_abs_Axc']=float(np.max(abs(x[:,3])))
  s['finite']=s['finite'] and bool(np.isfinite(x).all())
  s['usable']=bool(s['completed'] and s['finite'] and s.get('max_electron_error',1)<1e-4 and len(r)==round(960/dt))
 summary[d.name]=s
(out/'trajectory_status.json').write_text(json.dumps(summary,indent=2)+'\n')
energy=np.arange(.01,8.001,.01)
for g in ['g10','g40']:
 names=[f'td_{g}_{k}' for k in ['pump','plus','minus','ground']]
 if not all(summary.get(n,{}).get('usable') for n in names):continue
 pump,plus,minus,ground=[arrays[n] for n in names];t=pump[:,0];post=t>80;pre=~post
 for a in [plus,minus,ground]:np.testing.assert_allclose(a[:,0],t,atol=1e-8,rtol=0)
 np.testing.assert_allclose(plus[pre],pump[pre],atol=1e-10,rtol=1e-8)
 np.testing.assert_allclose(minus[pre],pump[pre],atol=1e-10,rtol=1e-8)
 np.testing.assert_allclose(plus[post,3]-pump[post,3],1e-4,atol=1e-13,rtol=0)
 np.testing.assert_allclose(minus[post,3]-pump[post,3],-1e-4,atol=1e-13,rtol=0)
 # response expects signed step and matter current. Time is relative to probe.
 for duration in [600,720,880]:
  mask=post&(t<=80+duration+1e-8)
  ep=response(t[mask]-80,(plus[mask,15]-minus[mask,15])/2,1e-4,energy)
  eg=response(t[mask]-80,ground[mask,15],1e-4,energy)
  np.savetxt(out/f'{g}_{duration}au.csv',np.column_stack([energy,eg.real,eg.imag,ep.real,ep.imag]),delimiter=',',header='eV,ground_re,ground_im,pumped_re,pumped_im')
print(json.dumps(summary,indent=2))
