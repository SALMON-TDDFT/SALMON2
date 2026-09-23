"""MLWF warm starts every ten RT steps; existing alpha is logged, not redefined."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wannier import read_u,geometry,initial_gauge,reconstruct,moments
from mlwf import overlap_mesh,functional,minimize,transport_gauge
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
work=Path('/private/tmp/salmon-si-time-wannier-dense');fs=.02418884326505;ang=.52917721067
r,c,k,b,L,h=geometry(root);dv=h**3;u0,_=read_u(root/'calculations/si_tdcdft_k4/gs/data_for_restart')
g0=np.load(out/'mlwf_initial.npz')['gauge'];steps=np.arange(0,1601,10);traces={};summary={}
for case in ('none','weak','strong'):
 status=json.loads((out/f'{case}_dense_status.json').read_text());assert status['complete']
 rt=np.loadtxt(work/case/'Si_rt_xc.data');g=g0.copy();records=[];previous_v=None;previous_centers=None;history=[]
 for step in steps:
  u=u0 if step==0 else read_u(work/case/f'checkpoint_rt_{step:06d}')[0]
  temporal_min=1.
  if previous_v is not None:
   g,sv=transport_gauge(u,previous_v,dv);temporal_min=float(sv.min())
  raw,neighbors,bv,weights=overlap_mesh(u,k,r,L,dv)
  g,log=minimize(g,raw,neighbors,bv,weights,maxiter=600,tolerance=1e-6)
  if not log['converged']:raise RuntimeError(f'MLWF not converged: {case}, {step}, {log["gradient_norm"]}')
  spread,grad,info=functional(g,raw,neighbors,bv,weights);v=u@g;centers=info['centers']
  if previous_v is not None:
   ov=np.einsum('krn,krm->nm',previous_v.conj(),v)*dv/len(k)
   diag=np.diag(abs(ov));label_ok=bool(np.all(np.argmax(abs(ov),axis=1)==np.arange(16)))
   center_delta=(centers-previous_centers+L/2)%L-L/2
   max_center_step=float(np.max(np.linalg.norm(center_delta,axis=1))*ang)
  else:diag=np.ones(16);label_ok=True;max_center_step=0.
  assert label_ok,(case,step,'orbital matching changed')
  alpha=1. if step==0 else float(rt[step-1,7])
  record=dict(step=int(step),time_fs=float(step*.08*fs),mean_spread_A2=float(spread/16*ang**2),
    mean_invariant_A2=float(info['omega_invariant']/16*ang**2),iterations=log['iterations'],gradient_norm=log['gradient_norm'],
    initial_spread_A2=float(log['initial_spread']/16*ang**2),min_orbital_overlap=float(diag.min()),max_center_step_A=max_center_step,
    existing_alpha=alpha,temporal_min_singular=temporal_min)
  records.append(record);history.append(centers);previous_v=v;previous_centers=centers
  if step%400==0:print(case,step,record['mean_spread_A2'],log['iterations'],flush=True)
 traces[case]=records
 np.savez_compressed(out/f'{case}_mlwf.npz',steps=steps,centers_bohr=np.array(history),final_gauge=g)
 (out/f'{case}_mlwf.json').write_text(json.dumps(records,indent=2)+'\n')
 summary[case]=dict(initial_mean_spread_A2=records[0]['mean_spread_A2'],final_mean_spread_A2=records[-1]['mean_spread_A2'],
    final_mean_invariant_A2=records[-1]['mean_invariant_A2'],max_iterations=max(z['iterations'] for z in records),
    mean_iterations=float(np.mean([z['iterations'] for z in records[1:]])),max_gradient_norm=max(z['gradient_norm'] for z in records),
    min_orbital_overlap=min(z['min_orbital_overlap'] for z in records),max_center_step_A=max(z['max_center_step_A'] for z in records),
    min_temporal_singular=min(z['temporal_min_singular'] for z in records),snapshots=len(records),all_converged=True,existing_alpha_final=records[-1]['existing_alpha'])
# Cold-start comparison at the final strong snapshot, to check sensitivity to warm starts.
u,_=read_u(work/'strong'/'checkpoint_rt_001600');raw,nb,bv,wt=overlap_mesh(u,k,r,L,dv)
g,log=minimize(g0,raw,nb,bv,wt,maxiter=1500,tolerance=1e-6)
summary['final_cold_start']=dict(converged=log['converged'],iterations=log['iterations'],mean_spread_A2=log['spread']/16*ang**2,
                              difference_from_warm_A2=log['spread']/16*ang**2-summary['strong']['final_mean_spread_A2'])
summary['initialization']='Previous U transported into current occupied basis via polar overlap, then spread minimized.'
summary['alpha_policy']='Original alpha logged every10 steps. No Wannier-derived alpha formula or feedback has been introduced.'
summary['sampling']=dict(dt_au=.08,stride=10,interval_fs=.8*fs,observation_fs=128*fs)
(out/'mlwf_metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,ax=plt.subplots(3,1,figsize=(9,9),sharex=True,layout='constrained');colors={'none':'0.4','weak':'tab:blue','strong':'tab:red'}
for case,records in traces.items():
 t=[z['time_fs'] for z in records]
 ax[0].plot(t,[z['mean_spread_A2'] for z in records],label=case,color=colors[case])
 ax[1].plot(t,[z['mean_invariant_A2'] for z in records],color=colors[case])
 ax[2].plot(t,[z['iterations'] for z in records],color=colors[case])
ax[0].set_ylabel('MLWF mean spread (A^2)');ax[0].legend();ax[1].set_ylabel('Gauge-invariant part (A^2)');ax[2].set_ylabel('Localization iterations');ax[2].set_xlabel('Time (fs)')
for a in ax:a.axvline(60*fs,color='0.5',ls=':');a.grid(alpha=.2)
fig.suptitle('Si 4x4x4 k: MLWF updated every10 RT steps\nPrevious U reused; alpha feedback unchanged')
fig.savefig(out/'mlwf_dynamics.png',dpi=150);fig.savefig(out/'mlwf_dynamics.pdf')
print(json.dumps(summary,indent=2),flush=True)
