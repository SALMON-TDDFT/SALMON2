"""Summarize completed bare-K trajectories; retain the absolute gate for a controlled comparison."""
from pathlib import Path
import json,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
work=root/'calculations/si_tdcdft_k4';fs=.02418884326505
metrics={};data={}
for name in ('weak','strong','strong_long','strong_matched_gate'):
 d=work/'bare_k'/name
 if not (d/'Si_rt_xc.data').exists(): continue
 log=(d/'outputfile').read_text()
 if 'end SALMON' not in log: continue
 x=np.loadtxt(d/'Si_rt_xc.data');r=np.loadtxt(d/'Si_rt.data')
 assert np.isfinite(x).all() and np.isfinite(r).all()
 np.testing.assert_allclose(x[:,0],r[:,0],rtol=0,atol=1e-8)
 complete='end SALMON' in log
 if not complete: continue
 assert len(x)==(12000 if name=='strong_long' else 1200)
 post=x[:,0]>62
 end=x[:,0]>max(62,x[-1,0]-3/fs)
 energy=np.loadtxt(d/'Si_rt_energy.data')
 norm=[]
 for line in log.splitlines():
  a=line.split()
  if len(a)!=7:continue
  try: step=int(a[0]);v=list(map(float,a[1:]))
  except ValueError:continue
  if step>0 and abs(v[0]-step*.08*fs)<1e-6:norm.append(v[4])
 expected=1/(1+4*np.pi*np.maximum(x[:,8],0)/.2**2)
 # Before the first valid field, alpha is initialized to one, consistent with K=0.
 np.testing.assert_allclose(x[:,7],expected,rtol=1e-10,atol=1e-12)
 metrics[name]=dict(final_window_mean_E=float(x[end,6].mean()),final_window_mean_J=float(r[end,15].mean()),electronic_energy_change_after_pulse=float(energy[-1,1]-energy[np.argmin(abs(energy[:,0]-60)),1]),steps=len(x),end_fs=float(x[-1,0]*fs),alpha_min=float(x[:,7].min()),alpha_max=float(x[:,7].max()),alpha_final=float(x[-1,7]),K_min=float(x[:,8].min()),K_max=float(x[:,8].max()),A_final=float(x[-1,3]),E_final=float(x[-1,6]),J_final=float(r[-1,15]),max_abs_A=float(abs(x[:,3]).max()),max_norm_error=float(np.max(abs(np.array(norm)-32))),post_offset_max=float(abs(x[post,6]-x[post,7]*x[post,11]).max()))
 data[name]=(x,r)
if 'weak' in data and 'strong' in data:
 wx,wr=data['weak'];sx,sr=data['strong'];scaled=wr[:,15]*np.sqrt(1e5)
 metrics['comparison']=dict(strong_vs_scaled_weak_current_relative_l2=float(np.linalg.norm(sr[:,15]-scaled)/np.linalg.norm(scaled)))
if 'weak' in data:
 wx,wr=data['weak'];a=wr[:,1:4]
 padded=np.vstack([np.zeros((1,3)),a,np.zeros((1,3))]);e=-(padded[2:]-padded[:-2])/.16
 denom=np.sum(a*a+(e/.2)**2,axis=1);num=np.sum(a*wr[:,13:16]-e*wx[:,9:12],axis=1)
 k=np.divide(num,denom,out=np.zeros_like(num),where=denom>0)
 valid=denom*1e5>(2e-5)**2;last=np.flatnonzero(valid)[-1]
 metrics['comparison']['scaled_linear_weak_offline_final_alpha']=float(1/(1+4*np.pi*max(k[last],0)/.04))
if 'strong_matched_gate' in data:
 mx,mr=data['strong_matched_gate'];wx,wr=data['weak']
 metrics['comparison']['matched_vs_absolute_current_relative_l2']=float(np.linalg.norm(mr[:,13:16]-data['strong'][1][:,13:16])/np.linalg.norm(data['strong'][1][:,13:16]))
 metrics['comparison']['matched_gate_current_relative_l2']=float(np.linalg.norm(mr[:,15]-wr[:,15]*np.sqrt(1e5))/np.linalg.norm(wr[:,15]*np.sqrt(1e5)))
if 'strong_long' in data:
 x,r=data['strong_long']
 np.testing.assert_allclose(x[:1200],data['strong'][0],rtol=0,atol=1e-12)
 old=np.loadtxt(work/'polarization/long/Si_rt.data');ox=np.loadtxt(work/'polarization/long/Si_rt_xc.data')
 metrics['comparison']['previous_alpha_final']=float(ox[-1,7])
 metrics['comparison']['previous_A_final']=float(ox[-1,3])
 metrics['comparison']['previous_E_final']=float(ox[-1,6])
 fig,ax=plt.subplots(3,1,figsize=(9,8),layout='constrained')
 for name in ('weak','strong','strong_matched_gate'):
  xx,rr=data[name];ax[0].plot(xx[:,0]*fs,xx[:,7],label=name)
 ax[0].set_ylabel('alpha');ax[0].set_xlabel('Time (fs)');ax[0].legend()
 ax[1].plot(x[:,0]*fs,x[:,3],label='alpha0=1, K0=0');ax[1].plot(ox[:,0]*fs,ox[:,3],label='Previous alpha0=0.2, K0=0.19591');ax[1].set_ylabel('XC A/c (a.u.)');ax[1].legend()
 ax[2].plot(r[:,0]*fs,r[:,15],label='alpha0=1, K0=0');ax[2].plot(old[:,0]*fs,old[:,15],alpha=.7,label='Previous');ax[2].set_ylabel('Number current (a.u.)');ax[2].set_xlabel('Time (fs)')
 for a in ax:a.axvline(60*fs,color='0.5',ls=':');a.grid(alpha=.2)
 fig.suptitle('Si 4x4x4 k: bare-K trial, instantaneous polarization closure')
 fig.savefig(out/'comparison.png',dpi=150);fig.savefig(out/'comparison.pdf')
(out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
print(json.dumps(metrics,indent=2))
