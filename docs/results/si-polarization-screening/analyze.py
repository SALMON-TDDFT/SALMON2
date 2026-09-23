"""Compare the old acceleration closure and polarization closure, with a half-dt check."""
from pathlib import Path
import json,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[3]
out=Path(__file__).resolve().parent
work=root/'calculations/si_tdcdft_k4'
fs=.02418884326505
paths={'old':work/'instant_long/strong_screened','polarization':work/'polarization/long',
       'fine':work/'polarization/fine'}
data={};metrics={}
for name,path in paths.items():
    x=np.loadtxt(path/'Si_rt_xc.data');j=np.loadtxt(path/'Si_rt.data')
    energy=np.loadtxt(path/'Si_rt_energy.data')
    assert len(x)==len(j) and np.isfinite(x).all() and np.isfinite(j).all()
    assert len(x)==(6000 if name=='fine' else 12000)
    assert 'end SALMON' in (path/'outputfile').read_text()
    np.testing.assert_allclose(x[:,0],j[:,0],atol=1e-8,rtol=0)
    post=x[:,0]>62;end=x[:,0]>max(62,x[-1,0]-3/fs)
    np.testing.assert_allclose(j[post,1:13],0,atol=1e-12,rtol=0)
    offset=x[post,6]-x[post,7]*x[post,11]
    norm=[]
    dt=.04 if name=='fine' else .08
    for line in (path/'outputfile').read_text().splitlines():
        a=line.split()
        if len(a)!=7: continue
        try:
            step=int(a[0]);v=list(map(float,a[1:]))
        except ValueError: continue
        if step>0 and abs(v[0]-step*dt*fs)<1e-6: norm.append(v[4])
    metrics[name]=dict(electronic_energy_change_after_pulse_au=float(energy[-1,1]-energy[np.argmin(abs(energy[:,0]-60)),1]),
                      final_time_fs=float(x[-1,0]*fs),A_final=float(x[-1,3]),
                      Exc_final=float(x[-1,6]),J_final=float(j[-1,15]),alpha_final=float(x[-1,7]),
                      offset_max_abs=float(np.max(abs(offset))),
                      max_norm_error=float(np.max(abs(np.array(norm)-32))),
                      post_Exc_rms=float(np.sqrt(np.mean(x[post,6]**2))),
                      final_window_mean_J=float(j[end,15].mean()),
                      final_window_rms_J=float(np.sqrt(np.mean(j[end,15]**2))),
                      final_window_mean_Exc=float(x[end,6].mean()),
                      final_window_A_slope_au=float(np.polyfit(x[end,0],x[end,3],1)[0]))
    data[name]=(x,j)
coarse,cj=data['polarization'];fine,fj=data['fine']
fine=fine[1::2];fj=fj[1::2];coarse=coarse[:len(fine)];cj=cj[:len(fine)]
np.testing.assert_allclose(coarse[:,0],fine[:,0],atol=1e-12,rtol=0)
metrics['timestep_comparison']=dict(overlap_time_fs=float(fine[-1,0]*fs),
    relative_current_l2=float(np.linalg.norm(fj[:,13:16]-cj[:,13:16])/np.linalg.norm(fj[:,13:16])),
    relative_postpulse_current_l2=float(np.linalg.norm(fj[fine[:,0]>60,13:16]-cj[fine[:,0]>60,13:16])/np.linalg.norm(fj[fine[:,0]>60,13:16])),
    max_A_difference=float(np.max(abs(fine[:,1:4]-coarse[:,1:4]))),
    max_Exc_difference=float(np.max(abs(fine[:,4:7]-coarse[:,4:7]))),
    coarse_alpha_at_overlap=float(coarse[-1,7]),fine_alpha_at_overlap=float(fine[-1,7]))
(out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
fig,axes=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
for name,label,color in [('old',"Previous: a''=alpha J",'tab:orange'),
                         ('polarization',"New: E_xc=alpha P",'tab:blue')]:
    x,j=data[name];t=x[:,0]*fs
    axes[0].plot(t,x[:,3],label=label,color=color)
    # Show post-pulse electric field separately so the pulse itself does not hide the residual.
    mask=x[:,0]>60
    axes[1].plot(t[mask],x[mask,6],color=color)
    axes[2].plot(t,j[:,15],color=color,lw=.8)
axes[0].legend();axes[0].set_ylabel('XC A/c z (a.u.)')
axes[1].set_ylabel('Post-pulse XC E z (a.u.)');axes[2].set_ylabel('Number current z (a.u.)')
for a in axes:
    a.axvline(60*fs,color='0.5',ls=':');a.axhline(0,color='0.7',lw=.6);a.grid(alpha=.2)
axes[2].set_xlabel('Time (fs)')
fig.suptitle('Si 4x4x4 k: polarization-consistent screening\nSame pulse and estimator; beta=gamma=0')
fig.savefig(out/'comparison.png',dpi=160);fig.savefig(out/'comparison.pdf')
print(json.dumps(metrics,indent=2))
