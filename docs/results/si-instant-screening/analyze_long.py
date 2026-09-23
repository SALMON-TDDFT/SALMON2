"""Analyze completed fixed/instantaneous long runs, including field drift and current envelopes."""
from pathlib import Path
import json, os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[3]
work=root/'calculations/si_tdcdft_k4/instant_long'
out=Path(__file__).resolve().parent
fs=.02418884326505
results={}; data={}
for name in ['strong_fixed','strong_screened']:
    d=work/name
    x=np.loadtxt(d/'Si_rt_xc.data');j=np.loadtxt(d/'Si_rt.data');energy=np.loadtxt(d/'Si_rt_energy.data')
    assert len(x)==len(j) and np.isfinite(x).all() and np.isfinite(j).all()
    completed='end SALMON' in (d/'outputfile').read_text()
    if name=='strong_screened': assert completed and len(x)==12000
    norm_samples=[]
    for line in (d/'outputfile').read_text().splitlines():
        cols=line.split()
        if len(cols)!=7: continue
        try:
            step=int(cols[0]); row=list(map(float,cols[1:]))
        except ValueError: continue
        if step>0 and abs(row[0]-step*.08*fs)<1e-6:
            norm_samples.append((row[0],row[4]))
    bad_norm=[t for t,n in norm_samples if abs(n-32)>1e-4]
    norm_loss_time=bad_norm[0] if bad_norm else None
    np.testing.assert_allclose(x[:,0],j[:,0],atol=1e-8,rtol=0)
    np.testing.assert_allclose(j[j[:,0]>60,1:13],0,atol=1e-12,rtol=0)
    t=x[:,0]*fs;post=t>60*fs
    windows=[]
    for left,right in [(60*fs,5),(5,10),(10,15),(15,20),(20,t[-1])]:
        right=min(right,float(t[-1]))
        if right<=left: continue
        mask=(t>=left)&(t<=right)
        if mask.sum()<2: continue
        tt=x[mask,0];a=x[mask,3];e=x[mask,6];current=j[mask,15]
        slope,offset=np.polyfit(tt,a,1)
        residual=a-(slope*tt+offset)
        windows.append(dict(from_fs=left,to_fs=right,mean_J=float(current.mean()),
                            rms_J=float(np.sqrt(np.mean(current**2))),
                            max_abs_J=float(np.max(abs(current))),mean_Exc=float(e.mean()),
                            rms_Exc=float(np.sqrt(np.mean(e**2))),
                            A_slope_au=float(slope),A_linear_residual_rms=float(np.sqrt(np.mean(residual**2)))))
    # First post-pulse sign crossing of XC A/c, linearly interpolated.
    indices=np.where((x[1:,3]*x[:-1,3]<0)&(t[:-1]>60*fs))[0]
    crossing=None
    if len(indices):
        i=indices[0];crossing=float(t[i]-x[i,3]*(t[i+1]-t[i])/(x[i+1,3]-x[i,3]))
    offset=x[post,6]-x[post,7]*x[post,11]
    results[name]=dict(post_Exc_minus_alphaP_mean=float(offset.mean()),
                      post_Exc_minus_alphaP_range=[float(offset.min()),float(offset.max())],
                      final_alphaP=float(x[-1,7]*x[-1,11]),completed=completed,rows=len(x),first_norm_error_over_1e4_electrons_fs=norm_loss_time,
                      max_printed_norm_error_electrons=max(abs(n-32) for _,n in norm_samples),final_time_fs=float(t[-1]),A_end=float(x[-1,3]),Exc_end=float(x[-1,6]),
                      J_end=float(j[-1,15]),alpha_end=float(x[-1,7]),first_postpulse_A_zero_fs=crossing,
                      max_post_abs_A=float(np.max(abs(x[post,3]))),
                      max_post_abs_J=float(np.max(abs(j[post,15]))),
                      post_alpha_range=[float(x[post,7].min()),float(x[post,7].max())],
                      electronic_energy_change_au=float(energy[-1,2]),windows=windows)
    data[name]=(t,x,j)
(out/'long_metrics.json').write_text(json.dumps(results,indent=2)+'\n')
fig,axes=plt.subplots(3,2,figsize=(11,8),sharex=True,layout='constrained')
for col,(name,label,color) in enumerate([('strong_fixed','Fixed alpha=0.2','0.35'),
                                        ('strong_screened','Instantaneous trial','tab:blue')]):
    t,x,j=data[name]
    axes[0,col].plot(t,x[:,3],color=color)
    axes[1,col].plot(t,x[:,6],color=color)
    axes[2,col].plot(t,j[:,15],color=color,lw=.8)
    axes[0,col].set_title(label+(' (stopped)' if not results[name]['completed'] else ''))
    for row in range(3):
        axes[row,col].axvline(60*fs,color='tab:red',ls=':',lw=1)
        axes[row,col].axhline(0,color='0.7',lw=.7)
        axes[row,col].grid(alpha=.2)
        axes[row,col].set_xlim(0,12000*.08*fs)
        loss=results[name]['first_norm_error_over_1e4_electrons_fs']
        if loss is not None: axes[row,col].axvspan(loss,t[-1],color='tab:orange',alpha=.2)
        axes[row,col].set_ylabel(['XC A/c z (a.u.)','XC electric field z (a.u.)','Number current z (a.u.)'][row])
    axes[2,col].set_xlabel('Time (fs)')
fig.suptitle('Si 4x4x4 k: post-pulse evolution to 23.22 fs\nSame 5.44 eV, 1e13 W/cm2 pulse; each panel uses its own vertical scale')
fig.savefig(out/'long_comparison.png',dpi=160)
fig.savefig(out/'long_comparison.pdf')
print(json.dumps(results,indent=2))
