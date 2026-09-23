"""Offline UEG-distance alpha on existing, unchanged TDCDFT trajectories."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
from pathlib import Path
import sys,json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;root=out.parents[2]
sys.path.insert(0,str(out.parent/'si-time-wannier'))
from wannier import read_u,geometry
from ueg import (global_indices,fourier_coefficients,fermi_reference,flow_reference,
                 mean_momentum,hs_distance,rank_lower_bound,alpha_from_distance)
r,c,k,b,L,h=geometry(root);nk=len(k);rank=16;grid=12;spacing=2*np.pi/(4*L)
indices,q=global_indices(k,grid,L,4)
assert len(np.unique(indices))==q[...,0].size
base=fermi_reference(q,nk*rank);base_k=base.ravel()[indices]
assert abs(base.sum()-nk*rank)<1e-9
np.testing.assert_allclose(mean_momentum(base,q),0,atol=1e-12)
u0,_=read_u(root/'calculations/si_tdcdft_k4/gs/data_for_restart')
c0=fourier_coefficients(u0,grid,h**3);d0=hs_distance(c0,base_k)
power0=np.sum(abs(c0)**2,axis=2)
pmean0=np.einsum('kg,kga->a',power0,q.reshape(-1,3)[indices])/power0.sum()
d0_modes={'flow':d0,'rest':d0,'spectral':hs_distance(c0,flow_reference(base,pmean0,spacing).ravel()[indices])}
summary={'initial_distance_one_spin_per_cell':d0,'initial_distances_by_reference':d0_modes,'reference':dict(number_one_spin_supercell=float(base.sum()),
    purity_ratio=float(np.sum(base**2)/base.sum()),fractional_states=int(np.count_nonzero((base>0)&(base<1))),
    per_k_number_min=float(base_k.sum(axis=1).min()),per_k_number_max=float(base_k.sum(axis=1).max()),
    fixed_rank_distance_lower_bound=rank_lower_bound(base_k,rank),
    fixed_rank_alpha_lower_bound=float(alpha_from_distance(rank_lower_bound(base_k,rank),d0))),
    'model':'alpha=.2*HS_squared_distance/initial; unclipped; offline only',
    'sampling':dict(stride=10,dt_au=.08,dt_fs=.08*.02418884326505)}
traces={};work=Path('/private/tmp/salmon-si-time-wannier-dense');started=time.monotonic()
for case in ('none','weak','strong'):
    source=out.parent/'si-time-wannier'
    rt=np.loadtxt(source/f'{case}_rt.data');xc=np.loadtxt(source/f'{case}_rt_xc.data')
    mlwf=json.loads((source/f'{case}_mlwf.json').read_text());records=[]
    for frame,step in enumerate(range(0,1601,10)):
        u=u0 if step==0 else read_u(work/case/f'checkpoint_rt_{step:06d}')[0]
        coeff=fourier_coefficients(u,grid,h**3)
        power=np.sum(abs(coeff)**2,axis=2);number=float(power.sum()/nk)
        if abs(number-rank)>1e-6:raise RuntimeError('Orbital normalization error')
        pmean=np.einsum('kg,kga->a',power,q.reshape(-1,3)[indices])/power.sum()
        a=np.zeros(3) if step==0 else rt[step-1,7:10]+xc[step-1,1:4]
        current=np.zeros(3) if step==0 else rt[step-1,13:16]
        target=current/(32/L**3)-a
        row=dict(step=step,time_fs=step*.08*.02418884326505,canonical_fft_mean=pmean.tolist(),
                 flow_target=target.tolist(),current_vs_spectral_velocity_error=float(np.linalg.norm(target-pmean)),
                 mlwf_spread_A2=mlwf[frame]['mean_spread_A2'],invariant_spread_A2=mlwf[frame]['mean_invariant_A2'],
                 density_relative_rms=float(np.std(np.mean(np.sum(abs(u)**2,axis=2),axis=0))/(rank/L**3)))
        for mode,p in [('flow',target),('spectral',pmean),('rest',np.zeros(3))]:
            f=flow_reference(base,p,spacing);fk=f.ravel()[indices]
            error=float(np.max(abs(mean_momentum(f,q)-p)))
            assert error<1e-10 and abs(f.sum()-nk*rank)<1e-8
            distance=hs_distance(coeff,fk)
            bound=rank_lower_bound(fk,rank)
            assert distance>=bound-1e-6
            row[mode]=dict(distance=distance,D=distance/d0_modes[mode],alpha=float(alpha_from_distance(distance,d0_modes[mode])),
                occupied_purity=float(np.sum(abs(coeff.conj().transpose(0,2,1)@coeff)**2)/nk),
                reference_purity_per_cell=float(np.sum(f*f)/nk),overlap_per_cell=float(np.sum(fk*power)/nk),
                reference_purity=float(np.sum(f*f)/(nk*rank)),flow_match_error=error,
                rank_bound_alpha=float(alpha_from_distance(bound,d0_modes[mode])))
        records.append(row)
        if step%400==0:print(case,step,row['flow']['alpha'],flush=True)
    traces[case]=records
    (out/f'{case}.json').write_text(json.dumps(records,indent=2)+'\n')
    summary[case]=dict(final={mode:records[-1][mode] for mode in ('flow','spectral','rest')},
        alpha_min=min(z['flow']['alpha'] for z in records),alpha_max=max(z['flow']['alpha'] for z in records),
        max_flow_vs_spectral_alpha_difference=max(abs(z['flow']['alpha']-z['spectral']['alpha']) for z in records),
        max_flow_vs_rest_alpha_difference=max(abs(z['flow']['alpha']-z['rest']['alpha']) for z in records),
        max_flow_match_error=max(z['flow']['flow_match_error'] for z in records),
        max_current_vs_spectral_velocity_error=max(z['current_vs_spectral_velocity_error'] for z in records),
        density_relative_rms_initial=records[0]['density_relative_rms'],density_relative_rms_final=records[-1]['density_relative_rms'])
summary['elapsed_seconds']=time.monotonic()-started
(out/'metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,ax=plt.subplots(3,1,figsize=(9,9),sharex=True,layout='constrained');colors={'none':'0.4','weak':'tab:blue','strong':'tab:red'}
for case,z in traces.items():
    t=[x['time_fs'] for x in z]
    ax[0].plot(t,[x['flow']['alpha'] for x in z],label=case,color=colors[case])
    ax[1].plot(t,[x['mlwf_spread_A2'] for x in z],color=colors[case])
    ax[2].plot(t,[x['flow']['reference_purity'] for x in z],color=colors[case])
ax[0].axhline(.2,color='0.5',ls='--');ax[0].set_ylabel('Proposed alpha = 0.2 D');ax[0].legend()
ax[1].set_ylabel('MLWF spread (A^2 / orbital)');ax[2].set_ylabel('UEG reference purity / N');ax[2].set_xlabel('Time (fs)')
for a in ax:a.axvline(60*.02418884326505,color='0.5',ls=':');a.grid(alpha=.2)
fig.suptitle('Si 4x4x4 k: UEG density-matrix distance\nUnclipped alpha; offline evaluation on previous trajectories')
fig.savefig(out/'distance.png',dpi=150);fig.savefig(out/'distance.pdf')
fig,ax=plt.subplots(figsize=(9,4),layout='constrained')
z=traces['strong'];t=[x['time_fs'] for x in z]
for mode,label in [('flow','Match SALMON current'),('spectral','Match spectral momentum'),('rest','Resting reference')]:
    ax.plot(t,[x[mode]['alpha'] for x in z],label=label)
ax.plot(t,[x['flow']['rank_bound_alpha'] for x in z],ls=':',label='Fixed-rank lower bound')
ax.set(xlabel='Time (fs)',ylabel='Proposed alpha',title='Strong pulse: reference convention sensitivity');ax.legend();ax.grid(alpha=.2)
fig.savefig(out/'reference_controls.png',dpi=150);fig.savefig(out/'reference_controls.pdf')
print(json.dumps(summary,indent=2),flush=True)
