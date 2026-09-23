"""TDELF diagnostics from unchanged Si RT wavefunction checkpoints."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
from pathlib import Path
import sys,json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;root=out.parents[2]
sys.path.insert(0,str(out.parent/'si-time-wannier'))
from wannier import geometry,read_u
from tdelf import gradients,fields
r,c,k,b,L,h=geometry(root);u0,_=read_u(root/'calculations/si_tdcdft_k4/gs/data_for_restart')
bond=np.min(np.linalg.norm((r[:,None,:]-b[None,:,:]+L/2)%L-L/2,axis=2),axis=1)<1.
assert bond.any()
bins=np.linspace(0,1,51);saved={};traces={};summary={};started=time.monotonic()
work=Path('/private/tmp/salmon-si-time-wannier-dense')
for case in ('none','weak','strong'):
    rows=[];hist=[];maps=[];densities=[]
    for step in range(0,1601,10):
        u=u0 if step==0 else read_u(work/case/f'checkpoint_rt_{step:06d}')[0]
        z=fields(u,gradients(u,k,L,12,zero_nyquist=True));n=z['n'];elf=z['elf'];w=n/n.sum()
        rows.append(dict(step=step,time_fs=step*.08*.02418884326505,
             mean_elf=float(w@elf),density_weighted_deviation_from_half=float(w@(elf-.5)**2),
             fraction_near_half=float(w[abs(elf-.5)<.05].sum()),fraction_above_075=float(w[elf>.75].sum()),
             bond_mean_elf=float(n[bond]@elf[bond]/n[bond].sum()),outside_bond_mean_elf=float(n[~bond]@elf[~bond]/n[~bond].sum()),
             max_current_correction=float(np.max(abs(elf-z['elf_without_current']))),
             mean_current_correction=float(w@(elf-z['elf_without_current'])),density_min=float(n.min()),
             min_curvature=float(z['curvature'].min()),max_elf=float(elf.max()),min_elf=float(elf.min())))
        hist.append(np.histogram(elf,bins=bins,weights=w)[0]);maps.append(elf);densities.append(n)
        if step in (0,800,1600):saved[f'{case}_{step}']=elf;saved[f'{case}_{step}_n']=n
        if step%400==0:print(case,step,rows[-1]['mean_elf'],rows[-1]['bond_mean_elf'],flush=True)
    traces[case]=rows
    (out/f'{case}.json').write_text(json.dumps(rows,indent=2)+'\n')
    np.savez_compressed(out/f'{case}_fields.npz',steps=np.arange(0,1601,10),elf=np.array(maps),one_spin_density=np.array(densities),histogram=np.array(hist),histogram_bins=bins)
    summary[case]=dict(initial=rows[0],final=rows[-1],max_current_correction=max(x['max_current_correction'] for x in rows),
        max_weighted_current_correction=max(x['mean_current_correction'] for x in rows),
        max_field_change_from_initial=float(np.max(abs(np.array(maps)-maps[0]))),snapshots=len(rows))
summary['elapsed_seconds']=time.monotonic()-started
summary['definition']='One-spin Burnus-Marques-Gross TDELF; orbital current correction; spectral Bloch gradients with zero odd-derivative Nyquist mode'
summary['bond_region']='Union of fixed spheres radius1 bohr around initial16 bond centers; electron weighted; no basin topology claim'
summary['interpretation']='Diagnostic only; no ELF-to-alpha rule and no proof of metallic transport'
(out/'metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,ax=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained');colors={'none':'0.4','weak':'tab:blue','strong':'tab:red'}
for case,rows in traces.items():
    t=[x['time_fs'] for x in rows]
    for a,key in zip(ax,['mean_elf','bond_mean_elf','fraction_near_half']):a.plot(t,[x[key] for x in rows],label=case,color=colors[case])
ax[0].set_ylabel('Electron-weighted mean ELF');ax[0].legend();ax[1].set_ylabel('Bond-region mean ELF');ax[2].set_ylabel('Electron fraction\n0.45 < ELF < 0.55');ax[2].set_xlabel('Time (fs)')
for a in ax:a.axvline(60*.02418884326505,color='0.5',ls=':');a.grid(alpha=.2)
fig.suptitle('Si 4x4x4 k: current-corrected time-dependent ELF\nExisting trajectories; localization diagnostic, not a metallicity test')
fig.savefig(out/'dynamics.png',dpi=150);fig.savefig(out/'dynamics.pdf')
fig,ax=plt.subplots(1,3,figsize=(12,4),layout='constrained')
for a,key,title in zip(ax,['none_0','strong_800','strong_1600'],['Initial','Strong: 1.548 fs','Strong: 3.096 fs']):
    cube=saved[key].reshape((12,)*3,order='F');plane=np.array([cube[i,i,:] for i in range(12)])
    im=a.imshow(plane.T,origin='lower',extent=[0,np.sqrt(2)*L*.529177,L*.529177*0,L*.529177],vmin=0,vmax=1,cmap='viridis',aspect='auto')
    a.set(title=title,xlabel='[110] distance (A)',ylabel='z (A)')
fig.colorbar(im,ax=ax,label='ELF');fig.suptitle('ELF on x=y plane; shared color scale')
fig.savefig(out/'maps.png',dpi=150);fig.savefig(out/'maps.pdf')
fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
for case,step,label in [('none',0,'Initial'),('strong',800,'Strong 1.548 fs'),('strong',1600,'Strong 3.096 fs')]:
    data=np.load(out/f'{case}_fields.npz');hh=data['histogram'][step//10]
    ax.stairs(hh,bins,label=label)
ax.axvline(.5,color='0.5',ls=':');ax.set(xlabel='ELF',ylabel='Electron fraction per bin',title='ELF distribution (bin width 0.02)');ax.legend()
fig.savefig(out/'histogram.png',dpi=150);fig.savefig(out/'histogram.pdf')
print(json.dumps(summary,indent=2),flush=True)
