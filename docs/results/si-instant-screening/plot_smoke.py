from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[3]
work=root/'calculations/si_tdcdft_k4/instant_smoke'
fig,ax=plt.subplots(3,1,figsize=(8,8),sharex=True,layout='constrained')
for name,label,color in [('strong_fixed','Fixed alpha=0.2','0.35'),('strong_screened','Instantaneous trial','tab:blue')]:
    j=np.loadtxt(work/name/'Si_rt.data');x=np.loadtxt(work/name/'Si_rt_xc.data')
    ax[0].plot(j[:,0]*.02418884326505,j[:,15],label=label,color=color)
    ax[1].plot(x[:,0]*.02418884326505,x[:,3],color=color)
    ax[2].plot(x[:,0]*.02418884326505,x[:,7],label=label,color=color)
for floor,color in [(1e-5,'tab:orange'),(4e-5,'tab:green')]:
    p=work/('strong_floor_'+str(floor))/'Si_rt_xc.data'
    if p.exists():
        x=np.loadtxt(p); ax[2].plot(x[:,0]*.02418884326505,x[:,7],ls='--',color=color,label=f'floor={floor:g}')
for a in ax:
    a.axvline(60*.02418884326505,ls=':',color='0.6');a.grid(alpha=.2)
ax[0].legend();ax[2].legend(fontsize=8)
ax[0].set_ylabel('Number current z (a.u.)');ax[1].set_ylabel('XC A/c z (a.u.)');ax[2].set_ylabel('Effective alpha')
ax[2].set_xlabel('Time (fs)');ax[2].set_ylim(-.01,.22)
fig.suptitle('Si 4x4x4 k points: short strong-pulse numerical trial\n5.44 eV, 1e13 W/cm2; no accuracy claim')
fig.savefig(Path(__file__).with_name('comparison.png'),dpi=160)
fig.savefig(Path(__file__).with_name('comparison.pdf'))
