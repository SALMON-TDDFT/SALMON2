import os
os.environ['MPLCONFIGDIR']='/private/tmp/salmon-mpl-cache'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
out=Path.cwd()/'docs/results/si-k4-pump-probe';a=np.loadtxt(out/'g40_880au.csv',delimiter=',');h=np.loadtxt(out/'g40_halfprobe_880au.csv',delimiter=',');f,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
ax[0].plot(a[:,0],a[:,2],label='Unpumped');ax[0].plot(h[:,0],h[:,2],label='Pumped (half probe)');ax[0].set(xlim=(1.5,6),xlabel='Energy (eV)',ylabel='Im epsilon');ax[0].legend()
for w in [600,720,880]:
 d=np.loadtxt(out/f'g40_{w}au.csv',delimiter=',');ax[1].plot(d[:,0],d[:,4],label=f'{w*.024188843:.2f} fs window')
ax[1].set(xlim=(1.5,6),xlabel='Energy (eV)',ylabel='Pumped Im epsilon');ax[1].legend()
for x in ax:
 values=np.concatenate([line.get_ydata()[(line.get_xdata()>=1.5)&(line.get_xdata()<=6)] for line in x.lines]);lo,hi=values.min(),values.max();x.set_ylim(lo-.08*(hi-lo),hi+.08*(hi-lo))
for x in ax:x.axhline(0,color='0.5',lw=.6);x.grid(alpha=.2)
f.suptitle('Si 4³ fixed-alpha TDCDFT: alpha=0.2, gamma=0.004; probe at 1.94 fs')
f.savefig(out/'tdcdft_g40_preliminary.png',dpi=150)
