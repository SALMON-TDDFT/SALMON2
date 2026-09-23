"""Reproduce the native Taylor/PT-CN comparison figure from completed runs."""
from pathlib import Path
import json,os
os.environ['MPLCONFIGDIR']='/private/tmp/hse-taylor-mpl'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
root=Path(__file__).resolve().parents[2]
p=root/'calculations/si_hse_native';out=root/'docs/results/si-hse-taylor'
fig,axes=plt.subplots(1,3,figsize=(13.8,4.0),layout='constrained')
colors={.32:'#be4935',.16:'#ca9227',.08:'#287ca5',.04:'#326b49'}
ref=np.loadtxt(p/'compare_hse_taylor4_dt0.04/Si_rt.data')
axes[0].plot(ref[:,0]*.024188843265857,ref[:,15]*1e6,color='black',lw=2,label='Taylor4, dt=0.04')
for dt in [.32,.16,.08]:
 a=np.loadtxt(p/f'compare_hse_ptcn_dt{dt:.2f}/Si_rt.data');axes[0].plot(a[:,0]*.024188843265857,a[:,15]*1e6,color=colors[dt],ls='--',lw=1.2,label=f'PT-CN, dt={dt}')
axes[0].set(xlabel='Time (fs)',ylabel='Current Jz (10^-6 a.u.)',title='Si impulse response');axes[0].legend(fontsize=8)
s=json.loads((out/'convergence.json').read_text())['results']
for mode,color,label in [('hse_ptcn','#287ca5','PT-CN'),('hse_taylor4','#222222','Taylor4 + PC')]:
 points=[]
 for name,v in s.items():
  if name.startswith('compare_'+mode+'_dt'):
   dt=float(name.rsplit('dt',1)[1]);err=v['current_relative_rms_error']*100
   if err>0:points.append((dt,err))
 points.sort();x,y=np.array(points).T;axes[1].loglog(x,y,'o-',color=color,label=label)
x=np.array([.04,.32]);axes[1].loglog(x,.065*(x/.16)**2,':',color='#999999',label='dt squared guide')
axes[1].set(xlabel='Time step (a.u.)',ylabel='Current RMS difference (%)',title='Against Taylor4, dt=0.04');axes[1].xaxis.set_minor_formatter(NullFormatter());axes[1].set_xticks([.04,.08,.16,.32], labels=['0.04','0.08','0.16','0.32']);axes[1].legend(fontsize=8);axes[1].grid(True,which='both',alpha=.2)
for name,label,col,ls in [('compare_long_hse_taylor4_dt0.08','Taylor4, dt=0.08','black','-'),('rt_stability_n8','PT-CN, dt=0.32',colors[.32],'--'),('compare_long_hse_ptcn_dt0.16','PT-CN, dt=0.16',colors[.16],':')]:
 a=np.loadtxt(p/name/'Si_rt.data');axes[2].plot(a[:,0]*.024188843265857,a[:,15]*1e6,label=label,color=col,ls=ls,lw=1.5)
axes[2].set(xlabel='Time (fs)',ylabel='Current Jz (10^-6 a.u.)',title='Longer pilot: 0.774 fs');axes[2].legend(fontsize=8)
fig.savefig(out/'comparison.png',dpi=190);fig.savefig(out/'comparison.pdf')
