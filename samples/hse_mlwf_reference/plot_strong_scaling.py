from pathlib import Path
import json,numpy as np,os
os.environ['MPLCONFIGDIR']='/private/tmp/hse-taylor-mpl'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2];out=root/'docs/results/si-hse-strong-scaling'
runs=json.loads((out/'runs.json').read_text());assert len(runs)==10
summary=[]
for n in (1,2,4,8,16):
 rows=[r for r in runs if r['ranks']==n];times=[r['seconds_per_step'] for r in rows]
 summary.append(dict(ranks=n,seconds_per_step=float(np.mean(times)),min_seconds_per_step=min(times),max_seconds_per_step=max(times),max_rank_rss_MiB=max(r['sampled_peak_rank_rss_KiB'] for r in rows)/1024,max_sum_rss_GiB=max(r['sampled_peak_sum_rss_KiB'] for r in rows)/1048576,orbital_relative_difference=max(r['orbital_relative_difference'] for r in rows)))
for r in summary:
 r['speedup']=summary[0]['seconds_per_step']/r['seconds_per_step'];r['parallel_efficiency']=r['speedup']/r['ranks']
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
x=np.array([r['ranks'] for r in summary]);t=np.array([r['seconds_per_step'] for r in summary]);e=np.array([[r['seconds_per_step']-r['min_seconds_per_step'] for r in summary],[r['max_seconds_per_step']-r['seconds_per_step'] for r in summary]])
fig,ax=plt.subplots(1,3,figsize=(12.8,3.7),layout='constrained')
ax[0].errorbar(x,t,yerr=e,fmt='o-',capsize=4,label='Measured (2 runs)');ax[0].plot(x,t[0]/x,'--',color='gray',label='Ideal 1/P');ax[0].set(ylabel='Seconds / RT step',title='Fixed Si: Taylor + ACE');ax[0].legend(fontsize=8)
ax[1].plot(x,[r['speedup'] for r in summary],'o-',label='Measured');ax[1].plot(x,x,'--',color='gray',label='Ideal');ax[1].set(ylabel='Speedup vs MPI1',title='Strong scaling');ax[1].legend(fontsize=8)
ax[2].plot(x,[r['max_rank_rss_MiB'] for r in summary],'o-',color='#287ca5',label='Largest rank');ax[2].set(ylabel='Largest rank RSS (MiB)',title='Sampled RSS maxima');a2=ax[2].twinx();a2.plot(x,[r['max_sum_rss_GiB'] for r in summary],'s--',color='#be4935',label='All ranks');a2.set_ylabel('Sum RSS (GiB)',color='#be4935');ax[2].legend(loc='upper center',fontsize=8);a2.legend(loc='center right',fontsize=8)
for a in ax:a.set_xscale('log',base=2);a.set_xticks(x,labels=[str(v) for v in x]);a.set_xlabel('MPI ranks');a.grid(alpha=.2)
fig.savefig(out/'scaling.png',dpi=180);fig.savefig(out/'scaling.pdf')
print(json.dumps(summary,indent=2))
