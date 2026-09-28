"""Render archived H2 SCF scaling best observed times and repeat ranges."""
import argparse,json,math,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/salmon-h2-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
d=json.loads(a.results.read_text());assert d['complete']
weak={tuple(s['shape']):s for s in d['summary'] if 'weak' in s['suites']}
strong=sorted([s for s in d['summary'] if 'strong' in s['suites']],key=lambda s:s['ranks'])
fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
series=[('n x 1 x 1',[(1,1,1),(2,1,1),(4,1,1),(8,1,1),(16,1,1)]),
        ('n x n x 1',[(1,1,1),(2,2,1),(4,4,1)]),('n x n x n',[(1,1,1),(2,2,2)])]
for ax,key,scale,title,ylabel in [(axes[0,0],'scf_seconds_per_iteration',1000,'Weak scaling: amortized SCF cost','ms / SCF iteration'),
        (axes[0,1],'peak_rank_bytes',1/2**20,'Weak scaling: process memory','Max rank peak RSS (MiB)')]:
 for label,shapes in series:
  rows=[weak[s][key] for s in shapes];stat='min' if key=='scf_seconds_per_iteration' else 'median';y=[r[stat]*scale for r in rows]
  ax.errorbar([math.prod(s) for s in shapes],y,yerr=[[(r[stat]-r['min'])*scale for r in rows],[(r['max']-r[stat])*scale for r in rows]],marker='o',capsize=3,label=label)
 ax.set(xscale='log',yscale='log',title=title,xlabel='Molecules = MPI ranks',ylabel=ylabel)
 ax.set_xticks([1,2,4,8,16],labels=['1','2','4','8','16']);ax.legend(fontsize=8);ax.grid(alpha=.25,which='both')
ax=axes[1,0];ranks=[s['ranks'] for s in strong]
for key,label in [('scf_max_seconds','Time to SCF convergence'),('scf_seconds_per_iteration','Amortized per iteration')]:
 base=strong[0][key]['min'];ax.plot(ranks,[base/s[key]['min'] for s in strong],'-o',label=label)
ax.plot(ranks,ranks,'--',color='gray',label='Ideal');ax.set(xscale='log',yscale='log',title='Strong scaling: 4 x 4 x 1',xlabel='MPI ranks',ylabel='Speedup vs MPI 1');ax.legend(fontsize=8)
ax.set_xticks(ranks,labels=list(map(str,ranks)));ax.grid(alpha=.25,which='both')
ax=axes[1,1];rows=[s['peak_rank_bytes'] for s in strong]
ax.errorbar(ranks,[r['median']/2**20 for r in rows],yerr=[[(r['median']-r['min'])/2**20 for r in rows],[(r['max']-r['median'])/2**20 for r in rows]],marker='o',capsize=3)
ax.set(xscale='log',title='Strong scaling: process memory',xlabel='MPI ranks',ylabel='Max rank peak RSS (MiB)');ax.set_xticks(ranks,labels=list(map(str,ranks)));ax.grid(alpha=.25)
fig.suptitle('H2 / PBEh(40), full support, single node\n3 runs: best times, median RSS, bars min-max; shared workstation / MPI 1 legacy',fontsize=12)
a.output.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(a.output.with_suffix('.png'),dpi=180);fig.savefig(a.output.with_suffix('.svg'))

svg=a.output.with_suffix('.svg')
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
