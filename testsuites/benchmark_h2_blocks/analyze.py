"""Summarize completed block measurements without hiding numerical errors."""
import json,math,statistics,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
src=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
d=json.loads(src.read_text())
partial='--partial' in sys.argv[3:]
if not d['complete'] and not partial:raise RuntimeError('incomplete measurement matrix')
summary=[]
for c in d['cases']:
    rows={m:[r for r in d['runs'] if r['shape']==c['shape'] and r['ranks']==c['ranks'] and r['mode']==m] for m in ('full','adaptive')}
    if any(len(rr)!=c.get('repeats',d['repeats']) for rr in rows.values()):
        if partial:continue
        raise RuntimeError('missing repeats')
    item=dict(c)
    for m,rr in rows.items():
        item[m]={k:dict(min=min(r[k] for r in rr),median=statistics.median(r[k] for r in rr),max=max(r[k] for r in rr)) for k in ('rt_max_seconds','peak_rank_bytes','post_impulse_energy_width_ha')}
    rr=rows['adaptive'];pairs=[p for r in rr for p in r['pair_counters']];supports=[s for r in rr for s in r['support_counters']]
    local=sum(p[0] for p in pairs);glob=sum(p[1] for p in pairs);points=sum(p[2] for p in pairs)
    item.update(speedup=item['full']['rt_max_seconds']['min']/item['adaptive']['rt_max_seconds']['min'],
      local_pairs=local,global_pairs=glob,local_fft_grid_ratio=points/local/(32**3*math.prod(c['shape'])) if local else None,
      max_radius=max((s[1] for s in supports),default=None),max_loss=max((s[2] for s in supports),default=None),
      active_updates=[len(r['support_counters']) for r in rr],retained_gauge_updates=[r['retained_gauge_updates'] for r in rr],
      max_current_difference=max(abs(a-b) for f,r in zip(rows['full'],rows['adaptive']) for fr,ar in zip(f['observables'],r['observables']) for a,b in zip(fr[13:16],ar[13:16])),
      max_energy_difference_ha=max(abs(fr[1]-ar[1]) for f,r in zip(rows['full'],rows['adaptive']) for fr,ar in zip(f['energies'],r['energies'])))
    summary.append(item)
# Compare each strong layout with MPI1, same mode and first repetition.
strong_differences={}
for mode in ('full','adaptive'):
    rows=[r for r in d['runs'] if 'strong' in r['suites'] and r['mode']==mode]
    if not rows:continue
    ref=next((r for r in rows if r['ranks']==1),None)
    if ref is None:continue
    strong_differences[mode]=dict(current=max(abs(x-y) for r in rows for a,b in zip(ref['observables'],r['observables']) for x,y in zip(a[13:16],b[13:16])),energy_ha=max(abs(a[1]-b[1]) for r in rows for a,b in zip(ref['energies'],r['energies'])))
(out/'summary.json').write_text(json.dumps(dict(complete=d['complete'],summary=summary,strong_differences=strong_differences,missing_cases=[c for c in d['cases'] if not any(r['shape']==c['shape'] and r['ranks']==c['ranks'] for r in summary)]),indent=2)+'\n')
lines=['Interim: remaining cases are pending.',''] if not d['complete'] else []
for suite in ('weak','strong'):
    lines+=['## '+suite,'','|Fragment array|H₂|MPI|Runs per mode|Full seconds: single or min [median,max]|99.9% seconds: single or min [median,max]|Speedup|Peak RSS full / .999 MiB|Radius bohr|Local FFT volume ratio|Energy width .999 Ha|','|---|---:|---:|---:|---|---|---:|---:|---:|---:|---:|']
    table_rows=[r for r in summary if suite in r['suites']]
    if suite=='strong':table_rows.sort(key=lambda r:r['ranks'])
    for r in table_rows:
        def t(m):
            x=r[m]['rt_max_seconds']
            if r.get('repeats',d['repeats'])==1:return f"{x['min']:.4f}"
            return f"{x['min']:.4f} [{x['median']:.4f},{x['max']:.4f}]"
        ratio='—' if r['local_fft_grid_ratio'] is None else f"{r['local_fft_grid_ratio']:.4f}"
        radius='—' if r['max_radius'] is None else f"{r['max_radius']:.3f}"
        lines.append('|'+ '×'.join(map(str,r['shape']))+f"|{8*math.prod(r['shape'])}|{r['ranks']}|{r.get('repeats',d['repeats'])}|{t('full')}|{t('adaptive')}|{r['speedup']:.3f}|{r['full']['peak_rank_bytes']['median']/2**20:.1f} / {r['adaptive']['peak_rank_bytes']['median']/2**20:.1f}|{radius}|{ratio}|{r['adaptive']['post_impulse_energy_width_ha']['max']:.3e}|")
    lines+=['']
(out/'tables.md').write_text('\n'.join(lines))
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for col,suite in enumerate(('weak','strong')):
    rr=[r for r in summary if suite in r['suites']]
    if suite=='strong':rr.sort(key=lambda r:r['ranks'])
    if not rr:
        for ax in axes[:,col]:ax.set_axis_off();ax.text(.5,.5,'Pending',ha='center',va='center')
        continue
    x=list(range(len(rr)));labels=[('x'.join(map(str,r['shape']))+'\n'+str(r['ranks'])+' MPI') if suite=='weak' else str(r['ranks']) for r in rr]
    for mode,color in [('full','#4968b1'),('adaptive','#d27724')]:
        ys=[r[mode]['rt_max_seconds']['min'] for r in rr]
        axes[0,col].plot(x,ys,'o-',label=mode,color=color)
        axes[0,col].vlines(x,ys,[r[mode]['rt_max_seconds']['max'] for r in rr],color=color,alpha=.5)
        axes[1,col].plot(x,[r[mode]['peak_rank_bytes']['median']/2**20 for r in rr],'o-',label=mode,color=color)
    for ax in axes[:,col]:
        ax.set_xticks(x,labels,fontsize=8);ax.grid(alpha=.2);ax.legend()
    axes[0,col].set_title(suite+' scaling: 8 H2 per fragment core');axes[0,col].set_ylabel('16-step RT time / s (single or minimum, range)');axes[0,col].set_yscale('log')
    axes[1,col].set_ylabel('Max-rank lifetime peak RSS / MiB');axes[1,col].set_xlabel('Fragment array / MPI ranks' if suite=='weak' else 'MPI ranks')
fig.savefig(out/'scaling.png',dpi=160);fig.savefig(out/'scaling.svg')
print('\n'.join(lines))
