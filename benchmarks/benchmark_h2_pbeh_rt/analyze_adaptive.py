"""Analyze paired RT measurements without treating completion as accuracy approval."""
import argparse, json, math, os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/salmon-h2-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('results',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();d=json.loads(a.results.read_text());assert d['complete']
    a.output.mkdir(parents=True,exist_ok=True)
    def runs(s,mode):
        return [r for r in d['runs'] if r['shape']==s['shape'] and r['ranks']==s['ranks'] and r['mode']==mode]
    def error(rows,reference,field,columns):
        return max(abs(row[c]-base[c]) for r in rows for row,base in zip(r[field],reference[field]) for c in columns)
    diagnostics=[]
    for s in d['summary']:
        if s['mode']!='adaptive':continue
        rows=runs(s,'adaptive');refs=runs(s,'full');ref=refs[0]
        counters=[v for r in rows for v in r['pair_counters']]
        lp=sum(v[0] for v in counters);gp=sum(v[1] for v in counters);pts=sum(v[2] for v in counters)
        ng=4096*math.prod(s['shape'])
        diagnostics.append(dict(shape=s['shape'],ranks=s['ranks'],
            active_refreshes=[len(r['support_counters']) for r in rows],
            full_reference_refreshes=[len(r['support_counters']) for r in refs],
            local_pairs=lp,global_pairs_during_active_support=gp,
            local_fft_volume_ratio=pts/(lp*ng) if lp else None,
            maximum_norm_loss=max((v[2] for r in rows for v in r['support_counters']),default=None),
            current_error_au=error(rows,ref,'observables',range(13,16)),
            reference_peak_current_au=max(abs(v[c]) for v in ref['observables'] for c in range(13,16)),
            energy_error_ha=error(rows,ref,'energies',[1]),
            repeat_current_error_au=error(rows,rows[0],'observables',range(13,16)),
            maximum_energy_width_ha=max(r['post_impulse_energy_width_ha'] for r in rows)))
    layout={}
    for mode in ('full','adaptive'):
        rows=[r for r in d['runs'] if 'strong' in r['suites'] and r['mode']==mode]
        ref=next(r for r in rows if r['ranks']==1)
        layout[mode]=dict(current_error_au=error(rows,ref,'observables',range(13,16)),
            energy_error_ha=error(rows,ref,'energies',[1]))
    (a.output/'diagnostics.json').write_text(json.dumps(dict(cases=diagnostics,strong_layout_reference_errors=layout),indent=2)+'\n')
    lines=[]
    for suite in ('weak','strong'):
        lines += ['## '+suite,'','|形状|MPI|全領域 秒 min [median,max]|99.9% 秒 min [median,max]|全領域/99.9%|RSS全領域→99.9% MiB|','|---|---:|---:|---:|---:|---:|']
        rows=[s for s in d['summary'] if s['mode']=='adaptive' and suite in s['suites']]
        if suite=='strong':rows.sort(key=lambda s:s['ranks'])
        for s in rows:
            f=next(t for t in d['summary'] if t['shape']==s['shape'] and t['ranks']==s['ranks'] and t['mode']=='full')
            def fmt(t):
                v=t['rt_max_seconds'];return f"{v['min']:.4f} [{v['median']:.4f},{v['max']:.4f}]"
            lines.append(f"|{'×'.join(map(str,s['shape']))}|{s['ranks']}|{fmt(f)}|{fmt(s)}|{f['rt_max_seconds']['min']/s['rt_max_seconds']['min']:.3f}|{f['peak_rank_bytes']['median']/2**20:.1f}→{s['peak_rank_bytes']['median']/2**20:.1f}|")
        lines.append('')
    lines+=['## diagnostics','','|形状|MPI|打切り有効更新数（各33回中）|局所FFT対数/有効時全対数（3試行計）|局所FFT格子比|最大電流差 a.u.|最大エネルギー幅 Ha|','|---|---:|---|---:|---:|---:|---:|']
    for v in diagnostics:
        ratio='—' if v['local_fft_volume_ratio'] is None else f"{v['local_fft_volume_ratio']:.4f}"
        lines.append(f"|{'×'.join(map(str,v['shape']))}|{v['ranks']}|{v['active_refreshes']}|{v['local_pairs']}/{v['local_pairs']+v['global_pairs_during_active_support']}|{ratio}|{v['current_error_au']:.4e}|{v['maximum_energy_width_ha']:.4e}|")
    (a.output/'tables.md').write_text('\n'.join(lines)+'\n')
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,suite in enumerate(('weak','strong')):
        for mode,label in [('full','Full support'),('adaptive','99.9% requested')]:
            rows=[s for s in d['summary'] if s['mode']==mode and suite in s['suites']]
            if suite=='strong':rows.sort(key=lambda s:s['ranks'])
            xs=list(range(len(rows)));labels=['x'.join(map(str,s['shape'])) if suite=='weak' else str(s['ranks']) for s in rows]
            for row,key,scale,stat in [(0,'rt_max_seconds',1,'min'),(1,'peak_rank_bytes',1/2**20,'median')]:
                vals=[s[key] for s in rows];ys=[v[stat]*scale for v in vals]
                axs[row,col].errorbar(xs,ys,yerr=[[(v[stat]-v['min'])*scale for v in vals],[(v['max']-v[stat])*scale for v in vals]],marker='o',capsize=3,label=label)
                axs[row,col].set_xticks(xs,labels,rotation=35 if suite=='weak' else 0)
                axs[row,col].set_ylabel('16-step RT time (s)' if row==0 else 'Max rank lifetime RSS (MiB)')
                axs[row,col].set_xlabel('Supercell (MPI ranks = molecules)' if suite=='weak' else 'MPI ranks (4x4x1)')
                axs[row,col].grid(alpha=.25);axs[row,col].legend(fontsize=8)
        axs[0,col].set_yscale('log');axs[0,col].set_title(suite.capitalize()+' scaling')
    fig.suptitle('H2 PBEh(40), impulse 16 steps: paired same-binary measurements\n3 repeats; minimum time / median RSS; bars min-max; adaptive support sometimes falls back')
    fig.savefig(a.output/'scaling.png',dpi=180);fig.savefig(a.output/'scaling.svg')
    svg=a.output/'scaling.svg';svg.write_text('\n'.join(s.rstrip() for s in svg.read_text().splitlines())+'\n')

if __name__=='__main__':main()
