"""Analyze completed scalar-only force-work diagnostics; no GS files read."""
import argparse
import json
from pathlib import Path
import numpy as np
from measure_nlcc_force_work import ordered_samples, summarize, simpson, write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('results',type=Path)
    args=parser.parse_args();root=args.results
    def rows(mode,grid):
        meta=json.loads((root/mode/'provenance.json').read_text())
        return ordered_samples(json.loads((root/mode/f'grid{grid}'/'samples.json').read_text()),meta['intervals'])
    report={'attribution':'Consistency residual is not an isolated NLCC force correction. Grid-dependent paths differ.'}
    paths={};rigid={};spots={}
    for g in (16,24):
        length=10.26/g;entry={'path_length_bohr':length,'start_atom_x_bohr':.12,'end_atom_x_bohr':.12+length}
        for mode in ('nlcc','no-nlcc-control'):
            r=rows(mode,g)
            if [v['index'] for v in r]!=list(range(17)):raise ValueError('incomplete path: '+mode+str(g))
            entry[mode]=summarize(r,length)
        paths[str(g)]=entry
        r=rows('rigid-nlcc',g)
        if [v['index'] for v in r]!=list(range(9)):raise ValueError('incomplete rigid control')
        rigid[str(g)]=dict(energy_peak_to_peak_hartree=float(np.ptp([v['energy_hartree'] for v in r])),
                          **summarize(r,length))
    for mode in ('nlcc','no-nlcc-control'):
        baseline={r['index']:r for r in rows(mode,24)}
        tight=rows('tight-'+mode,24)
        if [r['index'] for r in tight]!=[0,16]:raise ValueError('incomplete SCF spot check')
        differences=[t['energy_hartree']-baseline[t['index']]['energy_hartree'] for t in tight]
        spots[mode]=dict(max_abs_endpoint_energy_change_hartree=max(map(abs,differences)),
            delta_energy_change_hartree=differences[1]-differences[0],
            max_abs_endpoint_force_change_hartree_per_bohr=float(max(np.max(np.abs(
                np.asarray(t['forces_hartree_per_bohr'])-baseline[t['index']]['forces_hartree_per_bohr'])) for t in tight)),
            scope='endpoint SCF stability only; not a full force-integral tolerance study')
    lda={};points={}
    for xc in ('pz','libxc_pz'):
        a=rows('diagnostic-'+xc,16)[0];b=rows('diagnostic-'+xc+'-endpoint',16)[0]
        lda[xc]=dict(initial_energy_hartree=a['energy_hartree'],final_energy_hartree=b['energy_hartree'],
                     delta_energy_hartree=b['energy_hartree']-a['energy_hartree'])
        points[xc]=(a,b)
    lda['delta_energy_libxc_minus_builtin_hartree']=lda['libxc_pz']['delta_energy_hartree']-lda['pz']['delta_energy_hartree']
    lda['max_abs_endpoint_force_difference_hartree_per_bohr']=float(max(np.max(np.abs(
        np.array(a['forces_hartree_per_bohr'])-np.array(b['forces_hartree_per_bohr'])))
        for a,b in zip(points['pz'],points['libxc_pz'])))
    a=rows('diagnostic-no-core-pz',16)[0];b=rows('diagnostic-no-core-libxc_pz',16)[0]
    lda['no_core_initial_energy_difference_hartree']=b['energy_hartree']-a['energy_hartree']
    lda['scope']='Implementation comparison only: PBE-generated potential with PZ is not a physical XC-matched calculation.'
    report.update(paths=paths,rigid_translation=rigid,scf_endpoint_stability=spots,lda_implementation_diagnostic=lda)
    write_json(root/'analysis.json',report)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,8))
    for ax,g in zip(axs[0],(16,24)):
        for mode,color in (('nlcc','tab:blue'),('no-nlcc-control','tab:orange')):
            r=rows(mode,g);x=np.array([v['shift_bohr'] for v in r]);e=np.array([v['energy_hartree'] for v in r])
            f=np.array([v['path_force_hartree_per_bohr'] for v in r]);inds=np.arange(2,17,2)
            ax.plot(x,(e-e[0])*1000,color=color,label=mode+' energy')
            ax.plot(x[inds],[-simpson(f[:i+1],x[i])*1000 for i in inds],'.--',color=color,label=mode+' force work')
        ax.set(title=f'{g}³: one-atom displacement',xlabel='Displacement (bohr)',ylabel='Energy change (mHa)');ax.legend(fontsize=7)
    ax=axs[1,0]
    for g in (16,24):
        r=rows('rigid-nlcc',g);e=np.array([v['energy_hartree'] for v in r])
        ax.plot(np.linspace(0,1,len(r)),(e-e[0])*1000,'o-',label=f'{g}³')
    ax.set(title='Rigid translation: egg-box variation',xlabel='Shift / grid spacing',ylabel='Energy change (mHa)');ax.legend()
    ax=axs[1,1];labels=[];values=[];quad=[]
    for g in (16,24):
        for mode in ('nlcc','no-nlcc-control'):
            s=paths[str(g)][mode];labels.append(f'{g}³\n'+('NLCC' if mode=='nlcc' else 'core off'))
            values.append(s['consistency_residual_mHa_per_atom']);quad.append(s['quadrature_change_hartree']*1000/8)
    ax.bar(labels,values,yerr=quad,color=['tab:blue','tab:orange']*2,capsize=4)
    ax.set(title='ΔE + integral(F dx)',ylabel='Consistency residual (mHa/atom)')
    ax.axhline(0,color='black',lw=.7)
    fig.suptitle('Si8 force/energy consistency — diagnostic results')
    fig.text(.5,.015,'Different grids use different physical paths. Error bars: |I16 − I8| only. Core-off input is a diagnostic control.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.035,1,.96));fig.savefig(root/'comparison.png',dpi=180);plt.close(fig)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
