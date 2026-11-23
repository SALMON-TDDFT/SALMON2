"""Compare old, XC-only, and XC+force scalar diagnostics (no GS files)."""
import argparse,json
from pathlib import Path
import numpy as np
from measure_nlcc_force_work import ordered_samples,summarize,simpson,write_json

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('baseline',type=Path);p.add_argument('results',type=Path)
a=p.parse_args();root=a.results

def samples(base,mode,g):
    rows=ordered_samples(json.loads((base/mode/f'grid{g}/samples.json').read_text()),16)
    if [r['index'] for r in rows]!=list(range(17)):raise ValueError('incomplete path')
    return rows

report={'scope':'Unpolarized CPU NC; discrete derivatives validated. Full path force/energy residuals are not physically converged.'}
paths={}
for g in (16,24):
    datasets={'old':samples(a.baseline,'nlcc',g),'xc-only':samples(root,'xc-only',g),
              'xc-and-force':samples(root,'xc-and-force',g),'core-off-control':samples(a.baseline,'no-nlcc-control',g)}
    entry={name:summarize(rows,10.26/g) for name,rows in datasets.items()}
    xc=datasets['xc-only'];final=datasets['xc-and-force']
    e_delta=np.array([f['energy_hartree']-x['energy_hartree'] for x,f in zip(xc,final)])
    if np.max(np.abs(e_delta))!=0:raise AssertionError('force-only addition changed SCF energies')
    df=np.array([f['forces_hartree_per_bohr'] for f in final])-np.array([x['forces_hartree_per_bohr'] for x in xc])
    entry['core_force_component']=dict(max_abs_all_atoms_hartree_per_bohr=float(np.max(np.abs(df))),
        max_abs_moved_x_hartree_per_bohr=float(np.max(np.abs(df[:,0,0]))),
        moved_x_integral_hartree=simpson(df[:,0,0],10.26/g),
        max_scf_energy_change_hartree=float(np.max(np.abs(e_delta))))
    paths[str(g)]=entry
report['paths']=paths
lda={}
for xc in ('pz','libxc_pz'):
    rows=ordered_samples(json.loads((root/('lda-'+xc)/'grid16/samples.json').read_text()),16)
    if [r['index'] for r in rows]!=[0,16]:raise ValueError('incomplete LDA diagnostic')
    lda[xc]=dict(energy=[r['energy_hartree'] for r in rows],force=[r['forces_hartree_per_bohr'] for r in rows])
lda['deformation_energy_difference_hartree']=(lda['libxc_pz']['energy'][1]-lda['libxc_pz']['energy'][0])-(lda['pz']['energy'][1]-lda['pz']['energy'][0])
lda['max_abs_force_difference_hartree_per_bohr']=float(np.max(np.abs(np.array(lda['libxc_pz']['force'])-lda['pz']['force'])))
lda['scope']='Implementation-only PZ comparison using the same PBE-generated UPF, not an XC-matched physical prediction.'
report['lda_equivalence']=lda
report['mpi']=json.loads((root/'mpi-crosscheck/results.json').read_text())
report['no_core_regression']=json.loads((root/'no-core-regression/results.json').read_text())
write_json(root/'analysis.json',report)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axs=plt.subplots(2,2,figsize=(11,8))
g=24;xc=samples(root,'xc-only',g);final=samples(root,'xc-and-force',g)
x=np.array([r['shift_bohr'] for r in xc]);e=np.array([r['energy_hartree'] for r in xc]);inds=np.arange(2,17,2)
ax=axs[0,0];ax.plot(x,(e-e[0])*1000/8,'k-',label='Energy change (same for both)')
for label,rows,color in [('XC only',xc,'tab:orange'),('XC + core force',final,'tab:blue')]:
    f=np.array([r['path_force_hartree_per_bohr'] for r in rows])
    ax.plot(x[inds],[-simpson(f[:i+1],x[i])*1000/8 for i in inds],'o--',color=color,label=label+' work')
ax.set(title='24³: energy vs force work',xlabel='Displacement (bohr)',ylabel='mHa / atom (Si8)');ax.legend(fontsize=8)
ax=axs[0,1]
f=np.array([r['path_force_hartree_per_bohr'] for r in final])-np.array([r['path_force_hartree_per_bohr'] for r in xc])
ax.plot(x,f,'o-',color='tab:blue');ax.axhline(0,color='black',lw=.7)
ax.set(title='Explicit NLCC force on displaced atom',xlabel='Displacement (bohr)',ylabel='Added x force (Ha/bohr)')
ax=axs[1,0];locations=np.arange(2);width=.19
for offset,(name,color) in enumerate([('old','gray'),('xc-only','tab:orange'),('xc-and-force','tab:blue'),('core-off-control','tab:green')]):
    vals=[paths[str(g)][name]['consistency_residual_mHa_per_atom'] for g in (16,24)]
    q=[paths[str(g)][name]['quadrature_change_hartree']*1000/8 for g in (16,24)]
    ax.bar(locations+(offset-1.5)*width,vals,width,label=name,color=color,yerr=q,capsize=2)
ax.set_xticks(locations,['16³','24³']);ax.set(title='Full-path consistency residual remains',ylabel='ΔE + integral(F dx) (mHa/atom)');ax.axhline(0,color='black',lw=.7);ax.legend(fontsize=7)
ax=axs[1,1]
before=json.loads((a.baseline/'analysis.json').read_text())['lda_implementation_diagnostic']['delta_energy_libxc_minus_builtin_hartree']
ax.bar(['Before','Corrected'],[abs(before)*1000,abs(lda['deformation_energy_difference_hartree'])*1000],color=['gray','tab:blue'])
ax.set_yscale('log');ax.set(title='Builtin PZ vs Libxc-PZ energy change',ylabel='Absolute difference (mHa / Si8 cell)')
fig.suptitle('NLCC correction: discrete consistency and scalar Si8 diagnostics')
fig.text(.5,.015,'Different grids use different paths. Bars: 16 vs 8 interval changes, not error bounds. Core-off is a different Hamiltonian.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.035,1,.96));fig.savefig(root/'comparison.png',dpi=180);plt.close(fig)
print(json.dumps({'core_force':{g:v['core_force_component'] for g,v in paths.items()},'lda_deformation_difference':lda['deformation_energy_difference_hartree']},indent=2))
