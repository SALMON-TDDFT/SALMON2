from pathlib import Path
import json,re,math,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/work/si-laser-8x1x1-mesh12')
O=R.parent.parent/'SALMON2-pbeh40-rvv10/docs/reports/si-8x1x1-laser-5fs-2026-10-01';O.mkdir(parents=True,exist_ok=True)
C=[('PBE','pbe-laser'),('HSE WF 99.9% (old)','hse06-laser-omp-optimized'),('PBE0 WF 99.9%','pbe0-laser-latest'),('HSE WF+SR','hse06-sr-auto-laser-latest'),('HSE full EXX','hse06-full-exx-laser-latest')]
loads=json.loads((R/'external-load-events.json').read_text());ints=sorted((v['start_unix'],v['end_unix']) for v in loads if v.get('end_unix'));merged=[]
for a,b in ints:
 if merged and a<=merged[-1][1]:merged[-1][1]=max(b,merged[-1][1])
 else:merged.append([a,b])
def overlap(a,b):return sum(max(0,min(b,y)-max(a,x)) for x,y in merged)
D=[]
for label,name in C:
 p=R/name;s=json.loads((p/'status.json').read_text());log=(p/'output').read_text();j=np.loadtxt(p/'Si_rt.data');e=np.loadtxt(p/'Si_rt_energy.data')
 assert s['complete'] and s['returncode']==0 and len(j)==1477 and len(e)==1478 and np.isfinite(j).all() and np.isfinite(e).all() and 'end SALMON' in log
 ranks=[json.loads(q.read_text()) for q in sorted(p.glob('rank-????.json'))]
 if ranks:assert len(ranks)==8 and all(v['returncode']==0 for v in ranks)
 def timer(k):
  m=re.search('^'+re.escape(k)+r'\s+(\S+)\s+(\S+)',log,re.M);return float(m[2]) if m else None
 rt=timer('rt iterations');obs=json.loads((p/'step-observations.json').read_text());bins=[]
 for lo in np.arange(0,5,.25):
  z=[v for v in obs if lo<=v['time_fs']<lo+.25]
  if len(z)<2:continue
  a,b=z[0],z[-1];dt=b['elapsed_seconds']-a['elapsed_seconds'];dn=b['step']-a['step'];ov=overlap(s['started']+a['elapsed_seconds'],s['started']+b['elapsed_seconds'])
  if dn>0:bins.append(dict(fs=(a['time_fs']+b['time_fs'])/2,raw_seconds_step=dt/dn,clean_seconds_step=dt/dn if ov==0 else None,excluded_overlap_seconds=ov))
 # Diagnostics occur before their associated step row; use last diagnostic of each completed step.
 diag=[];state={};holds=0
 for line in log.splitlines():
  if 'retained accepted transported gauge' in line:holds+=1
  if line.startswith('EXX_ADAPTIVE fraction/'):
   v=list(map(float,line.split(':',1)[1].split()));state.update(radius=v[1],retained=1-v[2])
  if line.startswith('EXX_ADAPTIVE local/global'):
   v=list(map(int,line.split(':',1)[1].split()));state.update(local=v[0],global_pairs=v[1])
  if line.startswith('EXX_SR WF-local/'):
   v=list(map(int,line.split(':',1)[1].split()));state.update(wf_local=v[0],neighborhood=v[1])
  if line.startswith('EXX_SPATIAL refresh/iterations/'):
   v=line.split(':',1)[1].split();state.update(iterations=int(v[1]),mlwf_status=int(v[2]))
  t=line.split()
  if len(t)==7 and t[0].isdigit():
   try:fs=float(t[1]);norm=float(t[5])
   except ValueError:continue
   diag.append(dict(step=int(t[0]),fs=fs,norm_error=abs(norm-256),gauge_holds=holds,**state))
 row=dict(label=label,folder=name,binary_sha256=s['binary_sha256'],source_commit=s.get('source_commit'),wall_seconds=s['wall_seconds'],rt_timer_seconds=rt,rt_timer_seconds_step=rt/1477,reading_lda_seconds=timer('reading lda data'),rt_initialization_timer_seconds=timer('rt initialization'),wall_minus_rt_seconds=s['wall_seconds']-rt,rank_peak_rss_mib=[v/2**20 for v in s['rank_peak_rss_bytes']],development_overlap_seconds=overlap(s['started'],s['started']+s['wall_seconds']),ace_rejections=log.count('EXX_SUPPORT_ACE accepted: F'),gauge_holds=holds,max_norm_error=s['max_norm_error'],bins=bins,diagnostics=diag)
 D.append(dict(row=row,j=j,e=e))
 (O/(name+'-diagnostics.json')).write_text(json.dumps(diag,indent=2))
fig,ax=plt.subplots(3,1,figsize=(10,10),sharex=True)
for d in D:
 label=d['row']['label'];j=d['j'];e=d['e'];ax[0].plot(j[:,0]*.024188843265857,j[:,13],label=label);ax[1].plot(e[:,0]*.024188843265857,e[:,1]);ax[2].plot(e[:,0]*.024188843265857,e[:,2])
for a,y in zip(ax,['Jx (a.u.)','Total energy (Ha)','Energy change (Ha)']):a.set_ylabel(y);a.grid(alpha=.25)
ax[0].legend(fontsize=8);ax[-1].set_xlabel('Time (fs)');fig.tight_layout();fig.savefig(O/'current-energy.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(4,1,figsize=(10,12),sharex=True)
for d in D:
 row=d['row'];q=row['diagnostics'];t=[v['fs'] for v in q];label=row['label']
 if q and 'radius' in q[-1]:
  ax[0].plot(t,[v.get('radius',np.nan) for v in q],label=label);ax[1].plot(t,[v.get('retained',np.nan) for v in q]);ax[2].plot(t,[v.get('global_pairs',np.nan) for v in q]);ax[3].plot(t,[v['gauge_holds'] for v in q])
for a,y in zip(ax,['Maximum support R (bohr)','Minimum retained norm','Global FFT pairs','Cumulative gauge holds']):a.set_ylabel(y);a.grid(alpha=.25)
ax[0].legend(fontsize=8);ax[-1].set_xlabel('Time (fs)');fig.tight_layout();fig.savefig(O/'support-history.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(2,1,figsize=(10,8),sharex=True)
for d in D:
 row=d['row'];bins=row['bins'];line=ax[0].plot([v['fs'] for v in bins],[v['clean_seconds_step'] if v['clean_seconds_step'] is not None else np.nan for v in bins],'.-',label=row['label'])[0]
 dirty=[v for v in bins if v['clean_seconds_step'] is None];ax[0].scatter([v['fs'] for v in dirty],[v['raw_seconds_step'] for v in dirty],marker='x',color=line.get_color())
 q=row['diagnostics'];ax[1].plot([v['fs'] for v in q],[v.get('iterations',np.nan) for v in q])
ax[0].set_ylabel('Seconds/step (0.25 fs bins)');ax[0].set_yscale('log');ax[0].legend(fontsize=8);ax[1].set_ylabel('Last logged MLWF iterations');ax[1].set_xlabel('Time (fs)')
for a in ax:a.grid(alpha=.25)
fig.tight_layout();fig.savefig(O/'timing-iterations.png',dpi=160);plt.close(fig)
fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True);ref=D[-1];comparisons=[]
for d in [D[1],D[3]]:
 dj=d['j'][:,13]-ref['j'][:,13];de=d['e'][:,1]-ref['e'][:,1];dd=d['e'][:,2]-ref['e'][:,2]
 comparisons.append(dict(label=d['row']['label'],max_abs_Jx_difference=float(max(abs(dj))),Jx_relative_L2=float(np.linalg.norm(dj)/np.linalg.norm(ref['j'][:,13])),max_abs_energy_difference_Ha=float(max(abs(de))),initial_energy_difference_Ha=float(de[0]),final_energy_difference_Ha=float(de[-1]),max_abs_delta_energy_difference_Ha=float(max(abs(dd)))))
 ax[0].plot(d['j'][:,0]*.024188843265857,dj,label=d['row']['label']);ax[1].plot(d['e'][:,0]*.024188843265857,dd)
ax[0].legend();ax[0].set_ylabel('Jx - full (a.u.)');ax[1].set_ylabel('Delta E - full (Ha)');ax[1].set_xlabel('Time (fs)')
for a in ax:a.grid(alpha=.25)
fig.tight_layout();fig.savefig(O/'full-reference-errors.png',dpi=160);plt.close(fig)
# Detailed pair path / MLWF status / transverse currents, all retained in data.
fig,ax=plt.subplots(3,1,figsize=(10,9),sharex=True)
for d in D[1:]:
 q=d['row']['diagnostics'];t=[v['fs'] for v in q];label=d['row']['label'];ax[0].plot(t,[v.get('local',0) for v in q],label=label);ax[1].plot(t,[v.get('neighborhood',0) for v in q]);ax[2].plot(t,[v.get('mlwf_status',0) for v in q])
ax[0].legend(fontsize=8)
for a,y in zip(ax,['Local FFT pairs','SR neighborhood pairs','Last MLWF status']):a.set_ylabel(y);a.grid(alpha=.25)
ax[-1].set_xlabel('Time (fs)');fig.tight_layout();fig.savefig(O/'pair-paths.png',dpi=160);plt.close(fig)
summary=dict(conditions=[d['row'] for d in D],full_reference_comparisons=comparisons,load_intervals=loads)
(O/'summary.json').write_text(json.dumps(summary,indent=2));
if Path(__file__).resolve() != (O/'analyze.py').resolve():shutil.copy2(__file__,O/'analyze.py')
shutil.copy2(R/'external-load-events.json',O/'external-load-events.json')
print(json.dumps(dict(rows=[{k:v for k,v in d['row'].items() if k not in ['diagnostics','bins']} for d in D],comparisons=comparisons),indent=2))
