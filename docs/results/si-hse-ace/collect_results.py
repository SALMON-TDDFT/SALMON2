"""Collect small JSON from local checkpoints, excluding large orbital arrays."""
from pathlib import Path
import json
import numpy as np
root=Path(__file__).resolve().parents[3];p=root/'calculations/si_hse_reference';out=Path(__file__).resolve().parent
cases=['plus','halfdt','seeded','zero'];data={}
for case in cases:
 d=p/f'ace_{case}'
 data[case]={name:json.loads((d/f'{name}.json').read_text()) for name in ['result','trajectory']}
ref=json.loads((p/'rt_halfdt/trajectory.json').read_text())[-1]
jref=np.array(ref['current']); comparisons={}
for case in ['plus','halfdt']:
 row=data[case]['trajectory'][-1];u=np.load(p/f'ace_{case}/state.npz')['u'];v=np.load(p/'rt_halfdt/state.npz')['u']
 comparisons[case]=dict(current_relative_difference_vs_RK4_halfdt=float(np.linalg.norm(np.array(row['current'])-jref)/np.linalg.norm(jref)),
  orbital_relative_difference_vs_RK4_halfdt=float(np.linalg.norm(u-v)/np.linalg.norm(v)))
 # Orbital difference also measures harmless occupied gauge/phase error, not only density error.
 rho_u=np.sum(abs(u)**2,axis=(0,1));rho_v=np.sum(abs(v)**2,axis=(0,1))
 comparisons[case]['density_relative_difference_vs_RK4_halfdt']=float(np.linalg.norm(rho_u-rho_v)/np.linalg.norm(rho_v))
summary=dict(cases=data,comparisons=comparisons,reference_current=jref.tolist(),
 scope='Short ACE/midpoint pilots, not PT-CN or spectra; exact full exchange gate at every accepted midpoint')
(out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(comparisons,indent=2))
