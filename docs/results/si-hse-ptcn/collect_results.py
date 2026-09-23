"""Collect and compare equal-final-time PT-CN pilots with full RK4."""
from pathlib import Path
import json
import numpy as np
root=Path(__file__).resolve().parents[3];p=root/'calculations/si_hse_reference';out=Path(__file__).resolve().parent
cases={}
for name in ['dt032','dt016','dt008','rk4_032','zero032','midpoint032']:
 folder=p/f'ptcn_{name}'
 cases[name]={n:json.loads((folder/f'{n}.json').read_text()) for n in ['result','trajectory']}
v=np.load(p/'ptcn_rk4_032/state.npz')['u'];rho_ref=np.sum(abs(v)**2,axis=(0,1))
jref=np.array(cases['rk4_032']['trajectory'][-1]['current']);comparisons={}
for name in ['dt032','dt016','dt008','midpoint032']:
 u=np.load(p/f'ptcn_{name}/state.npz')['u'];rho=np.sum(abs(u)**2,axis=(0,1))
 j=np.array(cases[name]['trajectory'][-1]['current'])
 # Do not compare raw PT orbitals directly to Schrodinger orbitals: gauges differ.
 a=u.reshape(64,16,-1);b=v.reshape(64,16,-1);dv=.855**3
 cross=a.conj()@b.transpose(0,2,1)*dv
 left,_,right=np.linalg.svd(cross);rotation=left@right
 aligned=rotation.transpose(0,2,1)@a
 comparisons[name]=dict(current_relative_difference=float(np.linalg.norm(j-jref)/np.linalg.norm(jref)),
  current_absolute_difference=float(np.linalg.norm(j-jref)),density_relative_difference=float(np.linalg.norm(rho-rho_ref)/np.linalg.norm(rho_ref)),
  gauge_aligned_orbital_relative_difference=float(np.linalg.norm(aligned-b)/np.linalg.norm(b)),
  propagation_seconds=sum(row.get('step_wall_seconds',row['wall_seconds']) for row in cases[name]['trajectory']))
summary=dict(cases=cases,comparisons=comparisons,reference='Full self-consistent RK4, 4 steps of 0.08 au',final_time_au=.32,
 final_time_fs=.32*.024188843265857,scope='Short-pilot errors relative to stated finite-step RK4; not converged spectra')
(out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(comparisons,indent=2))
