"""Direct occupied Wannier propagation with a matched field-free control."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MPLCONFIGDIR','/private/tmp/salmon-mpl-cache')
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wannier import read_u,geometry,initial_gauge,reconstruct,moments
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
work=Path('/private/tmp/salmon-si-time-wannier');fs=.02418884326505;ang=.52917721067
r,c,k,b,L,h=geometry(root);dv=h**3;lengths=np.full(3,4*L);pos=(c[:,None,:]+r[None,:,:]).reshape(-1,3)
u0,_=read_u(root/'calculations/si_tdcdft_k4/gs/data_for_restart');g,s=initial_gauge(u0,k,r,b,L,dv)
np.savez_compressed(out/'initial_gauge.npz',gauge=g,singular_values=s,bond_centers=b)
summary=dict(gauge='Bond-Gaussian projection and symmetric orthonormalization; not MLWF',trial_sigma_bohr=1.2,
             min_projection_singular_value=float(s.min()),supercell_length_A=float(4*L*ang),cases={})
steps=np.arange(0,1601,200);traces={};clouds={};quality={}
for name in ('none','weak','strong'):
 status=json.loads((out/f'{name}_status.json').read_text());assert status['complete'] and status['exit_code']==0
 assert 'end SALMON' in (work/name/'outputfile').read_text()
 traces[name]=[];quality[name]=dict(max_orthogonality_error=0.,max_density_reconstruction_error=0.,max_norm_error=0.)
 for step in steps:
  u=u0 if step==0 else read_u(work/name/f'checkpoint_rt_{step:06d}')[0]
  free=u0 if step==0 else read_u(work/'none'/f'checkpoint_rt_{step:06d}')[0]
  orth=np.einsum('krn,krm->knm',u.conj(),u)*dv
  quality[name]['max_orthogonality_error']=max(quality[name]['max_orthogonality_error'],float(np.max(abs(orth-np.eye(16)))))
  overlap=np.einsum('krn,krm->knm',u0.conj(),free)*dv
  a,sv,z=np.linalg.svd(overlap);undo=(a@z).conj().transpose(0,2,1)
  result={}
  for mode,rotation in [('fixed_initial',g),('field_free_aligned',undo@g)]:
   v=u@rotation;w=reconstruct(v,k,r,c);rho=abs(w.reshape(-1,16))**2
   reconstructed=np.sum(abs(w)**2,axis=(0,2));target=np.mean(np.sum(abs(u)**2,axis=2),axis=0)
   err=float(np.max(abs(reconstructed-target)));quality[name]['max_density_reconstruction_error']=max(quality[name]['max_density_reconstruction_error'],err)
   assert err<1e-12
   mm=moments(rho,pos,lengths,dv)
   quality[name]['max_norm_error']=max(quality[name]['max_norm_error'],float(np.max(abs(mm['norm']-1))))
   assert np.max(abs(mm['norm']-1))<1e-6
   result[mode]=mm
   if mode=='field_free_aligned' and step in (0,800,1600):clouds[f'{name}_{step}']=rho[:,0].astype(np.float32)
  traces[name].append(result)
 nex=np.loadtxt(work/name/'Si_nex.data')
 (out/f'{name}_nex.data').write_text(''.join(line.rstrip()+'\n' for line in (work/name/'Si_nex.data').read_text().splitlines()))
 (out/f'{name}_rt.data').write_text(''.join(line.rstrip()+'\n' for line in (work/name/'Si_rt.data').read_text().splitlines()))
 (out/f'{name}_rt_xc.data').write_text(''.join(line.rstrip()+'\n' for line in (work/name/'Si_rt_xc.data').read_text().splitlines()))
 summary['cases'][name]=dict(quality=quality[name],excited_electrons_final=float(nex[-1,1]),excited_holes_final=float(nex[-1,2]),
    excited_electrons_cm3=float(nex[-1,1]/(L*ang*1e-8)**3),projection_missing_electrons=float(nex[-1,2]-nex[-1,1]))
 if name!='none':
  refdir=root/'calculations/si_tdcdft_k4/bare_k'/('weak' if name=='weak' else 'strong_long')
  reference=np.loadtxt(refdir/'Si_rt.data');current=np.loadtxt(work/name/'Si_rt.data');n=min(len(reference),len(current))
  err=float(np.max(abs(current[:n,13:16]-reference[:n,13:16])))
  summary['cases'][name]['max_current_difference_from_existing_run']=err
  assert err<1e-9
 np.savez_compressed(out/f'{name}_moments.npz',steps=steps,**{f'{mode}_{key}':np.array([item[mode][key] for item in traces[name]]) for mode in ('fixed_initial','field_free_aligned') for key in ('norm','center','spread','spread_axes','concentration','boundary_weight')})
 for mode in ('fixed_initial','field_free_aligned'):
  mm=traces[name][-1][mode];initial=traces[name][0][mode]
  delta=(mm['center']-initial['center']+lengths/2)%lengths-lengths/2
  summary['cases'][name][mode]=dict(initial_mean_spread_A2=float(initial['spread'].mean()*ang**2),final_mean_spread_A2=float(mm['spread'].mean()*ang**2),
     final_rms_radius_A=float(np.sqrt(mm['spread'].mean())*ang),final_mean_center_shift_A=(delta.mean(axis=0)*ang).tolist(),
     final_rms_center_shift_A=float(np.sqrt(np.mean(np.sum(delta**2,axis=1)))*ang),max_boundary_weight=float(mm['boundary_weight'].max()),min_circular_concentration=float(mm['concentration'].min()))
rows=[]
for name in traces:
 for step,item in zip(steps,traces[name]):
  for mode in item:
   mm=item[mode];rows.append([name,mode,int(step),step*.08*fs,float(mm['spread'].mean()*ang**2),float(mm['boundary_weight'].max())])
with (out/'moments.csv').open('w') as f:
 f.write('case,gauge,step,time_fs,mean_spread_A2,max_boundary_weight\n')
 for row in rows:f.write(','.join(map(str,row))+'\n')
np.savez_compressed(out/'representative_density.npz',positions_bohr=pos,**clouds)
(out/'metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,ax=plt.subplots(3,1,figsize=(9,10),sharex=True,layout='constrained')
colors={'none':'0.4','weak':'tab:blue','strong':'tab:red'}
for name in traces:
 for i,mode in enumerate(('fixed_initial','field_free_aligned')):
  ax[i].plot(steps*.08*fs,[item[mode]['spread'].mean()*ang**2 for item in traces[name]],'o-',ms=4,label=name,color=colors[name])
 nex=np.loadtxt(out/f'{name}_nex.data');ax[2].plot(nex[:,0]*fs,nex[:,1],label=name,color=colors[name])
ax[0].set_title('Direct propagation: fixed initial Wannier gauge');ax[1].set_title('Same field-free occupied rotation removed from all cases')
for a in ax[:2]:a.set_ylabel('Mean orbital spread (A^2)');a.legend()
ax[2].set_ylabel('Excited electrons / 8-atom cell');ax[2].set_xlabel('Time (fs)')
for a in ax:a.axvline(60*fs,color='0.5',ls=':');a.grid(alpha=.2)
fig.suptitle('Si 4x4x4 k: time-dependent projected Wannier orbitals\nSpread is gauge-dependent; not an exciton radius or screening coefficient')
fig.savefig(out/'dynamics.png',dpi=150);fig.savefig(out/'dynamics.pdf')
# One initial bond orbital and its final pump-induced density change, z integrated.
xy=(pos[:,:2]/h+.1).astype(int);ix,iy=xy[:,0],xy[:,1]
r0=clouds['none_0'];rf=clouds['none_1600'];rs=clouds['strong_1600'];images=[]
for rho in (r0,rf,rs,rs-rf):
 im=np.zeros((48,48));np.add.at(im,(ix,iy),rho*h);images.append(np.roll(np.roll(im,18,axis=0),18,axis=1))
fig,ax=plt.subplots(1,4,figsize=(14,3.5),layout='constrained')
for i,(a,im,title) in enumerate(zip(ax,images,['Initial','No pump, final','Strong pump, final','Strong - no pump'])):
 kwargs=dict(origin='lower',extent=[-18*h*ang,30*h*ang,-18*h*ang,30*h*ang])
 if i==3:kwargs.update(cmap='RdBu_r',vmin=-abs(im).max(),vmax=abs(im).max())
 else:kwargs.update(cmap='magma',vmin=0,vmax=max(z.max() for z in images[:3]))
 a.imshow(im.T,**kwargs);a.set_title(title);a.set_xlabel('x (A)');a.set_ylabel('y (A)')
fig.suptitle('Representative bond Wannier density, integrated along z; field-free-aligned gauge')
fig.savefig(out/'density.png',dpi=150);fig.savefig(out/'density.pdf')
print(json.dumps(summary,indent=2))
