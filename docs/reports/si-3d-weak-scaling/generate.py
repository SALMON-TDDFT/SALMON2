#!/usr/bin/env python3
"""Generate fixed-work-per-fragment silicon GS/RT inputs; no jobs are submitted."""
from pathlib import Path
import hashlib,json,re
root=Path(__file__).resolve().parent
template=root/'templates'
# Conventional cubic diamond basis, in units of a/4.
basis=((0,0,0),(1,1,1),(2,0,2),(0,2,2),(2,2,0),(3,1,3),(1,3,3),(3,3,1))
a=10.26
manifest={'lattice_bohr':a,'core_grid':[16]*3,'buffer_grid':[8]*3,
 'fragment_grid':[32]*3,'fragment_atoms':64,'fragment_states':256,
 'radius_bohr':9,'dt_au':0.02,'nt':16,'ace_interval':1,'u_interval':1,
 'radius_control':'&functional hse_lcfo_wf_radius, fixed bohr',
 'initial_norm_warning_threshold':0.999,
 'pseudo_source':'samples/dc_hse/si64-chain/Si_rps.dat',
 'algorithm_control':'&functional namelist only; no SALMON_LCFO environment controls',
 'seed_distributed':True,
 'status':'generated and statically validated; not run on Fugaku','cases':[]}
assert (root/'Si_rps.dat').is_file(), 'Bundled Si pseudopotential is required'
for n in (4,6,8,10):
 case=root/f'{n}x{n}x{n}';case.mkdir(exist_ok=True)
 nf=n**3;na=8*nf;ne=4*na;occupied=ne//2
 atoms=[]
 for z in range(n):
  for y in range(n):
   for x in range(n):
    for b in basis:
     atoms.append("'Si' "+' '.join(f'{a*(i+t/4):.10f}' for i,t in zip((x,y,z),b))+' 1')
 (case/'atom.dat').write_text('\n'.join(atoms)+'\n')
 for stage in ('gs','rt'):
  folder=case/stage;folder.mkdir(exist_ok=True)
  s=(template/f'{stage}-inputfile').read_text()
  s=s.replace("sysname='diamond64_hse'",f"sysname='si{na}_hse'")
  s=s.replace('num_fragment=8,1,1',f'num_fragment={n},{n},{n}')
  s=s.replace('num_rgrid_buffer=8,0,0','num_rgrid_buffer=8,8,8')
  s=s.replace('nproc_rgrid_tot=8,1,1',f'nproc_rgrid_tot={n},{n},{n}')
  s=s.replace('nstate_frag=64','nstate_frag=256')
  s=s.replace("lcfo_eigensolver='lapack'", "lcfo_eigensolver='chefsi'\n lcfo_diag_chefsi_filter_degree=60\n lcfo_diag_chefsi_max_cycle=200\n lcfo_diag_chefsi_residual_tolerance=1d-7")
  if stage=='rt':s=s.replace('nproc_rgrid=8,1,1',f'nproc_rgrid={n},{n},{n}')
  s=s.replace('al=53.76d0,6.72d0,6.72d0','al='+','.join(f'{a*n:.8f}d0' for _ in range(3)))
  s=re.sub(r'(?m)^ nstate=\d+',f' nstate={ne if stage=="gs" else occupied}',s)
  s=s.replace('nelec=256',f'nelec={ne}').replace('natom=64',f'natom={na}')
  s=s.replace("file_atom_coor='atom.dat'","file_atom_coor='../atom.dat'")
  s=s.replace("file_pseudo(1)='C_rps.dat'","file_pseudo(1)='../../Si_rps.dat'")
  s=s.replace('izatom(1)=6','izatom(1)=14').replace('lloc_ps(1)=1','lloc_ps(1)=2')
  s=s.replace('num_rgrid=128,16,16','num_rgrid='+','.join([str(16*n)]*3))
  # Algorithm choices belong to the input, not shell environment variables.
  s=s.replace(" yn_hse_wannier='n'\n",'')
  controls="\n yn_hse_lcfo_rt='{}'\n yn_hse_wannier='y'".format('y' if stage=='rt' else 'n')
  if stage=='rt':
   controls+="\n yn_hse_lcfo_direct_wf='y'\n yn_hse_lcfo_continuity='n'\n yn_hse_lcfo_fft_measure='n'"
   controls+="\n yn_hse_lcfo_seed_distributed='y'\n hse_lcfo_ace_interval=1\n hse_lcfo_u_interval=1"
   controls+="\n hse_lcfo_fft_batch=1\n hse_lcfo_wf_radius=9d0"
  s=s.replace("xc='hse06'","xc='hse06'"+controls)
  (folder/'inputfile').write_text(s)
 link=case/'rt/data_dcdft'
 if not link.is_symlink() and not link.exists():link.symlink_to('../gs/data_dcdft',target_is_directory=True)
 manifest['cases'].append({'case':case.name,'fragments':nf,'mpi_ranks_gs':nf,'mpi_ranks_rt':nf,
  'atoms':na,'electrons':ne,'gs_states':ne,'rt_states':occupied,'grid':[16*n]*3,
  'cell_bohr':[round(a*n,8)]*3,'one_dense_complex_matrix_GiB':16*occupied**2/2**30,
  'baseline_initial_raw_local_plus_raw_GiB_per_rank':12*16*occupied**2/2**30,
  'tiled_links_root_storage_GiB':6*16*occupied**2/2**30,
  'tiled_links_nonroot_storage_GiB':0,
  'tiled_links_scratch_GiB_per_rank':16*((4096+2*occupied)*min(64,occupied)+4096)/2**30,
  'baseline_root_gamma_link_arrays_lower_bound_GiB':18*16*occupied**2/2**30,
  'root_gamma_link_arrays_lower_bound_GiB':6*16*occupied**2/2**30})
manifest['sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
 for p in sorted(root.rglob('*')) if p.is_file() and p.name in ('inputfile','atom.dat','Si_rps.dat')}
(root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Generated 4 GS/RT input pairs and coordinates')
