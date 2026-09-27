#!/usr/bin/env python3
"""Static consistency check, not a substitute for a SALMON/Fugaku run."""
from pathlib import Path
from collections import Counter
import hashlib,json,math,re
root=Path(__file__).resolve().parent;m=json.loads((root/'manifest.json').read_text())
for c in m['cases']:
 n=int(c['case'].split('x')[0]);nf=n**3;case=root/c['case']
 coords=[list(map(float,l.split()[1:4])) for l in (case/'atom.dat').read_text().splitlines()]
 assert len(coords)==8*nf and len(set(map(tuple,coords)))==8*nf
 assert all(0<=x<6.72*n for p in coords for x in p)
 counts=Counter(tuple(int(math.floor((x+1e-9)/6.72)) for x in p) for p in coords)
 assert len(counts)==nf and set(counts.values())=={8}
 for stage in ('gs','rt'):
  folder=case/stage;text=(folder/'inputfile').read_text()
  def value(key):
   z=re.findall(r'^\s*'+re.escape(key)+r'\s*=\s*(.*?)\s*$',text,re.M)
   assert len(z)==1,(stage,key,z)
   return z[0]
  def vector(key):return list(map(int,value(key).split(',')))
  assert vector('num_fragment')==vector('nproc_rgrid_tot')==[n]*3
  assert vector('num_rgrid')==[16*n]*3 and vector('num_rgrid_buffer')==[8]*3
  assert vector('nproc_rgrid')==([1]*3 if stage=='gs' else [n]*3)
  assert int(value('nproc_ob'))==int(value('nproc_k'))==1
  assert int(value('natom'))==8*nf and int(value('nelec'))==32*nf
  assert int(value('nstate'))==(32*nf if stage=='gs' else 16*nf)
  assert int(value('nstate_frag'))==256 and value('lcfo_eigensolver')=="'chefsi'"
  for key in ('file_atom_coor','file_pseudo(1)'):assert (folder/value(key).strip("'")).is_file()
  if stage=='rt':assert value('dt')=='0.02d0' and value('nt')=='16'
 for k in ('one_dense_complex_matrix_GiB','baseline_initial_raw_local_plus_raw_GiB_per_rank'):assert c[k]>0
 print(c['case'],c['atoms'],'atoms /',nf,'MPI; initial link matrices alone',round(c['baseline_initial_raw_local_plus_raw_GiB_per_rank'],3),'GiB/rank')
for name,digest in m['sha256'].items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
print('Static geometry, input, path and hash checks passed; no numerical jobs run')
