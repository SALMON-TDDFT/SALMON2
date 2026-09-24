"""Prepare optional native Si parity fixtures from existing converged full-k checkpoints.
Does not run SALMON or alter the original checkpoints.
"""
from pathlib import Path
import numpy as np,struct,shutil,re,json
root=Path(__file__).resolve().parents[2];base=root/'calculations/si_hse_native';out=base/'symmetry_validation';out.mkdir(exist_ok=True)
k=np.loadtxt(base/'gs_reference/Si_k.data',skiprows=5,max_rows=64)[:,1:4]
ops=np.loadtxt(root/'samples/exercise_si_symmetry/sym.dat').reshape(-1,3,4)
seen=set();reps=[];weights=[]
for i,v in enumerate(k):
 key=tuple(np.rint(v*8).astype(int))
 if key in seen:continue
 orbit={tuple(np.rint(op[:,:3]@v*8).astype(int)) for op in ops}
 seen|=orbit;reps.append(i);weights.append(len(orbit)/len(k))
assert len(reps)==12
(out/'representatives.json').write_text(json.dumps(dict(indices_zero_based=reps,weights=weights),indent=2))
def records(p):
 b=p.read_bytes();r=[];i=0
 while i<len(b):
  n=struct.unpack_from('i',b,i)[0];r.append(b[i+4:i+4+n]);assert struct.unpack_from('i',b,i+4+n)[0]==n;i+=n+8
 return r
def write(p,r):p.write_bytes(b''.join(struct.pack('i',len(a))+a+struct.pack('i',len(a)) for a in r))
for model,src in [('hse',base/'reference_restart'),('td',root/'calculations/si_tdcdft_k4/gs/data_for_restart')]:
 d=out/(model+'_reduced_gs');shutil.copytree(src,d,dirs_exist_ok=True)
 r=records(d/'info.bin');no=struct.unpack('i',r[1])[0];r[0]=struct.pack('i',12);write(d/'info.bin',r)
 wf=np.fromfile(src/'wfn.bin',np.complex128).reshape(64,no,1728);wf[reps].tofile(d/'wfn.bin')
 occ=np.frombuffer(records(src/'occupation.bin')[0],np.float64).reshape(64,no);write(d/'occupation.bin',[occ[reps].tobytes()])
 for reduced in [False,True]:
  run=out/(model+('_reduced' if reduced else '_full'));run.mkdir(exist_ok=True)
  template=base/('impulse_hse06_fixedalpha' if model=='hse' else 'impulse_tdcdft_fixedalpha')/'inputfile'
  s=template.read_text();s=re.sub(r'directory_read_data=.*',"directory_read_data='"+str(d if reduced else src)+"/'",s)
  s=re.sub(r'(?m)^(\s*nt\s*=).*',r'\g<1> 4',s);s=s.replace('dt = 0.16','dt = 0.08');s=re.sub(r'nproc_k=\d+','nproc_k=4',s);s=s.replace('checkpoint_interval=1000','checkpoint_interval=2')
  if reduced:s=s.replace("yn_periodic = 'y'","yn_periodic = 'y'\n yn_symmetry='yyn'")
  s=s.replace('nproc_k=4','nproc_k=3') if reduced else s
  (run/'inputfile').write_text(s);shutil.copy(root/'samples/exercise_si_symmetry/sym.dat',run/'sym.dat')
