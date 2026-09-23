"""Atomic, self-contained restart archives for the constant-A PT-HSE driver."""
import hashlib,json,os,tempfile
from pathlib import Path
import numpy as np

FORMAT='SALMON_PT_HSE_RESTART_V1'

def fingerprint(paths):
 h=hashlib.sha256()
 for p in sorted(map(Path,paths),key=lambda p:p.name):
  h.update(p.name.encode()+b'\0')
  with p.open('rb') as f:
   while block:=f.read(1024*1024):h.update(block)
 return h.hexdigest()

def validate(u,g,m,rows):
 if m.get('format')!=FORMAT:raise ValueError('Unknown restart format')
 step=m.get('step')
 if not isinstance(step,int) or isinstance(step,bool) or step<0 or len(rows)!=step:raise ValueError('Checkpoint step/history mismatch')
 if not np.isfinite([m['dt'],m['amplitude'],m['initial_energy']]).all() or m['dt']<=0:raise ValueError('Invalid restart physics')
 if u.ndim<3 or g.shape!=(u.shape[0],u.shape[1],u.shape[1]) or not np.isfinite(u).all() or not np.isfinite(g).all():
  raise ValueError('Invalid restart arrays')
 if not np.allclose(g.conj().transpose(0,2,1)@g,np.eye(g.shape[-1]),atol=1e-9,rtol=0):raise ValueError('Nonunitary MLWF gauge')
 for i,row in enumerate(rows,1):
  if row['step']!=i or not np.isclose(row['time_au'],i*m['dt'],rtol=1e-12,atol=1e-14):raise ValueError('Inconsistent restart history')

def save_checkpoint(path,u,g,metadata,rows,previous=None):
 path=Path(path);validate(np.asarray(u),np.asarray(g),metadata,rows);temporary=None
 if previous is None:previous=np.array([a.reshape(u.shape[1],-1,order='F').T for a in u])@g
 if previous.shape!=(u.shape[0],int(np.prod(u.shape[2:])),u.shape[1]) or not np.isfinite(previous).all():raise ValueError('Invalid localization reference')
 try:
  with tempfile.NamedTemporaryFile(dir=path.parent,prefix='.restart-',suffix='.npz',delete=False) as f:
   temporary=Path(f.name)
   np.savez_compressed(f,u=u,gauge=g,localizer_previous=previous,metadata=json.dumps(metadata),trajectory=json.dumps(rows))
   f.flush();os.fsync(f.fileno())
  os.replace(temporary,path)
 finally:
  if temporary is not None and temporary.exists():temporary.unlink()

def load_checkpoint(path,expected):
 with np.load(path,allow_pickle=False) as d:
  u=d['u'].copy();g=d['gauge'].copy();previous=d['localizer_previous'].copy();m=json.loads(str(d['metadata']));rows=json.loads(str(d['trajectory']))
 validate(u,g,m,rows)
 for k,v in expected.items():
  if m.get(k)!=v:raise ValueError(f'Restart mismatch: {k}')
 if previous.shape!=(u.shape[0],int(np.prod(u.shape[2:])),u.shape[1]) or not np.isfinite(previous).all():raise ValueError('Invalid localization reference')
 return u,g,m,rows,previous

def write_json(path,value):
 path=Path(path);temporary=path.with_name('.'+path.name+'.tmp')
 temporary.write_text(json.dumps(value,indent=2)+'\n');os.replace(temporary,path)
