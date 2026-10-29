from pathlib import Path
import struct,json,hashlib,math
from continuation_checks import finite,check_factor_v3

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as q:
  for b in iter(lambda:q.read(1048576),b''):h.update(b)
 return h.hexdigest()
def records(p):
 b=p.read_bytes();o=0;out=[]
 while o<len(b):
  n=struct.unpack_from('<i',b,o)[0];assert n>=0 and o+n+8<=len(b)
  assert struct.unpack_from('<i',b,o+4+n)[0]==n
  out.append(b[o+4:o+4+n]);o+=n+8
 assert o==len(b)
 return out

def validate(cp,step):
 expected={'Vh_stock.bin','atomic_coor.txt','atomic_vel.txt','factor-mode.bin','hse_restart.bin','info.bin','occupation.bin','rho0_s.bin','rt-initial-energy.bin','symmetry_restart.bin','wfn.bin'}|{'factor-history-%08d.bin'%i for i in range(80)}
 assert {p.name for p in cp.iterdir()}==expected
 assert all(p.is_file() and not p.is_symlink() for p in cp.iterdir())
 manifest={p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in cp.iterdir()}
 assert [struct.unpack('<i',x)[0] for x in records(cp/'info.bin')]==[80,16,step,80,0]
 occ=records(cp/'occupation.bin');assert len(occ)==1 and len(occ[0])==80*16*8;finite(occ[0]);assert all(abs(x[0]-2)<1e-12 for x in struct.iter_unpack('<d',occ[0]))
 rr=records(cp/'factor-mode.bin');assert len(rr)==1 and struct.unpack('<3i',rr[0])==(1,step,2)
 rr=records(cp/'rt-initial-energy.bin');assert len(rr)==1 and len(rr[0])==16
 z=struct.unpack('<2id',rr[0]);assert z[:2]==(1,step) and math.isfinite(z[2])
 p=cp/'wfn.bin';assert p.stat().st_size==4096*16*80*16
 with p.open('rb') as q:
  for b in iter(lambda:q.read(1048576),b''):finite(b)
 for name,n in [('Vh_stock.bin',65536),('rho0_s.bin',32768)]:
  rr=records(cp/name);assert sum(map(len,rr))==n
  for b in rr:finite(b)
 rr=records(cp/'hse_restart.bin');assert [len(x) for x in rr]==[20,18,32,19424,4,260230]
 assert struct.unpack('<5i',rr[0])==(6,80,16,8,1);finite(rr[3]);assert struct.unpack('<i',rr[4])[0]==len(rr[5])
 rr=records(cp/'symmetry_restart.bin');assert [len(x) for x in rr]==[12,8704];finite(rr[1])
 for name in ['atomic_coor.txt','atomic_vel.txt']:
  numbers=[]
  for token in (cp/name).read_text().replace(',',' ').split():
   try:x=float(token.replace('D','E').replace('d','e'))
   except ValueError:continue
   assert math.isfinite(x);numbers.append(x)
  assert len(numbers)>=24
 for rank in range(80):
  z=check_factor_v3((cp/('factor-history-%08d.bin'%rank)).read_bytes(),step,16,16,80,rank,16,160)
  assert abs(z['dv']-(10.26/16)**3)<1e-12
 for p in cp.iterdir():assert sha(p)==manifest[p.name]['sha256']
 return {'pass':True,'step':step,'native_schema_version':3,'binary_payload_finite':True,'files':manifest,'automatic_resubmit':False}
