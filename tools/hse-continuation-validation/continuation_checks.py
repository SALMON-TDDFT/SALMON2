import struct, math, array, sys

def finite(b):
 assert len(b)%8==0
 a=array.array('d');a.frombytes(b)
 if sys.byteorder!='little':a.byteswap()
 assert all(math.isfinite(x) for x in a),'nonfinite payload'

def check_endpoints(steps,svds,nt,k,start=0):
 assert k>0 and nt>0 and 0<=start<nt
 expected=list(range(0 if start==0 else (start//k+1)*k,nt+1,k))
 assert steps==expected,('teacher schedule',len(steps),nt,k,start)
 assert svds==nt//k-start//k,('SVD count',svds,nt,k,start)
 return svds

def check_factor_v3(b,step,nr,no,nk,rank,k,warmup):
 ng=nr**3
 assert len(b)>=544
 h=struct.unpack_from('<20i',b,0)
 assert h[:2]==(3,step)
 assert h[2:17]==(1,1,1,nr,nr,nr,no,nk,1,nk,1,rank+1,1,k,warmup),h
 assert h[17:]==(ng,no,1)
 dv,condition,rows,ready=struct.unpack_from('<2d2i',b,80)
 assert math.isfinite(dv) and dv>0 and math.isfinite(condition) and rows==0 and ready in (0,1)
 assert struct.unpack_from('<3d',b,104)==(0.,0.,0.)
 assert struct.unpack_from('<5i',b,128)==(1,no,rank,0,0)
 assert all(not b[148+64*i:212+64*i].strip() for i in range(6))
 assert struct.unpack_from('<3i',b,532)==(1,1,1)
 off=544+ng*no*16
 finite(b[544:off])
 hh=struct.unpack_from('<9i',b,off)
 assert hh==(1,ng,no,3,k,3,step-step%k,step//k-2,ready),hh
 assert bool(ready)==(hh[7]>=8)
 off+=36
 finite(b[off:off+64]);off+=64
 assert len(b)==off+ng*no*3*16
 finite(b[off:])
 return {'version':3,'step':step,'last':hh[6],'teachers':hh[7],'ready':bool(ready),'dv':dv,'condition':condition,'bytes':len(b)}

def check_pulse(rows,dt=.03,omega=.0569614493722652368,tw=882.447392251770225,intensity=1e12,start=0):
 def ac(t):
  tt=t-tw/2
  return -5.338e-9*math.sqrt(intensity)/omega*math.cos(math.pi*tt/tw)**2*math.sin(omega*tt) if abs(tt)<tw/2 else 0.
 err=0.
 for i,row in enumerate(rows,start+1):
  assert len(row)==16 and all(math.isfinite(x) for x in row)
  assert abs(row[0]-i*dt)<1e-6
  av=ac(i*dt);ev=-(av-ac((i-1)*dt))/dt
  expected=[0.,0.,av,0.,0.,ev,0.,0.,av,0.,0.,ev]
  diff=max(abs(x-y) for x,y in zip(row[1:13],expected));assert diff<1e-10,('pulse phase',i,diff)
  err=max(err,diff)
 return err
