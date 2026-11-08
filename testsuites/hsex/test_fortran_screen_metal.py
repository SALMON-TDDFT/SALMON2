"""Smooth fractional two-band fixture: full metallic W pipeline, no tolerance relaxation."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/hsex'))
from prepare_w import prepare
from mpi_test_support import shared_fixture,run_mpi

def main():
 exe=Path(sys.argv[1]);nr=int(sys.argv[2]) if len(sys.argv)>2 else 4
 threads=int(sys.argv[3]) if len(sys.argv)>3 else 12
 with shared_fixture() as path:
  r=path/'data_for_restart';r.mkdir()
  n=np.array([4]*3);mesh=np.array([6]*3);nk=216;ng=64;nb=64;nv=1;boxes=1;eta=1e-13
  cell=np.eye(3)*9;sites=np.indices((6,)*3).reshape(3,-1)[::-1].T
  frac=(sites-2.5)/6;k=frac*2*np.pi/9
  # Occupation decreases linearly with band energy, yielding smooth screening.
  x=1e-7*np.random.default_rng(544).uniform(-1,1,nk)
  base=np.linspace(-1e-4,1e-4,ng);coeff=np.r_[-np.ones(ng//2),np.ones(ng//2)]
  eig=base[None,:]+x[:,None]*coeff
  pop=2/ng-100*eig
  gsocc=np.zeros((nk,ng));gsocc[:,0]=2
  with (r/'kli-v1.bin').open('wb') as f:
   f.write(b'KLI1');np.r_[n,1,nk,ng,1].astype('i4').tofile(f)
   for a in (cell,k,gsocc,np.zeros(3),np.zeros(ng),[0.]):np.asarray(a,dtype='f8').tofile(f)
  u=np.eye(ng,dtype=complex);states=u[:,:nb].T*np.sqrt(ng/729)
  with (r/'wfn.bin').open('wb') as f:
   for j in range(nk):(u.T*np.sqrt(ng/729)).tofile(f)
  for j in range(nk):
   with (path/f'hsex-pump-{j+1:06d}.bin').open('wb') as f:
    f.write(b'HSP1');np.r_[n,nk,j+1,nb,80].astype('i4').tofile(f)
    for a in (cell,k[j],np.zeros(3),eig[j,:nb],pop[j]):np.asarray(a,dtype='f8').tofile(f)
    states.tofile(f);np.diag(eig[j]).astype('c16').tofile(f)
  np.savetxt(path/'Si_eigen.data',np.c_[np.tile(np.arange(ng)+1,nk),eig.ravel(),gsocc.ravel()],fmt=['%d','%.17e','%.1f'])
  ref=path/'reference'
  # Exact diagonal resolvent provides an independent vectorized reference.
  # The original Python assembly, eta validation, head closure and I/O are used.
  import pump_screen
  from box_factors import integer_phases
  def diagonal_target(job,**kwargs):
   _,nn,kk,cc,ff,order,mm,shifts,shells,j,etas=job
   lookup={tuple(t%mm):i for i,t in enumerate(order)}
   values=[]
   for qi,g in zip(shifts,shells):
    src=lookup[tuple((order[j]-qi)%mm)]
    phase=integer_phases(nn,g+np.rint(ff[src]+qi/mm-ff[j]).astype(int))
    delta=eig[j]-eig[src]
    values.append(np.array([phase.conj().T@((.5*pop[src]*delta/(delta**2+ee**2))[:,None]*phase) for ee in etas]))
   return values,0.
  pump_screen.solve_pump_target=diagonal_target
  prepare(path,ref,mesh,boxes,.2,qmesh=mesh,pump=path,eta=eta)
  with (path/'pipeline.bin').open('wb') as f:
   np.r_[n,mesh,mesh,boxes,nv,nb].astype('i4').tofile(f)
   np.r_[.2,eta,cell.T.ravel()].astype('f8').tofile(f)
   for a in (k,eig[:,:nb],pop):np.asarray(a,dtype='f8').tofile(f)
   for j in range(nk):u[:,:nb].T.tofile(f)
   for j in range(nk):np.diag(eig[j]).astype('c16').tofile(f)
  # The 64-band all-fractional fixture also exercises compensated refinement.
  run_mpi(exe,nr,threads,path,'metal',timeout=600)
  actual=np.fromfile(path/'pipeline-output.bin','f8')
  offset=40+8*(9+3*nk+3+ng+1+3*nk)
  expected=np.frombuffer(ref.read_bytes(),'f8',count=nk,offset=offset)
  np.testing.assert_allclose(actual,expected,atol=2e-12,rtol=2e-11)
  print('FRACTIONAL_METAL_MAX_ERROR',np.max(abs(actual-expected)))
 print('FORTRAN_SCREEN_METAL_PASS')
if __name__=='__main__':main()
