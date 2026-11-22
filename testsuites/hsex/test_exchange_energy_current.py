from pathlib import Path
import os,shlex,subprocess,json,argparse,tempfile
import numpy as np
root=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description='Check actual HSEX current against fixed-orbital energy derivatives.')
parser.add_argument('--build',type=Path,required=True)
parser.add_argument('--ranks',type=int,nargs='+',default=[1,4])
parser.add_argument('--report',type=Path)
args=parser.parse_args()
build=args.build.resolve()
work=tempfile.TemporaryDirectory(prefix='hsex-energy-current-')
tmp=Path(work.name)
driver='''program check_current
use mpi
use structures, only:s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar
use salmon_global, only:num_kgrid,cval
use kli_hsex, only:hsex_init,hsex_stage,hsex_current,hsex_induced_energy
implicit none
type(s_dft_system)::sys
type(s_rgrid)::mg
type(s_parallel_info)::info
type(s_orbital)::orb
type(s_scalar)::xc
complex(8)::initial(64,2,4),excited(64,2,4)
real(8)::a(3),ac(3),j(3),energies(3,3,2),delta(3),e0
integer::rank,np,ierr,i,axis,sgn,d,provided
call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
if(provided<MPI_THREAD_FUNNELED)error stop 'MPI funneled threads unavailable'
call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
if(np>4)error stop 'test ranks must not exceed four'
sys%no=2;sys%nspin=1;sys%nk=4;sys%nion=1;sys%hvol=42d0/64
sys%primitive_a=0d0
do i=1,3
 sys%primitive_a(i,i)=2.5d0+.5d0*i
enddo
allocate(sys%vec_k(3,4),sys%Rion(3,1));sys%Rion=0d0
sys%vec_k=0d0
do i=1,4
 sys%vec_k(:,i)=[.03d0+(i-1)*acos(-1d0)/6,.02d0,-.01d0]
enddo
mg%is=1;mg%ie=4
info%ik_s=1+rank*4/np;info%ik_e=(rank+1)*4/np;info%numk=info%ik_e-info%ik_s+1
info%isize_ko=np;info%id_ko=rank;info%icomm_ko=MPI_COMM_WORLD
allocate(orb%zwf(4,4,4,1,2,info%ik_s:info%ik_e,1),xc%f(4,4,4));xc%f=0d0
open(11,file='orbitals.bin',form='unformatted',access='stream');read(11)initial,excited;close(11)
do i=info%ik_s,info%ik_e
 orb%zwf(:,:,:,1,:,i,1)=reshape(initial(:,:,i),[4,4,4,2])
enddo
num_kgrid=[4,1,1];cval=1d0
call hsex_init(sys,mg,info,orb,0d0,.true.,xc)
do i=info%ik_s,info%ik_e
 orb%zwf(:,:,:,1,:,i,1)=reshape(excited(:,:,i),[4,4,4,2])
enddo
a=[.13d0,.07d0,-.04d0];delta=[1d-3,1d-4,1d-5]
call hsex_stage(sys,mg,info,orb,a,2);j=0d0;call hsex_current(j);e0=hsex_induced_energy()
do d=1,3
 do axis=1,3
  do sgn=1,2
   ac=a;ac(axis)=ac(axis)+delta(d)*(2*sgn-3)
   call hsex_stage(sys,mg,info,orb,ac,2)
   energies(axis,d,sgn)=hsex_induced_energy()
  enddo
 enddo
enddo
if(rank==0)then
 open(11,file='energy-current.bin',form='unformatted',access='stream')
 write(11)j,e0,energies;close(11)
endif
call MPI_Finalize(ierr)
end program
'''
(tmp/'driver.f90').write_text(driver)
rng=np.random.default_rng(20261112)
initial=[];excited=[]
for _ in range(4):
 q=np.linalg.qr(rng.normal(size=(64,4))+1j*rng.normal(size=(64,4)))[0]*np.sqrt(64/42)
 initial.append(q[:,:2]);excited.append(np.cos(.19)*q[:,:2]+np.sin(.19)*q[:,2:])
initial=np.stack(initial,axis=2);excited=np.stack(excited,axis=2)
with (tmp/'orbitals.bin').open('wb') as f:
 initial.ravel(order='F').astype('c16').tofile(f);excited.ravel(order='F').astype('c16').tofile(f)
mesh=np.array([4,1,1]);indices=np.array([[i,0,0] for i in range(4)]);R=np.where(2*indices>=mesh,indices-mesh,indices)*np.array([3.,3.5,4.])
k=indices*2*np.pi/mesh/np.array([3.,3.5,4.])+np.array([.03,.02,-.01]);na=8
z=rng.normal(size=(4,na,na));w=np.empty_like(z)
for l in range(4):w[l]=.2+.01*(z[l]+z[(-l)%4].T)
with (tmp/'hsex-w.bin').open('wb') as f:
 f.write(b'HXT1');np.array([4,4,4,4,1,1,2,2,1],dtype='i4').tofile(f)
 for v in (np.diag([3.,3.5,4.]),k,np.zeros(3),np.zeros(64),[1.],R,w):np.asarray(v,dtype='f8').tofile(f)
# w is symmetric under pair transpose and R reversal; flatten each matrix in F order.
with (tmp/'hsex-w.bin').open('r+b') as f:
 f.seek(40+8*(9+12+3+64+1+12));np.concatenate([m.ravel(order='F') for m in w]).astype('f8').tofile(f)
original=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text());link_args=[];skip=False
for item in original[1:]:
 if skip:skip=False;continue
 if item=='-o':skip=True
 elif item.endswith('/main.f90.o'):continue
 elif item.endswith('.o'):link_args.append(str((build/'src'/item).resolve()))
 else:link_args.append(item)
checker=root/'tools/check_forbidden_eigensolvers.py'
subprocess.run(['python3',str(checker),'--source-root',str(root/'src')],check=True)
subprocess.run([original[0],'-fopenmp','-fcheck=all','-I'+str(build),'-c','driver.f90'],cwd=tmp,check=True)
for tree in (build,tmp):subprocess.run(['python3',str(checker),'--object-tree',str(tree)],check=True)
subprocess.run([original[0],*link_args,str(tmp/'driver.o'),'-o',str(tmp/'test')],cwd=tmp,check=True)
env=dict(os.environ,OMP_NUM_THREADS='12',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
records=[]
for ranks in args.ranks:
 subprocess.run(['mpiexec','-n',str(ranks),str(tmp/'test')],cwd=tmp,env=env,check=True)
 data=np.fromfile(tmp/'energy-current.bin','f8');j=data[:3];e0=data[3];e=data[4:].reshape((3,3,2),order='F')
 fd=(e[:,:,1]-e[:,:,0])/(42*2*np.array([1e-3,1e-4,1e-5]))
 error=np.max(abs(fd-j[:,None]),axis=0)
 record=dict(ranks=ranks,omp=12,J_number=j.tolist(),E_exchange=e0,deltas=[1e-3,1e-4,1e-5],energy_derivative=fd.tolist(),max_error=error.tolist(),relative_error=float(error[-1]/max(abs(j))),charge_current_convention='J_charge=-J_number; J_charge=-dE/dA/volume; volume=42')
 print(json.dumps(record),flush=True)
 assert error[-1]<1e-9 and error[-1]/max(abs(j))<1e-6
 records.append(record)
print('HSEX_EXCHANGE_ENERGY_CURRENT_PASS')

if args.report:args.report.write_text(json.dumps(records,indent=2)+'\n')
