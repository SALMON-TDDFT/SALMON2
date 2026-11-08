"""Exercise the actual HSEX stage/projected-current path on an even k mesh.

Requires an already configured/built MPI SALMON tree. Application objects and
this driver are audited before linking. The synthetic orbitals deliberately
have k-dependent amplitudes so the Nyquist translation is nonzero; this is an
operator identity test, not a normalized electronic-state calculation.

Requires the GNU/Makefiles build layout (link.txt and module files in build).
Example: python3 test_lattice_velocity.py --build /tmp/salmon-build --ranks 2
For isolated HSEX throughput, add --benchmark 20 --threads 1,4 --report PATH.
The benchmark runs three measured repeats after warmup, with 20 stage updates
and four Hamiltonian actions per update. It records stage/action wall times,
checksums, and the complete final action vector beside the JSON report.
Compare runs sequentially with OPENBLAS_NUM_THREADS=1 and no other workloads.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DRIVER = '''program lattice_velocity
use mpi
use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar
use salmon_global, only: num_kgrid,cval
use kli_hsex, only: hsex_init,hsex_stage,hsex_current_decomposed
implicit none
type(s_dft_system)::sys
type(s_rgrid)::mg
type(s_parallel_info)::info
type(s_orbital)::orb
type(s_scalar)::xc
real(8)::curr(3,1,1,4),pi
integer::i,j,rank,nproc,ierr,provided
call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
if(provided<MPI_THREAD_FUNNELED)error stop 'MPI funneled threads unavailable'
call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
if(nproc>4)error stop 'test ranks must not exceed four'
sys%no=1;sys%nspin=1;sys%nk=4;sys%nion=1;sys%hvol=1d0
sys%primitive_a=0d0
do i=1,3
  sys%primitive_a(i,i)=1d0
enddo
allocate(sys%vec_k(3,4),sys%Rion(3,1))
sys%vec_k=0d0;sys%Rion=0d0
pi=acos(-1d0)
do i=1,4
  sys%vec_k(1,i)=(i-1)*pi/2
enddo
mg%is=1;mg%ie=1
info%ik_s=1+rank*4/nproc;info%ik_e=(rank+1)*4/nproc
info%numk=info%ik_e-info%ik_s+1
info%isize_ko=nproc;info%id_ko=rank;info%icomm_ko=MPI_COMM_WORLD
allocate(orb%zwf(1,1,1,1,1,info%ik_s:info%ik_e,1))
do i=info%ik_s,info%ik_e
  orb%zwf(1,1,1,1,1,i,1)=sqrt(1d0+.2d0*(-1d0)**(i-1))
enddo
allocate(xc%f(1,1,1));xc%f=0d0
num_kgrid=[4,1,1];cval=1d0
call hsex_init(sys,mg,info,orb,0d0,.true.,xc)
call hsex_stage(sys,mg,info,orb,[.13d0,0d0,0d0],2)
curr=0d0
call hsex_current_decomposed(sys,mg,info,orb,curr)
if(rank==0)then
  open(11,file='velocity.bin',form='unformatted',access='stream')
  write(11)curr
  close(11)
endif
call MPI_Finalize(ierr)
end program lattice_velocity
'''


BENCHMARK = '''program hsex_benchmark
use mpi
use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar
use salmon_global, only: num_kgrid,cval
use kli_hsex, only: hsex_init,hsex_stage,hsex_add
implicit none
type(s_dft_system)::sys
type(s_rgrid)::mg
type(s_parallel_info)::info
type(s_orbital)::orb,hx
type(s_scalar)::xc
integer::i,j,g,x,y,z,rank,nproc,ierr,step,repeat,iterations,action,provided
real(8)::pi,elapsed,start,mark,tstage,tadd,times(3),maxtimes(3),cs(3),totalcs(3)
character(20)::argument
call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
if(provided<MPI_THREAD_FUNNELED)error stop 'MPI funneled threads unavailable'
call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
call get_command_argument(1,argument);read(argument,*)iterations
if(nproc>64)error stop 'test ranks must not exceed 64'
sys%no=4;sys%nspin=1;sys%nk=64;sys%nion=1;sys%hvol=1d0/1728
sys%primitive_a=0d0
do i=1,3
  sys%primitive_a(i,i)=1d0
enddo
allocate(sys%vec_k(3,64),sys%Rion(3,1))
sys%Rion=0d0;pi=acos(-1d0)
do i=1,64
  sys%vec_k(:,i)=[mod(i-1,4),mod((i-1)/4,4),(i-1)/16]*pi/2
enddo
mg%is=1;mg%ie=12
info%ik_s=1+rank*64/nproc;info%ik_e=(rank+1)*64/nproc
info%numk=info%ik_e-info%ik_s+1
info%isize_ko=nproc;info%id_ko=rank;info%icomm_ko=MPI_COMM_WORLD
allocate(orb%zwf(12,12,12,1,4,info%ik_s:info%ik_e,1))
allocate(hx%zwf(12,12,12,1,4,info%ik_s:info%ik_e,1))
do i=info%ik_s,info%ik_e
  do j=1,4
    do z=1,12
      do y=1,12
        do x=1,12
          g=x+12*(y-1)+144*(z-1)
          orb%zwf(x,y,z,1,j,i,1)=cmplx(sin(.03d0*g*j+.02d0*i),cos(.017d0*g*j-.03d0*i),8)
        enddo
      enddo
    enddo
  enddo
enddo
allocate(xc%f(12,12,12));xc%f=0d0
num_kgrid=[4,4,4];cval=1d0
call hsex_init(sys,mg,info,orb,0d0,.true.,xc)
! Warm up code paths and allocations outside measurements.
call hsex_stage(sys,mg,info,orb,[.13d0,0d0,0d0],2)
hx%zwf=0d0
call hsex_add(sys,mg,info,orb,hx)
if(rank==0)open(11,file='benchmark.bin',form='unformatted',access='stream')
do repeat=1,3
  tstage=0d0;tadd=0d0
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  start=MPI_Wtime()
  do step=1,iterations
    mark=MPI_Wtime()
    call hsex_stage(sys,mg,info,orb,[.13d0+.001d0*step,0d0,0d0],2)
    tstage=tstage+MPI_Wtime()-mark
    mark=MPI_Wtime()
    hx%zwf=0d0
    do action=1,4
      call hsex_add(sys,mg,info,orb,hx)
    enddo
    tadd=tadd+MPI_Wtime()-mark
  enddo
  elapsed=MPI_Wtime()-start
  times=[elapsed,tstage,tadd]
  call MPI_Allreduce(times,maxtimes,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  cs=[sum(real(hx%zwf,8)),sum(aimag(hx%zwf)),sum(abs(hx%zwf)**2)]
  call MPI_Allreduce(cs,totalcs,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(rank==0)write(11)maxtimes,totalcs
enddo
if(rank==0)close(11)
write(argument,'("action-",i0,".bin")')rank
open(11,file=trim(argument),form='unformatted',access='stream')
write(11)hx%zwf
close(11)
call MPI_Finalize(ierr)
end program hsex_benchmark
'''


def run(build, ranks, benchmark=0, threads=(1,), report=None):
    build = build.resolve()
    original = shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text())
    compiler = original[0]
    # Reuse exact application libraries and compiler flags; replace main only.
    args = []
    skip = False
    for item in original[1:]:
        if skip:
            skip = False
            continue
        if item == '-o':
            skip = True
        elif item.endswith('/main.f90.o'):
            continue
        elif item.endswith('.o'):
            args.append(str((build/'src'/item).resolve()))
        else:
            args.append(item)
    checker = ROOT/'tools/check_forbidden_eigensolvers.py'
    env = os.environ.copy()
    env.setdefault('OMP_NUM_THREADS', '1')
    env.setdefault('OPENBLAS_NUM_THREADS', '1')
    with tempfile.TemporaryDirectory(prefix='hsex-lattice-') as tmp:
        tmp = Path(tmp)
        (tmp/'driver.f90').write_text(BENCHMARK if benchmark else DRIVER)
        subprocess.run(['python3',str(checker),'--source-root',str(ROOT/'src')],check=True)
        subprocess.run([compiler,'-fopenmp','-fcheck=all','-I'+str(build),'-c','driver.f90'],cwd=tmp,check=True)
        for tree in (build,tmp):
            subprocess.run(['python3',str(checker),'--object-tree',str(tree)],check=True)
        subprocess.run([compiler,*args,str(tmp/'driver.o'),'-o',str(tmp/'test')],cwd=tmp,check=True)
        nk = 64 if benchmark else 4
        grid_size, boxes, bands = (12,6,4) if benchmark else (1,1,1)
        mesh = np.array([4,4,4] if benchmark else [4,1,1])
        index = np.indices(tuple(mesh[::-1])).reshape(3,-1)[::-1].T
        k = index*2*np.pi/mesh
        R = np.where(2*index >= mesh,index-mesh,index)
        with (tmp/'hsex-w.bin').open('wb') as f:
            f.write(b'HXT1')
            np.array([*[grid_size]*3,*mesh,boxes,bands,1],dtype='=i4').tofile(f)
            for value in (np.eye(3), k, np.zeros(3), np.zeros(grid_size**3), [1.], R):
                np.asarray(value,dtype='=f8').tofile(f)
            for r in R:
                np.full(boxes**6,.01/(1+np.dot(r,r)) if benchmark else 1.,dtype='=f8').tofile(f)
        if benchmark:
            digest = hashlib.sha256((tmp/'test').read_bytes()).hexdigest()
            reports = []
            for count in threads:
                assert hashlib.sha256((tmp/'test').read_bytes()).hexdigest() == digest
                env['OMP_NUM_THREADS'] = str(count)
                subprocess.run(['mpiexec','-n',str(ranks),str(tmp/'test'),str(benchmark)],
                               cwd=tmp,env=env,check=True)
                rows = np.fromfile(tmp/'benchmark.bin',dtype='=f8').reshape(3,6)
                record = dict(build=str(build),ranks=ranks,threads=count,iterations=benchmark,
                              binary_sha256=digest,seconds=rows[:,:3].tolist(),checksums=rows[:,3:].tolist())
                if report is not None:
                    action = np.concatenate([np.fromfile(tmp/f'action-{rank}.bin',dtype='=c16')
                                             for rank in range(ranks)])
                    action_path = report.with_suffix(f'.omp{count}.npy')
                    np.save(action_path,action)
                    record['action_file'] = str(action_path)
                reports.append(record)
                print(json.dumps(record),flush=True)
            if report is not None:
                report.write_text(json.dumps(reports,indent=2)+'\n')
            return
        subprocess.run(['mpiexec','-n',str(ranks),str(tmp/'test')],cwd=tmp,env=env,check=True)
        actual = np.fromfile(tmp/'velocity.bin',dtype='=f8').reshape(4,3)
        occupations = 1+.2*(-1.)**np.arange(4)
        A = np.array([.13,0,0])
        def density(a):
            return np.exp(-1j*R@(k+a).T)@occupations/4
        # Freeze the ORIGINAL translation kernel, before Hermitian projection.
        kernel = -(density(A)-density(np.zeros(3)))
        def energies(a):
            return occupations*(np.exp(1j*(k+a)@R.T)@kernel).real
        expected = np.column_stack([(energies(A+1e-6*e)-energies(A-1e-6*e))/2e-6 for e in np.eye(3)])
        assert np.max(abs(expected[:,0])) > .05
        np.testing.assert_allclose(actual,expected,rtol=1e-8,atol=2e-10)
        print(f'actual HSEX lattice velocity finite difference PASS: {ranks} MPI; max error {np.max(abs(actual-expected)):.3e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--ranks',type=int,default=1)
    parser.add_argument('--benchmark',type=int,default=0,help='iterations per repeat, three repeats')
    parser.add_argument('--threads',default='1',help='comma-separated OMP counts for benchmark')
    parser.add_argument('--report',type=Path)
    options = parser.parse_args()
    run(options.build,options.ranks,options.benchmark,tuple(map(int,options.threads.split(','))),options.report)
