! Focused primitive FCC C3 spinor-crystal-orbit action versus full-k exchange.
! Sources contain sqrt(f); independent targets deliberately differ from sources.
program primitive_orbits_probe
 use mpi
 use exx_symmetry
 use exx_k_exchange
 use exx_spinor
 use exx_k_orbits,only:apply_spinor_k_orbits
 implicit none
 type(s_exx_symmetry_map)::map
 type(exx_k_kernel)::op,refop
 integer,parameter::ng=64,nrep=4,nfull=8
 real(8)::ops(3,4,3),cell(3,3),reciprocal(3,3),ks(3,nrep),weights(nrep),h(3),pi,origin(3),errmax,globalerr
 complex(8)::source(2*ng,2,nrep),target(2*ng,3,nrep),ref(2*ng,3,nfull)
 complex(8),allocatable::es(:,:,:),et(:,:,:),action(:,:,:)
 integer::rank,np,e,i,j,k,s,p,first,nlocal,repfirst,nowned,q
 integer::starts(2),counts(2)
 call MPI_Init(e);call MPI_Comm_rank(MPI_COMM_WORLD,rank,e);call MPI_Comm_size(MPI_COMM_WORLD,np,e)
 if(np/=2)call MPI_Abort(MPI_COMM_WORLD,1,e)
 pi=acos(-1d0);cell=2d0
 do i=1,3
  cell(i,i)=0d0
 enddo
 reciprocal=pi/2
 do i=1,3
  reciprocal(i,i)=-pi/2
 enddo
 h=sqrt(8d0)/4
 ks(:,1)=matmul(reciprocal,[-.25d0,-.25d0,-.25d0])
 ks(:,2)=matmul(reciprocal,[.25d0,-.25d0,-.25d0])
 ks(:,3)=matmul(reciprocal,[.25d0,.25d0,-.25d0])
 ks(:,4)=matmul(reciprocal,[.25d0,.25d0,.25d0])
 weights=[1d0,3d0,3d0,1d0]/8;ops=0d0
 do i=1,3
  ops(i,i,1)=1d0
 enddo
 ops(1,2,2)=1;ops(2,3,2)=1;ops(3,1,2)=1
 ops(:,1:3,3)=matmul(ops(:,1:3,2),ops(:,1:3,2))
 origin=[.25d0,.5d0,0d0]
 do s=1,3
  ops(:,4,s)=origin-matmul(ops(:,1:3,s),origin)
 enddo
 call symmetry_init(map,[4,4,4],h,ks,weights,ops,ops,[2,2,2],e,cell=cell)
 if(e/=0)call MPI_Abort(MPI_COMM_WORLD,2,e)
 if(maxval(abs(map%cell-cell))>1d-13)call MPI_Abort(MPI_COMM_WORLD,3,e)
 do k=1,nrep;do j=1,2;do i=1,2*ng
  source(i,j,k)=sqrt(.2d0+.1d0*j+.05d0*k)*cmplx(sin(real(i+3*j+2*k,8)), &
   cos(real(2*i-j+k,8)),8)/7d0
 enddo;enddo;enddo
 do k=1,nrep;do j=1,3;do i=1,2*ng
  target(i,j,k)=cmplx(cos(real(i+2*j+k,8)),sin(real(3*i-j+k,8)),8)/8d0
 enddo;enddo;enddo
 call symmetry_expand(map,1,source,target,es,et,e)
 if(e/=0)call MPI_Abort(MPI_COMM_WORLD,4,e)
 call exx_k_kernel_init(refop,[4,4,4],[2,2,2],h,map%full_k,.11d0,2,e,cell=cell,reciprocal=reciprocal)
 if(e/=0)call MPI_Abort(MPI_COMM_WORLD,5,e)
 call exx_spinor_apply_distributed(refop,es,et,ref,[1],[nfull],0,MPI_COMM_SELF,transpose_local,e)
 if(e/=0)call MPI_Abort(MPI_COMM_WORLD,6,e)
 starts=[1,map%first(3)];counts=[map%first(3)-1,nfull-map%first(3)+1]
 repfirst=1+2*rank;nowned=2;first=starts(rank+1);nlocal=counts(rank+1)
 allocate(action(2*ng,3,nowned))
 call exx_k_kernel_init(op,[4,4,4],[2,2,2],h,map%full_k,.11d0,2,e, &
  phase_start=first,phase_count=nlocal,cell=cell,reciprocal=reciprocal)
 if(e/=0)call MPI_Abort(MPI_COMM_WORLD,7,e)
 call apply_spinor_k_orbits(op,map,target(:,:,repfirst:repfirst+nowned-1),action,starts,counts, &
  repfirst,rank,MPI_COMM_WORLD,transpose_tiles,fill,e)
 if(e/=0)call MPI_Abort(MPI_COMM_WORLD,8,e)
 errmax=0d0
 do j=1,nowned
  q=map%first(repfirst+j-1);errmax=max(errmax,maxval(abs(action(:,:,j)-ref(:,:,q))))
 enddo
 call MPI_Allreduce(errmax,globalerr,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,e)
 if(globalerr>2d-11)then
  print *,'PRIMITIVE_C3_SPINOR_ORBITS_FAIL',rank,errmax
  call MPI_Abort(MPI_COMM_WORLD,9,e)
 endif
 if(rank==0)print *,'PRIMITIVE_C3_SPINOR_ORBITS_PASS',globalerr
 call MPI_Finalize(e)
contains
 subroutine fill(j,points,density,spin_left,spin_right,columns)
  integer,intent(in)::j,points(:),spin_left,spin_right,columns(:)
  complex(8),intent(out)::density(:,:)
  complex(8)::left(size(points),size(es,2)),right(size(columns),size(es,2))
  integer::q,b
  q=first+j-1
  left=es(points+(spin_left-1)*ng,:,q);right=es(columns+(spin_right-1)*ng,:,q)
  do b=1,size(es,2)
   left(:,b)=left(:,b)*op%phase(points,j);right(:,b)=right(:,b)*op%phase(columns,j)
  enddo
  density=matmul(left,conjg(transpose(right)))
 end subroutine
 subroutine transpose_tiles(send,recv,count,comm)
  complex(8),intent(in)::send(:)
  complex(8),intent(out)::recv(:)
  integer,intent(in)::count,comm
  integer::e
  call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,comm,e)
 end subroutine
 subroutine transpose_local(send,recv,count,comm)
  complex(8),intent(in)::send(:)
  complex(8),intent(out)::recv(:)
  integer,intent(in)::count,comm
  recv=send
 end subroutine
end program
