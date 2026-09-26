program distributed_build_probe
 use mpi
 use lcfo_dist_rows
 use lcfo_dist_dense
 use hse_ace,only:hse_ace_state,hse_ace_build,hse_ace_apply
 implicit none
 type(s_lcfo_halo) :: halo
 type(hse_ace_state) :: ace,reference
 integer :: rank,np,ierr,nlocal,lo,j,k,p,n,ng,trial,status
 integer,allocatable :: counts(:),offsets(:),selected(:)
 complex(8),allocatable :: mixing(:,:), c(:,:),near(:,:),contribution(:,:),result(:,:),h(:,:),full(:,:),expected(:,:), &
  block(:,:),w(:,:),allw(:,:),rotation(:,:),previous(:,:),target(:,:,:),action(:,:,:),ref_action(:,:,:)
 real(8) :: error,minimum,pi
 call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 if(np/=2.and.np/=4)error stop 'probe requires MPI2 or MPI4'
 allocate(counts(np),offsets(np+1));counts=2;counts(1)=1;counts(np)=counts(np)+1
 if(np==2)counts=[3,5]
 offsets(1)=0
 do j=1,np;offsets(j+1)=offsets(j)+counts(j);enddo
 nlocal=counts(rank+1);lo=offsets(rank+1)
 if(rank==0)then
  selected=[7,1,4,2]
 else
  selected=[5,2,8,4,3]
 endif
 call lcfo_halo_init(halo,counts,selected,MPI_COMM_WORLD)
 allocate(full(8,3),h(8,8),block(size(selected),size(selected)))
 do j=1,3;do k=1,8;full(k,j)=cmplx(.03d0*k*j,.02d0*(k-j),8);enddo;enddo
 c=full(lo+1:lo+nlocal,:)
 call lcfo_halo_get(halo,c,near)
 if(maxval(abs(near-full(selected,:)))>1d-14)error stop 'requested rows do not match'
 do j=1,size(selected);do k=1,size(selected)
  block(k,j)=cmplx(.01d0*(j+k+rank),.02d0*(j-k),8)
 enddo;enddo
 h=0d0
 do j=1,size(selected);do k=1,size(selected);h(selected(k),selected(j))=block(k,j);enddo;enddo
 allocate(expected(8,8));call MPI_Allreduce(h,expected,64,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
 contribution=matmul(block,near)
 call lcfo_halo_sum(halo,contribution,result)
 w=matmul(expected,full)
 if(maxval(abs(result-w(lo+1:lo+nlocal,:)))>1d-13)error stop 'fragment Hermitian action mismatch'
 call lcfo_gather_root(c,counts,MPI_COMM_WORLD,previous)
 if(rank==0)then
  if(maxval(abs(previous-full))>1d-14)error stop 'initial gather mismatch'
 else
  if(size(previous)/=0)error stop 'nonroot holds global coefficients'
 endif
 call lcfo_halo_free(halo)
 deallocate(c,full,h,block,expected,w,previous)
 pi=acos(-1d0)
 do trial=1,2
  n=3;if(trial==2)n=129
  ng=2*n+1;counts=ng/np;counts(np)=counts(np)+mod(ng,np)
  offsets(1)=0
  do j=1,np;offsets(j+1)=offsets(j)+counts(j);enddo
  nlocal=counts(rank+1);lo=offsets(rank+1)
  allocate(full(ng,n),allw(ng,n),rotation(n,n),target(ng,2,1),action(nlocal,2,1),ref_action(ng,2,1))
  do j=1,n;do k=1,ng
   full(k,j)=exp(cmplx(0d0,2*pi*(k-1)*(j-1)/ng,8))/sqrt(dble(ng))
   allw(k,j)=-full(k,j)*(1d0+.2d0*k/ng)
  enddo;enddo
  c=full(lo+1:lo+nlocal,:);w=allw(lo+1:lo+nlocal,:);previous=c
  allocate(mixing(n,n))
  do j=1,n;do k=1,n
   mixing(k,j)=exp(cmplx(0d0,2*pi*(k-1)*(j-1)/n,8))/sqrt(dble(n))
  enddo;enddo
  previous=matmul(c,mixing)
  call lcfo_distributed_polar(c,previous,MPI_COMM_WORLD,rotation,minimum,status)
  if(status/=0.or.abs(minimum-1d0)>1d-10)error stop 'polar status/singular values'
  if(maxval(abs(matmul(c,rotation)-previous))>1d-10)error stop 'distributed polar mismatch'
  call lcfo_distributed_ace_build(ace,c,w,1d0,MPI_COMM_WORLD,status)
  if(status/=0)error stop 'distributed ACE rejected SPD metric'
  if(size(ace%factors,1)/=nlocal)error stop 'ACE factors are not local'
  call hse_ace_build(reference,reshape(full,[ng,n,1]),reshape(allw,[ng,n,1]),1d0,status)
  if(status/=0)error stop 'reference ACE failed'
  target(:,:,1)=full(:,1:2)
  call hse_ace_apply(reference,target,ref_action,status)
  block=matmul(conjg(transpose(ace%factors(:,:,1))),target(lo+1:lo+nlocal,:,1))
  allocate(expected(n,2));call MPI_Allreduce(block,expected,2*n,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  action(:,:,1)=-matmul(ace%factors(:,:,1),expected)
  error=maxval(abs(action(:,:,1)-ref_action(lo+1:lo+nlocal,:,1)))
  if(error>1d-10)error stop 'distributed ACE action mismatch'
  ! Use repeated identity rows so the remaining rows are well conditioned;
  ! a contiguous half of Fourier rows would be numerically rank deficient.
  c=0d0
  do k=1,nlocal;c(k,mod(lo+k-1,n)+1)=1d0;enddo
  ! One rank contributes zero; the other rows still span the occupied space.
  w=-c;if(rank==0)w=0d0
  call lcfo_distributed_ace_build(ace,c,w,1d0,MPI_COMM_WORLD,status)
  if(status/=0)error stop 'local zero exchange rejected globally SPD metric'
  ! Globally zero exchange is also valid.
  w=0d0
  call lcfo_distributed_ace_build(ace,c,w,1d0,MPI_COMM_WORLD,status)
  if(status/=0.or.maxval(abs(ace%factors))>0d0)error stop 'zero exchange'
  w=c
  call lcfo_distributed_ace_build(ace,c,w,1d0,MPI_COMM_WORLD,status)
  if(status==0)error stop 'positive exchange metric accepted'
  w=-c;w(:,n)=0d0
  call lcfo_distributed_ace_build(ace,c,w,1d0,MPI_COMM_WORLD,status)
  if(status==0)error stop 'singular exchange metric accepted'
  if(rank==0)write(*,*)'Distributed construction passed: orbitals/error =',n,error
  deallocate(mixing,full,allw,rotation,target,action,ref_action,c,w,previous,block,expected)
 enddo
 call MPI_Finalize(ierr)
end program
