program probe
 use mpi
 use hse_grid_exchange
 use hse_wannier
 implicit none
 type(s_hse_grid_exchange) :: op
 type(s_hse_wannier) :: ref
 integer :: rank,np,ierr,status,g,j,k,nlocal,ng,case_id
 integer,allocatable :: idx(:),badidx(:)
 complex(8),allocatable :: source(:,:),target(:,:),action(:,:),full_source(:,:),full_target(:,:,:),expected(:,:,:)
 complex(8) :: inner(2),total(2)
 real(8) :: h(3),err,globalerr,phase
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 ng=60;h=[.7d0,.9d0,1.1d0];nlocal=count([(mod(g-1,np)==rank,g=1,ng)])
 allocate(idx(nlocal),badidx(nlocal),source(nlocal,3),target(nlocal,2),action(nlocal,2))
 allocate(full_source(ng,3),full_target(ng,2,1),expected(ng,2,1))
 j=0
 do g=1,ng
  if(mod(g-1,np)==rank)then
   j=j+1;idx(j)=g
  endif
  do k=1,3
   phase=.17d0*g*k
   full_source(g,k)=cmplx(cos(phase),sin(phase*.73d0),8)/sqrt(real(ng,8))
  enddo
  ! Arbitrary oscillatory targets are not constrained to the occupied space.
  full_target(g,1,1)=cmplx(cos(2.8d0*g),sin(2.3d0*g),8)
  full_target(g,2,1)=cmplx(sin(1.9d0*g),cos(2.7d0*g),8)
 enddo
 target=full_target(idx,:,1)
 ! Reject a partition with a duplicated row and consequently a missing row.
 badidx=idx
 if(rank==0.and.size(idx)>1)badidx(1)=badidx(2)
 call grid_exchange_init(op,[5,4,3],h,.11d0,badidx,MPI_COMM_WORLD,2,.false.,status)
 if(status==0)error stop 'invalid partition accepted'
 call grid_exchange_init(op,[5,4,3],h,.11d0,idx,MPI_COMM_WORLD,2,.false.,status)
 if(status/=0)error stop 'grid init failed'
 call wannier_init(ref,[5,4,3],[1,1,1],h,reshape([0d0,0d0,0d0],[3,1]),.11d0,status)
 if(status/=0)error stop 'reference init failed'
 do case_id=1,2
  if(case_id==2)then
   ! Source support restriction only; targets remain unmasked.
   full_source(1:20,:)=0d0
  endif
  source=full_source(idx,:)
  call grid_exchange_set_source(op,source,status)
  if(status/=0)error stop 'source redistribution failed'
  if(size(op%kernel%source,2)/=count([(mod(k-1,np)==rank,k=1,3)]))error stop 'source replication'
  ref%source=full_source
  call wannier_apply(ref,full_target,expected,status)
  if(status/=0)error stop 'reference action failed'
  call grid_exchange_apply(op,target,action,status)
  if(status/=0)error stop 'grid action failed'
  err=maxval(abs(action-expected(idx,:,1)))
  call MPI_Allreduce(err,globalerr,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if(globalerr>2d-12)error stop 'grid/reference parity failed'
  inner(1)=sum(conjg(target(:,1))*action(:,2))*product(h)
  inner(2)=sum(conjg(action(:,1))*target(:,2))*product(h)
  call MPI_Allreduce(inner,total,2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(abs(total(1)-total(2))>2d-11)error stop 'Hermiticity failed'
  if(rank==0)write(*,*)'grid exchange parity/Hermiticity passed',case_id,globalerr,abs(total(1)-total(2))
 enddo
 ! One source exercises source-idle ranks in multi-rank tests.
 call grid_exchange_set_source(op,source(:,1:1),status)
 if(status/=0)error stop 'single source failed'
 ref%source=full_source(:,1:1)
 call wannier_apply(ref,full_target,expected,status)
 call grid_exchange_apply(op,target,action,status)
 if(status/=0.or.maxval(abs(action-expected(idx,:,1)))>2d-12)error stop 'idle source owner failed'
 call grid_exchange_free(op);call wannier_destroy(ref)
 call MPI_Finalize(ierr)
end program
