program probe
 use mpi
 use exx_spatial_local, only: s_exx_spatial_local,spatial_local_init,spatial_local_destroy, &
   s_exx_sr,sr_prepare,sr_apply,sr_destroy,sr_mask_kernel,spatial_local_apply
 implicit none
 type(s_exx_spatial_local) :: plan,full
 type(s_exx_sr) :: sr
 integer :: n(3),dims(3),coords(3),m(3),lo(3),groups(3),rank,peers,ierr,provided,a,b,color,status
 logical :: used
 complex(8),allocatable :: original_kernel(:),compact_source(:),near_action(:,:)
 integer :: x,y,z,g,j,p(3),q(3),d(3),idx,ng,trial
 real(8) :: error,all_error,h(3),tolerance
 real(8),allocatable :: multiplier(:),all_multiplier(:)
 complex(8),allocatable :: source(:),target(:,:),action(:,:),all_source(:),all_target(:),reference(:)
 character(64) :: arg
 call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,peers,ierr)
 call get_command_argument(1,arg);read(arg,*)dims
 if(product(dims)/=peers)error stop 'ranks'
 n=[24,8,8];h=1d0;ng=product(n);m=n/dims
 coords=[mod(rank,dims(1)),mod(rank/dims(1),dims(2)),rank/(dims(1)*dims(2))];lo=coords*m
 do a=1,3
  color=0
  do b=1,3
   if(b/=a)color=color*dims(b)+coords(b)
  enddo
  call MPI_Comm_split(MPI_COMM_WORLD,color,coords(a),groups(a),ierr)
 enddo
 allocate(multiplier(product(m)),all_multiplier(ng),source(product(m)),target(product(m),1), &
   action(product(m),1),all_source(ng),all_target(ng),reference(product(m)))
 do g=1,ng
  all_multiplier(g)=2d0+sin(.17d0*g)
  all_source(g)=cmplx(sin(.13d0*g),cos(.11d0*g),8)
  all_target(g)=cmplx(cos(.07d0*g),sin(.23d0*g),8)
 enddo
 g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
  g=g+1;p=lo+[x,y,z];idx=1+p(1)+n(1)*(p(2)+n(2)*p(3))
  multiplier(g)=all_multiplier(idx);source(g)=all_source(idx);target(g,1)=all_target(idx)
 enddo;enddo;enddo
 call spatial_local_init(full,n,[1,1,1],[0,0,0],[MPI_COMM_SELF,MPI_COMM_SELF,MPI_COMM_SELF], &
   all_multiplier,status)
 if(status/=0)error stop 'full kernel'
 call spatial_local_init(plan,n,dims,coords,groups,multiplier,status)
 if(status/=0)error stop 'distributed kernel'
 original_kernel=plan%kernel
 allocate(compact_source(size(source)),near_action(size(source),1))
 do trial=1,3
  plan%kernel=original_kernel
  tolerance=.1d0
  if(trial==2)tolerance=1d-3
  if(trial==3)tolerance=1d-30
  call sr_prepare(sr,plan,h,1d0,tolerance,MPI_COMM_WORLD,status)
  if(status/=0)error stop 'sr prepare'
  if(sr%ready)then
   call sr_apply(sr,MPI_COMM_WORLD,source,target,action)
   reference=0d0;g=0
   do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
    g=g+1;p=lo+[x,y,z]
    do j=1,ng
     q=[mod(j-1,n(1)),mod((j-1)/n(1),n(2)),(j-1)/(n(1)*n(2))]
     d=modulo(p-q+n/2,n)-n/2
     if(sum((d*h)**2)>sr%radius**2)cycle
     d=modulo(d,n);idx=1+d(1)+n(1)*(d(2)+n(2)*d(3))
     reference(g)=reference(g)-source(g)*full%kernel(idx)*conjg(all_source(j))*all_target(j)
    enddo
   enddo;enddo;enddo
   error=maxval(abs(reference-action(:,1)))
   call MPI_Allreduce(error,all_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
   if(all_error>1d-10)error stop 'SR direct convolution mismatch'
   if(rank==0)print *,'PASS SR',dims,trial,all_error,sr%radius,sr%length
   ! Evaluate identical cut kernels with broad and then compact, wrapped support.
   call sr_mask_kernel(sr,plan,h)
   compact_source=source;g=0
   do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
    g=g+1
    if(abs(modulo(lo(1)+x+1,n(1))-1)>1)compact_source(g)=0d0
   enddo;enddo;enddo
   call sr_apply(sr,MPI_COMM_WORLD,compact_source,target,near_action)
   call spatial_local_apply(plan,MPI_COMM_WORLD,compact_source,target,action,used,status)
   if(status/=0.or..not.used)error stop 'WF route not selected'
   error=maxval(abs(action-near_action))
   call MPI_Allreduce(error,all_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
   if(all_error>1d-10)error stop 'WF/near route mismatch'
   call spatial_local_apply(plan,MPI_COMM_WORLD,compact_source,target,action,used,status,fft_cost_limit=1d0)
   if(status/=0.or.used)error stop 'cost gate did not defer to neighborhood'
   if(rank==0)print *,'PASS SR route equivalence',dims,trial,all_error
  else
   if(rank==0)print *,'PASS SR full fallback',dims,trial
  endif
  call sr_destroy(sr)
 enddo
 call spatial_local_destroy(plan);call spatial_local_destroy(full)
 call MPI_Finalize(ierr)
end program
