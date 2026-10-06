! Compact source convolution on native Cartesian blocks. Only compact tiles
! are replicated; the inverse global kernel remains spatially distributed.
! All calls are collective over comm_r with matching target column counts.
module exx_spatial_local
 use iso_c_binding
 use iso_fortran_env, only: int64
 use exx_batch_backend, only: s_exx_batch_backend
 use communication, only: comm_summation,comm_get_max,comm_get_groupinfo,comm_create_group,comm_free_group,comm_bcast,comm_exchange
 use fftw_blocks, only: mesh_transform,block_layout
 use exx_local_fft, only: s_exx_local_fft,exx_local_prepare_compact,exx_local_apply,exx_local_destroy, &
                          compact_axis_size,compact_kernel_bounds,smooth_size,exx_local_cpu_prepare,exx_local_cpu_pair
 implicit none
 private
 include 'fftw3.f03'
 public :: s_exx_sr,sr_prepare,sr_apply,sr_destroy,sr_mask_kernel
 type s_exx_sr
  integer :: n(3)=0,m(3)=0,lo(3)=0,length(3)=0,origin(3)=0,rank=0,peers=1
  real(8) :: radius=0d0,tail_l1=0d0,tail_l2=0d0
  logical :: ready=.false.
  integer,allocatable :: origins(:,:)
  logical,allocatable :: neighbors(:)
  complex(8),allocatable :: work(:,:,:),filter(:,:,:),send(:,:,:),recv(:,:,:)
  type(c_ptr) :: forward=c_null_ptr,backward=c_null_ptr
 end type
 public :: s_exx_spatial_local,spatial_local_init,spatial_local_apply,spatial_local_destroy,local_batch_action
 abstract interface
  subroutine local_batch_action(padded,indices,filter,source,targets,action,status)
   implicit none
   integer,intent(in) :: padded(3),indices(:)
   complex(8),intent(in) :: filter(:,:,:),source(:),targets(:,:)
   complex(8),intent(out) :: action(:,:)
   integer,intent(out) :: status
  end subroutine
 end interface
 type s_exx_spatial_local
  integer :: n(3)=0,m(3)=0,lo(3)=0
  complex(8),allocatable :: kernel(:)
  type(s_exx_local_fft) :: fft
  class(s_exx_batch_backend),allocatable :: backend
  procedure(local_batch_action),pointer,nopass :: batch_action=>null()
  integer :: batch_size=8
 end type
contains
 subroutine sr_destroy(sr)
  implicit none
  type(s_exx_sr),intent(inout) :: sr
  if(c_associated(sr%forward))call fftw_destroy_plan(sr%forward)
  if(c_associated(sr%backward))call fftw_destroy_plan(sr%backward)
  sr%forward=c_null_ptr;sr%backward=c_null_ptr;sr%ready=.false.
  sr%radius=0d0;sr%tail_l1=0d0;sr%tail_l2=0d0;sr%length=0
  if(allocated(sr%work))deallocate(sr%work,sr%filter,sr%send,sr%recv,sr%origins,sr%neighbors)
 end subroutine

 subroutine sr_prepare(sr,plan,h,omega,tolerance,comm,status)
  implicit none
  type(s_exx_sr),intent(inout) :: sr
  type(s_exx_spatial_local),intent(in) :: plan
  real(8),intent(in) :: h(3),omega,tolerance
  integer,intent(in) :: comm
  integer,intent(out) :: status
  integer :: axis,j,g,x,y,z,p(3),d(3),q(3),halo(3),bad,other_origin(3)
  real(8) :: left,right,middle,tail(2),total(2)
  complex(8),allocatable :: buffer(:,:,:)
  call sr_destroy(sr)
  status=0
  if(omega<=0d0.or.tolerance<=0d0.or.tolerance>=1d0)then
   status=1;return
  endif
  left=0d0;right=1d0
  do while(erfc(right)>tolerance)
   right=2d0*right
  enddo
  do j=1,80
   middle=(left+right)/2d0
   if(erfc(middle)>tolerance)then
    left=middle
   else
    right=middle
   endif
  enddo
  sr%radius=right/omega
  sr%n=plan%n;sr%m=plan%m;sr%lo=plan%lo
  ! Clip before conversion to integer, including very small positive omega.
  halo=ceiling(min(sr%radius/h,real(sr%n,8)/2d0))
  do axis=1,3
   sr%length(axis)=min(sr%n(axis),smooth_size(min(sr%n(axis),sr%m(axis)+2*halo(axis))))
  enddo
  call comm_get_groupinfo(comm,sr%rank,sr%peers)
  if(all(sr%length==sr%n))return
  sr%origin=sr%lo-halo
  where(sr%length==sr%n)sr%origin=0
  allocate(sr%work(sr%length(1),sr%length(2),sr%length(3)),source=(0d0,0d0))
  allocate(sr%filter(sr%length(1),sr%length(2),sr%length(3)),buffer(sr%length(1),sr%length(2),sr%length(3)))
  allocate(sr%send(sr%m(1),sr%m(2),sr%m(3)),sr%recv(sr%m(1),sr%m(2),sr%m(3)))
  allocate(sr%origins(3,sr%peers),sr%neighbors(sr%peers))
  do j=0,sr%peers-1
   p=sr%lo
   call comm_bcast(p,comm,j)
   sr%origins(:,j+1)=p
  enddo
  do j=0,sr%peers-1
   other_origin=sr%origins(:,j+1)-halo
   where(sr%length==sr%n)other_origin=0
   sr%neighbors(j+1)=boxes_intersect(sr%origins(:,j+1),sr%origin).or.boxes_intersect(sr%lo,other_origin)
  enddo
  tail=0d0;g=0
  do z=0,sr%m(3)-1;do y=0,sr%m(2)-1;do x=0,sr%m(1)-1
   g=g+1;p=sr%lo+[x,y,z];d=modulo(p+sr%n/2,sr%n)-sr%n/2
   if(sum((d*h)**2)>sr%radius**2)then
    tail(1)=tail(1)+abs(plan%kernel(g));tail(2)=tail(2)+abs(plan%kernel(g))**2
   else
    q=modulo(d,sr%length)+1
    sr%work(q(1),q(2),q(3))=plan%kernel(g)
   endif
  enddo;enddo;enddo
  call comm_summation(tail,total,2,comm)
  sr%tail_l1=total(1);sr%tail_l2=sqrt(total(2))
  call comm_summation(sr%work,buffer,size(buffer),comm)
  sr%work=buffer
  sr%forward=fftw_plan_dft_3d(sr%length(3),sr%length(2),sr%length(1),sr%work,sr%work, &
    FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
  sr%backward=fftw_plan_dft_3d(sr%length(3),sr%length(2),sr%length(1),sr%work,sr%work, &
    FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
  bad=0
  if(.not.c_associated(sr%forward).or..not.c_associated(sr%backward))bad=1
  call comm_get_max(bad,comm)
  if(bad/=0)then
   call sr_destroy(sr);status=1;return
  endif
  call fftw_execute_dft(sr%forward,sr%work,sr%work)
  sr%filter=sr%work;sr%ready=.true.
 contains
  logical function boxes_intersect(block,origin) result(found)
   implicit none
   integer,intent(in) :: block(3),origin(3)
   integer :: a,k
   logical :: overlap
   found=.true.
   do a=1,3
    overlap=.false.
    do k=0,sr%m(a)-1
     if(modulo(block(a)+k-origin(a),sr%n(a))<sr%length(a))overlap=.true.
    enddo
    found=found.and.overlap
   enddo
  end function
 end subroutine

 ! Both FFT routes must see the same operator; selection changes only its evaluation.
 subroutine sr_mask_kernel(sr,plan,h)
  implicit none
  type(s_exx_sr),intent(in) :: sr
  type(s_exx_spatial_local),intent(inout) :: plan
  real(8),intent(in) :: h(3)
  integer :: x,y,z,g,d(3)
  if(.not.sr%ready)return
  g=0
  do z=0,plan%m(3)-1;do y=0,plan%m(2)-1;do x=0,plan%m(1)-1
   g=g+1;d=modulo(plan%lo+[x,y,z]+plan%n/2,plan%n)-plan%n/2
   if(sum((d*h)**2)>sr%radius**2)plan%kernel(g)=0d0
  enddo;enddo;enddo
 end subroutine

 subroutine sr_apply(sr,comm,source,targets,action)
  implicit none
  type(s_exx_sr),intent(inout) :: sr
  integer,intent(in) :: comm
  complex(8),intent(in) :: source(:),targets(:,:)
  complex(8),intent(out) :: action(:,:)
  integer :: j,round,partner,rounds,x,y,z,g,q(3)
  real(8) :: scale
  action=0d0
  rounds=1
  do while(rounds<sr%peers)
   rounds=2*rounds
  enddo
  scale=1d0/real(product(int(sr%length,int64)),8)
  do j=1,size(targets,2)
   sr%send=reshape(conjg(source)*targets(:,j),sr%m)
   sr%work=0d0
   call insert_block(sr%send,sr%lo)
   ! XOR rounds give symmetric Sendrecv partners without MPI calls in OMP.
   do round=1,rounds-1
    partner=ieor(sr%rank,round)
    if(partner>=sr%peers)cycle
    if(.not.sr%neighbors(partner+1))cycle
    call comm_exchange(sr%send,partner,sr%recv,partner,0,comm)
    call insert_block(sr%recv,sr%origins(:,partner+1))
   enddo
   call fftw_execute_dft(sr%forward,sr%work,sr%work)
   sr%work=sr%work*sr%filter
   call fftw_execute_dft(sr%backward,sr%work,sr%work)
   g=0
   do z=0,sr%m(3)-1;do y=0,sr%m(2)-1;do x=0,sr%m(1)-1
    g=g+1;q=modulo(sr%lo+[x,y,z]-sr%origin,sr%n)+1
    action(g,j)=-source(g)*sr%work(q(1),q(2),q(3))*scale
   enddo;enddo;enddo
  enddo
 contains
  subroutine insert_block(values,lo)
   implicit none
   complex(8),intent(in) :: values(:,:,:)
   integer,intent(in) :: lo(3)
   integer :: a,b,c,p(3)
   do c=0,sr%m(3)-1;do b=0,sr%m(2)-1;do a=0,sr%m(1)-1
    p=modulo(lo+[a,b,c]-sr%origin,sr%n)
    if(any(p>=sr%length))cycle
    sr%work(p(1)+1,p(2)+1,p(3)+1)=values(a+1,b+1,c+1)
   enddo;enddo;enddo
  end subroutine
 end subroutine
 subroutine spatial_local_destroy(plan,status)
  implicit none
  type(s_exx_spatial_local),intent(inout) :: plan
  integer,intent(out),optional :: status
  integer :: backend_status
  backend_status=0
  if(allocated(plan%backend))then
   call plan%backend%release(backend_status)
   deallocate(plan%backend)
  endif
  if(present(status))status=backend_status
  call exx_local_destroy(plan%fft)
  if(allocated(plan%kernel))deallocate(plan%kernel)
  plan%n=0;plan%m=0;plan%lo=0
  nullify(plan%batch_action);plan%batch_size=8
 end subroutine
 subroutine spatial_local_init(plan,n,dims,coords,comm,multiplier,status)
  implicit none
  type(s_exx_spatial_local),intent(inout) :: plan
  integer,intent(in) :: n(3),dims(:),coords(:),comm(:)
  real(8),intent(in) :: multiplier(:) ! Cartesian xyz block order
  integer,intent(out) :: status
  integer :: a
  complex(8),allocatable :: spectrum(:,:),realspace(:,:)
  call spatial_local_destroy(plan,status)
  status=merge(1,0,status/=0)
  do a=1,size(comm)
   call comm_get_max(status,comm(a))
  enddo
  if(status/=0)return
  allocate(spectrum(size(multiplier),1),realspace(size(multiplier),1))
  spectrum(:,1)=cmplx(multiplier,0d0,8)
  ! Includes global 1/N and the supplied G=0 value without modification.
  call mesh_transform(n,dims,coords,comm,spectrum,realspace,1,status)
  if(status/=0)return
  plan%n=n
  call block_layout(n,dims,coords,plan%m,plan%lo,status)
  plan%kernel=realspace(:,1)
 end subroutine
 subroutine spatial_local_apply(plan,comm_r,source,targets,action,used,status,pairs_executed,pair_fft_points,skip, &
                                fft_cost_limit)
  implicit none
  type(s_exx_spatial_local),intent(inout) :: plan
  integer,intent(in) :: comm_r
  complex(8),intent(in) :: source(:),targets(:,:)
  complex(8),intent(out) :: action(:,:)
  logical,intent(out) :: used
  integer,intent(out) :: status
  integer(int64),intent(out),optional :: pairs_executed,pair_fft_points
  logical,intent(in),optional :: skip(:)
  real(8),intent(in),optional :: fft_cost_limit
  integer,allocatable :: occupied(:,:),total(:,:),points(:)
  integer :: g,x,y,z,p(3),a,b,c,axis,j,gap,best,box(3),padded(3)
  integer :: rank,peers,group,member,bad,used_flag,lower(3),upper(3)
  integer(int64) :: pairs,fft_points
  real(8) :: counts(2),global_counts(2),fft_volume
  action=0d0;used=.false.;status=0;pairs=0_int64;fft_points=0_int64
  if(present(pairs_executed))pairs_executed=0_int64
  if(present(pair_fft_points))pair_fft_points=0_int64
  bad=0
  if(any(plan%n<1).or.any(plan%m<1))bad=1
  if(size(source)/=product(plan%m).or.size(targets,1)/=size(source))bad=1
  if(any(shape(action)/=shape(targets)))bad=1
  call comm_get_max(bad,comm_r)
  if(bad/=0)then
   status=1;return
  endif
  allocate(occupied(maxval(plan%n),3),total(maxval(plan%n),3));occupied=0
  g=0;member=0
  do z=0,plan%m(3)-1;do y=0,plan%m(2)-1;do x=0,plan%m(1)-1
   g=g+1
   if(source(g)==(0d0,0d0))cycle
   member=1;p=[x,y,z]+plan%lo
   do axis=1,3
    occupied(p(axis)+1,axis)=1
   enddo
  enddo;enddo;enddo
  call comm_summation(occupied,total,size(total),comm_r)
  if(.not.any(total/=0))then
   used=.true.;return
  endif
  do axis=1,3
   points=pack([(j-1,j=1,plan%n(axis))],total(1:plan%n(axis),axis)>0)
   best=0
   do j=1,size(points)
    gap=points(mod(j,size(points))+1)-points(j)
    if(j==size(points))gap=gap+plan%n(axis)
    best=max(best,gap)
   enddo
   box(axis)=plan%n(axis)-best+1;padded(axis)=compact_axis_size(plan%n(axis),box(axis))
   deallocate(points)
  enddo
  if(product(int(padded,int64))>=product(int(plan%n,int64)))return
  fft_volume=real(product(int(padded,int64)),8)
  if(present(fft_cost_limit))then
   if(fft_volume*log(max(2d0,fft_volume))>fft_cost_limit)return
  endif
  ! Kernel displacements may have owners outside the source support. Include
  ! those owners, but never allocate compact tiles on unrelated spatial ranks.
  call compact_kernel_bounds(plan%n,box,lower,upper)
  if(member==0)then
   kernel_owner: do c=lower(3),upper(3)
    do b=lower(2),upper(2);do a=lower(1),upper(1)
     p=modulo([a,b,c],plan%n)-plan%lo
     if(all(p>=0).and.all(p<plan%m))then
      member=1
      exit kernel_owner
     endif
    enddo;enddo
   enddo kernel_owner
  endif
  call comm_get_groupinfo(comm_r,rank,peers)
  group=comm_create_group(comm_r,member,rank)
  if(member==1)then
   call local_apply_group(plan,group,source,targets,action,used,status,pairs,fft_points,skip)
  else
   used=.true.
  endif
  call comm_free_group(group)
  call comm_get_max(status,comm_r)
  used_flag=merge(0,1,used);call comm_get_max(used_flag,comm_r);used=used_flag==0
  counts=real([pairs,fft_points],8)
  call comm_get_max(counts,global_counts,2,comm_r)
  if(present(pairs_executed))pairs_executed=int(global_counts(1),int64)
  if(present(pair_fft_points))pair_fft_points=int(global_counts(2),int64)
 end subroutine
 subroutine local_apply_group(plan,comm_r,source,targets,action,used,status,pairs_executed,pair_fft_points,skip)
!$ use omp_lib, only: omp_get_max_threads,omp_get_thread_num
  implicit none
  type(s_exx_spatial_local),intent(inout) :: plan
  integer,intent(in) :: comm_r
  complex(8),intent(in) :: source(:),targets(:,:)
  complex(8),intent(out) :: action(:,:)
  logical,intent(out) :: used
  logical,intent(in),optional :: skip(:)
  integer,intent(out) :: status
  integer(int64),intent(out),optional :: pairs_executed,pair_fft_points
  integer,allocatable :: occupied(:,:),total(:,:),axis_points(:),local_rows(:),tile_rows(:)
  complex(8),allocatable :: kernel_tile(:,:,:),kernel_sum(:,:,:),tile(:,:),tile_sum(:,:),result(:,:),result_sum(:,:)
  complex(8),allocatable :: source_tile(:),source_sum(:),batch_targets(:,:),batch_result(:,:)
  integer,allocatable :: owned(:)
  integer :: ng,nt,bad,axis,j,g,x,y,z,p(3),box(3),origin(3),padded(3),gap,best,start,a,b,c,r,rank,peers
  integer :: count,ntmax,first,nb,k,batch,nowned,capacity,first_owned,lower(3),upper(3),extent(3)
  integer :: workers,worker,executed
  logical :: batched,pair_used
  integer(int64) :: local_counts(2),global_counts(2)
  ! communication wrappers do not expose int64 reductions; sums are bounded
  ! here by the target count and accumulated through small real vectors.
  real(8) :: count_values(2),count_totals(2)
  status=1;used=.false.;action=0d0
  if(present(pairs_executed))pairs_executed=0_int64
  if(present(pair_fft_points))pair_fft_points=0_int64
  bad=0;ng=size(source);nt=size(targets,2)
  batched=associated(plan%batch_action).or.allocated(plan%backend)
  if(batched.and.plan%batch_size<1)bad=1
  if(associated(plan%batch_action).and.allocated(plan%backend))bad=1
  if(.not.allocated(plan%kernel))bad=1
  if(size(targets,1)/=ng.or.any(shape(action)/=shape(targets)))bad=1
  if(ng/=product(plan%m))bad=1
  if(present(skip))then
   if(size(skip)/=nt)bad=1
  endif
  ntmax=nt;call comm_get_max(ntmax,comm_r)
  if(ntmax/=nt)bad=1
  call comm_get_max(bad,comm_r)
  if(bad/=0)return
  status=0
  if(nt==0)then
   used=.true.;return
  endif
  allocate(occupied(maxval(plan%n),3),total(maxval(plan%n),3));occupied=0
  g=0
  do z=0,plan%m(3)-1;do y=0,plan%m(2)-1;do x=0,plan%m(1)-1
   g=g+1
   if(source(g)==(0d0,0d0))cycle
   p=[x,y,z]+plan%lo
   do axis=1,3
    occupied(p(axis)+1,axis)=1
   enddo
  enddo;enddo;enddo
  call comm_summation(occupied,total,size(total),comm_r)
  if(.not.any(total/=0))then
   used=.true.;return
  endif
  do axis=1,3
   axis_points=pack([(j-1,j=1,plan%n(axis))],total(1:plan%n(axis),axis)>0)
   best=-1;start=1
   do j=1,size(axis_points)
    if(j<size(axis_points))then
     gap=axis_points(j+1)-axis_points(j)
    else
     gap=axis_points(1)+plan%n(axis)-axis_points(j)
    endif
    if(gap>best)then
     best=gap;start=mod(j,size(axis_points))+1
    endif
   enddo
   origin(axis)=axis_points(start);box(axis)=plan%n(axis)-best+1
   padded(axis)=compact_axis_size(plan%n(axis),box(axis))
   deallocate(axis_points)
  enddo
  if(product(int(padded,int64))>=product(int(plan%n,int64)))return
  ! Gather only kernel displacements required by this compact source box.
  call compact_kernel_bounds(plan%n,box,lower,upper)
  extent=upper-lower+1
  allocate(kernel_tile(extent(1),extent(2),extent(3)),kernel_sum(extent(1),extent(2),extent(3)))
  kernel_tile=0d0
  do c=lower(3),upper(3);do b=lower(2),upper(2);do a=lower(1),upper(1)
   p=modulo([a,b,c],plan%n)-plan%lo
   if(any(p<0).or.any(p>=plan%m))cycle
   g=1+p(1)+plan%m(1)*(p(2)+plan%m(2)*p(3))
   kernel_tile(a-lower(1)+1,b-lower(2)+1,c-lower(3)+1)=plan%kernel(g)
  enddo;enddo;enddo
  call comm_summation(kernel_tile,kernel_sum,size(kernel_sum),comm_r)
  call exx_local_prepare_compact(plan%fft,plan%n,box,kernel_sum,status)
  bad=status;call comm_get_max(bad,comm_r)
  if(bad/=0)then
   status=1;return
  endif
  count=product(box)
  ! CPU workers process independent pairs with reusable private FFT buffers.
  ! Keep the batch bounded by the active worker count, not all target orbitals.
  call comm_get_groupinfo(comm_r,rank,peers)
  workers=1
!$ workers=omp_get_max_threads()
  if(batched)workers=1
  call comm_get_max(workers,comm_r)
  workers=min(workers,1+(nt-1)/peers)
  batch=int(min(int(nt,int64),int(peers,int64)*int(workers,int64)))
  if(batched) &
   batch=int(min(int(nt,int64),int(max(1,peers),int64)*int(plan%batch_size,int64)))
  if(.not.batched)then
   call exx_local_cpu_prepare(plan%fft,workers,bad)
   call comm_get_max(bad,comm_r)
   if(bad/=0)then
    status=1;return
   endif
  endif
  ! Communication wrappers take a default-integer element count.
  if(int(count,int64)*int(batch,int64)>int(huge(0),int64))then
   status=1;return
  endif
  allocate(tile(count,batch),tile_sum(count,batch),result(count,batch),result_sum(count,batch))
  allocate(source_tile(count),source_sum(count))
  allocate(local_rows(count),tile_rows(count))
  capacity=min(plan%batch_size,1+(nt-1)/max(1,peers))
  if(batched)allocate(owned(capacity))
  r=0
  do z=0,box(3)-1;do y=0,box(2)-1;do x=0,box(1)-1
   p=modulo(origin+[x,y,z],plan%n)-plan%lo
   if(any(p<0).or.any(p>=plan%m))cycle
   g=1+p(1)+plan%m(1)*(p(2)+plan%m(2)*p(3))
   j=1+x+box(1)*(y+box(2)*z)
   r=r+1;local_rows(r)=g;tile_rows(r)=j
  enddo;enddo;enddo
  ! The source is shared by every target batch: gather and upload it once.
  source_tile=0d0
  do j=1,r
   source_tile(tile_rows(j))=source(local_rows(j))
  enddo
  call comm_summation(source_tile,source_sum,count,comm_r)
  bad=0
  if(allocated(plan%backend).and.rank<nt)then
   call plan%backend%prepare(plan%fft%padded,plan%fft%indices,plan%fft%filter,source_sum,capacity,status)
   if(status/=0)bad=1
  endif
  call comm_get_max(bad,comm_r)
  if(bad/=0)then
   status=1;return
  endif
  local_counts=0_int64;bad=0
  do first=1,nt,batch
   nb=min(batch,nt-first+1);tile=0d0
   do j=1,r
    g=local_rows(j);k=tile_rows(j)
    tile(k,1:nb)=targets(g,first:first+nb-1)
   enddo
   call comm_summation(tile(:,1:nb),tile_sum(:,1:nb),count*nb,comm_r)
   result=0d0
   first_owned=1+modulo(rank-modulo(first-1,peers),peers)
   if(batched)then
    ! The owner rule is unchanged; each worker collects several of its pairs.
    nowned=0
    do j=first_owned,nb,peers
     if(present(skip))then
      if(skip(first+j-1))cycle
     endif
     if(.not.any(source_sum/=(0d0,0d0).and.tile_sum(:,j)/=(0d0,0d0)))cycle
     nowned=nowned+1;owned(nowned)=j
    enddo
    if(nowned>0)then
     allocate(batch_targets(count,nowned),batch_result(count,nowned))
     do k=1,nowned
      batch_targets(:,k)=tile_sum(:,owned(k))
     enddo
     if(allocated(plan%backend))then
      call plan%backend%apply(batch_targets,batch_result,status)
     else
      call plan%batch_action(plan%fft%padded,plan%fft%indices,plan%fft%filter, &
        source_sum,batch_targets,batch_result,status)
     endif
     if(status/=0)then
      bad=1
     else
      do k=1,nowned
       result(:,owned(k))=batch_result(:,k)
      enddo
      local_counts(1)=local_counts(1)+int(nowned,int64)
      local_counts(2)=local_counts(2)+int(nowned,int64)*int(plan%fft%fft_points,int64)
     endif
     deallocate(batch_targets,batch_result)
    endif
   else
   executed=0
!$omp parallel do default(none) schedule(static) num_threads(workers) if(workers>1) &
!$omp shared(first_owned,nb,peers,first,skip,source_sum,tile_sum,result,plan) &
!$omp private(j,worker,pair_used) reduction(+:executed)
   do j=first_owned,nb,peers
    if(present(skip))then
     if(skip(first+j-1))cycle
    endif
    worker=1
!$  worker=omp_get_thread_num()+1
    call exx_local_cpu_pair(plan%fft,worker,source_sum,tile_sum(:,j),result(:,j),pair_used)
    if(pair_used)executed=executed+1
   enddo
!$omp end parallel do
   local_counts(1)=local_counts(1)+int(executed,int64)
   local_counts(2)=local_counts(2)+int(executed,int64)*int(plan%fft%fft_points,int64)
   endif
   call comm_get_max(bad,comm_r)
   if(bad/=0)then
    status=1;return
   endif
   call comm_summation(result(:,1:nb),result_sum(:,1:nb),count*nb,comm_r)
   do j=1,r
    action(local_rows(j),first:first+nb-1)=result_sum(tile_rows(j),1:nb)
   enddo
  enddo
  count_values=real(local_counts,8)
  call comm_summation(count_values,count_totals,2,comm_r)
  global_counts=int(count_totals,int64)
  if(present(pairs_executed))pairs_executed=global_counts(1)
  if(present(pair_fft_points))pair_fft_points=global_counts(2)
  used=.true.;status=0
 end subroutine
end module
