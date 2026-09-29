! Compact source convolution on native Cartesian blocks. Only compact tiles
! are replicated; the inverse global kernel remains spatially distributed.
! All calls are collective over comm_r with matching target column counts.
module exx_spatial_local
 use iso_fortran_env, only: int64
 use exx_batch_backend, only: s_exx_batch_backend
 use communication, only: comm_summation,comm_get_max,comm_get_groupinfo,comm_create_group,comm_free_group
 use fftw_blocks, only: pencil_transform=>mesh_transform,block_layout
 use exx_local_fft, only: s_exx_local_fft,exx_local_prepare_compact,exx_local_apply,exx_local_destroy, &
                          compact_axis_size,compact_kernel_bounds
 implicit none
 private
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
  real(8),intent(in) :: multiplier(:) ! xyz block order (3D) or legacy Z pencil order (2D)
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
  call pencil_transform(n,dims,coords,comm,spectrum,realspace,1,status,spectral_z=.true.)
  if(status/=0)return
  plan%n=n
  call block_layout(n,dims,coords,plan%m,plan%lo,status)
  plan%kernel=realspace(:,1)
 end subroutine
 subroutine spatial_local_apply(plan,comm_r,source,targets,action,used,status,pairs_executed,pair_fft_points,skip)
  implicit none
  type(s_exx_spatial_local),intent(inout) :: plan
  integer,intent(in) :: comm_r
  complex(8),intent(in) :: source(:),targets(:,:)
  complex(8),intent(out) :: action(:,:)
  logical,intent(out) :: used
  integer,intent(out) :: status
  integer(int64),intent(out),optional :: pairs_executed,pair_fft_points
  logical,intent(in),optional :: skip(:)
  integer,allocatable :: occupied(:,:),total(:,:),points(:)
  integer :: g,x,y,z,p(3),a,b,c,axis,j,gap,best,box(3),padded(3)
  integer :: rank,peers,group,member,bad,used_flag,lower(3),upper(3)
  integer(int64) :: pairs,fft_points
  real(8) :: counts(2),global_counts(2)
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
  ! Kernel displacements may have owners outside the source support. Include
  ! those owners, but never allocate compact tiles on unrelated spatial ranks.
  call compact_kernel_bounds(plan%n,box,lower,upper)
  if(member==0)then
   do c=lower(3),upper(3);do b=lower(2),upper(2);do a=lower(1),upper(1)
    p=modulo([a,b,c],plan%n)-plan%lo
    if(all(p>=0).and.all(p<plan%m))member=1
   enddo;enddo;enddo
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
  complex(8),allocatable :: density(:),potential(:),source_tile(:),source_sum(:),batch_targets(:,:),batch_result(:,:)
  integer,allocatable :: owned(:)
  integer :: ng,nt,bad,axis,j,g,x,y,z,p(3),box(3),origin(3),padded(3),gap,best,start,a,b,c,r,rank,peers
  integer :: count,ntmax,first,nb,k,batch,nowned,capacity,lower(3),upper(3),extent(3)
  logical :: batched
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
  ! One selected pair per spatial worker; the target list has already been
  ! screened. No all-target compact action or fixed four-worker bottleneck.
  call comm_get_groupinfo(comm_r,rank,peers)
  batch=min(nt,max(1,peers))
  if(batched) &
   batch=int(min(int(nt,int64),int(max(1,peers),int64)*int(plan%batch_size,int64)))
  ! Communication wrappers take a default-integer element count.
  if(int(count,int64)*int(batch,int64)>int(huge(0),int64))then
   status=1;return
  endif
  allocate(tile(count,batch),tile_sum(count,batch),result(count,batch),result_sum(count,batch))
  allocate(source_tile(count),source_sum(count))
  allocate(local_rows(count),tile_rows(count),density(count),potential(count))
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
   if(batched)then
    ! The owner rule is unchanged; each worker collects several of its pairs.
    nowned=0
    do j=1,nb
     if(modulo(first+j-2,peers)/=rank)cycle
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
   do j=1,nb
    if(modulo(first+j-2,peers)/=rank)cycle
    if(present(skip))then
     if(skip(first+j-1))cycle
    endif
    density=conjg(source_sum)*tile_sum(:,j)
    if(all(density==(0d0,0d0)))cycle
    call exx_local_apply(plan%fft,density,potential,status)
    if(status/=0)then
     bad=1;exit
    endif
    result(:,j)=-source_sum*potential
    local_counts(1)=local_counts(1)+1_int64
    local_counts(2)=local_counts(2)+int(plan%fft%fft_points,int64)
   enddo
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
