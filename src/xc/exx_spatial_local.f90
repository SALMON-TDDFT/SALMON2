! Compact source convolution on x-complete y/z pencils. Only compact tiles
! are replicated; the inverse global kernel remains spatially distributed.
! All calls are collective over comm_r with matching target column counts.
module exx_spatial_local
 use iso_fortran_env, only: int64
 use communication, only: comm_summation,comm_get_max,comm_get_groupinfo
 use fftw_pencils, only: pencil_transform
 use exx_local_fft, only: s_exx_local_fft,exx_local_prepare_compact,exx_local_apply,exx_local_destroy,smooth_size
 implicit none
 private
 public :: s_exx_spatial_local,spatial_local_init,spatial_local_apply,spatial_local_destroy
 type s_exx_spatial_local
  integer :: n(3)=0,m(3)=0,lo(3)=0
  complex(8),allocatable :: kernel(:)
  type(s_exx_local_fft) :: fft
 end type
contains
 subroutine spatial_local_destroy(plan)
  implicit none
  type(s_exx_spatial_local),intent(inout) :: plan
  call exx_local_destroy(plan%fft)
  if(allocated(plan%kernel))deallocate(plan%kernel)
  plan%n=0;plan%m=0;plan%lo=0
 end subroutine
 subroutine spatial_local_init(plan,n,dims,coords,comm,multiplier,status)
  implicit none
  type(s_exx_spatial_local),intent(inout) :: plan
  integer,intent(in) :: n(3),dims(2),coords(2),comm(2)
  real(8),intent(in) :: multiplier(:) ! Z spectral order, (z,x,y)
  integer,intent(out) :: status
  complex(8),allocatable :: spectrum(:,:),realspace(:,:)
  call spatial_local_destroy(plan)
  allocate(spectrum(size(multiplier),1),realspace(size(multiplier),1))
  spectrum(:,1)=cmplx(multiplier,0d0,8)
  ! Includes global 1/N and the supplied G=0 value without modification.
  call pencil_transform(n,dims,coords,comm,spectrum,realspace,1,status,spectral_z=.true.)
  if(status/=0)return
  plan%n=n;plan%m=[n(1),n(2)/dims(1),n(3)/dims(2)]
  plan%lo=[0,coords(1)*plan%m(2),coords(2)*plan%m(3)]
  plan%kernel=realspace(:,1)
 end subroutine
 subroutine spatial_local_apply(plan,comm_r,source,targets,action,used,status,pairs_executed,pair_fft_points,skip)
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
  complex(8),allocatable :: density(:),potential(:)
  integer :: ng,nt,bad,axis,j,g,x,y,z,p(3),box(3),origin(3),padded(3),gap,best,start,a,b,c,r,rank,peers
  integer :: count,ntmax,first,nb,k,batch
  integer(int64) :: local_counts(2),global_counts(2)
  ! communication wrappers do not expose int64 reductions; sums are bounded
  ! here by the target count and accumulated through small real vectors.
  real(8) :: count_values(2),count_totals(2)
  status=1;used=.false.;action=0d0
  if(present(pairs_executed))pairs_executed=0_int64
  if(present(pair_fft_points))pair_fft_points=0_int64
  bad=0;ng=size(source);nt=size(targets,2)
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
   padded(axis)=smooth_size(2*box(axis)-1)
   deallocate(axis_points)
  enddo
  if(product(int(padded,int64))>=product(int(plan%n,int64)))return
  ! Gather only kernel displacements required by this compact source box.
  allocate(kernel_tile(2*box(1)-1,2*box(2)-1,2*box(3)-1),kernel_sum(2*box(1)-1,2*box(2)-1,2*box(3)-1))
  kernel_tile=0d0
  do c=1-box(3),box(3)-1;do b=1-box(2),box(2)-1;do a=1-box(1),box(1)-1
   p=modulo([a,b,c],plan%n)-plan%lo
   if(any(p<0).or.any(p>=plan%m))cycle
   g=1+p(1)+plan%m(1)*(p(2)+plan%m(2)*p(3))
   kernel_tile(a+box(1),b+box(2),c+box(3))=plan%kernel(g)
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
  allocate(tile(count,batch+1),tile_sum(count,batch+1),result(count,batch),result_sum(count,batch))
  allocate(local_rows(count),tile_rows(count),density(count),potential(count))
  r=0
  do z=0,box(3)-1;do y=0,box(2)-1;do x=0,box(1)-1
   p=modulo(origin+[x,y,z],plan%n)-plan%lo
   if(any(p<0).or.any(p>=plan%m))cycle
   g=1+p(1)+plan%m(1)*(p(2)+plan%m(2)*p(3))
   j=1+x+box(1)*(y+box(2)*z)
   r=r+1;local_rows(r)=g;tile_rows(r)=j
  enddo;enddo;enddo
  call comm_get_groupinfo(comm_r,rank,peers)
  local_counts=0_int64;bad=0
  do first=1,nt,batch
   nb=min(batch,nt-first+1);tile=0d0
   do j=1,r
    g=local_rows(j);k=tile_rows(j)
    tile(k,1)=source(g);tile(k,2:nb+1)=targets(g,first:first+nb-1)
   enddo
   call comm_summation(tile(:,1:nb+1),tile_sum(:,1:nb+1),count*(nb+1),comm_r)
   result=0d0
   do j=1,nb
    if(modulo(first+j-2,peers)/=rank)cycle
    if(present(skip))then
     if(skip(first+j-1))cycle
    endif
    density=conjg(tile_sum(:,1))*tile_sum(:,j+1)
    if(all(density==(0d0,0d0)))cycle
    call exx_local_apply(plan%fft,density,potential,status)
    if(status/=0)then
     bad=1;exit
    endif
    result(:,j)=-tile_sum(:,1)*potential
    local_counts(1)=local_counts(1)+1_int64
    local_counts(2)=local_counts(2)+int(plan%fft%fft_points,int64)
   enddo
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
