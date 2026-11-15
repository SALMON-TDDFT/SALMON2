! Sparse spatial-block envelopes for bounded EXX pair generation.
! Target columns/order and threshold agree within comm_r. Only block metadata
! is replicated; no target grid values and no source-by-target table are gathered.
module exx_pair_candidates
 use iso_fortran_env, only:int64
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 use communication, only:comm_bcast,comm_get_groupinfo,comm_get_max
 implicit none
 private
 public :: exx_pair_catalog,pair_catalog_build,pair_catalog_query,pair_source_box
 integer,parameter :: block_edge=8
 type exx_pair_catalog
  integer :: blocks(3)=0,targets=0,epoch=0
  real(8) :: floor=0d0
  integer(int64) :: entries=0,last_visited=0
  integer,allocatable :: offset(:),column(:),mark(:)
  real(8),allocatable :: amplitude(:)
 end type
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d
  end interface
contains
 subroutine pair_catalog_build(cat,n,lo,m,comm_r,target,threshold,status)
  implicit none
  type(exx_pair_catalog),intent(out) :: cat
  integer,intent(in) :: n(3),lo(3),m(3),comm_r
  complex(8),intent(in) :: target(:,:)
  real(8),intent(in) :: threshold
  integer,intent(out) :: status
  integer,allocatable :: block_of_row(:),cells(:),columns(:),counts(:),cursor(:),local_cell(:),local_column(:)
  real(8),allocatable :: maxima(:),values(:),local_value(:)
  integer :: bad,ncell,g,x,y,z,p(3),j,c,k,pass,local_count,total,rank,peers,owner,first,last
  status=1;bad=0
  if(any(n<1).or.any(lo<0).or.any(m<0).or.any(lo+m>n))bad=1
  if(size(target,1)/=product(m))bad=1
  if(.not.ieee_is_finite(threshold).or.threshold<0d0)bad=1
  if(.not.salmon_all_finite(real(target)).or..not.salmon_all_finite(aimag(target)))bad=1
  call comm_get_max(bad,comm_r)
  if(bad/=0)return
  cat%blocks=(n+block_edge-1)/block_edge;cat%targets=size(target,2);cat%floor=threshold
  ncell=product(cat%blocks)
  allocate(block_of_row(size(target,1)),maxima(ncell));g=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
   g=g+1;p=([x,y,z]+lo)/block_edge
   block_of_row(g)=1+p(1)+cat%blocks(1)*(p(2)+cat%blocks(2)*p(3))
  enddo;enddo;enddo
  ! Two passes retain only significant block envelopes, rather than allocating
  ! an ncell-by-targets temporary and compacting it after allocation.
  do pass=1,2
   local_count=0
   do j=1,cat%targets
    maxima=0d0
    do g=1,size(target,1)
     c=block_of_row(g);maxima(c)=max(maxima(c),abs(target(g,j)))
    enddo
    do c=1,ncell
     if(maxima(c)<=threshold)cycle
     local_count=local_count+1
     if(pass==2)then
      local_cell(local_count)=c;local_column(local_count)=j;local_value(local_count)=maxima(c)
     endif
    enddo
   enddo
   if(pass==1)allocate(local_cell(local_count),local_column(local_count),local_value(local_count))
  enddo
  call comm_get_groupinfo(comm_r,rank,peers)
  allocate(counts(0:peers-1));total=0
  do owner=0,peers-1
   k=local_count;call comm_bcast(k,comm_r,owner);counts(owner)=k
   if(k>huge(total)-total)bad=1
   if(bad==0)total=total+k
  enddo
  call comm_get_max(bad,comm_r)
  if(bad/=0)return
  allocate(cells(total),columns(total),values(total));first=1
  do owner=0,peers-1
   last=first+counts(owner)-1
   if(counts(owner)>0)then
    if(rank==owner)then
     cells(first:last)=local_cell;columns(first:last)=local_column;values(first:last)=local_value
    endif
    call comm_bcast(cells(first:last),comm_r,owner)
    call comm_bcast(columns(first:last),comm_r,owner)
    call comm_bcast(values(first:last),comm_r,owner)
   endif
   first=last+1
  enddo
  allocate(cat%offset(ncell+1),cat%column(total),cat%amplitude(total),cat%mark(cat%targets),cursor(ncell))
  cat%mark=0;cursor=0
  do k=1,total
   cursor(cells(k))=cursor(cells(k))+1
  enddo
  cat%offset(1)=1
  do c=1,ncell
   cat%offset(c+1)=cat%offset(c)+cursor(c)
  enddo
  cursor=cat%offset(:ncell)
  do k=1,total
   c=cells(k);j=cursor(c);cat%column(j)=columns(k);cat%amplitude(j)=values(k);cursor(c)=j+1
  enddo
  ! Duplicate block/column entries at rank boundaries are safe. Query marks
  ! remove duplicate candidates without any dense pair adjacency matrix.
  cat%entries=int(total,int64);status=0
 end subroutine

 subroutine pair_source_box(n,lo,m,comm_r,source,lower,upper,status)
  implicit none
  integer,intent(in) :: n(3),lo(3),m(3),comm_r
  complex(8),intent(in) :: source(:)
  integer,intent(out) :: lower(3),upper(3),status
  real(8) :: bounds(6),global_bounds(6)
  integer :: bad,g,x,y,z,p(3)
  status=1;bad=0;lower=n;upper=-1
  if(any(n<1).or.any(lo<0).or.any(m<0).or.any(lo+m>n).or.size(source)/=product(m))bad=1
  if(.not.salmon_all_finite(real(source)).or..not.salmon_all_finite(aimag(source)))bad=1
  call comm_get_max(bad,comm_r)
  if(bad/=0)return
  g=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
   g=g+1
   if(source(g)==(0d0,0d0))cycle
   p=[x,y,z]+lo;lower=min(lower,p);upper=max(upper,p)
  enddo;enddo;enddo
  bounds(:3)=-real(lower,8);bounds(4:)=real(upper,8)
  call comm_get_max(bounds,global_bounds,6,comm_r)
  lower=-nint(global_bounds(:3));upper=nint(global_bounds(4:));status=0
 end subroutine

 subroutine pair_catalog_query(cat,lower,upper,threshold,selected,nselected,status)
  implicit none
  type(exx_pair_catalog),intent(inout) :: cat
  integer,intent(in) :: lower(3),upper(3)
  real(8),intent(in) :: threshold
  integer,intent(out) :: selected(:),nselected,status
  integer :: low(3),high(3),x,y,z,c,k,j
  status=1;nselected=0;cat%last_visited=0
  if(.not.allocated(cat%offset).or.size(selected)<cat%targets)return
  if(.not.ieee_is_finite(threshold).or.threshold<cat%floor)return
  if(any(upper<lower))then
   status=0;return
  endif
  if(any(lower<0).or.any(upper/block_edge>=cat%blocks))return
  if(cat%epoch==huge(cat%epoch))then
   cat%mark=0;cat%epoch=0
  endif
  cat%epoch=cat%epoch+1;low=lower/block_edge;high=upper/block_edge
  do z=low(3),high(3);do y=low(2),high(2);do x=low(1),high(1)
   c=1+x+cat%blocks(1)*(y+cat%blocks(2)*z)
   do k=cat%offset(c),cat%offset(c+1)-1
    cat%last_visited=cat%last_visited+1_int64
    if(cat%amplitude(k)<=threshold)cycle
    j=cat%column(k)
    if(cat%mark(j)==cat%epoch)cycle
    cat%mark(j)=cat%epoch;nselected=nselected+1;selected(nselected)=j
   enddo
  enddo;enddo;enddo
  status=0
 end subroutine

  pure logical function finite_real_1d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:)
    real(8) :: value
    integer :: i
    finite=.false.
    do i=1,size(values,1)
      value=values(i)
      if(.not.ieee_is_finite(value))return
    enddo
    finite=.true.
  end function

  pure logical function finite_real_2d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:)
    real(8) :: value
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        value=values(i,j)
        if(.not.ieee_is_finite(value))return
      enddo
    enddo
    finite=.true.
  end function

end module
