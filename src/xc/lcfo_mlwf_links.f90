! Initial Gamma links: bounded column tiles, collective reduction only to root.
module lcfo_mlwf_links
 use iso_fortran_env,only:int64
 use communication,only:comm_summation,comm_get_groupinfo
 implicit none
 private
 public :: lcfo_initial_links
contains
 subroutine lcfo_initial_links(grid,position,lengths,dv,comm,raw,tile_width,scratch_elements)
  implicit none
  complex(8),intent(in),contiguous :: grid(:,:)
  real(8),intent(in) :: position(:,:),lengths(3),dv
  integer,intent(in) :: comm
  complex(8),allocatable,intent(out) :: raw(:,:,:,:)
  integer,optional,intent(in) :: tile_width
  integer(int64),optional,intent(out) :: scratch_elements
  complex(8),allocatable :: shifted(:,:),local(:,:),total(:,:),phase(:)
  integer :: n,ng,width,rank,np,a,first,count,j
  real(8) :: delta
  external :: zgemm
  n=size(grid,2);ng=size(grid,1);width=64
  if(present(tile_width))width=tile_width
  if(n<1.or.width<1.or.any(lengths<=0d0).or.dv<=0d0)error stop 'LCFO links: invalid dimensions'
  if(any(shape(position)/=[3,ng]))error stop 'LCFO links: invalid position shape'
  width=min(width,n)
  call comm_get_groupinfo(comm,rank,np)
  if(rank==0)then
   allocate(raw(n,n,6,1))
  else
   allocate(raw(0,0,0,0))
  endif
  allocate(shifted(ng,width),local(n,width),total(n,width),phase(ng))
  if(present(scratch_elements))scratch_elements= &
    (int(ng,int64)+2_int64*n)*width+ng
  do a=1,3
   delta=2*acos(-1d0)/lengths(a)
   phase=exp(cmplx(0d0,-delta*position(a,:),8))
   do first=1,n,width
    count=min(width,n-first+1)
    do j=1,count
     shifted(:,j)=grid(:,first+j-1)*phase
    enddo
    call zgemm('C','N',n,count,ng,cmplx(dv,0d0,8),grid,max(1,ng), &
      shifted,max(1,ng),(0d0,0d0),local,n)
    call comm_summation(local(:,1:count),total(:,1:count),n*count,comm,0)
    if(rank==0)then
     raw(:,first:first+count-1,a,1)=total(:,1:count)
     do j=1,count
      raw(first+j-1,:,a+3,1)=conjg(total(:,j))
     enddo
    endif
   enddo
  enddo
 end subroutine
end module
