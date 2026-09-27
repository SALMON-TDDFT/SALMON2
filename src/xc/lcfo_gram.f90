! Orthogonality validation using only the independent Hermitian entries.
module lcfo_gram
 implicit none
 private
 public :: lcfo_gram_error
contains
 subroutine lcfo_gram_error(coeff,comm,error)
  use iso_fortran_env, only:int64
  use communication, only:comm_summation
  use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
  implicit none
  complex(8),intent(in),contiguous :: coeff(:,:)
  integer,intent(in) :: comm
  real(8),intent(out) :: error
  complex(8),allocatable :: upper(:,:),packed(:,:),total(:,:)
  integer :: n,m,j,offset,npair
  integer(int64) :: count64
  external :: zherk
  n=size(coeff,2);m=size(coeff,1)
  if(n<1)error stop 'LCFO Gram: empty occupied space'
  count64=int(n,int64)*(int(n,int64)+1_int64)/2_int64
  if(count64>int(huge(npair),int64))error stop 'LCFO Gram: packed count exceeds communication limit'
  npair=int(count64)
  allocate(upper(n,n),packed(npair,1),total(npair,1))
  call zherk('U','C',n,m,1d0,coeff,max(1,m),0d0,upper,n)
  offset=0
  do j=1,n
   packed(offset+1:offset+j,1)=upper(1:j,j)
   offset=offset+j
  enddo
  call comm_summation(packed,total,npair,comm)
  ! A NaN can be hidden by MAXVAL; check every independent entry first.
  if(.not.all(ieee_is_finite(real(total))).or..not.all(ieee_is_finite(aimag(total))))then
   error=huge(1d0);return
  endif
  offset=0
  do j=1,n
   offset=offset+j;total(offset,1)=total(offset,1)-1d0
  enddo
  error=maxval(abs(total))
 end subroutine
end module
