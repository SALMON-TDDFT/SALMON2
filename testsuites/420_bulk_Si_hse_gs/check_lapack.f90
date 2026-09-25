! Check eigenvectors as well as eigenvalues: a failed fallback can otherwise
! appear to complete a short SCF/RT smoke test without converging the GS.
program check_hse_lapack
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: n=16
  complex(8) :: a(n,n),original(n,n),work(2*n-1),residual(n),metric(n,n)
  real(8) :: eigenvalues(n),rwork(3*n-2),error
  integer :: i,j,info

  do j=1,n
    do i=1,n
      a(i,j)=cmplx(sin(real(17*i+13*j,8)),cos(real(7*i-11*j,8)),8)
    end do
  end do
  a=a+transpose(conjg(a))
  original=a
  call zheev('V','U',n,a,n,eigenvalues,work,size(work),rwork,info)
  if(info/=0) error stop 'LAPACK ZHEEV failed'
  if(.not.all(ieee_is_finite(eigenvalues))) error stop 'Nonfinite LAPACK eigenvalues'
  if(.not.all(ieee_is_finite(real(a,8)))) error stop 'Nonfinite LAPACK eigenvectors'
  if(.not.all(ieee_is_finite(aimag(a)))) error stop 'Nonfinite LAPACK eigenvectors'
  error=0d0
  do i=1,n
    residual=matmul(original,a(:,i))-eigenvalues(i)*a(:,i)
    error=max(error,maxval(abs(residual))/maxval(abs(original)))
  end do
  ! The negated comparison also rejects NaN.
  if(.not.(error<1d-10)) error stop 'LAPACK eigenvector residual failed'
  metric=matmul(transpose(conjg(a)),a)
  do i=1,n
    metric(i,i)=metric(i,i)-1d0
  end do
  if(.not.(maxval(abs(metric))<1d-10)) error stop 'LAPACK orthonormality failed'
  print *, 'LAPACK eigenvectors passed; relative residual=',error
end program check_hse_lapack
