! Production-sized BLAS path: catches failures invisible on tiny matrices.
program probe
  use hse_ace
  implicit none
  integer,parameter :: ng=1728,no=16,nk=8
  type(hse_ace_state) :: ace
  complex(8),allocatable :: t(:,:,:),a(:,:,:),ref(:,:,:)
  integer :: i,j,k,ierr
  real(8) :: error
  allocate(ace%factors(ng,no,nk),t(ng,no,nk),a(ng,no,nk),ref(ng,no,nk))
  ace%dv=.3d0
  do k=1,nk;do j=1,no;do i=1,ng
    ace%factors(i,j,k)=cmplx(sin(.017d0*i*j+k),cos(.013d0*i*(j+k)),8)
    t(i,j,k)=cmplx(cos(.019d0*i*j+k),sin(.023d0*i*(j+k)),8)
  enddo;enddo;enddo
  do k=1,nk
    ref(:,:,k)=-ace%dv*matmul(ace%factors(:,:,k),matmul(transpose(conjg(ace%factors(:,:,k))),t(:,:,k)))
  enddo
  call hse_ace_apply(ace,t,a,ierr)
  if(ierr/=0)stop 1
  error=sqrt(sum(abs(a-ref)**2)/sum(abs(ref)**2))
  print *,error
  if(.not.(error<=1d-11))stop 2
end program
