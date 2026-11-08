program test_static_sternheimer
  use hsex_sternheimer, only: sternheimer_band
  implicit none
  complex(8) :: target(4,2), vertices(4,2), response(2,2), expected(2,2)
  complex(8) :: eigenvectors(4,4), overlap(4,2)
  real(8), parameter :: temperature=0.02d0
  real(8) :: levels(4), occupations(2), weight(4)
  integer :: i,j,status

  levels=[0d0,0.03d0,1d0,2d0]
  occupations=1d0/(1d0+exp(levels(:2)/temperature))
  eigenvectors=0d0
  eigenvectors(1,1)=1d0
  eigenvectors(2,2)=1d0
  eigenvectors(3:4,3)=[0.6d0,0.8d0]
  eigenvectors(3:4,4)=[-0.8d0,0.6d0]
  target=eigenvectors(:,:2)
  vertices(:,1)=[(1d0,0d0),(0d0,1d0),(0.3d0,0.4d0),(0.1d0,0.2d0)]
  vertices(:,2)=[(0d0,1d0),(0.5d0,0.1d0),(0.6d0,0d0),(0d0,0.4d0)]
  overlap=matmul(conjg(transpose(eigenvectors)),vertices)
  weight(1)=0.25d0/(4d0*temperature)
  weight(2)=0.25d0*(0.5d0-occupations(2))/levels(2)
  weight(3:4)=0.25d0/levels(3:4)
  do j=1,2
    do i=1,2
      expected(i,j)=sum(weight*conjg(overlap(:,i))*overlap(:,j))
    enddo
  enddo
  call sternheimer_band(action,target,levels(:2),occupations,0d0,0.5d0,temperature,vertices,response,status)
  if(status/=0)error stop 'Sternheimer solver rejected analytic fixture'
  if(maxval(abs(response-expected))>1d-12)error stop 'finite-temperature response differs from band sum'
  print *, 'STATIC_STERNHEIMER_MAX_ERROR',maxval(abs(response-expected))
  ! The explicit band sum must contain the empty-source orientation, too.
  call sternheimer_band(action,target,levels(:2),[1d0,0d0],levels(2),0d0,temperature,vertices,response,status)
  if(status/=0)error stop 'empty source was rejected'
  do j=1,2
    do i=1,2
      expected(i,j)=0.25d0/levels(2)*conjg(overlap(1,i))*overlap(1,j)
    enddo
  enddo
  if(maxval(abs(response-expected))>1d-12)error stop 'empty-source explicit transition is missing'
  print *, 'STATIC_STERNHEIMER_EMPTY_SOURCE_PASS'

contains
  subroutine action(x,y,ierr)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: ierr
    y=matmul(eigenvectors,levels*matmul(conjg(transpose(eigenvectors)),x))
    ierr=0
  end subroutine
end program
