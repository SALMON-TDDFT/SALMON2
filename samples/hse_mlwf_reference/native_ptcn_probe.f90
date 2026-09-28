program probe
  use hse_ptcn_core
  implicit none
  complex(8) :: u(3,1,1),x(3,1,1)
  real(8) :: residual
  integer :: ierr,builds,apps,j
  u(:,1,1)=[(1d0,0d0),(.2d0,.1d0),(-.1d0,.3d0)]
  u=u/sqrt(sum(abs(u)**2))
  call hse_ptcn_solve(u,.1d0,1d0,action,precondition,norm,x,residual,builds,apps,ierr)
  if(ierr/=0)error stop 'PT-CN did not converge'
  do j=1,3;write(*,'(2es25.16)')real(x(j,1,1)),aimag(x(j,1,1));enddo
contains
  subroutine action(a,b,refresh)
    complex(8),intent(in) :: a(:,:,:)
    complex(8),intent(out) :: b(:,:,:)
    logical,intent(in) :: refresh
    b(1,1,1)=.2d0*a(1,1,1)+.1d0*a(2,1,1)
    b(2,1,1)=.1d0*a(1,1,1)+.7d0*a(2,1,1)
    b(3,1,1)=1.1d0*a(3,1,1)
    b=b+.05d0*abs(a)**2*a
  end subroutine
  subroutine precondition(a,b)
    complex(8),intent(in) :: a(:,:,:)
    complex(8),intent(out) :: b(:,:,:)
    b=a
  end subroutine
  real(8) function norm(a)
    complex(8),intent(in) :: a(:,:,:)
    norm=sqrt(sum(abs(a)**2))
  end function
end program
