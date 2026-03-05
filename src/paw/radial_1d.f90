module radial_1d

  implicit none

  private
  public :: calc_hartree_radial_1d

contains

  subroutine calc_hartree_radial_1d( rad, rab, rho, Vh, Eh )
    implicit none
    real(8),intent(in) :: rad(:), rab(:)
    real(8),intent(in) :: rho(:)
    real(8),intent(out) :: Vh(:)
    real(8),intent(out) :: Eh
    integer :: nrad,i,j
    real(8) :: pi4,sum0,sum1

    Vh=0.0d0
    Eh=0.0d0
    nrad=size(rad)
    pi4=4.0d0*acos(-1.0d0)

    do i=1,nrad

       sum0=0.5d0*rad(1)**2*rho(1)*rab(1) &
           +0.5d0*rad(i)**2*rho(i)*rab(i)
       do j=2,i-1
          sum0=sum0+rad(j)**2*rho(j)*rab(j)
       end do
       sum0=sum0/rad(i)*pi4

       sum1=0.5d0*rad(i)*rho(i)*rab(i) &
           +0.5d0*rad(nrad)*rho(nrad)*rab(nrad)
       do j=i+1,nrad-1
          sum1=sum1+rad(j)*rho(j)*rab(j)
       end do
       sum1=sum1*pi4

       Vh(i) = sum1

       write(20,'(1x,4g20.10)') rad(i),sum0+sum1,sum0,sum1

    end do !i

  end subroutine calc_hartree_radial_1d

end module radial_1d
