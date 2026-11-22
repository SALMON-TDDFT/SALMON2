module hsex_frozen_kinetic
  implicit none
  private
  public :: frozen_kinetic_setup,frozen_kinetic_apply
contains
  subroutine frozen_kinetic_setup(n,cell,k,denominator,status)
    use hsex_screen_math, only: screen_geometry,screen_finite
    implicit none
    integer,intent(in) :: n(3)
    real(8),intent(in) :: cell(3,3),k(3)
    real(8),intent(out) :: denominator(:)
    integer,intent(out) :: status
    real(8) :: rec(3,3),vol,p(3)
    integer :: j,coord(3)
    status=1
    if(any(n<1).or.size(denominator)/=product(n))return
    if(.not.screen_finite(k))return
    call screen_geometry(cell,rec,vol,status)
    if(status/=0)return
    do j=0,product(n)-1
      coord=[mod(j,n(1)),mod(j/n(1),n(2)),j/(n(1)*n(2))]
      where(2*coord>=n)coord=coord-n
      p=matmul(rec,dble(coord))+k
      denominator(j+1)=1d0+.5d0*sum(p*p)
    enddo
    if(.not.screen_finite(denominator))status=1
  end subroutine
  subroutine frozen_kinetic_apply(n,denominator,x,y,status)
    use hsex_screen_math, only: screen_transform,screen_finite
    implicit none
    integer,intent(in) :: n(3)
    real(8),intent(in) :: denominator(:)
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    complex(8) :: transformed(1,size(x))
    status=1
    if(any(n<1).or.size(x)/=product(n).or.size(y)/=size(x).or.size(denominator)/=size(x))return
    if(.not.screen_finite(x).or..not.screen_finite(denominator))return
    if(any(denominator<=0d0))return
    ! Negative-sign forward FFT with 1/N; conjugation gives its inverse.
    transformed(1,:)=x;call screen_transform(transformed,n)
    transformed(1,:)=conjg(transformed(1,:)/denominator)
    call screen_transform(transformed,n)
    y=conjg(transformed(1,:))*product(n)
    if(.not.screen_finite(y))return
    status=0
  end subroutine
end module
