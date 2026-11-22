! Causal finite-delay spectral filters. Times are integral areas, in atomic units.
module hsex_frozen_window
  implicit none
  private
  public :: window_filter,window_filter_matrix,window_moments
  integer,parameter :: quadrature_order=64
contains
  subroutine window_moments(kind,area,first,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    integer,intent(in) :: kind
    real(8),intent(in) :: area
    real(8),intent(out) :: first
    integer,intent(out) :: status
    status=1;first=0d0
    if(kind<1.or.kind>2.or..not.screen_finite([area]))return
    if(area<sqrt(tiny(1d0)).or.area>sqrt(huge(1d0)))return
    first=area*area
    if(kind==1)first=first*(2d0/acos(-1d0))
    if(.not.screen_finite([first]))return
    status=0
  end subroutine

  subroutine legendre(nodes,weights,status)
    implicit none
    real(8),intent(out) :: nodes(quadrature_order),weights(quadrature_order)
    integer,intent(out) :: status
    real(8) :: z,previous,p0,p1,p2,derivative,change,weight
    integer :: i,j,iteration
    status=1
    do i=1,quadrature_order/2
      z=cos(acos(-1d0)*(i-.25d0)/(quadrature_order+.5d0))
      do iteration=1,100
        p0=1d0;p1=z
        do j=2,quadrature_order
          p2=((2*j-1)*z*p1-(j-1)*p0)/j;p0=p1;p1=p2
        enddo
        derivative=quadrature_order*(z*p1-p0)/(z*z-1d0)
        previous=z;change=p1/derivative;z=z-change
        if(abs(z-previous)<=4d0*epsilon(z))exit
      enddo
      if(iteration>100)return
      ! Re-evaluate the derivative at the converged root.
      p0=1d0;p1=z
      do j=2,quadrature_order
        p2=((2*j-1)*z*p1-(j-1)*p0)/j;p0=p1;p1=p2
      enddo
      derivative=quadrature_order*(z*p1-p0)/(z*z-1d0)
      weight=1d0/((1d0-z*z)*derivative*derivative)
      nodes(i)=.5d0*(1d0-z);nodes(quadrature_order+1-i)=.5d0*(1d0+z)
      weights(i)=weight;weights(quadrature_order+1-i)=weight
    enddo
    status=0
  end subroutine

  subroutine evaluate(kind,area,delta,nodes,weights,value,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    integer,intent(in) :: kind
    real(8),intent(in) :: area,delta,nodes(:),weights(:)
    complex(8),intent(out) :: value
    integer,intent(out) :: status
    real(8) :: x,dawson,total,term,previous,inverse,imaginary
    integer :: j
    value=0d0;status=1
    if(area>1d0)then
      if(abs(delta)>huge(1d0)/area)return
    endif
    if(kind==2)then
      value=1d0/(cmplx(0d0,1d0/area,8)-delta)
    else
      x=delta*area/sqrt(acos(-1d0))
      if(abs(x)<8d0)then
        dawson=0d0
        do j=1,size(nodes)
          dawson=dawson+weights(j)*exp(-x*x*(1d0-nodes(j)*nodes(j)))
        enddo
        dawson=x*dawson
        value=cmplx(-2d0*area/sqrt(acos(-1d0))*dawson,-area*exp(-x*x),8)
      else
        inverse=0d0
        if(abs(x)<1d8)inverse=.5d0/x/x
        total=1d0;term=1d0
        do j=1,100
          previous=term;term=term*(2*j-1)*inverse
          if(abs(term)>abs(previous))exit
          total=total+term
          if(abs(term)<=epsilon(total)*abs(total))exit
        enddo
        imaginary=0d0
        if(abs(x)<=sqrt(-log(tiny(1d0))))imaginary=-area*exp(-x*x)
        value=cmplx(-total/delta,imaginary,8)
      endif
    endif
    if(.not.screen_finite([value]))return
    status=0
  end subroutine

  subroutine window_filter(kind,area,delta,value,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    integer,intent(in) :: kind
    real(8),intent(in) :: area,delta
    complex(8),intent(out) :: value
    integer,intent(out) :: status
    real(8) :: nodes(quadrature_order),weights(quadrature_order),first
    value=0d0
    call window_moments(kind,area,first,status)
    if(status/=0)return
    status=1
    if(.not.screen_finite([delta]))return
    nodes=0d0;weights=0d0
    if(kind==1)then
      call legendre(nodes,weights,status)
      if(status/=0)return
    endif
    call evaluate(kind,area,delta,nodes,weights,value,status)
  end subroutine

  subroutine window_filter_matrix(kind,area,left,right,values,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    integer,intent(in) :: kind
    real(8),intent(in) :: area,left(:),right(:)
    complex(8),intent(out) :: values(:,:)
    integer,intent(out) :: status
    real(8) :: nodes(quadrature_order),weights(quadrature_order),first,delta
    integer :: i,j,code,invalid
    values=0d0
    call window_moments(kind,area,first,status)
    if(status/=0)return
    status=1
    if(size(values,1)/=size(left).or.size(values,2)/=size(right))return
    if(.not.screen_finite(left).or..not.screen_finite(right))return
    nodes=0d0;weights=0d0
    if(kind==1)then
      call legendre(nodes,weights,status)
      if(status/=0)return
    endif
    invalid=0
    !$omp parallel do collapse(2) default(shared) private(i,j,delta,code) reduction(max:invalid) schedule(static)
    do j=1,size(right)
      do i=1,size(left)
        delta=left(i)-right(j);code=1
        if(screen_finite([delta]))call evaluate(kind,area,delta,nodes,weights,values(i,j),code)
        if(code/=0)invalid=1
      enddo
    enddo
    !$omp end parallel do
    status=invalid
  end subroutine
end module
