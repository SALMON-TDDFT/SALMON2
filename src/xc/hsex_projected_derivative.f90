! A derivatives of the real-space box-projected nonlocal HSEX action.
module hsex_projected_derivative
  implicit none
  private
  public :: projected_derivative,rank1_second
contains
  subroutine rank1_second(indices,r,beta,norm,x,mass,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    integer,intent(in) :: indices(:)
    real(8),intent(in) :: r(:),norm
    complex(8),intent(in) :: beta(:),x(:)
    complex(8),intent(out) :: mass(:)
    integer,intent(out) :: status
    complex(8) :: p,rp,r2p
    integer :: j,m
    status=1;mass=0d0;m=size(indices)
    if(size(r)/=m.or.size(beta)/=m.or.size(mass)/=size(x))return
    if(any(indices<1).or.any(indices>size(x)))return
    if(.not.screen_finite(r).or..not.screen_finite(beta).or..not.screen_finite(x))return
    if(.not.screen_finite([norm]))return
    p=0d0;rp=0d0;r2p=0d0
    do j=1,m
      p=p+conjg(beta(j))*x(indices(j))
      rp=rp+r(j)*conjg(beta(j))*x(indices(j))
      r2p=r2p+r(j)**2*conjg(beta(j))*x(indices(j))
    enddo
    do j=1,m
      mass(indices(j))=mass(indices(j))-norm*beta(j)*(r(j)**2*p-2d0*r(j)*rp+r2p)
    enddo
    if(.not.screen_finite(mass))return
    status=0
  end subroutine

  subroutine projected_derivative(box,phase,r,sigma,first,hvol,x,velocity,status,second,mass)
    use hsex_screen_math, only: screen_finite
    implicit none
    integer,intent(in) :: box(:)
    complex(8),intent(in) :: phase(:),sigma(:,:),first(:,:),x(:)
    real(8),intent(in) :: r(:),hvol
    complex(8),intent(out) :: velocity(:)
    integer,intent(out) :: status
    complex(8),optional,intent(in) :: second(:,:)
    complex(8),optional,intent(out) :: mass(:)
    complex(8) :: p(size(sigma,1)),rp(size(sigma,1)),r2p(size(sigma,1)),hp(size(sigma,1)),dp(size(sigma,1))
    complex(8) :: ddp(size(sigma,1))
    integer :: ng,na,g
    status=1;velocity=0d0
    if(present(mass))mass=0d0
    ng=size(x);na=size(sigma,1)
    if(ng<1.or.na<1.or.size(box)/=ng.or.size(phase)/=ng.or.size(r)/=ng.or.size(velocity)/=ng)return
    if(size(sigma,2)/=na.or.any(shape(first)/=[na,na]).or.any(box<1).or.any(box>na))return
    if(.not.screen_finite(phase).or..not.screen_finite(r).or..not.screen_finite(x))return
    if(.not.screen_finite(sigma).or..not.screen_finite(first).or..not.screen_finite([hvol]))return
    if(hvol<=0d0.or.(present(second).neqv.present(mass)))return
    if(present(second))then
      if(any(shape(second)/=[na,na]).or.size(mass)/=ng)return
      if(.not.screen_finite(second))return
    endif
    p=0d0;rp=0d0;r2p=0d0
    do g=1,ng
      p(box(g))=p(box(g))+phase(g)*x(g)
      rp(box(g))=rp(box(g))+phase(g)*r(g)*x(g)
      if(present(second))r2p(box(g))=r2p(box(g))+phase(g)*r(g)**2*x(g)
    enddo
    hp=matmul(sigma,p)
    dp=matmul(first,p)+cmplx(0d0,1d0,8)*matmul(sigma,rp)
    if(present(second))ddp=matmul(second,p)+cmplx(0d0,2d0,8)*matmul(first,rp)-matmul(sigma,r2p)
    do g=1,ng
      ! The same 1/hvol as hsex_add is required for operator derivatives.
      velocity(g)=conjg(phase(g))*(dp(box(g))-cmplx(0d0,r(g),8)*hp(box(g)))/hvol
      if(present(mass))mass(g)=conjg(phase(g))*(ddp(box(g))- &
        cmplx(0d0,2d0*r(g),8)*dp(box(g))-r(g)**2*hp(box(g)))/hvol
    enddo
    if(.not.screen_finite(velocity))return
    if(present(mass))then
      if(.not.screen_finite(mass))return
    endif
    status=0
  end subroutine
end module
