! Frozen auxiliary-Pauli core in a declared ionic-rest reference frame.
! Scalar reference N,T and full reference J are supplied by a complete V2 PP.
! The velocity-gauge phase q=-Ac is a model prescription, not Dirac CKED.
module nc_mgga23_core
  implicit none
  private
  public :: core23_fields,core23_energy_deda
contains
  pure subroutine core23_fields(rho,tau,jref,ac,n,t,j,ierr)
    use,intrinsic :: ieee_arithmetic,only: ieee_is_finite
    implicit none
    real(8),intent(in)::rho,tau,jref(4,3),ac(3)
    real(8),intent(out)::n(4),t(4),j(4,3)
    integer,intent(out)::ierr
    real(8)::q(3)
    integer::a,c
    ierr=0;n=0d0;t=0d0;j=0d0
    if(.not.ieee_is_finite(rho).or..not.ieee_is_finite(tau)) then
      ierr=1
      return
    endif
    do c=1,3
      if(.not.ieee_is_finite(ac(c))) ierr=1
      do a=1,4
        if(.not.ieee_is_finite(jref(a,c))) ierr=1
      enddo
    enddo
    if(ierr/=0) return
    if(rho<0d0.or.tau<0d0.or.maxval(abs(jref(1,:)))>0d0) then
      ierr=2
      return
    endif
    if(rho==0d0.and.maxval(abs(jref))>0d0) then
      ierr=2
      return
    endif
    q=-ac;n(1)=rho;t(1)=tau
    do c=1,3
      j(:,c)=jref(:,c)+q(c)*n
      t=t+2d0*q(c)*jref(:,c)+q(c)**2*n
    enddo
    do a=1,4
      if(.not.ieee_is_finite(t(a))) ierr=3
      do c=1,3
        if(.not.ieee_is_finite(j(a,c))) ierr=3
      enddo
    enddo
    if(ierr/=0) then
      n=0d0;t=0d0;j=0d0
    endif
  end subroutine

  pure subroutine core23_energy_deda(n,jref,ac,dt,dj,deda)
    implicit none
    real(8),intent(in)::n(4),jref(4,3),ac(3),dt(4),dj(4,3)
    real(8),intent(out)::deda(3)
    integer::c
    do c=1,3
      deda(c)=-2d0*sum(dt*(jref(:,c)-ac(c)*n))-sum(dj(:,c)*n)
    enddo
  end subroutine
end module
