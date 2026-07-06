!
!  Copyright 2018-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!-----------------------------------------------------------------------------------------
module builtin_tbmbj
  use salmon_math, only: erfc_salmon
  implicit none
  
  real(8),parameter :: Pi=3.141592653589793d0

contains
  
  
  
  subroutine exc_cor_tbmbj(nl, rho, rho_s,  grho_s, lrho_s, tau_s, j_s, cval, eexc, vexc)
  !subroutine exc_cor_tbmbj(nl, rho, rho_s,  grho_s, lrho_s, tau_s, j_s, cval, eexc, vexc, Hxyz, aLxyz)
    implicit none
    integer, intent(in) :: nl
    real(8), intent(in) :: rho(nl), rho_s(nl)
    real(8), intent(in) :: grho_s(nl, 3), lrho_s(nl), tau_s(nl), j_s(nl, 3)
    real(8), intent(in) :: cval
    real(8), intent(out) :: eexc(nl), vexc(nl)
    !real(8), intent(in) :: Hxyz, aLxyz
    
    real(8),parameter :: alpha=-0.012d0,beta=1.023d0,gamma=0.80d0
    real(8),parameter :: rho_floor=1d-12 ! below this the mBJ evaluation is numerically unsafe
    real(8) :: c,tau_s_jrho,D_s_jrho,Q_s,rhs,x_s,b_s,Vx_BR,Vx_MBJ
    real(8) :: trho,rs,rhos,ec,dec_drhoa,dec_drhob
    integer :: i
    
    ! call rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)

    ! select case(functional)
    ! case('TBmBJ')
       if(cval > 1d-10) then
         c=cval ! use c-value given by input file
       else
         !c=sum(sqrt(grho_s(:,1)**2+grho_s(:,2)**2+grho_s(:,3)**2)/rho_s(:))*Hxyz/aLxyz
         ! nonpositive-density points (mixing transients) are excluded from the cell average
         c=sum(sqrt(grho_s(:,1)**2+grho_s(:,2)**2+grho_s(:,3)**2)/rho_s(:), mask=(rho_s(:) > 0d0)) / nl
         if(c >= 0d0 .and. c < huge(c)) then ! NaN/Inf-proof (contaminated transients)
           c=alpha+beta*sqrt(c)
         else
           c=1d0 ! fall back to the Becke-Johnson value
         endif
       endif
    ! case('BJ_PW')
    !    c=1d0
    ! end select
    
    do i=1,NL
      ! guard against nonpositive density (mixing transients); NaN also fails this test
      if(rho_s(i) > 0d0) then
        rhos=rho_s(i)
        tau_s_jrho=tau_s(i)-(j_s(i,1)**2+j_s(i,2)**2+j_s(i,3)**2)/rhos/2
        ! tau_s >= j_s^2/(2 rho_s) holds exactly (Cauchy-Schwarz); clamp negatives
        ! from floating-point cancellation at low density (would give sqrt(NaN))
        tau_s_jrho=max(tau_s_jrho,0d0)
        D_s_jrho=2*tau_s_jrho-0.25d0*(grho_s(i,1)**2+grho_s(i,2)**2+grho_s(i,3)**2)/rhos
        Q_s=(lrho_s(i)-2*gamma*D_s_jrho)/6d0
        rhs=2d0/3d0*pi**(2d0/3d0)*rhos**(5d0/3d0)/Q_s
        if(abs(rhs) > 0d0) then ! excludes rhs=0 and rhs=NaN (no bracket in BR_Newton)
          call BR_Newton(rhs,x_s)
          b_s=(x_s**3*exp(-x_s)/(8*pi*rhos))**(1d0/3d0)
          Vx_BR=-(1-exp(-x_s)-0.5d0*x_s*exp(-x_s))/b_s
          Vx_MBJ=c*Vx_BR+(3*c-2)/pi*sqrt(5d0/12d0)*sqrt(2*tau_s_jrho/rhos)
        else
          rhos=max(rhos,rho_floor) ! rhs=0 (underflow/huge Q_s) hung the old code; keep PWc in range
          Vx_MBJ=-(6d0*rhos/pi)**(1d0/3d0) ! LDA (Slater) exchange fallback
        endif
      else
        rhos=rho_floor
        Vx_MBJ=-(6d0*rhos/pi)**(1d0/3d0) ! LDA (Slater) exchange fallback
      endif

      ! total density guarded nonnegative (mixing transients); NaN also fails this test
      if(rho(i) > 0d0) then
        trho=rho(i)
      else
        trho=0d0
      endif
      rs=(3d0/(4*Pi*(trho+1d-10)))**(1d0/3d0)
      Eexc(i)=-.4582d0/rs
      Vexc(i)=Vx_MBJ
      call PWc(rhos,rhos,ec,dec_drhoa,dec_drhob)
      Vexc(i)=Vexc(i)+ec+trho*dec_drhoa
      Eexc(i)=Eexc(i)+ec
      Eexc(i)=Eexc(i)*(trho+1d-10)
    enddo
    return
  End Subroutine exc_cor_tbmbj



  SUBROUTINE BR_Newton(rhs,x_s)
    implicit none
    real(8),intent(IN) :: rhs
    real(8),intent(OUT) :: x_s
    integer,parameter :: maxiter=200 ! keep the search loops bounded for pathological rhs
    integer iter
    real(8) :: xmin,xmax,x,fx,dfx

  ! find xmax
    xmin=0d0
    x=1d0
    do iter=1,maxiter
      fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
      if(.not.(fx >= 0d0)) exit ! fx < 0 or fx is NaN
      x=x*2
    enddo
    xmax=x
  ! bi-section
    do iter=1,maxiter
      x=0.5d0*(xmin+xmax)
      fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
      if(fx < 0) then
        xmax=x
      else
        xmin=x
      endif
      if(.not.(xmax-xmin >= 1d-4)) exit ! converged (NaN-proof comparison)
    enddo
  ! Newton-Raphson
    do iter=1,5
      fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
      dfx=(1-2d0/3d0*x)*exp(-2d0/3d0*x)/rhs-1d0
      x=x-fx/dfx
    enddo
    x_s=x
    return
  end subroutine BR_Newton



  SUBROUTINE PWc(rhoa,rhob,ec,dec_drhoa,dec_drhob)
  implicit none
  real(8),parameter :: Pi=3.141592653589793d0
  real(8),intent(IN) :: rhoa,rhob
  real(8),intent(OUT) :: ec,dec_drhoa,dec_drhob
  real(8),parameter :: pec0=1d0,pec1=1d0,par=1d0
  real(8),parameter :: Aec0=0.031091d0,Aec1=0.015545d0,Aar=0.016887d0
  real(8),parameter :: a1ec0=0.21370d0,a1ec1=0.20548d0,a1ar=0.11125d0
  real(8),parameter :: b1ec0=7.5957d0,b1ec1=14.1189d0,b1ar=10.357d0
  real(8),parameter :: b2ec0=3.5876d0,b2ec1=6.1977d0,b2ar=3.6231d0
  real(8),parameter :: b3ec0=1.6382d0,b3ec1=3.3662d0,b3ar=0.88026d0
  real(8),parameter :: b4ec0=0.49294d0,b4ec1=0.62517d0,b4ar=0.49671d0
  real(8),parameter :: f2d=1.709921d0,eps=1d-15
  real(8) rs,sqrs,zeta,fzeta,dfzeta,Q0ec0,Q0ec1,Q0ar,Q1ec0,Q1ec1,Q1ar
  real(8) Qdec0,Qdec1,Qdar,Gec0,Gec1,Gar,dGec0,dGec1,dGar,decdrs,decdzeta

  rs=(3d0/(4d0*Pi*(rhoa+rhob)))**(1d0/3d0)
  sqrs=sqrt(rs)
  zeta=(rhoa-rhob)/(rhoa+rhob)
  if(abs(zeta-1d0) < eps) zeta=1d0-eps
  if(abs(zeta+1d0) < eps) zeta=-1d0+eps
  fzeta=((1+zeta)**(4d0/3d0)+(1-zeta)**(4d0/3d0)-2)/(2d0**(4d0/3d0)-2d0)
  dfzeta=4d0/3d0*((1+zeta)**(1d0/3d0)-(1-zeta)**(1d0/3d0))/(2d0**(4d0/3d0)-2d0)

  Q0ec0=-2*Aec0*(1+a1ec0*rs)
  Q0ec1=-2*Aec1*(1+a1ec1*rs)
  Q0ar =-2*Aar *(1+a1ar *rs)
  Q1ec0=2*Aec0*(b1ec0*sqrs+b2ec0*rs+b3ec0*rs*sqrs+b4ec0*rs**(pec0+1))
  Q1ec1=2*Aec1*(b1ec1*sqrs+b2ec1*rs+b3ec1*rs*sqrs+b4ec1*rs**(pec1+1))
  Q1ar =2*Aar *(b1ar *sqrs+b2ar *rs+b3ar *rs*sqrs+b4ar *rs**(par +1))
  Qdec0=Aec0*(b1ec0/sqrs+2*b2ec0+3*b3ec0*sqrs+2*(pec0+1)*b4ec0*rs**pec0)
  Qdec1=Aec1*(b1ec1/sqrs+2*b2ec1+3*b3ec1*sqrs+2*(pec1+1)*b4ec1*rs**pec1)
  Qdar =Aar *(b1ar /sqrs+2*b2ar +3*b3ar *sqrs+2*(par +1)*b4ar *rs**par )
  Gec0=Q0ec0*log(1+1d0/Q1ec0)
  Gec1=Q0ec1*log(1+1d0/Q1ec1)
  Gar =Q0ar *log(1+1d0/Q1ar )
  dGec0=-2*Aec0*a1ec0*log(1+1d0/Q1ec0)-Q0ec0*Qdec0/(Q1ec0**2+Q1ec0)
  dGec1=-2*Aec1*a1ec1*log(1+1d0/Q1ec1)-Q0ec1*Qdec1/(Q1ec1**2+Q1ec1)
  dGar =-2*Aar *a1ar *log(1+1d0/Q1ar )-Q0ar *Qdar /(Q1ar **2+Q1ar )

  ec=Gec0-Gar*fzeta/f2d*(1-zeta**4)+(Gec1-Gec0)*fzeta*zeta**4
  decdrs=dGec0*(1-fzeta*zeta**4)+dGec1*fzeta*zeta**4-dGar*fzeta/f2d*(1-zeta**4)
  decdzeta=dfzeta*(-Gar/f2d*(1-zeta**4)+(Gec1-Gec0)*zeta**4) &
&         +4*zeta**3*(fzeta*(Gec1-Gec0)+Gar*fzeta/f2d)

  dec_drhoa=(-rs/3*decdrs-(zeta-1)*decdzeta)/(rhoa+rhob)
  dec_drhob=(-rs/3*decdrs-(zeta+1)*decdzeta)/(rhoa+rhob)

  return
  end subroutine



end module builtin_tbmbj
