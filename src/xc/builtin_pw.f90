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
module builtin_pw
  implicit none
  
contains
  
  
  subroutine exc_cor_pw(nl, rho_s, exc, eexc, vexc)
    implicit none
    integer, intent(in) :: nl
    real(8), intent(in) :: rho_s(nl)
    real(8), intent(out) :: exc(nl), eexc(nl), vexc(nl)
    integer :: i
    real(8) :: trho, e_xc, de_xc_drho, de_xc_drhob

    ! call rho_j_tau(gs_rt, rho_s, tau_s, j_s, grho_s, lrho_s)
    ! rho_s=rho*0.5d0
    ! if(flag_nlcc)rho_s = rho_s + 0.5d0*rho_nlcc
    
#ifdef USE_OPENACC
!$acc kernels loop private(i,trho,e_xc,de_xc_drho)
#else
!$omp parallel do private(i,trho,e_xc,de_xc_drho)
#endif
    do i=1,NL
      trho=rho_s(i)
      call PWxc(trho,trho,e_xc,de_xc_drho,de_xc_drhob) ! de_xc_drhob is dummy
      exc(i)=e_xc
      Eexc(i)=e_xc*(2d0*trho)
      Vexc(i)=e_xc+(2d0*trho)*de_xc_drho
    enddo
#ifdef USE_OPENACC
!$acc end kernels
#endif
    return
  end subroutine exc_cor_pw
  
  SUBROUTINE PWxc(rhoa,rhob,exc,dexc_drhoa,dexc_drhob)
  !$acc routine seq
    implicit none
    real(8),parameter :: Pi=3.141592653589793d0
    real(8),intent(IN) :: rhoa,rhob
    real(8),intent(OUT) :: exc,dexc_drhoa,dexc_drhob
    real(8),parameter :: const = 0.75d0*(3d0/(2d0*pi))**(2d0/3d0)
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
    fzeta=((1d0+zeta)**(4d0/3d0)+(1d0-zeta)**(4d0/3d0)-2d0)/(2d0**(4d0/3d0)-2d0)
    dfzeta=4d0/3d0*((1d0+zeta)**(1d0/3d0)-(1d0-zeta)**(1d0/3d0))/(2d0**(4d0/3d0)-2d0)
    
!   Slater exchange
    exc = -const/rs                    !!! Warning: Spin-unpolarized (zeta=0)
    dexc_drhoa = exc/(3d0*(rhoa+rhob)) !!! Warning: Spin-unpolarized (zeta=0)
    dexc_drhob = dexc_drhoa            !!! Warning: Spin-unpolarized (zeta=0)

!   PW correlation [J. P. Perdew and Y. Wang., Phys. Rev. B 45, 13244 (1992)].
    Q0ec0=-2d0*Aec0*(1d0+a1ec0*rs)
    Q0ec1=-2d0*Aec1*(1d0+a1ec1*rs)
    Q0ar =-2d0*Aar *(1d0+a1ar *rs)
    Q1ec0=2d0*Aec0*(b1ec0*sqrs+b2ec0*rs+b3ec0*rs*sqrs+b4ec0*rs**(pec0+1d0))
    Q1ec1=2d0*Aec1*(b1ec1*sqrs+b2ec1*rs+b3ec1*rs*sqrs+b4ec1*rs**(pec1+1d0))
    Q1ar =2d0*Aar *(b1ar *sqrs+b2ar *rs+b3ar *rs*sqrs+b4ar *rs**(par +1d0))
    Qdec0=Aec0*(b1ec0/sqrs+2*b2ec0+3*b3ec0*sqrs+2d0*(pec0+1)*b4ec0*rs**pec0)
    Qdec1=Aec1*(b1ec1/sqrs+2*b2ec1+3*b3ec1*sqrs+2d0*(pec1+1)*b4ec1*rs**pec1)
    Qdar =Aar *(b1ar /sqrs+2*b2ar +3*b3ar *sqrs+2d0*(par +1)*b4ar *rs**par )
    Gec0=Q0ec0*log(1d0+1d0/Q1ec0)
    Gec1=Q0ec1*log(1d0+1d0/Q1ec1)
    Gar =Q0ar *log(1d0+1d0/Q1ar )
    dGec0=-2d0*Aec0*a1ec0*log(1d0+1d0/Q1ec0)-Q0ec0*Qdec0/(Q1ec0**2+Q1ec0)
    dGec1=-2d0*Aec1*a1ec1*log(1d0+1d0/Q1ec1)-Q0ec1*Qdec1/(Q1ec1**2+Q1ec1)
    dGar =-2d0*Aar *a1ar *log(1d0+1d0/Q1ar )-Q0ar *Qdar /(Q1ar **2+Q1ar )
    decdrs=dGec0*(1d0-fzeta*zeta**4)+dGec1*fzeta*zeta**4-dGar*fzeta/f2d*(1d0-zeta**4)
    decdzeta=dfzeta*(-Gar/f2d*(1d0-zeta**4)+(Gec1-Gec0)*zeta**4) &
  &         +4d0*zeta**3*(fzeta*(Gec1-Gec0)+Gar*fzeta/f2d)

    exc = exc + Gec0-Gar*fzeta/f2d*(1d0-zeta**4)+(Gec1-Gec0)*fzeta*zeta**4
    dexc_drhoa = dexc_drhoa + (-rs/3d0*decdrs-(zeta-1d0)*decdzeta)/(rhoa+rhob)
    dexc_drhob = dexc_drhob + (-rs/3d0*decdrs-(zeta+1d0)*decdzeta)/(rhoa+rhob)

    return
  end subroutine PWxc
  
end module builtin_pw
