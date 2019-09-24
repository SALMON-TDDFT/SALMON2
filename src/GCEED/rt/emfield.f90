!
!  Copyright 2019 SALMON developers
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
subroutine calc_emfields(nspin,curr_in)
  use math_constants, only : pi
  use salmon_global, only : ispin
  use scf_data, only : curr,itt,A_ind,dt,A_tot,A_ext,E_ext,E_ind,E_tot
  use inputoutput, only: trans_longi
  implicit none
  integer,intent(in) :: nspin
  real(8),intent(in) :: curr_in(3,nspin)

  curr(1:3,itt) = curr_in(1:3,1)
  if(ispin==1) curr(1:3,itt) = curr(1:3,itt) + curr_in(1:3,2)

  if(trans_longi=="lo")then
    A_ind(:,itt+1)=2.d0*A_ind(:,itt)-A_ind(:,itt-1)-4.d0*Pi*curr(:,itt)*dt**2
  else if(trans_longi=="tr")then
    A_ind(:,itt+1)=0.d0
  end if

  A_tot(:,itt+1)=A_ext(:,itt+1)+A_ind(:,itt+1)

  E_ext(:,itt)=-(A_ext(:,itt+1)-A_ext(:,itt-1))/(2.d0*dt)
  E_ind(:,itt)=-(A_ind(:,itt+1)-A_ind(:,itt-1))/(2.d0*dt)
  E_tot(:,itt)=-(A_tot(:,itt+1)-A_tot(:,itt-1))/(2.d0*dt)

end subroutine calc_emfields

subroutine calcAext
!$ use omp_lib
use em_field, only: calc_Ac_ext
use scf_data, only: dt,Miter_rt,itotNtime,A_ext
implicit none
integer :: itt
real(8) :: tt
do itt=Miter_rt+1,itotNtime+1
  tt = dt*dble(itt)
  call calc_Ac_ext(tt,A_ext(:,itt))
end do
return
end subroutine calcAext

subroutine calc_vecAc(vec_Ac,imode)
use scf_data
!$ use omp_lib
implicit none
real(8) :: vec_Ac(3)
integer :: imode
!
complex(8) :: vecA(3)

if(iSCFRT==1)then
  vecA=0.d0
else if(iSCFRT==2)then
  if(imode==1)then
    vecA(:)=A_ind(:,itt)
    if(epdir_re1(1)==1.d0)then
      vecA(1)=vecA(1)+A_ext(1,itt)
    else if(epdir_re1(2)==1.d0)then
      vecA(2)=vecA(2)+A_ext(2,itt)
    else if(epdir_re1(3)==1.d0)then
      vecA(3)=vecA(3)+A_ext(3,itt)
    end if
  else if(imode==2.or.imode==3)then
    vecA=0.d0
  else if(imode==4)then
    vecA(:)=A_ind(:,itt+1)
    if(epdir_re1(1)==1.d0)then
      vecA(1)=vecA(1)+A_ext(1,itt+1)
    else if(epdir_re1(2)==1.d0)then
      vecA(2)=vecA(2)+A_ext(2,itt+1)
    else if(epdir_re1(3)==1.d0)then
      vecA(3)=vecA(3)+A_ext(3,itt+1)
    end if
  end if
end if
vec_Ac = vecA

end subroutine calc_vecAc

subroutine initA(Ntime)
use scf_data
implicit none
integer :: Ntime
integer :: t_max

if(IC_rt==0)then
  t_max=Ntime
else
  t_max=Ntime+Miter_rt
end if

allocate( curr(3,0:t_max) )
allocate( curr_ion(3,0:t_max) )
allocate( sumcurr(3,0:t_max) )
allocate( A_ext(3,0:t_max+1) )
allocate( A_ind(3,0:t_max+1) )
allocate( A_tot(3,0:t_max+1) )
allocate( E_ext(3,0:t_max) )
allocate( E_ind(3,0:t_max) )
allocate( E_tot(3,0:t_max) )
allocate( rE_ind(3,0:t_max) )
curr=0.d0
curr_ion=0.d0
sumcurr=0.d0
A_ext=0.d0
A_ind=0.d0
A_tot=0.d0
E_ind=0.d0
rE_ind=0.d0

E_ext=0.d0
E_ind(:,0)=0.d0
E_tot(:,0)=0.d0

end subroutine initA

subroutine calc_env_trigon(ipulse,tenv_trigon)
  !$ use omp_lib
  use inputoutput
  use scf_data

  implicit none
  integer,intent(in)  :: ipulse  ! 1: first pulse, 2: second pulse
  real(8),intent(out) :: tenv_trigon

  real(8) :: alpha,beta
  real(8) :: theta1,theta2,theta1_0,theta2_0,tenv_trigon_0
  character(16) :: tae_shape

  if(ipulse==1)then
    tae_shape=ae_shape1
    ! cos(theta1)**2
    theta1=Pi/tw1*(dble(itt)*dt-0.5d0*tw1)
    alpha=Pi/tw1
    ! cos(theta2)
    theta2=omega1*(dble(itt)*dt-0.5d0*tw1)+phi_cep1*2d0*pi
    beta=omega1
    ! for iperiodic=3 .and. ae_shape1='Ecos2'
    theta1_0=Pi/tw1*(0-0.5d0*tw1)
    theta2_0=omega1*(0-0.5d0*tw1)+phi_cep1*2d0*pi
  else if(ipulse==2)then
    tae_shape=ae_shape2
    ! cos(theta1)**2
    theta1=Pi/tw2*(dble(itt)*dt-t1_t2-0.5d0*tw1)
    alpha=Pi/tw2
    ! cos(theta2)
    theta2=omega2*(dble(itt)*dt-t1_t2-0.5d0*tw1)+phi_cep2*2d0*pi
    beta=omega2
    ! for iperiodic=3 .and. ae_shape2='Ecos2'
    theta1_0=Pi/tw2*(0-0.5d0*tw1)
    theta2_0=omega2*(0-0.5d0*tw1)+phi_cep2*2d0*pi
  end if

  select case(iperiodic)
  case(0)
    if(tae_shape=='Ecos2')then
      tenv_trigon=cos(theta1)**2*cos(theta2)
    else if(tae_shape=='Acos2')then
      tenv_trigon=-(-alpha*sin(2.d0*theta1)*cos(theta2)   &
                    -beta*cos(theta1)**2*sin(theta2))/beta
    end if
  case(3)
    if(tae_shape=='Ecos2')then
      tenv_trigon_0=sin(theta2_0)/(2.d0*beta)   &
                   +sin(2.d0*theta1_0+theta2_0)/(4.d0*(2.d0*alpha+beta))   &
                   +sin(2.d0*theta1_0-theta2_0)/(4.d0*(2.d0*alpha-beta))
      tenv_trigon  =-(sin(theta2)/(2.d0*beta)   &
                   +sin(2.d0*theta1+theta2)/(4.d0*(2.d0*alpha+beta))   &
                   +sin(2.d0*theta1-theta2)/(4.d0*(2.d0*alpha-beta))   &
                   -tenv_trigon_0)
    else if(tae_shape=='Acos2')then
      tenv_trigon=(1/beta)*cos(theta1)**2*cos(theta2)
    end if
  end select

  return

end subroutine calc_env_trigon

subroutine calcEstatic(lg, ng, info, sVh, srg_ng)
use salmon_communication, only: comm_summation
use scf_data
use new_world_sub
use structures, only: s_rgrid, s_orbital_parallel, s_scalar, s_sendrecv_grid
use sendrecv_grid, only: update_overlap_real8
implicit none
type(s_rgrid),intent(in)            :: lg, ng
type(s_orbital_parallel),intent(in) :: info
type(s_scalar),intent(in) :: sVh
type(s_sendrecv_grid),intent(inout) :: srg_ng

integer :: ist,ix,iy,iz
real(8) :: Vh_wk(ng%is(1)-Ndh:ng%ie(1)+Ndh,   &
                 ng%is(2)-Ndh:ng%ie(2)+Ndh,   &
                 ng%is(3)-Ndh:ng%ie(3)+Ndh)
real(8) :: Ex_static2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
real(8) :: Ey_static2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
real(8) :: Ez_static2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))


!$OMP parallel do private(iz,iy,ix)
do iz=ng%is(3)-Ndh,ng%ie(3)+Ndh
do iy=ng%is(2)-Ndh,ng%ie(2)+Ndh
do ix=ng%is(1)-Ndh,ng%ie(1)+Ndh
  Vh_wk(ix,iy,iz)=0.d0
end do
end do
end do

!$OMP parallel do private(iz,iy,ix)
do iz=ng%is(3),ng%ie(3)
do iy=ng%is(2),ng%ie(2)
do ix=ng%is(1),ng%ie(1)
  Vh_wk(ix,iy,iz) = sVh%f(ix,iy,iz)
end do
end do
end do

call update_overlap_real8(srg_ng, ng, Vh_wk)

if(ng%is(1)==lg%is(1))then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
    do ix=1,Ndh
      Vh_wk(ng%is(1)-ix,iy,iz) = Vh_wk(ng%is(1),iy,iz)
    end do
  end do
  end do
end if

if(ng%ie(1)==lg%ie(1))then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
    do ix=1,Ndh
      Vh_wk(ng%ie(1)+ix,iy,iz) = Vh_wk(ng%ie(1),iy,iz)
    end do
  end do
  end do
end if

if(ng%is(2)==lg%is(2))then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do ix=ng%is(1),ng%ie(1)
    do iy=1,Ndh
      Vh_wk(ix,ng%is(2)-iy,iz) = Vh_wk(ix,ng%is(2),iz)
    end do
  end do
  end do
end if

if(ng%ie(2)==lg%ie(2))then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do ix=ng%is(1),ng%ie(1)
    do iy=1,Ndh
      Vh_wk(ix,ng%ie(2)+iy,iz) = Vh_wk(ix,ng%ie(2),iz)
    end do
  end do
  end do
end if

if(ng%is(3)==lg%is(3))then
!$OMP parallel do private(iz,iy,ix)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    do iz=1,Ndh
      Vh_wk(ix,iy,ng%is(3)-iz) = Vh_wk(ix,iy,ng%is(3))
    end do
  end do
  end do
end if

if(ng%ie(3)==lg%ie(3))then
!$OMP parallel do private(iz,iy,ix)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    do iz=1,Ndh
      Vh_wk(ix,iy,ng%ie(3)+iz) = Vh_wk(ix,iy,ng%ie(3))
    end do
  end do
  end do
end if

!$OMP parallel do private(iz,iy,ix)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  Ex_static2(ix,iy,iz)=0.d0
  Ey_static2(ix,iy,iz)=0.d0
  Ez_static2(ix,iy,iz)=0.d0
end do
end do
end do

!$OMP parallel do private(iz,iy,ix,ist)
do iz=ng%is(3),ng%ie(3)
do iy=ng%is(2),ng%ie(2)
do ix=ng%is(1),ng%ie(1)
  do ist=1,Ndh
    Ex_static2(ix,iy,iz)=Ex_static2(ix,iy,iz)               &
       -(bNmat(ist,Ndh)*(Vh_wk(ix+ist,iy,iz)-Vh_wk(ix-ist,iy,iz)))/Hgs(1)
    Ey_static2(ix,iy,iz)=Ey_static2(ix,iy,iz)               &
       -(bNmat(ist,Ndh)*(Vh_wk(ix,iy+ist,iz)-Vh_wk(ix,iy-ist,iz)))/Hgs(2)
    Ez_static2(ix,iy,iz)=Ez_static2(ix,iy,iz)               &
       -(bNmat(ist,Ndh)*(Vh_wk(ix,iy,iz+ist)-Vh_wk(ix,iy,iz-ist)))/Hgs(3)
  end do
end do
end do
end do

call comm_summation(Ex_static2,Ex_static,mg_num(1)*mg_num(2)*mg_num(3),info%icomm_ko)
call comm_summation(Ey_static2,Ey_static,mg_num(1)*mg_num(2)*mg_num(3),info%icomm_ko)
call comm_summation(Ez_static2,Ez_static,mg_num(1)*mg_num(2)*mg_num(3),info%icomm_ko)

end subroutine calcEstatic
