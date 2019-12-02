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
subroutine calc_emfields(nspin,rt,curr_in)
  use structures, only : s_rt
  use math_constants, only : pi
  use salmon_global, only : ispin
  use scf_data, only : itt,dt
  use inputoutput, only: trans_longi
  implicit none
  type(s_rt),intent(inout) :: rt
  integer,intent(in) :: nspin
  real(8),intent(in) :: curr_in(3,2)  !curr_in(3,nspin)??

  rt%curr(1:3,itt) = curr_in(1:3,1)
  if(ispin==1) rt%curr(1:3,itt) = rt%curr(1:3,itt) + curr_in(1:3,2)  !<-- only if nspin==2??

  if(trans_longi=="lo")then
    rt%Ac_ind(:,itt+1)=2d0*rt%Ac_ind(:,itt) -rt%Ac_ind(:,itt-1) -4d0*Pi*rt%curr(:,itt)*dt**2
  else if(trans_longi=="tr")then
    rt%Ac_ind(:,itt+1)=0d0
  end if

  rt%Ac_tot(:,itt+1) = rt%Ac_ext(:,itt+1) + rt%Ac_ind(:,itt+1)

  rt%E_ext(:,itt) = -(rt%Ac_ext(:,itt+1) - rt%Ac_ext(:,itt-1))/(2d0*dt)
  rt%E_ind(:,itt) = -(rt%Ac_ind(:,itt+1) - rt%Ac_ind(:,itt-1))/(2d0*dt)
  rt%E_tot(:,itt) = -(rt%Ac_tot(:,itt+1) - rt%Ac_tot(:,itt-1))/(2d0*dt)

end subroutine calc_emfields

subroutine calc_Aext(Mit,rt)
!$ use omp_lib
use structures, only : s_rt
use em_field, only: calc_Ac_ext
use scf_data, only: dt,itotNtime
use math_constants, only: pi
implicit none
integer :: itt,Mit
  type(s_rt),intent(inout) :: rt
real(8) :: tt
do itt=Mit+1,itotNtime+1
   tt = dt*dble(itt)
   call calc_Ac_ext(tt,rt%Ac_ext(:,itt))
end do
return
end subroutine calc_Aext

subroutine init_A(Ntime,Mit,rt)
use structures, only : s_rt
use scf_data
use salmon_global, only: yn_restart
implicit none
type(s_rt),intent(inout) :: rt
integer :: Ntime,Mit
integer :: t_max

if(yn_restart /= 'y')then
  t_max = Ntime
else
  t_max = Ntime + Mit
end if

allocate( rt%curr( 3,0:t_max) )
allocate( rt%E_ext(3,0:t_max) )
allocate( rt%E_ind(3,0:t_max) )
allocate( rt%E_tot(3,0:t_max) )
allocate( rt%Ac_ext(3,0:t_max+1) )
allocate( rt%Ac_ind(3,0:t_max+1) )
allocate( rt%Ac_tot(3,0:t_max+1) )

rt%curr =0d0
rt%E_ext=0d0
rt%E_ind=0d0
rt%E_tot=0d0

rt%Ac_ext   =0d0
rt%Ac_ind   =0d0
rt%Ac_tot   =0d0

end subroutine init_A

subroutine calc_env_trigon(ipulse,tenv_trigon)
  !$ use omp_lib
  use inputoutput
  use scf_data
  use math_constants, only: pi

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
