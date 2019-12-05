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
subroutine calc_emfields(itt,nspin,rt,curr_in)
  use structures, only : s_rt
  use math_constants, only : pi
  use salmon_global, only : ispin
  use scf_data, only : dt
  use inputoutput, only: trans_longi
  implicit none
  integer,intent(in) :: itt
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

subroutine calc_Aext(Mit,itotNtime,rt)
!$ use omp_lib
use structures, only : s_rt
use em_field, only: calc_Ac_ext
use scf_data, only: dt
use math_constants, only: pi
implicit none
integer :: itt,Mit,itotNtime
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
