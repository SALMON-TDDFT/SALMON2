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
  use scf_data, only : curr,itt,iflag_indA,A_ind,dt,A_tot,A_ext,E_ext,E_ind,E_tot
  implicit none
  integer,intent(in) :: nspin
  real(8),intent(in) :: curr_in(3,nspin)

  curr(1:3,itt) = curr_in(1:3,1)
  if(ispin==1) curr(1:3,itt) = curr(1:3,itt) + curr_in(1:3,2)

  if(iflag_indA==1)then
    A_ind(:,itt+1)=2.d0*A_ind(:,itt)-A_ind(:,itt-1)-4.d0*Pi*curr(:,itt)*dt**2
  else if(iflag_indA==0)then
    A_ind(:,itt+1)=0.d0
  end if

  A_tot(:,itt+1)=A_ext(:,itt+1)+A_ind(:,itt+1)

  E_ext(:,itt)=-(A_ext(:,itt+1)-A_ext(:,itt-1))/(2.d0*dt)
  E_ind(:,itt)=-(A_ind(:,itt+1)-A_ind(:,itt-1))/(2.d0*dt)
  E_tot(:,itt)=-(A_tot(:,itt+1)-A_tot(:,itt-1))/(2.d0*dt)

end subroutine calc_emfields

subroutine calc_current_ion(lg,system,pp,j_ion)
  use structures
  use salmon_global, only: MI,Kion
  use scf_data, only: Hvol
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_dft_system) :: system
  type(s_pp_info) :: pp
  integer :: ia
  real(8) :: j_ion(3)

  !AY memo
  !current of ion: defined by positive charge-->minus sign
  !This is NOT matter current NOR electric current.... strange definition....
  !This is defined so as to the total electric current = -(curr + curr_ion)
  !Should change this ion current but if you change, 
  !please change all part in ARTED, multiscale ..... 
  j_ion(:)=0d0
  do ia=1,MI
     j_ion(:) = j_ion(:) - pp%Zps(Kion(ia)) * system%Velocity(:,ia)
  enddo
  j_ion(:) = j_ion(:)/(dble(lg%num(1)*lg%num(2)*lg%num(3))*Hvol)

end subroutine calc_current_ion
