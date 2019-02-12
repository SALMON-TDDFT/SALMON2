!
!  Copyright 2017 SALMON developers
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
subroutine taylor(mg,info,tzpsi_in,tzpsi_out,htpsi)
use structures, only: s_rgrid, s_wf_info
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
implicit none

type(s_rgrid),intent(in) :: mg
type(s_wf_info),intent(in) :: info
integer :: nn,ix,iy,iz
complex(8) :: tzpsi_in(mg%is_array(1):mg%ie_array(1),  &
                       mg%is_array(2):mg%ie_array(2),  &
                       mg%is_array(3):mg%ie_array(3),1:info%numo,info%ik_s:info%ik_e)
complex(8) :: htpsi(mg%is_array(1):mg%ie_array(1),  &
                    mg%is_array(1):mg%ie_array(1),  &
                    mg%is_array(3):mg%ie_array(3),1:info%numo,info%ik_s:info%ik_e)
complex(8) :: tzpsi_out(mg%is_array(1):mg%ie_array(1),  &
                        mg%is_array(2):mg%ie_array(2),  &
                        mg%is_array(3):mg%ie_array(3),1:info%numo,info%ik_s:info%ik_e)

iwk_size=2
call make_iwksta_iwkend

if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rhobox(ix,iy,iz) = 0.d0
  end do
  end do
  end do
else if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rhobox_s(ix,iy,iz,1) = 0.d0
    rhobox_s(ix,iy,iz,2) = 0.d0
  end do
  end do
  end do
end if

if(iperiodic==0.and.ihpsieff==1)then
  if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
    end do
    end do
    end do
  else if(ilsda==1)then 
!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
      Vlocal2(ix,iy,iz,2)=Vlocal(ix,iy,iz,2)+Vbox(ix,iy,iz)
    end do
    end do
    end do
  end if
end if


  do nn=1,N_hamil
    if(iperiodic==0.and.ihpsieff==1)then
      if(mod(nn,2)==1)then
        call hpsi_groupob(tzpsi_in,htpsi,tzpsi_out,Vlocal2,nn,1)
      else
        call hpsi_groupob(htpsi,tzpsi_in,tzpsi_out,Vlocal2,nn,1)
      end if
    else
      if(mod(nn,2)==1)then
        call hpsi_groupob(tzpsi_in,htpsi,tzpsi_out,Vlocal,nn,1)
      else
        call hpsi_groupob(htpsi,tzpsi_in,tzpsi_out,Vlocal,nn,1)
      end if
    end if
  end do


end subroutine taylor

