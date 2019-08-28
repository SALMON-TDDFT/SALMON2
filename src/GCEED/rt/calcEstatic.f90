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
subroutine calcEstatic(ng, info, sVh, srg_ng)
use salmon_communication, only: comm_summation
use scf_data
use new_world_sub
use structures, only: s_rgrid, s_orbital_parallel, s_scalar, s_sendrecv_grid
use sendrecv_grid, only: update_overlap_real8
implicit none
type(s_rgrid),intent(in)            :: ng
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

if(ng%is(1)==lg_sta(1))then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
    do ix=1,Ndh
      Vh_wk(ng%is(1)-ix,iy,iz) = Vh_wk(ng%is(1),iy,iz)
    end do
  end do
  end do
end if

if(ng%ie(1)==lg_end(1))then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
    do ix=1,Ndh
      Vh_wk(ng%ie(1)+ix,iy,iz) = Vh_wk(ng%ie(1),iy,iz)
    end do
  end do
  end do
end if

if(ng%is(2)==lg_sta(2))then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do ix=ng%is(1),ng%ie(1)
    do iy=1,Ndh
      Vh_wk(ix,ng%is(2)-iy,iz) = Vh_wk(ix,ng%is(2),iz)
    end do
  end do
  end do
end if

if(ng%ie(2)==lg_end(2))then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do ix=ng%is(1),ng%ie(1)
    do iy=1,Ndh
      Vh_wk(ix,ng%ie(2)+iy,iz) = Vh_wk(ix,ng%ie(2),iz)
    end do
  end do
  end do
end if

if(ng%is(3)==lg_sta(3))then
!$OMP parallel do private(iz,iy,ix)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    do iz=1,Ndh
      Vh_wk(ix,iy,ng%is(3)-iz) = Vh_wk(ix,iy,ng%is(3))
    end do
  end do
  end do
end if

if(ng%ie(3)==lg_end(3))then
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

