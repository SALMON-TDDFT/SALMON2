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
!=======================================================================
!=======================================================================

subroutine set_vonf(mg,lg,Hgs,rt)
!$ use omp_lib
  use structures, only: s_rgrid,s_rt
  use salmon_communication, only: comm_is_root, comm_summation
  use salmon_global
  implicit none
  real(8),intent(in) :: Hgs(3)
  type(s_rgrid) :: mg,lg
  type(s_rt) :: rt
  integer :: i
  integer :: ix,iy,iz,iix,iiy,iiz
  integer :: max_icell
  real(8) :: rr

  if(iperiodic==0)then
    max_icell=0
  else if(iperiodic==3)then
    max_icell=2
  end if

  rt%vonf%f=0.d0
  rt%eonf(1)%f=0.d0
  rt%eonf(2)%f=0.d0
  rt%eonf(3)%f=0.d0

  do i=1,num_dipole_source
    do iiz=-max_icell,max_icell
    do iiy=-max_icell,max_icell
    do iix=-max_icell,max_icell
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rr=sqrt((lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1)))**2 &
               +(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2)))**2 &
               +(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3)))**2)
        if(rr>rad_dipole_source)then
          rt%vonf%f(ix,iy,iz)=rt%vonf%f(ix,iy,iz)  &
                        -(vec_dipole_source(1,i)*(lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1))) &
                         +vec_dipole_source(2,i)*(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2))) &
                         +vec_dipole_source(3,i)*(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3))))/rr**3
          rt%eonf(1)%f(ix,iy,iz)=rt%eonf(1)%f(ix,iy,iz)  &
                    +(3.d0*(vec_dipole_source(1,i)*(lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1)))/rr &
                           +vec_dipole_source(2,i)*(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2)))/rr &
                           +vec_dipole_source(3,i)*(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3)))/rr)* &
                              (lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1)))/rr &
                        -vec_dipole_source(1,i)) /rr**3
          rt%eonf(2)%f(ix,iy,iz)=rt%eonf(2)%f(ix,iy,iz)  &
                    +(3.d0*(vec_dipole_source(1,i)*(lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1)))/rr &
                           +vec_dipole_source(2,i)*(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2)))/rr &
                           +vec_dipole_source(3,i)*(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3)))/rr)* &
                              (lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2)))/rr &
                        -vec_dipole_source(2,i)) /rr**3
          rt%eonf(3)%f(ix,iy,iz)=rt%eonf(3)%f(ix,iy,iz)  &
                    +(3.d0*(vec_dipole_source(1,i)*(lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1)))/rr &
                           +vec_dipole_source(2,i)*(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2)))/rr &
                           +vec_dipole_source(3,i)*(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3)))/rr)* &
                              (lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3)))/rr &
                          -vec_dipole_source(3,i)) /rr**3
        else
          rt%vonf%f(ix,iy,iz)=rt%vonf%f(ix,iy,iz)  &
                        -(vec_dipole_source(1,i)*(lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1))) &
                         +vec_dipole_source(2,i)*(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2))) &
                         +vec_dipole_source(3,i)*(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3))))&
                                                                                                        /rad_dipole_source**3
          rt%eonf(1)%f(ix,iy,iz)=rt%eonf(1)%f(ix,iy,iz)  &
                         +vec_dipole_source(1,i)/rad_dipole_source**3
          rt%eonf(2)%f(ix,iy,iz)=rt%eonf(2)%f(ix,iy,iz)  &
                         +vec_dipole_source(2,i)/rad_dipole_source**3
          rt%eonf(3)%f(ix,iy,iz)=rt%eonf(3)%f(ix,iy,iz)  &
                         +vec_dipole_source(3,i)/rad_dipole_source**3
        end if
      end do
      end do
      end do
    end do
    end do
    end do
  end do

  return

end subroutine set_vonf

