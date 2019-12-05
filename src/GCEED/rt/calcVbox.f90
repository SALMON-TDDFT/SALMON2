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

SUBROUTINE calcVbox(mg,lg,itt_t,system,Vbox)
  use structures, only: s_rgrid, s_dft_system, s_scalar
  use salmon_communication, only: comm_is_root
  use misc_routines, only: get_wtime
  use em_field, only: calc_E_ext
  use inputoutput
  use scf_data
  
  implicit none
  
  type(s_rgrid),intent(in) :: mg,lg
  integer :: itt_t
  type(s_dft_system),intent(inout) :: system
  type(s_scalar),intent(inout)     :: Vbox
  integer :: ix,iy,iz,jj
  integer :: ix_sta_Vbox(3),ix_end_Vbox(3)
  integer :: ipulse
  real(8) :: env_trigon_1,env_trigon_2
  integer,parameter :: Nd = 4



  if(iperiodic==0)then
    if(yn_md=='y' .or. yn_out_rvf_rt=='y')then
      do jj=1,3
        if(lg%is(jj)==mg%is(jj))then
          ix_sta_Vbox(jj)=mg%is(jj)
        else
          ix_sta_Vbox(jj)=mg%is(jj)-Nd
        end if
        if(lg%ie(jj)==mg%ie(jj))then
          ix_end_Vbox(jj)=mg%ie(jj)
        else
          ix_end_Vbox(jj)=mg%ie(jj)+Nd
        end if
      end do
    else
      ix_sta_Vbox(1:3)=mg%is(1:3)
      ix_end_Vbox(1:3)=mg%ie(1:3)
    end if
  end if
 
  !$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    Vbox%f(ix,iy,iz)=0.d0
  end do
  end do
  end do

  system%vec_E_ext(1:3)=0.d0
   
  if(iperiodic==0)then
    if(ae_shape1 == 'impulse')then
      continue
    else
        if(dt*dble(itt_t) <= tw1)then
          ipulse=1
          call calc_E_ext(ipulse,dt*dble(itt_t),env_trigon_1,'y')
        !$OMP parallel do collapse(2) private(ix,iy,iz)
          do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
          do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
          do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
            Vbox%f(ix,iy,iz)=Vbox%f(ix,iy,iz)  &
                             +( epdir_re1(1)*lg%coordinate(ix,1)   &
                               +epdir_re1(2)*lg%coordinate(iy,2)   &
                               +epdir_re1(3)*lg%coordinate(iz,3) )*env_trigon_1  &
                             +( epdir_im1(1)*lg%coordinate(ix,1)   &
                               +epdir_im1(2)*lg%coordinate(iy,2)   &
                               +epdir_im1(3)*lg%coordinate(iz,3) )*env_trigon_1
          end do
          end do
          end do
          system%vec_E_ext(1:3)=system%vec_E_ext(1:3)+E_amplitude1*epdir_re1(1:3)*env_trigon_1   &
                                                       +E_amplitude1*epdir_im1(1:3)*env_trigon_1
        end if
        if(abs(dt*dble(itt_t)-0.5d0*tw1-t1_t2) < 0.5d0*tw2)then
          ipulse=2
          call calc_E_ext(ipulse,dt*dble(itt_t),env_trigon_2,'y')
          !$OMP parallel do collapse(2) private(ix,iy,iz)
          do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
          do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
          do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
            Vbox%f(ix,iy,iz)=Vbox%f(ix,iy,iz)   &
                             +( epdir_re2(1)*lg%coordinate(ix,1)   &
                               +epdir_re2(2)*lg%coordinate(iy,2)   &
                               +epdir_re2(3)*lg%coordinate(iz,3) )*env_trigon_2  &
                             +( epdir_im2(1)*lg%coordinate(ix,1)   &
                               +epdir_im2(2)*lg%coordinate(iy,2)   &
                               +epdir_im2(3)*lg%coordinate(iz,3) )*env_trigon_2
          end do
          end do
          end do
          system%vec_E_ext(1:3)=system%vec_E_ext(1:3)+E_amplitude2*epdir_re2(1:3)*env_trigon_2   &
                                                       +E_amplitude2*epdir_im2(1:3)*env_trigon_2
        end if
    end if
  end if
  system%vec_E(1:3)=system%vec_E_ext(1:3) 
  system%vec_Ac(1:3)=system%vec_Ac(1:3)-system%vec_E(1:3)*dt
  system%vec_Ac_ext(1:3)=system%vec_Ac(1:3) 
   
  if(num_dipole_source>=1)then
    if(dt*dble(itt_t) <= tw1)then
      ipulse=1
      call calc_E_ext(ipulse,dt*dble(itt_t),env_trigon_1,'n')
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        Vbox%f(ix,iy,iz)=Vbox%f(ix,iy,iz)+vonf_sd(ix,iy,iz)*env_trigon_1
      end do
      end do
      end do
    end if
    if(abs(dt*dble(itt_t)-0.5d0*tw1-t1_t2) < 0.5d0*tw2)then
      ipulse=2
      call calc_E_ext(ipulse,dt*dble(itt_t),env_trigon_2,'n')
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        Vbox%f(ix,iy,iz)=Vbox%f(ix,iy,iz)+vonf_sd(ix,iy,iz)*env_trigon_2
      end do
      end do
      end do
    end if
  end if

  return
  
END SUBROUTINE calcVbox
