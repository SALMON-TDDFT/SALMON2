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

SUBROUTINE calcVbox(itt_t)
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use misc_routines, only: get_wtime
  use inputoutput
  use scf_data
  
  implicit none
  integer :: itt_t
  integer :: ix,iy,iz,jj
  integer :: ix_sta_Vbox(3),ix_end_Vbox(3)
  integer :: ipulse
  real(8) :: env_trigon_1,env_trigon_2



  if(iperiodic==0)then
    if(alocal_laser=='y')then
      do jj=1,3
        if(mg_sta(jj)-Nd<ilasbound_sta(jj))then
          ix_sta_Vbox(jj)=mg_sta(jj)-Nd
        else
          ix_sta_Vbox(jj)=ilasbound_sta(jj)
        end if
        if(mg_end(jj)+Nd>ilasbound_sta(jj))then
          ix_end_Vbox(jj)=mg_end(jj)+Nd
        else
          ix_end_Vbox(jj)=ilasbound_end(jj)
        end if
      end do
    else
      if(icalcforce==1.or.iflag_md==1)then
        do jj=1,3
          if(lg_sta(jj)==mg_sta(jj))then
            ix_sta_Vbox(jj)=mg_sta(jj)
          else
            ix_sta_Vbox(jj)=mg_sta(jj)-Nd
          end if
          if(lg_end(jj)==mg_end(jj))then
            ix_end_Vbox(jj)=mg_end(jj)
          else
            ix_end_Vbox(jj)=mg_end(jj)+Nd
          end if
        end do
      else
        ix_sta_Vbox(1:3)=mg_sta(1:3)
        ix_end_Vbox(1:3)=mg_end(1:3)
      end if
    end if
  end if
 
  !$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    Vbox(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
  if(iperiodic==0)then
    if(ae_shape1 == 'impulse')then
      continue
    else
        if(dt*dble(itt_t) <= pulse_tw1)then
          ipulse=1
          call calc_env_trigon(ipulse,env_trigon_1)
        !$OMP parallel do collapse(2) private(ix,iy,iz)
          do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
          do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
          do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
            Vbox(ix,iy,iz)=Vbox(ix,iy,iz)+  &
                           amplitude1*(epdir_re1(1)*(gridcoo(ix,1)-rlaser_center(1))+   &
                                       epdir_re1(2)*(gridcoo(iy,2)-rlaser_center(2))+   &
                                       epdir_re1(3)*(gridcoo(iz,3)-rlaser_center(3)))*env_trigon_1  &
                          +amplitude1*(epdir_im1(1)*(gridcoo(ix,1)-rlaser_center(1))+   &
                                       epdir_im1(2)*(gridcoo(iy,2)-rlaser_center(2))+   &
                                       epdir_im1(3)*(gridcoo(iz,3)-rlaser_center(3)))*env_trigon_1
          end do
          end do
          end do
        end if
        if(abs(dt*dble(itt_t)-0.5d0*pulse_tw1-t1_t2) < 0.5d0*pulse_tw2)then
          ipulse=2
          call calc_env_trigon(ipulse,env_trigon_2)
          !$OMP parallel do collapse(2) private(ix,iy,iz)
          do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
          do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
          do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
            Vbox(ix,iy,iz)=Vbox(ix,iy,iz)   &
                          +amplitude2*(epdir_re2(1)*(gridcoo(ix,1)-rlaser_center(1))+   &
                                       epdir_re2(2)*(gridcoo(iy,2)-rlaser_center(2))+   &
                                       epdir_re2(3)*(gridcoo(iz,3)-rlaser_center(3)))*env_trigon_2  &
                          +amplitude2*(epdir_im2(1)*(gridcoo(ix,1)-rlaser_center(1))+   &
                                       epdir_im2(2)*(gridcoo(iy,2)-rlaser_center(2))+   &
                                       epdir_im2(3)*(gridcoo(iz,3)-rlaser_center(3)))*env_trigon_2
          end do
          end do
          end do
        end if
    end if
  end if
   
  if(nump>=1)then
    if(dt*dble(itt_t) <= pulse_tw1)then
      ipulse=1
      call calc_env_trigon(ipulse,env_trigon_1)
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        Vbox(ix,iy,iz)=Vbox(ix,iy,iz)+vonf_sd(ix,iy,iz)*env_trigon_1
      end do
      end do
      end do
    end if
    if(abs(dt*dble(itt_t)-0.5d0*pulse_tw1-t1_t2) < 0.5d0*pulse_tw2)then
      ipulse=2
      call calc_env_trigon(ipulse,env_trigon_2)
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        Vbox(ix,iy,iz)=Vbox(ix,iy,iz)+vonf_sd(ix,iy,iz)*env_trigon_2
      end do
      end do
      end do
    end if
  end if

  return
  
END SUBROUTINE calcVbox
