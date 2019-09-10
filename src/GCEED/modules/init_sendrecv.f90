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
MODULE init_sendrecv_sub

implicit none 

CONTAINS

!=======================================================================
!=======================================================================

subroutine create_sendrecv_neig_mg(neig_mg, ob_para_info, iperiodic)
  use structures, only: s_orbital_parallel
  use sendrecv_grid, only: init_sendrecv_grid
  use salmon_communication, only: comm_proc_null
  use salmon_global, only: nproc_domain_orbital
  implicit none
  integer, intent(out) :: neig_mg(1:2, 1:3)
  type(s_orbital_parallel), intent(in) :: ob_para_info
  integer, intent(in) :: iperiodic
  !
  integer :: imr(3)
  integer :: nproc_d_o(3)

  imr = ob_para_info%imr
  nproc_d_o = nproc_domain_orbital

  neig_mg(1, 1)=ob_para_info%id_r+1
  neig_mg(2, 1)=ob_para_info%id_r-1
  select case(iperiodic)
  case(0)
    if(imr(1)==nproc_d_o(1)-1) neig_mg(1, 1)=comm_proc_null
    if(imr(1)==0) neig_mg(2, 1)=comm_proc_null
  case(3)
    if(imr(1)==nproc_d_o(1)-1) then
      neig_mg(1, 1)=ob_para_info%id_r-(nproc_d_o(1)-1)
    end if
    if(imr(1)==0) then
      neig_mg(2, 1)=ob_para_info%id_r+(nproc_d_o(1)-1)
    end if
  end select
  
  neig_mg(1, 2)=ob_para_info%id_r+nproc_d_o(1)
  neig_mg(2, 2)=ob_para_info%id_r-nproc_d_o(1)
  select case(iperiodic)
  case(0)
    if(imr(2)==nproc_d_o(2)-1) neig_mg(1, 2)=comm_proc_null
    if(imr(2)==0) neig_mg(2, 2)=comm_proc_null
  case(3)
    if(imr(2)==nproc_d_o(2)-1) then
      neig_mg(1, 2)=ob_para_info%id_r-(nproc_d_o(2)-1)*nproc_d_o(1)
    end if
    if(imr(2)==0) then
      neig_mg(2, 2)=ob_para_info%id_r+(nproc_d_o(2)-1)*nproc_d_o(1)
    end if
  end select
  
  neig_mg(1, 3)=ob_para_info%id_r+nproc_d_o(1)*nproc_d_o(2)
  neig_mg(2, 3)=ob_para_info%id_r-nproc_d_o(1)*nproc_d_o(2)
  select case(iperiodic)
  case(0)
    if(imr(3)==nproc_d_o(3)-1) neig_mg(1, 3)=comm_proc_null
    if(imr(3)==0) neig_mg(2, 3)=comm_proc_null
  case(3)
    if(imr(3)==nproc_d_o(3)-1) then
      neig_mg(1, 3)=ob_para_info%id_r-(nproc_d_o(3)-1)*nproc_d_o(1)*nproc_d_o(2)
    end if
    if(imr(3)==0) then
      neig_mg(2, 3)=ob_para_info%id_r+(nproc_d_o(3)-1)*nproc_d_o(1)*nproc_d_o(2)
    end if
  end select
  
  return
end subroutine create_sendrecv_neig_mg


subroutine create_sendrecv_neig_ng(neig_ng, info_field, iperiodic)
  use structures, only: s_field_parallel
  use sendrecv_grid, only: init_sendrecv_grid
  use salmon_communication, only: comm_proc_null
  use salmon_global, only: process_allocation,nproc_domain_orbital,nproc_domain_general,nproc_ob,nproc_k
  implicit none
  integer, intent(out) :: neig_ng(1:2, 1:3)
  type(s_field_parallel), intent(in) :: info_field
  integer, intent(in) :: iperiodic
  !
  integer :: myrank,imr(3),imrs(3)
  integer :: nproc_d_o(3),nproc_d_g(3),nproc_d_o_mul,nproc_d_g_dm(3),nproc_d_g_mul_dm

  myrank = info_field%id_all
  imr = info_field%imr
  imrs = info_field%imrs

  nproc_d_o = nproc_domain_orbital
  nproc_d_g = nproc_domain_general
  nproc_d_o_mul = nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)
  nproc_d_g_dm = nproc_d_g/nproc_d_o
  nproc_d_g_mul_dm = nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)

  if(process_allocation=='orbital_sequential')then
    if(imr(1)==nproc_d_o(1)-1.and.imrs(1)==nproc_d_g_dm(1)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 1)=comm_proc_null
      case(3)
        neig_ng(1, 1)=myrank+nproc_d_g_mul_dm-nproc_d_g_dm(1)+1-nproc_d_g_mul_dm*nproc_d_o(1)
      end select
    else if(imrs(1)==nproc_d_g_dm(1)-1) then
      neig_ng(1, 1)=myrank+nproc_d_g_mul_dm-nproc_d_g_dm(1)+1
    else
      neig_ng(1, 1)=myrank+1
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(1)==nproc_d_o(1)-1.and.imrs(1)==nproc_d_g_dm(1)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 1)=comm_proc_null
      case(3)
        neig_ng(1, 1)=myrank+nproc_d_o_mul+1-nproc_d_o_mul*nproc_d_g_dm(1)-nproc_d_o(1)
      end select
    else if(imrs(1)==nproc_d_g_dm(1)-1) then
      neig_ng(1, 1)=myrank+nproc_d_o_mul+1-nproc_d_o_mul*nproc_d_g_dm(1)
    else
      neig_ng(1, 1)=myrank+nproc_d_o_mul
    end if
  end if
  
  if(process_allocation=='orbital_sequential')then
    if(imr(1)==0.and.imrs(1)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 1)=comm_proc_null
      case(3)
        neig_ng(2, 1)=myrank-nproc_d_g_mul_dm+nproc_d_g_dm(1)-1+nproc_d_g_mul_dm*nproc_d_o(1)
      end select
    else if(imrs(1)==0) then
      neig_ng(2, 1)=myrank-nproc_d_g_mul_dm+nproc_d_g_dm(1)-1
    else
      neig_ng(2, 1)=myrank-1
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(1)==0.and.imrs(1)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 1)=comm_proc_null
      case(3)
        neig_ng(2, 1)=myrank-nproc_d_o_mul-1+nproc_d_o_mul*nproc_d_g_dm(1)+nproc_d_o(1)
      end select
    else if(imrs(1)==0) then
      neig_ng(2, 1)=myrank-nproc_d_o_mul-1+nproc_d_o_mul*nproc_d_g_dm(1)
    else
      neig_ng(2, 1)=myrank-nproc_d_o_mul
    end if
  end if
  
  if(process_allocation=='orbital_sequential')then
    if(imr(2)==nproc_d_o(2)-1.and.imrs(2)==nproc_d_g_dm(2)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 2)=comm_proc_null
      case(3)
        neig_ng(1, 2)=myrank+nproc_d_g_mul_dm*nproc_d_o(1)    &
                                      -(nproc_d_g_dm(2)-1)*nproc_d_g_dm(1)-nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2)
      end select 
    else if(imrs(2)==nproc_d_g_dm(2)-1) then
      neig_ng(1, 2)=myrank+nproc_d_g_mul_dm*nproc_d_o(1)    &
                                      -(nproc_d_g_dm(2)-1)*nproc_d_g_dm(1)
    else
      neig_ng(1, 2)=myrank+nproc_d_g_dm(1)
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(2)==nproc_d_o(2)-1.and.imrs(2)==nproc_d_g_dm(2)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 2)=comm_proc_null
      case(3)
        neig_ng(1, 2)=myrank+nproc_d_o_mul*nproc_d_g_dm(1)  &
              +nproc_d_o(1)-nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)-nproc_d_o(1)*nproc_d_o(2)
      end select
    else if(imrs(2)==nproc_d_g_dm(2)-1) then
      neig_ng(1, 2)=myrank+nproc_d_o_mul*nproc_d_g_dm(1)  &
              +nproc_d_o(1)-nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    else
      neig_ng(1, 2)=myrank+nproc_d_o_mul*nproc_d_g_dm(1)
    end if
  end if
  
  if(process_allocation=='orbital_sequential')then
    if(imr(2)==0.and.imrs(2)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 2)=comm_proc_null
      case(3)
        neig_ng(2, 2)=myrank-nproc_d_g_mul_dm*nproc_d_o(1)    &
                                      +(nproc_d_g_dm(2)-1)*nproc_d_g_dm(1)+nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2)
      end select
    else if(imrs(2)==0) then
      neig_ng(2, 2)=myrank-nproc_d_g_mul_dm*nproc_d_o(1)    &
                                      +(nproc_d_g_dm(2)-1)*nproc_d_g_dm(1)
    else
      neig_ng(2, 2)=myrank-nproc_d_g_dm(1)
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(2)==0.and.imrs(2)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 2)=comm_proc_null
      case(3)
        neig_ng(2, 2)=myrank-nproc_d_o_mul*nproc_d_g_dm(1)  &
              -nproc_d_o(1)+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)+nproc_d_o(1)*nproc_d_o(2)
      end select
    else if(imrs(2)==0) then
      neig_ng(2, 2)=myrank-nproc_d_o_mul*nproc_d_g_dm(1)  &
                -nproc_d_o(1)+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    else
      neig_ng(2, 2)=myrank-nproc_d_o_mul*nproc_d_g_dm(1)
    end if
  end if

  if(process_allocation=='orbital_sequential')then
    if(imr(3)==nproc_d_o(3)-1.and.imrs(3)==nproc_d_g_dm(3)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 3)=comm_proc_null
      case(3)
        neig_ng(1, 3)=myrank+nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2) &
                                      -(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2) &
                                      -nproc_d_g_mul_dm*nproc_d_o_mul
      end select
    else if(imrs(3)==nproc_d_g_dm(3)-1) then
      neig_ng(1, 3)=myrank+nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2)   &
                                      -(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    else
      neig_ng(1, 3)=myrank+nproc_d_g_dm(1)*nproc_d_g_dm(2)
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(3)==nproc_d_o(3)-1.and.imrs(3)==nproc_d_g_dm(3)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 3)=comm_proc_null
      case(3)
        neig_ng(1, 3)=myrank+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2) &
                +nproc_d_o(1)*nproc_d_o(2)   &
                -nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)-nproc_d_o_mul
      end select
    else if(imrs(3)==nproc_d_g_dm(3)-1) then
      neig_ng(1, 3)=myrank+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
                +nproc_d_o(1)*nproc_d_o(2)   &
                -nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)
    else
      neig_ng(1, 3)=myrank+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    end if
  end if
  
  if(process_allocation=='orbital_sequential')then
    if(imr(3)==0.and.imrs(3)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 3)=comm_proc_null
      case(3)
        neig_ng(2, 3)=myrank-nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2) &
                                      +(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2) &
                                      +nproc_d_g_mul_dm*nproc_d_o_mul
      end select
    else if(imrs(3)==0) then
      neig_ng(2, 3)=myrank-nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2)   &
                                      +(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    else
      neig_ng(2, 3)=myrank-nproc_d_g_dm(1)*nproc_d_g_dm(2)
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(3)==0.and.imrs(3)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 3)=comm_proc_null
      case(3)
        neig_ng(2, 3)=myrank-nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2) &
                                      +(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2) &
                                      +nproc_d_g_mul_dm*nproc_d_o_mul
      end select
    else if(imrs(3)==0) then
      neig_ng(2, 3)=myrank-nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
                -nproc_d_o(1)*nproc_d_o(2)   &
                +nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)
    else
      neig_ng(2, 3)=myrank-nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    end if
  end if
  
  return
end subroutine create_sendrecv_neig_ng

END MODULE init_sendrecv_sub

