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

  use scf_data, only: imr, nproc_d_o
  ! Note: the above will be moved to s_parallel in near future

  implicit none
  integer, intent(out) :: neig_mg(1:2, 1:3)
  type(s_orbital_parallel), intent(in) :: ob_para_info
  integer, intent(in) :: iperiodic

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


subroutine create_sendrecv_neig_ng(neig_ng, ob_para_info, iperiodic)
  use structures, only: s_orbital_parallel
  use sendrecv_grid, only: init_sendrecv_grid
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_proc_null
  use inputoutput, only: process_allocation

  use scf_data, only: imr, imrs, nproc_d_o, nproc_d_g_dm, nproc_d_g_mul_dm, nproc_d_o_mul

  ! Note: the above will be moved to s_parallel in near future

  implicit none
  integer, intent(out) :: neig_ng(1:2, 1:3)
  type(s_orbital_parallel), intent(in) :: ob_para_info
  integer, intent(in) :: iperiodic
  

  if(process_allocation=='orbital_sequential')then
    if(imr(1)==nproc_d_o(1)-1.and.imrs(1)==nproc_d_g_dm(1)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 1)=comm_proc_null
      case(3)
        neig_ng(1, 1)=nproc_id_global+nproc_d_g_mul_dm-nproc_d_g_dm(1)+1-nproc_d_g_mul_dm*nproc_d_o(1)
      end select
    else if(imrs(1)==nproc_d_g_dm(1)-1) then
      neig_ng(1, 1)=nproc_id_global+nproc_d_g_mul_dm-nproc_d_g_dm(1)+1
    else
      neig_ng(1, 1)=nproc_id_global+1
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(1)==nproc_d_o(1)-1.and.imrs(1)==nproc_d_g_dm(1)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 1)=comm_proc_null
      case(3)
        neig_ng(1, 1)=nproc_id_global+nproc_d_o_mul+1-nproc_d_o_mul*nproc_d_g_dm(1)-nproc_d_o(1)
      end select
    else if(imrs(1)==nproc_d_g_dm(1)-1) then
      neig_ng(1, 1)=nproc_id_global+nproc_d_o_mul+1-nproc_d_o_mul*nproc_d_g_dm(1)
    else
      neig_ng(1, 1)=nproc_id_global+nproc_d_o_mul
    end if
  end if
  
  if(process_allocation=='orbital_sequential')then
    if(imr(1)==0.and.imrs(1)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 1)=comm_proc_null
      case(3)
        neig_ng(2, 1)=nproc_id_global-nproc_d_g_mul_dm+nproc_d_g_dm(1)-1+nproc_d_g_mul_dm*nproc_d_o(1)
      end select
    else if(imrs(1)==0) then
      neig_ng(2, 1)=nproc_id_global-nproc_d_g_mul_dm+nproc_d_g_dm(1)-1
    else
      neig_ng(2, 1)=nproc_id_global-1
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(1)==0.and.imrs(1)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 1)=comm_proc_null
      case(3)
        neig_ng(2, 1)=nproc_id_global-nproc_d_o_mul-1+nproc_d_o_mul*nproc_d_g_dm(1)+nproc_d_o(1)
      end select
    else if(imrs(1)==0) then
      neig_ng(2, 1)=nproc_id_global-nproc_d_o_mul-1+nproc_d_o_mul*nproc_d_g_dm(1)
    else
      neig_ng(2, 1)=nproc_id_global-nproc_d_o_mul
    end if
  end if
  
  if(process_allocation=='orbital_sequential')then
    if(imr(2)==nproc_d_o(2)-1.and.imrs(2)==nproc_d_g_dm(2)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 2)=comm_proc_null
      case(3)
        neig_ng(1, 2)=nproc_id_global+nproc_d_g_mul_dm*nproc_d_o(1)    &
                                      -(nproc_d_g_dm(2)-1)*nproc_d_g_dm(1)-nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2)
      end select 
    else if(imrs(2)==nproc_d_g_dm(2)-1) then
      neig_ng(1, 2)=nproc_id_global+nproc_d_g_mul_dm*nproc_d_o(1)    &
                                      -(nproc_d_g_dm(2)-1)*nproc_d_g_dm(1)
    else
      neig_ng(1, 2)=nproc_id_global+nproc_d_g_dm(1)
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(2)==nproc_d_o(2)-1.and.imrs(2)==nproc_d_g_dm(2)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 2)=comm_proc_null
      case(3)
        neig_ng(1, 2)=nproc_id_global+nproc_d_o_mul*nproc_d_g_dm(1)  &
              +nproc_d_o(1)-nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)-nproc_d_o(1)*nproc_d_o(2)
      end select
    else if(imrs(2)==nproc_d_g_dm(2)-1) then
      neig_ng(1, 2)=nproc_id_global+nproc_d_o_mul*nproc_d_g_dm(1)  &
              +nproc_d_o(1)-nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    else
      neig_ng(1, 2)=nproc_id_global+nproc_d_o_mul*nproc_d_g_dm(1)
    end if
  end if
  
  if(process_allocation=='orbital_sequential')then
    if(imr(2)==0.and.imrs(2)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 2)=comm_proc_null
      case(3)
        neig_ng(2, 2)=nproc_id_global-nproc_d_g_mul_dm*nproc_d_o(1)    &
                                      +(nproc_d_g_dm(2)-1)*nproc_d_g_dm(1)+nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2)
      end select
    else if(imrs(2)==0) then
      neig_ng(2, 2)=nproc_id_global-nproc_d_g_mul_dm*nproc_d_o(1)    &
                                      +(nproc_d_g_dm(2)-1)*nproc_d_g_dm(1)
    else
      neig_ng(2, 2)=nproc_id_global-nproc_d_g_dm(1)
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(2)==0.and.imrs(2)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 2)=comm_proc_null
      case(3)
        neig_ng(2, 2)=nproc_id_global-nproc_d_o_mul*nproc_d_g_dm(1)  &
              -nproc_d_o(1)+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)+nproc_d_o(1)*nproc_d_o(2)
      end select
    else if(imrs(2)==0) then
      neig_ng(2, 2)=nproc_id_global-nproc_d_o_mul*nproc_d_g_dm(1)  &
                -nproc_d_o(1)+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    else
      neig_ng(2, 2)=nproc_id_global-nproc_d_o_mul*nproc_d_g_dm(1)
    end if
  end if

  if(process_allocation=='orbital_sequential')then
    if(imr(3)==nproc_d_o(3)-1.and.imrs(3)==nproc_d_g_dm(3)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 3)=comm_proc_null
      case(3)
        neig_ng(1, 3)=nproc_id_global+nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2) &
                                      -(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2) &
                                      -nproc_d_g_mul_dm*nproc_d_o_mul
      end select
    else if(imrs(3)==nproc_d_g_dm(3)-1) then
      neig_ng(1, 3)=nproc_id_global+nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2)   &
                                      -(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    else
      neig_ng(1, 3)=nproc_id_global+nproc_d_g_dm(1)*nproc_d_g_dm(2)
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(3)==nproc_d_o(3)-1.and.imrs(3)==nproc_d_g_dm(3)-1) then
      select case(iperiodic)
      case(0)
        neig_ng(1, 3)=comm_proc_null
      case(3)
        neig_ng(1, 3)=nproc_id_global+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2) &
                +nproc_d_o(1)*nproc_d_o(2)   &
                -nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)-nproc_d_o_mul
      end select
    else if(imrs(3)==nproc_d_g_dm(3)-1) then
      neig_ng(1, 3)=nproc_id_global+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
                +nproc_d_o(1)*nproc_d_o(2)   &
                -nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)
    else
      neig_ng(1, 3)=nproc_id_global+nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    end if
  end if
  
  if(process_allocation=='orbital_sequential')then
    if(imr(3)==0.and.imrs(3)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 3)=comm_proc_null
      case(3)
        neig_ng(2, 3)=nproc_id_global-nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2) &
                                      +(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2) &
                                      +nproc_d_g_mul_dm*nproc_d_o_mul
      end select
    else if(imrs(3)==0) then
      neig_ng(2, 3)=nproc_id_global-nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2)   &
                                      +(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    else
      neig_ng(2, 3)=nproc_id_global-nproc_d_g_dm(1)*nproc_d_g_dm(2)
    end if
  else if(process_allocation=='grid_sequential')then
    if(imr(3)==0.and.imrs(3)==0) then
      select case(iperiodic)
      case(0)
        neig_ng(2, 3)=comm_proc_null
      case(3)
        neig_ng(2, 3)=nproc_id_global-nproc_d_g_mul_dm*nproc_d_o(1)*nproc_d_o(2) &
                                      +(nproc_d_g_dm(3)-1)*nproc_d_g_dm(1)*nproc_d_g_dm(2) &
                                      +nproc_d_g_mul_dm*nproc_d_o_mul
      end select
    else if(imrs(3)==0) then
      neig_ng(2, 3)=nproc_id_global-nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
                -nproc_d_o(1)*nproc_d_o(2)   &
                +nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)
    else
      neig_ng(2, 3)=nproc_id_global-nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
    end if
  end if
  
  return
end subroutine create_sendrecv_neig_ng


!=====================================================================
subroutine init_itype
use salmon_parallel, only: nproc_size_global, nproc_id_global
use scf_data, only: inum_Mxin, inum_Mxin_s, nd

implicit none

integer :: ilap_s
integer :: ibox
integer :: isize1(3),isize2(3)
integer :: isubsize1(3),isubsize2(3)
integer :: istart1(3),istart2(3)
integer :: inum_Mxin2(3,0:nproc_size_global-1)

do ilap_s=1,2

  if(ilap_s==1)then
    inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin(1:3,0:nproc_size_global-1)
  else if(ilap_s==2)then
    inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin_s(1:3,0:nproc_size_global-1)
  end if

!send from idw to iup

  isize1(1:3)=inum_Mxin2(1:3,nproc_id_global)+2*Nd
  isize2(1:3)=inum_Mxin2(1:3,nproc_id_global)+2*Nd

  isubsize1(1)=Nd
  isubsize1(2)=inum_Mxin2(2,nproc_id_global)
  isubsize1(3)=inum_Mxin2(3,nproc_id_global)
  istart1(1)=inum_Mxin2(1,nproc_id_global)
  istart1(2)=Nd
  istart1(3)=Nd
  isubsize2(1)=Nd
  isubsize2(2)=inum_Mxin2(2,nproc_id_global)
  isubsize2(3)=inum_Mxin2(3,nproc_id_global)
  istart2(1)=0
  istart2(2)=Nd
  istart2(3)=Nd

  if(ilap_s==1)then
    ibox=1
  else if(ilap_s==2)then
    ibox=13
  end if

!send from iup to idw

  istart1(1)=Nd
  istart2(1)=Nd+inum_Mxin2(1,nproc_id_global)

  if(ilap_s==1)then
    ibox=3
  else if(ilap_s==2)then
    ibox=15
  end if

!send from jdw to jup

  isize1(1:3)=inum_Mxin2(1:3,nproc_id_global)+2*Nd
  isubsize1(1)=inum_Mxin2(1,nproc_id_global)
  isubsize1(2)=Nd
  isubsize1(3)=inum_Mxin2(3,nproc_id_global)
  istart1(1)=Nd
  istart1(2)=inum_Mxin2(2,nproc_id_global)
  istart1(3)=Nd
  isize2(1:3)=inum_Mxin2(1:3,nproc_id_global)+2*Nd
  isubsize2(1)=inum_Mxin2(1,nproc_id_global)
  isubsize2(2)=Nd
  isubsize2(3)=inum_Mxin2(3,nproc_id_global)
  istart2(1)=Nd
  istart2(2)=0
  istart2(3)=Nd

  if(ilap_s==1)then
    ibox=5
  else if(ilap_s==2)then
    ibox=17
  end if

!send from jup to jdw
  
  istart1(2)=Nd
  istart2(2)=Nd+inum_Mxin2(2,nproc_id_global)

  if(ilap_s==1)then
    ibox=7
  else if(ilap_s==2)then
    ibox=19
  end if

!send from kup to kdw

  isubsize1(1)=inum_Mxin2(1,nproc_id_global)
  isubsize1(2)=inum_Mxin2(2,nproc_id_global)
  isubsize1(3)=Nd
  istart1(1)=Nd
  istart1(2)=Nd
  istart1(3)=inum_Mxin2(3,nproc_id_global)
  isubsize2(1)=inum_Mxin2(1,nproc_id_global)
  isubsize2(2)=inum_Mxin2(2,nproc_id_global)
  isubsize2(3)=Nd
  istart2(1)=Nd
  istart2(2)=Nd
  istart2(3)=0

  if(ilap_s==1)then
    ibox=9
  else if(ilap_s==2)then
    ibox=21
  end if

!send from kup to kdw

  istart1(3)=Nd
  istart2(3)=Nd+inum_Mxin2(3,nproc_id_global)

  if(ilap_s==1)then
    ibox=11
  else if(ilap_s==2)then
    ibox=23
  end if

end do

end subroutine init_itype

!======================================================================
subroutine init_sendrecv_matrix
use salmon_parallel, only: nproc_size_global
use scf_data, only: inum_Mxin_s, inum_Mxin
integer :: inum_Mxin2(3,0:nproc_size_global-1)

inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin_s(1:3,0:nproc_size_global-1)

inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin(1:3,0:nproc_size_global-1)

inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin_s(1:3,0:nproc_size_global-1)

end subroutine init_sendrecv_matrix

END MODULE init_sendrecv_sub

