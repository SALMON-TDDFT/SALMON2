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
MODULE new_world_sub

use scf_data

implicit none

integer,allocatable :: icorr_polenum(:)

! FFTE routine
integer :: iquot
integer :: i11,i12,i13,i14,i15

CONTAINS

!=======================================================================
subroutine make_new_world(info,info_field)
use inputoutput, only: process_allocation
use structures, only: s_orbital_parallel, s_field_parallel
use salmon_parallel
use salmon_communication, only: comm_create_group, comm_get_groupinfo, &
                                comm_summation
use misc_routines, only: get_wtime
implicit none
type(s_orbital_parallel),intent(inout) :: info
type(s_field_parallel),intent(inout) :: info_field
integer :: ii
integer :: i1,i2,i3,i4,i5
integer :: ix,iy,iz
integer :: ixs,iys,izs
integer :: ibox
integer :: icolor,ikey

!new_world for comm_kgrid
if(process_allocation=='orbital_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
      if(nproc_id_global==ibox)then
        icolor=i5*nproc_ob+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
        ikey=i4
      end if
    end do
    end do
  end do
  end do
  end do
else if(process_allocation=='grid_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
      if(nproc_id_global==ibox)then
        icolor=i5+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_k
        ikey=i4
      end if
    end do
    end do
  end do
  end do
  end do
end if

info%icomm_o = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(info%icomm_o, nproc_id_kgrid, nproc_size_kgrid)

!new_world for comm_korbital
if(process_allocation=='orbital_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
      if(nproc_id_global==ibox)then
        icolor=i5*nproc_ob+i4
        ikey=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
      end if
    end do
    end do
  end do
  end do
  end do
else if(process_allocation=='grid_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
      if(nproc_id_global==ibox)then
        icolor=i5*nproc_ob+i4
        ikey=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
      end if
    end do
    end do
  end do
  end do
  end do
end if

info%icomm_r = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(info%icomm_r, info%id_r, info%isize_r)

!new_world for comm_k
if(process_allocation=='orbital_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i5
          ikey=i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=2*i5+0
            ikey=i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob_spin(1)
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=2*i5+1
            ikey=i4-nproc_ob_spin(1)+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob_spin(2)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
else if(process_allocation=='grid_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i5
          ikey=i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=2*i5+0
            ikey=i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob_spin(1)
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=2*i5+1
            ikey=i4-nproc_ob_spin(1)+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob_spin(2)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
end if
 
info%icomm_ro = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(info%icomm_ro, info%id_ro, info%isize_ro)

!new_world for comm_grid
if(process_allocation=='orbital_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
          ikey=i5*nproc_ob+i4
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
            ikey=i5*nproc_ob_spin(1)+i4
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
            ikey=i5*nproc_ob_spin(2)+i4-nproc_ob_spin(1)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
else if(process_allocation=='grid_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
          ikey=i5*nproc_ob+i4
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
            ikey=i5*nproc_ob_spin(1)+i4
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)
            ikey=i5*nproc_ob_spin(2)+i4-nproc_ob_spin(1)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
end if
 
info%icomm_ko = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(info%icomm_ko, info%id_ko, info%isize_ko)

!new_world for comm_orbitalgrid
if(process_allocation=='orbital_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))*nproc_ob*nproc_k
      if(nproc_id_global==ibox)then
        icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)+i4*nproc_d_o_mul
        ikey=i5
      end if
    end do
    end do
  end do
  end do
  end do
else if(process_allocation=='grid_sequential')then
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2))+(i5*nproc_ob+i4)*nproc_d_o_mul
      if(nproc_id_global==ibox)then
        icolor=i1+i2*nproc_d_o(1)+i3*nproc_d_o(1)*nproc_d_o(2)+i4*nproc_d_o_mul
        ikey=i5
      end if
    end do
    end do
  end do
  end do
  end do
end if
 
info%icomm_k = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(info%icomm_k, info%id_k, info%isize_k)

!new_world for comm_mesh_s

nproc_d_g_mul_dm=nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)

if(process_allocation=='orbital_sequential')then
  do i4=0,nproc_size_global/nproc_d_o_mul/nproc_d_g_mul_dm-1
  do i3=0,nproc_d_g_dm(3)-1
  do i2=0,nproc_d_g_dm(2)-1
  do i1=0,nproc_d_g_dm(1)-1
    do ii=0,nproc_d_o_mul-1
      ibox=i1+i2*nproc_d_g_dm(1)   &
             +i3*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
             +i4*nproc_d_g_mul_dm  &
             +ii*nproc_size_global/nproc_d_o_mul
      if(nproc_id_global==ibox)then
        icolor=i4
        ikey=i1+i2*nproc_d_g_dm(1)   &
               +i3*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
               +ii*nproc_d_g_mul_dm
      end if
    end do
  end do
  end do
  end do
  end do
else if(process_allocation=='grid_sequential')then
  do i4=0,nproc_size_global/nproc_d_o_mul/nproc_d_g_mul_dm-1
  do i3=0,nproc_d_g_dm(3)-1
  do i2=0,nproc_d_g_dm(2)-1
  do i1=0,nproc_d_g_dm(1)-1
    do ii=0,nproc_d_o_mul-1
      ibox=ii+i1*nproc_d_o_mul+i2*nproc_d_o_mul*nproc_d_g_dm(1)   &
            +i3*nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
            +i4*nproc_d_o_mul*nproc_d_g_mul_dm
      if(nproc_id_global==ibox)then
        icolor=i4
        ikey=ii+i1*nproc_d_o_mul+i2*nproc_d_o_mul*nproc_d_g_dm(1)   &
              +i3*nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)
      end if
    end do
  end do
  end do
  end do
  end do
end if

nproc_group_h = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_h, nproc_id_h, nproc_size_h)

if(process_allocation=='orbital_sequential')then
  do iz=0,nproc_d_o(3)-1
  do iy=0,nproc_d_o(2)-1
  do ix=0,nproc_d_o(1)-1
    do i4=0,nproc_size_global/nproc_d_o_mul/nproc_d_g_mul_dm-1
    do izs=0,nproc_d_g_dm(3)-1
    do iys=0,nproc_d_g_dm(2)-1
    do ixs=0,nproc_d_g_dm(1)-1
      ibox=ixs+iys*nproc_d_g_dm(1)   &
              +izs*nproc_d_g_dm(1)*nproc_d_g_dm(2)   &
              +i4*nproc_d_g_mul_dm    &
              +ix*nproc_size_global/nproc_d_o_mul    &
              +iy*nproc_size_global/nproc_d_o_mul*nproc_d_o(1)   &
              +iz*nproc_size_global/nproc_d_o_mul*nproc_d_o(1)*nproc_d_o(2)
      if(nproc_id_global==ibox)then
        imr(1)=ix
        imr(2)=iy
        imr(3)=iz
        imrs(1)=ixs
        imrs(2)=iys
        imrs(3)=izs
        igroup=i4
      end if
    end do 
    end do 
    end do 
  end do 
  end do 
  end do
  end do
else if(process_allocation=='grid_sequential')then
  do i4=0,nproc_size_global/nproc_d_o_mul/nproc_d_g_mul_dm-1
  do izs=0,nproc_d_g_dm(3)-1
  do iys=0,nproc_d_g_dm(2)-1
  do ixs=0,nproc_d_g_dm(1)-1
    do iz=0,nproc_d_o(3)-1
    do iy=0,nproc_d_o(2)-1
    do ix=0,nproc_d_o(1)-1
      ibox=ix+iy*nproc_d_o(1)+iz*nproc_d_o(1)*nproc_d_o(2)  &
             +ixs*nproc_d_o_mul  &
             +iys*nproc_d_o_mul*nproc_d_g_dm(1)  &
             +izs*nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
             +i4*nproc_d_o_mul*nproc_d_g_mul_dm
      if(nproc_id_global==ibox)then
        imr(1)=ix
        imr(2)=iy
        imr(3)=iz
        imrs(1)=ixs
        imrs(2)=iys
        imrs(3)=izs
        igroup=i4
      end if
    end do 
    end do 
    end do 
  end do 
  end do 
  end do
  end do
end if

if(process_allocation=='orbital_sequential')then
  icolor=imrs(2)+imrs(3)*nproc_d_g_dm(2)   &
                +igroup*nproc_d_g_dm(2)*nproc_d_g_dm(3)   &
                +imr(2)*nproc_size_global/nproc_d_o_mul/nproc_d_g_dm(1)   &
                +imr(3)*nproc_size_global/nproc_d_o_mul/nproc_d_g_dm(1)*nproc_d_o(2)
  ikey=imrs(1)+imr(1)*nproc_d_g_dm(1)
else if(process_allocation=='grid_sequential')then
  icolor=imr(2)+imr(3)*nproc_d_o(2)   &
               +imrs(2)*nproc_d_o(2)*nproc_d_o(3)   &
               +imrs(3)*nproc_d_o(2)*nproc_d_o(3)*nproc_d_g_dm(2)  &
               +nproc_id_global/(nproc_d_o_mul*nproc_d_g_mul_dm)*nproc_d_o(2)  &
               *nproc_d_o(3)*nproc_d_g_dm(2)*nproc_d_g_dm(3)
!  ikey=imr(1)+imrs(1)*nproc_d_o(1)
  ikey=imrs(1)+imr(1)*nproc_d_g_dm(1)
end if

info_field%icomm(1) = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(info_field%icomm(1), info_field%id(1), info_field%isize(1))

if(process_allocation=='orbital_sequential')then
  icolor=imrs(1)+imrs(3)*nproc_d_g_dm(1)   &
                +igroup*nproc_d_g_dm(1)*nproc_d_g_dm(3)   &
                +imr(1)*nproc_size_global/nproc_d_o_mul/nproc_d_g_dm(2)   &
                +imr(3)*nproc_size_global/nproc_d_o_mul/nproc_d_g_dm(2)*nproc_d_o(1)
  ikey=imrs(2)+imr(2)*nproc_d_g_dm(2)
else if(process_allocation=='grid_sequential')then
  icolor=imr(1)+imr(3)*nproc_d_o(1)   &
               +imrs(1)*nproc_d_o(1)*nproc_d_o(3)   &
               +imrs(3)*nproc_d_o(1)*nproc_d_o(3)*nproc_d_g_dm(1)  &
               +nproc_id_global/(nproc_d_o_mul*nproc_d_g_mul_dm)*nproc_d_o(1)*nproc_d_o(3)  &
               *nproc_d_g_dm(1)*nproc_d_g_dm(3)
!  ikey=imr(2)+imrs(2)*nproc_d_o(2)
  ikey=imrs(2)+imr(2)*nproc_d_g_dm(2)
end if

info_field%icomm(2) = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(info_field%icomm(2), info_field%id(2), info_field%isize(2))

if(process_allocation=='orbital_sequential')then
  icolor=imrs(1)+imrs(2)*nproc_d_g_dm(1)   &
                +igroup*nproc_d_g_dm(1)*nproc_d_g_dm(2)   &
                +imr(1)*nproc_size_global/nproc_d_o_mul/nproc_d_g_dm(3)   &
                +imr(2)*nproc_size_global/nproc_d_o_mul/nproc_d_g_dm(3)*nproc_d_o(1)
  ikey=imrs(3)+imr(3)*nproc_d_g_dm(3)
else if(process_allocation=='grid_sequential')then
  icolor=imr(1)+imr(2)*nproc_d_o(1)   &
               +imrs(1)*nproc_d_o(1)*nproc_d_o(2)   &
               +imrs(2)*nproc_d_o(1)*nproc_d_o(2)*nproc_d_g_dm(1)  &
               +nproc_id_global/(nproc_d_o_mul*nproc_d_g_mul_dm)*nproc_d_o(1)*nproc_d_o(2)  &
               *nproc_d_g_dm(1)*nproc_d_g_dm(2)
!  ikey=imr(3)+imrs(3)*nproc_d_o(3)
  ikey=imrs(3)+imr(3)*nproc_d_g_dm(3)
end if

info_field%icomm(3) = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(info_field%icomm(3), info_field%id(3), info_field%isize(3))

if(process_allocation=='orbital_sequential')then
  do ii=0,nproc_d_o_mul-1
    do i4=0,nproc_size_global/nproc_d_o_mul/nproc_d_g_mul_dm-1
    do i3=0,nproc_d_g_dm(3)-1
    do i2=0,nproc_d_g_dm(2)-1
    do i1=0,nproc_d_g_dm(1)-1
      ibox=i1+i2*nproc_d_g_dm(1)   &
             +i3*nproc_d_g_dm(1)*nproc_d_g_dm(2)   &
             +i4*nproc_d_g_mul_dm   &
             +ii*nproc_size_global/nproc_d_o_mul
      if(nproc_id_global==ibox)then
        icolor=i4+ii*nproc_size_global/nproc_d_o_mul/nproc_d_g_mul_dm
        ikey=i1+i2*nproc_d_g_dm(1)+i3*nproc_d_g_dm(1)*nproc_d_g_dm(2)
      end if
    end do
    end do
    end do
    end do
  end do
else if(process_allocation=='grid_sequential')then
  do i4=0,nproc_size_global/nproc_d_o_mul/nproc_d_g_mul_dm-1
  do i3=0,nproc_d_g_dm(3)-1
  do i2=0,nproc_d_g_dm(2)-1
  do i1=0,nproc_d_g_dm(1)-1
    do ii=0,nproc_d_o_mul-1
      ibox=ii+i1*nproc_d_o_mul+i2*nproc_d_o_mul*nproc_d_g_dm(1)   &
            +i3*nproc_d_o_mul*nproc_d_g_dm(1)*nproc_d_g_dm(2)  &
            +i4*nproc_d_o_mul*nproc_d_g_mul_dm
      if(nproc_id_global==ibox)then
        icolor=ii+i4*nproc_d_o_mul
        ikey=i1+i2*nproc_d_g_dm(1)+i3*nproc_d_g_dm(1)*nproc_d_g_dm(2)
      end if
    end do
  end do
  end do
  end do
  end do
end if

nproc_group_korbital_vhxc = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_korbital_vhxc, nproc_id_korbital_vhxc, nproc_size_korbital_vhxc)

!call FACTOR(nproc,LNPU)
!NPUZ=(2**(LNPU(1)/2))*(3**(LNPU(2)/2))*(5**(LNPU(3)/2))
!NPUY=nproc/NPUZ

!if(iflag_hartree==4)then

! communicators for FFTE routine
  NPUW=nproc_d_g_dm(1)*nproc_d_o(1)
  NPUY=nproc_d_g_dm(2)*nproc_d_o(2)
  NPUZ=nproc_d_g_dm(3)*nproc_d_o(3)

  icolor=info_field%id(3)+info_field%id(1)*NPUZ
  ikey=info_field%id(2)
  nproc_group_icommy = comm_create_group(nproc_group_global, icolor, ikey)
  call comm_get_groupinfo(nproc_group_icommy, nproc_id_icommy, nproc_size_icommy)

  icolor=info_field%id(2)+info_field%id(1)*NPUY
  ikey=info_field%id(3)
  nproc_group_icommz = comm_create_group(nproc_group_global, icolor, ikey)
  call comm_get_groupinfo(nproc_group_icommz, nproc_id_icommz, nproc_size_icommz)

  icolor=info_field%id(2)+info_field%id(3)*NPUY
  ikey=info_field%id(1)
  nproc_group_icommw = comm_create_group(nproc_group_global, icolor, ikey)
  call comm_get_groupinfo(nproc_group_icommw, nproc_id_icommw, nproc_size_icommw)

  iquot=nproc_id_global/(NPUY*NPUZ)
  
  i11=mod(nproc_id_global,nproc_d_o(2)*nproc_d_o(3))
  i12=i11/nproc_d_o(2)
  i13=i12*nproc_d_o(3)
  i14=nproc_id_global/(NPUY*nproc_d_o(3))
  icolor=i13+i14+iquot*NPUZ
  
  i11=mod(nproc_id_global,nproc_d_o(2))
  i12=nproc_id_global/(nproc_d_o(2)*nproc_d_o(3))
  ikey=i11*NPUY/nproc_d_o(2)+mod(i12,NPUY/nproc_d_o(2))
  
  nproc_group_icommy_copy = comm_create_group(nproc_group_global, icolor, ikey)
  call comm_get_groupinfo(nproc_group_icommy_copy, nproc_id_icommy_copy, nproc_size_icommy_copy)

!end if

  nproc_group_tdks = nproc_group_global
  nproc_id_tdks    = nproc_id_global
  nproc_size_tdks  = nproc_size_global

end subroutine make_new_world

!=====================================================================
subroutine make_corr_pole(ng,poisson)
use structures, only: s_rgrid,s_poisson
implicit none
type(s_rgrid), intent(in) :: ng
type(s_poisson),intent(inout) :: poisson
integer :: a,i
integer :: ix,iy,iz
integer :: ibox
integer :: j1,j2,j3
integer,allocatable :: ista_Mxin_pole(:,:)
integer,allocatable :: iend_Mxin_pole(:,:)
integer,allocatable :: inum_Mxin_pole(:,:)
integer,allocatable :: iflag_pole(:)
integer :: amin,amax
real(8) :: rmin,r
real(8),allocatable :: Rion2(:,:)
integer,allocatable :: nearatomnum(:,:,:)
integer,allocatable :: inv_icorr_polenum(:)
integer :: maxval_ig_num

if(layout_multipole==2)then

  if(iflag_ps==1)then
    amax=MI
    allocate(Rion2(3,MI))
    Rion2(:,:)=Rion(:,:)
  end if

  allocate(poisson%ig_num(1:amax))
  allocate(nearatomnum(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
  poisson%ig_num=0
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rmin=1.d6
    do a=1,amax
      r=sqrt( (gridcoo(ix,1)-Rion2(1,a))**2      &
            + (gridcoo(iy,2)-Rion2(2,a))**2      &
            + (gridcoo(iz,3)-Rion2(3,a))**2 )
      if ( r < rmin ) then
        rmin=r ; amin=a
      end if
    end do
    poisson%ig_num(amin)=poisson%ig_num(amin)+1
    nearatomnum(ix,iy,iz)=amin
  end do
  end do
  end do

  allocate(poisson%ipole_tbl(1:amax))
  allocate(inv_icorr_polenum(1:amax))
  poisson%ipole_tbl=0
  inv_icorr_polenum=0
  ibox=0
  do a=1,amax
    if(poisson%ig_num(a)>=1)then
      ibox=ibox+1
      poisson%ipole_tbl(ibox)=a
      inv_icorr_polenum(a)=ibox
    end if
  end do
  poisson%npole_partial=ibox

  maxval_ig_num=maxval(poisson%ig_num(:)) 
  allocate(poisson%ig(3,maxval(poisson%ig_num(:)),poisson%npole_partial))

  poisson%ig_num(:)=0

  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    ibox=inv_icorr_polenum(nearatomnum(ix,iy,iz))
    poisson%ig_num(ibox)=poisson%ig_num(ibox)+1
    poisson%ig(1,poisson%ig_num(ibox),ibox)=ix
    poisson%ig(2,poisson%ig_num(ibox),ibox)=iy
    poisson%ig(3,poisson%ig_num(ibox),ibox)=iz
  end do
  end do
  end do

  if(iflag_ps==1)then
    deallocate(Rion2)
  end if
  deallocate(nearatomnum)
  deallocate(inv_icorr_polenum)

else if(layout_multipole==3)then

  allocate(ista_Mxin_pole(3,0:poisson%npole_total-1))
  allocate(iend_Mxin_pole(3,0:poisson%npole_total-1))
  allocate(inum_Mxin_pole(3,0:poisson%npole_total-1))
  allocate(iflag_pole(1:poisson%npole_total))

  do j3=0,num_multipole_xyz(3)-1
  do j2=0,num_multipole_xyz(2)-1
  do j1=0,num_multipole_xyz(1)-1
    ibox = j1 + num_multipole_xyz(1)*j2 + num_multipole_xyz(1)*num_multipole_xyz(2)*j3 
    ista_Mxin_pole(1,ibox)=j1*lg_num(1)/num_multipole_xyz(1)+lg_sta(1)
    iend_Mxin_pole(1,ibox)=(j1+1)*lg_num(1)/num_multipole_xyz(1)+lg_sta(1)-1
    ista_Mxin_pole(2,ibox)=j2*lg_num(2)/num_multipole_xyz(2)+lg_sta(2)
    iend_Mxin_pole(2,ibox)=(j2+1)*lg_num(2)/num_multipole_xyz(2)+lg_sta(2)-1
    ista_Mxin_pole(3,ibox)=j3*lg_num(3)/num_multipole_xyz(3)+lg_sta(3)
    iend_Mxin_pole(3,ibox)=(j3+1)*lg_num(3)/num_multipole_xyz(3)+lg_sta(3)-1
  end do
  end do
  end do

  iflag_pole=0

  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    do i=1,poisson%npole_total
      if(ista_Mxin_pole(3,i-1)<=iz.and.iend_Mxin_pole(3,i-1)>=iz.and.   &
         ista_Mxin_pole(2,i-1)<=iy.and.iend_Mxin_pole(2,i-1)>=iy.and.   &
         ista_Mxin_pole(1,i-1)<=ix.and.iend_Mxin_pole(1,i-1)>=ix)then
        iflag_pole(i)=1
      end if
    end do
  end do
  end do
  end do

  poisson%npole_partial=0
  do i=1,poisson%npole_total
    if(iflag_pole(i)==1)then
      poisson%npole_partial=poisson%npole_partial+1
    end if
  end do

  allocate(poisson%ipole_tbl(1:poisson%npole_partial))
  allocate(poisson%ig_num(1:poisson%npole_partial))

  ibox=1
  do i=1,poisson%npole_total
    if(iflag_pole(i)==1)then
      poisson%ipole_tbl(ibox)=i
      ibox=ibox+1
    end if
  end do

  poisson%ig_num=0

  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    do i=1,poisson%npole_partial
      if(ista_Mxin_pole(3,poisson%ipole_tbl(i)-1)<=iz.and.iend_Mxin_pole(3,poisson%ipole_tbl(i)-1)>=iz.and.   &
         ista_Mxin_pole(2,poisson%ipole_tbl(i)-1)<=iy.and.iend_Mxin_pole(2,poisson%ipole_tbl(i)-1)>=iy.and.   &
         ista_Mxin_pole(1,poisson%ipole_tbl(i)-1)<=ix.and.iend_Mxin_pole(1,poisson%ipole_tbl(i)-1)>=ix)then
        poisson%ig_num(i)=poisson%ig_num(i)+1
      end if
    end do
  end do
  end do
  end do

  maxval_ig_num=maxval(poisson%ig_num(:)) 
  allocate(poisson%ig(3,maxval(poisson%ig_num(:)),poisson%npole_partial))
 
  poisson%ig_num=0

  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    do i=1,poisson%npole_partial
      if(ista_Mxin_pole(3,poisson%ipole_tbl(i)-1)<=iz.and.iend_Mxin_pole(3,poisson%ipole_tbl(i)-1)>=iz.and.   &
         ista_Mxin_pole(2,poisson%ipole_tbl(i)-1)<=iy.and.iend_Mxin_pole(2,poisson%ipole_tbl(i)-1)>=iy.and.   &
         ista_Mxin_pole(1,poisson%ipole_tbl(i)-1)<=ix.and.iend_Mxin_pole(1,poisson%ipole_tbl(i)-1)>=ix)then
        poisson%ig_num(i)=poisson%ig_num(i)+1
        poisson%ig(1,poisson%ig_num(i),i)=ix
        poisson%ig(2,poisson%ig_num(i),i)=iy
        poisson%ig(3,poisson%ig_num(i),i)=iz
      end if
    end do
  end do
  end do
  end do

  deallocate(ista_Mxin_pole,iend_Mxin_pole,inum_Mxin_pole)
  deallocate(iflag_pole)

end if

end subroutine make_corr_pole

!=====================================================================
subroutine set_ig_bound(ng,poisson)
use structures, only: s_rgrid,s_poisson
use salmon_parallel
implicit none
type(s_rgrid), intent(in) :: ng
type(s_poisson),intent(inout) :: poisson
integer :: ix,iy,iz
integer :: ibox
integer :: icount

ibox=inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)  &
     *inum_Mxin_s(3,nproc_id_global)/minval(inum_Mxin_s(1:3,nproc_id_global))*2*Ndh

allocate( poisson%ig_bound(3,ibox,3) )

icount=0
do iz=ng%is(3),ng%ie(3)
do iy=ng%is(2),ng%ie(2)
do ix=lg_sta(1)-Ndh,lg_sta(1)-1
  icount=icount+1
  poisson%ig_bound(1,icount,1)=ix
  poisson%ig_bound(2,icount,1)=iy
  poisson%ig_bound(3,icount,1)=iz
end do
end do
end do
do iz=ng%is(3),ng%ie(3)
do iy=ng%is(2),ng%ie(2)
do ix=lg_end(1)+1,lg_end(1)+Ndh
  icount=icount+1
  poisson%ig_bound(1,icount,1)=ix
  poisson%ig_bound(2,icount,1)=iy
  poisson%ig_bound(3,icount,1)=iz
end do
end do
end do
icount=0
do iz=ng%is(3),ng%ie(3)
do iy=lg_sta(2)-Ndh,lg_sta(2)-1
do ix=ng%is(1),ng%ie(1)
  icount=icount+1
  poisson%ig_bound(1,icount,2)=ix
  poisson%ig_bound(2,icount,2)=iy
  poisson%ig_bound(3,icount,2)=iz
end do
end do
end do
do iz=ng%is(3),ng%ie(3)
do iy=lg_end(2)+1,lg_end(2)+Ndh
do ix=ng%is(1),ng%ie(1)
  icount=icount+1
  poisson%ig_bound(1,icount,2)=ix
  poisson%ig_bound(2,icount,2)=iy
  poisson%ig_bound(3,icount,2)=iz
end do
end do
end do
icount=0
do iz=lg_sta(3)-Ndh,lg_sta(3)-1
do iy=ng%is(2),ng%ie(2)
do ix=ng%is(1),ng%ie(1)
  icount=icount+1
  poisson%ig_bound(1,icount,3)=ix
  poisson%ig_bound(2,icount,3)=iy
  poisson%ig_bound(3,icount,3)=iz
end do
end do
end do
do iz=lg_end(3)+1,lg_end(3)+Ndh
do iy=ng%is(2),ng%ie(2)
do ix=ng%is(1),ng%ie(1)
  icount=icount+1
  poisson%ig_bound(1,icount,3)=ix
  poisson%ig_bound(2,icount,3)=iy
  poisson%ig_bound(3,icount,3)=iz
end do
end do
end do

return
end subroutine set_ig_bound

!=====================================================================
!======================================================================
subroutine allgatherv_vlocal(ng,info,nspin,Vh,Vpsl,Vxc,Vlocal)
use structures, only: s_rgrid, s_orbital_parallel, s_scalar
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_allgatherv
use timer
implicit none
type(s_rgrid),           intent(in) :: ng
type(s_orbital_parallel),intent(in) :: info
integer       ,intent(in) :: nspin
type(s_scalar),intent(in) :: Vh,Vpsl,Vxc(nspin)
type(s_scalar)            :: Vlocal(nspin)
!
integer :: i
integer :: i1,i2,i3
integer :: ix,iy,iz
integer :: ibox,ibox2,ibox3
real(8),allocatable :: matbox11(:),matbox12(:)
integer :: iscnt
integer,allocatable :: ircnt(:)
integer,allocatable :: idisp(:)
integer :: is

allocate(ircnt(0:nproc_d_g_mul_dm-1))
allocate(idisp(0:nproc_d_g_mul_dm-1))

allocate (matbox11(0:(inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)*inum_Mxin_s(3,nproc_id_global))-1))
allocate (matbox12(0:(mg_num(1)*mg_num(2)*mg_num(3))-1))

iscnt=inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)*inum_Mxin_s(3,nproc_id_global)
if(process_allocation=='orbital_sequential')then
  do i=0,nproc_d_g_mul_dm-1
    ibox=(nproc_id_global/nproc_d_g_mul_dm)*nproc_d_g_mul_dm+i
    ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
  end do
else if(process_allocation=='grid_sequential')then
  do i=0,nproc_d_g_mul_dm-1
    ibox=mod(nproc_id_global,nproc_d_o_mul)+i*nproc_d_o_mul
    ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
  end do
end if

idisp(0)=0
do i=1,nproc_d_g_mul_dm-1
  idisp(i)=idisp(i-1)+ircnt(i-1)
end do

do is=1,nspin
!$OMP parallel do private(ibox3,ix,iy,iz)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    ibox3=ix-ng%is(1)+(iy-ng%is(2))*inum_Mxin_s(1,nproc_id_global)   &
                  +(iz-ng%is(3))*inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)
    matbox11(ibox3) = Vpsl%f(ix,iy,iz) + Vh%f(ix,iy,iz) + Vxc(is)%f(ix,iy,iz)
  end do
  end do
  end do

  call timer_begin(LOG_ALLGATHERV_TOTAL)
  call comm_allgatherv(matbox11,matbox12,ircnt,idisp,info%icomm_ko)
  call timer_end(LOG_ALLGATHERV_TOTAL)

  if(process_allocation=='orbital_sequential')then
!$OMP parallel do private(i1,i2,i3,ibox,ibox2) collapse(3)
    do i3=0,nproc_d_g_dm(3)-1
    do i2=0,nproc_d_g_dm(2)-1
    do i1=0,nproc_d_g_dm(1)-1
      ibox=(nproc_id_global/nproc_d_g_mul_dm)*nproc_d_g_mul_dm    &
            +(i1+i2*nproc_d_g_dm(1)+i3*nproc_d_g_dm(1)*nproc_d_g_dm(2))
      ibox2=i1+i2*nproc_d_g_dm(1)+i3*nproc_d_g_dm(1)*nproc_d_g_dm(2)

      call copyVlocal(matbox12(idisp(ibox2):  &
                      (idisp(ibox2)+inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox))-1),  &
                      ibox,Vlocal(is)%f)

    end do
    end do
    end do
  else if(process_allocation=='grid_sequential')then
!$OMP parallel do private(i1,i2,i3,ibox,ibox2) collapse(3)
    do i3=0,nproc_d_g_dm(3)-1
    do i2=0,nproc_d_g_dm(2)-1
    do i1=0,nproc_d_g_dm(1)-1
      ibox=mod(nproc_id_global,nproc_d_o_mul)    &
          +(i1+i2*nproc_d_g_dm(1)+i3*nproc_d_g_dm(1)*nproc_d_g_dm(2))*nproc_d_o_mul
      ibox2=i1+i2*nproc_d_g_dm(1)+i3*nproc_d_g_dm(1)*nproc_d_g_dm(2)

      call copyVlocal(matbox12(idisp(ibox2):  &
                      (idisp(ibox2)+inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox))-1),  &
                      ibox,Vlocal(is)%f)

    end do
    end do
    end do
  end if

end do ! is=1,nspin

deallocate (ircnt,idisp)
deallocate (matbox11)
deallocate (matbox12)

CONTAINS
  subroutine copyVlocal(matbox,ibox,V)
    implicit none
    integer,intent(in) :: ibox
    real(8),intent(in) :: matbox(ista_Mxin_s(1,ibox):iend_Mxin_s(1,ibox),     &
                        ista_Mxin_s(2,ibox):iend_Mxin_s(2,ibox),     &
                        ista_Mxin_s(3,ibox):iend_Mxin_s(3,ibox))
    real(8) :: V(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
    V( ista_Mxin_s(1,ibox):iend_Mxin_s(1,ibox),     &
            ista_Mxin_s(2,ibox):iend_Mxin_s(2,ibox),     &
            ista_Mxin_s(3,ibox):iend_Mxin_s(3,ibox)) = &
             matbox(ista_Mxin_s(1,ibox):iend_Mxin_s(1,ibox),     &
                      ista_Mxin_s(2,ibox):iend_Mxin_s(2,ibox),     &
                      ista_Mxin_s(3,ibox):iend_Mxin_s(3,ibox))

    return
  end subroutine copyVlocal
end subroutine allgatherv_vlocal

subroutine wrapper_allgatherv_vlocal(ng,info) ! --> remove (future works)
  use structures
  implicit none
  type(s_rgrid),           intent(in) :: ng
  type(s_orbital_parallel),intent(in) :: info
  type(s_scalar) :: sVh,sVpsl
  type(s_scalar),allocatable :: sVlocal(:),sVxc(:)
  integer :: nspin,jspin

  nspin = 1
  if(ispin==1) nspin = 2
  allocate(sVlocal(nspin),sVxc(nspin))
  do jspin=1,nspin
    allocate(sVlocal(jspin)%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
    allocate(sVxc(jspin)%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  end do
  allocate(sVh%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  allocate(sVpsl%f(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
  sVpsl%f = Vpsl
  sVh%f = Vh
  if(ilsda == 1) then
    do jspin=1,nspin
      sVxc(jspin)%f = Vxc_s(:,:,:,jspin)
    end do
  else
    sVxc(1)%f = Vxc
  end if

  call allgatherv_vlocal(ng,info,nspin,sVh,sVpsl,sVxc,sVlocal)

  do jspin=1,nspin
    Vlocal(:,:,:,jspin) = sVlocal(jspin)%f
    call deallocate_scalar(sVxc(jspin))
    call deallocate_scalar(sVlocal(jspin))
  end do
  call deallocate_scalar(sVh)
  call deallocate_scalar(sVpsl)
  return
end subroutine wrapper_allgatherv_vlocal

END MODULE new_world_sub

