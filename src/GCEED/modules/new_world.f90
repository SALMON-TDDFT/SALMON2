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

integer :: npuy,npuz

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
call comm_get_groupinfo(info%icomm_o, info%id_o, info%isize_o)

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

! communicators for FFTE routine
  npuy = nproc_d_g_dm(2)*nproc_d_o(2)
  npuz = nproc_d_g_dm(3)*nproc_d_o(3)

  icolor=info_field%id(3)+info_field%id(1)*npuz
  ikey=info_field%id(2)
  info_field%icomm_ffte(2) = comm_create_group(nproc_group_global, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm_ffte(2), info_field%id_ffte(2), info_field%isize_ffte(2))

  icolor=info_field%id(2)+info_field%id(1)*npuy
  ikey=info_field%id(3)
  info_field%icomm_ffte(3) = comm_create_group(nproc_group_global, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm_ffte(3), info_field%id_ffte(3), info_field%isize_ffte(3))

  icolor=info_field%id(2)+info_field%id(3)*npuy
  ikey=info_field%id(1)
  info_field%icomm_ffte(1) = comm_create_group(nproc_group_global, icolor, ikey)
  call comm_get_groupinfo(info_field%icomm_ffte(1), info_field%id_ffte(1), info_field%isize_ffte(1))

  nproc_group_tdks = nproc_group_global
  nproc_id_tdks    = nproc_id_global
  nproc_size_tdks  = nproc_size_global

end subroutine make_new_world

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

