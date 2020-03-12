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
MODULE init_communicator
  implicit none

CONTAINS

!===================================================================================================================================

subroutine init_communicator_dft(comm,pinfo,info,info_field)
  use salmon_global, only: process_allocation
  use structures, only: s_orbital_parallel, s_field_parallel, s_process_info
  use communication, only: comm_create_group, comm_get_groupinfo, &
                           comm_is_root, comm_summation, comm_create_group_byid
  implicit none
  integer,      intent(in) :: comm
  type(s_process_info), intent(in) :: pinfo
  type(s_orbital_parallel) :: info
  type(s_field_parallel)   :: info_field
  !
  integer :: myrank,nproc
  integer :: nproc_k,nproc_ob
  integer :: nproc_d_o(3)
  integer :: i1,i2,i3,i4,i5,ix,iy,iz,nl
  integer,allocatable :: iranklists(:)

  call comm_get_groupinfo(comm, myrank, nproc)

  nproc_k   = pinfo%npk
  nproc_ob  = pinfo%nporbital
  nproc_d_o = pinfo%nprgrid

  allocate(iranklists(nproc))

! info
  info%icomm_rko = comm
  info%id_rko    = myrank
  info%isize_rko = nproc

  allocate(info%imap(0:nproc_d_o(1)-1, &
                     0:nproc_d_o(2)-1, &
                     0:nproc_d_o(3)-1, &
                     0:nproc_ob-1, &
                     0:nproc_k-1))

! communicator r,o,k,ro,ko
  nl = -1
  if (process_allocation == 'grid_sequential') then
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      nl = nl + 1
      info%imap(i1,i2,i3,i4,i5) = nl
      if (nl == myrank) then
        info%iaddress = [i1,i2,i3,i4,i5]
      end if
    end do
    end do
    end do
    end do
    end do
  else if (process_allocation == 'orbital_sequential') then
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      nl = nl + 1
      info%imap(i1,i2,i3,i4,i5) = nl
      if (nl == myrank) then
        info%iaddress = [i1,i2,i3,i4,i5]
      end if
    end do
    end do
    end do
    end do
    end do
  else
    stop 'undefined: process_allocation'
  end if

  if (nl /= info%isize_rko-1) &
    stop '[FATAL ERROR] init_communicator_dft'

  ! icomm_r
  i5 = info%iaddress(5)
  i4 = info%iaddress(4)
  nl = 0
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    nl = nl + 1
    iranklists(nl) = info%imap(i1,i2,i3,i4,i5)
  end do
  end do
  end do

  info%icomm_r = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_r, info%id_r, info%isize_r)

  ! icomm_o
  i5 = info%iaddress(5)
  i3 = info%iaddress(3)
  i2 = info%iaddress(2)
  i1 = info%iaddress(1)
  nl = 0
  do i4=0,nproc_ob-1
    nl = nl + 1
    iranklists(nl) = info%imap(i1,i2,i3,i4,i5)
  end do

  info%icomm_o = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_o, info%id_o, info%isize_o)

  ! icomm_k
  i4 = info%iaddress(4)
  i3 = info%iaddress(3)
  i2 = info%iaddress(2)
  i1 = info%iaddress(1)
  nl = 0
  do i5=0,nproc_k-1
    nl = nl + 1
    iranklists(nl) = info%imap(i1,i2,i3,i4,i5)
  end do

  info%icomm_k = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_k, info%id_k, info%isize_k)

  ! icomm_ro
  i5 = info%iaddress(5)
  nl = 0
  do i4=0,nproc_ob-1
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    nl = nl + 1
    iranklists(nl) = info%imap(i1,i2,i3,i4,i5)
  end do
  end do
  end do
  end do

  info%icomm_ro = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_ro, info%id_ro, info%isize_ro)

  ! icomm_ko
  i3 = info%iaddress(3)
  i2 = info%iaddress(2)
  i1 = info%iaddress(1)
  nl = 0
  do i5=0,nproc_k-1
  do i4=0,nproc_ob-1
    nl = nl + 1
    iranklists(nl) = info%imap(i1,i2,i3,i4,i5)
  end do
  end do

  info%icomm_ko = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_ko, info%id_ko, info%isize_ko)


! info_field
  info_field%icomm_all = info%icomm_r
  info_field%id_all    = info%id_r
  info_field%isize_all = info%isize_r

  allocate(info_field%imap(0:nproc_d_o(1)-1, &
                           0:nproc_d_o(2)-1, &
                           0:nproc_d_o(3)-1))

  i5 = info%iaddress(5)
  i4 = info%iaddress(4)
  info_field%imap(:,:,:)   = info%imap(:,:,:,i4,i5)
  info_field%iaddress(1:3) = info%iaddress(1:3)

! x-dir summation
  iz = info_field%iaddress(3)
  iy = info_field%iaddress(2)
  nl = 0
  do ix=0,nproc_d_o(1)-1
    nl = nl + 1
    iranklists(nl) = info_field%imap(ix,iy,iz)
  end do
  info_field%icomm(1) = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info_field%icomm(1), info_field%id(1), info_field%isize(1))

! y-dir summation
  iz = info_field%iaddress(3)
  ix = info_field%iaddress(1)
  nl = 0
  do iy=0,nproc_d_o(2)-1
    nl = nl + 1
    iranklists(nl) = info_field%imap(ix,iy,iz)
  end do
  info_field%icomm(2) = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info_field%icomm(2), info_field%id(2), info_field%isize(2))

! z-dir summation
  iy = info_field%iaddress(2)
  ix = info_field%iaddress(1)
  nl = 0
  do iz=0,nproc_d_o(3)-1
    nl = nl + 1
    iranklists(nl) = info_field%imap(ix,iy,iz)
  end do
  info_field%icomm(3) = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info_field%icomm(3), info_field%id(3), info_field%isize(3))

! xy-dir summation (for singlescale FDTD)
  iz = info_field%iaddress(3)
  nl = 0
  do iy=0,nproc_d_o(2)-1
  do ix=0,nproc_d_o(1)-1
    nl = nl + 1
    iranklists(nl) = info_field%imap(ix,iy,iz)
  end do
  end do
  info_field%icomm_xy = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info_field%icomm_xy, info_field%id_xy, info_field%isize_xy)

end subroutine init_communicator_dft

END MODULE init_communicator

