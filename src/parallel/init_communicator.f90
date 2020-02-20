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
  use salmon_global, only: ispin,process_allocation
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
  integer :: nproc_d_o(3),nproc_d_g(3),nproc_d_o_mul,nproc_d_g_dm(3),nproc_d_g_mul_dm,nproc_ob_spin(2)
  integer :: i1,i2,i3,i4,i5,ix,iy,iz,ixs,iys,izs,nsize,ilocal_rank
  integer :: ibox
  integer :: mesh(3),network(3),mesh_step(3)
  integer :: icolor_r,icolor_o,icolor_k,icolor_ro,icolor_ko
  integer :: ikey_ro,ikey_ko,nl
  integer,allocatable :: iranklists(:)

  call comm_get_groupinfo(comm, myrank, nproc)

  nproc_k          = pinfo%npk
  nproc_ob         = pinfo%nporbital
  nproc_ob_spin    = pinfo%nporbital_spin
  nproc_d_o        = pinfo%npdomain_orbital
  nproc_d_g        = pinfo%npdomain_general
  nproc_d_g_dm     = pinfo%npdomain_general_dm
  nproc_d_o_mul    = product(nproc_d_o)
  nproc_d_g_mul_dm = product(nproc_d_g_dm)

  call check_network_form


! info
  info%icomm_rko = comm
  info%id_rko    = myrank
  info%isize_rko = nproc

! communicator r,o,k,ro,ko
  do i5=0,nproc_k-1
  do i4=0,nproc_ob-1
  do i3=0,nproc_d_o(3)-1
  do i2=0,nproc_d_o(2)-1
  do i1=0,nproc_d_o(1)-1
    if (process_allocation=='grid_sequential') then
      ibox = (i5 * nproc_ob + i4) * product(nproc_d_o) &
           + (i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1)
    else if (process_allocation=='orbital_sequential') then
      ibox = (i5 * nproc_ob + i4) &
           + (i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1) * nproc_ob * nproc_k
    end if

    if(myrank == ibox)then
      info%iaddress = [i1, i2, i3, i4, i5]

      ! icomm_r
      icolor_r = i5 * nproc_ob + i4

      ! icomm_o
      icolor_o = i5 &
               + (i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1) * nproc_k

      ! icomm_k
      icolor_k = i4 * product(nproc_d_o) &
               + (i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1)

      ! icomm_ro, ko
      if(ispin==0)then
        icolor_ro = i5
        ikey_ro   = i4 &
                  + (i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1) * nproc_ob
        ikey_ko   = i5 * nproc_ob + i4
      else if(i4<nproc_ob_spin(1))then
        icolor_ro = 2 * i5 + 0
        ikey_ro   = i4 &
                  + (i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1) * nproc_ob_spin(1)
        ikey_ko   = i5 * nproc_ob_spin(1) + i4
      else
        icolor_ro = 2 * i5 + 1
        ikey_ro   = i4 - nproc_ob_spin(1) &
                  + (i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1) * nproc_ob_spin(2)
        ikey_ko   = i5 * nproc_ob_spin(2) + i4 - nproc_ob_spin(1)
      end if
      icolor_ko = (i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1)
    end if
  end do
  end do
  end do
  end do
  end do

  info%icomm_r = comm_create_group(comm, icolor_r, myrank)
  call comm_get_groupinfo(info%icomm_r, info%id_r, info%isize_r)

  info%icomm_o = comm_create_group(comm, icolor_o, myrank)
  call comm_get_groupinfo(info%icomm_o, info%id_o, info%isize_o)

  info%icomm_k = comm_create_group(comm, icolor_k, myrank)
  call comm_get_groupinfo(info%icomm_k, info%id_k, info%isize_k)

  info%icomm_ro = comm_create_group(comm, icolor_ro, ikey_ro)
  call comm_get_groupinfo(info%icomm_ro, info%id_ro, info%isize_ro)

  info%icomm_ko = comm_create_group(comm, icolor_ko, ikey_ko)
  call comm_get_groupinfo(info%icomm_ko, info%id_ko, info%isize_ko)

  allocate(info%imap(0:nproc_d_o(1)-1, &
                     0:nproc_d_o(2)-1, &
                     0:nproc_d_o(3)-1, &
                     0:nproc_ob-1, &
                     0:nproc_k-1))
  info%imap = 0
  info%imap(info%iaddress(1), &
            info%iaddress(2), &
            info%iaddress(3), &
            info%iaddress(4), &
            info%iaddress(5)) = myrank
  call comm_summation(info%imap, comm)


! info_field
  info_field%icomm_all = comm
  info_field%id_all    = myrank
  info_field%isize_all = nproc

  allocate(info_field%imap(0:nproc_d_g(1)-1, &
                           0:nproc_d_g(2)-1, &
                           0:nproc_d_g(3)-1))

  allocate(iranklists(product(nproc_d_g)))

  ! packing same id_r
  ilocal_rank = 0
  do iz=0,nproc_d_g(3)-1,nproc_d_g_dm(3)
  do iy=0,nproc_d_g(2)-1,nproc_d_g_dm(2)
  do ix=0,nproc_d_g(1)-1,nproc_d_g_dm(1)
    ! search same id_r...
    nsize = 0
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
    do i3=0,nproc_d_o(3)-1
    do i2=0,nproc_d_o(2)-1
    do i1=0,nproc_d_o(1)-1
      ! `process_allocation` is not related (maybe...)
      ibox = i3 * nproc_d_o(1) * nproc_d_o(2) + i2 * nproc_d_o(1) + i1
      if (ibox == ilocal_rank) then
        nsize = nsize + 1
        iranklists(nsize) = info%imap(i1,i2,i3,i4,i5)
      end if
    end do
    end do
    end do
    end do
    end do

    nsize = 0
    do izs=iz,iz+nproc_d_g_dm(3)-1
    do iys=iy,iy+nproc_d_g_dm(2)-1
    do ixs=ix,ix+nproc_d_g_dm(1)-1
      nsize = nsize + 1
      info_field%imap(ixs,iys,izs) = iranklists(nsize)
      if (iranklists(nsize) == myrank) then
        info_field%iaddress = [ixs,iys,izs]
      end if
    end do
    end do
    end do
    ilocal_rank = ilocal_rank + 1
  end do
  end do
  end do

! x-dir summation
  iz = info_field%iaddress(3)
  iy = info_field%iaddress(2)
  nl = 0
  do ix=0,nproc_d_g(1)-1
    nl = nl + 1
    iranklists(nl) = info_field%imap(ix,iy,iz)
  end do
  info_field%icomm(1) = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info_field%icomm(1), info_field%id(1), info_field%isize(1))

! y-dir summation
  iz = info_field%iaddress(3)
  ix = info_field%iaddress(1)
  nl = 0
  do iy=0,nproc_d_g(2)-1
    nl = nl + 1
    iranklists(nl) = info_field%imap(ix,iy,iz)
  end do
  info_field%icomm(2) = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info_field%icomm(2), info_field%id(2), info_field%isize(2))

! z-dir summation
  iy = info_field%iaddress(2)
  ix = info_field%iaddress(1)
  nl = 0
  do iz=0,nproc_d_g(3)-1
    nl = nl + 1
    iranklists(nl) = info_field%imap(ix,iy,iz)
  end do
  info_field%icomm(3) = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info_field%icomm(3), info_field%id(3), info_field%isize(3))


! for allgatherv_vlocal
  info_field%icomm_v  = info%icomm_ko
  info_field%ngo(1:3) = nproc_d_g_dm(1:3)
  info_field%ngo_xyz  = nproc_d_g_mul_dm
  info_field%nproc_o  = nproc_d_o_mul


! communicators for FFTE routine
  info_field%icomm_ffte(1:3) = info_field%icomm(1:3)

  call comm_get_groupinfo(info_field%icomm_ffte(1), info_field%id_ffte(1), info_field%isize_ffte(1))
  call comm_get_groupinfo(info_field%icomm_ffte(2), info_field%id_ffte(2), info_field%isize_ffte(2))
  call comm_get_groupinfo(info_field%icomm_ffte(3), info_field%id_ffte(3), info_field%isize_ffte(3))

contains
  subroutine check_network_form
    implicit none

    mesh(:)      = nproc_d_o(:)
    mesh_step(:) = nproc_d_g(:) / nproc_d_o(:)
    network(:)   = mesh_step(:) * nproc_d_o(:)

    if (comm_is_root(myrank)) then
      print *, '================ MPI-Info ================='
      print *, 'we assume the network is fat-tree.'
      print *, 'we allocate the process as grid sequential...'

      print '(A,5I4)', ' requested network (5D):',mesh,nproc_ob,nproc_k
      print '(A,3I4)', ' assumed   network (3D):',network

      if (product(mesh) * nproc_ob * nproc_k /= nproc) then
        print *, 'invalid requested network' ; stop
      end if

      if (product(network) /= nproc) then
        print *, 'invalid assumed network' ; stop
      end if
    end if
  end subroutine check_network_form
end subroutine init_communicator_dft

END MODULE init_communicator

