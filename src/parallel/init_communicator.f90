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

subroutine init_communicator_dft(comm,pinfo,info)
  use salmon_global, only: process_allocation
  use structures, only: s_process_info, s_parallel_info
  use communication, only: comm_create_group, comm_get_groupinfo, &
                           comm_is_root, comm_summation, comm_create_group_byid
  implicit none
  integer             ,intent(in) :: comm
  type(s_process_info),intent(in) :: pinfo
  type(s_parallel_info)           :: info
  !
  integer :: myrank,nproc
  integer :: nproc_k,nproc_ob
  integer :: nproc_d_o(3)
  integer :: i1,i2,i3,i4,i5,ix,iy,iz,nl
  integer,allocatable :: iranklists(:)

#ifdef __FUJITSU
  integer :: iret
#endif

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
#ifdef __FUJITSU
  call tofu_network_oriented_mapping(iret)
  if (iret < 0) then
    if (comm_is_root(info%id_rko)) then
      print *, 'Requested process shape does not match Tofu network shape.'
      print *, 'We fallback to usual mapping...'
    end if
#endif
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
#ifdef __FUJITSU
  end if
#endif

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

  i5 = info%iaddress(5)
  i4 = info%iaddress(4)

! x-dir summation
  iz = info%iaddress(3)
  iy = info%iaddress(2)
  nl = 0
  do ix=0,nproc_d_o(1)-1
    nl = nl + 1
    iranklists(nl) = info%imap(ix,iy,iz,i4,i5)
  end do
  info%icomm_x = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_x, info%id_x, info%isize_x)

! y-dir summation
  iz = info%iaddress(3)
  ix = info%iaddress(1)
  nl = 0
  do iy=0,nproc_d_o(2)-1
    nl = nl + 1
    iranklists(nl) = info%imap(ix,iy,iz,i4,i5)
  end do
  info%icomm_y = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_y, info%id_y, info%isize_y)

! z-dir summation
  iy = info%iaddress(2)
  ix = info%iaddress(1)
  nl = 0
  do iz=0,nproc_d_o(3)-1
    nl = nl + 1
    iranklists(nl) = info%imap(ix,iy,iz,i4,i5)
  end do
  info%icomm_z = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_z, info%id_z, info%isize_z)

! xy-dir summation (for singlescale FDTD)
  iz = info%iaddress(3)
  nl = 0
  do iy=0,nproc_d_o(2)-1
  do ix=0,nproc_d_o(1)-1
    nl = nl + 1
    iranklists(nl) = info%imap(ix,iy,iz,i4,i5)
  end do
  end do
  info%icomm_xy = comm_create_group_byid(comm, iranklists(1:nl))
  call comm_get_groupinfo(info%icomm_xy, info%id_xy, info%isize_xy)

#ifdef __FUJITSU
contains
  subroutine tofu_network_oriented_mapping(iret)
    use mpi_ext
    implicit none
    integer,intent(out) :: iret
    integer,parameter :: maxppn = 16
    integer :: n,pw(3)
    integer :: tofu_dim, tofu_shape(3), nprocs_per_node, icoords(3)
    integer :: pshape(3), iaddress(5)
    integer :: iranks(maxppn), outppn
    integer :: ierr

    iret = -1

    if (process_allocation /= 'grid_sequential') then
      print *, 'tofu_network_oriented_mapping: support grid_sequential only...'
      return
    end if

    call FJMPI_Topology_get_dimension(tofu_dim, ierr)
    if (ierr /= MPI_SUCCESS) &
      stop 'FJMPI_Topology_get_dimension: error'

    call FJMPI_Topology_get_shape(tofu_shape(1), tofu_shape(2), tofu_shape(3), ierr)
    if (ierr /= MPI_SUCCESS) &
      stop 'FJMPI_Topology_get_shape: error'

    call FJMPI_Topology_get_coords(info%icomm_rko, info%id_rko, FJMPI_LOGICAL, tofu_dim, icoords, ierr)
    if (ierr /= MPI_SUCCESS) &
      stop 'FJMPI_Topology_get_coords: error'

    nprocs_per_node = info%isize_rko / product(tofu_shape(1:tofu_dim))
    if (mod(info%isize_rko, product(tofu_shape(1:tofu_dim))) /= 0) &
      stop 'logical error: calcluate # of process/node'

    call FJMPI_Topology_get_ranks(info%icomm_rko, FJMPI_LOGICAL, icoords, maxppn, outppn, iranks, ierr)
    if (ierr /= MPI_SUCCESS) &
      stop 'FJMPI_Topology_get_ranks: error'
    if (outppn /= nprocs_per_node) &
      stop 'logical error: get # of process/node'

    select case(tofu_dim)

    case default
      if (comm_is_root(info%id_rko)) then
        print *, 'tofu_network_oriented_mapping: unsupported Tofu dimension,', tofu_dim
        return
      end if

    ! 3-dimensional network
    case(3)
      if (comm_is_root(info%id_rko)) then
        print *, '======================================================='
        print *, 'SALMON executes on 3-dimensional Tofu-network system...'
        print *, 'nproc_rgrid      :', nproc_d_o
        print *, 'nproc_ob*nproc_k :', nproc_ob*nproc_k
        print *, 'Tofu shape       :', tofu_shape
        print *, '# of process/node:', nprocs_per_node
        print *, ''
        print *, '[PROCESS MAPPING RULE]'
        print *, '  Requested network shape = (PX, PY, PZ, PW)'
        print *, '  Tofu-network shape      = (TX, TY, TZ)'
        print *, '  Process shape           = (TX*NPN, TY, TZ)'
        print *, ''
        print *, '  (PX, PY, PZ)      = nproc_rgrid'
        print *, '  PW  = PW1*PW2*PW3 = nproc_ob*nproc_k'
        print *, '  PW1 = TX*NPN / PX'
        print *, '  PW2 = TY     / PY'
        print *, '  PW3 = TZ     / PZ'
        print *, '  TX  = PX*PW1 / NPN'
        print *, '  TY  = PY*PW2'
        print *, '  TZ  = PZ*PW3'
        print *, '  NPN = # of process/node'
      end if

      pw(1) = tofu_shape(1) * nprocs_per_node / nproc_d_o(1)
      pw(2) = tofu_shape(2)                   / nproc_d_o(2)
      pw(3) = tofu_shape(3)                   / nproc_d_o(3)

      pshape(1:3) = nproc_d_o(1:3) * pw(1:3)

      if (comm_is_root(info%id_rko)) then
        print *, '  (PW1, PW2, PW3) =', pw(1:3)
        print *, '  (PX*PW1, PY*PW2, PZ*PW3) =', pshape(1:3)
        print *, '======================================================='
      end if

      if (product(pw) /= nproc_ob*nproc_k) then
        if (comm_is_root(info%id_rko)) then
          print *, 'product(pw) /= nproc_ob*nproc_k'
          print *, '(PW1, PW2, PW3) =', pw
        end if
        return
      end if

      info%imap = -1
      nl = -1
      do iz=0,pshape(3)-1
      do iy=0,pshape(2)-1
      do ix=0,pshape(1)-1
        nl = nl + 1

        iaddress(1) = mod(ix, nproc_d_o(1))
        iaddress(2) = mod(iy, nproc_d_o(2))
        iaddress(3) = mod(iz, nproc_d_o(3))

        n = ((ix + nproc_d_o(1)) / nproc_d_o(1) - 1) &
          + ((iy + nproc_d_o(2)) / nproc_d_o(2) - 1) * pw(1) &
          + ((iz + nproc_d_o(3)) / nproc_d_o(3) - 1) * pw(1) * pw(2)
        iaddress(4) = mod(n, nproc_ob)
        iaddress(5) = n / nproc_ob

        if (info%imap(iaddress(1),iaddress(2),iaddress(3),iaddress(4),iaddress(5)) /= -1) then
          if (comm_is_root(info%id_rko)) then
            print *, nl,'conflict:',iaddress,' pad =',ix,iy,iz
          end if
          return
        end if

        info%imap(iaddress(1),iaddress(2),iaddress(3),iaddress(4),iaddress(5)) = nl
        if (nl == myrank) then
          info%iaddress = iaddress
        end if
      end do
      end do
      end do

      if (minval(info%imap) < 0 .or. maxval(info%imap) >= info%isize_rko) then
        if (comm_is_root(info%id_rko)) then
          print *, 'tofu_network_oriented_mapping: logical error.'
        end if
        return
      end if

      if (nl /= info%isize_rko-1) then
        if (comm_is_root(info%id_rko)) then
          print *, 'tofu_network_oriented_mapping: Tofu process mapping.'
        end if
        return
      end if

    end select

    iret = 0
    return
  end subroutine tofu_network_oriented_mapping
#endif
end subroutine init_communicator_dft

END MODULE init_communicator

