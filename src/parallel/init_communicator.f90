!
!  Copyright 2019-2020 SALMON developers
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

subroutine init_communicator_dft(comm,info)
  use salmon_global, only: process_allocation, method_poisson, yn_ffte
  use structures, only: s_parallel_info
  use communication, only: comm_create_group, comm_get_groupinfo, &
                           comm_is_root, comm_summation, comm_create_group_byid
  implicit none
  integer             ,intent(in) :: comm
  type(s_parallel_info)           :: info
  !
  integer :: myrank,nproc
  integer :: nproc_k,nproc_ob
  integer :: nproc_d_o(3)
  integer :: i1,i2,i3,i4,i5,ix,iy,iz,nl,io1,io2,io3,io4
  integer,allocatable :: iranklists(:)

#ifdef __FUJITSU
  integer :: iret
#endif

  call comm_get_groupinfo(comm, myrank, nproc)

  nproc_k   = info%npk
  nproc_ob  = info%nporbital
  nproc_d_o = info%nprgrid

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

! begin: preparation of communicator for isolated_ffte
  if(method_poisson=='ft'.and.yn_ffte=='y')then

! condition of parallel calulations
    if(mod(nproc_ob,4)==0) then

! begin: allocation of imap_isolated_ffte
      allocate(info%imap_isolated_ffte(0:nproc_d_o(1)-1, &
                                       0:nproc_d_o(2)-1, &
                                       0:nproc_d_o(3)-1, &
                                       0:1, &
                                       0:1, &
                                       0:nproc_ob/4-1, &
                                       0:nproc_k-1))
! end: allocation of imap_isolated_ffte
    
! begin: definition of iaddress_isolated_ffte
      nl = -1
      if (process_allocation == 'grid_sequential') then
        do i5=0,nproc_k-1
        do io4=0,nproc_ob/4-1
        do io3=0,1 ! z-dir
        do io2=0,1 ! y-dir
        do i3=0,nproc_d_o(3)-1
        do i2=0,nproc_d_o(2)-1
        do i1=0,nproc_d_o(1)-1
          nl = nl + 1
          info%imap_isolated_ffte(i1,i2,i3,io2,io3,io4,i5) = nl
          if (nl == myrank) then
            info%iaddress_isolated_ffte = [i1,i2,i3,io2,io3,io4,i5]
          end if
        end do
        end do
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
        do io4=0,nproc_ob/4-1
        do io3=0,1 ! z-dir
        do io2=0,1 ! y-dir
          nl = nl + 1
          info%imap_isolated_ffte(i1,i2,i3,io2,io3,io4,i5) = nl
          if (nl == myrank) then
            info%iaddress_isolated_ffte = [i1,i2,i3,io2,io3,io4,i5]
          end if
        end do
        end do
        end do
        end do
        end do
        end do
        end do
      else
        stop 'undefined: process_allocation'
      end if
! end: definition of iaddress_isolated_ffte

! begin: definition of communicator
! y-dir 
      i5 = info%iaddress_isolated_ffte(7)
      io4 = info%iaddress_isolated_ffte(6)
      io3 = info%iaddress_isolated_ffte(5)
      iz = info%iaddress_isolated_ffte(3)
      ix = info%iaddress_isolated_ffte(1)
      nl = 0
      do io2=0,1
      do iy=0,nproc_d_o(2)-1
        nl = nl + 1
        iranklists(nl) = info%imap_isolated_ffte(ix,iy,iz,io2,io3,io4,i5)
      end do
      end do
      info%icomm_y_isolated_ffte = comm_create_group_byid(comm, iranklists(1:nl))
      call comm_get_groupinfo(info%icomm_y_isolated_ffte, info%id_y_isolated_ffte, info%isize_y_isolated_ffte)

! z-dir 
      i5 = info%iaddress_isolated_ffte(7)
      io4 = info%iaddress_isolated_ffte(6)
      io2 = info%iaddress_isolated_ffte(4)
      iy = info%iaddress_isolated_ffte(2)
      ix = info%iaddress_isolated_ffte(1)
      nl = 0
      do io3=0,1
      do iz=0,nproc_d_o(3)-1
        nl = nl + 1
        iranklists(nl) = info%imap_isolated_ffte(ix,iy,iz,io2,io3,io4,i5)
      end do
      end do
      info%icomm_z_isolated_ffte = comm_create_group_byid(comm, iranklists(1:nl))
      call comm_get_groupinfo(info%icomm_z_isolated_ffte, info%id_z_isolated_ffte, info%isize_z_isolated_ffte)
    
! orbital-dir 
      i5 = info%iaddress_isolated_ffte(7)
      io4 = info%iaddress_isolated_ffte(6)
      iz = info%iaddress_isolated_ffte(3)
      iy = info%iaddress_isolated_ffte(2)
      ix = info%iaddress_isolated_ffte(1)
      nl = 0
      do io3=0,1
      do io2=0,1
        nl = nl + 1
        iranklists(nl) = info%imap_isolated_ffte(ix,iy,iz,io2,io3,io4,i5)
      end do
      end do
      info%icomm_o_isolated_ffte = comm_create_group_byid(comm, iranklists(1:nl))
      call comm_get_groupinfo(info%icomm_o_isolated_ffte, info%id_o_isolated_ffte, info%isize_o_isolated_ffte)

! condition of single calulations
    else
! x-dir, y-dir, z-dir, o-dir
      i5 = info%iaddress(5)
      i4 = info%iaddress(4)
      iz = info%iaddress(3)
      iy = info%iaddress(2)
      ix = info%iaddress(1)
      nl = 1
      iranklists(nl) = info%imap(ix,iy,iz,i4,i5)
      info%icomm_x_isolated_ffte = comm_create_group_byid(comm, iranklists(1:nl))
      call comm_get_groupinfo(info%icomm_x_isolated_ffte, info%id_x_isolated_ffte, info%isize_x_isolated_ffte)
      info%icomm_y_isolated_ffte = comm_create_group_byid(comm, iranklists(1:nl))
      call comm_get_groupinfo(info%icomm_y_isolated_ffte, info%id_y_isolated_ffte, info%isize_y_isolated_ffte)
      info%icomm_z_isolated_ffte = comm_create_group_byid(comm, iranklists(1:nl))
      call comm_get_groupinfo(info%icomm_z_isolated_ffte, info%id_z_isolated_ffte, info%isize_z_isolated_ffte)
      info%icomm_o_isolated_ffte = comm_create_group_byid(comm, iranklists(1:nl))
      call comm_get_groupinfo(info%icomm_o_isolated_ffte, info%id_o_isolated_ffte, info%isize_o_isolated_ffte)
    end if
  end if
! end: preparation of communicator for isolated_ffte

#ifdef __FUJITSU
contains
  subroutine tofu_network_oriented_mapping(iret)
    use mpi_ext
    use salmon_global, only: nx_m,ny_m,nz_m, theory
    implicit none
    integer,intent(out) :: iret
    integer,parameter :: maxppn = 16
    integer :: n,pw(3)
    integer :: tofu_dim, tofu_shape(3), nprocs_per_node, icoords(3)
    integer :: pshape(3), iaddress(5)
    integer :: iranks(maxppn), outppn
    integer :: ierr

    iret = -1

    call FJMPI_Topology_get_dimension(tofu_dim, ierr)
    if (ierr /= MPI_SUCCESS) &
      stop 'FJMPI_Topology_get_dimension: error'

    call FJMPI_Topology_get_shape(tofu_shape(1), tofu_shape(2), tofu_shape(3), ierr)
    if (ierr /= MPI_SUCCESS) &
      stop 'FJMPI_Topology_get_shape: error'

    call FJMPI_Topology_get_coords(info%icomm_rko, info%id_rko, FJMPI_LOGICAL, tofu_dim, icoords, ierr)
    if (ierr /= MPI_SUCCESS) &
      stop 'FJMPI_Topology_get_coords: error'


    if(theory=='multi_scale_maxwell_tddft') then
       nprocs_per_node = info%isize_rko *nx_m*ny_m*nz_m / product(tofu_shape(1:tofu_dim))
      if (mod(info%isize_rko*nx_m*ny_m*nz_m, product(tofu_shape(1:tofu_dim))) /= 0) &
          stop 'logical error (ms): calcluate # of process/node'
    else
       nprocs_per_node = info%isize_rko / product(tofu_shape(1:tofu_dim))
      if (mod(info%isize_rko, product(tofu_shape(1:tofu_dim))) /= 0) &
          stop 'logical error: calcluate # of process/node'
    endif

    call FJMPI_Topology_get_ranks(info%icomm_rko, FJMPI_LOGICAL, icoords, maxppn, outppn, iranks, ierr)
    if (ierr /= MPI_SUCCESS) &
      stop 'FJMPI_Topology_get_ranks: error'
    if (outppn /= nprocs_per_node) &
      stop 'logical error: get # of process/node'

    select case(tofu_dim)

    case default
      if (comm_is_root(info%id_rko)) then
        print *, 'tofu_network_oriented_mapping: unsupported Tofu dimension,', tofu_dim
      end if
      return

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
        print *, '  Process shape           = (TX*PPN, TY, TZ)'
        print *, ''
        print *, '  (PX, PY, PZ) = nproc_rgrid'
        print *, '  PW           = nproc_ob*nproc_k'
        print *, '  PPN          = # of process/node'
        if (process_allocation == 'grid_sequential') then
          print *, '  PW  = PW1*PW2*PW3'
          print *, '  PW1 = TX*PPN / PX'
          print *, '  PW2 = TY     / PY'
          print *, '  PW3 = TZ     / PZ'
          print *, '  TX  = PX*PW1 / PPN'
          print *, '  TY  = PY*PW2'
          print *, '  TZ  = PZ*PW3'
        else if (process_allocation == 'orbital_sequential') then
          print *, '  PX  = PX1*PX2*PX3'
          print *, '  PX1 = TX*PPN / PW'
          print *, '  PX2 = TY     / PY'
          print *, '  PX3 = TZ     / PZ'
          print *, '  TX  = PW*PX1 / PPN'
          print *, '  TY  = PY*PX2'
          print *, '  TZ  = PZ*PX3'
        end if
      end if

      if (process_allocation == 'grid_sequential') then
        pw(1) = tofu_shape(1) * nprocs_per_node / nproc_d_o(1)
      else if (process_allocation == 'orbital_sequential') then
        pw(1) = tofu_shape(1) * nprocs_per_node / (nproc_ob*nproc_k)
      end if
      pw(2) = tofu_shape(2) / nproc_d_o(2)
      pw(3) = tofu_shape(3) / nproc_d_o(3)

      if (process_allocation == 'grid_sequential') then
        pshape(1:3) = nproc_d_o(1:3) * pw(1:3)
      else if (process_allocation == 'orbital_sequential') then
        pshape(1)   = nproc_ob*nproc_k * pw(1)
        pshape(2:3) = nproc_d_o(2:3) * pw(2:3)
      end if

      if (comm_is_root(info%id_rko)) then
        if (process_allocation == 'grid_sequential') then
          print *, '  (PW1, PW2, PW3) =', pw(1:3)
          print *, '  (PX*PW1, PY*PW2, PZ*PW3) =', pshape(1:3)
        else if (process_allocation == 'orbital_sequential') then
          print *, '  (PX1, PX2, PX3) =', pw(1:3)
          print *, '  (PW*PX1, PY*PX2, PZ*PX3) =', pshape(1:3)
        end if
        print *, '======================================================='
      end if

      if (process_allocation == 'grid_sequential') then
        if (product(pw) /= nproc_ob*nproc_k) then
          if (comm_is_root(info%id_rko)) then
            print *, 'product(PW[1-3]) /= PW'
            print *, '(PW1, PW2, PW3) =', pw
          end if
          return
        end if
      else if (process_allocation == 'orbital_sequential') then
        if (product(pw) /= nproc_d_o(1)) then
          if (comm_is_root(info%id_rko)) then
            print *, 'product(PX[1-3]) /= PX'
            print *, '(PX1, PX2, PX3) =', pw
          end if
          return
        end if
      end if

      info%imap = -1
      nl = -1
      do iz=0,pshape(3)-1
      do iy=0,pshape(2)-1
      do ix=0,pshape(1)-1
        nl = nl + 1

        iaddress(2) = mod(iy, nproc_d_o(2))
        iaddress(3) = mod(iz, nproc_d_o(3))

        if (process_allocation == 'grid_sequential') then
          iaddress(1) = mod(ix, nproc_d_o(1))

          n = ((ix + nproc_d_o(1)) / nproc_d_o(1) - 1) &
            + ((iy + nproc_d_o(2)) / nproc_d_o(2) - 1) * pw(1) &
            + ((iz + nproc_d_o(3)) / nproc_d_o(3) - 1) * pw(1) * pw(2)

          iaddress(4) = mod(n, nproc_ob)
          iaddress(5) = n / nproc_ob
        else if (process_allocation == 'orbital_sequential') then
          n = ((ix + nproc_ob*nproc_k) / nproc_ob*nproc_k - 1) &
            + ((iy + nproc_d_o(2))     / nproc_d_o(2)     - 1) * pw(1) &
            + ((iz + nproc_d_o(3))     / nproc_d_o(3)     - 1) * pw(1) * pw(2)
          iaddress(1) = mod(n, nproc_d_o(1))

          n = mod(ix, nproc_ob*nproc_k)
          iaddress(4) = mod(n, nproc_ob)
          iaddress(5) = n / nproc_ob
        end if

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

