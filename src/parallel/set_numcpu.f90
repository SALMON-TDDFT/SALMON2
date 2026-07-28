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
module set_numcpu
  implicit none

  integer,parameter,public :: iprefer_k_distribution        =  1
  integer,parameter,public :: iprefer_orbital_distribution  =  2
  integer,parameter,public :: iprefer_domain_distribution   =  3

contains

  function check_numcpu(icomm1, info) result(iok)
  use structures, only: s_parallel_info
  use communication, only: comm_get_groupinfo
  implicit none
  integer, intent(in) :: icomm1 ! Communicator for a single grid system
  type(s_parallel_info),intent(inout) :: info
  logical :: iok

  integer :: nproc_k, nproc_ob
  integer,dimension(3) :: nproc_d_o
  integer :: nproc_total_wf
  integer :: nproc_id_comm1, nproc_size_comm1

  call comm_get_groupinfo(icomm1, nproc_id_comm1, nproc_size_comm1)

  nproc_k   = info%npk
  nproc_ob  = info%nporbital
  nproc_d_o = info%nprgrid

  nproc_total_wf = nproc_k * nproc_ob * product(nproc_d_o(:))

  iok = .true.

  if(nproc_total_wf/=nproc_size_comm1)then
    print "(A)",   'number of MPI process is not correct'
    print "(A,I7)", '  nproc_k * nproc_ob * product(nproc_domain_orbital) =', nproc_total_wf
    print "(A,I7)", '  number of MPI process                              =', nproc_size_comm1
    iok = .false.
  end if
end function check_numcpu

subroutine set_numcpu_general(iprefer_dist,numk,numo,icomm1,info,rspace_allowed)
  use structures, only: s_parallel_info
  use communication, only: comm_get_groupinfo
  implicit none
  integer,intent(in)                :: iprefer_dist,numk,numo,icomm1
  type(s_parallel_info),intent(out) :: info
  ! rspace_allowed=.false. forbids real-space (domain) parallelization, used for
  ! non-orthogonal lattices. Optional; defaults to .true. (previous behaviour).
  logical,intent(in),optional       :: rspace_allowed
  logical :: allow_rspace
  logical :: if_found
  integer :: nproc_size_comm1, nproc_id_comm1

  integer :: ip
  integer :: nproc_k,nproc_ob
  integer :: nproc_d_o(3)
  integer :: nproc_size_comm1_tmp

  integer :: ii,icount
  integer :: num_factor2
  integer :: num_factor3
  integer :: num_factor5

  integer :: nk,no

  allow_rspace = .true.
  if (present(rspace_allowed)) allow_rspace = rspace_allowed

  call comm_get_groupinfo(icomm1, nproc_id_comm1, nproc_size_comm1)
  nproc_size_comm1_tmp=nproc_size_comm1

  nk = max(numk,1)
  no = max(numo,1)

  nproc_k   = 1
  nproc_ob  = 1
  nproc_d_o = 1

  if (.not. allow_rspace) then
    ! Non-orthogonal lattice: r-space (domain) parallelization is unavailable, so
    ! every rank must be placed in k x orbital only (nprgrid = 1). The greedy loop
    ! below can leave an unplaceable residual and abort even when a valid k x orbital
    ! split exists, so search the rank-count divisors for one (preferring a clean
    ! division of nk / no and the requested distribution axis) before giving up.
    call set_kob_only(iprefer_dist, nproc_size_comm1, nk, no, nproc_k, nproc_ob, if_found)
    if (.not. if_found) then
      stop "automatic process distribution: a non-orthogonal lattice cannot use &
           &r-space parallelization, and no nproc_k*nproc_ob split fits the rank &
           &count; set nproc_k and nproc_ob explicitly so their product equals the &
           &number of MPI processes."
    end if
    nproc_d_o = 1
    info%npk          = nproc_k
    info%nporbital    = nproc_ob
    info%nprgrid(1:3) = nproc_d_o(1:3)
    return
  end if

  do ip=iprefer_dist,iprefer_domain_distribution

    select case(ip)

    ! k-point
    case(iprefer_k_distribution)
      if (1 < nk .and. nk < nproc_size_comm1_tmp) then
        if (mod(nproc_size_comm1_tmp,nk) == 0) then
          nproc_k               = nproc_size_comm1_tmp / nk
          nproc_size_comm1_tmp = nproc_size_comm1_tmp / nproc_k
        end if
      else
        nproc_k               = min(nproc_size_comm1_tmp,nk)
        nproc_size_comm1_tmp = nproc_size_comm1_tmp / nproc_k
      end if

    ! orbital
    case(iprefer_orbital_distribution)
      if (1 < no .and. no < nproc_size_comm1_tmp) then
        if (mod(nproc_size_comm1_tmp,no) == 0) then
          nproc_ob              = nproc_size_comm1_tmp / no
          nproc_size_comm1_tmp = nproc_size_comm1_tmp / nproc_ob
        end if
      else
        nproc_ob              = min(nproc_size_comm1_tmp,no)
        nproc_size_comm1_tmp = nproc_size_comm1_tmp / nproc_ob
      end if

    ! rgrid
    case(iprefer_domain_distribution)
      if(.not. allow_rspace)then
        ! Non-orthogonal lattice: real-space (domain) parallelization is not
        ! supported. Any processes left after k/orbital distribution cannot be
        ! placed without it, so this is a user configuration error.
        if(nproc_size_comm1_tmp /= 1) &
          stop "automatic process distribution: a non-orthogonal lattice cannot use &
               &r-space parallelization; set nproc_k*nproc_ob = (number of MPI processes)."
        nproc_d_o = 1
        cycle
      end if
      num_factor2=0
      do ii=1,26
        if(mod(nproc_size_comm1_tmp,2)==0)then
          num_factor2=num_factor2+1
          nproc_size_comm1_tmp=nproc_size_comm1_tmp/2
        end if
      end do

      num_factor3=0
      do ii=1,17
        if(mod(nproc_size_comm1_tmp,3)==0)then
          num_factor3=num_factor3+1
          nproc_size_comm1_tmp=nproc_size_comm1_tmp/3
        end if
      end do

      num_factor5=0
      do ii=1,11
        if(mod(nproc_size_comm1_tmp,5)==0)then
          num_factor5=num_factor5+1
          nproc_size_comm1_tmp=nproc_size_comm1_tmp/5
        end if
      end do

      if(nproc_size_comm1_tmp/=1)then
        stop "In automatic process distribution, prime factors for number of processes must be combination of 2, 3 or 5."
      end if

      nproc_d_o = 1
      icount = 0

      do ii=1,num_factor5
        icount=icount+1
        if(mod(icount,3)==1)then
          nproc_d_o(3)=nproc_d_o(3)*5
        else if(mod(icount,3)==2)then
          nproc_d_o(2)=nproc_d_o(2)*5
        else
          nproc_d_o(1)=nproc_d_o(1)*5
        end if
      end do

      do ii=1,num_factor3
        icount=icount+1
        if(mod(icount,3)==1)then
          nproc_d_o(3)=nproc_d_o(3)*3
        else if(mod(icount,3)==2)then
          nproc_d_o(2)=nproc_d_o(2)*3
        else
          nproc_d_o(1)=nproc_d_o(1)*3
        end if
      end do

      do ii=1,num_factor2
        icount=icount+1
        if(mod(icount,3)==1)then
          nproc_d_o(3)=nproc_d_o(3)*2
        else if(mod(icount,3)==2)then
          nproc_d_o(2)=nproc_d_o(2)*2
        else
          nproc_d_o(1)=nproc_d_o(1)*2
        end if
      end do

    end select
  end do

  info%npk         = nproc_k
  info%nporbital   = nproc_ob
  info%nprgrid(1:3) = nproc_d_o(1:3)

end subroutine set_numcpu_general

subroutine set_kob_only(iprefer_dist, nproc, nk, no, nproc_k, nproc_ob, if_found)
  ! Distribute `nproc` ranks over k x orbital only (no r-space parallelization,
  ! for non-orthogonal lattices): find nproc_k*nproc_ob = nproc with nproc_k<=nk
  ! and nproc_ob<=no. Pass 1 also requires nk and no to divide evenly; pass 2
  ! relaxes that to the <= bounds. Among candidates the requested axis (k or
  ! orbital) is maximized. if_found=.false. if no such split exists.
  implicit none
  integer,intent(in)  :: iprefer_dist, nproc, nk, no
  integer,intent(out) :: nproc_k, nproc_ob
  logical,intent(out) :: if_found
  integer :: p, pk, po, pass
  logical :: prefer_k, clean, ok

  prefer_k = (iprefer_dist == iprefer_k_distribution)
  if_found = .false.
  nproc_k  = 1
  nproc_ob = 1

  do pass = 1, 2
    clean = (pass == 1)
    do p = 1, nproc
      if (mod(nproc,p) /= 0) cycle
      if (prefer_k) then
        pk = p ;          po = nproc/p
      else
        po = p ;          pk = nproc/p
      end if
      if (pk > nk .or. po > no) cycle
      if (clean) then
        if (mod(nk,pk) /= 0 .or. mod(no,po) /= 0) cycle
      end if
      if (prefer_k) then
        ok = (.not. if_found) .or. (pk > nproc_k)
      else
        ok = (.not. if_found) .or. (po > nproc_ob)
      end if
      if (ok) then
        nproc_k  = pk
        nproc_ob = po
        if_found = .true.
      end if
    end do
    if (if_found) exit
  end do
end subroutine set_kob_only

end module set_numcpu
