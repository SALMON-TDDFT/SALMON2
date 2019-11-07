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
module set_numcpu
  implicit none

  integer,parameter,public :: iprefer_k_distribution        =  1
  integer,parameter,public :: iprefer_orbital_distribution  =  2
  integer,parameter,public :: iprefer_domain_distribution   =  3

contains

  function check_numcpu(pinfo) result(iok)
  use structures, only: s_process_info
  use salmon_parallel, only: nproc_size_global
  implicit none
  type(s_process_info),intent(inout) :: pinfo
  logical :: iok

  integer :: j
  integer :: nproc_k, nproc_ob
  integer,dimension(3) :: nproc_d_o,nproc_d_g
  integer :: nproc_total_wf,nproc_total_g

  nproc_k   = pinfo%npk
  nproc_ob  = pinfo%nporbital
  nproc_d_o = pinfo%npdomain_orbital
  nproc_d_g = pinfo%npdomain_general

  nproc_total_wf = nproc_k * nproc_ob * product(nproc_d_o(:))
  nproc_total_g  = product(nproc_d_g(:))

  iok = .true.

  if(nproc_total_wf/=nproc_size_global)then
    print "(A)",   'number of MPI process is not correct'
    print "(A,I7)", '  nproc_k * nproc_ob * product(nproc_domain_orbital) =', nproc_total_wf
    print "(A,I7)", '  number of MPI process                              =', nproc_size_global
    iok = .false.
  end if

  if(nproc_total_g/=nproc_size_global)then
    print "(A)",    'product of nproc_domain_general is not correct'
    print "(A,3I5)", '  product(nproc_domain_general) =', nproc_total_g
    print "(A,3I5)", '  number of MPI process         =', nproc_size_global
    iok = .false.
  end if

  do j=1,3
    if(nproc_d_g(j)<nproc_d_o(j))then
      print "('nproc_domain_general(',I1,') is smaller than nproc_domain_orbital(',I1,')')",j,j
      print "('  nproc_domain_general(',I1,') = ',I5)", j, nproc_d_g(j)
      print "('  nproc_domain_orbital(',I1,') = ',I5)", j, nproc_d_o(j)
      iok = .false.
    end if

    if(mod(nproc_d_g(j),nproc_d_o(j))/=0)then
      print "('nproc_domain_general(',I1,') is not mutiple of nproc_domain_orbital(',I1,')')",j,j
      print "('  nproc_domain_general(',I1,') = ',I5)", j, nproc_d_g(j)
      print "('  nproc_domain_orbital(',I1,') = ',I5)", j, nproc_d_o(j)
      iok = .false.
    end if
  end do
end function check_numcpu

subroutine set_numcpu_general(iprefer_dist,numk,numo,pinfo)
  use structures, only: s_process_info
  use salmon_parallel, only: nproc_size_global
  implicit none
  integer,intent(in)               :: iprefer_dist,numk,numo
  type(s_process_info),intent(out) :: pinfo

  integer :: ip
  integer :: nproc_k,nproc_ob
  integer :: nproc_d_o(3),nproc_d_g(3)
  integer :: nproc_size_global_tmp

  integer :: ii,icount
  integer :: ir_num_factor2  ! ir means ireduced
  integer :: num_factor2
  integer :: num_factor3
  integer :: num_factor5

  integer :: nk,no

  nproc_size_global_tmp=nproc_size_global

  nk = max(numk,1)
  no = max(numo,1)

  nproc_k   = 1
  nproc_ob  = 1
  nproc_d_o = 1

  do ip=iprefer_dist,iprefer_domain_distribution

    select case(ip)

    ! k-point
    case(iprefer_k_distribution)
      if (1 < nk .and. nk < nproc_size_global_tmp) then
        if (mod(nproc_size_global_tmp,nk) == 0) then
          nproc_k               = nproc_size_global_tmp / nk
          nproc_size_global_tmp = nproc_size_global_tmp / nproc_k
        end if
      else
        nproc_k               = min(nproc_size_global_tmp,nk)
        nproc_size_global_tmp = nproc_size_global_tmp / nproc_k
      end if

    ! orbital
    case(iprefer_orbital_distribution)
      if (1 < no .and. no < nproc_size_global_tmp) then
        if (mod(nproc_size_global_tmp,no) == 0) then
          nproc_ob              = nproc_size_global_tmp / no
          nproc_size_global_tmp = nproc_size_global_tmp / nproc_ob
        end if
      else
        nproc_ob              = min(nproc_size_global_tmp,no)
        nproc_size_global_tmp = nproc_size_global_tmp / nproc_ob
      end if

    ! rgrid
    case(iprefer_domain_distribution)
      num_factor2=0
      do ii=1,26
        if(mod(nproc_size_global_tmp,2)==0)then
          num_factor2=num_factor2+1
          nproc_size_global_tmp=nproc_size_global_tmp/2
        end if
      end do

      num_factor3=0
      do ii=1,17
        if(mod(nproc_size_global_tmp,3)==0)then
          num_factor3=num_factor3+1
          nproc_size_global_tmp=nproc_size_global_tmp/3
        end if
      end do

      num_factor5=0
      do ii=1,11
        if(mod(nproc_size_global_tmp,5)==0)then
          num_factor5=num_factor5+1
          nproc_size_global_tmp=nproc_size_global_tmp/5
        end if
      end do

      if(nproc_size_global_tmp/=1)then
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

  ! general
  if (iprefer_dist == iprefer_domain_distribution) then
    nproc_d_g = nproc_d_o
  else
    nproc_size_global_tmp = nproc_size_global

    num_factor2=0
    do ii=1,26
      if(mod(nproc_size_global_tmp,2)==0)then
        num_factor2=num_factor2+1
        nproc_size_global_tmp=nproc_size_global_tmp/2
      end if
    end do

    num_factor3=0
    do ii=1,17
      if(mod(nproc_size_global_tmp,3)==0)then
        num_factor3=num_factor3+1
        nproc_size_global_tmp=nproc_size_global_tmp/3
      end if
    end do

    num_factor5=0
    do ii=1,11
      if(mod(nproc_size_global_tmp,5)==0)then
        num_factor5=num_factor5+1
        nproc_size_global_tmp=nproc_size_global_tmp/5
      end if
    end do

    if(nproc_size_global_tmp/=1)then
      stop "In automatic process distribution, prime factors for number of processes must be combination of 2, 3 or 5."
    end if

    if(num_factor2>=3)then
      nproc_d_g(1)=2
      nproc_d_g(2)=2
      nproc_d_g(3)=2
      icount=0
      ir_num_factor2=num_factor2-3
    else if(num_factor2==2)then
      nproc_d_g(1)=1
      nproc_d_g(2)=2
      nproc_d_g(3)=2
      icount=2
      ir_num_factor2=num_factor2-2
    else if(num_factor2==1)then
      nproc_d_g(1)=1
      nproc_d_g(2)=1
      nproc_d_g(3)=2
      icount=1
      ir_num_factor2=num_factor2-1
    else
      nproc_d_g(1)=1
      nproc_d_g(2)=1
      nproc_d_g(3)=1
      icount=0
      ir_num_factor2=num_factor2
    end if

    do ii=1,num_factor5
      icount=icount+1
      if(mod(icount,3)==1)then
        nproc_d_g(3)=nproc_d_g(3)*5
      else if(mod(icount,3)==2)then
        nproc_d_g(2)=nproc_d_g(2)*5
      else
        nproc_d_g(1)=nproc_d_g(1)*5
      end if
    end do

    do ii=1,num_factor3
      icount=icount+1
      if(mod(icount,3)==1)then
        nproc_d_g(3)=nproc_d_g(3)*3
      else if(mod(icount,3)==2)then
        nproc_d_g(2)=nproc_d_g(2)*3
      else
        nproc_d_g(1)=nproc_d_g(1)*3
      end if
    end do

    do ii=1,ir_num_factor2
      icount=icount+1
      if(mod(icount,3)==1)then
        nproc_d_g(3)=nproc_d_g(3)*2
      else if(mod(icount,3)==2)then
        nproc_d_g(2)=nproc_d_g(2)*2
      else
        nproc_d_g(1)=nproc_d_g(1)*2
      end if
    end do
  end if

  pinfo%npk                      = nproc_k
  pinfo%nporbital                = nproc_ob
  pinfo%npdomain_orbital(1:3)    = nproc_d_o(1:3)
  pinfo%npdomain_general(1:3)    = nproc_d_g(1:3)
  pinfo%npdomain_general_dm(1:3) = pinfo%npdomain_general(1:3)/pinfo%npdomain_orbital(1:3)

end subroutine set_numcpu_general

end module set_numcpu
