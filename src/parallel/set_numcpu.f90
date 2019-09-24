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

contains

subroutine check_numcpu(nproc_d_o,nproc_d_g,nproc_d_g_dm)
  use salmon_global, only: nproc_k,nproc_ob
  use salmon_parallel, only: nproc_size_global
  implicit none
  integer,intent(in) :: nproc_d_o(3)
  integer,intent(in) :: nproc_d_g(3)
  integer,intent(out) :: nproc_d_g_dm(3)
  integer :: j

  if(nproc_k*nproc_ob*nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)/=nproc_size_global)then
    write(*,*) "inumcpu_check error!"
    write(*,*) "number of cpu is not correct!"
    stop
  end if
  do j=1,3
    if(nproc_d_g(j)<nproc_d_o(j))then
      write(*,*) "inumcpu_check error!"
      write(*,*) "nproc_domain_general is smaller than nproc_domain_orbital."
      stop
    end if
  end do
  if(nproc_d_g(1)*nproc_d_g(2)*nproc_d_g(3)>nproc_size_global)then
    write(*,*) "inumcpu_check error!"
    write(*,*) "product of nproc_domain_general is larger than nproc."
    stop
  end if
  if(mod(nproc_d_g(1),nproc_d_o(1))/=0)then
    write(*,*) "inumcpu_check error!"
    write(*,*) "nproc_domain_general(1) is not mutiple of nproc_domain_orbital(1)."
    stop
  end if
  if(mod(nproc_d_g(2),nproc_d_o(2))/=0)then
    write(*,*) "inumcpu_check error!"
    write(*,*) "nproc_domain_general(2) is not mutiple of nproc_domain_orbital(2)."
    stop
  end if
  if(mod(nproc_d_g(3),nproc_d_o(3))/=0)then
    write(*,*) "inumcpu_check error!"
    write(*,*) "nproc_domain_general(3) is not mutiple of nproc_domain_orbital(3)."
    stop
  end if
  nproc_d_g_dm(1:3)=nproc_d_g(1:3)/nproc_d_o(1:3)

end subroutine check_numcpu

subroutine set_numcpu_gs(nproc_d_o,nproc_d_g,nproc_d_g_dm)
  use salmon_global, only: nproc_k,nproc_ob,nproc_domain_orbital,nproc_domain_general
  use salmon_parallel, only: nproc_size_global
  implicit none
  integer,intent(out) :: nproc_d_o(3)
  integer,intent(out) :: nproc_d_g(3)
  integer,intent(out) :: nproc_d_g_dm(3)
  integer :: ii
  integer :: nproc_size_global_tmp
  integer :: nproc_d_o_tmp(3)
  
  integer :: num_factor2
  integer :: num_factor3
  integer :: num_factor5

  integer :: icount
 
  nproc_size_global_tmp=nproc_size_global
  
  ! this code treats the situation that nproc_size_global is less than or equal to 48,828,125
  
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
    stop "In automatic process assignment, prime factors for number of processes must be combination of 2, 3 or 5."
  end if

  nproc_d_o_tmp(1:3)=1
 
  icount=0

  do ii=1,num_factor5
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_d_o_tmp(3)=nproc_d_o_tmp(3)*5
    else if(mod(icount,3)==2)then
      nproc_d_o_tmp(2)=nproc_d_o_tmp(2)*5
    else
      nproc_d_o_tmp(1)=nproc_d_o_tmp(1)*5
    end if
  end do

  do ii=1,num_factor3
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_d_o_tmp(3)=nproc_d_o_tmp(3)*3
    else if(mod(icount,3)==2)then
      nproc_d_o_tmp(2)=nproc_d_o_tmp(2)*3
    else
      nproc_d_o_tmp(1)=nproc_d_o_tmp(1)*3
    end if
  end do

  do ii=1,num_factor2
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_d_o_tmp(3)=nproc_d_o_tmp(3)*2
    else if(mod(icount,3)==2)then
      nproc_d_o_tmp(2)=nproc_d_o_tmp(2)*2
    else
      nproc_d_o_tmp(1)=nproc_d_o_tmp(1)*2
    end if
  end do

  nproc_k=1
  nproc_ob=1
  nproc_d_o(1:3)=nproc_d_o_tmp(1:3)
  nproc_d_g(1:3)=nproc_d_o_tmp(1:3)
  nproc_d_g_dm(1:3)=nproc_d_g(1:3)/nproc_d_o(1:3)

  nproc_domain_orbital = nproc_d_o
  nproc_domain_general = nproc_d_g
  
end subroutine set_numcpu_gs

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

subroutine set_numcpu_rt(nproc_d_o,nproc_d_g,nproc_d_g_dm)
  use salmon_global, only: nproc_k,nproc_ob,nproc_domain_orbital,nproc_domain_general
  use salmon_parallel, only: nproc_size_global
  implicit none
  integer,intent(out) :: nproc_d_o(3)
  integer,intent(out) :: nproc_d_g(3)
  integer,intent(out) :: nproc_d_g_dm(3)
  integer :: ii
  integer :: nproc_size_global_tmp
  integer :: nproc_ob_tmp
  integer :: nproc_d_o_tmp(3)
  integer :: nproc_d_g_tmp(3)
  
  integer :: num_factor2
  integer :: ir_num_factor2  ! ir means ireduced
  integer :: num_factor3
  integer :: num_factor5

  integer :: icount
 
  nproc_size_global_tmp=nproc_size_global
  
  ! this code treats the situation that nproc_size_global is less than or equal to 48,828,125

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
    stop "In automatic process assignment, prime factors for number of processes must be combination of 2, 3 or 5."
  end if


  if(num_factor2>=3)then
    nproc_ob_tmp=nproc_size_global/8
    nproc_d_o_tmp(1)=2 
    nproc_d_o_tmp(2)=2 
    nproc_d_o_tmp(3)=2 
    icount=0
    ir_num_factor2=num_factor2-3
  else if(num_factor2==2)then
    nproc_ob_tmp=nproc_size_global/4
    nproc_d_o_tmp(1)=1 
    nproc_d_o_tmp(2)=2 
    nproc_d_o_tmp(3)=2 
    icount=2
    ir_num_factor2=num_factor2-2
  else if(num_factor2==1)then
    nproc_ob_tmp=nproc_size_global/2
    nproc_d_o_tmp(1)=1 
    nproc_d_o_tmp(2)=1 
    nproc_d_o_tmp(3)=2 
    icount=1
    ir_num_factor2=num_factor2-1
  else
    nproc_ob_tmp=nproc_size_global
    nproc_d_o_tmp(1)=1 
    nproc_d_o_tmp(2)=1 
    nproc_d_o_tmp(3)=1 
    icount=0
    ir_num_factor2=num_factor2
  end if

  nproc_d_g_tmp(1:3)=nproc_d_o_tmp(1:3)

  do ii=1,num_factor5
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_d_g_tmp(3)=nproc_d_g_tmp(3)*5
    else if(mod(icount,3)==2)then
      nproc_d_g_tmp(2)=nproc_d_g_tmp(2)*5
    else
      nproc_d_g_tmp(1)=nproc_d_g_tmp(1)*5
    end if
  end do

  do ii=1,num_factor3
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_d_g_tmp(3)=nproc_d_g_tmp(3)*3
    else if(mod(icount,3)==2)then
      nproc_d_g_tmp(2)=nproc_d_g_tmp(2)*3
    else
      nproc_d_g_tmp(1)=nproc_d_g_tmp(1)*3
    end if
  end do

  do ii=1,ir_num_factor2
    icount=icount+1
    if(mod(icount,3)==1)then
      nproc_d_g_tmp(3)=nproc_d_g_tmp(3)*2
    else if(mod(icount,3)==2)then
      nproc_d_g_tmp(2)=nproc_d_g_tmp(2)*2
    else
      nproc_d_g_tmp(1)=nproc_d_g_tmp(1)*2
    end if
  end do

  nproc_k=1
  nproc_ob=nproc_ob_tmp
  nproc_d_o(1:3)=nproc_d_o_tmp(1:3)
  nproc_d_g(1:3)=nproc_d_g_tmp(1:3)
  nproc_d_g_dm(1:3)=nproc_d_g(1:3)/nproc_d_o(1:3)

  nproc_domain_orbital = nproc_d_o
  nproc_domain_general = nproc_d_g

end subroutine set_numcpu_rt



end module set_numcpu
