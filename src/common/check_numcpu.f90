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
module check_numcpu_sub
  implicit none

contains

subroutine check_numcpu(nproc_d_o,nproc_d_g,nproc_d_g_dm)
  use inputoutput, only: nproc_k,nproc_ob
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

end module check_numcpu_sub
