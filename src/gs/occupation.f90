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
module occupation
  implicit none

contains

SUBROUTINE ne2mu(energy,system,ilevel_print)
  use structures
  use occupation_kernel,only:solve_spectrum_occupations
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  use salmon_global, only: nelec, nelec_spin, temperature, yn_spinorbit
  implicit none
  type(s_dft_energy),intent(in) :: energy
  type(s_dft_system)            :: system
  !
  integer :: ilevel_print
  integer :: jspin,io,ik,nspin,nk,no,p5,p1,p2
  real(8) :: nein,muout,electron_count
  real(8),allocatable :: solved_occupations(:,:,:)
  logical :: occupation_ok
  character(256) :: occupation_message

  nspin = system%nspin
  nk = system%nk
  no = system%no
  
  nein = -1d0
  if(nelec/=0) then
    nein = dble(nelec)
  else if(nspin==2 .and. sum(nelec_spin(:))>0) then
    nein = dble(nelec_spin(1)+nelec_spin(2))
  end if

  call solve_spectrum_occupations(energy%esp,system%wtk,nein,temperature,yn_spinorbit=='y',&
    solved_occupations,muout,electron_count,occupation_ok,occupation_message)
  if(.not.occupation_ok)then
    if(comm_is_root(nproc_id_global))write(0,'(2a)')'Occupation solve failed: ',trim(occupation_message)
    error stop 'authoritative occupation kernel failed'
  endif
  system%rocc=solved_occupations
  system%mu = muout

  if(ilevel_print.ge.3 .and. comm_is_root(nproc_id_global)) then
     write(*,*)
     write(*,*) 'Fractional Occupation Numbers :'
     write(*,*)
     do ik=1,nk
     do jspin=1,nspin
        if(ik<=10)then
          print *, ' iik = ',ik
          do p5=1,(no+4)/5
             p1=5*(p5-1)+1
             p2=5*p5 ; if ( p2 > no ) p2=no
             write(*,'(1x,5(i6,f8.4,1x))')  (io,system%rocc(io,ik,jspin),io=p1,p2)
          end do
        endif
    end do
    end do
    write(*,*)
    write(*,'(a,f15.8)') ' Fermi level (a.u.)  = ',muout
    write(*,'(a,f15.8)') ' Number of Electrons = ',electron_count
    write(*,*)
  end if

  return
END SUBROUTINE ne2mu

end module Occupation
