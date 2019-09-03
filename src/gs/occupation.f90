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
module occupation
  implicit none

contains

SUBROUTINE ne2mu(energy,system,info)
  use structures
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use salmon_global, only: nelec, temperature
  implicit none
  type(s_dft_energy),intent(in) :: energy
  type(s_dft_system)            :: system
  type(s_orbital_parallel)      :: info
  !
  integer :: jspin,io,ik,nspin,nk,no,iter
  real(8) :: nein,muout
  real(8) :: mu1,mu2,mu3,ne1,ne3,ne3o,diff_ne
  integer :: ii,p5,p1,p2
  integer :: nc
  real(8) :: diff_mu,muo
  real(8) :: diff_ne2

  nspin = system%nspin
  nk = system%nk
  no = system%no

  nein = dble(nelec)

  mu1 = minval(energy%esp)
  mu2 = maxval(energy%esp)
  nc=0

  ITERATION: do iter=1,100000
    diff_ne = 100.d0
    diff_ne2 = 100.d0

    ne1  = 0d0
    ne3  = 0d0
    ne3o = 0d0

    diff_mu =100.d0
    muo=0.0d0

    call mu2ne(mu1,ne1)

    do ii=1,1000
      if ( ii .eq. 1000 ) then
        if ( nc .le. 50)  then
           nc= nc + 1
           mu1 = minval(energy%esp) - 0.2d0*dble(nc)
           mu2 = maxval(energy%esp) + 0.2d0*dble(nc)
           cycle ITERATION
        else
           print *,'=================================='
           print *,'Const Ne does not converged!!!!!!!'
           print *,'=================================='
           exit ITERATION
        endif
      endif
      if ( diff_ne < 1d-10  .and. diff_mu < 1d-9  &
         .and. diff_ne2 < 1.d-9 ) exit ITERATION
      mu3 = mu1 + (mu2-mu1)/2.d0
      ne3=0
      call mu2ne(mu3,ne3)
      diff_ne = abs((ne3-ne3o)/ne3)
      diff_ne2 = abs(nein-ne3)
      if ( (ne1-nein)*(ne3-nein) > 0 ) then
        mu1=mu3
        ne1=ne3
      else
        mu2=mu3
      end if

      ne3o = ne3
      diff_mu = mu3 - muo
      muo  = mu3
    end do
  end do ITERATION

  muout = mu3
  call mu2ne(muout,ne3)

  if(comm_is_root(nproc_id_global))then
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
             write(*,'(1x,5(i3,f8.4,1x))')  (io,system%rocc(io,ik,jspin),io=p1,p2)
          end do
        endif
    end do
    end do
    write(*,*)
    write(*,'(a,f15.8)') ' Fermi level (a.u.)  = ',muout
    write(*,'(a,f15.8)') ' Number of Electrons = ',ne3
    write(*,*)
  end if

  do jspin=1,system%nspin
    do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        info%occ(io,ik,jspin,1) = system%rocc(io,ik,jspin)*system%wtk(ik)
      end do
    end do
  end do

  return
contains

  SUBROUTINE mu2ne(muin,neout)
    implicit none
    real(8), intent(in)  :: muin
    real(8), intent(out) :: neout
    !
    real(8) :: fact
    neout=0d0

    do jspin=1,nspin
    do ik=1,nk
    do io=1,no
       fact = (energy%esp(io,ik,jspin)-muin)/temperature
       if(fact.ge.40.d0) then
          system%rocc(io,ik,jspin) = 0.d0
       else
          system%rocc(io,ik,jspin) = 2d0/( 1.d0 + exp( (energy%esp(io,ik,jspin)-muin)/temperature ) )
       endif
       neout = neout + system%rocc(io,ik,jspin)*system%wtk(ik)
    end do
    end do
    end do

  END SUBROUTINE mu2ne

END SUBROUTINE ne2mu

end module Occupation
