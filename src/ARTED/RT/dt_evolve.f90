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
!This file is "dt_evolve.f90"
!This file contain one subroutine.
!Subroutine dt_evolve
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine dt_evolve_KB(iter)
  use global_variables, only: propagator,kAc,kAc0,Ac_tot,zu_t
  implicit none
  integer, intent(in) :: iter

  select case(propagator)
    case('middlepoint')
      call default_propagator
    case('etrs')
      call etrs_propagator
    case default
      call err_finalize('invalid propagator')
  end select

contains
  subroutine default_propagator
    implicit none
    integer :: ixyz
    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+0.5d0*(Ac_tot(iter,ixyz) + Ac_tot(iter+1,ixyz) )
    enddo
    call dt_evolve_omp_KB(zu_t)
  end subroutine

  subroutine etrs_propagator
    use global_variables, only: kAc_new
    implicit none
    integer :: ixyz
    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter,ixyz)
      kAc_new(:,ixyz)=kAc0(:,ixyz)+Ac_tot(iter+1,ixyz)
    enddo
    call dt_evolve_etrs_omp_KB(zu_t)
  end subroutine
end subroutine

subroutine dt_evolve_KB_MS(imacro)
  use global_variables, only: propagator,kAc,kAc0,Ac_new_m,ac_m,zu_m
  implicit none
  integer, intent(in) :: imacro

  select case(propagator)
    case('middlepoint')
      call default_propagator
    case('etrs')
      call etrs_propagator
    case default
      call err_finalize('invalid propagator')
  end select

contains
  subroutine default_propagator
    implicit none
    integer :: ixyz
    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+(Ac_new_m(ixyz,imacro)+Ac_m(ixyz,imacro))/2d0
    enddo
    call dt_evolve_omp_KB_MS(zu_m(:,:,:,imacro))
  end subroutine

  subroutine etrs_propagator
    use global_variables, only: kAc_new
    implicit none
    integer :: ixyz
    do ixyz=1,3
      kAc(:,ixyz)=kAc0(:,ixyz)+Ac_m(ixyz,imacro)
      kAc_new(:,ixyz)=kAc0(:,ixyz)+Ac_new_m(ixyz,imacro)
    enddo
    call dt_evolve_etrs_omp_KB(zu_m(:,:,:,imacro))
  end subroutine
end subroutine

! ---------------------------------------------

Subroutine dt_evolve_omp_KB(zu)
  use Global_Variables
  use projector
  use timer
  use opt_variables
  implicit none
  complex(8),intent(inout) :: zu(NL,NBoccmax,NK_s:NK_e)
  integer    :: ik,ib,ikb
  integer    :: i

  call timer_begin(LOG_CALC_TIME_PROPAGATION)

!Constructing nonlocal part
  call update_projector(kac)

! yabana
  select case(functional)
  case('VS98','TPSS','TBmBJ','BJ_PW','tbmbj','bj_pw')

!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu_GS(:,ib,ik)=zu(:,ib,ik)
  end do
    
  Vloc_t=Vloc

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass_t=tmass
    tjr_t=tjr
    tjr2_t=tjr2
  end if

  call hamiltonian_arted(zu,.false.)

  call psi_rho_RT(zu)
  if(yn_fix_func=='n')then
     call Hartree
     call Exc_Cor(calc_mode_rt,NBoccmax,zu)
  endif

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  Vloc=0.5d0*(Vloc+Vloc_t)

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass=0.5d0*(tmass+tmass_t)
    tjr=0.5d0*(tjr+tjr_t)
    tjr2=0.5d0*(tjr2+tjr2_t)
  end if


!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu(:,ib,ik)=zu_GS(:,ib,ik)
  end do

  end select
! yabana

  call hamiltonian_arted(zu,.true.)

  call psi_rho_RT(zu)

  if(yn_fix_func=='n')then
  call Hartree

! yabana
  call Exc_Cor(calc_mode_rt,NBoccmax,zu)
! yabana
  endif

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  call timer_end(LOG_CALC_TIME_PROPAGATION)

  return
End Subroutine dt_evolve_omp_KB
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine dt_evolve_etrs_omp_KB(zu)
  use Global_Variables
  use projector
  use timer
  use opt_variables
  implicit none
  complex(8),intent(inout) :: zu(NL,NBoccmax,NK_s:NK_e)
  integer    :: ik,ib,ikb
  integer    :: i
  real(8)    :: dt_t

  call timer_begin(LOG_CALC_TIME_PROPAGATION)

  dt_t = dt; dt = 0.5d0*dt

!Constructing nonlocal part
  call update_projector(kac)

  call hamiltonian_arted(zu,.false.)

  Vloc_t=Vloc
  Vloc_new(:) = 3d0*Vloc(:) - 3d0*Vloc_old(:,1) + Vloc_old(:,2)
  Vloc_old(:,2) = Vloc_old(:,1)
  Vloc_old(:,1) = Vloc(:)
  Vloc(:) = Vloc_new(:)

  kAc=kAc_new

!Constructing nonlocal part
  call update_projector(kac)

!== predictor-corrector ==
  select case(functional)
  case('VS98','TPSS','TBmBJ','BJ_PW','tbmbj','bj_pw')

!$omp parallel do private(ik,ib)
     do ikb=1,NKB
        ik=ik_table(ikb) ; ib=ib_table(ikb)
        zu_GS(:,ib,ik)=zu(:,ib,ik)
     end do

     call hamiltonian_arted(zu,.false.)

     call psi_rho_RT(zu)

     if(yn_fix_func=='n')then
     call Hartree

     call Exc_Cor(calc_mode_rt,NBoccmax,zu)
     endif

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

!$omp parallel do private(ik,ib)
     do ikb=1,NKB
        ik=ik_table(ikb) ; ib=ib_table(ikb)
        zu(:,ib,ik)=zu_GS(:,ib,ik)
     end do

  end select


  call hamiltonian_arted(zu,.true.)

  call psi_rho_RT(zu)

  if(yn_fix_func=='n')then
  call Hartree

  call Exc_Cor(calc_mode_rt,NBoccmax,zu)
  endif

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  dt = dt_t
  call timer_end(LOG_CALC_TIME_PROPAGATION)

  return
End Subroutine dt_evolve_etrs_omp_KB
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine dt_evolve_omp_KB_MS(zu)
  use Global_Variables
  use projector
  use timer
  use opt_variables
  implicit none
  complex(8),intent(inout) :: zu(NL,NBoccmax,NK_s:NK_e)
  integer    :: ik,ib,ikb
  integer    :: i

  call timer_begin(LOG_CALC_TIME_PROPAGATION)

!Constructing nonlocal part ! sato
  call update_projector(kac)


! yabana
  select case(functional)
  case('VS98','TPSS','TBmBJ','BJ_PW','tbmbj','bj_pw')


!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu_GS(:,ib,ik)=zu(:,ib,ik)
  end do

  Vloc_t=Vloc

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass_t=tmass
    tjr_t=tjr
    tjr2_t=tjr2
  end if

  call hamiltonian_arted(zu,.false.)

  call psi_rho_RT(zu)
  if(yn_fix_func=='n')then
     call Hartree
     call Exc_Cor(calc_mode_rt,NBoccmax,zu)
  endif

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  Vloc=0.5d0*(Vloc+Vloc_t)

  if(functional == 'VS98' .or. functional == 'TPSS')then
    tmass=0.5d0*(tmass+tmass_t)
    tjr=0.5d0*(tjr+tjr_t)
    tjr2=0.5d0*(tjr2+tjr2_t)
  end if

!$omp parallel do private(ik,ib)
  do ikb=1,NKB
    ik=ik_table(ikb) ; ib=ib_table(ikb)
    zu(:,ib,ik)=zu_GS(:,ib,ik)
  end do

  end select
! yabana

  call hamiltonian_arted(zu,.true.)

  call psi_rho_RT(zu)

  if(yn_fix_func=='n')then
  call Hartree

! yabana
  call Exc_Cor(calc_mode_rt,NBoccmax,zu)
! yabana
  endif

!$omp parallel do
  do i=1,NL
    Vloc(i)=Vh(i)+Vpsl(i)+Vexc(i)
  end do

  call timer_end(LOG_CALC_TIME_PROPAGATION)

  return
End Subroutine dt_evolve_omp_KB_MS
