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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

! for single-scale Maxwell-TDDFT method
module fdtd_coulomb_gauge
  implicit none

contains

subroutine fdtd_singlescale(itt,lg,mg,system,info,rho,Vh,j_e,srg_scalar,Ac,div_Ac,fw)
  use structures
  use math_constants,only : zi,pi
  use phys_constants, only: cspeed_au
  use salmon_global, only: dt,yn_gbp
  use sendrecv_grid, only: update_overlap_real8
  use stencil_sub, only: calc_gradient_field
  use communication, only: comm_is_root, comm_summation
  use inputoutput, only: t_unit_time
  use timer
  implicit none
  integer                 ,intent(in) :: itt
  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info)   ,intent(in) :: info
  type(s_scalar)          ,intent(in) :: rho,Vh ! electron number density & Hartree potential
  type(s_vector)          ,intent(in) :: j_e    ! electron number current density (without rho*A/c)
  type(s_sendrecv_grid)               :: srg_scalar
  type(s_vector)                      :: Ac     ! A/c, A: vector potential, c: speed of light
  type(s_scalar)                      :: div_Ac ! div(A/c)
  type(s_singlescale)                 :: fw     ! FDTD working arrays, etc.
  !
  integer,parameter :: mstep=100
  integer :: ix,iy,iz,i1,ii,krd(3,3),lcs(3,3,3),dr(3)
  real(8) :: Hvol,hgs(3),dt_m,tm,coef,lap_A,Energy_em,diff_A,coef2 &
  & ,e_em,e_em_wrk,e_joule,e_joule_wrk,e_poynting(2),e_poynting_wrk(2),rho_t
  real(8),dimension(3) :: out_curr,out_Aext,out_Ab1,out_Ab2,wrk,wrk2,wrk3,wrk4,vec_je,Aext0,Aext1,Aext0_old,Aext1_old
  real(8) :: e_poy1,e_poy2,rtmp1(6),rtmp2(6)

  call timer_begin(LOG_SS_FDTD_CALC)

  krd = 0
  krd(1,1) = 1; krd(2,2) = 1; krd(3,3) = 1

  lcs = 0
  lcs(1,2,3) = 1; lcs(3,1,2) = 1; lcs(2,3,1) = 1
  lcs(1,3,2) = -1; lcs(2,1,3) = -1; lcs(3,2,1) = -1

  hgs = system%hgs
  Hvol = system%Hvol
  dt_m = dt / dble(mstep)

!-----------------------------------------------------------------------------------------------------------------------------------

  fw%box = 0d0
  wrk = 0d0
!$OMP parallel do collapse(2) private(ix,iy,iz,vec_je,rho_t) reduction(+:wrk)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    vec_je = ( j_e%v(1:3,ix,iy,iz) + fw%vec_je_old(1:3,ix,iy,iz) )*0.5d0 ! j_e(t) = ( j_e(t+dt/2) + j_e(t-dt/2) )/2
    rho_t  = ( rho%f(ix,iy,iz)     + fw%rho_old(ix,iy,iz)        )*0.5d0 ! rho(t) = ( rho(t+dt/2) + rho(t-dt/2) )/2
    fw%curr(ix,iy,iz,1:3) = vec_je + rho_t * fw%vec_Ac_m(1,ix,iy,iz,1:3) ! curr(t): electron number current density
    wrk = wrk + fw%curr(ix,iy,iz,1:3)

    fw%box(ix,iy,iz) = Vh%f(ix,iy,iz) ! Vh(t+dt/2)
  end do
  end do
  end do
  wrk = wrk/dble(lg%num(1)*lg%num(2)*lg%num(3))
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  call comm_summation(wrk,out_curr,3,info%icomm_r)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

! gradient of d(Vh)/dt (Vh: Hartree potential)
  call timer_begin(LOG_SS_FDTD_COMM)
  if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, fw%box1)
  call timer_end(LOG_SS_FDTD_COMM)

  call timer_begin(LOG_SS_FDTD_CALC)
  call calc_gradient_field(mg,fw%coef_nab,fw%box1,fw%grad_Vh) ! grad[Vh(t+dt/2)]

  !$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    fw%current4pi(ix,iy,iz,1:3) = ( fw%grad_Vh(1:3,ix,iy,iz) - fw%grad_Vh_old(1:3,ix,iy,iz) ) /dt & ! d(grad(Vh))/dt
                                & - 4d0*pi * fw%curr(ix,iy,iz,1:3) ! 4*pi* (electric current density)
  end do
  end do
  end do
  call timer_end(LOG_SS_FDTD_CALC)

!-----------------------------------------------------------------------------------------------------------------------------------

! FDTD loop: Ac(t) --> Ac(t+dt)

  if(yn_gbp=='y') then
    call fdtd_gbp
  else
    call fdtd
  end if

!-----------------------------------------------------------------------------------------------------------------------------------

! divergence & rotation of Ac(t+dt)

  call timer_begin(LOG_SS_FDTD_CALC)
  fw%div_Ac = 0d0
  fw%rot_Ac = 0d0
  do i1=1,3
  
!$OMP parallel do collapse(2) private(ix,iy,iz)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      fw%box(ix,iy,iz) = fw%vec_Ac_m(1,ix,iy,iz,i1)
    end do
    end do
    end do
    call timer_begin(LOG_SS_FDTD_COMM)
    call update_overlap_real8(fw%srg_eg, fw%eg, fw%box)
    call timer_end(LOG_SS_FDTD_COMM)
    if(mg%is(3)==lg%is(3))then
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        fw%box(ix,iy,lg%is(3)-1) = fw%vec_Ac_boundary_bottom(ix,iy,i1)
      end do
      end do
    end if
    if(mg%ie(3)==lg%ie(3))then
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        fw%box(ix,iy,lg%ie(3)+1) = fw%vec_Ac_boundary_top(ix,iy,i1)
      end do
      end do
    end if

  !$OMP parallel do collapse(2) private(ix,iy,iz,wrk,dr,diff_A)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
    ! rot(A)
      wrk(1) = ( fw%box(ix+1,iy,iz) - fw%box(ix-1,iy,iz) ) / ( 2d0* Hgs(1) )
      wrk(2) = ( fw%box(ix,iy+1,iz) - fw%box(ix,iy-1,iz) ) / ( 2d0* Hgs(2) )
      wrk(3) = ( fw%box(ix,iy,iz+1) - fw%box(ix,iy,iz-1) ) / ( 2d0* Hgs(3) )
      fw%rot_Ac(:,ix,iy,iz) = fw%rot_Ac(:,ix,iy,iz) + lcs(:,1,i1) * wrk(1) + lcs(:,2,i1) * wrk(2) + lcs(:,3,i1) * wrk(3)
    ! div(A)
      dr = krd(:,i1)
      diff_A = ( fw%box(ix+dr(1),iy+dr(2),iz+dr(3)) - fw%box(ix-dr(1),iy-dr(2),iz-dr(3)) ) / ( 2d0* Hgs(i1) )
      fw%div_Ac(ix,iy,iz) = fw%div_Ac(ix,iy,iz) + diff_A
    end do
    end do
    end do
    
  end do

!-----------------------------------------------------------------------------------------------------------------------------------

! Ac & div_Ac for TDDFT

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)

    Ac%v(1:3,ix,iy,iz) = ( fw%vec_Ac_m(1,ix,iy,iz,1:3) + fw%vec_Ac_old(1:3,ix,iy,iz) ) * 0.5d0 ! Ac(t+dt/2) = ( A(t+dt) + A(t) )/2
    div_Ac%f(ix,iy,iz)   = ( fw%div_Ac(ix,iy,iz) + fw%div_Ac_old(ix,iy,iz) ) * 0.5d0 ! div ( A(t+dt) + A(t) )/2

  end do
  end do
  end do

!-----------------------------------------------------------------------------------------------------------------------------------

! Poynting vector, energy, etc.

  coef = cspeed_au / (4d0*pi)
  coef2 = Hvol / (8d0*pi)
  e_em_wrk = 0d0 ! for Electro-Magnetic energy
  e_joule_wrk = 0d0 ! for Joule dissipated power

!$OMP parallel do collapse(2) private(ix,iy,iz,wrk,wrk2,wrk3,wrk4) reduction(+:e_em_wrk,e_joule_wrk)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)

    wrk4 = ( fw%vec_Ac_m(1,ix,iy,iz,:) - fw%vec_Ac_old(:,ix,iy,iz) ) / dt ! (A(t+dt)-A(t))/dt
    wrk  = - (-1d0)*fw%grad_Vh(:,ix,iy,iz) - wrk4 ! E
    wrk2 = cspeed_au * fw%rot_Ac(:,ix,iy,iz)      ! B
    wrk3 = - fw%curr(ix,iy,iz,1:3)                ! j
    fw%poynting_vector(:,ix,iy,iz) = coef * ( lcs(:,1,2) * wrk(1) * wrk2(2) + lcs(:,1,3) * wrk(1) * wrk2(3) &
                                            + lcs(:,2,1) * wrk(2) * wrk2(1) + lcs(:,2,3) * wrk(2) * wrk2(3) &
                                            + lcs(:,3,1) * wrk(3) * wrk2(1) + lcs(:,3,2) * wrk(3) * wrk2(2) ) ! E x B
    e_em_wrk = e_em_wrk + coef2 * ( sum(wrk**2) + sum(wrk2**2) ) ! ( E^2 + B^2 )/(8*pi)
    e_joule_wrk = e_joule_wrk + sum(wrk3*wrk) * Hvol ! j*E

  end do
  end do
  end do
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  call comm_summation(e_em_wrk,e_em,info%icomm_r)
  call comm_summation(e_joule_wrk,e_joule,info%icomm_r)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

  call timer_begin(LOG_SS_FDTD_CALC)

!-----------------------------------------------------------------------------------------------------------------------------------

  ! integral(A) @ z = 0 (bottom boundary)
  wrk = 0d0
  if (mg%is(3) == lg%is(3)) then
!$omp parallel do private(ix,iy,iz) reduction(+:wrk)
    do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        wrk(1) = wrk(1) + Ac%v(1,ix,iy,mg%is(3))
        wrk(2) = wrk(2) + Ac%v(2,ix,iy,mg%is(3))
        wrk(3) = wrk(3) + Ac%v(3,ix,iy,mg%is(3))
      end do
    end do
  end if

  ! integral(A) @ z = az (top boundary)
  wrk3 = 0d0
  if (mg%ie(3) == lg%ie(3)) then
!$omp parallel do private(ix,iy,iz) reduction(+:wrk3)
    do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        wrk3(1) = wrk3(1) + Ac%v(1,ix,iy,mg%ie(3))
        wrk3(2) = wrk3(2) + Ac%v(2,ix,iy,mg%ie(3))
        wrk3(3) = wrk3(3) + Ac%v(3,ix,iy,mg%ie(3))
      end do
    end do
  end if
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  rtmp1 = [wrk, wrk3]
  call comm_summation(rtmp1,rtmp2,6,info%icomm_r)
  wrk (1:3) = rtmp2(1:3)
  wrk3(1:3) = rtmp2(4:6)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

  call timer_begin(LOG_SS_FDTD_CALC)
  out_Ab1 = wrk  / dble(lg%num(1)*lg%num(2))
  out_Ab2 = wrk3 / dble(lg%num(1)*lg%num(2))

! Surface integral of the poynting vector S
  e_poy1 = 0d0
  e_poy2 = 0d0
  coef = Hgs(1)*Hgs(2)
  if(mg%is(3)==lg%is(3)) then ! integral(S) @ z = 0 (bottom boundary)
!$omp parallel do private(ix,iy,iz) reduction(+:e_poy1)
    do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        e_poy1 = e_poy1 + fw%poynting_vector(3,ix,iy,lg%is(3)) * coef
      end do
    end do
  end if
  if(mg%ie(3)==lg%ie(3)) then ! integral(S) @ z = az (top boundary)
!$omp parallel do private(ix,iy,iz) reduction(+:e_poy2)
    do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        e_poy2 = e_poy2 + fw%poynting_vector(3,ix,iy,lg%ie(3)) * coef
      end do
    end do
  end if
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  e_poynting_wrk = [e_poy1, e_poy2]
  call comm_summation(e_poynting_wrk,e_poynting,2,info%icomm_r)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

  call timer_begin(LOG_SS_FDTD_CALC)
  Energy_em = e_em
  fw%Energy_joule = fw%Energy_joule + dt*e_joule
  fw%Energy_poynting = fw%Energy_poynting + dt*e_poynting

  if(comm_is_root(info%id_rko)) write(fw%fh_rt_micro,'(99(1X,E23.15E3))') &
    dble(itt)*dt*t_unit_time%conv,out_Ab1,out_Ab2,out_Aext,out_curr,fw%E_electron,fw%Energy_poynting,Energy_em,fw%Energy_joule
    
!-----------------------------------------------------------------------------------------------------------------------------------

! FIXME: is following result unused?
! for spatial distribution of excitation energy
!  coef = Hgs(1)*Hgs(2)
!  fw%integral_poynting_tmp = 0d0
!!$omp parallel do collapse(2) private(iz,iy,ix)
!  do iy=mg%is(2),mg%ie(2)
!  do ix=mg%is(1),mg%ie(1)
!    do iz=mg%is(3),mg%ie(3)
!      fw%integral_poynting_tmp(iz) = fw%integral_poynting_tmp(iz) + fw%poynting_vector(3,ix,iy,iz) * coef
!    end do
!  end do
!  end do
!  call timer_end(LOG_SS_FDTD_CALC)
!
!  call timer_begin(LOG_SS_FDTD_COMM_COLL)
!  call comm_summation(fw%integral_poynting_tmp,fw%integral_poynting_tmp2,lg_num(3),info%icomm_r)
!  call timer_end(LOG_SS_FDTD_COMM_COLL)
!
!  call timer_begin(LOG_SS_FDTD_CALC)
!  fw%integral_poynting = fw%integral_poynting + dt * fw%integral_poynting_tmp2

!  fw%Ac_zt_t = 0d0
!! for the vector potential Ax(z,t)
!  do iz=mg%is(3),mg%ie(3)
!!$omp parallel do collapse(2) private(iy,ix)
!    do iy=mg%is(2),mg%ie(2)
!    do ix=mg%is(1),mg%ie(1)
!        fw%Ac_zt_t(1,iz) = fw%Ac_zt_t(1,iz) + Ac%v(1,ix,iy,iz) / (lg%num(1)*lg%num(2))
!        fw%Ac_zt_t(2,iz) = fw%Ac_zt_t(2,iz) + Ac%v(2,ix,iy,iz) / (lg%num(1)*lg%num(2))
!        fw%Ac_zt_t(3,iz) = fw%Ac_zt_t(3,iz) + Ac%v(3,ix,iy,iz) / (lg%num(1)*lg%num(2))
!    end do
!    end do
!  end do
!  call timer_end(LOG_SS_FDTD_CALC)
!
!  call timer_begin(LOG_SS_FDTD_COMM_COLL)
!  call comm_summation(fw%Ac_zt_t,fw%Ac_zt,size(fw%Ac_zt),info%icomm_r)
!  call timer_end(LOG_SS_FDTD_COMM_COLL)
!
!  if(comm_is_root(info%id_rko)) then
!    do iz=lg%is(3),lg%ie(3)
!      write(fw%fh_Ac_zt,fmt='(99(1X,E23.15E3))',advance='no') dble(iz)*hgs(3),fw%Ac_zt(1,iz),fw%Ac_zt(2,iz),fw%Ac_zt(3,iz)
!    end do
!    write(fw%fh_Ac_zt,'()')
!  end if

!-----------------------------------------------------------------------------------------------------------------------------------

! stock old Ac, grad_Vh, j_e, & rho

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    fw%vec_Ac_old(:,ix,iy,iz)    = fw%vec_Ac_m(1,ix,iy,iz,1:3) ! Ac(t+dt) --> Ac(t) of next step
    fw%div_Ac_old(ix,iy,iz)      = fw%div_Ac(ix,iy,iz)      ! div Ac(t+dt) --> div Ac(t) of next step
    fw%grad_Vh_old(1:3,ix,iy,iz) = fw%grad_Vh(1:3,ix,iy,iz) ! grad[Vh(t-dt/2)] of next step
    fw%vec_je_old(1:3,ix,iy,iz)  = j_e%v(1:3,ix,iy,iz) ! j_e(t-dt/2) of next step
    fw%rho_old(ix,iy,iz)         = rho%f(ix,iy,iz)     ! rho(t-dt/2) of next step
  end do
  end do
  end do
  
  call timer_end(LOG_SS_FDTD_CALC)

!-----------------------------------------------------------------------------------------------------------------------------------

  return

contains

!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine pulse(t,r,A_ext)
    use em_field, only: calc_Ac_ext
    implicit none
    real(8),intent(in)  :: t,r
    real(8),intent(out) :: A_ext(3)
    !
    real(8) :: tt

    tt = t - r/cspeed_au
    call calc_Ac_ext(tt,A_ext)

    return
  end subroutine pulse
  
!-----------------------------------------------------------------------------------------------------------------------------------
  
  subroutine fdtd
    implicit none

    call timer_begin(LOG_SS_FDTD_CALC)
    do ii=1,mstep

    !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)

        fw%vec_Ac_m(-1,ix,iy,iz,1:3) = fw%vec_Ac_m(0,ix,iy,iz,1:3)
        fw%vec_Ac_m(0 ,ix,iy,iz,1:3) = fw%vec_Ac_m(1,ix,iy,iz,1:3)

      end do
      end do
      end do

      do i1=1,3

    !$OMP parallel do collapse(2) private(ix,iy,iz)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          fw%box(ix,iy,iz) = fw%vec_Ac_m(0,ix,iy,iz,i1)
        end do
        end do
        end do
        call timer_end(LOG_SS_FDTD_CALC)

        call timer_begin(LOG_SS_FDTD_COMM)
        call update_overlap_real8(fw%srg_eg, fw%eg, fw%box)
        call timer_end(LOG_SS_FDTD_COMM)

        call timer_begin(LOG_SS_FDTD_CALC)
        if(mg%is(3)==lg%is(3))then
    !$OMP parallel do collapse(2) private(ix,iy,iz)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            fw%box(ix,iy,lg%is(3)-1) = fw%vec_Ac_boundary_bottom(ix,iy,i1)
          end do
          end do
        end if
        if(mg%ie(3)==lg%ie(3))then
    !$OMP parallel do collapse(2) private(ix,iy,iz)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            fw%box(ix,iy,lg%ie(3)+1) = fw%vec_Ac_boundary_top(ix,iy,i1)
          end do
          end do
        end if

    !$OMP parallel do collapse(2) private(ix,iy,iz,lap_A)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          lap_A = ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix-1,iy,iz) + fw%box(ix+1,iy,iz) ) / Hgs(1)**2 &
                + ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix,iy-1,iz) + fw%box(ix,iy+1,iz) ) / Hgs(2)**2 &
                + ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix,iy,iz-1) + fw%box(ix,iy,iz+1) ) / Hgs(3)**2
          fw%vec_Ac_m(1,ix,iy,iz,i1) = ( cspeed_au * dt_m )**2 * lap_A &
                                    + 2.d0* fw%box(ix,iy,iz) - fw%vec_Ac_m(-1,ix,iy,iz,i1) &
                                    + dt_m**2 * fw%current4pi(ix,iy,iz,i1)
        end do
        end do
        end do

      end do ! i1 (spacial )

    ! external field
      tm = ( dble(itt-1) + dble(ii-1)/dble(mstep) ) *dt
      call pulse(tm,0d0,   Aext0_old)
      call pulse(tm,Hgs(3),Aext1_old)
      call pulse(tm+dt_m,0d0,   Aext0)
      call pulse(tm+dt_m,Hgs(3),Aext1)
      out_Aext = Aext1

    ! z axis: Mur absorbing boundary condition
      coef = ( cspeed_au * dt_m - Hgs(3) ) / ( cspeed_au * dt_m + Hgs(3) )
      if(mg%is(3)==lg%is(3))then
    !$OMP parallel do collapse(2) private(ix,iy,iz)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
        ! absorbing boundary condition with the incident field vec_Ac_ext
          fw%vec_Ac_boundary_bottom(ix,iy,1:3) = Aext0 &
                                          + ( fw%vec_Ac_m(0,ix,iy,lg%is(3),1:3) - Aext1_old )  &
                                          + coef* ( ( fw%vec_Ac_m(1,ix,iy,lg%is(3),1:3) - Aext1 ) &
                                                  - ( fw%vec_Ac_boundary_bottom_old(ix,iy,1:3) - Aext0_old ) )
        end do
        end do
      end if
      if(mg%ie(3)==lg%ie(3))then
    !$OMP parallel do collapse(2) private(ix,iy,iz)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          fw%vec_Ac_boundary_top(ix,iy,1:3) = fw%vec_Ac_m(0,ix,iy,lg%ie(3),1:3)   &
                                      + coef* ( fw%vec_Ac_m(1,ix,iy,lg%ie(3),1:3) - fw%vec_Ac_boundary_top_old(ix,iy,1:3) )
        end do
        end do
      end if

    !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        fw%vec_Ac_boundary_bottom_old(ix,iy,1:3) = fw%vec_Ac_boundary_bottom(ix,iy,1:3)
        fw%vec_Ac_boundary_top_old   (ix,iy,1:3) = fw%vec_Ac_boundary_top   (ix,iy,1:3)
      end do
      end do

    end do ! ii=1,mstep

    call timer_end(LOG_SS_FDTD_CALC)

    return
  end subroutine fdtd
  
!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine fdtd_gbp
    implicit none

    call timer_begin(LOG_SS_FDTD_CALC)
    fw%tmp_zt = 0d0
    do iz=mg%is(3),mg%ie(3)
  !$omp parallel do collapse(2) private(iy,ix)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        fw%tmp_zt(iz,1:3) = fw%tmp_zt(iz,1:3) + fw%current4pi(ix,iy,iz,1:3) / (lg%num(1)*lg%num(2))
      end do
      end do
    end do
    call timer_end(LOG_SS_FDTD_CALC)

    call timer_begin(LOG_SS_FDTD_COMM_COLL)
    call comm_summation(fw%tmp_zt,fw%curr4pi_zt,size(fw%curr4pi_zt),info%icomm_r)
    call timer_end(LOG_SS_FDTD_COMM_COLL)

    call timer_begin(LOG_SS_FDTD_CALC)
    do ii=1,mstep

    !$OMP parallel do
      do iz=lg%is(3),lg%ie(3)

        fw%Ac_zt_m(iz,-1,1:3) = fw%Ac_zt_m(iz,0,1:3)
        fw%Ac_zt_m(iz,0 ,1:3) = fw%Ac_zt_m(iz,1,1:3)

      end do

      fw%Ac_zt_m(lg%is(3)-1,0,1:3) = fw%Ac_zt_boundary_bottom(1:3)
      fw%Ac_zt_m(lg%ie(3)+1,0,1:3) = fw%Ac_zt_boundary_top   (1:3)

!$OMP parallel do collapse(2) private(i1,iz,lap_A)
      do i1=1,3
        do iz=lg%is(3),lg%ie(3)
          lap_A = ( - 2d0* fw%Ac_zt_m(iz,0,i1) + fw%Ac_zt_m(iz-1,0,i1) + fw%Ac_zt_m(iz+1,0,i1) ) / Hgs(3)**2
          fw%Ac_zt_m(iz,1,i1) = ( cspeed_au * dt_m )**2 * lap_A &
                                    + 2.d0* fw%Ac_zt_m(iz,0,i1) - fw%Ac_zt_m(iz,-1,i1) &
                                    + dt_m**2 * fw%curr4pi_zt(iz,i1)
        end do
      end do ! i1 (spacial )

    ! external field
      tm = ( dble(itt-1) + dble(ii-1)/dble(mstep) ) *dt
      call pulse(tm,0d0,   Aext0_old)
      call pulse(tm,Hgs(3),Aext1_old)
      call pulse(tm+dt_m,0d0,   Aext0)
      call pulse(tm+dt_m,Hgs(3),Aext1)
      out_Aext = Aext1

    ! z axis: Mur absorbing boundary condition
      coef = ( cspeed_au * dt_m - Hgs(3) ) / ( cspeed_au * dt_m + Hgs(3) )

    ! absorbing boundary condition with the incident field vec_Ac_ext
      fw%Ac_zt_boundary_bottom = Aext0 &
                                      + ( fw%Ac_zt_m(lg%is(3),0,1:3) - Aext1_old )  &
                                      + coef* ( ( fw%Ac_zt_m(lg%is(3),1,1:3) - Aext1 ) &
                                              - ( fw%Ac_zt_boundary_bottom_old(1:3) - Aext0_old ) )

      fw%Ac_zt_boundary_top = fw%Ac_zt_m(lg%ie(3),0,1:3)   &
                                  + coef* ( fw%Ac_zt_m(lg%ie(3),1,1:3) - fw%Ac_zt_boundary_top_old(1:3) )

      fw%Ac_zt_boundary_bottom_old(1:3) = fw%Ac_zt_boundary_bottom(1:3)
      fw%Ac_zt_boundary_top_old   (1:3) = fw%Ac_zt_boundary_top   (1:3)

    end do ! ii=1,mstep

  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)

      fw%vec_Ac_m(1,ix,iy,iz,1:3) = fw%Ac_zt_m(iz,1,1:3) + fw%Ac_fourier(ix,iy,iz,1:3)

    end do
    end do
    end do

    if(mg%is(3)==lg%is(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        fw%vec_Ac_boundary_bottom(ix,iy,1:3) = fw%Ac_zt_boundary_bottom
      end do
      end do
    end if

    if(mg%ie(3)==lg%ie(3))then
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        fw%vec_Ac_boundary_top(ix,iy,1:3) = fw%Ac_zt_boundary_top
      end do
      end do
    end if
    call timer_end(LOG_SS_FDTD_CALC)

    return
  end subroutine fdtd_gbp

end subroutine fdtd_singlescale

!===================================================================================================================================

subroutine fourier_singlescale(lg,mg,info,fg,rho,j_e,Vh,poisson,singlescale)
  use structures
  use math_constants,only : zi,pi
  use salmon_global,only: dt,yn_gbp_fourier0
  use phys_constants, only: cspeed_au
  use communication, only: comm_summation,comm_bcast
  implicit none
  type(s_rgrid)          ,intent(in) :: lg
  type(s_rgrid)          ,intent(in) :: mg
  type(s_parallel_info)  ,intent(in) :: info
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)         ,intent(in) :: rho
  type(s_vector)         ,intent(in) :: j_e ! electron number current density (without rho*A/c)
  type(s_scalar)                     :: Vh
  type(s_poisson)                    :: poisson
  type(s_singlescale)                :: singlescale
  !
  integer :: ix,iy,iz
  integer :: iiy,iiz,iix
  real(8) :: inv_lgnum3,vec_je(3),rho_t
  integer :: i
  real(8) :: sin,cos
  complex(8) :: j0,F0,c0,s0,F_old,c_old,s_old,Ac
  
  if(info%isize_x < 4) stop "isize(1) must be > 3"
  
  i = mod(info%id_x,4) ! i=0,1,2,3

  inv_lgnum3=1.d0/(lg%num(1)*lg%num(2)*lg%num(3))

  singlescale%b_ffte=0.d0
  !$OMP parallel do collapse(2) private(iiz,iiy,ix,iy,iz,vec_je,rho_t)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    do ix=mg%is(1),mg%ie(1)
      singlescale%b_ffte(ix,iy,iz,0) = cmplx(rho%f(ix,iiy,iiz)) ! charge density rho
      vec_je = ( j_e%v(1:3,ix,iiy,iiz) + singlescale%vec_je_old(1:3,ix,iiy,iiz) )*0.5d0 ! j(t) = ( j(t+dt/2) + j(t-dt/2) )/2
      rho_t  = ( rho%f(ix,iiy,iiz) + singlescale%rho_old(ix,iiy,iiz) )*0.5d0 ! rho(t) = ( rho(t+dt/2) + rho(t-dt/2) )/2
      vec_je = vec_je + rho_t * singlescale%vec_Ac_m(1,ix,iiy,iiz,1:3) ! electron number current density
      singlescale%b_ffte(ix,iy,iz,1) = cmplx(vec_je(1)) ! current density j_x
      singlescale%b_ffte(ix,iy,iz,2) = cmplx(vec_je(2)) ! current density j_y
      singlescale%b_ffte(ix,iy,iz,3) = cmplx(vec_je(3)) ! current density j_z
    end do
  end do
  end do
  
  call comm_summation(singlescale%b_ffte,singlescale%a_ffte,size(singlescale%a_ffte),info%icomm_x)

  CALL PZFFT3DV_MOD(singlescale%a_ffte(:,:,:,i),singlescale%b_ffte(:,:,:,i),lg%num(1),lg%num(2),lg%num(3),   &
                    info%isize_y,info%isize_z,-1, &
                    info%icomm_y,info%icomm_z)

  call comm_bcast(singlescale%b_ffte(:,:,:,0),info%icomm_x, 0)

! Poisson eq.: singlescale%b_ffte(ix,iy,iz,0)=rho(G) --> poisson%b_ffte(ix,iy,iz)=Vh(G)
  poisson%zrhoG_ele=0d0
  !$omp parallel do collapse(2) default(none) &
  !$omp             private(iz,iy,ix,iiy,iiz,iix) &
  !$omp             shared(mg,lg,poisson,singlescale,inv_lgnum3,fg)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    do ix=1,mg%num(1)
      iix=ix+mg%is(1)-1
      poisson%zrhoG_ele(iix,iiy,iiz) = singlescale%b_ffte(iix,iy,iz,0)*inv_lgnum3
    end do
    do ix=1,lg%num(1)
      poisson%b_ffte(ix,iy,iz) = singlescale%b_ffte(ix,iy,iz,0)*fg%coef(ix,iiy,iiz)
    end do
  end do
  end do
  !$omp end parallel do

! Maxwell eq.: singlescale%b_ffte(ix,iy,iz,i)=j(G,t) --> singlescale%b_ffte(ix,iy,iz,i)=Ac(G,t+dt)
  if(i/=0) then
    !$omp parallel do collapse(2) private(iz,iy,ix,iiz,iiy,j0,F0,c0,s0,F_old,c_old,s_old,Ac,sin,cos)
    do iz=1,mg%num(3)
    do iy=1,mg%num(2)
    do ix=1,lg%num(1)
      iiz=iz+mg%is(3)-1
      iiy=iy+mg%is(2)-1
    ! j(transverse) = j - (1/(4*pi))* d(grad(phi))/dt
      j0 = singlescale%b_ffte(ix,iy,iz,i) &
      & - (1d0/(4d0*pi))* fg%coef_nabla(ix,iiy,iiz,i) * ( poisson%b_ffte(ix,iy,iz) - singlescale%Vh_ffte_old(ix,iy,iz) )/dt
    ! F_{i} = - (4*pi/(c*G)**2)* j(t)
      F0 = - (1/cspeed_au**2)*fg%coef(ix,iiy,iiz) * j0
    ! old variables
      F_old = singlescale%zf_old(ix,iy,iz,i)
      c_old = singlescale%zc_old(ix,iy,iz,i)
      s_old = singlescale%zs_old(ix,iy,iz,i)
    ! cos(c*G*dt), sin(c*G*dt)
      cos = fg%cos_cGdt(ix,iiy,iiz)
      sin = fg%sin_cGdt(ix,iiy,iiz)
    ! c_{i} = c_{i-1}* cos(c*G*dt) + s_{i-1}* sin(c*G*dt) - (F_{i}-F_{i-1})
      c0 = c_old* cos + s_old* sin - (F0-F_old)
    ! s_{i} = - c_{i-1}* sin(c*G*dt) + s_{i-1}* cos(c*G*dt)
      s0 = - c_old* sin + s_old* cos
    ! Ac(t+dt) = c_{i}* cos(c*G*dt) + s_{i}* sin(c*G*dt) + F_{i}
      Ac = c0* cos + s0* sin + F0
    ! Ac(gx=gy=0) --> 0
      singlescale%b_ffte(ix,iy,iz,i) = Ac * fg%coef_gxgy0(ix,iiy,iiz)
    ! F_old,c_old,s_old for next step
      singlescale%zf_old(ix,iy,iz,i) = F0
      singlescale%zc_old(ix,iy,iz,i) = c0
      singlescale%zs_old(ix,iy,iz,i) = s0
    end do
    end do
    end do
    !$omp end parallel do
  end if
  
  !$omp parallel do collapse(2) private(iz,iy,ix,iiy,iiz,iix)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
  do ix=1,lg%num(1)
    singlescale%b_ffte(ix,iy,iz,0)    = poisson%b_ffte(ix,iy,iz)
    singlescale%Vh_ffte_old(ix,iy,iz) = poisson%b_ffte(ix,iy,iz)
  end do
  end do
  end do
  !$omp end parallel do

  CALL PZFFT3DV_MOD(singlescale%b_ffte(:,:,:,i),singlescale%a_ffte(:,:,:,i),lg%num(1),lg%num(2),lg%num(3), &
                    info%isize_y,info%isize_z,1, &
                    info%icomm_y,info%icomm_z)

  call comm_bcast(singlescale%a_ffte(:,:,:,0),info%icomm_x, 0)
  call comm_bcast(singlescale%a_ffte(:,:,:,1),info%icomm_x, 1)
  call comm_bcast(singlescale%a_ffte(:,:,:,2),info%icomm_x, 2)
  call comm_bcast(singlescale%a_ffte(:,:,:,3),info%icomm_x, 3)
  
  !$OMP parallel do private(iiz,iiy,ix,iy,iz)
  do iz=1,mg%num(3)
  do iy=1,mg%num(2)
    iiz=iz+mg%is(3)-1
    iiy=iy+mg%is(2)-1
    Vh%f(mg%is(1):mg%ie(1),iiy,iiz) = singlescale%a_ffte(mg%is(1):mg%ie(1),iy,iz,0)
    singlescale%Ac_fourier(mg%is(1):mg%ie(1),iiy,iiz,1:3) = singlescale%a_ffte(mg%is(1):mg%ie(1),iy,iz,1:3)
  end do
  end do
  
  if(yn_gbp_fourier0=='y') then
    singlescale%Ac_fourier = 0d0
  end if

  return
end subroutine fourier_singlescale

!===================================================================================================================================

subroutine init_singlescale(mg,lg,info,hgs,rho,Vh,srg_scalar,fw,Ac,div_Ac)
  use structures
  use sendrecv_grid, only: update_overlap_real8
  use stencil_sub, only: calc_gradient_field
  use salmon_global, only: sysname,base_directory,yn_restart,yn_ffte,yn_gbp
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  use initialization_sub, only: set_bn
  use filesystem, only: open_filehandle
  use inputoutput, only: t_unit_time
  use checkpoint_restart_sub, only: restart_singlescale
  use sendrecv_grid, only: init_sendrecv_grid
  implicit none
  type(s_rgrid)         ,intent(in) :: lg,mg
  type(s_parallel_info) ,intent(in) :: info
  real(8)               ,intent(in) :: hgs(3)
  type(s_scalar)        ,intent(in) :: rho,Vh ! electron number density & Hartree potential
  type(s_sendrecv_grid)             :: srg_scalar
  type(s_singlescale)               :: fw
  type(s_vector)                    :: Ac
  type(s_scalar)                    :: div_Ac
  !
  character(100) :: filename
  integer :: ii,jj,ix,iy,iz
  real(8) :: bnmat(4,4)
  
  call allocate_scalar(mg,div_Ac)
  call allocate_vector(mg,Ac)

  fw%Energy_poynting = 0d0
  fw%Energy_joule = 0d0

  call set_bn(bnmat)
  do jj=1,3
    do ii=1,4
      fw%coef_nab(ii,jj) = bnmat(ii,4)/hgs(jj)
    end do
  end do

  allocate( fw%vec_Ac_old (3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) )
  allocate( fw%curr         (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),3) )
  allocate( fw%vec_je_old (3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) )
  allocate( fw%rho_old      (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) )
  allocate( fw%current4pi   (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),3) )
  allocate( fw%grad_Vh    (3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) )
  allocate( fw%grad_Vh_old(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) )

!1st element: time step (-1->m-1, 0->m, 1->m+1)
!5th element: components of A vector (1->Ax, 2->Ay, 3->Az)
  allocate( fw%vec_Ac_m(-1:1,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:3) )

  allocate( fw%vec_Ac_boundary_bottom(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3) &
           ,fw%vec_Ac_boundary_bottom_old(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3) &
           ,fw%vec_Ac_boundary_top(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3) &
           ,fw%vec_Ac_boundary_top_old(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3) )

  allocate(fw%integral_poynting(lg%is(3):lg%ie(3)))

  allocate(fw%rot_Ac    (3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) &
   & ,fw%poynting_vector(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) &
        & ,fw%div_Ac      (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) &
        & ,fw%div_Ac_old  (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) &
        & ,fw%integral_poynting_tmp(lg%num(3)),fw%integral_poynting_tmp2(lg%num(3)) &
!        & ,fw%Ac_zt(3,lg%is(3):lg%ie(3)) &
        & ,fw%tmp_zt(lg%is(3):lg%ie(3),3))

  fw%vec_Ac_old = 0d0
  fw%vec_Ac_m = 0d0
  fw%vec_Ac_boundary_bottom = 0d0
  fw%vec_Ac_boundary_bottom_old = 0d0
  fw%vec_Ac_boundary_top = 0d0
  fw%vec_Ac_boundary_top_old = 0d0
  fw%div_Ac_old = 0d0
  fw%integral_poynting = 0d0
  fw%curr = 0d0
  fw%vec_je_old = 0d0
!  fw%Ac_zt = 0d0
  fw%tmp_zt = 0d0
  
! gbp
  if(yn_gbp=='y') then
    allocate( fw%curr4pi_zt(lg%is(3):lg%ie(3),3) )
    allocate(fw%Ac_zt_m(lg%is(3)-1:lg%ie(3)+1,-1:1,1:3))
    allocate(fw%Ac_fourier(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),3))
    allocate(fw%a_ffte(lg%num(1),mg%num(2),mg%num(3),0:3),fw%b_ffte(lg%num(1),mg%num(2),mg%num(3),0:3))
    allocate(fw%Vh_ffte_old(lg%num(1),mg%num(2),mg%num(3)))
    allocate(fw%zf_old(lg%num(1),mg%num(2),mg%num(3),0:3), &
           & fw%zc_old(lg%num(1),mg%num(2),mg%num(3),0:3), &
           & fw%zs_old(lg%num(1),mg%num(2),mg%num(3),0:3))
    fw%curr4pi_zt = 0d0
    fw%Ac_zt_m = 0d0
    fw%Ac_zt_boundary_bottom = 0d0
    fw%Ac_zt_boundary_top = 0d0
    fw%Ac_zt_boundary_bottom_old = 0d0
    fw%Ac_zt_boundary_top_old = 0d0
    fw%Ac_fourier = 0d0
    fw%zf_old = 0d0
    fw%zc_old = 0d0
    fw%zs_old = 0d0
    if(yn_ffte=='y') then
    ! Vh --> fw%Vh_ffte_old
      call calc_Vh_ffte
    end if
  end if

  if(comm_is_root(nproc_id_global)) then
    write(filename,"(2A,'_rt_micro.data')") trim(base_directory),trim(SYSname)
    fw%fh_rt_micro = open_filehandle(filename)
    write(fw%fh_rt_micro, '("#",99(1X,I0,":",A,"[",A,"]"))') &
      & 1, "time", trim(t_unit_time%name), &
      & 2, "Ac_tot_x(z=0)", "a.u.", &
      & 3, "Ac_tot_y(z=0)", "a.u.", &
      & 4, "Ac_tot_z(z=0)", "a.u.", &
      & 5, "Ac_tot_x(z=L)", "a.u.", &
      & 6, "Ac_tot_y(z=L)", "a.u.", &
      & 7, "Ac_tot_z(z=L)", "a.u.", &
      & 8, "Ac_ext_x",      "a.u.", &
      & 9, "Ac_ext_y",      "a.u.", &
     & 10, "Ac_ext_z",      "a.u.", &
     & 11, "J_x(w/o Ac)",   "a.u.", &
     & 12, "J_y(w/o Ac)",   "a.u.", &
     & 13, "J_z(w/o Ac)",   "a.u.", &
     & 14, "E_electron",    "a.u.", &
     & 15, "E_poynting(z=0)", "a.u.", &
     & 16, "E_poynting(z=L)", "a.u.", &
     & 17, "E_em",            "a.u.", &
     & 18, "E_joule",         "a.u."

!  ! for spatial distribution of excitation energy
!    write(filename,"(2A,'_excitation.data')") trim(base_directory),trim(SYSname)
!    fw%fh_excitation = open_filehandle(filename)

!  ! for the vector potential Ac(z,t)
!    write(filename,"(2A,'_Ac_zt.data')") trim(base_directory),trim(SYSname)
!    fw%fh_Ac_zt = open_filehandle(filename)
  end if
  
  allocate(fw%box1(mg%is_array(1):mg%ie_array(1), &
  & mg%is_array(2):mg%ie_array(2), &
  & mg%is_array(3):mg%ie_array(3)))
  fw%box1 = 0d0

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    fw%box1(ix,iy,iz) = Vh%f(ix,iy,iz)
    fw%rho_old(ix,iy,iz) = rho%f(ix,iy,iz)
  end do
  end do
  end do
  if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, fw%box1)
  call calc_gradient_field(mg,fw%coef_nab,fw%box1,fw%grad_Vh_old)
  
  if(yn_restart=='y') then
    call restart_singlescale(info%icomm_rko,lg,mg,fw,Ac,div_Ac)
  end if
  
! specialized in FDTD timestep
  fw%eg%nd = 1
  fw%eg%is = mg%is
  fw%eg%ie = mg%ie
  fw%eg%is_array = mg%is - 1
  fw%eg%ie_array = mg%ie + 1
  call init_sendrecv_grid(fw%srg_eg, fw%eg, 1, srg_scalar%icomm, srg_scalar%neig)
  allocate(fw%box(fw%eg%is_array(1):fw%eg%ie_array(1), &
                & fw%eg%is_array(2):fw%eg%ie_array(2), &
                & fw%eg%is_array(3):fw%eg%ie_array(3)))

  return
  
contains

  subroutine calc_Vh_ffte
    use communication, only: comm_summation
    implicit none
    integer :: ix,iy,iz
    integer :: iiy,iiz,iix

    fw%b_ffte = 0d0
  !$OMP parallel do private(iiz,iiy,ix)
    do iz=1,mg%num(3)
    do iy=1,mg%num(2)
      iiz=iz+mg%is(3)-1
      iiy=iy+mg%is(2)-1
      fw%b_ffte(mg%is(1):mg%ie(1),iy,iz,0) = Vh%f(mg%is(1):mg%ie(1),iiy,iiz)
    end do
    end do
    call comm_summation(fw%b_ffte,fw%a_ffte,size(fw%a_ffte),info%icomm_x)

    CALL PZFFT3DV_MOD(fw%a_ffte(:,:,:,0),fw%Vh_ffte_old,lg%num(1),lg%num(2),lg%num(3),   &
                      info%isize_y,info%isize_z,-1, &
                      info%icomm_y,info%icomm_z)
  
    return
  end subroutine calc_Vh_ffte
  
end subroutine init_singlescale

!===================================================================================================================================

! line integral \int_{{\bf r}_0 \rightarrow {\bf r}_1} {\bf A}({\bf r}) \cdot d{\bf x}
! path: r0 --> r1 = (ix1*Hx,iy1*Hy,iz1*Hz)
subroutine line_integral(integral,r0,A,nx,ny,nz,ix1,iy1,iz1,Hx,Hy,Hz &
                        ,A_lerp,line,wrk,n_max,index)
  implicit none
  integer,intent(in)  :: nx,ny,nz,ix1,iy1,iz1,n_max
  real(8),intent(in)  :: r0(3),A(3,0:nx-1,0:ny-1,0:nz-1),Hx,Hy,Hz
  real(8),intent(out) :: integral
  integer             :: index(n_max)
  real(8)             :: A_lerp(3,n_max),line(3,n_max),wrk(n_max)
  !
  integer :: ix,iy,iz,ix0,iy0,iz0,i,j,k,n,ixp,iyp,izp
  real(8) :: r1(3),p(3),q(3),r(3),h(3),l(3),t

# define IMOD(ix,iy,iz) modulo(ix,nx),modulo(iy,ny),modulo(iz,nz)

  r1 = (/ ix1*Hx, iy1*Hy, iz1*Hz /)

  ix0 = floor( r0(1)/Hx )
  iy0 = floor( r0(2)/Hy )
  iz0 = floor( r0(3)/Hz )

  h = (/ Hx, Hy, Hz /)

! trilinear interpolation @ r0
  ix = ix0
  iy = iy0
  iz = iz0
  ixp = ix + 1
  iyp = iy + 1
  izp = iz + 1
  p = (/ ix*Hx, iy*Hy, iz*Hz /)
  r = ( r0 - p )/h
  A_lerp(:,1) = A(:,IMOD(ix,iy,iz)) * ( 1d0 - r(1) ) * ( 1d0 - r(2) ) * ( 1d0 - r(3) ) &
              + A(:,IMOD(ixp,iy,iz)) * r(1) * ( 1d0 - r(2) ) * ( 1d0 - r(3) ) &
              + A(:,IMOD(ix,iyp,iz)) * ( 1d0 - r(1) ) * r(2) * ( 1d0 - r(3) ) &
              + A(:,IMOD(ix,iy,izp)) * ( 1d0 - r(1) ) * ( 1d0 - r(2) ) * r(3) &
              + A(:,IMOD(ix,iyp,izp)) * ( 1d0 - r(1) ) * r(2) * r(3) &
              + A(:,IMOD(ixp,iy,izp)) * r(1) * ( 1d0 - r(2) ) * r(3) &
              + A(:,IMOD(ixp,iyp,iz)) * r(1) * r(2) * ( 1d0 - r(3) ) &
              + A(:,IMOD(ixp,iyp,izp)) * r(1) * r(2) * r(3)
  line(:,1) = 0d0
  wrk(1) = 0d0

! bilinear interpolation @ intersections
  i = 1
  do ix=min(ix0,ix1),max(ix0,ix1)
    if(ix0==ix1) cycle
    t = ( ix*Hx - r0(1) ) / ( r1(1) - r0(1) )
    if(t<0d0 .or. 1d0<t) cycle
    l = r0 + t * ( r1 - r0 )
    iy = floor( l(2)/Hy )
    iz = floor( l(3)/Hz )
    iyp = iy + 1
    izp = iz + 1
    q = (/ ix*Hx, iy*Hy, iz*Hz /)
    r = ( l - q )/h
    i = i + 1
    if(i>n_max) stop "n_max is too small"
    A_lerp(:,i) = A(:,IMOD(ix,iy,iz)) * ( 1d0 - r(2) ) * ( 1d0 - r(3) ) &
                + A(:,IMOD(ix,iyp,iz)) * r(2) * ( 1d0 - r(3) ) &
                + A(:,IMOD(ix,iy,izp)) * ( 1d0 - r(2) ) * r(3) &
                + A(:,IMOD(ix,iyp,izp)) * r(2) * r(3)
    line(:,i) = l - r0
    wrk(i) = line(1,i)**2 + line(2,i)**2 + line(3,i)**2
  end do
  do iy=min(iy0,iy1),max(iy0,iy1)
    if(iy0==iy1) cycle
    t = ( iy*Hy - r0(2) ) / ( r1(2) - r0(2) )
    if(t<0d0 .or. 1d0<t) cycle
    l = r0 + t * ( r1 - r0 )
    ix = floor( l(1)/Hx )
    iz = floor( l(3)/Hz )
    ixp = ix + 1
    izp = iz + 1
    q = (/ ix*Hx, iy*Hy, iz*Hz /)
    r = ( l - q )/h
    i = i + 1
    if(i>n_max) stop "n_max is too small"
    A_lerp(:,i) = A(:,IMOD(ix,iy,iz)) * ( 1d0 - r(1) ) * ( 1d0 - r(3) ) &
                + A(:,IMOD(ixp,iy,iz)) * r(1) * ( 1d0 - r(3) ) &
                + A(:,IMOD(ix,iy,izp)) * ( 1d0 - r(1) ) * r(3) &
                + A(:,IMOD(ixp,iy,izp)) * r(1) * r(3)
    line(:,i) = l - r0
    wrk(i) = line(1,i)**2 + line(2,i)**2 + line(3,i)**2
  end do
  do iz=min(iz0,iz1),max(iz0,iz1)
    if(iz0==iz1) cycle
    t = ( iz*Hz - r0(3) ) / ( r1(3) - r0(3) )
    if(t<0d0 .or. 1d0<t) cycle
    l = r0 + t * ( r1 - r0 )
    ix = floor( l(1)/Hx )
    iy = floor( l(2)/Hy )
    ixp = ix + 1
    iyp = iy + 1
    q = (/ ix*Hx, iy*Hy, iz*Hz /)
    r = ( l - q )/h
    i = i + 1
    if(i>n_max) stop "n_max is too small"
    A_lerp(:,i) = A(:,IMOD(ix,iy,iz)) * ( 1d0 - r(1) ) * ( 1d0 - r(2) ) &
                + A(:,IMOD(ixp,iy,iz)) * r(1) * ( 1d0 - r(2) ) &
                + A(:,IMOD(ix,iyp,iz)) * ( 1d0 - r(1) ) * r(2) &
                + A(:,IMOD(ixp,iyp,iz)) * r(1) * r(2)
    line(:,i) = l - r0
    wrk(i) = line(1,i)**2 + line(2,i)**2 + line(3,i)**2
  end do
  if(i==1) then
    i = i + 1
    if(i>n_max) stop "n_max is too small"
    A_lerp(:,i) = A(:,IMOD(ix0,iy0,iz0))
    line(:,i) = p - r0
    wrk(i) = line(1,i)**2 + line(2,i)**2 + line(3,i)**2
  end if
  n = i

  call heapsort(n,wrk,index)

  integral = 0d0
  do i=2,n
    j = index(i)
    k = index(i-1)
    integral = integral + 0.5d0* ( A_lerp(1,j) + A_lerp(1,k) ) &
                        * ( line(1,j) - line(1,k) ) &
                        + 0.5d0* ( A_lerp(2,j) + A_lerp(2,k) ) &
                        * ( line(2,j) - line(2,k) ) &
                        + 0.5d0* ( A_lerp(3,j) + A_lerp(3,k) ) &
                        * ( line(3,j) - line(3,k) )
  end do

  return
contains

  subroutine heapsort(n,array,index)
    implicit none
    integer,intent(in)    :: n
    real(8),intent(inout) :: array(n)
    integer,intent(out)   :: index(n)
    !
    integer :: i,k,j,l,m
    real(8) :: t

    if(n < 2) return

    do i=1,n
      index(i) = i
    end do

    l = n/2 + 1
    k = n
    do while(k /= 1)
      if(l > 1)then
        l = l - 1
        t = array(l)
        m = index(l)
      else
        t = array(k)
        m = index(k)
        array(k) = array(1)
        index(k) = index(1)
        k = k - 1
        if(k == 1) then
          array(1) = t
          index(1) = m
          exit
        end if
      end if
      i = l
      j = l + l
      do while(j <= k)
        if(j < k) then
          if(array(j) < array(j+1)) j = j + 1
        endif
        if (t < array(j))then
          array(i) = array(j)
          index(i) = index(j)
          i = j
          j = j + j
        else
          j = k + 1
        end if
      end do
      array(i) = t
      index(i) = m
    end do

    return
  end subroutine heapsort
end subroutine line_integral

end module fdtd_coulomb_gauge
