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
  type ls_singlescale
    integer :: fh_rt_micro,fh_excitation,fh_Ac_zt
    real(8) :: E_electron,Energy_poynting(2),coef_nab(4,3),bnmat(4,4)
    real(8),allocatable :: vec_Ac(:,:,:,:),vec_Ac_stock(:,:,:,:),Vh_n(:,:,:),curr1_m(:,:,:,:),vec_Ac_m(:,:,:,:,:) &
    & ,vec_Ac_boundary_bottom(:,:,:),vec_Ac_boundary_bottom_old(:,:,:),vec_Ac_boundary_top(:,:,:),vec_Ac_boundary_top_old(:,:,:) &
    & ,integral_poynting(:),Ac_zt(:,:)
    real(8),allocatable :: box(:,:,:),grad_Vh(:,:,:,:),gradient_V(:,:,:,:),rotation_A(:,:,:,:),poynting_vector(:,:,:,:) &
    & ,divergence_A(:,:,:),vbox(:,:,:,:),lgbox1(:,:,:),lgbox2(:,:,:),integral_poynting_tmp(:),integral_poynting_tmp2(:) &
    & ,vec_Ac_ext(:,:,:,:),vec_Ac_ext_old(:,:,:,:)
  end type ls_singlescale

contains

subroutine fdtd_singlescale(itt,lg,mg,ng,hgs,rho,Vh,j_e,srg_ng,Ac,div_Ac,fw)
  use structures
  use math_constants,only : zi,pi
  use phys_constants, only: cspeed_au
  use salmon_global, only: dt
  use sendrecv_grid, only: update_overlap_real8
  use salmon_parallel, only: nproc_id_global, nproc_group_global
  use salmon_communication, only: comm_is_root, comm_summation
  use inputoutput, only: t_unit_time
  implicit none
  integer       ,intent(in) :: itt
  type(s_rgrid) ,intent(in) :: lg,mg,ng
  real(8)       ,intent(in) :: hgs(3)
  type(s_scalar),intent(in) :: rho,Vh ! electron number density & Hartree potential
  type(s_vector),intent(in) :: j_e    ! electron number current density
  type(s_sendrecv_grid)     :: srg_ng
  type(s_vector)            :: Ac     ! A/c, A: vector potential, c: speed of light
  type(s_scalar)            :: div_Ac ! div(A/c)
  type(ls_singlescale)      :: fw     ! FDTD working arrays, etc.
  !
  integer,parameter :: mstep=100
  integer,parameter :: Nd = 4
  integer,dimension(3) :: ng_sta,ng_end,ng_num,mg_sta,mg_end,mg_num,lg_sta,lg_end,lg_num
  integer :: ix,iy,iz,i1,ii,krd(3,3),lcs(3,3,3),dr(3),comm

  real(8) :: Hvol,dt_m,tm,coef,lap_A,Energy_em,Energy_joule,diff_A,coef2 &
  & ,e_em,e_em_wrk,e_joule,e_joule_wrk,e_poynting(2),e_poynting_wrk(2)
  real(8),dimension(3) :: out_curr,out_Aext,out_Ab1,out_Ab2,wrk,wrk2,wrk4

  comm = nproc_group_global ! for comm_summation: ng --> lg

  krd = 0
  krd(1,1) = 1; krd(2,2) = 1; krd(3,3) = 1

  lcs = 0
  lcs(1,2,3) = 1; lcs(3,1,2) = 1; lcs(2,3,1) = 1
  lcs(1,3,2) = -1; lcs(2,1,3) = -1; lcs(3,2,1) = -1

  Hvol = hgs(1)*hgs(2)*hgs(3)
  dt_m = dt / dble(mstep)

  ng_sta = ng%is
  ng_end = ng%ie
  ng_num = ng%num
  mg_sta = mg%is
  mg_end = mg%ie
  mg_num = mg%num
  lg_sta = lg%is
  lg_end = lg%ie
  lg_num = lg%num

  if(.not.allocated(fw%vec_Ac)) then
    call init(ng_sta,ng_end,lg_sta,lg_end,hgs,fw)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      fw%Vh_n(ix,iy,iz) = Vh%f(ix,iy,iz)
    end do
    end do
    end do
  end if

!-----------------------------------------------------------------------------------------------------------------------------------

! calculate grad_Vh: gradient of d(Vh)/dt (Vh: Hartree potential)

  !$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    fw%box(ix,iy,iz) = ( Vh%f(ix,iy,iz) - fw%Vh_n(ix,iy,iz) ) /dt ! t differential ! Vh = V_H(t), Vh_n = V_H(t-dt)
  end do
  end do
  end do

  call update_overlap_real8(srg_ng, ng, fw%box)
  call calc_gradient(ng_sta,ng_end,fw%coef_nab,fw%box,fw%grad_Vh) ! grad_Vh: grad( d(Vh)/dt )

!-----------------------------------------------------------------------------------------------------------------------------------

! for electric field E

  !$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    fw%box(ix,iy,iz) = - Vh%f(ix,iy,iz) ! Vh_wk: scalar potential, Vh: Hartree potential
    fw%Vh_n(ix,iy,iz) = Vh%f(ix,iy,iz)  ! old hartree potential
  end do
  end do
  end do
  call update_overlap_real8(srg_ng, ng, fw%box)
  call calc_gradient(ng_sta,ng_end,fw%coef_nab,fw%box,fw%gradient_V)

!-----------------------------------------------------------------------------------------------------------------------------------

  fw%vbox = 0d0
  wrk = 0d0
!$OMP parallel do collapse(2) private(ix,iy,iz) reduction(+:wrk)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    fw%curr1_m(ix,iy,iz,1:3) = j_e%v(1:3,ix,iy,iz) + rho%f(ix,iy,iz) * fw%vec_Ac_m(1,ix,iy,iz,1:3) ! j(t)
    wrk = wrk + j_e%v(1:3,ix,iy,iz) ! definition of out_curr. j_e%v --> fw%curr1_m ?
  end do
  end do
  end do
  wrk = wrk/dble(lg_num(1)*lg_num(2)*lg_num(3))
  call comm_summation(wrk,out_curr,3,comm)

!-----------------------------------------------------------------------------------------------------------------------------------

  ! FDTD loop: A(t) --> A(t+dt)

  fw%rotation_A = 0d0
  fw%divergence_A = 0d0
  do ii=1,mstep

  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)

      fw%vec_Ac_m(-1,ix,iy,iz,1:3) = fw%vec_Ac_m(0,ix,iy,iz,1:3)
      fw%vec_Ac_m(0 ,ix,iy,iz,1:3) = fw%vec_Ac_m(1,ix,iy,iz,1:3)

    end do
    end do
    end do

    do i1=1,3

  !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        fw%box(ix,iy,iz) = fw%vec_Ac_m(0,ix,iy,iz,i1)
      end do
      end do
      end do

      call update_overlap_real8(srg_ng, ng, fw%box)
      if(ng_sta(3)==lg_sta(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          fw%box(ix,iy,lg_sta(3)-1) = fw%vec_Ac_boundary_bottom(ix,iy,i1)
        end do
        end do
      end if
      if(ng_end(3)==lg_end(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          fw%box(ix,iy,lg_end(3)+1) = fw%vec_Ac_boundary_top(ix,iy,i1)
        end do
        end do
      end if

  !$OMP parallel do collapse(2) private(ix,iy,iz,lap_A)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          lap_A = ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix-1,iy,iz) + fw%box(ix+1,iy,iz) ) / Hgs(1)**2 &
                + ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix,iy-1,iz) + fw%box(ix,iy+1,iz) ) / Hgs(2)**2 &
                + ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix,iy,iz-1) + fw%box(ix,iy,iz+1) ) / Hgs(3)**2
          fw%vec_Ac_m(1,ix,iy,iz,i1) = ( cspeed_au * dt_m )**2 * lap_A &
                                    + 2.d0* fw%box(ix,iy,iz) - fw%vec_Ac_m(-1,ix,iy,iz,i1) &
                                    + dt_m**2 * ( fw%grad_Vh(i1,ix,iy,iz) - 4d0*pi * fw%curr1_m(ix,iy,iz,i1) )
        end do
        end do
        end do

  !   rotation & divergence of A
      if(ii==mstep/2) then
  !$OMP parallel do collapse(2) private(ix,iy,iz,wrk,dr,diff_A)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
        ! rot(A)
          wrk(1) = ( fw%box(ix+1,iy,iz) - fw%box(ix-1,iy,iz) ) / ( 2d0* Hgs(1) )
          wrk(2) = ( fw%box(ix,iy+1,iz) - fw%box(ix,iy-1,iz) ) / ( 2d0* Hgs(2) )
          wrk(3) = ( fw%box(ix,iy,iz+1) - fw%box(ix,iy,iz-1) ) / ( 2d0* Hgs(3) )
          fw%rotation_A(:,ix,iy,iz) = fw%rotation_A(:,ix,iy,iz) + lcs(:,1,i1) * wrk(1) + lcs(:,2,i1) * wrk(2) + lcs(:,3,i1) * wrk(3)

        ! div(A)
          dr = krd(:,i1)
          diff_A = ( fw%box(ix+dr(1),iy+dr(2),iz+dr(3)) - fw%box(ix-dr(1),iy-dr(2),iz-dr(3)) ) / ( 2d0* Hgs(i1) )
          fw%divergence_A(ix,iy,iz) = fw%divergence_A(ix,iy,iz) + diff_A
        end do
        end do
        end do
      end if

    end do ! i1 (spacial )

  ! external field
    tm = ( dble(itt-1) + dble(ii-1)/dble(mstep) ) *dt
    call pulse(tm,0d0,wrk)
  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      fw%vec_Ac_ext_old(ix,iy,0,1:3) = wrk
    end do
    end do
    call pulse(tm,Hgs(3),wrk)
  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      fw%vec_Ac_ext_old(ix,iy,1,1:3) = wrk
    end do
    end do
    call pulse(tm+dt_m,0d0,wrk)
  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      fw%vec_Ac_ext(ix,iy,0,1:3) = wrk
    end do
    end do
    call pulse(tm+dt_m,Hgs(3),wrk)
  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      fw%vec_Ac_ext(ix,iy,1,1:3) = wrk
    end do
    end do

    out_Aext = wrk

  ! z axis: Mur absorbing boundary condition
    coef = ( cspeed_au * dt_m - Hgs(3) ) / ( cspeed_au * dt_m + Hgs(3) )
    if(ng_sta(3)==lg_sta(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
      ! absorbing boundary condition with the incident field vec_Ac_ext
        fw%vec_Ac_boundary_bottom(ix,iy,1:3) = fw%vec_Ac_ext(ix,iy,0,1:3) &
                                        + ( fw%vec_Ac_m(0,ix,iy,lg_sta(3),1:3) - fw%vec_Ac_ext_old(ix,iy,1,1:3) )  &
                                        + coef* ( ( fw%vec_Ac_m(1,ix,iy,lg_sta(3),1:3) - fw%vec_Ac_ext(ix,iy,1,1:3) ) &
                                                - ( fw%vec_Ac_boundary_bottom_old(ix,iy,1:3) - fw%vec_Ac_ext_old(ix,iy,0,1:3) ) )
      end do
      end do
    end if
    if(ng_end(3)==lg_end(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        fw%vec_Ac_boundary_top(ix,iy,1:3) = fw%vec_Ac_m(0,ix,iy,lg_end(3),1:3)   &
                                    + coef* ( fw%vec_Ac_m(1,ix,iy,lg_end(3),1:3) - fw%vec_Ac_boundary_top_old(ix,iy,1:3) )
      end do
      end do
    end if

  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      fw%vec_Ac_boundary_bottom_old(ix,iy,1:3) = fw%vec_Ac_boundary_bottom(ix,iy,1:3)
      fw%vec_Ac_boundary_top_old   (ix,iy,1:3) = fw%vec_Ac_boundary_top   (ix,iy,1:3)
    end do
    end do

  end do ! ii=1,mstep

!-----------------------------------------------------------------------------------------------------------------------------------

  fw%vbox = 0d0

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    fw%vbox(1:3,ix,iy,iz) = ( fw%vec_Ac_m(1,ix,iy,iz,1:3) + fw%vec_Ac_stock(1:3,ix,iy,iz) ) * 0.5d0 ! ( A(t+dt) + A(t) )/2
  end do
  end do
  end do

  call comm_summation(fw%vbox,fw%vec_Ac,3*lg_num(1)*lg_num(2)*lg_num(3),comm)

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    Ac%v(:,ix,iy,iz) = fw%vec_Ac(:,ix,iy,iz)
  end do
  end do
  end do

!-----------------------------------------------------------------------------------------------------------------------------------

  fw%lgbox1 = 0d0

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    fw%lgbox1(ix,iy,iz) = fw%divergence_A(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(fw%lgbox1,fw%lgbox2,lg_num(1)*lg_num(2)*lg_num(3),comm)

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    div_Ac%f(ix,iy,iz) = fw%lgbox2(ix,iy,iz)
  end do
  end do
  end do

!-----------------------------------------------------------------------------------------------------------------------------------

  ! integral(A) @ z = 0 (bottom boundary)
  wrk = 0d0
  wrk(1) = sum(fw%vec_Ac(1,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3)))
  wrk(2) = sum(fw%vec_Ac(2,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3)))
  wrk(3) = sum(fw%vec_Ac(3,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3)))
  out_Ab1 = wrk / dble(lg_num(1)*lg_num(2))

  ! integral(A) @ z = az (top boundary)
  wrk = 0d0
  wrk(1) = sum(fw%vec_Ac(1,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_end(3)))
  wrk(2) = sum(fw%vec_Ac(2,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_end(3)))
  wrk(3) = sum(fw%vec_Ac(3,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_end(3)))
  out_Ab2 = wrk / dble(lg_num(1)*lg_num(2))

  ! Electro-Magnetic energy & Joule dissipated power
  coef = cspeed_au / (4d0*pi)
  coef2 = Hvol / (8d0*pi)
  e_em_wrk = 0d0
  e_joule_wrk = 0d0
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    wrk4 = ( fw%vec_Ac_m(1,ix,iy,iz,:) - fw%vec_Ac_stock(:,ix,iy,iz) ) / dt ! (A(t+dt)-A(t))/dt
    wrk  = - fw%gradient_V(:,ix,iy,iz) - wrk4    ! E
    wrk2 = cspeed_au * fw%rotation_A(:,ix,iy,iz) ! B
    fw%poynting_vector(:,ix,iy,iz) = coef * ( lcs(:,1,2) * wrk(1) * wrk2(2) + lcs(:,1,3) * wrk(1) * wrk2(3) &
                                            + lcs(:,2,1) * wrk(2) * wrk2(1) + lcs(:,2,3) * wrk(2) * wrk2(3) &
                                            + lcs(:,3,1) * wrk(3) * wrk2(1) + lcs(:,3,2) * wrk(3) * wrk2(2) ) ! E x B
    e_em_wrk = e_em_wrk + coef2 * ( sum(wrk**2) + sum(wrk2**2) ) ! ( E^2 + B^2 )/(8*pi)
  end do
  end do
  end do
  call comm_summation(e_em_wrk,e_em,comm)
  call comm_summation(e_joule_wrk,e_joule,comm)

! Surface integral of the poynting vector S
  coef = Hgs(1)*Hgs(2)
  e_poynting_wrk = 0d0
  if(ng_sta(3)==lg_sta(3)) then ! integral(S) @ z = 0 (bottom boundary)
    e_poynting_wrk(1) = sum( fw%poynting_vector(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),lg_sta(3)) ) * coef
  end if
  if(ng_end(3)==lg_end(3)) then ! integral(S) @ z = az (top boundary)
    e_poynting_wrk(2) = sum( fw%poynting_vector(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),lg_end(3)) ) * coef
  end if
  call comm_summation(e_poynting_wrk,e_poynting,2,comm)

  Energy_em = e_em
!  Energy_joule = Energy_joule + dt*e_joule
  fw%Energy_poynting = fw%Energy_poynting + dt*e_poynting

  if(comm_is_root(nproc_id_global)) write(fw%fh_rt_micro,'(99(1X,E23.15E3))') &
  dble(itt)*dt*t_unit_time%conv,out_Ab1,out_Ab2,out_Aext,out_curr,fw%E_electron,fw%Energy_poynting,Energy_em,Energy_joule

! for spatial distribution of excitation energy
  coef = Hgs(1)*Hgs(2)
  fw%integral_poynting_tmp = 0d0
  do iz=ng_sta(3),ng_end(3)
    fw%integral_poynting_tmp(iz) = sum( fw%poynting_vector(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),iz) ) * coef
  end do
  call comm_summation(fw%integral_poynting_tmp,fw%integral_poynting_tmp2,lg_num(3),comm)
  fw%integral_poynting = fw%integral_poynting + dt * fw%integral_poynting_tmp2

! for the vector potential Ax(z,t)
  do iz=lg_sta(3),lg_end(3)
    fw%Ac_zt(:,iz) = sum( fw%vec_Ac(:,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),iz) )/(lg_num(1)*lg_num(2))
  end do
  if(comm_is_root(nproc_id_global)) then
    do iz=lg_sta(3),lg_end(3)
      write(fw%fh_Ac_zt,fmt='(99(1X,E23.15E3))',advance='no') dble(iz)*hgs(3),fw%Ac_zt(1,iz),fw%Ac_zt(2,iz),fw%Ac_zt(3,iz)
    end do
    write(fw%fh_Ac_zt,'()')
  end if

!-----------------------------------------------------------------------------------------------------------------------------------

  ! stock vec_Ac

  !$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    fw%vec_Ac_stock(1,ix,iy,iz) = fw%vec_Ac_m(1,ix,iy,iz,1) ! A(t+dt) --> A(t)
    fw%vec_Ac_stock(2,ix,iy,iz) = fw%vec_Ac_m(1,ix,iy,iz,2)
    fw%vec_Ac_stock(3,ix,iy,iz) = fw%vec_Ac_m(1,ix,iy,iz,3)
  end do
  end do
  end do

  return

  contains

# define DX(dt) (ix+(dt)),iy,iz
# define DY(dt) ix,(iy+(dt)),iz
# define DZ(dt) ix,iy,(iz+(dt))

  subroutine calc_gradient(is,ie,nabt,box,grad)
    implicit none
    integer      ,intent(in) :: is(3),ie(3)
    real(8)      ,intent(in) :: nabt(4,3)
    real(8)      ,intent(in) :: box(is(1)-Nd:ie(1)+Nd,is(2)-Nd:ie(2)+Nd,is(3)-Nd:ie(3)+Nd)
    real(8)                  :: grad(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
    !
    integer :: ix,iy,iz
    real(8) :: w(3)
  !$OMP parallel
  !$OMP do private(iz,iy,ix,w)
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      w(1) =  nabt(1,1)*(box(DX(1)) - box(DX(-1))) &
           & +nabt(2,1)*(box(DX(2)) - box(DX(-2))) &
           & +nabt(3,1)*(box(DX(3)) - box(DX(-3))) &
           & +nabt(4,1)*(box(DX(4)) - box(DX(-4)))
      w(2) =  nabt(1,2)*(box(DY(1)) - box(DY(-1))) &
           & +nabt(2,2)*(box(DY(2)) - box(DY(-2))) &
           & +nabt(3,2)*(box(DY(3)) - box(DY(-3))) &
           & +nabt(4,2)*(box(DY(4)) - box(DY(-4)))
      w(3) =  nabt(1,3)*(box(DZ(1)) - box(DZ(-1))) &
           & +nabt(2,3)*(box(DZ(2)) - box(DZ(-2))) &
           & +nabt(3,3)*(box(DZ(3)) - box(DZ(-3))) &
           & +nabt(4,3)*(box(DZ(4)) - box(DZ(-4)))
      grad(:,ix,iy,iz) = w
    end do
    end do
    end do
  !$OMP end do
  !$OMP end parallel
  end subroutine calc_gradient

  subroutine pulse(t,r,A_ext)
    use salmon_global, only: pulse_tw1,omega1,phi_CEP1,amplitude1,rlaser_int_wcm2_1,epdir_re1,ae_shape1
    implicit none
    real(8),intent(in)  :: t,r
    real(8),intent(out) :: A_ext(3)
    !
    real(8) :: wrk,theta1,theta2,rr,amplitude

    if(amplitude1/=0d0) then
      amplitude = amplitude1
    else
      amplitude = sqrt(rlaser_int_wcm2_1)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    end if

    rr = r - cspeed_au * t

    rr = rr + 0.5d0*pulse_tw1*cspeed_au

  !  theta1 = Pi/pulse_tw1*(tt-0.5d0*pulse_tw1)
    theta1 = pi * rr / (pulse_tw1*cspeed_au)
  !  theta2 = omega1*(tt-0.5d0*pulse_tw1)+phi_cep1*2d0*pi
    theta2 = omega1 * rr / cspeed_au + phi_CEP1*2d0*pi

    if(ae_shape1=='Acos2')then
      wrk = cos(theta1)**2 * aimag(exp(zI*theta2)) / omega1
    end if

    A_ext = 0d0
    if(abs(rr) < 0.5d0*pulse_tw1*cspeed_au)then
      A_ext(1:3) = amplitude * epdir_re1(1:3) * wrk
    end if

    return
  end subroutine pulse

  subroutine init(ng_sta,ng_end,lg_sta,lg_end,hgs,fw)
    use salmon_global, only: sysname,directory
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    use salmon_initialization, only: setbn
    use salmon_file, only: open_filehandle
    use inputoutput, only: t_unit_time
    implicit none
    integer,dimension(3),intent(in) :: ng_sta,ng_end,lg_sta,lg_end
    real(8)  ,intent(in) :: hgs(3)
    type(ls_singlescale) :: fw
    !
    integer,parameter :: Nd=4
    character(100) :: filename
    integer :: lg_num(3),ii,jj
    lg_num = lg_end - lg_sta + 1

    fw%Energy_poynting = 0d0

    call setbn(fw%bnmat)
    do jj=1,3
      do ii=1,4
        fw%coef_nab(ii,jj) = fw%bnmat(ii,4)/hgs(jj)
      end do
    end do

    allocate( fw%vec_Ac(3,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) )
    allocate( fw%vec_Ac_stock(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) )
    allocate( fw%Vh_n      (  ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) )
    allocate( fw%curr1_m   (  ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),3))

  !1st element: time step (-1->m-1, 0->m, 1->m+1)
  !5th element: components of A vector (1->Ax, 2->Ay, 3->Az)
    allocate( fw%vec_Ac_m(-1:1,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:3) )

    allocate( fw%vec_Ac_boundary_bottom(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),1:3) &
             ,fw%vec_Ac_boundary_bottom_old(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),1:3) &
             ,fw%vec_Ac_boundary_top(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),1:3) &
             ,fw%vec_Ac_boundary_top_old(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),1:3) )

    allocate(fw%integral_poynting(lg_sta(3):lg_end(3)))

    fw%vec_Ac = 0d0
    fw%vec_Ac_stock = 0d0
    fw%vec_Ac_m = 0d0
    fw%vec_Ac_boundary_bottom = 0d0
    fw%vec_Ac_boundary_bottom_old = 0d0
    fw%vec_Ac_boundary_top = 0d0
    fw%vec_Ac_boundary_top_old = 0d0
    fw%integral_poynting = 0d0

  ! temporary arrays
    allocate(fw%box(ng_sta(1)-Nd:ng_end(1)+Nd,ng_sta(2)-Nd:ng_end(2)+Nd,ng_sta(3)-Nd:ng_end(3)+Nd) &
          & ,fw%grad_Vh   (3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
          & ,fw%gradient_V(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
          & ,fw%rotation_A(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
     & ,fw%poynting_vector(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
          & ,fw%divergence_A(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) &
          & ,fw%vbox      (3,lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) &
          & ,fw%lgbox1      (lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) &
          & ,fw%lgbox2      (lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)) &
          & ,fw%integral_poynting_tmp(lg_num(3)),fw%integral_poynting_tmp2(lg_num(3)) &
          & ,fw%vec_Ac_ext    (ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),0:1,1:3) &
          & ,fw%vec_Ac_ext_old(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),0:1,1:3) &
          & ,fw%Ac_zt(3,lg_sta(3):lg_end(3)))

    fw%grad_Vh = 0d0
    fw%gradient_V = 0d0

    if(comm_is_root(nproc_id_global)) then
      write(filename,"(2A,'_rt_micro.data')") trim(directory),trim(SYSname)
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

    ! for spatial distribution of excitation energy
      write(filename,"(2A,'_excitation.data')") trim(directory),trim(SYSname)
      fw%fh_excitation = open_filehandle(filename)

    ! for the vector potential Ac(z,t)
      write(filename,"(2A,'_Ac_zt.data')") trim(directory),trim(SYSname)
      fw%fh_Ac_zt = open_filehandle(filename)
    end if

  end subroutine init

end subroutine fdtd_singlescale

end module fdtd_coulomb_gauge
