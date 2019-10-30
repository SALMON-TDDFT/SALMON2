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
!-----------------------------------------------------------------------------------------
subroutine weyl_calc(fs,ff,fw)
  use inputoutput, only: fdtddim, dt_em
  use structures,     only: s_fdtd_system, s_fdtd_field
  use salmon_maxwell, only: ls_fdtd_work
  use phys_constants,    only: cspeed_au
  use math_constants,    only: pi
  use pack_unpack, only: copy_data
  implicit none
  type(s_fdtd_system) :: fs
  type(s_fdtd_field)  :: ff
  type(ls_fdtd_work)  :: fw

  ! Update electromagnetic field
  select case(trim(fdtddim))
  case('1d', '1D')
    call dt_evolve_Ac_1d()
    call calc_vec_h_1d()
  case('2d', '2D')
    call dt_evolve_Ac_2d()
    call calc_vec_h_2d()
  case('3d', '3D')
    call dt_evolve_Ac_3d()
    call calc_vec_h_3d()
  ! case('3d_cylindal', '2dc', '2DC')
  !   call dt_evolve_Ac_3d_cylindal()
  !   call calc_vec_h_3d_cylindal()
  end select
  call calc_vec_e()

  call calc_energy_density()
  
  return

contains

  subroutine dt_evolve_Ac_1d()
    implicit none
    integer :: i1, i2, i3
    real(8) :: rot2_Ac(3) ! rot rot Ac
    real(8) :: r_inv_h(3)
    i2 = fs%mg%is(2)
    i3 = fs%mg%is(3)
    r_inv_h(1) = 1d0 / fs%hgs(1)
    !$omp parallel do default(shared) private(i1, rot2_Ac)
    do i1 = fs%mg%is(1), fs%mg%ie(1)
      rot2_Ac(1) = 0d0
      rot2_Ac(2:3) = -( &
        &      + ff%vec_Ac%v(2:3, i1+1, i2, i3) &
        & -2d0 * ff%vec_Ac%v(2:3, i1,   i2, i3) &
        &      + ff%vec_Ac%v(2:3, i1-1, i2, i3) &
        & ) * r_inv_h(1) ** 2
      fw%vec_Ac_tmp%v(:,i1, i2, i3) = (2 * ff%vec_Ac%v(:,i1, i2, i3) - ff%vec_Ac_old%v(:,i1, i2, i3) &
        & + ff%vec_j_em%v(:,i1, i2, i3) * 4.0 * pi * (dt_em**2) - rot2_Ac(:) * (cspeed_au * dt_em)**2 )
    end do
    !$omp end parallel do
    ! Temporal data transfer: Ac_tmp -> Ac -> Ac_old
    call copy_data(ff%vec_Ac%v, ff%vec_Ac_old%v)
    call copy_data(fw%vec_Ac_tmp%v, ff%vec_Ac%v)
    return
  end subroutine dt_evolve_Ac_1d


  subroutine dt_evolve_Ac_2d()
    implicit none
    integer :: i1, i2, i3
    real(8) :: rot2_Ac(3) ! rot rot Ac
    real(8) :: r_inv_h(3)
    r_inv_h(1:2) = 1d0 / fs%hgs(1:2)
    i3 = fs%mg%is(3)
    !$omp parallel do collapse(2) default(shared) private(i1, i2, rot2_Ac)
    do i2 = fs%mg%is(2), fs%mg%ie(2)
      do i1 = fs%mg%is(1), fs%mg%ie(1)
        rot2_Ac(1) = +(-1.00d0 * r_inv_h(2)**2) * ff%vec_Ac%v(1, i1, i2-1, i3) &
          & +(+2.00d0 * r_inv_h(2)**2) * ff%vec_Ac%v(1, i1, i2, i3) &
          & +(-1.00d0 * r_inv_h(2)**2) * ff%vec_Ac%v(1, i1, i2+1, i3) &
          & +(+0.25d0 * r_inv_h(1) * r_inv_h(2)) * ff%vec_Ac%v(2, i1-1, i2-1, i3) &
          & +(-0.25d0 * r_inv_h(1) * r_inv_h(2)) * ff%vec_Ac%v(2, i1-1, i2+1, i3) &
          & +(-0.25d0 * r_inv_h(1) * r_inv_h(2)) * ff%vec_Ac%v(2, i1+1, i2-1, i3) &
          & +(+0.25d0 * r_inv_h(1) * r_inv_h(2)) * ff%vec_Ac%v(2, i1+1, i2+1, i3)
        rot2_Ac(2) = +(+0.25d0 * r_inv_h(1) * r_inv_h(2)) * ff%vec_Ac%v(1, i1-1, i2-1, i3) &
          & +(-0.25d0 * r_inv_h(1) * r_inv_h(2)) * ff%vec_Ac%v(1, i1-1, i2+1, i3) &
          & +(-0.25d0 * r_inv_h(1) * r_inv_h(2)) * ff%vec_Ac%v(1, i1+1, i2-1, i3) &
          & +(+0.25d0 * r_inv_h(1) * r_inv_h(2)) * ff%vec_Ac%v(1, i1+1, i2+1, i3) &
          & +(-1.00d0 * r_inv_h(1)**2) * ff%vec_Ac%v(2, i1-1, i2, i3) &
          & +(+2.00d0 * r_inv_h(1)**2) * ff%vec_Ac%v(2, i1, i2, i3) &
          & +(-1.00d0 * r_inv_h(1)**2) * ff%vec_Ac%v(2, i1+1, i2, i3)
        rot2_Ac(3) = +(-1.00d0 * r_inv_h(1)**2) * ff%vec_Ac%v(3, i1-1, i2, i3) &
          & +(-1.00d0 * r_inv_h(2)**2) * ff%vec_Ac%v(3, i1, i2-1, i3) &
          & +(+2.00d0 * r_inv_h(2)**2 +2.00d0 * r_inv_h(1)**2) * ff%vec_Ac%v(3, i1, i2, i3) &
          & +(-1.00d0 * r_inv_h(2)**2) * ff%vec_Ac%v(3, i1, i2+1, i3) &
          & +(-1.00d0 * r_inv_h(1)**2) * ff%vec_Ac%v(3, i1+1, i2, i3)
        fw%vec_Ac_tmp%v(:,i1, i2, i3) = (2 * ff%vec_Ac%v(:,i1, i2, i3) - ff%vec_Ac_old%v(:,i1, i2, i3) &
          & + ff%vec_j_em%v(:,i1, i2, i3) * 4.0 * pi* (dt_em**2) - rot2_Ac(:)* (cspeed_au * dt_em)**2 )
      end do
    end do
    !$omp end parallel do
    ! Temporal data transfer: Ac_tmp -> Ac -> Ac_old
    call copy_data(ff%vec_Ac%v, ff%vec_Ac_old%v)
    call copy_data(fw%vec_Ac_tmp%v, ff%vec_Ac%v)
    return
  end subroutine dt_evolve_Ac_2d


  subroutine dt_evolve_Ac_3d()
    implicit none
    integer :: i1, i2, i3
    real(8) :: rot2_Ac(3) ! rot rot Ac
    real(8) :: r_inv_h(3)
    r_inv_h(:) = 1.00 / fs%hgs(:)
    !$omp parallel do collapse(3) default(shared) private(i1,i2,i3,rot2_Ac)
    do i3 = fs%mg%is(3), fs%mg%ie(3)
      do i2 = fs%mg%is(2), fs%mg%ie(2)
        do i1 = fs%mg%is(1), fs%mg%ie(1)
          ! Calculate Rot Rot A
          rot2_Ac(1) = - (r_inv_h(2)**2) * ff%vec_Ac%v(1, i1+0, i2-1, i3+0) &
            & - (r_inv_h(3)**2) * ff%vec_Ac%v(1, i1+0, i2+0, i3-1) &
            & + (2d0* (r_inv_h(2)**2 + r_inv_h(3)**2)) * ff%vec_Ac%v(1, i1+0, i2+0, i3+0) &
            & - (r_inv_h(3)**2) * ff%vec_Ac%v(1, i1+0, i2+0, i3+1) &
            & - (r_inv_h(2)**2) * ff%vec_Ac%v(1, i1+0, i2+1, i3+0) &
            & + (r_inv_h(1) * r_inv_h(2) * 0.25d0) * ff%vec_Ac%v(2, i1-1, i2-1, i3+0) &
            & - (r_inv_h(1) * r_inv_h(2) * 0.25d0) * ff%vec_Ac%v(2, i1-1, i2+1, i3+0) &
            & - (r_inv_h(1) * r_inv_h(2) * 0.25d0) * ff%vec_Ac%v(2, i1+1, i2-1, i3+0) &
            & + (r_inv_h(1) * r_inv_h(2) * 0.25d0) * ff%vec_Ac%v(2, i1+1, i2+1, i3+0) &
            & + (r_inv_h(1) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(3, i1-1, i2+0, i3-1) &
            & - (r_inv_h(1) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(3, i1-1, i2+0, i3+1) &
            & - (r_inv_h(1) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(3, i1+1, i2+0, i3-1) &
            & + (r_inv_h(1) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(3, i1+1, i2+0, i3+1)
          rot2_Ac(2) = + (r_inv_h(1) * r_inv_h(2) * 0.25d0) * ff%vec_Ac%v(1, i1-1, i2-1, i3+0) &
            & - (r_inv_h(1) * r_inv_h(2) * 0.25d0) * ff%vec_Ac%v(1, i1-1, i2+1, i3+0) &
            & - (r_inv_h(1) * r_inv_h(2) * 0.25d0) * ff%vec_Ac%v(1, i1+1, i2-1, i3+0) &
            & + (r_inv_h(1) * r_inv_h(2) * 0.25d0) * ff%vec_Ac%v(1, i1+1, i2+1, i3+0) &
            & - (r_inv_h(1)**2) * ff%vec_Ac%v(2, i1-1, i2+0, i3+0) &
            & - (r_inv_h(3)**2) * ff%vec_Ac%v(2, i1+0, i2+0, i3-1) &
            & + (2d0* (r_inv_h(1)**2 + r_inv_h(3)**2)) * ff%vec_Ac%v(2, i1+0, i2+0, i3+0) &
            & - (r_inv_h(3)**2) * ff%vec_Ac%v(2, i1+0, i2+0, i3+1) &
            & - (r_inv_h(1)**2) * ff%vec_Ac%v(2, i1+1, i2+0, i3+0) &
            & + (r_inv_h(2) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(3, i1+0, i2-1, i3-1) &
            & - (r_inv_h(2) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(3, i1+0, i2-1, i3+1) &
            & - (r_inv_h(2) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(3, i1+0, i2+1, i3-1) &
            & + (r_inv_h(2) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(3, i1+0, i2+1, i3+1)
          rot2_Ac(3) = + (r_inv_h(1) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(1, i1-1, i2+0, i3-1) &
            & - (r_inv_h(1) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(1, i1-1, i2+0, i3+1) &
            & - (r_inv_h(1) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(1, i1+1, i2+0, i3-1) &
            & + (r_inv_h(1) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(1, i1+1, i2+0, i3+1) &
            & + (r_inv_h(2) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(2, i1+0, i2-1, i3-1) &
            & - (r_inv_h(2) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(2, i1+0, i2-1, i3+1) &
            & - (r_inv_h(2) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(2, i1+0, i2+1, i3-1) &
            & + (r_inv_h(2) * r_inv_h(3) * 0.25d0) * ff%vec_Ac%v(2, i1+0, i2+1, i3+1) &
            & - (r_inv_h(1)**2) * ff%vec_Ac%v(3, i1-1, i2+0, i3+0) &
            & - (r_inv_h(2)**2) * ff%vec_Ac%v(3, i1+0, i2-1, i3+0) &
            & + (2d0 * (r_inv_h(1)**2 + r_inv_h(2)**2)) * ff%vec_Ac%v(3, i1+0, i2+0, i3+0) &
            & - (r_inv_h(2)**2) * ff%vec_Ac%v(3, i1+0, i2+1, i3+0) &
            & - (r_inv_h(1)**2) * ff%vec_Ac%v(3, i1+1, i2+0, i3+0)
          fw%vec_Ac_tmp%v(:,i1, i2, i3) = &
            & + (2 * ff%vec_Ac%v(:,i1, i2, i3) - ff%vec_Ac_old%v(:,i1, i2, i3) &
            & + ff%vec_j_em%v(:,i1, i2, i3) * 4.0 * pi* (dt_em**2) - rot2_Ac(:) * (cspeed_au * dt_em)**2)
        end do
      end do
    end do
    !$omp end parallel do
    ! Temporal data transfer: Ac_tmp -> Ac -> Ac_old
    call copy_data(ff%vec_Ac%v, ff%vec_Ac_old%v)
    call copy_data(fw%vec_Ac_tmp%v, ff%vec_Ac%v)
    return
  end subroutine dt_evolve_Ac_3d


  subroutine dt_evolve_Ac_3d_cylindal()
    implicit none
    integer :: i1, i2
    integer :: i3
    real(8) :: Y, rot2_Ac(3) ! rot rot Ac
    real(8) :: r_inv_h(3)
    r_inv_h(:) = 1d0 / fs%hgs(:)
    i3 = fs%mg%is(3)
    !$omp parallel do collapse(2) default(shared) private(i1, i2, rot2_Ac,Y)
    do i2 = fs%mg%is(2), fs%mg%ie(2)
      do i1 = fs%mg%is(1), fs%mg%ie(1)
        Y = i2 * fs%hgs(2) + fs%origin(2)
        rot2_Ac(1) = +(+0.50d0* (1.00d0 * r_inv_h(2))* (1.00d0/Y)-(1.00d0 * r_inv_h(2)**2))*ff%vec_Ac%v(1,i1+0,i2-1,i3) &
          & +2.00d0* (1.00d0 * r_inv_h(2)**2)*ff%vec_Ac%v(1,i1+0,i2+0,i3) &
          & +(-0.50d0* (1.00d0 * r_inv_h(2))* (1.00d0/Y)-(1.00d0 * r_inv_h(2)**2))*ff%vec_Ac%v(1,i1+0,i2+1,i3) &
          & +0.25d0* (1.00d0 * r_inv_h(1) * r_inv_h(2))*ff%vec_Ac%v(2,i1-1,i2-1,i3) &
          & -0.50d0* (1.00d0 * r_inv_h(1))* (1.00d0/Y)*ff%vec_Ac%v(2,i1-1,i2+0,i3) &
          & -0.25d0* (1.00d0 * r_inv_h(1) * r_inv_h(2))*ff%vec_Ac%v(2,i1-1,i2+1,i3) &
          & -0.25d0* (1.00d0 * r_inv_h(1) * r_inv_h(2))*ff%vec_Ac%v(2,i1+1,i2-1,i3) &
          & +0.50d0* (1.00d0 * r_inv_h(1))* (1.00d0/Y)*ff%vec_Ac%v(2,i1+1,i2+0,i3) &
          & +0.25d0* (1.00d0 * r_inv_h(1) * r_inv_h(2))*ff%vec_Ac%v(2,i1+1,i2+1,i3)
        rot2_Ac(2) = +0.25d0* (1.00d0 * r_inv_h(1) * r_inv_h(2))*ff%vec_Ac%v(1,i1-1,i2-1,i3) &
          & -0.25d0* (1.00d0 * r_inv_h(1) * r_inv_h(2))*ff%vec_Ac%v(1,i1-1,i2+1,i3) &
          & -0.25d0* (1.00d0 * r_inv_h(1) * r_inv_h(2))*ff%vec_Ac%v(1,i1+1,i2-1,i3) &
          & +0.25d0* (1.00d0 * r_inv_h(1) * r_inv_h(2))*ff%vec_Ac%v(1,i1+1,i2+1,i3) &
          & -(1.00d0 * r_inv_h(1)**2)*ff%vec_Ac%v(2,i1-1,i2+0,i3) &
          & +2.00d0* (1.00d0 * r_inv_h(1)**2)*ff%vec_Ac%v(2,i1+0,i2+0,i3) &
          & -(1.00d0 * r_inv_h(1)**2)*ff%vec_Ac%v(2,i1+1,i2+0,i3)
        rot2_Ac(3) = -(1.00d0 * r_inv_h(1)**2)*ff%vec_Ac%v(3,i1-1,i2+0,i3) &
          & +(+0.50d0* (1.00d0 * r_inv_h(2))* (1.00d0/Y)-(1.00d0 * r_inv_h(2)**2))*ff%vec_Ac%v(3,i1+0,i2-1,i3) &
          & +(+(1.00d0/Y**2)+2.00d0* (1.00d0 * r_inv_h(1)**2)+2.00d0* (1.00d0 * r_inv_h(2)**2))*ff%vec_Ac%v(3,i1+0,i2+0,i3) &
          & +(-0.50d0* (1.00d0 * r_inv_h(2))* (1.00d0/Y)-(1.00d0 * r_inv_h(2)**2))*ff%vec_Ac%v(3,i1+0,i2+1,i3) &
          & -(1.00d0 * r_inv_h(1)**2)*ff%vec_Ac%v(3,i1+1,i2+0,i3)
        fw%vec_Ac_tmp%v(:,i1, i2, i3) = (2 * ff%vec_Ac%v(:,i1, i2, i3) - ff%vec_Ac_old%v(:,i1, i2, i3) &
          & + ff%vec_j_em%v(:,i1, i2, i3) * 4.0 * pi* (dt_em**2) - rot2_Ac(:)* (cspeed_au * dt_em)**2)
      end do
    end do
    !$omp end parallel do
    ! Temporal data transfer: Ac_tmp -> Ac -> Ac_old
    call copy_data(ff%vec_Ac%v, ff%vec_Ac_old%v)
    call copy_data(fw%vec_Ac_tmp%v, ff%vec_Ac%v)
    return
  end subroutine dt_evolve_Ac_3d_cylindal


  subroutine calc_vec_e()
    implicit none
    integer ::  i1, i2, i3
    real(8) :: r_inv_dt
    r_inv_dt = 1d0 / dt_em
    !$omp parallel do collapse(3) default(shared) private(i1, i2, i3)
    do i3 = fs%mg%is(3), fs%mg%ie(3)
      do i2 = fs%mg%is(2), fs%mg%ie(2)
        do i1 = fs%mg%is(1), fs%mg%ie(1)
          ! Calculate: E = - diff(Ac, t)
          ff%vec_e%v(:, i1, i2, i3) = - (ff%vec_Ac%v(:, i1, i2, i3) - ff%vec_Ac_old%v(:, i1, i2, i3)) * r_inv_dt
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine calc_vec_e


  subroutine calc_vec_h_1d()
    implicit none
    integer :: i1, i2, i3
    real(8) :: rot_Ac(3)
    real(8) :: r_inv_h(3)
    i2 = fs%mg%is(2)
    i3 = fs%mg%is(3)
    r_inv_h(1) = 1d0 / fs%hgs(1)
    !$omp parallel do default(shared) private(i1, rot_Ac)
    do i1 = fs%mg%is(1), fs%mg%ie(1)
      rot_Ac(1) = 0.0d0
      rot_Ac(2) = - (ff%vec_Ac%v(3, i1+1, i2, i3) - ff%vec_Ac%v(3, i1-1, i2, i3)) * (0.5d0 * r_inv_h(1))
      rot_Ac(3) = + (ff%vec_Ac%v(2, i1+1, i2, i3) - ff%vec_Ac%v(2, i1-1, i2, i3)) * (0.5d0 * r_inv_h(1))
      ff%vec_h%v(:, i1, i2, i3) = rot_Ac(:) * cspeed_au
    end do
    !$omp end parallel do
    return
  end subroutine calc_vec_h_1d


  subroutine calc_vec_h_2d()
    implicit none
    integer :: i1, i2, i3
    real(8) :: rot_Ac(3)
    real(8) :: r_inv_h(3)
    i3 = fs%mg%is(3)
    r_inv_h(1:2) = 1d0 / fs%hgs(1:2)
    !$omp parallel do collapse(2) default(shared) private(i1, i2, rot_Ac)
    do i2 = fs%mg%is(2), fs%mg%ie(2)
      do i1 = fs%mg%is(1), fs%mg%ie(1)
        rot_Ac(1) = +(+ff%vec_Ac%v(3, i1+0, i2+1, i3+0)-ff%vec_Ac%v(3, i1+0, i2-1, i3+0)) * (0.5d0 * r_inv_h(2))
        rot_Ac(2) = -(+ff%vec_Ac%v(3, i1+1, i2+0, i3+0)-ff%vec_Ac%v(3, i1-1, i2+0, i3+0)) * (0.5d0 * r_inv_h(1))
        rot_Ac(3) = +(+ff%vec_Ac%v(2, i1+1, i2+0, i3+0)-ff%vec_Ac%v(2, i1-1, i2+0, i3+0)) * (0.5d0 * r_inv_h(1)) &
               & -(+ff%vec_Ac%v(1, i1+0, i2+1, i3+0)-ff%vec_Ac%v(1, i1+0, i2-1, i3+0)) * (0.5d0 * r_inv_h(2))
        ff%vec_h%v(:, i1, i2, i3) = rot_Ac(:) * cspeed_au
      end do
    end do
    !$omp end parallel do
    return
  end subroutine calc_vec_h_2d


  subroutine calc_vec_h_3d()
    implicit none
    integer :: i1, i2, i3
    real(8) :: rot_Ac(3)
    real(8) :: r_inv_h(3)
    r_inv_h(:) = 1d0 / fs%hgs(:)
    !$omp parallel do collapse(3) default(shared) private(i1, i2, i3, rot_Ac)
    do i3 = fs%mg%is(3), fs%mg%ie(3)
      do i2 = fs%mg%is(2), fs%mg%ie(2)
        do i1 = fs%mg%is(1), fs%mg%ie(1)
          rot_Ac(1) = +(+ff%vec_Ac%v(3, i1+0, i2+1, i3+0)-ff%vec_Ac%v(3, i1+0, i2-1, i3+0)) * (0.5d0 * r_inv_h(2)) &
            & -(+ff%vec_Ac%v(2, i1+0, i2+0, i3+1)-ff%vec_Ac%v(2, i1+0, i2+0, i3-1)) * (0.5d0 * r_inv_h(3))
          rot_Ac(2) = +(+ff%vec_Ac%v(1, i1+0, i2+0, i3+1)-ff%vec_Ac%v(1, i1+0, i2+0, i3-1)) * (0.5d0 * r_inv_h(3)) &
            & -(+ff%vec_Ac%v(3, i1+1, i2+0, i3+0)-ff%vec_Ac%v(3, i1-1, i2+0, i3+0)) * (0.5d0 * r_inv_h(1))
          rot_Ac(3) = +(+ff%vec_Ac%v(2, i1+1, i2+0, i3+0)-ff%vec_Ac%v(2, i1-1, i2+0, i3+0)) * (0.5d0 * r_inv_h(1)) &
            & -(+ff%vec_Ac%v(1, i1+0, i2+1, i3+0)-ff%vec_Ac%v(1, i1+0, i2-1, i3+0)) * (0.5d0 * r_inv_h(2))
          ff%vec_h%v(:, i1, i2, i3) = rot_Ac * cspeed_au
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine calc_vec_h_3d


  subroutine calc_vec_h_3d_cylindal()
    implicit none
    integer :: i1, i2, i3
    real(8) :: Y, rot_Ac(3)
    real(8) :: r_inv_h(3)
    i3 = fs%mg%is(3)
    r_inv_h(1:2) = 1d0 / fs%hgs(1:2)
    !$omp parallel do collapse(2) default(shared) private(i1, i2, Y, rot_Ac)
    do i2 = fs%mg%is(2), fs%mg%ie(2)
      do i1 = fs%mg%is(1), fs%mg%ie(1)
        Y = i2 * fs%hgs(2) + fs%origin(2)
        rot_Ac(1) = -0.50d0*(r_inv_h(2))*ff%vec_Ac%v(3,i1+0,i2-1,i3) &
          & +1.00d0*(1.00d0/Y)*ff%vec_Ac%v(3,i1+0,i2+0,i3) &
          & +0.50d0*(r_inv_h(2))*ff%vec_Ac%v(3,i1+0,i2+1,i3)
        rot_Ac(2) = +0.50d0*(r_inv_h(1))*ff%vec_Ac%v(3,i1-1,i2+0,i3) &
          & -0.50d0*(r_inv_h(1))*ff%vec_Ac%v(3,i1+1,i2+0,i3)
        rot_Ac(3) = +0.50d0*(r_inv_h(2))*ff%vec_Ac%v(1,i1+0,i2-1,i3) &
          & -0.50d0*(r_inv_h(2))*ff%vec_Ac%v(1,i1+0,i2+1,i3) &
          & -0.50d0*(r_inv_h(1))*ff%vec_Ac%v(2,i1-1,i2+0,i3) &
          & +0.50d0*(r_inv_h(1))*ff%vec_Ac%v(2,i1+1,i2+0,i3)
        ff%vec_h%v(:, i1, i2, i3) = rot_Ac(:) * cspeed_au
      end do
    end do
    !$omp end parallel do
    return
  end subroutine calc_vec_h_3d_cylindal


  subroutine calc_energy_density()
    implicit none
    integer :: i1, i2, i3
    !$omp parallel do collapse(3) default(shared) private(i1,i2,i3)
    do i3 = fs%mg%is(3), fs%mg%ie(3)
      do i2 = fs%mg%is(2), fs%mg%ie(2)
        do i1 = fs%mg%is(1), fs%mg%ie(1)
          fw%edensity_emfield%f(i1, i2, i3) = (1d0 / (8d0 * pi)) * ( &
            & norm2(ff%vec_e%v(:, i1, i2, i3)) + norm2(ff%vec_h%v(:, i1, i2, i3)) &
            & )
          fw%edensity_absorb%f(i1, i2, i3) = fw%edensity_absorb(i1, i2, i3) + ( &
            & dot_product(ff%vec_e%v(:, i1, i2, i3), ff%vec_j_em%v(:, i1, i2, i3)) &
            & ) * dt_em
        end do
      end do
    end do
    return
  end subroutine calc_energy_density
    

end subroutine weyl_calc
