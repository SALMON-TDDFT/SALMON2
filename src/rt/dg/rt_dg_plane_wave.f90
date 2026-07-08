!
!  Copyright 2019-2026 SALMON developers
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
!=======================================================================
! Plane-wave mixing utilities for DG-Fragment RT-TDDFT
!=======================================================================

module rt_dg_plane_wave
  use salmon_global, only: dg_mixed_z_local_prop_backend, dg_mixed_z_neighbor_env_shell
  use rt_dg_fragment_types, only: s_dg_fragment_rt, complex_matrix_block_info
  use rt_dg_wpw_reduced_mixedz_diag, only: wpw_reduced_mixedz_operator_stats, &
    wpw_reduced_canon_mixedz_current_coeff_stats, wpw_reduced_canon_mixedz_bridge_stats, &
    wpw_reduced_canon_pz_block_operator_stats, wpw_reduced_canon_p_projection_stats, &
    wpw_reduced_prod_p_basis_save_stats
  use rt_dg_wpw_window, only: wpw_window_buffer_axis, wpw_window_transition_width_axis, &
    wpw_fragment_box_size, wpw_fragment_box_bounds, wpw_raw_window_at_grid, &
    wpw_normalized_window_at_grid, wpw_face_neighbor_fragment, wpw_local_is_neighbor_pair, &
    map_global_to_phi_box_coord_pw
  use rt_dg_wpw_linalg, only: build_wpw_sorth_reduced_neighbor_block, &
    apply_s_orthogonal_transform, build_symmetric_orthogonal_hamiltonian, &
    build_hermitian_pseudoinverse, build_hermitian_inverse_sqrt, build_hermitian_inverse, &
    build_wpw_reduced_c_can_reference, build_wpw_reduced_raw_back_hybrid_from_inverse, &
    build_wpw_reduced_raw_back_hybrid_matrix, zheev_with_query, zhegv_with_query, &
    rectangular_singular_minmax, wpw_local_herm_max
  use structures, only: s_dft_system, s_scalar
  implicit none

  private
  public :: init_plane_wave_basis
  public :: compute_fragment_pw_overlap
  public :: compute_fragment_pw_gradient_from_overlap
  public :: compute_fragment_pw_position_overlap
  public :: compute_mean_potential
  public :: compute_fragment_pw_hamiltonian
  public :: build_mixed_hamiltonian
  public :: build_local_fragment_pw_row_list
  public :: prepare_local_fragment_pw_blocks
  public :: diagonalize_mixed_basis
  public :: assemble_mixed_hamiltonian_dense
  public :: prepare_mixed_basis_startup
  public :: initialize_wpw_windows
  public :: ensure_wpw_local_pp_blocks
  public :: build_wpw_sorth_reduced_neighbor_block
  public :: initialize_wpw_reduced_self_projection
  public :: diagnose_wpw_reduced_embed_local
  public :: diagnose_wpw_reduced_embed_prodcoef
  public :: diagnose_wpw_reduced_density
  public :: diagnose_wpw_mixed_neighbor_hamiltonian
  public :: apply_wpw_reduced_density_to_production
  public :: apply_wpw_reduced_pz_to_production
  public :: wpw_normalized_window_at_grid
  public :: wpw_face_neighbor_fragment
  public :: map_global_to_phi_box_coord_pw
  public :: reconstruct_psi_from_C_can
  public :: build_wpw_reduced_c_can_reference
  public :: build_wpw_reduced_raw_back_hybrid_from_inverse
  public :: build_wpw_reduced_raw_back_hybrid_matrix

contains

  logical function plane_wave_trace_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_PW_TRACE', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function plane_wave_trace_enabled

  logical function wpw_trace_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_TRACE', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_trace_enabled

  logical function wpw_kinetic_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_KINETIC_DIAG', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_kinetic_diag_enabled

  logical function wpw_kinetic_use_weak() result(enabled)
    implicit none
    character(len=32) :: env_mode
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_mode = ''
      call get_environment_variable('SALMON_DG_WPW_KINETIC', env_mode, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_mode(1:env_len))))
        case ('weak', 'WEAK', 'Weak', 'wpw', 'WPW')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_kinetic_use_weak

  logical function wpw_kinetic_unit_test_enabled() result(enabled)
    implicit none
    character(len=32) :: env_mode
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_mode = ''
      call get_environment_variable('SALMON_DG_WPW_UNIT_WINDOW_TEST', env_mode, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_mode(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_kinetic_unit_test_enabled

  logical function wpw_local_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_LOCAL_DIAG', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_local_diag_enabled

  logical function wpw_interface_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_INTERFACE_DIAG', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_interface_diag_enabled

  logical function wpw_pp_blocks_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_PP_BLOCKS', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_pp_blocks_enabled

  logical function wpw_pp_prop_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_PP_PROP_DIAG', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_pp_prop_diag_enabled

  logical function wpw_mixed_block_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_MIXED_BLOCK_DIAG', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_mixed_block_diag_enabled

  logical function wpw_mixed_h_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_MIXED_H_DIAG', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_mixed_h_diag_enabled

  logical function wpw_mixed_neighbor_h_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_MIXED_NEIGHBOR_H_DIAG', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_mixed_neighbor_h_diag_enabled

  logical function wpw_neighbor_sorth_candidate_enabled() result(enabled)
    implicit none
    character(len=32) :: env_trace
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env_trace = ''
      call get_environment_variable('SALMON_DG_WPW_NEIGHBOR_SORTH', env_trace, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_trace(1:env_len))))
        case ('1', 'y', 'Y', 'yes', 'YES', 'true', 'TRUE', 'on', 'ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_neighbor_sorth_candidate_enabled

  real(8) function wpw_neighbor_sorth_tol() result(tol)
    implicit none
    character(len=32) :: env_value
    integer :: env_len, env_stat
    real(8) :: parsed
    logical, save :: initialized = .false.
    real(8), save :: cached_tol = 1.0d-8

    if (.not. initialized) then
      env_value = ''
      call get_environment_variable('SALMON_DG_WPW_NEIGHBOR_SORTH_TOL', env_value, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=env_stat) parsed
        if (env_stat == 0 .and. parsed >= 0.0d0) cached_tol = parsed
      end if
      initialized = .true.
    end if
    tol = cached_tol
  end function wpw_neighbor_sorth_tol

  real(8) function wpw_interface_penalty_factor() result(factor)
    implicit none
    character(len=64) :: env_value
    integer :: env_len, env_stat, ios
    real(8), save :: cached_factor = 0.3d0
    logical, save :: initialized = .false.

    if (.not. initialized) then
      env_value = ''
      call get_environment_variable('SALMON_DG_WPW_INTERFACE_PENALTY', env_value, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=ios) cached_factor
        if (ios /= 0) cached_factor = 0.3d0
      end if
      cached_factor = max(0.0d0, cached_factor)
      initialized = .true.
    end if
    factor = cached_factor
  end function wpw_interface_penalty_factor



  integer function wpw_neighbor_filter_drop_g2_override() result(drop_g2)
    implicit none
    character(len=64) :: env_value
    integer :: env_len, env_stat, ios
    integer, save :: cached_drop_g2 = huge(1)
    logical, save :: initialized = .false.

    if (.not. initialized) then
      env_value = ''
      call get_environment_variable('SALMON_DG_WPW_NEIGHBOR_FILTER_DROP_G2', env_value, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=ios) cached_drop_g2
        if (ios /= 0) cached_drop_g2 = huge(1)
      end if
      initialized = .true.
    end if
    drop_g2 = cached_drop_g2
  end function wpw_neighbor_filter_drop_g2_override

  integer function wpw_neighbor_env_shell_from_backend() result(shell)
    shell = max(1, dg_mixed_z_neighbor_env_shell)
    select case (trim(adjustl(dg_mixed_z_local_prop_backend)))
    case ('neighbor_env_delta_shell2','neighbor_delta_shell2','delta_shell2')
      shell = max(shell, 2)
    end select
  end function wpw_neighbor_env_shell_from_backend

  subroutine collect_wpw_fragment_shell_local(dg_frag, ifrag_center, shell_max, pfrag_ids, pfrag_axis, pfrag_side, n_pfrag)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_center, shell_max
    integer, allocatable, intent(out) :: pfrag_ids(:), pfrag_axis(:), pfrag_side(:)
    integer, intent(out) :: n_pfrag
    integer :: max_frag, axis, side, level, start_idx, end_idx, idx, jfrag, base_frag
    integer :: radius, dx, dy, dz
    integer, allocatable :: queue(:), queue_axis(:), queue_side(:)

    max_frag = max(1, product(dg_frag%num_fragment(1:3)))
    allocate(queue(max_frag), queue_axis(max_frag), queue_side(max_frag))
    queue(:) = 0
    queue_axis(:) = 0
    queue_side(:) = 0
    n_pfrag = 1
    queue(1) = ifrag_center
    if (shell_max >= 2) then
      radius = max(1, shell_max - 1)
      do dx = -radius, radius
        do dy = -radius, radius
          do dz = -radius, radius
            jfrag = wpw_periodic_offset_fragment_id(dg_frag, ifrag_center, dx, dy, dz)
            if (jfrag <= 0) cycle
            if (.not. any(queue(1:n_pfrag) == jfrag) .and. n_pfrag < max_frag) then
              n_pfrag = n_pfrag + 1
              queue(n_pfrag) = jfrag
            end if
          end do
        end do
      end do
      allocate(pfrag_ids(n_pfrag), pfrag_axis(n_pfrag), pfrag_side(n_pfrag))
      pfrag_ids(1:n_pfrag) = queue(1:n_pfrag)
      pfrag_axis(1:n_pfrag) = queue_axis(1:n_pfrag)
      pfrag_side(1:n_pfrag) = queue_side(1:n_pfrag)
      deallocate(queue, queue_axis, queue_side)
      return
    end if
    start_idx = 1
    end_idx = 1
    do level = 1, max(0, shell_max)
      idx = start_idx
      start_idx = end_idx + 1
      do while (idx <= end_idx)
        base_frag = queue(idx)
        do axis = 1, 3
          do side = -1, 1, 2
            jfrag = wpw_periodic_face_neighbor_id(dg_frag, base_frag, axis, side)
            if (jfrag <= 0) cycle
            if (.not. any(queue(1:n_pfrag) == jfrag) .and. n_pfrag < max_frag) then
              n_pfrag = n_pfrag + 1
              queue(n_pfrag) = jfrag
              if (base_frag == ifrag_center) then
                queue_axis(n_pfrag) = axis
                queue_side(n_pfrag) = side
              end if
            end if
          end do
        end do
        idx = idx + 1
      end do
      end_idx = n_pfrag
    end do
    allocate(pfrag_ids(n_pfrag), pfrag_axis(n_pfrag), pfrag_side(n_pfrag))
    pfrag_ids(1:n_pfrag) = queue(1:n_pfrag)
    pfrag_axis(1:n_pfrag) = queue_axis(1:n_pfrag)
    pfrag_side(1:n_pfrag) = queue_side(1:n_pfrag)
    deallocate(queue, queue_axis, queue_side)
  end subroutine collect_wpw_fragment_shell_local

  integer function wpw_periodic_face_neighbor_id(dg_frag, ifrag, axis, side) result(jfrag)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, side
    integer :: coords(3), rem, nfrag_total

    jfrag = 0
    nfrag_total = max(1, product(dg_frag%num_fragment(1:3)))
    if (ifrag < 1 .or. ifrag > nfrag_total) return
    if (axis < 1 .or. axis > 3) return
    coords(1) = (ifrag - 1) / max(1, dg_frag%num_fragment(2) * dg_frag%num_fragment(3)) + 1
    rem = modulo(ifrag - 1, max(1, dg_frag%num_fragment(2) * dg_frag%num_fragment(3)))
    coords(2) = rem / max(1, dg_frag%num_fragment(3)) + 1
    coords(3) = modulo(rem, max(1, dg_frag%num_fragment(3))) + 1
    coords(axis) = modulo(coords(axis) - 1 + side + max(1, dg_frag%num_fragment(axis)), &
      max(1, dg_frag%num_fragment(axis))) + 1
    jfrag = ((coords(1) - 1) * max(1, dg_frag%num_fragment(2)) + coords(2) - 1) * &
      max(1, dg_frag%num_fragment(3)) + coords(3)
  end function wpw_periodic_face_neighbor_id

  integer function wpw_periodic_offset_fragment_id(dg_frag, ifrag, dx, dy, dz) result(jfrag)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, dx, dy, dz
    integer :: coords(3), rem, nfrag_total

    jfrag = 0
    nfrag_total = max(1, product(dg_frag%num_fragment(1:3)))
    if (ifrag < 1 .or. ifrag > nfrag_total) return
    coords(1) = (ifrag - 1) / max(1, dg_frag%num_fragment(2) * dg_frag%num_fragment(3)) + 1
    rem = modulo(ifrag - 1, max(1, dg_frag%num_fragment(2) * dg_frag%num_fragment(3)))
    coords(2) = rem / max(1, dg_frag%num_fragment(3)) + 1
    coords(3) = modulo(rem, max(1, dg_frag%num_fragment(3))) + 1
    coords(1) = modulo(coords(1) - 1 + dx + max(1, dg_frag%num_fragment(1)), &
      max(1, dg_frag%num_fragment(1))) + 1
    coords(2) = modulo(coords(2) - 1 + dy + max(1, dg_frag%num_fragment(2)), &
      max(1, dg_frag%num_fragment(2))) + 1
    coords(3) = modulo(coords(3) - 1 + dz + max(1, dg_frag%num_fragment(3)), &
      max(1, dg_frag%num_fragment(3))) + 1
    jfrag = ((coords(1) - 1) * max(1, dg_frag%num_fragment(2)) + coords(2) - 1) * &
      max(1, dg_frag%num_fragment(3)) + coords(3)
  end function wpw_periodic_offset_fragment_id

  integer function fragment_local_row_index_pw(dg_frag, ig, ispin, nrow) result(iloc)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ig, ispin, nrow

    iloc = 0
    if (ig < 1 .or. ig > dg_frag%n_mat_max) return
    if (allocated(dg_frag%coef_global_to_local)) then
      if (ispin < 1 .or. ispin > size(dg_frag%coef_global_to_local, 2)) return
      iloc = dg_frag%coef_global_to_local(ig, ispin)
    else
      iloc = ig
    end if
    if (iloc < 1 .or. iloc > nrow) iloc = 0
  end function fragment_local_row_index_pw

  logical function cubic_k_order_less(a, b) result(is_less)
    implicit none
    integer, intent(in) :: a(3), b(3)
    integer :: ia, ib

    is_less = .false.
    ia = abs(a(1)) + abs(a(2)) + abs(a(3))
    ib = abs(b(1)) + abs(b(2)) + abs(b(3))
    if (ia /= ib) then
      is_less = ia < ib
    else if (abs(a(3)) /= abs(b(3))) then
      is_less = abs(a(3)) < abs(b(3))
    else if (abs(a(2)) /= abs(b(2))) then
      is_less = abs(a(2)) < abs(b(2))
    else if (abs(a(1)) /= abs(b(1))) then
      is_less = abs(a(1)) < abs(b(1))
    else if (a(3) /= b(3)) then
      is_less = a(3) < b(3)
    else if (a(2) /= b(2)) then
      is_less = a(2) < b(2)
    else if (a(1) /= b(1)) then
      is_less = a(1) < b(1)
    end if
  end function cubic_k_order_less

  subroutine assemble_mixed_hamiltonian_dense(dg_frag, ispin, H_frag_pw, mat)
    use rt_dg_fragment_ops, only: copy_matrix_blocks_to_complex_dense
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: H_frag_pw(:,:,:)
    complex(8), intent(inout) :: mat(:,:)

    integer :: n_frag, n_pw, ipw, i

    n_frag = dg_frag%n_mat_max
    n_pw = dg_frag%n_plane_waves
    mat(:, :) = (0.0d0, 0.0d0)

    if (allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)) then
      mat(1:n_frag, 1:n_frag) = dg_frag%H_mat_c(1:n_frag, 1:n_frag, ispin)
    else if (allocated(dg_frag%H_mat_blocks)) then
      call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, mat(1:n_frag, 1:n_frag))
    else if (allocated(dg_frag%H_mat)) then
      mat(1:n_frag, 1:n_frag) = cmplx(dg_frag%H_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
    end if

    if (n_pw <= 0) return

    if (allocated(dg_frag%H_mat_pw)) then
      mat(n_frag+1:n_frag+n_pw, n_frag+1:n_frag+n_pw) = dg_frag%H_mat_pw(1:n_pw, 1:n_pw, ispin)
    else
      do ipw = 1, n_pw
        i = n_frag + ipw
        if (allocated(dg_frag%H_mat_pw_diag)) mat(i, i) = dg_frag%H_mat_pw_diag(ipw, ispin)
      end do
    end if
    mat(1:n_frag, n_frag+1:n_frag+n_pw) = H_frag_pw(1:n_frag, 1:n_pw, ispin)
    mat(n_frag+1:n_frag+n_pw, 1:n_frag) = conjg(transpose(H_frag_pw(1:n_frag, 1:n_pw, ispin)))
  end subroutine assemble_mixed_hamiltonian_dense

  !=======================================================================
  ! Initialize plane wave basis for mixing with fragment basis
  !=======================================================================
  subroutine init_plane_wave_basis(dg_frag, system, lg, info)
    use structures
    use salmon_global, only: yn_plane_wave_basis, n_plane_waves_dg, k_cutoff_plane_wave
    use inputoutput, only: t_unit_energy
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: lg
    type(s_parallel_info),  intent(in)    :: info

    integer :: ikx, iky, ikz, ipw, nk(3)
    real(8) :: Lbox(3), k_vec(3), k_norm
    real(8) :: energy_cutoff_pw
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    integer, allocatable :: k_indices(:,:)
    integer :: n_pw_candidate, n_selected, n_shell_selected

    ! Check if plane wave basis is enabled
    if (yn_plane_wave_basis /= 'y') then
      dg_frag%use_plane_wave_basis = .false.
      dg_frag%n_plane_waves = 0
      return
    end if

    dg_frag%use_plane_wave_basis = .true.
    energy_cutoff_pw = max(0.0d0, k_cutoff_plane_wave)
    dg_frag%k_cutoff_pw = sqrt(2.0d0 * energy_cutoff_pw)

    ! Calculate box size from grid
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))

    if (comm_is_root(info%id_rko)) then
      write(*,*)
      write(*,*) "=== Initializing plane wave basis for DG-Fragment ==="
      write(*,'(1x,a,3f12.6)') "Box size [a.u.]: ", Lbox(1:3)
      write(*,'(1x,a,f10.4,1x,a)') "PW cutoff energy [", energy_cutoff_pw * t_unit_energy%conv, trim(t_unit_energy%name) // "]"
      write(*,'(1x,a,f10.4,a)') "k cutoff [a.u.^-1]: ", dg_frag%k_cutoff_pw, ""
    end if

    ! Estimate number of k-points within cutoff sphere
    nk(1:3) = ceiling(dg_frag%k_cutoff_pw * Lbox(1:3) / (2.0d0*pi)) + 1
    n_pw_candidate = (2*nk(1)+1) * (2*nk(2)+1) * (2*nk(3)+1)

    ! Allocate temporary array for k-point selection
    allocate(k_indices(3, n_pw_candidate))
    n_selected = 0

    ! Select k-points within cutoff sphere
    do ikz = -nk(3), nk(3)
      do iky = -nk(2), nk(2)
        do ikx = -nk(1), nk(1)
          if (ikx == 0 .and. iky == 0 .and. ikz == 0) cycle

          k_vec(1) = 2.0d0*pi/Lbox(1) * dble(ikx)
          k_vec(2) = 2.0d0*pi/Lbox(2) * dble(iky)
          k_vec(3) = 2.0d0*pi/Lbox(3) * dble(ikz)
          k_norm = sqrt(sum(k_vec**2))

          if (k_norm <= dg_frag%k_cutoff_pw) then
            n_selected = n_selected + 1
            if (n_selected <= n_pw_candidate) then
              k_indices(1:3, n_selected) = [ikx, iky, ikz]
            end if
          end if
        end do
      end do
    end do

    ! Sort k-points by cubic shell.  A shell is all integer G vectors
    ! with the same ikx^2+iky^2+ikz^2; never cut a shell halfway.
    block
      integer, allocatable :: shell_keys(:)
      integer :: ii, jj, itemp(3)
      integer :: itemp_key
      integer :: shell_start, shell_end, shell_count
      integer :: prev_total, next_total, target

      allocate(shell_keys(n_selected))

      do ipw = 1, n_selected
        ikx = k_indices(1, ipw)
        iky = k_indices(2, ipw)
        ikz = k_indices(3, ipw)
        shell_keys(ipw) = ikx*ikx + iky*iky + ikz*ikz
      end do

      do ii = 2, n_selected
        do jj = ii, 2, -1
          if (shell_keys(jj) < shell_keys(jj-1) .or. &
              (shell_keys(jj) == shell_keys(jj-1) .and. cubic_k_order_less(k_indices(:,jj), k_indices(:,jj-1)))) then
            itemp_key = shell_keys(jj)
            shell_keys(jj) = shell_keys(jj-1)
            shell_keys(jj-1) = itemp_key
            itemp = k_indices(:, jj)
            k_indices(:, jj) = k_indices(:, jj-1)
            k_indices(:, jj-1) = itemp
          else
            exit
          end if
        end do
      end do

      target = max(0, n_plane_waves_dg)
      n_shell_selected = 0
      shell_start = 1
      do while (shell_start <= n_selected)
        shell_end = shell_start
        do while (shell_end < n_selected .and. shell_keys(shell_end+1) == shell_keys(shell_start))
          shell_end = shell_end + 1
        end do
        shell_count = shell_end - shell_start + 1
        prev_total = n_shell_selected
        next_total = shell_end
        if (target <= 0) then
          n_shell_selected = 0
          exit
        end if
        if (next_total <= target) then
          n_shell_selected = next_total
        else
          if (prev_total <= 0) then
            n_shell_selected = next_total
          else if (abs(next_total - target) < abs(prev_total - target)) then
            n_shell_selected = next_total
          end if
          exit
        end if
        shell_start = shell_end + 1
      end do

      deallocate(shell_keys)
    end block

    dg_frag%n_plane_waves = min(n_selected, n_shell_selected)

    if (comm_is_root(info%id_rko)) then
      write(*,'(1x,a,i0)') "k-points within cutoff: ", n_selected
      write(*,'(1x,a,i0)') "Requested plane waves: ", n_plane_waves_dg
      write(*,'(1x,a,i0,a)') "Using plane waves: ", dg_frag%n_plane_waves, " (complete cubic shells)"
    end if

    allocate(dg_frag%k_pw(3, dg_frag%n_plane_waves))
    allocate(dg_frag%coef_pw(dg_frag%n_plane_waves, dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%coef_pw = (0.0d0, 0.0d0)

    do ipw = 1, dg_frag%n_plane_waves
      ikx = k_indices(1, ipw)
      iky = k_indices(2, ipw)
      ikz = k_indices(3, ipw)
      dg_frag%k_pw(1, ipw) = 2.0d0*pi/Lbox(1) * dble(ikx)
      dg_frag%k_pw(2, ipw) = 2.0d0*pi/Lbox(2) * dble(iky)
      dg_frag%k_pw(3, ipw) = 2.0d0*pi/Lbox(3) * dble(ikz)
    end do

    if (comm_is_root(info%id_rko)) then
      if (dg_frag%n_plane_waves > 0) then
        write(*,'(1x,a,f10.4,1x,a)') "Lowest PW energy: ", &
             0.5d0*sum(dg_frag%k_pw(:,1)**2) * t_unit_energy%conv, trim(t_unit_energy%name)
        if (dg_frag%n_plane_waves > 1) then
          write(*,'(1x,a,f10.4,1x,a)') "Highest PW energy: ", &
               0.5d0*sum(dg_frag%k_pw(:,dg_frag%n_plane_waves)**2) * t_unit_energy%conv, trim(t_unit_energy%name)
        end if
      else
        write(*,'(1x,a)') "No PW mode selected within cutoff/limit."
      end if
    end if

    deallocate(k_indices)

    if (comm_is_root(info%id_rko)) then
      write(*,*) "Plane wave basis initialization complete"
      write(*,*)
    end if

  end subroutine init_plane_wave_basis

  !=======================================================================
  ! Initialize Windowed Plane Wave amplitude windows on fragment+buffer grids.
  ! WPW-1 only stores chi_A and grad(chi_A); kinetic terms are added later
  ! through the weak form 1/2 int grad(phi_i)^* . grad(phi_j).
  !=======================================================================
  subroutine initialize_wpw_windows(dg_frag, info)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_parallel_info),  intent(in)    :: info

    integer :: nlocal, ifrag, ilocal, jfrag
    integer :: nx_max, ny_max, nz_max
    integer :: nx, ny, nz, ix, iy, iz
    integer :: gx, gy, gz
    integer :: nfrag
    real(8), allocatable :: q(:), gq(:,:)
    real(8) :: qsum, sqrt_qsum, qg_sum(3)
    real(8) :: chi2_sum, chi2_min, chi2_max, chi2_maxdev
    real(8) :: raw_chi2_sum, raw_chi2_min, raw_chi2_max, raw_chi2_maxdev
    real(8) :: q_local, gq_local(3)
    real(8), parameter :: tiny_q = 1.0d-28

    dg_frag%has_wpw_window = .false.
    dg_frag%wpw_partition_sum_chi2_min = 0.0d0
    dg_frag%wpw_partition_sum_chi2_max = 0.0d0
    dg_frag%wpw_partition_sum_chi2_maxdev = 0.0d0

    if (.not. dg_frag%use_plane_wave_basis .and. .not. wpw_trace_enabled()) return
    if (.not. allocated(dg_frag%nxyz_domain)) return
    if (.not. allocated(dg_frag%ixyz_frag)) return
    if (dg_frag%ifrag_end < dg_frag%ifrag_start) return

    nlocal = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    nfrag = dg_frag%n_frag
    nx_max = 1
    ny_max = 1
    nz_max = 1
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      call wpw_fragment_box_size(dg_frag, ifrag, nx, ny, nz)
      nx_max = max(nx_max, nx)
      ny_max = max(ny_max, ny)
      nz_max = max(nz_max, nz)
    end do

    if (allocated(dg_frag%wpw_window_box_lo)) deallocate(dg_frag%wpw_window_box_lo)
    if (allocated(dg_frag%wpw_window_box_hi)) deallocate(dg_frag%wpw_window_box_hi)
    if (allocated(dg_frag%wpw_chi)) deallocate(dg_frag%wpw_chi)
    if (allocated(dg_frag%wpw_grad_chi)) deallocate(dg_frag%wpw_grad_chi)

    allocate(dg_frag%wpw_window_box_lo(3, nlocal), dg_frag%wpw_window_box_hi(3, nlocal))
    allocate(dg_frag%wpw_chi(nx_max, ny_max, nz_max, nlocal))
    allocate(dg_frag%wpw_grad_chi(3, nx_max, ny_max, nz_max, nlocal))
    dg_frag%wpw_chi = 0.0d0
    dg_frag%wpw_grad_chi = 0.0d0

    allocate(q(nfrag), gq(3, nfrag))
    chi2_min = huge(1.0d0)
    chi2_max = -huge(1.0d0)
    chi2_maxdev = 0.0d0
    raw_chi2_min = huge(1.0d0)
    raw_chi2_max = -huge(1.0d0)
    raw_chi2_maxdev = 0.0d0

    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      ilocal = ifrag - dg_frag%ifrag_start + 1
      call wpw_fragment_box_bounds(dg_frag, ifrag, dg_frag%wpw_window_box_lo(:, ilocal), &
        dg_frag%wpw_window_box_hi(:, ilocal))
      nx = dg_frag%wpw_window_box_hi(1, ilocal) - dg_frag%wpw_window_box_lo(1, ilocal) + 1
      ny = dg_frag%wpw_window_box_hi(2, ilocal) - dg_frag%wpw_window_box_lo(2, ilocal) + 1
      nz = dg_frag%wpw_window_box_hi(3, ilocal) - dg_frag%wpw_window_box_lo(3, ilocal) + 1

      do iz = 1, nz
        gz = dg_frag%wpw_window_box_lo(3, ilocal) + iz - 1
        do iy = 1, ny
          gy = dg_frag%wpw_window_box_lo(2, ilocal) + iy - 1
          do ix = 1, nx
            gx = dg_frag%wpw_window_box_lo(1, ilocal) + ix - 1
            q = 0.0d0
            gq = 0.0d0
            do jfrag = 1, nfrag
              call wpw_raw_window_at_grid(dg_frag, jfrag, gx, gy, gz, q(jfrag), gq(:, jfrag))
            end do
            qsum = sum(q(1:nfrag)**2)
            raw_chi2_sum = qsum
            raw_chi2_min = min(raw_chi2_min, raw_chi2_sum)
            raw_chi2_max = max(raw_chi2_max, raw_chi2_sum)
            raw_chi2_maxdev = max(raw_chi2_maxdev, abs(raw_chi2_sum - 1.0d0))
            if (qsum <= tiny_q) cycle

            sqrt_qsum = sqrt(qsum)
            qg_sum(1:3) = 0.0d0
            do jfrag = 1, nfrag
              qg_sum(1:3) = qg_sum(1:3) + q(jfrag) * gq(1:3, jfrag)
            end do
            q_local = q(ifrag)
            gq_local(1:3) = gq(1:3, ifrag)
            dg_frag%wpw_chi(ix, iy, iz, ilocal) = q_local / sqrt_qsum
            dg_frag%wpw_grad_chi(1:3, ix, iy, iz, ilocal) = &
              gq_local(1:3) / sqrt_qsum - q_local * qg_sum(1:3) / (qsum * sqrt_qsum)

            chi2_sum = sum((q(1:nfrag) / sqrt_qsum)**2)
            chi2_min = min(chi2_min, chi2_sum)
            chi2_max = max(chi2_max, chi2_sum)
            chi2_maxdev = max(chi2_maxdev, abs(chi2_sum - 1.0d0))
          end do
        end do
      end do
    end do

    deallocate(q, gq)
    if (chi2_min == huge(1.0d0)) chi2_min = 0.0d0
    if (chi2_max == -huge(1.0d0)) chi2_max = 0.0d0
    if (raw_chi2_min == huge(1.0d0)) raw_chi2_min = 0.0d0
    if (raw_chi2_max == -huge(1.0d0)) raw_chi2_max = 0.0d0
    dg_frag%wpw_partition_sum_chi2_min = chi2_min
    dg_frag%wpw_partition_sum_chi2_max = chi2_max
    dg_frag%wpw_partition_sum_chi2_maxdev = chi2_maxdev
    dg_frag%has_wpw_window = .true.

    if (wpw_trace_enabled()) then
      write(*,'(1x,a,i0,a,i0,a,3es12.4,a,3es12.4,a,3(1x,i0))') &
        "[DG-WPW] rank=", info%id_r, " nlocal=", nlocal, &
        " chi2(min,max,maxdev)=", chi2_min, chi2_max, chi2_maxdev, &
        " raw_chi2(min,max,maxdev)=", raw_chi2_min, raw_chi2_max, raw_chi2_maxdev, &
        " wpw_buffer=", wpw_window_buffer_axis(dg_frag, 1), wpw_window_buffer_axis(dg_frag, 2), &
        wpw_window_buffer_axis(dg_frag, 3)
    end if

    if (wpw_local_diag_enabled() .or. wpw_pp_blocks_enabled() .or. wpw_pp_prop_diag_enabled() .or. &
        wpw_mixed_block_diag_enabled() .or. wpw_mixed_h_diag_enabled() .or. wpw_mixed_neighbor_h_diag_enabled()) &
      call diagnose_wpw_local_ag_volume(dg_frag)

  end subroutine initialize_wpw_windows





  !=======================================================================
  ! Prepare mixed-basis operator blocks at startup without dense EVP.
  ! Fragment initial states are kept as-is and PW coefficients start from zero.
  !=======================================================================
  subroutine prepare_mixed_basis_startup(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot, mg, stencil)
    use structures
    use rt_dg_fragment_ops, only: ensure_gradient_basis_cache
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    real(8), intent(in) :: Ac_tot(3)
    type(s_rgrid), optional, intent(in) :: mg
    type(s_stencil), optional, intent(in) :: stencil

    integer :: n_frag, n_pw, n_tot, ispin, irow
    complex(8), allocatable :: S_frag_pw(:,:,:), H_frag_pw(:,:,:), P_frag_pw(:,:,:,:)

    n_frag = dg_frag%n_mat_max
    if (allocated(dg_frag%coef)) n_frag = size(dg_frag%coef, 1)
    n_pw = dg_frag%n_plane_waves
    if (.not. dg_frag%use_plane_wave_basis .or. n_pw <= 0) return

    allocate(S_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(H_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(P_frag_pw(3, n_frag, n_pw, dg_frag%nspin))

    call compute_fragment_pw_overlap(dg_frag, S_frag_pw)
    call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_frag_pw)
    call compute_fragment_pw_gradient_from_overlap(dg_frag, S_frag_pw, P_frag_pw)
    if ((wpw_mixed_h_diag_enabled() .or. wpw_mixed_neighbor_h_diag_enabled()) .and. &
        present(mg) .and. present(stencil)) then
      call ensure_gradient_basis_cache(dg_frag, mg, stencil)
    end if

    if (.not. allocated(dg_frag%S_mat_frag_pw)) then
      allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    else if (size(dg_frag%S_mat_frag_pw,1) /= n_frag .or. size(dg_frag%S_mat_frag_pw,2) /= n_pw .or. &
             size(dg_frag%S_mat_frag_pw,3) /= dg_frag%nspin) then
      deallocate(dg_frag%S_mat_frag_pw)
      allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    end if
    dg_frag%S_mat_frag_pw(:, :, :) = S_frag_pw(:, :, :)

    if (.not. allocated(dg_frag%H_mat_frag_pw)) then
      allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    else if (size(dg_frag%H_mat_frag_pw,1) /= n_frag .or. size(dg_frag%H_mat_frag_pw,2) /= n_pw .or. &
             size(dg_frag%H_mat_frag_pw,3) /= dg_frag%nspin) then
      deallocate(dg_frag%H_mat_frag_pw)
      allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    end if
    dg_frag%H_mat_frag_pw(:, :, :) = H_frag_pw(:, :, :)

    if (.not. allocated(dg_frag%P_mat_frag_pw)) then
      allocate(dg_frag%P_mat_frag_pw(3, n_frag, n_pw, dg_frag%nspin))
    else if (size(dg_frag%P_mat_frag_pw,1) /= 3 .or. size(dg_frag%P_mat_frag_pw,2) /= n_frag .or. &
             size(dg_frag%P_mat_frag_pw,3) /= n_pw .or. size(dg_frag%P_mat_frag_pw,4) /= dg_frag%nspin) then
      deallocate(dg_frag%P_mat_frag_pw)
      allocate(dg_frag%P_mat_frag_pw(3, n_frag, n_pw, dg_frag%nspin))
    end if
    dg_frag%P_mat_frag_pw(:, :, :, :) = P_frag_pw(:, :, :, :)

    call build_mixed_hamiltonian(dg_frag, dg_frag%lg, Vh, Vxc, Vpsl, Ac_tot, S_frag_pw, H_frag_pw)

    n_tot = n_frag + n_pw
    if (.not. allocated(dg_frag%mixed_basis_dim)) then
      allocate(dg_frag%mixed_basis_dim(dg_frag%nspin))
    else if (size(dg_frag%mixed_basis_dim) /= dg_frag%nspin) then
      deallocate(dg_frag%mixed_basis_dim)
      allocate(dg_frag%mixed_basis_dim(dg_frag%nspin))
    end if
    dg_frag%mixed_basis_dim(1:dg_frag%nspin) = n_tot

    if (.not. allocated(dg_frag%mixed_transform)) then
      allocate(dg_frag%mixed_transform(n_tot, n_tot, dg_frag%nspin))
    else if (size(dg_frag%mixed_transform, 1) /= n_tot .or. &
             size(dg_frag%mixed_transform, 2) /= n_tot .or. &
             size(dg_frag%mixed_transform, 3) /= dg_frag%nspin) then
      deallocate(dg_frag%mixed_transform)
      allocate(dg_frag%mixed_transform(n_tot, n_tot, dg_frag%nspin))
    end if
    dg_frag%mixed_transform(:, :, :) = (0.0d0, 0.0d0)
    do ispin = 1, dg_frag%nspin
      do irow = 1, n_tot
        dg_frag%mixed_transform(irow, irow, ispin) = (1.0d0, 0.0d0)
      end do
    end do

    if (.not. allocated(dg_frag%coef_mix)) then
      allocate(dg_frag%coef_mix(n_tot, dg_frag%nstate_tot, dg_frag%nspin))
    else if (size(dg_frag%coef_mix, 1) /= n_tot .or. &
             size(dg_frag%coef_mix, 2) /= dg_frag%nstate_tot .or. &
             size(dg_frag%coef_mix, 3) /= dg_frag%nspin) then
      deallocate(dg_frag%coef_mix)
      allocate(dg_frag%coef_mix(n_tot, dg_frag%nstate_tot, dg_frag%nspin))
    end if
    dg_frag%coef_mix(:, :, :) = (0.0d0, 0.0d0)
    do ispin = 1, dg_frag%nspin
      if (allocated(dg_frag%coef)) then
        dg_frag%coef_mix(1:min(n_frag, size(dg_frag%coef, 1)), 1:min(dg_frag%nstate_tot, size(dg_frag%coef, 2)), ispin) = &
          dg_frag%coef(1:min(n_frag, size(dg_frag%coef, 1)), 1:min(dg_frag%nstate_tot, size(dg_frag%coef, 2)), ispin)
      end if
    end do
    dg_frag%mixed_basis_ready = .true.
    dg_frag%mixed_basis_identity_raw = .true.

    if (allocated(dg_frag%coef_pw)) dg_frag%coef_pw(:, :, :) = (0.0d0, 0.0d0)

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "  [init] Prepared mixed FP/PP operators without dense diagonalization"
      write(*,'(1x,a)') "  [init] Mixed density uses raw FP+PW identity transform"
      write(*,'(1x,a)') "  [init] PW coefficients start from zero and couple in during propagation"
    end if

    deallocate(S_frag_pw, H_frag_pw, P_frag_pw)
  end subroutine prepare_mixed_basis_startup

  !=======================================================================
  ! Compute overlap matrix between fragment basis and plane waves
  !=======================================================================
  subroutine compute_fragment_pw_overlap(dg_frag, S_complex)
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(out) :: S_complex(:,:,:)  ! (n_mat_max, n_pw, nspin)

    integer :: env_len, env_stat
    character(len=64) :: env_sfp_mode
    logical :: use_fft_gspace

    env_sfp_mode = ''
    call get_environment_variable('SALMON_DG_SFP_MODE', env_sfp_mode, length=env_len, status=env_stat)
    use_fft_gspace = .false.
    if (env_stat == 0 .and. env_len > 0) then
      if (env_sfp_mode(1:1) == 'f' .or. env_sfp_mode(1:1) == 'F' .or. env_sfp_mode(1:1) == 'g' .or. env_sfp_mode(1:1) == 'G') then
        use_fft_gspace = .true.
      end if
    end if

    if (use_fft_gspace) then
      call compute_fragment_pw_overlap_fft_gspace(dg_frag, S_complex)
      if (.not. allocated(dg_frag%S_mat_frag_pw_g)) then
        allocate(dg_frag%S_mat_frag_pw_g(size(S_complex,1), size(S_complex,2), size(S_complex,3)))
      else if (size(dg_frag%S_mat_frag_pw_g,1) /= size(S_complex,1) .or. size(dg_frag%S_mat_frag_pw_g,2) /= size(S_complex,2) .or. &
               size(dg_frag%S_mat_frag_pw_g,3) /= size(S_complex,3)) then
        deallocate(dg_frag%S_mat_frag_pw_g)
        allocate(dg_frag%S_mat_frag_pw_g(size(S_complex,1), size(S_complex,2), size(S_complex,3)))
      end if
      dg_frag%S_mat_frag_pw_g(:, :, :) = S_complex(:, :, :)
    else
      call compute_fragment_pw_overlap_realspace(dg_frag, S_complex)
    end if

  end subroutine compute_fragment_pw_overlap

  subroutine compute_fragment_pw_overlap_fft_gspace(dg_frag, S_complex)
    use structures, only: s_parallel_info
    use communication, only: comm_bcast, COMM_GROUP_NULL
    use rt_dg_fragment_parallel, only: init_fragment_poisson_info, finalize_fragment_parallel
    use salmon_global, only: nproc_rgrid
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(out) :: S_complex(:,:,:)  ! (n_mat_max, n_pw, nspin)

    integer :: ifrag, i_local, ispin, io, ig, iloc, ipw
    integer :: ix, iy, iz, gx, gy, gz, bx, by, bz
    integer :: nx, ny, nz, gx0, gy0, gz0
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: env_len, env_stat
    integer :: ik_mode(3)
    integer, allocatable :: k_fft_slot(:,:), k_fft_slot_neg(:,:)
    integer :: nfft1, nfft2, nfft3
    real(8) :: vol_elem, Lbox(3), sqrt_V, fft_scale, s_fp_norm
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    complex(8) :: raw_coeff
    complex(8), allocatable :: fft_in(:,:,:), fft_out(:,:,:)
    complex(8), allocatable :: frag_block(:,:,:)
    type(s_parallel_info) :: info_fft
    logical :: use_complex_basis, owns_ifrag, participates_ifrag_fft, fft_info_ready
    logical, save :: fp_domain_initialized = .false.
    logical, save :: use_buffered_domain = .false.
    character(len=64), save :: fp_domain_mode = 'core'
    character(len=64) :: env_fp_domain

    if (.not. dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return
    if (dg_frag%n_plane_waves <= 0) return

    vol_elem = product(dg_frag%hgs(1:3))
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))
    sqrt_V = sqrt(product(Lbox))
    fft_scale = vol_elem / sqrt_V
    nfft1 = dg_frag%lgnum_total(1)
    nfft2 = dg_frag%lgnum_total(2)
    nfft3 = dg_frag%lgnum_total(3)
    use_complex_basis = allocated(dg_frag%phi_frag_c)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_ub3 = ubound(dg_frag%phi_frag, 3)

    if (.not. fp_domain_initialized) then
      env_fp_domain = ''
      call get_environment_variable('SALMON_DG_FP_INTEGRAL_DOMAIN', env_fp_domain, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        fp_domain_mode = adjustl(env_fp_domain(1:env_len))
      else
        ! Match the real-space density materialization domain by default.
        ! Buffered FP integrals can still be requested explicitly for tests.
        fp_domain_mode = 'core'
      end if
      select case (fp_domain_mode(1:1))
      case ('b','B','p','P','f','F','1')
        use_buffered_domain = .true.
      case default
        use_buffered_domain = .false.
      end select
      fp_domain_initialized = .true.
      if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
        write(*,'(1x,a,a)') '[FP-DOMAIN] S_fp integral domain mode = ', trim(fp_domain_mode)
      end if
    end if

    allocate(k_fft_slot(3, dg_frag%n_plane_waves))
    allocate(k_fft_slot_neg(3, dg_frag%n_plane_waves))
    do ipw = 1, dg_frag%n_plane_waves
      ik_mode(1) = nint(dg_frag%k_pw(1, ipw) * Lbox(1) / (2.0d0 * pi))
      ik_mode(2) = nint(dg_frag%k_pw(2, ipw) * Lbox(2) / (2.0d0 * pi))
      ik_mode(3) = nint(dg_frag%k_pw(3, ipw) * Lbox(3) / (2.0d0 * pi))
      k_fft_slot(1, ipw) = map_fft_mode_to_slot_pw(ik_mode(1), nfft1)
      k_fft_slot(2, ipw) = map_fft_mode_to_slot_pw(ik_mode(2), nfft2)
      k_fft_slot(3, ipw) = map_fft_mode_to_slot_pw(ik_mode(3), nfft3)
      k_fft_slot_neg(1, ipw) = map_fft_mode_to_slot_pw(-ik_mode(1), nfft1)
      k_fft_slot_neg(2, ipw) = map_fft_mode_to_slot_pw(-ik_mode(2), nfft2)
      k_fft_slot_neg(3, ipw) = map_fft_mode_to_slot_pw(-ik_mode(3), nfft3)
    end do

    allocate(fft_in(nfft1, nfft2, nfft3))
    allocate(fft_out(nfft1, nfft2, nfft3))
    allocate(frag_block(dg_frag%nstate_frag, dg_frag%n_plane_waves, dg_frag%nspin))

    fft_info_ready = .false.
    if (dg_frag%icomm_frag /= COMM_GROUP_NULL) then
      info_fft%nprgrid(1:3) = nproc_rgrid(1:3)
      call init_fragment_poisson_info(dg_frag%icomm_frag, info_fft)
      fft_info_ready = .true.
    end if

    fft_in(:, :, :) = (0.0d0, 0.0d0)
    fft_out(:, :, :) = (0.0d0, 0.0d0)

    S_complex(:, :, :) = (0.0d0, 0.0d0)

    do ifrag = 1, dg_frag%n_frag
      owns_ifrag = (ifrag >= dg_frag%ifrag_start .and. ifrag <= dg_frag%ifrag_end)
      participates_ifrag_fft = fft_info_ready .and. dg_frag%ifrag_group == ifrag
      frag_block(:, :, :) = (0.0d0, 0.0d0)

      if (use_buffered_domain .and. allocated(dg_frag%frag_buf_lo) .and. allocated(dg_frag%frag_buf_hi)) then
        gx0 = dg_frag%frag_buf_lo(1, ifrag)
        gy0 = dg_frag%frag_buf_lo(2, ifrag)
        gz0 = dg_frag%frag_buf_lo(3, ifrag)
        nx = dg_frag%frag_buf_hi(1, ifrag) - dg_frag%frag_buf_lo(1, ifrag) + 1
        ny = dg_frag%frag_buf_hi(2, ifrag) - dg_frag%frag_buf_lo(2, ifrag) + 1
        nz = dg_frag%frag_buf_hi(3, ifrag) - dg_frag%frag_buf_lo(3, ifrag) + 1
      else
        if (use_buffered_domain) then
          if (dg_frag%id == 0) then
            write(*,'(1x,a)') '[FP-DOMAIN] buffered FP domain requested but fragment buffer bounds are unavailable'
          end if
          stop "DG-Fragment RT: missing fragment buffer bounds for buffered FP domain"
        end if
        gx0 = dg_frag%ixyz_frag(1, ifrag)
        gy0 = dg_frag%ixyz_frag(2, ifrag)
        gz0 = dg_frag%ixyz_frag(3, ifrag)
        nx = dg_frag%nxyz_domain(1, ifrag)
        ny = dg_frag%nxyz_domain(2, ifrag)
        nz = dg_frag%nxyz_domain(3, ifrag)
      end if

      if (owns_ifrag) i_local = ifrag - dg_frag%ifrag_start + 1
      if (participates_ifrag_fft) then
        fft_in(:, :, :) = (0.0d0, 0.0d0)
        fft_out(:, :, :) = (0.0d0, 0.0d0)
        call PZFFT3DV_MOD(fft_in, fft_out, nfft1, nfft2, nfft3, info_fft%isize_y, info_fft%isize_z, 0, &
                          info_fft%icomm_y, info_fft%icomm_z)
      end if

      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%n_basis(ifrag, ispin)
          fft_in(:, :, :) = (0.0d0, 0.0d0)

          if (owns_ifrag) then
            if (use_complex_basis) then
              do iz = 1, nz
                gz = gz0 + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, nfft3)
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = 1, ny
                  gy = gy0 + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, nfft2)
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  do ix = 1, nx
                    gx = gx0 + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, nfft1)
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    fft_in(map_global_to_fft_slot_pw(gx, nfft1), map_global_to_fft_slot_pw(gy, nfft2), &
                           map_global_to_fft_slot_pw(gz, nfft3)) = &
                         fft_in(map_global_to_fft_slot_pw(gx, nfft1), map_global_to_fft_slot_pw(gy, nfft2), &
                                map_global_to_fft_slot_pw(gz, nfft3)) + conjg(dg_frag%phi_frag_c(bx, by, bz, io, i_local))
                  end do
                end do
              end do
            else
              do iz = 1, nz
                gz = gz0 + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, nfft3)
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = 1, ny
                  gy = gy0 + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, nfft2)
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  do ix = 1, nx
                    gx = gx0 + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, nfft1)
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    fft_in(map_global_to_fft_slot_pw(gx, nfft1), map_global_to_fft_slot_pw(gy, nfft2), &
                           map_global_to_fft_slot_pw(gz, nfft3)) = &
                         fft_in(map_global_to_fft_slot_pw(gx, nfft1), map_global_to_fft_slot_pw(gy, nfft2), &
                                map_global_to_fft_slot_pw(gz, nfft3)) + cmplx(dg_frag%phi_frag(bx, by, bz, io, i_local), 0.0d0, kind=8)
                  end do
                end do
              end do
            end if
          end if

          if (participates_ifrag_fft) then
            call PZFFT3DV_MOD(fft_in, fft_out, nfft1, nfft2, nfft3, info_fft%isize_y, info_fft%isize_z, -1, &
                              info_fft%icomm_y, info_fft%icomm_z)
          end if

          if (owns_ifrag) then
            do ipw = 1, dg_frag%n_plane_waves
              if (use_complex_basis) then
                raw_coeff = fft_out(k_fft_slot_neg(1, ipw), k_fft_slot_neg(2, ipw), k_fft_slot_neg(3, ipw))
              else
                raw_coeff = fft_out(k_fft_slot(1, ipw), k_fft_slot(2, ipw), k_fft_slot(3, ipw))
              end if
              frag_block(io, ipw, ispin) = raw_coeff * fft_scale
            end do
          end if
        end do
      end do

      call comm_bcast(frag_block, dg_frag%icomm, dg_frag%id_array(ifrag))
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%n_basis(ifrag, ispin)
          ig = dg_frag%index_basis(io, ifrag, ispin)
          iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(S_complex, 1))
          if (iloc <= 0) cycle
          S_complex(iloc, 1:dg_frag%n_plane_waves, ispin) = frag_block(io, 1:dg_frag%n_plane_waves, ispin)
        end do
      end do
    end do

    s_fp_norm = sqrt(sum(abs(S_complex(1:size(S_complex, 1), 1:dg_frag%n_plane_waves, 1:dg_frag%nspin))**2))
    if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
      write(*,'(1x,a,1x,1pe14.6)') '[FP-DOMAIN] ||S_fp||_F =', s_fp_norm
    end if

    if (fft_info_ready) call finalize_fragment_parallel(info_fft)

    deallocate(fft_in, fft_out, frag_block, k_fft_slot, k_fft_slot_neg)
  end subroutine compute_fragment_pw_overlap_fft_gspace

  subroutine compute_fragment_pw_overlap_realspace(dg_frag, S_complex)
    use structures
    use communication, only: comm_bcast, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(out) :: S_complex(:,:,:)  ! (n_mat_max, n_pw, nspin)

    integer :: ipw, ifrag, i_local, io, ig, iloc, ispin, ix, iy, iz
    integer :: nx, ny, nz, nx_max, ny_max, nz_max
    integer :: gx, gy, gz, gx0, gy0, gz0, bx, by, bz
    integer :: loc_s(3), loc_e(3)
    integer :: ipw_s, ipw_e
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: env_len, env_stat
    real(8) :: k_vec(3), Lbox(3), sqrt_V, inv_sqrt_V
    real(8) :: vol_elem, s_fp_norm
    complex(8) :: pw_val, overlap_local, phase_yz
    complex(8) :: step_x, step_y, step_z
    complex(8) :: phase_x0, phase_y0, phase_z0
    complex(8), allocatable :: frag_block(:,:,:), frag_block_sum(:,:,:)
    complex(8), allocatable :: phase_x(:), phase_y(:), phase_z(:)
    logical :: use_complex_basis
    logical :: owns_pw_col
    logical, save :: fp_domain_initialized = .false.
    logical, save :: use_buffered_domain = .false.
    character(len=64), save :: fp_domain_mode = 'core'
    character(len=64) :: env_fp_domain

    if (.not. dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return

    if (dg_frag%ifrag_start > dg_frag%ifrag_end) return

    vol_elem = product(dg_frag%hgs(1:3))
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))
    sqrt_V = sqrt(product(Lbox))
    inv_sqrt_V = 1.0d0 / sqrt_V
    use_complex_basis = allocated(dg_frag%phi_frag_c)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_ub3 = ubound(dg_frag%phi_frag, 3)

    if (.not. fp_domain_initialized) then
      env_fp_domain = ''
      call get_environment_variable('SALMON_DG_FP_INTEGRAL_DOMAIN', env_fp_domain, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        fp_domain_mode = adjustl(env_fp_domain(1:env_len))
      else
        ! Match the real-space density materialization domain by default.
        ! Buffered FP integrals can still be requested explicitly for tests.
        fp_domain_mode = 'core'
      end if
      select case (fp_domain_mode(1:1))
      case ('b','B','p','P','f','F','1')
        use_buffered_domain = .true.
      case default
        use_buffered_domain = .false.
      end select
      fp_domain_initialized = .true.
      if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
        write(*,'(1x,a,a)') '[FP-DOMAIN] S_fp integral domain mode = ', trim(fp_domain_mode)
      end if
    end if

    if (use_buffered_domain .and. allocated(dg_frag%frag_buf_lo) .and. allocated(dg_frag%frag_buf_hi)) then
      nx_max = max(1, maxval(dg_frag%frag_buf_hi(1, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(1, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
      ny_max = max(1, maxval(dg_frag%frag_buf_hi(2, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(2, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
      nz_max = max(1, maxval(dg_frag%frag_buf_hi(3, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(3, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
    else
      if (use_buffered_domain) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a)') '[FP-DOMAIN] buffered FP domain requested but fragment buffer bounds are unavailable'
        end if
        stop "DG-Fragment RT: missing fragment buffer bounds for buffered FP domain"
      end if
      nx_max = maxval(dg_frag%nxyz_domain(1, dg_frag%ifrag_start:dg_frag%ifrag_end))
      ny_max = maxval(dg_frag%nxyz_domain(2, dg_frag%ifrag_start:dg_frag%ifrag_end))
      nz_max = maxval(dg_frag%nxyz_domain(3, dg_frag%ifrag_start:dg_frag%ifrag_end))
    end if
    allocate(phase_x(nx_max), phase_y(ny_max), phase_z(nz_max))

    S_complex = (0.0d0, 0.0d0)
    call get_fragment_pw_column_range(dg_frag, dg_frag%n_plane_waves, ipw_s, ipw_e)

    do ispin = 1, dg_frag%nspin
      do ipw = 1, dg_frag%n_plane_waves
        owns_pw_col = (.not. dg_frag%parallel_mode_orbital) .or. (ipw >= ipw_s .and. ipw <= ipw_e)
        if (.not. owns_pw_col) cycle
        k_vec(1:3) = dg_frag%k_pw(1:3, ipw)

        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          if (use_buffered_domain .and. allocated(dg_frag%frag_buf_lo) .and. allocated(dg_frag%frag_buf_hi)) then
            gx0 = dg_frag%frag_buf_lo(1, ifrag)
            gy0 = dg_frag%frag_buf_lo(2, ifrag)
            gz0 = dg_frag%frag_buf_lo(3, ifrag)
            nx = dg_frag%frag_buf_hi(1, ifrag) - dg_frag%frag_buf_lo(1, ifrag) + 1
            ny = dg_frag%frag_buf_hi(2, ifrag) - dg_frag%frag_buf_lo(2, ifrag) + 1
            nz = dg_frag%frag_buf_hi(3, ifrag) - dg_frag%frag_buf_lo(3, ifrag) + 1
          else
            nx = dg_frag%nxyz_domain(1, ifrag)
            ny = dg_frag%nxyz_domain(2, ifrag)
            nz = dg_frag%nxyz_domain(3, ifrag)
            gx0 = dg_frag%ixyz_frag(1, ifrag)
            gy0 = dg_frag%ixyz_frag(2, ifrag)
            gz0 = dg_frag%ixyz_frag(3, ifrag)
          end if

          step_x = exp(cmplx(0.0d0, k_vec(1) * dg_frag%hgs(1), kind=8))
          step_y = exp(cmplx(0.0d0, k_vec(2) * dg_frag%hgs(2), kind=8))
          step_z = exp(cmplx(0.0d0, k_vec(3) * dg_frag%hgs(3), kind=8))
          phase_x0 = exp(cmplx(0.0d0, k_vec(1) * dble(gx0) * dg_frag%hgs(1), kind=8))
          phase_y0 = exp(cmplx(0.0d0, k_vec(2) * dble(gy0) * dg_frag%hgs(2), kind=8))
          phase_z0 = exp(cmplx(0.0d0, k_vec(3) * dble(gz0) * dg_frag%hgs(3), kind=8))
          phase_x(1) = phase_x0
          do ix = 2, nx
            phase_x(ix) = phase_x(ix-1) * step_x
          end do
          phase_y(1) = phase_y0
          do iy = 2, ny
            phase_y(iy) = phase_y(iy-1) * step_y
          end do
          phase_z(1) = phase_z0
          do iz = 2, nz
            phase_z(iz) = phase_z(iz-1) * step_z
          end do

          if (dg_frag%parallel_mode_orbital) then
            ! Orbital mode distributes fragment-PW matrix work by PW columns.
            ! The fragment box is replicated; it is not split into real-space boxes.
            loc_s(:) = [1, 1, 1]
            loc_e(:) = [nx, ny, nz]
          else
            call get_fragment_subgroup_box_range_pw(dg_frag, [nx, ny, nz], loc_s, loc_e)
          end if
          if (any(loc_s(:) > loc_e(:))) cycle

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig = dg_frag%index_basis(io, ifrag, ispin)
            iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(S_complex, 1))
            if (iloc <= 0) cycle
            overlap_local = (0.0d0, 0.0d0)

            if (use_complex_basis) then
              do iz = loc_s(3), loc_e(3)
                gz = gz0 + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = loc_s(2), loc_e(2)
                  gy = gy0 + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  phase_yz = phase_y(iy) * phase_z(iz)
                  do ix = loc_s(1), loc_e(1)
                    gx = gx0 + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    pw_val = phase_x(ix) * phase_yz * inv_sqrt_V
                    overlap_local = overlap_local + &
                         conjg(dg_frag%phi_frag_c(bx,by,bz,io,i_local)) * pw_val * vol_elem
                  end do
                end do
              end do
            else
              ! Real fragment orbitals have conjg(phi)=phi, so <phi|PW>
              ! must use the same +ik phase as density reconstruction.
              do iz = loc_s(3), loc_e(3)
                gz = gz0 + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = loc_s(2), loc_e(2)
                  gy = gy0 + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  phase_yz = phase_y(iy) * phase_z(iz)
                  do ix = loc_s(1), loc_e(1)
                    gx = gx0 + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    pw_val = phase_x(ix) * phase_yz * inv_sqrt_V
                    overlap_local = overlap_local + &
                         dg_frag%phi_frag(bx,by,bz,io,i_local) * pw_val * vol_elem
                  end do
                end do
              end do
            end if

            S_complex(iloc, ipw, ispin) = S_complex(iloc, ipw, ispin) + overlap_local
          end do
        end do
      end do
    end do

    if (allocated(phase_x)) deallocate(phase_x)
    if (allocated(phase_y)) deallocate(phase_y)
    if (allocated(phase_z)) deallocate(phase_z)

    allocate(frag_block(dg_frag%nstate_frag, dg_frag%n_plane_waves, dg_frag%nspin))
    allocate(frag_block_sum(dg_frag%nstate_frag, dg_frag%n_plane_waves, dg_frag%nspin))
    do ifrag = 1, dg_frag%n_frag
      frag_block(:, :, :) = (0.0d0, 0.0d0)
      if (ifrag >= dg_frag%ifrag_start .and. ifrag <= dg_frag%ifrag_end) then
        i_local = ifrag - dg_frag%ifrag_start + 1
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig = dg_frag%index_basis(io, ifrag, ispin)
            iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(S_complex, 1))
            if (iloc <= 0) cycle
            frag_block(io, 1:dg_frag%n_plane_waves, ispin) = S_complex(iloc, 1:dg_frag%n_plane_waves, ispin)
          end do
        end do
        ! Fragment-PW overlaps are linear integrals.  Legacy mode reduces
        ! real-space pieces; orbital mode reduces disjoint PW-column pieces.
        if (dg_frag%isize_frag > 1) then
          call comm_summation(frag_block, frag_block_sum, size(frag_block), dg_frag%icomm_frag)
          frag_block(:, :, :) = frag_block_sum(:, :, :)
        end if
      end if
      call comm_bcast(frag_block, dg_frag%icomm, dg_frag%id_array(ifrag))
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%n_basis(ifrag, ispin)
          ig = dg_frag%index_basis(io, ifrag, ispin)
          iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(S_complex, 1))
          if (iloc <= 0) cycle
          S_complex(iloc, 1:dg_frag%n_plane_waves, ispin) = frag_block(io, 1:dg_frag%n_plane_waves, ispin)
        end do
      end do
    end do
    deallocate(frag_block, frag_block_sum)

    s_fp_norm = sqrt(sum(abs(S_complex(1:size(S_complex, 1), 1:dg_frag%n_plane_waves, 1:dg_frag%nspin))**2))
    if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
      write(*,'(1x,a,1x,1pe14.6)') '[FP-DOMAIN] ||S_fp||_F =', s_fp_norm
    end if

  end subroutine compute_fragment_pw_overlap_realspace

  subroutine compute_fragment_pw_position_overlap(dg_frag, R_complex)
    use structures
    use communication, only: comm_bcast, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(out) :: R_complex(:,:,:,:)  ! (3, n_mat_max, n_pw, nspin)

    integer :: ipw, ifrag, i_local, io, ig, iloc, ispin, ix, iy, iz, idir
    integer :: nx, ny, nz, nx_max, ny_max, nz_max
    integer :: gx, gy, gz, gx0, gy0, gz0, bx, by, bz
    integer :: loc_s(3), loc_e(3)
    integer :: ipw_s, ipw_e
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: env_len, env_stat
    real(8) :: k_vec(3), Lbox(3), sqrt_V, inv_sqrt_V
    real(8) :: vol_elem, rcoord(3), r_norm
    complex(8) :: pw_val, phase_yz, phi_conj
    complex(8) :: step_x, step_y, step_z
    complex(8) :: phase_x0, phase_y0, phase_z0
    complex(8), allocatable :: frag_block(:,:,:,:), frag_block_sum(:,:,:,:)
    complex(8), allocatable :: phase_x(:), phase_y(:), phase_z(:)
    logical :: use_complex_basis
    logical :: owns_pw_col
    logical, save :: fp_domain_initialized = .false.
    logical, save :: use_buffered_domain = .false.
    character(len=64), save :: fp_domain_mode = 'core'
    character(len=64) :: env_fp_domain

    R_complex(:, :, :, :) = (0.0d0, 0.0d0)
    if (.not. dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return
    if (dg_frag%ifrag_start > dg_frag%ifrag_end) return

    vol_elem = product(dg_frag%hgs(1:3))
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))
    sqrt_V = sqrt(product(Lbox))
    inv_sqrt_V = 1.0d0 / sqrt_V
    use_complex_basis = allocated(dg_frag%phi_frag_c)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_ub3 = ubound(dg_frag%phi_frag, 3)

    if (.not. fp_domain_initialized) then
      env_fp_domain = ''
      call get_environment_variable('SALMON_DG_FP_INTEGRAL_DOMAIN', env_fp_domain, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        fp_domain_mode = adjustl(env_fp_domain(1:env_len))
      else
        fp_domain_mode = 'core'
      end if
      select case (fp_domain_mode(1:1))
      case ('b','B','p','P','f','F','1')
        use_buffered_domain = .true.
      case default
        use_buffered_domain = .false.
      end select
      fp_domain_initialized = .true.
      if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
        write(*,'(1x,a,a)') '[FP-DOMAIN] R_fp integral domain mode = ', trim(fp_domain_mode)
      end if
    end if

    if (use_buffered_domain .and. allocated(dg_frag%frag_buf_lo) .and. allocated(dg_frag%frag_buf_hi)) then
      nx_max = max(1, maxval(dg_frag%frag_buf_hi(1, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(1, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
      ny_max = max(1, maxval(dg_frag%frag_buf_hi(2, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(2, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
      nz_max = max(1, maxval(dg_frag%frag_buf_hi(3, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(3, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
    else
      if (use_buffered_domain) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a)') '[FP-DOMAIN] buffered R_fp domain requested but fragment buffer bounds are unavailable'
        end if
        stop "DG-Fragment RT: missing fragment buffer bounds for buffered R_fp domain"
      end if
      nx_max = maxval(dg_frag%nxyz_domain(1, dg_frag%ifrag_start:dg_frag%ifrag_end))
      ny_max = maxval(dg_frag%nxyz_domain(2, dg_frag%ifrag_start:dg_frag%ifrag_end))
      nz_max = maxval(dg_frag%nxyz_domain(3, dg_frag%ifrag_start:dg_frag%ifrag_end))
    end if
    allocate(phase_x(nx_max), phase_y(ny_max), phase_z(nz_max))

    call get_fragment_pw_column_range(dg_frag, dg_frag%n_plane_waves, ipw_s, ipw_e)

    do ispin = 1, dg_frag%nspin
      do ipw = 1, dg_frag%n_plane_waves
        owns_pw_col = (.not. dg_frag%parallel_mode_orbital) .or. (ipw >= ipw_s .and. ipw <= ipw_e)
        if (.not. owns_pw_col) cycle
        k_vec(1:3) = dg_frag%k_pw(1:3, ipw)

        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          if (use_buffered_domain .and. allocated(dg_frag%frag_buf_lo) .and. allocated(dg_frag%frag_buf_hi)) then
            gx0 = dg_frag%frag_buf_lo(1, ifrag)
            gy0 = dg_frag%frag_buf_lo(2, ifrag)
            gz0 = dg_frag%frag_buf_lo(3, ifrag)
            nx = dg_frag%frag_buf_hi(1, ifrag) - dg_frag%frag_buf_lo(1, ifrag) + 1
            ny = dg_frag%frag_buf_hi(2, ifrag) - dg_frag%frag_buf_lo(2, ifrag) + 1
            nz = dg_frag%frag_buf_hi(3, ifrag) - dg_frag%frag_buf_lo(3, ifrag) + 1
          else
            nx = dg_frag%nxyz_domain(1, ifrag)
            ny = dg_frag%nxyz_domain(2, ifrag)
            nz = dg_frag%nxyz_domain(3, ifrag)
            gx0 = dg_frag%ixyz_frag(1, ifrag)
            gy0 = dg_frag%ixyz_frag(2, ifrag)
            gz0 = dg_frag%ixyz_frag(3, ifrag)
          end if

          step_x = exp(cmplx(0.0d0, k_vec(1) * dg_frag%hgs(1), kind=8))
          step_y = exp(cmplx(0.0d0, k_vec(2) * dg_frag%hgs(2), kind=8))
          step_z = exp(cmplx(0.0d0, k_vec(3) * dg_frag%hgs(3), kind=8))
          phase_x0 = exp(cmplx(0.0d0, k_vec(1) * dble(gx0) * dg_frag%hgs(1), kind=8))
          phase_y0 = exp(cmplx(0.0d0, k_vec(2) * dble(gy0) * dg_frag%hgs(2), kind=8))
          phase_z0 = exp(cmplx(0.0d0, k_vec(3) * dble(gz0) * dg_frag%hgs(3), kind=8))
          phase_x(1) = phase_x0
          do ix = 2, nx
            phase_x(ix) = phase_x(ix-1) * step_x
          end do
          phase_y(1) = phase_y0
          do iy = 2, ny
            phase_y(iy) = phase_y(iy-1) * step_y
          end do
          phase_z(1) = phase_z0
          do iz = 2, nz
            phase_z(iz) = phase_z(iz-1) * step_z
          end do

          if (dg_frag%parallel_mode_orbital) then
            loc_s(:) = [1, 1, 1]
            loc_e(:) = [nx, ny, nz]
          else
            call get_fragment_subgroup_box_range_pw(dg_frag, [nx, ny, nz], loc_s, loc_e)
          end if
          if (any(loc_s(:) > loc_e(:))) cycle

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig = dg_frag%index_basis(io, ifrag, ispin)
            iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(R_complex, 2))
            if (iloc <= 0) cycle
            if (use_complex_basis) then
              do iz = loc_s(3), loc_e(3)
                gz = gz0 + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = loc_s(2), loc_e(2)
                  gy = gy0 + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  phase_yz = phase_y(iy) * phase_z(iz)
                  do ix = loc_s(1), loc_e(1)
                    gx = gx0 + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    pw_val = phase_x(ix) * phase_yz * inv_sqrt_V
                    phi_conj = conjg(dg_frag%phi_frag_c(bx,by,bz,io,i_local))
                    rcoord(1) = dble(gx) * dg_frag%hgs(1)
                    rcoord(2) = dble(gy) * dg_frag%hgs(2)
                    rcoord(3) = dble(gz) * dg_frag%hgs(3)
                    do idir = 1, 3
                      R_complex(idir,iloc,ipw,ispin) = R_complex(idir,iloc,ipw,ispin) + &
                        phi_conj * rcoord(idir) * pw_val * vol_elem
                    end do
                  end do
                end do
              end do
            else
              do iz = loc_s(3), loc_e(3)
                gz = gz0 + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = loc_s(2), loc_e(2)
                  gy = gy0 + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  phase_yz = phase_y(iy) * phase_z(iz)
                  do ix = loc_s(1), loc_e(1)
                    gx = gx0 + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    pw_val = phase_x(ix) * phase_yz * inv_sqrt_V
                    phi_conj = cmplx(dg_frag%phi_frag(bx,by,bz,io,i_local), 0.0d0, kind=8)
                    rcoord(1) = dble(gx) * dg_frag%hgs(1)
                    rcoord(2) = dble(gy) * dg_frag%hgs(2)
                    rcoord(3) = dble(gz) * dg_frag%hgs(3)
                    do idir = 1, 3
                      R_complex(idir,iloc,ipw,ispin) = R_complex(idir,iloc,ipw,ispin) + &
                        phi_conj * rcoord(idir) * pw_val * vol_elem
                    end do
                  end do
                end do
              end do
            end if
          end do
        end do
      end do
    end do

    if (allocated(phase_x)) deallocate(phase_x)
    if (allocated(phase_y)) deallocate(phase_y)
    if (allocated(phase_z)) deallocate(phase_z)

    allocate(frag_block(3,dg_frag%nstate_frag,dg_frag%n_plane_waves,dg_frag%nspin))
    allocate(frag_block_sum(3,dg_frag%nstate_frag,dg_frag%n_plane_waves,dg_frag%nspin))
    do ifrag = 1, dg_frag%n_frag
      frag_block(:, :, :, :) = (0.0d0, 0.0d0)
      if (ifrag >= dg_frag%ifrag_start .and. ifrag <= dg_frag%ifrag_end) then
        i_local = ifrag - dg_frag%ifrag_start + 1
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig = dg_frag%index_basis(io, ifrag, ispin)
            iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(R_complex, 2))
            if (iloc <= 0) cycle
            frag_block(1:3,io,1:dg_frag%n_plane_waves,ispin) = R_complex(1:3,iloc,1:dg_frag%n_plane_waves,ispin)
          end do
        end do
        if (dg_frag%isize_frag > 1) then
          call comm_summation(frag_block, frag_block_sum, size(frag_block), dg_frag%icomm_frag)
          frag_block(:, :, :, :) = frag_block_sum(:, :, :, :)
        end if
      end if
      call comm_bcast(frag_block, dg_frag%icomm, dg_frag%id_array(ifrag))
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%n_basis(ifrag, ispin)
          ig = dg_frag%index_basis(io, ifrag, ispin)
          iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(R_complex, 2))
          if (iloc <= 0) cycle
          R_complex(1:3,iloc,1:dg_frag%n_plane_waves,ispin) = frag_block(1:3,io,1:dg_frag%n_plane_waves,ispin)
        end do
      end do
    end do
    deallocate(frag_block, frag_block_sum)

    r_norm = sqrt(sum(abs(R_complex(1:3,1:size(R_complex,2),1:dg_frag%n_plane_waves,1:dg_frag%nspin))**2))
    if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
      write(*,'(1x,a,1x,1pe14.6)') '[FP-DOMAIN] ||R_fp||_F =', r_norm
    end if
  end subroutine compute_fragment_pw_position_overlap

  subroutine compute_wpw_kinetic_weak(dg_frag, use_unit_window, T_pw)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    logical, intent(in) :: use_unit_window
    complex(8), intent(out) :: T_pw(:,:)

    integer :: n_pw, ifrag, ilocal, ix, iy, iz, ipw, jpw
    integer :: nx, ny, nz, gx, gy, gz
    real(8) :: vol_elem, box_volume, phase_arg, chi, grad_chi(3)
    complex(8), allocatable :: grad_phi(:,:), local_block(:,:), global_block(:,:)
    complex(8) :: phase
    complex(8), parameter :: zi = (0.0d0, 1.0d0)

    n_pw = dg_frag%n_plane_waves
    T_pw(:, :) = (0.0d0, 0.0d0)
    if (n_pw <= 0) return
    if (.not. use_unit_window .and. .not. dg_frag%has_wpw_window) return

    vol_elem = product(dg_frag%hgs(1:3))
    box_volume = product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    allocate(grad_phi(3, n_pw), local_block(n_pw, n_pw), global_block(n_pw, n_pw))
    local_block(:, :) = (0.0d0, 0.0d0)

    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      ilocal = ifrag - dg_frag%ifrag_start + 1
      if (use_unit_window) then
        nx = dg_frag%nxyz_domain(1, ifrag)
        ny = dg_frag%nxyz_domain(2, ifrag)
        nz = dg_frag%nxyz_domain(3, ifrag)
      else
        nx = dg_frag%wpw_window_box_hi(1, ilocal) - dg_frag%wpw_window_box_lo(1, ilocal) + 1
        ny = dg_frag%wpw_window_box_hi(2, ilocal) - dg_frag%wpw_window_box_lo(2, ilocal) + 1
        nz = dg_frag%wpw_window_box_hi(3, ilocal) - dg_frag%wpw_window_box_lo(3, ilocal) + 1
      end if

      do iz = 1, nz
        do iy = 1, ny
          do ix = 1, nx
            if (use_unit_window) then
              gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
              gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
              gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
              chi = 1.0d0
              grad_chi(1:3) = 0.0d0
            else
              gx = dg_frag%wpw_window_box_lo(1, ilocal) + ix - 1
              gy = dg_frag%wpw_window_box_lo(2, ilocal) + iy - 1
              gz = dg_frag%wpw_window_box_lo(3, ilocal) + iz - 1
              chi = dg_frag%wpw_chi(ix, iy, iz, ilocal)
              grad_chi(1:3) = dg_frag%wpw_grad_chi(1:3, ix, iy, iz, ilocal)
            end if
            if (abs(chi) <= 1.0d-30 .and. maxval(abs(grad_chi(1:3))) <= 1.0d-30) cycle

            do ipw = 1, n_pw
              phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                          dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                          dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
              phase = exp(cmplx(0.0d0, phase_arg, kind=8))
              grad_phi(1:3, ipw) = (cmplx(grad_chi(1:3), 0.0d0, kind=8) + &
                zi * dg_frag%k_pw(1:3, ipw) * chi) * phase
            end do

            do jpw = 1, n_pw
              do ipw = 1, n_pw
                local_block(ipw, jpw) = local_block(ipw, jpw) + 0.5d0 * &
                  sum(conjg(grad_phi(1:3, ipw)) * grad_phi(1:3, jpw)) * vol_elem
              end do
            end do
          end do
        end do
      end do
    end do

    call comm_summation(local_block, global_block, n_pw * n_pw, dg_frag%icomm)
    T_pw(1:n_pw, 1:n_pw) = global_block(1:n_pw, 1:n_pw) / box_volume

    do jpw = 1, n_pw
      do ipw = jpw + 1, n_pw
        T_pw(ipw, jpw) = 0.5d0 * (T_pw(ipw, jpw) + conjg(T_pw(jpw, ipw)))
        T_pw(jpw, ipw) = conjg(T_pw(ipw, jpw))
      end do
      T_pw(jpw, jpw) = cmplx(real(T_pw(jpw, jpw), 8), 0.0d0, kind=8)
    end do

    deallocate(grad_phi, local_block, global_block)
  end subroutine compute_wpw_kinetic_weak

  subroutine compute_wpw_overlap(dg_frag, use_unit_window, S_pw)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    logical, intent(in) :: use_unit_window
    complex(8), intent(out) :: S_pw(:,:)

    integer :: n_pw, ifrag, ilocal, ix, iy, iz, ipw, jpw
    integer :: nx, ny, nz, gx, gy, gz
    real(8) :: vol_elem, box_volume, phase_arg, chi
    complex(8), allocatable :: phi(:), local_block(:,:), global_block(:,:)

    n_pw = dg_frag%n_plane_waves
    S_pw(:, :) = (0.0d0, 0.0d0)
    if (n_pw <= 0) return
    if (.not. use_unit_window .and. .not. dg_frag%has_wpw_window) return

    vol_elem = product(dg_frag%hgs(1:3))
    box_volume = product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    allocate(phi(n_pw), local_block(n_pw, n_pw), global_block(n_pw, n_pw))
    local_block(:, :) = (0.0d0, 0.0d0)

    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      ilocal = ifrag - dg_frag%ifrag_start + 1
      if (use_unit_window) then
        nx = dg_frag%nxyz_domain(1, ifrag)
        ny = dg_frag%nxyz_domain(2, ifrag)
        nz = dg_frag%nxyz_domain(3, ifrag)
      else
        nx = dg_frag%wpw_window_box_hi(1, ilocal) - dg_frag%wpw_window_box_lo(1, ilocal) + 1
        ny = dg_frag%wpw_window_box_hi(2, ilocal) - dg_frag%wpw_window_box_lo(2, ilocal) + 1
        nz = dg_frag%wpw_window_box_hi(3, ilocal) - dg_frag%wpw_window_box_lo(3, ilocal) + 1
      end if

      do iz = 1, nz
        do iy = 1, ny
          do ix = 1, nx
            if (use_unit_window) then
              gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
              gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
              gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
              chi = 1.0d0
            else
              gx = dg_frag%wpw_window_box_lo(1, ilocal) + ix - 1
              gy = dg_frag%wpw_window_box_lo(2, ilocal) + iy - 1
              gz = dg_frag%wpw_window_box_lo(3, ilocal) + iz - 1
              chi = dg_frag%wpw_chi(ix, iy, iz, ilocal)
            end if
            if (abs(chi) <= 1.0d-30) cycle

            do ipw = 1, n_pw
              phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                          dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                          dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
              phi(ipw) = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
            end do

            do jpw = 1, n_pw
              do ipw = 1, n_pw
                local_block(ipw, jpw) = local_block(ipw, jpw) + conjg(phi(ipw)) * phi(jpw) * vol_elem
              end do
            end do
          end do
        end do
      end do
    end do

    call comm_summation(local_block, global_block, n_pw * n_pw, dg_frag%icomm)
    S_pw(1:n_pw, 1:n_pw) = global_block(1:n_pw, 1:n_pw) / box_volume

    do jpw = 1, n_pw
      do ipw = jpw + 1, n_pw
        S_pw(ipw, jpw) = 0.5d0 * (S_pw(ipw, jpw) + conjg(S_pw(jpw, ipw)))
        S_pw(jpw, ipw) = conjg(S_pw(ipw, jpw))
      end do
      S_pw(jpw, jpw) = cmplx(real(S_pw(jpw, jpw), 8), 0.0d0, kind=8)
    end do

    deallocate(phi, local_block, global_block)
  end subroutine compute_wpw_overlap

  subroutine diagnose_wpw_kinetic(dg_frag, T_pw, S_pw, use_unit_window)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    complex(8), intent(in) :: T_pw(:,:)
    complex(8), intent(in) :: S_pw(:,:)
    logical, intent(in) :: use_unit_window

    integer :: n_pw, ipw, jpw, info, lwork
    real(8) :: herm_max, s_herm_max, min_eval, max_eval, diag_ref_max
    real(8) :: trace_t, max_abs_t, s_min, s_max, s_cond, gen_min, gen_max
    complex(8), allocatable :: work_mat(:,:), work_s(:,:), work_c(:)
    real(8), allocatable :: eval(:), eval_s(:), rwork(:)
    external :: zheev, zhegv

    n_pw = min(size(T_pw, 1), size(T_pw, 2), dg_frag%n_plane_waves)
    if (n_pw <= 0) return
    herm_max = 0.0d0
    s_herm_max = 0.0d0
    diag_ref_max = 0.0d0
    trace_t = 0.0d0
    max_abs_t = 0.0d0
    do jpw = 1, n_pw
      trace_t = trace_t + real(T_pw(jpw, jpw), 8)
      do ipw = 1, n_pw
        herm_max = max(herm_max, abs(T_pw(ipw, jpw) - conjg(T_pw(jpw, ipw))))
        s_herm_max = max(s_herm_max, abs(S_pw(ipw, jpw) - conjg(S_pw(jpw, ipw))))
        max_abs_t = max(max_abs_t, abs(T_pw(ipw, jpw)))
      end do
      diag_ref_max = max(diag_ref_max, abs(T_pw(jpw, jpw) - 0.5d0 * sum(dg_frag%k_pw(:, jpw)**2)))
    end do

    allocate(work_mat(n_pw, n_pw), work_s(n_pw, n_pw), eval(n_pw), eval_s(n_pw), rwork(max(1, 3*n_pw-2)))
    work_mat(:, :) = T_pw(1:n_pw, 1:n_pw)
    lwork = -1
    allocate(work_c(1))
    call ZHEEV('N', 'U', n_pw, work_mat, n_pw, eval, work_c, lwork, rwork, info)
    lwork = max(1, int(real(work_c(1))))
    deallocate(work_c)
    allocate(work_c(lwork))
    work_mat(:, :) = T_pw(1:n_pw, 1:n_pw)
    call ZHEEV('N', 'U', n_pw, work_mat, n_pw, eval, work_c, lwork, rwork, info)
    if (info == 0) then
      min_eval = eval(1)
      max_eval = eval(n_pw)
    else
      min_eval = huge(1.0d0)
      max_eval = -huge(1.0d0)
    end if

    work_s(:, :) = S_pw(1:n_pw, 1:n_pw)
    work_c(1) = (0.0d0, 0.0d0)
    lwork = -1
    call ZHEEV('N', 'U', n_pw, work_s, n_pw, eval_s, work_c, lwork, rwork, info)
    lwork = max(1, int(real(work_c(1))))
    if (size(work_c) /= lwork) then
      deallocate(work_c)
      allocate(work_c(lwork))
    end if
    work_s(:, :) = S_pw(1:n_pw, 1:n_pw)
    call ZHEEV('N', 'U', n_pw, work_s, n_pw, eval_s, work_c, lwork, rwork, info)
    if (info == 0) then
      s_min = eval_s(1)
      s_max = eval_s(n_pw)
      if (s_min > 0.0d0) then
        s_cond = s_max / s_min
      else
        s_cond = huge(1.0d0)
      end if
    else
      s_min = huge(1.0d0)
      s_max = -huge(1.0d0)
      s_cond = huge(1.0d0)
    end if

    work_mat(:, :) = T_pw(1:n_pw, 1:n_pw)
    work_s(:, :) = S_pw(1:n_pw, 1:n_pw)
    lwork = -1
    work_c(1) = (0.0d0, 0.0d0)
    call ZHEGV(1, 'N', 'U', n_pw, work_mat, n_pw, work_s, n_pw, eval, work_c, lwork, rwork, info)
    lwork = max(1, int(real(work_c(1))))
    if (size(work_c) /= lwork) then
      deallocate(work_c)
      allocate(work_c(lwork))
    end if
    work_mat(:, :) = T_pw(1:n_pw, 1:n_pw)
    work_s(:, :) = S_pw(1:n_pw, 1:n_pw)
    call ZHEGV(1, 'N', 'U', n_pw, work_mat, n_pw, work_s, n_pw, eval, work_c, lwork, rwork, info)
    if (info == 0) then
      gen_min = eval(1)
      gen_max = eval(n_pw)
    else
      gen_min = huge(1.0d0)
      gen_max = -huge(1.0d0)
    end if
    deallocate(work_mat, work_s, eval, eval_s, rwork, work_c)

    if (dg_frag%id == 0) then
      write(*,'(1x,a,l1,13(a,1pe12.4))') &
        '[DG-WPW-KINETIC] unit_window=', use_unit_window, &
        ' herm=', herm_max, ' eval_min=', min_eval, ' eval_max=', max_eval, &
        ' trace=', trace_t, ' max_abs=', max_abs_t, ' diag_ref_max=', diag_ref_max, &
        ' S_herm=', s_herm_max, ' S_min=', s_min, ' S_max=', s_max, ' S_cond=', s_cond, &
        ' gen_eval_min=', gen_min, ' gen_eval_max=', gen_max, &
        ' chi2_maxdev=', dg_frag%wpw_partition_sum_chi2_maxdev
    end if
  end subroutine diagnose_wpw_kinetic

  subroutine diagnose_wpw_local_ag_volume(dg_frag, force_store_blocks, quiet)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    logical, intent(in), optional :: force_store_blocks, quiet

    integer :: n_frag, n_pw, n_tot
    integer :: ifrag, jfrag, ipw, jpw, idx, jdx
    integer :: ix, iy, iz, gx, gy, gz
    integer :: nx, ny, nz
    integer :: info, lwork, rank_def
    real(8) :: vol_elem, box_volume, phase_arg
    real(8) :: qsum, sqrt_qsum, qg_sum(3)
    real(8) :: s_herm, t_herm, s_min, s_max, s_cond, gen_min, gen_max
    real(8) :: s_self_max, s_neigh_max, s_non_max
    real(8) :: t_self_max, t_neigh_max, t_non_max
    real(8) :: if_self_max, if_neigh_max, if_non_max
    real(8) :: if_herm, if_gen_min, if_gen_max, if_shift_min, if_shift_max, if_penalty
    real(8), allocatable :: q(:), gq(:,:), chi(:), grad_chi(:,:)
    real(8), allocatable :: eval(:), eval_s(:), rwork(:)
    complex(8), allocatable :: phi(:), grad_phi(:,:)
    complex(8), allocatable :: local_s(:,:), local_t(:,:), S_ag(:,:), T_ag(:,:), T_if(:,:)
    complex(8), allocatable :: work_mat(:,:), work_s(:,:), work_c(:)
    complex(8) :: phase
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    real(8), parameter :: tiny_q = 1.0d-28
    logical :: force_store_mode, quiet_mode
    external :: zheev, zhegv

    force_store_mode = .false.
    if (present(force_store_blocks)) force_store_mode = force_store_blocks
    quiet_mode = .false.
    if (present(quiet)) quiet_mode = quiet

    n_frag = dg_frag%n_frag
    n_pw = dg_frag%n_plane_waves
    n_tot = n_frag * n_pw
    if (n_frag <= 0 .or. n_pw <= 0) return
    if (.not. dg_frag%has_wpw_window) return

    vol_elem = product(dg_frag%hgs(1:3))
    box_volume = product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))

    allocate(q(n_frag), gq(3, n_frag), chi(n_frag), grad_chi(3, n_frag))
    allocate(phi(n_tot), grad_phi(3, n_tot))
    allocate(local_s(n_tot, n_tot), local_t(n_tot, n_tot), S_ag(n_tot, n_tot), T_ag(n_tot, n_tot), T_if(n_tot, n_tot))
    local_s(:, :) = (0.0d0, 0.0d0)
    local_t(:, :) = (0.0d0, 0.0d0)

    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      nx = dg_frag%nxyz_domain(1, ifrag)
      ny = dg_frag%nxyz_domain(2, ifrag)
      nz = dg_frag%nxyz_domain(3, ifrag)
      do iz = 1, nz
        gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
        do iy = 1, ny
          gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
          do ix = 1, nx
            gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1

            q = 0.0d0
            gq = 0.0d0
            do jfrag = 1, n_frag
              call wpw_raw_window_at_grid(dg_frag, jfrag, gx, gy, gz, q(jfrag), gq(:, jfrag))
            end do
            qsum = sum(q(1:n_frag)**2)
            if (qsum <= tiny_q) cycle
            sqrt_qsum = sqrt(qsum)
            qg_sum(1:3) = 0.0d0
            do jfrag = 1, n_frag
              qg_sum(1:3) = qg_sum(1:3) + q(jfrag) * gq(1:3, jfrag)
            end do
            do jfrag = 1, n_frag
              chi(jfrag) = q(jfrag) / sqrt_qsum
              grad_chi(1:3, jfrag) = gq(1:3, jfrag) / sqrt_qsum - &
                q(jfrag) * qg_sum(1:3) / (qsum * sqrt_qsum)
            end do

            do jfrag = 1, n_frag
              do ipw = 1, n_pw
                idx = (jfrag - 1) * n_pw + ipw
                phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                            dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                            dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                phase = exp(cmplx(0.0d0, phase_arg, kind=8))
                phi(idx) = chi(jfrag) * phase
                grad_phi(1:3, idx) = (cmplx(grad_chi(1:3, jfrag), 0.0d0, kind=8) + &
                  zi * dg_frag%k_pw(1:3, ipw) * chi(jfrag)) * phase
              end do
            end do

            do jdx = 1, n_tot
              do idx = 1, n_tot
                local_s(idx, jdx) = local_s(idx, jdx) + conjg(phi(idx)) * phi(jdx) * vol_elem
                local_t(idx, jdx) = local_t(idx, jdx) + 0.5d0 * &
                  sum(conjg(grad_phi(1:3, idx)) * grad_phi(1:3, jdx)) * vol_elem
              end do
            end do
          end do
        end do
      end do
    end do

    call comm_summation(local_s, S_ag, n_tot * n_tot, dg_frag%icomm)
    call comm_summation(local_t, T_ag, n_tot * n_tot, dg_frag%icomm)
    S_ag(:, :) = S_ag(:, :) / box_volume
    T_ag(:, :) = T_ag(:, :) / box_volume

    call hermitize_matrix(S_ag, n_tot)
    call hermitize_matrix(T_ag, n_tot)
    call wpw_local_herm_max(S_ag, n_tot, s_herm)
    call wpw_local_herm_max(T_ag, n_tot, t_herm)
    call wpw_local_block_magnitudes(dg_frag, S_ag, s_self_max, s_neigh_max, s_non_max)
    call wpw_local_block_magnitudes(dg_frag, T_ag, t_self_max, t_neigh_max, t_non_max)

    T_if(:, :) = (0.0d0, 0.0d0)
    if (force_store_mode .or. wpw_interface_diag_enabled() .or. wpw_pp_blocks_enabled() .or. &
        wpw_pp_prop_diag_enabled() .or. wpw_mixed_block_diag_enabled()) then
      if_penalty = wpw_interface_penalty_factor()
      call compute_wpw_local_pp_interface(dg_frag, T_if)
      call hermitize_matrix(T_if, n_tot)
      call wpw_local_herm_max(T_if, n_tot, if_herm)
      call wpw_local_block_magnitudes(dg_frag, T_if, if_self_max, if_neigh_max, if_non_max)
    else
      if_herm = 0.0d0
      if_penalty = 0.0d0
      if_self_max = 0.0d0
      if_neigh_max = 0.0d0
      if_non_max = 0.0d0
    end if

    allocate(work_mat(n_tot, n_tot), work_s(n_tot, n_tot), eval(n_tot), eval_s(n_tot), &
      rwork(max(1, 3*n_tot-2)))
    allocate(work_c(1))

    work_s(:, :) = S_ag(:, :)
    lwork = -1
    call ZHEEV('N', 'U', n_tot, work_s, n_tot, eval_s, work_c, lwork, rwork, info)
    lwork = max(1, int(real(work_c(1))))
    deallocate(work_c)
    allocate(work_c(lwork))
    work_s(:, :) = S_ag(:, :)
    call ZHEEV('N', 'U', n_tot, work_s, n_tot, eval_s, work_c, lwork, rwork, info)
    if (info == 0) then
      s_min = eval_s(1)
      s_max = eval_s(n_tot)
      rank_def = count(eval_s(1:n_tot) < 1.0d-8)
      if (s_min > 0.0d0) then
        s_cond = s_max / s_min
      else
        s_cond = huge(1.0d0)
      end if
    else
      s_min = huge(1.0d0)
      s_max = -huge(1.0d0)
      s_cond = huge(1.0d0)
      rank_def = n_tot
    end if

    work_mat(:, :) = T_ag(:, :)
    work_s(:, :) = S_ag(:, :)
    lwork = -1
    work_c(1) = (0.0d0, 0.0d0)
    call ZHEGV(1, 'N', 'U', n_tot, work_mat, n_tot, work_s, n_tot, eval, work_c, lwork, rwork, info)
    lwork = max(1, int(real(work_c(1))))
    if (size(work_c) /= lwork) then
      deallocate(work_c)
      allocate(work_c(lwork))
    end if
    work_mat(:, :) = T_ag(:, :)
    work_s(:, :) = S_ag(:, :)
    call ZHEGV(1, 'N', 'U', n_tot, work_mat, n_tot, work_s, n_tot, eval, work_c, lwork, rwork, info)
    if (info == 0) then
      gen_min = eval(1)
      gen_max = eval(n_tot)
    else
      gen_min = huge(1.0d0)
      gen_max = -huge(1.0d0)
    end if

    if (wpw_interface_diag_enabled()) then
      work_mat(:, :) = T_ag(:, :) + T_if(:, :)
      work_s(:, :) = S_ag(:, :)
      lwork = -1
      work_c(1) = (0.0d0, 0.0d0)
      call ZHEGV(1, 'N', 'U', n_tot, work_mat, n_tot, work_s, n_tot, eval, work_c, lwork, rwork, info)
      lwork = max(1, int(real(work_c(1))))
      if (size(work_c) /= lwork) then
        deallocate(work_c)
        allocate(work_c(lwork))
      end if
      work_mat(:, :) = T_ag(:, :) + T_if(:, :)
      work_s(:, :) = S_ag(:, :)
      call ZHEGV(1, 'N', 'U', n_tot, work_mat, n_tot, work_s, n_tot, eval, work_c, lwork, rwork, info)
      if (info == 0) then
        if_gen_min = eval(1)
        if_gen_max = eval(n_tot)
      else
        if_gen_min = huge(1.0d0)
        if_gen_max = -huge(1.0d0)
      end if
      if_shift_min = if_gen_min - gen_min
      if_shift_max = if_gen_max - gen_max
    else
      if_gen_min = 0.0d0
      if_gen_max = 0.0d0
      if_shift_min = 0.0d0
      if_shift_max = 0.0d0
    end if

    if (force_store_mode .or. wpw_pp_blocks_enabled() .or. wpw_pp_prop_diag_enabled() .or. &
        wpw_mixed_block_diag_enabled() .or. wpw_mixed_h_diag_enabled() .or. wpw_mixed_neighbor_h_diag_enabled()) &
      call store_wpw_local_pp_blocks(dg_frag, S_ag, T_ag, T_if, quiet_mode)
    if (wpw_pp_prop_diag_enabled()) call diagnose_wpw_pp_only_propagation(dg_frag)
    if (wpw_mixed_block_diag_enabled()) call diagnose_wpw_mixed_self_overlap(dg_frag)

    if (dg_frag%id == 0 .and. .not. quiet_mode) then
      write(*,'(1x,a,i0,a,i0,a,i0,13(a,1pe12.4),a,i0)') &
        '[DG-WPW-LOCAL] dim=', n_tot, ' nfrag=', n_frag, ' npw=', n_pw, &
        ' S_herm=', s_herm, ' T_herm=', t_herm, &
        ' S_eval_min=', s_min, ' S_eval_max=', s_max, ' S_cond=', s_cond, &
        ' gen_eval_min=', gen_min, ' gen_eval_max=', gen_max, &
        ' S_self_block_max=', s_self_max, ' S_neighbor_block_max=', s_neigh_max, &
        ' S_nonneighbor_block_max=', s_non_max, &
        ' T_self_block_max=', t_self_max, ' T_neighbor_block_max=', t_neigh_max, &
        ' T_nonneighbor_block_max=', t_non_max, &
        ' rank_deficiency=', rank_def
      if (wpw_interface_diag_enabled()) then
        write(*,'(1x,a,11(a,1pe12.4))') &
          '[DG-WPW-INTERFACE] PP_only volume_separate=T', &
          ' penalty=', if_penalty, ' herm=', if_herm, &
          ' self_block_max=', if_self_max, ' neighbor_block_max=', if_neigh_max, &
          ' nonneighbor_block_max=', if_non_max, &
          ' neighbor_over_self=', if_neigh_max / max(if_self_max, 1.0d-300), &
          ' gen_eval_min=', if_gen_min, ' gen_eval_max=', if_gen_max, &
          ' shift_min=', if_shift_min, ' shift_max=', if_shift_max, &
          ' interface_over_volume=', max(if_self_max, if_neigh_max) / max(t_self_max, 1.0d-300)
      end if
    end if

    deallocate(q, gq, chi, grad_chi, phi, grad_phi, local_s, local_t, S_ag, T_ag, T_if)
    deallocate(work_mat, work_s, eval, eval_s, rwork, work_c)
  end subroutine diagnose_wpw_local_ag_volume

  subroutine ensure_wpw_local_pp_blocks(dg_frag, quiet)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    logical, intent(in), optional :: quiet
    logical :: quiet_mode

    quiet_mode = .false.
    if (present(quiet)) quiet_mode = quiet
    if (dg_frag%wpw_pp_blocks_ready .and. allocated(dg_frag%wpw_S_pp_blocks) .and. &
        allocated(dg_frag%wpw_T_pp_volume_blocks) .and. &
        allocated(dg_frag%wpw_T_pp_interface_blocks)) return
    call diagnose_wpw_local_ag_volume(dg_frag, .true., quiet_mode)
  end subroutine ensure_wpw_local_pp_blocks

  subroutine store_wpw_local_pp_blocks(dg_frag, S_ag, T_ag, T_if, quiet)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(in) :: S_ag(:,:), T_ag(:,:), T_if(:,:)
    logical, intent(in), optional :: quiet

    integer :: n_pw, nfrag, nblocks
    integer :: ifrag, jfrag, ipw, jpw, ispin, iblk, idx, jdx
    logical :: quiet_mode

    quiet_mode = .false.
    if (present(quiet)) quiet_mode = quiet
    n_pw = dg_frag%n_plane_waves
    nfrag = dg_frag%n_frag
    dg_frag%wpw_pp_blocks_ready = .false.
    if (n_pw <= 0 .or. nfrag <= 0) return
    if (size(S_ag, 1) < n_pw*nfrag .or. size(T_ag, 1) < n_pw*nfrag .or. size(T_if, 1) < n_pw*nfrag) return

    if (allocated(dg_frag%wpw_S_pp_blocks)) deallocate(dg_frag%wpw_S_pp_blocks)
    if (allocated(dg_frag%wpw_T_pp_volume_blocks)) deallocate(dg_frag%wpw_T_pp_volume_blocks)
    if (allocated(dg_frag%wpw_T_pp_interface_blocks)) deallocate(dg_frag%wpw_T_pp_interface_blocks)

    nblocks = 0
    do jfrag = 1, nfrag
      do ifrag = 1, nfrag
        if (.not. wpw_local_is_neighbor_pair(dg_frag, ifrag, jfrag)) cycle
        nblocks = nblocks + 1
      end do
    end do
    if (nblocks <= 0) return

    allocate(dg_frag%wpw_S_pp_blocks(nblocks))
    allocate(dg_frag%wpw_T_pp_volume_blocks(nblocks))
    allocate(dg_frag%wpw_T_pp_interface_blocks(nblocks))

    iblk = 0
    do jfrag = 1, nfrag
      do ifrag = 1, nfrag
        if (.not. wpw_local_is_neighbor_pair(dg_frag, ifrag, jfrag)) cycle
        iblk = iblk + 1
        call allocate_wpw_pp_block(dg_frag%wpw_S_pp_blocks(iblk), ifrag, jfrag, n_pw, dg_frag%nspin)
        call allocate_wpw_pp_block(dg_frag%wpw_T_pp_volume_blocks(iblk), ifrag, jfrag, n_pw, dg_frag%nspin)
        call allocate_wpw_pp_block(dg_frag%wpw_T_pp_interface_blocks(iblk), ifrag, jfrag, n_pw, dg_frag%nspin)
        do ispin = 1, dg_frag%nspin
          do jpw = 1, n_pw
            jdx = (jfrag - 1) * n_pw + jpw
            do ipw = 1, n_pw
              idx = (ifrag - 1) * n_pw + ipw
              dg_frag%wpw_S_pp_blocks(iblk)%val(ipw, jpw, ispin) = S_ag(idx, jdx)
              dg_frag%wpw_T_pp_volume_blocks(iblk)%val(ipw, jpw, ispin) = T_ag(idx, jdx)
              dg_frag%wpw_T_pp_interface_blocks(iblk)%val(ipw, jpw, ispin) = T_if(idx, jdx)
            end do
          end do
        end do
      end do
    end do
    dg_frag%wpw_pp_blocks_ready = .true.

    if (dg_frag%id == 0 .and. .not. quiet_mode) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') &
        '[DG-WPW-BLOCKS] PP fragment-local blocks stored: nblocks=', nblocks, &
        ' nfrag=', nfrag, ' npw_per_frag=', n_pw, &
        ' penalty=', wpw_interface_penalty_factor()
    end if
  end subroutine store_wpw_local_pp_blocks

  subroutine allocate_wpw_pp_block(block, ifrag, jfrag, n_pw, nspin)
    implicit none
    type(complex_matrix_block_info), intent(inout) :: block
    integer, intent(in) :: ifrag, jfrag, n_pw, nspin

    block%ifrag_row = ifrag
    block%ifrag_col = jfrag
    block%nrow_max = n_pw
    block%ncol_max = n_pw
    if (allocated(block%val)) deallocate(block%val)
    allocate(block%val(n_pw, n_pw, nspin))
    block%val(:, :, :) = (0.0d0, 0.0d0)
  end subroutine allocate_wpw_pp_block

  subroutine diagnose_wpw_pp_only_propagation(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag

    integer :: iblk, n_pw, nblocks_self, ispin, info, lwork, i, j
    real(8) :: s_herm, h_herm, s_herm_max, h_herm_max
    real(8) :: eval_min, eval_max, eval_min_all, eval_max_all
    real(8) :: norm0, norm1, norm_dev, norm_dev_max, dt
    complex(8), allocatable :: S(:,:), H(:,:), C(:,:), work(:)
    complex(8), allocatable :: c0(:), c1(:), amp(:), tmp(:)
    real(8), allocatable :: eval(:), rwork(:)
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    external :: zhegv, zheev

    if (.not. dg_frag%wpw_pp_blocks_ready) return
    if (.not. allocated(dg_frag%wpw_S_pp_blocks)) return
    if (.not. allocated(dg_frag%wpw_T_pp_volume_blocks)) return
    if (.not. allocated(dg_frag%wpw_T_pp_interface_blocks)) return

    n_pw = dg_frag%n_plane_waves
    if (n_pw <= 0) return
    dt = 0.1d0

    allocate(S(n_pw, n_pw), H(n_pw, n_pw), C(n_pw, n_pw))
    allocate(c0(n_pw), c1(n_pw), amp(n_pw), tmp(n_pw), eval(n_pw), rwork(max(1, 3*n_pw - 2)))

    nblocks_self = 0
    s_herm_max = 0.0d0
    h_herm_max = 0.0d0
    eval_min_all = huge(1.0d0)
    eval_max_all = -huge(1.0d0)
    norm_dev_max = 0.0d0

    do iblk = 1, size(dg_frag%wpw_S_pp_blocks)
      if (dg_frag%wpw_S_pp_blocks(iblk)%ifrag_row /= dg_frag%wpw_S_pp_blocks(iblk)%ifrag_col) cycle
      nblocks_self = nblocks_self + 1
      do ispin = 1, dg_frag%nspin
        S(:, :) = dg_frag%wpw_S_pp_blocks(iblk)%val(1:n_pw, 1:n_pw, ispin)
        H(:, :) = dg_frag%wpw_T_pp_volume_blocks(iblk)%val(1:n_pw, 1:n_pw, ispin) + &
                  dg_frag%wpw_T_pp_interface_blocks(iblk)%val(1:n_pw, 1:n_pw, ispin)
        call wpw_local_herm_max(S, n_pw, s_herm)
        call wpw_local_herm_max(H, n_pw, h_herm)
        s_herm_max = max(s_herm_max, s_herm)
        h_herm_max = max(h_herm_max, h_herm)

        lwork = -1
        allocate(work(1))
        C(:, :) = H(:, :)
        call ZHEGV(1, 'V', 'U', n_pw, C, n_pw, S, n_pw, eval, work, lwork, rwork, info)
        lwork = max(1, int(real(work(1), kind=8)))
        deallocate(work)
        allocate(work(lwork))
        S(:, :) = dg_frag%wpw_S_pp_blocks(iblk)%val(1:n_pw, 1:n_pw, ispin)
        C(:, :) = H(:, :)
        call ZHEGV(1, 'V', 'U', n_pw, C, n_pw, S, n_pw, eval, work, lwork, rwork, info)
        deallocate(work)
        if (info /= 0) cycle
        eval_min_all = min(eval_min_all, eval(1))
        eval_max_all = max(eval_max_all, eval(n_pw))

        S(:, :) = dg_frag%wpw_S_pp_blocks(iblk)%val(1:n_pw, 1:n_pw, ispin)
        c0(:) = (0.0d0, 0.0d0)
        c0(1) = (1.0d0, 0.0d0)
        tmp(:) = matmul(S, c0)
        norm0 = real(dot_product(c0, tmp), kind=8)
        do i = 1, n_pw
          amp(i) = dot_product(C(:, i), tmp(:))
          amp(i) = amp(i) * exp(-zi * eval(i) * dt)
        end do
        c1(:) = (0.0d0, 0.0d0)
        do i = 1, n_pw
          do j = 1, n_pw
            c1(j) = c1(j) + C(j, i) * amp(i)
          end do
        end do
        tmp(:) = matmul(S, c1)
        norm1 = real(dot_product(c1, tmp), kind=8)
        norm_dev = abs(norm1 - norm0) / max(abs(norm0), 1.0d-300)
        norm_dev_max = max(norm_dev_max, norm_dev)
      end do
    end do

    if (eval_min_all == huge(1.0d0)) eval_min_all = 0.0d0
    if (eval_max_all == -huge(1.0d0)) eval_max_all = 0.0d0
    if (dg_frag%id == 0) then
      write(*,'(1x,a,a,i0,a,i0,6(a,1pe12.4))') &
        '[DG-WPW-PP-PROP] self block generalized exp diagnostic:', &
        ' nself=', nblocks_self, ' npw=', n_pw, &
        ' S_herm=', s_herm_max, ' H_herm=', h_herm_max, &
        ' eval_min=', eval_min_all, ' eval_max=', eval_max_all, &
        ' dt=', dt, ' S_norm_relerr=', norm_dev_max
    end if

    deallocate(S, H, C, c0, c1, amp, tmp, eval, rwork)
  end subroutine diagnose_wpw_pp_only_propagation

  subroutine diagnose_wpw_mixed_self_overlap(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag

    integer :: i_local, ifrag, iw, ib, ipw, ix, iy, iz, gx, gy, gz, idir
    integer :: n_w, n_pw, n_mix, nbf, p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: bx, by, bz, iblk, info, lwork, rank_def
    real(8) :: vol_elem, box_volume, chi, grad_chi(3), phase_arg
    real(8) :: s_herm, s_min, s_max, s_cond, swp_norm, swp_norm_max
    complex(8) :: phase
    complex(8), allocatable :: S_mix(:,:), S_wp(:,:), S_pp(:,:), wval(:), work(:)
    real(8), allocatable :: eval(:), rwork(:)
    external :: zheev

    if (.not. dg_frag%has_global_wannier_local_basis) return
    if (.not. allocated(dg_frag%global_wannier_local_coef)) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return
    if (.not. allocated(dg_frag%wpw_S_pp_blocks)) return
    if (.not. allocated(dg_frag%phi_frag) .and. .not. allocated(dg_frag%phi_frag_c)) return

    n_pw = dg_frag%n_plane_waves
    if (n_pw <= 0) return
    vol_elem = product(dg_frag%hgs(1:3))
    box_volume = product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_ub3 = ubound(dg_frag%phi_frag, 3)

    swp_norm_max = 0.0d0
    do i_local = 1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1
      ifrag = dg_frag%ifrag_start + i_local - 1
      if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
      n_w = dg_frag%global_wannier_local_nkeep(i_local)
      if (n_w <= 0) cycle
      n_mix = n_w + n_pw
      allocate(S_mix(n_mix, n_mix), S_wp(n_w, n_pw), S_pp(n_pw, n_pw), wval(n_w))
      S_mix(:, :) = (0.0d0, 0.0d0)
      S_wp(:, :) = (0.0d0, 0.0d0)
      S_pp(:, :) = (0.0d0, 0.0d0)
      do iw = 1, n_w
        S_mix(iw, iw) = (1.0d0, 0.0d0)
      end do

      iblk = find_wpw_pp_block(dg_frag%wpw_S_pp_blocks, ifrag, ifrag)
      if (iblk <= 0) then
        deallocate(S_mix, S_wp, S_pp, wval)
        cycle
      end if
      S_pp(1:n_pw, 1:n_pw) = dg_frag%wpw_S_pp_blocks(iblk)%val(1:n_pw, 1:n_pw, 1)

      nbf = min(dg_frag%n_basis(ifrag, 1), size(dg_frag%global_wannier_local_coef, 1), size(dg_frag%phi_frag, 4))
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
        bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
        if (bz < p_lb3 .or. bz > p_ub3) cycle
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
          by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          if (by < p_lb2 .or. by > p_ub2) cycle
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
            bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            if (bx < p_lb1 .or. bx > p_ub1) cycle
            wval(:) = (0.0d0, 0.0d0)
            do iw = 1, n_w
              do ib = 1, nbf
                if (allocated(dg_frag%phi_frag_c)) then
                  wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, 1, i_local) * &
                    dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                else
                  wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, 1, i_local) * &
                    cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                end if
              end do
            end do
            call wpw_normalized_window_at_grid(dg_frag, ifrag, gx, gy, gz, chi, grad_chi)
            do ipw = 1, n_pw
              phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                          dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                          dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
              phase = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
              do iw = 1, n_w
                S_wp(iw, ipw) = S_wp(iw, ipw) + conjg(wval(iw)) * phase * vol_elem / box_volume
              end do
            end do
          end do
        end do
      end do

      S_mix(1:n_w, n_w+1:n_mix) = S_wp(1:n_w, 1:n_pw)
      S_mix(n_w+1:n_mix, 1:n_w) = conjg(transpose(S_wp(1:n_w, 1:n_pw)))
      S_mix(n_w+1:n_mix, n_w+1:n_mix) = S_pp(1:n_pw, 1:n_pw)
      call hermitize_matrix(S_mix, n_mix)
      call wpw_local_herm_max(S_mix, n_mix, s_herm)
      allocate(eval(n_mix), rwork(max(1, 3*n_mix - 2)), work(1))
      lwork = -1
      call ZHEEV('N', 'U', n_mix, S_mix, n_mix, eval, work, lwork, rwork, info)
      lwork = max(1, int(real(work(1), kind=8)))
      deallocate(work)
      allocate(work(lwork))
      call ZHEEV('N', 'U', n_mix, S_mix, n_mix, eval, work, lwork, rwork, info)
      if (info == 0) then
        s_min = eval(1)
        s_max = eval(n_mix)
        rank_def = count(eval(1:n_mix) < 1.0d-8)
        if (s_min > 0.0d0) then
          s_cond = s_max / s_min
        else
          s_cond = huge(1.0d0)
        end if
      else
        s_min = 0.0d0
        s_max = 0.0d0
        s_cond = huge(1.0d0)
        rank_def = n_mix
      end if
      swp_norm = sqrt(sum(abs(S_wp(1:n_w, 1:n_pw))**2))
      swp_norm_max = max(swp_norm_max, swp_norm)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a,i0,a,i0,5(a,1pe12.4),a,i0)') &
          '[DG-WPW-MIXED-S] self overlap block:', &
          ' ifrag=', ifrag, ' nw=', n_w, ' npw=', n_pw, &
          ' ||S_WP||=', swp_norm, ' S_herm=', s_herm, &
          ' S_eval_min=', s_min, ' S_eval_max=', s_max, &
          ' S_cond=', s_cond, ' S_rank_def=', rank_def
      end if
      deallocate(S_mix, S_wp, S_pp, wval, eval, rwork, work)
    end do
    if (dg_frag%id == 0) then
      write(*,'(1x,a,1pe12.4)') '[DG-WPW-MIXED-S] max ||S_WP|| over local fragments=', swp_norm_max
    end if
  end subroutine diagnose_wpw_mixed_self_overlap

  subroutine diagnose_wpw_mixed_self_hamiltonian(dg_frag, Vh, Vxc, Vpsl)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl

    integer :: i_local, ifrag, iw, ib, ipw, ix, iy, iz, gx, gy, gz, idir
    integer :: n_w, n_pw, n_mix, nbf, iblk, ispin, info, lwork, rank_def
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: v_lb1, v_lb2, v_lb3, v_ub1, v_ub2, v_ub3
    integer :: bx, by, bz, vx, vy, vz
    real(8) :: vol_elem, box_volume, chi, grad_chi(3), phase_arg, v_total, Lbox(3)
    real(8) :: h_herm, hwp_norm, hwp_norm_max, hwp_local_norm, twp_norm
    real(8) :: eval_shift_zero_max, eval_shift_zero_local, eval_shift_local_max, eval_shift_local
    real(8) :: eval_min, eval_max, eval0_min, eval0_max, evalloc_min, evalloc_max, s_min, s_max, s_cond
    logical :: have_grad_cache
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    complex(8) :: phase, grad_p(3)
    complex(8), allocatable :: S_mix(:,:), H_mix(:,:), H_zero(:,:), H_local(:,:), S_work(:,:), H_work(:,:)
    complex(8), allocatable :: S_wp(:,:), H_wp(:,:), H_wp_local(:,:), T_wp(:,:), S_pp(:,:), H_pp(:,:)
    complex(8), allocatable :: wval(:), grad_wval(:,:), work(:)
    real(8), allocatable :: eval(:), eval0(:), evalloc(:), rwork(:)
    external :: zhegv

    if (.not. dg_frag%has_global_wannier_local_basis) return
    if (.not. allocated(dg_frag%global_wannier_local_coef)) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return
    if (.not. allocated(dg_frag%wpw_S_pp_blocks)) return
    if (.not. allocated(dg_frag%wpw_T_pp_volume_blocks)) return
    if (.not. allocated(dg_frag%wpw_T_pp_interface_blocks)) return
    if (.not. allocated(dg_frag%phi_frag) .and. .not. allocated(dg_frag%phi_frag_c)) return

    n_pw = dg_frag%n_plane_waves
    if (n_pw <= 0) return
    vol_elem = product(dg_frag%hgs(1:3))
    box_volume = product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))
    if (allocated(dg_frag%phi_frag_c)) then
      p_lb1 = lbound(dg_frag%phi_frag_c, 1); p_ub1 = ubound(dg_frag%phi_frag_c, 1)
      p_lb2 = lbound(dg_frag%phi_frag_c, 2); p_ub2 = ubound(dg_frag%phi_frag_c, 2)
      p_lb3 = lbound(dg_frag%phi_frag_c, 3); p_ub3 = ubound(dg_frag%phi_frag_c, 3)
    else
      p_lb1 = lbound(dg_frag%phi_frag, 1); p_ub1 = ubound(dg_frag%phi_frag, 1)
      p_lb2 = lbound(dg_frag%phi_frag, 2); p_ub2 = ubound(dg_frag%phi_frag, 2)
      p_lb3 = lbound(dg_frag%phi_frag, 3); p_ub3 = ubound(dg_frag%phi_frag, 3)
    end if
    v_lb1 = lbound(Vpsl%f, 1); v_ub1 = ubound(Vpsl%f, 1)
    v_lb2 = lbound(Vpsl%f, 2); v_ub2 = ubound(Vpsl%f, 2)
    v_lb3 = lbound(Vpsl%f, 3); v_ub3 = ubound(Vpsl%f, 3)

    have_grad_cache = dg_frag%gradient_basis_cache_valid .and. allocated(dg_frag%gradient_basis_cache)
    hwp_norm_max = 0.0d0
    eval_shift_zero_max = 0.0d0
    eval_shift_local_max = 0.0d0
    do i_local = 1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1
      ifrag = dg_frag%ifrag_start + i_local - 1
      if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
      n_w = dg_frag%global_wannier_local_nkeep(i_local)
      if (n_w <= 0) cycle
      iblk = find_wpw_pp_block(dg_frag%wpw_S_pp_blocks, ifrag, ifrag)
      if (iblk <= 0) cycle
      n_mix = n_w + n_pw

      allocate(S_mix(n_mix, n_mix), H_mix(n_mix, n_mix), H_zero(n_mix, n_mix), H_local(n_mix, n_mix))
      allocate(S_work(n_mix, n_mix), H_work(n_mix, n_mix))
      allocate(S_wp(n_w, n_pw), H_wp(n_w, n_pw), H_wp_local(n_w, n_pw), T_wp(n_w, n_pw))
      allocate(S_pp(n_pw, n_pw), H_pp(n_pw, n_pw), wval(n_w), grad_wval(3, n_w))
      S_mix(:, :) = (0.0d0, 0.0d0)
      H_mix(:, :) = (0.0d0, 0.0d0)
      H_zero(:, :) = (0.0d0, 0.0d0)
      H_local(:, :) = (0.0d0, 0.0d0)
      S_wp(:, :) = (0.0d0, 0.0d0)
      H_wp(:, :) = (0.0d0, 0.0d0)
      H_wp_local(:, :) = (0.0d0, 0.0d0)
      T_wp(:, :) = (0.0d0, 0.0d0)
      S_pp(:, :) = dg_frag%wpw_S_pp_blocks(iblk)%val(1:n_pw, 1:n_pw, 1)
      H_pp(:, :) = dg_frag%wpw_T_pp_volume_blocks(iblk)%val(1:n_pw, 1:n_pw, 1) + &
                   dg_frag%wpw_T_pp_interface_blocks(iblk)%val(1:n_pw, 1:n_pw, 1)
      do iw = 1, n_w
        S_mix(iw, iw) = (1.0d0, 0.0d0)
      end do

      nbf = min(dg_frag%n_basis(ifrag, 1), size(dg_frag%global_wannier_local_coef, 1))
      if (allocated(dg_frag%phi_frag_c)) then
        nbf = min(nbf, size(dg_frag%phi_frag_c, 4))
      else
        nbf = min(nbf, size(dg_frag%phi_frag, 4))
      end if
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
        bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
        vz = map_global_to_phi_box_coord_pw(gz, v_lb3, v_ub3, dg_frag%lgnum_total(3))
        if (bz < p_lb3 .or. bz > p_ub3 .or. vz < v_lb3 .or. vz > v_ub3) cycle
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
          by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          vy = map_global_to_phi_box_coord_pw(gy, v_lb2, v_ub2, dg_frag%lgnum_total(2))
          if (by < p_lb2 .or. by > p_ub2 .or. vy < v_lb2 .or. vy > v_ub2) cycle
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
            bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            vx = map_global_to_phi_box_coord_pw(gx, v_lb1, v_ub1, dg_frag%lgnum_total(1))
            if (bx < p_lb1 .or. bx > p_ub1 .or. vx < v_lb1 .or. vx > v_ub1) cycle
            wval(:) = (0.0d0, 0.0d0)
            grad_wval(:, :) = (0.0d0, 0.0d0)
            do iw = 1, n_w
              do ib = 1, nbf
                if (allocated(dg_frag%phi_frag_c)) then
                  wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, 1, i_local) * &
                    dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                else
                  wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, 1, i_local) * &
                    cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                end if
                if (have_grad_cache) then
                  do idir = 1, 3
                    grad_wval(idir, iw) = grad_wval(idir, iw) + &
                      dg_frag%global_wannier_local_coef(ib, iw, 1, i_local) * &
                      cmplx(dg_frag%gradient_basis_cache(ix, iy, iz, idir, ib, i_local), 0.0d0, kind=8)
                  end do
                end if
              end do
            end do
            v_total = Vpsl%f(vx, vy, vz) + Vh%f(vx, vy, vz) + Vxc(1)%f(vx, vy, vz)
            call wpw_normalized_window_at_grid(dg_frag, ifrag, gx, gy, gz, chi, grad_chi)
            do ipw = 1, n_pw
              phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                          dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                          dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
              phase = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
              grad_p(1:3) = (cmplx(grad_chi(1:3), 0.0d0, kind=8) + zi * dg_frag%k_pw(1:3, ipw) * chi) * &
                exp(cmplx(0.0d0, phase_arg, kind=8))
              do iw = 1, n_w
                S_wp(iw, ipw) = S_wp(iw, ipw) + conjg(wval(iw)) * phase * vol_elem / box_volume
                H_wp_local(iw, ipw) = H_wp_local(iw, ipw) + conjg(wval(iw)) * v_total * phase * vol_elem / box_volume
                if (have_grad_cache) then
                  T_wp(iw, ipw) = T_wp(iw, ipw) + 0.5d0 * &
                    sum(conjg(grad_wval(1:3, iw)) * grad_p(1:3)) * vol_elem / box_volume
                end if
              end do
            end do
          end do
        end do
      end do
      H_wp(1:n_w, 1:n_pw) = H_wp_local(1:n_w, 1:n_pw) + T_wp(1:n_w, 1:n_pw)

      S_mix(1:n_w, n_w+1:n_mix) = S_wp(1:n_w, 1:n_pw)
      S_mix(n_w+1:n_mix, 1:n_w) = conjg(transpose(S_wp(1:n_w, 1:n_pw)))
      S_mix(n_w+1:n_mix, n_w+1:n_mix) = S_pp(1:n_pw, 1:n_pw)
      H_mix(1:n_w, n_w+1:n_mix) = H_wp(1:n_w, 1:n_pw)
      H_mix(n_w+1:n_mix, 1:n_w) = conjg(transpose(H_wp(1:n_w, 1:n_pw)))
      H_mix(n_w+1:n_mix, n_w+1:n_mix) = H_pp(1:n_pw, 1:n_pw)
      H_local(1:n_w, n_w+1:n_mix) = H_wp_local(1:n_w, 1:n_pw)
      H_local(n_w+1:n_mix, 1:n_w) = conjg(transpose(H_wp_local(1:n_w, 1:n_pw)))
      H_local(n_w+1:n_mix, n_w+1:n_mix) = H_pp(1:n_pw, 1:n_pw)
      H_zero(n_w+1:n_mix, n_w+1:n_mix) = H_pp(1:n_pw, 1:n_pw)
      call hermitize_matrix(S_mix, n_mix)
      call hermitize_matrix(H_mix, n_mix)
      call hermitize_matrix(H_local, n_mix)
      call hermitize_matrix(H_zero, n_mix)
      call wpw_local_herm_max(H_mix, n_mix, h_herm)
      hwp_norm = sqrt(sum(abs(H_wp(1:n_w, 1:n_pw))**2))
      hwp_local_norm = sqrt(sum(abs(H_wp_local(1:n_w, 1:n_pw))**2))
      twp_norm = sqrt(sum(abs(T_wp(1:n_w, 1:n_pw))**2))
      hwp_norm_max = max(hwp_norm_max, hwp_norm)

      allocate(eval(n_mix), eval0(n_mix), evalloc(n_mix), rwork(max(1, 3*n_mix - 2)), work(1))
      H_work(:, :) = H_mix(:, :)
      S_work(:, :) = S_mix(:, :)
      lwork = -1
      call ZHEGV(1, 'N', 'U', n_mix, H_work, n_mix, S_work, n_mix, eval, work, lwork, rwork, info)
      lwork = max(1, int(real(work(1), kind=8)))
      deallocate(work)
      allocate(work(lwork))
      H_work(:, :) = H_mix(:, :)
      S_work(:, :) = S_mix(:, :)
      call ZHEGV(1, 'N', 'U', n_mix, H_work, n_mix, S_work, n_mix, eval, work, lwork, rwork, info)
      if (info == 0) then
        eval_min = eval(1)
        eval_max = eval(n_mix)
      else
        eval_min = 0.0d0
        eval_max = 0.0d0
      end if
      H_work(:, :) = H_zero(:, :)
      S_work(:, :) = S_mix(:, :)
      call ZHEGV(1, 'N', 'U', n_mix, H_work, n_mix, S_work, n_mix, eval0, work, lwork, rwork, info)
      if (info == 0) then
        eval0_min = eval0(1)
        eval0_max = eval0(n_mix)
        eval_shift_zero_local = maxval(abs(eval(1:n_mix) - eval0(1:n_mix)))
      else
        eval0_min = 0.0d0
        eval0_max = 0.0d0
        eval_shift_zero_local = 0.0d0
      end if
      H_work(:, :) = H_local(:, :)
      S_work(:, :) = S_mix(:, :)
      call ZHEGV(1, 'N', 'U', n_mix, H_work, n_mix, S_work, n_mix, evalloc, work, lwork, rwork, info)
      if (info == 0) then
        evalloc_min = evalloc(1)
        evalloc_max = evalloc(n_mix)
        eval_shift_local = maxval(abs(eval(1:n_mix) - evalloc(1:n_mix)))
      else
        evalloc_min = 0.0d0
        evalloc_max = 0.0d0
        eval_shift_local = 0.0d0
      end if
      eval_shift_zero_max = max(eval_shift_zero_max, eval_shift_zero_local)
      eval_shift_local_max = max(eval_shift_local_max, eval_shift_local)

      deallocate(work)
      allocate(work(1))
      H_work(:, :) = S_mix(:, :)
      lwork = -1
      call ZHEEV('N', 'U', n_mix, H_work, n_mix, eval, work, lwork, rwork, info)
      lwork = max(1, int(real(work(1), kind=8)))
      deallocate(work)
      allocate(work(lwork))
      H_work(:, :) = S_mix(:, :)
      call ZHEEV('N', 'U', n_mix, H_work, n_mix, eval, work, lwork, rwork, info)
      if (info == 0) then
        s_min = eval(1)
        s_max = eval(n_mix)
        rank_def = count(eval(1:n_mix) < 1.0d-8)
        if (s_min > 0.0d0) then
          s_cond = s_max / s_min
        else
          s_cond = huge(1.0d0)
        end if
      else
        s_min = 0.0d0
        s_max = 0.0d0
        s_cond = huge(1.0d0)
        rank_def = n_mix
      end if

      if (dg_frag%id == 0) then
        write(*,'(1x,a,a,i0,a,i0,a,i0,10(a,1pe12.4),a,i0)') &
          '[DG-WPW-MIXED-H] self local+weak-kinetic WP block:', &
          ' ifrag=', ifrag, ' nw=', n_w, ' npw=', n_pw, &
          ' ||H_WP_local||=', hwp_local_norm, ' ||T_WP||=', twp_norm, &
          ' ||H_WP||=', hwp_norm, ' H_herm=', h_herm, &
          ' eval_min=', eval_min, ' eval_max=', eval_max, &
          ' evalloc_min=', evalloc_min, ' evalloc_max=', evalloc_max, &
          ' eval_shift_vs_local=', eval_shift_local, ' S_cond=', s_cond, ' S_rank_def=', rank_def
      end if
      deallocate(S_mix, H_mix, H_zero, H_local, S_work, H_work, S_wp, H_wp, H_wp_local, T_wp)
      deallocate(S_pp, H_pp, wval, grad_wval, eval, eval0, evalloc, rwork, work)
    end do

    if (dg_frag%id == 0) then
      write(*,'(1x,a,2(a,1pe12.4))') &
        '[DG-WPW-MIXED-H] max over local fragments:', &
        ' ||H_WP||=', hwp_norm_max, ' eval_shift_vs_zero_max=', eval_shift_zero_max
      write(*,'(1x,a,1pe12.4,a,l1)') &
        '[DG-WPW-MIXED-H] max eval_shift_vs_local=', eval_shift_local_max, &
        ' grad_cache=', have_grad_cache
    end if
  end subroutine diagnose_wpw_mixed_self_hamiltonian

  subroutine diagnose_wpw_mixed_neighbor_hamiltonian(dg_frag, Vh, Vxc, Vpsl, force_sorth, quiet)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    logical, intent(in), optional :: force_sorth, quiet

    integer :: i_local, ifrag, n_pfrag, pidx, qidx, iblk
    integer :: iw, jw, ib, ipw, ix, iy, iz, gx, gy, gz, idir, n_w, n_pw, n_mix, n_self
    integer :: nbf, info, lwork, rank_def
    integer :: n_ab, n_local_frag, wpw_red_max_dim
    integer :: n_red, n_keep, n_drop, info_red
    integer :: neighbor_shell
    integer :: max_g2, drop_g2, ipw_g2, ik_mode(3)
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: v_lb1, v_lb2, v_lb3, v_ub1, v_ub2, v_ub3
    integer :: bx, by, bz, vx, vy, vz, row0, col0
    real(8) :: vol_elem, box_volume, chi, grad_chi(3), phase_arg, v_total, Lbox(3)
    real(8) :: h_herm, s_herm, hwp_self_norm, hwp_neigh_norm, hpw_neigh_norm, eval_shift_self
    real(8) :: swp_self_norm, swp_neigh_norm, swp_ab_sv_min, swp_ab_sv_max
    real(8) :: hwp_neigh_perp_norm, swp_neigh_perp_norm, tgrad_neigh_norm, tig_neigh_norm
    real(8) :: hlocal_neigh_norm, hwp_perp_ratio, swp_perp_ratio
    real(8) :: tsym_neigh_norm, tdiff_sym_norm, tcurrent_neigh_norm
    real(8) :: tgrad_face_norm(3), gradchi_a_face_norm(3), gradchi_b_face_norm(3)
    real(8) :: hwp_neigh_filtered_norm, swp_neigh_filtered_norm, tgrad_filtered_norm, tig_filtered_norm
    real(8) :: eval_min, eval_max, eval_self_min, eval_self_max, s_min, s_max, s_cond
    real(8) :: eval_filt_min, eval_filt_max, eval_shift_filt_self, eval_shift_raw_filt
    real(8) :: neigh_norm_max, shift_max
    real(8) :: tol_red, lambda_min, lambda_max, s_sn_after, snn_i_err, sss_min, sss_max
    real(8) :: red_s_min, red_s_max, red_h_herm
    real(8) :: red_lambda_min_all, red_lambda_max_all, red_s_sn_max, red_snn_i_max, red_h_herm_max
    real(8) :: red_s_min_all, red_s_max_all
    logical :: have_grad_cache, do_sorth, quiet_mode
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    complex(8) :: phase, grad_p(3), grad_sym_w(3)
    complex(8), allocatable :: S_mix(:,:), H_mix(:,:), S_filt(:,:), H_filt(:,:)
    complex(8), allocatable :: S_self(:,:), H_self(:,:), S_work(:,:), H_work(:,:)
    complex(8), allocatable :: S_self_ref(:,:), H_self_ref(:,:)
    complex(8), allocatable :: H_red(:,:), S_red(:,:), T_red(:,:)
    complex(8), allocatable :: S_wp(:,:,:), H_wp(:,:,:), H_wp_local(:,:,:), T_wp(:,:,:)
    complex(8), allocatable :: T_wp_gradchi(:,:,:), T_wp_igchi(:,:,:), T_wp_sym(:,:,:), H_ww(:,:)
    complex(8), allocatable :: T_wp_gradchi_face(:,:,:)
    complex(8), allocatable :: S_ab(:,:), H_ab(:,:), H_ab_perp(:,:), S_ab_perp(:,:)
    complex(8), allocatable :: T_grad_ab(:,:), T_ig_ab(:,:), T_sym_ab(:,:), H_local_ab(:,:)
    integer, allocatable :: pfrag_ids(:), pfrag_axis(:), pfrag_side(:)
    real(8) :: chi_a, grad_chi_a(3)
    complex(8), allocatable :: wval(:), grad_wval(:,:), work(:)
    real(8), allocatable :: eval(:), eval_raw(:), eval_filt(:), eval_self(:), eval_red(:), rwork(:)
    external :: zhegv, zheev

    do_sorth = wpw_neighbor_sorth_candidate_enabled()
    if (present(force_sorth)) do_sorth = do_sorth .or. force_sorth
    quiet_mode = .false.
    if (present(quiet)) quiet_mode = quiet

    if (.not. dg_frag%has_global_wannier_local_basis) return
    if (.not. allocated(dg_frag%global_wannier_local_coef)) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return
    if (.not. allocated(dg_frag%wpw_S_pp_blocks)) return
    if (.not. allocated(dg_frag%wpw_T_pp_volume_blocks)) return
    if (.not. allocated(dg_frag%wpw_T_pp_interface_blocks)) return
    if (.not. allocated(dg_frag%phi_frag) .and. .not. allocated(dg_frag%phi_frag_c)) return

    n_pw = dg_frag%n_plane_waves
    if (n_pw <= 0) return
    have_grad_cache = dg_frag%gradient_basis_cache_valid .and. allocated(dg_frag%gradient_basis_cache)
    vol_elem = product(dg_frag%hgs(1:3))
    box_volume = product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))
    if (allocated(dg_frag%phi_frag_c)) then
      p_lb1 = lbound(dg_frag%phi_frag_c, 1); p_ub1 = ubound(dg_frag%phi_frag_c, 1)
      p_lb2 = lbound(dg_frag%phi_frag_c, 2); p_ub2 = ubound(dg_frag%phi_frag_c, 2)
      p_lb3 = lbound(dg_frag%phi_frag_c, 3); p_ub3 = ubound(dg_frag%phi_frag_c, 3)
    else
      p_lb1 = lbound(dg_frag%phi_frag, 1); p_ub1 = ubound(dg_frag%phi_frag, 1)
      p_lb2 = lbound(dg_frag%phi_frag, 2); p_ub2 = ubound(dg_frag%phi_frag, 2)
      p_lb3 = lbound(dg_frag%phi_frag, 3); p_ub3 = ubound(dg_frag%phi_frag, 3)
    end if
    v_lb1 = lbound(Vpsl%f, 1); v_ub1 = ubound(Vpsl%f, 1)
    v_lb2 = lbound(Vpsl%f, 2); v_ub2 = ubound(Vpsl%f, 2)
    v_lb3 = lbound(Vpsl%f, 3); v_ub3 = ubound(Vpsl%f, 3)

    neigh_norm_max = 0.0d0
    shift_max = 0.0d0
    red_lambda_min_all = huge(1.0d0)
    red_lambda_max_all = 0.0d0
    red_s_sn_max = 0.0d0
    red_snn_i_max = 0.0d0
    red_h_herm_max = 0.0d0
    red_s_min_all = huge(1.0d0)
    red_s_max_all = 0.0d0
    tol_red = wpw_neighbor_sorth_tol()
    n_local_frag = max(0, dg_frag%ifrag_end - dg_frag%ifrag_start + 1)
    neighbor_shell = wpw_neighbor_env_shell_from_backend()
    wpw_red_max_dim = 0
    do i_local = 1, n_local_frag
      if (i_local >= 1 .and. i_local <= size(dg_frag%global_wannier_local_nkeep)) then
        wpw_red_max_dim = max(wpw_red_max_dim, dg_frag%global_wannier_local_nkeep(i_local) + &
          max(1, product(dg_frag%num_fragment(1:3))) * n_pw)
      end if
    end do
    if (do_sorth .and. n_local_frag > 0 .and. wpw_red_max_dim > 0) then
      if (allocated(dg_frag%wpw_reduced_dim)) deallocate(dg_frag%wpw_reduced_dim)
      if (allocated(dg_frag%wpw_reduced_nself)) deallocate(dg_frag%wpw_reduced_nself)
      if (allocated(dg_frag%wpw_reduced_nkeep)) deallocate(dg_frag%wpw_reduced_nkeep)
      if (allocated(dg_frag%wpw_reduced_ndrop)) deallocate(dg_frag%wpw_reduced_ndrop)
      if (allocated(dg_frag%wpw_reduced_H)) deallocate(dg_frag%wpw_reduced_H)
      if (allocated(dg_frag%wpw_reduced_S)) deallocate(dg_frag%wpw_reduced_S)
      if (allocated(dg_frag%wpw_reduced_transform)) deallocate(dg_frag%wpw_reduced_transform)
      if (allocated(dg_frag%wpw_reduced_Hraw_build)) deallocate(dg_frag%wpw_reduced_Hraw_build)
      if (allocated(dg_frag%wpw_reduced_Sraw_build)) deallocate(dg_frag%wpw_reduced_Sraw_build)
      if (allocated(dg_frag%wpw_reduced_nraw)) deallocate(dg_frag%wpw_reduced_nraw)
      if (allocated(dg_frag%wpw_reproject_prev_coef)) deallocate(dg_frag%wpw_reproject_prev_coef)
      allocate(dg_frag%wpw_reduced_dim(n_local_frag), dg_frag%wpw_reduced_nself(n_local_frag))
      allocate(dg_frag%wpw_reduced_nkeep(n_local_frag), dg_frag%wpw_reduced_ndrop(n_local_frag))
      allocate(dg_frag%wpw_reduced_nraw(n_local_frag))
      allocate(dg_frag%wpw_reduced_H(wpw_red_max_dim, wpw_red_max_dim, dg_frag%nspin, n_local_frag))
      allocate(dg_frag%wpw_reduced_S(wpw_red_max_dim, wpw_red_max_dim, dg_frag%nspin, n_local_frag))
      allocate(dg_frag%wpw_reduced_transform(wpw_red_max_dim, wpw_red_max_dim, n_local_frag))
      allocate(dg_frag%wpw_reduced_Hraw_build(wpw_red_max_dim, wpw_red_max_dim, n_local_frag))
      allocate(dg_frag%wpw_reduced_Sraw_build(wpw_red_max_dim, wpw_red_max_dim, n_local_frag))
      dg_frag%wpw_reduced_dim(:) = 0
      dg_frag%wpw_reduced_nself(:) = 0
      dg_frag%wpw_reduced_nkeep(:) = 0
      dg_frag%wpw_reduced_ndrop(:) = 0
      dg_frag%wpw_reduced_nraw(:) = 0
      dg_frag%wpw_reduced_H(:, :, :, :) = (0.0d0, 0.0d0)
      dg_frag%wpw_reduced_S(:, :, :, :) = (0.0d0, 0.0d0)
      dg_frag%wpw_reduced_transform(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%wpw_reduced_Hraw_build(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%wpw_reduced_Sraw_build(:, :, :) = (0.0d0, 0.0d0)
      dg_frag%wpw_reduced_ready = .false.
      dg_frag%wpw_reduced_shell = 0
      dg_frag%wpw_reduced_max_dim = wpw_red_max_dim
      dg_frag%wpw_reproject_prev_valid = .false.
    end if
    do i_local = 1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1
      ifrag = dg_frag%ifrag_start + i_local - 1
      if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
      n_w = dg_frag%global_wannier_local_nkeep(i_local)
      if (n_w <= 0) cycle

      call collect_wpw_fragment_shell_local(dg_frag, ifrag, neighbor_shell, pfrag_ids, pfrag_axis, pfrag_side, n_pfrag)
      if (n_pfrag <= 1) then
        if (allocated(pfrag_ids)) deallocate(pfrag_ids)
        if (allocated(pfrag_axis)) deallocate(pfrag_axis)
        if (allocated(pfrag_side)) deallocate(pfrag_side)
        cycle
      end if

      n_mix = n_w + n_pfrag * n_pw
      n_self = n_w + n_pw
      allocate(S_mix(n_mix, n_mix), H_mix(n_mix, n_mix), S_filt(n_mix, n_mix), H_filt(n_mix, n_mix))
      allocate(S_work(n_mix, n_mix), H_work(n_mix, n_mix))
      allocate(S_self(n_self, n_self), H_self(n_self, n_self), S_self_ref(n_self, n_self), H_self_ref(n_self, n_self))
      allocate(S_wp(n_w, n_pw, n_pfrag), H_wp(n_w, n_pw, n_pfrag))
      allocate(H_wp_local(n_w, n_pw, n_pfrag), T_wp(n_w, n_pw, n_pfrag))
      allocate(T_wp_gradchi(n_w, n_pw, n_pfrag), T_wp_igchi(n_w, n_pw, n_pfrag))
      allocate(T_wp_sym(n_w, n_pw, n_pfrag), T_wp_gradchi_face(n_w, n_pw, 3), H_ww(n_w, n_w))
      allocate(wval(n_w), grad_wval(3, n_w))
      S_mix(:, :) = (0.0d0, 0.0d0)
      H_mix(:, :) = (0.0d0, 0.0d0)
      S_self(:, :) = (0.0d0, 0.0d0)
      H_self(:, :) = (0.0d0, 0.0d0)
      S_wp(:, :, :) = (0.0d0, 0.0d0)
      H_wp(:, :, :) = (0.0d0, 0.0d0)
      H_wp_local(:, :, :) = (0.0d0, 0.0d0)
      T_wp(:, :, :) = (0.0d0, 0.0d0)
      T_wp_gradchi(:, :, :) = (0.0d0, 0.0d0)
      T_wp_igchi(:, :, :) = (0.0d0, 0.0d0)
      T_wp_sym(:, :, :) = (0.0d0, 0.0d0)
      T_wp_gradchi_face(:, :, :) = (0.0d0, 0.0d0)
      H_ww(:, :) = (0.0d0, 0.0d0)
      gradchi_a_face_norm(:) = 0.0d0
      gradchi_b_face_norm(:) = 0.0d0
      do iw = 1, n_w
        S_mix(iw, iw) = (1.0d0, 0.0d0)
        S_self(iw, iw) = (1.0d0, 0.0d0)
      end do

      do pidx = 1, n_pfrag
        row0 = n_w + (pidx - 1) * n_pw
        do qidx = 1, n_pfrag
          col0 = n_w + (qidx - 1) * n_pw
          iblk = find_wpw_pp_block(dg_frag%wpw_S_pp_blocks, pfrag_ids(pidx), pfrag_ids(qidx))
          if (iblk <= 0) cycle
          S_mix(row0+1:row0+n_pw, col0+1:col0+n_pw) = dg_frag%wpw_S_pp_blocks(iblk)%val(1:n_pw, 1:n_pw, 1)
          H_mix(row0+1:row0+n_pw, col0+1:col0+n_pw) = &
            dg_frag%wpw_T_pp_volume_blocks(iblk)%val(1:n_pw, 1:n_pw, 1) + &
            dg_frag%wpw_T_pp_interface_blocks(iblk)%val(1:n_pw, 1:n_pw, 1)
          if (pidx == 1 .and. qidx == 1) then
            S_self(n_w+1:n_self, n_w+1:n_self) = S_mix(row0+1:row0+n_pw, col0+1:col0+n_pw)
            H_self(n_w+1:n_self, n_w+1:n_self) = H_mix(row0+1:row0+n_pw, col0+1:col0+n_pw)
          end if
        end do
      end do

      nbf = min(dg_frag%n_basis(ifrag, 1), size(dg_frag%global_wannier_local_coef, 1))
      if (allocated(dg_frag%phi_frag_c)) then
        nbf = min(nbf, size(dg_frag%phi_frag_c, 4))
      else
        nbf = min(nbf, size(dg_frag%phi_frag, 4))
      end if
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
        bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
        vz = map_global_to_phi_box_coord_pw(gz, v_lb3, v_ub3, dg_frag%lgnum_total(3))
        if (bz < p_lb3 .or. bz > p_ub3 .or. vz < v_lb3 .or. vz > v_ub3) cycle
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
          by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
          vy = map_global_to_phi_box_coord_pw(gy, v_lb2, v_ub2, dg_frag%lgnum_total(2))
          if (by < p_lb2 .or. by > p_ub2 .or. vy < v_lb2 .or. vy > v_ub2) cycle
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
            bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
            vx = map_global_to_phi_box_coord_pw(gx, v_lb1, v_ub1, dg_frag%lgnum_total(1))
            if (bx < p_lb1 .or. bx > p_ub1 .or. vx < v_lb1 .or. vx > v_ub1) cycle
            wval(:) = (0.0d0, 0.0d0)
            grad_wval(:, :) = (0.0d0, 0.0d0)
            do iw = 1, n_w
              do ib = 1, nbf
                if (allocated(dg_frag%phi_frag_c)) then
                  wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, 1, i_local) * &
                    dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                else
                  wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, 1, i_local) * &
                    cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                end if
                if (have_grad_cache) then
                  do idir = 1, 3
                    grad_wval(idir, iw) = grad_wval(idir, iw) + &
                      dg_frag%global_wannier_local_coef(ib, iw, 1, i_local) * &
                      cmplx(dg_frag%gradient_basis_cache(ix, iy, iz, idir, ib, i_local), 0.0d0, kind=8)
                  end do
                end if
              end do
            end do
            v_total = Vpsl%f(vx, vy, vz) + Vh%f(vx, vy, vz) + Vxc(1)%f(vx, vy, vz)
            call wpw_normalized_window_at_grid(dg_frag, ifrag, gx, gy, gz, chi_a, grad_chi_a)
            do jw = 1, n_w
              do iw = 1, n_w
                H_ww(iw, jw) = H_ww(iw, jw) + conjg(wval(iw)) * v_total * wval(jw) * vol_elem / box_volume
                if (have_grad_cache) then
                  H_ww(iw, jw) = H_ww(iw, jw) + 0.5d0 * &
                    sum(conjg(grad_wval(1:3, iw)) * grad_wval(1:3, jw)) * vol_elem / box_volume
                end if
              end do
            end do
            do pidx = 1, n_pfrag
              call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
              if (pidx > 1 .and. pfrag_axis(pidx) >= 1 .and. pfrag_axis(pidx) <= 3) then
                gradchi_a_face_norm(pfrag_axis(pidx)) = gradchi_a_face_norm(pfrag_axis(pidx)) + &
                  sum(grad_chi_a(1:3)**2) * vol_elem / box_volume
                gradchi_b_face_norm(pfrag_axis(pidx)) = gradchi_b_face_norm(pfrag_axis(pidx)) + &
                  sum(grad_chi(1:3)**2) * vol_elem / box_volume
              end if
              do ipw = 1, n_pw
                phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                            dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                            dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                phase = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                grad_p(1:3) = (cmplx(grad_chi(1:3), 0.0d0, kind=8) + zi * dg_frag%k_pw(1:3, ipw) * chi) * &
                  exp(cmplx(0.0d0, phase_arg, kind=8))
                do iw = 1, n_w
                  S_wp(iw, ipw, pidx) = S_wp(iw, ipw, pidx) + conjg(wval(iw)) * phase * vol_elem / box_volume
                  H_wp_local(iw, ipw, pidx) = H_wp_local(iw, ipw, pidx) + &
                    conjg(wval(iw)) * v_total * phase * vol_elem / box_volume
                  if (have_grad_cache) then
                    grad_sym_w(1:3) = chi_a * grad_wval(1:3, iw) + &
                      cmplx(grad_chi_a(1:3), 0.0d0, kind=8) * wval(iw)
                    T_wp_sym(iw, ipw, pidx) = T_wp_sym(iw, ipw, pidx) + 0.5d0 * &
                      sum(conjg(grad_sym_w(1:3)) * grad_p(1:3)) * vol_elem / box_volume
                    T_wp_gradchi(iw, ipw, pidx) = T_wp_gradchi(iw, ipw, pidx) + 0.5d0 * &
                      sum(conjg(grad_wval(1:3, iw)) * cmplx(grad_chi(1:3), 0.0d0, kind=8) * &
                      exp(cmplx(0.0d0, phase_arg, kind=8))) * vol_elem / box_volume
                    T_wp_igchi(iw, ipw, pidx) = T_wp_igchi(iw, ipw, pidx) + 0.5d0 * &
                      sum(conjg(grad_wval(1:3, iw)) * zi * dg_frag%k_pw(1:3, ipw) * chi * &
                      exp(cmplx(0.0d0, phase_arg, kind=8))) * vol_elem / box_volume
                    T_wp(iw, ipw, pidx) = T_wp(iw, ipw, pidx) + 0.5d0 * &
                      sum(conjg(grad_wval(1:3, iw)) * grad_p(1:3)) * vol_elem / box_volume
                    if (pidx > 1 .and. pfrag_axis(pidx) >= 1 .and. pfrag_axis(pidx) <= 3) then
                      T_wp_gradchi_face(iw, ipw, pfrag_axis(pidx)) = &
                        T_wp_gradchi_face(iw, ipw, pfrag_axis(pidx)) + 0.5d0 * &
                        sum(conjg(grad_wval(1:3, iw)) * cmplx(grad_chi(1:3), 0.0d0, kind=8) * &
                        exp(cmplx(0.0d0, phase_arg, kind=8))) * vol_elem / box_volume
                    end if
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
      H_wp(:, :, :) = H_wp_local(:, :, :) + T_wp(:, :, :)
      call hermitize_matrix(H_ww, n_w)

      do pidx = 1, n_pfrag
        row0 = n_w + (pidx - 1) * n_pw
        S_mix(1:n_w, row0+1:row0+n_pw) = S_wp(1:n_w, 1:n_pw, pidx)
        S_mix(row0+1:row0+n_pw, 1:n_w) = conjg(transpose(S_wp(1:n_w, 1:n_pw, pidx)))
        H_mix(1:n_w, row0+1:row0+n_pw) = H_wp(1:n_w, 1:n_pw, pidx)
        H_mix(row0+1:row0+n_pw, 1:n_w) = conjg(transpose(H_wp(1:n_w, 1:n_pw, pidx)))
        if (pidx == 1) then
          S_self(1:n_w, n_w+1:n_self) = S_wp(1:n_w, 1:n_pw, pidx)
          S_self(n_w+1:n_self, 1:n_w) = conjg(transpose(S_wp(1:n_w, 1:n_pw, pidx)))
          H_self(1:n_w, n_w+1:n_self) = H_wp(1:n_w, 1:n_pw, pidx)
          H_self(n_w+1:n_self, 1:n_w) = conjg(transpose(H_wp(1:n_w, 1:n_pw, pidx)))
        end if
      end do

      call hermitize_matrix(S_mix, n_mix)
      call hermitize_matrix(H_mix, n_mix)
      call hermitize_matrix(S_self, n_self)
      call hermitize_matrix(H_self, n_self)
      S_self_ref(:, :) = S_self(:, :)
      H_self_ref(:, :) = H_self(:, :)
      call wpw_local_herm_max(H_mix, n_mix, h_herm)
      call wpw_local_herm_max(S_mix, n_mix, s_herm)
      hwp_self_norm = sqrt(sum(abs(H_wp(1:n_w, 1:n_pw, 1))**2))
      hwp_neigh_norm = sqrt(sum(abs(H_wp(1:n_w, 1:n_pw, 2:n_pfrag))**2))
      hpw_neigh_norm = hwp_neigh_norm
      swp_self_norm = sqrt(sum(abs(S_wp(1:n_w, 1:n_pw, 1))**2))
      swp_neigh_norm = sqrt(sum(abs(S_wp(1:n_w, 1:n_pw, 2:n_pfrag))**2))
      hlocal_neigh_norm = sqrt(sum(abs(H_wp_local(1:n_w, 1:n_pw, 2:n_pfrag))**2))
      tgrad_neigh_norm = sqrt(sum(abs(T_wp_gradchi(1:n_w, 1:n_pw, 2:n_pfrag))**2))
      tig_neigh_norm = sqrt(sum(abs(T_wp_igchi(1:n_w, 1:n_pw, 2:n_pfrag))**2))
      max_g2 = 0
      do ipw = 1, n_pw
        ik_mode(1) = nint(dg_frag%k_pw(1, ipw) * Lbox(1) / (2.0d0 * acos(-1.0d0)))
        ik_mode(2) = nint(dg_frag%k_pw(2, ipw) * Lbox(2) / (2.0d0 * acos(-1.0d0)))
        ik_mode(3) = nint(dg_frag%k_pw(3, ipw) * Lbox(3) / (2.0d0 * acos(-1.0d0)))
        max_g2 = max(max_g2, sum(ik_mode(1:3)**2))
      end do
      drop_g2 = wpw_neighbor_filter_drop_g2_override()
      if (drop_g2 == huge(1)) drop_g2 = max_g2
      hwp_neigh_filtered_norm = 0.0d0
      swp_neigh_filtered_norm = 0.0d0
      tgrad_filtered_norm = 0.0d0
      tig_filtered_norm = 0.0d0
      do pidx = 2, n_pfrag
        do ipw = 1, n_pw
          ik_mode(1) = nint(dg_frag%k_pw(1, ipw) * Lbox(1) / (2.0d0 * acos(-1.0d0)))
          ik_mode(2) = nint(dg_frag%k_pw(2, ipw) * Lbox(2) / (2.0d0 * acos(-1.0d0)))
          ik_mode(3) = nint(dg_frag%k_pw(3, ipw) * Lbox(3) / (2.0d0 * acos(-1.0d0)))
          ipw_g2 = sum(ik_mode(1:3)**2)
          if (ipw_g2 >= drop_g2) cycle
          hwp_neigh_filtered_norm = hwp_neigh_filtered_norm + sum(abs(H_wp(1:n_w, ipw, pidx))**2)
          swp_neigh_filtered_norm = swp_neigh_filtered_norm + sum(abs(S_wp(1:n_w, ipw, pidx))**2)
          tgrad_filtered_norm = tgrad_filtered_norm + sum(abs(T_wp_gradchi(1:n_w, ipw, pidx))**2)
          tig_filtered_norm = tig_filtered_norm + sum(abs(T_wp_igchi(1:n_w, ipw, pidx))**2)
        end do
      end do
      hwp_neigh_filtered_norm = sqrt(hwp_neigh_filtered_norm)
      swp_neigh_filtered_norm = sqrt(swp_neigh_filtered_norm)
      tgrad_filtered_norm = sqrt(tgrad_filtered_norm)
      tig_filtered_norm = sqrt(tig_filtered_norm)
      n_ab = (n_pfrag - 1) * n_pw
      allocate(S_ab(n_w, n_ab), H_ab(n_w, n_ab), H_ab_perp(n_w, n_ab), S_ab_perp(n_w, n_ab))
      allocate(T_grad_ab(n_w, n_ab), T_ig_ab(n_w, n_ab), T_sym_ab(n_w, n_ab), H_local_ab(n_w, n_ab))
      do pidx = 2, n_pfrag
        col0 = (pidx - 2) * n_pw
        S_ab(1:n_w, col0+1:col0+n_pw) = S_wp(1:n_w, 1:n_pw, pidx)
        H_ab(1:n_w, col0+1:col0+n_pw) = H_wp(1:n_w, 1:n_pw, pidx)
        H_local_ab(1:n_w, col0+1:col0+n_pw) = H_wp_local(1:n_w, 1:n_pw, pidx)
        T_grad_ab(1:n_w, col0+1:col0+n_pw) = T_wp_gradchi(1:n_w, 1:n_pw, pidx)
        T_ig_ab(1:n_w, col0+1:col0+n_pw) = T_wp_igchi(1:n_w, 1:n_pw, pidx)
        T_sym_ab(1:n_w, col0+1:col0+n_pw) = T_wp_sym(1:n_w, 1:n_pw, pidx)
      end do
      S_ab_perp(:, :) = (0.0d0, 0.0d0)
      H_ab_perp(:, :) = H_ab(:, :) - matmul(H_ww, S_ab)
      swp_neigh_perp_norm = sqrt(sum(abs(S_ab_perp)**2))
      hwp_neigh_perp_norm = sqrt(sum(abs(H_ab_perp)**2))
      call rectangular_singular_minmax(S_ab, n_w, n_ab, swp_ab_sv_min, swp_ab_sv_max)
      tcurrent_neigh_norm = sqrt(sum(abs(T_grad_ab + T_ig_ab)**2))
      tsym_neigh_norm = sqrt(sum(abs(T_sym_ab)**2))
      tdiff_sym_norm = sqrt(sum(abs((T_grad_ab + T_ig_ab) - T_sym_ab)**2))
      do idir = 1, 3
        tgrad_face_norm(idir) = sqrt(sum(abs(T_wp_gradchi_face(1:n_w, 1:n_pw, idir))**2))
        gradchi_a_face_norm(idir) = sqrt(max(0.0d0, gradchi_a_face_norm(idir)))
        gradchi_b_face_norm(idir) = sqrt(max(0.0d0, gradchi_b_face_norm(idir)))
      end do
      hwp_perp_ratio = hwp_neigh_perp_norm / max(hwp_neigh_norm, 1.0d-300)
      swp_perp_ratio = swp_neigh_perp_norm / max(swp_neigh_norm, 1.0d-300)
      neigh_norm_max = max(neigh_norm_max, hwp_neigh_norm)

      S_filt(:, :) = S_mix(:, :)
      H_filt(:, :) = H_mix(:, :)
      do pidx = 2, n_pfrag
        row0 = n_w + (pidx - 1) * n_pw
        do ipw = 1, n_pw
          ik_mode(1) = nint(dg_frag%k_pw(1, ipw) * Lbox(1) / (2.0d0 * acos(-1.0d0)))
          ik_mode(2) = nint(dg_frag%k_pw(2, ipw) * Lbox(2) / (2.0d0 * acos(-1.0d0)))
          ik_mode(3) = nint(dg_frag%k_pw(3, ipw) * Lbox(3) / (2.0d0 * acos(-1.0d0)))
          ipw_g2 = sum(ik_mode(1:3)**2)
          if (ipw_g2 < drop_g2) cycle
          S_filt(1:n_w, row0+ipw) = (0.0d0, 0.0d0)
          S_filt(row0+ipw, 1:n_w) = (0.0d0, 0.0d0)
          H_filt(1:n_w, row0+ipw) = (0.0d0, 0.0d0)
          H_filt(row0+ipw, 1:n_w) = (0.0d0, 0.0d0)
        end do
      end do
      call hermitize_matrix(S_filt, n_mix)
      call hermitize_matrix(H_filt, n_mix)

      if (dg_frag%id == 0 .and. .not. quiet_mode) then
        call diagnose_self_like_eigen_tracking('raw', ifrag, H_self_ref, S_self_ref, H_mix, S_mix, n_self, n_mix)
        call diagnose_self_like_eigen_tracking('filtered', ifrag, H_self_ref, S_self_ref, H_filt, S_filt, n_self, n_mix)
        call diagnose_overlap_schur_effects('raw', ifrag, H_self_ref, S_self_ref, H_mix, S_mix, n_self, n_mix)
        if (do_sorth) then
          call build_s_orthogonal_neighbor_variant('sorth_neighbor', ifrag, H_mix, S_mix, n_self, n_mix, H_work, S_work)
          call diagnose_self_like_eigen_tracking('sorth_neighbor', ifrag, H_self_ref, S_self_ref, H_work, S_work, n_self, n_mix)
          call diagnose_overlap_schur_effects('sorth_neighbor', ifrag, H_self_ref, S_self_ref, H_work, S_work, n_self, n_mix)
          call diagnose_s_orthogonal_reduced_neighbor_variant('sorth_neighbor_reduced', ifrag, &
            H_self_ref, S_self_ref, H_mix, S_mix, n_self, n_mix)
        end if

        S_work(:, :) = S_mix(:, :)
        H_work(:, :) = (0.0d0, 0.0d0)
        H_work(1:n_self, 1:n_self) = H_self_ref(1:n_self, 1:n_self)
        if (n_mix > n_self) then
          H_work(n_self+1:n_mix, n_self+1:n_mix) = H_mix(n_self+1:n_mix, n_self+1:n_mix)
        end if
        call hermitize_matrix(H_work, n_mix)
        call diagnose_self_like_eigen_tracking('overlap_only', ifrag, H_self_ref, S_self_ref, H_work, S_work, n_self, n_mix)
        call diagnose_overlap_schur_effects('overlap_only', ifrag, H_self_ref, S_self_ref, H_work, S_work, n_self, n_mix)

        S_work(:, :) = S_mix(:, :)
        H_work(:, :) = H_mix(:, :)
        H_work(1:n_w, n_w+1:n_mix) = (0.0d0, 0.0d0)
        H_work(n_w+1:n_mix, 1:n_w) = (0.0d0, 0.0d0)
        call hermitize_matrix(H_work, n_mix)
        call diagnose_self_like_eigen_tracking('h_pp_only', ifrag, H_self_ref, S_self_ref, H_work, S_work, n_self, n_mix)

        S_work(:, :) = S_mix(:, :)
        if (n_mix > n_self) then
          S_work(1:n_self, n_self+1:n_mix) = (0.0d0, 0.0d0)
          S_work(n_self+1:n_mix, 1:n_self) = (0.0d0, 0.0d0)
        end if
        call hermitize_matrix(S_work, n_mix)
        H_work(:, :) = H_mix(:, :)
        call diagnose_self_like_eigen_tracking('h_coupling_only', ifrag, H_self_ref, S_self_ref, H_work, S_work, n_self, n_mix)
      end if

      if (do_sorth .and. allocated(dg_frag%wpw_reduced_H)) then
        call build_wpw_sorth_reduced_neighbor_block(H_mix, S_mix, n_self, n_mix, tol_red, H_red, S_red, &
          n_red, n_keep, n_drop, lambda_min, lambda_max, s_sn_after, snn_i_err, sss_min, sss_max, info_red, T_red)
        if (info_red == 0 .and. n_red <= dg_frag%wpw_reduced_max_dim) then
          dg_frag%wpw_reduced_dim(i_local) = n_red
          dg_frag%wpw_reduced_nself(i_local) = n_self
          dg_frag%wpw_reduced_nkeep(i_local) = n_keep
          dg_frag%wpw_reduced_ndrop(i_local) = n_drop
          dg_frag%wpw_reduced_nraw(i_local) = n_mix
          if (allocated(dg_frag%wpw_reduced_Hraw_build)) &
            dg_frag%wpw_reduced_Hraw_build(1:n_mix, 1:n_mix, i_local) = H_mix(1:n_mix, 1:n_mix)
          if (allocated(dg_frag%wpw_reduced_Sraw_build)) &
            dg_frag%wpw_reduced_Sraw_build(1:n_mix, 1:n_mix, i_local) = S_mix(1:n_mix, 1:n_mix)
          dg_frag%wpw_reduced_H(1:n_red, 1:n_red, 1:dg_frag%nspin, i_local) = &
            spread(H_red(1:n_red, 1:n_red), 3, dg_frag%nspin)
          dg_frag%wpw_reduced_S(1:n_red, 1:n_red, 1:dg_frag%nspin, i_local) = &
            spread(S_red(1:n_red, 1:n_red), 3, dg_frag%nspin)
          if (allocated(T_red)) dg_frag%wpw_reduced_transform(1:n_mix, 1:n_red, i_local) = T_red(1:n_mix, 1:n_red)
          call wpw_local_herm_max(H_red, n_red, red_h_herm)
          allocate(eval_red(n_red), work(1), rwork(max(1, 3*n_red - 2)))
          H_work(1:n_red, 1:n_red) = S_red(1:n_red, 1:n_red)
          lwork = -1
          call ZHEEV('N', 'U', n_red, H_work, n_red, eval_red, work, lwork, rwork, info)
          lwork = max(1, int(real(work(1), kind=8)))
          deallocate(work)
          allocate(work(lwork))
          H_work(1:n_red, 1:n_red) = S_red(1:n_red, 1:n_red)
          call ZHEEV('N', 'U', n_red, H_work, n_red, eval_red, work, lwork, rwork, info)
          if (info == 0) then
            red_s_min = eval_red(1)
            red_s_max = eval_red(n_red)
          else
            red_s_min = 0.0d0
            red_s_max = 0.0d0
          end if
          deallocate(eval_red, work, rwork)
          red_lambda_min_all = min(red_lambda_min_all, lambda_min)
          red_lambda_max_all = max(red_lambda_max_all, lambda_max)
          red_s_sn_max = max(red_s_sn_max, s_sn_after)
          red_snn_i_max = max(red_snn_i_max, snn_i_err)
          red_h_herm_max = max(red_h_herm_max, red_h_herm)
          red_s_min_all = min(red_s_min_all, red_s_min)
          red_s_max_all = max(red_s_max_all, red_s_max)
          if (red_s_min <= 0.0d0 .or. s_sn_after > 1.0d-8 .or. snn_i_err > 1.0d-8) then
            dg_frag%wpw_reduced_dim(i_local) = 0
          end if
          if (dg_frag%id == 0 .and. .not. quiet_mode) then
            write(*,'(1x,a,a,i0,5(a,i0),8(a,1pe12.4),a,i0)') &
              '[DG-WPW-RED-INIT] fragment reduced H/S:', &
              ' ifrag=', ifrag, ' dim=', n_red, ' nself=', n_self, &
              ' nkeep=', n_keep, ' ndrop=', n_drop, ' info=', info_red, &
              ' ||S_SN_after||=', s_sn_after, ' ||S_NN-I||=', snn_i_err, &
              ' lambda_keep_min=', lambda_min, ' lambda_keep_max=', lambda_max, &
              ' H_herm=', red_h_herm, ' S_eval_min=', red_s_min, &
              ' S_eval_max=', red_s_max, ' SSS_min=', sss_min
          end if
        else
          if (dg_frag%id == 0 .and. .not. quiet_mode) then
            write(*,'(1x,a,a,i0,4(a,i0),a,1pe12.4)') &
              '[DG-WPW-RED-INIT][WARN] reduced builder failed:', &
              ' ifrag=', ifrag, ' nself=', n_self, ' next=', n_mix, &
              ' nred=', n_red, ' info=', info_red, ' tol=', tol_red
          end if
        end if
        if (allocated(H_red)) deallocate(H_red)
        if (allocated(S_red)) deallocate(S_red)
        if (allocated(T_red)) deallocate(T_red)
      end if

      allocate(eval(n_mix), eval_raw(n_mix), eval_filt(n_mix), eval_self(n_self), rwork(max(1, 3*n_mix - 2)), work(1))
      H_work(:, :) = H_mix(:, :)
      S_work(:, :) = S_mix(:, :)
      lwork = -1
      call ZHEGV(1, 'N', 'U', n_mix, H_work, n_mix, S_work, n_mix, eval, work, lwork, rwork, info)
      lwork = max(1, int(real(work(1), kind=8)))
      deallocate(work)
      allocate(work(lwork))
      H_work(:, :) = H_mix(:, :)
      S_work(:, :) = S_mix(:, :)
      call ZHEGV(1, 'N', 'U', n_mix, H_work, n_mix, S_work, n_mix, eval, work, lwork, rwork, info)
      if (info == 0) then
        eval_min = eval(1)
        eval_max = eval(n_mix)
        eval_raw(1:n_mix) = eval(1:n_mix)
      else
        eval_min = 0.0d0
        eval_max = 0.0d0
        eval_raw(1:n_mix) = 0.0d0
      end if
      H_work(:, :) = H_filt(:, :)
      S_work(:, :) = S_filt(:, :)
      call ZHEGV(1, 'N', 'U', n_mix, H_work, n_mix, S_work, n_mix, eval, work, lwork, rwork, info)
      if (info == 0) then
        eval_filt_min = eval(1)
        eval_filt_max = eval(n_mix)
        eval_filt(1:n_mix) = eval(1:n_mix)
        eval_shift_raw_filt = maxval(abs(eval_raw(1:n_mix) - eval_filt(1:n_mix)))
      else
        eval_filt_min = 0.0d0
        eval_filt_max = 0.0d0
        eval_filt(1:n_mix) = 0.0d0
        eval_shift_raw_filt = 0.0d0
      end if
      deallocate(work)
      allocate(work(1))
      lwork = -1
      call ZHEGV(1, 'N', 'U', n_self, H_self, n_self, S_self, n_self, eval_self, work, lwork, rwork, info)
      lwork = max(1, int(real(work(1), kind=8)))
      deallocate(work)
      allocate(work(lwork))
      H_self(:, :) = H_self_ref(:, :)
      S_self(:, :) = S_self_ref(:, :)
      call ZHEGV(1, 'N', 'U', n_self, H_self, n_self, S_self, n_self, eval_self, work, lwork, rwork, info)
      if (info == 0) then
        eval_self_min = eval_self(1)
        eval_self_max = eval_self(n_self)
        eval_shift_self = maxval(abs(eval_raw(1:n_self) - eval_self(1:n_self)))
        eval_shift_filt_self = maxval(abs(eval_filt(1:n_self) - eval_self(1:n_self)))
      else
        eval_self_min = 0.0d0
        eval_self_max = 0.0d0
        eval_shift_self = 0.0d0
        eval_shift_filt_self = 0.0d0
      end if
      shift_max = max(shift_max, eval_shift_self)

      deallocate(work)
      allocate(work(1))
      H_work(:, :) = S_mix(:, :)
      lwork = -1
      call ZHEEV('N', 'U', n_mix, H_work, n_mix, eval, work, lwork, rwork, info)
      lwork = max(1, int(real(work(1), kind=8)))
      deallocate(work)
      allocate(work(lwork))
      H_work(:, :) = S_mix(:, :)
      call ZHEEV('N', 'U', n_mix, H_work, n_mix, eval, work, lwork, rwork, info)
      if (info == 0) then
        s_min = eval(1)
        s_max = eval(n_mix)
        rank_def = count(eval(1:n_mix) < 1.0d-8)
        if (s_min > 0.0d0) then
          s_cond = s_max / s_min
        else
          s_cond = huge(1.0d0)
        end if
      else
        s_min = 0.0d0
        s_max = 0.0d0
        s_cond = huge(1.0d0)
        rank_def = n_mix
      end if

      if (dg_frag%id == 0 .and. .not. quiet_mode) then
        write(*,'(1x,a,a,i0,a,i0,a,i0,a,i0,a,i0,36(a,1pe12.4),a,i0)') &
          '[DG-WPW-MIXED-NEIGHBOR-H] face-neighbor WP block:', &
          ' ifrag=', ifrag, ' nw=', n_w, ' npw=', n_pw, ' npfrag=', n_pfrag, ' filtered_drop_G2=', drop_g2, &
          ' ||S_WP_self||=', swp_self_norm, &
          ' ||S_WP_AB||=', swp_neigh_norm, &
          ' S_AB_over_self=', swp_neigh_norm / max(swp_self_norm, 1.0d-300), &
          ' S_AB_sv_min=', swp_ab_sv_min, ' S_AB_sv_max=', swp_ab_sv_max, &
          ' ||S_WP_AB_perp||=', swp_neigh_perp_norm, &
          ' S_perp_over_raw=', swp_perp_ratio, &
          ' ||H_WP_self||=', hwp_self_norm, &
          ' ||H_WP_AB||=', hwp_neigh_norm, ' ||H_PW_AB||=', hpw_neigh_norm, &
          ' AB_over_self=', hwp_neigh_norm / max(hwp_self_norm, 1.0d-300), &
          ' ||H_WP_AB_perp||=', hwp_neigh_perp_norm, &
          ' H_perp_over_raw=', hwp_perp_ratio, &
          ' ||H_WP_AB_local||=', hlocal_neigh_norm, &
          ' ||T_WP_AB_gradchi||=', tgrad_neigh_norm, &
          ' ||T_WP_AB_iGchi||=', tig_neigh_norm, &
          ' ||S_WP_AB_filtered||=', swp_neigh_filtered_norm, &
          ' ||H_WP_AB_filtered||=', hwp_neigh_filtered_norm, &
          ' ||T_gradchi_AB_filtered||=', tgrad_filtered_norm, &
          ' ||T_iGchi_AB_filtered||=', tig_filtered_norm, &
          ' ||T_AB_current||=', tcurrent_neigh_norm, &
          ' ||T_AB_sym||=', tsym_neigh_norm, &
          ' ||T_AB_current_minus_sym||=', tdiff_sym_norm, &
          ' H_herm=', h_herm, ' S_herm=', s_herm, &
          ' eval_shift_vs_self=', eval_shift_self, ' eval_min=', eval_min, &
          ' eval_shift_vs_self_filtered=', eval_shift_filt_self, &
          ' eval_shift_raw_to_filtered=', eval_shift_raw_filt, &
          ' eval_filtered_min=', eval_filt_min, ' eval_filtered_max=', eval_filt_max, &
          ' eval_self_min=', eval_self_min, ' eval_self_max=', eval_self_max, &
          ' eval_max=', eval_max, ' S_min=', s_min, ' S_cond=', s_cond, ' S_rank_def=', rank_def
        write(*,'(1x,a,a,i0,9(a,1pe12.4))') &
          '[DG-WPW-MIXED-NEIGHBOR-H] gradchi face diagnostic:', &
          ' ifrag=', ifrag, &
          ' T_gradchi_face_x=', tgrad_face_norm(1), &
          ' T_gradchi_face_y=', tgrad_face_norm(2), &
          ' T_gradchi_face_z=', tgrad_face_norm(3), &
          ' gradchi_A_face_x=', gradchi_a_face_norm(1), &
          ' gradchi_A_face_y=', gradchi_a_face_norm(2), &
          ' gradchi_A_face_z=', gradchi_a_face_norm(3), &
          ' gradchi_B_face_x=', gradchi_b_face_norm(1), &
          ' gradchi_B_face_y=', gradchi_b_face_norm(2), &
          ' gradchi_B_face_z=', gradchi_b_face_norm(3)
      end if
      deallocate(S_ab, H_ab, H_ab_perp, S_ab_perp, T_grad_ab, T_ig_ab, T_sym_ab, H_local_ab)
      deallocate(S_mix, H_mix, S_filt, H_filt, S_self, H_self, S_self_ref, H_self_ref, S_work, H_work)
      deallocate(S_wp, H_wp, H_wp_local, T_wp)
      deallocate(T_wp_gradchi, T_wp_igchi, T_wp_sym, T_wp_gradchi_face, H_ww)
      deallocate(wval, grad_wval, eval, eval_raw, eval_filt, eval_self, rwork, work)
      if (allocated(pfrag_ids)) deallocate(pfrag_ids)
      if (allocated(pfrag_axis)) deallocate(pfrag_axis)
      if (allocated(pfrag_side)) deallocate(pfrag_side)
    end do

    if (do_sorth .and. allocated(dg_frag%wpw_reduced_dim)) then
      dg_frag%wpw_reduced_ready = any(dg_frag%wpw_reduced_dim(:) > 0)
      if (dg_frag%wpw_reduced_ready) then
        dg_frag%wpw_reduced_shell = neighbor_shell
      else
        dg_frag%wpw_reduced_shell = 0
      end if
      if (.not. dg_frag%wpw_reduced_ready) then
        red_lambda_min_all = 0.0d0
        red_s_min_all = 0.0d0
      end if
      if (dg_frag%id == 0 .and. .not. quiet_mode) then
        write(*,'(1x,a,5(a,i0),7(a,1pe12.4),a,l1,2(a,a))') &
          '[DG-WPW-RED-INIT] summary:', &
          ' local_fragments=', n_local_frag, &
          ' neighbor_shell=', neighbor_shell, &
          ' dim_min=', minval(dg_frag%wpw_reduced_dim), &
          ' dim_max=', maxval(dg_frag%wpw_reduced_dim), &
          ' max_dim=', dg_frag%wpw_reduced_max_dim, &
          ' lambda_keep_min_global=', red_lambda_min_all, &
          ' lambda_keep_max_global=', red_lambda_max_all, &
          ' max||S_SN_after||=', red_s_sn_max, &
          ' max||S_NN-I||=', red_snn_i_max, &
          ' max_H_herm=', red_h_herm_max, &
          ' S_eval_min_global=', red_s_min_all, &
          ' S_eval_max_global=', red_s_max_all, &
          ' ready=', dg_frag%wpw_reduced_ready, &
          ' backend=', trim(dg_mixed_z_local_prop_backend), &
          ' shell_source=', 'namelist_or_backend_alias'
      end if
    end if

    if (dg_frag%id == 0 .and. .not. quiet_mode) then
      write(*,'(1x,a,2(a,1pe12.4),a,l1)') &
        '[DG-WPW-MIXED-NEIGHBOR-H] max over local fragments:', &
        ' ||H_WP_AB||=', neigh_norm_max, ' eval_shift_vs_self=', shift_max, &
        ' grad_cache=', have_grad_cache
    end if
  end subroutine diagnose_wpw_mixed_neighbor_hamiltonian

  subroutine diagnose_self_like_eigen_tracking(label, ifrag, H_self_in, S_self_in, H_ext_in, S_ext_in, n_self, n_ext)
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(in) :: ifrag, n_self, n_ext
    complex(8), intent(in) :: H_self_in(n_self, n_self), S_self_in(n_self, n_self)
    complex(8), intent(in) :: H_ext_in(n_ext, n_ext), S_ext_in(n_ext, n_ext)

    integer :: info, lwork, n, m, best_m, c05, c08, c09, c09c05, c09c08, top_n
    real(8) :: best_w, ov_abs, shift, w_min, w_max, w_avg
    real(8) :: ext_self_norm, comp_min, comp_max, comp_avg, best_comp
    real(8) :: shift05_max, shift08_max, shift09_max, shift_all_max
    complex(8) :: ov
    complex(8), allocatable :: Hs(:,:), Ss(:,:), He(:,:), Se(:,:), work(:), svec(:)
    integer, allocatable :: best_idx(:)
    real(8), allocatable :: es(:), ee(:), rwork(:), best_weight(:), best_component(:), best_shift(:)
    external :: zhegv

    if (n_self <= 0 .or. n_ext < n_self) return
    allocate(Hs(n_self,n_self), Ss(n_self,n_self), He(n_ext,n_ext), Se(n_ext,n_ext))
    allocate(es(n_self), ee(n_ext), rwork(max(1, 3*n_ext - 2)), work(1), svec(n_self))
    allocate(best_idx(n_self), best_weight(n_self), best_component(n_self), best_shift(n_self))
    Hs(:, :) = H_self_in(:, :)
    Ss(:, :) = S_self_in(:, :)
    lwork = -1
    call ZHEGV(1, 'V', 'U', n_self, Hs, n_self, Ss, n_self, es, work, lwork, rwork, info)
    lwork = max(1, int(real(work(1), kind=8)))
    deallocate(work)
    allocate(work(lwork))
    Hs(:, :) = H_self_in(:, :)
    Ss(:, :) = S_self_in(:, :)
    call ZHEGV(1, 'V', 'U', n_self, Hs, n_self, Ss, n_self, es, work, lwork, rwork, info)
    if (info /= 0) then
      write(*,'(1x,a,a,a,i0,a,i0)') '[DG-WPW-MIXED-NEIGHBOR-TRACK] ', trim(label), &
        ' ifrag=', ifrag, ' self_zhegv_info=', info
      deallocate(Hs, Ss, He, Se, es, ee, rwork, work, svec)
      deallocate(best_idx, best_weight, best_component, best_shift)
      return
    end if

    if (size(rwork) < max(1, 3*n_ext - 2)) then
      deallocate(rwork)
      allocate(rwork(max(1, 3*n_ext - 2)))
    end if
    if (size(work) < 1) then
      if (allocated(work)) deallocate(work)
      allocate(work(1))
    end if
    He(:, :) = H_ext_in(:, :)
    Se(:, :) = S_ext_in(:, :)
    lwork = -1
    call ZHEGV(1, 'V', 'U', n_ext, He, n_ext, Se, n_ext, ee, work, lwork, rwork, info)
    lwork = max(1, int(real(work(1), kind=8)))
    deallocate(work)
    allocate(work(lwork))
    He(:, :) = H_ext_in(:, :)
    Se(:, :) = S_ext_in(:, :)
    call ZHEGV(1, 'V', 'U', n_ext, He, n_ext, Se, n_ext, ee, work, lwork, rwork, info)
    if (info /= 0) then
      write(*,'(1x,a,a,a,i0,a,i0)') '[DG-WPW-MIXED-NEIGHBOR-TRACK] ', trim(label), &
        ' ifrag=', ifrag, ' ext_zhegv_info=', info
      deallocate(Hs, Ss, He, Se, es, ee, rwork, work, svec)
      deallocate(best_idx, best_weight, best_component, best_shift)
      return
    end if

    w_min = huge(1.0d0)
    w_max = 0.0d0
    w_avg = 0.0d0
    comp_min = huge(1.0d0)
    comp_max = 0.0d0
    comp_avg = 0.0d0
    shift_all_max = 0.0d0
    shift05_max = 0.0d0
    shift08_max = 0.0d0
    shift09_max = 0.0d0
    c05 = 0
    c08 = 0
    c09 = 0
    c09c05 = 0
    c09c08 = 0
    do n = 1, n_self
      best_w = -1.0d0
      best_m = 1
      best_comp = 0.0d0
      do m = 1, n_ext
        svec(1:n_self) = matmul(S_self_in(1:n_self, 1:n_self), He(1:n_self, m))
        ext_self_norm = real(sum(conjg(He(1:n_self, m)) * svec(1:n_self)), kind=8)
        ext_self_norm = max(ext_self_norm, 1.0d-300)
        ov = sum(conjg(Hs(1:n_self, n)) * svec(1:n_self))
        ov_abs = abs(ov)**2 / ext_self_norm
        if (ov_abs > best_w) then
          best_w = ov_abs
          best_m = m
          best_comp = ext_self_norm
        end if
      end do
      shift = ee(best_m) - es(n)
      best_idx(n) = best_m
      best_weight(n) = best_w
      best_component(n) = best_comp
      best_shift(n) = shift
      w_min = min(w_min, best_w)
      w_max = max(w_max, best_w)
      w_avg = w_avg + best_w
      comp_min = min(comp_min, best_comp)
      comp_max = max(comp_max, best_comp)
      comp_avg = comp_avg + best_comp
      shift_all_max = max(shift_all_max, abs(shift))
      if (best_w >= 0.5d0) then
        c05 = c05 + 1
        shift05_max = max(shift05_max, abs(shift))
      end if
      if (best_w >= 0.8d0) then
        c08 = c08 + 1
        shift08_max = max(shift08_max, abs(shift))
      end if
      if (best_w >= 0.9d0) then
        c09 = c09 + 1
        shift09_max = max(shift09_max, abs(shift))
      end if
      if (best_w >= 0.9d0 .and. best_comp >= 0.5d0) c09c05 = c09c05 + 1
      if (best_w >= 0.9d0 .and. best_comp >= 0.8d0) c09c08 = c09c08 + 1
    end do
    w_avg = w_avg / dble(max(1, n_self))
    if (w_min == huge(1.0d0)) w_min = 0.0d0
    comp_avg = comp_avg / dble(max(1, n_self))
    if (comp_min == huge(1.0d0)) comp_min = 0.0d0

    write(*,'(1x,a,a,a,i0,a,i0,a,i0,10(a,1pe12.4),5(a,i0))') &
      '[DG-WPW-MIXED-NEIGHBOR-TRACK] ', trim(label), &
      ' ifrag=', ifrag, ' nself=', n_self, ' next=', n_ext, &
      ' weight_min=', w_min, ' weight_avg=', w_avg, ' weight_max=', w_max, &
      ' self_component_min=', comp_min, ' self_component_avg=', comp_avg, ' self_component_max=', comp_max, &
      ' shift_all_max=', shift_all_max, &
      ' shift_w05_max=', shift05_max, ' shift_w08_max=', shift08_max, ' shift_w09_max=', shift09_max, &
      ' count_w05=', c05, ' count_w08=', c08, ' count_w09=', c09, &
      ' count_w09c05=', c09c05, ' count_w09c08=', c09c08

    top_n = wpw_neighbor_track_top_n()
    call print_self_like_tracking_outliers(label, ifrag, 'all', top_n, n_self, es, ee, &
      best_idx, best_weight, best_component, best_shift, -1.0d0, -1.0d0)
    call print_self_like_tracking_outliers(label, ifrag, 'w05', top_n, n_self, es, ee, &
      best_idx, best_weight, best_component, best_shift, 0.5d0, -1.0d0)
    call print_self_like_tracking_outliers(label, ifrag, 'w08', top_n, n_self, es, ee, &
      best_idx, best_weight, best_component, best_shift, 0.8d0, -1.0d0)
    call print_self_like_tracking_outliers(label, ifrag, 'w09', top_n, n_self, es, ee, &
      best_idx, best_weight, best_component, best_shift, 0.9d0, -1.0d0)
    call print_self_like_tracking_outliers(label, ifrag, 'w09c05', top_n, n_self, es, ee, &
      best_idx, best_weight, best_component, best_shift, 0.9d0, 0.5d0)
    call print_self_like_tracking_outliers(label, ifrag, 'w09c08', top_n, n_self, es, ee, &
      best_idx, best_weight, best_component, best_shift, 0.9d0, 0.8d0)

    deallocate(Hs, Ss, He, Se, es, ee, rwork, work, svec)
    deallocate(best_idx, best_weight, best_component, best_shift)
  end subroutine diagnose_self_like_eigen_tracking

  integer function wpw_neighbor_track_top_n() result(top_n)
    implicit none
    character(len=32) :: env_value
    integer :: env_len, env_stat, parsed
    logical, save :: initialized = .false.
    integer, save :: cached_top_n = 8

    if (.not. initialized) then
      env_value = ''
      call get_environment_variable('SALMON_DG_WPW_TRACK_TOP_N', env_value, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_value(1:env_len), *, iostat=env_stat) parsed
        if (env_stat == 0) cached_top_n = max(0, parsed)
      end if
      initialized = .true.
    end if
    top_n = cached_top_n
  end function wpw_neighbor_track_top_n

  subroutine print_self_like_tracking_outliers(label, ifrag, tag, top_n, n_self, es, ee, &
      best_idx, best_weight, best_component, best_shift, weight_threshold, component_threshold)
    implicit none
    character(len=*), intent(in) :: label, tag
    integer, intent(in) :: ifrag, top_n, n_self
    integer, intent(in) :: best_idx(n_self)
    real(8), intent(in) :: es(n_self), ee(:), best_weight(n_self), best_component(n_self), best_shift(n_self)
    real(8), intent(in) :: weight_threshold, component_threshold

    integer :: rank, n, best_n, count_selected, m_ext
    real(8) :: best_abs, abs_shift
    logical, allocatable :: used(:)

    if (top_n <= 0 .or. n_self <= 0) return
    allocate(used(n_self))
    used(:) = .false.
    count_selected = 0

    do rank = 1, min(top_n, n_self)
      best_n = 0
      best_abs = -1.0d0
      do n = 1, n_self
        if (used(n)) cycle
        if (weight_threshold >= 0.0d0 .and. best_weight(n) < weight_threshold) cycle
        if (component_threshold >= 0.0d0 .and. best_component(n) < component_threshold) cycle
        abs_shift = abs(best_shift(n))
        if (abs_shift > best_abs) then
          best_abs = abs_shift
          best_n = n
        end if
      end do
      if (best_n <= 0) exit
      used(best_n) = .true.
      count_selected = count_selected + 1
      m_ext = best_idx(best_n)
      if (m_ext < 1 .or. m_ext > size(ee)) cycle
      write(*,'(1x,a,a,a,a,a,i0,a,i0,a,i0,a,i0,6(a,1pe12.4))') &
        '[DG-WPW-MIXED-NEIGHBOR-TRACK-OUTLIER] ', trim(label), &
        ' tag=', trim(tag), ' rank=', count_selected, &
        ' ifrag=', ifrag, ' n_self=', best_n, ' m_ext=', m_ext, &
        ' E_self=', es(best_n), ' E_ext=', ee(m_ext), &
        ' shift=', best_shift(best_n), ' abs_shift=', abs(best_shift(best_n)), &
        ' weight=', best_weight(best_n), ' self_component=', best_component(best_n)
    end do

    if (count_selected == 0) then
      write(*,'(1x,a,a,a,a,a,i0,a,i0)') &
        '[DG-WPW-MIXED-NEIGHBOR-TRACK-OUTLIER] ', trim(label), &
        ' tag=', trim(tag), ' ifrag=', ifrag, ' count=0'
    end if
    deallocate(used)
  end subroutine print_self_like_tracking_outliers

  subroutine diagnose_overlap_schur_effects(label, ifrag, H_self_in, S_self_in, H_ext_in, S_ext_in, n_self, n_ext)
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(in) :: ifrag, n_self, n_ext
    complex(8), intent(in) :: H_self_in(n_self, n_self), S_self_in(n_self, n_self)
    complex(8), intent(in) :: H_ext_in(n_ext, n_ext), S_ext_in(n_ext, n_ext)

    integer :: n_neigh, info, top_n, n, m, a, rank, best_n, best_m, best_a, selected_modes(3)
    integer :: c09c08
    real(8) :: s_ext_min, s_ext_max, s_self_min, s_self_max, s_nn_min, s_nn_max
    real(8) :: sn_norm, sn_schur_norm, delta_norm, delta_eval_min, delta_eval_max
    real(8) :: best_w, best_comp, ext_self_norm, ov_abs, shift, best_abs, contrib, best_contrib
    real(8) :: en, denom
    complex(8) :: ov, vna
    complex(8), allocatable :: Horth_ext(:,:), Horth_self(:,:), delta_h(:,:), delta_work(:,:)
    complex(8), allocatable :: Hs(:,:), Ss(:,:), He(:,:), Se(:,:), Hn(:,:), Sn(:,:)
    complex(8), allocatable :: ys(:,:), svec(:), sn_schur(:,:), s_nn_inv(:,:)
    complex(8), allocatable :: work_vec(:)
    integer, allocatable :: best_idx(:), outlier_idx(:), outlier_m(:)
    real(8), allocatable :: es(:), ee(:), en_neigh(:), eval_delta(:), rwork(:)
    real(8), allocatable :: best_weight(:), best_component(:), best_shift(:), outlier_abs(:)
    logical, allocatable :: used(:)
    external :: zhegv, zheev

    if (n_self <= 0 .or. n_ext <= n_self) return
    n_neigh = n_ext - n_self
    top_n = wpw_neighbor_track_top_n()
    if (top_n <= 0) top_n = 3

    allocate(Horth_ext(n_ext,n_ext), Horth_self(n_self,n_self), delta_h(n_self,n_self), delta_work(n_self,n_self))
    call build_symmetric_orthogonal_hamiltonian(H_ext_in, S_ext_in, n_ext, Horth_ext, info, s_ext_min, s_ext_max)
    if (info /= 0) then
      write(*,'(1x,a,a,a,i0,a,i0)') '[DG-WPW-MIXED-SCHUR] ', trim(label), &
        ' ifrag=', ifrag, ' ext_orth_info=', info
      deallocate(Horth_ext, Horth_self, delta_h, delta_work)
      return
    end if
    call build_symmetric_orthogonal_hamiltonian(H_self_in, S_self_in, n_self, Horth_self, info, s_self_min, s_self_max)
    if (info /= 0) then
      write(*,'(1x,a,a,a,i0,a,i0)') '[DG-WPW-MIXED-SCHUR] ', trim(label), &
        ' ifrag=', ifrag, ' self_orth_info=', info
      deallocate(Horth_ext, Horth_self, delta_h, delta_work)
      return
    end if

    delta_h(:, :) = Horth_ext(1:n_self, 1:n_self) - Horth_self(:, :)
    call hermitize_matrix(delta_h, n_self)
    delta_norm = sqrt(sum(abs(delta_h)**2))
    allocate(eval_delta(n_self), rwork(max(1, 3*max(n_ext, n_self) - 2)), work_vec(1))
    delta_work(:, :) = delta_h(:, :)
    call zheev_with_query(delta_work, n_self, eval_delta, info)
    if (info == 0) then
      delta_eval_min = eval_delta(1)
      delta_eval_max = eval_delta(n_self)
    else
      delta_eval_min = 0.0d0
      delta_eval_max = 0.0d0
    end if

    sn_norm = sqrt(sum(abs(S_ext_in(1:n_self, n_self+1:n_ext))**2))
    allocate(s_nn_inv(n_neigh,n_neigh), sn_schur(n_self,n_self))
    call build_hermitian_inverse(S_ext_in(n_self+1:n_ext, n_self+1:n_ext), n_neigh, s_nn_inv, info, s_nn_min, s_nn_max)
    if (info == 0) then
      sn_schur(:, :) = matmul(matmul(S_ext_in(1:n_self, n_self+1:n_ext), s_nn_inv), &
        S_ext_in(n_self+1:n_ext, 1:n_self))
      sn_schur_norm = sqrt(sum(abs(sn_schur)**2))
    else
      sn_schur_norm = 0.0d0
      s_nn_min = 0.0d0
      s_nn_max = 0.0d0
    end if

    allocate(Hs(n_self,n_self), Ss(n_self,n_self), He(n_ext,n_ext), Se(n_ext,n_ext))
    allocate(Hn(n_neigh,n_neigh), Sn(n_neigh,n_neigh), ys(n_self,n_self), svec(n_self))
    allocate(es(n_self), ee(n_ext), en_neigh(n_neigh), best_idx(n_self), best_weight(n_self))
    allocate(best_component(n_self), best_shift(n_self), outlier_idx(top_n), outlier_m(top_n), outlier_abs(top_n))

    Hs(:, :) = H_self_in(:, :)
    Ss(:, :) = S_self_in(:, :)
    call zhegv_with_query(Hs, Ss, n_self, es, info)
    if (info /= 0) then
      write(*,'(1x,a,a,a,i0,a,i0)') '[DG-WPW-MIXED-SCHUR] ', trim(label), &
        ' ifrag=', ifrag, ' self_zhegv_info=', info
      deallocate(Horth_ext, Horth_self, delta_h, delta_work, eval_delta, rwork, work_vec)
      deallocate(s_nn_inv, sn_schur, Hs, Ss, He, Se, Hn, Sn, ys, svec, es, ee, en_neigh)
      deallocate(best_idx, best_weight, best_component, best_shift, outlier_idx, outlier_m, outlier_abs)
      return
    end if
    He(:, :) = H_ext_in(:, :)
    Se(:, :) = S_ext_in(:, :)
    call zhegv_with_query(He, Se, n_ext, ee, info)
    if (info /= 0) then
      write(*,'(1x,a,a,a,i0,a,i0)') '[DG-WPW-MIXED-SCHUR] ', trim(label), &
        ' ifrag=', ifrag, ' ext_zhegv_info=', info
      deallocate(Horth_ext, Horth_self, delta_h, delta_work, eval_delta, rwork, work_vec)
      deallocate(s_nn_inv, sn_schur, Hs, Ss, He, Se, Hn, Sn, ys, svec, es, ee, en_neigh)
      deallocate(best_idx, best_weight, best_component, best_shift, outlier_idx, outlier_m, outlier_abs)
      return
    end if

    ys(:, :) = Horth_self(:, :)
    call zheev_with_query(ys, n_self, es, info)

    c09c08 = 0
    do n = 1, n_self
      best_w = -1.0d0
      best_m = 1
      best_comp = 0.0d0
      do m = 1, n_ext
        svec(1:n_self) = matmul(S_self_in(1:n_self, 1:n_self), He(1:n_self, m))
        ext_self_norm = real(sum(conjg(He(1:n_self, m)) * svec(1:n_self)), kind=8)
        ext_self_norm = max(ext_self_norm, 1.0d-300)
        ov = sum(conjg(Hs(1:n_self, n)) * svec(1:n_self))
        ov_abs = abs(ov)**2 / ext_self_norm
        if (ov_abs > best_w) then
          best_w = ov_abs
          best_m = m
          best_comp = ext_self_norm
        end if
      end do
      best_idx(n) = best_m
      best_weight(n) = best_w
      best_component(n) = best_comp
      best_shift(n) = ee(best_m) - es(n)
      if (best_w >= 0.9d0 .and. best_comp >= 0.8d0) c09c08 = c09c08 + 1
    end do

    write(*,'(1x,a,a,a,i0,9(a,1pe12.4),a,i0)') &
      '[DG-WPW-MIXED-SCHUR] ', trim(label), ' ifrag=', ifrag, &
      ' ||S_SN||=', sn_norm, &
      ' ||S_SN_SNNinv_S_NS||=', sn_schur_norm, &
      ' S_ext_min=', s_ext_min, ' S_ext_max=', s_ext_max, &
      ' S_self_min=', s_self_min, ' S_self_max=', s_self_max, &
      ' S_NN_min=', s_nn_min, ' S_NN_max=', s_nn_max, &
      ' ||DeltaH_SS||=', delta_norm, &
      ' count_w09c08=', c09c08
    write(*,'(1x,a,a,a,i0,2(a,1pe12.4))') &
      '[DG-WPW-MIXED-SCHUR] ', trim(label), ' ifrag=', ifrag, &
      ' DeltaH_eval_min=', delta_eval_min, ' DeltaH_eval_max=', delta_eval_max

    call collect_top_self_like_outliers(n_self, top_n, best_weight, best_component, best_shift, best_idx, &
      outlier_idx, outlier_m, outlier_abs)
    do rank = 1, top_n
      n = outlier_idx(rank)
      if (n <= 0) cycle
      en = es(n)
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,6(a,1pe12.4))') &
        '[DG-WPW-MIXED-SCHUR-OUTLIER] ', trim(label), &
        ' ifrag=', ifrag, ' rank=', rank, ' n_self=', n, ' m_ext=', outlier_m(rank), &
        ' E_self=', es(n), ' shift=', best_shift(n), ' weight=', best_weight(n), &
        ' self_component=', best_component(n), &
        ' dH_expect=', real(sum(conjg(ys(1:n_self, n)) * matmul(delta_h, ys(1:n_self, n))), kind=8), &
        ' abs_shift=', abs(best_shift(n))
    end do

    Hn(:, :) = H_ext_in(n_self+1:n_ext, n_self+1:n_ext)
    Sn(:, :) = S_ext_in(n_self+1:n_ext, n_self+1:n_ext)
    call zhegv_with_query(Hn, Sn, n_neigh, en_neigh, info)
    if (info == 0) then
      do rank = 1, top_n
        n = outlier_idx(rank)
        if (n <= 0) cycle
        en = es(n)
        selected_modes(:) = 0
        do a = 1, min(3, n_neigh)
          best_a = 0
          best_contrib = -1.0d0
          do m = 1, n_neigh
            if (a > 1) then
              if (any(selected_modes(1:a-1) == m)) cycle
            end if
            vna = sum(conjg(Hs(1:n_self, n)) * matmul( &
              H_ext_in(1:n_self, n_self+1:n_ext) - en * S_ext_in(1:n_self, n_self+1:n_ext), Hn(1:n_neigh, m)))
            denom = max(abs(en - en_neigh(m)), 1.0d-12)
            contrib = abs(vna)**2 / denom
            if (contrib > best_contrib) then
              best_contrib = contrib
              best_a = m
            end if
          end do
          if (best_a > 0) then
            selected_modes(a) = best_a
            vna = sum(conjg(Hs(1:n_self, n)) * matmul( &
              H_ext_in(1:n_self, n_self+1:n_ext) - en * S_ext_in(1:n_self, n_self+1:n_ext), Hn(1:n_neigh, best_a)))
            write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,5(a,1pe12.4))') &
              '[DG-WPW-MIXED-SCHUR-MODE] ', trim(label), &
              ' ifrag=', ifrag, ' outlier_rank=', rank, ' n_self=', n, ' mode=', best_a, &
              ' E_self=', en, ' E_neigh=', en_neigh(best_a), &
              ' abs_V=', abs(vna), ' contribution=', best_contrib, &
              ' denom_abs=', abs(en - en_neigh(best_a))
          end if
        end do
      end do
    else
      write(*,'(1x,a,a,a,i0,a,i0)') '[DG-WPW-MIXED-SCHUR] ', trim(label), &
        ' ifrag=', ifrag, ' neighbor_zhegv_info=', info
    end if

    deallocate(Horth_ext, Horth_self, delta_h, delta_work, eval_delta, rwork, work_vec)
    deallocate(s_nn_inv, sn_schur, Hs, Ss, He, Se, Hn, Sn, ys, svec, es, ee, en_neigh)
    deallocate(best_idx, best_weight, best_component, best_shift, outlier_idx, outlier_m, outlier_abs)
  end subroutine diagnose_overlap_schur_effects

  subroutine collect_top_self_like_outliers(n_self, top_n, best_weight, best_component, best_shift, best_idx, &
      outlier_idx, outlier_m, outlier_abs)
    implicit none
    integer, intent(in) :: n_self, top_n, best_idx(n_self)
    real(8), intent(in) :: best_weight(n_self), best_component(n_self), best_shift(n_self)
    integer, intent(out) :: outlier_idx(top_n), outlier_m(top_n)
    real(8), intent(out) :: outlier_abs(top_n)

    integer :: rank, n, best_n
    real(8) :: best_abs
    logical, allocatable :: used(:)

    outlier_idx(:) = 0
    outlier_m(:) = 0
    outlier_abs(:) = 0.0d0
    allocate(used(n_self))
    used(:) = .false.
    do rank = 1, top_n
      best_n = 0
      best_abs = -1.0d0
      do n = 1, n_self
        if (used(n)) cycle
        if (best_weight(n) < 0.9d0 .or. best_component(n) < 0.8d0) cycle
        if (abs(best_shift(n)) > best_abs) then
          best_abs = abs(best_shift(n))
          best_n = n
        end if
      end do
      if (best_n <= 0) exit
      used(best_n) = .true.
      outlier_idx(rank) = best_n
      outlier_m(rank) = best_idx(best_n)
      outlier_abs(rank) = best_abs
    end do
    deallocate(used)
  end subroutine collect_top_self_like_outliers

  subroutine build_s_orthogonal_neighbor_variant(label, ifrag, H_in, S_in, n_self, n_ext, H_out, S_out)
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(in) :: ifrag, n_self, n_ext
    complex(8), intent(in) :: H_in(n_ext,n_ext), S_in(n_ext,n_ext)
    complex(8), intent(out) :: H_out(n_ext,n_ext), S_out(n_ext,n_ext)

    integer :: n_neigh, info, n_drop, n_keep
    real(8) :: sn_before, sn_after, snn_min, snn_max, sss_min, sss_max, tol
    complex(8), allocatable :: sss_inv(:,:), xmat(:,:), work_nn(:,:)
    real(8), allocatable :: eval_nn(:)

    H_out(:, :) = H_in(:, :)
    S_out(:, :) = S_in(:, :)
    if (n_ext <= n_self .or. n_self <= 0) return
    n_neigh = n_ext - n_self
    allocate(sss_inv(n_self,n_self), xmat(n_self,n_neigh), work_nn(n_neigh,n_neigh), eval_nn(n_neigh))

    call build_hermitian_inverse(S_in(1:n_self, 1:n_self), n_self, sss_inv, info, sss_min, sss_max)
    if (info /= 0) then
      write(*,'(1x,a,a,a,i0,a,i0)') '[DG-WPW-MIXED-SORTH] ', trim(label), &
        ' ifrag=', ifrag, ' SSS_inverse_info=', info
      deallocate(sss_inv, xmat, work_nn, eval_nn)
      return
    end if

    xmat(:, :) = matmul(sss_inv, S_in(1:n_self, n_self+1:n_ext))
    sn_before = sqrt(sum(abs(S_in(1:n_self, n_self+1:n_ext))**2))

    call apply_s_orthogonal_transform(H_in, xmat, n_self, n_ext, H_out)
    call apply_s_orthogonal_transform(S_in, xmat, n_self, n_ext, S_out)
    call hermitize_matrix(H_out, n_ext)
    call hermitize_matrix(S_out, n_ext)

    sn_after = sqrt(sum(abs(S_out(1:n_self, n_self+1:n_ext))**2))
    work_nn(:, :) = S_out(n_self+1:n_ext, n_self+1:n_ext)
    call zheev_with_query(work_nn, n_neigh, eval_nn, info)
    if (info == 0) then
      snn_min = eval_nn(1)
      snn_max = eval_nn(n_neigh)
      tol = wpw_neighbor_sorth_tol()
      n_drop = count(eval_nn(1:n_neigh) <= tol)
      n_keep = n_neigh - n_drop
    else
      snn_min = 0.0d0
      snn_max = 0.0d0
      tol = wpw_neighbor_sorth_tol()
      n_drop = n_neigh
      n_keep = 0
    end if

    write(*,'(1x,a,a,a,i0,7(a,1pe12.4),3(a,i0))') &
      '[DG-WPW-MIXED-SORTH] ', trim(label), ' ifrag=', ifrag, &
      ' ||S_SN_before||=', sn_before, &
      ' ||S_SN_after||=', sn_after, &
      ' SSS_min=', sss_min, ' SSS_max=', sss_max, &
      ' SNNprime_min=', snn_min, ' SNNprime_max=', snn_max, &
      ' SNNprime_tol=', tol, &
      ' n_keep=', n_keep, ' n_drop=', n_drop, ' info=', info

    deallocate(sss_inv, xmat, work_nn, eval_nn)
  end subroutine build_s_orthogonal_neighbor_variant

  logical function wpw_reduced_density_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_DENSITY_DIAG', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_density_diag_enabled

  logical function wpw_reduced_canonical_reproject_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANONICAL_REPROJECT_EXPERIMENTAL', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_canonical_reproject_enabled

  logical function wpw_reduced_canon_pack_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANON_PACK_DIAG', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_canon_pack_diag_enabled

  logical function wpw_reduced_canon_recon_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANON_RECON_DIAG', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_canon_recon_diag_enabled

  logical function wpw_reduced_canon_obs_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANON_OBS_DIAG', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_canon_obs_diag_enabled

  logical function wpw_reduced_canon_hook_dryrun_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANON_HOOK_DRYRUN', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_canon_hook_dryrun_enabled

  logical function wpw_reduced_canon_use_density_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANON_USE_DENSITY', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_canon_use_density_enabled

  logical function wpw_reduced_canon_use_pz_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANON_USE_PZ', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_canon_use_pz_enabled

  logical function wpw_reduced_canon_mixedz_pz_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANON_MIXEDZ_PZ_DIAG', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_canon_mixedz_pz_diag_enabled

  logical function wpw_reduced_drift_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_DRIFT_DIAG', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_drift_diag_enabled

  integer function wpw_reduced_drift_interval() result(interval)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat, parsed
    logical, save :: initialized = .false.
    integer, save :: cached_interval = 100

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_DRIFT_INTERVAL', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env(1:env_len), *, iostat=env_stat) parsed
        if (env_stat == 0 .and. parsed > 0) cached_interval = parsed
      end if
      initialized = .true.
    end if
    interval = cached_interval
  end function wpw_reduced_drift_interval

  logical function wpw_reduced_propagated_debug_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_PROPAGATED_DEBUG', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_propagated_debug_enabled

  logical function wpw_reduced_state_prop_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_STATE_PROP_DIAG', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_state_prop_diag_enabled

  logical function wpw_reduced_sample_u_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_SAMPLE_U_DIAG', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      if (.not. cached_enabled) then
        env = ''
        call get_environment_variable('SALMON_DG_WPW_REDUCED_PRODOP_ACTION_DIAG', env, length=env_len, status=env_stat)
        if (env_stat == 0 .and. env_len > 0) then
          select case (adjustl(trim(env(1:env_len))))
          case ('1','y','Y','yes','YES','true','TRUE','on','ON')
            cached_enabled = .true.
          case default
            cached_enabled = .false.
          end select
        end if
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_sample_u_diag_enabled

  logical function wpw_reduced_embed_local_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_EMBED_LOCAL_DIAG', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_embed_local_diag_enabled

  logical function wpw_reduced_embed_prodcoef_diag_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_EMBED_PRODCOEF_DIAG', env, &
        length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_embed_prodcoef_diag_enabled

  real(8) function wpw_reduced_sample_u_tol() result(tol)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    real(8) :: parsed
    logical, save :: initialized = .false.
    real(8), save :: cached_tol = 1.0d-10

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_SAMPLE_U_TOL', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env(1:env_len), *, iostat=env_stat) parsed
        if (env_stat == 0 .and. parsed > 0.0d0 .and. parsed < 1.0d0) cached_tol = parsed
      end if
      initialized = .true.
    end if
    tol = cached_tol
  end function wpw_reduced_sample_u_tol

  integer function wpw_reduced_keep_n() result(keep_n)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat, parsed
    logical, save :: initialized = .false.
    integer, save :: cached_keep_n = -1

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_KEEP_N', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env(1:env_len), *, iostat=env_stat) parsed
        if (env_stat == 0) cached_keep_n = max(0, parsed)
      end if
      initialized = .true.
    end if
    keep_n = cached_keep_n
  end function wpw_reduced_keep_n

  logical function wpw_reduced_pz_series_enabled() result(enabled)
    implicit none
    character(len=32) :: env
    integer :: env_len, env_stat
    logical, save :: initialized = .false.
    logical, save :: cached_enabled = .false.

    if (.not. initialized) then
      env = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_PZ_SERIES', env, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          cached_enabled = .true.
        case default
          cached_enabled = .false.
        end select
      end if
      initialized = .true.
    end if
    enabled = cached_enabled
  end function wpw_reduced_pz_series_enabled

  character(len=128) function wpw_reduced_pz_filename(keep_n) result(filename)
    implicit none
    integer, intent(in) :: keep_n

    if (keep_n < 0) then
      filename = 'dg_wpw_reduced_pz_cmp_keep_all.data'
    else
      write(filename,'(a,i0,a)') 'dg_wpw_reduced_pz_cmp_keep_', keep_n, '.data'
    end if
  end function wpw_reduced_pz_filename

  subroutine append_wpw_reduced_pz_cmp(step, time_au, keep_n, pz_prod, pz_reduced)
    implicit none
    integer, intent(in) :: step, keep_n
    real(8), intent(in) :: time_au, pz_prod, pz_reduced

    integer :: iunit, io_stat
    logical :: exists
    character(len=128) :: filename
    real(8) :: rel_diff

    filename = wpw_reduced_pz_filename(keep_n)
    inquire(file=trim(filename), exist=exists)
    if (exists) then
      open(newunit=iunit, file=trim(filename), status='old', position='append', action='write', iostat=io_stat)
    else
      open(newunit=iunit, file=trim(filename), status='replace', action='write', iostat=io_stat)
      if (io_stat == 0) then
        write(iunit,'(a)') '# WPW reduced diagnostic Pz comparison'
        write(iunit,'(a)') '# 1:step 2:time[a.u.] 3:keep_n 4:Pz_prod 5:Pz_reduced 6:Pz_diff 7:rel_Pz_diff'
      end if
    end if
    if (io_stat /= 0) return
    rel_diff = (pz_prod - pz_reduced) / max(abs(pz_prod), 1.0d-300)
    write(iunit,'(i10,1x,es24.16,1x,i10,4(1x,es24.16))') &
      step, time_au, keep_n, pz_prod, pz_reduced, pz_prod - pz_reduced, rel_diff
    close(iunit)
  end subroutine append_wpw_reduced_pz_cmp

  subroutine initialize_wpw_reduced_self_projection(dg_frag, state_s, state_e, did_project)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: state_s, state_e
    logical, intent(out) :: did_project

    integer :: i_local, ifrag, ispin, n_pw, n_w, nself, nred, nraw, nneigh, nbf
    integer :: ist, iw, ib, ipw, ix, iy, iz, gx, gy, gz, bx, by, bz
    integer :: n_pfrag, pfrag_ids(7), pidx, row0, axis, side, jfrag
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: global_idx, local_idx, info_inv
    real(8) :: vol_weight, phase_arg, chi, grad_chi(3), smin, smax
    real(8) :: norm_before, norm_after, norm_loss_sum, norm_loss_max
    complex(8) :: phase, psi
    complex(8), allocatable :: rhs(:), cproj(:), s_inv(:,:), s_tmp(:), cw0(:), pval(:)
    complex(8), allocatable :: basis_raw(:), wval(:), S_density(:,:), S_density_inv(:,:), rhs_density(:), c_raw(:)
    logical :: have_sdensity

    did_project = .false.
    if (.not. dg_frag%wpw_reduced_ready) return
    if (.not. allocated(dg_frag%wpw_reduced_S)) return
    if (.not. allocated(dg_frag%wpw_reduced_transform)) return
    if (.not. allocated(dg_frag%wpw_reduced_nraw)) return
    if (.not. allocated(dg_frag%coef_wpw_self)) return
    if (.not. allocated(dg_frag%coef_wpw_neighbor_reduced)) return
    if (.not. allocated(dg_frag%global_wannier_local_coef)) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return
    if (.not. allocated(dg_frag%phi_frag) .and. .not. allocated(dg_frag%phi_frag_c)) return

    n_pw = dg_frag%n_plane_waves
    if (n_pw <= 0) return
    vol_weight = product(dg_frag%hgs(1:3)) / &
      product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    if (allocated(dg_frag%phi_frag_c)) then
      p_lb1 = lbound(dg_frag%phi_frag_c, 1); p_ub1 = ubound(dg_frag%phi_frag_c, 1)
      p_lb2 = lbound(dg_frag%phi_frag_c, 2); p_ub2 = ubound(dg_frag%phi_frag_c, 2)
      p_lb3 = lbound(dg_frag%phi_frag_c, 3); p_ub3 = ubound(dg_frag%phi_frag_c, 3)
    else
      p_lb1 = lbound(dg_frag%phi_frag, 1); p_ub1 = ubound(dg_frag%phi_frag, 1)
      p_lb2 = lbound(dg_frag%phi_frag, 2); p_ub2 = ubound(dg_frag%phi_frag, 2)
      p_lb3 = lbound(dg_frag%phi_frag, 3); p_ub3 = ubound(dg_frag%phi_frag, 3)
    end if

    norm_loss_sum = 0.0d0
    norm_loss_max = 0.0d0
    do i_local = 1, size(dg_frag%wpw_reduced_dim)
      ifrag = dg_frag%ifrag_start + i_local - 1
      nred = dg_frag%wpw_reduced_dim(i_local)
      nself = dg_frag%wpw_reduced_nself(i_local)
      nraw = dg_frag%wpw_reduced_nraw(i_local)
      if (nred <= 0 .or. nself <= 0 .or. nraw <= 0) cycle
      if (i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
      n_w = dg_frag%global_wannier_local_nkeep(i_local)
      if (n_w <= 0 .or. nself /= n_w + n_pw) cycle
      nneigh = nred - nself
      n_pfrag = 1
      pfrag_ids(1) = ifrag
      do axis = 1, 3
        do side = -1, 1, 2
          jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle
          if (.not. any(pfrag_ids(1:n_pfrag) == jfrag) .and. n_pfrag < size(pfrag_ids)) then
            n_pfrag = n_pfrag + 1
            pfrag_ids(n_pfrag) = jfrag
          end if
        end do
      end do
      if (nraw /= n_w + n_pfrag * n_pw) cycle
      nbf = min(dg_frag%n_basis(ifrag, 1), size(dg_frag%global_wannier_local_coef, 1))
      if (allocated(dg_frag%phi_frag_c)) then
        nbf = min(nbf, size(dg_frag%phi_frag_c, 4))
      else
        nbf = min(nbf, size(dg_frag%phi_frag, 4))
      end if
      if (nbf <= 0) cycle

      allocate(rhs(nraw), cproj(nred), s_inv(nred,nred), s_tmp(nred), cw0(n_w), pval(n_pw))
      allocate(basis_raw(nraw), wval(n_w), S_density(nraw,nraw), S_density_inv(nraw,nraw), rhs_density(nraw), c_raw(nraw))
      do ispin = 1, dg_frag%nspin
        call build_hermitian_inverse(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
          nred, s_inv, info_inv, smin, smax)
        if (info_inv /= 0) cycle
        have_sdensity = .false.
        do ist = state_s, min(state_e, size(dg_frag%coef_wpw_self, 2), size(dg_frag%coef, 2))
          rhs(:) = (0.0d0, 0.0d0)
          rhs_density(:) = (0.0d0, 0.0d0)
          if (.not. have_sdensity) S_density(:, :) = (0.0d0, 0.0d0)
          cw0(:) = (0.0d0, 0.0d0)
          do iw = 1, n_w
            do ib = 1, nbf
              global_idx = dg_frag%index_basis(ib, ifrag, ispin)
              local_idx = 0
              if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) &
                local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
              if (local_idx <= 0 .or. local_idx > size(dg_frag%coef, 1)) cycle
              cw0(iw) = cw0(iw) + conjg(dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local)) * &
                dg_frag%coef(local_idx, ist, ispin)
            end do
          end do
          ! Hybrid metric convention: W-W is orthonormal in the reduced builder.
          ! Do not use the real-space W-W integral for coefficient projection.
          rhs(1:n_w) = cw0(1:n_w)
          do iz = 1, dg_frag%nxyz_domain(3, ifrag)
            gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
            bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bz < p_lb3 .or. bz > p_ub3) cycle
            do iy = 1, dg_frag%nxyz_domain(2, ifrag)
              gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
              by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              if (by < p_lb2 .or. by > p_ub2) cycle
              do ix = 1, dg_frag%nxyz_domain(1, ifrag)
                gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                if (bx < p_lb1 .or. bx > p_ub1) cycle
                psi = (0.0d0, 0.0d0)
                wval(:) = (0.0d0, 0.0d0)
                do ib = 1, nbf
                  global_idx = dg_frag%index_basis(ib, ifrag, ispin)
                  local_idx = 0
                  if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) &
                    local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
                  if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
                    if (allocated(dg_frag%phi_frag_c)) then
                      psi = psi + dg_frag%coef(local_idx, ist, ispin) * dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                    else
                      psi = psi + dg_frag%coef(local_idx, ist, ispin) * &
                        cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                    end if
                  end if
                  do iw = 1, n_w
                    if (allocated(dg_frag%phi_frag_c)) then
                      wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                        dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                    else
                      wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                        cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                    end if
                  end do
                end do
                basis_raw(:) = (0.0d0, 0.0d0)
                basis_raw(1:n_w) = wval(1:n_w)
                do pidx = 1, n_pfrag
                  row0 = n_w + (pidx - 1) * n_pw
                  call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                  do ipw = 1, n_pw
                    phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                    pval(ipw) = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                  end do
                  do ipw = 1, n_pw
                    rhs(row0 + ipw) = rhs(row0 + ipw) + conjg(pval(ipw)) * psi * vol_weight
                    basis_raw(row0 + ipw) = pval(ipw)
                  end do
                end do
                do iw = 1, nraw
                  rhs_density(iw) = rhs_density(iw) + conjg(basis_raw(iw)) * psi * vol_weight
                  if (.not. have_sdensity) then
                    do ipw = 1, nraw
                      S_density(iw, ipw) = S_density(iw, ipw) + conjg(basis_raw(iw)) * basis_raw(ipw) * vol_weight
                    end do
                  end if
                end do
              end do
            end do
          end do
          if (.not. have_sdensity) then
            call hermitize_matrix(S_density, nraw)
            call build_hermitian_pseudoinverse(S_density, nraw, 1.0d-8, S_density_inv, info_inv, smin, smax, ipw)
            if (info_inv /= 0) cycle
            have_sdensity = .true.
          end if
          if (info_inv == 0 .and. allocated(dg_frag%wpw_reduced_Sraw_build)) then
            c_raw(:) = matmul(S_density_inv, rhs_density(:))
            rhs(:) = matmul(dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local), c_raw(:))
          end if
          cproj(:) = matmul(s_inv, matmul( &
            conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), rhs(1:nraw)))
          s_tmp(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), cproj(:))
          norm_after = real(sum(conjg(cproj(:)) * s_tmp(:)), kind=8)
          norm_before = norm_after
          norm_loss_sum = norm_loss_sum + abs(norm_before - norm_after)
          norm_loss_max = max(norm_loss_max, abs(norm_before - norm_after))
          dg_frag%coef_wpw_self(1:nself, ist, ispin, i_local) = cproj(1:nself)
          if (nneigh > 0) dg_frag%coef_wpw_neighbor_reduced(1:nneigh, ist, ispin, i_local) = &
            cproj(nself+1:nred)
          did_project = .true.
        end do
      end do
      deallocate(rhs, cproj, s_inv, s_tmp, cw0, pval)
      deallocate(basis_raw, wval, S_density, S_density_inv, rhs_density, c_raw)
    end do
    if (did_project .and. dg_frag%id == 0) then
      write(*,'(1x,a,2(a,1pe12.4))') '[DG-WPW-RED-INIT-PROJ] S-metric self projection applied:', &
        ' residual_sum=', norm_loss_sum, ' residual_max=', norm_loss_max
    end if
  end subroutine initialize_wpw_reduced_self_projection

  subroutine diagnose_wpw_reduced_embed_local(dg_frag, istep)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in), optional :: istep

    integer :: step_use, i_local, ifrag, ispin, nred, nraw, nself, n_w, n_pw, n_pmix
    integer :: n_pfrag, pfrag_ids(7), axis, side, jfrag, nbf
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: ix, iy, iz, gx, gy, gz, bx, by, bz
    integer :: iw, ib, ipw, pidx, row0, jred, info_inv
    real(8) :: vol_weight, phase_arg, chi, grad_chi(3)
    real(8) :: local_norm, local_charge, coef_diff_norm, basis_leakage
    real(8) :: eval_min, eval_max
    complex(8) :: phase, psi
    complex(8), allocatable :: raw_vec(:), basis_raw(:), wval(:)
    complex(8), allocatable :: Sred_inv(:,:), c_back(:), c_diff(:), tmp_red(:)
    complex(8), allocatable :: raw_back(:), raw_resid(:), tmp_raw(:)
    logical :: found, bad
    logical, save :: warned_unready = .false.

    if (.not. wpw_reduced_embed_local_diag_enabled()) return
    step_use = -1
    if (present(istep)) step_use = istep
    if (.not. dg_frag%wpw_reduced_ready .or. .not. allocated(dg_frag%wpw_reduced_dim) .or. &
        .not. allocated(dg_frag%wpw_reduced_nself) .or. .not. allocated(dg_frag%wpw_reduced_nraw) .or. &
        .not. allocated(dg_frag%wpw_reduced_transform) .or. .not. allocated(dg_frag%wpw_reduced_Sraw_build) .or. &
        .not. allocated(dg_frag%wpw_reduced_S) .or. .not. allocated(dg_frag%global_wannier_local_coef) .or. &
        .not. allocated(dg_frag%global_wannier_local_nkeep) .or. &
        (.not. allocated(dg_frag%phi_frag) .and. .not. allocated(dg_frag%phi_frag_c))) then
      if (dg_frag%id == 0 .and. .not. warned_unready) then
        write(*,'(1x,a,6(a,i0),5(a,1pe12.4),a,l1,a)') &
          '[DG-WPW-RED-DIAG-EMBED-LOCAL]', &
          ' step=', step_use, ' rank=', dg_frag%id, ' frag=', 0, ' spin=', 0, ' nred=', 0, ' j=', 0, &
          ' local_norm=', 0.0d0, ' local_charge=', 0.0d0, &
          ' reproj_coef_diff_Snorm=', 0.0d0, ' rel_reproj_coef_diff=', 0.0d0, &
          ' basis_leakage_Snorm=', 0.0d0, ' bad=', .true., ' reason=reduced_not_ready'
        warned_unready = .true.
      end if
      return
    end if

    n_pw = dg_frag%n_plane_waves
    if (n_pw <= 0) return
    vol_weight = product(dg_frag%hgs(1:3)) / &
      product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    if (allocated(dg_frag%phi_frag_c)) then
      p_lb1 = lbound(dg_frag%phi_frag_c, 1); p_ub1 = ubound(dg_frag%phi_frag_c, 1)
      p_lb2 = lbound(dg_frag%phi_frag_c, 2); p_ub2 = ubound(dg_frag%phi_frag_c, 2)
      p_lb3 = lbound(dg_frag%phi_frag_c, 3); p_ub3 = ubound(dg_frag%phi_frag_c, 3)
    else
      p_lb1 = lbound(dg_frag%phi_frag, 1); p_ub1 = ubound(dg_frag%phi_frag, 1)
      p_lb2 = lbound(dg_frag%phi_frag, 2); p_ub2 = ubound(dg_frag%phi_frag, 2)
      p_lb3 = lbound(dg_frag%phi_frag, 3); p_ub3 = ubound(dg_frag%phi_frag, 3)
    end if

    do i_local = 1, size(dg_frag%wpw_reduced_dim)
      ifrag = dg_frag%ifrag_start + i_local - 1
      nred = dg_frag%wpw_reduced_dim(i_local)
      nself = dg_frag%wpw_reduced_nself(i_local)
      nraw = dg_frag%wpw_reduced_nraw(i_local)
      n_w = dg_frag%global_wannier_local_nkeep(i_local)
      n_pmix = nraw - n_w
      if (nred <= 0 .or. nraw <= 0 .or. nself <= 0 .or. n_w <= 0) cycle
      if (n_pmix <= 0) cycle
      if (nself /= n_w + n_pw) cycle
      n_pfrag = 1
      pfrag_ids(1) = ifrag
      do axis = 1, 3
        do side = -1, 1, 2
          jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle
          found = any(pfrag_ids(1:n_pfrag) == jfrag)
          if (.not. found .and. n_pfrag < size(pfrag_ids)) then
            n_pfrag = n_pfrag + 1
            pfrag_ids(n_pfrag) = jfrag
          end if
        end do
      end do
      if (nraw /= n_w + n_pfrag * n_pw) cycle
      n_pmix = nraw - n_w
      if (n_pmix <= 0) cycle
      nbf = min(dg_frag%n_basis(ifrag, 1), size(dg_frag%global_wannier_local_coef, 1))
      if (allocated(dg_frag%phi_frag_c)) then
        nbf = min(nbf, size(dg_frag%phi_frag_c, 4))
      else
        nbf = min(nbf, size(dg_frag%phi_frag, 4))
      end if
      if (nbf <= 0) cycle

      allocate(raw_vec(nraw), basis_raw(nraw), wval(n_w))
      allocate(Sred_inv(nred,nred), c_back(nred), c_diff(nred), tmp_red(nred))
      allocate(raw_back(nraw), raw_resid(nraw), tmp_raw(nraw))
      do ispin = 1, dg_frag%nspin
        call build_hermitian_inverse(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
          nred, Sred_inv, info_inv, eval_min, eval_max)
        do jred = 1, nred
          bad = .false.
          raw_vec(1:nraw) = dg_frag%wpw_reduced_transform(1:nraw, jred, i_local)
          local_norm = real(dg_frag%wpw_reduced_S(jred, jred, ispin, i_local), kind=8)
          local_charge = 0.0d0
          if (info_inv == 0) then
            tmp_raw(:) = matmul(dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local), raw_vec(:))
            c_back(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), tmp_raw))
            c_diff(:) = c_back(:)
            c_diff(jred) = c_diff(jred) - (1.0d0, 0.0d0)
            tmp_red(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), c_diff(:))
            coef_diff_norm = sqrt(max(0.0d0, real(sum(conjg(c_diff(:)) * tmp_red(:)), kind=8)))
            raw_back(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local), c_back(:))
            raw_resid(:) = raw_vec(:) - raw_back(:)
            tmp_raw(:) = matmul(dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local), raw_resid(:))
            basis_leakage = sqrt(max(0.0d0, real(sum(conjg(raw_resid(:)) * tmp_raw(:)), kind=8)))
          else
            coef_diff_norm = huge(1.0d0)
            basis_leakage = huge(1.0d0)
            bad = .true.
          end if

          do iz = 1, dg_frag%nxyz_domain(3, ifrag)
            gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
            bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bz < p_lb3 .or. bz > p_ub3) cycle
            do iy = 1, dg_frag%nxyz_domain(2, ifrag)
              gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
              by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              if (by < p_lb2 .or. by > p_ub2) cycle
              do ix = 1, dg_frag%nxyz_domain(1, ifrag)
                gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                if (bx < p_lb1 .or. bx > p_ub1) cycle
                basis_raw(:) = (0.0d0, 0.0d0)
                wval(:) = (0.0d0, 0.0d0)
                do iw = 1, n_w
                  do ib = 1, nbf
                    if (allocated(dg_frag%phi_frag_c)) then
                      wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                        dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                    else
                      wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                        cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                    end if
                  end do
                end do
                basis_raw(1:n_w) = wval(1:n_w)
                do pidx = 1, n_pfrag
                  row0 = n_w + (pidx - 1) * n_pw
                  call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                  do ipw = 1, n_pw
                    phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                    phase = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                    basis_raw(row0 + ipw) = phase
                  end do
                end do
                psi = sum(raw_vec(1:nraw) * basis_raw(1:nraw))
                local_charge = local_charge + abs(psi)**2 * vol_weight
              end do
            end do
          end do
          if (local_norm /= local_norm .or. local_charge /= local_charge .or. &
              coef_diff_norm /= coef_diff_norm .or. basis_leakage /= basis_leakage) bad = .true.
          write(*,'(1x,a,6(a,i0),5(a,1pe12.4),a,l1)') &
            '[DG-WPW-RED-DIAG-EMBED-LOCAL]', &
            ' step=', step_use, ' rank=', dg_frag%id, ' frag=', ifrag, ' spin=', ispin, &
            ' nred=', nred, ' j=', jred, &
            ' local_norm=', local_norm, &
            ' local_charge=', local_charge, &
            ' reproj_coef_diff_Snorm=', coef_diff_norm, &
            ' rel_reproj_coef_diff=', coef_diff_norm / max(sqrt(max(0.0d0, local_norm)), 1.0d-300), &
            ' basis_leakage_Snorm=', basis_leakage, &
            ' bad=', bad
        end do
      end do
      deallocate(raw_vec, basis_raw, wval, Sred_inv, c_back, c_diff, tmp_red, raw_back, raw_resid, tmp_raw)
    end do
  end subroutine diagnose_wpw_reduced_embed_local

  subroutine reconstruct_canonical_W_from_coef(dg_frag, ispin, i_local, nbf, n_w, cW, bx, by, bz, psi_grid)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin, i_local, nbf, n_w, bx, by, bz
    complex(8), intent(in) :: cW(n_w)
    complex(8), intent(out) :: psi_grid

    integer :: iw, ib
    complex(8) :: wval

    psi_grid = (0.0d0, 0.0d0)
    do iw = 1, n_w
      wval = (0.0d0, 0.0d0)
      do ib = 1, nbf
        if (allocated(dg_frag%phi_frag_c)) then
          wval = wval + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
            dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
        else
          wval = wval + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
            cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
        end if
      end do
      psi_grid = psi_grid + cW(iw) * wval
    end do
  end subroutine reconstruct_canonical_W_from_coef

  subroutine reconstruct_psi_from_C_can(dg_frag, ispin, i_local, nbf, n_w, n_pw, n_pfrag, pfrag_ids, &
      c_can, gx, gy, gz, bx, by, bz, psi_grid)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin, i_local, nbf, n_w, n_pw, n_pfrag
    integer, intent(in) :: pfrag_ids(n_pfrag)
    integer, intent(in) :: gx, gy, gz, bx, by, bz
    complex(8), intent(in) :: c_can(n_w + n_pw * n_pfrag)
    complex(8), intent(out) :: psi_grid

    integer :: ipw, pidx, row0
    real(8) :: phase_arg, chi, grad_chi(3)
    complex(8) :: phase, psi_w

    call reconstruct_canonical_W_from_coef(dg_frag, ispin, i_local, nbf, n_w, c_can(1:n_w), &
      bx, by, bz, psi_w)
    psi_grid = psi_w
    do pidx = 1, n_pfrag
      row0 = n_w + (pidx - 1) * n_pw
      call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
      do ipw = 1, n_pw
        phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                    dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                    dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
        phase = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
        psi_grid = psi_grid + c_can(row0 + ipw) * phase
      end do
    end do
  end subroutine reconstruct_psi_from_C_can

  subroutine diagnose_wpw_reduced_embed_prodcoef(dg_frag, istep)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in), optional :: istep

    integer :: step_use, i_local, ifrag, ispin, nred, nraw, nself, n_w, n_pw, n_pmix
    integer :: n_pfrag, pfrag_ids(7), axis, side, jfrag, nbf
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: ix, iy, iz, gx, gy, gz, bx, by, bz
    integer :: iw, ib, jb, ipw, jpw, pidx, row0, jred, dominant_raw_component
    integer :: p_neighbor_owner_frag, p_neighbor_source_frag, pidx_dom
    integer :: info_prod, info_raw, info_red, info_bpw, info_mixed_p, info_W, info_Wgrid, nkeep_tmp, nkeep_Wgrid
    real(8) :: vol_weight, phase_arg, chi, grad_chi(3)
    real(8) :: prodcoef_norm, local_norm, back_local_norm
    real(8) :: norm_in, norm_after_fragment_projection, norm_residual_after_fragment
    real(8) :: norm_after_bpw_perp_projection, norm_residual_after_bpw_perp
    real(8) :: rel_residual_after_fragment, rel_residual_after_bpw_perp
    real(8) :: norm_after_mixed_p_projection, norm_residual_after_mixed_p
    real(8) :: norm_after_neighbor_p_projection, norm_residual_after_neighbor_p
    real(8) :: norm_after_all_available_projection, norm_residual_after_all_available
    real(8) :: rel_residual_after_mixed_p, rel_residual_after_neighbor_p, rel_residual_after_all_available
    real(8) :: mixed_p_coef_norm, neighbor_p_coef_norm
    real(8) :: reduced_reproj_norm, reduced_overlap_with_input, missing_reduced_norm
    real(8) :: raw_W_norm, raw_P_self_norm, raw_P_neighbor_norm
    real(8) :: raw_BPW_perp_norm, raw_mixed_P_norm
    real(8) :: prod_W_norm, rawW_grid_norm, prodW_grid_norm
    real(8) :: rawW_to_rawW_diff, prodW_to_rawW_diff
    real(8) :: rawW_prodW_overlap, rawW_prodW_rel_diff
    real(8) :: raw_W_norm_euclid, raw_W_norm_reduced_metric
    real(8) :: S_W_diag_min, S_W_diag_max, S_W_diag_avg
    real(8) :: gridS_W_diag_min, gridS_W_diag_max, gridS_W_diag_avg
    real(8) :: diag_ratio_min, diag_ratio_max
    real(8) :: scale_to_match_S_min, scale_to_match_S_max, scale_to_match_S_avg
    real(8) :: plain_dot_diff, dV_dot_diff, Sinv_project_diff
    real(8) :: conjugate_basis_diff, transpose_transform_diff
    real(8) :: best_roundtrip_diff, best_norm_ratio
    real(8) :: normW_grid_norm, normW_to_rawW_diff, normW_red_reproj_diff
    real(8) :: normW_reduced_reproj_norm, normW_reduced_overlap, normW_missing_norm
    real(8) :: normW_basis_leakage, normW_grid_diag_min, normW_grid_diag_max, normW_grid_diag_avg
    real(8) :: canonicalW_grid_norm, canonicalW_to_rawW_diff
    real(8) :: canonicalW_red_reproj_diff, canonicalW_reduced_reproj_norm
    real(8) :: canonicalW_reduced_overlap, canonicalW_missing_norm, canonicalW_basis_leakage
    real(8) :: canonicalWP_red_reproj_diff, canonicalWP_reduced_reproj_norm
    real(8) :: canonicalWP_reduced_overlap, canonicalWP_missing_norm, canonicalWP_basis_leakage
    real(8) :: rawPself_to_rawPself_diff, rawPneighbor_to_rawPneighbor_diff
    real(8) :: Pself_grid_norm, Pneighbor_grid_norm
    real(8) :: eval_min_Wgrid, eval_max_Wgrid
    real(8) :: rel_prod_roundtrip, rel_reduced_reproj
    real(8) :: prod_roundtrip_diff, red_reproj_diff, basis_leakage, eval_min, eval_max
    complex(8) :: phase, psi_local, psi_frag, psi_bpw, psi_back, psi_resid, psi_rawW, psi_prodW, psi_normW
    complex(8) :: psi_Pself, psi_Pneighbor
    complex(8), allocatable :: raw_vec(:), raw_from_prod(:), raw_resid(:), raw_back(:), tmp_raw(:)
    complex(8), allocatable :: c_red_back(:), c_red_diff(:), tmp_red(:), Sred_inv(:,:)
    complex(8), allocatable :: basis_raw(:), wval(:), phi_val(:), pw_val(:), bpw_val(:)
    complex(8), allocatable :: S_prod(:,:), S_prod_inv(:,:), rhs_prod(:), coef_prod(:)
    complex(8), allocatable :: S_density(:,:), S_density_inv(:,:), rhs_density(:), rhs_hybrid(:)
    complex(8), allocatable :: S_frag_pw(:,:), C_frag_pw(:,:)
    complex(8), allocatable :: S_bpw(:,:), S_bpw_inv(:,:), rhs_bpw(:), coef_bpw(:)
    complex(8), allocatable :: S_mixed_p(:,:), S_mixed_p_inv(:,:), rhs_mixed_p(:), coef_mixed_p(:), tmp_mixed_p(:)
    complex(8), allocatable :: S_WW(:,:), S_WW_inv(:,:), raw_W_coef(:), prod_W_coef(:)
    complex(8), allocatable :: rhs_rawW(:), rhs_prodW(:), back_W_coef(:), diff_W_coef(:), tmp_W(:)
    complex(8), allocatable :: rhs_plainW(:), rhs_conjW(:), rhs_transposeW(:)
    complex(8), allocatable :: S_W_grid(:,:), S_W_grid_invsqrt(:,:), Wnorm_metric(:,:)
    complex(8), allocatable :: wnorm_val(:), rhs_normW(:), back_normW_coef(:), diff_normW_coef(:)
    complex(8), allocatable :: rhs_hybrid_normW(:), c_red_back_normW(:), c_red_diff_normW(:), tmp_red_normW(:)
    complex(8), allocatable :: rhs_hybrid_canonicalW(:), c_red_back_canonicalW(:)
    complex(8), allocatable :: c_red_diff_canonicalW(:), tmp_red_canonicalW(:)
    complex(8), allocatable :: rhs_hybrid_canonicalWP(:), c_red_back_canonicalWP(:)
    complex(8), allocatable :: c_red_diff_canonicalWP(:), tmp_red_canonicalWP(:)
    complex(8), allocatable :: raw_back_canonicalWP(:), raw_diff_canonicalWP(:), tmp_raw_canonicalWP(:)
    logical :: found, bad, ready, basis_order_match, metric_match, normalization_match
    logical :: neighbor_mapping_match
    logical :: wbasis_object_match, wbasis_index_range_match, wbasis_normalization_match, wbasis_bad
    character(len=64) :: reason, best_variant, wbasis_reason
    logical, save :: warned_unready = .false.

    if (.not. wpw_reduced_embed_prodcoef_diag_enabled()) return
    step_use = -1
    if (present(istep)) step_use = istep
    ready = dg_frag%wpw_reduced_ready .and. allocated(dg_frag%wpw_reduced_dim) .and. &
      allocated(dg_frag%wpw_reduced_nself) .and. allocated(dg_frag%wpw_reduced_nraw) .and. &
      allocated(dg_frag%wpw_reduced_transform) .and. allocated(dg_frag%wpw_reduced_Sraw_build) .and. &
      allocated(dg_frag%wpw_reduced_S) .and. allocated(dg_frag%global_wannier_local_coef) .and. &
      allocated(dg_frag%global_wannier_local_nkeep) .and. &
      (allocated(dg_frag%phi_frag) .or. allocated(dg_frag%phi_frag_c))
    if (.not. ready) then
      if (dg_frag%id == 0 .and. .not. warned_unready) then
        write(*,'(1x,a,6(a,i0),8(a,1pe12.4),a,l1,a)') &
          '[DG-WPW-RED-DIAG-EMBED-PRODCOEF]', &
          ' step=', step_use, ' rank=', dg_frag%id, ' frag=', 0, ' spin=', 0, ' nred=', 0, ' j=', 0, &
          ' prodcoef_norm=', 0.0d0, ' local_norm=', 0.0d0, ' back_local_norm=', 0.0d0, &
          ' prod_roundtrip_diff_Snorm=', 0.0d0, ' rel_prod_roundtrip_diff=', 0.0d0, &
          ' reduced_reproj_diff_Snorm=', 0.0d0, ' rel_reduced_reproj_diff=', 0.0d0, &
          ' basis_leakage_Snorm=', 0.0d0, ' bad=', .true., ' reason=reduced_not_ready'
        warned_unready = .true.
      end if
      return
    end if

    n_pw = dg_frag%n_plane_waves
    if (n_pw <= 0) return
    vol_weight = product(dg_frag%hgs(1:3)) / &
      product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    if (allocated(dg_frag%phi_frag_c)) then
      p_lb1 = lbound(dg_frag%phi_frag_c, 1); p_ub1 = ubound(dg_frag%phi_frag_c, 1)
      p_lb2 = lbound(dg_frag%phi_frag_c, 2); p_ub2 = ubound(dg_frag%phi_frag_c, 2)
      p_lb3 = lbound(dg_frag%phi_frag_c, 3); p_ub3 = ubound(dg_frag%phi_frag_c, 3)
    else
      p_lb1 = lbound(dg_frag%phi_frag, 1); p_ub1 = ubound(dg_frag%phi_frag, 1)
      p_lb2 = lbound(dg_frag%phi_frag, 2); p_ub2 = ubound(dg_frag%phi_frag, 2)
      p_lb3 = lbound(dg_frag%phi_frag, 3); p_ub3 = ubound(dg_frag%phi_frag, 3)
    end if

    do i_local = 1, size(dg_frag%wpw_reduced_dim)
      ifrag = dg_frag%ifrag_start + i_local - 1
      nred = dg_frag%wpw_reduced_dim(i_local)
      nself = dg_frag%wpw_reduced_nself(i_local)
      nraw = dg_frag%wpw_reduced_nraw(i_local)
      n_w = dg_frag%global_wannier_local_nkeep(i_local)
      if (nred <= 0 .or. nraw <= 0 .or. nself <= 0 .or. n_w <= 0) cycle
      n_pmix = nraw - n_w
      if (n_pmix <= 0) cycle
      if (nself /= n_w + n_pw) cycle
      n_pfrag = 1
      pfrag_ids(1) = ifrag
      do axis = 1, 3
        do side = -1, 1, 2
          jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle
          found = any(pfrag_ids(1:n_pfrag) == jfrag)
          if (.not. found .and. n_pfrag < size(pfrag_ids)) then
            n_pfrag = n_pfrag + 1
            pfrag_ids(n_pfrag) = jfrag
          end if
        end do
      end do
      if (nraw /= n_w + n_pfrag * n_pw) cycle
      nbf = min(dg_frag%n_basis(ifrag, 1), size(dg_frag%global_wannier_local_coef, 1))
      if (allocated(dg_frag%phi_frag_c)) then
        nbf = min(nbf, size(dg_frag%phi_frag_c, 4))
      else
        nbf = min(nbf, size(dg_frag%phi_frag, 4))
      end if
      if (nbf <= 0) cycle

      allocate(raw_vec(nraw), raw_from_prod(nraw), raw_resid(nraw), raw_back(nraw), tmp_raw(nraw))
      allocate(c_red_back(nred), c_red_diff(nred), tmp_red(nred), Sred_inv(nred,nred))
      allocate(basis_raw(nraw), wval(n_w), phi_val(nbf), pw_val(n_pw), bpw_val(n_pw))
      allocate(S_prod(nbf,nbf), S_prod_inv(nbf,nbf), rhs_prod(nbf), coef_prod(nbf))
      allocate(S_density(nraw,nraw), S_density_inv(nraw,nraw), rhs_density(nraw), rhs_hybrid(nraw))
      allocate(S_frag_pw(nbf,n_pw), C_frag_pw(nbf,n_pw))
      allocate(S_bpw(n_pw,n_pw), S_bpw_inv(n_pw,n_pw), rhs_bpw(n_pw), coef_bpw(n_pw))
      allocate(S_mixed_p(n_pmix,n_pmix), S_mixed_p_inv(n_pmix,n_pmix), rhs_mixed_p(n_pmix), &
        coef_mixed_p(n_pmix), tmp_mixed_p(n_pmix))
      allocate(S_WW(n_w,n_w), S_WW_inv(n_w,n_w), raw_W_coef(n_w), prod_W_coef(n_w), &
        rhs_rawW(n_w), rhs_prodW(n_w), back_W_coef(n_w), diff_W_coef(n_w), tmp_W(n_w))
      allocate(rhs_plainW(n_w), rhs_conjW(n_w), rhs_transposeW(n_w))
      allocate(S_W_grid(n_w,n_w), S_W_grid_invsqrt(n_w,n_w), Wnorm_metric(n_w,n_w))
      allocate(wnorm_val(n_w), rhs_normW(n_w), back_normW_coef(n_w), diff_normW_coef(n_w))
      allocate(rhs_hybrid_normW(nraw), c_red_back_normW(nred), c_red_diff_normW(nred), tmp_red_normW(nred))
      allocate(rhs_hybrid_canonicalW(nraw), c_red_back_canonicalW(nred), c_red_diff_canonicalW(nred), &
        tmp_red_canonicalW(nred))
      allocate(rhs_hybrid_canonicalWP(nraw), c_red_back_canonicalWP(nred), c_red_diff_canonicalWP(nred), &
        tmp_red_canonicalWP(nred), raw_back_canonicalWP(nraw), raw_diff_canonicalWP(nraw), &
        tmp_raw_canonicalWP(nraw))
      do ispin = 1, dg_frag%nspin
        S_prod(:, :) = (0.0d0, 0.0d0)
        S_density(:, :) = (0.0d0, 0.0d0)
        S_frag_pw(:, :) = (0.0d0, 0.0d0)
        do iz = 1, dg_frag%nxyz_domain(3, ifrag)
          gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
          bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
          if (bz < p_lb3 .or. bz > p_ub3) cycle
          do iy = 1, dg_frag%nxyz_domain(2, ifrag)
            gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
            by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            if (by < p_lb2 .or. by > p_ub2) cycle
            do ix = 1, dg_frag%nxyz_domain(1, ifrag)
              gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
              bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              if (bx < p_lb1 .or. bx > p_ub1) cycle
              phi_val(:) = (0.0d0, 0.0d0)
              do ib = 1, nbf
                if (allocated(dg_frag%phi_frag_c)) then
                  phi_val(ib) = dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                else
                  phi_val(ib) = cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                end if
              end do
              do ipw = 1, n_pw
                phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                            dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                            dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                pw_val(ipw) = exp(cmplx(0.0d0, phase_arg, kind=8))
              end do
              do ib = 1, nbf
                do jb = 1, nbf
                  S_prod(ib, jb) = S_prod(ib, jb) + conjg(phi_val(ib)) * phi_val(jb) * vol_weight
                end do
                do ipw = 1, n_pw
                  S_frag_pw(ib, ipw) = S_frag_pw(ib, ipw) + conjg(phi_val(ib)) * pw_val(ipw) * vol_weight
                end do
              end do
              basis_raw(:) = (0.0d0, 0.0d0)
              wval(:) = (0.0d0, 0.0d0)
              do iw = 1, n_w
                do ib = 1, nbf
                  if (allocated(dg_frag%phi_frag_c)) then
                    wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                      dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                  else
                    wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                      cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                  end if
                end do
              end do
              basis_raw(1:n_w) = wval(1:n_w)
              do pidx = 1, n_pfrag
                row0 = n_w + (pidx - 1) * n_pw
                call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                do ipw = 1, n_pw
                  phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                              dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                              dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                  basis_raw(row0 + ipw) = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                end do
              end do
              do iw = 1, nraw
                do ipw = 1, nraw
                  S_density(iw, ipw) = S_density(iw, ipw) + conjg(basis_raw(iw)) * basis_raw(ipw) * vol_weight
                end do
              end do
            end do
          end do
        end do
        call hermitize_matrix(S_prod, nbf)
        call hermitize_matrix(S_density, nraw)
        call build_hermitian_pseudoinverse(S_prod, nbf, 1.0d-8, S_prod_inv, info_prod, eval_min, eval_max, &
          nkeep_tmp)
        C_frag_pw(:, :) = matmul(S_prod_inv, S_frag_pw)
        S_bpw(:, :) = (0.0d0, 0.0d0)
        do iz = 1, dg_frag%nxyz_domain(3, ifrag)
          gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
          bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
          if (bz < p_lb3 .or. bz > p_ub3) cycle
          do iy = 1, dg_frag%nxyz_domain(2, ifrag)
            gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
            by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            if (by < p_lb2 .or. by > p_ub2) cycle
            do ix = 1, dg_frag%nxyz_domain(1, ifrag)
              gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
              bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              if (bx < p_lb1 .or. bx > p_ub1) cycle
              phi_val(:) = (0.0d0, 0.0d0)
              do ib = 1, nbf
                if (allocated(dg_frag%phi_frag_c)) then
                  phi_val(ib) = dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                else
                  phi_val(ib) = cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                end if
              end do
              do ipw = 1, n_pw
                phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                            dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                            dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                bpw_val(ipw) = exp(cmplx(0.0d0, phase_arg, kind=8)) - sum(C_frag_pw(1:nbf, ipw) * phi_val(1:nbf))
              end do
              do ipw = 1, n_pw
                do jpw = 1, n_pw
                  S_bpw(ipw, jpw) = S_bpw(ipw, jpw) + conjg(bpw_val(ipw)) * bpw_val(jpw) * vol_weight
                end do
              end do
            end do
          end do
        end do
        call hermitize_matrix(S_bpw, n_pw)
        call build_hermitian_pseudoinverse(S_bpw, n_pw, 1.0d-8, S_bpw_inv, info_bpw, eval_min, eval_max, &
          nkeep_tmp)
        info_raw = 0
        S_mixed_p(:, :) = dg_frag%wpw_reduced_Sraw_build(n_w+1:nraw, n_w+1:nraw, i_local)
        call hermitize_matrix(S_mixed_p, n_pmix)
        call build_hermitian_pseudoinverse(S_mixed_p, n_pmix, 1.0d-8, S_mixed_p_inv, info_mixed_p, &
          eval_min, eval_max, nkeep_tmp)
        S_WW(:, :) = dg_frag%wpw_reduced_Sraw_build(1:n_w, 1:n_w, i_local)
        call hermitize_matrix(S_WW, n_w)
        call build_hermitian_pseudoinverse(S_WW, n_w, 1.0d-8, S_WW_inv, info_W, eval_min, eval_max, &
          nkeep_tmp)
        S_W_grid(:, :) = S_density(1:n_w, 1:n_w)
        call hermitize_matrix(S_W_grid, n_w)
        call build_hermitian_inverse_sqrt(S_W_grid, n_w, 1.0d-8, S_W_grid_invsqrt, info_Wgrid, &
          eval_min_Wgrid, eval_max_Wgrid, nkeep_tmp)
        nkeep_Wgrid = nkeep_tmp
        if (info_Wgrid == 0) then
          Wnorm_metric(:, :) = matmul(conjg(transpose(S_W_grid_invsqrt)), matmul(S_W_grid, S_W_grid_invsqrt))
          call hermitize_matrix(Wnorm_metric, n_w)
        else
          Wnorm_metric(:, :) = (0.0d0, 0.0d0)
        end if
        S_W_diag_min = huge(1.0d0)
        S_W_diag_max = -huge(1.0d0)
        S_W_diag_avg = 0.0d0
        gridS_W_diag_min = huge(1.0d0)
        gridS_W_diag_max = -huge(1.0d0)
        gridS_W_diag_avg = 0.0d0
        normW_grid_diag_min = huge(1.0d0)
        normW_grid_diag_max = -huge(1.0d0)
        normW_grid_diag_avg = 0.0d0
        diag_ratio_min = huge(1.0d0)
        diag_ratio_max = -huge(1.0d0)
        scale_to_match_S_min = huge(1.0d0)
        scale_to_match_S_max = -huge(1.0d0)
        scale_to_match_S_avg = 0.0d0
        do iw = 1, n_w
          S_W_diag_min = min(S_W_diag_min, real(S_WW(iw, iw), kind=8))
          S_W_diag_max = max(S_W_diag_max, real(S_WW(iw, iw), kind=8))
          S_W_diag_avg = S_W_diag_avg + real(S_WW(iw, iw), kind=8)
          gridS_W_diag_min = min(gridS_W_diag_min, real(S_density(iw, iw), kind=8))
          gridS_W_diag_max = max(gridS_W_diag_max, real(S_density(iw, iw), kind=8))
          gridS_W_diag_avg = gridS_W_diag_avg + real(S_density(iw, iw), kind=8)
          normW_grid_diag_min = min(normW_grid_diag_min, real(Wnorm_metric(iw, iw), kind=8))
          normW_grid_diag_max = max(normW_grid_diag_max, real(Wnorm_metric(iw, iw), kind=8))
          normW_grid_diag_avg = normW_grid_diag_avg + real(Wnorm_metric(iw, iw), kind=8)
          diag_ratio_min = min(diag_ratio_min, real(S_density(iw, iw), kind=8) / &
            max(abs(real(S_WW(iw, iw), kind=8)), 1.0d-300))
          diag_ratio_max = max(diag_ratio_max, real(S_density(iw, iw), kind=8) / &
            max(abs(real(S_WW(iw, iw), kind=8)), 1.0d-300))
          scale_to_match_S_min = min(scale_to_match_S_min, sqrt(abs(real(S_WW(iw, iw), kind=8)) / &
            max(abs(real(S_density(iw, iw), kind=8)), 1.0d-300)))
          scale_to_match_S_max = max(scale_to_match_S_max, sqrt(abs(real(S_WW(iw, iw), kind=8)) / &
            max(abs(real(S_density(iw, iw), kind=8)), 1.0d-300)))
          scale_to_match_S_avg = scale_to_match_S_avg + sqrt(abs(real(S_WW(iw, iw), kind=8)) / &
            max(abs(real(S_density(iw, iw), kind=8)), 1.0d-300))
        end do
        S_W_diag_avg = S_W_diag_avg / max(1, n_w)
        gridS_W_diag_avg = gridS_W_diag_avg / max(1, n_w)
        normW_grid_diag_avg = normW_grid_diag_avg / max(1, n_w)
        scale_to_match_S_avg = scale_to_match_S_avg / max(1, n_w)
        wbasis_index_range_match = (n_w <= size(dg_frag%global_wannier_local_coef, 2)) .and. &
          (ispin <= size(dg_frag%global_wannier_local_coef, 3)) .and. &
          (i_local <= size(dg_frag%global_wannier_local_coef, 4)) .and. &
          (nbf <= size(dg_frag%global_wannier_local_coef, 1))
        wbasis_object_match = max(abs(diag_ratio_min - 1.0d0), abs(diag_ratio_max - 1.0d0)) < 1.0d-8
        wbasis_normalization_match = (abs(scale_to_match_S_avg - 1.0d0) < 1.0d-8) .and. &
          (abs(scale_to_match_S_max - scale_to_match_S_min) < 1.0d-8)
        wbasis_bad = (.not. wbasis_index_range_match) .or. &
          S_W_diag_min /= S_W_diag_min .or. S_W_diag_max /= S_W_diag_max .or. &
          gridS_W_diag_min /= gridS_W_diag_min .or. gridS_W_diag_max /= gridS_W_diag_max .or. &
          scale_to_match_S_min /= scale_to_match_S_min .or. scale_to_match_S_max /= scale_to_match_S_max
        wbasis_reason = 'ok'
        if (.not. wbasis_index_range_match) then
          wbasis_reason = 'index_range_mismatch'
        else if (.not. wbasis_object_match .and. wbasis_normalization_match) then
          wbasis_reason = 'basis_object_mismatch'
        else if (.not. wbasis_normalization_match) then
          wbasis_reason = 'normalization_mismatch'
        end if
        write(*,'(1x,a,5(a,i0),9(a,1pe12.4),4(a,l1),a,a)') &
          '[DG-WPW-RED-DIAG-WBASIS-IDENTITY]', &
          ' step=', step_use, ' frag=', ifrag, ' spin=', ispin, ' nW=', n_w, &
          ' nW_grid_keep=', nkeep_Wgrid, &
          ' S_W_diag_min=', S_W_diag_min, &
          ' S_W_diag_max=', S_W_diag_max, &
          ' S_W_diag_avg=', S_W_diag_avg, &
          ' gridW_norm_min=', gridS_W_diag_min, &
          ' gridW_norm_max=', gridS_W_diag_max, &
          ' gridW_norm_avg=', gridS_W_diag_avg, &
          ' scale_to_match_S_min=', scale_to_match_S_min, &
          ' scale_to_match_S_max=', scale_to_match_S_max, &
          ' scale_to_match_S_avg=', scale_to_match_S_avg, &
          ' basis_object_match=', wbasis_object_match, &
          ' index_range_match=', wbasis_index_range_match, &
          ' normalization_match=', wbasis_normalization_match, &
          ' bad=', wbasis_bad, ' reason=', trim(wbasis_reason)
        call build_hermitian_inverse(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
          nred, Sred_inv, info_red, eval_min, eval_max)
        do jred = 1, nred
          bad = .false.
          reason = 'ok'
          raw_vec(1:nraw) = dg_frag%wpw_reduced_transform(1:nraw, jred, i_local)
          dominant_raw_component = maxloc(abs(raw_vec(1:nraw)), dim=1)
          raw_W_norm = 0.0d0
          raw_P_self_norm = 0.0d0
          raw_P_neighbor_norm = 0.0d0
          raw_BPW_perp_norm = 0.0d0
          raw_mixed_P_norm = 0.0d0
          prod_W_norm = 0.0d0
          rawW_grid_norm = 0.0d0
          prodW_grid_norm = 0.0d0
          rawW_to_rawW_diff = 0.0d0
          prodW_to_rawW_diff = 0.0d0
          rawW_prodW_overlap = 0.0d0
          rawW_prodW_rel_diff = 0.0d0
          raw_W_norm_euclid = 0.0d0
          raw_W_norm_reduced_metric = raw_W_norm
          plain_dot_diff = 0.0d0
          dV_dot_diff = 0.0d0
          Sinv_project_diff = 0.0d0
          conjugate_basis_diff = 0.0d0
          transpose_transform_diff = 0.0d0
          best_roundtrip_diff = huge(1.0d0)
          best_norm_ratio = 0.0d0
          best_variant = 'none'
          normW_grid_norm = 0.0d0
          normW_to_rawW_diff = huge(1.0d0)
          normW_red_reproj_diff = huge(1.0d0)
          normW_reduced_reproj_norm = 0.0d0
          normW_reduced_overlap = 0.0d0
          normW_missing_norm = huge(1.0d0)
          normW_basis_leakage = huge(1.0d0)
          canonicalW_grid_norm = 0.0d0
          canonicalW_to_rawW_diff = 0.0d0
          canonicalW_red_reproj_diff = huge(1.0d0)
          canonicalW_reduced_reproj_norm = 0.0d0
          canonicalW_reduced_overlap = 0.0d0
          canonicalW_missing_norm = huge(1.0d0)
          canonicalW_basis_leakage = huge(1.0d0)
          canonicalWP_red_reproj_diff = huge(1.0d0)
          canonicalWP_reduced_reproj_norm = 0.0d0
          canonicalWP_reduced_overlap = 0.0d0
          canonicalWP_missing_norm = huge(1.0d0)
          canonicalWP_basis_leakage = huge(1.0d0)
          rawPself_to_rawPself_diff = huge(1.0d0)
          rawPneighbor_to_rawPneighbor_diff = huge(1.0d0)
          Pself_grid_norm = 0.0d0
          Pneighbor_grid_norm = 0.0d0
          p_neighbor_owner_frag = ifrag
          p_neighbor_source_frag = 0
          neighbor_mapping_match = .true.
          if (dominant_raw_component > n_w) then
            pidx_dom = (dominant_raw_component - n_w - 1) / n_pw + 1
            if (pidx_dom >= 1 .and. pidx_dom <= n_pfrag) p_neighbor_source_frag = pfrag_ids(pidx_dom)
            if (dominant_raw_component > n_w + n_pw) then
              neighbor_mapping_match = (p_neighbor_source_frag > 0 .and. p_neighbor_source_frag /= ifrag)
            end if
          end if
          raw_W_coef(:) = raw_vec(1:n_w)
          raw_W_norm_euclid = real(sum(conjg(raw_W_coef(:)) * raw_W_coef(:)), kind=8)
          prod_W_coef(:) = (0.0d0, 0.0d0)
          rhs_rawW(:) = (0.0d0, 0.0d0)
          rhs_prodW(:) = (0.0d0, 0.0d0)
          rhs_plainW(:) = (0.0d0, 0.0d0)
          rhs_conjW(:) = (0.0d0, 0.0d0)
          rhs_transposeW(:) = (0.0d0, 0.0d0)
          rhs_normW(:) = (0.0d0, 0.0d0)
          basis_order_match = .false.
          metric_match = .false.
          normalization_match = .false.
          tmp_raw(:) = (0.0d0, 0.0d0)
          if (n_w > 0) then
            tmp_raw(1:n_w) = matmul(dg_frag%wpw_reduced_Sraw_build(1:n_w, 1:n_w, i_local), raw_vec(1:n_w))
            raw_W_norm = real(sum(conjg(raw_vec(1:n_w)) * tmp_raw(1:n_w)), kind=8)
            raw_W_norm_reduced_metric = raw_W_norm
            canonicalW_grid_norm = raw_W_norm
          end if
          if (n_pw > 0 .and. n_w + n_pw <= nraw) then
            tmp_raw(n_w+1:n_w+n_pw) = matmul( &
              dg_frag%wpw_reduced_Sraw_build(n_w+1:n_w+n_pw, n_w+1:n_w+n_pw, i_local), &
              raw_vec(n_w+1:n_w+n_pw))
            raw_P_self_norm = real(sum(conjg(raw_vec(n_w+1:n_w+n_pw)) * tmp_raw(n_w+1:n_w+n_pw)), kind=8)
          end if
          if (n_w + n_pw < nraw) then
            tmp_raw(n_w+n_pw+1:nraw) = matmul( &
              dg_frag%wpw_reduced_Sraw_build(n_w+n_pw+1:nraw, n_w+n_pw+1:nraw, i_local), &
              raw_vec(n_w+n_pw+1:nraw))
            raw_P_neighbor_norm = real(sum(conjg(raw_vec(n_w+n_pw+1:nraw)) * &
              tmp_raw(n_w+n_pw+1:nraw)), kind=8)
          end if
          local_norm = real(dg_frag%wpw_reduced_S(jred, jred, ispin, i_local), kind=8)
          prodcoef_norm = 0.0d0
          back_local_norm = 0.0d0
          norm_in = 0.0d0
          norm_after_fragment_projection = 0.0d0
          norm_residual_after_fragment = 0.0d0
          norm_after_bpw_perp_projection = 0.0d0
          norm_residual_after_bpw_perp = 0.0d0
          norm_after_mixed_p_projection = 0.0d0
          norm_residual_after_mixed_p = 0.0d0
          norm_after_neighbor_p_projection = 0.0d0
          norm_residual_after_neighbor_p = 0.0d0
          norm_after_all_available_projection = 0.0d0
          norm_residual_after_all_available = 0.0d0
          rel_residual_after_fragment = 0.0d0
          rel_residual_after_bpw_perp = 0.0d0
          rel_residual_after_mixed_p = 0.0d0
          rel_residual_after_neighbor_p = 0.0d0
          rel_residual_after_all_available = 0.0d0
          mixed_p_coef_norm = 0.0d0
          neighbor_p_coef_norm = 0.0d0
          reduced_reproj_norm = 0.0d0
          reduced_overlap_with_input = 0.0d0
          missing_reduced_norm = 0.0d0
          prod_roundtrip_diff = huge(1.0d0)
          red_reproj_diff = huge(1.0d0)
          basis_leakage = huge(1.0d0)
          rel_prod_roundtrip = 0.0d0
          rel_reduced_reproj = 0.0d0
          if (info_prod /= 0 .or. info_bpw /= 0 .or. info_raw /= 0 .or. info_red /= 0 .or. &
              info_W /= 0 .or. info_Wgrid /= 0) then
            bad = .true.
            reason = 'metric_inverse_failed'
          else
            rhs_prod(:) = (0.0d0, 0.0d0)
            local_norm = 0.0d0
            do iz = 1, dg_frag%nxyz_domain(3, ifrag)
              gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
              bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              if (bz < p_lb3 .or. bz > p_ub3) cycle
              do iy = 1, dg_frag%nxyz_domain(2, ifrag)
                gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
                by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                if (by < p_lb2 .or. by > p_ub2) cycle
                do ix = 1, dg_frag%nxyz_domain(1, ifrag)
                  gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                  bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                  if (bx < p_lb1 .or. bx > p_ub1) cycle
                  basis_raw(:) = (0.0d0, 0.0d0)
                  wval(:) = (0.0d0, 0.0d0)
                  do iw = 1, n_w
                    do ib = 1, nbf
                      if (allocated(dg_frag%phi_frag_c)) then
                        wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                          dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                      else
                        wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                          cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                      end if
                    end do
                  end do
                  basis_raw(1:n_w) = wval(1:n_w)
                  wnorm_val(:) = (0.0d0, 0.0d0)
                  do iw = 1, n_w
                    do jb = 1, n_w
                      wnorm_val(jb) = wnorm_val(jb) + wval(iw) * S_W_grid_invsqrt(iw, jb)
                    end do
                  end do
                  do pidx = 1, n_pfrag
                    row0 = n_w + (pidx - 1) * n_pw
                    call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                    do ipw = 1, n_pw
                      phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                  dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                  dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                      basis_raw(row0 + ipw) = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                    end do
                  end do
                  psi_local = sum(raw_vec(1:nraw) * basis_raw(1:nraw))
                  local_norm = local_norm + abs(psi_local)**2 * vol_weight
                  do ib = 1, nbf
                    if (allocated(dg_frag%phi_frag_c)) then
                      rhs_prod(ib) = rhs_prod(ib) + conjg(dg_frag%phi_frag_c(bx, by, bz, ib, i_local)) * &
                        psi_local * vol_weight
                    else
                      rhs_prod(ib) = rhs_prod(ib) + &
                        dg_frag%phi_frag(bx, by, bz, ib, i_local) * psi_local * vol_weight
                    end if
                  end do
                end do
              end do
            end do
            norm_in = local_norm
            coef_prod(:) = matmul(S_prod_inv, rhs_prod)
            prodcoef_norm = real(sum(conjg(coef_prod(:)) * matmul(S_prod, coef_prod(:))), kind=8)
            prod_W_coef(:) = matmul(conjg(transpose( &
              dg_frag%global_wannier_local_coef(1:nbf, 1:n_w, ispin, i_local))), rhs_prod(1:nbf))
            prod_W_norm = real(sum(conjg(prod_W_coef(:)) * matmul(S_WW, prod_W_coef(:))), kind=8)
            diff_W_coef(:) = raw_W_coef(:) - prod_W_coef(:)
            rawW_prodW_rel_diff = sqrt(max(0.0d0, real(sum(conjg(diff_W_coef(:)) * &
              matmul(S_WW, diff_W_coef(:))), kind=8))) / max(sqrt(max(0.0d0, raw_W_norm)), 1.0d-300)
            rawW_prodW_overlap = abs(sum(conjg(raw_W_coef(:)) * matmul(S_WW, prod_W_coef(:))))
            rhs_bpw(:) = (0.0d0, 0.0d0)
            rhs_density(:) = (0.0d0, 0.0d0)
            rhs_hybrid(:) = (0.0d0, 0.0d0)
            do iz = 1, dg_frag%nxyz_domain(3, ifrag)
              gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
              bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              if (bz < p_lb3 .or. bz > p_ub3) cycle
              do iy = 1, dg_frag%nxyz_domain(2, ifrag)
                gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
                by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                if (by < p_lb2 .or. by > p_ub2) cycle
                do ix = 1, dg_frag%nxyz_domain(1, ifrag)
                  gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                  bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                  if (bx < p_lb1 .or. bx > p_ub1) cycle
                  phi_val(:) = (0.0d0, 0.0d0)
                  do ib = 1, nbf
                    if (allocated(dg_frag%phi_frag_c)) then
                      phi_val(ib) = dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                    else
                      phi_val(ib) = cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                    end if
                  end do
                  do ipw = 1, n_pw
                    phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                    bpw_val(ipw) = exp(cmplx(0.0d0, phase_arg, kind=8)) - sum(C_frag_pw(1:nbf, ipw) * phi_val(1:nbf))
                  end do
                  basis_raw(:) = (0.0d0, 0.0d0)
                  wval(:) = (0.0d0, 0.0d0)
                  do iw = 1, n_w
                    do ib = 1, nbf
                      if (allocated(dg_frag%phi_frag_c)) then
                        wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                          dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                      else
                        wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                          cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                      end if
                    end do
                  end do
                  basis_raw(1:n_w) = wval(1:n_w)
                  do pidx = 1, n_pfrag
                    row0 = n_w + (pidx - 1) * n_pw
                    call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                    do ipw = 1, n_pw
                      phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                  dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                  dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                      basis_raw(row0 + ipw) = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                    end do
                  end do
                  psi_local = sum(raw_vec(1:nraw) * basis_raw(1:nraw))
                  psi_Pself = (0.0d0, 0.0d0)
                  psi_Pneighbor = (0.0d0, 0.0d0)
                  if (n_pw > 0) psi_Pself = sum(raw_vec(n_w+1:n_w+n_pw) * basis_raw(n_w+1:n_w+n_pw))
                  if (n_w + n_pw < nraw) psi_Pneighbor = &
                    sum(raw_vec(n_w+n_pw+1:nraw) * basis_raw(n_w+n_pw+1:nraw))
                  call reconstruct_canonical_W_from_coef(dg_frag, ispin, i_local, nbf, n_w, raw_W_coef, &
                    bx, by, bz, psi_rawW)
                  psi_prodW = sum(prod_W_coef(1:n_w) * basis_raw(1:n_w))
                  psi_normW = sum(raw_W_coef(1:n_w) * wnorm_val(1:n_w))
                  rawW_grid_norm = rawW_grid_norm + abs(psi_rawW)**2 * vol_weight
                  prodW_grid_norm = prodW_grid_norm + abs(psi_prodW)**2 * vol_weight
                  normW_grid_norm = normW_grid_norm + abs(psi_normW)**2 * vol_weight
                  Pself_grid_norm = Pself_grid_norm + abs(psi_Pself)**2 * vol_weight
                  Pneighbor_grid_norm = Pneighbor_grid_norm + abs(psi_Pneighbor)**2 * vol_weight
                  rhs_rawW(1:n_w) = rhs_rawW(1:n_w) + conjg(basis_raw(1:n_w)) * psi_rawW * vol_weight
                  rhs_plainW(1:n_w) = rhs_plainW(1:n_w) + conjg(basis_raw(1:n_w)) * psi_rawW
                  rhs_conjW(1:n_w) = rhs_conjW(1:n_w) + basis_raw(1:n_w) * psi_rawW * vol_weight
                  rhs_transposeW(1:n_w) = rhs_transposeW(1:n_w) + basis_raw(1:n_w) * conjg(psi_rawW) * vol_weight
                  rhs_prodW(1:n_w) = rhs_prodW(1:n_w) + conjg(basis_raw(1:n_w)) * psi_prodW * vol_weight
                  rhs_normW(1:n_w) = rhs_normW(1:n_w) + conjg(wnorm_val(1:n_w)) * psi_normW * vol_weight
	                  psi_frag = sum(coef_prod(1:nbf) * phi_val(1:nbf))
                  psi_resid = psi_local - psi_frag
                  norm_after_fragment_projection = norm_after_fragment_projection + abs(psi_frag)**2 * vol_weight
                  norm_residual_after_fragment = norm_residual_after_fragment + abs(psi_resid)**2 * vol_weight
                  rhs_bpw(1:n_pw) = rhs_bpw(1:n_pw) + conjg(bpw_val(1:n_pw)) * psi_resid * vol_weight
                end do
              end do
            end do
            coef_bpw(:) = matmul(S_bpw_inv, rhs_bpw)
            raw_BPW_perp_norm = real(sum(conjg(coef_bpw(:)) * matmul(S_bpw, coef_bpw(:))), kind=8)
            prodcoef_norm = prodcoef_norm + raw_BPW_perp_norm
            rhs_mixed_p(:) = (0.0d0, 0.0d0)
            do iz = 1, dg_frag%nxyz_domain(3, ifrag)
              gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
              bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              if (bz < p_lb3 .or. bz > p_ub3) cycle
              do iy = 1, dg_frag%nxyz_domain(2, ifrag)
                gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
                by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                if (by < p_lb2 .or. by > p_ub2) cycle
                do ix = 1, dg_frag%nxyz_domain(1, ifrag)
                  gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                  bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                  if (bx < p_lb1 .or. bx > p_ub1) cycle
                  phi_val(:) = (0.0d0, 0.0d0)
                  do ib = 1, nbf
                    if (allocated(dg_frag%phi_frag_c)) then
                      phi_val(ib) = dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                    else
                      phi_val(ib) = cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                    end if
                  end do
                  do ipw = 1, n_pw
                    phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                    bpw_val(ipw) = exp(cmplx(0.0d0, phase_arg, kind=8)) - sum(C_frag_pw(1:nbf, ipw) * phi_val(1:nbf))
                  end do
                  basis_raw(:) = (0.0d0, 0.0d0)
                  wval(:) = (0.0d0, 0.0d0)
                  do iw = 1, n_w
                    do ib = 1, nbf
                      if (allocated(dg_frag%phi_frag_c)) then
                        wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                          dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                      else
                        wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                          cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                      end if
                    end do
                  end do
                  basis_raw(1:n_w) = wval(1:n_w)
                  do pidx = 1, n_pfrag
                    row0 = n_w + (pidx - 1) * n_pw
                    call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                    do ipw = 1, n_pw
                      phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                  dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                  dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                      basis_raw(row0 + ipw) = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                    end do
                  end do
                  psi_local = sum(raw_vec(1:nraw) * basis_raw(1:nraw))
                  psi_frag = sum(coef_prod(1:nbf) * phi_val(1:nbf))
                  psi_bpw = sum(coef_bpw(1:n_pw) * bpw_val(1:n_pw))
                  psi_resid = psi_local - psi_frag - psi_bpw
                  rhs_mixed_p(1:n_pmix) = rhs_mixed_p(1:n_pmix) + &
                    conjg(basis_raw(n_w+1:nraw)) * psi_resid * vol_weight
                end do
              end do
            end do
            if (info_mixed_p == 0) then
              coef_mixed_p(:) = matmul(S_mixed_p_inv, rhs_mixed_p)
            else
              coef_mixed_p(:) = (0.0d0, 0.0d0)
            end if
            mixed_p_coef_norm = real(sum(conjg(coef_mixed_p(:)) * matmul(S_mixed_p, coef_mixed_p(:))), kind=8)
            raw_mixed_P_norm = mixed_p_coef_norm
            tmp_mixed_p(:) = (0.0d0, 0.0d0)
            if (n_pmix > n_pw) tmp_mixed_p(n_pw+1:n_pmix) = coef_mixed_p(n_pw+1:n_pmix)
            neighbor_p_coef_norm = real(sum(conjg(tmp_mixed_p(:)) * matmul(S_mixed_p, tmp_mixed_p(:))), kind=8)
            prodcoef_norm = prodcoef_norm + mixed_p_coef_norm
            rhs_density(:) = (0.0d0, 0.0d0)
            rhs_hybrid(:) = (0.0d0, 0.0d0)
            rhs_prod(:) = (0.0d0, 0.0d0)
            do iz = 1, dg_frag%nxyz_domain(3, ifrag)
              gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
              bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
              if (bz < p_lb3 .or. bz > p_ub3) cycle
              do iy = 1, dg_frag%nxyz_domain(2, ifrag)
                gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
                by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                if (by < p_lb2 .or. by > p_ub2) cycle
                do ix = 1, dg_frag%nxyz_domain(1, ifrag)
                  gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                  bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                  if (bx < p_lb1 .or. bx > p_ub1) cycle
                  phi_val(:) = (0.0d0, 0.0d0)
                  do ib = 1, nbf
                    if (allocated(dg_frag%phi_frag_c)) then
                      phi_val(ib) = dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                    else
                      phi_val(ib) = cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                    end if
                  end do
                  do ipw = 1, n_pw
                    phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                    bpw_val(ipw) = exp(cmplx(0.0d0, phase_arg, kind=8)) - sum(C_frag_pw(1:nbf, ipw) * phi_val(1:nbf))
                  end do
                  basis_raw(:) = (0.0d0, 0.0d0)
                  wval(:) = (0.0d0, 0.0d0)
                  do iw = 1, n_w
                    do ib = 1, nbf
                      if (allocated(dg_frag%phi_frag_c)) then
                        wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                          dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                      else
                        wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                          cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                      end if
                    end do
                  end do
                  basis_raw(1:n_w) = wval(1:n_w)
                  do pidx = 1, n_pfrag
                    row0 = n_w + (pidx - 1) * n_pw
                    call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                    do ipw = 1, n_pw
                      phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                  dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                  dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                      basis_raw(row0 + ipw) = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                    end do
                  end do
                  psi_local = sum(raw_vec(1:nraw) * basis_raw(1:nraw))
                  psi_frag = sum(coef_prod(1:nbf) * phi_val(1:nbf))
                  psi_bpw = sum(coef_bpw(1:n_pw) * bpw_val(1:n_pw))
                  psi_back = psi_frag + psi_bpw
                  norm_after_bpw_perp_projection = norm_after_bpw_perp_projection + abs(psi_back)**2 * vol_weight
                  norm_residual_after_bpw_perp = norm_residual_after_bpw_perp + abs(psi_local - psi_back)**2 * vol_weight
                  if (n_pw > 0) psi_back = psi_back + &
                    sum(coef_mixed_p(1:min(n_pw, n_pmix)) * basis_raw(n_w+1:n_w+min(n_pw, n_pmix)))
                  norm_after_mixed_p_projection = norm_after_mixed_p_projection + abs(psi_back)**2 * vol_weight
                  norm_residual_after_mixed_p = norm_residual_after_mixed_p + abs(psi_local - psi_back)**2 * vol_weight
                  if (n_pmix > n_pw) psi_back = psi_back + &
                    sum(coef_mixed_p(n_pw+1:n_pmix) * basis_raw(n_w+n_pw+1:nraw))
                  norm_after_neighbor_p_projection = norm_after_neighbor_p_projection + abs(psi_back)**2 * vol_weight
                  norm_residual_after_neighbor_p = norm_residual_after_neighbor_p + abs(psi_local - psi_back)**2 * vol_weight
                  back_local_norm = back_local_norm + abs(psi_back)**2 * vol_weight
                  norm_after_all_available_projection = norm_after_all_available_projection + abs(psi_back)**2 * vol_weight
                  norm_residual_after_all_available = norm_residual_after_all_available + abs(psi_local - psi_back)**2 * vol_weight
                  rhs_density(1:nraw) = rhs_density(1:nraw) + conjg(basis_raw(1:nraw)) * psi_back * vol_weight
                  rhs_hybrid(n_w+1:nraw) = rhs_hybrid(n_w+1:nraw) + conjg(basis_raw(n_w+1:nraw)) * psi_back * vol_weight
                  rhs_prod(1:nbf) = rhs_prod(1:nbf) + conjg(phi_val(1:nbf)) * psi_back * vol_weight
                end do
              end do
            end do
            rhs_hybrid(1:n_w) = matmul(conjg(transpose( &
              dg_frag%global_wannier_local_coef(1:nbf, 1:n_w, ispin, i_local))), rhs_prod(1:nbf))
            rhs_hybrid_canonicalW(:) = rhs_hybrid(:)
            rhs_hybrid_canonicalW(1:n_w) = raw_W_coef(1:n_w)
            back_W_coef(:) = matmul(S_WW_inv, rhs_rawW(:))
            diff_W_coef(:) = back_W_coef(:) - raw_W_coef(:)
            tmp_W(:) = matmul(S_WW, diff_W_coef(:))
            rawW_to_rawW_diff = sqrt(max(0.0d0, real(sum(conjg(diff_W_coef(:)) * tmp_W(:)), kind=8)))
            Sinv_project_diff = rawW_to_rawW_diff
            diff_W_coef(:) = rhs_plainW(:) - raw_W_coef(:)
            tmp_W(:) = matmul(S_WW, diff_W_coef(:))
            plain_dot_diff = sqrt(max(0.0d0, real(sum(conjg(diff_W_coef(:)) * tmp_W(:)), kind=8)))
            diff_W_coef(:) = rhs_rawW(:) - raw_W_coef(:)
            tmp_W(:) = matmul(S_WW, diff_W_coef(:))
            dV_dot_diff = sqrt(max(0.0d0, real(sum(conjg(diff_W_coef(:)) * tmp_W(:)), kind=8)))
            back_W_coef(:) = matmul(S_WW_inv, rhs_conjW(:))
            diff_W_coef(:) = back_W_coef(:) - raw_W_coef(:)
            tmp_W(:) = matmul(S_WW, diff_W_coef(:))
            conjugate_basis_diff = sqrt(max(0.0d0, real(sum(conjg(diff_W_coef(:)) * tmp_W(:)), kind=8)))
            back_W_coef(:) = conjg(matmul(S_WW_inv, rhs_transposeW(:)))
            diff_W_coef(:) = back_W_coef(:) - raw_W_coef(:)
            tmp_W(:) = matmul(S_WW, diff_W_coef(:))
            transpose_transform_diff = sqrt(max(0.0d0, real(sum(conjg(diff_W_coef(:)) * tmp_W(:)), kind=8)))
            back_normW_coef(:) = rhs_normW(:)
            diff_normW_coef(:) = back_normW_coef(:) - raw_W_coef(:)
            normW_to_rawW_diff = sqrt(max(0.0d0, real(sum(conjg(diff_normW_coef(:)) * diff_normW_coef(:)), kind=8)))
            best_variant = 'plain_dot'
            best_roundtrip_diff = plain_dot_diff
            if (dV_dot_diff < best_roundtrip_diff) then
              best_variant = 'dV_dot'
              best_roundtrip_diff = dV_dot_diff
            end if
            if (Sinv_project_diff < best_roundtrip_diff) then
              best_variant = 'Sinv_project'
              best_roundtrip_diff = Sinv_project_diff
            end if
            if (conjugate_basis_diff < best_roundtrip_diff) then
              best_variant = 'conjugate_basis'
              best_roundtrip_diff = conjugate_basis_diff
            end if
            if (transpose_transform_diff < best_roundtrip_diff) then
              best_variant = 'transpose_transform'
              best_roundtrip_diff = transpose_transform_diff
            end if
            best_norm_ratio = rawW_grid_norm / max(abs(raw_W_norm), 1.0d-300)
            back_W_coef(:) = matmul(S_WW_inv, rhs_prodW(:))
            diff_W_coef(:) = back_W_coef(:) - raw_W_coef(:)
            tmp_W(:) = matmul(S_WW, diff_W_coef(:))
            prodW_to_rawW_diff = sqrt(max(0.0d0, real(sum(conjg(diff_W_coef(:)) * tmp_W(:)), kind=8)))
            basis_order_match = rawW_to_rawW_diff <= 1.0d-8
            metric_match = prodW_to_rawW_diff <= 1.0d-6
            normalization_match = abs(rawW_grid_norm - prodW_grid_norm) <= &
              1.0d-8 * max(1.0d0, abs(rawW_grid_norm))
            rel_residual_after_fragment = sqrt(max(0.0d0, norm_residual_after_fragment)) / &
              max(sqrt(max(0.0d0, norm_in)), 1.0d-300)
            rel_residual_after_bpw_perp = sqrt(max(0.0d0, norm_residual_after_bpw_perp)) / &
              max(sqrt(max(0.0d0, norm_in)), 1.0d-300)
            rel_residual_after_mixed_p = sqrt(max(0.0d0, norm_residual_after_mixed_p)) / &
              max(sqrt(max(0.0d0, norm_in)), 1.0d-300)
            rel_residual_after_neighbor_p = sqrt(max(0.0d0, norm_residual_after_neighbor_p)) / &
              max(sqrt(max(0.0d0, norm_in)), 1.0d-300)
            rel_residual_after_all_available = sqrt(max(0.0d0, norm_residual_after_all_available)) / &
              max(sqrt(max(0.0d0, norm_in)), 1.0d-300)
            prod_roundtrip_diff = sqrt(max(0.0d0, norm_residual_after_all_available))
            c_red_back(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), rhs_hybrid(:)))
            tmp_red(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), c_red_back(:))
            reduced_reproj_norm = real(sum(conjg(c_red_back(:)) * tmp_red(:)), kind=8)
            reduced_overlap_with_input = abs(tmp_red(jred))
            c_red_diff(:) = c_red_back(:)
            c_red_diff(jred) = c_red_diff(jred) - (1.0d0, 0.0d0)
            tmp_red(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), c_red_diff(:))
            red_reproj_diff = sqrt(max(0.0d0, real(sum(conjg(c_red_diff(:)) * tmp_red(:)), kind=8)))
            missing_reduced_norm = red_reproj_diff
            raw_back(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local), c_red_back(:))
            tmp_raw(:) = matmul(dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local), raw_back(:))
            tmp_red(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), tmp_raw(:)))
            tmp_red(:) = tmp_red(:) - c_red_back(:)
            tmp_red(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), tmp_red(:))
            basis_leakage = sqrt(max(0.0d0, real(sum(conjg(tmp_red(:)) * matmul(Sred_inv, tmp_red(:))), kind=8)))
            rhs_hybrid_normW(:) = rhs_hybrid(:)
            rhs_hybrid_normW(1:n_w) = rhs_normW(1:n_w)
            c_red_back_normW(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), &
              rhs_hybrid_normW(:)))
            tmp_red_normW(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
              c_red_back_normW(:))
            normW_reduced_reproj_norm = real(sum(conjg(c_red_back_normW(:)) * tmp_red_normW(:)), kind=8)
            normW_reduced_overlap = abs(tmp_red_normW(jred))
            c_red_diff_normW(:) = c_red_back_normW(:)
            c_red_diff_normW(jred) = c_red_diff_normW(jred) - (1.0d0, 0.0d0)
            tmp_red_normW(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
              c_red_diff_normW(:))
            normW_red_reproj_diff = sqrt(max(0.0d0, real(sum(conjg(c_red_diff_normW(:)) * &
              tmp_red_normW(:)), kind=8)))
            normW_missing_norm = normW_red_reproj_diff
            raw_back(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local), c_red_back_normW(:))
            tmp_raw(:) = matmul(dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local), raw_back(:))
            tmp_red_normW(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), tmp_raw(:)))
            tmp_red_normW(:) = tmp_red_normW(:) - c_red_back_normW(:)
            tmp_red_normW(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), tmp_red_normW(:))
            normW_basis_leakage = sqrt(max(0.0d0, real(sum(conjg(tmp_red_normW(:)) * &
              matmul(Sred_inv, tmp_red_normW(:))), kind=8)))
            c_red_back_canonicalW(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), &
              rhs_hybrid_canonicalW(:)))
            tmp_red_canonicalW(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
              c_red_back_canonicalW(:))
            canonicalW_reduced_reproj_norm = real(sum(conjg(c_red_back_canonicalW(:)) * tmp_red_canonicalW(:)), kind=8)
            canonicalW_reduced_overlap = abs(tmp_red_canonicalW(jred))
            c_red_diff_canonicalW(:) = c_red_back_canonicalW(:)
            c_red_diff_canonicalW(jred) = c_red_diff_canonicalW(jred) - (1.0d0, 0.0d0)
            tmp_red_canonicalW(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
              c_red_diff_canonicalW(:))
            canonicalW_red_reproj_diff = sqrt(max(0.0d0, real(sum(conjg(c_red_diff_canonicalW(:)) * &
              tmp_red_canonicalW(:)), kind=8)))
            canonicalW_missing_norm = canonicalW_red_reproj_diff
            raw_back(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local), c_red_back_canonicalW(:))
            tmp_raw(:) = matmul(dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local), raw_back(:))
            tmp_red_canonicalW(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), tmp_raw(:)))
            tmp_red_canonicalW(:) = tmp_red_canonicalW(:) - c_red_back_canonicalW(:)
            tmp_red_canonicalW(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
              tmp_red_canonicalW(:))
            canonicalW_basis_leakage = sqrt(max(0.0d0, real(sum(conjg(tmp_red_canonicalW(:)) * &
              matmul(Sred_inv, tmp_red_canonicalW(:))), kind=8)))
            rhs_hybrid_canonicalWP(:) = matmul(dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local), &
              raw_vec(1:nraw))
            c_red_back_canonicalWP(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), &
              rhs_hybrid_canonicalWP(:)))
            tmp_red_canonicalWP(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
              c_red_back_canonicalWP(:))
            canonicalWP_reduced_reproj_norm = real(sum(conjg(c_red_back_canonicalWP(:)) * &
              tmp_red_canonicalWP(:)), kind=8)
            canonicalWP_reduced_overlap = abs(tmp_red_canonicalWP(jred))
            c_red_diff_canonicalWP(:) = c_red_back_canonicalWP(:)
            c_red_diff_canonicalWP(jred) = c_red_diff_canonicalWP(jred) - (1.0d0, 0.0d0)
            tmp_red_canonicalWP(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
              c_red_diff_canonicalWP(:))
            canonicalWP_red_reproj_diff = sqrt(max(0.0d0, real(sum(conjg(c_red_diff_canonicalWP(:)) * &
              tmp_red_canonicalWP(:)), kind=8)))
            canonicalWP_missing_norm = canonicalWP_red_reproj_diff
            raw_back_canonicalWP(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local), &
              c_red_back_canonicalWP(:))
            raw_diff_canonicalWP(:) = raw_back_canonicalWP(:) - raw_vec(:)
            if (n_pw > 0 .and. n_w + n_pw <= nraw) then
              tmp_raw_canonicalWP(n_w+1:n_w+n_pw) = matmul( &
                dg_frag%wpw_reduced_Sraw_build(n_w+1:n_w+n_pw, n_w+1:n_w+n_pw, i_local), &
                raw_diff_canonicalWP(n_w+1:n_w+n_pw))
              rawPself_to_rawPself_diff = sqrt(max(0.0d0, real(sum( &
                conjg(raw_diff_canonicalWP(n_w+1:n_w+n_pw)) * tmp_raw_canonicalWP(n_w+1:n_w+n_pw)), kind=8)))
            end if
            if (n_w + n_pw < nraw) then
              tmp_raw_canonicalWP(n_w+n_pw+1:nraw) = matmul( &
                dg_frag%wpw_reduced_Sraw_build(n_w+n_pw+1:nraw, n_w+n_pw+1:nraw, i_local), &
                raw_diff_canonicalWP(n_w+n_pw+1:nraw))
              rawPneighbor_to_rawPneighbor_diff = sqrt(max(0.0d0, real(sum( &
                conjg(raw_diff_canonicalWP(n_w+n_pw+1:nraw)) * tmp_raw_canonicalWP(n_w+n_pw+1:nraw)), kind=8)))
            end if
            tmp_raw_canonicalWP(:) = matmul(dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local), &
              raw_back_canonicalWP(:))
            tmp_red_canonicalWP(:) = matmul(Sred_inv, &
              matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), &
              tmp_raw_canonicalWP(:)))
            tmp_red_canonicalWP(:) = tmp_red_canonicalWP(:) - c_red_back_canonicalWP(:)
            tmp_red_canonicalWP(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
              tmp_red_canonicalWP(:))
            canonicalWP_basis_leakage = sqrt(max(0.0d0, real(sum(conjg(tmp_red_canonicalWP(:)) * &
              matmul(Sred_inv, tmp_red_canonicalWP(:))), kind=8)))
            if (prodcoef_norm /= prodcoef_norm .or. local_norm /= local_norm .or. &
                back_local_norm /= back_local_norm .or. prod_roundtrip_diff /= prod_roundtrip_diff .or. &
                red_reproj_diff /= red_reproj_diff .or. basis_leakage /= basis_leakage .or. &
                norm_in /= norm_in .or. norm_after_fragment_projection /= norm_after_fragment_projection .or. &
                norm_residual_after_fragment /= norm_residual_after_fragment .or. &
                norm_after_bpw_perp_projection /= norm_after_bpw_perp_projection .or. &
                norm_residual_after_bpw_perp /= norm_residual_after_bpw_perp .or. &
                norm_after_mixed_p_projection /= norm_after_mixed_p_projection .or. &
                norm_residual_after_mixed_p /= norm_residual_after_mixed_p .or. &
                norm_after_neighbor_p_projection /= norm_after_neighbor_p_projection .or. &
                norm_residual_after_neighbor_p /= norm_residual_after_neighbor_p .or. &
                norm_after_all_available_projection /= norm_after_all_available_projection .or. &
                norm_residual_after_all_available /= norm_residual_after_all_available .or. &
                rel_residual_after_fragment /= rel_residual_after_fragment .or. &
                rel_residual_after_bpw_perp /= rel_residual_after_bpw_perp .or. &
                rel_residual_after_mixed_p /= rel_residual_after_mixed_p .or. &
                rel_residual_after_neighbor_p /= rel_residual_after_neighbor_p .or. &
                rel_residual_after_all_available /= rel_residual_after_all_available .or. &
                mixed_p_coef_norm /= mixed_p_coef_norm .or. neighbor_p_coef_norm /= neighbor_p_coef_norm .or. &
                reduced_reproj_norm /= reduced_reproj_norm .or. &
                reduced_overlap_with_input /= reduced_overlap_with_input .or. &
                missing_reduced_norm /= missing_reduced_norm .or. raw_W_norm /= raw_W_norm .or. &
                raw_P_self_norm /= raw_P_self_norm .or. raw_P_neighbor_norm /= raw_P_neighbor_norm .or. &
                raw_BPW_perp_norm /= raw_BPW_perp_norm .or. raw_mixed_P_norm /= raw_mixed_P_norm .or. &
                prod_W_norm /= prod_W_norm .or. rawW_grid_norm /= rawW_grid_norm .or. &
                prodW_grid_norm /= prodW_grid_norm .or. rawW_to_rawW_diff /= rawW_to_rawW_diff .or. &
                prodW_to_rawW_diff /= prodW_to_rawW_diff .or. &
                rawW_prodW_overlap /= rawW_prodW_overlap .or. &
                rawW_prodW_rel_diff /= rawW_prodW_rel_diff .or. &
                raw_W_norm_euclid /= raw_W_norm_euclid .or. &
                raw_W_norm_reduced_metric /= raw_W_norm_reduced_metric .or. &
                S_W_diag_min /= S_W_diag_min .or. S_W_diag_max /= S_W_diag_max .or. &
                gridS_W_diag_min /= gridS_W_diag_min .or. gridS_W_diag_max /= gridS_W_diag_max .or. &
                diag_ratio_min /= diag_ratio_min .or. diag_ratio_max /= diag_ratio_max .or. &
                plain_dot_diff /= plain_dot_diff .or. dV_dot_diff /= dV_dot_diff .or. &
                Sinv_project_diff /= Sinv_project_diff .or. &
                conjugate_basis_diff /= conjugate_basis_diff .or. &
                transpose_transform_diff /= transpose_transform_diff .or. &
                best_roundtrip_diff /= best_roundtrip_diff .or. best_norm_ratio /= best_norm_ratio .or. &
                normW_grid_norm /= normW_grid_norm .or. normW_to_rawW_diff /= normW_to_rawW_diff .or. &
                normW_red_reproj_diff /= normW_red_reproj_diff .or. &
                normW_reduced_reproj_norm /= normW_reduced_reproj_norm .or. &
                normW_reduced_overlap /= normW_reduced_overlap .or. normW_missing_norm /= normW_missing_norm .or. &
                normW_basis_leakage /= normW_basis_leakage .or. &
                canonicalW_grid_norm /= canonicalW_grid_norm .or. &
                canonicalW_to_rawW_diff /= canonicalW_to_rawW_diff .or. &
                canonicalW_red_reproj_diff /= canonicalW_red_reproj_diff .or. &
                canonicalW_reduced_reproj_norm /= canonicalW_reduced_reproj_norm .or. &
                canonicalW_reduced_overlap /= canonicalW_reduced_overlap .or. &
                canonicalW_missing_norm /= canonicalW_missing_norm .or. &
                canonicalW_basis_leakage /= canonicalW_basis_leakage .or. &
                canonicalWP_red_reproj_diff /= canonicalWP_red_reproj_diff .or. &
                canonicalWP_reduced_reproj_norm /= canonicalWP_reduced_reproj_norm .or. &
                canonicalWP_reduced_overlap /= canonicalWP_reduced_overlap .or. &
                canonicalWP_missing_norm /= canonicalWP_missing_norm .or. &
                canonicalWP_basis_leakage /= canonicalWP_basis_leakage .or. &
                rawPself_to_rawPself_diff /= rawPself_to_rawPself_diff .or. &
                rawPneighbor_to_rawPneighbor_diff /= rawPneighbor_to_rawPneighbor_diff .or. &
                Pself_grid_norm /= Pself_grid_norm .or. Pneighbor_grid_norm /= Pneighbor_grid_norm .or. &
                normW_grid_diag_min /= normW_grid_diag_min .or. normW_grid_diag_max /= normW_grid_diag_max) then
              bad = .true.
              reason = 'nan_or_inf'
            else
              reason = 'ok'
            end if
          end if
          if (.not. bad) then
            rel_prod_roundtrip = prod_roundtrip_diff / max(sqrt(max(0.0d0, local_norm)), 1.0d-300)
            rel_reduced_reproj = red_reproj_diff / max(sqrt(max(0.0d0, local_norm)), 1.0d-300)
          end if
          write(*,'(1x,a,6(a,i0),34(a,1pe12.4),a,a,a,l1,a,a)') &
            '[DG-WPW-RED-DIAG-EMBED-PRODCOEF]', &
            ' step=', step_use, ' frag=', ifrag, ' spin=', ispin, ' nred=', nred, ' j=', jred, &
            ' dominant_raw_component=', dominant_raw_component, &
            ' prodcoef_norm=', prodcoef_norm, &
            ' local_norm=', local_norm, &
            ' back_local_norm=', back_local_norm, &
            ' norm_in=', norm_in, &
            ' norm_after_fragment_projection=', norm_after_fragment_projection, &
            ' norm_residual_after_fragment=', norm_residual_after_fragment, &
            ' norm_after_bpw_perp_projection=', norm_after_bpw_perp_projection, &
            ' norm_residual_after_bpw_perp=', norm_residual_after_bpw_perp, &
            ' rel_residual_after_fragment=', rel_residual_after_fragment, &
            ' rel_residual_after_bpw_perp=', rel_residual_after_bpw_perp, &
            ' norm_after_mixed_p_projection=', norm_after_mixed_p_projection, &
            ' norm_residual_after_mixed_p=', norm_residual_after_mixed_p, &
            ' rel_residual_after_mixed_p=', rel_residual_after_mixed_p, &
            ' norm_after_neighbor_p_projection=', norm_after_neighbor_p_projection, &
            ' norm_residual_after_neighbor_p=', norm_residual_after_neighbor_p, &
            ' rel_residual_after_neighbor_p=', rel_residual_after_neighbor_p, &
            ' norm_after_all_available_projection=', norm_after_all_available_projection, &
            ' norm_residual_after_all_available=', norm_residual_after_all_available, &
            ' rel_residual_after_all_available=', rel_residual_after_all_available, &
            ' mixed_p_coef_norm=', mixed_p_coef_norm, &
            ' neighbor_p_coef_norm=', neighbor_p_coef_norm, &
            ' prod_roundtrip_diff_Snorm=', prod_roundtrip_diff, &
            ' rel_prod_roundtrip_diff=', rel_prod_roundtrip, &
            ' reduced_reproj_diff_Snorm=', red_reproj_diff, &
            ' rel_reduced_reproj_diff=', rel_reduced_reproj, &
            ' reduced_reproj_norm=', reduced_reproj_norm, &
            ' reduced_overlap_with_input=', reduced_overlap_with_input, &
            ' missing_reduced_norm=', missing_reduced_norm, &
            ' raw_W_norm=', raw_W_norm, &
            ' raw_P_self_norm=', raw_P_self_norm, &
            ' raw_P_neighbor_norm=', raw_P_neighbor_norm, &
            ' raw_BPW_perp_norm=', raw_BPW_perp_norm, &
            ' raw_mixed_P_norm=', raw_mixed_P_norm, &
            ' basis_leakage_Snorm=', basis_leakage, &
            ' route=', 'prod_fragment_plus_bpw_perp_plus_mixed_p', ' bad=', bad, ' reason=', trim(reason)
          write(*,'(1x,a,5(a,i0),15(a,1pe12.4),3(a,l1),a,a,a,l1,a,a)') &
            '[DG-WPW-RED-DIAG-EMBED-PRODCOEF]', &
            ' step=', step_use, ' frag=', ifrag, ' spin=', ispin, ' j=', jred, &
            ' nW_grid_keep=', nkeep_Wgrid, &
            ' raw_W_norm=', raw_W_norm, &
            ' gridW_norm_diag=', normW_grid_norm, &
            ' normW_grid_diag_min=', normW_grid_diag_min, &
            ' normW_grid_diag_max=', normW_grid_diag_max, &
            ' normW_grid_diag_avg=', normW_grid_diag_avg, &
            ' rawW_to_rawW_diff=', normW_to_rawW_diff, &
            ' reduced_reproj_diff_Snorm=', normW_red_reproj_diff, &
            ' reduced_reproj_norm=', normW_reduced_reproj_norm, &
            ' reduced_overlap_with_input=', normW_reduced_overlap, &
            ' missing_reduced_norm=', normW_missing_norm, &
            ' basis_leakage_Snorm=', normW_basis_leakage, &
            ' eval_min_Wgrid=', eval_min_Wgrid, &
            ' eval_max_Wgrid=', eval_max_Wgrid, &
            ' original_reduced_reproj_diff_Snorm=', red_reproj_diff, &
            ' original_missing_reduced_norm=', missing_reduced_norm, &
            ' normalization_match=', normW_to_rawW_diff <= 1.0d-8, &
            ' basis_order_match=', normW_to_rawW_diff <= 1.0d-8, &
            ' metric_match=', normW_basis_leakage <= 1.0d-8, &
            ' route=', 'prodcoef_with_reduced_normalized_W', ' bad=', bad, ' reason=', trim(reason)
          write(*,'(1x,a,4(a,i0),12(a,1pe12.4),3(a,l1),a,a,a,l1,a,a)') &
            '[DG-WPW-RED-DIAG-EMBED-PRODCOEF]', &
            ' step=', step_use, ' frag=', ifrag, ' spin=', ispin, ' j=', jred, &
            ' raw_W_norm=', raw_W_norm, &
            ' rawW_grid_norm=', canonicalW_grid_norm, &
            ' rawW_to_rawW_diff=', canonicalW_to_rawW_diff, &
            ' reduced_reproj_diff_Snorm=', canonicalW_red_reproj_diff, &
            ' reduced_reproj_norm=', canonicalW_reduced_reproj_norm, &
            ' reduced_overlap_with_input=', canonicalW_reduced_overlap, &
            ' missing_reduced_norm=', canonicalW_missing_norm, &
            ' basis_leakage_Snorm=', canonicalW_basis_leakage, &
            ' original_grid_rawW_grid_norm=', rawW_grid_norm, &
            ' original_grid_rawW_to_rawW_diff=', rawW_to_rawW_diff, &
            ' original_reduced_reproj_diff_Snorm=', red_reproj_diff, &
            ' original_missing_reduced_norm=', missing_reduced_norm, &
            ' normalization_match=', canonicalW_to_rawW_diff <= 1.0d-12, &
            ' basis_order_match=', canonicalW_to_rawW_diff <= 1.0d-12, &
            ' metric_match=', canonicalW_basis_leakage <= 1.0d-8, &
            ' route=', 'prodcoef_with_canonical_W', ' bad=', bad, ' reason=', trim(reason)
          write(*,'(1x,a,8(a,i0),16(a,1pe12.4),4(a,l1),a,a,a,l1,a,a)') &
            '[DG-WPW-RED-DIAG-EMBED-PRODCOEF]', &
            ' step=', step_use, ' frag=', ifrag, ' spin=', ispin, ' nred=', nred, ' j=', jred, &
            ' dominant_raw_component=', dominant_raw_component, &
            ' P_neighbor_owner_frag=', p_neighbor_owner_frag, &
            ' P_neighbor_source_frag=', p_neighbor_source_frag, &
            ' raw_W_norm=', raw_W_norm, &
            ' raw_P_self_norm=', raw_P_self_norm, &
            ' raw_P_neighbor_norm=', raw_P_neighbor_norm, &
            ' Pself_grid_norm=', Pself_grid_norm, &
            ' Pneighbor_grid_norm=', Pneighbor_grid_norm, &
            ' rawPself_to_rawPself_diff=', rawPself_to_rawPself_diff, &
            ' rawPneighbor_to_rawPneighbor_diff=', rawPneighbor_to_rawPneighbor_diff, &
            ' reduced_reproj_diff_Snorm=', canonicalWP_red_reproj_diff, &
            ' reduced_reproj_norm=', canonicalWP_reduced_reproj_norm, &
            ' reduced_overlap_with_input=', canonicalWP_reduced_overlap, &
            ' missing_reduced_norm=', canonicalWP_missing_norm, &
            ' basis_leakage_Snorm=', canonicalWP_basis_leakage, &
            ' original_canonicalW_reduced_reproj_diff_Snorm=', canonicalW_red_reproj_diff, &
            ' original_grid_reduced_reproj_diff_Snorm=', red_reproj_diff, &
            ' original_Pself_dominated_hint=', raw_P_self_norm, &
            ' original_Pneighbor_dominated_hint=', raw_P_neighbor_norm, &
            ' normalization_match=', canonicalWP_red_reproj_diff <= 1.0d-8, &
            ' basis_order_match=', canonicalWP_red_reproj_diff <= 1.0d-8, &
            ' metric_match=', canonicalWP_basis_leakage <= 1.0d-8, &
            ' neighbor_mapping_match=', neighbor_mapping_match, &
            ' route=', 'prodcoef_with_canonical_WP', ' bad=', bad, ' reason=', trim(reason)
          if (.not. bad .and. red_reproj_diff > 1.0d-1) then
            write(*,'(1x,a,6(a,i0),12(a,1pe12.4),a,a)') &
              '[DG-WPW-RED-DIAG-EMBED-PRODCOEF-DETAIL]', &
              ' step=', step_use, ' frag=', ifrag, ' spin=', ispin, ' nred=', nred, ' j=', jred, &
              ' dominant_raw_component=', dominant_raw_component, &
              ' input_norm=', local_norm, &
              ' recovered_norm=', reduced_reproj_norm, &
              ' overlap=', reduced_overlap_with_input, &
              ' missing_norm=', missing_reduced_norm, &
              ' raw_W_norm=', raw_W_norm, &
              ' raw_P_self_norm=', raw_P_self_norm, &
              ' raw_P_neighbor_norm=', raw_P_neighbor_norm, &
              ' raw_BPW_perp_norm=', raw_BPW_perp_norm, &
              ' raw_mixed_P_norm=', raw_mixed_P_norm, &
              ' norm_residual_after_bpw_perp=', norm_residual_after_bpw_perp, &
              ' norm_residual_after_all_available=', norm_residual_after_all_available, &
              ' basis_leakage_Snorm=', basis_leakage, &
              ' route=', 'prod_fragment_plus_bpw_perp_plus_mixed_p'
          end if
          if (.not. bad .and. (red_reproj_diff > 1.0d-1 .or. raw_W_norm > 9.0d-1)) then
            write(*,'(1x,a,4(a,i0),10(a,1pe12.4),4(a,l1),a,a,a,a)') &
              '[DG-WPW-RED-DIAG-WBLOCK-CONVENTION]', &
              ' step=', step_use, ' frag=', ifrag, ' spin=', ispin, ' j=', jred, &
              ' raw_W_norm=', raw_W_norm, &
              ' prod_W_norm=', prod_W_norm, &
              ' rawW_grid_norm=', rawW_grid_norm, &
              ' prodW_grid_norm=', prodW_grid_norm, &
              ' rawW_to_rawW_diff=', rawW_to_rawW_diff, &
              ' prodW_to_rawW_diff=', prodW_to_rawW_diff, &
              ' rawW_prodW_overlap=', rawW_prodW_overlap, &
              ' rawW_prodW_rel_diff=', rawW_prodW_rel_diff, &
              ' reduced_reproj_diff_Snorm=', red_reproj_diff, &
              ' missing_reduced_norm=', missing_reduced_norm, &
              ' basis_order_match=', basis_order_match, &
              ' metric_match=', metric_match, &
              ' normalization_match=', normalization_match, &
              ' bad=', bad, &
              ' route=', 'prod_W_block_convention_check', ' reason=', trim(reason)
            write(*,'(1x,a,4(a,i0),20(a,1pe12.4),a,a,a,l1,a,a)') &
              '[DG-WPW-RED-DIAG-WBLOCK-METRIC-AUDIT]', &
              ' step=', step_use, ' frag=', ifrag, ' spin=', ispin, ' j=', jred, &
              ' raw_W_norm_euclid=', raw_W_norm_euclid, &
              ' raw_W_norm_S=', raw_W_norm, &
              ' raw_W_norm_reduced_metric=', raw_W_norm_reduced_metric, &
              ' rawW_grid_norm=', rawW_grid_norm, &
              ' S_W_diag_min=', S_W_diag_min, &
              ' S_W_diag_max=', S_W_diag_max, &
              ' S_W_diag_avg=', S_W_diag_avg, &
              ' gridS_W_diag_min=', gridS_W_diag_min, &
              ' gridS_W_diag_max=', gridS_W_diag_max, &
              ' gridS_W_diag_avg=', gridS_W_diag_avg, &
              ' diag_ratio_min=', diag_ratio_min, &
              ' diag_ratio_max=', diag_ratio_max, &
              ' plain_dot_diff=', plain_dot_diff, &
              ' dV_dot_diff=', dV_dot_diff, &
              ' Sinv_project_diff=', Sinv_project_diff, &
              ' conjugate_basis_diff=', conjugate_basis_diff, &
              ' transpose_transform_diff=', transpose_transform_diff, &
              ' best_roundtrip_diff=', best_roundtrip_diff, &
              ' best_norm_ratio=', best_norm_ratio, &
              ' missing_reduced_norm=', missing_reduced_norm, &
              ' best_variant=', trim(best_variant), &
              ' bad=', bad, ' reason=', trim(reason)
          end if
        end do
      end do
      deallocate(raw_vec, raw_from_prod, raw_resid, raw_back, tmp_raw)
      deallocate(c_red_back, c_red_diff, tmp_red, Sred_inv)
      deallocate(basis_raw, wval, phi_val, pw_val, bpw_val, S_prod, S_prod_inv, rhs_prod, coef_prod)
      deallocate(S_density, S_density_inv, rhs_density, rhs_hybrid)
      deallocate(S_frag_pw, C_frag_pw, S_bpw, S_bpw_inv, rhs_bpw, coef_bpw)
      deallocate(S_mixed_p, S_mixed_p_inv, rhs_mixed_p, coef_mixed_p, tmp_mixed_p)
      deallocate(S_WW, S_WW_inv, raw_W_coef, prod_W_coef, rhs_rawW, rhs_prodW, back_W_coef, diff_W_coef, tmp_W)
      deallocate(rhs_plainW, rhs_conjW, rhs_transposeW)
      deallocate(S_W_grid, S_W_grid_invsqrt, Wnorm_metric)
      deallocate(wnorm_val, rhs_normW, back_normW_coef, diff_normW_coef)
      deallocate(rhs_hybrid_normW, c_red_back_normW, c_red_diff_normW, tmp_red_normW)
      deallocate(rhs_hybrid_canonicalW, c_red_back_canonicalW, c_red_diff_canonicalW, tmp_red_canonicalW)
      deallocate(rhs_hybrid_canonicalWP, c_red_back_canonicalWP, c_red_diff_canonicalWP, tmp_red_canonicalWP)
      deallocate(raw_back_canonicalWP, raw_diff_canonicalWP, tmp_raw_canonicalWP)
    end do
  end subroutine diagnose_wpw_reduced_embed_prodcoef

  subroutine diagnose_wpw_reduced_density(dg_frag, system, rho_prod, istep, nstep, coef_bad, dt, &
      rho_s_prod, replace_density, density_replaced, propagated_density_source, reproject_stage, &
      prodop_route, prodop_field_included, prodop_mixed_z_included, prodop_global_flux_included, &
      prodop_kick_applied, prodop_predictor_corrector_included, canonical_density_source, &
      canonical_pz_source, pz_can_raw, pz_prod_raw, pz_bad)
    use communication, only: comm_summation, comm_get_max
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(inout), optional :: rho_prod
    integer, intent(in), optional :: istep, nstep
    logical, intent(in), optional :: coef_bad
    real(8), intent(in), optional :: dt
    type(s_scalar), intent(inout), optional :: rho_s_prod(system%nspin)
    logical, intent(in), optional :: replace_density
    logical, intent(out), optional :: density_replaced
    logical, intent(in), optional :: propagated_density_source
    character(*), intent(in), optional :: reproject_stage
    character(*), intent(in), optional :: prodop_route
    logical, intent(in), optional :: prodop_field_included, prodop_mixed_z_included, prodop_global_flux_included
    logical, intent(in), optional :: prodop_kick_applied, prodop_predictor_corrector_included
    logical, intent(in), optional :: canonical_density_source
    logical, intent(in), optional :: canonical_pz_source
    real(8), intent(out), optional :: pz_can_raw, pz_prod_raw
    logical, intent(out), optional :: pz_bad

    integer :: i_local, ifrag, ispin, n_pw, n_w, nself, nred, nraw, nneigh, n_pfrag
    integer :: keep_n_env, keep_neighbor, nred_eff
    integer :: pfrag_ids(7), axis, side, jfrag, pidx, ipw, iw, jw, ib, nbf, ist, nocc_spin
    integer :: ix, iy, iz, gx, gy, gz, bx, by, bz, row0
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: nstate_store, nstate_metric, state_idx, i_state, top_iter, top_idx
    integer :: info_raw, info_build, raw_keep_count, nkeep_raw, nkeep_build, top_raw_state
    integer :: diag_step, diag_nstep, drift_interval
    integer :: canon_pack_count, canon_pack_bad, canon_recon_count, canon_recon_bad
    real(8) :: vol_weight, phase_arg, occ, chi, grad_chi(3), z_coord
    real(8) :: diag_time, dt_use
    real(8) :: q_self, q_neigh, q_cross, q_total
    real(8) :: q_reproject, rho_reproject_grid
    real(8) :: q_w, q_p, q_wp_cross
    real(8) :: rho_self_grid, rho_reduced_grid, delta_rho_grid
    real(8) :: rho_prod_grid, delta_prod_grid, rho_can_grid, canon_density_grid_diff
    real(8) :: rho_can_replace_grid
    real(8) :: local_max_delta, local_max_rho_self
    real(8) :: local_max_prod_delta, local_max_rho_prod
    real(8) :: local_max_coef_diff, local_max_hred_herm, local_max_sred_herm
    real(8) :: local_max_sred_eval, local_min_sred_eval, local_max_sred_cond
    real(8) :: local_max_sample_cond, local_max_sample_rank, local_min_sample_rank
    real(8) :: local_max_canon_density_diff
    real(8) :: max_in(20), max_out(20), rms_delta, delta_over_self, rms_prod_delta, prod_delta_over_prod
    real(8) :: hvol_weight, phys_scale
    real(8) :: occ_state, self_state_hvol, red_state_hvol, norm_state, loss_state
    real(8) :: sum_loss_self, sum_loss_red, max_loss_self, max_loss_red, avg_loss_self, avg_loss_red
    real(8) :: min_norm_red, max_norm_red, cond_norm_red, top_abs_loss
    real(8) :: frag_sum(5), local_sum(90), global_sum(90)
    real(8) :: trace_w, trace_self, trace_red
    real(8) :: loss_w, loss_self_case, loss_red_case, loss_raw_case
    real(8) :: trace_raw, max_loss_w_case, max_loss_self_case, max_loss_red_case, max_loss_raw_case
    real(8) :: sraw_min, sraw_max, sbuild_min, sbuild_max, cond_raw
    real(8) :: local_max_sdiff, local_max_sref, sdiff_entry_count, sdiff_rms, order_mismatch_score
    real(8) :: sbuild_trace, sdensity_trace, sbuild_norm2, sdensity_norm2, sdiff_norm2
    real(8) :: bsb_build_norm2, bsb_density_norm2, bsb_build_rms, bsb_density_rms
    real(8) :: bsb_hybrid_norm2, hybrid_norm2, hybrid_trace, hybrid_diff_norm2
    real(8) :: norm_current_sred, norm_reproject_sred, diff_sred
    real(8) :: sred_min_tmp, sred_max_tmp, hred_herm_tmp, sred_herm_tmp
    real(8) :: sample_smin, sample_smax, sample_resid_norm, sample_next_norm, sample_pred_norm
    real(8) :: sample_cond, sample_tol
    real(8) :: canon_norm, canon_roundtrip, canon_w_diff, canon_pself_diff, canon_pneighbor_diff
    real(8) :: canon_norm_sum, canon_roundtrip_max, canon_w_diff_max, canon_pself_diff_max
    real(8) :: canon_pneighbor_diff_max
    real(8) :: canon_recon_max_abs_diff, canon_recon_norm_diff_grid, canon_recon_prod_norm_grid
    real(8) :: canon_recon_can_norm_grid, canon_recon_density_int_diff, canon_recon_rms_accum
    real(8) :: phase2e_can_w_norm, phase2e_can_pself_norm, phase2e_can_pneighbor_norm
    real(8) :: phase2e_can_total_norm
    complex(8) :: sample_overlap
    complex(8) :: overlap_sred
    complex(8) :: phase_static
    integer :: top_w_state, top_self_state, top_red_state
    logical :: prod_grid_ok
    complex(8) :: phase, psi_w, psi_p, psi_self, psi_neigh, psi_total, psi_prod, psi_can, psi_diff
    complex(8), allocatable :: c_red(:), raw_coef(:,:), raw_coef_reproj(:,:), basis_raw(:), wval(:)
    complex(8), allocatable :: phi_val(:), S_prod(:,:), coef_prod(:), rhs_prod(:), rhs_canonical(:)
    complex(8), allocatable :: S_raw(:,:), rhs_raw(:,:), S_raw_inv(:,:), S_red_inv(:,:)
    complex(8), allocatable :: S_build(:,:), S_hybrid(:,:), S_build_inv(:,:)
    complex(8), allocatable :: bsb_build(:,:), bsb_density(:,:), bsb_hybrid(:,:)
    complex(8), allocatable :: c_raw(:), tmp_raw(:), tmp_raw_build(:)
    complex(8), allocatable :: c_can(:), c_can_back(:), c_can_diff(:), tmp_can(:), c_red_can(:), tmp_red_can(:)
    complex(8), allocatable :: c_can_store(:,:)
    complex(8), allocatable :: c_red_proj(:), c_red_build(:), c_red_hybrid(:), tmp_red(:), tmp_red_build(:), tmp_red_hybrid(:)
    complex(8), allocatable :: raw_back(:), raw_back_build(:), raw_back_hybrid(:), raw_resid(:), raw_resid_build(:), raw_resid_hybrid(:)
    complex(8), allocatable :: c_prev(:), c_pred(:), c_pred_resid(:), tmp_prev(:), tmp_pred(:), eig_amp(:)
    complex(8), allocatable :: c_reproj_curr(:,:), gram_sample(:,:), gram_inv_sample(:,:), sample_rhs(:), sample_weight(:)
    complex(8), allocatable :: sample_pred(:), sample_resid(:), sample_tmp(:)
    real(8), allocatable :: state_local(:,:), state_global(:,:), state_loss_work(:)
    logical :: found, bad_density, bad_reproject
    logical :: do_heavy_diag, do_drift_log, do_reproject, do_projection_diag, drift_enabled, density_enabled
    logical :: use_canonical_reproject, canon_pack_diag, canon_recon_diag, canon_obs_diag, canon_hook_dryrun
    logical :: canon_use_density, canon_use_pz
    logical :: canon_grid_diag
    logical :: propagated_debug, state_prop_diag, sample_u_diag, sample_u_requested
    logical :: do_pz_series, replace_density_active, replace_density_zeroed
    logical :: canon_replace_density_active, canon_replace_density_zeroed
    logical :: replace_density_from_propagated
    logical :: state_prop_before, state_prop_after, state_prop_can_static, prev_realloc
    logical :: prodop_field_flag, prodop_mixed_z_flag, prodop_global_flux_flag
    logical :: prodop_kick_flag, prodop_predictor_corrector_flag
    logical :: phase2e_source_from_prodcoef, phase2e_source_from_reproject_coef
    logical :: phase2e_source_from_density_reconstruction
    character(len=32) :: prodop_route_label
    logical, save :: heavy_diag_done = .false.

    if (present(density_replaced)) density_replaced = .false.
    if (present(pz_can_raw)) pz_can_raw = 0.0d0
    if (present(pz_prod_raw)) pz_prod_raw = 0.0d0
    if (present(pz_bad)) pz_bad = .true.
    density_enabled = wpw_reduced_density_diag_enabled()
    use_canonical_reproject = wpw_reduced_canonical_reproject_enabled()
    canon_pack_diag = wpw_reduced_canon_pack_diag_enabled()
    canon_recon_diag = wpw_reduced_canon_recon_diag_enabled()
    canon_obs_diag = wpw_reduced_canon_obs_diag_enabled()
    canon_hook_dryrun = wpw_reduced_canon_hook_dryrun_enabled()
    canon_use_density = .false.
    if (present(canonical_density_source)) canon_use_density = canonical_density_source
    canon_use_pz = .false.
    if (present(canonical_pz_source)) canon_use_pz = canonical_pz_source
    canon_replace_density_active = canon_use_density .and. present(rho_prod) .and. present(rho_s_prod)
    canon_replace_density_zeroed = .false.
    canon_grid_diag = canon_recon_diag .or. canon_obs_diag .or. canon_hook_dryrun .or. canon_use_pz .or. &
      canon_replace_density_active
    drift_enabled = wpw_reduced_drift_diag_enabled()
    propagated_debug = wpw_reduced_propagated_debug_enabled()
    state_prop_diag = wpw_reduced_state_prop_diag_enabled()
    sample_u_requested = wpw_reduced_sample_u_diag_enabled()
    sample_u_diag = sample_u_requested
    sample_tol = wpw_reduced_sample_u_tol()
    keep_n_env = wpw_reduced_keep_n()
    do_pz_series = wpw_reduced_pz_series_enabled()
    replace_density_active = .false.
    if (present(replace_density)) replace_density_active = replace_density
    replace_density_active = replace_density_active .and. present(rho_prod) .and. present(rho_s_prod)
    replace_density_zeroed = .false.
    replace_density_from_propagated = .false.
    if (present(propagated_density_source)) replace_density_from_propagated = propagated_density_source
    if (.not. density_enabled .and. .not. drift_enabled .and. .not. replace_density_active .and. &
        .not. canon_replace_density_active .and. .not. canon_use_pz .and. &
        .not. canon_pack_diag .and. .not. canon_grid_diag .and. &
        .not. state_prop_diag .and. .not. sample_u_requested) return
    if (.not. dg_frag%wpw_reduced_ready) return
    if (.not. dg_frag%wpw_reduced_coef_initialized) return
    if (.not. allocated(dg_frag%wpw_reduced_transform)) return
    if (.not. allocated(dg_frag%wpw_reduced_nraw)) return
    if (.not. allocated(dg_frag%coef_wpw_self)) return
    if (.not. allocated(dg_frag%coef_wpw_neighbor_reduced)) return
    if (.not. allocated(dg_frag%global_wannier_local_coef)) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return
    if (.not. allocated(dg_frag%phi_frag) .and. .not. allocated(dg_frag%phi_frag_c)) return

    diag_step = -1
    diag_nstep = -1
    if (present(istep)) diag_step = istep
    if (present(nstep)) diag_nstep = nstep
    diag_time = 0.0d0
    dt_use = 0.0d0
    if (present(dt)) dt_use = dt
    if (diag_step >= 0) diag_time = dble(diag_step) * dt_use
    drift_interval = wpw_reduced_drift_interval()
    do_heavy_diag = density_enabled .and. propagated_debug .and. .not. heavy_diag_done
    do_drift_log = drift_enabled
    if (do_drift_log .and. diag_step > 0) then
      do_drift_log = (diag_step == 1 .or. diag_step == 2 .or. diag_step == 10 .or. diag_step == 100 .or. &
        (diag_nstep > 0 .and. diag_step == diag_nstep) .or. mod(diag_step, drift_interval) == 0)
    end if
    state_prop_before = .false.
    state_prop_after = .false.
    if (present(reproject_stage)) then
      state_prop_before = trim(reproject_stage) == 'before-production'
      state_prop_after = trim(reproject_stage) == 'after-production'
    end if
    prodop_route_label = 'unknown'
    if (present(prodop_route)) prodop_route_label = trim(prodop_route)
    prodop_field_flag = .false.
    prodop_mixed_z_flag = .false.
    prodop_global_flux_flag = .false.
    prodop_kick_flag = .false.
    prodop_predictor_corrector_flag = .false.
    if (present(prodop_field_included)) prodop_field_flag = prodop_field_included
    if (present(prodop_mixed_z_included)) prodop_mixed_z_flag = prodop_mixed_z_included
    if (present(prodop_global_flux_included)) prodop_global_flux_flag = prodop_global_flux_included
    if (present(prodop_kick_applied)) prodop_kick_flag = prodop_kick_applied
    if (present(prodop_predictor_corrector_included)) &
      prodop_predictor_corrector_flag = prodop_predictor_corrector_included
    state_prop_diag = (state_prop_diag .or. sample_u_requested) .and. (state_prop_before .or. state_prop_after)
    sample_u_diag = sample_u_requested .and. state_prop_after
    do_reproject = do_drift_log .or. do_heavy_diag .or. do_pz_series .or. &
      (replace_density_active .and. .not. replace_density_from_propagated) .or. state_prop_diag .or. sample_u_diag
    do_projection_diag = do_heavy_diag .or. do_reproject .or. canon_pack_diag .or. canon_grid_diag
    if (.not. do_heavy_diag .and. .not. do_drift_log .and. .not. do_pz_series .and. &
        .not. replace_density_active .and. .not. canon_replace_density_active .and. .not. canon_use_pz .and. &
        .not. state_prop_diag .and. .not. sample_u_diag .and. &
        .not. canon_pack_diag .and. .not. canon_grid_diag) return

    n_pw = dg_frag%n_plane_waves
    if (n_pw <= 0) return
    nstate_store = max(1, size(dg_frag%coef_wpw_self, 2))
    nstate_metric = nstate_store * max(1, dg_frag%nspin)
    if (state_prop_diag) then
      prev_realloc = .false.
      if (.not. allocated(dg_frag%wpw_reproject_prev_coef)) then
        prev_realloc = .true.
      else if (size(dg_frag%wpw_reproject_prev_coef, 1) /= dg_frag%wpw_reduced_max_dim .or. &
               size(dg_frag%wpw_reproject_prev_coef, 2) /= nstate_store .or. &
               size(dg_frag%wpw_reproject_prev_coef, 3) /= dg_frag%nspin .or. &
               size(dg_frag%wpw_reproject_prev_coef, 4) /= size(dg_frag%wpw_reduced_dim)) then
        deallocate(dg_frag%wpw_reproject_prev_coef)
        prev_realloc = .true.
      end if
      if (prev_realloc) then
        allocate(dg_frag%wpw_reproject_prev_coef(dg_frag%wpw_reduced_max_dim, nstate_store, &
          dg_frag%nspin, size(dg_frag%wpw_reduced_dim)))
        dg_frag%wpw_reproject_prev_coef(:, :, :, :) = (0.0d0, 0.0d0)
        dg_frag%wpw_reproject_prev_valid = .false.
      end if
    end if
    allocate(state_local(9, nstate_metric), state_global(9, nstate_metric), state_loss_work(nstate_metric))
    state_local(:, :) = 0.0d0
    state_global(:, :) = 0.0d0
    state_loss_work(:) = -1.0d300
    vol_weight = product(dg_frag%hgs(1:3)) / &
      product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    hvol_weight = product(dg_frag%hgs(1:3))
    if (allocated(dg_frag%phi_frag_c)) then
      p_lb1 = lbound(dg_frag%phi_frag_c, 1); p_ub1 = ubound(dg_frag%phi_frag_c, 1)
      p_lb2 = lbound(dg_frag%phi_frag_c, 2); p_ub2 = ubound(dg_frag%phi_frag_c, 2)
      p_lb3 = lbound(dg_frag%phi_frag_c, 3); p_ub3 = ubound(dg_frag%phi_frag_c, 3)
    else
      p_lb1 = lbound(dg_frag%phi_frag, 1); p_ub1 = ubound(dg_frag%phi_frag, 1)
      p_lb2 = lbound(dg_frag%phi_frag, 2); p_ub2 = ubound(dg_frag%phi_frag, 2)
      p_lb3 = lbound(dg_frag%phi_frag, 3); p_ub3 = ubound(dg_frag%phi_frag, 3)
    end if

    local_sum(:) = 0.0d0
    local_max_delta = 0.0d0
    local_max_rho_self = 0.0d0
    local_max_prod_delta = 0.0d0
    local_max_rho_prod = 0.0d0
    local_max_coef_diff = 0.0d0
    local_max_hred_herm = 0.0d0
    local_max_sred_herm = 0.0d0
    local_max_sred_eval = 0.0d0
    local_min_sred_eval = huge(1.0d0)
    local_max_sred_cond = 0.0d0
    local_max_sample_cond = 0.0d0
    local_max_sample_rank = 0.0d0
    local_min_sample_rank = huge(1.0d0)
    local_max_canon_density_diff = 0.0d0
    local_max_sdiff = 0.0d0
    local_max_sref = 0.0d0
    raw_keep_count = 0
    cond_raw = 0.0d0
    bad_density = .false.
    bad_reproject = .false.
    if (present(coef_bad)) then
      if (coef_bad) local_sum(53) = 1.0d0
    end if
    do i_local = 1, size(dg_frag%wpw_reduced_dim)
      ifrag = dg_frag%ifrag_start + i_local - 1
      nred = dg_frag%wpw_reduced_dim(i_local)
      nself = dg_frag%wpw_reduced_nself(i_local)
      nraw = dg_frag%wpw_reduced_nraw(i_local)
      n_w = dg_frag%global_wannier_local_nkeep(i_local)
      if (nred <= 0 .or. nself <= 0 .or. nraw <= 0 .or. n_w <= 0) cycle
      if (nself /= n_w + n_pw) cycle
      nneigh = nred - nself
      if (keep_n_env >= 0) then
        keep_neighbor = min(max(0, keep_n_env), max(0, nneigh))
      else
        keep_neighbor = max(0, nneigh)
      end if
      nred_eff = nself + keep_neighbor

      n_pfrag = 1
      pfrag_ids(1) = ifrag
      do axis = 1, 3
        do side = -1, 1, 2
          jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle
          found = any(pfrag_ids(1:n_pfrag) == jfrag)
          if (.not. found .and. n_pfrag < size(pfrag_ids)) then
            n_pfrag = n_pfrag + 1
            pfrag_ids(n_pfrag) = jfrag
          end if
        end do
      end do
      if (nraw /= n_w + n_pfrag * n_pw) cycle
      nbf = min(dg_frag%n_basis(ifrag, 1), size(dg_frag%global_wannier_local_coef, 1))
      if (allocated(dg_frag%phi_frag_c)) then
        nbf = min(nbf, size(dg_frag%phi_frag_c, 4))
      else
        nbf = min(nbf, size(dg_frag%phi_frag, 4))
      end if

      allocate(c_red(nred), raw_coef(nraw, max(1, size(dg_frag%coef_wpw_self, 2))))
      if (do_reproject) allocate(raw_coef_reproj(nraw, max(1, size(dg_frag%coef_wpw_self, 2))))
      allocate(basis_raw(nraw), wval(n_w))
      if (do_projection_diag .and. (use_canonical_reproject .or. canon_pack_diag .or. canon_grid_diag)) then
        allocate(phi_val(nbf), S_prod(nbf,nbf), coef_prod(nbf), rhs_prod(nbf), rhs_canonical(nraw))
        S_prod(:, :) = (0.0d0, 0.0d0)
      end if
      local_sum(18) = local_sum(18) + dble(n_w)
      local_sum(19) = local_sum(19) + dble(nself)
      local_sum(20) = local_sum(20) + dble(nred)
      local_sum(59) = local_sum(59) + dble(nred_eff)
      frag_sum(:) = 0.0d0
      do ispin = 1, dg_frag%nspin
        nocc_spin = 0
        if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
        nocc_spin = min(nocc_spin, size(dg_frag%coef_wpw_self, 2))
        if (allocated(system%rocc)) nocc_spin = min(nocc_spin, size(system%rocc, 1))
        if (nocc_spin <= 0) cycle
        if (do_projection_diag) then
          allocate(S_raw(nraw,nraw), rhs_raw(nraw,nocc_spin), S_raw_inv(nraw,nraw), S_red_inv(nred,nred))
          allocate(S_build(nraw,nraw), S_hybrid(nraw,nraw), S_build_inv(nraw,nraw))
          allocate(bsb_build(nred,nred), bsb_density(nred,nred), bsb_hybrid(nred,nred))
          allocate(c_raw(nraw), tmp_raw(nraw), tmp_raw_build(nraw))
          if (canon_pack_diag .or. canon_grid_diag) allocate(c_can(nraw), c_can_back(nraw), c_can_diff(nraw), tmp_can(nraw), &
            c_red_can(nred), tmp_red_can(nred))
          if (canon_grid_diag) allocate(c_can_store(nraw,nocc_spin))
          allocate(c_red_proj(nred), c_red_build(nred), c_red_hybrid(nred))
          allocate(tmp_red(nred), tmp_red_build(nred), tmp_red_hybrid(nred))
          allocate(raw_back(nraw), raw_back_build(nraw), raw_back_hybrid(nraw))
          allocate(raw_resid(nraw), raw_resid_build(nraw), raw_resid_hybrid(nraw))
          if (sample_u_diag) then
            allocate(c_reproj_curr(nred,nocc_spin), gram_sample(nocc_spin,nocc_spin), &
              gram_inv_sample(nocc_spin,nocc_spin))
            allocate(sample_rhs(nocc_spin), sample_weight(nocc_spin))
            allocate(sample_pred(nred), sample_resid(nred), sample_tmp(nred))
            c_reproj_curr(:, :) = (0.0d0, 0.0d0)
          end if
          if (state_prop_diag) then
            allocate(c_prev(nred), c_pred(nred), c_pred_resid(nred), tmp_prev(nred), tmp_pred(nred), eig_amp(nred))
          end if
          S_raw(:, :) = (0.0d0, 0.0d0)
          S_build(:, :) = (0.0d0, 0.0d0)
          S_hybrid(:, :) = (0.0d0, 0.0d0)
          if ((use_canonical_reproject .or. canon_pack_diag .or. canon_grid_diag) .and. allocated(S_prod)) &
            S_prod(:, :) = (0.0d0, 0.0d0)
          if (allocated(c_can_store)) c_can_store(:, :) = (0.0d0, 0.0d0)
          if (allocated(dg_frag%wpw_reduced_Sraw_build)) then
            S_build(:, :) = dg_frag%wpw_reduced_Sraw_build(1:nraw, 1:nraw, i_local)
          end if
          rhs_raw(:, :) = (0.0d0, 0.0d0)
        end if
        raw_coef(:, :) = (0.0d0, 0.0d0)
        if (do_reproject) raw_coef_reproj(:, :) = (0.0d0, 0.0d0)
        do ist = 1, nocc_spin
          c_red(:) = (0.0d0, 0.0d0)
          c_red(1:nself) = dg_frag%coef_wpw_self(1:nself, ist, ispin, i_local)
          if (nneigh > 0) c_red(nself+1:nred) = &
            dg_frag%coef_wpw_neighbor_reduced(1:nneigh, ist, ispin, i_local)
          raw_coef(1:nraw, ist) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local), c_red(1:nred))
        end do

        do iz = 1, dg_frag%nxyz_domain(3, ifrag)
          gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
          bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
          if (bz < p_lb3 .or. bz > p_ub3) cycle
          do iy = 1, dg_frag%nxyz_domain(2, ifrag)
            gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
            by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
            if (by < p_lb2 .or. by > p_ub2) cycle
            do ix = 1, dg_frag%nxyz_domain(1, ifrag)
              gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
              bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
              if (bx < p_lb1 .or. bx > p_ub1) cycle
              basis_raw(:) = (0.0d0, 0.0d0)
              wval(:) = (0.0d0, 0.0d0)
              do iw = 1, n_w
                do ib = 1, nbf
                  if (allocated(dg_frag%phi_frag_c)) then
                    wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                      dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                  else
                    wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                      cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                  end if
                end do
              end do
              basis_raw(1:n_w) = wval(1:n_w)
              do pidx = 1, n_pfrag
                row0 = n_w + (pidx - 1) * n_pw
                call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                do ipw = 1, n_pw
                  phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                              dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                              dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                  phase = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                  basis_raw(row0 + ipw) = phase
                end do
              end do
              if (do_projection_diag) then
                do iw = 1, nraw
                  do ipw = 1, nraw
                    S_raw(iw, ipw) = S_raw(iw, ipw) + conjg(basis_raw(iw)) * basis_raw(ipw) * vol_weight
                  end do
                end do
                if (use_canonical_reproject .or. canon_pack_diag .or. canon_grid_diag) then
                  phi_val(:) = (0.0d0, 0.0d0)
                  do ib = 1, nbf
                    if (allocated(dg_frag%phi_frag_c)) then
                      phi_val(ib) = dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                    else
                      phi_val(ib) = cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                    end if
                  end do
                  do iw = 1, nbf
                    do ipw = 1, nbf
                      S_prod(iw, ipw) = S_prod(iw, ipw) + conjg(phi_val(iw)) * phi_val(ipw) * vol_weight
                    end do
                  end do
                end if
              end if
              do ist = 1, nocc_spin
                if (ist == 1) then
                  rho_self_grid = 0.0d0
                  rho_reduced_grid = 0.0d0
                end if
                occ = 1.0d0
                if (allocated(system%rocc)) occ = max(0.0d0, system%rocc(ist, 1, ispin))
                if (occ <= 0.0d0) cycle
                state_idx = (ispin - 1) * nstate_store + ist
                if (ix == 1 .and. iy == 1 .and. iz == 1) then
                  local_sum(15) = local_sum(15) + occ
                  local_sum(16) = local_sum(16) + 1.0d0
                  state_local(1, state_idx) = state_local(1, state_idx) + occ
                end if
                psi_w = sum(raw_coef(1:n_w, ist) * basis_raw(1:n_w))
                if (n_pw > 0) then
                  psi_p = sum(raw_coef(n_w+1:nself, ist) * basis_raw(n_w+1:nself))
                else
                  psi_p = (0.0d0, 0.0d0)
                end if
                psi_self = psi_w + psi_p
                if (nraw > nself) then
                  psi_neigh = sum(raw_coef(nself+1:nraw, ist) * basis_raw(nself+1:nraw))
                else
                  psi_neigh = (0.0d0, 0.0d0)
                end if
                psi_total = psi_self + psi_neigh
                psi_prod = (0.0d0, 0.0d0)
                do ib = 1, nbf
                  if (dg_frag%index_basis(ib, ifrag, ispin) <= 0) cycle
                  if (dg_frag%index_basis(ib, ifrag, ispin) > dg_frag%n_mat_max) cycle
                  if (dg_frag%coef_global_to_local(dg_frag%index_basis(ib, ifrag, ispin), ispin) <= 0) cycle
                  if (dg_frag%coef_global_to_local(dg_frag%index_basis(ib, ifrag, ispin), ispin) > size(dg_frag%coef, 1)) cycle
                  if (allocated(dg_frag%phi_frag_c)) then
                    psi_prod = psi_prod + dg_frag%coef(dg_frag%coef_global_to_local( &
                      dg_frag%index_basis(ib, ifrag, ispin), ispin), ist, ispin) * &
                      dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                  else
                    psi_prod = psi_prod + dg_frag%coef(dg_frag%coef_global_to_local( &
                      dg_frag%index_basis(ib, ifrag, ispin), ispin), ist, ispin) * &
                      cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                  end if
                end do
                if (do_projection_diag) &
                  rhs_raw(1:nraw, ist) = rhs_raw(1:nraw, ist) + conjg(basis_raw(1:nraw)) * psi_prod * vol_weight
                q_self = occ * abs(psi_self)**2 * vol_weight
                q_neigh = occ * abs(psi_neigh)**2 * vol_weight
                q_total = occ * abs(psi_total)**2 * vol_weight
                q_cross = q_total - q_self - q_neigh
                q_w = occ * abs(psi_w)**2 * vol_weight
                q_p = occ * abs(psi_p)**2 * vol_weight
                q_wp_cross = q_self - q_w - q_p
                rho_self_grid = rho_self_grid + occ * abs(psi_self)**2
                rho_reduced_grid = rho_reduced_grid + occ * abs(psi_total)**2
                frag_sum(1) = frag_sum(1) + q_total
                frag_sum(2) = frag_sum(2) + q_self
                frag_sum(3) = frag_sum(3) + q_neigh
                frag_sum(4) = frag_sum(4) + q_cross
                frag_sum(5) = frag_sum(5) + occ * vol_weight
                local_sum(1) = local_sum(1) + q_total
                local_sum(2) = local_sum(2) + q_self
                local_sum(3) = local_sum(3) + q_neigh
                local_sum(4) = local_sum(4) + q_cross
                local_sum(5) = local_sum(5) + occ * vol_weight
                state_local(2, state_idx) = state_local(2, state_idx) + q_self
                state_local(3, state_idx) = state_local(3, state_idx) + q_total
                state_local(4, state_idx) = state_local(4, state_idx) + q_neigh
                state_local(5, state_idx) = state_local(5, state_idx) + q_cross
                state_local(6, state_idx) = state_local(6, state_idx) + q_w
                state_local(7, state_idx) = state_local(7, state_idx) + q_p
                state_local(8, state_idx) = state_local(8, state_idx) + q_wp_cross
                if (q_total /= q_total .or. q_self /= q_self .or. q_neigh /= q_neigh) bad_density = .true.
              end do
              if (replace_density_active .and. replace_density_from_propagated) then
                if (.not. replace_density_zeroed) then
                  rho_prod%f(:, :, :) = 0.0d0
                  do iw = 1, system%nspin
                    rho_s_prod(iw)%f(:, :, :) = 0.0d0
                  end do
                  replace_density_zeroed = .true.
                end if
                rho_s_prod(ispin)%f(gx, gy, gz) = rho_s_prod(ispin)%f(gx, gy, gz) + rho_reduced_grid
                rho_prod%f(gx, gy, gz) = rho_prod%f(gx, gy, gz) + rho_reduced_grid
              end if
              delta_rho_grid = rho_reduced_grid - rho_self_grid
              local_sum(7) = local_sum(7) + delta_rho_grid * vol_weight
              local_sum(8) = local_sum(8) + delta_rho_grid**2 * vol_weight
              local_max_delta = max(local_max_delta, abs(delta_rho_grid))
              local_max_rho_self = max(local_max_rho_self, abs(rho_self_grid))
              if (present(rho_prod)) then
                rho_prod_grid = rho_prod%f(gx, gy, gz)
                prod_grid_ok = (rho_prod_grid == rho_prod_grid) .and. abs(rho_prod_grid) < huge(1.0d0)
                if (prod_grid_ok) then
                  if (mod(dg_frag%lgnum_total(3), 2) == 0) then
                    z_coord = dble(gz) - 0.5d0
                  else
                    z_coord = dble(gz)
                  end if
                  delta_prod_grid = rho_reduced_grid - rho_prod_grid
                  local_sum(9) = local_sum(9) + rho_prod_grid * vol_weight
                  local_sum(10) = local_sum(10) + delta_prod_grid * vol_weight
                  local_sum(11) = local_sum(11) + delta_prod_grid**2 * vol_weight
                  local_sum(12) = local_sum(12) + vol_weight
                  local_sum(60) = local_sum(60) + rho_prod_grid * z_coord * vol_weight
                  local_max_prod_delta = max(local_max_prod_delta, abs(delta_prod_grid))
                  local_max_rho_prod = max(local_max_rho_prod, abs(rho_prod_grid))
                  if (delta_prod_grid /= delta_prod_grid) bad_density = .true.
                else
                  local_sum(13) = local_sum(13) + vol_weight
                  bad_density = .true.
                end if
              end if
              if (ispin == 1) local_sum(14) = local_sum(14) + vol_weight
              if (ispin == 1) local_sum(17) = local_sum(17) + hvol_weight
              if (rho_reduced_grid /= rho_reduced_grid .or. delta_rho_grid /= delta_rho_grid) bad_density = .true.
            end do
          end do
        end do
        if (do_projection_diag) then
          call hermitize_matrix(S_raw, nraw)
          if (use_canonical_reproject .or. canon_pack_diag .or. canon_grid_diag) call hermitize_matrix(S_prod, nbf)
          if (allocated(dg_frag%wpw_reduced_Sraw_build)) then
          call hermitize_matrix(S_build, nraw)
          S_hybrid(:, :) = S_build(:, :)
          call hermitize_matrix(S_hybrid, nraw)
          call build_hermitian_pseudoinverse(S_build, nraw, 1.0d-8, S_build_inv, &
            info_build, sbuild_min, sbuild_max, nkeep_build)
          sdiff_norm2 = sum(abs(S_raw(1:nraw,1:nraw) - S_build(1:nraw,1:nraw))**2)
          sbuild_norm2 = sum(abs(S_build(1:nraw,1:nraw))**2)
          sdensity_norm2 = sum(abs(S_raw(1:nraw,1:nraw))**2)
          sbuild_trace = real(sum([(S_build(iw,iw), iw=1,nraw)]), kind=8)
          sdensity_trace = real(sum([(S_raw(iw,iw), iw=1,nraw)]), kind=8)
          if (do_heavy_diag) then
            local_sum(26) = local_sum(26) + 1.0d0
            local_sum(27) = local_sum(27) + sdiff_norm2
            local_sum(28) = local_sum(28) + sbuild_trace
            local_sum(29) = local_sum(29) + sdensity_trace
            local_sum(30) = local_sum(30) + sbuild_norm2
            local_sum(31) = local_sum(31) + sdensity_norm2
            local_sum(32) = local_sum(32) + dble(max(0, nkeep_build))
            local_sum(33) = local_sum(33) + dble(max(0, nraw))
            local_sum(39) = local_sum(39) + dble(nraw) * dble(nraw)
            local_max_sdiff = max(local_max_sdiff, maxval(abs(S_raw(1:nraw,1:nraw) - S_build(1:nraw,1:nraw))))
            local_max_sref = max(local_max_sref, maxval(abs(S_raw(1:nraw,1:nraw))), &
              maxval(abs(S_build(1:nraw,1:nraw))))
          end if

          if (do_heavy_diag) then
            bsb_build(:, :) = matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), &
              matmul(S_build(1:nraw,1:nraw), dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local)))
            bsb_density(:, :) = matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), &
              matmul(S_raw(1:nraw,1:nraw), dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local)))
            bsb_build_norm2 = sum(abs(bsb_build(1:nred,1:nred) - &
              dg_frag%wpw_reduced_S(1:nred,1:nred,ispin,i_local))**2)
            bsb_density_norm2 = sum(abs(bsb_density(1:nred,1:nred) - &
              dg_frag%wpw_reduced_S(1:nred,1:nred,ispin,i_local))**2)
            bsb_hybrid(:, :) = matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local))), &
              matmul(S_hybrid(1:nraw,1:nraw), dg_frag%wpw_reduced_transform(1:nraw, 1:nred, i_local)))
            bsb_hybrid_norm2 = sum(abs(bsb_hybrid(1:nred,1:nred) - &
              dg_frag%wpw_reduced_S(1:nred,1:nred,ispin,i_local))**2)
            hybrid_diff_norm2 = sum(abs(S_hybrid(1:nraw,1:nraw) - S_build(1:nraw,1:nraw))**2)
            hybrid_norm2 = sum(abs(S_hybrid(1:nraw,1:nraw))**2)
            hybrid_trace = real(sum([(S_hybrid(iw,iw), iw=1,nraw)]), kind=8)
            local_sum(37) = local_sum(37) + bsb_build_norm2
            local_sum(38) = local_sum(38) + bsb_density_norm2
            local_sum(46) = local_sum(46) + hybrid_trace
            local_sum(47) = local_sum(47) + hybrid_diff_norm2
            local_sum(48) = local_sum(48) + hybrid_norm2
            local_sum(49) = local_sum(49) + bsb_hybrid_norm2
            do iw = 1, n_w
              local_sum(40) = local_sum(40) + real(S_build(iw, iw), kind=8)
              local_sum(43) = local_sum(43) + real(S_raw(iw, iw), kind=8)
            end do
            do iw = n_w + 1, nself
              local_sum(41) = local_sum(41) + real(S_build(iw, iw), kind=8)
              local_sum(44) = local_sum(44) + real(S_raw(iw, iw), kind=8)
            end do
            do iw = nself + 1, nraw
              local_sum(42) = local_sum(42) + real(S_build(iw, iw), kind=8)
              local_sum(45) = local_sum(45) + real(S_raw(iw, iw), kind=8)
            end do
          end if
        else
          info_build = -1
          nkeep_build = 0
          sbuild_min = 0.0d0
          sbuild_max = 0.0d0
          S_build_inv(:, :) = (0.0d0, 0.0d0)
        end if
        call build_hermitian_pseudoinverse(S_raw, nraw, 1.0d-8, S_raw_inv, info_raw, sraw_min, sraw_max, nkeep_raw)
        if (info_raw == 0) then
          call build_hermitian_inverse(dg_frag%wpw_reduced_S(1:nred_eff, 1:nred_eff, ispin, i_local), &
            nred_eff, S_red_inv(1:nred_eff, 1:nred_eff), info_raw, sraw_min, sraw_max)
          if (info_raw == 0) then
            if (do_heavy_diag) then
              call wpw_local_herm_max(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), nred, sred_herm_tmp)
              call wpw_local_herm_max(dg_frag%wpw_reduced_H(1:nred, 1:nred, ispin, i_local), nred, hred_herm_tmp)
              local_max_sred_herm = max(local_max_sred_herm, sred_herm_tmp)
              local_max_hred_herm = max(local_max_hred_herm, hred_herm_tmp)
              sred_min_tmp = sraw_min
              sred_max_tmp = sraw_max
              local_min_sred_eval = min(local_min_sred_eval, sred_min_tmp)
              local_max_sred_eval = max(local_max_sred_eval, sred_max_tmp)
              if (sred_max_tmp > 0.0d0) &
                local_max_sred_cond = max(local_max_sred_cond, sred_max_tmp / max(sred_min_tmp, 1.0d-300))
            end if
            if (do_heavy_diag) then
              raw_keep_count = raw_keep_count + 1
              local_sum(21) = local_sum(21) + dble(nkeep_raw)
              local_sum(22) = local_sum(22) + 1.0d0
            end if
            if (sraw_max > 0.0d0) then
              cond_raw = max(cond_raw, sraw_max / max(sraw_min, 1.0d-8 * sraw_max, 1.0d-14))
            else
              cond_raw = max(cond_raw, huge(1.0d0))
            end if
            canon_norm_sum = 0.0d0
            canon_roundtrip_max = 0.0d0
            canon_w_diff_max = 0.0d0
            canon_pself_diff_max = 0.0d0
            canon_pneighbor_diff_max = 0.0d0
            canon_pack_bad = 0
            canon_pack_count = 0
            canon_recon_max_abs_diff = 0.0d0
            canon_recon_norm_diff_grid = 0.0d0
            canon_recon_prod_norm_grid = 0.0d0
            canon_recon_can_norm_grid = 0.0d0
            canon_recon_density_int_diff = 0.0d0
            canon_recon_rms_accum = 0.0d0
            canon_recon_count = 0
            canon_recon_bad = 0
            phase2e_can_w_norm = 0.0d0
            phase2e_can_pself_norm = 0.0d0
            phase2e_can_pneighbor_norm = 0.0d0
            phase2e_can_total_norm = 0.0d0
            phase2e_source_from_prodcoef = .false.
            phase2e_source_from_reproject_coef = .false.
            phase2e_source_from_density_reconstruction = .false.
            do ist = 1, nocc_spin
              occ = 1.0d0
              if (allocated(system%rocc)) occ = max(0.0d0, system%rocc(ist, 1, ispin))
              if (occ <= 0.0d0) cycle
              state_idx = (ispin - 1) * nstate_store + ist
              if ((use_canonical_reproject .or. canon_pack_diag .or. canon_grid_diag) .and. &
                  info_build == 0 .and. allocated(rhs_canonical)) then
                coef_prod(:) = (0.0d0, 0.0d0)
                do ib = 1, nbf
                  if (dg_frag%index_basis(ib, ifrag, ispin) <= 0) cycle
                  if (dg_frag%index_basis(ib, ifrag, ispin) > dg_frag%n_mat_max) cycle
                  if (dg_frag%coef_global_to_local(dg_frag%index_basis(ib, ifrag, ispin), ispin) <= 0) cycle
                  if (dg_frag%coef_global_to_local(dg_frag%index_basis(ib, ifrag, ispin), ispin) > size(dg_frag%coef, 1)) cycle
                  coef_prod(ib) = dg_frag%coef(dg_frag%coef_global_to_local( &
                    dg_frag%index_basis(ib, ifrag, ispin), ispin), ist, ispin)
                end do
                ! Experimental partial-canonical route: only the W block is
                ! sourced from production coefficients. Keep this off by
                ! default; the production path needs a full W/P/neighbor
                ! canonical packed object before density/current use.
                rhs_prod(:) = matmul(S_prod, coef_prod)
                rhs_canonical(:) = rhs_raw(:, ist)
                rhs_canonical(1:n_w) = matmul(conjg(transpose( &
                  dg_frag%global_wannier_local_coef(1:nbf, 1:n_w, ispin, i_local))), rhs_prod(:))
                c_can(:) = matmul(S_build_inv, rhs_canonical)
                if (canon_pack_diag) then
                  tmp_can(:) = matmul(S_build, c_can)
                  canon_norm = real(sum(conjg(c_can(:)) * tmp_can(:)), kind=8)
                  c_red_can(:) = (0.0d0, 0.0d0)
                  c_red_can(1:nred_eff) = matmul(S_red_inv(1:nred_eff,1:nred_eff), matmul( &
                    conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred_eff, i_local))), tmp_can(:)))
                  tmp_red_can(:) = (0.0d0, 0.0d0)
                  tmp_red_can(1:nred_eff) = matmul(dg_frag%wpw_reduced_S(1:nred_eff, 1:nred_eff, ispin, i_local), &
                    c_red_can(1:nred_eff))
                  c_can_back(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred_eff, i_local), &
                    c_red_can(1:nred_eff))
                  c_can_diff(:) = c_can(:) - c_can_back(:)
                  tmp_can(:) = matmul(S_build, c_can_diff)
                  canon_roundtrip = sqrt(max(0.0d0, real(sum(conjg(c_can_diff(:)) * tmp_can(:)), kind=8)))
                  if (n_w > 0) then
                    tmp_can(1:n_w) = matmul(S_build(1:n_w, 1:n_w), c_can_diff(1:n_w))
                    canon_w_diff = sqrt(max(0.0d0, real(sum(conjg(c_can_diff(1:n_w)) * tmp_can(1:n_w)), kind=8)))
                  else
                    canon_w_diff = 0.0d0
                  end if
                  if (n_pw > 0 .and. n_w + n_pw <= nraw) then
                    tmp_can(n_w+1:n_w+n_pw) = matmul(S_build(n_w+1:n_w+n_pw, n_w+1:n_w+n_pw), &
                      c_can_diff(n_w+1:n_w+n_pw))
                    canon_pself_diff = sqrt(max(0.0d0, real(sum(conjg(c_can_diff(n_w+1:n_w+n_pw)) * &
                      tmp_can(n_w+1:n_w+n_pw)), kind=8)))
                  else
                    canon_pself_diff = 0.0d0
                  end if
                  if (n_w + n_pw < nraw) then
                    tmp_can(n_w+n_pw+1:nraw) = matmul(S_build(n_w+n_pw+1:nraw, n_w+n_pw+1:nraw), &
                      c_can_diff(n_w+n_pw+1:nraw))
                    canon_pneighbor_diff = sqrt(max(0.0d0, real(sum(conjg(c_can_diff(n_w+n_pw+1:nraw)) * &
                      tmp_can(n_w+n_pw+1:nraw)), kind=8)))
                  else
                    canon_pneighbor_diff = 0.0d0
                  end if
                  canon_pack_count = canon_pack_count + 1
                  canon_norm_sum = canon_norm_sum + canon_norm
                  canon_roundtrip_max = max(canon_roundtrip_max, canon_roundtrip)
                  canon_w_diff_max = max(canon_w_diff_max, canon_w_diff)
                  canon_pself_diff_max = max(canon_pself_diff_max, canon_pself_diff)
                  canon_pneighbor_diff_max = max(canon_pneighbor_diff_max, canon_pneighbor_diff)
                  if (canon_norm /= canon_norm .or. canon_roundtrip /= canon_roundtrip .or. &
                      canon_w_diff /= canon_w_diff .or. canon_pself_diff /= canon_pself_diff .or. &
                      canon_pneighbor_diff /= canon_pneighbor_diff) canon_pack_bad = canon_pack_bad + 1
                end if
                if (use_canonical_reproject) then
                  c_raw(:) = c_can(:)
                else
                  c_raw(:) = matmul(S_raw_inv, rhs_raw(:, ist))
                end if
              else
                c_raw(:) = matmul(S_raw_inv, rhs_raw(:, ist))
              end if
              tmp_raw(:) = matmul(S_raw, c_raw)
              if (do_heavy_diag) state_local(9, state_idx) = state_local(9, state_idx) + &
                occ * real(sum(conjg(c_raw(:)) * tmp_raw(:)), kind=8)
              c_red_proj(:) = (0.0d0, 0.0d0)
              c_red_proj(1:nred_eff) = matmul(S_red_inv(1:nred_eff,1:nred_eff), &
                matmul(conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred_eff, i_local))), tmp_raw))
              tmp_red(:) = (0.0d0, 0.0d0)
              tmp_red(1:nred_eff) = matmul(dg_frag%wpw_reduced_S(1:nred_eff, 1:nred_eff, ispin, i_local), &
                c_red_proj(1:nred_eff))
              raw_back(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred_eff, i_local), &
                c_red_proj(1:nred_eff))
              raw_resid(:) = c_raw(:) - raw_back(:)
              if (do_heavy_diag) then
                local_sum(23) = local_sum(23) + occ * real(sum(conjg(c_raw(:)) * tmp_raw(:)), kind=8)
                local_sum(24) = local_sum(24) + occ * real(sum(conjg(c_red_proj(:)) * tmp_red(:)), kind=8)
                local_sum(25) = local_sum(25) + occ * real(sum(conjg(raw_resid(:)) * matmul(S_raw, raw_resid(:))), kind=8)
              end if
              if (info_build == 0) then
                tmp_raw_build(:) = matmul(S_build, c_raw)
                c_red_build(:) = (0.0d0, 0.0d0)
                c_red_build(1:nred_eff) = matmul(S_red_inv(1:nred_eff,1:nred_eff), matmul( &
                  conjg(transpose(dg_frag%wpw_reduced_transform(1:nraw, 1:nred_eff, i_local))), tmp_raw_build))
                tmp_red_build(:) = (0.0d0, 0.0d0)
                tmp_red_build(1:nred_eff) = matmul(dg_frag%wpw_reduced_S(1:nred_eff, 1:nred_eff, ispin, i_local), &
                  c_red_build(1:nred_eff))
                raw_back_build(:) = matmul(dg_frag%wpw_reduced_transform(1:nraw, 1:nred_eff, i_local), &
                  c_red_build(1:nred_eff))
                raw_resid_build(:) = c_raw(:) - raw_back_build(:)
                if (do_heavy_diag) then
                  local_sum(34) = local_sum(34) + occ * real(sum(conjg(c_raw(:)) * tmp_raw_build(:)), kind=8)
                  local_sum(35) = local_sum(35) + occ * real(sum(conjg(c_red_build(:)) * tmp_red_build(:)), kind=8)
                  local_sum(36) = local_sum(36) + occ * real(sum(conjg(raw_resid_build(:)) * &
                    matmul(S_build, raw_resid_build(:))), kind=8)
                end if
                c_red_hybrid(:) = (0.0d0, 0.0d0)
                raw_back_hybrid(:) = (0.0d0, 0.0d0)
                if (use_canonical_reproject) then
                  call build_wpw_reduced_c_can_reference(S_hybrid, S_red_inv(1:nred_eff,1:nred_eff), &
                    dg_frag%wpw_reduced_transform(1:nraw, 1:nred_eff, i_local), &
                    c_raw, nraw, nred_eff, raw_back_hybrid, c_red_hybrid, info_build)
                else
                  call build_wpw_reduced_raw_back_hybrid_from_inverse(S_raw_inv, rhs_raw(1:nraw, ist), &
                    S_hybrid, S_red_inv(1:nred_eff,1:nred_eff), &
                    dg_frag%wpw_reduced_transform(1:nraw, 1:nred_eff, i_local), &
                    nraw, nred_eff, raw_back_hybrid, c_red_hybrid, c_raw, info_build)
                end if
                tmp_red_hybrid(:) = (0.0d0, 0.0d0)
                tmp_red_hybrid(1:nred_eff) = matmul(dg_frag%wpw_reduced_S(1:nred_eff, 1:nred_eff, ispin, i_local), &
                  c_red_hybrid(1:nred_eff))
                if (sample_u_diag .and. nred_eff == nred) c_reproj_curr(1:nred, ist) = c_red_hybrid(1:nred)
                if (do_reproject) then
                  c_red(:) = (0.0d0, 0.0d0)
                  c_red(1:nself) = dg_frag%coef_wpw_self(1:nself, ist, ispin, i_local)
                  if (nneigh > 0) c_red(nself+1:nred) = &
                    dg_frag%coef_wpw_neighbor_reduced(1:nneigh, ist, ispin, i_local)
                  tmp_red(:) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), c_red(:))
                  norm_current_sred = real(sum(conjg(c_red(:)) * tmp_red(:)), kind=8)
                  norm_reproject_sred = real(sum(conjg(c_red_hybrid(:)) * tmp_red_hybrid(:)), kind=8)
                  overlap_sred = sum(conjg(c_red(:)) * tmp_red_hybrid(:))
                  diff_sred = max(0.0d0, norm_current_sred + norm_reproject_sred - &
                    2.0d0 * real(overlap_sred, kind=8))
                  local_sum(56) = local_sum(56) + occ * norm_current_sred
                  local_sum(57) = local_sum(57) + occ * norm_reproject_sred
                  local_sum(58) = local_sum(58) + occ * (norm_current_sred - norm_reproject_sred)
                  local_sum(62) = local_sum(62) + occ * diff_sred
                  local_sum(63) = local_sum(63) + occ * real(overlap_sred, kind=8)
                  local_sum(64) = local_sum(64) + occ * aimag(overlap_sred)
                  local_sum(65) = local_sum(65) + occ * abs(overlap_sred)
                  local_max_coef_diff = max(local_max_coef_diff, maxval(abs(c_red(1:nred) - c_red_hybrid(1:nred))))
                end if
                if (state_prop_diag .and. nred_eff == nred) then
                  if (state_prop_before) then
                    dg_frag%wpw_reproject_prev_coef(1:nred, ist, ispin, i_local) = c_red_hybrid(1:nred)
                    local_sum(66) = local_sum(66) + occ * norm_reproject_sred
                    local_sum(70) = local_sum(70) + occ
                  else if (state_prop_after .and. dg_frag%wpw_reproject_prev_valid .and. &
                           allocated(dg_frag%wpw_reduced_eval) .and. allocated(dg_frag%wpw_reduced_evec)) then
                    c_prev(1:nred) = dg_frag%wpw_reproject_prev_coef(1:nred, ist, ispin, i_local)
                    tmp_prev(1:nred) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), c_prev(1:nred))
                    eig_amp(1:nred) = matmul(conjg(transpose(dg_frag%wpw_reduced_evec(1:nred, 1:nred, ispin, i_local))), &
                      tmp_prev(1:nred))
                    do iw = 1, nred
                      phase_static = exp(cmplx(0.0d0, -dg_frag%wpw_reduced_eval(iw, ispin, i_local) * dt_use, kind=8))
                      eig_amp(iw) = phase_static * eig_amp(iw)
                    end do
                    c_pred(1:nred) = matmul(dg_frag%wpw_reduced_evec(1:nred, 1:nred, ispin, i_local), eig_amp(1:nred))
                    c_pred_resid(1:nred) = c_red_hybrid(1:nred) - c_pred(1:nred)
                    tmp_pred(1:nred) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), c_pred_resid(1:nred))
                    local_sum(66) = local_sum(66) + occ * norm_reproject_sred
                    local_sum(67) = local_sum(67) + occ * &
                      real(sum(conjg(c_pred_resid(1:nred)) * tmp_pred(1:nred)), kind=8)
                    local_sum(69) = local_sum(69) + occ * real(sum(conjg(c_pred(1:nred)) * &
                      matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), c_pred(1:nred))), kind=8)
                    local_sum(70) = local_sum(70) + occ
                  end if
                end if
                if (canon_grid_diag .and. allocated(c_can_store)) c_can_store(1:nraw, ist) = raw_back_hybrid(1:nraw)
                if (canon_grid_diag .and. allocated(c_can_store)) then
                  phase2e_source_from_reproject_coef = .true.
                  phase2e_source_from_density_reconstruction = .true.
                  phase2e_source_from_prodcoef = use_canonical_reproject .or. canon_pack_diag .or. canon_grid_diag
                end if
                raw_resid_hybrid(:) = c_raw(:) - raw_back_hybrid(:)
                if (state_prop_diag .and. state_prop_after .and. dg_frag%wpw_reproject_prev_valid .and. &
                    nred_eff == nred) then
                  local_sum(68) = local_sum(68) + occ * real(sum(conjg(raw_resid_hybrid(:)) * &
                    matmul(S_hybrid, raw_resid_hybrid(:))), kind=8)
                end if
                if (do_reproject) raw_coef_reproj(1:nraw, ist) = raw_back_hybrid(1:nraw)
                if (do_heavy_diag) then
                  local_sum(50) = local_sum(50) + occ * real(sum(conjg(c_raw(:)) * matmul(S_hybrid, c_raw)), kind=8)
                  local_sum(51) = local_sum(51) + occ * real(sum(conjg(c_red_hybrid(:)) * tmp_red_hybrid(:)), kind=8)
                  local_sum(52) = local_sum(52) + occ * real(sum(conjg(raw_resid_hybrid(:)) * &
                    matmul(S_hybrid, raw_resid_hybrid(:))), kind=8)
                end if
              else
                if (do_reproject) raw_coef_reproj(1:nraw, ist) = raw_back(1:nraw)
              end if
            end do
            if (canon_pack_diag .and. dg_frag%id == dg_frag%id_array(ifrag)) then
              write(*,'(1x,a,7(a,i0),7(a,1pe12.4),2(a,i0),a,l1,a,a)') &
                '[DG-WPW-RED-CANON-PACK-DIAG]', &
                ' step=', diag_step, ' frag=', ifrag, ' spin=', ispin, &
                ' n_state=', canon_pack_count, ' n_W=', n_w, ' n_P_self=', n_pw, &
                ' n_P_neighbor=', max(0, nraw - nself), &
                ' norm_C_can_avg=', canon_norm_sum / dble(max(1, canon_pack_count)), &
                ' roundtrip_diff_Snorm=', canon_roundtrip_max, &
                ' W_roundtrip_diff=', canon_w_diff_max, &
                ' Pself_roundtrip_diff=', canon_pself_diff_max, &
                ' Pneighbor_roundtrip_diff=', canon_pneighbor_diff_max, &
                ' metric_norm=', sbuild_max, &
                ' metric_min=', sbuild_min, &
                ' neighbor_mapping_false=', 0, ' bad_count=', canon_pack_bad, &
                ' bad=', canon_pack_bad /= 0, ' metric_tag=S_build'
            end if
            if (canon_grid_diag .and. allocated(c_can_store)) then
              do ist = 1, nocc_spin
                occ = 1.0d0
                if (allocated(system%rocc)) occ = max(0.0d0, system%rocc(ist, 1, ispin))
                if (occ <= 0.0d0) cycle
                tmp_can(:) = matmul(S_build, c_can_store(1:nraw, ist))
                phase2e_can_total_norm = phase2e_can_total_norm + occ * &
                  real(sum(conjg(c_can_store(1:nraw, ist)) * tmp_can(:)), kind=8)
                if (n_w > 0) then
                  tmp_can(1:n_w) = matmul(S_build(1:n_w, 1:n_w), c_can_store(1:n_w, ist))
                  phase2e_can_w_norm = phase2e_can_w_norm + occ * &
                    real(sum(conjg(c_can_store(1:n_w, ist)) * tmp_can(1:n_w)), kind=8)
                end if
                if (nself > n_w) then
                  tmp_can(n_w+1:nself) = matmul(S_build(n_w+1:nself, n_w+1:nself), &
                    c_can_store(n_w+1:nself, ist))
                  phase2e_can_pself_norm = phase2e_can_pself_norm + occ * &
                    real(sum(conjg(c_can_store(n_w+1:nself, ist)) * tmp_can(n_w+1:nself)), kind=8)
                end if
                if (nraw > nself) then
                  tmp_can(nself+1:nraw) = matmul(S_build(nself+1:nraw, nself+1:nraw), &
                    c_can_store(nself+1:nraw, ist))
                  phase2e_can_pneighbor_norm = phase2e_can_pneighbor_norm + occ * &
                    real(sum(conjg(c_can_store(nself+1:nraw, ist)) * tmp_can(nself+1:nraw)), kind=8)
                end if
              end do
              if (dg_frag%id == dg_frag%id_array(ifrag)) then
                write(*,'(1x,a,4(a,i0),4(a,1pe16.8),7(a,l1),4a)') &
                  '[DG-MIXEDZ-PHASE2E-CAN-SOURCE-DIAG]', &
                  ' step=', diag_step, &
                  ' frag=', ifrag, &
                  ' spin=', ispin, &
                  ' n_state=', nocc_spin, &
                  ' C_can_W_norm=', phase2e_can_w_norm, &
                  ' C_can_Pself_norm=', phase2e_can_pself_norm, &
                  ' C_can_Pneighbor_norm=', phase2e_can_pneighbor_norm, &
                  ' C_can_total_norm=', phase2e_can_total_norm, &
                  ' source_from_prodcoef=', phase2e_source_from_prodcoef, &
                  ' source_from_reproject_coef=', phase2e_source_from_reproject_coef, &
                  ' source_from_density_reconstruction=', phase2e_source_from_density_reconstruction, &
                  ' canon_grid_diag=', canon_grid_diag, &
                  ' canon_pack_diag=', canon_pack_diag, &
                  ' bad=', phase2e_can_total_norm /= phase2e_can_total_norm, &
                  ' production_replacement=', canon_replace_density_active, &
                  ' source_array_name=', 'c_can_store/raw_back_hybrid', &
                  ' code_path=', 'Phase2e-S_hybrid-reproject-C_can'
              end if
            end if
            if (canon_grid_diag .and. allocated(c_can_store)) then
              if (canon_replace_density_active .and. .not. canon_replace_density_zeroed) then
                rho_prod%f(:, :, :) = 0.0d0
                do iw = 1, system%nspin
                  rho_s_prod(iw)%f(:, :, :) = 0.0d0
                end do
                canon_replace_density_zeroed = .true.
              end if
              do iz = 1, dg_frag%nxyz_domain(3, ifrag)
                gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = 1, dg_frag%nxyz_domain(2, ifrag)
                  gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  do ix = 1, dg_frag%nxyz_domain(1, ifrag)
                    gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    rho_can_replace_grid = 0.0d0
                    do ist = 1, nocc_spin
                      occ = 1.0d0
                      if (allocated(system%rocc)) occ = max(0.0d0, system%rocc(ist, 1, ispin))
                      if (occ <= 0.0d0) cycle
                      psi_prod = (0.0d0, 0.0d0)
                      do ib = 1, nbf
                        if (dg_frag%index_basis(ib, ifrag, ispin) <= 0) cycle
                        if (dg_frag%index_basis(ib, ifrag, ispin) > dg_frag%n_mat_max) cycle
                        if (dg_frag%coef_global_to_local(dg_frag%index_basis(ib, ifrag, ispin), ispin) <= 0) cycle
                        if (dg_frag%coef_global_to_local(dg_frag%index_basis(ib, ifrag, ispin), ispin) > &
                            size(dg_frag%coef, 1)) cycle
                        if (allocated(dg_frag%phi_frag_c)) then
                          psi_prod = psi_prod + dg_frag%coef(dg_frag%coef_global_to_local( &
                            dg_frag%index_basis(ib, ifrag, ispin), ispin), ist, ispin) * &
                            dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                        else
                          psi_prod = psi_prod + dg_frag%coef(dg_frag%coef_global_to_local( &
                            dg_frag%index_basis(ib, ifrag, ispin), ispin), ist, ispin) * &
                            cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                        end if
                      end do
                      call reconstruct_psi_from_C_can(dg_frag, ispin, i_local, nbf, n_w, n_pw, n_pfrag, &
                        pfrag_ids(1:n_pfrag), c_can_store(1:nraw, ist), gx, gy, gz, bx, by, bz, psi_can)
                      psi_diff = psi_can - psi_prod
                      rho_can_grid = occ * abs(psi_can)**2
                      rho_prod_grid = occ * abs(psi_prod)**2
                      canon_density_grid_diff = rho_can_grid - rho_prod_grid
                      rho_can_replace_grid = rho_can_replace_grid + rho_can_grid
                      canon_recon_max_abs_diff = max(canon_recon_max_abs_diff, abs(psi_diff))
                      canon_recon_norm_diff_grid = canon_recon_norm_diff_grid + occ * abs(psi_diff)**2 * vol_weight
                      canon_recon_prod_norm_grid = canon_recon_prod_norm_grid + occ * abs(psi_prod)**2 * vol_weight
                      canon_recon_can_norm_grid = canon_recon_can_norm_grid + occ * abs(psi_can)**2 * vol_weight
                      canon_recon_density_int_diff = canon_recon_density_int_diff + canon_density_grid_diff * vol_weight
                      canon_recon_rms_accum = canon_recon_rms_accum + abs(psi_diff)**2
                      canon_recon_count = canon_recon_count + 1
                      if (canon_obs_diag .or. canon_hook_dryrun .or. canon_replace_density_active .or. canon_use_pz) then
                        if (mod(dg_frag%lgnum_total(3), 2) == 0) then
                          z_coord = dble(gz) - 0.5d0
                        else
                          z_coord = dble(gz)
                        end if
                        local_sum(80) = local_sum(80) + rho_prod_grid * vol_weight
                        local_sum(81) = local_sum(81) + rho_can_grid * vol_weight
                        local_sum(82) = local_sum(82) + canon_density_grid_diff * vol_weight
                        local_sum(83) = local_sum(83) + canon_density_grid_diff**2 * vol_weight
                        local_sum(84) = local_sum(84) + occ * abs(psi_diff)**2 * vol_weight
                        local_sum(85) = local_sum(85) + rho_prod_grid * z_coord * vol_weight
                        local_sum(86) = local_sum(86) + rho_can_grid * z_coord * vol_weight
                        local_sum(87) = local_sum(87) + vol_weight
                        local_max_canon_density_diff = max(local_max_canon_density_diff, abs(canon_density_grid_diff))
                      end if
                      if (abs(psi_diff) /= abs(psi_diff) .or. abs(psi_can) /= abs(psi_can) .or. &
                          abs(psi_prod) /= abs(psi_prod) .or. &
                          canon_density_grid_diff /= canon_density_grid_diff) canon_recon_bad = canon_recon_bad + 1
                    end do
                    if (canon_replace_density_active) then
                      rho_s_prod(ispin)%f(gx, gy, gz) = rho_s_prod(ispin)%f(gx, gy, gz) + rho_can_replace_grid
                      rho_prod%f(gx, gy, gz) = rho_prod%f(gx, gy, gz) + rho_can_replace_grid
                    end if
                  end do
                end do
              end do
              if (canon_obs_diag .or. canon_hook_dryrun .or. canon_replace_density_active .or. canon_use_pz) &
                local_sum(88) = local_sum(88) + dble(canon_recon_bad)
              if (dg_frag%id == dg_frag%id_array(ifrag)) then
                write(*,'(1x,a,6(a,i0),8(a,1pe12.4),a,i0,a,l1,a,a)') &
                  '[DG-WPW-RED-CANON-RECON-DIAG]', &
                  ' step=', diag_step, ' frag=', ifrag, ' spin=', ispin, &
                  ' n_state=', nocc_spin, ' n_W=', n_w, ' n_P_total=', nraw - n_w, &
                  ' max_abs_diff=', canon_recon_max_abs_diff, &
                  ' rms_diff=', sqrt(canon_recon_rms_accum / dble(max(1, canon_recon_count))), &
                  ' norm_diff_grid=', sqrt(max(0.0d0, canon_recon_norm_diff_grid)), &
                  ' rel_norm_diff_grid=', sqrt(max(0.0d0, canon_recon_norm_diff_grid)) / &
                    max(sqrt(max(0.0d0, canon_recon_prod_norm_grid)), 1.0d-300), &
                  ' density_integral_diff=', canon_recon_density_int_diff, &
                  ' prod_norm_grid=', canon_recon_prod_norm_grid, &
                  ' can_norm_grid=', canon_recon_can_norm_grid, &
                  ' count=', dble(canon_recon_count), &
                  ' bad_count=', canon_recon_bad, ' bad=', canon_recon_bad /= 0, &
                  ' route=C_can_to_grid metric_tag=S_build'
              end if
            end if
            if (sample_u_diag .and. dg_frag%wpw_reproject_prev_valid .and. nred_eff == nred) then
              gram_sample(:, :) = (0.0d0, 0.0d0)
              do jw = 1, nocc_spin
                tmp_pred(1:nred) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
                  dg_frag%wpw_reproject_prev_coef(1:nred, jw, ispin, i_local))
                do iw = 1, nocc_spin
                  gram_sample(iw, jw) = sum(conjg(dg_frag%wpw_reproject_prev_coef(1:nred, iw, ispin, i_local)) * &
                    tmp_pred(1:nred))
                end do
              end do
              call hermitize_matrix(gram_sample, nocc_spin)
              call build_hermitian_pseudoinverse(gram_sample, nocc_spin, sample_tol, gram_inv_sample, &
                info_build, sample_smin, sample_smax, nkeep_build)
              if (info_build == 0) then
                if (sample_smax > 0.0d0) then
                  sample_cond = sample_smax / max(sample_smin, sample_tol * sample_smax, 1.0d-300)
                else
                  sample_cond = 0.0d0
                end if
                local_max_sample_cond = max(local_max_sample_cond, sample_cond)
                local_max_sample_rank = max(local_max_sample_rank, dble(nkeep_build))
                local_min_sample_rank = min(local_min_sample_rank, dble(nkeep_build))
                do ist = 1, nocc_spin
                  occ = 1.0d0
                  if (allocated(system%rocc)) occ = max(0.0d0, system%rocc(ist, 1, ispin))
                  if (occ <= 0.0d0) cycle
                  sample_rhs(1:nocc_spin) = gram_sample(1:nocc_spin, ist)
                  sample_weight(1:nocc_spin) = matmul(gram_inv_sample(1:nocc_spin,1:nocc_spin), sample_rhs(1:nocc_spin))
                  sample_pred(1:nred) = matmul(c_reproj_curr(1:nred,1:nocc_spin), sample_weight(1:nocc_spin))
                  sample_resid(1:nred) = c_reproj_curr(1:nred, ist) - sample_pred(1:nred)
                  sample_tmp(1:nred) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
                    c_reproj_curr(1:nred, ist))
                  sample_next_norm = real(sum(conjg(c_reproj_curr(1:nred, ist)) * sample_tmp(1:nred)), kind=8)
                  sample_overlap = sum(conjg(sample_pred(1:nred)) * sample_tmp(1:nred))
                  sample_tmp(1:nred) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
                    sample_pred(1:nred))
                  sample_pred_norm = real(sum(conjg(sample_pred(1:nred)) * sample_tmp(1:nred)), kind=8)
                  sample_tmp(1:nred) = matmul(dg_frag%wpw_reduced_S(1:nred, 1:nred, ispin, i_local), &
                    sample_resid(1:nred))
                  sample_resid_norm = real(sum(conjg(sample_resid(1:nred)) * sample_tmp(1:nred)), kind=8)
                  local_sum(71) = local_sum(71) + occ * sample_next_norm
                  local_sum(72) = local_sum(72) + occ * sample_resid_norm
                  local_sum(73) = local_sum(73) + occ
                  local_sum(78) = local_sum(78) + occ * sample_pred_norm
                  local_sum(79) = local_sum(79) + occ * real(sample_overlap, kind=8)
                end do
                local_sum(74) = local_sum(74) + dble(nkeep_build)
                local_sum(75) = local_sum(75) + 1.0d0
                local_sum(76) = local_sum(76) + sample_cond
              else
                local_sum(77) = local_sum(77) + 1.0d0
              end if
            end if
          end if
        end if
        end if
        if (do_reproject) then
          if (replace_density_active .and. .not. replace_density_zeroed) then
            rho_prod%f(:, :, :) = 0.0d0
            do iw = 1, system%nspin
              rho_s_prod(iw)%f(:, :, :) = 0.0d0
            end do
            replace_density_zeroed = .true.
          end if
          do iz = 1, dg_frag%nxyz_domain(3, ifrag)
            gz = dg_frag%ixyz_frag(3, ifrag) + iz - 1
            bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
            if (bz < p_lb3 .or. bz > p_ub3) cycle
            do iy = 1, dg_frag%nxyz_domain(2, ifrag)
              gy = dg_frag%ixyz_frag(2, ifrag) + iy - 1
              by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
              if (by < p_lb2 .or. by > p_ub2) cycle
              do ix = 1, dg_frag%nxyz_domain(1, ifrag)
                gx = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                if (bx < p_lb1 .or. bx > p_ub1) cycle
                basis_raw(:) = (0.0d0, 0.0d0)
                wval(:) = (0.0d0, 0.0d0)
                do iw = 1, n_w
                  do ib = 1, nbf
                    if (allocated(dg_frag%phi_frag_c)) then
                      wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                        dg_frag%phi_frag_c(bx, by, bz, ib, i_local)
                    else
                      wval(iw) = wval(iw) + dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local) * &
                        cmplx(dg_frag%phi_frag(bx, by, bz, ib, i_local), 0.0d0, kind=8)
                    end if
                  end do
                end do
                basis_raw(1:n_w) = wval(1:n_w)
                do pidx = 1, n_pfrag
                  row0 = n_w + (pidx - 1) * n_pw
                  call wpw_normalized_window_at_grid(dg_frag, pfrag_ids(pidx), gx, gy, gz, chi, grad_chi)
                  do ipw = 1, n_pw
                    phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                                dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                                dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                    basis_raw(row0 + ipw) = chi * exp(cmplx(0.0d0, phase_arg, kind=8))
                  end do
                end do
                do ist = 1, nocc_spin
                  if (ist == 1) rho_reproject_grid = 0.0d0
                  occ = 1.0d0
                  if (allocated(system%rocc)) occ = max(0.0d0, system%rocc(ist, 1, ispin))
                  if (occ <= 0.0d0) cycle
                  psi_total = sum(raw_coef_reproj(1:nraw, ist) * basis_raw(1:nraw))
                  q_reproject = occ * abs(psi_total)**2 * vol_weight
                  rho_reproject_grid = rho_reproject_grid + occ * abs(psi_total)**2
                  local_sum(54) = local_sum(54) + q_reproject
                  if (mod(dg_frag%lgnum_total(3), 2) == 0) then
                    z_coord = dble(gz) - 0.5d0
                  else
                    z_coord = dble(gz)
                  end if
                  local_sum(61) = local_sum(61) + occ * abs(psi_total)**2 * z_coord * vol_weight
                  if (q_reproject /= q_reproject) bad_reproject = .true.
                end do
                if (replace_density_active .and. .not. replace_density_from_propagated) then
                  rho_s_prod(ispin)%f(gx, gy, gz) = rho_s_prod(ispin)%f(gx, gy, gz) + rho_reproject_grid
                  rho_prod%f(gx, gy, gz) = rho_prod%f(gx, gy, gz) + rho_reproject_grid
                end if
              end do
            end do
          end do
        end if
        if (do_projection_diag) then
          deallocate(S_raw, rhs_raw, S_raw_inv, S_red_inv)
          deallocate(S_build, S_hybrid, S_build_inv, bsb_build, bsb_density, bsb_hybrid)
          deallocate(c_raw, tmp_raw, tmp_raw_build)
          if (allocated(c_can)) deallocate(c_can, c_can_back, c_can_diff, tmp_can, c_red_can, tmp_red_can)
          if (allocated(c_can_store)) deallocate(c_can_store)
          deallocate(c_red_proj, c_red_build, c_red_hybrid, tmp_red, tmp_red_build, tmp_red_hybrid)
          deallocate(raw_back, raw_back_build, raw_back_hybrid, raw_resid, raw_resid_build, raw_resid_hybrid)
          if (allocated(c_prev)) deallocate(c_prev, c_pred, c_pred_resid, tmp_prev, tmp_pred, eig_amp)
          if (allocated(c_reproj_curr)) then
            deallocate(c_reproj_curr, gram_sample, gram_inv_sample)
            deallocate(sample_rhs, sample_weight, sample_pred, sample_resid, sample_tmp)
          end if
        end if
      end do
      if (do_heavy_diag .and. dg_frag%id == dg_frag%id_array(ifrag)) then
        write(*,'(1x,a,i0,5(a,1pe12.4),a,l1)') &
          '[DG-WPW-RED-DENSITY] fragment experimental:', ifrag, &
          ' charge=', frag_sum(1), ' self=', frag_sum(2), &
          ' neighbor=', frag_sum(3), ' cross=', frag_sum(4), &
          ' grid_occ_weight=', frag_sum(5), ' bad=', bad_density
      end if
      if (allocated(raw_coef_reproj)) deallocate(raw_coef_reproj)
      if (allocated(phi_val)) deallocate(phi_val, S_prod, coef_prod, rhs_prod, rhs_canonical)
      deallocate(c_red, raw_coef, basis_raw, wval)
    end do

    if (state_prop_diag .and. state_prop_before) dg_frag%wpw_reproject_prev_valid = .true.
    if (bad_density) local_sum(6) = 1.0d0
    if (bad_reproject) local_sum(55) = 1.0d0
    call comm_summation(local_sum, global_sum, size(local_sum), dg_frag%icomm)
    call comm_summation(state_local, state_global, size(state_local), dg_frag%icomm)
    max_in(:) = 0.0d0
    max_in(1) = local_max_delta
    max_in(2) = local_max_rho_self
    max_in(3) = local_max_prod_delta
    max_in(4) = local_max_rho_prod
    max_in(5) = cond_raw
    max_in(6) = local_max_sdiff
    max_in(7) = local_max_sref
    max_in(8) = local_max_coef_diff
    max_in(9) = local_max_hred_herm
    max_in(10) = local_max_sred_herm
    max_in(11) = local_max_sred_eval
    max_in(12) = local_max_sred_cond
    if (local_min_sred_eval < huge(1.0d0)) then
      max_in(13) = -local_min_sred_eval
    else
      max_in(13) = 0.0d0
    end if
    max_in(14) = local_max_sample_cond
    max_in(15) = local_max_sample_rank
    if (local_min_sample_rank < huge(1.0d0)) then
      max_in(16) = -local_min_sample_rank
    else
      max_in(16) = 0.0d0
    end if
    max_in(17) = local_max_canon_density_diff
    call comm_get_max(max_in, max_out, 20, dg_frag%icomm)
    if (dg_frag%id == 0) then
      rms_delta = sqrt(max(0.0d0, global_sum(8)))
      if (max_out(2) > 0.0d0) then
        delta_over_self = max_out(1) / max_out(2)
      else
        delta_over_self = 0.0d0
      end if
      if (global_sum(14) > 0.0d0) then
        phys_scale = global_sum(17) / global_sum(14)
      else
        phys_scale = 0.0d0
      end if
      if (do_drift_log) then
        write(*,'(1x,a,2(a,i0),4(a,1pe12.4),2(a,l1))') &
          '[DG-WPW-RED-DIAG-DENSITY-DRIFT]', &
          ' step=', diag_step, ' nstep=', diag_nstep, &
          ' int_prod_hvol=', global_sum(9) * phys_scale, &
          ' int_reduced_hvol=', global_sum(54) * phys_scale, &
          ' int_prod_minus_reduced_hvol=', (global_sum(9) - global_sum(54)) * phys_scale, &
          ' trace_occ_norm=', global_sum(54) * phys_scale, &
          ' bad_coef=', global_sum(53) > 0.5d0, &
          ' bad_density=', (global_sum(6) > 0.5d0) .or. (global_sum(55) > 0.5d0)
        write(*,'(1x,a,4(a,i0),5(a,1pe12.4),2(a,l1))') &
          '[DG-WPW-RED-DIAG-COMPRESS-SWEEP]', &
          ' step=', diag_step, ' keep_n=', keep_n_env, &
          ' total_nred=', nint(global_sum(20)), ' effective_nred=', nint(global_sum(59)), &
          ' int_prod_hvol=', global_sum(9) * phys_scale, &
          ' int_reduced_hvol=', global_sum(54) * phys_scale, &
          ' int_diff=', (global_sum(9) - global_sum(54)) * phys_scale, &
          ' rel_int_diff=', (global_sum(9) - global_sum(54)) / max(abs(global_sum(9)), 1.0d-300), &
          ' trace_occ_norm_reduced=', global_sum(54) * phys_scale, &
          ' bad_coef=', global_sum(53) > 0.5d0, &
          ' bad_density=', (global_sum(6) > 0.5d0) .or. (global_sum(55) > 0.5d0)
        write(*,'(1x,a,2(a,i0),5(a,1pe12.4))') &
          '[DG-WPW-RED-DIAG-PZ-CMP]', &
          ' step=', diag_step, ' keep_n=', keep_n_env, &
          ' Pz_prod=', global_sum(60) * phys_scale, &
          ' Pz_reduced=', global_sum(61) * phys_scale, &
          ' Pz_diff=', (global_sum(60) - global_sum(61)) * phys_scale, &
          ' rel_Pz_diff=', (global_sum(60) - global_sum(61)) / max(abs(global_sum(60)), 1.0d-300), &
          ' hvol_over_norm_weight=', phys_scale
        if (propagated_debug) then
          write(*,'(1x,a,2(a,i0),8(a,1pe12.4),3(a,l1))') &
            '[DG-WPW-RED-DIAG-DENSITY-PROPAGATED-DEBUG]', &
            ' step=', diag_step, ' nstep=', diag_nstep, &
            ' int_prod_hvol=', global_sum(9) * phys_scale, &
            ' int_reduced_hvol_current=', global_sum(1) * phys_scale, &
            ' diff_current=', (global_sum(9) - global_sum(1)) * phys_scale, &
            ' trace_occ_norm_current=', global_sum(1) * phys_scale, &
            ' int_reduced_hvol_reproject=', global_sum(54) * phys_scale, &
            ' diff_reproject=', (global_sum(9) - global_sum(54)) * phys_scale, &
            ' trace_occ_norm_reproject=', global_sum(54) * phys_scale, &
            ' hvol_over_norm_weight=', phys_scale, &
            ' bad_coef=', global_sum(53) > 0.5d0, &
            ' bad_density=', global_sum(6) > 0.5d0, &
            ' reproject=', do_reproject
          write(*,'(1x,a,2(a,i0),11(a,1pe12.4),a,l1)') &
            '[DG-WPW-RED-DIAG-COEF-NORM-DRIFT]', &
            ' step=', diag_step, ' nstep=', diag_nstep, &
            ' norm_current_Sred=', global_sum(56), &
            ' norm_reproject_Sred=', global_sum(57), &
            ' norm_diff_current_minus_reproject=', global_sum(58), &
            ' rel_norm_diff=', global_sum(58) / max(abs(global_sum(57)), 1.0d-300), &
            ' coef_diff_Snorm=', sqrt(max(0.0d0, global_sum(62))), &
            ' rel_coef_diff_Snorm=', sqrt(max(0.0d0, global_sum(62))) / &
              max(sqrt(max(0.0d0, global_sum(57))), 1.0d-300), &
            ' overlap_S_real=', global_sum(63), &
            ' overlap_S_imag=', global_sum(64), &
            ' overlap_S_abs_sum=', global_sum(65), &
            ' overlap_S_normed=', global_sum(63) / &
              max(sqrt(max(0.0d0, global_sum(56) * global_sum(57))), 1.0d-300), &
            ' max_abs_coef_diff=', max_out(8), &
            ' bad_coef=', global_sum(53) > 0.5d0
        end if
        if (state_prop_diag) then
          write(*,'(1x,a,2(a,i0),a,a,8(a,1pe12.4),4(a,l1))') &
            '[DG-WPW-RED-DIAG-STATE-PROP]', &
            ' step=', diag_step, ' nstep=', diag_nstep, &
            ' stage=', trim(reproject_stage), &
            ' norm_reproj=', global_sum(66), &
            ' norm_static_pred=', global_sum(69), &
            ' static_residual_Snorm=', sqrt(max(0.0d0, global_sum(67))), &
            ' rel_static_residual=', sqrt(max(0.0d0, global_sum(67))) / &
              max(sqrt(max(0.0d0, global_sum(66))), 1.0d-300), &
            ' basis_leakage_Snorm=', sqrt(max(0.0d0, global_sum(68))), &
            ' rel_basis_leakage=', sqrt(max(0.0d0, global_sum(68))) / &
              max(sqrt(max(0.0d0, global_sum(66))), 1.0d-300), &
            ' occ_weight=', global_sum(70), &
            ' dt=', dt_use, &
            ' prev_valid=', dg_frag%wpw_reproject_prev_valid, &
            ' after_stage=', state_prop_after, &
            ' static_available=', state_prop_after .and. dg_frag%wpw_reproject_prev_valid .and. &
              allocated(dg_frag%wpw_reduced_eval) .and. allocated(dg_frag%wpw_reduced_evec), &
            ' bad_coef=', global_sum(53) > 0.5d0
        end if
        if (sample_u_diag) then
          write(*,'(1x,a,2(a,i0),7(a,1pe12.4),4(a,i0),a,1pe12.4,a,l1)') &
            '[DG-WPW-RED-DIAG-SAMPLE-U]', &
            ' step=', diag_step, ' nstep=', diag_nstep, &
            ' norm_reproj_next=', global_sum(71), &
            ' sample_residual_Snorm=', sqrt(max(0.0d0, global_sum(72))), &
            ' rel_sample_residual=', sqrt(max(0.0d0, global_sum(72))) / &
              max(sqrt(max(0.0d0, global_sum(71))), 1.0d-300), &
            ' occ_weight=', global_sum(73), &
            ' avg_sample_rank=', global_sum(74) / max(global_sum(75), 1.0d0), &
            ' avg_sample_cond=', global_sum(76) / max(global_sum(75), 1.0d0), &
            ' sample_tol=', sample_tol, &
            ' failed_blocks=', nint(global_sum(77)), &
            ' fitted_blocks=', nint(global_sum(75)), &
            ' min_sample_rank=', max(0, nint(-max_out(16))), &
            ' max_sample_rank=', nint(max_out(15)), &
            ' max_sample_cond=', max_out(14), &
            ' full_reduced_basis_action=', .false.
          write(*,'(1x,a,2(a,i0),a,a,5(a,l1),8(a,1pe12.4),2(a,l1))') &
            '[DG-WPW-RED-DIAG-PRODOP-ACTION]', &
            ' step=', diag_step, ' nred=', nint(global_sum(59)), &
            ' route=', trim(prodop_route_label), &
            ' field_on=', prodop_field_flag, &
            ' kick_applied=', prodop_kick_flag, &
            ' mixed_z_included=', prodop_mixed_z_flag, &
            ' global_flux_included=', prodop_global_flux_flag, &
            ' predictor_corrector_included=', prodop_predictor_corrector_flag, &
            ' rel_action_residual=', sqrt(max(0.0d0, global_sum(72))) / &
              max(sqrt(max(0.0d0, global_sum(71))), 1.0d-300), &
            ' rel_coef_diff_Snorm=', sqrt(max(0.0d0, global_sum(72))) / &
              max(sqrt(max(0.0d0, global_sum(71))), 1.0d-300), &
            ' coef_diff_Snorm=', sqrt(max(0.0d0, global_sum(72))), &
            ' overlap_S_normed=', global_sum(79) / &
              max(sqrt(max(0.0d0, global_sum(78) * global_sum(71))), 1.0d-300), &
            ' Pz_diff=', (global_sum(60) - global_sum(61)) * phys_scale, &
            ' sample_tol=', sample_tol, &
            ' avg_sample_rank=', global_sum(74) / max(global_sum(75), 1.0d0), &
            ' max_sample_cond=', max_out(14), &
            ' sampled_state_action=', .true., &
            ' full_reduced_basis_action=', .false.
        end if
      end if
      if (canon_obs_diag) then
        write(*,'(1x,a,2(a,i0),13(a,1pe12.4),2(a,l1))') &
          '[DG-WPW-RED-CANON-OBS-DIAG]', &
          ' step=', diag_step, ' nstep=', diag_nstep, &
          ' density_integral_prod=', global_sum(80) * phys_scale, &
          ' density_integral_can=', global_sum(81) * phys_scale, &
          ' density_integral_diff=', global_sum(82) * phys_scale, &
          ' density_max_abs_diff=', max_out(17), &
          ' density_rms_diff=', sqrt(max(0.0d0, global_sum(83))), &
          ' psi_norm_diff_grid=', sqrt(max(0.0d0, global_sum(84))), &
          ' Pz_prod=', global_sum(85) * phys_scale, &
          ' Pz_can=', global_sum(86) * phys_scale, &
          ' Pz_diff=', (global_sum(86) - global_sum(85)) * phys_scale, &
          ' rel_Pz_diff=', (global_sum(86) - global_sum(85)) / max(abs(global_sum(85)), 1.0d-300), &
          ' current_prod=', 0.0d0, &
          ' current_can=', 0.0d0, &
          ' current_diff=', 0.0d0, &
          ' current_available=', .false., &
          ' bad=', (global_sum(88) > 0.5d0)
      end if
      if (canon_hook_dryrun) then
        write(*,'(1x,a,2(a,i0),13(a,1pe12.4),3(a,l1))') &
          '[DG-WPW-RED-CANON-HOOK-DRYRUN]', &
          ' step=', diag_step, ' nstep=', diag_nstep, &
          ' density_integral_prod=', global_sum(80) * phys_scale, &
          ' density_integral_can=', global_sum(81) * phys_scale, &
          ' density_integral_diff=', global_sum(82) * phys_scale, &
          ' density_max_abs_diff=', max_out(17), &
          ' density_rms_diff=', sqrt(max(0.0d0, global_sum(83))), &
          ' psi_norm_diff_grid=', sqrt(max(0.0d0, global_sum(84))), &
          ' Pz_prod=', global_sum(85) * phys_scale, &
          ' Pz_can=', global_sum(86) * phys_scale, &
          ' Pz_diff=', (global_sum(86) - global_sum(85)) * phys_scale, &
          ' rel_Pz_diff=', (global_sum(86) - global_sum(85)) / max(abs(global_sum(85)), 1.0d-300), &
          ' current_prod=', 0.0d0, &
          ' current_can=', 0.0d0, &
          ' current_diff=', 0.0d0, &
          ' current_available=', .false., &
          ' density_replaced=', .false., &
          ' bad=', (global_sum(88) > 0.5d0)
      end if
      if (canon_replace_density_active) then
        write(*,'(1x,a,2(a,i0),6(a,1pe12.4),2(a,l1))') &
          '[DG-WPW-RED-CANON-DENSITY-PRODUCTION]', &
          ' step=', diag_step, ' nstep=', diag_nstep, &
          ' density_integral_prod_before=', global_sum(80) * phys_scale, &
          ' density_integral_can=', global_sum(81) * phys_scale, &
          ' density_integral_diff=', global_sum(82) * phys_scale, &
          ' density_max_abs_diff=', max_out(17), &
          ' density_rms_diff=', sqrt(max(0.0d0, global_sum(83))), &
          ' psi_norm_diff_grid=', sqrt(max(0.0d0, global_sum(84))), &
          ' density_replaced=', (global_sum(88) <= 0.5d0), &
          ' bad=', (global_sum(88) > 0.5d0)
      end if
      if (canon_use_pz) then
        write(*,'(1x,a,2(a,i0),5(a,1pe12.4),2(a,l1))') &
          '[DG-WPW-RED-CANON-PZ-PRODUCTION]', &
          ' step=', diag_step, ' nstep=', diag_nstep, &
          ' Pz_prod_before=', global_sum(85) * phys_scale, &
          ' Pz_can=', global_sum(86) * phys_scale, &
          ' Pz_diff=', (global_sum(86) - global_sum(85)) * phys_scale, &
          ' rel_Pz_diff=', (global_sum(86) - global_sum(85)) / max(abs(global_sum(85)), 1.0d-300), &
          ' psi_norm_diff_grid=', sqrt(max(0.0d0, global_sum(84))), &
          ' Pz_candidate_available=', (global_sum(88) <= 0.5d0), &
          ' bad=', (global_sum(88) > 0.5d0)
      end if
      if (present(pz_prod_raw)) pz_prod_raw = global_sum(85) * phys_scale
      if (present(pz_can_raw)) pz_can_raw = global_sum(86) * phys_scale
      if (present(pz_bad)) pz_bad = global_sum(88) > 0.5d0
      if (do_pz_series .or. do_drift_log) then
        call append_wpw_reduced_pz_cmp(diag_step, diag_time, keep_n_env, &
          global_sum(60) * phys_scale, global_sum(61) * phys_scale)
      end if
    end if
    if (present(pz_prod_raw)) pz_prod_raw = global_sum(85) * phys_scale
    if (present(pz_can_raw)) pz_can_raw = global_sum(86) * phys_scale
    if (present(pz_bad)) pz_bad = global_sum(88) > 0.5d0
    if (.not. do_heavy_diag) then
      deallocate(state_local, state_global, state_loss_work)
      return
    end if
    if (dg_frag%id == 0) then
      if (do_heavy_diag) then
        write(*,'(1x,a,5(a,1pe12.4))') &
          '[DG-WPW-RED-DIAG-HS-CHECK]', &
          ' herm_Hred=', max_out(9), &
          ' herm_Sred=', max_out(10), &
          ' min_eval_Sred=', -max_out(13), &
          ' max_eval_Sred=', max_out(11), &
          ' cond_Sred=', max_out(12)
      end if
      write(*,'(1x,a,5(a,1pe12.4),a,l1)') &
        '[DG-WPW-RED-DENSITY] total experimental:', &
        ' charge=', global_sum(1), ' self=', global_sum(2), &
        ' neighbor=', global_sum(3), ' cross=', global_sum(4), &
        ' grid_occ_weight=', global_sum(5), ' bad=', global_sum(6) > 0.5d0
      write(*,'(1x,a,7(a,1pe12.4),a,l1)') &
        '[DG-WPW-RED-DENSITY-DELTA] total experimental:', &
        ' int_self_only=', global_sum(2), ' int_reduced=', global_sum(1), &
        ' int_delta=', global_sum(7), ' max_abs_delta=', max_out(1), &
        ' rms_delta=', rms_delta, ' max_abs_delta_over_self=', delta_over_self, &
        ' max_abs_rho_self=', max_out(2), ' bad=', global_sum(6) > 0.5d0
      if (present(rho_prod)) then
        if (global_sum(14) > 0.0d0) then
          phys_scale = global_sum(17) / global_sum(14)
        else
          phys_scale = 0.0d0
        end if
        rms_prod_delta = sqrt(max(0.0d0, global_sum(11)))
        if (max_out(4) > 0.0d0) then
          prod_delta_over_prod = max_out(3) / max_out(4)
        else
          prod_delta_over_prod = 0.0d0
        end if
        write(*,'(1x,a,7(a,1pe12.4),a,l1)') &
          '[DG-WPW-RED-DENSITY-PROD-CMP] total experimental:', &
          ' int_rho_prod_current=', global_sum(9), ' int_rho_reduced_diag=', global_sum(1), &
          ' int_diff=', global_sum(10), ' max_abs_grid_diff=', max_out(3), &
          ' rms_grid_diff=', rms_prod_delta, ' max_abs_grid_diff_over_prod=', prod_delta_over_prod, &
          ' max_abs_rho_prod=', max_out(4), ' bad=', global_sum(6) > 0.5d0
        write(*,'(1x,a,6(a,1pe12.4),a,l1)') &
          '[DG-WPW-RED-DENSITY-SPLIT] total experimental:', &
          ' int_prod_total=', global_sum(9), &
          ' int_diag_self_only=', global_sum(2), &
          ' int_diag_neighbor=', global_sum(3), &
          ' int_diag_cross=', global_sum(4), &
          ' int_diag_reduced=', global_sum(1), &
          ' int_prod_minus_self=', global_sum(9) - global_sum(2), &
          ' prod_split_available=', .false.
        write(*,'(1x,a,6(a,1pe12.4),a,i0,a,l1,a,l1)') &
          '[DG-WPW-RED-DENSITY-WEIGHT] total experimental:', &
          ' occ_weight_sum_diag=', global_sum(5), &
          ' prod_valid_grid_weight=', global_sum(12), &
          ' prod_bad_grid_weight=', global_sum(13), &
          ' diag_grid_weight=', global_sum(14), &
          ' int_prod_minus_reduced=', global_sum(9) - global_sum(1), &
          ' int_prod_minus_self=', global_sum(9) - global_sum(2), &
          ' nspin=', dg_frag%nspin, &
          ' occ_weight_sum_prod_available=', .false., &
          ' bad=', global_sum(6) > 0.5d0
        write(*,'(1x,a,8(a,1pe12.4),a,l1)') &
          '[DG-WPW-RED-DENSITY-PHYS] total experimental:', &
          ' hvol_over_norm_weight=', phys_scale, &
          ' int_prod_hvol=', global_sum(9) * phys_scale, &
          ' int_diag_self_hvol=', global_sum(2) * phys_scale, &
          ' int_diag_neighbor_hvol=', global_sum(3) * phys_scale, &
          ' int_diag_cross_hvol=', global_sum(4) * phys_scale, &
          ' int_diag_reduced_hvol=', global_sum(1) * phys_scale, &
          ' int_prod_minus_reduced_hvol=', (global_sum(9) - global_sum(1)) * phys_scale, &
          ' rocc_sum_diag_states=', global_sum(15) / dble(max(1, dg_frag%n_frag)), &
          ' bad=', global_sum(6) > 0.5d0
        write(*,'(1x,a,5(a,1pe12.4),a,l1)') &
          '[DG-WPW-RED-DENSITY-PROD-TRACE] total experimental:', &
          ' int_prod_total_norm=', global_sum(9), &
          ' int_prod_self_basis_like=', global_sum(2), &
          ' int_prod_self_minus_diag_self=', global_sum(9) - global_sum(2), &
          ' int_diag_self_only=', global_sum(2), &
          ' diag_state_count_weight=', global_sum(16) / dble(max(1, dg_frag%n_frag)), &
          ' bad=', global_sum(6) > 0.5d0
        sum_loss_self = 0.0d0
        sum_loss_red = 0.0d0
        max_loss_self = 0.0d0
        max_loss_red = 0.0d0
        min_norm_red = huge(1.0d0)
        max_norm_red = 0.0d0
        do i_state = 1, nstate_metric
          occ_state = state_global(1, i_state) / dble(max(1, dg_frag%n_frag))
          if (occ_state <= 0.0d0) cycle
          self_state_hvol = state_global(2, i_state) * phys_scale
          red_state_hvol = state_global(3, i_state) * phys_scale
          sum_loss_self = sum_loss_self + (occ_state - self_state_hvol)
          sum_loss_red = sum_loss_red + (occ_state - red_state_hvol)
          max_loss_self = max(max_loss_self, abs(occ_state - self_state_hvol))
          max_loss_red = max(max_loss_red, abs(occ_state - red_state_hvol))
          norm_state = red_state_hvol / max(occ_state, 1.0d-300)
          min_norm_red = min(min_norm_red, norm_state)
          max_norm_red = max(max_norm_red, norm_state)
        end do
        if (global_sum(15) > 0.0d0) then
          avg_loss_self = sum_loss_self / (global_sum(16) / dble(max(1, dg_frag%n_frag)))
          avg_loss_red = sum_loss_red / (global_sum(16) / dble(max(1, dg_frag%n_frag)))
        else
          avg_loss_self = 0.0d0
          avg_loss_red = 0.0d0
        end if
        if (min_norm_red < huge(1.0d0) .and. min_norm_red > 0.0d0) then
          cond_norm_red = max_norm_red / min_norm_red
        else
          min_norm_red = 0.0d0
          cond_norm_red = 0.0d0
        end if
        write(*,'(1x,a,10(a,1pe12.4),a,l1)') &
          '[DG-WPW-RED-BASIS-METRIC] total experimental:', &
          ' trace_S_WP_occ=', global_sum(1) * phys_scale, &
          ' trace_S_WP_unocc=', 0.0d0, &
          ' min_eig_S_WP=', min_norm_red, &
          ' max_eig_S_WP=', max_norm_red, &
          ' cond_S_WP=', cond_norm_red, &
          ' sum_occ_norm_loss=', sum_loss_red, &
          ' max_occ_norm_loss=', max_loss_red, &
          ' avg_occ_norm_loss=', avg_loss_red, &
          ' sum_self_norm_loss=', sum_loss_self, &
          ' max_self_norm_loss=', max_loss_self, &
          ' bad=', global_sum(6) > 0.5d0
        if (global_sum(26) > 0.0d0) then
          if (global_sum(39) > 0.0d0) then
            sdiff_rms = sqrt(max(0.0d0, global_sum(27) / global_sum(39)))
          else
            sdiff_rms = 0.0d0
          end if
          order_mismatch_score = max_out(6) / max(max_out(7), 1.0d-300)
          bsb_build_rms = sqrt(max(0.0d0, global_sum(37)))
          bsb_density_rms = sqrt(max(0.0d0, global_sum(38)))
          write(*,'(1x,a,12(a,1pe12.4),2(a,i0),a,l1)') &
            '[DG-WPW-RED-SRAW-CMP]', &
            ' nraw_build=', global_sum(33), &
            ' nraw_density=', global_sum(33), &
            ' max_abs_Sdiff=', max_out(6), &
            ' rms_Sdiff=', sdiff_rms, &
            ' trace_S_build=', global_sum(28), &
            ' trace_S_density=', global_sum(29), &
            ' norm_S_build=', sqrt(max(0.0d0, global_sum(30))), &
            ' norm_S_density=', sqrt(max(0.0d0, global_sum(31))), &
            ' order_mismatch_score=', order_mismatch_score, &
            ' BSB_build_frob_diff=', bsb_build_rms, &
            ' BSB_density_frob_diff=', bsb_density_rms, &
            ' ncmp=', global_sum(26), &
            ' rank_build=', nint(global_sum(32)), &
            ' rank_density=', nint(global_sum(21)), &
            ' available=', .true.
          write(*,'(1x,a,9(a,1pe12.4),a,l1)') &
            '[DG-WPW-RED-SRAW-DIAG-BLOCK]', &
            ' build_trace_W=', global_sum(40), &
            ' build_trace_P=', global_sum(41), &
            ' build_trace_neighbor=', global_sum(42), &
            ' density_trace_W=', global_sum(43), &
            ' density_trace_P=', global_sum(44), &
            ' density_trace_neighbor=', global_sum(45), &
            ' ratio_W=', global_sum(43) / max(global_sum(40), 1.0d-300), &
            ' ratio_P=', global_sum(44) / max(global_sum(41), 1.0d-300), &
            ' ratio_neighbor=', global_sum(45) / max(global_sum(42), 1.0d-300), &
            ' available=', .true.
          write(*,'(1x,a,4(a,1pe12.4),a,l1)') &
            '[DG-WPW-RED-SRAW-HYBRID-CMP]', &
            ' trace_hybrid=', global_sum(46), &
            ' max_abs_diff_vs_build=', sqrt(max(0.0d0, global_sum(47))), &
            ' rms_diff_vs_build=', sqrt(max(0.0d0, global_sum(47) / max(global_sum(39), 1.0d0))), &
            ' BSB_hybrid_frob_diff=', sqrt(max(0.0d0, global_sum(49))), &
            ' available=', .true.
        else
          write(*,'(1x,a,12(a,1pe12.4),2(a,i0),a,l1)') &
            '[DG-WPW-RED-SRAW-CMP]', &
            ' nraw_build=', 0.0d0, ' nraw_density=', 0.0d0, &
            ' max_abs_Sdiff=', 0.0d0, ' rms_Sdiff=', 0.0d0, &
            ' trace_S_build=', 0.0d0, ' trace_S_density=', 0.0d0, &
            ' norm_S_build=', 0.0d0, ' norm_S_density=', 0.0d0, &
            ' order_mismatch_score=', 0.0d0, &
            ' BSB_build_frob_diff=', 0.0d0, ' BSB_density_frob_diff=', 0.0d0, &
            ' ncmp=', 0.0d0, ' rank_build=', 0, ' rank_density=', 0, &
            ' available=', .false.
          write(*,'(1x,a,9(a,1pe12.4),a,l1)') &
            '[DG-WPW-RED-SRAW-DIAG-BLOCK]', &
            ' build_trace_W=', 0.0d0, ' build_trace_P=', 0.0d0, &
            ' build_trace_neighbor=', 0.0d0, ' density_trace_W=', 0.0d0, &
            ' density_trace_P=', 0.0d0, ' density_trace_neighbor=', 0.0d0, &
            ' ratio_W=', 0.0d0, ' ratio_P=', 0.0d0, &
            ' ratio_neighbor=', 0.0d0, ' available=', .false.
          write(*,'(1x,a,4(a,1pe12.4),a,l1)') &
            '[DG-WPW-RED-SRAW-HYBRID-CMP]', &
            ' trace_hybrid=', 0.0d0, ' max_abs_diff_vs_build=', 0.0d0, &
            ' rms_diff_vs_build=', 0.0d0, ' BSB_hybrid_frob_diff=', 0.0d0, &
            ' available=', .false.
        end if
        trace_w = 0.0d0
        trace_self = global_sum(2) * phys_scale
        trace_red = global_sum(1) * phys_scale
        trace_raw = 0.0d0
        max_loss_w_case = 0.0d0
        max_loss_self_case = 0.0d0
        max_loss_red_case = 0.0d0
        max_loss_raw_case = 0.0d0
        top_w_state = 0
        top_self_state = 0
        top_red_state = 0
        top_raw_state = 0
        do i_state = 1, nstate_metric
          occ_state = state_global(1, i_state) / dble(max(1, dg_frag%n_frag))
          if (occ_state <= 0.0d0) cycle
          trace_w = trace_w + state_global(6, i_state) * phys_scale
          trace_raw = trace_raw + state_global(9, i_state) * phys_scale
          loss_w = occ_state - state_global(6, i_state) * phys_scale
          loss_self_case = occ_state - state_global(2, i_state) * phys_scale
          loss_red_case = occ_state - state_global(3, i_state) * phys_scale
          loss_raw_case = occ_state - state_global(9, i_state) * phys_scale
          if (abs(loss_w) > max_loss_w_case) then
            max_loss_w_case = abs(loss_w)
            top_w_state = i_state
          end if
          if (abs(loss_self_case) > max_loss_self_case) then
            max_loss_self_case = abs(loss_self_case)
            top_self_state = i_state
          end if
          if (abs(loss_red_case) > max_loss_red_case) then
            max_loss_red_case = abs(loss_red_case)
            top_red_state = i_state
          end if
          if (abs(loss_raw_case) > max_loss_raw_case) then
            max_loss_raw_case = abs(loss_raw_case)
            top_raw_state = i_state
          end if
        end do
        write(*,'(1x,a,a,3(a,1pe12.4),a,i0,a,i0,a,l1)') &
          '[DG-WPW-RED-INIT-PROJ-COMPARE]', ' case=self_W_only', &
          ' trace_occ_norm=', trace_w, ' occ_norm_loss=', global_sum(15) / dble(max(1, dg_frag%n_frag)) - trace_w, &
          ' max_state_loss=', max_loss_w_case, ' top_loss_state=', top_w_state, &
          ' basis_count=', nint(global_sum(18)), ' available=', .true.
        write(*,'(1x,a,a,4(a,1pe12.4),2(a,i0),a,l1)') &
          '[DG-WPW-RED-INIT-PROJ-COMPARE]', ' case=self_WP', &
          ' trace_occ_norm=', trace_self, ' occ_norm_loss=', global_sum(15) / dble(max(1, dg_frag%n_frag)) - trace_self, &
          ' max_state_loss=', max_loss_self_case, ' cond_S=', cond_norm_red, &
          ' top_loss_state=', top_self_state, ' basis_count=', nint(global_sum(19)), ' available=', .true.
        if (global_sum(22) > 0.0d0) then
          write(*,'(1x,a,a,4(a,1pe12.4),2(a,i0),a,l1)') &
            '[DG-WPW-RED-INIT-PROJ-COMPARE]', ' case=self_WP_raw_neighbor', &
            ' trace_occ_norm=', trace_raw, ' occ_norm_loss=', global_sum(15) / dble(max(1, dg_frag%n_frag)) - trace_raw, &
            ' max_state_loss=', max_loss_raw_case, ' cond_S=', max_out(5), &
            ' top_loss_state=', top_raw_state, ' basis_count=', nint(global_sum(21)), ' available=', .true.
        else
          write(*,'(1x,a,a,3(a,1pe12.4),a,i0,a,i0,a,l1)') &
            '[DG-WPW-RED-INIT-PROJ-COMPARE]', ' case=self_WP_raw_neighbor', &
            ' trace_occ_norm=', 0.0d0, ' occ_norm_loss=', 0.0d0, ' max_state_loss=', 0.0d0, &
            ' top_loss_state=', 0, ' basis_count=', 0, ' available=', .false.
        end if
        if (global_sum(22) > 0.0d0) then
          write(*,'(1x,a,8(a,1pe12.4),2(a,i0),a,l1)') &
            '[DG-WPW-RED-RAW-TO-RED-PROJ]', &
            ' raw_charge_recovered=', global_sum(23) * phys_scale, &
            ' red_projected_charge=', global_sum(24) * phys_scale, &
            ' raw_to_red_loss=', (global_sum(23) - global_sum(24)) * phys_scale, &
            ' norm_raw_S=', global_sum(23) * phys_scale, &
            ' norm_residual_S=', global_sum(25) * phys_scale, &
            ' rel_residual_S=', global_sum(25) / max(global_sum(23), 1.0d-300), &
            ' sigma_min_kept=', 0.0d0, ' sigma_max=', 0.0d0, &
            ' rank_raw=', nint(global_sum(21)), ' rank_red=', nint(global_sum(20)), &
            ' available=', .true.
          if (global_sum(34) > 0.0d0) then
            write(*,'(1x,a,6(a,1pe12.4),2(a,i0),a,l1)') &
              '[DG-WPW-RED-RAW-TO-RED-PROJ-BUILD-S]', &
              ' raw_charge_recovered=', global_sum(34) * phys_scale, &
              ' red_projected_charge=', global_sum(35) * phys_scale, &
              ' raw_to_red_loss=', (global_sum(34) - global_sum(35)) * phys_scale, &
              ' norm_raw_S=', global_sum(34) * phys_scale, &
              ' norm_residual_S=', global_sum(36) * phys_scale, &
              ' rel_residual_S=', global_sum(36) / max(global_sum(34), 1.0d-300), &
              ' rank_raw=', nint(global_sum(32)), ' rank_red=', nint(global_sum(20)), &
              ' available=', .true.
            write(*,'(1x,a,6(a,1pe12.4),2(a,i0),a,l1)') &
              '[DG-WPW-RED-RAW-TO-RED-PROJ-HYBRID-S]', &
              ' raw_charge_recovered=', global_sum(50) * phys_scale, &
              ' red_projected_charge=', global_sum(51) * phys_scale, &
              ' raw_to_red_loss=', (global_sum(50) - global_sum(51)) * phys_scale, &
              ' norm_raw_S=', global_sum(50) * phys_scale, &
              ' norm_residual_S=', global_sum(52) * phys_scale, &
              ' rel_residual_S=', global_sum(52) / max(global_sum(50), 1.0d-300), &
              ' rank_raw=', nint(global_sum(32)), ' rank_red=', nint(global_sum(20)), &
              ' available=', .true.
          else
            write(*,'(1x,a,6(a,1pe12.4),2(a,i0),a,l1)') &
              '[DG-WPW-RED-RAW-TO-RED-PROJ-BUILD-S]', &
              ' raw_charge_recovered=', 0.0d0, ' red_projected_charge=', 0.0d0, &
              ' raw_to_red_loss=', 0.0d0, ' norm_raw_S=', 0.0d0, &
              ' norm_residual_S=', 0.0d0, ' rel_residual_S=', 0.0d0, &
              ' rank_raw=', 0, ' rank_red=', 0, ' available=', .false.
            write(*,'(1x,a,6(a,1pe12.4),2(a,i0),a,l1)') &
              '[DG-WPW-RED-RAW-TO-RED-PROJ-HYBRID-S]', &
              ' raw_charge_recovered=', 0.0d0, ' red_projected_charge=', 0.0d0, &
              ' raw_to_red_loss=', 0.0d0, ' norm_raw_S=', 0.0d0, &
              ' norm_residual_S=', 0.0d0, ' rel_residual_S=', 0.0d0, &
              ' rank_raw=', 0, ' rank_red=', 0, ' available=', .false.
          end if
        else
          write(*,'(1x,a,8(a,1pe12.4),2(a,i0),a,l1)') &
            '[DG-WPW-RED-RAW-TO-RED-PROJ]', &
            ' raw_charge_recovered=', 0.0d0, ' red_projected_charge=', 0.0d0, &
            ' raw_to_red_loss=', 0.0d0, ' norm_raw_S=', 0.0d0, &
            ' norm_residual_S=', 0.0d0, ' rel_residual_S=', 0.0d0, &
            ' sigma_min_kept=', 0.0d0, ' sigma_max=', 0.0d0, &
            ' rank_raw=', 0, ' rank_red=', 0, ' available=', .false.
          write(*,'(1x,a,6(a,1pe12.4),2(a,i0),a,l1)') &
            '[DG-WPW-RED-RAW-TO-RED-PROJ-BUILD-S]', &
            ' raw_charge_recovered=', 0.0d0, ' red_projected_charge=', 0.0d0, &
            ' raw_to_red_loss=', 0.0d0, ' norm_raw_S=', 0.0d0, &
            ' norm_residual_S=', 0.0d0, ' rel_residual_S=', 0.0d0, &
            ' rank_raw=', 0, ' rank_red=', 0, ' available=', .false.
          write(*,'(1x,a,6(a,1pe12.4),2(a,i0),a,l1)') &
            '[DG-WPW-RED-RAW-TO-RED-PROJ-HYBRID-S]', &
            ' raw_charge_recovered=', 0.0d0, ' red_projected_charge=', 0.0d0, &
            ' raw_to_red_loss=', 0.0d0, ' norm_raw_S=', 0.0d0, &
            ' norm_residual_S=', 0.0d0, ' rel_residual_S=', 0.0d0, &
            ' rank_raw=', 0, ' rank_red=', 0, ' available=', .false.
        end if
        write(*,'(1x,a,a,4(a,1pe12.4),2(a,i0),a,l1)') &
          '[DG-WPW-RED-INIT-PROJ-COMPARE]', ' case=self_WP_reduced_neighbor', &
          ' trace_occ_norm=', trace_red, ' occ_norm_loss=', global_sum(15) / dble(max(1, dg_frag%n_frag)) - trace_red, &
          ' max_state_loss=', max_loss_red_case, ' cond_S=', cond_norm_red, &
          ' top_loss_state=', top_red_state, ' basis_count=', nint(global_sum(20)), ' available=', .true.
        state_loss_work(:) = -1.0d300
        do i_state = 1, nstate_metric
          occ_state = state_global(1, i_state) / dble(max(1, dg_frag%n_frag))
          if (occ_state <= 0.0d0) cycle
          red_state_hvol = state_global(3, i_state) * phys_scale
          state_loss_work(i_state) = abs(occ_state - red_state_hvol)
        end do
        do top_iter = 1, min(5, nstate_metric)
          top_idx = maxloc(state_loss_work, dim=1)
          top_abs_loss = state_loss_work(top_idx)
          if (top_abs_loss < 0.0d0) exit
          occ_state = state_global(1, top_idx) / dble(max(1, dg_frag%n_frag))
          self_state_hvol = state_global(2, top_idx) * phys_scale
          red_state_hvol = state_global(3, top_idx) * phys_scale
          loss_state = occ_state - red_state_hvol
          norm_state = red_state_hvol / max(occ_state, 1.0d-300)
          write(*,'(1x,a,2(a,i0),7(a,1pe12.4))') &
            '[DG-WPW-RED-STATE-NORM] top_loss:', &
            ' rank=', top_iter, ' state=', top_idx, &
            ' occ=', occ_state, &
            ' norm_self_hvol=', self_state_hvol, &
            ' norm_red_hvol=', red_state_hvol, &
            ' loss_red=', loss_state, &
            ' abs_loss_red=', abs(loss_state), &
            ' norm_red_over_occ=', norm_state, &
            ' neighbor_cross_hvol=', (state_global(4, top_idx) + state_global(5, top_idx)) * phys_scale
          write(*,'(1x,a,2(a,i0),12(a,1pe12.4),2(a,l1))') &
            '[DG-WPW-RED-STATE-PROJ] top_loss:', &
            ' rank=', top_iter, ' state=', top_idx, &
            ' occ=', occ_state, &
            ' norm_prod=', occ_state, &
            ' norm_self_W=', state_global(6, top_idx) * phys_scale, &
            ' norm_self_P=', state_global(7, top_idx) * phys_scale, &
            ' norm_self_WP=', self_state_hvol, &
            ' norm_neighbor_raw=', 0.0d0, &
            ' norm_neighbor_reduced=', state_global(4, top_idx) * phys_scale, &
            ' norm_total_raw_WP=', 0.0d0, &
            ' norm_total_reduced_WP=', red_state_hvol, &
            ' loss_self_WP=', occ_state - self_state_hvol, &
            ' max_fragment_weight=', 0.0d0, &
            ' boundary_weight=', 0.0d0, &
            ' neighbor_raw_available=', .false., &
            ' fragment_weight_available=', .false.
          state_loss_work(top_idx) = -1.0d300
        end do
      end if
    end if
    if (do_heavy_diag) heavy_diag_done = .true.
    if (present(density_replaced)) density_replaced = &
      (replace_density_active .and. replace_density_zeroed .and. .not. bad_density .and. .not. bad_reproject) .or. &
      (canon_replace_density_active .and. canon_replace_density_zeroed .and. &
       global_sum(88) <= 0.5d0 .and. .not. bad_density)
    deallocate(state_local, state_global, state_loss_work)
  end subroutine diagnose_wpw_reduced_density

  subroutine apply_wpw_reduced_density_to_production(dg_frag, system, rho, rho_s, istep)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(inout) :: rho
    type(s_scalar), intent(inout) :: rho_s(system%nspin)
    integer, intent(in), optional :: istep
    character(len=32) :: env_density, env_current, env_source, env_reprojected, env_propagated, env_canon_density
    integer :: env_len, env_stat, step_use
    logical :: use_density, replaced, use_propagated_source, use_canon_density
    logical, save :: checked = .false.
    logical, save :: enabled = .false.
    logical, save :: source_propagated = .false.
    logical, save :: canonical_density_enabled = .false.
    logical, save :: warned_current = .false.
    logical, save :: warned_failed = .false.
    logical, save :: warned_propagated_invalid = .false.
    logical, save :: logged_enabled = .false.

    if (.not. checked) then
      env_source = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_DENSITY_SOURCE', &
        env_source, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_source(1:env_len))))
        case ('propagated','PROPAGATED','prop','PROP')
          enabled = .true.
          source_propagated = .true.
        case ('reprojected','REPROJECTED','reproject','REPROJECT')
          enabled = .true.
          source_propagated = .false.
        end select
      end if
      env_propagated = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_USE_PROPAGATED_DENSITY', &
        env_propagated, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_propagated(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
          source_propagated = .true.
        end select
      end if
      env_reprojected = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_USE_REPROJECTED_DENSITY', &
        env_reprojected, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_reprojected(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
          source_propagated = .false.
        end select
      end if
      env_density = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_USE_COEF_FOR_DENSITY', &
        env_density, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_density(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
          if (trim(env_source) == '' .and. trim(env_reprojected) == '' .and. trim(env_propagated) == '') &
            source_propagated = .false.
        end select
      end if
      env_canon_density = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_CANON_USE_DENSITY', &
        env_canon_density, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        select case (adjustl(trim(env_canon_density(1:env_len))))
        case ('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
          canonical_density_enabled = .true.
          source_propagated = .false.
        end select
      end if
      checked = .true.
    end if
    use_density = enabled
    use_canon_density = canonical_density_enabled
    use_propagated_source = source_propagated
    if (.not. use_density) return
    if (.not. logged_enabled .and. dg_frag%id == 0) then
      if (use_canon_density) then
        write(*,'(1x,a)') '[DG-WPW-RED-CANON-DENSITY-USE] enabled density-only current=off source=canonical'
      else if (use_propagated_source) then
        write(*,'(1x,a)') '[DG-WPW-RED-DIAG-DENSITY-USE] enabled source=propagated current=off'
      else
        write(*,'(1x,a)') '[DG-WPW-RED-DIAG-DENSITY-USE] enabled source=reprojected current=off'
      end if
      logged_enabled = .true.
    end if
    if (use_propagated_source .and. .not. warned_propagated_invalid .and. dg_frag%id == 0) then
      write(*,'(1x,a)') &
        '[DG-WPW-RED-DIAG-DENSITY-USE] WARNING: propagated WPW reduced density source is experimental and not valid for production density/current because auxiliary propagation does not match the production DG propagator.'
      write(*,'(1x,a)') &
        '[DG-WPW-RED-DIAG-DENSITY-USE] Use source=reprojected for the validated raw-state -> reduced density path.'
      warned_propagated_invalid = .true.
    end if
    if (use_propagated_source .and. .not. use_canon_density) return

    env_current = ''
    call get_environment_variable('SALMON_DG_WPW_REDUCED_USE_COEF_FOR_CURRENT', &
      env_current, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_current(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        if (.not. warned_current .and. dg_frag%id == 0) then
          write(*,'(1x,a)') &
            '[DG-WPW-RED-DIAG-DENSITY-USE] current replacement is not implemented; density-only replacement remains active.'
          warned_current = .true.
        end if
      end select
    end if

    step_use = -1
    if (present(istep)) step_use = istep
    replaced = .false.
    if (use_canon_density) then
      call diagnose_wpw_reduced_density(dg_frag, system, rho, step_use, -1, .false., 0.0d0, &
        rho_s, replace_density=.false., density_replaced=replaced, propagated_density_source=.false., &
        canonical_density_source=.true.)
    else
      call diagnose_wpw_reduced_density(dg_frag, system, rho, step_use, -1, .false., 0.0d0, &
        rho_s, .true., replaced, use_propagated_source)
    end if
    if (.not. replaced) then
      if (.not. warned_failed .and. dg_frag%id == 0) then
        if (use_canon_density) then
          write(*,'(1x,a)') &
            '[DG-WPW-RED-CANON-DENSITY-USE] requested canonical density replacement was unavailable or failed validation at this call.'
        else if (use_propagated_source) then
          write(*,'(1x,a)') &
            '[DG-WPW-RED-DIAG-DENSITY-USE] requested propagated density replacement was unavailable or failed validation at this call.'
        else
          write(*,'(1x,a)') &
            '[DG-WPW-RED-DIAG-DENSITY-USE] requested reprojected density replacement was unavailable or failed validation at this call.'
        end if
        warned_failed = .true.
      end if
    end if
  end subroutine apply_wpw_reduced_density_to_production

  subroutine apply_wpw_reduced_pz_to_production(dg_frag, system, istep)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    integer, intent(in), optional :: istep
    integer :: step_use
    real(8) :: pz_can_raw, pz_prod_raw, pz_can_density, pz_prod_density, pz_before
    real(8) :: pz_candidate, pz_output_delta, pz_replace_tol
    real(8) :: z_herm_diff, z_trace
    real(8) :: bridge_roundtrip_diff, bridge_t_norm, bridge_zcan_herm, bridge_pz_can
    real(8) :: bridge_pz_diff, bridge_rel_pz_diff
    real(8) :: block_zww_norm, block_zwp_norm, block_zpp_norm, block_zherm
    real(8) :: block_pz_can, block_pz_diff, block_rel_pz_diff
    real(8) :: canon_mixedz_pz_can, canon_mixedz_pz_diff, canon_mixedz_rel_pz_diff
    real(8) :: canon_mixedz_zherm
    real(8) :: pproj_norm_pself, pproj_norm_pneighbor, pproj_norm_pcan
    real(8) :: pproj_proj_pself, pproj_proj_pneighbor, pproj_proj_pcan
    real(8) :: pproj_leakage, pproj_rel_leakage, pproj_overlap_self_neighbor
    real(8) :: pproj_pz_can, pproj_pz_diff, pproj_rel_pz_diff
    real(8) :: pbasis_herm, pbasis_smin, pbasis_smax, pbasis_cond
    integer :: z_basis_dim, bridge_dim_can, bridge_dim_mixed, bridge_rank_est
    integer :: block_dim_w, block_dim_pself, block_dim_pneighbor
    integer :: pproj_dim_prod, pproj_dim_self, pproj_dim_neighbor, pproj_rank
    integer :: pbasis_dim_prod, pbasis_raw_before, pbasis_after_perp, pbasis_after_qmat
    logical :: replace_ok
    logical :: use_pz, pz_bad, mixedz_diag
    logical :: z_neighbor_terms, bridge_built, bridge_match, bridge_bad, block_bad
    logical :: canon_mixedz_built, canon_mixedz_match, canon_mixedz_bad
    logical :: pproj_built, pproj_bad
    logical :: pbasis_saved, ptransform_saved, pmetric_saved, pbasis_bad
    logical, save :: logged_enabled = .false.
    logical, save :: warned_unavailable = .false.
    logical, save :: interval_checked = .false.
    integer, save :: mixedz_diag_interval = 1
    integer :: env_len, env_stat
    character(len=32) :: env_interval
    character(len=32) :: mixed_z_h_route

    use_pz = wpw_reduced_canon_use_pz_enabled()
    mixedz_diag = wpw_reduced_canon_mixedz_pz_diag_enabled()
    if (.not. use_pz .and. .not. mixedz_diag) return
    if (.not. logged_enabled .and. dg_frag%id == 0) then
      if (use_pz) then
        write(*,'(1x,a)') '[DG-WPW-RED-CANON-PZ-USE] enabled z-only current=off source=canonical'
      else
        write(*,'(1x,a)') '[DG-WPW-RED-CANON-PZ-MIXEDZ-DIAG] enabled diagnostic-only'
      end if
      logged_enabled = .true.
    end if

    step_use = -1
    if (present(istep)) step_use = istep
    if (.not. interval_checked) then
      env_interval = ''
      call get_environment_variable('SALMON_DG_WPW_REDUCED_MIXEDZ_DIAG_INTERVAL', &
        env_interval, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_interval(1:env_len), *, iostat=env_stat) mixedz_diag_interval
        if (env_stat /= 0 .or. mixedz_diag_interval < 1) mixedz_diag_interval = 1
      end if
      interval_checked = .true.
    end if
    if (mixedz_diag .and. .not. use_pz) then
      if (step_use > 1 .and. mod(step_use, mixedz_diag_interval) /= 0) return
    end if
    pz_can_raw = 0.0d0
    pz_prod_raw = 0.0d0
    pz_bad = .true.
    call diagnose_wpw_reduced_density(dg_frag, system, istep=step_use, nstep=-1, coef_bad=.false., dt=0.0d0, &
      canonical_pz_source=.true., pz_can_raw=pz_can_raw, pz_prod_raw=pz_prod_raw, pz_bad=pz_bad)
    if (pz_bad) then
      if (.not. warned_unavailable .and. dg_frag%id == 0) then
        write(*,'(1x,a)') &
          '[DG-WPW-RED-CANON-PZ-USE] requested canonical Pz replacement was unavailable or failed validation at this call.'
        warned_unavailable = .true.
      end if
      return
    end if

    if (system%ngrid > 0 .and. system%hvol > 0.0d0) then
      pz_can_density = pz_can_raw / (dble(system%ngrid) * system%hvol)
      pz_prod_density = pz_prod_raw / (dble(system%ngrid) * system%hvol)
    else
      pz_can_density = 0.0d0
      pz_prod_density = 0.0d0
    end if
    pz_before = dg_frag%polarization_lg(3)
    pz_candidate = pz_can_density - dg_frag%polarization_lg_ref(3)
    pz_output_delta = pz_candidate - pz_before
    pz_replace_tol = 1.0d-8 * max(1.0d0, abs(pz_before))
    replace_ok = abs(pz_output_delta) <= pz_replace_tol
    canon_mixedz_pz_can = 0.0d0
    canon_mixedz_pz_diff = -dg_frag%dipole_lg_raw(3)
    canon_mixedz_rel_pz_diff = -1.0d0
    canon_mixedz_zherm = huge(1.0d0)
    canon_mixedz_built = .false.
    canon_mixedz_match = .false.
    canon_mixedz_bad = .true.
    if (dg_frag%mixed_wannier_bpw_final_uses_h_evec) then
      mixed_z_h_route = 'h_p_default'
    else
      mixed_z_h_route = 'legacy_h_vec'
    end if
    if (mixedz_diag) then
      call wpw_reduced_canon_mixedz_current_coeff_stats(dg_frag, system, dg_frag%dipole_lg_raw(3), &
        canon_mixedz_pz_can, canon_mixedz_pz_diff, canon_mixedz_rel_pz_diff, canon_mixedz_zherm, &
        canon_mixedz_built, canon_mixedz_match, canon_mixedz_bad)
    end if
    if (dg_frag%id == 0 .and. mixedz_diag) then
      call wpw_reduced_mixedz_operator_stats(dg_frag, z_herm_diff, z_trace, z_basis_dim, z_neighbor_terms)
      call wpw_reduced_canon_mixedz_bridge_stats(dg_frag, dg_frag%dipole_lg_raw(3), bridge_dim_can, &
        bridge_dim_mixed, bridge_t_norm, bridge_rank_est, bridge_roundtrip_diff, bridge_zcan_herm, &
        bridge_pz_can, bridge_pz_diff, bridge_rel_pz_diff, bridge_built, bridge_match, bridge_bad)
      call wpw_reduced_canon_pz_block_operator_stats(dg_frag, dg_frag%dipole_lg_raw(3), &
        block_dim_w, block_dim_pself, block_dim_pneighbor, block_zww_norm, block_zwp_norm, block_zpp_norm, &
        block_zherm, block_pz_can, block_pz_diff, block_rel_pz_diff, block_bad)
      call wpw_reduced_canon_p_projection_stats(dg_frag, dg_frag%dipole_lg_raw(3), &
        pproj_dim_prod, pproj_dim_self, pproj_dim_neighbor, pproj_rank, pproj_norm_pself, &
        pproj_norm_pneighbor, pproj_norm_pcan, pproj_proj_pself, pproj_proj_pneighbor, &
        pproj_proj_pcan, pproj_leakage, pproj_rel_leakage, pproj_overlap_self_neighbor, &
        pproj_pz_can, pproj_pz_diff, pproj_rel_pz_diff, pproj_built, pproj_bad)
      call wpw_reduced_prod_p_basis_save_stats(dg_frag, pbasis_dim_prod, pbasis_raw_before, &
        pbasis_after_perp, pbasis_after_qmat, pbasis_herm, pbasis_smin, pbasis_smax, pbasis_cond, &
        pbasis_saved, ptransform_saved, pmetric_saved, pbasis_bad)
      write(*,'(1x,a,4(a,i0),23(a,1pe12.4),a,l1)') &
        '[DG-WPW-RED-PROD-P-QMAT-METRIC-DIAG]', &
        ' step=', step_use, &
        ' dim_P_raw=', dg_frag%mixed_wannier_bpw_praw_dim, &
        ' dim_P_perp=', dg_frag%mixed_wannier_bpw_np, &
        ' dim_P_prod=', dg_frag%mixed_wannier_bpw_np, &
        ' S_raw_herm_diff=', dg_frag%mixed_wannier_bpw_sraw_herm_diff, &
        ' S_perp_herm_diff=', dg_frag%mixed_wannier_bpw_sperp_herm_diff, &
        ' qmat_metric_herm_diff=', dg_frag%mixed_wannier_bpw_qmat_metric_herm_diff, &
        ' qmat_metric_min_eval=', dg_frag%mixed_wannier_bpw_qmat_metric_min_eval, &
        ' qmat_metric_max_eval=', dg_frag%mixed_wannier_bpw_qmat_metric_max_eval, &
        ' qmat_metric_cond=', dg_frag%mixed_wannier_bpw_qmat_metric_cond, &
        ' qmat_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_qmat_metric_diff_from_i, &
        ' qmat_left_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_qleft_metric_diff_from_i, &
        ' final_metric_herm_diff=', dg_frag%mixed_wannier_bpw_final_metric_herm_diff, &
        ' final_metric_min_eval=', dg_frag%mixed_wannier_bpw_final_metric_min_eval, &
        ' final_metric_max_eval=', dg_frag%mixed_wannier_bpw_final_metric_max_eval, &
        ' final_metric_cond=', dg_frag%mixed_wannier_bpw_final_metric_cond, &
        ' final_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_final_metric_diff_from_i, &
        ' transform_metric_herm_diff=', dg_frag%mixed_wannier_bpw_transform_metric_herm_diff, &
        ' transform_metric_min_eval=', dg_frag%mixed_wannier_bpw_transform_metric_min_eval, &
        ' transform_metric_max_eval=', dg_frag%mixed_wannier_bpw_transform_metric_max_eval, &
        ' transform_metric_cond=', dg_frag%mixed_wannier_bpw_transform_metric_cond, &
        ' transform_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_transform_metric_diff_from_i, &
        ' transform_metric_diff_from_saved_metric=', dg_frag%mixed_wannier_bpw_transform_metric_diff_saved, &
        ' qmat_column_norm_min=', dg_frag%mixed_wannier_bpw_qmat_col_norm_min, &
        ' qmat_column_norm_max=', dg_frag%mixed_wannier_bpw_qmat_col_norm_max, &
        ' qmat_row_norm_min=', dg_frag%mixed_wannier_bpw_qmat_row_norm_min, &
        ' qmat_row_norm_max=', dg_frag%mixed_wannier_bpw_qmat_row_norm_max, &
        ' bad=', pbasis_bad
      write(*,'(1x,a,4(a,i0),8(a,1pe12.4),4(a,l1))') &
        '[DG-WPW-RED-CORRECTED-MIXEDZ-PZ-DIAG]', &
        ' step=', step_use, &
        ' dim_P_raw=', dg_frag%mixed_wannier_bpw_praw_dim, &
        ' dim_P_prod=', dg_frag%mixed_wannier_bpw_np, &
        ' basis_dim=', z_basis_dim, &
        ' Pz_prod_mixedZ=', dg_frag%dipole_lg_raw(3), &
        ' Pz_can_mixedZ=', canon_mixedz_pz_can, &
        ' Pz_diff=', canon_mixedz_pz_diff, &
        ' rel_Pz_diff=', canon_mixedz_rel_pz_diff, &
        ' qmat_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_qmat_metric_diff_from_i, &
        ' S_Pprod_cond=', pbasis_cond, &
        ' Zop_herm_diff=', z_herm_diff, &
        ' transform_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_transform_metric_diff_from_i, &
        ' Pprod_basis_saved=', pbasis_saved, &
        ' canonical_operator_built=', canon_mixedz_built, &
        ' convention_match=', canon_mixedz_match, &
        ' bad=', pbasis_bad .or. canon_mixedz_bad
      write(*,'(1x,a,1(a,i0),15(a,1pe12.4),2(a,l1),a,a,2(a,l1))') &
        '[DG-WPW-RED-PROD-P-FINAL-METRIC-DIAG]', &
        ' step=', step_use, &
        ' h_input_herm_diff=', dg_frag%mixed_wannier_bpw_h_input_herm_diff, &
        ' h_evec_unitarity_diff=', dg_frag%mixed_wannier_bpw_h_evec_unitarity_diff, &
        ' h_input_evec_diff=', dg_frag%mixed_wannier_bpw_h_input_evec_diff, &
        ' qmat_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_qmat_metric_diff_from_i, &
        ' final_metric_herm_diff=', dg_frag%mixed_wannier_bpw_final_metric_herm_diff, &
        ' final_metric_min_eval=', dg_frag%mixed_wannier_bpw_final_metric_min_eval, &
        ' final_metric_max_eval=', dg_frag%mixed_wannier_bpw_final_metric_max_eval, &
        ' final_metric_cond=', dg_frag%mixed_wannier_bpw_final_metric_cond, &
        ' final_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_final_metric_diff_from_i, &
        ' h_evec_metric_herm_diff=', dg_frag%mixed_wannier_bpw_transform_metric_herm_diff, &
        ' h_evec_metric_min_eval=', dg_frag%mixed_wannier_bpw_transform_metric_min_eval, &
        ' h_evec_metric_max_eval=', dg_frag%mixed_wannier_bpw_transform_metric_max_eval, &
        ' h_evec_metric_cond=', dg_frag%mixed_wannier_bpw_transform_metric_cond, &
        ' h_evec_metric_diff_from_I=', dg_frag%mixed_wannier_bpw_transform_metric_diff_from_i, &
        ' h_evec_metric_diff_from_saved=', dg_frag%mixed_wannier_bpw_transform_metric_diff_saved, &
        ' eigen_zheev_input_preserved=', .true., &
        ' eigen_zheev_evec_in_third_arg=', .true., &
        ' mixed_z_h_route=', trim(mixed_z_h_route), &
        ' legacy_route=', .not. dg_frag%mixed_wannier_bpw_final_uses_h_evec, &
        ' bad=', pbasis_bad
      write(*,'(1x,a,5(a,i0),4(a,1pe12.4),4(a,l1))') &
        '[DG-WPW-RED-PROD-P-BASIS-SAVE-DIAG]', &
        ' step=', step_use, &
        ' dim_P_prod=', pbasis_dim_prod, &
        ' dim_P_raw_before_perp=', pbasis_raw_before, &
        ' dim_P_after_perp=', pbasis_after_perp, &
        ' dim_P_after_qmat=', pbasis_after_qmat, &
        ' S_Pprod_herm_diff=', pbasis_herm, &
        ' S_Pprod_min_eval=', pbasis_smin, &
        ' S_Pprod_max_eval=', pbasis_smax, &
        ' S_Pprod_cond=', pbasis_cond, &
        ' Pprod_basis_saved=', pbasis_saved, &
        ' transform_saved=', ptransform_saved, &
        ' metric_saved=', pmetric_saved, &
        ' bad=', pbasis_bad
      write(*,'(1x,a,4(a,i0),7(a,1pe12.4),4(a,l1))') &
        '[DG-WPW-RED-CANON-MIXEDZ-BRIDGE-DIAG]', &
        ' step=', step_use, &
        ' T_dim_can=', bridge_dim_can, &
        ' T_dim_mixed=', bridge_dim_mixed, &
        ' T_rank_est=', bridge_rank_est, &
        ' T_norm=', bridge_t_norm, &
        ' bridge_roundtrip_diff=', bridge_roundtrip_diff, &
        ' Zcan_herm_diff=', bridge_zcan_herm, &
        ' Pz_prod_mixedZ=', dg_frag%dipole_lg_raw(3), &
        ' Pz_can_mixedZ=', bridge_pz_can, &
        ' Pz_diff=', bridge_pz_diff, &
        ' rel_Pz_diff=', bridge_rel_pz_diff, &
        ' canonical_operator_built=', bridge_built, &
        ' convention_match=', bridge_match, &
        ' Pz_replaced=', .false., &
        ' bad=', bridge_bad
      write(*,'(1x,a,4(a,i0),8(a,1pe12.4),4(a,l1))') &
        '[DG-WPW-RED-CANON-PZ-BLOCK-OP-DIAG]', &
        ' step=', step_use, &
        ' dim_W=', block_dim_w, &
        ' dim_P_self=', block_dim_pself, &
        ' dim_P_neighbor=', block_dim_pneighbor, &
        ' Z_WW_norm=', block_zww_norm, &
        ' Z_WP_norm=', block_zwp_norm, &
        ' Z_PP_norm=', block_zpp_norm, &
        ' Z_herm_diff=', block_zherm, &
        ' Pz_prod_mixedZ=', dg_frag%dipole_lg_raw(3), &
        ' Pz_can_block=', block_pz_can, &
        ' Pz_diff=', block_pz_diff, &
        ' rel_Pz_diff=', block_rel_pz_diff, &
        ' canonical_operator_built=', .false., &
        ' convention_match=', .false., &
        ' Pz_replaced=', .false., &
        ' bad=', block_bad
      write(*,'(1x,a,5(a,i0),13(a,1pe12.4),3(a,l1))') &
        '[DG-WPW-RED-CANON-P-PROJ-DIAG]', &
        ' step=', step_use, &
        ' dim_P_prod=', pproj_dim_prod, &
        ' dim_P_self=', pproj_dim_self, &
        ' dim_P_neighbor=', pproj_dim_neighbor, &
        ' rank_Pprod_Pcan=', pproj_rank, &
        ' norm_Pself=', pproj_norm_pself, &
        ' norm_Pneighbor=', pproj_norm_pneighbor, &
        ' norm_Pcan=', pproj_norm_pcan, &
        ' proj_norm_Pself_to_Pprod=', pproj_proj_pself, &
        ' proj_norm_Pneighbor_to_Pprod=', pproj_proj_pneighbor, &
        ' proj_norm_Pcan_to_Pprod=', pproj_proj_pcan, &
        ' leakage_norm_Pcan_from_Pprod=', pproj_leakage, &
        ' rel_leakage_Pcan_from_Pprod=', pproj_rel_leakage, &
        ' overlap_Pself_Pneighbor_norm=', pproj_overlap_self_neighbor, &
        ' Pz_prod_mixedZ=', dg_frag%dipole_lg_raw(3), &
        ' Pz_can_projectedP=', pproj_pz_can, &
        ' Pz_diff_projectedP=', pproj_pz_diff, &
        ' rel_Pz_diff_projectedP=', pproj_rel_pz_diff, &
        ' projector_built=', pproj_built, &
        ' Pz_replaced=', .false., &
        ' bad=', pproj_bad
      write(*,'(1x,a,1(a,i0),6(a,1pe12.4),1(a,i0),5(a,l1))') &
        '[DG-WPW-RED-CANON-MIXEDZ-OP-DIAG]', &
        ' step=', step_use, &
        ' Pz_prod_mixedZ=', dg_frag%dipole_lg_raw(3), &
        ' Pz_can_mixedZ=', canon_mixedz_pz_can, &
        ' Pz_diff=', canon_mixedz_pz_diff, &
        ' rel_Pz_diff=', canon_mixedz_rel_pz_diff, &
        ' Zop_herm_diff=', max(z_herm_diff, canon_mixedz_zherm), &
        ' Zop_trace=', z_trace, &
        ' basis_dim=', z_basis_dim, &
        ' neighbor_terms_included=', z_neighbor_terms, &
        ' canonical_operator_built=', canon_mixedz_built, &
        ' convention_match=', canon_mixedz_match, &
        ' Pz_replaced=', .false., &
        ' bad=', canon_mixedz_bad
      write(*,'(1x,a,1(a,i0),8(a,1pe12.4),4(a,l1),a)') &
        '[DG-WPW-RED-CANON-PZ-MIXEDZ-DIAG]', &
        ' step=', step_use, &
        ' Pz_mixedZ_raw=', dg_frag%dipole_lg_raw(3), &
        ' Pz_mixedZ_output=', pz_before, &
        ' Pz_grid_raw=', pz_can_raw, &
        ' Pz_grid_density=', pz_can_density, &
        ' Pz_grid_candidate=', pz_candidate, &
        ' mixedZ_minus_grid_raw=', dg_frag%dipole_lg_raw(3) - pz_can_raw, &
        ' output_delta_if_grid_used=', pz_output_delta, &
        ' rel_grid_raw_diff=', (pz_can_raw - pz_prod_raw) / max(abs(pz_prod_raw), 1.0d-300), &
        ' mixedz_operator_available=', .true., &
        ' canonical_mixedz_operator_available=', canon_mixedz_built, &
        ' convention_match=', canon_mixedz_match, &
        ' replacement_allowed=', replace_ok, &
        ' reason=canonical_mixedZ_available_grid_Pz_not_used_for_replacement'
    end if
    if (replace_ok) then
      dg_frag%dipole_lg_raw(3) = pz_can_raw
      dg_frag%polarization_lg(3) = pz_candidate
    end if
    if (dg_frag%id == 0) then
      if (.not. replace_ok) then
        write(*,'(1x,a)') &
          '[DG-WPW-RED-CANON-PZ-USE] canonical grid Pz is not in the production polarization convention; replacement skipped.'
      end if
      write(*,'(1x,a,1(a,i0),9(a,1pe12.4),3(a,l1))') &
        '[DG-WPW-RED-CANON-PZ-REPLACE]', &
        ' step=', step_use, &
        ' Pz_prod_before_raw=', pz_prod_raw, &
        ' Pz_can_raw=', pz_can_raw, &
        ' Pz_raw_diff=', pz_can_raw - pz_prod_raw, &
        ' rel_Pz_raw_diff=', (pz_can_raw - pz_prod_raw) / max(abs(pz_prod_raw), 1.0d-300), &
        ' Pz_prod_density=', pz_prod_density, &
        ' Pz_can_density=', pz_can_density, &
        ' Pz_output_before=', pz_before, &
        ' Pz_candidate=', pz_candidate, &
        ' dPz_output_delta=', pz_output_delta, &
        ' Pz_replaced=', replace_ok, &
        ' convention_match=', replace_ok, &
        ' current_replaced=', .false.
    end if
  end subroutine apply_wpw_reduced_pz_to_production








  subroutine diagnose_s_orthogonal_reduced_neighbor_variant(label, ifrag, H_self_ref, S_self_ref, H_in, S_in, n_self, n_ext)
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(in) :: ifrag, n_self, n_ext
    complex(8), intent(in) :: H_self_ref(n_self,n_self), S_self_ref(n_self,n_self)
    complex(8), intent(in) :: H_in(n_ext,n_ext), S_in(n_ext,n_ext)

    integer :: n_keep, n_drop, n_red, info
    real(8) :: tol, sn_after, snn_i_err, lambda_min, lambda_max, sss_min, sss_max
    complex(8), allocatable :: H_red(:,:), S_red(:,:)

    tol = wpw_neighbor_sorth_tol()
    call build_wpw_sorth_reduced_neighbor_block(H_in, S_in, n_self, n_ext, tol, H_red, S_red, &
      n_red, n_keep, n_drop, lambda_min, lambda_max, sn_after, snn_i_err, sss_min, sss_max, info)
    if (info /= 0) then
      write(*,'(1x,a,a,a,i0,a,1pe12.4,3(a,i0))') '[DG-WPW-MIXED-SORTH-REDUCED] ', trim(label), &
        ' ifrag=', ifrag, ' tol=', tol, ' n_keep=', n_keep, ' n_drop=', n_drop, ' info=', info
      if (allocated(H_red)) deallocate(H_red)
      if (allocated(S_red)) deallocate(S_red)
      return
    end if

    write(*,'(1x,a,a,a,i0,7(a,1pe12.4),4(a,i0))') &
      '[DG-WPW-MIXED-SORTH-REDUCED] ', trim(label), ' ifrag=', ifrag, &
      ' tol=', tol, ' lambda_keep_min=', lambda_min, ' lambda_keep_max=', lambda_max, &
      ' ||S_SN_after||=', sn_after, ' ||S_NN_reduced-I||=', snn_i_err, &
      ' SSS_min=', sss_min, ' SSS_max=', sss_max, &
      ' n_ext=', n_ext, ' n_red=', n_red, ' n_keep=', n_keep, ' n_drop=', n_drop

    call diagnose_self_like_eigen_tracking(label, ifrag, H_self_ref, S_self_ref, H_red, S_red, n_self, n_red)
    call diagnose_overlap_schur_effects(label, ifrag, H_self_ref, S_self_ref, H_red, S_red, n_self, n_red)

    deallocate(H_red, S_red)
  end subroutine diagnose_s_orthogonal_reduced_neighbor_variant













  integer function find_wpw_pp_block(blocks, ifrag, jfrag) result(iblk)
    implicit none
    type(complex_matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ifrag, jfrag
    integer :: i

    iblk = 0
    do i = 1, size(blocks)
      if (blocks(i)%ifrag_row == ifrag .and. blocks(i)%ifrag_col == jfrag) then
        iblk = i
        return
      end if
    end do
  end function find_wpw_pp_block

  subroutine compute_wpw_local_pp_interface(dg_frag, T_if)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    complex(8), intent(out) :: T_if(:,:)

    integer :: n_frag, n_pw, n_tot
    integer :: ifrag, jfrag, axis, side, ipw, jpw, idx, jdx
    integer :: ix, iy, iz, gx, gy, gz, grow(3), gcol(3)
    integer :: loop_lo(3), loop_hi(3)
    real(8) :: area_weight, alpha
    complex(8), allocatable :: local_if(:,:), global_if(:,:)
    complex(8) :: v_l, dnv_l, u_l, dnu_l, u_r, dnu_r, term_self, term_cross
    real(8) :: surface_penalty_factor

    n_frag = dg_frag%n_frag
    n_pw = dg_frag%n_plane_waves
    n_tot = n_frag * n_pw
    T_if(:, :) = (0.0d0, 0.0d0)
    if (n_frag <= 0 .or. n_pw <= 0) return
    surface_penalty_factor = wpw_interface_penalty_factor()

    allocate(local_if(n_tot, n_tot), global_if(n_tot, n_tot))
    local_if(:, :) = (0.0d0, 0.0d0)

    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      do axis = 1, 3
        if (dg_frag%num_fragment(axis) <= 1) cycle
        if (dg_frag%hgs(axis) <= 0.0d0) cycle
        area_weight = product(dg_frag%hgs(1:3)) / dg_frag%hgs(axis)
        alpha = surface_penalty_factor / dg_frag%hgs(axis)
        do side = -1, 1, 2
          jfrag = wpw_face_neighbor_fragment(dg_frag, ifrag, axis, side)
          if (jfrag <= 0 .or. jfrag == ifrag) cycle

          loop_lo(:) = dg_frag%ixyz_frag(:, ifrag)
          loop_hi(:) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
          if (side > 0) then
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag) + dg_frag%nxyz_domain(axis, ifrag) - 1
            loop_hi(axis) = loop_lo(axis)
          else
            loop_lo(axis) = dg_frag%ixyz_frag(axis, ifrag)
            loop_hi(axis) = loop_lo(axis)
          end if

          do iz = loop_lo(3), loop_hi(3)
            do iy = loop_lo(2), loop_hi(2)
              do ix = loop_lo(1), loop_hi(1)
                grow = [ix, iy, iz]
                gcol = grow
                if (side > 0) then
                  gcol(axis) = dg_frag%ixyz_frag(axis, jfrag)
                else
                  gcol(axis) = dg_frag%ixyz_frag(axis, jfrag) + dg_frag%nxyz_domain(axis, jfrag) - 1
                end if

                do jpw = 1, n_pw
                  jdx = (ifrag - 1) * n_pw + jpw
                  call wpw_basis_trace_value(dg_frag, ifrag, jpw, grow, side, axis, u_l, dnu_l)
                  do ipw = 1, n_pw
                    idx = (ifrag - 1) * n_pw + ipw
                    call wpw_basis_trace_value(dg_frag, ifrag, ipw, grow, side, axis, v_l, dnv_l)
                    term_self = (-0.25d0 * conjg(v_l) * dnu_l - 0.25d0 * conjg(dnv_l) * u_l + &
                      alpha * conjg(v_l) * u_l) * area_weight
                    local_if(idx, jdx) = local_if(idx, jdx) + term_self
                  end do
                end do

                do jpw = 1, n_pw
                  jdx = (jfrag - 1) * n_pw + jpw
                  call wpw_basis_trace_value(dg_frag, jfrag, jpw, gcol, side, axis, u_r, dnu_r)
                  do ipw = 1, n_pw
                    idx = (ifrag - 1) * n_pw + ipw
                    call wpw_basis_trace_value(dg_frag, ifrag, ipw, grow, side, axis, v_l, dnv_l)
                    term_cross = (-0.25d0 * conjg(v_l) * dnu_r + 0.25d0 * conjg(dnv_l) * u_r - &
                      alpha * conjg(v_l) * u_r) * area_weight
                    local_if(idx, jdx) = local_if(idx, jdx) + term_cross
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    call comm_summation(local_if, global_if, n_tot * n_tot, dg_frag%icomm)
    T_if(:, :) = global_if(:, :) / product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    deallocate(local_if, global_if)
  end subroutine compute_wpw_local_pp_interface

  subroutine wpw_basis_trace_value(dg_frag, ifrag, ipw, g, side, axis, val, dn)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, ipw, g(3), side, axis
    complex(8), intent(out) :: val, dn
    real(8) :: chi, grad_chi(3), phase_arg
    complex(8) :: phase, grad_phi(3)
    complex(8), parameter :: zi = (0.0d0, 1.0d0)

    call wpw_normalized_window_at_grid(dg_frag, ifrag, g(1), g(2), g(3), chi, grad_chi)
    phase_arg = dg_frag%k_pw(1, ipw) * dble(g(1)) * dg_frag%hgs(1) + &
                dg_frag%k_pw(2, ipw) * dble(g(2)) * dg_frag%hgs(2) + &
                dg_frag%k_pw(3, ipw) * dble(g(3)) * dg_frag%hgs(3)
    phase = exp(cmplx(0.0d0, phase_arg, kind=8))
    val = chi * phase
    grad_phi(1:3) = (cmplx(grad_chi(1:3), 0.0d0, kind=8) + zi * dg_frag%k_pw(1:3, ipw) * chi) * phase
    dn = dble(side) * grad_phi(axis)
  end subroutine wpw_basis_trace_value



  subroutine hermitize_matrix(mat, n)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(inout) :: mat(n, n)
    integer :: i, j

    do j = 1, n
      do i = j + 1, n
        mat(i, j) = 0.5d0 * (mat(i, j) + conjg(mat(j, i)))
        mat(j, i) = conjg(mat(i, j))
      end do
      mat(j, j) = cmplx(real(mat(j, j), 8), 0.0d0, kind=8)
    end do
  end subroutine hermitize_matrix


  subroutine wpw_local_block_magnitudes(dg_frag, mat, self_max, neighbor_max, nonneighbor_max)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    complex(8), intent(in) :: mat(:,:)
    real(8), intent(out) :: self_max, neighbor_max, nonneighbor_max
    integer :: n_pw, ifrag, jfrag, ipw, jpw, idx, jdx
    real(8) :: block_max

    n_pw = dg_frag%n_plane_waves
    self_max = 0.0d0
    neighbor_max = 0.0d0
    nonneighbor_max = 0.0d0
    if (n_pw <= 0) return

    do jfrag = 1, dg_frag%n_frag
      do ifrag = 1, dg_frag%n_frag
        block_max = 0.0d0
        do jpw = 1, n_pw
          jdx = (jfrag - 1) * n_pw + jpw
          do ipw = 1, n_pw
            idx = (ifrag - 1) * n_pw + ipw
            block_max = max(block_max, abs(mat(idx, jdx)))
          end do
        end do
        if (ifrag == jfrag) then
          self_max = max(self_max, block_max)
        else if (wpw_local_is_neighbor_pair(dg_frag, ifrag, jfrag)) then
          neighbor_max = max(neighbor_max, block_max)
        else
          nonneighbor_max = max(nonneighbor_max, block_max)
        end if
      end do
    end do
  end subroutine wpw_local_block_magnitudes



  !=======================================================================
  ! Compute mean potential energy for plane wave basis
  !=======================================================================
  subroutine compute_mean_potential(dg_frag, Vh, Vxc, Vpsl, V_mean)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in)  :: dg_frag
    type(s_scalar),         intent(in)  :: Vh, Vxc(:), Vpsl
    real(8),                intent(out) :: V_mean(:)  ! (nspin)

    integer :: ix, iy, iz, ispin
    real(8) :: vol_elem, total_volume, integral_local
    real(8), allocatable :: integral_all(:)

    allocate(integral_all(dg_frag%nspin))
    vol_elem = product(dg_frag%hgs(1:3))
    total_volume = product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))

    V_mean = 0.0d0
    integral_all = 0.0d0

    do ispin = 1, dg_frag%nspin
      integral_local = 0.0d0
!$omp parallel do collapse(3) reduction(+:integral_local) schedule(static)
      do iz = lbound(Vpsl%f, 3), ubound(Vpsl%f, 3)
        do iy = lbound(Vpsl%f, 2), ubound(Vpsl%f, 2)
          do ix = lbound(Vpsl%f, 1), ubound(Vpsl%f, 1)
            integral_local = integral_local + &
                 (Vpsl%f(ix,iy,iz) + Vh%f(ix,iy,iz) + Vxc(ispin)%f(ix,iy,iz)) * vol_elem
          end do
        end do
      end do
!$omp end parallel do
      integral_all(ispin) = integral_local
    end do

    call comm_summation(integral_all, V_mean, dg_frag%nspin, dg_frag%icomm)
    V_mean(1:dg_frag%nspin) = V_mean(1:dg_frag%nspin) / total_volume

    deallocate(integral_all)

  end subroutine compute_mean_potential

  !=======================================================================
  ! Compute fragment-plane wave Hamiltonian matrix elements
  !=======================================================================
  subroutine compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_complex)
    use structures
    use communication, only: comm_bcast, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_scalar),         intent(in)    :: Vh, Vxc(:), Vpsl
    complex(8),             intent(out)   :: H_complex(:,:,:)  ! (n_mat_max, n_pw, nspin)

    integer :: ipw, ifrag, i_local, io, ig, iloc, ispin, ix, iy, iz
    integer :: nx, ny, nz, nx_max, ny_max, nz_max
    integer :: gx, gy, gz, gx0, gy0, gz0, bx, by, bz, vx, vy, vz
    integer :: loc_s(3), loc_e(3)
    integer :: ipw_s, ipw_e
    integer :: p_lb1, p_lb2, p_lb3, p_ub1, p_ub2, p_ub3
    integer :: v_lb1, v_lb2, v_lb3, v_ub1, v_ub2, v_ub3
    integer :: env_len, env_stat
    real(8) :: k_vec(3), Lbox(3), sqrt_V, inv_sqrt_V
    real(8) :: vol_elem, k_squared, V_total, h_fp_norm
    complex(8) :: pw_val, pw_laplacian, hamiltonian_local, phase_yz
    complex(8) :: step_x, step_y, step_z
    complex(8) :: phase_x0, phase_y0, phase_z0
    complex(8), allocatable :: frag_block(:,:,:), frag_block_sum(:,:,:)
    complex(8), allocatable :: phase_x(:), phase_y(:), phase_z(:)
    real(8), allocatable :: V_box_cache(:,:,:,:)
    logical :: use_complex_basis
    logical :: owns_pw_col
    logical, allocatable :: V_box_ready(:)
    logical, save :: fp_domain_initialized = .false.
    logical, save :: use_buffered_domain = .false.
    character(len=64), save :: fp_domain_mode = 'core'
    character(len=64) :: env_fp_domain

    if (.not. dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return

    if (dg_frag%ifrag_start > dg_frag%ifrag_end) return

    vol_elem = product(dg_frag%hgs(1:3))
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))
    sqrt_V = sqrt(product(Lbox))
    inv_sqrt_V = 1.0d0 / sqrt_V
    use_complex_basis = allocated(dg_frag%phi_frag_c)
    p_lb1 = lbound(dg_frag%phi_frag, 1)
    p_lb2 = lbound(dg_frag%phi_frag, 2)
    p_lb3 = lbound(dg_frag%phi_frag, 3)
    p_ub1 = ubound(dg_frag%phi_frag, 1)
    p_ub2 = ubound(dg_frag%phi_frag, 2)
    p_ub3 = ubound(dg_frag%phi_frag, 3)
    v_lb1 = lbound(Vpsl%f, 1)
    v_lb2 = lbound(Vpsl%f, 2)
    v_lb3 = lbound(Vpsl%f, 3)
    v_ub1 = ubound(Vpsl%f, 1)
    v_ub2 = ubound(Vpsl%f, 2)
    v_ub3 = ubound(Vpsl%f, 3)

    if (.not. fp_domain_initialized) then
      env_fp_domain = ''
      call get_environment_variable('SALMON_DG_FP_INTEGRAL_DOMAIN', env_fp_domain, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        fp_domain_mode = adjustl(env_fp_domain(1:env_len))
      else
        ! Match the real-space density materialization domain by default.
        ! Buffered FP integrals can still be requested explicitly for tests.
        fp_domain_mode = 'core'
      end if
      select case (fp_domain_mode(1:1))
      case ('b','B','p','P','f','F','1')
        use_buffered_domain = .true.
      case default
        use_buffered_domain = .false.
      end select
      fp_domain_initialized = .true.
      if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
        write(*,'(1x,a,a)') '[FP-DOMAIN] H_fp integral domain mode = ', trim(fp_domain_mode)
      end if
    end if

    if (use_buffered_domain .and. allocated(dg_frag%frag_buf_lo) .and. allocated(dg_frag%frag_buf_hi)) then
      nx_max = max(1, maxval(dg_frag%frag_buf_hi(1, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(1, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
      ny_max = max(1, maxval(dg_frag%frag_buf_hi(2, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(2, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
      nz_max = max(1, maxval(dg_frag%frag_buf_hi(3, dg_frag%ifrag_start:dg_frag%ifrag_end) - &
                           dg_frag%frag_buf_lo(3, dg_frag%ifrag_start:dg_frag%ifrag_end) + 1))
    else
      if (use_buffered_domain) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a)') '[FP-DOMAIN] buffered FP domain requested but fragment buffer bounds are unavailable'
        end if
        stop "DG-Fragment RT: missing fragment buffer bounds for buffered FP domain"
      end if
      nx_max = maxval(dg_frag%nxyz_domain(1, dg_frag%ifrag_start:dg_frag%ifrag_end))
      ny_max = maxval(dg_frag%nxyz_domain(2, dg_frag%ifrag_start:dg_frag%ifrag_end))
      nz_max = maxval(dg_frag%nxyz_domain(3, dg_frag%ifrag_start:dg_frag%ifrag_end))
    end if
    allocate(phase_x(nx_max), phase_y(ny_max), phase_z(nz_max))
    if (dg_frag%parallel_mode_orbital) then
      allocate(V_box_cache(nx_max, ny_max, nz_max, max(1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1)))
      allocate(V_box_ready(max(1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1)))
      V_box_cache(:, :, :, :) = 0.0d0
      V_box_ready(:) = .false.
    end if

    H_complex = (0.0d0, 0.0d0)
    call get_fragment_pw_column_range(dg_frag, dg_frag%n_plane_waves, ipw_s, ipw_e)

    do ispin = 1, dg_frag%nspin
      if (allocated(V_box_ready)) V_box_ready(:) = .false.
      do ipw = 1, dg_frag%n_plane_waves
        owns_pw_col = (.not. dg_frag%parallel_mode_orbital) .or. (ipw >= ipw_s .and. ipw <= ipw_e)
        if (.not. owns_pw_col) cycle
        k_vec(1:3) = dg_frag%k_pw(1:3, ipw)
        k_squared = sum(k_vec**2)

        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          if (use_buffered_domain .and. allocated(dg_frag%frag_buf_lo) .and. allocated(dg_frag%frag_buf_hi)) then
            gx0 = dg_frag%frag_buf_lo(1, ifrag)
            gy0 = dg_frag%frag_buf_lo(2, ifrag)
            gz0 = dg_frag%frag_buf_lo(3, ifrag)
            nx = dg_frag%frag_buf_hi(1, ifrag) - dg_frag%frag_buf_lo(1, ifrag) + 1
            ny = dg_frag%frag_buf_hi(2, ifrag) - dg_frag%frag_buf_lo(2, ifrag) + 1
            nz = dg_frag%frag_buf_hi(3, ifrag) - dg_frag%frag_buf_lo(3, ifrag) + 1
          else
            nx = dg_frag%nxyz_domain(1, ifrag)
            ny = dg_frag%nxyz_domain(2, ifrag)
            nz = dg_frag%nxyz_domain(3, ifrag)
            gx0 = dg_frag%ixyz_frag(1, ifrag)
            gy0 = dg_frag%ixyz_frag(2, ifrag)
            gz0 = dg_frag%ixyz_frag(3, ifrag)
          end if

          step_x = exp(cmplx(0.0d0, k_vec(1) * dg_frag%hgs(1), kind=8))
          step_y = exp(cmplx(0.0d0, k_vec(2) * dg_frag%hgs(2), kind=8))
          step_z = exp(cmplx(0.0d0, k_vec(3) * dg_frag%hgs(3), kind=8))
          phase_x0 = exp(cmplx(0.0d0, k_vec(1) * dble(gx0) * dg_frag%hgs(1), kind=8))
          phase_y0 = exp(cmplx(0.0d0, k_vec(2) * dble(gy0) * dg_frag%hgs(2), kind=8))
          phase_z0 = exp(cmplx(0.0d0, k_vec(3) * dble(gz0) * dg_frag%hgs(3), kind=8))
          phase_x(1) = phase_x0
          do ix = 2, nx
            phase_x(ix) = phase_x(ix-1) * step_x
          end do
          phase_y(1) = phase_y0
          do iy = 2, ny
            phase_y(iy) = phase_y(iy-1) * step_y
          end do
          phase_z(1) = phase_z0
          do iz = 2, nz
            phase_z(iz) = phase_z(iz-1) * step_z
          end do

          if (dg_frag%parallel_mode_orbital) then
            ! H_fp uses the same PW-column ownership as S_fp.  The scalar
            ! potential is independent of the PW column.  Assemble it once per
            ! fragment/spin and reuse it for all PW columns owned by this rank.
            loc_s(:) = [1, 1, 1]
            loc_e(:) = [nx, ny, nz]
            if (.not. V_box_ready(i_local)) then
              call build_fragment_pw_total_potential_box(dg_frag, ifrag, Vh, Vxc(ispin), Vpsl, gx0, gy0, gz0, &
                V_box_cache(1:nx, 1:ny, 1:nz, i_local))
              V_box_ready(i_local) = .true.
            end if
          else
            call get_fragment_subgroup_box_range_pw(dg_frag, [nx, ny, nz], loc_s, loc_e)
          end if
          if (any(loc_s(:) > loc_e(:))) cycle

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig = dg_frag%index_basis(io, ifrag, ispin)
            iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(H_complex, 1))
            if (iloc <= 0) cycle
            hamiltonian_local = (0.0d0, 0.0d0)

            if (use_complex_basis) then
              do iz = loc_s(3), loc_e(3)
                gz = gz0 + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = loc_s(2), loc_e(2)
                  gy = gy0 + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  phase_yz = phase_y(iy) * phase_z(iz)
                  do ix = loc_s(1), loc_e(1)
                    gx = gx0 + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    pw_val = phase_x(ix) * phase_yz * inv_sqrt_V
                    pw_laplacian = (k_squared / 2.0d0) * pw_val
                    if (dg_frag%parallel_mode_orbital) then
                      V_total = V_box_cache(ix, iy, iz, i_local)
                    else
                      vx = map_global_to_phi_box_coord_pw(gx, v_lb1, v_ub1, dg_frag%lgnum_total(1))
                      if (vx < v_lb1 .or. vx > v_ub1) cycle
                      vy = map_global_to_phi_box_coord_pw(gy, v_lb2, v_ub2, dg_frag%lgnum_total(2))
                      if (vy < v_lb2 .or. vy > v_ub2) cycle
                      vz = map_global_to_phi_box_coord_pw(gz, v_lb3, v_ub3, dg_frag%lgnum_total(3))
                      if (vz < v_lb3 .or. vz > v_ub3) cycle
                      V_total = Vpsl%f(vx, vy, vz) + Vh%f(vx, vy, vz) + Vxc(ispin)%f(vx, vy, vz)
                    end if

                    hamiltonian_local = hamiltonian_local + &
                         conjg(dg_frag%phi_frag_c(bx,by,bz,io,i_local)) * &
                         (pw_laplacian + V_total * pw_val) * vol_elem
                  end do
                end do
              end do
            else
              ! Keep H_fp in the same PW phase convention as S_fp and density.
              do iz = loc_s(3), loc_e(3)
                gz = gz0 + iz - 1
                bz = map_global_to_phi_box_coord_pw(gz, p_lb3, p_ub3, dg_frag%lgnum_total(3))
                if (bz < p_lb3 .or. bz > p_ub3) cycle
                do iy = loc_s(2), loc_e(2)
                  gy = gy0 + iy - 1
                  by = map_global_to_phi_box_coord_pw(gy, p_lb2, p_ub2, dg_frag%lgnum_total(2))
                  if (by < p_lb2 .or. by > p_ub2) cycle
                  phase_yz = phase_y(iy) * phase_z(iz)
                  do ix = loc_s(1), loc_e(1)
                    gx = gx0 + ix - 1
                    bx = map_global_to_phi_box_coord_pw(gx, p_lb1, p_ub1, dg_frag%lgnum_total(1))
                    if (bx < p_lb1 .or. bx > p_ub1) cycle
                    pw_val = phase_x(ix) * phase_yz * inv_sqrt_V
                    pw_laplacian = (k_squared / 2.0d0) * pw_val
                    if (dg_frag%parallel_mode_orbital) then
                      V_total = V_box_cache(ix, iy, iz, i_local)
                    else
                      vx = map_global_to_phi_box_coord_pw(gx, v_lb1, v_ub1, dg_frag%lgnum_total(1))
                      if (vx < v_lb1 .or. vx > v_ub1) cycle
                      vy = map_global_to_phi_box_coord_pw(gy, v_lb2, v_ub2, dg_frag%lgnum_total(2))
                      if (vy < v_lb2 .or. vy > v_ub2) cycle
                      vz = map_global_to_phi_box_coord_pw(gz, v_lb3, v_ub3, dg_frag%lgnum_total(3))
                      if (vz < v_lb3 .or. vz > v_ub3) cycle
                      V_total = Vpsl%f(vx, vy, vz) + Vh%f(vx, vy, vz) + Vxc(ispin)%f(vx, vy, vz)
                    end if

                    hamiltonian_local = hamiltonian_local + &
                      dg_frag%phi_frag(bx,by,bz,io,i_local) * &
                      (pw_laplacian + V_total * pw_val) * vol_elem
                  end do
                end do
              end do
            end if

            H_complex(iloc, ipw, ispin) = H_complex(iloc, ipw, ispin) + hamiltonian_local
          end do
        end do
      end do
    end do

    if (allocated(phase_x)) deallocate(phase_x)
    if (allocated(phase_y)) deallocate(phase_y)
    if (allocated(phase_z)) deallocate(phase_z)
    if (allocated(V_box_cache)) deallocate(V_box_cache)
    if (allocated(V_box_ready)) deallocate(V_box_ready)

    allocate(frag_block(dg_frag%nstate_frag, dg_frag%n_plane_waves, dg_frag%nspin))
    allocate(frag_block_sum(dg_frag%nstate_frag, dg_frag%n_plane_waves, dg_frag%nspin))
    do ifrag = 1, dg_frag%n_frag
      frag_block(:, :, :) = (0.0d0, 0.0d0)
      if (ifrag >= dg_frag%ifrag_start .and. ifrag <= dg_frag%ifrag_end) then
        i_local = ifrag - dg_frag%ifrag_start + 1
        do ispin = 1, dg_frag%nspin
          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig = dg_frag%index_basis(io, ifrag, ispin)
            iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(H_complex, 1))
            if (iloc <= 0) cycle
            frag_block(io, 1:dg_frag%n_plane_waves, ispin) = H_complex(iloc, 1:dg_frag%n_plane_waves, ispin)
          end do
        end do
        ! H_fp follows the same ownership as S_fp: real-space mode reduces spatial
        ! pieces, orbital mode reduces disjoint PW-column pieces.
        if (dg_frag%isize_frag > 1) then
          call comm_summation(frag_block, frag_block_sum, size(frag_block), dg_frag%icomm_frag)
          frag_block(:, :, :) = frag_block_sum(:, :, :)
        end if
      end if
      call comm_bcast(frag_block, dg_frag%icomm, dg_frag%id_array(ifrag))
      do ispin = 1, dg_frag%nspin
        do io = 1, dg_frag%n_basis(ifrag, ispin)
          ig = dg_frag%index_basis(io, ifrag, ispin)
          iloc = fragment_local_row_index_pw(dg_frag, ig, ispin, size(H_complex, 1))
          if (iloc <= 0) cycle
          H_complex(iloc, 1:dg_frag%n_plane_waves, ispin) = frag_block(io, 1:dg_frag%n_plane_waves, ispin)
        end do
      end do
    end do
    deallocate(frag_block, frag_block_sum)

    h_fp_norm = sqrt(sum(abs(H_complex(1:size(H_complex, 1), 1:dg_frag%n_plane_waves, 1:dg_frag%nspin))**2))
    if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
      write(*,'(1x,a,1x,1pe14.6)') '[FP-DOMAIN] ||H_fp||_F =', h_fp_norm
    end if

  end subroutine compute_fragment_pw_hamiltonian

  subroutine compute_fragment_pw_gradient_from_overlap(dg_frag, S_frag_pw, P_frag_pw)
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(in) :: S_frag_pw(:,:,:)
    complex(8), intent(out) :: P_frag_pw(:,:,:,:)

    integer :: ispin, ipw, idir
    integer :: env_len, env_stat
    character(len=64) :: env_mfp_mode
    logical :: use_fft_mfp
    complex(8), parameter :: zi = (0.0d0, 1.0d0)

    P_frag_pw(:, :, :, :) = (0.0d0, 0.0d0)
    if (.not. dg_frag%use_plane_wave_basis) return
    if (dg_frag%n_plane_waves <= 0) return

    env_mfp_mode = ''
    call get_environment_variable('SALMON_DG_MFP_MODE', env_mfp_mode, length=env_len, status=env_stat)
    use_fft_mfp = .false.
    if (env_stat == 0 .and. env_len > 0) then
      if (env_mfp_mode(1:1) == 'f' .or. env_mfp_mode(1:1) == 'F' .or. env_mfp_mode(1:1) == 'g' .or. env_mfp_mode(1:1) == 'G') then
        use_fft_mfp = .true.
      end if
    end if

    do ispin = 1, dg_frag%nspin
      do ipw = 1, dg_frag%n_plane_waves
        do idir = 1, 3
          P_frag_pw(idir, :, ipw, ispin) = zi * dg_frag%k_pw(idir, ipw) * S_frag_pw(:, ipw, ispin)
        end do
      end do
    end do

    if (use_fft_mfp) then
      if (.not. allocated(dg_frag%P_mat_frag_pw_g)) then
        allocate(dg_frag%P_mat_frag_pw_g(size(P_frag_pw,1), size(P_frag_pw,2), size(P_frag_pw,3), size(P_frag_pw,4)))
      else if (size(dg_frag%P_mat_frag_pw_g,1) /= size(P_frag_pw,1) .or. size(dg_frag%P_mat_frag_pw_g,2) /= size(P_frag_pw,2) .or. &
               size(dg_frag%P_mat_frag_pw_g,3) /= size(P_frag_pw,3) .or. size(dg_frag%P_mat_frag_pw_g,4) /= size(P_frag_pw,4)) then
        deallocate(dg_frag%P_mat_frag_pw_g)
        allocate(dg_frag%P_mat_frag_pw_g(size(P_frag_pw,1), size(P_frag_pw,2), size(P_frag_pw,3), size(P_frag_pw,4)))
      end if
      dg_frag%P_mat_frag_pw_g(:, :, :, :) = P_frag_pw(:, :, :, :)
    end if
  end subroutine compute_fragment_pw_gradient_from_overlap

  subroutine build_local_fragment_pw_row_list(dg_frag, ispin)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin

    integer :: i, nrow, nvalid, nsrc, ispin_eff, ispin_loop, nspin_src
    integer, allocatable :: row_ids(:)
    logical, allocatable :: row_seen(:)

    if (allocated(dg_frag%fp_local_row_ids)) deallocate(dg_frag%fp_local_row_ids)

    ispin_eff = ispin
    nvalid = 0

    if (ispin_eff <= 0) then
      if (allocated(dg_frag%local_coef_global_ids)) then
        nspin_src = min(max(1, dg_frag%nspin), size(dg_frag%local_coef_global_ids, 2))
        allocate(row_ids(max(0, dg_frag%n_mat_max)))
        allocate(row_seen(max(0, dg_frag%n_mat_max)))
        row_seen(:) = .false.
        do ispin_loop = 1, nspin_src
          nsrc = size(dg_frag%local_coef_global_ids, 1)
          if (allocated(dg_frag%local_coef_count)) then
            if (ispin_loop <= size(dg_frag%local_coef_count)) then
              nsrc = min(nsrc, max(0, dg_frag%local_coef_count(ispin_loop)))
            end if
          end if
          do i = 1, nsrc
            nrow = dg_frag%local_coef_global_ids(i, ispin_loop)
            if (nrow < 1 .or. nrow > dg_frag%n_mat_max) cycle
            if (allocated(dg_frag%coef_owner)) then
              if (nrow <= size(dg_frag%coef_owner, 1) .and. ispin_loop <= size(dg_frag%coef_owner, 2)) then
                if (dg_frag%coef_owner(nrow, ispin_loop) /= dg_frag%id) cycle
              end if
            end if
            if (row_seen(nrow)) cycle
            nvalid = nvalid + 1
            row_ids(nvalid) = nrow
            row_seen(nrow) = .true.
          end do
        end do
        allocate(dg_frag%fp_local_row_ids(nvalid))
        if (nvalid > 0) dg_frag%fp_local_row_ids(1:nvalid) = row_ids(1:nvalid)
        deallocate(row_ids, row_seen)
        return
      end if

      if (allocated(dg_frag%coef_owner)) then
        nspin_src = min(max(1, dg_frag%nspin), size(dg_frag%coef_owner, 2))
        nsrc = min(dg_frag%n_mat_max, size(dg_frag%coef_owner, 1))
        allocate(row_ids(max(0, nsrc)))
        allocate(row_seen(max(0, dg_frag%n_mat_max)))
        row_seen(:) = .false.
        do ispin_loop = 1, nspin_src
          do i = 1, nsrc
            if (dg_frag%coef_owner(i, ispin_loop) /= dg_frag%id) cycle
            if (row_seen(i)) cycle
            nvalid = nvalid + 1
            row_ids(nvalid) = i
            row_seen(i) = .true.
          end do
        end do
        allocate(dg_frag%fp_local_row_ids(nvalid))
        if (nvalid > 0) dg_frag%fp_local_row_ids(1:nvalid) = row_ids(1:nvalid)
        deallocate(row_ids, row_seen)
        return
      end if

      nsrc = 0
      if (allocated(dg_frag%coef)) nsrc = min(size(dg_frag%coef, 1), dg_frag%n_mat_max)
      allocate(dg_frag%fp_local_row_ids(nsrc))
      do i = 1, nsrc
        dg_frag%fp_local_row_ids(i) = i
      end do
      return
    end if

    ispin_eff = max(1, min(ispin_eff, max(1, dg_frag%nspin)))

    if (allocated(dg_frag%local_coef_global_ids)) then
      if (ispin_eff <= size(dg_frag%local_coef_global_ids, 2)) then
        nsrc = size(dg_frag%local_coef_global_ids, 1)
        if (allocated(dg_frag%local_coef_count)) then
          if (ispin_eff <= size(dg_frag%local_coef_count)) then
            nsrc = min(nsrc, max(0, dg_frag%local_coef_count(ispin_eff)))
          end if
        end if
        allocate(row_ids(max(0, nsrc)))
        do i = 1, nsrc
          nrow = dg_frag%local_coef_global_ids(i, ispin_eff)
          if (nrow < 1 .or. nrow > dg_frag%n_mat_max) cycle
          if (allocated(dg_frag%coef_owner)) then
            if (nrow <= size(dg_frag%coef_owner, 1) .and. ispin_eff <= size(dg_frag%coef_owner, 2)) then
              if (dg_frag%coef_owner(nrow, ispin_eff) /= dg_frag%id) cycle
            end if
          end if
          nvalid = nvalid + 1
          row_ids(nvalid) = nrow
        end do
        allocate(dg_frag%fp_local_row_ids(nvalid))
        if (nvalid > 0) dg_frag%fp_local_row_ids(1:nvalid) = row_ids(1:nvalid)
        deallocate(row_ids)
        return
      end if
    end if

    if (allocated(dg_frag%coef_owner)) then
      if (ispin_eff <= size(dg_frag%coef_owner, 2)) then
        nsrc = min(dg_frag%n_mat_max, size(dg_frag%coef_owner, 1))
        allocate(row_ids(nsrc))
        do i = 1, nsrc
          if (dg_frag%coef_owner(i, ispin_eff) /= dg_frag%id) cycle
          nvalid = nvalid + 1
          row_ids(nvalid) = i
        end do
        allocate(dg_frag%fp_local_row_ids(nvalid))
        if (nvalid > 0) dg_frag%fp_local_row_ids(1:nvalid) = row_ids(1:nvalid)
        deallocate(row_ids)
        return
      end if
    end if

    nsrc = 0
    if (allocated(dg_frag%coef)) nsrc = min(size(dg_frag%coef, 1), dg_frag%n_mat_max)
    allocate(dg_frag%fp_local_row_ids(nsrc))
    do i = 1, nsrc
      dg_frag%fp_local_row_ids(i) = i
    end do
  end subroutine build_local_fragment_pw_row_list

  subroutine prepare_local_fragment_pw_blocks(dg_frag, ispin)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in), optional :: ispin

    integer :: ispin_eff, nrow_local, npw_local, n_pw
    integer :: idir, irow, ipw, ispin_loop, global_row, global_pw, source_row
    logical :: have_S_full, have_H_full, have_P_full
    complex(8), parameter :: zi = (0.0d0, 1.0d0)

    ispin_eff = 0
    if (present(ispin)) ispin_eff = ispin
    call build_local_fragment_pw_row_list(dg_frag, ispin_eff)

    n_pw = max(0, dg_frag%n_plane_waves)
    if (allocated(dg_frag%fp_local_pw_ids)) deallocate(dg_frag%fp_local_pw_ids)
    allocate(dg_frag%fp_local_pw_ids(n_pw))
    ! Transitional Task-2 layout: keep every PW row requested until Task 3
    ! adds a distributed PW-row list for Y_P.
    do ipw = 1, n_pw
      dg_frag%fp_local_pw_ids(ipw) = ipw
    end do

    nrow_local = size(dg_frag%fp_local_row_ids)
    npw_local = size(dg_frag%fp_local_pw_ids)

    if (.not. allocated(dg_frag%S_mat_frag_pw_local)) then
      allocate(dg_frag%S_mat_frag_pw_local(nrow_local, npw_local, dg_frag%nspin))
    else if (size(dg_frag%S_mat_frag_pw_local, 1) /= nrow_local .or. &
             size(dg_frag%S_mat_frag_pw_local, 2) /= npw_local .or. &
             size(dg_frag%S_mat_frag_pw_local, 3) /= dg_frag%nspin) then
      deallocate(dg_frag%S_mat_frag_pw_local)
      allocate(dg_frag%S_mat_frag_pw_local(nrow_local, npw_local, dg_frag%nspin))
    end if

    if (.not. allocated(dg_frag%H_mat_frag_pw_local)) then
      allocate(dg_frag%H_mat_frag_pw_local(nrow_local, npw_local, dg_frag%nspin))
    else if (size(dg_frag%H_mat_frag_pw_local, 1) /= nrow_local .or. &
             size(dg_frag%H_mat_frag_pw_local, 2) /= npw_local .or. &
             size(dg_frag%H_mat_frag_pw_local, 3) /= dg_frag%nspin) then
      deallocate(dg_frag%H_mat_frag_pw_local)
      allocate(dg_frag%H_mat_frag_pw_local(nrow_local, npw_local, dg_frag%nspin))
    end if

    if (.not. allocated(dg_frag%P_mat_frag_pw_local)) then
      allocate(dg_frag%P_mat_frag_pw_local(3, nrow_local, npw_local, dg_frag%nspin))
    else if (size(dg_frag%P_mat_frag_pw_local, 1) /= 3 .or. &
             size(dg_frag%P_mat_frag_pw_local, 2) /= nrow_local .or. &
             size(dg_frag%P_mat_frag_pw_local, 3) /= npw_local .or. &
             size(dg_frag%P_mat_frag_pw_local, 4) /= dg_frag%nspin) then
      deallocate(dg_frag%P_mat_frag_pw_local)
      allocate(dg_frag%P_mat_frag_pw_local(3, nrow_local, npw_local, dg_frag%nspin))
    end if

    dg_frag%S_mat_frag_pw_local(:, :, :) = (0.0d0, 0.0d0)
    dg_frag%H_mat_frag_pw_local(:, :, :) = (0.0d0, 0.0d0)
    dg_frag%P_mat_frag_pw_local(:, :, :, :) = (0.0d0, 0.0d0)

    have_S_full = allocated(dg_frag%S_mat_frag_pw)
    have_H_full = allocated(dg_frag%H_mat_frag_pw)
    have_P_full = allocated(dg_frag%P_mat_frag_pw)
    do ispin_loop = 1, dg_frag%nspin
      do ipw = 1, npw_local
        global_pw = dg_frag%fp_local_pw_ids(ipw)
        if (global_pw < 1) cycle
        do irow = 1, nrow_local
          global_row = dg_frag%fp_local_row_ids(irow)
          if (global_row < 1) cycle
          if (have_S_full) then
            source_row = fragment_local_row_index_pw(dg_frag, global_row, ispin_loop, size(dg_frag%S_mat_frag_pw, 1))
            if (source_row > 0 .and. &
                global_pw <= size(dg_frag%S_mat_frag_pw, 2) .and. &
                ispin_loop <= size(dg_frag%S_mat_frag_pw, 3)) then
              dg_frag%S_mat_frag_pw_local(irow, ipw, ispin_loop) = &
                dg_frag%S_mat_frag_pw(source_row, global_pw, ispin_loop)
            end if
          end if
          if (have_H_full) then
            source_row = fragment_local_row_index_pw(dg_frag, global_row, ispin_loop, size(dg_frag%H_mat_frag_pw, 1))
            if (source_row > 0 .and. &
                global_pw <= size(dg_frag%H_mat_frag_pw, 2) .and. &
                ispin_loop <= size(dg_frag%H_mat_frag_pw, 3)) then
              dg_frag%H_mat_frag_pw_local(irow, ipw, ispin_loop) = &
                dg_frag%H_mat_frag_pw(source_row, global_pw, ispin_loop)
            end if
          end if
          if (have_P_full) then
            source_row = fragment_local_row_index_pw(dg_frag, global_row, ispin_loop, size(dg_frag%P_mat_frag_pw, 2))
            if (source_row > 0 .and. &
                global_pw <= size(dg_frag%P_mat_frag_pw, 3) .and. &
                ispin_loop <= size(dg_frag%P_mat_frag_pw, 4)) then
              dg_frag%P_mat_frag_pw_local(1:3, irow, ipw, ispin_loop) = &
                dg_frag%P_mat_frag_pw(1:3, source_row, global_pw, ispin_loop)
            end if
          else if (allocated(dg_frag%k_pw)) then
            if (global_pw <= size(dg_frag%k_pw, 2)) then
              do idir = 1, 3
                dg_frag%P_mat_frag_pw_local(idir, irow, ipw, ispin_loop) = &
                  zi * dg_frag%k_pw(idir, global_pw) * dg_frag%S_mat_frag_pw_local(irow, ipw, ispin_loop)
              end do
            end if
          end if
        end do
      end do
    end do
  end subroutine prepare_local_fragment_pw_blocks

  subroutine compute_pw_hamiltonian_local_potential(dg_frag, Vh, Vxc, Vpsl, H_pw)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    complex(8), intent(out) :: H_pw(:,:,:)  ! (n_pw, n_pw, nspin)

    integer :: n_pw, ispin, ifrag, ix, iy, iz, ipw, jpw
    integer :: nx, ny, nz, gx0, gy0, gz0, gx, gy, gz
    integer :: vx, vy, vz
    integer :: v_lb1, v_lb2, v_lb3, v_ub1, v_ub2, v_ub3
    real(8) :: vol_elem, box_volume, phase_arg, v_total
    complex(8), allocatable :: pw_phase(:), local_block(:,:), global_block(:,:)
    complex(8), allocatable :: T_wpw(:,:), S_wpw(:,:)
    logical :: use_wpw_weak, do_wpw_diag, use_unit_window

    n_pw = dg_frag%n_plane_waves
    H_pw(:, :, :) = (0.0d0, 0.0d0)
    if (n_pw <= 0) return
    use_wpw_weak = wpw_kinetic_use_weak()
    do_wpw_diag = wpw_kinetic_diag_enabled()
    use_unit_window = wpw_kinetic_unit_test_enabled()

    vol_elem = product(dg_frag%hgs(1:3))
    box_volume = product(dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3)))
    v_lb1 = lbound(Vpsl%f, 1)
    v_lb2 = lbound(Vpsl%f, 2)
    v_lb3 = lbound(Vpsl%f, 3)
    v_ub1 = ubound(Vpsl%f, 1)
    v_ub2 = ubound(Vpsl%f, 2)
    v_ub3 = ubound(Vpsl%f, 3)

    allocate(pw_phase(n_pw), local_block(n_pw, n_pw), global_block(n_pw, n_pw))

    do ispin = 1, dg_frag%nspin
      local_block(:, :) = (0.0d0, 0.0d0)

      if (dg_frag%ifrag_start <= dg_frag%ifrag_end) then
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          nx = dg_frag%nxyz_domain(1, ifrag)
          ny = dg_frag%nxyz_domain(2, ifrag)
          nz = dg_frag%nxyz_domain(3, ifrag)
          gx0 = dg_frag%ixyz_frag(1, ifrag)
          gy0 = dg_frag%ixyz_frag(2, ifrag)
          gz0 = dg_frag%ixyz_frag(3, ifrag)

          do iz = 1, nz
            gz = gz0 + iz - 1
            vz = map_global_to_phi_box_coord_pw(gz, v_lb3, v_ub3, dg_frag%lgnum_total(3))
            if (vz < v_lb3 .or. vz > v_ub3) cycle
            do iy = 1, ny
              gy = gy0 + iy - 1
              vy = map_global_to_phi_box_coord_pw(gy, v_lb2, v_ub2, dg_frag%lgnum_total(2))
              if (vy < v_lb2 .or. vy > v_ub2) cycle
              do ix = 1, nx
                gx = gx0 + ix - 1
                vx = map_global_to_phi_box_coord_pw(gx, v_lb1, v_ub1, dg_frag%lgnum_total(1))
                if (vx < v_lb1 .or. vx > v_ub1) cycle

                v_total = Vpsl%f(vx, vy, vz) + Vh%f(vx, vy, vz) + Vxc(ispin)%f(vx, vy, vz)
                do ipw = 1, n_pw
                  phase_arg = dg_frag%k_pw(1, ipw) * dble(gx) * dg_frag%hgs(1) + &
                              dg_frag%k_pw(2, ipw) * dble(gy) * dg_frag%hgs(2) + &
                              dg_frag%k_pw(3, ipw) * dble(gz) * dg_frag%hgs(3)
                  pw_phase(ipw) = exp(cmplx(0.0d0, phase_arg, kind=8))
                end do

                do jpw = 1, n_pw
                  do ipw = 1, n_pw
                    local_block(ipw, jpw) = local_block(ipw, jpw) + &
                      conjg(pw_phase(ipw)) * pw_phase(jpw) * v_total * vol_elem
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

      call comm_summation(local_block, global_block, n_pw * n_pw, dg_frag%icomm)
      H_pw(1:n_pw, 1:n_pw, ispin) = global_block(1:n_pw, 1:n_pw) / box_volume

      if (use_wpw_weak .or. do_wpw_diag) then
        if (.not. allocated(T_wpw)) then
          allocate(T_wpw(n_pw, n_pw), S_wpw(n_pw, n_pw))
          call compute_wpw_kinetic_weak(dg_frag, use_unit_window, T_wpw)
          if (do_wpw_diag) then
            call compute_wpw_overlap(dg_frag, use_unit_window, S_wpw)
            call diagnose_wpw_kinetic(dg_frag, T_wpw, S_wpw, use_unit_window)
          end if
        end if
      end if

      if (use_wpw_weak) then
        H_pw(1:n_pw, 1:n_pw, ispin) = H_pw(1:n_pw, 1:n_pw, ispin) + T_wpw(1:n_pw, 1:n_pw)
      else
        do ipw = 1, n_pw
          H_pw(ipw, ipw, ispin) = H_pw(ipw, ipw, ispin) + 0.5d0 * sum(dg_frag%k_pw(:, ipw)**2)
        end do
      end if

      do jpw = 1, n_pw
        do ipw = jpw + 1, n_pw
          H_pw(ipw, jpw, ispin) = 0.5d0 * (H_pw(ipw, jpw, ispin) + conjg(H_pw(jpw, ipw, ispin)))
          H_pw(jpw, ipw, ispin) = conjg(H_pw(ipw, jpw, ispin))
        end do
      end do
    end do

    if (allocated(T_wpw)) deallocate(T_wpw)
    if (allocated(S_wpw)) deallocate(S_wpw)
    deallocate(pw_phase, local_block, global_block)
  end subroutine compute_pw_hamiltonian_local_potential

  !=======================================================================
  ! Build mixed Hamiltonian matrix with plane wave basis
  !=======================================================================
  subroutine build_mixed_hamiltonian(dg_frag, lg, Vh, Vxc, Vpsl, Ac_tot, &
                                      S_frag_pw, H_frag_pw)
    use structures
    use communication, only: comm_is_root, comm_summation
    use rt_dg_fragment_ops, only: apply_complex_matrix_blocks_batch
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: lg
    type(s_scalar),         intent(in)    :: Vh, Vxc(:), Vpsl
    real(8),                intent(in)    :: Ac_tot(3)
    complex(8),             intent(in)    :: S_frag_pw(:,:,:)
    complex(8),             intent(inout) :: H_frag_pw(:,:,:)

    integer :: ispin, ipw, n_frag, n_pw
    complex(8), allocatable :: H_pw(:,:,:)
    complex(8), allocatable :: H_nl_fp(:,:), H_nl_pp(:,:), H_nl_pp_local(:,:)
    logical :: include_nl_mixed

    if (.not. dg_frag%use_plane_wave_basis) return

    n_frag = min(size(S_frag_pw, 1), size(H_frag_pw, 1))
    n_pw = min(dg_frag%n_plane_waves, size(S_frag_pw, 2), size(H_frag_pw, 2))

    include_nl_mixed = dg_frag%has_nl_cache .and. allocated(dg_frag%H_nl_blocks) .and. n_frag > 0 .and. n_pw > 0
    if (include_nl_mixed) then
      allocate(H_nl_fp(n_frag, n_pw), H_nl_pp(n_pw, n_pw), H_nl_pp_local(n_pw, n_pw))
      H_nl_fp(:, :) = (0.0d0, 0.0d0)
      H_nl_pp(:, :) = (0.0d0, 0.0d0)
      H_nl_pp_local(:, :) = (0.0d0, 0.0d0)
    end if

    if (.not. allocated(dg_frag%H_mat_pw_diag)) then
      allocate(dg_frag%H_mat_pw_diag(n_pw, dg_frag%nspin))
    else if (size(dg_frag%H_mat_pw_diag,1) /= n_pw .or. size(dg_frag%H_mat_pw_diag,2) /= dg_frag%nspin) then
      deallocate(dg_frag%H_mat_pw_diag)
      allocate(dg_frag%H_mat_pw_diag(n_pw, dg_frag%nspin))
    end if

    if (n_pw > 0) then
      if (.not. allocated(dg_frag%H_mat_pw)) then
        allocate(dg_frag%H_mat_pw(n_pw, n_pw, dg_frag%nspin))
      else if (size(dg_frag%H_mat_pw,1) /= n_pw .or. size(dg_frag%H_mat_pw,2) /= n_pw .or. &
               size(dg_frag%H_mat_pw,3) /= dg_frag%nspin) then
        deallocate(dg_frag%H_mat_pw)
        allocate(dg_frag%H_mat_pw(n_pw, n_pw, dg_frag%nspin))
      end if
      allocate(H_pw(n_pw, n_pw, dg_frag%nspin))
      call compute_pw_hamiltonian_local_potential(dg_frag, Vh, Vxc, Vpsl, H_pw)
      if (include_nl_mixed) then
        do ispin = 1, dg_frag%nspin
          H_nl_fp(:, :) = (0.0d0, 0.0d0)
          if (allocated(dg_frag%H_nl_local_block_ids)) then
            call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
                 S_frag_pw(1:n_frag, 1:n_pw, ispin), H_nl_fp(:, :), dg_frag%H_nl_local_block_ids)
          else
            call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, &
                 S_frag_pw(1:n_frag, 1:n_pw, ispin), H_nl_fp(:, :))
          end if
          H_frag_pw(1:n_frag, 1:n_pw, ispin) = H_frag_pw(1:n_frag, 1:n_pw, ispin) + H_nl_fp(:, :)
          H_nl_pp_local(:, :) = matmul(conjg(transpose(S_frag_pw(1:n_frag, 1:n_pw, ispin))), H_nl_fp(:, :))
          if (n_frag < dg_frag%n_mat_max) then
            call comm_summation(H_nl_pp_local, H_nl_pp, n_pw * n_pw, dg_frag%icomm)
          else
            H_nl_pp(:, :) = H_nl_pp_local(:, :)
          end if
          H_pw(1:n_pw, 1:n_pw, ispin) = H_pw(1:n_pw, 1:n_pw, ispin) + H_nl_pp(:, :)
        end do
      end if
      dg_frag%H_mat_pw(:, :, :) = H_pw(:, :, :)
      if (wpw_mixed_h_diag_enabled()) call diagnose_wpw_mixed_self_hamiltonian(dg_frag, Vh, Vxc, Vpsl)
      if (wpw_mixed_neighbor_h_diag_enabled()) call diagnose_wpw_mixed_neighbor_hamiltonian(dg_frag, Vh, Vxc, Vpsl)
    end if

    if (allocated(dg_frag%H_mat_frag_pw)) then
      if (size(dg_frag%H_mat_frag_pw,1) == n_frag .and. size(dg_frag%H_mat_frag_pw,2) == n_pw .and. &
          size(dg_frag%H_mat_frag_pw,3) == dg_frag%nspin) then
        dg_frag%H_mat_frag_pw(1:n_frag, 1:n_pw, 1:dg_frag%nspin) = H_frag_pw(1:n_frag, 1:n_pw, 1:dg_frag%nspin)
      end if
    end if

    do ispin = 1, dg_frag%nspin
      if (n_pw <= 0) cycle

      do ipw = 1, n_pw
        dg_frag%H_mat_pw_diag(ipw, ispin) = dg_frag%H_mat_pw(ipw, ipw, ispin)
      end do
    end do

    if (allocated(H_pw)) deallocate(H_pw)
    if (allocated(H_nl_fp)) deallocate(H_nl_fp)
    if (allocated(H_nl_pp)) deallocate(H_nl_pp)
    if (allocated(H_nl_pp_local)) deallocate(H_nl_pp_local)

    if (plane_wave_trace_enabled() .and. comm_is_root(dg_frag%id) .and. n_pw > 0) then
      write(*,'(1x,a)') '[HPP-MODE] full PW-PW potential enabled'
    end if

  end subroutine build_mixed_hamiltonian

  !=======================================================================
  ! Diagonalize mixed basis using generalized eigenproblem:
  !   H c = eps S c
  !=======================================================================
  subroutine diagonalize_mixed_basis(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    real(8), intent(in) :: Ac_tot(3)

    call prepare_mixed_basis_startup(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    dg_frag%mixed_basis_ready = dg_frag%mixed_basis_identity_raw
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "Mixed basis dense EVP skipped"
      write(*,'(1x,a)') "Mixed basis startup uses block/raw PW coupling path"
    end if
  end subroutine diagonalize_mixed_basis

  !=======================================================================
  ! Project PW components onto the fragment overlap-orthogonal complement:
  !   PW_perp = PW - Frag * (Sff^{-1} Sfp)
  ! and apply the same transform to fragment-PW Hamiltonian coupling.
  !=======================================================================
  subroutine project_pw_to_fragment_orthogonal_complement(dg_frag, eps_abs, eps_rel, S_frag_pw, H_frag_pw, P_frag_pw, protect_low_g)
    use communication, only: comm_is_root
    use rt_dg_fragment_ops, only: copy_matrix_blocks_to_complex_dense, copy_momentum_blocks_to_complex_dense
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: eps_abs, eps_rel
    complex(8), intent(inout) :: S_frag_pw(:,:,:)
    complex(8), intent(inout) :: H_frag_pw(:,:,:)
    complex(8), intent(inout) :: P_frag_pw(:,:,:,:)
    logical, intent(in) :: protect_low_g(:)

    integer :: ispin, n_frag, n_pw, info, lwork
    integer :: i, j, idir, n_keep_s, n_pw_dom, n_pw_nontriv
    integer :: n_proj, n_active_s
    real(8) :: eps_s, lam_max, lam_min_keep
    real(8) :: pw_eff_before, pw_eff_after
    real(8), parameter :: pw_dom_thresh = 0.5d0, pw_nontriv_thresh = 0.1d0
    complex(8), allocatable :: Sff(:,:), Sff_raw(:,:), Hff(:,:), Cproj(:,:), tmp_c(:,:), work_c(:)
    complex(8), allocatable :: Pff(:,:,:)
    complex(8), allocatable :: Pff_tmp(:,:)
    complex(8), allocatable :: Sfp_proj(:,:), Hfp_proj(:,:), Pfp_proj(:,:,:)
    integer, allocatable :: proj_idx(:)
    real(8), allocatable :: eval_sff(:), rwork(:)

    external :: zheev

    n_frag = dg_frag%n_mat_max
    n_pw = dg_frag%n_plane_waves
    if (n_frag <= 0 .or. n_pw <= 0) return

    allocate(Sff(n_frag, n_frag), Sff_raw(n_frag, n_frag), Hff(n_frag, n_frag))
    allocate(Pff(3, n_frag, n_frag))
    allocate(Pff_tmp(n_frag, n_frag))
    allocate(Cproj(n_frag, n_pw), tmp_c(n_frag, n_pw))
    allocate(eval_sff(n_frag), rwork(max(1, 3*n_frag-2)))

    do ispin = 1, dg_frag%nspin
      call get_fragment_overlap_dense(dg_frag, ispin, Sff_raw)
      call get_fragment_hamiltonian_dense(dg_frag, ispin, Hff)
      Pff(:, :, :) = (0.0d0, 0.0d0)
      do idir = 1, 3
        Pff_tmp(:, :) = (0.0d0, 0.0d0)
        select case (idir)
        case (1)
          call copy_momentum_blocks_to_complex_dense(dg_frag, ispin, [1.0d0, 0.0d0, 0.0d0], Pff_tmp)
        case (2)
          call copy_momentum_blocks_to_complex_dense(dg_frag, ispin, [0.0d0, 1.0d0, 0.0d0], Pff_tmp)
        case default
          call copy_momentum_blocks_to_complex_dense(dg_frag, ispin, [0.0d0, 0.0d0, 1.0d0], Pff_tmp)
        end select
        Pff(idir, :, :) = Pff_tmp(:, :)
      end do

      call evaluate_mixed_pw_effective_dim(Sff_raw, S_frag_pw(:, :, ispin), eps_abs, eps_rel, pw_eff_before, n_keep_s, n_pw_dom, n_pw_nontriv)

      Sff(:, :) = Sff_raw(:, :)
      lwork = -1
      allocate(work_c(1))
      call ZHEEV('V', 'U', n_frag, Sff, n_frag, eval_sff, work_c, lwork, rwork, info)
      lwork = int(real(work_c(1), kind=8)) + 1
      deallocate(work_c)
      allocate(work_c(lwork))
      call ZHEEV('V', 'U', n_frag, Sff, n_frag, eval_sff, work_c, lwork, rwork, info)
      deallocate(work_c)
      if (info /= 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,i0)') 'WARN: PW_perp projection skipped (Sff diagonalization failed), spin=', ispin, ' info=', info
        end if
        cycle
      end if

      lam_max = eval_sff(n_frag)
      eps_s = max(eps_abs, eps_rel * max(lam_max, 1.0d0))
      n_active_s = count(eval_sff(1:n_frag) >= eps_s)
      if (n_active_s <= 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,1pe11.3,a)') 'PW_perp projection skipped for spin=', ispin, ' (no Sff mode above eps_s=', eps_s, ')'
        end if
        cycle
      end if

      n_proj = count(.not. protect_low_g(1:n_pw))
      if (n_proj <= 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a)') 'PW_perp projection skipped for spin=', ispin, ' (all PW columns protected)'
        end if
        cycle
      end if

      allocate(proj_idx(n_proj), Sfp_proj(n_frag, n_proj), Hfp_proj(n_frag, n_proj), Pfp_proj(3, n_frag, n_proj))
      j = 0
      do i = 1, n_pw
        if (protect_low_g(i)) cycle
        j = j + 1
        proj_idx(j) = i
      end do
      Sfp_proj(:, :) = S_frag_pw(:, proj_idx(:), ispin)
      Hfp_proj(:, :) = H_frag_pw(:, proj_idx(:), ispin)
      Pfp_proj(:, :, :) = P_frag_pw(:, :, proj_idx(:), ispin)

      tmp_c(:, 1:n_proj) = matmul(conjg(transpose(Sff)), Sfp_proj)
      do i = 1, n_frag
        if (eval_sff(i) >= eps_s) then
          tmp_c(i, 1:n_proj) = tmp_c(i, 1:n_proj) / eval_sff(i)
        else
          tmp_c(i, 1:n_proj) = (0.0d0, 0.0d0)
        end if
      end do
      Cproj(:, 1:n_proj) = matmul(Sff, tmp_c(:, 1:n_proj))

      Sfp_proj(:, :) = Sfp_proj(:, :) - matmul(Sff_raw, Cproj(:, 1:n_proj))
      Hfp_proj(:, :) = Hfp_proj(:, :) - matmul(Hff, Cproj(:, 1:n_proj))
      do idir = 1, 3
        Pfp_proj(idir, :, :) = Pfp_proj(idir, :, :) - matmul(Pff(idir, :, :), Cproj(:, 1:n_proj))
      end do
      S_frag_pw(:, proj_idx(:), ispin) = Sfp_proj(:, :)
      H_frag_pw(:, proj_idx(:), ispin) = Hfp_proj(:, :)
      P_frag_pw(:, :, proj_idx(:), ispin) = Pfp_proj(:, :, :)

      deallocate(proj_idx, Sfp_proj, Hfp_proj, Pfp_proj)

      call evaluate_mixed_pw_effective_dim(Sff_raw, S_frag_pw(:, :, ispin), eps_abs, eps_rel, pw_eff_after, n_keep_s, n_pw_dom, n_pw_nontriv)

      if (comm_is_root(dg_frag%id)) then
         lam_min_keep = minval(eval_sff, mask=(eval_sff >= eps_s))
         write(*,'(1x,a,i0,a,1pe11.3,a,1pe11.3,a,1pe11.3,a,i0)') 'PW_perp projection (spin=', ispin, '): pw_eff_before=', &
           pw_eff_before, ' pw_eff_after=', pw_eff_after, ' lam_min_keep=', lam_min_keep, ' n_active_s=', n_active_s
      end if
    end do

    deallocate(Sff, Sff_raw, Hff, Pff, Pff_tmp, Cproj, tmp_c, eval_sff, rwork)
  end subroutine project_pw_to_fragment_orthogonal_complement

  subroutine select_protected_low_g_modes(dg_frag, protect_low_g, n_protect, k_thr)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    logical, allocatable, intent(out) :: protect_low_g(:)
    integer, intent(out) :: n_protect
    real(8), intent(out) :: k_thr

    integer :: n_pw, ifrag, idir, ipw, env_len, env_stat, info
    real(8) :: l_frag_max, k_norm
    real(8) :: protect_scale
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    character(len=64) :: env_scale

    n_pw = dg_frag%n_plane_waves
    allocate(protect_low_g(max(0, n_pw)))
    if (n_pw <= 0) then
      n_protect = 0
      k_thr = 0.0d0
      return
    end if
    protect_low_g(:) = .false.

    protect_scale = 0.0d0
    env_scale = ''
    call get_environment_variable('SALMON_DG_PW_PROTECT_SCALE', env_scale, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_scale(1:env_len), *, iostat=info) protect_scale
      if (info /= 0) protect_scale = 0.0d0
    end if
    if (protect_scale <= 0.0d0) then
      n_protect = 0
      k_thr = 0.0d0
      return
    end if

    l_frag_max = 0.0d0
    if (allocated(dg_frag%nxyz_domain)) then
      do ifrag = 1, min(dg_frag%n_frag, size(dg_frag%nxyz_domain, 2))
        do idir = 1, 3
          l_frag_max = max(l_frag_max, dg_frag%hgs(idir) * dble(max(1, dg_frag%nxyz_domain(idir, ifrag))))
        end do
      end do
    end if
    if (l_frag_max <= 0.0d0) then
      n_protect = 0
      k_thr = 0.0d0
      return
    end if

    k_thr = protect_scale * (2.0d0 * pi / l_frag_max)
    do ipw = 1, n_pw
      k_norm = sqrt(sum(dg_frag%k_pw(:, ipw)**2))
      if (k_norm <= k_thr) protect_low_g(ipw) = .true.
    end do
    n_protect = count(protect_low_g)
  end subroutine select_protected_low_g_modes

  subroutine get_fragment_overlap_dense(dg_frag, ispin, Sff)
    use rt_dg_fragment_ops, only: copy_matrix_blocks_to_complex_dense
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(out) :: Sff(:,:)

    integer :: n_frag, i, j

    n_frag = dg_frag%n_mat_max
    Sff(:, :) = (0.0d0, 0.0d0)
    if (allocated(dg_frag%S_mat_c)) then
      Sff(1:n_frag, 1:n_frag) = dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin)
    else if (allocated(dg_frag%S_mat_blocks)) then
      call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%S_mat_blocks, ispin, Sff)
    else if (allocated(dg_frag%S_mat)) then
      Sff(1:n_frag, 1:n_frag) = cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
    else
      do i = 1, n_frag
        Sff(i, i) = (1.0d0, 0.0d0)
      end do
    end if
    do j = 1, n_frag
      Sff(j, j) = cmplx(real(Sff(j, j), kind=8), 0.0d0, kind=8)
      do i = j + 1, n_frag
        Sff(i, j) = conjg(Sff(j, i))
      end do
    end do
  end subroutine get_fragment_overlap_dense

  subroutine get_fragment_hamiltonian_dense(dg_frag, ispin, Hff)
    use rt_dg_fragment_ops, only: copy_matrix_blocks_to_complex_dense
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(out) :: Hff(:,:)

    integer :: n_frag

    n_frag = dg_frag%n_mat_max
    Hff(:, :) = (0.0d0, 0.0d0)
    if (allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)) then
      Hff(1:n_frag, 1:n_frag) = dg_frag%H_mat_c(1:n_frag, 1:n_frag, ispin)
    else if (allocated(dg_frag%H_mat_blocks)) then
      call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, Hff)
    else if (allocated(dg_frag%H_mat)) then
      Hff(1:n_frag, 1:n_frag) = cmplx(dg_frag%H_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
    end if
  end subroutine get_fragment_hamiltonian_dense

  subroutine evaluate_mixed_pw_effective_dim(Sff, Sfp, eps_abs, eps_rel, pw_eff, n_keep_s, n_pw_dom, n_pw_nontriv)
    implicit none
    complex(8), intent(in) :: Sff(:,:), Sfp(:,:)
    real(8), intent(in) :: eps_abs, eps_rel
    real(8), intent(out) :: pw_eff
    integer, intent(out) :: n_keep_s, n_pw_dom, n_pw_nontriv

    integer :: n_frag, n_pw, n_total, i, j, info, lwork, k
    real(8) :: eps_s, tau_s, pw_weight
    real(8), parameter :: pw_dom_thresh = 0.5d0, pw_nontriv_thresh = 0.1d0
    complex(8), allocatable :: Smix(:,:), work_c(:)
    real(8), allocatable :: eval_s(:), rwork(:)
    integer, allocatable :: keep_idx(:)

    external :: zheev

    n_frag = size(Sff, 1)
    n_pw = size(Sfp, 2)
    n_total = n_frag + n_pw

    pw_eff = 0.0d0
    n_keep_s = 0
    n_pw_dom = 0
    n_pw_nontriv = 0
    if (n_total <= 0) return

    allocate(Smix(n_total, n_total), eval_s(n_total), rwork(max(1, 3*n_total-2)))
    Smix(:, :) = (0.0d0, 0.0d0)
    Smix(1:n_frag, 1:n_frag) = Sff(:, :)
    if (n_pw > 0) then
      Smix(1:n_frag, n_frag+1:n_total) = Sfp(:, :)
      Smix(n_frag+1:n_total, 1:n_frag) = conjg(transpose(Sfp(:, :)))
      do i = n_frag + 1, n_total
        Smix(i, i) = (1.0d0, 0.0d0)
      end do
    end if

    lwork = -1
    allocate(work_c(1))
    call ZHEEV('V', 'U', n_total, Smix, n_total, eval_s, work_c, lwork, rwork, info)
    lwork = int(real(work_c(1), kind=8)) + 1
    deallocate(work_c)
    allocate(work_c(lwork))
    call ZHEEV('V', 'U', n_total, Smix, n_total, eval_s, work_c, lwork, rwork, info)
    deallocate(work_c)
    if (info /= 0) then
      deallocate(Smix, eval_s, rwork)
      return
    end if

    eps_s = max(eps_abs, eps_rel * max(eval_s(n_total), 1.0d0))
    tau_s = eps_s
    do i = 1, n_total
      if (eval_s(i) >= tau_s) n_keep_s = n_keep_s + 1
    end do
    if (n_keep_s <= 0) n_keep_s = 1

    allocate(keep_idx(n_keep_s))
    k = 0
    do i = 1, n_total
      if (eval_s(i) >= tau_s) then
        k = k + 1
        if (k <= n_keep_s) keep_idx(k) = i
      end if
    end do
    if (k < n_keep_s) then
      do i = k + 1, n_keep_s
        keep_idx(i) = n_total - n_keep_s + i
      end do
    end if

    if (n_pw > 0) then
      do j = 1, n_keep_s
        i = keep_idx(j)
        pw_weight = sum(abs(Smix(n_frag+1:n_total, i))**2)
        pw_eff = pw_eff + pw_weight
        if (pw_weight >= pw_dom_thresh) n_pw_dom = n_pw_dom + 1
        if (pw_weight >= pw_nontriv_thresh) n_pw_nontriv = n_pw_nontriv + 1
      end do
    end if

    deallocate(keep_idx, Smix, eval_s, rwork)
  end subroutine evaluate_mixed_pw_effective_dim

  subroutine compute_pw_weight_proxy_from_overlap(S_frag_pw, pw_proxy)
    implicit none
    complex(8), intent(in) :: S_frag_pw(:,:,:)
    real(8), intent(out) :: pw_proxy(:)

    integer :: n_frag, n_pw, nspin, ispin, ipw
    real(8) :: frag_overlap_norm

    n_frag = size(S_frag_pw, 1)
    n_pw = size(S_frag_pw, 2)
    nspin = size(S_frag_pw, 3)
    pw_proxy(:) = 0.0d0
    if (n_pw <= 0) return
    do ipw = 1, n_pw
      frag_overlap_norm = 0.0d0
      do ispin = 1, nspin
        frag_overlap_norm = frag_overlap_norm + sum(abs(S_frag_pw(1:n_frag, ipw, ispin))**2)
      end do
      pw_proxy(ipw) = 1.0d0 / (1.0d0 + frag_overlap_norm)
    end do
  end subroutine compute_pw_weight_proxy_from_overlap

  subroutine build_keep_sets_from_jpvt(n_pw, jpvt, n_keep_base, pw_proxy, n_pw_weight_protect, &
                                       keep_jpvt, keep_pw_protect, n_keep_pw_protect)
    implicit none
    integer, intent(in) :: n_pw, jpvt(:), n_keep_base, n_pw_weight_protect
    real(8), intent(in) :: pw_proxy(:)
    integer, allocatable, intent(out) :: keep_jpvt(:), keep_pw_protect(:)
    integer, intent(out) :: n_keep_pw_protect

    integer :: i, idx, n_valid_base, n_extra, best_idx, n_added
    real(8) :: best_val
    logical, allocatable :: selected(:)

    n_valid_base = max(0, min(n_keep_base, n_pw))
    allocate(keep_jpvt(max(1, n_valid_base)))
    allocate(keep_pw_protect(max(1, n_pw)))
    allocate(selected(n_pw))
    selected(:) = .false.

    n_valid_base = 0
    do i = 1, min(size(jpvt), n_pw)
      if (n_valid_base >= n_keep_base) exit
      idx = jpvt(i)
      if (idx < 1 .or. idx > n_pw) cycle
      if (selected(idx)) cycle
      n_valid_base = n_valid_base + 1
      keep_jpvt(n_valid_base) = idx
      keep_pw_protect(n_valid_base) = idx
      selected(idx) = .true.
    end do

    n_extra = min(max(0, n_pw_weight_protect), n_pw - n_valid_base)
    n_added = 0
    do i = 1, n_extra
      best_idx = 0
      best_val = -1.0d0
      do idx = 1, n_pw
        if (selected(idx)) cycle
        if (pw_proxy(idx) > best_val) then
          best_val = pw_proxy(idx)
          best_idx = idx
        end if
      end do
      if (best_idx <= 0) exit
      n_added = n_added + 1
      keep_pw_protect(n_valid_base + n_added) = best_idx
      selected(best_idx) = .true.
    end do

    n_keep_pw_protect = n_valid_base + n_added
    if (n_valid_base <= 0) keep_jpvt(1) = 1
    deallocate(selected)
  end subroutine build_keep_sets_from_jpvt

  !=======================================================================
  ! Re-check the final PW keep set after protection has been merged in.
  ! Protection is intentionally soft here: a protected PW mode is retained
  ! only if it adds an independent direction to the fragment-overlap space.
  ! The modified Gram-Schmidt pass is deterministic because it follows the
  ! existing keep order and does not depend on nearly tied QR pivots.
  !=======================================================================
  subroutine filter_pw_keep_by_real_rank(S_frag_pw, keep_idx, n_keep, tau_rel, n_dropped)
    implicit none
    complex(8), intent(in) :: S_frag_pw(:,:,:)
    integer, intent(inout) :: keep_idx(:)
    integer, intent(inout) :: n_keep
    real(8), intent(in) :: tau_rel
    integer, intent(out) :: n_dropped

    integer :: n_frag, n_pw, nspin, m_real, j, k, idx, n_new
    real(8) :: max_col_norm, tau_abs, norm_col, norm_res, proj
    real(8), allocatable :: q(:,:), col(:)
    integer, allocatable :: keep_new(:)

    n_dropped = 0
    if (n_keep <= 1) return

    n_frag = size(S_frag_pw, 1)
    n_pw = size(S_frag_pw, 2)
    nspin = size(S_frag_pw, 3)
    if (n_frag <= 0 .or. n_pw <= 0 .or. nspin <= 0) return

    m_real = 2 * n_frag * nspin
    allocate(q(m_real, n_keep), col(m_real), keep_new(n_keep))

    max_col_norm = 0.0d0
    do j = 1, n_keep
      idx = keep_idx(j)
      if (idx < 1 .or. idx > n_pw) cycle
      call fill_real_pw_column(idx, col)
      norm_col = sqrt(sum(col(:) * col(:)))
      max_col_norm = max(max_col_norm, norm_col)
    end do

    tau_abs = max(1.0d-14, tau_rel * max(max_col_norm, 1.0d0))
    q(:, :) = 0.0d0
    n_new = 0
    do j = 1, n_keep
      idx = keep_idx(j)
      if (idx < 1 .or. idx > n_pw) cycle
      call fill_real_pw_column(idx, col)
      norm_col = sqrt(sum(col(:) * col(:)))
      if (norm_col < tau_abs) cycle

      do k = 1, n_new
        proj = dot_product(q(:, k), col(:))
        col(:) = col(:) - proj * q(:, k)
      end do
      norm_res = sqrt(sum(col(:) * col(:)))

      if (norm_res >= tau_abs) then
        n_new = n_new + 1
        q(:, n_new) = col(:) / norm_res
        keep_new(n_new) = idx
      end if
    end do

    if (n_new <= 0) then
      do j = 1, n_keep
        idx = keep_idx(j)
        if (idx < 1 .or. idx > n_pw) cycle
        n_new = 1
        keep_new(1) = idx
        exit
      end do
    end if

    if (n_new > 0) keep_idx(1:n_new) = keep_new(1:n_new)
    n_dropped = max(0, n_keep - n_new)
    n_keep = n_new

    deallocate(q, col, keep_new)

  contains

    subroutine fill_real_pw_column(idx_col, vec)
      implicit none
      integer, intent(in) :: idx_col
      real(8), intent(out) :: vec(:)

      integer :: ispin_l, ifrag_l, row_l

      row_l = 0
      do ispin_l = 1, nspin
        do ifrag_l = 1, n_frag
          row_l = row_l + 1
          vec(row_l) = real(S_frag_pw(ifrag_l, idx_col, ispin_l), kind=8)
        end do
        do ifrag_l = 1, n_frag
          row_l = row_l + 1
          vec(row_l) = aimag(S_frag_pw(ifrag_l, idx_col, ispin_l))
        end do
      end do
    end subroutine fill_real_pw_column
  end subroutine filter_pw_keep_by_real_rank

  !=======================================================================
  ! Ensure the final fragment+PW metric is full-rank for every spin channel.
  ! When protected PW modes recreate a near-null mixed-S direction, the PW
  ! component with the largest weight in the discarded metric eigenvectors is
  ! removed. Ties are resolved by dropping the later keep entry so the earlier
  ! QR/proxy ordering remains stable.
  !=======================================================================
  subroutine filter_pw_keep_by_mixed_metric(dg_frag, S_frag_pw, keep_idx, n_keep, n_preferred, eps_abs, eps_rel, n_dropped)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    complex(8), intent(in) :: S_frag_pw(:,:,:)
    integer, intent(inout) :: keep_idx(:)
    integer, intent(inout) :: n_keep
    integer, intent(in) :: n_preferred
    real(8), intent(in) :: eps_abs, eps_rel
    integer, intent(out) :: n_dropped

    integer :: n_frag, n_pw, nspin, n_total, n_total_max
    integer :: iter, ispin, i, j, idx, info, lwork, drop_pos, preferred_limit
    real(8) :: tau_s, best_score
    logical :: full_rank
    complex(8), allocatable :: Sff(:,:), Smix(:,:), work_c(:)
    real(8), allocatable :: eval_s(:), rwork(:), drop_score(:)

    external :: zheev

    n_dropped = 0
    if (n_keep <= 0) return

    n_frag = size(S_frag_pw, 1)
    n_pw = size(S_frag_pw, 2)
    nspin = size(S_frag_pw, 3)
    if (n_frag <= 0 .or. n_pw <= 0 .or. nspin <= 0) return

    n_total_max = n_frag + n_keep
    allocate(Sff(n_frag, n_frag), Smix(n_total_max, n_total_max))
    allocate(eval_s(n_total_max), rwork(max(1, 3*n_total_max - 2)), drop_score(n_keep))

    do iter = 1, n_pw
      if (n_keep <= 0) exit
      n_total = n_frag + n_keep
      drop_score(1:n_keep) = 0.0d0
      full_rank = .true.

      do ispin = 1, nspin
        call get_fragment_overlap_dense(dg_frag, ispin, Sff)
        Smix(1:n_total, 1:n_total) = (0.0d0, 0.0d0)
        Smix(1:n_frag, 1:n_frag) = Sff(:, :)
        do j = 1, n_keep
          idx = keep_idx(j)
          if (idx < 1 .or. idx > n_pw) cycle
          do i = 1, n_frag
            Smix(i, n_frag + j) = S_frag_pw(i, idx, ispin)
            Smix(n_frag + j, i) = conjg(Smix(i, n_frag + j))
          end do
          Smix(n_frag + j, n_frag + j) = (1.0d0, 0.0d0)
        end do
        do i = 1, n_total
          Smix(i, i) = cmplx(real(Smix(i, i), kind=8), 0.0d0, kind=8)
        end do

        lwork = -1
        allocate(work_c(1))
        call ZHEEV('V', 'U', n_total, Smix, n_total_max, eval_s, work_c, lwork, rwork, info)
        lwork = max(1, int(real(work_c(1), kind=8)) + 1)
        deallocate(work_c)
        allocate(work_c(lwork))
        call ZHEEV('V', 'U', n_total, Smix, n_total_max, eval_s, work_c, lwork, rwork, info)
        deallocate(work_c)

        if (info /= 0) then
          full_rank = .false.
          drop_score(n_keep) = drop_score(n_keep) + 1.0d0
          cycle
        end if

        tau_s = max(eps_abs, eps_rel * max(eval_s(n_total), 1.0d0))
        do i = 1, n_total
          if (eval_s(i) >= tau_s) cycle
          full_rank = .false.
          do j = 1, n_keep
            drop_score(j) = drop_score(j) + abs(Smix(n_frag + j, i))**2
          end do
        end do
      end do

      if (full_rank) exit

      preferred_limit = max(0, min(n_preferred, n_keep))
      drop_pos = 0
      best_score = -1.0d0
      do j = preferred_limit + 1, n_keep
        if (drop_score(j) > best_score + 1.0d-18) then
          best_score = drop_score(j)
          drop_pos = j
        end if
      end do
      if (drop_pos <= 0) then
        drop_pos = n_keep
        best_score = drop_score(drop_pos)
        do j = 1, preferred_limit
          if (drop_score(j) > best_score + 1.0d-18) then
            best_score = drop_score(j)
            drop_pos = j
          end if
        end do
      end if

      if (drop_pos < n_keep) keep_idx(drop_pos:n_keep-1) = keep_idx(drop_pos+1:n_keep)
      n_keep = n_keep - 1
      n_dropped = n_dropped + 1
    end do

    deallocate(Sff, Smix, eval_s, rwork, drop_score)
  end subroutine filter_pw_keep_by_mixed_metric

  !=======================================================================
  ! Keep only rank-revealed independent PW columns.
  !=======================================================================
  subroutine compact_plane_wave_basis(dg_frag, S_frag_pw, H_frag_pw, P_frag_pw, keep_idx, n_keep)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), allocatable, intent(inout) :: S_frag_pw(:,:,:)
    complex(8), allocatable, intent(inout) :: H_frag_pw(:,:,:)
    complex(8), allocatable, intent(inout) :: P_frag_pw(:,:,:,:)
    integer, intent(in) :: keep_idx(:)
    integer, intent(in) :: n_keep

    integer :: n_frag, n_pw_old, i
    real(8), allocatable :: k_new(:,:)
    complex(8), allocatable :: coef_new(:,:,:), Sfp_new(:,:,:), Hfp_new(:,:,:), Pfp_new(:,:,:,:)
    complex(8), allocatable :: Hpp_diag_new(:,:), Hpp_full_new(:,:,:)

    n_frag = dg_frag%n_mat_max
    n_pw_old = dg_frag%n_plane_waves
    if (n_keep <= 0 .or. n_keep >= n_pw_old) return

    allocate(k_new(3, n_keep), coef_new(n_keep, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(Sfp_new(n_frag, n_keep, dg_frag%nspin), Hfp_new(n_frag, n_keep, dg_frag%nspin), &
             Pfp_new(3, n_frag, n_keep, dg_frag%nspin), Hpp_diag_new(n_keep, dg_frag%nspin), &
             Hpp_full_new(n_keep, n_keep, dg_frag%nspin))
    coef_new(:, :, :) = (0.0d0, 0.0d0)
    Hpp_diag_new(:, :) = (0.0d0, 0.0d0)
    Hpp_full_new(:, :, :) = (0.0d0, 0.0d0)
    do i = 1, n_keep
      k_new(:, i) = dg_frag%k_pw(:, keep_idx(i))
      coef_new(i, :, :) = dg_frag%coef_pw(keep_idx(i), :, :)
      Sfp_new(:, i, :) = S_frag_pw(:, keep_idx(i), :)
      Hfp_new(:, i, :) = H_frag_pw(:, keep_idx(i), :)
      Pfp_new(:, :, i, :) = P_frag_pw(:, :, keep_idx(i), :)
      if (allocated(dg_frag%H_mat_pw_diag)) Hpp_diag_new(i, :) = dg_frag%H_mat_pw_diag(keep_idx(i), :)
    end do
    if (allocated(dg_frag%H_mat_pw)) then
      do i = 1, n_keep
        Hpp_full_new(i, :, :) = dg_frag%H_mat_pw(keep_idx(i), keep_idx(1:n_keep), :)
      end do
    end if

    deallocate(dg_frag%k_pw, dg_frag%coef_pw)
    if (allocated(dg_frag%coef_pw_owner)) deallocate(dg_frag%coef_pw_owner)
    if (allocated(dg_frag%coef_pw_full_cache)) deallocate(dg_frag%coef_pw_full_cache)
    dg_frag%coef_pw_full_cache_nstate = 0
    if (allocated(dg_frag%S_mat_frag_pw)) deallocate(dg_frag%S_mat_frag_pw)
    if (allocated(dg_frag%H_mat_frag_pw)) deallocate(dg_frag%H_mat_frag_pw)
    if (allocated(dg_frag%P_mat_frag_pw)) deallocate(dg_frag%P_mat_frag_pw)
    if (allocated(dg_frag%fp_local_row_ids)) deallocate(dg_frag%fp_local_row_ids)
    if (allocated(dg_frag%fp_local_pw_ids)) deallocate(dg_frag%fp_local_pw_ids)
    if (allocated(dg_frag%S_mat_frag_pw_local)) deallocate(dg_frag%S_mat_frag_pw_local)
    if (allocated(dg_frag%H_mat_frag_pw_local)) deallocate(dg_frag%H_mat_frag_pw_local)
    if (allocated(dg_frag%P_mat_frag_pw_local)) deallocate(dg_frag%P_mat_frag_pw_local)
    if (allocated(dg_frag%H_mat_pw_diag)) deallocate(dg_frag%H_mat_pw_diag)
    if (allocated(dg_frag%H_mat_pw)) deallocate(dg_frag%H_mat_pw)

    dg_frag%n_plane_waves = n_keep
    dg_frag%owned_coef_pw_start = 0
    dg_frag%owned_coef_pw_end = -1
    allocate(dg_frag%k_pw(3, n_keep), dg_frag%coef_pw(n_keep, dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%k_pw(:, :) = k_new(:, :)
    dg_frag%coef_pw(:, :, :) = coef_new(:, :, :)
    allocate(dg_frag%coef_pw_owner(dg_frag%n_plane_waves))
    do i = 1, dg_frag%n_plane_waves
      dg_frag%coef_pw_owner(i) = min(dg_frag%isize - 1, ((i - 1) * dg_frag%isize) / dg_frag%n_plane_waves)
    end do
    do i = 1, dg_frag%n_plane_waves
      if (dg_frag%coef_pw_owner(i) /= dg_frag%id) cycle
      if (dg_frag%owned_coef_pw_start == 0) dg_frag%owned_coef_pw_start = i
      dg_frag%owned_coef_pw_end = i
    end do
    allocate(dg_frag%S_mat_frag_pw(n_frag, n_keep, dg_frag%nspin))
    allocate(dg_frag%H_mat_frag_pw(n_frag, n_keep, dg_frag%nspin), dg_frag%P_mat_frag_pw(3, n_frag, n_keep, dg_frag%nspin))
    allocate(dg_frag%H_mat_pw_diag(n_keep, dg_frag%nspin), dg_frag%H_mat_pw(n_keep, n_keep, dg_frag%nspin))
    dg_frag%S_mat_frag_pw(:, :, :) = Sfp_new(:, :, :)
    dg_frag%H_mat_frag_pw(:, :, :) = Hfp_new(:, :, :)
    dg_frag%P_mat_frag_pw(:, :, :, :) = Pfp_new(:, :, :, :)
    dg_frag%H_mat_pw_diag(:, :) = Hpp_diag_new(:, :)
    dg_frag%H_mat_pw(:, :, :) = Hpp_full_new(:, :, :)

    if (allocated(S_frag_pw)) deallocate(S_frag_pw)
    if (allocated(H_frag_pw)) deallocate(H_frag_pw)
    if (allocated(P_frag_pw)) deallocate(P_frag_pw)
    allocate(S_frag_pw(n_frag, n_keep, dg_frag%nspin), H_frag_pw(n_frag, n_keep, dg_frag%nspin), &
             P_frag_pw(3, n_frag, n_keep, dg_frag%nspin))
    S_frag_pw(:, :, :) = Sfp_new(:, :, :)
    H_frag_pw(:, :, :) = Hfp_new(:, :, :)
    P_frag_pw(:, :, :, :) = Pfp_new(:, :, :, :)

    deallocate(k_new, coef_new, Sfp_new, Hfp_new, Pfp_new, Hpp_diag_new, Hpp_full_new)
  end subroutine compact_plane_wave_basis

  subroutine build_fragment_pw_total_potential_box(dg_frag, ifrag, Vh, Vxc_spin, Vpsl, gx0, gy0, gz0, V_box)
    use structures, only: s_scalar
    use communication, only: comm_summation, COMM_GROUP_NULL
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    type(s_scalar), intent(in) :: Vh, Vxc_spin, Vpsl
    integer, intent(in) :: gx0, gy0, gz0
    real(8), intent(inout) :: V_box(:,:,:)

    integer :: ix, iy, iz, gx, gy, gz, vx, vy, vz
    integer :: v_lb1, v_lb2, v_lb3, v_ub1, v_ub2, v_ub3
    integer :: npt
    real(8), allocatable :: V_sum(:,:,:)

    V_box(:, :, :) = 0.0d0
    v_lb1 = lbound(Vpsl%f, 1)
    v_lb2 = lbound(Vpsl%f, 2)
    v_lb3 = lbound(Vpsl%f, 3)
    v_ub1 = ubound(Vpsl%f, 1)
    v_ub2 = ubound(Vpsl%f, 2)
    v_ub3 = ubound(Vpsl%f, 3)

    ! H_fp is column-partitioned in orbital mode.  Each column owner needs the
    ! scalar potential over the whole fragment box, so the parent-grid pieces
    ! held by the subgroup ranks are gathered once into this local box.
!$omp parallel do private(iy,ix,gz,gy,gx,vz,vy,vx) schedule(static)
    do iz = lbound(V_box, 3), ubound(V_box, 3)
      gz = gz0 + iz - 1
      vz = map_global_to_phi_box_coord_pw(gz, v_lb3, v_ub3, dg_frag%lgnum_total(3))
      do iy = lbound(V_box, 2), ubound(V_box, 2)
        gy = gy0 + iy - 1
        vy = map_global_to_phi_box_coord_pw(gy, v_lb2, v_ub2, dg_frag%lgnum_total(2))
!$omp simd private(gx,vx)
        do ix = lbound(V_box, 1), ubound(V_box, 1)
          gx = gx0 + ix - 1
          vx = map_global_to_phi_box_coord_pw(gx, v_lb1, v_ub1, dg_frag%lgnum_total(1))
          if (vx < v_lb1 .or. vx > v_ub1 .or. &
              vy < v_lb2 .or. vy > v_ub2 .or. &
              vz < v_lb3 .or. vz > v_ub3) cycle
          V_box(ix, iy, iz) = Vpsl%f(vx, vy, vz) + Vh%f(vx, vy, vz) + Vxc_spin%f(vx, vy, vz)
        end do
      end do
    end do
!$omp end parallel do

    if (dg_frag%isize_frag > 1 .and. dg_frag%icomm_frag /= COMM_GROUP_NULL) then
      npt = size(V_box)
      allocate(V_sum(lbound(V_box,1):ubound(V_box,1), &
                     lbound(V_box,2):ubound(V_box,2), &
                     lbound(V_box,3):ubound(V_box,3)))
      call comm_summation(V_box, V_sum, npt, dg_frag%icomm_frag)
      V_box(:, :, :) = V_sum(:, :, :)
      deallocate(V_sum)
    end if
  end subroutine build_fragment_pw_total_potential_box

  subroutine get_fragment_pw_column_range(dg_frag, ncol, col_s, col_e)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ncol
    integer, intent(out) :: col_s, col_e

    integer :: base, extra, rank_in_frag, nworker

    if (ncol <= 0) then
      col_s = 1
      col_e = 0
      return
    end if
    if (.not. dg_frag%parallel_mode_orbital .or. dg_frag%isize_frag <= 1) then
      col_s = 1
      col_e = ncol
      return
    end if

    nworker = max(1, dg_frag%isize_frag)
    rank_in_frag = max(0, min(dg_frag%id_frag, nworker - 1))
    base = ncol / nworker
    extra = mod(ncol, nworker)
    if (rank_in_frag < extra) then
      col_s = rank_in_frag * (base + 1) + 1
      col_e = col_s + base
    else
      col_s = extra * (base + 1) + (rank_in_frag - extra) * base + 1
      col_e = col_s + base - 1
    end if
    if (col_s > ncol) then
      col_s = 1
      col_e = 0
    else
      col_e = min(col_e, ncol)
    end if
  end subroutine get_fragment_pw_column_range

  subroutine get_fragment_subgroup_box_range_pw(dg_frag, nbox, loc_s, loc_e)
    use salmon_global, only: nproc_rgrid
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: nbox(3)
    integer, intent(out) :: loc_s(3), loc_e(3)

    integer :: axis, nproc_axis(3), coords(3), nsize

    nproc_axis(:) = max(1, nproc_rgrid(:))
    coords(1) = mod(dg_frag%id_frag, nproc_axis(1))
    coords(2) = mod(dg_frag%id_frag / nproc_axis(1), nproc_axis(2))
    coords(3) = dg_frag%id_frag / max(1, nproc_axis(1) * nproc_axis(2))

    do axis = 1, 3
      if (nbox(axis) <= 0) then
        loc_s(axis) = 1
        loc_e(axis) = 0
      else
        nsize = (nbox(axis) + nproc_axis(axis) - 1) / nproc_axis(axis)
        loc_s(axis) = 1 + nsize * coords(axis)
        loc_e(axis) = min(nbox(axis), loc_s(axis) + nsize - 1)
      end if
    end do
  end subroutine get_fragment_subgroup_box_range_pw

  !=======================================================================
  ! Map global periodic grid index to this rank's phi-box coordinate.
  ! If the point is not represented in the local box, caller must skip it.
  !=======================================================================

  pure integer function map_global_to_fft_slot_pw(ig, lgnum) result(slot_idx)
    implicit none
    integer, intent(in) :: ig, lgnum

    slot_idx = modulo(ig, lgnum) + 1
  end function map_global_to_fft_slot_pw

  pure integer function map_fft_mode_to_slot_pw(ik, lgnum) result(slot_idx)
    implicit none
    integer, intent(in) :: ik, lgnum

    slot_idx = modulo(ik, lgnum) + 1
  end function map_fft_mode_to_slot_pw

end module rt_dg_plane_wave
