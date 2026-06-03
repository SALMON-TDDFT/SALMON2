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

#include "config.h"
module rt_dg_initial_state
  use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info, complex_matrix_block_info
  use rt_dg_fragment_ops, only: apply_matrix_blocks_batch, apply_complex_matrix_blocks_batch, &
                                apply_overlap_operator_batch_orbital_fragment_self, &
                                apply_overlap_operator_batch, solve_overlap_operator_batch
  implicit none

  private
  public :: measure_fragment_initial_surface_residual
  public :: diagonalize_initial_dg_full_distributed
  public :: relax_initial_occupied_subspace_block_sparse
  public :: solve_fragment_generalized_eigen
  public :: build_fd_occupations
  public :: build_face_cluster_fragment_list

contains

  subroutine dg_scalapack_index_1d(ig, block_size, nproc_axis, proc, loc)
    implicit none
    integer, intent(in) :: ig, block_size, nproc_axis
    integer, intent(out) :: proc, loc
    integer :: iblock, inblock, local_block

    if (ig <= 0 .or. block_size <= 0 .or. nproc_axis <= 0) then
      proc = 0
      loc = 0
      return
    end if
    iblock = (ig - 1) / block_size
    inblock = mod(ig - 1, block_size) + 1
    proc = mod(iblock, nproc_axis)
    local_block = iblock / nproc_axis
    loc = local_block * block_size + inblock
  end subroutine dg_scalapack_index_1d

  subroutine dg_scalapack_global_index_1d(loc, block_size, nproc_axis, myproc_axis, ig)
    implicit none
    integer, intent(in) :: loc, block_size, nproc_axis, myproc_axis
    integer, intent(out) :: ig
    integer :: local_block, inblock, global_block

    if (loc <= 0 .or. block_size <= 0 .or. nproc_axis <= 0) then
      ig = 0
      return
    end if
    local_block = (loc - 1) / block_size
    inblock = mod(loc - 1, block_size) + 1
    global_block = local_block * nproc_axis + myproc_axis
    ig = global_block * block_size + inblock
  end subroutine dg_scalapack_global_index_1d

  subroutine get_initial_state_spin_occ_info(dg_frag, system, ispin, occ_weight, nocc, state_cap)
    use structures
    use salmon_global, only: nelec, nelec_spin
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    integer, intent(in) :: ispin
    real(8), intent(out) :: occ_weight(:)
    integer, intent(out) :: nocc
    integer, intent(in), optional :: state_cap

    integer :: cap, io, nocc_eff

    occ_weight(:) = 0.0d0
    cap = min(size(occ_weight), dg_frag%nstate_tot)
    if (present(state_cap)) cap = min(cap, max(0, state_cap))

    if (cap <= 0) then
      nocc = 0
      return
    end if

    if (allocated(system%rocc)) then
      if (ispin <= size(system%rocc, 3)) then
        occ_weight(1:min(cap, size(system%rocc, 1))) = &
          max(0.0d0, system%rocc(1:min(cap, size(system%rocc, 1)), 1, ispin))
      end if
    else
      if (dg_frag%nspin == 1) then
        nocc_eff = min(cap, int(nelec / 2.0d0 + 1.0d-12))
        if (nocc_eff > 0) occ_weight(1:nocc_eff) = 2.0d0
      else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
        nocc_eff = min(cap, nelec_spin(ispin))
        if (nocc_eff > 0) occ_weight(1:nocc_eff) = 1.0d0
      else
        nocc_eff = min(cap, int(nelec / 2.0d0 + 1.0d-12))
        if (nocc_eff > 0) occ_weight(1:nocc_eff) = 1.0d0
      end if
    end if

    nocc = 0
    do io = 1, cap
      if (occ_weight(io) > 0.0d0) nocc = io
    end do
  end subroutine get_initial_state_spin_occ_info

  subroutine measure_fragment_initial_surface_residual(dg_frag, full_rel, core_rel, surface_rel)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    real(8), intent(out) :: full_rel, core_rel, surface_rel

    integer :: ispin, ncheck, io
    real(8) :: sums_local(6), sums_global(6)
    real(8), allocatable :: eps(:), eps_core(:)
    real(8), allocatable :: num_local(:), den_local(:), num_global(:), den_global(:)
    real(8), allocatable :: num_core_local(:), num_core_global(:)
    complex(8), allocatable :: coef_chk(:, :), hc(:, :), hcorec(:, :), sc(:, :), res(:, :)

    full_rel = huge(1.0d0)
    core_rel = huge(1.0d0)
    surface_rel = huge(1.0d0)
    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%H_mat_core_blocks)) return
    if (.not. allocated(dg_frag%S_mat_blocks)) return
    if (.not. allocated(dg_frag%coef)) return
    if (.not. allocated(dg_frag%esp)) return

    ispin = 1
    ncheck = min(dg_frag%nstate_tot, size(dg_frag%coef, 2), size(dg_frag%esp, 1))
    if (allocated(dg_frag%nocc_spin)) then
      if (ispin <= size(dg_frag%nocc_spin)) ncheck = min(ncheck, max(0, dg_frag%nocc_spin(ispin)))
    end if
    if (ncheck <= 0) return

    allocate(coef_chk(size(dg_frag%coef, 1), ncheck))
    allocate(hc(size(dg_frag%coef, 1), ncheck))
    allocate(hcorec(size(dg_frag%coef, 1), ncheck))
    allocate(sc(size(dg_frag%coef, 1), ncheck))
    allocate(res(size(dg_frag%coef, 1), ncheck))
    allocate(eps(ncheck), eps_core(ncheck))
    allocate(num_local(ncheck), den_local(ncheck), num_global(ncheck), den_global(ncheck))
    allocate(num_core_local(ncheck), num_core_global(ncheck))

    coef_chk(:, :) = dg_frag%coef(:, 1:ncheck, ispin)
    hc(:, :) = (0.0d0, 0.0d0)
    hcorec(:, :) = (0.0d0, 0.0d0)
    sc(:, :) = (0.0d0, 0.0d0)

    if (allocated(dg_frag%H_local_block_ids)) then
      call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_chk, hc, dg_frag%H_local_block_ids)
      call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_core_blocks, ispin, coef_chk, hcorec, dg_frag%H_local_block_ids)
    else
      call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, coef_chk, hc)
      call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_core_blocks, ispin, coef_chk, hcorec)
    end if
    if (allocated(dg_frag%H_nl_blocks)) then
      if (allocated(dg_frag%H_nl_local_block_ids)) then
        call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, coef_chk, hc, &
                                               dg_frag%H_nl_local_block_ids)
        call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, coef_chk, hcorec, &
                                               dg_frag%H_nl_local_block_ids)
      else
        call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, coef_chk, hc)
        call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, coef_chk, hcorec)
      end if
    end if
    call apply_overlap_operator_batch_orbital_fragment_self(dg_frag, ispin, coef_chk, sc, .true.)

    do io = 1, ncheck
      num_local(io) = real(sum(conjg(coef_chk(:, io)) * hc(:, io)), kind=8)
      num_core_local(io) = real(sum(conjg(coef_chk(:, io)) * hcorec(:, io)), kind=8)
      den_local(io) = real(sum(conjg(coef_chk(:, io)) * sc(:, io)), kind=8)
    end do
    call comm_summation(num_local, num_global, ncheck, dg_frag%icomm)
    call comm_summation(num_core_local, num_core_global, ncheck, dg_frag%icomm)
    call comm_summation(den_local, den_global, ncheck, dg_frag%icomm)
    do io = 1, ncheck
      eps(io) = num_global(io) / max(den_global(io), 1.0d-300)
      eps_core(io) = num_core_global(io) / max(den_global(io), 1.0d-300)
    end do

    res(:, :) = hc(:, :)
    do io = 1, ncheck
      res(:, io) = res(:, io) - eps(io) * sc(:, io)
    end do
    sums_local(1) = sum(abs(res(:, :))**2)
    sums_local(2) = sum(abs(hc(:, :))**2)
    sums_local(3) = sum(abs(sc(:, :))**2)
    res(:, :) = hcorec(:, :)
    do io = 1, ncheck
      res(:, io) = res(:, io) - eps_core(io) * sc(:, io)
    end do
    sums_local(4) = sum(abs(res(:, :))**2)
    sums_local(5) = sum(abs(hc(:, :) - hcorec(:, :))**2)
    sums_local(6) = sum(abs(hcorec(:, :))**2)

    call comm_summation(sums_local, sums_global, 6, dg_frag%icomm)
    full_rel = sqrt(max(0.0d0, sums_global(1)) / max(sums_global(2), 1.0d-300))
    core_rel = sqrt(max(0.0d0, sums_global(4)) / max(sums_global(6), 1.0d-300))
    surface_rel = sqrt(max(0.0d0, sums_global(5)) / max(sums_global(2), 1.0d-300))

    deallocate(coef_chk, hc, hcorec, sc, res, eps, eps_core, num_local, den_local, num_global, den_global, &
               num_core_local, num_core_global)
  end subroutine measure_fragment_initial_surface_residual

  subroutine diagonalize_initial_dg_full_distributed(dg_frag, system, did_solve)
    use structures
    use communication, only: comm_is_root, comm_summation, comm_get_min, comm_get_max
    use salmon_global, only: yn_scalapack, yn_eigenexa
#ifdef USE_EIGENEXA
    use eigen_libs_mod
#endif
#ifdef USE_SCALAPACK
    use mpi, only: MPI_COMM_WORLD, MPI_UNDEFINED, MPI_INTEGER, MPI_DOUBLE_PRECISION, &
                   MPI_Comm_group, MPI_Group_translate_ranks, MPI_Group_free, MPI_Alltoall, MPI_Alltoallv
#endif
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    logical, intent(out) :: did_solve

    character(len=16) :: env_value_main
    integer :: env_len_main, env_status_main
    integer :: full_eig_max_n, env_read_status
    logical :: use_eigenexa_main

#ifdef USE_SCALAPACK
    integer :: ispin, n, nkeep, nloc_row, nloc_col, nloc_vec_col, i, j
    integer :: nprow, npcol, myrow, mycol, ictxt, ierr, iam, nprocs
    integer :: proc_row, proc_col, loc_row, loc_col
    integer :: mb, nb, nmap, idx, group_src, group_world
    integer :: lwork, liwork, m_found, nz_found
    integer :: n_ortho
    integer :: desca(9), descb(9)
    integer, allocatable :: gridmap_world(:, :), rank_src(:), rank_world(:)
    integer, allocatable :: iwork(:), ifail(:), iclustr(:)
    integer :: iwork_query(1)
    real(8), allocatable :: h_div(:, :), s_div(:, :), eval(:), work(:)
    real(8), allocatable :: y_div(:, :)
    real(8), allocatable :: gap(:)
    real(8) :: work_query(3)
    real(8) :: scalapack_alpha, abstol, vl, vu, orfac, scale_chol
    real(8) :: eig_cut, eig_gap, eig_gap_tol
    real(8) :: s_diag_err, s_offdiag_max, s_frob_err, s_ortho_tol
    real(8) :: scalapack_est_gb_per_rank, scalapack_max_est_gb_per_rank, scalapack_max_gb_per_rank
    real(8) :: scalapack_dense_bytes, scalapack_coef_bytes
    complex(8), allocatable :: coef_part(:, :), coef_sum(:, :)
    integer :: nlocal
    integer :: NUMROC
#endif
#ifdef USE_EIGENEXA
    integer, parameter :: panel_ex = 32
    integer :: ispin_ex, n_ex, nkeep_ex, nx_ex, ny_ex, nlocal_ex
    integer :: nnod_ex, x_nnod_ex, y_nnod_ex, inod_ex, x_inod_ex, y_inod_ex
    integer :: ix_s_ex, ix_e_ex, iy_s_ex, iy_e_ex, ix_loc_ex, iy_loc_ex
    integer :: ig_ex, jg_ex, idx_ex, ib_ex, ie_ex, nb_ex, ip_ex
    integer :: iloc_ex, state_ex, col_ex, ie_col_ex, nb_col_ex, gid_ex
    integer :: ix_panel_s_ex, ix_panel_e_ex, nb_left_ex
    integer :: ifrag_s_ex, nbf_s_ex, iblk_s_ex, nbf_max_ex
    integer, allocatable :: ix_panel_gid_ex(:)
    integer, allocatable :: ifrag_of_gid_ex(:), io_of_gid_ex(:)
    real(8), allocatable :: s_div_ex(:, :), u_div_ex(:, :), lambda_ex(:)
    real(8), allocatable :: s_block_all_ex(:, :, :)
    real(8), allocatable :: horth_div_ex(:, :), y_div_ex(:, :)
    real(8), allocatable :: w_left_panel_ex(:, :), w_panel_ex(:, :), hw_part_ex(:, :), hw_panel_ex(:, :)
    real(8), allocatable :: y_panel_ex(:, :)
    real(8), allocatable :: eval_ex(:)
    real(8) :: s_absmax_local_ex, s_absmax_global_ex, s_val_abs_ex
    real(8) :: s_diag_min_local_ex, s_diag_min_global_ex
    real(8) :: s_diag_max_local_ex, s_diag_max_global_ex
    real(8) :: s_diag_count_local_ex, s_diag_count_global_ex
    real(8) :: s_nonzero_local_ex, s_nonzero_global_ex
    real(8) :: s_asym_local_ex, s_asym_global_ex
    real(8) :: eig_min_ex, eig_max_ex
    integer :: missing_gid_ex, bad_diag_ex
    complex(8), allocatable :: x_panel_ex(:, :), hx_panel_ex(:, :)
#endif

    did_solve = .false.
    use_eigenexa_main = (yn_eigenexa == 'y')
    env_value_main = ''
    call get_environment_variable('SALMON_DG_FULL_EIGENEXA', env_value_main, &
                                  length=env_len_main, status=env_status_main)
    if (env_status_main == 0 .and. env_len_main > 0) then
      select case (env_value_main(1:1))
      case ('0', 'n', 'N', 'f', 'F')
        use_eigenexa_main = .false.
      case ('1', 'y', 'Y', 't', 'T')
        use_eigenexa_main = .true.
      end select
    end if
    full_eig_max_n = 32768
    env_value_main = ''
    call get_environment_variable('SALMON_DG_FULL_EIG_MAX_N', env_value_main, &
                                  length=env_len_main, status=env_status_main)
    if (env_status_main == 0 .and. env_len_main > 0) then
      read(env_value_main(1:env_len_main), *, iostat=env_read_status) full_eig_max_n
      if (env_read_status /= 0) full_eig_max_n = 32768
    end if
    if (use_eigenexa_main) then
#ifdef USE_EIGENEXA
      if (.not. allocated(dg_frag%H_mat_blocks)) return
      if (.not. allocated(dg_frag%S_mat_blocks)) return
      if (.not. allocated(dg_frag%coef)) return
      if (.not. allocated(dg_frag%esp)) return
      if (dg_frag%nspin /= 1) return

      ispin_ex = 1
      n_ex = dg_frag%n_mat(ispin_ex)
      if (n_ex <= 0) return
      if (n_ex > full_eig_max_n) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,i0)') '[WARN] DG-DIST-EIG full dense DG eigensolver skipped: n=', &
            n_ex, ' max_n=', full_eig_max_n
          write(*,'(1x,a)') '[DG-DIST-EIG] trying ScaLAPACK distributed generalized eigensolver instead.'
          flush(6)
        end if
        goto 900
      end if
      nkeep_ex = min(dg_frag%nstate_tot, size(dg_frag%coef, 2), size(dg_frag%esp, 1), n_ex)
      if (nkeep_ex <= 0) return

      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,i0)') '[DG-DIST-EIG] EigenExa full DG generalized solve setup: n=', &
          n_ex, ' nkeep=', nkeep_ex
        flush(6)
      end if

      allocate(ifrag_of_gid_ex(n_ex), io_of_gid_ex(n_ex))
      ifrag_of_gid_ex(:) = 0
      io_of_gid_ex(:) = 0
      do ig_ex = 1, dg_frag%n_frag
        do jg_ex = 1, dg_frag%n_basis(ig_ex, ispin_ex)
          idx_ex = dg_frag%index_basis(jg_ex, ig_ex, ispin_ex)
          if (idx_ex >= 1 .and. idx_ex <= n_ex) then
            ifrag_of_gid_ex(idx_ex) = ig_ex
            io_of_gid_ex(idx_ex) = jg_ex
          end if
        end do
      end do
      missing_gid_ex = 0
      do ig_ex = 1, n_ex
        if (ifrag_of_gid_ex(ig_ex) <= 0 .or. io_of_gid_ex(ig_ex) <= 0) &
          missing_gid_ex = missing_gid_ex + 1
      end do
      call comm_get_max(missing_gid_ex, dg_frag%icomm)
      if (missing_gid_ex > 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0)') '[WARN] DG-DIST-EIG full DG index_basis map is incomplete: missing_gid_max=', &
            missing_gid_ex
          flush(6)
        end if
        deallocate(ifrag_of_gid_ex, io_of_gid_ex)
        return
      end if

      nbf_max_ex = maxval(dg_frag%n_basis(:, ispin_ex))
      allocate(s_block_all_ex(max(1, nbf_max_ex), max(1, nbf_max_ex), dg_frag%n_frag))
      s_block_all_ex(:, :, :) = 0.0d0
      if (allocated(dg_frag%S_block_map) .and. allocated(dg_frag%S_mat_blocks)) then
        do ifrag_s_ex = 1, dg_frag%n_frag
          if (dg_frag%ifrag_group /= ifrag_s_ex) cycle
          if (.not. dg_frag%is_frag_root) cycle
          iblk_s_ex = dg_frag%S_block_map(ifrag_s_ex, ifrag_s_ex)
          if (iblk_s_ex <= 0 .or. iblk_s_ex > size(dg_frag%S_mat_blocks)) cycle
          if (.not. allocated(dg_frag%S_mat_blocks(iblk_s_ex)%val)) cycle
          nbf_s_ex = min(dg_frag%n_basis(ifrag_s_ex, ispin_ex), &
                         size(dg_frag%S_mat_blocks(iblk_s_ex)%val, 1), &
                         size(dg_frag%S_mat_blocks(iblk_s_ex)%val, 2), nbf_max_ex)
          if (nbf_s_ex <= 0) cycle
          s_block_all_ex(1:nbf_s_ex, 1:nbf_s_ex, ifrag_s_ex) = &
            dg_frag%S_mat_blocks(iblk_s_ex)%val(1:nbf_s_ex, 1:nbf_s_ex, ispin_ex)
        end do
      end if
      call comm_summation(s_block_all_ex, dg_frag%icomm)

      call eigen_init(dg_frag%icomm)
      call eigen_get_matdims(n_ex, nx_ex, ny_ex)
      call eigen_get_procs(nnod_ex, x_nnod_ex, y_nnod_ex)
      call eigen_get_id(inod_ex, x_inod_ex, y_inod_ex)
      ix_s_ex = eigen_loop_start(1, x_nnod_ex, x_inod_ex)
      ix_e_ex = eigen_loop_end(n_ex, x_nnod_ex, x_inod_ex)
      iy_s_ex = eigen_loop_start(1, y_nnod_ex, y_inod_ex)
      iy_e_ex = eigen_loop_end(n_ex, y_nnod_ex, y_inod_ex)
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,8(a,i0))') '[DG-DIST-EIG-MAT] EigenExa layout', &
          ' nx=', nx_ex, ' ny=', ny_ex, ' nnod=', nnod_ex, ' x_nnod=', x_nnod_ex, &
          ' y_nnod=', y_nnod_ex, ' ix_s=', ix_s_ex, ' ix_e=', ix_e_ex, ' iy_s=', iy_s_ex
        write(*,'(1x,a,i0)') '[DG-DIST-EIG-MAT] EigenExa local y end: iy_e=', iy_e_ex
        flush(6)
      end if

      allocate(s_div_ex(nx_ex, ny_ex), u_div_ex(nx_ex, ny_ex), lambda_ex(n_ex))
      s_div_ex(:, :) = 0.0d0
      s_absmax_local_ex = 0.0d0
      s_diag_min_local_ex = huge(1.0d0)
      s_diag_max_local_ex = -huge(1.0d0)
      s_diag_count_local_ex = 0.0d0
      s_nonzero_local_ex = 0.0d0
      s_asym_local_ex = 0.0d0
      do iy_loc_ex = iy_s_ex, iy_e_ex
        jg_ex = eigen_translate_l2g(iy_loc_ex, y_nnod_ex, y_inod_ex)
        if (jg_ex > n_ex) cycle
        do ix_loc_ex = ix_s_ex, ix_e_ex
          ig_ex = eigen_translate_l2g(ix_loc_ex, x_nnod_ex, x_inod_ex)
          if (ig_ex > n_ex) cycle
          s_div_ex(ix_loc_ex, iy_loc_ex) = dg_full_s_element(ig_ex, jg_ex, ispin_ex, &
                                                            ifrag_of_gid_ex, io_of_gid_ex)
          s_val_abs_ex = abs(s_div_ex(ix_loc_ex, iy_loc_ex))
          s_absmax_local_ex = max(s_absmax_local_ex, s_val_abs_ex)
          s_asym_local_ex = max(s_asym_local_ex, &
            abs(s_div_ex(ix_loc_ex, iy_loc_ex) - dg_full_s_element(jg_ex, ig_ex, ispin_ex, &
                                                                   ifrag_of_gid_ex, io_of_gid_ex)))
          if (s_val_abs_ex > 0.0d0) s_nonzero_local_ex = s_nonzero_local_ex + 1.0d0
          if (ig_ex == jg_ex) then
            s_diag_count_local_ex = s_diag_count_local_ex + 1.0d0
            s_diag_min_local_ex = min(s_diag_min_local_ex, s_div_ex(ix_loc_ex, iy_loc_ex))
            s_diag_max_local_ex = max(s_diag_max_local_ex, s_div_ex(ix_loc_ex, iy_loc_ex))
          end if
        end do
      end do
      s_absmax_global_ex = -s_absmax_local_ex
      call comm_get_min(s_absmax_global_ex, dg_frag%icomm)
      s_absmax_global_ex = -s_absmax_global_ex
      s_diag_min_global_ex = s_diag_min_local_ex
      call comm_get_min(s_diag_min_global_ex, dg_frag%icomm)
      s_diag_max_global_ex = -s_diag_max_local_ex
      call comm_get_min(s_diag_max_global_ex, dg_frag%icomm)
      s_diag_max_global_ex = -s_diag_max_global_ex
      s_asym_global_ex = -s_asym_local_ex
      call comm_get_min(s_asym_global_ex, dg_frag%icomm)
      s_asym_global_ex = -s_asym_global_ex
      call comm_summation(s_diag_count_local_ex, s_diag_count_global_ex, dg_frag%icomm)
      call comm_summation(s_nonzero_local_ex, s_nonzero_global_ex, dg_frag%icomm)
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,6(a,1pe13.5))') '[DG-DIST-EIG-MAT] S stats before EigenExa', &
          ' absmax=', s_absmax_global_ex, ' diag_min=', s_diag_min_global_ex, &
          ' diag_max=', s_diag_max_global_ex, ' diag_count=', s_diag_count_global_ex, &
          ' nonzero_count=', s_nonzero_global_ex, ' asym_max=', s_asym_global_ex
        flush(6)
      end if
      bad_diag_ex = 0
      if (s_diag_count_global_ex < dble(n_ex)) bad_diag_ex = 1
      if (s_diag_min_global_ex <= 0.0d0) bad_diag_ex = 1
      if (s_absmax_global_ex <= 0.0d0) bad_diag_ex = 1
      if (s_asym_global_ex > 1.0d-8 * max(1.0d0, s_absmax_global_ex)) bad_diag_ex = 1
      if (bad_diag_ex /= 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a)') '[WARN] DG-DIST-EIG invalid S matrix setup before EigenExa; using fallback.'
          flush(6)
        end if
        call eigen_free()
        deallocate(s_div_ex, u_div_ex, lambda_ex, ifrag_of_gid_ex, io_of_gid_ex, s_block_all_ex)
        return
      end if

      lambda_ex(:) = huge(1.0d0)
      call eigen_sx(n_ex, n_ex, s_div_ex, nx_ex, lambda_ex, u_div_ex, nx_ex)
      eig_min_ex = minval(lambda_ex(1:n_ex))
      eig_max_ex = maxval(lambda_ex(1:n_ex))
      if (abs(eig_min_ex) >= 1.0d100 .or. abs(eig_max_ex) >= 1.0d100) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,2(a,1pe13.5))') '[WARN] DG-DIST-EIG EigenExa S eigenvalues look unset; using fallback.', &
            ' lambda_min=', eig_min_ex, ' lambda_max=', eig_max_ex
          flush(6)
        end if
        call eigen_free()
        deallocate(s_div_ex, u_div_ex, lambda_ex, ifrag_of_gid_ex, io_of_gid_ex, s_block_all_ex)
        return
      end if
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,2(a,1pe13.5))') '[DG-DIST-EIG] EigenExa S diagonalization done', &
          ' lambda_min=', eig_min_ex, ' lambda_max=', eig_max_ex
        flush(6)
      end if
      if (eig_min_ex <= 1.0d-12) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,1pe13.5)') '[WARN] DG-DIST-EIG S is nearly singular; lambda_min=', &
            eig_min_ex
          flush(6)
        end if
        call eigen_free()
        deallocate(s_div_ex, u_div_ex, lambda_ex, ifrag_of_gid_ex, io_of_gid_ex, s_block_all_ex)
        return
      end if

      allocate(horth_div_ex(nx_ex, ny_ex), y_div_ex(nx_ex, ny_ex), eval_ex(n_ex))
      horth_div_ex(:, :) = 0.0d0
      nlocal_ex = size(dg_frag%coef, 1)
      allocate(x_panel_ex(nlocal_ex, panel_ex), hx_panel_ex(nlocal_ex, panel_ex))
      allocate(ix_panel_gid_ex(panel_ex), w_left_panel_ex(n_ex, panel_ex))
      allocate(w_panel_ex(n_ex, panel_ex), hw_part_ex(n_ex, panel_ex), hw_panel_ex(n_ex, panel_ex))
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,i0)') '[DG-DIST-EIG] building W^T H W by block-sparse H panels: n=', &
          n_ex, ' panel=', panel_ex
        flush(6)
      end if

      do ib_ex = 1, n_ex, panel_ex
        ie_ex = min(n_ex, ib_ex + panel_ex - 1)
        nb_ex = ie_ex - ib_ex + 1
        call gather_w_range(ib_ex, ie_ex, w_panel_ex)
        x_panel_ex(:, 1:nb_ex) = (0.0d0, 0.0d0)
        do iloc_ex = 1, nlocal_ex
          if (iloc_ex > size(dg_frag%local_coef_global_ids, 1)) cycle
          gid_ex = dg_frag%local_coef_global_ids(iloc_ex, ispin_ex)
          if (gid_ex < 1 .or. gid_ex > n_ex) cycle
          x_panel_ex(iloc_ex, 1:nb_ex) = cmplx(w_panel_ex(gid_ex, 1:nb_ex), 0.0d0, kind=8)
        end do
        hx_panel_ex(:, 1:nb_ex) = (0.0d0, 0.0d0)
        call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin_ex, &
                                       x_panel_ex(:, 1:nb_ex), hx_panel_ex(:, 1:nb_ex))
        if (allocated(dg_frag%H_nl_blocks)) then
          call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin_ex, &
                                                 x_panel_ex(:, 1:nb_ex), hx_panel_ex(:, 1:nb_ex))
        end if

        hw_part_ex(:, 1:nb_ex) = 0.0d0
        do iloc_ex = 1, nlocal_ex
          if (iloc_ex > size(dg_frag%local_coef_global_ids, 1)) cycle
          gid_ex = dg_frag%local_coef_global_ids(iloc_ex, ispin_ex)
          if (gid_ex < 1 .or. gid_ex > n_ex) cycle
          hw_part_ex(gid_ex, 1:nb_ex) = real(hx_panel_ex(iloc_ex, 1:nb_ex), kind=8)
        end do
        call comm_summation(hw_part_ex(1:n_ex, 1:nb_ex), hw_panel_ex(1:n_ex, 1:nb_ex), n_ex * nb_ex, dg_frag%icomm)

        do ix_panel_s_ex = ix_s_ex, ix_e_ex, panel_ex
          ix_panel_e_ex = min(ix_e_ex, ix_panel_s_ex + panel_ex - 1)
          nb_left_ex = ix_panel_e_ex - ix_panel_s_ex + 1
          do ip_ex = 1, nb_left_ex
            ix_panel_gid_ex(ip_ex) = eigen_translate_l2g(ix_panel_s_ex + ip_ex - 1, x_nnod_ex, x_inod_ex)
          end do
          call gather_w_columns(ix_panel_gid_ex(1:nb_left_ex), nb_left_ex, w_left_panel_ex)
          do iy_loc_ex = iy_s_ex, iy_e_ex
            jg_ex = eigen_translate_l2g(iy_loc_ex, y_nnod_ex, y_inod_ex)
            if (jg_ex < ib_ex .or. jg_ex > ie_ex) cycle
            do ix_loc_ex = ix_panel_s_ex, ix_panel_e_ex
              ig_ex = eigen_translate_l2g(ix_loc_ex, x_nnod_ex, x_inod_ex)
              if (ig_ex > n_ex) cycle
              ip_ex = ix_loc_ex - ix_panel_s_ex + 1
              horth_div_ex(ix_loc_ex, iy_loc_ex) = &
                sum(w_left_panel_ex(1:n_ex, ip_ex) * hw_panel_ex(1:n_ex, jg_ex - ib_ex + 1))
            end do
          end do
        end do
      end do

      s_absmax_local_ex = 0.0d0
      s_diag_min_local_ex = huge(1.0d0)
      s_diag_max_local_ex = -huge(1.0d0)
      s_diag_count_local_ex = 0.0d0
      s_nonzero_local_ex = 0.0d0
      do iy_loc_ex = iy_s_ex, iy_e_ex
        jg_ex = eigen_translate_l2g(iy_loc_ex, y_nnod_ex, y_inod_ex)
        if (jg_ex > n_ex) cycle
        do ix_loc_ex = ix_s_ex, ix_e_ex
          ig_ex = eigen_translate_l2g(ix_loc_ex, x_nnod_ex, x_inod_ex)
          if (ig_ex > n_ex) cycle
          s_val_abs_ex = abs(horth_div_ex(ix_loc_ex, iy_loc_ex))
          s_absmax_local_ex = max(s_absmax_local_ex, s_val_abs_ex)
          if (s_val_abs_ex > 0.0d0) s_nonzero_local_ex = s_nonzero_local_ex + 1.0d0
          if (ig_ex == jg_ex) then
            s_diag_count_local_ex = s_diag_count_local_ex + 1.0d0
            s_diag_min_local_ex = min(s_diag_min_local_ex, horth_div_ex(ix_loc_ex, iy_loc_ex))
            s_diag_max_local_ex = max(s_diag_max_local_ex, horth_div_ex(ix_loc_ex, iy_loc_ex))
          end if
        end do
      end do
      s_absmax_global_ex = -s_absmax_local_ex
      call comm_get_min(s_absmax_global_ex, dg_frag%icomm)
      s_absmax_global_ex = -s_absmax_global_ex
      s_diag_min_global_ex = s_diag_min_local_ex
      call comm_get_min(s_diag_min_global_ex, dg_frag%icomm)
      s_diag_max_global_ex = -s_diag_max_local_ex
      call comm_get_min(s_diag_max_global_ex, dg_frag%icomm)
      s_diag_max_global_ex = -s_diag_max_global_ex
      call comm_summation(s_diag_count_local_ex, s_diag_count_global_ex, dg_frag%icomm)
      call comm_summation(s_nonzero_local_ex, s_nonzero_global_ex, dg_frag%icomm)
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,5(a,1pe13.5))') '[DG-DIST-EIG-MAT] Horth stats before EigenExa', &
          ' absmax=', s_absmax_global_ex, ' diag_min=', s_diag_min_global_ex, &
          ' diag_max=', s_diag_max_global_ex, ' diag_count=', s_diag_count_global_ex, &
          ' nonzero_count=', s_nonzero_global_ex
        flush(6)
      end if
      if (s_diag_count_global_ex < dble(n_ex) .or. s_absmax_global_ex <= 0.0d0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a)') '[WARN] DG-DIST-EIG invalid Horth matrix setup before EigenExa; using fallback.'
          flush(6)
        end if
        call eigen_free()
        deallocate(s_div_ex, u_div_ex, lambda_ex, ifrag_of_gid_ex, io_of_gid_ex, ix_panel_gid_ex, s_block_all_ex)
        deallocate(horth_div_ex, y_div_ex, eval_ex, w_left_panel_ex, w_panel_ex, hw_part_ex, hw_panel_ex)
        deallocate(x_panel_ex, hx_panel_ex)
        return
      end if

      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[DG-DIST-EIG] diagonalizing orthonormalized H with EigenExa'
        flush(6)
      end if
      eval_ex(:) = huge(1.0d0)
      call eigen_sx(n_ex, n_ex, horth_div_ex, nx_ex, eval_ex, y_div_ex, nx_ex)
      eig_min_ex = minval(eval_ex(1:n_ex))
      eig_max_ex = maxval(eval_ex(1:n_ex))
      if (abs(eig_min_ex) >= 1.0d100 .or. abs(eig_max_ex) >= 1.0d100) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,2(a,1pe13.5))') '[WARN] DG-DIST-EIG EigenExa H eigenvalues look unset; using fallback.', &
            ' eig_min=', eig_min_ex, ' eig_max=', eig_max_ex
          flush(6)
        end if
        call eigen_free()
        deallocate(s_div_ex, u_div_ex, lambda_ex, ifrag_of_gid_ex, io_of_gid_ex, ix_panel_gid_ex, s_block_all_ex)
        deallocate(horth_div_ex, y_div_ex, eval_ex, w_left_panel_ex, w_panel_ex, hw_part_ex, hw_panel_ex)
        deallocate(x_panel_ex, hx_panel_ex)
        return
      end if
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,2(a,1pe13.5))') '[DG-DIST-EIG] EigenExa H diagonalization done', &
          ' eig_min=', eig_min_ex, ' eig_max=', eig_max_ex
        flush(6)
      end if

      allocate(y_panel_ex(n_ex, panel_ex))
      dg_frag%coef(:, 1:nkeep_ex, ispin_ex) = (0.0d0, 0.0d0)
      do ib_ex = 1, nkeep_ex, panel_ex
        ie_ex = min(nkeep_ex, ib_ex + panel_ex - 1)
        nb_ex = ie_ex - ib_ex + 1
        call gather_y_range(ib_ex, ie_ex, y_panel_ex)
        do col_ex = 1, n_ex, panel_ex
          ie_col_ex = min(n_ex, col_ex + panel_ex - 1)
          nb_col_ex = ie_col_ex - col_ex + 1
          call gather_w_range(col_ex, ie_col_ex, w_panel_ex)
          do state_ex = 1, nb_ex
            do iloc_ex = 1, nlocal_ex
              if (iloc_ex > size(dg_frag%local_coef_global_ids, 1)) cycle
              gid_ex = dg_frag%local_coef_global_ids(iloc_ex, ispin_ex)
              if (gid_ex < 1 .or. gid_ex > n_ex) cycle
              dg_frag%coef(iloc_ex, ib_ex + state_ex - 1, ispin_ex) = &
                dg_frag%coef(iloc_ex, ib_ex + state_ex - 1, ispin_ex) + &
                cmplx(sum(w_panel_ex(gid_ex, 1:nb_col_ex) * &
                          y_panel_ex(col_ex:ie_col_ex, state_ex)), 0.0d0, kind=8)
            end do
          end do
        end do
      end do
      if (allocated(dg_frag%coef_work)) dg_frag%coef_work(:, 1:nkeep_ex, ispin_ex) = dg_frag%coef(:, 1:nkeep_ex, ispin_ex)
      if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, 1:nkeep_ex, ispin_ex) = dg_frag%coef(:, 1:nkeep_ex, ispin_ex)
      dg_frag%esp(1:nkeep_ex, ispin_ex) = eval_ex(1:nkeep_ex)
      did_solve = .true.
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,2(a,1pe13.5))') '[DG-DIST-EIG] EigenExa generalized solve done', &
          ' eig_min=', dg_frag%esp(1, ispin_ex), ' eig_keep=', dg_frag%esp(nkeep_ex, ispin_ex)
        flush(6)
      end if
      call eigen_free()
      deallocate(s_div_ex, u_div_ex, lambda_ex, ifrag_of_gid_ex, io_of_gid_ex, ix_panel_gid_ex, s_block_all_ex)
      deallocate(horth_div_ex, y_div_ex, eval_ex, w_left_panel_ex, w_panel_ex, hw_part_ex, hw_panel_ex)
      deallocate(x_panel_ex, hx_panel_ex, y_panel_ex)
      return
#else
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[DG-DIST-EIG] EigenExa full DG path requested/default but this binary was built without USE_EIGENEXA;'
        write(*,'(1x,a)') '[DG-DIST-EIG] trying ScaLAPACK distributed generalized eigensolver instead.'
        flush(6)
      end if
      goto 900
#endif
    end if
900 continue

#ifndef USE_SCALAPACK
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') '[DG-DIST-EIG] this binary was built without USE_SCALAPACK; using block-sparse relaxation.'
      flush(6)
    end if
    return
#else
    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%S_mat_blocks)) return
    if (.not. allocated(dg_frag%coef)) return
    if (.not. allocated(dg_frag%esp)) return
    if (dg_frag%nspin /= 1) return

    ispin = 1
    n = dg_frag%n_mat(ispin)
    if (n <= 0) return
    nkeep = min(dg_frag%nstate_tot, size(dg_frag%coef, 2), size(dg_frag%esp, 1), n)
    if (nkeep <= 0) return

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0)') '[DG-DIST-EIG] ScaLAPACK PDPOTRF/PDSYGST/PDSYEVX partial DG generalized solve: n=', n, &
        ' nkeep=', nkeep
      flush(6)
    end if

    npcol = int(sqrt(dble(max(1, dg_frag%isize))))
    do
      if (npcol <= 1) exit
      if (mod(dg_frag%isize, npcol) == 0) exit
      npcol = npcol - 1
    end do
    nprow = max(1, dg_frag%isize / max(1, npcol))
    if (nprow > npcol) then
      i = nprow
      nprow = npcol
      npcol = i
    end if
    if (nprow * npcol /= dg_frag%isize) stop 'DG-DIST-EIG: invalid ScaLAPACK process grid'
    ! A 1x1 block-cyclic layout makes ScaLAPACK factorizations communication
    ! dominated at n ~ 1e5.  Use square tiles so PDPOTRF/PDSYGST/PDSYEVX can
    ! run mostly through BLAS3 kernels.
    mb = 64
    nb = 64

    call BLACS_PINFO(iam, nprocs)
    if (nprocs < 1) call BLACS_SETUP(iam, nprocs)
    call BLACS_GET(0, 0, ictxt)

    nmap = nprow * npcol
    allocate(rank_src(nmap), rank_world(nmap), gridmap_world(nprow, npcol))
    idx = 0
    do j = 1, npcol
      do i = 1, nprow
        idx = idx + 1
        rank_src(idx) = idx - 1
      end do
    end do
    call MPI_Comm_group(dg_frag%icomm, group_src, ierr)
    call MPI_Comm_group(MPI_COMM_WORLD, group_world, ierr)
    call MPI_Group_translate_ranks(group_src, nmap, rank_src, group_world, rank_world, ierr)
    call MPI_Group_free(group_world, ierr)
    call MPI_Group_free(group_src, ierr)
    idx = 0
    do j = 1, npcol
      do i = 1, nprow
        idx = idx + 1
        if (rank_world(idx) == MPI_UNDEFINED) stop 'DG-DIST-EIG: rank translation failed'
        gridmap_world(i, j) = rank_world(idx)
      end do
    end do
    deallocate(rank_src, rank_world)

    call BLACS_GRIDMAP(ictxt, gridmap_world, nprow, nprow, npcol)
    deallocate(gridmap_world)
    call BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)
    nloc_row = NUMROC(n, mb, myrow, 0, nprow)
    nloc_col = NUMROC(n, nb, mycol, 0, npcol)
    call DESCINIT(desca, n, n, mb, nb, 0, 0, ictxt, max(1, nloc_row), ierr)
    ! PDSYEVX documents Z/DESCZ as a global (N,N) matrix.  Even when RANGE
    ! selects only the occupied window, implementations may validate DESCZ
    ! against DESCA before the number of selected vectors is known.
    nloc_vec_col = nloc_col
    call DESCINIT(descb, n, n, mb, nb, 0, 0, ictxt, max(1, nloc_row), ierr)
    scalapack_max_gb_per_rank = 24.0d0
    env_value_main = ''
    call get_environment_variable('SALMON_DG_SCALAPACK_MAX_GB_PER_RANK', env_value_main, &
                                  length=env_len_main, status=env_status_main)
    if (env_status_main == 0 .and. env_len_main > 0) then
      read(env_value_main(1:env_len_main), *, iostat=env_read_status) scalapack_max_gb_per_rank
      if (env_read_status /= 0) scalapack_max_gb_per_rank = 24.0d0
    end if
    ! Main local dense arrays are H, S, Z/Y, plus one dense-equivalent margin
    ! for ScaLAPACK work buffers.  Coefficient gather buffers are separate.
    scalapack_dense_bytes = 4.0d0 * 8.0d0 * dble(max(1, nloc_row)) * dble(max(1, nloc_col))
    scalapack_coef_bytes = 4.0d0 * 8.0d0 * dble(size(dg_frag%coef, 1)) * dble(max(1, nkeep))
    scalapack_est_gb_per_rank = (scalapack_dense_bytes + scalapack_coef_bytes) / 1024.0d0**3
    scalapack_max_est_gb_per_rank = -scalapack_est_gb_per_rank
    call comm_get_min(scalapack_max_est_gb_per_rank, dg_frag%icomm)
    scalapack_max_est_gb_per_rank = -scalapack_max_est_gb_per_rank
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,5(a,i0),2(a,1pe13.5))') '[DG-DIST-EIG-MEM] ScaLAPACK estimate', &
        ' n=', n, ' ranks=', dg_frag%isize, ' nloc_row=', nloc_row, ' nloc_col=', nloc_col, &
        ' nkeep=', nkeep, ' max_GB_per_rank_est=', scalapack_max_est_gb_per_rank, &
        ' max_GB_per_rank=', scalapack_max_gb_per_rank
      flush(6)
    end if
    if (scalapack_max_gb_per_rank > 0.0d0 .and. scalapack_max_est_gb_per_rank > scalapack_max_gb_per_rank) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,2(a,1pe13.5),a)') '[WARN] DG-DIST-EIG skipping ScaLAPACK full dense solve:', &
          ' estimated_GB_per_rank=', scalapack_max_est_gb_per_rank, &
          ' limit_GB_per_rank=', scalapack_max_gb_per_rank, &
          ' (set SALMON_DG_SCALAPACK_MAX_GB_PER_RANK<=0 to force)'
        flush(6)
      end if
      call BLACS_GRIDEXIT(ictxt)
      return
    end if
    allocate(h_div(max(1, nloc_row), max(1, nloc_col)))
    allocate(s_div(max(1, nloc_row), max(1, nloc_col)))
    allocate(eval(n))
    h_div(:, :) = 0.0d0
    s_div(:, :) = 0.0d0
    call assemble_scalapack_hs_from_blocks(h_div, s_div, ispin)

    allocate(y_div(max(1, nloc_row), max(1, nloc_vec_col)))

    call PDPOTRF('L', n, s_div, 1, 1, desca, ierr)
    if (ierr /= 0) then
      if (comm_is_root(dg_frag%id)) write(*,'(1x,a,i0)') '[WARN] DG ScaLAPACK PDPOTRF(S) failed: info=', ierr
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, y_div)
      return
    end if
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') '[DG-DIST-EIG] ScaLAPACK S Cholesky factorization done'
      flush(6)
    end if

    call PDSYGST(1, 'L', n, h_div, 1, 1, desca, s_div, 1, 1, desca, scale_chol, ierr)
    if (ierr /= 0) then
      if (comm_is_root(dg_frag%id)) write(*,'(1x,a,i0)') '[WARN] DG ScaLAPACK PDSYGST(H,S) failed: info=', ierr
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, y_div)
      return
    end if
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe13.5)') '[DG-DIST-EIG] ScaLAPACK H reduced to standard form: scale=', scale_chol
      flush(6)
    end if

    allocate(ifail(n))
    allocate(iclustr(2 * max(1, nprow * npcol)))
    allocate(gap(max(1, nprow * npcol)))
    vl = 0.0d0
    vu = 0.0d0
    ! PDSYEVX checks that scalar inputs are bitwise-consistent across the
    ! BLACS grid by broadcasting rank-0 values through WORK(1:3).  Avoid
    ! per-rank PDLAMCH differences here; use a deterministic positive floor.
    abstol = 1.0d-12
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe13.5)') '[DG-DIST-EIG] ScaLAPACK PDSYEVX abstol=', abstol
      flush(6)
    end if
    orfac = 1.0d-3
    lwork = -1
    liwork = -1
    work_query(:) = 0.0d0
    call PDSYEVX('N', 'A', 'L', n, h_div, 1, 1, desca, vl, vu, 1, n, abstol, &
                 m_found, nz_found, eval, orfac, y_div, 1, 1, descb, work_query, lwork, &
                 iwork_query, liwork, ifail, iclustr, gap, ierr)
    if (ierr /= 0) then
      if (comm_is_root(dg_frag%id)) write(*,'(1x,a,i0)') &
        '[WARN] DG ScaLAPACK PDSYEVX(Hstd eigenvalues) workspace query failed: info=', ierr
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, y_div, ifail, iclustr, gap)
      return
    end if
    lwork = max(3, int(work_query(1)))
    liwork = max(1, iwork_query(1))
    allocate(work(lwork))
    allocate(iwork(liwork))

    y_div(:, :) = 0.0d0
    call PDSYEVX('N', 'A', 'L', n, h_div, 1, 1, desca, vl, vu, 1, n, abstol, &
                 m_found, nz_found, eval, orfac, y_div, 1, 1, descb, work, lwork, &
                 iwork, liwork, ifail, iclustr, gap, ierr)
    if (ierr /= 0 .or. m_found < n) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,3(a,i0))') '[WARN] DG ScaLAPACK PDSYEVX(Hstd eigenvalues) failed/incomplete:', &
          ' info=', ierr, ' m=', m_found, ' nz=', nz_found
        if (ierr /= 0) then
          write(*,'(1x,a,1pe13.5)') '[WARN] DG ScaLAPACK PDSYEVX min_gap=', minval(gap)
          write(*,'(1x,a,4(i0,1x))') '[WARN] DG ScaLAPACK PDSYEVX ifail sample=', &
            ifail(1), ifail(min(2, n)), ifail(min(3, n)), ifail(min(4, n))
        end if
      end if
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, work, iwork, y_div, ifail, iclustr, gap)
      return
    end if
    if (nkeep >= n) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,2(a,i0))') '[WARN] DG ScaLAPACK value-window path requires nkeep < n:', &
          ' nkeep=', nkeep, ' n=', n
      end if
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, work, iwork, y_div, ifail, iclustr, gap)
      return
    end if
    eig_cut = eval(nkeep)
    eig_gap = eval(nkeep + 1) - eval(nkeep)
    eig_gap_tol = max(1.0d-10, 1.0d-8 * max(1.0d0, abs(eig_cut)))
    if (eig_gap <= eig_gap_tol) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,3(a,1pe13.5))') '[WARN] DG ScaLAPACK eigenvalue cutoff is degenerate; cannot use RANGE=V safely:', &
          ' eig_cut=', eig_cut, ' gap=', eig_gap, ' tol=', eig_gap_tol
      end if
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, work, iwork, y_div, ifail, iclustr, gap)
      return
    end if
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,3(a,1pe13.5))') '[DG-DIST-EIG] ScaLAPACK Hstd eigenvalue window selected', &
        ' eig_min=', eval(1), ' eig_keep=', eig_cut, ' gap=', eig_gap
      flush(6)
    end if

    deallocate(work, iwork)
    h_div(:, :) = 0.0d0
    call assemble_scalapack_hs_from_blocks(h_div, s_div, ispin, .false., .false.)
    call PDSYGST(1, 'L', n, h_div, 1, 1, desca, s_div, 1, 1, desca, scale_chol, ierr)
    if (ierr /= 0) then
      if (comm_is_root(dg_frag%id)) write(*,'(1x,a,i0)') '[WARN] DG ScaLAPACK PDSYGST(H,S) retry failed: info=', ierr
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, y_div, ifail, iclustr, gap)
      return
    end if

    vl = eval(1) - max(1.0d0, abs(eval(1))) * 1.0d-8
    vu = eig_cut + 0.5d0 * eig_gap
    lwork = -1
    liwork = -1
    work_query(:) = 0.0d0
    call PDSYEVX('V', 'V', 'L', n, h_div, 1, 1, desca, vl, vu, 1, n, abstol, &
                 m_found, nz_found, eval, orfac, y_div, 1, 1, descb, work_query, lwork, &
                 iwork_query, liwork, ifail, iclustr, gap, ierr)
    if (ierr /= 0) then
      if (comm_is_root(dg_frag%id)) write(*,'(1x,a,i0)') &
        '[WARN] DG ScaLAPACK PDSYEVX(Hstd vectors) workspace query failed: info=', ierr
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, y_div, ifail, iclustr, gap)
      return
    end if
    lwork = max(3, int(work_query(1)))
    liwork = max(1, iwork_query(1))
    allocate(work(lwork))
    allocate(iwork(liwork))

    y_div(:, :) = 0.0d0
    call PDSYEVX('V', 'V', 'L', n, h_div, 1, 1, desca, vl, vu, 1, n, abstol, &
                 m_found, nz_found, eval, orfac, y_div, 1, 1, descb, work, lwork, &
                 iwork, liwork, ifail, iclustr, gap, ierr)
    if (ierr /= 0 .or. m_found < nkeep .or. nz_found < nkeep) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,3(a,i0),2(a,1pe13.5))') '[WARN] DG ScaLAPACK PDSYEVX(Hstd vectors) failed/incomplete:', &
          ' info=', ierr, ' m=', m_found, ' nz=', nz_found, ' vl=', vl, ' vu=', vu
        if (ierr /= 0) then
          write(*,'(1x,a,1pe13.5)') '[WARN] DG ScaLAPACK PDSYEVX min_gap=', minval(gap)
          write(*,'(1x,a,4(i0,1x))') '[WARN] DG ScaLAPACK PDSYEVX ifail sample=', &
            ifail(1), ifail(min(2, n)), ifail(min(3, n)), ifail(min(4, n))
        end if
      end if
      call BLACS_GRIDEXIT(ictxt)
      deallocate(h_div, s_div, eval, work, iwork, y_div, ifail, iclustr, gap)
      return
    end if
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,2(a,1pe13.5),2(a,i0))') '[DG-DIST-EIG] ScaLAPACK partial Hstd diagonalization done', &
        ' eig_min=', eval(1), ' eig_keep=', eval(nkeep), ' m=', m_found, ' nz=', nz_found
      flush(6)
    end if

    ! Convert standard-problem eigenvectors y back to generalized vectors:
    ! S = L L^T and x = L^{-T} y.
    scalapack_alpha = 1.0d0
    deallocate(h_div)
    call PDTRSM('L', 'L', 'T', 'N', n, nkeep, scalapack_alpha, s_div, 1, 1, desca, &
                y_div, 1, 1, descb)
    deallocate(s_div)

    nlocal = size(dg_frag%coef, 1)
    allocate(coef_part(nlocal, nkeep))
    coef_part(:, :) = (0.0d0, 0.0d0)
    do i = 1, nlocal
      if (i > size(dg_frag%local_coef_global_ids, 1)) cycle
      idx = dg_frag%local_coef_global_ids(i, ispin)
      if (idx < 1 .or. idx > n) cycle
      call dg_scalapack_index_1d(idx, mb, nprow, proc_row, loc_row)
      if (proc_row /= myrow) cycle
      do j = 1, nkeep
        call dg_scalapack_index_1d(j, nb, npcol, proc_col, loc_col)
        if (proc_col == mycol .and. loc_row <= nloc_row .and. loc_col <= nloc_vec_col) then
          coef_part(i, j) = cmplx(y_div(loc_row, loc_col), 0.0d0, kind=8)
        end if
      end do
    end do
    deallocate(y_div)
    allocate(coef_sum(nlocal, nkeep))
    call comm_summation(coef_part, coef_sum, nlocal * nkeep, dg_frag%icomm)
    n_ortho = nkeep
    call measure_s_orthogonality_for_coef(coef_sum, ispin, n_ortho, s_diag_err, s_offdiag_max, s_frob_err)
    s_ortho_tol = 1.0d-7
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,3(a,1pe13.5))') '[DG-DIST-EIG-SORTHO] ncheck=', n_ortho, &
        ' diag_err=', s_diag_err, ' offdiag_max=', s_offdiag_max, ' frob_err=', s_frob_err
      flush(6)
    end if
    if (s_diag_err > s_ortho_tol .or. s_offdiag_max > s_ortho_tol) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,3(a,1pe13.5))') '[WARN] DG ScaLAPACK eigenvectors failed S-orthogonality check:', &
          ' diag_err=', s_diag_err, ' offdiag_max=', s_offdiag_max, ' tol=', s_ortho_tol
      end if
      call BLACS_GRIDEXIT(ictxt)
      deallocate(eval, work, iwork, ifail, iclustr, gap, coef_part, coef_sum)
      return
    end if
    dg_frag%coef(:, 1:nkeep, ispin) = coef_sum(:, 1:nkeep)
    if (allocated(dg_frag%coef_work)) dg_frag%coef_work(:, 1:nkeep, ispin) = dg_frag%coef(:, 1:nkeep, ispin)
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, 1:nkeep, ispin) = dg_frag%coef(:, 1:nkeep, ispin)
    dg_frag%esp(1:nkeep, ispin) = eval(1:nkeep)

    call BLACS_GRIDEXIT(ictxt)
    deallocate(eval, work, iwork, ifail, iclustr, gap, coef_part, coef_sum)
    did_solve = .true.
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,2(a,1pe13.5))') '[DG-DIST-EIG] done', ' eig_min=', dg_frag%esp(1, ispin), &
        ' eig_keep=', dg_frag%esp(nkeep, ispin)
      flush(6)
    end if
#endif
  contains
#ifdef USE_SCALAPACK
    subroutine assemble_scalapack_hs_from_blocks(hloc, sloc, ispin_in, include_s_in, do_log_in)
      implicit none
      real(8), intent(inout) :: hloc(:, :), sloc(:, :)
      integer, intent(in) :: ispin_in
      logical, intent(in), optional :: include_s_in, do_log_in

      integer :: nproc_pack, total_send, total_recv, ierr_pack
      integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
      integer, allocatable :: send_counts_meta(:), recv_counts_meta(:), send_displs_meta(:), recv_displs_meta(:)
      integer, allocatable :: send_meta(:), recv_meta(:), fill_counts(:)
      real(8), allocatable :: send_vals(:), recv_vals(:)
      integer :: p, k, row_l, col_l, kind_l
      real(8) :: h_abs_local, s_abs_local, h_abs_global, s_abs_global
      real(8) :: s_diag_min_local, s_diag_min_global, s_diag_max_local, s_diag_max_global
      real(8) :: s_diag_count_local, s_diag_count_global
      logical :: include_s, do_log

      include_s = .true.
      do_log = .true.
      if (present(include_s_in)) include_s = include_s_in
      if (present(do_log_in)) do_log = do_log_in
      nproc_pack = dg_frag%isize
      allocate(send_counts(nproc_pack), recv_counts(nproc_pack), send_displs(nproc_pack), recv_displs(nproc_pack))
      allocate(send_counts_meta(nproc_pack), recv_counts_meta(nproc_pack))
      allocate(send_displs_meta(nproc_pack), recv_displs_meta(nproc_pack), fill_counts(nproc_pack))
      send_counts(:) = 0

      if (dg_frag%is_frag_root) then
        call count_real_blocks(dg_frag%H_mat_blocks, 1, send_counts)
        if (allocated(dg_frag%H_nl_blocks)) call count_complex_blocks(dg_frag%H_nl_blocks, 1, send_counts)
        if (include_s) call count_real_blocks(dg_frag%S_mat_blocks, 2, send_counts)
      end if

      call MPI_Alltoall(send_counts, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, dg_frag%icomm, ierr_pack)
      send_displs(1) = 0
      recv_displs(1) = 0
      do p = 2, nproc_pack
        send_displs(p) = send_displs(p - 1) + send_counts(p - 1)
        recv_displs(p) = recv_displs(p - 1) + recv_counts(p - 1)
      end do
      total_send = send_displs(nproc_pack) + send_counts(nproc_pack)
      total_recv = recv_displs(nproc_pack) + recv_counts(nproc_pack)
      send_counts_meta(:) = 3 * send_counts(:)
      recv_counts_meta(:) = 3 * recv_counts(:)
      send_displs_meta(:) = 3 * send_displs(:)
      recv_displs_meta(:) = 3 * recv_displs(:)

      allocate(send_meta(max(1, 3 * total_send)), recv_meta(max(1, 3 * total_recv)))
      allocate(send_vals(max(1, total_send)), recv_vals(max(1, total_recv)))
      fill_counts(:) = 0
      if (dg_frag%is_frag_root) then
        call pack_real_blocks(dg_frag%H_mat_blocks, 1, send_displs, fill_counts, send_meta, send_vals)
        if (allocated(dg_frag%H_nl_blocks)) &
          call pack_complex_blocks(dg_frag%H_nl_blocks, 1, send_displs, fill_counts, send_meta, send_vals)
        if (include_s) call pack_real_blocks(dg_frag%S_mat_blocks, 2, send_displs, fill_counts, send_meta, send_vals)
      end if

      call MPI_Alltoallv(send_meta, send_counts_meta, send_displs_meta, MPI_INTEGER, &
                         recv_meta, recv_counts_meta, recv_displs_meta, MPI_INTEGER, dg_frag%icomm, ierr_pack)
      call MPI_Alltoallv(send_vals, send_counts, send_displs, MPI_DOUBLE_PRECISION, &
                         recv_vals, recv_counts, recv_displs, MPI_DOUBLE_PRECISION, dg_frag%icomm, ierr_pack)

      do k = 1, total_recv
        row_l = recv_meta(3 * k - 2)
        col_l = recv_meta(3 * k - 1)
        kind_l = recv_meta(3 * k)
        if (row_l < 1 .or. row_l > size(hloc, 1)) cycle
        if (col_l < 1 .or. col_l > size(hloc, 2)) cycle
        if (kind_l == 1) then
          hloc(row_l, col_l) = hloc(row_l, col_l) + recv_vals(k)
        else if (kind_l == 2 .and. include_s) then
          sloc(row_l, col_l) = sloc(row_l, col_l) + recv_vals(k)
        end if
      end do

      if (do_log) then
        h_abs_local = maxval(abs(hloc))
        s_abs_local = maxval(abs(sloc))
        h_abs_global = -h_abs_local
        call comm_get_min(h_abs_global, dg_frag%icomm)
        h_abs_global = -h_abs_global
        s_abs_global = -s_abs_local
        call comm_get_min(s_abs_global, dg_frag%icomm)
        s_abs_global = -s_abs_global
        call scalapack_s_diag_stats(sloc, s_diag_min_local, s_diag_max_local, s_diag_count_local)
        s_diag_min_global = s_diag_min_local
        call comm_get_min(s_diag_min_global, dg_frag%icomm)
        s_diag_max_global = -s_diag_max_local
        call comm_get_min(s_diag_max_global, dg_frag%icomm)
        s_diag_max_global = -s_diag_max_global
        call comm_summation(s_diag_count_local, s_diag_count_global, dg_frag%icomm)
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,5(a,1pe13.5))') '[DG-DIST-EIG-MAT] ScaLAPACK H/S assembled', &
            ' H_absmax=', h_abs_global, ' S_absmax=', s_abs_global, &
            ' S_diag_min=', s_diag_min_global, ' S_diag_max=', s_diag_max_global, &
            ' S_diag_count=', s_diag_count_global
          flush(6)
        end if
      end if

      deallocate(send_counts, recv_counts, send_displs, recv_displs)
      deallocate(send_counts_meta, recv_counts_meta, send_displs_meta, recv_displs_meta)
      deallocate(send_meta, recv_meta, send_vals, recv_vals, fill_counts)
    end subroutine assemble_scalapack_hs_from_blocks

    subroutine count_real_blocks(blocks, kind, counts)
      implicit none
      type(matrix_block_info), intent(in) :: blocks(:)
      integer, intent(in) :: kind
      integer, intent(inout) :: counts(:)
      integer :: iblk, ii, jj, ifr, jfr, nrow_b, ncol_b, ig, jg, target, lrow, lcol
      real(8) :: val

      do iblk = 1, size(blocks)
        if (.not. allocated(blocks(iblk)%val)) cycle
        ifr = blocks(iblk)%ifrag_row
        jfr = blocks(iblk)%ifrag_col
        if (ifr < 1 .or. ifr > dg_frag%n_frag) cycle
        if (jfr < 1 .or. jfr > dg_frag%n_frag) cycle
        nrow_b = min(dg_frag%n_basis(ifr, ispin), size(blocks(iblk)%val, 1))
        ncol_b = min(dg_frag%n_basis(jfr, ispin), size(blocks(iblk)%val, 2))
        do jj = 1, ncol_b
          jg = dg_frag%index_basis(jj, jfr, ispin)
          if (jg < 1 .or. jg > n) cycle
          do ii = 1, nrow_b
            ig = dg_frag%index_basis(ii, ifr, ispin)
            if (ig < jg .or. ig < 1 .or. ig > n) cycle
            val = blocks(iblk)%val(ii, jj, ispin)
            if (abs(val) <= 0.0d0) cycle
            call scalapack_owner_blockcyclic(ig, jg, target, lrow, lcol)
            counts(target + 1) = counts(target + 1) + 1
          end do
        end do
      end do
    end subroutine count_real_blocks

    subroutine count_complex_blocks(blocks, kind, counts)
      implicit none
      type(complex_matrix_block_info), intent(in) :: blocks(:)
      integer, intent(in) :: kind
      integer, intent(inout) :: counts(:)
      integer :: iblk, ii, jj, ifr, jfr, nrow_b, ncol_b, ig, jg, target, lrow, lcol
      real(8) :: val

      do iblk = 1, size(blocks)
        if (.not. allocated(blocks(iblk)%val)) cycle
        ifr = blocks(iblk)%ifrag_row
        jfr = blocks(iblk)%ifrag_col
        if (ifr < 1 .or. ifr > dg_frag%n_frag) cycle
        if (jfr < 1 .or. jfr > dg_frag%n_frag) cycle
        nrow_b = min(dg_frag%n_basis(ifr, ispin), size(blocks(iblk)%val, 1))
        ncol_b = min(dg_frag%n_basis(jfr, ispin), size(blocks(iblk)%val, 2))
        do jj = 1, ncol_b
          jg = dg_frag%index_basis(jj, jfr, ispin)
          if (jg < 1 .or. jg > n) cycle
          do ii = 1, nrow_b
            ig = dg_frag%index_basis(ii, ifr, ispin)
            if (ig < jg .or. ig < 1 .or. ig > n) cycle
            val = real(blocks(iblk)%val(ii, jj, ispin), kind=8)
            if (abs(val) <= 0.0d0) cycle
            call scalapack_owner_blockcyclic(ig, jg, target, lrow, lcol)
            counts(target + 1) = counts(target + 1) + 1
          end do
        end do
      end do
    end subroutine count_complex_blocks

    subroutine pack_real_blocks(blocks, kind, displs, filled, meta, vals)
      implicit none
      type(matrix_block_info), intent(in) :: blocks(:)
      integer, intent(in) :: kind, displs(:)
      integer, intent(inout) :: filled(:), meta(:)
      real(8), intent(inout) :: vals(:)
      integer :: iblk, ii, jj, ifr, jfr, nrow_b, ncol_b, ig, jg, target, lrow, lcol, pos
      real(8) :: val

      do iblk = 1, size(blocks)
        if (.not. allocated(blocks(iblk)%val)) cycle
        ifr = blocks(iblk)%ifrag_row
        jfr = blocks(iblk)%ifrag_col
        if (ifr < 1 .or. ifr > dg_frag%n_frag) cycle
        if (jfr < 1 .or. jfr > dg_frag%n_frag) cycle
        nrow_b = min(dg_frag%n_basis(ifr, ispin), size(blocks(iblk)%val, 1))
        ncol_b = min(dg_frag%n_basis(jfr, ispin), size(blocks(iblk)%val, 2))
        do jj = 1, ncol_b
          jg = dg_frag%index_basis(jj, jfr, ispin)
          if (jg < 1 .or. jg > n) cycle
          do ii = 1, nrow_b
            ig = dg_frag%index_basis(ii, ifr, ispin)
            if (ig < jg .or. ig < 1 .or. ig > n) cycle
            val = blocks(iblk)%val(ii, jj, ispin)
            if (abs(val) <= 0.0d0) cycle
            call scalapack_owner_blockcyclic(ig, jg, target, lrow, lcol)
            pos = displs(target + 1) + filled(target + 1) + 1
            meta(3 * pos - 2) = lrow
            meta(3 * pos - 1) = lcol
            meta(3 * pos) = kind
            vals(pos) = val
            filled(target + 1) = filled(target + 1) + 1
          end do
        end do
      end do
    end subroutine pack_real_blocks

    subroutine pack_complex_blocks(blocks, kind, displs, filled, meta, vals)
      implicit none
      type(complex_matrix_block_info), intent(in) :: blocks(:)
      integer, intent(in) :: kind, displs(:)
      integer, intent(inout) :: filled(:), meta(:)
      real(8), intent(inout) :: vals(:)
      integer :: iblk, ii, jj, ifr, jfr, nrow_b, ncol_b, ig, jg, target, lrow, lcol, pos
      real(8) :: val

      do iblk = 1, size(blocks)
        if (.not. allocated(blocks(iblk)%val)) cycle
        ifr = blocks(iblk)%ifrag_row
        jfr = blocks(iblk)%ifrag_col
        if (ifr < 1 .or. ifr > dg_frag%n_frag) cycle
        if (jfr < 1 .or. jfr > dg_frag%n_frag) cycle
        nrow_b = min(dg_frag%n_basis(ifr, ispin), size(blocks(iblk)%val, 1))
        ncol_b = min(dg_frag%n_basis(jfr, ispin), size(blocks(iblk)%val, 2))
        do jj = 1, ncol_b
          jg = dg_frag%index_basis(jj, jfr, ispin)
          if (jg < 1 .or. jg > n) cycle
          do ii = 1, nrow_b
            ig = dg_frag%index_basis(ii, ifr, ispin)
            if (ig < jg .or. ig < 1 .or. ig > n) cycle
            val = real(blocks(iblk)%val(ii, jj, ispin), kind=8)
            if (abs(val) <= 0.0d0) cycle
            call scalapack_owner_blockcyclic(ig, jg, target, lrow, lcol)
            pos = displs(target + 1) + filled(target + 1) + 1
            meta(3 * pos - 2) = lrow
            meta(3 * pos - 1) = lcol
            meta(3 * pos) = kind
            vals(pos) = val
            filled(target + 1) = filled(target + 1) + 1
          end do
        end do
      end do
    end subroutine pack_complex_blocks

    subroutine scalapack_owner_blockcyclic(ig, jg, target, lrow, lcol)
      implicit none
      integer, intent(in) :: ig, jg
      integer, intent(out) :: target, lrow, lcol
      integer :: prow, pcol

      call dg_scalapack_index_1d(ig, mb, nprow, prow, lrow)
      call dg_scalapack_index_1d(jg, nb, npcol, pcol, lcol)
      target = pcol * nprow + prow
    end subroutine scalapack_owner_blockcyclic

    subroutine scalapack_s_diag_stats(sloc, dmin, dmax, dcount)
      implicit none
      real(8), intent(in) :: sloc(:, :)
      real(8), intent(out) :: dmin, dmax, dcount
      integer :: iloc, jloc, ig, jg

      dmin = huge(1.0d0)
      dmax = -huge(1.0d0)
      dcount = 0.0d0
      do jloc = 1, size(sloc, 2)
        call dg_scalapack_global_index_1d(jloc, nb, npcol, mycol, jg)
        if (jg > n) cycle
        do iloc = 1, size(sloc, 1)
          call dg_scalapack_global_index_1d(iloc, mb, nprow, myrow, ig)
          if (ig /= jg .or. ig > n) cycle
          dmin = min(dmin, sloc(iloc, jloc))
          dmax = max(dmax, sloc(iloc, jloc))
          dcount = dcount + 1.0d0
        end do
      end do
    end subroutine scalapack_s_diag_stats

    subroutine measure_s_orthogonality_for_coef(coef_candidate, ispin_in, ncheck_in, diag_err, offdiag_max, frob_err)
      implicit none
      complex(8), intent(in) :: coef_candidate(:, :)
      integer, intent(in) :: ispin_in, ncheck_in
      real(8), intent(out) :: diag_err, offdiag_max, frob_err

      integer, parameter :: panel_ortho = 32
      integer :: ncheck, nloc, jb, je, nb_panel, ii, jj, jglob
      complex(8), allocatable :: cblk(:, :), sblk(:, :)
      complex(8), allocatable :: gram_local(:, :), gram_global(:, :)
      complex(8) :: target
      real(8) :: err

      ncheck = min(ncheck_in, size(coef_candidate, 2))
      nloc = size(coef_candidate, 1)
      diag_err = huge(1.0d0)
      offdiag_max = huge(1.0d0)
      frob_err = huge(1.0d0)
      if (ncheck <= 0 .or. nloc <= 0) return

      diag_err = 0.0d0
      offdiag_max = 0.0d0
      frob_err = 0.0d0
      allocate(cblk(nloc, panel_ortho), sblk(nloc, panel_ortho))
      allocate(gram_local(ncheck, panel_ortho), gram_global(ncheck, panel_ortho))

      do jb = 1, ncheck, panel_ortho
        je = min(ncheck, jb + panel_ortho - 1)
        nb_panel = je - jb + 1
        cblk(:, 1:nb_panel) = coef_candidate(:, jb:je)
        sblk(:, 1:nb_panel) = (0.0d0, 0.0d0)
        call apply_overlap_operator_batch_orbital_fragment_self(dg_frag, ispin_in, &
          cblk(:, 1:nb_panel), sblk(:, 1:nb_panel), .true.)

        gram_local(:, 1:nb_panel) = matmul(conjg(transpose(coef_candidate(:, 1:ncheck))), sblk(:, 1:nb_panel))
        call comm_summation(gram_local(:, 1:nb_panel), gram_global(:, 1:nb_panel), ncheck * nb_panel, dg_frag%icomm)

        do jj = 1, nb_panel
          jglob = jb + jj - 1
          do ii = 1, ncheck
            target = (0.0d0, 0.0d0)
            if (ii == jglob) target = (1.0d0, 0.0d0)
            err = abs(gram_global(ii, jj) - target)
            frob_err = frob_err + err * err
            if (ii == jglob) then
              diag_err = max(diag_err, err)
            else
              offdiag_max = max(offdiag_max, err)
            end if
          end do
        end do
      end do
      frob_err = sqrt(max(0.0d0, frob_err))

      deallocate(cblk, sblk, gram_local, gram_global)
    end subroutine measure_s_orthogonality_for_coef
#endif

    real(8) function dg_full_h_element(ig, jg, ispin_in, ifrag_of_gid, io_of_gid) result(val)
      integer, intent(in) :: ig, jg, ispin_in
      integer, intent(in) :: ifrag_of_gid(:), io_of_gid(:)
      integer :: fi, fj, ii, jj, iblk

      val = 0.0d0
      fi = ifrag_of_gid(ig); fj = ifrag_of_gid(jg)
      ii = io_of_gid(ig); jj = io_of_gid(jg)
      if (fi <= 0 .or. fj <= 0 .or. ii <= 0 .or. jj <= 0) return
      iblk = 0
      if (allocated(dg_frag%H_block_map)) iblk = dg_frag%H_block_map(fi, fj)
      if (iblk > 0 .and. iblk <= size(dg_frag%H_mat_blocks)) then
        val = val + dg_frag%H_mat_blocks(iblk)%val(ii, jj, ispin_in)
      else if (allocated(dg_frag%H_block_map)) then
        iblk = dg_frag%H_block_map(fj, fi)
        if (iblk > 0 .and. iblk <= size(dg_frag%H_mat_blocks)) &
          val = val + dg_frag%H_mat_blocks(iblk)%val(jj, ii, ispin_in)
      end if
      iblk = 0
      if (allocated(dg_frag%H_nl_block_map)) iblk = dg_frag%H_nl_block_map(fi, fj)
      if (iblk > 0 .and. allocated(dg_frag%H_nl_blocks) .and. iblk <= size(dg_frag%H_nl_blocks)) then
        val = val + real(dg_frag%H_nl_blocks(iblk)%val(ii, jj, ispin_in), kind=8)
      else if (allocated(dg_frag%H_nl_block_map)) then
        iblk = dg_frag%H_nl_block_map(fj, fi)
        if (iblk > 0 .and. allocated(dg_frag%H_nl_blocks) .and. iblk <= size(dg_frag%H_nl_blocks)) &
          val = val + real(conjg(dg_frag%H_nl_blocks(iblk)%val(jj, ii, ispin_in)), kind=8)
      end if
    end function dg_full_h_element

#ifdef USE_EIGENEXA
    subroutine gather_w_range(col_first, col_last, wcols)
      integer, intent(in) :: col_first, col_last
      real(8), intent(out) :: wcols(:, :)
      integer :: cols(panel_ex), ncols_local, cc

      ncols_local = col_last - col_first + 1
      do cc = 1, ncols_local
        cols(cc) = col_first + cc - 1
      end do
      call gather_w_columns(cols(1:ncols_local), ncols_local, wcols)
    end subroutine gather_w_range

    subroutine gather_w_columns(cols, ncols, wcols)
      integer, intent(in) :: ncols
      integer, intent(in) :: cols(ncols)
      real(8), intent(out) :: wcols(:, :)
      real(8), allocatable :: part(:, :)
      integer :: lx, ly, row_g, col_g, pos

      wcols(:, 1:ncols) = 0.0d0
      allocate(part(n_ex, ncols))
      part(:, :) = 0.0d0
      do ly = iy_s_ex, iy_e_ex
        col_g = eigen_translate_l2g(ly, y_nnod_ex, y_inod_ex)
        pos = column_position(col_g, cols, ncols)
        if (pos <= 0) cycle
        do lx = ix_s_ex, ix_e_ex
          row_g = eigen_translate_l2g(lx, x_nnod_ex, x_inod_ex)
          if (row_g < 1 .or. row_g > n_ex) cycle
          part(row_g, pos) = u_div_ex(lx, ly) / sqrt(lambda_ex(col_g))
        end do
      end do
      call comm_summation(part, wcols(:, 1:ncols), n_ex * ncols, dg_frag%icomm)
      deallocate(part)
    end subroutine gather_w_columns

    subroutine gather_y_range(col_first, col_last, ycols)
      integer, intent(in) :: col_first, col_last
      real(8), intent(out) :: ycols(:, :)
      real(8), allocatable :: part(:, :)
      integer :: ncols_local, lx, ly, row_g, col_g, pos

      ncols_local = col_last - col_first + 1
      ycols(:, 1:ncols_local) = 0.0d0
      allocate(part(n_ex, ncols_local))
      part(:, :) = 0.0d0
      do ly = iy_s_ex, iy_e_ex
        col_g = eigen_translate_l2g(ly, y_nnod_ex, y_inod_ex)
        if (col_g < col_first .or. col_g > col_last) cycle
        pos = col_g - col_first + 1
        do lx = ix_s_ex, ix_e_ex
          row_g = eigen_translate_l2g(lx, x_nnod_ex, x_inod_ex)
          if (row_g < 1 .or. row_g > n_ex) cycle
          part(row_g, pos) = y_div_ex(lx, ly)
        end do
      end do
      call comm_summation(part, ycols(:, 1:ncols_local), n_ex * ncols_local, dg_frag%icomm)
      deallocate(part)
    end subroutine gather_y_range

    integer function column_position(col, cols, ncols) result(pos)
      integer, intent(in) :: col, ncols
      integer, intent(in) :: cols(ncols)
      integer :: kk

      pos = 0
      do kk = 1, ncols
        if (cols(kk) == col) then
          pos = kk
          return
        end if
      end do
    end function column_position

#endif

    real(8) function dg_full_s_element(ig, jg, ispin_in, ifrag_of_gid, io_of_gid) result(val)
      integer, intent(in) :: ig, jg, ispin_in
      integer, intent(in) :: ifrag_of_gid(:), io_of_gid(:)
      integer :: fi, fj, ii, jj, iblk

      val = 0.0d0
      fi = ifrag_of_gid(ig); fj = ifrag_of_gid(jg)
      ii = io_of_gid(ig); jj = io_of_gid(jg)
      if (fi <= 0 .or. fj <= 0 .or. ii <= 0 .or. jj <= 0) return
      if (fi /= fj) return
#ifdef USE_EIGENEXA
      if (allocated(s_block_all_ex)) then
        if (ii <= size(s_block_all_ex, 1) .and. jj <= size(s_block_all_ex, 2) .and. &
            fi <= size(s_block_all_ex, 3)) then
          val = s_block_all_ex(ii, jj, fi)
          return
        end if
      end if
#endif
      iblk = 0
      if (allocated(dg_frag%S_block_map)) iblk = dg_frag%S_block_map(fi, fj)
      if (allocated(dg_frag%S_mat_blocks)) then
        if (iblk > 0 .and. iblk <= size(dg_frag%S_mat_blocks)) then
          val = dg_frag%S_mat_blocks(iblk)%val(ii, jj, ispin_in)
        else if (allocated(dg_frag%S_block_map)) then
          iblk = dg_frag%S_block_map(fj, fi)
          if (iblk > 0 .and. iblk <= size(dg_frag%S_mat_blocks)) &
            val = dg_frag%S_mat_blocks(iblk)%val(jj, ii, ispin_in)
        end if
      end if
    end function dg_full_s_element
  end subroutine diagonalize_initial_dg_full_distributed

  subroutine relax_initial_occupied_subspace_block_sparse(dg_frag, system)
    use structures
    use communication, only: comm_is_root, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system

    integer, parameter :: max_iter = 4
    character(len=16) :: env_relax
    integer :: iter, ispin, nocc, nrelax, io
    integer :: max_relax_occ, lobpcg_block, ib, ie, nb, env_len, env_status, env_read_status
    logical :: ortho_ok
    real(8) :: rel_res_global, hnorm_global, rnorm_global, best_rel
    real(8), allocatable :: eps(:), eps_best(:), num_local(:), den_local(:), num_global(:), den_global(:)
    real(8), allocatable :: occ_weight(:), sums_local(:), sums_global(:)
    complex(8), allocatable :: c(:,:), hc(:,:), sc(:,:), res(:,:), corr(:,:), pdir(:,:), coef_best(:,:)
    complex(8), allocatable :: x_old(:,:), x_new(:,:), p_new(:,:)

    if (.not. allocated(dg_frag%coef)) return
    if (.not. allocated(dg_frag%esp)) return
    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%S_mat_blocks)) return
    if (dg_frag%nspin /= 1) return

    ispin = 1
    nocc = min(dg_frag%nstate_tot, size(dg_frag%coef, 2), size(dg_frag%esp, 1))
    if (allocated(dg_frag%nocc_spin)) then
      if (ispin <= size(dg_frag%nocc_spin)) nocc = min(nocc, max(0, dg_frag%nocc_spin(ispin)))
    end if
    if (nocc <= 0) return
    allocate(occ_weight(nocc))
    call get_initial_state_spin_occ_info(dg_frag, system, ispin, occ_weight, nocc, nocc)
    if (nocc <= 0) then
      deallocate(occ_weight)
      return
    end if
    max_relax_occ = 4096
    env_relax = ''
    call get_environment_variable('SALMON_DG_BLOCK_SPARSE_RELAX_MAX_OCC', env_relax, &
                                  length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_relax(1:env_len), *, iostat=env_read_status) max_relax_occ
      if (env_read_status /= 0) max_relax_occ = 4096
    end if
    lobpcg_block = 64
    env_relax = ''
    call get_environment_variable('SALMON_DG_LOBPCG_BLOCK', env_relax, &
                                  length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_relax(1:env_len), *, iostat=env_read_status) lobpcg_block
      if (env_read_status /= 0) lobpcg_block = 64
    end if
    lobpcg_block = max(1, min(lobpcg_block, max_relax_occ))
    nrelax = min(nocc, max(1, max_relax_occ))
    if (nrelax <= 0) then
      deallocate(occ_weight)
      return
    end if

    allocate(c(size(dg_frag%coef, 1), nrelax), hc(size(dg_frag%coef, 1), nrelax), &
             sc(size(dg_frag%coef, 1), nrelax), res(size(dg_frag%coef, 1), nrelax), &
             corr(size(dg_frag%coef, 1), nrelax), pdir(size(dg_frag%coef, 1), nrelax), &
             coef_best(size(dg_frag%coef, 1), nrelax))
    allocate(eps(nrelax), eps_best(nrelax), num_local(nrelax), den_local(nrelax), &
             num_global(nrelax), den_global(nrelax), sums_local(2), sums_global(2))
    allocate(x_old(size(dg_frag%coef, 1), lobpcg_block), x_new(size(dg_frag%coef, 1), lobpcg_block), &
             p_new(size(dg_frag%coef, 1), lobpcg_block))

    call s_orthonormalize_coef_columns(dg_frag, ispin, nrelax, ortho_ok)
    if (.not. ortho_ok) then
      deallocate(c, hc, sc, res, corr, coef_best, eps, eps_best, num_local, den_local, num_global, den_global, &
                 occ_weight, sums_local, sums_global)
      return
    end if
    best_rel = huge(1.0d0)
    coef_best(:, :) = dg_frag%coef(:, 1:nrelax, ispin)
    eps_best(:) = dg_frag%esp(1:nrelax, ispin)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,1pe13.5)') '[DG-BLOCK-SPARSE-LOBPCG] active occupied block=', &
        nrelax, ' of nocc=', nocc, ' lobpcg_block=', lobpcg_block, ' occ=', occ_weight(1)
      flush(6)
    end if
    pdir(:, :) = (0.0d0, 0.0d0)

    do iter = 0, max_iter
      c(:, :) = dg_frag%coef(:, 1:nrelax, ispin)
      hc(:, :) = (0.0d0, 0.0d0)
      sc(:, :) = (0.0d0, 0.0d0)
      if (allocated(dg_frag%H_local_block_ids)) then
        call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, c, hc, dg_frag%H_local_block_ids)
      else
        call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, c, hc)
      end if
      if (allocated(dg_frag%H_nl_blocks)) then
        if (allocated(dg_frag%H_nl_local_block_ids)) then
          call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, c, hc, &
                                                 dg_frag%H_nl_local_block_ids)
        else
          call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, c, hc)
        end if
      end if
      call apply_overlap_operator_batch(dg_frag, ispin, c, sc, .true.)

      do io = 1, nrelax
        num_local(io) = real(sum(conjg(c(:, io)) * hc(:, io)), kind=8)
        den_local(io) = real(sum(conjg(c(:, io)) * sc(:, io)), kind=8)
      end do
      call comm_summation(num_local, num_global, nrelax, dg_frag%icomm)
      call comm_summation(den_local, den_global, nrelax, dg_frag%icomm)
      do io = 1, nrelax
        eps(io) = num_global(io) / max(den_global(io), 1.0d-300)
      end do

      res(:, :) = hc(:, :)
      do io = 1, nrelax
        res(:, io) = res(:, io) - eps(io) * sc(:, io)
      end do
      call remove_occupied_internal_residual(dg_frag, ispin, nrelax, res)
      sums_local(1) = sum(abs(res(:, :))**2)
      sums_local(2) = sum(abs(hc(:, :))**2)
      call comm_summation(sums_local, sums_global, 2, dg_frag%icomm)
      rnorm_global = sqrt(max(0.0d0, sums_global(1)))
      hnorm_global = sqrt(max(sums_global(2), 1.0d-300))
      rel_res_global = rnorm_global / hnorm_global
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,3(a,1pe13.5))') '[DG-BLOCK-SPARSE-LOBPCG] iter=', iter, &
          ' rel=', rel_res_global, ' rnorm=', rnorm_global, ' hnorm=', hnorm_global
        flush(6)
      end if

      if (rel_res_global <= best_rel) then
        best_rel = rel_res_global
        coef_best(:, :) = dg_frag%coef(:, 1:nrelax, ispin)
        eps_best(:) = eps(:)
      else
        dg_frag%coef(:, 1:nrelax, ispin) = coef_best(:, :)
        dg_frag%esp(1:nrelax, ispin) = eps_best(:)
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,2(a,1pe13.5))') '[WARN] DG block-sparse relaxation rejected: iter=', &
            iter, ' rel=', rel_res_global, ' best=', best_rel
          flush(6)
        end if
        exit
      end if
      if (iter == max_iter) exit

      corr(:, :) = (0.0d0, 0.0d0)
      call solve_overlap_operator_batch(dg_frag, ispin, res, corr, .true.)
      corr(:, :) = -corr(:, :)
      call remove_occupied_internal_residual(dg_frag, ispin, nrelax, corr)

      do ib = 1, nrelax, lobpcg_block
        ie = min(nrelax, ib + lobpcg_block - 1)
        nb = ie - ib + 1
        x_old(:, 1:nb) = dg_frag%coef(:, ib:ie, ispin)
        call lobpcg_ritz_update_block(dg_frag, ispin, ib, ie, iter, &
                                      corr(:, ib:ie), pdir(:, ib:ie), eps(ib:ie), &
                                      x_new(:, 1:nb), ortho_ok)
        if (.not. ortho_ok) exit
        p_new(:, 1:nb) = x_new(:, 1:nb) - x_old(:, 1:nb)
        dg_frag%coef(:, ib:ie, ispin) = x_new(:, 1:nb)
        pdir(:, ib:ie) = p_new(:, 1:nb)
      end do
      if (.not. ortho_ok) then
        dg_frag%coef(:, 1:nrelax, ispin) = coef_best(:, :)
        dg_frag%esp(1:nrelax, ispin) = eps_best(:)
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,1pe13.5)') '[WARN] DG block-sparse LOBPCG stopped after Ritz failure: iter=', &
            iter, ' best=', best_rel
          flush(6)
        end if
        exit
      end if
      call s_orthonormalize_coef_columns(dg_frag, ispin, nrelax, ortho_ok)
      if (.not. ortho_ok) then
        dg_frag%coef(:, 1:nrelax, ispin) = coef_best(:, :)
        dg_frag%esp(1:nrelax, ispin) = eps_best(:)
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,1pe13.5)') '[WARN] DG block-sparse relaxation stopped after S-orthonormalization failure: iter=', &
            iter, ' best=', best_rel
          flush(6)
        end if
        exit
      end if
      dg_frag%esp(1:nrelax, ispin) = eps(:)
    end do
    dg_frag%esp(1:nrelax, ispin) = eps_best(:)

    if (allocated(dg_frag%coef_work)) dg_frag%coef_work(:, 1:nrelax, ispin) = dg_frag%coef(:, 1:nrelax, ispin)
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, 1:nrelax, ispin) = dg_frag%coef(:, 1:nrelax, ispin)

    deallocate(c, hc, sc, res, corr, pdir, coef_best, x_old, x_new, p_new, eps, eps_best, &
               num_local, den_local, num_global, den_global, occ_weight, sums_local, sums_global)
  end subroutine relax_initial_occupied_subspace_block_sparse

  subroutine lobpcg_ritz_update_block(dg_frag, ispin, first_col, last_col, iter, w_in, p_in, eig_out, x_out, success)
    use communication, only: comm_summation
    use eigen_subdiag_sub, only: eigen_zheev
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin, first_col, last_col, iter
    complex(8), intent(in) :: w_in(:, :), p_in(:, :)
    real(8), intent(out) :: eig_out(:)
    complex(8), intent(out) :: x_out(:, :)
    logical, intent(out) :: success

    integer :: nrow, nb, ntrial, nvalid, i, j
    real(8), allocatable :: eval(:)
    complex(8), allocatable :: v(:,:), hv(:,:), hsmall(:,:), hsmall_global(:,:), evec(:,:)

    success = .false.
    nb = last_col - first_col + 1
    if (nb <= 0) return
    nrow = size(dg_frag%coef, 1)
    if (nrow <= 0) return
    ntrial = 2 * nb
    if (iter > 0) ntrial = 3 * nb

    allocate(v(nrow, ntrial))
    v(:, 1:nb) = dg_frag%coef(:, first_col:last_col, ispin)
    v(:, nb + 1:2 * nb) = w_in(:, 1:nb)
    if (iter > 0) v(:, 2 * nb + 1:3 * nb) = p_in(:, 1:nb)

    call s_orthonormalize_trial_matrix(dg_frag, ispin, v, ntrial, nvalid, success)
    if (.not. success .or. nvalid < nb) then
      deallocate(v)
      return
    end if

    allocate(hv(nrow, nvalid), hsmall(nvalid, nvalid), hsmall_global(nvalid, nvalid))
    hv(:, :) = (0.0d0, 0.0d0)
    if (allocated(dg_frag%H_local_block_ids)) then
      call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, v(:, 1:nvalid), hv, &
                                     dg_frag%H_local_block_ids)
    else
      call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, v(:, 1:nvalid), hv)
    end if
    if (allocated(dg_frag%H_nl_blocks)) then
      if (allocated(dg_frag%H_nl_local_block_ids)) then
        call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, v(:, 1:nvalid), hv, &
                                               dg_frag%H_nl_local_block_ids)
      else
        call apply_complex_matrix_blocks_batch(dg_frag, dg_frag%H_nl_blocks, ispin, v(:, 1:nvalid), hv)
      end if
    end if

    hsmall(:, :) = matmul(conjg(transpose(v(:, 1:nvalid))), hv(:, 1:nvalid))
    call comm_summation(hsmall, hsmall_global, nvalid * nvalid, dg_frag%icomm)
    do j = 1, nvalid
      do i = j + 1, nvalid
        hsmall_global(i, j) = 0.5d0 * (hsmall_global(i, j) + conjg(hsmall_global(j, i)))
        hsmall_global(j, i) = conjg(hsmall_global(i, j))
      end do
      hsmall_global(j, j) = cmplx(real(hsmall_global(j, j), kind=8), 0.0d0, kind=8)
    end do

    allocate(eval(nvalid), evec(nvalid, nvalid))
    hsmall(:, :) = hsmall_global(:, :)
    call eigen_zheev(hsmall, eval, evec)
    if (any(eval /= eval) .or. any(evec /= evec) .or. any(abs(evec) > huge(1.0d0))) then
      deallocate(v, hv, hsmall, hsmall_global, eval, evec)
      return
    end if
    x_out(:, 1:nb) = matmul(v(:, 1:nvalid), evec(1:nvalid, 1:nb))
    eig_out(1:nb) = eval(1:nb)
    success = .true.
    deallocate(v, hv, hsmall, hsmall_global, eval, evec)
  end subroutine lobpcg_ritz_update_block

  subroutine s_orthonormalize_trial_matrix(dg_frag, ispin, v, ntrial, nvalid, success)
    use communication, only: comm_summation
    use eigen_subdiag_sub, only: eigen_zheev
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin, ntrial
    complex(8), intent(inout) :: v(:, :)
    integer, intent(out) :: nvalid
    logical, intent(out) :: success

    integer :: nrow, i, j, k
    real(8), parameter :: norm_floor = 1.0d-24
    real(8), allocatable :: eval(:)
    complex(8), allocatable :: sv(:,:), gram(:,:), gram_global(:,:), evec(:,:), transform(:,:)

    success = .false.
    nvalid = 0
    if (ntrial <= 0) return
    nrow = size(v, 1)
    if (nrow <= 0) return

    allocate(sv(nrow, ntrial), gram(ntrial, ntrial), gram_global(ntrial, ntrial))
    sv(:, :) = (0.0d0, 0.0d0)
    call apply_overlap_operator_batch(dg_frag, ispin, v(:, 1:ntrial), sv, .true.)
    gram(:, :) = matmul(conjg(transpose(v(:, 1:ntrial))), sv(:, 1:ntrial))
    call comm_summation(gram, gram_global, ntrial * ntrial, dg_frag%icomm)
    do j = 1, ntrial
      do i = j + 1, ntrial
        gram_global(i, j) = 0.5d0 * (gram_global(i, j) + conjg(gram_global(j, i)))
        gram_global(j, i) = conjg(gram_global(i, j))
      end do
      gram_global(j, j) = cmplx(real(gram_global(j, j), kind=8), 0.0d0, kind=8)
    end do

    allocate(eval(ntrial), evec(ntrial, ntrial))
    gram(:, :) = gram_global(:, :)
    call eigen_zheev(gram, eval, evec)
    if (any(eval /= eval) .or. any(evec /= evec) .or. any(abs(evec) > huge(1.0d0))) then
      deallocate(sv, gram, gram_global, eval, evec)
      return
    end if
    nvalid = count(eval(1:ntrial) > norm_floor)
    if (nvalid <= 0) then
      deallocate(sv, gram, gram_global, eval, evec)
      return
    end if

    allocate(transform(ntrial, nvalid))
    k = 0
    do i = 1, ntrial
      if (eval(i) <= norm_floor) cycle
      k = k + 1
      transform(:, k) = evec(:, i) / sqrt(eval(i))
    end do
    v(:, 1:nvalid) = matmul(v(:, 1:ntrial), transform(:, 1:nvalid))
    success = .true.
    deallocate(sv, gram, gram_global, eval, evec, transform)
  end subroutine s_orthonormalize_trial_matrix

  subroutine remove_occupied_internal_residual(dg_frag, ispin, nocc, res)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin, nocc
    complex(8), intent(inout) :: res(:, :)

    integer, parameter :: panel_size = 32, state_chunk = 512
    integer :: ipanel, panel_end, bsz, istate, state_end, nts, nrow
    complex(8), allocatable :: q(:,:), sq(:,:), proj_local(:,:), proj_global(:,:)

    if (nocc <= 0) return
    if (.not. allocated(dg_frag%coef)) return
    nrow = size(dg_frag%coef, 1)
    if (nrow <= 0) return

    allocate(q(nrow, panel_size), sq(nrow, panel_size))
    allocate(proj_local(panel_size, state_chunk), proj_global(panel_size, state_chunk))

    do ipanel = 1, nocc, panel_size
      panel_end = min(nocc, ipanel + panel_size - 1)
      bsz = panel_end - ipanel + 1
      q(:, :) = (0.0d0, 0.0d0)
      sq(:, :) = (0.0d0, 0.0d0)
      q(:, 1:bsz) = dg_frag%coef(:, ipanel:panel_end, ispin)
      call apply_overlap_operator_batch(dg_frag, ispin, q(:, 1:bsz), sq(:, 1:bsz), .true.)

      do istate = 1, nocc, state_chunk
        state_end = min(nocc, istate + state_chunk - 1)
        nts = state_end - istate + 1
        proj_local(:, :) = (0.0d0, 0.0d0)
        proj_local(1:bsz, 1:nts) = matmul(conjg(transpose(q(:, 1:bsz))), res(:, istate:state_end))
        call comm_summation(proj_local(1:bsz, 1:nts), proj_global(1:bsz, 1:nts), &
                            bsz * nts, dg_frag%icomm)
        res(:, istate:state_end) = res(:, istate:state_end) - &
          matmul(sq(:, 1:bsz), proj_global(1:bsz, 1:nts))
      end do
    end do

    deallocate(q, sq, proj_local, proj_global)
  end subroutine remove_occupied_internal_residual

  subroutine s_orthonormalize_coef_columns(dg_frag, ispin, nocc, success)
    use communication, only: comm_summation
    use eigen_subdiag_sub, only: eigen_zheev
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin, nocc
    logical, intent(out) :: success

    integer, parameter :: panel_size = 32, trail_chunk = 512
    integer :: ipanel, panel_end, bsz, itrail, trail_end, nts, nrow, i, j
    real(8), parameter :: norm_floor = 1.0d-24
    real(8), allocatable :: eval(:)
    complex(8), allocatable :: q(:,:), sq(:,:), gram_local(:,:), gram_global(:,:)
    complex(8), allocatable :: gram_panel(:,:), evec(:,:), evec_scaled(:,:), inv_sqrt(:,:)
    complex(8), allocatable :: proj_local(:,:), proj_global(:,:)

    success = .false.
    if (nocc <= 0) return
    if (.not. allocated(dg_frag%coef)) return
    nrow = size(dg_frag%coef, 1)
    if (nrow <= 0) return

    allocate(q(nrow, panel_size), sq(nrow, panel_size))
    allocate(gram_local(panel_size, panel_size), gram_global(panel_size, panel_size))
    allocate(proj_local(panel_size, trail_chunk), proj_global(panel_size, trail_chunk))

    do ipanel = 1, nocc, panel_size
      panel_end = min(nocc, ipanel + panel_size - 1)
      bsz = panel_end - ipanel + 1

      q(:, :) = (0.0d0, 0.0d0)
      sq(:, :) = (0.0d0, 0.0d0)
      q(:, 1:bsz) = dg_frag%coef(:, ipanel:panel_end, ispin)
      call apply_overlap_operator_batch(dg_frag, ispin, q(:, 1:bsz), sq(:, 1:bsz), .true.)

      gram_local(:, :) = (0.0d0, 0.0d0)
      do j = 1, bsz
        do i = 1, bsz
          gram_local(i, j) = sum(conjg(q(:, i)) * sq(:, j))
        end do
      end do
      call comm_summation(gram_local(1:bsz, 1:bsz), gram_global(1:bsz, 1:bsz), bsz * bsz, dg_frag%icomm)
      do j = 1, bsz
        do i = j + 1, bsz
          gram_global(i, j) = 0.5d0 * (gram_global(i, j) + conjg(gram_global(j, i)))
          gram_global(j, i) = conjg(gram_global(i, j))
        end do
        gram_global(j, j) = cmplx(real(gram_global(j, j), kind=8), 0.0d0, kind=8)
      end do

      allocate(gram_panel(bsz, bsz), evec(bsz, bsz), evec_scaled(bsz, bsz), inv_sqrt(bsz, bsz), eval(bsz))
      gram_panel(:, :) = gram_global(1:bsz, 1:bsz)
      call eigen_zheev(gram_panel, eval, evec)
      if (minval(eval(1:bsz)) <= norm_floor) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,1pe13.5)') &
          '[WARN] S-panel orthonormalization found null occupied subspace: rank=', dg_frag%id, &
          ' ispin=', ispin, ' first_col=', ipanel, ' min_norm=', minval(eval(1:bsz))
        deallocate(gram_panel, evec, evec_scaled, inv_sqrt, eval)
        deallocate(q, sq, gram_local, gram_global, proj_local, proj_global)
        return
      end if
      evec_scaled(:, :) = evec(:, :)
      do i = 1, bsz
        evec_scaled(:, i) = evec_scaled(:, i) / sqrt(eval(i))
      end do
      inv_sqrt(:, :) = matmul(evec_scaled(:, :), conjg(transpose(evec(:, :))))
      q(:, 1:bsz) = matmul(q(:, 1:bsz), inv_sqrt(:, :))
      dg_frag%coef(:, ipanel:panel_end, ispin) = q(:, 1:bsz)
      deallocate(gram_panel, evec, evec_scaled, inv_sqrt, eval)

      sq(:, :) = (0.0d0, 0.0d0)
      call apply_overlap_operator_batch(dg_frag, ispin, q(:, 1:bsz), sq(:, 1:bsz), .true.)
      do itrail = panel_end + 1, nocc, trail_chunk
        trail_end = min(nocc, itrail + trail_chunk - 1)
        nts = trail_end - itrail + 1
        proj_local(:, :) = (0.0d0, 0.0d0)
        proj_local(1:bsz, 1:nts) = matmul(conjg(transpose(sq(:, 1:bsz))), &
                                           dg_frag%coef(:, itrail:trail_end, ispin))
        call comm_summation(proj_local(1:bsz, 1:nts), proj_global(1:bsz, 1:nts), &
                            bsz * nts, dg_frag%icomm)
        dg_frag%coef(:, itrail:trail_end, ispin) = dg_frag%coef(:, itrail:trail_end, ispin) - &
          matmul(q(:, 1:bsz), proj_global(1:bsz, 1:nts))
      end do
    end do

    deallocate(q, sq, gram_local, gram_global, proj_local, proj_global)
    success = .true.
  end subroutine s_orthonormalize_coef_columns

  subroutine build_face_cluster_fragment_list(dg_frag, ifrag, cfrag, ncf)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    integer, intent(out) :: cfrag(:)
    integer, intent(out) :: ncf

    integer :: axis, side, jfrag, i
    logical :: exists

    ncf = 0
    call append_cluster_fragment(ifrag)
    do axis = 1, 3
      if (dg_frag%num_fragment(axis) <= 1) cycle
      do side = -1, 1, 2
        jfrag = face_neighbor_fragment_rt(dg_frag, ifrag, axis, side)
        if (jfrag > 0) call append_cluster_fragment(jfrag)
      end do
    end do

  contains
    subroutine append_cluster_fragment(jfrag_in)
      integer, intent(in) :: jfrag_in

      if (jfrag_in < 1 .or. jfrag_in > dg_frag%n_frag) return
      exists = .false.
      do i = 1, ncf
        if (cfrag(i) == jfrag_in) exists = .true.
      end do
      if (exists) return
      if (ncf >= size(cfrag)) return
      ncf = ncf + 1
      cfrag(ncf) = jfrag_in
    end subroutine append_cluster_fragment
  end subroutine build_face_cluster_fragment_list

  integer function face_neighbor_fragment_rt(dg_frag, ifrag, axis, side) result(jfrag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, axis, side
    integer :: coords(3), rem

    jfrag = 0
    if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
    if (axis < 1 .or. axis > 3) return
    if (dg_frag%num_fragment(axis) <= 1) return
    coords(1) = (ifrag - 1) / (dg_frag%num_fragment(2) * dg_frag%num_fragment(3)) + 1
    rem = modulo(ifrag - 1, dg_frag%num_fragment(2) * dg_frag%num_fragment(3))
    coords(2) = rem / dg_frag%num_fragment(3) + 1
    coords(3) = modulo(rem, dg_frag%num_fragment(3)) + 1
    coords(axis) = modulo(coords(axis) - 1 + side + dg_frag%num_fragment(axis), &
                          dg_frag%num_fragment(axis)) + 1
    jfrag = ((coords(1) - 1) * dg_frag%num_fragment(2) + coords(2) - 1) * &
            dg_frag%num_fragment(3) + coords(3)
  end function face_neighbor_fragment_rt

  subroutine solve_fragment_generalized_eigen(H_in, S_in, n, nkeep, eig, C_keep, nkeep_valid)
    use eigen_subdiag_sub, only: eigen_dsyev
    implicit none
    integer, intent(in) :: n, nkeep
    real(8), intent(in) :: H_in(n, n), S_in(n, n)
    real(8), intent(out) :: eig(nkeep), C_keep(n, nkeep)
    integer, intent(out) :: nkeep_valid

    integer :: i, j, k, p, q, npos
    integer :: nskip_bad
    real(8) :: s_tol, smax
    real(8), parameter :: coeff_abs_limit = 1.0d8
    integer, allocatable :: pos(:)
    real(8), allocatable :: H_sym(:,:), S_sym(:,:), lambda(:), U(:,:), W(:,:), A(:,:), Z(:,:), eig_all(:)
    real(8), allocatable :: C_all(:,:)

    eig(:) = 0.0d0
    C_keep(:, :) = 0.0d0
    nkeep_valid = 0
    allocate(H_sym(n, n), S_sym(n, n), lambda(n), U(n, n))
    do j = 1, n
      do i = 1, n
        H_sym(i, j) = 0.5d0 * (H_in(i, j) + H_in(j, i))
        S_sym(i, j) = 0.5d0 * (S_in(i, j) + S_in(j, i))
      end do
    end do
    if (any(H_sym /= H_sym) .or. any(abs(H_sym) > huge(1.0d0)) .or. &
        any(S_sym /= S_sym) .or. any(abs(S_sym) > huge(1.0d0))) then
      write(*,'(1x,a)') "[FATAL] NaN/Inf in fragment-local generalized eigen input H/S"
      stop "DG-Fragment RT: invalid local eigen input"
    end if

    call eigen_dsyev(S_sym, lambda, U)
    if (any(lambda /= lambda) .or. any(abs(lambda) > huge(1.0d0)) .or. &
        any(U /= U) .or. any(abs(U) > huge(1.0d0))) then
      write(*,'(1x,a)') "[FATAL] NaN/Inf in fragment-local overlap diagonalization"
      stop "DG-Fragment RT: invalid local overlap eigensystem"
    end if
    smax = maxval(lambda(1:n))
    s_tol = max(1.0d-8, 1.0d-8 * max(1.0d0, smax))
    npos = count(lambda(1:n) > s_tol)
    if (npos <= 0) then
      deallocate(H_sym, S_sym, lambda, U)
      return
    end if

    allocate(pos(npos))
    p = 0
    do k = 1, n
      if (lambda(k) <= s_tol) cycle
      p = p + 1
      pos(p) = k
    end do

    allocate(W(n, npos), A(npos, npos), Z(npos, npos), eig_all(npos))
    do p = 1, npos
      W(1:n, p) = U(1:n, pos(p)) / sqrt(lambda(pos(p)))
    end do
    if (any(W /= W) .or. any(abs(W) > coeff_abs_limit)) then
      write(*,'(1x,a,1pe13.5,a,1pe13.5)') &
        "[FATAL] ill-conditioned fragment-local overlap transform: smax=", smax, " stol=", s_tol
      stop "DG-Fragment RT: ill-conditioned local overlap"
    end if
    A = matmul(transpose(W), matmul(H_sym, W))
    do q = 1, npos
      do p = q + 1, npos
        A(p, q) = 0.5d0 * (A(p, q) + A(q, p))
        A(q, p) = A(p, q)
      end do
    end do
    if (any(A /= A) .or. any(abs(A) > huge(1.0d0))) then
      write(*,'(1x,a)') "[FATAL] NaN/Inf in fragment-local orthonormal Hamiltonian"
      stop "DG-Fragment RT: invalid local orthonormal Hamiltonian"
    end if
    call eigen_dsyev(A, eig_all, Z)
    if (any(eig_all /= eig_all) .or. any(abs(eig_all) > huge(1.0d0)) .or. &
        any(Z /= Z) .or. any(abs(Z) > huge(1.0d0))) then
      write(*,'(1x,a)') "[FATAL] NaN/Inf in fragment-local Hamiltonian diagonalization"
      stop "DG-Fragment RT: invalid local Hamiltonian eigensystem"
    end if
    allocate(C_all(n, npos))
    C_all(1:n, 1:npos) = matmul(W, Z(1:npos, 1:npos))
    nskip_bad = 0
    do k = 1, npos
      if (any(C_all(1:n, k) /= C_all(1:n, k)) .or. any(abs(C_all(1:n, k)) > coeff_abs_limit)) then
        nskip_bad = nskip_bad + 1
        cycle
      end if
      nkeep_valid = nkeep_valid + 1
      C_keep(1:n, nkeep_valid) = C_all(1:n, k)
      eig(nkeep_valid) = eig_all(k)
      if (nkeep_valid >= nkeep) exit
    end do
    if (nskip_bad > 0) then
      write(*,'(1x,a,i0,a,i0,a,1pe13.5)') "[WARN] skipped ill-conditioned local eigenvectors: skipped=", &
        nskip_bad, " kept=", nkeep_valid, " coeff_limit=", coeff_abs_limit
    end if
    deallocate(H_sym, S_sym, lambda, U, W, pos, A, Z, eig_all, C_all)
  end subroutine solve_fragment_generalized_eigen

  subroutine build_fd_occupations(eig, nkeep, max_keep, nfrag, nspin, ne_target, temp_in, occ, mu, neout)
    implicit none
    integer, intent(in) :: max_keep, nfrag, nspin
    integer, intent(in) :: nkeep(nfrag, nspin)
    real(8), intent(in) :: eig(max_keep, nfrag, nspin), ne_target, temp_in
    real(8), intent(out) :: occ(max_keep, nfrag, nspin), mu, neout

    integer :: iter, ifrag, ispin, io
    real(8) :: emin, emax, mu1, mu2, mu3, ne3, temp_eff, wspin

    if (nspin == 1) then
      wspin = 2.0d0
    else
      wspin = 1.0d0
    end if
    temp_eff = max(temp_in, 1.0d-12)
    emin = huge(1.0d0)
    emax = -huge(1.0d0)
    do ispin = 1, nspin
      do ifrag = 1, nfrag
        do io = 1, nkeep(ifrag, ispin)
          emin = min(emin, eig(io, ifrag, ispin))
          emax = max(emax, eig(io, ifrag, ispin))
        end do
      end do
    end do
    if (emin > emax) then
      occ(:, :, :) = 0.0d0
      mu = 0.0d0
      neout = 0.0d0
      return
    end if

    mu1 = emin - max(1.0d0, 100.0d0 * temp_eff)
    mu2 = emax + max(1.0d0, 100.0d0 * temp_eff)
    if (fd_electron_count(eig, nkeep, max_keep, nfrag, nspin, mu2, temp_eff, wspin) < ne_target - 1.0d-8) then
      write(*,'(1x,a,1pe13.5,a,1pe13.5)') "[FATAL] fragment-local eigen cap cannot hold requested electrons: capacity=", &
        fd_electron_count(eig, nkeep, max_keep, nfrag, nspin, mu2, temp_eff, wspin), " target=", ne_target
      stop "DG-Fragment RT: local eigen cap electron capacity too small"
    end if
    do iter = 1, 240
      mu3 = 0.5d0 * (mu1 + mu2)
      ne3 = fd_electron_count(eig, nkeep, max_keep, nfrag, nspin, mu3, temp_eff, wspin)
      if (ne3 < ne_target) then
        mu1 = mu3
      else
        mu2 = mu3
      end if
      if (abs(ne3 - ne_target) < 1.0d-10) exit
    end do
    mu = 0.5d0 * (mu1 + mu2)
    occ(:, :, :) = 0.0d0
    neout = 0.0d0
    do ispin = 1, nspin
      do ifrag = 1, nfrag
        do io = 1, nkeep(ifrag, ispin)
          occ(io, ifrag, ispin) = fd_occ(eig(io, ifrag, ispin), mu, temp_eff, wspin)
          neout = neout + occ(io, ifrag, ispin)
        end do
      end do
    end do
  end subroutine build_fd_occupations

  real(8) function fd_electron_count(eig, nkeep, max_keep, nfrag, nspin, mu, temp_eff, wspin) result(neout)
    implicit none
    integer, intent(in) :: max_keep, nfrag, nspin
    integer, intent(in) :: nkeep(nfrag, nspin)
    real(8), intent(in) :: eig(max_keep, nfrag, nspin), mu, temp_eff, wspin
    integer :: ifrag, ispin, io

    neout = 0.0d0
    do ispin = 1, nspin
      do ifrag = 1, nfrag
        do io = 1, nkeep(ifrag, ispin)
          neout = neout + fd_occ(eig(io, ifrag, ispin), mu, temp_eff, wspin)
        end do
      end do
    end do
  end function fd_electron_count

  real(8) function fd_occ(eps, mu, temp_eff, wspin) result(occ)
    implicit none
    real(8), intent(in) :: eps, mu, temp_eff, wspin
    real(8) :: fact

    fact = (eps - mu) / temp_eff
    if (fact >= 40.0d0) then
      occ = 0.0d0
    else if (fact <= -40.0d0) then
      occ = wspin
    else
      occ = wspin / (1.0d0 + exp(fact))
    end if
  end function fd_occ

end module rt_dg_initial_state
