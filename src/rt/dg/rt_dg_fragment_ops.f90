module rt_dg_fragment_ops
  use communication, only: comm_bcast, comm_get_max, comm_is_root, comm_summation, COMM_GROUP_NULL
  use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info, complex_matrix_block_info, &
                                  invalidate_coef_exchange_cache
  implicit none

  private
  public :: ensure_nonlocal_pp_matrix_A
  public :: ensure_overlap_prop_available
  public :: calculate_microscopic_current_dg
  public :: calculate_macroscopic_current_dg
  public :: calculate_nonlocal_current_dg
  public :: calculate_local_wannier_polarization_dg
  public :: ensure_gradient_basis_cache
  public :: apply_gradient_to_basis
  public :: apply_momentum_blocks
  public :: ensure_nonlocal_projector_overlap_cache
  public :: apply_nonlocal_projector_overlap_batch
  public :: apply_matrix_blocks
  public :: apply_matrix_blocks_batch
  public :: apply_complex_matrix_blocks_batch
  public :: apply_mixed_hamiltonian
  public :: apply_mixed_hamiltonian_local_rows
  public :: rebuild_local_h_block_ids
  public :: copy_matrix_blocks_to_complex_dense
  public :: copy_momentum_blocks_to_complex_dense
  public :: symmetrize_real_matrix_blocks
  public :: mixed_fp_coupling_active
  public :: apply_overlap_operator
  public :: apply_overlap_operator_batch
  public :: solve_overlap_operator_batch
  public :: solve_overlap_operator_batch_local
  public :: copy_overlap_operator_to_dense
  public :: pack_owned_coef
  public :: fetch_remote_coef_rows
  public :: pack_owned_coef_pw
  public :: fetch_remote_coef_pw_rows
  public :: refresh_pw_coef_cache
  public :: gather_fragment_coef_view
  public :: gather_full_coef_view
  public :: zero_nonowned_coefficients
  public :: zero_nonlocal_h_matrix_blocks

contains

  integer function dg_state_owner_offset(nstate_tot, nworker, state_col) result(owner_offset)
    implicit none
    integer, intent(in) :: nstate_tot, nworker, state_col
    integer :: worker, base, extra, state_s, state_e

    owner_offset = 0
    if (nworker <= 1 .or. nstate_tot <= 0) return
    base = nstate_tot / nworker
    extra = mod(nstate_tot, nworker)
    do worker = 0, nworker - 1
      if (worker < extra) then
        state_s = worker * (base + 1) + 1
        state_e = state_s + base
      else
        state_s = extra * (base + 1) + (worker - extra) * base + 1
        state_e = state_s + base - 1
      end if
      if (state_col >= state_s .and. state_col <= state_e) then
        owner_offset = worker
        return
      end if
    end do
    owner_offset = max(0, min(nworker - 1, nworker - 1))
  end function dg_state_owner_offset

  subroutine ensure_gradient_basis_cache(dg_frag, mg, stencil)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil

    integer :: ifrag_count, i_local, ifrag, jo
    integer :: nx_max, ny_max, nz_max
    integer :: nxyz(3)

    if (.not. dg_frag%has_real_space_basis) return
    if (dg_frag%gradient_basis_cache_valid) return

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) return

    nx_max = 0
    ny_max = 0
    nz_max = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      nx_max = max(nx_max, dg_frag%nxyz_domain(1, ifrag))
      ny_max = max(ny_max, dg_frag%nxyz_domain(2, ifrag))
      nz_max = max(nz_max, dg_frag%nxyz_domain(3, ifrag))
    end do

    if (allocated(dg_frag%gradient_basis_cache)) then
      if (size(dg_frag%gradient_basis_cache, 1) /= nx_max .or. &
          size(dg_frag%gradient_basis_cache, 2) /= ny_max .or. &
          size(dg_frag%gradient_basis_cache, 3) /= nz_max .or. &
          size(dg_frag%gradient_basis_cache, 5) /= dg_frag%nstate_frag .or. &
          size(dg_frag%gradient_basis_cache, 6) /= ifrag_count) then
        deallocate(dg_frag%gradient_basis_cache)
      end if
    end if
    if (.not. allocated(dg_frag%gradient_basis_cache)) then
      allocate(dg_frag%gradient_basis_cache(nx_max, ny_max, nz_max, 3, dg_frag%nstate_frag, ifrag_count))
    end if
    dg_frag%gradient_basis_cache(:, :, :, :, :, :) = 0.0d0

    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      nxyz(:) = dg_frag%nxyz_domain(:, ifrag)
      do jo = 1, dg_frag%n_basis(ifrag, 1)
        ! spin-independent real-space basis; use max active basis count across spins below
      end do
      do jo = 1, maxval(dg_frag%n_basis(ifrag, 1:dg_frag%nspin))
        call apply_gradient_to_basis(dg_frag, i_local, jo, mg, stencil, &
          dg_frag%gradient_basis_cache(1:nxyz(1), 1:nxyz(2), 1:nxyz(3), 1:3, jo, i_local))
      end do
    end do

    dg_frag%gradient_basis_cache_valid = .true.
  end subroutine ensure_gradient_basis_cache

  real(8) function mixed_fp_maxabs(dg_frag, ispin) result(maxabs_fp)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin

    maxabs_fp = 0.0d0
    if (.not. dg_frag%use_plane_wave_basis) return
    if (.not. allocated(dg_frag%S_mat_frag_pw)) return
    if (ispin < 1 .or. ispin > size(dg_frag%S_mat_frag_pw, 3)) return
    if (size(dg_frag%S_mat_frag_pw, 1) <= 0 .or. size(dg_frag%S_mat_frag_pw, 2) <= 0) return
    maxabs_fp = maxval(abs(dg_frag%S_mat_frag_pw(:, :, ispin)))
  end function mixed_fp_maxabs

  logical function mixed_fp_coupling_active(dg_frag, ispin) result(active)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    real(8), parameter :: fp_tol = 1.0d-12

    active = (mixed_fp_maxabs(dg_frag, ispin) > fp_tol)
  end function mixed_fp_coupling_active

  subroutine symmetrize_real_matrix_blocks(dg_frag, blocks)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(inout) :: blocks(:)
    integer :: ispin, iblk, nbf, io, jo

    !$omp parallel do collapse(2) private(ispin,iblk,nbf,jo,io) schedule(static)
    do ispin = 1, dg_frag%nspin
      do iblk = 1, size(blocks)
        if (blocks(iblk)%ifrag_row /= blocks(iblk)%ifrag_col) cycle
        nbf = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        if (nbf <= 0) cycle
        do jo = 1, nbf
          do io = jo + 1, nbf
            blocks(iblk)%val(io, jo, ispin) = 0.5d0 * (blocks(iblk)%val(io, jo, ispin) + blocks(iblk)%val(jo, io, ispin))
            blocks(iblk)%val(jo, io, ispin) = blocks(iblk)%val(io, jo, ispin)
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine symmetrize_real_matrix_blocks

  logical function is_runtime_neighbor_axis(lg, s1, n1, s2, n2) result(ok)
    implicit none
    integer, intent(in) :: lg, s1, n1, s2, n2
    integer :: e1, e2, s1_next, s2_next

    e1 = s1 + n1 - 1
    e2 = s2 + n2 - 1
    s1_next = modulo(e1, lg) + 1
    s2_next = modulo(e2, lg) + 1
    ok = ((s1 == s2) .and. (n1 == n2)) .or. (s1 == s2_next) .or. (s2 == s1_next)
  end function is_runtime_neighbor_axis

  logical function is_runtime_neighbor_pair(dg_frag, ifrag_row, ifrag_col) result(is_pair)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row, ifrag_col
    integer :: axis
    logical :: axis_ok(3)

    is_pair = .false.
    if (ifrag_row == ifrag_col) then
      is_pair = .true.
      return
    end if
    if (allocated(dg_frag%runtime_neighbor_pair_cache)) then
      if (ifrag_row >= 1 .and. ifrag_row <= size(dg_frag%runtime_neighbor_pair_cache, 1) .and. &
          ifrag_col >= 1 .and. ifrag_col <= size(dg_frag%runtime_neighbor_pair_cache, 2)) then
        is_pair = dg_frag%runtime_neighbor_pair_cache(ifrag_row, ifrag_col)
        return
      end if
    end if

    do axis = 1, 3
      axis_ok(axis) = is_runtime_neighbor_axis(dg_frag%lgnum_total(axis), &
        dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
        dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col))
    end do

    is_pair = all(axis_ok)
  end function is_runtime_neighbor_pair

  subroutine ensure_runtime_neighbor_pair_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag_row, ifrag_col, axis
    logical :: axis_ok(3)
    logical, save :: option_initialized = .false.
    logical, save :: enable_pair_cache = .false.
    character(16) :: env_value
    integer :: env_status

    if (allocated(dg_frag%runtime_neighbor_pair_cache)) return
    if (.not. option_initialized) then
      env_value = ''
      call get_environment_variable('SALMON_DG_CACHE_RUNTIME_NEIGHBOR_PAIRS', env_value, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_value)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_pair_cache = .true.
        end select
      end if
      option_initialized = .true.
    end if
    if (.not. enable_pair_cache) return

    allocate(dg_frag%runtime_neighbor_pair_cache(dg_frag%n_frag, dg_frag%n_frag))
    dg_frag%runtime_neighbor_pair_cache(:, :) = .false.
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (ifrag_row == ifrag_col) then
          dg_frag%runtime_neighbor_pair_cache(ifrag_row, ifrag_col) = .true.
        else
          do axis = 1, 3
            axis_ok(axis) = is_runtime_neighbor_axis(dg_frag%lgnum_total(axis), &
              dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
              dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col))
          end do
          dg_frag%runtime_neighbor_pair_cache(ifrag_row, ifrag_col) = all(axis_ok)
        end if
      end do
    end do
  end subroutine ensure_runtime_neighbor_pair_cache

  subroutine ensure_runtime_neighbor_list_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag, jfrag, nnei, max_nnei

    if (allocated(dg_frag%runtime_neighbor_frag_count) .and. &
        allocated(dg_frag%runtime_neighbor_frag_ids)) then
      if (size(dg_frag%runtime_neighbor_frag_count) == dg_frag%n_frag .and. &
          size(dg_frag%runtime_neighbor_frag_ids, 2) == dg_frag%n_frag) return
      deallocate(dg_frag%runtime_neighbor_frag_count, dg_frag%runtime_neighbor_frag_ids)
    end if

    allocate(dg_frag%runtime_neighbor_frag_count(max(1, dg_frag%n_frag)))
    dg_frag%runtime_neighbor_frag_count(:) = 0
    max_nnei = 0
    do ifrag = 1, dg_frag%n_frag
      nnei = 0
      do jfrag = 1, dg_frag%n_frag
        if (is_runtime_neighbor_pair(dg_frag, ifrag, jfrag)) nnei = nnei + 1
      end do
      dg_frag%runtime_neighbor_frag_count(ifrag) = nnei
      max_nnei = max(max_nnei, nnei)
    end do

    allocate(dg_frag%runtime_neighbor_frag_ids(max(1, max_nnei), max(1, dg_frag%n_frag)))
    dg_frag%runtime_neighbor_frag_ids(:, :) = 0
    do ifrag = 1, dg_frag%n_frag
      nnei = 0
      do jfrag = 1, dg_frag%n_frag
        if (.not. is_runtime_neighbor_pair(dg_frag, ifrag, jfrag)) cycle
        nnei = nnei + 1
        dg_frag%runtime_neighbor_frag_ids(nnei, ifrag) = jfrag
      end do
    end do
  end subroutine ensure_runtime_neighbor_list_cache

  subroutine ensure_coef_row_fragment_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ispin, ifrag, ib, gid

    if (allocated(dg_frag%coef_row_fragment)) then
      if (size(dg_frag%coef_row_fragment, 1) == dg_frag%n_mat_max .and. &
          size(dg_frag%coef_row_fragment, 2) == dg_frag%nspin) return
      deallocate(dg_frag%coef_row_fragment)
    end if
    allocate(dg_frag%coef_row_fragment(max(1, dg_frag%n_mat_max), max(1, dg_frag%nspin)))
    dg_frag%coef_row_fragment(:, :) = 0
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. allocated(dg_frag%n_basis)) return

    do ispin = 1, dg_frag%nspin
      do ifrag = 1, dg_frag%n_frag
        do ib = 1, min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
          gid = dg_frag%index_basis(ib, ifrag, ispin)
          if (gid < 1 .or. gid > dg_frag%n_mat_max) cycle
          dg_frag%coef_row_fragment(gid, ispin) = ifrag
        end do
      end do
    end do
  end subroutine ensure_coef_row_fragment_cache

  subroutine ensure_coef_exchange_peer_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    logical, allocatable :: peer_needed(:), owned_frag(:)
    integer :: ispin, ifrag_owned, ifrag_req, ib, gid, requester, npeer, peer, inei

    if (dg_frag%coef_exchange_peer_cache_valid .and. allocated(dg_frag%coef_exchange_peer_ranks) .and. &
        allocated(dg_frag%coef_exchange_peer_count)) return
    if (allocated(dg_frag%coef_exchange_peer_ranks)) deallocate(dg_frag%coef_exchange_peer_ranks)
    if (allocated(dg_frag%coef_exchange_peer_count)) deallocate(dg_frag%coef_exchange_peer_count)
    allocate(dg_frag%coef_exchange_peer_ranks(max(1, dg_frag%isize), max(1, dg_frag%nspin)))
    allocate(dg_frag%coef_exchange_peer_count(max(1, dg_frag%nspin)))
    dg_frag%coef_exchange_peer_ranks(:, :) = -1
    dg_frag%coef_exchange_peer_count(:) = 0

    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%index_basis)) return
    call ensure_runtime_neighbor_pair_cache(dg_frag)
    call ensure_runtime_neighbor_list_cache(dg_frag)

    allocate(peer_needed(0:max(0, dg_frag%isize-1)))
    allocate(owned_frag(max(1, dg_frag%n_frag)))
    do ispin = 1, dg_frag%nspin
      peer_needed(:) = .false.
      owned_frag(:) = .false.

      do ifrag_owned = 1, dg_frag%n_frag
        do ib = 1, min(dg_frag%n_basis(ifrag_owned, ispin), size(dg_frag%index_basis, 1))
          gid = dg_frag%index_basis(ib, ifrag_owned, ispin)
          if (gid < 1 .or. gid > size(dg_frag%coef_owner, 1)) cycle
          if (dg_frag%coef_owner(gid, ispin) == dg_frag%id) then
            owned_frag(ifrag_owned) = .true.
            exit
          end if
        end do
      end do

      do ifrag_owned = 1, dg_frag%n_frag
        if (.not. owned_frag(ifrag_owned)) cycle
        do inei = 1, dg_frag%runtime_neighbor_frag_count(ifrag_owned)
          ifrag_req = dg_frag%runtime_neighbor_frag_ids(inei, ifrag_owned)
          if (ifrag_req < 1 .or. ifrag_req > dg_frag%n_frag) cycle
          do ib = 1, min(dg_frag%n_basis(ifrag_req, ispin), size(dg_frag%index_basis, 1))
            gid = dg_frag%index_basis(ib, ifrag_req, ispin)
            if (gid < 1 .or. gid > size(dg_frag%coef_owner, 1)) cycle
            requester = dg_frag%coef_owner(gid, ispin)
            if (requester < 0 .or. requester >= dg_frag%isize) cycle
            if (requester == dg_frag%id) cycle
            peer_needed(requester) = .true.
          end do
        end do
      end do

      npeer = 0
      do peer = 0, dg_frag%isize - 1
        if (.not. peer_needed(peer)) cycle
        npeer = npeer + 1
        dg_frag%coef_exchange_peer_ranks(npeer, ispin) = peer
      end do
      dg_frag%coef_exchange_peer_count(ispin) = npeer
    end do
    deallocate(peer_needed, owned_frag)
    dg_frag%coef_exchange_peer_cache_valid = .true.
  end subroutine ensure_coef_exchange_peer_cache

  subroutine ensure_coef_allowed_request_frag_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    logical, allocatable :: local_frag(:)
    integer :: i, ifrag_local, ifrag_req, ispin, ib, gid, inei

    if (allocated(dg_frag%coef_allowed_request_frag)) then
      if (size(dg_frag%coef_allowed_request_frag) == dg_frag%n_frag) return
      deallocate(dg_frag%coef_allowed_request_frag)
    end if
    allocate(dg_frag%coef_allowed_request_frag(max(1, dg_frag%n_frag)))
    dg_frag%coef_allowed_request_frag(:) = .false.
    call ensure_runtime_neighbor_pair_cache(dg_frag)
    call ensure_runtime_neighbor_list_cache(dg_frag)

    allocate(local_frag(max(1, dg_frag%n_frag)))
    local_frag(:) = .false.
    if (allocated(dg_frag%H_local_rows)) then
      do i = 1, size(dg_frag%H_local_rows)
        ifrag_local = dg_frag%H_local_rows(i)
        if (ifrag_local >= 1 .and. ifrag_local <= dg_frag%n_frag) local_frag(ifrag_local) = .true.
      end do
    else if (allocated(dg_frag%coef_owner) .and. allocated(dg_frag%index_basis)) then
      do ispin = 1, dg_frag%nspin
        do ifrag_local = 1, dg_frag%n_frag
          do ib = 1, min(dg_frag%n_basis(ifrag_local, ispin), size(dg_frag%index_basis, 1))
            gid = dg_frag%index_basis(ib, ifrag_local, ispin)
            if (gid < 1 .or. gid > size(dg_frag%coef_owner, 1)) cycle
            if (dg_frag%coef_owner(gid, ispin) == dg_frag%id) then
              local_frag(ifrag_local) = .true.
              exit
            end if
          end do
        end do
      end do
    end if
    if (.not. any(local_frag)) then
      do ifrag_local = max(1, dg_frag%ifrag_start), min(dg_frag%n_frag, dg_frag%ifrag_end)
        local_frag(ifrag_local) = .true.
      end do
    end if

    do ifrag_local = 1, dg_frag%n_frag
      if (.not. local_frag(ifrag_local)) cycle
      do inei = 1, dg_frag%runtime_neighbor_frag_count(ifrag_local)
        ifrag_req = dg_frag%runtime_neighbor_frag_ids(inei, ifrag_local)
        if (ifrag_req >= 1 .and. ifrag_req <= dg_frag%n_frag) dg_frag%coef_allowed_request_frag(ifrag_req) = .true.
      end do
    end do
    deallocate(local_frag)
  end subroutine ensure_coef_allowed_request_frag_cache

  subroutine ensure_nonlocal_projector_overlap_cache(dg_frag, mg, ppg, nspin, hvol, Ac_tot)
    use structures
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    integer,                intent(in)    :: nspin
    real(8),                intent(in)    :: hvol
    real(8),                intent(in)    :: Ac_tot(3)

    integer :: ifrag_count, i_local, ifrag, ispin, io, ilma, ia, j, i_halo
    integer :: ix, iy, iz, lx, ly, lz, nbf_row
    integer :: is(3), ie(3)
    real(8) :: x, y, z, phase, delta_A
    complex(8) :: basis_val, overlap_i, overlap_r_i(3)
    complex(8), allocatable :: proj_local(:,:), proj_sum(:,:), proj_r_local(:,:,:), proj_r_sum(:,:,:)

    if (ppg%Nlma <= 0) return
    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) return

    delta_A = maxval(abs(Ac_tot - dg_frag%Ac_nl_projector_cache))
    if (dg_frag%has_nl_projector_cache .and. allocated(dg_frag%nl_projector_overlap)) then
      if (size(dg_frag%nl_projector_overlap, 1) == dg_frag%nstate_frag .and. &
          size(dg_frag%nl_projector_overlap, 2) == ppg%Nlma .and. &
          size(dg_frag%nl_projector_overlap, 3) == nspin .and. &
          size(dg_frag%nl_projector_overlap, 4) == ifrag_count .and. &
          allocated(dg_frag%nl_projector_overlap_halo) .and. &
          allocated(dg_frag%nl_projector_r_overlap) .and. &
          allocated(dg_frag%nl_projector_r_overlap_halo) .and. &
          size(dg_frag%nl_projector_overlap_halo, 1) == dg_frag%nstate_frag .and. &
          size(dg_frag%nl_projector_overlap_halo, 2) == ppg%Nlma .and. &
          size(dg_frag%nl_projector_overlap_halo, 3) == nspin .and. &
          size(dg_frag%nl_projector_overlap_halo, 4) == max(1, dg_frag%n_halo) .and. &
          dg_frag%nl_projector_cache_nlma == ppg%Nlma .and. &
          delta_A <= dg_frag%Ac_nl_cache_tol) return
    end if

    if (allocated(dg_frag%nl_projector_overlap)) deallocate(dg_frag%nl_projector_overlap)
    if (allocated(dg_frag%nl_projector_overlap_halo)) deallocate(dg_frag%nl_projector_overlap_halo)
    if (allocated(dg_frag%nl_projector_r_overlap)) deallocate(dg_frag%nl_projector_r_overlap)
    if (allocated(dg_frag%nl_projector_r_overlap_halo)) deallocate(dg_frag%nl_projector_r_overlap_halo)
    allocate(dg_frag%nl_projector_overlap(dg_frag%nstate_frag, ppg%Nlma, nspin, ifrag_count))
    allocate(dg_frag%nl_projector_overlap_halo(dg_frag%nstate_frag, ppg%Nlma, nspin, max(1, dg_frag%n_halo)))
    allocate(dg_frag%nl_projector_r_overlap(3, dg_frag%nstate_frag, ppg%Nlma, nspin, ifrag_count))
    allocate(dg_frag%nl_projector_r_overlap_halo(3, dg_frag%nstate_frag, ppg%Nlma, nspin, max(1, dg_frag%n_halo)))
    dg_frag%nl_projector_overlap(:, :, :, :) = (0.0d0, 0.0d0)
    dg_frag%nl_projector_overlap_halo(:, :, :, :) = (0.0d0, 0.0d0)
    dg_frag%nl_projector_r_overlap(:, :, :, :, :) = (0.0d0, 0.0d0)
    dg_frag%nl_projector_r_overlap_halo(:, :, :, :, :) = (0.0d0, 0.0d0)
    allocate(proj_local(dg_frag%nstate_frag, ppg%Nlma))
    allocate(proj_sum(dg_frag%nstate_frag, ppg%Nlma))
    allocate(proj_r_local(3, dg_frag%nstate_frag, ppg%Nlma))
    allocate(proj_r_sum(3, dg_frag%nstate_frag, ppg%Nlma))

    is = mg%is
    ie = mg%ie
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      do ispin = 1, nspin
        nbf_row = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
        proj_local(:, :) = (0.0d0, 0.0d0)
        proj_sum(:, :) = (0.0d0, 0.0d0)
        proj_r_local(:, :, :) = (0.0d0, 0.0d0)
        proj_r_sum(:, :, :) = (0.0d0, 0.0d0)
        if (nbf_row > 0) then
          !$omp parallel do collapse(2) &
          !$omp& private(ilma, io, ia, j, ix, iy, iz, lx, ly, lz, x, y, z, phase, basis_val, overlap_i, overlap_r_i)
          do ilma = 1, ppg%Nlma
            do io = 1, nbf_row
              ia = ppg%ia_tbl(ilma)
              overlap_i = (0.0d0, 0.0d0)
              overlap_r_i(:) = (0.0d0, 0.0d0)
              do j = 1, ppg%mps(ia)
                ix = ppg%jxyz(1, j, ia)
                iy = ppg%jxyz(2, j, ia)
                iz = ppg%jxyz(3, j, ia)
                if (ix < is(1) .or. ix > ie(1) .or. &
                    iy < is(2) .or. iy > ie(2) .or. &
                    iz < is(3) .or. iz > ie(3)) cycle
                lx = map_global_to_phi_box_coord(ix, lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1), &
                                                 dg_frag%lgnum_total(1))
                ly = map_global_to_phi_box_coord(iy, lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2), &
                                                 dg_frag%lgnum_total(2))
                lz = map_global_to_phi_box_coord(iz, lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3), &
                                                 dg_frag%lgnum_total(3))
                if (lx == 0 .or. ly == 0 .or. lz == 0) cycle
                x = ppg%rxyz(1, j, ia)
                y = ppg%rxyz(2, j, ia)
                z = ppg%rxyz(3, j, ia)
                phase = Ac_tot(1) * x + Ac_tot(2) * y + Ac_tot(3) * z
                if (allocated(dg_frag%phi_frag_c)) then
                  basis_val = dg_frag%phi_frag_c(lx, ly, lz, io, i_local)
                else
                  basis_val = cmplx(dg_frag%phi_frag(lx, ly, lz, io, i_local), 0.0d0, kind=8)
                end if
                overlap_i = overlap_i + basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
                overlap_r_i(1) = overlap_r_i(1) + x * basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
                overlap_r_i(2) = overlap_r_i(2) + y * basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
                overlap_r_i(3) = overlap_r_i(3) + z * basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
              end do
              proj_local(io, ilma) = overlap_i
              proj_r_local(1:3, io, ilma) = overlap_r_i(1:3)
            end do
          end do
          !$omp end parallel do
        end if
        call comm_summation(proj_local, proj_sum, dg_frag%nstate_frag * ppg%Nlma, dg_frag%icomm_frag)
        call comm_summation(proj_r_local, proj_r_sum, 3 * dg_frag%nstate_frag * ppg%Nlma, dg_frag%icomm_frag)
        dg_frag%nl_projector_overlap(:, :, ispin, i_local) = proj_sum(:, :)
        dg_frag%nl_projector_r_overlap(:, :, :, ispin, i_local) = proj_r_sum(:, :, :)
      end do
    end do

    do i_halo = 1, dg_frag%n_halo
      do ispin = 1, nspin
        nbf_row = min(dg_frag%n_basis(dg_frag%halo(i_halo)%ifrag_src, ispin), dg_frag%nstate_frag)
        proj_local(:, :) = (0.0d0, 0.0d0)
        proj_sum(:, :) = (0.0d0, 0.0d0)
        proj_r_local(:, :, :) = (0.0d0, 0.0d0)
        proj_r_sum(:, :, :) = (0.0d0, 0.0d0)
        if (nbf_row > 0) then
          !$omp parallel do collapse(2) &
          !$omp& private(ilma, io, ia, j, ix, iy, iz, lx, ly, lz, x, y, z, phase, basis_val, overlap_i, overlap_r_i)
          do ilma = 1, ppg%Nlma
            do io = 1, nbf_row
              ia = ppg%ia_tbl(ilma)
              overlap_i = (0.0d0, 0.0d0)
              overlap_r_i(:) = (0.0d0, 0.0d0)
              do j = 1, ppg%mps(ia)
                ix = ppg%jxyz(1, j, ia)
                iy = ppg%jxyz(2, j, ia)
                iz = ppg%jxyz(3, j, ia)
                if (ix < is(1) .or. ix > ie(1) .or. &
                    iy < is(2) .or. iy > ie(2) .or. &
                    iz < is(3) .or. iz > ie(3)) cycle
                lx = map_global_to_halo_recv_buf_coord(dg_frag, dg_frag%halo(i_halo), 1, ix)
                ly = map_global_to_halo_recv_buf_coord(dg_frag, dg_frag%halo(i_halo), 2, iy)
                lz = map_global_to_halo_recv_buf_coord(dg_frag, dg_frag%halo(i_halo), 3, iz)
                if (lx < 1 .or. lx > dg_frag%halo(i_halo)%length(1)) cycle
                if (ly < 1 .or. ly > dg_frag%halo(i_halo)%length(2)) cycle
                if (lz < 1 .or. lz > dg_frag%halo(i_halo)%length(3)) cycle
                x = ppg%rxyz(1, j, ia)
                y = ppg%rxyz(2, j, ia)
                z = ppg%rxyz(3, j, ia)
                phase = Ac_tot(1) * x + Ac_tot(2) * y + Ac_tot(3) * z
                if (allocated(dg_frag%halo(i_halo)%buf_recv_c)) then
                  basis_val = dg_frag%halo(i_halo)%buf_recv_c(lx, ly, lz, io, 1)
                else
                  basis_val = cmplx(dg_frag%halo(i_halo)%buf_recv(lx, ly, lz, io, 1), 0.0d0, kind=8)
                end if
                overlap_i = overlap_i + basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
                overlap_r_i(1) = overlap_r_i(1) + x * basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
                overlap_r_i(2) = overlap_r_i(2) + y * basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
                overlap_r_i(3) = overlap_r_i(3) + z * basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
              end do
              proj_local(io, ilma) = overlap_i
              proj_r_local(1:3, io, ilma) = overlap_r_i(1:3)
            end do
          end do
          !$omp end parallel do
        end if
        call comm_summation(proj_local, proj_sum, dg_frag%nstate_frag * ppg%Nlma, dg_frag%icomm_frag)
        call comm_summation(proj_r_local, proj_r_sum, 3 * dg_frag%nstate_frag * ppg%Nlma, dg_frag%icomm_frag)
        dg_frag%nl_projector_overlap_halo(:, :, ispin, i_halo) = proj_sum(:, :)
        dg_frag%nl_projector_r_overlap_halo(:, :, :, ispin, i_halo) = proj_r_sum(:, :, :)
      end do
    end do

    deallocate(proj_r_local, proj_r_sum, proj_local, proj_sum)
    dg_frag%Ac_nl_projector_cache = Ac_tot
    dg_frag%nl_projector_cache_nlma = ppg%Nlma
    dg_frag%has_nl_projector_cache = .true.
  end subroutine ensure_nonlocal_projector_overlap_cache

  subroutine apply_nonlocal_projector_overlap_batch(dg_frag, mg, ppg, hvol, Ac_tot, ispin, x, y, row_frag_ids)
    use structures
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid),          intent(in) :: mg
    type(s_pp_grid),        intent(in) :: ppg
    real(8),                intent(in) :: hvol
    real(8),                intent(in) :: Ac_tot(3)
    integer,                intent(in) :: ispin
    complex(8),             intent(in) :: x(:, :)
    complex(8),             intent(inout) :: y(:, :)
    integer, optional,      intent(in) :: row_frag_ids(:)

    integer :: iblk_idx, iblk, irow_frag, ifrag, jfrag, i_local, nbf_row, nbf_col, nstate
    integer :: io, jo, ilma, istate, gid_i, gid_j, halo_slot, i_halo
    integer :: valid_row_count, valid_col_count
    integer, allocatable :: valid_row_basis(:), valid_row_gid(:)
    integer, allocatable :: valid_col_basis(:), valid_col_gid(:)
    complex(8), allocatable :: tmp_lma(:,:), uVcol(:,:)
    complex(8) :: pjl, pil, accum

    if (ppg%Nlma <= 0) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. dg_frag%has_nl_projector_cache) return
    if (.not. allocated(dg_frag%nl_projector_overlap)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return
    nstate = min(size(x, 2), size(y, 2))
    if (nstate <= 0) return

    if (.not. allocated(dg_frag%H_mat_blocks)) return

    allocate(valid_row_basis(max(1, size(dg_frag%index_basis, 1))))
    allocate(valid_row_gid(max(1, size(dg_frag%index_basis, 1))))
    allocate(valid_col_basis(max(1, size(dg_frag%index_basis, 1))))
    allocate(valid_col_gid(max(1, size(dg_frag%index_basis, 1))))
    allocate(tmp_lma(ppg%Nlma, nstate))
    allocate(uVcol(max(1, size(dg_frag%index_basis, 1)), ppg%Nlma))

    if (present(row_frag_ids)) then
      do irow_frag = 1, size(row_frag_ids)
        ifrag = row_frag_ids(irow_frag)
        call apply_one_fragment(ifrag)
      end do
    else
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        call apply_one_fragment(ifrag)
      end do
    end if

    deallocate(uVcol, tmp_lma, valid_col_gid, valid_col_basis, valid_row_gid, valid_row_basis)

  contains
    subroutine apply_one_fragment(ifrag_in)
      integer, intent(in) :: ifrag_in

      if (ifrag_in < dg_frag%ifrag_start .or. ifrag_in > dg_frag%ifrag_end) return
      i_local = ifrag_in - dg_frag%ifrag_start + 1
      if (i_local < 1 .or. i_local > size(dg_frag%nl_projector_overlap, 4)) return
      nbf_row = min(dg_frag%n_basis(ifrag_in, ispin), size(dg_frag%index_basis, 1), &
                size(dg_frag%nl_projector_overlap, 1))
      if (nbf_row <= 0) return

      valid_row_count = 0
      do io = 1, nbf_row
        gid_i = dg_frag%index_basis(io, ifrag_in, ispin)
        if (gid_i < 1 .or. gid_i > size(y, 1)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_basis(valid_row_count) = io
        valid_row_gid(valid_row_count) = gid_i
      end do
      if (valid_row_count <= 0) return

      if (allocated(dg_frag%H_local_block_ids)) then
        do iblk_idx = 1, size(dg_frag%H_local_block_ids)
          iblk = dg_frag%H_local_block_ids(iblk_idx)
          if (iblk < 1 .or. iblk > size(dg_frag%H_mat_blocks)) cycle
          if (dg_frag%H_mat_blocks(iblk)%ifrag_row /= ifrag_in) cycle
          call apply_one_column_fragment(dg_frag%H_mat_blocks(iblk)%ifrag_col)
        end do
      else
        do iblk = 1, size(dg_frag%H_mat_blocks)
          if (dg_frag%H_mat_blocks(iblk)%ifrag_row /= ifrag_in) cycle
          call apply_one_column_fragment(dg_frag%H_mat_blocks(iblk)%ifrag_col)
        end do
      end if
    end subroutine apply_one_fragment

    subroutine apply_one_column_fragment(jfrag_in)
      integer, intent(in) :: jfrag_in

      if (jfrag_in < 1 .or. jfrag_in > dg_frag%n_frag) return
      nbf_col = min(dg_frag%n_basis(jfrag_in, ispin), size(dg_frag%index_basis, 1), size(uVcol, 1))
      if (nbf_col <= 0) return

      valid_col_count = 0
      do jo = 1, nbf_col
        gid_j = dg_frag%index_basis(jo, jfrag_in, ispin)
        if (gid_j < 1 .or. gid_j > size(x, 1)) cycle
        valid_col_count = valid_col_count + 1
        valid_col_basis(valid_col_count) = jo
        valid_col_gid(valid_col_count) = gid_j
      end do
      if (valid_col_count <= 0) return

      uVcol(:, :) = (0.0d0, 0.0d0)
      if (jfrag_in == ifrag) then
        uVcol(1:nbf_col, 1:ppg%Nlma) = dg_frag%nl_projector_overlap(1:nbf_col, 1:ppg%Nlma, ispin, i_local)
      else
        halo_slot = 0
        do i_halo = 1, dg_frag%n_halo
          if (dg_frag%halo(i_halo)%ifrag_dst == ifrag .and. dg_frag%halo(i_halo)%ifrag_src == jfrag_in) then
            halo_slot = i_halo
            exit
          end if
        end do
        if (halo_slot <= 0) return
        if (.not. allocated(dg_frag%nl_projector_overlap_halo)) return
        if (halo_slot > size(dg_frag%nl_projector_overlap_halo, 4)) return
        uVcol(1:nbf_col, 1:ppg%Nlma) = dg_frag%nl_projector_overlap_halo(1:nbf_col, 1:ppg%Nlma, ispin, halo_slot)
      end if

      tmp_lma(:, :) = (0.0d0, 0.0d0)
!$omp parallel do private(ilma, istate, jo, gid_j, pjl) schedule(static)
      do istate = 1, nstate
        do ilma = 1, ppg%Nlma
          do jo = 1, valid_col_count
            gid_j = valid_col_gid(jo)
            pjl = uVcol(valid_col_basis(jo), ilma)
            tmp_lma(ilma, istate) = tmp_lma(ilma, istate) + pjl * x(gid_j, istate)
          end do
        end do
      end do
!$omp end parallel do

!$omp parallel do private(istate, io, gid_i, ilma, pil, accum) schedule(static)
      do istate = 1, nstate
        do io = 1, valid_row_count
          gid_i = valid_row_gid(io)
          accum = (0.0d0, 0.0d0)
          do ilma = 1, ppg%Nlma
            pil = dg_frag%nl_projector_overlap(valid_row_basis(io), ilma, ispin, i_local)
            accum = accum + conjg(pil) * ppg%rinv_uvu(ilma) * tmp_lma(ilma, istate)
          end do
          if (real(accum, kind=8) /= real(accum, kind=8) .or. &
              aimag(accum) /= aimag(accum)) then
            stop "NaN in DG nonlocal projector apply"
          end if
          y(gid_i, istate) = y(gid_i, istate) + accum
        end do
      end do
!$omp end parallel do
    end subroutine apply_one_column_fragment

  end subroutine apply_nonlocal_projector_overlap_batch

  subroutine build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, nspin, hvol, Ac_tot, use_micro_A, Ac_micro, &
       H_nl_blocks, H_nl_block_map, use_projector_cache)
    use structures
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in) :: mg
    type(s_pp_grid),        intent(in) :: ppg
    integer,                intent(in) :: nspin
    real(8),                intent(in) :: hvol
    real(8),                intent(in) :: Ac_tot(3)
    logical,                intent(in) :: use_micro_A
    real(8),                intent(in), optional :: Ac_micro(:, :, :, :)
    type(complex_matrix_block_info), intent(inout) :: H_nl_blocks(:)
    integer,                intent(in) :: H_nl_block_map(:, :)
    logical, optional,      intent(in) :: use_projector_cache

    integer :: ifrag, jfrag, ispin, io, jo, i_local, ilma, ia, j, ix, iy, iz
    integer :: nbf_row, nbf_col, iblk, halo_slot, i_halo
    integer :: lx, ly, lz
    integer :: is(3), ie(3)
    real(8) :: x, y, z, phase
    real(8) :: A_local(3)
    complex(8), allocatable :: uVrow(:,:), uVcol(:,:)
    complex(8) :: overlap_i, overlap_j, nlpp_contrib
    complex(8) :: basis_val
    logical :: use_cached_projectors

    if (ppg%Nlma == 0) return
    use_cached_projectors = .false.
    if (present(use_projector_cache)) use_cached_projectors = use_projector_cache
    use_cached_projectors = use_cached_projectors .and. (.not. use_micro_A) .and. &
                            dg_frag%has_nl_projector_cache .and. allocated(dg_frag%nl_projector_overlap)

    is = mg%is
    ie = mg%ie
    do iblk = 1, size(H_nl_blocks)
      H_nl_blocks(iblk)%val(:, :, :) = (0.0d0, 0.0d0)
    end do

    allocate(uVrow(dg_frag%nstate_frag, ppg%Nlma))
    allocate(uVcol(dg_frag%nstate_frag, ppg%Nlma))

    if (use_cached_projectors .and. .not. dg_frag%is_frag_root) then
      deallocate(uVrow, uVcol)
      return
    end if

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1

      do ispin = 1, nspin
        nbf_row = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
        if (nbf_row <= 0) cycle
        uVrow(:, :) = (0.0d0, 0.0d0)
        if (use_cached_projectors) then
          uVrow(1:nbf_row, 1:ppg%Nlma) = dg_frag%nl_projector_overlap(1:nbf_row, 1:ppg%Nlma, ispin, i_local)
        else
          !$omp parallel do collapse(2) &
          !$omp& private(ilma, io, ia, j, ix, iy, iz, lx, ly, lz, x, y, z, phase, A_local, basis_val, overlap_i)
          do ilma = 1, ppg%Nlma
            do io = 1, nbf_row
              ia = ppg%ia_tbl(ilma)
              overlap_i = (0.0d0, 0.0d0)
              do j = 1, ppg%mps(ia)
                ix = ppg%jxyz(1, j, ia)
                iy = ppg%jxyz(2, j, ia)
                iz = ppg%jxyz(3, j, ia)

                if (ix >= is(1) .and. ix <= ie(1) .and. &
                    iy >= is(2) .and. iy <= ie(2) .and. &
                    iz >= is(3) .and. iz <= ie(3)) then
                  lx = map_global_to_phi_box_coord(ix, lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1), &
                                                   dg_frag%lgnum_total(1))
                  ly = map_global_to_phi_box_coord(iy, lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2), &
                                                   dg_frag%lgnum_total(2))
                  lz = map_global_to_phi_box_coord(iz, lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3), &
                                                   dg_frag%lgnum_total(3))
                  if (lx == 0 .or. ly == 0 .or. lz == 0) cycle
                  x = ppg%rxyz(1, j, ia)
                  y = ppg%rxyz(2, j, ia)
                  z = ppg%rxyz(3, j, ia)
                  if (use_micro_A .and. present(Ac_micro)) then
                    A_local(1:3) = Ac_micro(1:3, ix, iy, iz)
                  else
                    A_local(1:3) = Ac_tot(1:3)
                  end if
                  phase = A_local(1) * x + A_local(2) * y + A_local(3) * z
                  if (allocated(dg_frag%phi_frag_c)) then
                    basis_val = dg_frag%phi_frag_c(lx, ly, lz, io, i_local)
                  else
                    basis_val = cmplx(dg_frag%phi_frag(lx, ly, lz, io, i_local), 0.0d0, kind=8)
                  end if
                  overlap_i = overlap_i + basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
                end if
              end do
              uVrow(io, ilma) = overlap_i
            end do
          end do
          !$omp end parallel do
        end if

        do iblk = 1, size(H_nl_blocks)
          if (H_nl_blocks(iblk)%ifrag_row /= ifrag) cycle
          jfrag = H_nl_blocks(iblk)%ifrag_col
          if (jfrag < 1 .or. jfrag > dg_frag%n_frag) cycle
          nbf_col = min(dg_frag%n_basis(jfrag, ispin), dg_frag%nstate_frag)
          if (nbf_col <= 0) cycle

          uVcol(:, :) = (0.0d0, 0.0d0)
          if (jfrag == ifrag) then
            uVcol(1:nbf_col, 1:ppg%Nlma) = uVrow(1:nbf_col, 1:ppg%Nlma)
          else
            halo_slot = 0
            do i_halo = 1, dg_frag%n_halo
              if (dg_frag%halo(i_halo)%ifrag_dst == ifrag .and. dg_frag%halo(i_halo)%ifrag_src == jfrag) then
                halo_slot = i_halo
                exit
              end if
            end do
            if (halo_slot <= 0) cycle

            !$omp parallel do collapse(2) &
            !$omp& private(ilma, jo, ia, j, ix, iy, iz, lx, ly, lz, x, y, z, phase, A_local, basis_val, overlap_j)
            do ilma = 1, ppg%Nlma
              do jo = 1, nbf_col
                ia = ppg%ia_tbl(ilma)
                overlap_j = (0.0d0, 0.0d0)
                do j = 1, ppg%mps(ia)
                  ix = ppg%jxyz(1, j, ia)
                  iy = ppg%jxyz(2, j, ia)
                  iz = ppg%jxyz(3, j, ia)
                  if (ix < is(1) .or. ix > ie(1) .or. iy < is(2) .or. iy > ie(2) .or. iz < is(3) .or. iz > ie(3)) cycle
                  lx = map_global_to_halo_recv_buf_coord(dg_frag, dg_frag%halo(halo_slot), 1, ix)
                  ly = map_global_to_halo_recv_buf_coord(dg_frag, dg_frag%halo(halo_slot), 2, iy)
                  lz = map_global_to_halo_recv_buf_coord(dg_frag, dg_frag%halo(halo_slot), 3, iz)
                  if (lx < 1 .or. lx > dg_frag%halo(halo_slot)%length(1)) cycle
                  if (ly < 1 .or. ly > dg_frag%halo(halo_slot)%length(2)) cycle
                  if (lz < 1 .or. lz > dg_frag%halo(halo_slot)%length(3)) cycle
                  x = ppg%rxyz(1, j, ia)
                  y = ppg%rxyz(2, j, ia)
                  z = ppg%rxyz(3, j, ia)
                  if (use_micro_A .and. present(Ac_micro)) then
                    A_local(1:3) = Ac_micro(1:3, ix, iy, iz)
                  else
                    A_local(1:3) = Ac_tot(1:3)
                  end if
                  phase = A_local(1) * x + A_local(2) * y + A_local(3) * z
                  if (allocated(dg_frag%halo(halo_slot)%buf_recv_c)) then
                    basis_val = dg_frag%halo(halo_slot)%buf_recv_c(lx, ly, lz, jo, 1)
                  else
                    basis_val = cmplx(dg_frag%halo(halo_slot)%buf_recv(lx, ly, lz, jo, 1), 0.0d0, kind=8)
                  end if
                  overlap_j = overlap_j + basis_val * ppg%uV(j, ilma) * exp(-zi * phase) * hvol
                end do
                uVcol(jo, ilma) = overlap_j
              end do
            end do
            !$omp end parallel do
          end if

          !$omp parallel do collapse(2) private(jo, io, ilma, nlpp_contrib, overlap_i, overlap_j)
          do jo = 1, nbf_col
            do io = 1, nbf_row
              nlpp_contrib = (0.0d0, 0.0d0)
              do ilma = 1, ppg%Nlma
                overlap_i = uVrow(io, ilma)
                overlap_j = uVcol(jo, ilma)
                nlpp_contrib = nlpp_contrib + conjg(overlap_i) * overlap_j * ppg%rinv_uvu(ilma)
              end do
              H_nl_blocks(iblk)%val(io, jo, ispin) = H_nl_blocks(iblk)%val(io, jo, ispin) + nlpp_contrib
            end do
          end do
          !$omp end parallel do
        end do
      end do
    end do

    deallocate(uVrow, uVcol)
  end subroutine build_nonlocal_pp_matrix_A_blocks

  !=======================================================================
  ! Ensure cached non-local PP matrix for current A(t)
  !=======================================================================
  subroutine ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_tot, require_dense_cache)
    use structures
    use salmon_global, only: ae_shape1, ae_shape2, theory
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(in)    :: Ac_tot(3)
    logical, optional,      intent(in)    :: require_dense_cache

    real(8) :: delta_A
    real(8) :: zero_A(3)
    logical :: reuse_allowed
    logical :: use_micro_A
    logical :: self_only_nl
    integer :: i, iblk
    type(complex_matrix_block_info), allocatable :: H_nl_ref_blocks(:)
    integer, allocatable :: H_nl_ref_block_map(:, :)
    logical, parameter :: enable_hmat_nl_progress = .false.

    if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        hmat-nl trace: stage=entry"
      flush(6)
    end if

    if (ppg%Nlma == 0 .or. .not. allocated(ppg%uV)) then
      if (allocated(dg_frag%H_nl_blocks)) then
        do i = 1, size(dg_frag%H_nl_blocks)
          if (allocated(dg_frag%H_nl_blocks(i)%val)) deallocate(dg_frag%H_nl_blocks(i)%val)
        end do
        deallocate(dg_frag%H_nl_blocks)
      end if
      if (allocated(dg_frag%H_nl_block_map)) deallocate(dg_frag%H_nl_block_map)
      if (allocated(dg_frag%H_nl_local_block_ids)) deallocate(dg_frag%H_nl_local_block_ids)
      dg_frag%n_H_nl_blocks = 0
      dg_frag%has_nl_cache = .false.
      return
    end if

    use_micro_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v))
    if (present(require_dense_cache)) then
      if (require_dense_cache) then
        write(*,'(1x,a)') "[FATAL] dense nonlocal-PP cache is disabled in DG block-sparse RT."
        stop "DG-Fragment RT: dense nonlocal cache requested"
      end if
    end if

    reuse_allowed = (.not. use_micro_A) .and. &
                    (trim(ae_shape1) == 'impulse' .or. trim(ae_shape1) == 'none') .and. &
                    (trim(ae_shape2) == 'impulse' .or. trim(ae_shape2) == 'none')
    self_only_nl = (dg_frag%dc_lcfo_seed_basis_cleaned .and. .not. dg_frag%identity_seed_coefficients)

    delta_A = maxval(abs(Ac_tot - dg_frag%Ac_nl_cache))
    if (.not. dg_frag%has_nl_cache .or. (.not. reuse_allowed) .or. delta_A > dg_frag%Ac_nl_cache_tol) then
      if (self_only_nl) then
        call init_complex_matrix_self_blocks_runtime(dg_frag, dg_frag%H_nl_blocks, dg_frag%H_nl_block_map)
        dg_frag%n_H_nl_blocks = size(dg_frag%H_nl_blocks)
      else if (.not. allocated(dg_frag%H_nl_blocks) .or. .not. allocated(dg_frag%H_nl_block_map)) then
        call init_complex_matrix_blocks_runtime(dg_frag, dg_frag%H_nl_blocks, dg_frag%H_nl_block_map)
        dg_frag%n_H_nl_blocks = size(dg_frag%H_nl_blocks)
      end if
      if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        hmat-nl trace: stage=after-init-blocks"
        flush(6)
      end if
      if (dg_frag%H_blocks_include_nonlocal) then
        zero_A(:) = 0.0d0
        call init_complex_matrix_self_blocks_runtime(dg_frag, H_nl_ref_blocks, H_nl_ref_block_map)
        if (use_micro_A) then
          call build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
               .true., system%Ac_micro%v, dg_frag%H_nl_blocks, dg_frag%H_nl_block_map)
        else
          call build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
               .false., H_nl_blocks=dg_frag%H_nl_blocks, H_nl_block_map=dg_frag%H_nl_block_map)
        end if
        call build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, system%nspin, system%hvol, zero_A, &
             .false., H_nl_blocks=H_nl_ref_blocks, H_nl_block_map=H_nl_ref_block_map)
        do iblk = 1, min(size(dg_frag%H_nl_blocks), size(H_nl_ref_blocks))
          if (allocated(dg_frag%H_nl_blocks(iblk)%val) .and. allocated(H_nl_ref_blocks(iblk)%val)) then
            dg_frag%H_nl_blocks(iblk)%val(:, :, :) = dg_frag%H_nl_blocks(iblk)%val(:, :, :) - &
                                                     H_nl_ref_blocks(iblk)%val(:, :, :)
          end if
        end do
        do iblk = 1, size(H_nl_ref_blocks)
          if (allocated(H_nl_ref_blocks(iblk)%val)) deallocate(H_nl_ref_blocks(iblk)%val)
        end do
        if (allocated(H_nl_ref_blocks)) deallocate(H_nl_ref_blocks)
        if (allocated(H_nl_ref_block_map)) deallocate(H_nl_ref_block_map)
      else if (use_micro_A) then
        call build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
             .true., system%Ac_micro%v, dg_frag%H_nl_blocks, dg_frag%H_nl_block_map)
      else if (self_only_nl) then
        call ensure_nonlocal_projector_overlap_cache(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot)
        call build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
             .false., H_nl_blocks=dg_frag%H_nl_blocks, H_nl_block_map=dg_frag%H_nl_block_map, &
             use_projector_cache=.true.)
      else
        call build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
             .false., H_nl_blocks=dg_frag%H_nl_blocks, H_nl_block_map=dg_frag%H_nl_block_map)
      end if
      if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        hmat-nl trace: stage=after-build-blocks"
        flush(6)
        write(*,'(1x,a)') "        hmat-nl trace: stage=before-block-reduce"
        flush(6)
      end if
      call reduce_complex_matrix_blocks_runtime(dg_frag, dg_frag%H_nl_blocks, "hmat-nl", dg_frag%icomm)
      call rebuild_local_nl_block_ids(dg_frag)
      if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        hmat-nl trace: stage=after-block-reduce"
        flush(6)
      end if
      dg_frag%Ac_nl_cache = Ac_tot
      dg_frag%has_nl_cache = .true.
    end if

  end subroutine ensure_nonlocal_pp_matrix_A

  integer function map_global_to_phi_box_coord(ig, lb, ub, lgtot) result(iloc)
    implicit none
    integer, intent(in) :: ig, lb, ub, lgtot

    iloc = modulo(ig - 1, lgtot) + 1
    if (iloc < lb) then
      iloc = iloc + ((lb - iloc + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc > ub) then
      iloc = iloc - ((iloc - ub + lgtot - 1) / lgtot) * lgtot
    end if
    if (iloc < lb .or. iloc > ub) iloc = 0
  end function map_global_to_phi_box_coord

  integer function map_global_to_halo_recv_buf_coord(dg_frag, halo, axis, ig) result(ibuf)
    use rt_dg_fragment_types, only: halo_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(halo_info), intent(in) :: halo
    integer, intent(in) :: axis, ig

    integer :: local_idx

    local_idx = map_global_to_phi_box_coord(ig, lbound(dg_frag%phi_frag, axis), ubound(dg_frag%phi_frag, axis), &
                                            dg_frag%lgnum_total(axis))
    if (local_idx == 0) then
      ibuf = 0
      return
    end if

    if (local_idx < halo%recv_lo(axis) .or. local_idx > halo%recv_hi(axis)) then
      ibuf = 0
      return
    end if

    ibuf = local_idx - halo%recv_lo(axis) + 1
  end function map_global_to_halo_recv_buf_coord

  subroutine ensure_overlap_prop_available(dg_frag, n_use)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: n_use

    if (n_use <= 0) return
  end subroutine ensure_overlap_prop_available

  subroutine apply_mixed_hamiltonian(dg_frag, ispin, x, y)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:, :)
    complex(8), intent(inout) :: y(:, :)

    integer :: n_frag, n_pw, n_tot, ipw

    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    y(:, :) = (0.0d0, 0.0d0)

    if (n_frag > 0) then
      if (allocated(dg_frag%H_mat_blocks)) then
        call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, x(1:n_frag, :), y(1:n_frag, :))
      else if (allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)) then
        y(1:n_frag, :) = matmul(dg_frag%H_mat_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag, :))
      else if (allocated(dg_frag%H_mat)) then
        y(1:n_frag, :) = matmul(cmplx(dg_frag%H_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag, :))
      end if
    end if

    if (n_pw > 0 .and. allocated(dg_frag%H_mat_frag_pw)) then
      y(1:n_frag, :) = y(1:n_frag, :) + matmul(dg_frag%H_mat_frag_pw(1:n_frag, 1:n_pw, ispin), x(n_frag+1:n_tot, :))
      y(n_frag+1:n_tot, :) = y(n_frag+1:n_tot, :) + &
        matmul(conjg(transpose(dg_frag%H_mat_frag_pw(1:n_frag, 1:n_pw, ispin))), x(1:n_frag, :))
    end if

    if (n_pw > 0) then
      if (allocated(dg_frag%H_mat_pw)) then
        y(n_frag+1:n_tot, :) = y(n_frag+1:n_tot, :) + &
          matmul(dg_frag%H_mat_pw(1:n_pw, 1:n_pw, ispin), x(n_frag+1:n_tot, :))
      else if (allocated(dg_frag%H_mat_pw_diag)) then
        do ipw = 1, n_pw
          y(n_frag+ipw, :) = y(n_frag+ipw, :) + dg_frag%H_mat_pw_diag(ipw, ispin) * x(n_frag+ipw, :)
        end do
      else
        do ipw = 1, n_pw
          y(n_frag+ipw, :) = y(n_frag+ipw, :) + 0.5d0 * sum(dg_frag%k_pw(:, ipw)**2) * x(n_frag+ipw, :)
        end do
      end if
    end if
  end subroutine apply_mixed_hamiltonian

  subroutine apply_mixed_hamiltonian_local_rows(dg_frag, ispin, frag_row_ids, pw_row_ids, &
                                                coef_frag, coef_pw, y_frag, y_pw)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: frag_row_ids(:), pw_row_ids(:)
    complex(8), intent(in) :: coef_frag(:, :), coef_pw(:, :)
    complex(8), intent(inout) :: y_frag(:, :), y_pw(:, :)

    integer :: irow_local, ipw_local, jpw_local, frag_row, pw_row, pw_col, frag_pos, pw_pos, col_pos
    integer :: istate, nstate
    integer, allocatable :: frag_pos_map(:), pw_pos_map(:)
    complex(8) :: hfp
    logical :: pw_row_owned

    if (.not. allocated(dg_frag%H_mat_frag_pw_local)) return
    if (.not. allocated(dg_frag%fp_local_row_ids)) return
    if (.not. allocated(dg_frag%fp_local_pw_ids)) return
    if (ispin < 1 .or. ispin > size(dg_frag%H_mat_frag_pw_local, 3)) return

    nstate = min(size(y_frag, 2), size(y_pw, 2), size(coef_frag, 2), size(coef_pw, 2))
    if (nstate <= 0) return
    if (size(frag_row_ids) > 0) then
      allocate(frag_pos_map(max(1, dg_frag%n_mat_max)))
      frag_pos_map(:) = 0
      do frag_pos = 1, min(size(frag_row_ids), size(y_frag, 1), size(coef_frag, 1))
        frag_row = frag_row_ids(frag_pos)
        if (frag_row < 1 .or. frag_row > size(frag_pos_map)) cycle
        frag_pos_map(frag_row) = frag_pos
      end do
    else
      allocate(frag_pos_map(1))
      frag_pos_map(:) = 0
    end if

    if (size(pw_row_ids) > 0) then
      allocate(pw_pos_map(max(1, dg_frag%n_plane_waves)))
      pw_pos_map(:) = 0
      do pw_pos = 1, min(size(pw_row_ids), size(y_pw, 1), size(coef_pw, 1))
        pw_row = pw_row_ids(pw_pos)
        if (pw_row < 1 .or. pw_row > size(pw_pos_map)) cycle
        pw_pos_map(pw_row) = pw_pos
      end do
    else
      allocate(pw_pos_map(1))
      pw_pos_map(:) = 0
    end if

    do ipw_local = 1, min(size(dg_frag%fp_local_pw_ids), size(dg_frag%H_mat_frag_pw_local, 2))
      pw_row = dg_frag%fp_local_pw_ids(ipw_local)
      if (pw_row < 1 .or. pw_row > size(pw_pos_map)) cycle
      pw_pos = pw_pos_map(pw_row)
      if (pw_pos <= 0) cycle
      pw_row_owned = .true.
      if (allocated(dg_frag%coef_pw_owner)) then
        pw_row_owned = (pw_row <= size(dg_frag%coef_pw_owner) .and. dg_frag%coef_pw_owner(pw_row) == dg_frag%id)
      end if
      if (pw_row_owned .and. allocated(dg_frag%H_mat_pw)) then
        if (pw_row <= size(dg_frag%H_mat_pw, 1) .and. ispin <= size(dg_frag%H_mat_pw, 3)) then
          do jpw_local = 1, size(pw_row_ids)
            pw_col = pw_row_ids(jpw_local)
            if (pw_col < 1 .or. pw_col > size(pw_pos_map)) cycle
            col_pos = pw_pos_map(pw_col)
            if (col_pos <= 0) cycle
            if (pw_col > size(dg_frag%H_mat_pw, 2)) cycle
            hfp = dg_frag%H_mat_pw(pw_row, pw_col, ispin)
            if (abs(hfp) == 0.0d0) cycle
            do istate = 1, nstate
              y_pw(pw_pos, istate) = y_pw(pw_pos, istate) + hfp * coef_pw(col_pos, istate)
            end do
          end do
        end if
      else if (pw_row_owned .and. allocated(dg_frag%H_mat_pw_diag)) then
        if (pw_row <= size(dg_frag%H_mat_pw_diag, 1) .and. ispin <= size(dg_frag%H_mat_pw_diag, 2)) then
          do istate = 1, nstate
            y_pw(pw_pos, istate) = y_pw(pw_pos, istate) + dg_frag%H_mat_pw_diag(pw_row, ispin) * coef_pw(pw_pos, istate)
          end do
        end if
      else if (pw_row_owned .and. allocated(dg_frag%k_pw)) then
        if (pw_row <= size(dg_frag%k_pw, 2)) then
          do istate = 1, nstate
            y_pw(pw_pos, istate) = y_pw(pw_pos, istate) + 0.5d0 * sum(dg_frag%k_pw(:, pw_row)**2) * coef_pw(pw_pos, istate)
          end do
        end if
      end if
      do irow_local = 1, min(size(dg_frag%fp_local_row_ids), size(dg_frag%H_mat_frag_pw_local, 1))
        frag_row = dg_frag%fp_local_row_ids(irow_local)
        if (frag_row < 1 .or. frag_row > size(frag_pos_map)) cycle
        frag_pos = frag_pos_map(frag_row)
        if (frag_pos <= 0) cycle
        hfp = dg_frag%H_mat_frag_pw_local(irow_local, ipw_local, ispin)
        if (abs(hfp) == 0.0d0) cycle
        do istate = 1, nstate
          y_frag(frag_pos, istate) = y_frag(frag_pos, istate) + hfp * coef_pw(pw_pos, istate)
        end do
        if (.not. pw_row_owned) cycle
        do istate = 1, nstate
          y_pw(pw_pos, istate) = y_pw(pw_pos, istate) + conjg(hfp) * coef_frag(frag_pos, istate)
        end do
      end do
    end do

    deallocate(frag_pos_map, pw_pos_map)
  end subroutine apply_mixed_hamiltonian_local_rows

  subroutine apply_overlap_operator(dg_frag, ispin, x, y, use_prop)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:)
    complex(8), intent(inout) :: y(:)
    logical, intent(in) :: use_prop

    integer :: n_frag, n_pw, n_tot

    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    if (size(x, 1) < n_tot .or. size(y, 1) < n_tot) then
      n_pw = 0
      n_tot = n_frag
    end if
    y(:) = (0.0d0, 0.0d0)

    if (dg_frag%dc_lcfo_seed_basis_cleaned .and. .not. dg_frag%identity_seed_coefficients .and. n_pw == 0) then
      y(1:n_frag) = x(1:n_frag)
      return
    end if

    if (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin)) then
      if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks)) then
        call apply_matrix_blocks(dg_frag, dg_frag%S_mat_prop_blocks, ispin, x(1:n_frag), y(1:n_frag))
      else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_blocks)) then
        call apply_matrix_blocks(dg_frag, dg_frag%S_mat_blocks, ispin, x(1:n_frag), y(1:n_frag))
      else if (use_prop .and. allocated(dg_frag%S_mat_prop_c)) then
        y(1:n_frag) = matmul(dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag))
      else if (use_prop .and. allocated(dg_frag%S_mat_prop)) then
        y(1:n_frag) = matmul(cmplx(dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag))
      else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_c)) then
        y(1:n_frag) = matmul(dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag))
      else if (allocated(dg_frag%S_mat)) then
        y(1:n_frag) = matmul(cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag))
      end if
      y(1:n_frag) = y(1:n_frag) + matmul(dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin), x(n_frag+1:n_tot))
      y(n_frag+1:n_tot) = x(n_frag+1:n_tot) + matmul(conjg(transpose(dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin))), x(1:n_frag))
    else if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks) .and. n_pw == 0) then
      call apply_matrix_blocks(dg_frag, dg_frag%S_mat_prop_blocks, ispin, x(1:n_frag), y(1:n_frag))
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_blocks) .and. n_pw == 0) then
      call apply_matrix_blocks(dg_frag, dg_frag%S_mat_blocks, ispin, x(1:n_frag), y(1:n_frag))
    else if (use_prop .and. allocated(dg_frag%S_mat_prop_c)) then
      y(1:n_frag) = matmul(dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag))
    else if (use_prop .and. allocated(dg_frag%S_mat_prop)) then
      y(1:n_frag) = matmul(cmplx(dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag))
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_c)) then
      y(1:n_frag) = matmul(dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag))
    else if (allocated(dg_frag%S_mat)) then
      y(1:n_frag) = matmul(cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag))
    end if
    if (n_pw > 0) then
      y(n_frag+1:n_tot) = x(n_frag+1:n_tot)
    end if
  end subroutine apply_overlap_operator

  subroutine apply_overlap_operator_batch(dg_frag, ispin, x, y, use_prop)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:, :)
    complex(8), intent(inout) :: y(:, :)
    logical, intent(in) :: use_prop

    integer :: n_frag, n_pw, n_tot

    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    if (size(x, 1) < n_tot .or. size(y, 1) < n_tot) then
      n_pw = 0
      n_tot = n_frag
    end if
    y(:, :) = (0.0d0, 0.0d0)

    if (dg_frag%dc_lcfo_seed_basis_cleaned .and. .not. dg_frag%identity_seed_coefficients .and. n_pw == 0) then
      y(1:n_frag, :) = x(1:n_frag, :)
      return
    end if

    if (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin)) then
      if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks)) then
        call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_prop_blocks, ispin, x(1:n_frag, :), y(1:n_frag, :))
      else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_blocks)) then
        call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, x(1:n_frag, :), y(1:n_frag, :))
      else if (use_prop .and. allocated(dg_frag%S_mat_prop_c)) then
        y(1:n_frag, :) = matmul(dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag, :))
      else if (use_prop .and. allocated(dg_frag%S_mat_prop)) then
        y(1:n_frag, :) = matmul(cmplx(dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag, :))
      else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_c)) then
        y(1:n_frag, :) = matmul(dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag, :))
      else if (allocated(dg_frag%S_mat)) then
        y(1:n_frag, :) = matmul(cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag, :))
      end if
      y(1:n_frag, :) = y(1:n_frag, :) + matmul(dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin), x(n_frag+1:n_tot, :))
      y(n_frag+1:n_tot, :) = x(n_frag+1:n_tot, :) + &
        matmul(conjg(transpose(dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin))), x(1:n_frag, :))
    else if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks) .and. n_pw == 0) then
      call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_prop_blocks, ispin, x(1:n_frag, :), y(1:n_frag, :))
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_blocks) .and. n_pw == 0) then
      call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, x(1:n_frag, :), y(1:n_frag, :))
    else if (use_prop .and. allocated(dg_frag%S_mat_prop_c)) then
      y(1:n_frag, :) = matmul(dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag, :))
    else if (use_prop .and. allocated(dg_frag%S_mat_prop)) then
      y(1:n_frag, :) = matmul(cmplx(dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag, :))
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_c)) then
      y(1:n_frag, :) = matmul(dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin), x(1:n_frag, :))
    else if (allocated(dg_frag%S_mat)) then
      y(1:n_frag, :) = matmul(cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8), x(1:n_frag, :))
    end if
    if (n_pw > 0) then
      y(n_frag+1:n_tot, :) = x(n_frag+1:n_tot, :)
    end if
  end subroutine apply_overlap_operator_batch

  subroutine build_overlap_operator_diagonal(dg_frag, ispin, use_prop, diag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    logical, intent(in) :: use_prop
    real(8), intent(out) :: diag(:)

    integer :: n_frag, n_pw, n_tot
    integer :: ifrag, iblk, ib, ig, nbf, nblock_diag

    diag(:) = 0.0d0
    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    if (size(diag) < n_tot) then
      n_pw = 0
      n_tot = n_frag
    end if
    if (size(diag) < n_tot) return

    if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks)) then
      do ifrag = 1, dg_frag%n_frag
        iblk = find_matrix_block_runtime(dg_frag%S_block_map, ifrag, ifrag)
        if (iblk <= 0 .or. iblk > size(dg_frag%S_mat_prop_blocks)) cycle
        nbf = dg_frag%n_basis(ifrag, ispin)
        nblock_diag = min(nbf, size(dg_frag%S_mat_prop_blocks(iblk)%val, 1), size(dg_frag%S_mat_prop_blocks(iblk)%val, 2))
        do ib = 1, nblock_diag
          ig = dg_frag%index_basis(ib, ifrag, ispin)
          if (ig < 1 .or. ig > size(diag)) cycle
          diag(ig) = real(dg_frag%S_mat_prop_blocks(iblk)%val(ib, ib, ispin), kind=8)
        end do
      end do
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_blocks)) then
      do ifrag = 1, dg_frag%n_frag
        iblk = find_matrix_block_runtime(dg_frag%S_block_map, ifrag, ifrag)
        if (iblk <= 0 .or. iblk > size(dg_frag%S_mat_blocks)) cycle
        nbf = dg_frag%n_basis(ifrag, ispin)
        nblock_diag = min(nbf, size(dg_frag%S_mat_blocks(iblk)%val, 1), size(dg_frag%S_mat_blocks(iblk)%val, 2))
        do ib = 1, nblock_diag
          ig = dg_frag%index_basis(ib, ifrag, ispin)
          if (ig < 1 .or. ig > size(diag)) cycle
          diag(ig) = real(dg_frag%S_mat_blocks(iblk)%val(ib, ib, ispin), kind=8)
        end do
      end do
    else if (use_prop .and. allocated(dg_frag%S_mat_prop_c)) then
      diag(1:n_frag) = real([(dg_frag%S_mat_prop_c(ib, ib, ispin), ib=1,n_frag)], kind=8)
    else if (use_prop .and. allocated(dg_frag%S_mat_prop)) then
      diag(1:n_frag) = [(dg_frag%S_mat_prop(ib, ib, ispin), ib=1,n_frag)]
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_c)) then
      diag(1:n_frag) = real([(dg_frag%S_mat_c(ib, ib, ispin), ib=1,n_frag)], kind=8)
    else if (allocated(dg_frag%S_mat)) then
      diag(1:n_frag) = [(dg_frag%S_mat(ib, ib, ispin), ib=1,n_frag)]
    end if

    if (n_pw > 0) diag(n_frag+1:n_tot) = 1.0d0
  end subroutine build_overlap_operator_diagonal

  subroutine solve_overlap_operator_batch(dg_frag, ispin, rhs, sol, use_prop)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: rhs(:, :)
    complex(8), intent(out) :: sol(:, :)
    logical, intent(in) :: use_prop

    integer :: n_dim, n_rhs, icol, iter, max_iter
    integer :: irow
    real(8), parameter :: diag_floor = 1.0d-10, tol_rel = 1.0d-10
    complex(8), allocatable :: r(:,:), z(:,:), p(:,:), p_full(:,:), ap(:,:)
    real(8), allocatable :: diag(:), tol_abs2(:), rho(:), rho_new(:), denom(:)
    real(8), allocatable :: dot_local(:), dot_global(:), res_local(:), res_global(:)
    complex(8), allocatable :: alpha(:), beta(:)
    logical, allocatable :: active(:)
    logical, allocatable :: owned(:)
    logical :: distributed_rows

    n_dim = size(rhs, 1)
    n_rhs = size(rhs, 2)
    sol(:, :) = (0.0d0, 0.0d0)
    if (n_dim <= 0 .or. n_rhs <= 0) return

    allocate(diag(n_dim), r(n_dim, n_rhs), z(n_dim, n_rhs), p(n_dim, n_rhs), p_full(n_dim, n_rhs), ap(n_dim, n_rhs))
    allocate(tol_abs2(n_rhs), rho(n_rhs), rho_new(n_rhs), denom(n_rhs))
    allocate(dot_local(n_rhs), dot_global(n_rhs), res_local(n_rhs), res_global(n_rhs))
    allocate(alpha(n_rhs), beta(n_rhs), active(n_rhs))
    allocate(owned(n_dim))

    distributed_rows = allocated(dg_frag%coef_owner)
    if (distributed_rows) then
      distributed_rows = (ispin >= 1 .and. ispin <= size(dg_frag%coef_owner, 2))
    end if

    owned(:) = .true.
    if (distributed_rows) then
      owned(:) = .false.
      do irow = 1, min(n_dim, size(dg_frag%coef_owner, 1))
        owned(irow) = (dg_frag%coef_owner(irow, ispin) == dg_frag%id)
      end do
    end if

    call build_overlap_operator_diagonal(dg_frag, ispin, use_prop, diag)
    where (diag < diag_floor) diag = diag_floor

    max_iter = max(20, min(6 * max(1, n_dim), 400))

    sol(:, :) = (0.0d0, 0.0d0)
    do icol = 1, n_rhs
      do irow = 1, n_dim
        if (owned(irow)) sol(irow, icol) = rhs(irow, icol)
      end do
    end do

    if (distributed_rows) then
      call comm_summation(sol, p_full, n_dim * n_rhs, dg_frag%icomm)
    else
      p_full(:, :) = sol(:, :)
    end if
    call apply_overlap_operator_batch(dg_frag, ispin, p_full, ap, use_prop)
    r(:, :) = (0.0d0, 0.0d0)
    dot_local(:) = 0.0d0
    res_local(:) = 0.0d0
    do icol = 1, n_rhs
      do irow = 1, n_dim
        if (.not. owned(irow)) cycle
        r(irow, icol) = rhs(irow, icol) - ap(irow, icol)
        dot_local(icol) = dot_local(icol) + real(conjg(rhs(irow, icol)) * rhs(irow, icol), kind=8)
        res_local(icol) = res_local(icol) + real(conjg(r(irow, icol)) * r(irow, icol), kind=8)
      end do
    end do
    if (distributed_rows) then
      call comm_summation(dot_local, dot_global, n_rhs, dg_frag%icomm)
      call comm_summation(res_local, res_global, n_rhs, dg_frag%icomm)
    else
      dot_global(:) = dot_local(:)
      res_global(:) = res_local(:)
    end if

    do icol = 1, n_rhs
      tol_abs2(icol) = max(1.0d-24, (tol_rel * max(sqrt(max(0.0d0, dot_global(icol))), 1.0d0))**2)
      active(icol) = (res_global(icol) > tol_abs2(icol))
    end do

    if (.not. any(active)) then
      deallocate(diag, r, z, p, p_full, ap, tol_abs2, rho, rho_new, denom)
      deallocate(dot_local, dot_global, res_local, res_global, alpha, beta, active, owned)
      return
    end if

    z(:, :) = (0.0d0, 0.0d0)
    p(:, :) = (0.0d0, 0.0d0)
    dot_local(:) = 0.0d0
    do icol = 1, n_rhs
      if (.not. active(icol)) cycle
      do irow = 1, n_dim
        if (.not. owned(irow)) cycle
        z(irow, icol) = r(irow, icol) / diag(irow)
        p(irow, icol) = z(irow, icol)
        dot_local(icol) = dot_local(icol) + real(conjg(r(irow, icol)) * z(irow, icol), kind=8)
      end do
    end do
    if (distributed_rows) then
      call comm_summation(dot_local, dot_global, n_rhs, dg_frag%icomm)
    else
      dot_global(:) = dot_local(:)
    end if
    rho(:) = dot_global(:)
    do icol = 1, n_rhs
      if (.not. active(icol)) cycle
      if (rho(icol) <= 0.0d0) then
        active(icol) = .false.
      end if
    end do

    do iter = 1, max_iter
      if (.not. any(active)) exit
      do icol = 1, n_rhs
        if (.not. active(icol)) p(:, icol) = (0.0d0, 0.0d0)
      end do

      if (distributed_rows) then
        call comm_summation(p, p_full, n_dim * n_rhs, dg_frag%icomm)
      else
        p_full(:, :) = p(:, :)
      end if
      call apply_overlap_operator_batch(dg_frag, ispin, p_full, ap, use_prop)

      dot_local(:) = 0.0d0
      do icol = 1, n_rhs
        if (.not. active(icol)) cycle
        do irow = 1, n_dim
          if (.not. owned(irow)) cycle
          dot_local(icol) = dot_local(icol) + real(conjg(p(irow, icol)) * ap(irow, icol), kind=8)
        end do
      end do
      if (distributed_rows) then
        call comm_summation(dot_local, dot_global, n_rhs, dg_frag%icomm)
      else
        dot_global(:) = dot_local(:)
      end if
      denom(:) = dot_global(:)

      res_local(:) = 0.0d0
      do icol = 1, n_rhs
        if (.not. active(icol)) cycle
        if (abs(denom(icol)) <= 1.0d-30) then
          active(icol) = .false.
          cycle
        end if
        alpha(icol) = cmplx(rho(icol) / denom(icol), 0.0d0, kind=8)
        do irow = 1, n_dim
          if (.not. owned(irow)) cycle
          sol(irow, icol) = sol(irow, icol) + alpha(icol) * p(irow, icol)
          r(irow, icol) = r(irow, icol) - alpha(icol) * ap(irow, icol)
          res_local(icol) = res_local(icol) + real(conjg(r(irow, icol)) * r(irow, icol), kind=8)
        end do
      end do
      if (distributed_rows) then
        call comm_summation(res_local, res_global, n_rhs, dg_frag%icomm)
      else
        res_global(:) = res_local(:)
      end if

      dot_local(:) = 0.0d0
      rho_new(:) = 0.0d0
      do icol = 1, n_rhs
        if (.not. active(icol)) cycle
        if (res_global(icol) <= tol_abs2(icol)) then
          active(icol) = .false.
          p(:, icol) = (0.0d0, 0.0d0)
        else
          do irow = 1, n_dim
            if (.not. owned(irow)) cycle
            z(irow, icol) = r(irow, icol) / diag(irow)
            dot_local(icol) = dot_local(icol) + real(conjg(r(irow, icol)) * z(irow, icol), kind=8)
          end do
        end if
      end do
      if (distributed_rows) then
        call comm_summation(dot_local, dot_global, n_rhs, dg_frag%icomm)
      else
        dot_global(:) = dot_local(:)
      end if
      rho_new(:) = dot_global(:)

      do icol = 1, n_rhs
        if (.not. active(icol)) cycle
        if (rho(icol) <= 0.0d0) then
          active(icol) = .false.
          cycle
        end if
        if (rho_new(icol) <= 0.0d0) then
          active(icol) = .false.
          cycle
        end if
        beta(icol) = cmplx(rho_new(icol) / rho(icol), 0.0d0, kind=8)
        do irow = 1, n_dim
          if (.not. owned(irow)) cycle
          p(irow, icol) = z(irow, icol) + beta(icol) * p(irow, icol)
        end do
        rho(icol) = rho_new(icol)
      end do
    end do

    deallocate(diag, r, z, p, p_full, ap, tol_abs2, rho, rho_new, denom)
    deallocate(dot_local, dot_global, res_local, res_global, alpha, beta, active, owned)
  end subroutine solve_overlap_operator_batch

  subroutine solve_overlap_operator_batch_local(dg_frag, ispin, rhs, sol, use_prop)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: rhs(:, :)
    complex(8), intent(out) :: sol(:, :)
    logical, intent(in) :: use_prop

    integer :: n_dim, n_rhs, n_loc, ifrag, iblk, ii, jj, ig_i, ig_j, info, n_frag, n_pw, ipw
    integer :: idx_ii, idx_jj, nrow, ncol, valid_row_count, valid_col_count
    integer, allocatable :: row_ids(:), local_pos(:), ipiv(:)
    logical, allocatable :: frag_selected(:)
    complex(8), allocatable :: mat_loc(:,:), rhs_loc(:,:)
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))

    n_dim = size(rhs, 1)
    n_rhs = size(rhs, 2)
    sol(:, :) = (0.0d0, 0.0d0)
    if (n_dim <= 0 .or. n_rhs <= 0) return
    if (.not. allocated(dg_frag%H_local_rows)) return
    if (size(dg_frag%H_local_rows) <= 0) return
    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    if (n_dim < n_frag) return

    allocate(frag_selected(dg_frag%n_frag))
    frag_selected = .false.
    do ifrag = 1, size(dg_frag%H_local_rows)
      if (dg_frag%H_local_rows(ifrag) < 1 .or. dg_frag%H_local_rows(ifrag) > dg_frag%n_frag) cycle
      frag_selected(dg_frag%H_local_rows(ifrag)) = .true.
    end do

    if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks)) then
      do iblk = 1, size(dg_frag%S_mat_prop_blocks)
        if (frag_selected(dg_frag%S_mat_prop_blocks(iblk)%ifrag_row) .or. &
            frag_selected(dg_frag%S_mat_prop_blocks(iblk)%ifrag_col)) then
          frag_selected(dg_frag%S_mat_prop_blocks(iblk)%ifrag_row) = .true.
          frag_selected(dg_frag%S_mat_prop_blocks(iblk)%ifrag_col) = .true.
        end if
      end do
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_blocks)) then
      do iblk = 1, size(dg_frag%S_mat_blocks)
        if (frag_selected(dg_frag%S_mat_blocks(iblk)%ifrag_row) .or. &
            frag_selected(dg_frag%S_mat_blocks(iblk)%ifrag_col)) then
          frag_selected(dg_frag%S_mat_blocks(iblk)%ifrag_row) = .true.
          frag_selected(dg_frag%S_mat_blocks(iblk)%ifrag_col) = .true.
        end if
      end do
    else
      deallocate(frag_selected)
      return
    end if

    allocate(local_pos(n_dim))
    local_pos = 0
    n_loc = 0
    do ifrag = 1, dg_frag%n_frag
      if (.not. frag_selected(ifrag)) cycle
      do ii = 1, dg_frag%n_basis(ifrag, ispin)
        ig_i = dg_frag%index_basis(ii, ifrag, ispin)
        if (ig_i < 1 .or. ig_i > n_dim) cycle
        if (local_pos(ig_i) /= 0) cycle
        n_loc = n_loc + 1
        local_pos(ig_i) = n_loc
      end do
    end do
    if (n_loc <= 0) then
      deallocate(frag_selected, local_pos)
      return
    end if

    allocate(row_ids(n_loc))
    do ig_i = 1, n_dim
      if (local_pos(ig_i) > 0) row_ids(local_pos(ig_i)) = ig_i
    end do

    allocate(mat_loc(n_loc, n_loc), rhs_loc(n_loc, n_rhs), ipiv(n_loc))
    mat_loc(:, :) = (0.0d0, 0.0d0)
    rhs_loc(:, :) = rhs(row_ids(:), :)

    if (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin)) then
      do ipw = 1, n_pw
        do jj = 1, n_rhs
          rhs_loc(:, jj) = rhs_loc(:, jj) - dg_frag%S_mat_frag_pw(row_ids(:), ipw, ispin) * rhs(n_frag + ipw, jj)
        end do
      end do
      sol(n_frag+1:n_frag+n_pw, :) = rhs(n_frag+1:n_frag+n_pw, :)
    end if

    if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks)) then
!$omp parallel do schedule(static) &
!$omp& private(iblk, nrow, ncol, valid_row_count, valid_col_count, ii, jj, idx_ii, idx_jj, &
!$omp&         ig_i, ig_j, row_gid, col_gid, valid_row_ids, valid_col_ids)
      do iblk = 1, size(dg_frag%S_mat_prop_blocks)
        if (.not. frag_selected(dg_frag%S_mat_prop_blocks(iblk)%ifrag_row)) cycle
        if (.not. frag_selected(dg_frag%S_mat_prop_blocks(iblk)%ifrag_col)) cycle
        nrow = dg_frag%n_basis(dg_frag%S_mat_prop_blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(dg_frag%S_mat_prop_blocks(iblk)%ifrag_col, ispin)
        valid_row_count = 0
        do ii = 1, nrow
          row_gid(ii) = dg_frag%index_basis(ii, dg_frag%S_mat_prop_blocks(iblk)%ifrag_row, ispin)
          if (row_gid(ii) < 1 .or. row_gid(ii) > n_dim) cycle
          if (local_pos(row_gid(ii)) <= 0) cycle
          valid_row_count = valid_row_count + 1
          valid_row_ids(valid_row_count) = ii
        end do
        if (valid_row_count <= 0) cycle
        valid_col_count = 0
        do jj = 1, ncol
          col_gid(jj) = dg_frag%index_basis(jj, dg_frag%S_mat_prop_blocks(iblk)%ifrag_col, ispin)
          if (col_gid(jj) < 1 .or. col_gid(jj) > n_dim) cycle
          if (local_pos(col_gid(jj)) <= 0) cycle
          valid_col_count = valid_col_count + 1
          valid_col_ids(valid_col_count) = jj
        end do
        if (valid_col_count <= 0) cycle
        do idx_jj = 1, valid_col_count
          jj = valid_col_ids(idx_jj)
          ig_j = col_gid(jj)
!$omp simd private(ii,ig_i)
          do idx_ii = 1, valid_row_count
            ii = valid_row_ids(idx_ii)
            ig_i = row_gid(ii)
            mat_loc(local_pos(ig_i), local_pos(ig_j)) = dg_frag%S_mat_prop_blocks(iblk)%val(ii, jj, ispin)
          end do
        end do
      end do
!$omp end parallel do
    else
!$omp parallel do schedule(static) &
!$omp& private(iblk, nrow, ncol, valid_row_count, valid_col_count, ii, jj, idx_ii, idx_jj, &
!$omp&         ig_i, ig_j, row_gid, col_gid, valid_row_ids, valid_col_ids)
      do iblk = 1, size(dg_frag%S_mat_blocks)
        if (.not. frag_selected(dg_frag%S_mat_blocks(iblk)%ifrag_row)) cycle
        if (.not. frag_selected(dg_frag%S_mat_blocks(iblk)%ifrag_col)) cycle
        nrow = dg_frag%n_basis(dg_frag%S_mat_blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(dg_frag%S_mat_blocks(iblk)%ifrag_col, ispin)
        valid_row_count = 0
        do ii = 1, nrow
          row_gid(ii) = dg_frag%index_basis(ii, dg_frag%S_mat_blocks(iblk)%ifrag_row, ispin)
          if (row_gid(ii) < 1 .or. row_gid(ii) > n_dim) cycle
          if (local_pos(row_gid(ii)) <= 0) cycle
          valid_row_count = valid_row_count + 1
          valid_row_ids(valid_row_count) = ii
        end do
        if (valid_row_count <= 0) cycle
        valid_col_count = 0
        do jj = 1, ncol
          col_gid(jj) = dg_frag%index_basis(jj, dg_frag%S_mat_blocks(iblk)%ifrag_col, ispin)
          if (col_gid(jj) < 1 .or. col_gid(jj) > n_dim) cycle
          if (local_pos(col_gid(jj)) <= 0) cycle
          valid_col_count = valid_col_count + 1
          valid_col_ids(valid_col_count) = jj
        end do
        if (valid_col_count <= 0) cycle
        do idx_jj = 1, valid_col_count
          jj = valid_col_ids(idx_jj)
          ig_j = col_gid(jj)
!$omp simd private(ii,ig_i)
          do idx_ii = 1, valid_row_count
            ii = valid_row_ids(idx_ii)
            ig_i = row_gid(ii)
            mat_loc(local_pos(ig_i), local_pos(ig_j)) = cmplx(dg_frag%S_mat_blocks(iblk)%val(ii, jj, ispin), 0.0d0, kind=8)
          end do
        end do
      end do
!$omp end parallel do
    end if

    call zgesv(n_loc, n_rhs, mat_loc, n_loc, ipiv, rhs_loc, n_loc, info)
    if (info == 0) then
      sol(row_ids(:), :) = rhs_loc(:, :)
    else
      sol(:, :) = rhs(:, :)
      call solve_overlap_operator_batch(dg_frag, ispin, rhs, sol, use_prop)
    end if

    deallocate(frag_selected, local_pos, row_ids, mat_loc, rhs_loc, ipiv)
  end subroutine solve_overlap_operator_batch_local

  subroutine copy_overlap_operator_to_dense(dg_frag, ispin, use_prop, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    logical, intent(in) :: use_prop
    complex(8), intent(inout) :: mat(:, :)

    integer :: n_frag, n_pw, n_tot, ipw

    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    mat(:, :) = (0.0d0, 0.0d0)

    if (n_pw > 0 .and. mixed_fp_coupling_active(dg_frag, ispin)) then
      if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks)) then
        call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%S_mat_prop_blocks, ispin, mat(1:n_frag, 1:n_frag))
      else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_blocks)) then
        call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%S_mat_blocks, ispin, mat(1:n_frag, 1:n_frag))
      else if (use_prop .and. allocated(dg_frag%S_mat_prop_c)) then
        mat(1:n_frag, 1:n_frag) = dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin)
      else if (use_prop .and. allocated(dg_frag%S_mat_prop)) then
        mat(1:n_frag, 1:n_frag) = cmplx(dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
      else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_c)) then
        mat(1:n_frag, 1:n_frag) = dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin)
      else if (allocated(dg_frag%S_mat)) then
        mat(1:n_frag, 1:n_frag) = cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
      end if
      mat(1:n_frag, n_frag+1:n_tot) = dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin)
      mat(n_frag+1:n_tot, 1:n_frag) = conjg(transpose(dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin)))
      do ipw = 1, n_pw
        mat(n_frag+ipw, n_frag+ipw) = (1.0d0, 0.0d0)
      end do
    else if (use_prop .and. allocated(dg_frag%S_mat_prop_blocks) .and. n_pw == 0) then
      call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%S_mat_prop_blocks, ispin, mat(1:n_frag, 1:n_frag))
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_blocks) .and. n_pw == 0) then
      call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%S_mat_blocks, ispin, mat(1:n_frag, 1:n_frag))
    else if (use_prop .and. allocated(dg_frag%S_mat_prop_c)) then
      mat(1:n_frag, 1:n_frag) = dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin)
    else if (use_prop .and. allocated(dg_frag%S_mat_prop)) then
      mat(1:n_frag, 1:n_frag) = cmplx(dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
    else if ((.not. use_prop) .and. allocated(dg_frag%S_mat_c)) then
      mat(1:n_frag, 1:n_frag) = dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin)
    else if (allocated(dg_frag%S_mat)) then
      mat(1:n_frag, 1:n_frag) = cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
    end if
    if (n_pw > 0) then
      do ipw = 1, n_pw
        mat(n_frag+ipw, n_frag+ipw) = (1.0d0, 0.0d0)
      end do
    end if
  end subroutine copy_overlap_operator_to_dense

  subroutine copy_momentum_blocks_to_complex_dense(dg_frag, ispin, scale_vec, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    real(8), intent(in) :: scale_vec(3)
    complex(8), intent(inout) :: mat(:, :)

    integer :: iblk, idir, ib, jb, row_idx, col_idx, n_frag
    integer :: nrow, ncol, idx_ib, idx_jb, valid_row_count, valid_col_count
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))

    mat(:, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%momentum_blocks)) return
    if (.not. allocated(dg_frag%index_basis)) return
    n_frag = dg_frag%n_mat_max
!$omp parallel do schedule(static) &
!$omp& private(iblk, nrow, ncol, valid_row_count, valid_col_count, ib, jb, idir, idx_ib, idx_jb, &
!$omp&         row_idx, col_idx, row_gid, col_gid, valid_row_ids, valid_col_ids)
    do iblk = 1, dg_frag%n_momentum_blocks
      if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
      nrow = dg_frag%n_basis(dg_frag%momentum_blocks(iblk)%ifrag_row, ispin)
      ncol = dg_frag%n_basis(dg_frag%momentum_blocks(iblk)%ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) cycle
      valid_row_count = 0
      do ib = 1, nrow
        row_gid(ib) = dg_frag%index_basis(ib, dg_frag%momentum_blocks(iblk)%ifrag_row, ispin)
        if (row_gid(ib) < 1 .or. row_gid(ib) > n_frag) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ib
      end do
      valid_col_count = 0
      do jb = 1, ncol
        col_gid(jb) = dg_frag%index_basis(jb, dg_frag%momentum_blocks(iblk)%ifrag_col, ispin)
        if (col_gid(jb) < 1 .or. col_gid(jb) > n_frag) cycle
        valid_col_count = valid_col_count + 1
        valid_col_ids(valid_col_count) = jb
      end do
      do idir = 1, 3
        if (abs(scale_vec(idir)) <= 0.0d0) cycle
        do idx_jb = 1, valid_col_count
          jb = valid_col_ids(idx_jb)
          col_idx = col_gid(jb)
!$omp simd private(ib,row_idx)
          do idx_ib = 1, valid_row_count
            ib = valid_row_ids(idx_ib)
            row_idx = row_gid(ib)
            mat(row_idx, col_idx) = mat(row_idx, col_idx) + &
              cmplx(scale_vec(idir) * dg_frag%momentum_blocks(iblk)%val(idir, ib, jb, ispin), 0.0d0, kind=8)
          end do
        end do
      end do
    end do
!$omp end parallel do
  end subroutine copy_momentum_blocks_to_complex_dense

  subroutine init_matrix_blocks_runtime(dg_frag, blocks, block_map)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)

    integer :: ifrag_row, ifrag_col, iblk, n_blocks, nrow_max, ncol_max, inei

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)
    call ensure_runtime_neighbor_pair_cache(dg_frag)
    call ensure_runtime_neighbor_list_cache(dg_frag)

    n_blocks = sum(dg_frag%runtime_neighbor_frag_count(1:max(1, dg_frag%n_frag)))
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_row = 1, dg_frag%n_frag
      do inei = 1, dg_frag%runtime_neighbor_frag_count(ifrag_row)
        ifrag_col = dg_frag%runtime_neighbor_frag_ids(inei, ifrag_row)
        if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
        iblk = iblk + 1
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        block_map(ifrag_row, ifrag_col) = iblk
        blocks(iblk)%ifrag_row = ifrag_row
        blocks(iblk)%ifrag_col = ifrag_col
        blocks(iblk)%nrow_max = nrow_max
        blocks(iblk)%ncol_max = ncol_max
        allocate(blocks(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_matrix_blocks_runtime

  subroutine reduce_matrix_blocks_runtime(dg_frag, blocks, label, icomm_reduce)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(inout) :: blocks(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce

    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: chunk_begin, chunk_count, offset_flat
    logical, parameter :: trace_block_reduce = .false.

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    if (trace_block_reduce .and. comm_is_root(dg_frag%id) .and. trim(label) /= "hmat-nl") then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size <= 0) return
    allocate(send_block(max_block_size), recv_block(max_block_size))

    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = blocks(iblk)%val(ii, jj, ispin)
            offset_flat = offset_flat + 1
          end do
        end do

        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
            offset_flat = offset_flat + 1
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (trace_block_reduce .and. comm_is_root(dg_frag%id) .and. trim(label) /= "hmat-nl") then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_matrix_blocks_runtime

  subroutine init_complex_matrix_blocks_runtime(dg_frag, blocks, block_map)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(complex_matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer :: ifrag_row, ifrag_col, iblk, n_blocks, nrow_max, ncol_max, inei

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)
    call ensure_runtime_neighbor_pair_cache(dg_frag)
    call ensure_runtime_neighbor_list_cache(dg_frag)

    n_blocks = sum(dg_frag%runtime_neighbor_frag_count(1:max(1, dg_frag%n_frag)))
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_row = 1, dg_frag%n_frag
      do inei = 1, dg_frag%runtime_neighbor_frag_count(ifrag_row)
        ifrag_col = dg_frag%runtime_neighbor_frag_ids(inei, ifrag_row)
        if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
        iblk = iblk + 1
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        block_map(ifrag_row, ifrag_col) = iblk
        blocks(iblk)%ifrag_row = ifrag_row
        blocks(iblk)%ifrag_col = ifrag_col
        blocks(iblk)%nrow_max = nrow_max
        blocks(iblk)%ncol_max = ncol_max
        allocate(blocks(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks(iblk)%val(:, :, :) = (0.0d0, 0.0d0)
      end do
    end do
  end subroutine init_complex_matrix_blocks_runtime

  subroutine init_complex_matrix_self_blocks_runtime(dg_frag, blocks, block_map)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(complex_matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer :: ifrag, iblk, nrow_max
    logical :: valid

    valid = allocated(blocks) .and. allocated(block_map)
    if (valid) valid = (size(blocks) == dg_frag%n_frag)
    if (valid) valid = (size(block_map, 1) == dg_frag%n_frag .and. size(block_map, 2) == dg_frag%n_frag)
    if (valid) then
      do ifrag = 1, dg_frag%n_frag
        if (block_map(ifrag, ifrag) /= ifrag) valid = .false.
      end do
    end if
    if (.not. valid) then
      if (allocated(blocks)) then
        do iblk = 1, size(blocks)
          if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
        end do
        deallocate(blocks)
      end if
      if (allocated(block_map)) deallocate(block_map)
      allocate(blocks(max(1, dg_frag%n_frag)))
      allocate(block_map(max(1, dg_frag%n_frag), max(1, dg_frag%n_frag)))
      block_map(:, :) = 0
      do ifrag = 1, dg_frag%n_frag
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag, 1:dg_frag%nspin)))
        block_map(ifrag, ifrag) = ifrag
        blocks(ifrag)%ifrag_row = ifrag
        blocks(ifrag)%ifrag_col = ifrag
        blocks(ifrag)%nrow_max = nrow_max
        blocks(ifrag)%ncol_max = nrow_max
        allocate(blocks(ifrag)%val(nrow_max, nrow_max, dg_frag%nspin))
      end do
    end if

    do iblk = 1, size(blocks)
      if (allocated(blocks(iblk)%val)) blocks(iblk)%val(:, :, :) = (0.0d0, 0.0d0)
    end do
  end subroutine init_complex_matrix_self_blocks_runtime

  subroutine reduce_complex_matrix_blocks_runtime(dg_frag, blocks, label, icomm_reduce)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(complex_matrix_block_info), intent(inout) :: blocks(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, ispin, ii, jj
    integer :: nrow, ncol, block_size, packed_size, max_block_size, total_active_size
    integer :: chunk_begin, chunk_count, offset_flat
    logical, parameter :: trace_block_reduce = .false.

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    if (trace_block_reduce .and. comm_is_root(dg_frag%id) .and. trim(label) /= "hmat-nl") then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size <= 0) return
    allocate(send_block(2 * max_block_size), recv_block(2 * max_block_size))

    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        packed_size = 2 * block_size

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = real(blocks(iblk)%val(ii, jj, ispin), kind=8)
            send_block(block_size + offset_flat) = aimag(blocks(iblk)%val(ii, jj, ispin))
            offset_flat = offset_flat + 1
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= packed_size)
          chunk_count = min(reduce_chunk_size, packed_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks(iblk)%val(ii, jj, ispin) = cmplx(recv_block(offset_flat), recv_block(block_size + offset_flat), kind=8)
            offset_flat = offset_flat + 1
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (trace_block_reduce .and. comm_is_root(dg_frag%id) .and. trim(label) /= "hmat-nl") then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_complex_matrix_blocks_runtime

  integer function find_matrix_block_runtime(block_map, ifrag_row, ifrag_col) result(iblk)
    implicit none
    integer, intent(in) :: block_map(:, :)
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (ifrag_row < 1 .or. ifrag_row > size(block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(block_map, 2)) return
    iblk = block_map(ifrag_row, ifrag_col)
  end function find_matrix_block_runtime

  subroutine pack_owned_coef(dg_frag, ispin, row_ids, packed)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: packed(:, :)

    integer :: irow, global_row, local_idx

    packed(:, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    if (.not. allocated(dg_frag%coef)) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return

!$omp parallel do private(irow, global_row, local_idx) schedule(static)
    do irow = 1, min(size(row_ids), size(packed, 1))
      global_row = row_ids(irow)
      if (global_row < 1 .or. global_row > size(dg_frag%coef_owner, 1)) cycle
      if (dg_frag%coef_owner(global_row, ispin) /= dg_frag%id) cycle
      local_idx = dg_frag%coef_global_to_local(global_row, ispin)
      if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
      packed(irow, 1:size(packed, 2)) = dg_frag%coef(local_idx, 1:size(packed, 2), ispin)
    end do
!$omp end parallel do
  end subroutine pack_owned_coef

  subroutine fetch_remote_coef_rows(dg_frag, ispin, row_ids, fetched, col_start, col_end, collective_counts)
    use mpi, only: MPI_Isend, MPI_Irecv, MPI_Waitall, MPI_Alltoallv, MPI_DOUBLE_COMPLEX, MPI_INTEGER, &
                   MPI_SUCCESS, MPI_STATUSES_IGNORE
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: fetched(:, :)
    integer, intent(in), optional :: col_start, col_end
    logical, intent(in), optional :: collective_counts

    integer :: irow, global_row, owner_rank, local_idx, row_frag
    integer :: nrows_req, nstate_req, owner, peer, k, ist, pos, slot_idx, ipeer
    integer :: cstart, cend, total_req_send, total_req_recv, total_data_send, total_data_recv
    integer :: state_owner_offset, state_owner_end_offset, local_cstart
    integer :: ierr
    integer, allocatable, save :: req_send_counts(:), req_recv_counts(:)
    integer, allocatable, save :: req_send_displs(:), req_recv_displs(:)
    integer, allocatable, save :: data_send_counts(:), data_recv_counts(:)
    integer, allocatable, save :: data_send_displs(:), data_recv_displs(:)
    integer, allocatable, save :: req_send_fill(:), req_send_slot(:), req_send_rows(:), req_recv_rows(:)
    integer, allocatable, save :: peer_ranks(:), requests(:)
    logical, allocatable, save :: peer_needed(:)
    complex(8), allocatable, save :: data_send_buf(:), data_recv_buf(:)
    integer :: npeer, nreq
    logical :: use_cached_plan
    logical :: use_alltoallv_data
    logical :: use_collective_counts
    logical, save :: fetch_options_initialized = .false.
    logical, save :: enable_fetch_alltoallv = .false.
    character(16) :: env_fetch_alltoallv
    integer :: env_status
    integer, save :: fetch_work_isize = -1
    integer, save :: fetch_work_peer_cap = 0
    integer, save :: fetch_work_request_cap = 0
    integer, save :: fetch_work_req_send_cap = 0
    integer, save :: fetch_work_req_recv_cap = 0
    integer, save :: fetch_work_data_send_cap = 0
    integer, save :: fetch_work_data_recv_cap = 0
    logical, save :: fetch_plan_valid = .false.
    integer, save :: fetch_plan_ispin = -1
    integer, save :: fetch_plan_isize = -1
    integer, save :: fetch_plan_nrows = -1
    integer, allocatable, save :: fetch_plan_row_ids(:)
    integer, allocatable, save :: fetch_plan_peer_ranks(:)
    integer, allocatable, save :: fetch_plan_req_send_counts(:)
    integer, allocatable, save :: fetch_plan_req_recv_counts(:)
    integer, allocatable, save :: fetch_plan_req_send_displs(:)
    integer, allocatable, save :: fetch_plan_req_recv_displs(:)
    integer, allocatable, save :: fetch_plan_req_send_slot(:)
    integer, allocatable, save :: fetch_plan_req_send_rows(:)
    integer, allocatable, save :: fetch_plan_req_recv_rows(:)
    integer, parameter :: tag_count = 8411
    integer, parameter :: tag_rows  = 8412
    integer, parameter :: tag_data  = 8413
    logical, parameter :: enable_coef_gather_trace = .false.

    if (.not. fetch_options_initialized) then
      env_fetch_alltoallv = ''
      call get_environment_variable('SALMON_DG_FETCH_ALLTOALLV', env_fetch_alltoallv, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_fetch_alltoallv)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enable_fetch_alltoallv = .true.
        end select
      end if
      fetch_options_initialized = .true.
    end if

    fetched(:, :) = (0.0d0, 0.0d0)
    use_collective_counts = .false.
    if (present(collective_counts)) use_collective_counts = collective_counts
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    if (.not. allocated(dg_frag%coef)) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (dg_frag%isize <= 0) return

    nrows_req = min(size(row_ids), size(fetched, 1))
    nstate_req = size(fetched, 2)
    if (nrows_req <= 0 .or. nstate_req <= 0) return
    cstart = 1
    cend = nstate_req
    if (present(col_start)) cstart = col_start
    if (present(col_end)) cend = col_end
    if (dg_frag%coef_state_block_mode) then
      if (cstart < 1 .or. cend < cstart .or. cend > dg_frag%nstate_tot) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] invalid global state range in state-block fetch: rank=", &
          dg_frag%id, " cstart=", cstart, " cend=", cend, " nstate_tot=", dg_frag%nstate_tot
        stop "DG-Fragment RT: invalid state-block coefficient range"
      end if
      state_owner_offset = dg_state_owner_offset(dg_frag%nstate_tot, max(1, dg_frag%isize_frag), cstart)
      state_owner_end_offset = dg_state_owner_offset(dg_frag%nstate_tot, max(1, dg_frag%isize_frag), cend)
      if (state_owner_offset /= state_owner_end_offset) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] state-block fetch crosses state owners: rank=", &
          dg_frag%id, " cstart=", cstart, " cend=", cend, " owner_s=", state_owner_offset, &
          " owner_e=", state_owner_end_offset
        stop "DG-Fragment RT: state-block coefficient fetch crosses owners"
      end if
    else if (cstart < 1 .or. cend < cstart .or. cend > size(dg_frag%coef, 2)) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] invalid state range in fetch_remote_coef_rows: rank=", dg_frag%id, &
        " cstart=", cstart, " cend=", cend, " coef_cols=", size(dg_frag%coef, 2)
      stop "DG-Fragment RT: invalid state range in fetch_remote_coef_rows"
    end if
    if (cend - cstart + 1 /= nstate_req) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] fetch_remote_coef_rows state-range mismatch: rank=", dg_frag%id, &
        " nstate_req=", nstate_req, " cstart=", cstart, " cend=", cend
      stop "DG-Fragment RT: state range/shape mismatch in fetch_remote_coef_rows"
    end if

    call ensure_runtime_neighbor_pair_cache(dg_frag)
    call ensure_coef_row_fragment_cache(dg_frag)
    call ensure_coef_exchange_peer_cache(dg_frag)
    call ensure_coef_allowed_request_frag_cache(dg_frag)
    if (fetch_work_isize /= dg_frag%isize) then
      if (allocated(req_send_counts)) deallocate(req_send_counts, req_recv_counts)
      if (allocated(req_send_displs)) deallocate(req_send_displs, req_recv_displs)
      if (allocated(data_send_counts)) deallocate(data_send_counts, data_recv_counts)
      if (allocated(data_send_displs)) deallocate(data_send_displs, data_recv_displs)
      if (allocated(req_send_fill)) deallocate(req_send_fill)
      if (allocated(peer_needed)) deallocate(peer_needed)
      if (allocated(peer_ranks)) deallocate(peer_ranks)
      if (allocated(requests)) deallocate(requests)
      if (allocated(req_send_slot)) deallocate(req_send_slot, req_send_rows)
      if (allocated(req_recv_rows)) deallocate(req_recv_rows)
      if (allocated(data_send_buf)) deallocate(data_send_buf)
      if (allocated(data_recv_buf)) deallocate(data_recv_buf)
      fetch_work_isize = dg_frag%isize
      fetch_work_peer_cap = 0
      fetch_work_request_cap = 0
      fetch_work_req_send_cap = 0
      fetch_work_req_recv_cap = 0
      fetch_work_data_send_cap = 0
      fetch_work_data_recv_cap = 0
    end if
    if (.not. allocated(req_send_counts)) allocate(req_send_counts(0:dg_frag%isize-1), req_recv_counts(0:dg_frag%isize-1))
    if (.not. allocated(req_send_displs)) allocate(req_send_displs(0:dg_frag%isize-1), req_recv_displs(0:dg_frag%isize-1))
    if (.not. allocated(data_send_counts)) allocate(data_send_counts(0:dg_frag%isize-1), data_recv_counts(0:dg_frag%isize-1))
    if (.not. allocated(data_send_displs)) allocate(data_send_displs(0:dg_frag%isize-1), data_recv_displs(0:dg_frag%isize-1))
    if (.not. allocated(req_send_fill)) allocate(req_send_fill(0:dg_frag%isize-1))
    if (.not. allocated(peer_needed)) allocate(peer_needed(0:dg_frag%isize-1))
    req_send_counts(:) = 0
    req_recv_counts(:) = 0
    req_send_displs(:) = 0
    req_recv_displs(:) = 0
    data_send_counts(:) = 0
    data_recv_counts(:) = 0
    data_send_displs(:) = 0
    data_recv_displs(:) = 0
    req_send_fill(:) = 0
    peer_needed(:) = .false.
    if ((.not. dg_frag%coef_state_block_mode) .and. &
        allocated(dg_frag%coef_exchange_peer_count) .and. allocated(dg_frag%coef_exchange_peer_ranks)) then
      if (ispin <= size(dg_frag%coef_exchange_peer_count)) then
        do ipeer = 1, dg_frag%coef_exchange_peer_count(ispin)
          peer = dg_frag%coef_exchange_peer_ranks(ipeer, ispin)
          if (peer >= 0 .and. peer < dg_frag%isize .and. peer /= dg_frag%id) peer_needed(peer) = .true.
        end do
      end if
    end if

    if (enable_coef_gather_trace .and. dg_frag%id == 0 .and. ispin == 1) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        coef gather trace: rank=", dg_frag%id, &
        " ispin=", ispin, " nrows=", nrows_req, " nstate=", nstate_req
      flush(6)
    end if

    do irow = 1, nrows_req
      global_row = row_ids(irow)
      if (global_row < 1 .or. global_row > size(dg_frag%coef_owner, 1)) cycle
      row_frag = 0
      if (allocated(dg_frag%coef_row_fragment)) then
        if (global_row <= size(dg_frag%coef_row_fragment, 1) .and. ispin <= size(dg_frag%coef_row_fragment, 2)) then
          row_frag = dg_frag%coef_row_fragment(global_row, ispin)
        end if
      end if
      if (row_frag < 1 .or. row_frag > dg_frag%n_frag) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] coefficient row has no fragment owner: rank=", dg_frag%id, &
          " ispin=", ispin, " row=", global_row
        stop "DG-Fragment RT: coefficient row fragment map is invalid"
      end if
      if (.not. dg_frag%coef_allowed_request_frag(row_frag)) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] nonlocal coefficient request is outside runtime-neighbor stencil: rank=", &
          dg_frag%id, " ispin=", ispin, " row=", global_row, " row_frag=", row_frag
        stop "DG-Fragment RT: non-neighbor coefficient request"
      end if
      if (dg_frag%coef_state_block_mode) then
        owner_rank = dg_frag%id_array(row_frag) + state_owner_offset
      else
        owner_rank = dg_frag%coef_owner(global_row, ispin)
      end if
      if (owner_rank < 0) cycle
      if (owner_rank >= dg_frag%isize) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] invalid coef owner rank: rank=", dg_frag%id, &
          " ispin=", ispin, " row=", global_row, " owner_rank=", owner_rank, " isize=", dg_frag%isize
        stop "DG-Fragment RT: invalid coef owner rank in fetch_remote_coef_rows"
      end if
      if (owner_rank == dg_frag%id) then
        local_idx = dg_frag%coef_global_to_local(global_row, ispin)
        if (local_idx >= 1 .and. local_idx <= size(dg_frag%coef, 1)) then
          if (dg_frag%coef_state_block_mode) then
            local_cstart = cstart - dg_frag%coef_state_start + 1
            if (local_cstart >= 1 .and. local_cstart + nstate_req - 1 <= size(dg_frag%coef, 2)) then
              fetched(irow, 1:nstate_req) = dg_frag%coef(local_idx, local_cstart:local_cstart+nstate_req-1, ispin)
            end if
          else
            fetched(irow, 1:nstate_req) = dg_frag%coef(local_idx, cstart:cend, ispin)
          end if
        end if
        cycle
      end if
      req_send_counts(owner_rank) = req_send_counts(owner_rank) + 1
      peer_needed(owner_rank) = .true.
    end do
    if (dg_frag%coef_state_block_mode .and. .not. use_collective_counts) then
      do row_frag = 1, dg_frag%n_frag
        if (allocated(dg_frag%coef_allowed_request_frag)) then
          if (.not. dg_frag%coef_allowed_request_frag(row_frag)) cycle
        else if (.not. is_runtime_neighbor_pair(dg_frag, dg_frag%ifrag_group, row_frag)) then
          cycle
        end if
        owner_rank = dg_frag%id_array(row_frag) + state_owner_offset
        if (owner_rank >= 0 .and. owner_rank < dg_frag%isize .and. owner_rank /= dg_frag%id) then
          peer_needed(owner_rank) = .true.
        end if
      end do
    end if

    use_cached_plan = (.not. dg_frag%coef_state_block_mode) .and. fetch_plan_valid .and. fetch_plan_ispin == ispin .and. &
                      fetch_plan_isize == dg_frag%isize .and. fetch_plan_nrows == nrows_req .and. &
                      allocated(fetch_plan_row_ids)
    if (use_cached_plan) then
      if (size(fetch_plan_row_ids) /= nrows_req) then
        use_cached_plan = .false.
      else if (any(fetch_plan_row_ids(1:nrows_req) /= row_ids(1:nrows_req))) then
        use_cached_plan = .false.
      end if
    end if

    if (use_cached_plan) then
      npeer = size(fetch_plan_peer_ranks)
      if (fetch_work_peer_cap < max(1, npeer)) then
        if (allocated(peer_ranks)) deallocate(peer_ranks)
        fetch_work_peer_cap = max(1, npeer)
        allocate(peer_ranks(fetch_work_peer_cap))
      end if
      if (npeer > 0) peer_ranks(1:npeer) = fetch_plan_peer_ranks(1:npeer)
      req_send_counts(:) = fetch_plan_req_send_counts(:)
      req_recv_counts(:) = fetch_plan_req_recv_counts(:)
      req_send_displs(:) = fetch_plan_req_send_displs(:)
      req_recv_displs(:) = fetch_plan_req_recv_displs(:)
      total_req_send = 0
      total_req_recv = 0
      if (dg_frag%isize > 0) then
        total_req_send = req_send_displs(dg_frag%isize-1) + req_send_counts(dg_frag%isize-1)
        total_req_recv = req_recv_displs(dg_frag%isize-1) + req_recv_counts(dg_frag%isize-1)
      end if
      if (fetch_work_req_send_cap < max(1, total_req_send)) then
        if (allocated(req_send_slot)) deallocate(req_send_slot, req_send_rows)
        fetch_work_req_send_cap = max(1, total_req_send)
        allocate(req_send_rows(fetch_work_req_send_cap), req_send_slot(fetch_work_req_send_cap))
      end if
      if (fetch_work_req_recv_cap < max(1, total_req_recv)) then
        if (allocated(req_recv_rows)) deallocate(req_recv_rows)
        fetch_work_req_recv_cap = max(1, total_req_recv)
        allocate(req_recv_rows(fetch_work_req_recv_cap))
      end if
      if (total_req_send > 0) then
        req_send_rows(1:total_req_send) = fetch_plan_req_send_rows(1:total_req_send)
        req_send_slot(1:total_req_send) = fetch_plan_req_send_slot(1:total_req_send)
      end if
      if (total_req_recv > 0) req_recv_rows(1:total_req_recv) = fetch_plan_req_recv_rows(1:total_req_recv)
    else
      if (dg_frag%coef_state_block_mode .and. use_collective_counts) then
        do peer = 0, dg_frag%isize-1
          peer_needed(peer) = (peer /= dg_frag%id)
        end do
      end if

      npeer = count(peer_needed)
      if (npeer <= 0) then
        return
      end if
      if (fetch_work_peer_cap < max(1, npeer)) then
        if (allocated(peer_ranks)) deallocate(peer_ranks)
        fetch_work_peer_cap = max(1, npeer)
        allocate(peer_ranks(fetch_work_peer_cap))
      end if
      npeer = 0
      do peer = 0, dg_frag%isize-1
        if (.not. peer_needed(peer)) cycle
        npeer = npeer + 1
        peer_ranks(npeer) = peer
      end do

      if (fetch_work_request_cap < max(1, 2*npeer)) then
        if (allocated(requests)) deallocate(requests)
        fetch_work_request_cap = max(1, 2*npeer)
        allocate(requests(fetch_work_request_cap))
      end if
      nreq = 0
      do ipeer = 1, npeer
        peer = peer_ranks(ipeer)
        nreq = nreq + 1
        call MPI_Irecv(req_recv_counts(peer), 1, MPI_INTEGER, peer, tag_count, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: request-count recv failed in fetch_remote_coef_rows"
      end do
      do ipeer = 1, npeer
        peer = peer_ranks(ipeer)
        nreq = nreq + 1
        call MPI_Isend(req_send_counts(peer), 1, MPI_INTEGER, peer, tag_count, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: request-count send failed in fetch_remote_coef_rows"
      end do
      if (nreq > 0) then
        call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: request-count wait failed in fetch_remote_coef_rows"
      end if

      req_send_displs(0) = 0
      req_recv_displs(0) = 0
      data_send_displs(0) = 0
      data_recv_displs(0) = 0
      do owner = 1, dg_frag%isize-1
        req_send_displs(owner) = req_send_displs(owner-1) + req_send_counts(owner-1)
        req_recv_displs(owner) = req_recv_displs(owner-1) + req_recv_counts(owner-1)
      end do
      total_req_send = req_send_displs(dg_frag%isize-1) + req_send_counts(dg_frag%isize-1)
      total_req_recv = req_recv_displs(dg_frag%isize-1) + req_recv_counts(dg_frag%isize-1)

      if (fetch_work_req_send_cap < max(1, total_req_send)) then
        if (allocated(req_send_slot)) deallocate(req_send_slot, req_send_rows)
        fetch_work_req_send_cap = max(1, total_req_send)
        allocate(req_send_rows(fetch_work_req_send_cap), req_send_slot(fetch_work_req_send_cap))
      end if
      if (fetch_work_req_recv_cap < max(1, total_req_recv)) then
        if (allocated(req_recv_rows)) deallocate(req_recv_rows)
        fetch_work_req_recv_cap = max(1, total_req_recv)
        allocate(req_recv_rows(fetch_work_req_recv_cap))
      end if
      req_send_rows(1:max(1, total_req_send)) = 0
      req_send_slot(1:max(1, total_req_send)) = 0
      req_recv_rows(1:max(1, total_req_recv)) = 0

      do irow = 1, nrows_req
        global_row = row_ids(irow)
        if (global_row < 1 .or. global_row > size(dg_frag%coef_owner, 1)) cycle
        if (dg_frag%coef_state_block_mode) then
          row_frag = dg_frag%coef_row_fragment(global_row, ispin)
          owner_rank = dg_frag%id_array(row_frag) + state_owner_offset
        else
          owner_rank = dg_frag%coef_owner(global_row, ispin)
        end if
        if (owner_rank < 0 .or. owner_rank >= dg_frag%isize) cycle
        if (owner_rank == dg_frag%id) cycle
        pos = req_send_displs(owner_rank) + req_send_fill(owner_rank) + 1
        req_send_fill(owner_rank) = req_send_fill(owner_rank) + 1
        req_send_slot(pos) = irow
        req_send_rows(pos) = global_row
      end do

      nreq = 0
      do ipeer = 1, npeer
        peer = peer_ranks(ipeer)
        if (req_recv_counts(peer) <= 0) cycle
        nreq = nreq + 1
        call MPI_Irecv(req_recv_rows(req_recv_displs(peer)+1), req_recv_counts(peer), MPI_INTEGER, &
                       peer, tag_rows, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: request-row recv failed in fetch_remote_coef_rows"
      end do
      do ipeer = 1, npeer
        peer = peer_ranks(ipeer)
        if (req_send_counts(peer) <= 0) cycle
        nreq = nreq + 1
        call MPI_Isend(req_send_rows(req_send_displs(peer)+1), req_send_counts(peer), MPI_INTEGER, &
                       peer, tag_rows, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: request-row send failed in fetch_remote_coef_rows"
      end do
      if (nreq > 0) then
        call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: request-row wait failed in fetch_remote_coef_rows"
      end if

      if (allocated(fetch_plan_row_ids)) deallocate(fetch_plan_row_ids)
      if (allocated(fetch_plan_peer_ranks)) deallocate(fetch_plan_peer_ranks)
      if (allocated(fetch_plan_req_send_counts)) deallocate(fetch_plan_req_send_counts)
      if (allocated(fetch_plan_req_recv_counts)) deallocate(fetch_plan_req_recv_counts)
      if (allocated(fetch_plan_req_send_displs)) deallocate(fetch_plan_req_send_displs)
      if (allocated(fetch_plan_req_recv_displs)) deallocate(fetch_plan_req_recv_displs)
      if (allocated(fetch_plan_req_send_slot)) deallocate(fetch_plan_req_send_slot)
      if (allocated(fetch_plan_req_send_rows)) deallocate(fetch_plan_req_send_rows)
      if (allocated(fetch_plan_req_recv_rows)) deallocate(fetch_plan_req_recv_rows)
      allocate(fetch_plan_row_ids(nrows_req))
      allocate(fetch_plan_peer_ranks(npeer))
      allocate(fetch_plan_req_send_counts(0:dg_frag%isize-1), fetch_plan_req_recv_counts(0:dg_frag%isize-1))
      allocate(fetch_plan_req_send_displs(0:dg_frag%isize-1), fetch_plan_req_recv_displs(0:dg_frag%isize-1))
      allocate(fetch_plan_req_send_slot(max(1, total_req_send)), fetch_plan_req_send_rows(max(1, total_req_send)))
      allocate(fetch_plan_req_recv_rows(max(1, total_req_recv)))
      fetch_plan_row_ids(:) = row_ids(1:nrows_req)
      if (npeer > 0) fetch_plan_peer_ranks(:) = peer_ranks(:)
      fetch_plan_req_send_counts(:) = req_send_counts(:)
      fetch_plan_req_recv_counts(:) = req_recv_counts(:)
      fetch_plan_req_send_displs(:) = req_send_displs(:)
      fetch_plan_req_recv_displs(:) = req_recv_displs(:)
      if (total_req_send > 0) then
        fetch_plan_req_send_slot(1:total_req_send) = req_send_slot(1:total_req_send)
        fetch_plan_req_send_rows(1:total_req_send) = req_send_rows(1:total_req_send)
      end if
      if (total_req_recv > 0) fetch_plan_req_recv_rows(1:total_req_recv) = req_recv_rows(1:total_req_recv)
      fetch_plan_ispin = ispin
      fetch_plan_isize = dg_frag%isize
      fetch_plan_nrows = nrows_req
      fetch_plan_valid = .true.
    end if

    do peer = 0, dg_frag%isize-1
      data_send_counts(peer) = req_recv_counts(peer) * nstate_req
      data_recv_counts(peer) = req_send_counts(peer) * nstate_req
    end do
    do peer = 1, dg_frag%isize-1
      data_send_displs(peer) = data_send_displs(peer-1) + data_send_counts(peer-1)
      data_recv_displs(peer) = data_recv_displs(peer-1) + data_recv_counts(peer-1)
    end do
    total_data_send = data_send_displs(dg_frag%isize-1) + data_send_counts(dg_frag%isize-1)
    total_data_recv = data_recv_displs(dg_frag%isize-1) + data_recv_counts(dg_frag%isize-1)

    if (fetch_work_data_send_cap < max(1, total_data_send)) then
      if (allocated(data_send_buf)) deallocate(data_send_buf)
      fetch_work_data_send_cap = max(1, total_data_send)
      allocate(data_send_buf(fetch_work_data_send_cap))
    end if
    if (fetch_work_data_recv_cap < max(1, total_data_recv)) then
      if (allocated(data_recv_buf)) deallocate(data_recv_buf)
      fetch_work_data_recv_cap = max(1, total_data_recv)
      allocate(data_recv_buf(fetch_work_data_recv_cap))
    end if
    data_send_buf(1:max(1, total_data_send)) = (0.0d0, 0.0d0)
    data_recv_buf(1:max(1, total_data_recv)) = (0.0d0, 0.0d0)
    if (fetch_work_request_cap < max(1, 2*npeer)) then
      if (allocated(requests)) deallocate(requests)
      fetch_work_request_cap = max(1, 2*npeer)
      allocate(requests(fetch_work_request_cap))
    end if

!$omp parallel do private(peer, k, irow, global_row, local_idx, pos, ist) schedule(static)
    do peer = 0, dg_frag%isize-1
      do k = 1, req_recv_counts(peer)
        irow = req_recv_displs(peer) + k
        global_row = req_recv_rows(irow)
        if (global_row < 1 .or. global_row > size(dg_frag%coef_global_to_local, 1)) cycle
        local_idx = dg_frag%coef_global_to_local(global_row, ispin)
        if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
        pos = data_send_displs(peer) + (k - 1) * nstate_req
        do ist = 1, nstate_req
          if (dg_frag%coef_state_block_mode) then
            local_cstart = cstart - dg_frag%coef_state_start + 1
            if (local_cstart + ist - 1 >= 1 .and. local_cstart + ist - 1 <= size(dg_frag%coef, 2)) then
              data_send_buf(pos + ist) = dg_frag%coef(local_idx, local_cstart + ist - 1, ispin)
            end if
          else
            data_send_buf(pos + ist) = dg_frag%coef(local_idx, cstart + ist - 1, ispin)
          end if
        end do
      end do
    end do
!$omp end parallel do

    use_alltoallv_data = enable_fetch_alltoallv .and. use_cached_plan
    if (use_alltoallv_data) then
      call MPI_Alltoallv(data_send_buf, data_send_counts, data_send_displs, MPI_DOUBLE_COMPLEX, &
                         data_recv_buf, data_recv_counts, data_recv_displs, MPI_DOUBLE_COMPLEX, &
                         dg_frag%icomm, ierr)
      if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: coefficient alltoallv failed in fetch_remote_coef_rows"
    else
      nreq = 0
      do ipeer = 1, npeer
        peer = peer_ranks(ipeer)
        if (data_recv_counts(peer) <= 0) cycle
        nreq = nreq + 1
        call MPI_Irecv(data_recv_buf(data_recv_displs(peer)+1), data_recv_counts(peer), MPI_DOUBLE_COMPLEX, &
                       peer, tag_data, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: coefficient recv failed in fetch_remote_coef_rows"
      end do
      do ipeer = 1, npeer
        peer = peer_ranks(ipeer)
        if (data_send_counts(peer) <= 0) cycle
        nreq = nreq + 1
        call MPI_Isend(data_send_buf(data_send_displs(peer)+1), data_send_counts(peer), MPI_DOUBLE_COMPLEX, &
                       peer, tag_data, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: coefficient send failed in fetch_remote_coef_rows"
      end do
      if (nreq > 0) then
        call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
        if (ierr /= MPI_SUCCESS) stop "DG-Fragment RT: coefficient wait failed in fetch_remote_coef_rows"
      end if
    end if

!$omp parallel do private(owner, k, irow, slot_idx, pos, ist) schedule(static)
    do owner = 0, dg_frag%isize-1
      do k = 1, req_send_counts(owner)
        irow = req_send_displs(owner) + k
        slot_idx = req_send_slot(irow)
        if (slot_idx < 1 .or. slot_idx > nrows_req) cycle
        pos = data_recv_displs(owner) + (k - 1) * nstate_req
        do ist = 1, nstate_req
          fetched(slot_idx, ist) = data_recv_buf(pos + ist)
        end do
      end do
    end do
!$omp end parallel do

  end subroutine fetch_remote_coef_rows

  subroutine pack_owned_coef_pw(dg_frag, row_ids, packed)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: packed(:, :, :)

    integer :: irow, pw_row

    packed(:, :, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_pw_owner)) return
    if (.not. allocated(dg_frag%coef_pw)) return

!$omp parallel do private(irow, pw_row) schedule(static)
    do irow = 1, min(size(row_ids), size(packed, 1))
      pw_row = row_ids(irow)
      if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) cycle
      if (dg_frag%coef_pw_owner(pw_row) /= dg_frag%id) cycle
      packed(irow, 1:size(packed, 2), 1:size(packed, 3)) = &
        dg_frag%coef_pw(pw_row, 1:size(packed, 2), 1:size(packed, 3))
    end do
!$omp end parallel do
  end subroutine pack_owned_coef_pw

  subroutine fetch_remote_coef_pw_rows(dg_frag, row_ids, fetched, col_start, col_end, ispin_req)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: fetched(:, :, :)
    integer, intent(in), optional :: col_start, col_end
    integer, intent(in), optional :: ispin_req

    integer :: irow, pw_row, owner_rank
    integer :: nrows_req, nstate_req, nspin_req, owner, k, cnt, max_owner_rows
    integer :: cstart, cend, owner_chunk_rows, chunk0, chunk1, chunk_cnt, slot_idx
    integer :: spin_sel
    integer(8) :: target_bytes, bytes_per_pw_row
    integer, allocatable :: owner_counts(:), owner_offsets(:), owner_fill(:)
    integer, allocatable :: owner_slot(:), owner_row(:)
    complex(8), allocatable :: packed_rows(:, :, :)
    integer, parameter :: pw_bcast_target_mb = 64

    fetched(:, :, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_pw_owner)) return
    if (.not. allocated(dg_frag%coef_pw)) return
    if (dg_frag%isize <= 0) return

    nrows_req = min(size(row_ids), size(fetched, 1))
    nstate_req = size(fetched, 2)
    nspin_req = size(fetched, 3)
    if (nrows_req <= 0 .or. nstate_req <= 0 .or. nspin_req <= 0) return
    cstart = 1
    cend = nstate_req
    if (present(col_start)) cstart = col_start
    if (present(col_end)) cend = col_end
    if (cstart < 1 .or. cend < cstart .or. cend > size(dg_frag%coef_pw, 2)) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] invalid state range in fetch_remote_coef_pw_rows: rank=", dg_frag%id, &
        " cstart=", cstart, " cend=", cend, " coef_pw_cols=", size(dg_frag%coef_pw, 2)
      stop "DG-Fragment RT: invalid state range in fetch_remote_coef_pw_rows"
    end if
    if (cend - cstart + 1 /= nstate_req) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] fetch_remote_coef_pw_rows state-range mismatch: rank=", dg_frag%id, &
        " nstate_req=", nstate_req, " cstart=", cstart, " cend=", cend
      stop "DG-Fragment RT: state range/shape mismatch in fetch_remote_coef_pw_rows"
    end if
    spin_sel = 1
    if (present(ispin_req)) then
      spin_sel = ispin_req
      if (spin_sel < 1 .or. spin_sel > dg_frag%nspin) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] invalid ispin_req in fetch_remote_coef_pw_rows: rank=", dg_frag%id, &
          " ispin_req=", spin_sel, " nspin=", dg_frag%nspin
        stop "DG-Fragment RT: invalid ispin_req in fetch_remote_coef_pw_rows"
      end if
    end if

    allocate(owner_counts(0:dg_frag%isize-1), owner_offsets(0:dg_frag%isize-1), owner_fill(0:dg_frag%isize-1))
    allocate(owner_slot(nrows_req), owner_row(nrows_req))
    owner_counts(:) = 0
    owner_offsets(:) = 0
    owner_fill(:) = 0
    owner_slot(:) = 0
    owner_row(:) = 0

    do irow = 1, nrows_req
      pw_row = row_ids(irow)
      if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) cycle
      owner_rank = dg_frag%coef_pw_owner(pw_row)
      if (owner_rank < 0) cycle
      if (owner_rank >= dg_frag%isize) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] invalid coef_pw owner rank: rank=", dg_frag%id, &
          " row=", pw_row, " owner_rank=", owner_rank, " isize=", dg_frag%isize
        stop "DG-Fragment RT: invalid coef_pw owner rank in fetch_remote_coef_pw_rows"
      end if
      owner_counts(owner_rank) = owner_counts(owner_rank) + 1
    end do

    max_owner_rows = 0
    do owner = 0, dg_frag%isize-1
      if (owner_counts(owner) > max_owner_rows) max_owner_rows = owner_counts(owner)
    end do
    if (max_owner_rows <= 0) then
      deallocate(owner_counts, owner_offsets, owner_fill, owner_slot, owner_row)
      return
    end if

    owner_offsets(0) = 1
    do owner = 1, dg_frag%isize-1
      owner_offsets(owner) = owner_offsets(owner-1) + owner_counts(owner-1)
    end do

    do irow = 1, nrows_req
      pw_row = row_ids(irow)
      if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) cycle
      owner_rank = dg_frag%coef_pw_owner(pw_row)
      if (owner_rank < 0 .or. owner_rank >= dg_frag%isize) cycle
      k = owner_offsets(owner_rank) + owner_fill(owner_rank)
      owner_fill(owner_rank) = owner_fill(owner_rank) + 1
      owner_slot(k) = irow
      owner_row(k) = pw_row
    end do

    target_bytes = int(pw_bcast_target_mb, kind=8) * 1024_8 * 1024_8
    bytes_per_pw_row = 16_8 * int(max(1, nstate_req), kind=8) * int(max(1, nspin_req), kind=8)
    owner_chunk_rows = max(1, min(max_owner_rows, int(max(1_8, target_bytes / max(1_8, bytes_per_pw_row)))))
    do owner = 0, dg_frag%isize-1
      cnt = owner_counts(owner)
      if (cnt <= 0) cycle
      do chunk0 = 1, cnt, owner_chunk_rows
        chunk1 = min(cnt, chunk0 + owner_chunk_rows - 1)
        chunk_cnt = chunk1 - chunk0 + 1
        allocate(packed_rows(chunk_cnt, nstate_req, nspin_req))
        packed_rows(:, :, :) = (0.0d0, 0.0d0)
        if (owner == dg_frag%id) then
          do k = 1, chunk_cnt
            irow = owner_offsets(owner) + (chunk0 + k - 2)
            pw_row = owner_row(irow)
            if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) cycle
            if (present(ispin_req)) then
              packed_rows(k, 1:nstate_req, 1) = dg_frag%coef_pw(pw_row, cstart:cend, spin_sel)
            else
              packed_rows(k, 1:nstate_req, 1:nspin_req) = dg_frag%coef_pw(pw_row, cstart:cend, 1:nspin_req)
            end if
          end do
        end if
        call comm_bcast(packed_rows, dg_frag%icomm, owner)
        do k = 1, chunk_cnt
          irow = owner_offsets(owner) + (chunk0 + k - 2)
          slot_idx = owner_slot(irow)
          if (slot_idx < 1 .or. slot_idx > nrows_req) cycle
          fetched(slot_idx, 1:nstate_req, 1:nspin_req) = packed_rows(k, 1:nstate_req, 1:nspin_req)
        end do
        deallocate(packed_rows)
      end do
    end do

    deallocate(owner_counts, owner_offsets, owner_fill, owner_slot, owner_row)
  end subroutine fetch_remote_coef_pw_rows

  subroutine refresh_pw_coef_cache(dg_frag, nstate_use)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in), optional :: nstate_use

    integer :: i, n_pw, nstate_cache
    integer :: env_status, env_len
    integer, allocatable :: pw_row_ids(:)
    real(8) :: cache_diff, cache_norm
    logical :: use_sum_cache, compare_cache, trace_cache
    character(32) :: env_cache_mode, env_cache_trace
    complex(8), allocatable :: coef_pw_reduce(:,:,:)

    if (.not. dg_frag%use_plane_wave_basis) then
      if (allocated(dg_frag%coef_pw_full_cache)) deallocate(dg_frag%coef_pw_full_cache)
      dg_frag%coef_pw_full_cache_nstate = 0
      return
    end if
    if (.not. allocated(dg_frag%coef_pw) .or. .not. allocated(dg_frag%coef_pw_owner)) then
      if (allocated(dg_frag%coef_pw_full_cache)) deallocate(dg_frag%coef_pw_full_cache)
      dg_frag%coef_pw_full_cache_nstate = 0
      return
    end if

    n_pw = dg_frag%n_plane_waves
    nstate_cache = dg_frag%nstate_tot
    if (present(nstate_use)) nstate_cache = min(max(0, nstate_use), dg_frag%nstate_tot)
    if (n_pw <= 0) then
      if (allocated(dg_frag%coef_pw_full_cache)) deallocate(dg_frag%coef_pw_full_cache)
      dg_frag%coef_pw_full_cache_nstate = 0
      return
    end if
    if (nstate_cache <= 0) then
      if (allocated(dg_frag%coef_pw_full_cache)) deallocate(dg_frag%coef_pw_full_cache)
      dg_frag%coef_pw_full_cache_nstate = 0
      return
    end if

    if (.not. allocated(dg_frag%coef_pw_full_cache)) then
      allocate(dg_frag%coef_pw_full_cache(n_pw, nstate_cache, dg_frag%nspin))
    else if (size(dg_frag%coef_pw_full_cache, 1) /= n_pw .or. &
             size(dg_frag%coef_pw_full_cache, 2) /= nstate_cache .or. &
             size(dg_frag%coef_pw_full_cache, 3) /= dg_frag%nspin) then
      deallocate(dg_frag%coef_pw_full_cache)
      allocate(dg_frag%coef_pw_full_cache(n_pw, nstate_cache, dg_frag%nspin))
    end if

    env_cache_mode = ''
    call get_environment_variable('SALMON_DG_PW_CACHE_MODE', env_cache_mode, length=env_len, status=env_status)
    use_sum_cache = .true.
    compare_cache = .false.
    if (env_status == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_cache_mode(1:env_len))))
      case ('bcast','BCAST','broadcast','BROADCAST')
        use_sum_cache = .false.
      case ('sum','SUM','allreduce','ALLREDUCE','reduce','REDUCE')
        use_sum_cache = .true.
      case ('compare','COMPARE','check','CHECK')
        compare_cache = .true.
      end select
    end if
    env_cache_trace = ''
    call get_environment_variable('SALMON_DG_PW_CACHE_TRACE', env_cache_trace, length=env_len, status=env_status)
    trace_cache = .false.
    if (env_status == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_cache_trace(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        trace_cache = .true.
      end select
    end if

    dg_frag%coef_pw_full_cache(:, :, :) = (0.0d0, 0.0d0)
    if (use_sum_cache .or. compare_cache) then
      allocate(coef_pw_reduce(n_pw, nstate_cache, dg_frag%nspin))
      coef_pw_reduce(:, :, :) = (0.0d0, 0.0d0)
      do i = 1, n_pw
        if (i > size(dg_frag%coef_pw, 1)) cycle
        if (dg_frag%coef_pw_owner(i) /= dg_frag%id) cycle
        coef_pw_reduce(i, 1:nstate_cache, 1:dg_frag%nspin) = &
          dg_frag%coef_pw(i, 1:nstate_cache, 1:dg_frag%nspin)
      end do
      if (dg_frag%isize > 1 .and. dg_frag%icomm /= COMM_GROUP_NULL) then
        call comm_summation(coef_pw_reduce, dg_frag%coef_pw_full_cache, &
                            n_pw * nstate_cache * dg_frag%nspin, dg_frag%icomm)
      else
        dg_frag%coef_pw_full_cache(:, :, :) = coef_pw_reduce(:, :, :)
      end if
      if (compare_cache) then
        allocate(pw_row_ids(n_pw))
        do i = 1, n_pw
          pw_row_ids(i) = i
        end do
        coef_pw_reduce(:, :, :) = dg_frag%coef_pw_full_cache(:, :, :)
        dg_frag%coef_pw_full_cache(:, :, :) = (0.0d0, 0.0d0)
        call fetch_remote_coef_pw_rows(dg_frag, pw_row_ids, dg_frag%coef_pw_full_cache)
        cache_diff = maxval(abs(dg_frag%coef_pw_full_cache - coef_pw_reduce))
        cache_norm = max(1.0d-300, maxval(abs(dg_frag%coef_pw_full_cache)))
        if (dg_frag%id == 0 .or. trace_cache) then
          write(*,'(1x,a,i0,3(a,1pe12.4))') &
            '[DG-PW-CACHE-COMPARE] rank=', dg_frag%id, &
            ' max_abs_diff=', cache_diff, ' ref_norm=', cache_norm, &
            ' rel_diff=', cache_diff / cache_norm
          flush(6)
        end if
        deallocate(pw_row_ids)
      end if
      deallocate(coef_pw_reduce)
    else
      allocate(pw_row_ids(n_pw))
      do i = 1, n_pw
        pw_row_ids(i) = i
      end do
      call fetch_remote_coef_pw_rows(dg_frag, pw_row_ids, dg_frag%coef_pw_full_cache)
      deallocate(pw_row_ids)
    end if

    if (trace_cache) then
      cache_norm = maxval(abs(dg_frag%coef_pw_full_cache))
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1,a,l1,a,1pe12.4)') &
        '[DG-PW-CACHE] rank=', dg_frag%id, ' n_pw=', n_pw, &
        ' nstate=', nstate_cache, ' nspin=', dg_frag%nspin, &
        ' sum_mode=', use_sum_cache, ' compare=', compare_cache, &
        ' max_abs=', cache_norm
      flush(6)
    end if
    dg_frag%coef_pw_full_cache_nstate = nstate_cache
  end subroutine refresh_pw_coef_cache

  subroutine gather_fragment_coef_view(dg_frag, ispin, n_frag_rows, nstate_use, coef_frag, state_start, state_end)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: n_frag_rows
    integer, intent(in) :: nstate_use
    complex(8), allocatable, intent(inout) :: coef_frag(:,:)
    integer, intent(in), optional :: state_start, state_end

    integer :: i, cstart, cend
    integer, allocatable, save :: frag_row_ids(:)

    if (n_frag_rows < 0 .or. nstate_use < 0) then
      stop "DG gather_fragment_coef_view negative dimensions"
    end if
    if (n_frag_rows > dg_frag%n_mat_max) then
      stop "DG gather_fragment_coef_view invalid fragment row count"
    end if
    if (nstate_use > size(dg_frag%coef, 2)) then
      stop "DG gather_fragment_coef_view invalid state count"
    end if

    cstart = 1
    cend = nstate_use
    if (present(state_start)) cstart = state_start
    if (present(state_end)) cend = state_end
    if (cstart < 1 .or. cend < cstart .or. cend > size(dg_frag%coef, 2)) then
      stop "DG gather_fragment_coef_view invalid state range"
    end if
    if (cend - cstart + 1 /= nstate_use) then
      stop "DG gather_fragment_coef_view state range/shape mismatch"
    end if

    if (.not. allocated(coef_frag)) then
      allocate(coef_frag(max(0, n_frag_rows), max(0, nstate_use)))
    else if (size(coef_frag, 1) /= max(0, n_frag_rows) .or. size(coef_frag, 2) /= max(0, nstate_use)) then
      deallocate(coef_frag)
      allocate(coef_frag(max(0, n_frag_rows), max(0, nstate_use)))
    end if
    coef_frag(:, :) = (0.0d0, 0.0d0)
    if (n_frag_rows <= 0 .or. nstate_use <= 0) return

    if (.not. allocated(frag_row_ids)) then
      allocate(frag_row_ids(n_frag_rows))
    else if (size(frag_row_ids) /= n_frag_rows) then
      deallocate(frag_row_ids)
      allocate(frag_row_ids(n_frag_rows))
    end if
    do i = 1, n_frag_rows
      frag_row_ids(i) = i
    end do
    call fetch_remote_coef_rows(dg_frag, ispin, frag_row_ids, coef_frag, cstart, cend)
  end subroutine gather_fragment_coef_view

  subroutine gather_full_coef_view(dg_frag, ispin, n_frag_rows, nstate_use, coef_frag, coef_pw, state_start, state_end)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: n_frag_rows
    integer, intent(in) :: nstate_use
    complex(8), allocatable, intent(inout) :: coef_frag(:,:)
    complex(8), allocatable, intent(inout) :: coef_pw(:,:)
    integer, intent(in), optional :: state_start, state_end

    integer :: i, n_pw, cstart, cend
    integer, allocatable, save :: frag_row_ids(:), pw_row_ids(:)
    complex(8), allocatable, save :: coef_pw_all(:,:,:)
    integer :: ispin_eff
    logical, parameter :: enable_coef_gather_trace = .false.

    if (enable_coef_gather_trace .and. dg_frag%id == 0 .and. ispin == 1) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        coef gather entry: rank=", dg_frag%id, &
        " ispin=", ispin, " n_frag_rows=", n_frag_rows, " nstate_use=", nstate_use
      flush(6)
    end if
    if (n_frag_rows < 0 .or. nstate_use < 0) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] negative gather dimensions: rank=", dg_frag%id, &
        " n_frag_rows=", n_frag_rows, " nstate_use=", nstate_use
      stop "DG gather_full_coef_view negative dimensions"
    end if
    if (n_frag_rows > dg_frag%n_mat_max) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] gather rows exceed coef extent: rank=", dg_frag%id, &
        " n_frag_rows=", n_frag_rows, " n_mat_max=", dg_frag%n_mat_max
      stop "DG gather_full_coef_view invalid fragment row count"
    end if
    if (n_frag_rows >= dg_frag%n_mat_max .and. dg_frag%n_mat_max > 0) then
      write(*,'(1x,a,i0,a,i0)') "[FATAL] all-row coefficient gather is disabled: rank=", dg_frag%id, &
        " n_mat_max=", dg_frag%n_mat_max
      stop "DG-Fragment RT: all-row coefficient gather is disabled"
    end if
    if (nstate_use > size(dg_frag%coef, 2)) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] gather states exceed coef extent: rank=", dg_frag%id, &
        " nstate_use=", nstate_use, " coef_cols=", size(dg_frag%coef, 2)
      stop "DG gather_full_coef_view invalid state count"
    end if
    cstart = 1
    cend = nstate_use
    if (present(state_start)) cstart = state_start
    if (present(state_end)) cend = state_end
    if (cstart < 1 .or. cend < cstart .or. cend > size(dg_frag%coef, 2)) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] gather state-range invalid: rank=", dg_frag%id, &
        " cstart=", cstart, " cend=", cend, " coef_cols=", size(dg_frag%coef, 2)
      stop "DG gather_full_coef_view invalid state range"
    end if
    if (cend - cstart + 1 /= nstate_use) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] gather state-range/shape mismatch: rank=", dg_frag%id, &
        " nstate_use=", nstate_use, " cstart=", cstart, " cend=", cend
      stop "DG gather_full_coef_view state range/shape mismatch"
    end if

    if (.not. allocated(coef_frag)) then
      allocate(coef_frag(max(0, n_frag_rows), max(0, nstate_use)))
    else if (size(coef_frag, 1) /= max(0, n_frag_rows) .or. size(coef_frag, 2) /= max(0, nstate_use)) then
      deallocate(coef_frag)
      allocate(coef_frag(max(0, n_frag_rows), max(0, nstate_use)))
    end if
    coef_frag(:, :) = (0.0d0, 0.0d0)
    if (n_frag_rows > 0 .and. nstate_use > 0) then
      if (.not. allocated(frag_row_ids)) then
        allocate(frag_row_ids(n_frag_rows))
      else if (size(frag_row_ids) /= n_frag_rows) then
        deallocate(frag_row_ids)
        allocate(frag_row_ids(n_frag_rows))
      end if
      do i = 1, n_frag_rows
        frag_row_ids(i) = i
      end do
      call fetch_remote_coef_rows(dg_frag, ispin, frag_row_ids, coef_frag, cstart, cend)
    end if

    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    if (.not. allocated(coef_pw)) then
      allocate(coef_pw(max(0, n_pw), max(0, nstate_use)))
    else if (size(coef_pw, 1) /= max(0, n_pw) .or. size(coef_pw, 2) /= max(0, nstate_use)) then
      deallocate(coef_pw)
      allocate(coef_pw(max(0, n_pw), max(0, nstate_use)))
    end if
    coef_pw(:, :) = (0.0d0, 0.0d0)
    if (n_pw > 0 .and. nstate_use > 0) then
      if (.not. allocated(pw_row_ids)) then
        allocate(pw_row_ids(n_pw))
      else if (size(pw_row_ids) /= n_pw) then
        deallocate(pw_row_ids)
        allocate(pw_row_ids(n_pw))
      end if
      do i = 1, n_pw
        pw_row_ids(i) = i
      end do
      ispin_eff = min(max(ispin, 1), dg_frag%nspin)
      if (.not. allocated(coef_pw_all)) then
        allocate(coef_pw_all(n_pw, nstate_use, 1))
      else if (size(coef_pw_all, 1) /= n_pw .or. size(coef_pw_all, 2) /= nstate_use .or. &
               size(coef_pw_all, 3) /= 1) then
        deallocate(coef_pw_all)
        allocate(coef_pw_all(n_pw, nstate_use, 1))
      end if
      coef_pw_all(:, :, :) = (0.0d0, 0.0d0)
      call fetch_remote_coef_pw_rows(dg_frag, pw_row_ids, coef_pw_all, cstart, cend, ispin_eff)
      coef_pw(:, :) = coef_pw_all(:, :, 1)
    end if
  end subroutine gather_full_coef_view

  subroutine zero_nonowned_coefficients(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ispin, i, global_idx

    if (allocated(dg_frag%coef_owner) .and. allocated(dg_frag%coef)) then
      do ispin = 1, min(size(dg_frag%coef_owner, 2), size(dg_frag%coef, 3))
        if (allocated(dg_frag%local_coef_global_ids)) then
          do i = 1, size(dg_frag%coef, 1)
            global_idx = dg_frag%local_coef_global_ids(i, ispin)
            if (global_idx >= 1 .and. global_idx <= size(dg_frag%coef_owner, 1)) then
              if (dg_frag%coef_owner(global_idx, ispin) == dg_frag%id) cycle
            end if
            dg_frag%coef(i, :, ispin) = (0.0d0, 0.0d0)
          end do
        else
          do i = 1, min(size(dg_frag%coef_owner, 1), size(dg_frag%coef, 1))
            if (dg_frag%coef_owner(i, ispin) == dg_frag%id) cycle
            dg_frag%coef(i, :, ispin) = (0.0d0, 0.0d0)
          end do
        end if
      end do
    end if

    if (allocated(dg_frag%coef_pw_owner) .and. allocated(dg_frag%coef_pw)) then
      do i = 1, min(size(dg_frag%coef_pw_owner), size(dg_frag%coef_pw, 1))
        if (dg_frag%coef_pw_owner(i) == dg_frag%id) cycle
        dg_frag%coef_pw(i, :, :) = (0.0d0, 0.0d0)
      end do
    end if
  end subroutine zero_nonowned_coefficients

  subroutine rebuild_local_h_block_ids(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    logical, allocatable :: row_is_local(:)
    integer :: ifrag, iblk, io, ispin, n_local_rows, n_local_blocks
    integer :: global_idx

    if (dg_frag%H_local_block_ids_valid .and. allocated(dg_frag%H_local_block_ids) .and. allocated(dg_frag%H_local_rows)) return
    call invalidate_coef_exchange_cache(dg_frag)
    if (allocated(dg_frag%H_local_block_ids)) deallocate(dg_frag%H_local_block_ids)
    if (allocated(dg_frag%H_local_rows)) deallocate(dg_frag%H_local_rows)
    dg_frag%H_local_block_ids_valid = .false.
    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return

    allocate(row_is_local(dg_frag%n_frag))
    row_is_local = .false.
    n_local_rows = 0
    do ifrag = 1, dg_frag%n_frag
      do ispin = 1, dg_frag%nspin
        if (dg_frag%n_basis(ifrag, ispin) <= 0) cycle
        do io = 1, min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
          global_idx = dg_frag%index_basis(io, ifrag, ispin)
          if (global_idx < 1 .or. global_idx > size(dg_frag%coef_owner, 1)) cycle
          if (dg_frag%coef_owner(global_idx, ispin) /= dg_frag%id) cycle
          row_is_local(ifrag) = .true.
          exit
        end do
        if (row_is_local(ifrag)) exit
      end do
      if (row_is_local(ifrag)) n_local_rows = n_local_rows + 1
    end do

    if (n_local_rows > 0) then
      allocate(dg_frag%H_local_rows(n_local_rows))
      n_local_rows = 0
      do ifrag = 1, dg_frag%n_frag
        if (.not. row_is_local(ifrag)) cycle
        n_local_rows = n_local_rows + 1
        dg_frag%H_local_rows(n_local_rows) = ifrag
      end do
    end if

    n_local_blocks = 0
    do iblk = 1, size(dg_frag%H_mat_blocks)
      if (dg_frag%H_mat_blocks(iblk)%ifrag_row < 1 .or. dg_frag%H_mat_blocks(iblk)%ifrag_row > dg_frag%n_frag) cycle
      if (.not. row_is_local(dg_frag%H_mat_blocks(iblk)%ifrag_row)) cycle
      n_local_blocks = n_local_blocks + 1
    end do

    if (n_local_blocks > 0) then
      allocate(dg_frag%H_local_block_ids(n_local_blocks))
      n_local_blocks = 0
      do iblk = 1, size(dg_frag%H_mat_blocks)
        if (dg_frag%H_mat_blocks(iblk)%ifrag_row < 1 .or. dg_frag%H_mat_blocks(iblk)%ifrag_row > dg_frag%n_frag) cycle
        if (.not. row_is_local(dg_frag%H_mat_blocks(iblk)%ifrag_row)) cycle
        n_local_blocks = n_local_blocks + 1
        dg_frag%H_local_block_ids(n_local_blocks) = iblk
      end do
    end if

    deallocate(row_is_local)
    dg_frag%H_local_block_ids_valid = .true.
  end subroutine rebuild_local_h_block_ids

  subroutine rebuild_local_nl_block_ids(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    logical, allocatable :: row_is_local(:)
    integer :: ifrag, iblk, io, ispin, n_local_blocks
    integer :: global_idx

    if (allocated(dg_frag%H_nl_local_block_ids)) deallocate(dg_frag%H_nl_local_block_ids)
    if (.not. allocated(dg_frag%H_nl_blocks)) return
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%n_basis)) return
    if (.not. allocated(dg_frag%index_basis)) return

    allocate(row_is_local(dg_frag%n_frag))
    row_is_local = .false.
    if (dg_frag%coef_state_block_mode) then
      if (dg_frag%ifrag_group >= 1 .and. dg_frag%ifrag_group <= dg_frag%n_frag) then
        row_is_local(dg_frag%ifrag_group) = .true.
      end if
    else
      do ifrag = 1, dg_frag%n_frag
        do ispin = 1, dg_frag%nspin
          if (dg_frag%n_basis(ifrag, ispin) <= 0) cycle
          do io = 1, min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            if (global_idx < 1 .or. global_idx > size(dg_frag%coef_owner, 1)) cycle
            if (dg_frag%coef_owner(global_idx, ispin) /= dg_frag%id) cycle
            row_is_local(ifrag) = .true.
            exit
          end do
          if (row_is_local(ifrag)) exit
        end do
      end do
    end if

    n_local_blocks = 0
    do iblk = 1, size(dg_frag%H_nl_blocks)
      if (dg_frag%H_nl_blocks(iblk)%ifrag_row < 1 .or. dg_frag%H_nl_blocks(iblk)%ifrag_row > dg_frag%n_frag) cycle
      if (.not. row_is_local(dg_frag%H_nl_blocks(iblk)%ifrag_row)) cycle
      n_local_blocks = n_local_blocks + 1
    end do

    if (n_local_blocks > 0) then
      allocate(dg_frag%H_nl_local_block_ids(n_local_blocks))
      n_local_blocks = 0
      do iblk = 1, size(dg_frag%H_nl_blocks)
        if (dg_frag%H_nl_blocks(iblk)%ifrag_row < 1 .or. dg_frag%H_nl_blocks(iblk)%ifrag_row > dg_frag%n_frag) cycle
        if (.not. row_is_local(dg_frag%H_nl_blocks(iblk)%ifrag_row)) cycle
        n_local_blocks = n_local_blocks + 1
        dg_frag%H_nl_local_block_ids(n_local_blocks) = iblk
      end do
    end if

    deallocate(row_is_local)
  end subroutine rebuild_local_nl_block_ids

  subroutine zero_nonlocal_h_matrix_blocks(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: iblk

    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%H_local_block_ids)) then
!$omp parallel do private(iblk) schedule(static)
      do iblk = 1, size(dg_frag%H_mat_blocks)
        dg_frag%H_mat_blocks(iblk)%val(:, :, :) = 0.0d0
        if (allocated(dg_frag%H_mat_kinetic_blocks) .and. iblk <= size(dg_frag%H_mat_kinetic_blocks)) then
          dg_frag%H_mat_kinetic_blocks(iblk)%val(:, :, :) = 0.0d0
        end if
      end do
!$omp end parallel do
      return
    end if

!$omp parallel do private(iblk) schedule(static)
    do iblk = 1, size(dg_frag%H_mat_blocks)
      if (any(dg_frag%H_local_block_ids == iblk)) cycle
      dg_frag%H_mat_blocks(iblk)%val(:, :, :) = 0.0d0
      if (allocated(dg_frag%H_mat_kinetic_blocks) .and. iblk <= size(dg_frag%H_mat_kinetic_blocks)) then
        dg_frag%H_mat_kinetic_blocks(iblk)%val(:, :, :) = 0.0d0
      end if
    end do
!$omp end parallel do
  end subroutine zero_nonlocal_h_matrix_blocks

  subroutine apply_matrix_blocks(dg_frag, blocks, ispin, x, y)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:)
    complex(8), intent(inout) :: y(:)

    integer :: iblk, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ig_i, ig_j, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))
    complex(8) :: xj

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return

    do iblk = 1, size(blocks)
      ifrag_row = blocks(iblk)%ifrag_row
      ifrag_col = blocks(iblk)%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
      nrow = dg_frag%n_basis(ifrag_row, ispin)
      ncol = dg_frag%n_basis(ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) cycle

      valid_row_count = 0
      do ii = 1, nrow
        row_gid(ii) = dg_frag%index_basis(ii, ifrag_row, ispin)
        if (row_gid(ii) < 1 .or. row_gid(ii) > size(y)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ii
      end do
      if (valid_row_count <= 0) cycle
      valid_col_count = 0
      do jj = 1, ncol
        col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (col_gid(jj) < 1 .or. col_gid(jj) > size(x)) cycle
        valid_col_count = valid_col_count + 1
        valid_col_ids(valid_col_count) = jj
      end do
      if (valid_col_count <= 0) cycle

      do idx_jj = 1, valid_col_count
        jj = valid_col_ids(idx_jj)
        ig_j = col_gid(jj)
        xj = x(ig_j)
!$omp simd private(ii,ig_i)
        do idx_ii = 1, valid_row_count
          ii = valid_row_ids(idx_ii)
          ig_i = row_gid(ii)
          y(ig_i) = y(ig_i) + blocks(iblk)%val(ii, jj, ispin) * xj
        end do
      end do
    end do
  end subroutine apply_matrix_blocks

  subroutine apply_matrix_blocks_batch(dg_frag, blocks, ispin, x, y, block_ids)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:, :)
    complex(8), intent(inout) :: y(:, :)
    integer, intent(in), optional :: block_ids(:)

    integer, parameter :: block_apply_target_mb = 128
    integer :: iblk, iblk_idx, nstate, max_basis, state_chunk
    integer(8) :: target_bytes, bytes_per_state
    complex(8), allocatable :: mat_block(:, :), x_block(:, :), y_block(:, :)

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return
    nstate = min(size(x, 2), size(y, 2))
    if (nstate <= 0) return
    max_basis = size(dg_frag%index_basis, 1)
    if (allocated(dg_frag%n_basis)) then
      max_basis = max(1, min(max_basis, maxval(dg_frag%n_basis(:, ispin))))
    end if
    target_bytes = int(block_apply_target_mb, kind=8) * 1024_8 * 1024_8
    bytes_per_state = int(max(1, 2 * max_basis), kind=8) * 16_8
    state_chunk = max(1, min(nstate, int(max(1_8, target_bytes / bytes_per_state))))
    allocate(mat_block(max_basis, max_basis), x_block(max_basis, state_chunk), y_block(max_basis, state_chunk))

    if (present(block_ids)) then
      do iblk_idx = 1, size(block_ids)
        iblk = block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(blocks)) cycle
        call apply_one_real_block(iblk)
      end do
    else
      do iblk = 1, size(blocks)
        call apply_one_real_block(iblk)
      end do
    end if
    deallocate(mat_block, x_block, y_block)

  contains
    subroutine apply_one_real_block(iblk_in)
      implicit none
      integer, intent(in) :: iblk_in
      integer :: ifrag_row, ifrag_col, nrow, ncol
      integer :: ii, jj, idx_ii, idx_jj, ig_i, ig_j, valid_row_count, valid_col_count
      integer :: state0, state1, nstate_part
      integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
      integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))

      ifrag_row = blocks(iblk_in)%ifrag_row
      ifrag_col = blocks(iblk_in)%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) return
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) return
      nrow = dg_frag%n_basis(ifrag_row, ispin)
      ncol = dg_frag%n_basis(ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) return

      valid_row_count = 0
      do ii = 1, nrow
        row_gid(ii) = dg_frag%index_basis(ii, ifrag_row, ispin)
        if (row_gid(ii) < 1 .or. row_gid(ii) > size(y, 1)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ii
      end do
      valid_col_count = 0
      do jj = 1, ncol
        col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (col_gid(jj) < 1 .or. col_gid(jj) > size(x, 1)) cycle
        valid_col_count = valid_col_count + 1
        valid_col_ids(valid_col_count) = jj
      end do
      if (valid_row_count <= 0 .or. valid_col_count <= 0) return

      do idx_jj = 1, valid_col_count
        jj = valid_col_ids(idx_jj)
        do idx_ii = 1, valid_row_count
          ii = valid_row_ids(idx_ii)
          mat_block(idx_ii, idx_jj) = cmplx(blocks(iblk_in)%val(ii, jj, ispin), 0.0d0, kind=8)
        end do
      end do

      do state0 = 1, nstate, state_chunk
        state1 = min(nstate, state0 + state_chunk - 1)
        nstate_part = state1 - state0 + 1
        do idx_jj = 1, valid_col_count
          jj = valid_col_ids(idx_jj)
          ig_j = col_gid(jj)
          x_block(idx_jj, 1:nstate_part) = x(ig_j, state0:state1)
        end do
        call zgemm('N', 'N', valid_row_count, nstate_part, valid_col_count, (1.0d0, 0.0d0), &
                   mat_block, max_basis, x_block, max_basis, (0.0d0, 0.0d0), y_block, max_basis)
        do idx_ii = 1, valid_row_count
          ii = valid_row_ids(idx_ii)
          ig_i = row_gid(ii)
          y(ig_i, state0:state1) = y(ig_i, state0:state1) + y_block(idx_ii, 1:nstate_part)
        end do
      end do
    end subroutine apply_one_real_block
  end subroutine apply_matrix_blocks_batch

  subroutine apply_complex_matrix_blocks_batch(dg_frag, blocks, ispin, x, y, block_ids)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(complex_matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:, :)
    complex(8), intent(inout) :: y(:, :)
    integer, intent(in), optional :: block_ids(:)

    integer, parameter :: block_apply_target_mb = 128
    integer :: iblk, iblk_idx, nstate, max_basis, state_chunk
    integer(8) :: target_bytes, bytes_per_state
    complex(8), allocatable :: mat_block(:, :), x_block(:, :), y_block(:, :)

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return
    nstate = min(size(x, 2), size(y, 2))
    if (nstate <= 0) return
    max_basis = size(dg_frag%index_basis, 1)
    if (allocated(dg_frag%n_basis)) then
      max_basis = max(1, min(max_basis, maxval(dg_frag%n_basis(:, ispin))))
    end if
    target_bytes = int(block_apply_target_mb, kind=8) * 1024_8 * 1024_8
    bytes_per_state = int(max(1, 2 * max_basis), kind=8) * 16_8
    state_chunk = max(1, min(nstate, int(max(1_8, target_bytes / bytes_per_state))))
    allocate(mat_block(max_basis, max_basis), x_block(max_basis, state_chunk), y_block(max_basis, state_chunk))

    if (present(block_ids)) then
      do iblk_idx = 1, size(block_ids)
        iblk = block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(blocks)) cycle
        call apply_one_complex_block(iblk)
      end do
    else
      do iblk = 1, size(blocks)
        call apply_one_complex_block(iblk)
      end do
    end if
    deallocate(mat_block, x_block, y_block)

  contains
    subroutine apply_one_complex_block(iblk_in)
      implicit none
      integer, intent(in) :: iblk_in
      integer :: ifrag_row, ifrag_col, nrow, ncol
      integer :: ii, jj, idx_ii, idx_jj, ig_i, ig_j, valid_row_count, valid_col_count
      integer :: state0, state1, nstate_part
      integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
      integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))

      ifrag_row = blocks(iblk_in)%ifrag_row
      ifrag_col = blocks(iblk_in)%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) return
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) return
      nrow = dg_frag%n_basis(ifrag_row, ispin)
      ncol = dg_frag%n_basis(ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) return

      valid_row_count = 0
      do ii = 1, nrow
        row_gid(ii) = dg_frag%index_basis(ii, ifrag_row, ispin)
        if (row_gid(ii) < 1 .or. row_gid(ii) > size(y, 1)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ii
      end do
      valid_col_count = 0
      do jj = 1, ncol
        col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (col_gid(jj) < 1 .or. col_gid(jj) > size(x, 1)) cycle
        valid_col_count = valid_col_count + 1
        valid_col_ids(valid_col_count) = jj
      end do
      if (valid_row_count <= 0 .or. valid_col_count <= 0) return

      do idx_jj = 1, valid_col_count
        jj = valid_col_ids(idx_jj)
        do idx_ii = 1, valid_row_count
          ii = valid_row_ids(idx_ii)
          mat_block(idx_ii, idx_jj) = blocks(iblk_in)%val(ii, jj, ispin)
        end do
      end do

      do state0 = 1, nstate, state_chunk
        state1 = min(nstate, state0 + state_chunk - 1)
        nstate_part = state1 - state0 + 1
        do idx_jj = 1, valid_col_count
          jj = valid_col_ids(idx_jj)
          ig_j = col_gid(jj)
          x_block(idx_jj, 1:nstate_part) = x(ig_j, state0:state1)
        end do
        call zgemm('N', 'N', valid_row_count, nstate_part, valid_col_count, (1.0d0, 0.0d0), &
                   mat_block, max_basis, x_block, max_basis, (0.0d0, 0.0d0), y_block, max_basis)
        do idx_ii = 1, valid_row_count
          ii = valid_row_ids(idx_ii)
          ig_i = row_gid(ii)
          y(ig_i, state0:state1) = y(ig_i, state0:state1) + y_block(idx_ii, 1:nstate_part)
        end do
      end do
    end subroutine apply_one_complex_block
  end subroutine apply_complex_matrix_blocks_batch

  subroutine copy_matrix_blocks_to_complex_dense(dg_frag, blocks, ispin, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    complex(8), intent(inout) :: mat(:, :)

    integer :: iblk, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ig_i, ig_j, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return

    do iblk = 1, size(blocks)
      ifrag_row = blocks(iblk)%ifrag_row
      ifrag_col = blocks(iblk)%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
      nrow = dg_frag%n_basis(ifrag_row, ispin)
      ncol = dg_frag%n_basis(ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) cycle

      valid_row_count = 0
      do ii = 1, nrow
        row_gid(ii) = dg_frag%index_basis(ii, ifrag_row, ispin)
        if (row_gid(ii) < 1 .or. row_gid(ii) > size(mat, 1)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ii
      end do
      if (valid_row_count <= 0) cycle

      valid_col_count = 0
      do jj = 1, ncol
        col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (col_gid(jj) < 1 .or. col_gid(jj) > size(mat, 2)) cycle
        valid_col_count = valid_col_count + 1
        valid_col_ids(valid_col_count) = jj
      end do
      if (valid_col_count <= 0) cycle

      do idx_jj = 1, valid_col_count
        jj = valid_col_ids(idx_jj)
        ig_j = col_gid(jj)
!$omp simd private(ii,ig_i)
        do idx_ii = 1, valid_row_count
          ii = valid_row_ids(idx_ii)
          ig_i = row_gid(ii)
          mat(ig_i, ig_j) = cmplx(blocks(iblk)%val(ii, jj, ispin), 0.0d0, kind=8)
        end do
      end do
    end do
  end subroutine copy_matrix_blocks_to_complex_dense

  !=======================================================================
  ! Apply gradient operator to a basis function using finite differences
  !=======================================================================
  subroutine apply_gradient_to_basis(dg_frag, i_local, jo, mg, stencil, grad_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    real(8),                intent(out) :: grad_phi(:,:,:,:)

    integer :: lx, ly, lz, ifrag
    integer :: ixg, iyg, izg
    integer :: ix0, iy0, iz0
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3
    integer :: ixp1, ixm1, ixp2, ixm2, ixp3, ixm3, ixp4, ixm4
    integer :: iyp1, iym1, iyp2, iym2, iyp3, iym3, iyp4, iym4
    integer :: izp1, izm1, izp2, izm2, izp3, izm3, izp4, izm4
    real(8) :: gx, gy, gz
    real(8) :: nabt(4,3)
    integer :: ndom(3)
    logical :: use_direct_halo_fd

    nabt = stencil%coef_nab
    ifrag = dg_frag%ifrag_start + i_local - 1
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)
    use_direct_halo_fd = all(dg_frag%nxyz_buffer(1:3) >= 4)

    grad_phi = 0.0d0

    if (use_direct_halo_fd) then
      !$omp parallel do collapse(2) private(lx, ly, lz, ixg, iyg, izg, ix0, iy0, iz0, gx, gy, gz, &
      !$omp& ixp1, ixm1, ixp2, ixm2, ixp3, ixm3, ixp4, ixm4, &
      !$omp& iyp1, iym1, iyp2, iym2, iyp3, iym3, iyp4, iym4, &
      !$omp& izp1, izm1, izp2, izm2, izp3, izm3, izp4, izm4) schedule(static)
      do lz = 1, ndom(3)
        do ly = 1, ndom(2)
          !$omp simd private(gx, gy, gz)
          do lx = 1, ndom(1)
            ixg = modulo(dg_frag%ixyz_frag(1, ifrag) + lx - 2, dg_frag%lgnum_total(1)) + 1
            iyg = modulo(dg_frag%ixyz_frag(2, ifrag) + ly - 2, dg_frag%lgnum_total(2)) + 1
            izg = modulo(dg_frag%ixyz_frag(3, ifrag) + lz - 2, dg_frag%lgnum_total(3)) + 1
            ix0 = map_global_to_phi_box_coord(ixg, phi_lb1, phi_ub1, dg_frag%lgnum_total(1))
            iy0 = map_global_to_phi_box_coord(iyg, phi_lb2, phi_ub2, dg_frag%lgnum_total(2))
            iz0 = map_global_to_phi_box_coord(izg, phi_lb3, phi_ub3, dg_frag%lgnum_total(3))

            ixp1 = ix0 + 1
            ixm1 = ix0 - 1
            ixp2 = ix0 + 2
            ixm2 = ix0 - 2
            ixp3 = ix0 + 3
            ixm3 = ix0 - 3
            ixp4 = ix0 + 4
            ixm4 = ix0 - 4
            iyp1 = iy0 + 1
            iym1 = iy0 - 1
            iyp2 = iy0 + 2
            iym2 = iy0 - 2
            iyp3 = iy0 + 3
            iym3 = iy0 - 3
            iyp4 = iy0 + 4
            iym4 = iy0 - 4
            izp1 = iz0 + 1
            izm1 = iz0 - 1
            izp2 = iz0 + 2
            izm2 = iz0 - 2
            izp3 = iz0 + 3
            izm3 = iz0 - 3
            izp4 = iz0 + 4
            izm4 = iz0 - 4

            gx = nabt(1,1) * (dg_frag%phi_frag(ixp1, iy0, iz0, jo, i_local) - dg_frag%phi_frag(ixm1, iy0, iz0, jo, i_local)) + &
                 nabt(2,1) * (dg_frag%phi_frag(ixp2, iy0, iz0, jo, i_local) - dg_frag%phi_frag(ixm2, iy0, iz0, jo, i_local)) + &
                 nabt(3,1) * (dg_frag%phi_frag(ixp3, iy0, iz0, jo, i_local) - dg_frag%phi_frag(ixm3, iy0, iz0, jo, i_local)) + &
                 nabt(4,1) * (dg_frag%phi_frag(ixp4, iy0, iz0, jo, i_local) - dg_frag%phi_frag(ixm4, iy0, iz0, jo, i_local))

            gy = nabt(1,2) * (dg_frag%phi_frag(ix0, iyp1, iz0, jo, i_local) - dg_frag%phi_frag(ix0, iym1, iz0, jo, i_local)) + &
                 nabt(2,2) * (dg_frag%phi_frag(ix0, iyp2, iz0, jo, i_local) - dg_frag%phi_frag(ix0, iym2, iz0, jo, i_local)) + &
                 nabt(3,2) * (dg_frag%phi_frag(ix0, iyp3, iz0, jo, i_local) - dg_frag%phi_frag(ix0, iym3, iz0, jo, i_local)) + &
                 nabt(4,2) * (dg_frag%phi_frag(ix0, iyp4, iz0, jo, i_local) - dg_frag%phi_frag(ix0, iym4, iz0, jo, i_local))

            gz = nabt(1,3) * (dg_frag%phi_frag(ix0, iy0, izp1, jo, i_local) - dg_frag%phi_frag(ix0, iy0, izm1, jo, i_local)) + &
                 nabt(2,3) * (dg_frag%phi_frag(ix0, iy0, izp2, jo, i_local) - dg_frag%phi_frag(ix0, iy0, izm2, jo, i_local)) + &
                 nabt(3,3) * (dg_frag%phi_frag(ix0, iy0, izp3, jo, i_local) - dg_frag%phi_frag(ix0, iy0, izm3, jo, i_local)) + &
                 nabt(4,3) * (dg_frag%phi_frag(ix0, iy0, izp4, jo, i_local) - dg_frag%phi_frag(ix0, iy0, izm4, jo, i_local))

            grad_phi(lx, ly, lz, 1) = gx
            grad_phi(lx, ly, lz, 2) = gy
            grad_phi(lx, ly, lz, 3) = gz
          end do
        end do
      end do
      !$omp end parallel do
    else
      stop "DG-Fragment RT: apply_gradient_to_basis requires nxyz_buffer>=4"
    end if

  end subroutine apply_gradient_to_basis

  subroutine apply_momentum_blocks(dg_frag, ispin, scale_vec, x, y, row_frag_ids)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: ispin
    real(8),                intent(in) :: scale_vec(3)
    complex(8),             intent(in) :: x(:,:)
    complex(8),             intent(inout) :: y(:,:)
    integer, optional,      intent(in) :: row_frag_ids(:)

    integer :: iblk, idir, ib, jb, row_idx, col_idx, nstate, istate
    integer :: active_dir_count, valid_row_count, valid_col_count, idx_dir, idx_ib, idx_jb
    integer :: ifrag_row, ifrag_col, nrow, ncol
    integer :: active_dirs(3), valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    real(8) :: scale

    if (.not. allocated(dg_frag%momentum_blocks)) return
    nstate = min(size(x, 2), size(y, 2))
    active_dir_count = 0
    do idir = 1, 3
      if (abs(scale_vec(idir)) < 1.0d-30) cycle
      active_dir_count = active_dir_count + 1
      active_dirs(active_dir_count) = idir
    end do
    if (active_dir_count <= 0) return

!$omp parallel do schedule(static) &
!$omp& private(istate, iblk, idir, scale, jb, col_idx, ib, row_idx, valid_row_count, valid_col_count, &
!$omp&         idx_dir, idx_ib, idx_jb, ifrag_row, ifrag_col, nrow, ncol, valid_row_ids, valid_col_ids, &
!$omp&         row_gid, col_gid)
    do istate = 1, nstate
      do iblk = 1, dg_frag%n_momentum_blocks
        if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle

        ifrag_row = dg_frag%momentum_blocks(iblk)%ifrag_row
        ifrag_col = dg_frag%momentum_blocks(iblk)%ifrag_col
        if (present(row_frag_ids)) then
          if (.not. any(row_frag_ids == ifrag_row)) cycle
        end if
        nrow = dg_frag%n_basis(ifrag_row, ispin)
        ncol = dg_frag%n_basis(ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle

        valid_row_count = 0
        do ib = 1, nrow
          row_gid(ib) = dg_frag%index_basis(ib, ifrag_row, ispin)
          if (row_gid(ib) < 1 .or. row_gid(ib) > size(y, 1)) cycle
          valid_row_count = valid_row_count + 1
          valid_row_ids(valid_row_count) = ib
        end do
        if (valid_row_count <= 0) cycle

        valid_col_count = 0
        do jb = 1, ncol
          col_gid(jb) = dg_frag%index_basis(jb, ifrag_col, ispin)
          if (col_gid(jb) < 1 .or. col_gid(jb) > size(x, 1)) cycle
          valid_col_count = valid_col_count + 1
          valid_col_ids(valid_col_count) = jb
        end do
        if (valid_col_count <= 0) cycle

        do idx_dir = 1, active_dir_count
          idir = active_dirs(idx_dir)
          scale = scale_vec(idir)
          do idx_jb = 1, valid_col_count
            jb = valid_col_ids(idx_jb)
            col_idx = col_gid(jb)
!$omp simd private(ib,row_idx)
            do idx_ib = 1, valid_row_count
              ib = valid_row_ids(idx_ib)
              row_idx = row_gid(ib)
              y(row_idx, istate) = y(row_idx, istate) + &
                scale * dg_frag%momentum_blocks(iblk)%val(idir, ib, jb, ispin) * x(col_idx, istate)
            end do
          end do
        end do
      end do
    end do
!$omp end parallel do

  end subroutine apply_momentum_blocks

  subroutine calculate_macroscopic_current_from_momentum_blocks(dg_frag, system, current_raw)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(out)   :: current_raw(3)

    integer :: ispin, iblk, ifrag_row, ifrag_col, nrow, ncol, io, jo, idir
    integer :: nocc_spin, occ0, nbatch, state_work, occ_first, occ_last, local_col, local_idx
    integer :: gid_i, gid_j, idx_row, n_global, n_needed, pos_i, pos_j
    integer :: max_basis, n_output, iout, n_active_blocks, iactive
    integer :: io_block, jo_block, nrow_active, ncol_active
    integer, allocatable, save :: needed_ids(:), needed_pos(:), output_pos(:)
    integer, allocatable, save :: block_ids(:), block_nrow(:), block_ncol(:)
    integer, allocatable, save :: block_row_lid(:,:), block_row_pos(:,:)
    integer, allocatable, save :: block_col_lid(:,:), block_col_pos(:,:)
    logical, allocatable, save :: row_needed(:), row_output(:)
    complex(8), allocatable, save :: coef_work(:,:), pc_work(:,:,:)
    complex(8), allocatable, save :: coef_col_work(:,:), pc_block_work(:,:), mom_block_work(:,:)
    complex(8), allocatable, save :: pc_self_work(:,:,:), pc_cross_work(:,:,:)
    real(8) :: curr_local(3), curr_sum(3), pair_x, pair_y, pair_z, occ_factor
    real(8) :: pair_self_x, pair_self_y, pair_self_z, pair_cross_x, pair_cross_y, pair_cross_z
    real(8) :: curr_x, curr_y, curr_z, curr_self_x, curr_self_y, curr_self_z, curr_cross_x, curr_cross_y, curr_cross_z
    real(8) :: curr_self_local(3), curr_cross_local(3), curr_self_sum(3), curr_cross_sum(3)
    real(8) :: block_val
    logical :: have_owned_rows
    logical, save :: trace_decomp_initialized = .false.
    logical, save :: trace_decomp = .false.
    character(64) :: env_trace_decomp
    integer :: env_status
    integer, parameter :: current_coef_cache_target_mb = 64
    integer(8) :: target_bytes, bytes_per_state

    current_raw(:) = 0.0d0
    curr_local(:) = 0.0d0
    curr_sum(:) = 0.0d0
    curr_self_local(:) = 0.0d0
    curr_cross_local(:) = 0.0d0
    dg_frag%current_momentum_self_raw(:) = 0.0d0
    dg_frag%current_momentum_cross_raw(:) = 0.0d0
    dg_frag%current_momentum_decomp_ready = .false.
    if (.not. trace_decomp_initialized) then
      env_trace_decomp = ''
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_PROP_OBS_SERIES_CURRENT_DETAIL', &
        env_trace_decomp, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_trace_decomp)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          trace_decomp = .true.
        end select
      end if
      env_trace_decomp = ''
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_CURRENT_OP_DIAG', env_trace_decomp, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_trace_decomp)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          trace_decomp = .true.
        end select
      end if
      env_trace_decomp = ''
      call get_environment_variable('SALMON_DG_MIXED_Z_LOCAL_CURRENT_OP_ALL', env_trace_decomp, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_trace_decomp)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          trace_decomp = .true.
        end select
      end if
      trace_decomp_initialized = .true.
    end if
    if (.not. allocated(dg_frag%momentum_blocks)) then
      call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
      current_raw(:) = curr_sum(:)
      return
    end if
    if (dg_frag%dc_lcfo_seed_basis_cleaned .and. .not. dg_frag%momentum_blocks_include_dg_flux) then
      stop "DG current requires covariant Flux contribution in DG velocity blocks"
    end if
    if (.not. allocated(dg_frag%index_basis)) then
      call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
      current_raw(:) = curr_sum(:)
      return
    end if
    if (.not. allocated(dg_frag%coef_global_to_local)) then
      call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
      current_raw(:) = curr_sum(:)
      return
    end if

    n_global = size(dg_frag%coef_owner, 1)
    if (.not. allocated(row_needed)) then
      allocate(row_needed(n_global), row_output(n_global), needed_pos(n_global))
    else if (size(row_needed, 1) /= n_global) then
      deallocate(row_needed, row_output, needed_pos)
      allocate(row_needed(n_global), row_output(n_global), needed_pos(n_global))
    end if
    do ispin = 1, min(dg_frag%nspin, system%nspin)
      nocc_spin = 0
      if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
      if (dg_frag%coef_state_block_mode) then
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(system%rocc, 1))
      else
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(dg_frag%coef, 2), size(system%rocc, 1))
      end if
      if (nocc_spin <= 0) cycle
      occ_first = 1
      occ_last = nocc_spin
      if (dg_frag%coef_state_block_mode) then
        occ_first = max(1, dg_frag%coef_state_start)
        occ_last = min(nocc_spin, dg_frag%coef_state_end)
        if (occ_first > occ_last) cycle
      end if

      row_needed(:) = .false.
      row_output(:) = .false.
      do iblk = 1, dg_frag%n_momentum_blocks
        if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
        ifrag_row = dg_frag%momentum_blocks(iblk)%ifrag_row
        ifrag_col = dg_frag%momentum_blocks(iblk)%ifrag_col
        if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
        if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
        nrow = min(dg_frag%n_basis(ifrag_row, ispin), size(dg_frag%momentum_blocks(iblk)%val, 2), &
                   size(dg_frag%index_basis, 1))
        ncol = min(dg_frag%n_basis(ifrag_col, ispin), size(dg_frag%momentum_blocks(iblk)%val, 3), &
                   size(dg_frag%index_basis, 1))
        if (nrow <= 0 .or. ncol <= 0) cycle
        have_owned_rows = .false.
        do io = 1, nrow
          gid_i = dg_frag%index_basis(io, ifrag_row, ispin)
          if (gid_i < 1 .or. gid_i > n_global) cycle
          if (dg_frag%coef_state_block_mode) then
            if (dg_frag%coef_global_to_local(gid_i, ispin) <= 0) cycle
          else
            if (dg_frag%coef_owner(gid_i, ispin) /= dg_frag%id) cycle
          end if
          have_owned_rows = .true.
          row_output(gid_i) = .true.
          row_needed(gid_i) = .true.
        end do
        if (.not. have_owned_rows) cycle
        do jo = 1, ncol
          gid_j = dg_frag%index_basis(jo, ifrag_col, ispin)
          if (gid_j < 1 .or. gid_j > n_global) cycle
          row_needed(gid_j) = .true.
        end do
      end do

      n_needed = count(row_needed)
      if (n_needed <= 0) cycle
      if (.not. allocated(needed_ids)) then
        allocate(needed_ids(n_needed))
      else if (size(needed_ids, 1) /= n_needed) then
        deallocate(needed_ids)
        allocate(needed_ids(n_needed))
      end if
      needed_pos(:) = 0
      n_needed = 0
      do gid_i = 1, n_global
        if (.not. row_needed(gid_i)) cycle
        n_needed = n_needed + 1
        needed_ids(n_needed) = gid_i
        needed_pos(gid_i) = n_needed
      end do

      n_output = count(row_output)
      if (n_output <= 0) cycle
      if (.not. allocated(output_pos)) then
        allocate(output_pos(n_output))
      else if (size(output_pos, 1) /= n_output) then
        deallocate(output_pos)
        allocate(output_pos(n_output))
      end if
      n_output = 0
      do gid_i = 1, n_global
        if (.not. row_output(gid_i)) cycle
        pos_i = needed_pos(gid_i)
        if (pos_i <= 0) cycle
        n_output = n_output + 1
        output_pos(n_output) = pos_i
      end do
      if (n_output <= 0) cycle

      max_basis = size(dg_frag%index_basis, 1)
      if (dg_frag%n_momentum_blocks <= 0) cycle
      if (allocated(block_ids)) then
        if (size(block_ids, 1) /= dg_frag%n_momentum_blocks .or. &
            size(block_row_lid, 1) /= max_basis .or. &
            size(block_row_lid, 2) /= dg_frag%n_momentum_blocks) then
          deallocate(block_ids, block_nrow, block_ncol)
          deallocate(block_row_lid, block_row_pos, block_col_lid, block_col_pos)
        end if
      end if
      if (.not. allocated(block_ids)) then
        allocate(block_ids(dg_frag%n_momentum_blocks), block_nrow(dg_frag%n_momentum_blocks), &
                 block_ncol(dg_frag%n_momentum_blocks))
        allocate(block_row_lid(max_basis, dg_frag%n_momentum_blocks), &
                 block_row_pos(max_basis, dg_frag%n_momentum_blocks))
        allocate(block_col_lid(max_basis, dg_frag%n_momentum_blocks), &
                 block_col_pos(max_basis, dg_frag%n_momentum_blocks))
      end if
      block_nrow(:) = 0
      block_ncol(:) = 0
      n_active_blocks = 0
      do iblk = 1, dg_frag%n_momentum_blocks
        if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
        ifrag_row = dg_frag%momentum_blocks(iblk)%ifrag_row
        ifrag_col = dg_frag%momentum_blocks(iblk)%ifrag_col
        if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
        if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
        nrow = min(dg_frag%n_basis(ifrag_row, ispin), size(dg_frag%momentum_blocks(iblk)%val, 2), max_basis)
        ncol = min(dg_frag%n_basis(ifrag_col, ispin), size(dg_frag%momentum_blocks(iblk)%val, 3), max_basis)
        if (nrow <= 0 .or. ncol <= 0) cycle

        nrow_active = 0
        do io = 1, nrow
          gid_i = dg_frag%index_basis(io, ifrag_row, ispin)
          if (gid_i < 1 .or. gid_i > n_global) cycle
          if (.not. row_output(gid_i)) cycle
          pos_i = needed_pos(gid_i)
          if (pos_i <= 0) cycle
          nrow_active = nrow_active + 1
          block_row_lid(nrow_active, iblk) = io
          block_row_pos(nrow_active, iblk) = pos_i
        end do
        if (nrow_active <= 0) cycle

        ncol_active = 0
        do jo = 1, ncol
          gid_j = dg_frag%index_basis(jo, ifrag_col, ispin)
          if (gid_j < 1 .or. gid_j > n_global) cycle
          pos_j = needed_pos(gid_j)
          if (pos_j <= 0) cycle
          ncol_active = ncol_active + 1
          block_col_lid(ncol_active, iblk) = jo
          block_col_pos(ncol_active, iblk) = pos_j
        end do
        if (ncol_active <= 0) cycle

        n_active_blocks = n_active_blocks + 1
        block_ids(n_active_blocks) = iblk
        block_nrow(iblk) = nrow_active
        block_ncol(iblk) = ncol_active
      end do
      if (n_active_blocks <= 0) cycle

      target_bytes = int(current_coef_cache_target_mb, kind=8) * 1024_8 * 1024_8
      bytes_per_state = 16_8 * int(max(1, n_needed), kind=8) * 4_8
      state_work = max(1, min(occ_last - occ_first + 1, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))
      if (.not. allocated(coef_work)) then
        allocate(coef_work(n_needed, state_work), pc_work(n_needed, state_work, 3))
      else if (size(coef_work, 1) /= n_needed .or. size(coef_work, 2) /= state_work) then
        deallocate(coef_work, pc_work)
        allocate(coef_work(n_needed, state_work), pc_work(n_needed, state_work, 3))
      end if
      if (.not. allocated(coef_col_work)) then
        allocate(coef_col_work(max_basis, state_work), pc_block_work(max_basis, state_work), &
                 mom_block_work(max_basis, max_basis))
      else if (size(coef_col_work, 1) /= max_basis .or. size(coef_col_work, 2) /= state_work) then
        deallocate(coef_col_work, pc_block_work, mom_block_work)
        allocate(coef_col_work(max_basis, state_work), pc_block_work(max_basis, state_work), &
                 mom_block_work(max_basis, max_basis))
      end if
      if (trace_decomp) then
        if (.not. allocated(pc_self_work)) then
          allocate(pc_self_work(n_needed, state_work, 3), pc_cross_work(n_needed, state_work, 3))
        else if (size(pc_self_work, 1) /= n_needed .or. size(pc_self_work, 2) /= state_work) then
          deallocate(pc_self_work, pc_cross_work)
          allocate(pc_self_work(n_needed, state_work, 3), pc_cross_work(n_needed, state_work, 3))
        end if
      end if

      do occ0 = occ_first, occ_last, state_work
        nbatch = min(state_work, occ_last - occ0 + 1)
        coef_work(:, :) = (0.0d0, 0.0d0)
        pc_work(:, :, :) = (0.0d0, 0.0d0)
        if (trace_decomp) then
          pc_self_work(:, :, :) = (0.0d0, 0.0d0)
          pc_cross_work(:, :, :) = (0.0d0, 0.0d0)
        end if
        call fetch_remote_coef_rows(dg_frag, ispin, needed_ids, coef_work(1:n_needed, 1:nbatch), &
                                    occ0, occ0 + nbatch - 1)

        if (trace_decomp) then
          do iactive = 1, n_active_blocks
            iblk = block_ids(iactive)
            ifrag_row = dg_frag%momentum_blocks(iblk)%ifrag_row
            ifrag_col = dg_frag%momentum_blocks(iblk)%ifrag_col
            nrow = block_nrow(iblk)
            ncol = block_ncol(iblk)
            do io = 1, nrow
              io_block = block_row_lid(io, iblk)
              pos_i = block_row_pos(io, iblk)
              do jo = 1, ncol
                jo_block = block_col_lid(jo, iblk)
                pos_j = block_col_pos(jo, iblk)
                do idir = 1, 3
                  block_val = dg_frag%momentum_blocks(iblk)%val(idir, io_block, jo_block, ispin)
                  if (block_val == 0.0d0) cycle
!$omp simd
                  do idx_row = 1, nbatch
                    pc_work(pos_i, idx_row, idir) = pc_work(pos_i, idx_row, idir) + &
                      block_val * coef_work(pos_j, idx_row)
                    if (ifrag_row == ifrag_col) then
                      pc_self_work(pos_i, idx_row, idir) = pc_self_work(pos_i, idx_row, idir) + &
                        block_val * coef_work(pos_j, idx_row)
                    else
                      pc_cross_work(pos_i, idx_row, idir) = pc_cross_work(pos_i, idx_row, idir) + &
                        block_val * coef_work(pos_j, idx_row)
                    end if
                  end do
                end do
              end do
            end do
          end do
        else
          do iactive = 1, n_active_blocks
            iblk = block_ids(iactive)
            nrow = block_nrow(iblk)
            ncol = block_ncol(iblk)
            coef_col_work(1:ncol,1:nbatch) = (0.0d0, 0.0d0)
            do jo = 1, ncol
              pos_j = block_col_pos(jo, iblk)
              coef_col_work(jo,1:nbatch) = coef_work(pos_j,1:nbatch)
            end do
            do idir = 1, 3
              do jo = 1, ncol
                jo_block = block_col_lid(jo, iblk)
                do io = 1, nrow
                  io_block = block_row_lid(io, iblk)
                  mom_block_work(io, jo) = cmplx(dg_frag%momentum_blocks(iblk)%val(idir, io_block, jo_block, ispin), &
                    0.0d0, kind=8)
                end do
              end do
              call zgemm('N', 'N', nrow, nbatch, ncol, (1.0d0, 0.0d0), &
                mom_block_work, max(1,max_basis), coef_col_work, max(1,max_basis), &
                (0.0d0, 0.0d0), pc_block_work, max(1,max_basis))
              if (dg_frag%mixed_z_perf_count_enabled) then
                dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 1_8
              end if
              do io = 1, nrow
                pos_i = block_row_pos(io, iblk)
                pc_work(pos_i,1:nbatch,idir) = pc_work(pos_i,1:nbatch,idir) + pc_block_work(io,1:nbatch)
              end do
            end do
          end do
        end if

        curr_x = 0.0d0
        curr_y = 0.0d0
        curr_z = 0.0d0
        curr_self_x = 0.0d0
        curr_self_y = 0.0d0
        curr_self_z = 0.0d0
        curr_cross_x = 0.0d0
        curr_cross_y = 0.0d0
        curr_cross_z = 0.0d0
!$omp parallel do private(iout, pos_i, idx_row, occ_factor, pair_x, pair_y, pair_z) &
!$omp& private(pair_self_x,pair_self_y,pair_self_z,pair_cross_x,pair_cross_y,pair_cross_z) &
!$omp& reduction(+:curr_x,curr_y,curr_z,curr_self_x,curr_self_y,curr_self_z,curr_cross_x,curr_cross_y,curr_cross_z) &
!$omp& schedule(static)
        do iout = 1, n_output
          pos_i = output_pos(iout)
          pair_x = 0.0d0
          pair_y = 0.0d0
          pair_z = 0.0d0
          pair_self_x = 0.0d0
          pair_self_y = 0.0d0
          pair_self_z = 0.0d0
          pair_cross_x = 0.0d0
          pair_cross_y = 0.0d0
          pair_cross_z = 0.0d0
!$omp simd reduction(+:pair_x,pair_y,pair_z,pair_self_x,pair_self_y,pair_self_z,pair_cross_x,pair_cross_y,pair_cross_z)
          do idx_row = 1, nbatch
            occ_factor = system%rocc(occ0 + idx_row - 1, 1, ispin)
            pair_x = pair_x + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_work(pos_i, idx_row, 1))
            pair_y = pair_y + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_work(pos_i, idx_row, 2))
            pair_z = pair_z + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_work(pos_i, idx_row, 3))
            if (trace_decomp) then
              pair_self_x = pair_self_x + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_self_work(pos_i, idx_row, 1))
              pair_self_y = pair_self_y + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_self_work(pos_i, idx_row, 2))
              pair_self_z = pair_self_z + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_self_work(pos_i, idx_row, 3))
              pair_cross_x = pair_cross_x + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_cross_work(pos_i, idx_row, 1))
              pair_cross_y = pair_cross_y + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_cross_work(pos_i, idx_row, 2))
              pair_cross_z = pair_cross_z + occ_factor * aimag(conjg(coef_work(pos_i, idx_row)) * pc_cross_work(pos_i, idx_row, 3))
            end if
          end do
          curr_x = curr_x + pair_x
          curr_y = curr_y + pair_y
          curr_z = curr_z + pair_z
          curr_self_x = curr_self_x + pair_self_x
          curr_self_y = curr_self_y + pair_self_y
          curr_self_z = curr_self_z + pair_self_z
          curr_cross_x = curr_cross_x + pair_cross_x
          curr_cross_y = curr_cross_y + pair_cross_y
          curr_cross_z = curr_cross_z + pair_cross_z
        end do
!$omp end parallel do
        if (curr_x /= curr_x .or. curr_y /= curr_y .or. curr_z /= curr_z) then
          stop "DG-Fragment RT: NaN in block current density matrix"
        end if
        curr_local(1) = curr_local(1) + curr_x
        curr_local(2) = curr_local(2) + curr_y
        curr_local(3) = curr_local(3) + curr_z
        if (trace_decomp) then
          curr_self_local(:) = curr_self_local(:) + [curr_self_x, curr_self_y, curr_self_z]
          curr_cross_local(:) = curr_cross_local(:) + [curr_cross_x, curr_cross_y, curr_cross_z]
        end if
      end do
    end do

    call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
    if (trace_decomp) then
      call comm_summation(curr_self_local, curr_self_sum, 3, dg_frag%icomm)
      call comm_summation(curr_cross_local, curr_cross_sum, 3, dg_frag%icomm)
      dg_frag%current_momentum_self_raw(:) = curr_self_sum(:)
      dg_frag%current_momentum_cross_raw(:) = curr_cross_sum(:)
      dg_frag%current_momentum_decomp_ready = .true.
    end if
    if (curr_sum(1) /= curr_sum(1) .or. curr_sum(2) /= curr_sum(2) .or. curr_sum(3) /= curr_sum(3)) then
      stop "DG-Fragment RT: NaN in block macroscopic current reduction"
    end if
    current_raw(:) = curr_sum(:)
  end subroutine calculate_macroscopic_current_from_momentum_blocks


  subroutine calculate_local_wannier_polarization_dg(dg_frag, system, polarization_raw)
    use structures, only: s_dft_system
    use communication, only: comm_summation
    use salmon_global, only: yn_dg_mixed_z
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(out)   :: polarization_raw(3)

    integer, parameter :: pol_coef_cache_target_mb = 64
    integer :: ispin, ifrag, i_local, ib, iw, jw, idir
    integer :: iblk_idx, iblk, ifrag_row, ifrag_col, io, jo, nrow, ncol
    integer :: row_gid, col_gid, row_pos, col_pos
    integer :: istate, occ0, nbatch, state_work
    integer :: nocc_spin, occ_first, occ_last
    integer :: n_global, n_needed, gid, pos, nkeep, nkeep_max, nbasis
    integer :: local_idx, local_col
    integer, allocatable :: needed_ids(:), needed_pos(:)
    logical, allocatable :: row_needed(:)
    complex(8), allocatable :: coef_work(:,:), cw(:,:)
    real(8) :: pol_local(3), pol_sum(3), occ, contrib
    integer(8) :: target_bytes, bytes_per_state
    logical :: use_buffer_wannier

    polarization_raw(:) = 0.0d0
    pol_local(:) = 0.0d0
    pol_sum(:) = 0.0d0
    use_buffer_wannier = dg_frag%buffer_wannier_flux_seed_applied .and. &
      dg_frag%has_buffer_periodic_wannier_basis .and. &
      allocated(dg_frag%buffer_wannier_coef) .and. allocated(dg_frag%buffer_wannier_v)
    if (yn_dg_mixed_z == 'y' .and. dg_frag%has_mixed_wannier_bpw_position .and. &
        allocated(dg_frag%mixed_wannier_bpw_z) .and. allocated(dg_frag%mixed_wannier_bpw_pcoef)) then
      call calculate_global_mixed_wannier_bpw_polarization_dg(dg_frag, system, polarization_raw)
      return
    end if
    if (.not. use_buffer_wannier) then
      if (dg_frag%has_global_wannier_basis .and. allocated(dg_frag%global_wannier_coef) .and. &
          allocated(dg_frag%global_wannier_center)) then
        call calculate_global_wannier_center_polarization_dg(dg_frag, system, polarization_raw)
        return
      end if
      if (.not. dg_frag%has_local_wannier_basis) return
      if (.not. allocated(dg_frag%local_wannier_coef)) return
      if (.not. allocated(dg_frag%local_wannier_r)) return
    end if
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. allocated(dg_frag%n_basis)) return

    if (.not. dg_frag%coef_state_block_mode .and. .not. dg_frag%is_frag_root) then
      call comm_summation(pol_local, pol_sum, 3, dg_frag%icomm)
      polarization_raw(:) = pol_sum(:)
      return
    end if

    n_global = size(dg_frag%coef_owner, 1)
    allocate(row_needed(n_global), needed_pos(n_global))
    do ispin = 1, min(dg_frag%nspin, system%nspin)
      nocc_spin = 0
      if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
      if (dg_frag%coef_state_block_mode) then
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(system%rocc, 1))
      else
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(dg_frag%coef, 2), size(system%rocc, 1))
      end if
      if (nocc_spin <= 0) cycle
      occ_first = 1
      occ_last = nocc_spin
      if (dg_frag%coef_state_block_mode) then
        occ_first = max(1, dg_frag%coef_state_start)
        occ_last = min(nocc_spin, dg_frag%coef_state_end)
        if (occ_first > occ_last) cycle
      end if

      row_needed(:) = .false.
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = ifrag - dg_frag%ifrag_start + 1
        if (use_buffer_wannier) then
          if (i_local < 1 .or. i_local > size(dg_frag%buffer_wannier_nkeep)) cycle
          nkeep = dg_frag%buffer_wannier_nkeep(i_local)
        else
          if (i_local < 1 .or. i_local > size(dg_frag%local_wannier_nkeep)) cycle
          nkeep = dg_frag%local_wannier_nkeep(i_local)
        end if
        if (nkeep <= 0) cycle
        if (use_buffer_wannier) then
          nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                       size(dg_frag%buffer_wannier_coef, 1))
        else
          nbasis = min(dg_frag%local_wannier_nbasis(i_local), dg_frag%n_basis(ifrag, ispin), &
                       size(dg_frag%index_basis, 1))
        end if
        do ib = 1, nbasis
          gid = dg_frag%index_basis(ib, ifrag, ispin)
          if (gid >= 1 .and. gid <= n_global) row_needed(gid) = .true.
        end do
      end do
      if (use_buffer_wannier .and. dg_frag%buffer_wannier_xi_flux_available .and. &
          allocated(dg_frag%buffer_wannier_xi_flux_blocks) .and. &
          allocated(dg_frag%buffer_wannier_xi_flux_local_block_ids)) then
        do iblk_idx = 1, size(dg_frag%buffer_wannier_xi_flux_local_block_ids)
          iblk = dg_frag%buffer_wannier_xi_flux_local_block_ids(iblk_idx)
          if (iblk < 1 .or. iblk > size(dg_frag%buffer_wannier_xi_flux_blocks)) cycle
          ifrag_row = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_row
          ifrag_col = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_col
          if (ifrag_row < dg_frag%ifrag_start .or. ifrag_row > dg_frag%ifrag_end) cycle
          if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
          if (.not. allocated(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val)) cycle
          nrow = min(dg_frag%n_basis(ifrag_row, ispin), &
                     size(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val, 2), &
                     size(dg_frag%index_basis, 1))
          ncol = min(dg_frag%n_basis(ifrag_col, ispin), &
                     size(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val, 3), &
                     size(dg_frag%index_basis, 1))
          do io = 1, nrow
            row_gid = dg_frag%index_basis(io, ifrag_row, ispin)
            if (row_gid >= 1 .and. row_gid <= n_global) row_needed(row_gid) = .true.
          end do
          do jo = 1, ncol
            col_gid = dg_frag%index_basis(jo, ifrag_col, ispin)
            if (col_gid >= 1 .and. col_gid <= n_global) row_needed(col_gid) = .true.
          end do
        end do
      end if
      n_needed = count(row_needed)
      if (n_needed <= 0) cycle
      if (allocated(needed_ids)) deallocate(needed_ids)
      allocate(needed_ids(n_needed))
      needed_pos(:) = 0
      n_needed = 0
      do gid = 1, n_global
        if (.not. row_needed(gid)) cycle
        n_needed = n_needed + 1
        needed_ids(n_needed) = gid
        needed_pos(gid) = n_needed
      end do

      target_bytes = int(pol_coef_cache_target_mb, kind=8) * 1024_8 * 1024_8
      bytes_per_state = 16_8 * int(max(1, n_needed), kind=8) * 4_8
      state_work = max(1, min(occ_last - occ_first + 1, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))
      if (use_buffer_wannier) then
        nkeep_max = max(1, maxval(dg_frag%buffer_wannier_nkeep))
      else
        nkeep_max = max(1, maxval(dg_frag%local_wannier_nkeep))
      end if
      if (allocated(coef_work)) deallocate(coef_work)
      if (allocated(cw)) deallocate(cw)
      allocate(coef_work(n_needed, state_work))
      allocate(cw(nkeep_max, state_work))

      do occ0 = occ_first, occ_last, state_work
        nbatch = min(state_work, occ_last - occ0 + 1)
        coef_work(:, :) = (0.0d0, 0.0d0)
        if (dg_frag%coef_state_block_mode) then
          do pos = 1, n_needed
            local_idx = dg_frag%coef_global_to_local(needed_ids(pos), ispin)
            if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
            do istate = 1, nbatch
              local_col = occ0 + istate - 1 - dg_frag%coef_state_start + 1
              if (local_col < 1 .or. local_col > size(dg_frag%coef, 2)) cycle
              coef_work(pos, istate) = dg_frag%coef(local_idx, local_col, ispin)
            end do
          end do
        else
          call fetch_remote_coef_rows(dg_frag, ispin, needed_ids, coef_work(1:n_needed, 1:nbatch), &
                                      occ0, occ0 + nbatch - 1)
        end if

        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = ifrag - dg_frag%ifrag_start + 1
          if (use_buffer_wannier) then
            if (i_local < 1 .or. i_local > size(dg_frag%buffer_wannier_nkeep)) cycle
            nkeep = dg_frag%buffer_wannier_nkeep(i_local)
          else
            if (i_local < 1 .or. i_local > size(dg_frag%local_wannier_nkeep)) cycle
            nkeep = dg_frag%local_wannier_nkeep(i_local)
          end if
          if (nkeep <= 0) cycle
          if (use_buffer_wannier) then
            nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                         size(dg_frag%buffer_wannier_coef, 1))
          else
            nbasis = min(dg_frag%local_wannier_nbasis(i_local), dg_frag%n_basis(ifrag, ispin), &
                         size(dg_frag%index_basis, 1))
          end if
          cw(1:nkeep, 1:nbatch) = (0.0d0, 0.0d0)
          do iw = 1, nkeep
            do ib = 1, nbasis
              gid = dg_frag%index_basis(ib, ifrag, ispin)
              if (gid < 1 .or. gid > n_global) cycle
              pos = needed_pos(gid)
              if (pos <= 0) cycle
!$omp simd
              do istate = 1, nbatch
                if (use_buffer_wannier) then
                  cw(iw, istate) = cw(iw, istate) + &
                    dg_frag%buffer_wannier_coef(ib, iw, ispin, i_local) * coef_work(pos, istate)
                else
                  cw(iw, istate) = cw(iw, istate) + &
                    dg_frag%local_wannier_coef(ib, iw, ispin, i_local) * coef_work(pos, istate)
                end if
              end do
            end do
          end do
          do idir = 1, 3
            contrib = 0.0d0
            do istate = 1, nbatch
              occ = max(0.0d0, system%rocc(occ0 + istate - 1, 1, ispin))
              if (occ <= 0.0d0) cycle
              do iw = 1, nkeep
                do jw = 1, nkeep
                  if (use_buffer_wannier) then
                    if (iw == jw .and. allocated(dg_frag%buffer_wannier_frag_center)) then
                      contrib = contrib - occ * real(conjg(cw(iw, istate)) * &
                        (dg_frag%buffer_wannier_v(idir, iw, jw, i_local) + &
                         dg_frag%buffer_wannier_frag_center(idir, i_local)) * cw(jw, istate), 8)
                    else
                      contrib = contrib - occ * real(conjg(cw(iw, istate)) * &
                        dg_frag%buffer_wannier_v(idir, iw, jw, i_local) * cw(jw, istate), 8)
                    end if
                  else
                    contrib = contrib - occ * real(conjg(cw(iw, istate)) * &
                      dg_frag%local_wannier_r(idir, iw, jw, ispin, i_local) * cw(jw, istate), 8)
                  end if
                end do
              end do
            end do
            pol_local(idir) = pol_local(idir) + contrib
          end do
        end do
        if (use_buffer_wannier .and. dg_frag%buffer_wannier_xi_flux_available .and. &
            allocated(dg_frag%buffer_wannier_xi_flux_blocks) .and. &
            allocated(dg_frag%buffer_wannier_xi_flux_local_block_ids)) then
          do iblk_idx = 1, size(dg_frag%buffer_wannier_xi_flux_local_block_ids)
            iblk = dg_frag%buffer_wannier_xi_flux_local_block_ids(iblk_idx)
            if (iblk < 1 .or. iblk > size(dg_frag%buffer_wannier_xi_flux_blocks)) cycle
            ifrag_row = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_row
            ifrag_col = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_col
            if (ifrag_row < dg_frag%ifrag_start .or. ifrag_row > dg_frag%ifrag_end) cycle
            if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
            if (.not. allocated(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val)) cycle
            nrow = min(dg_frag%n_basis(ifrag_row, ispin), &
                       size(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val, 2), &
                       size(dg_frag%index_basis, 1))
            ncol = min(dg_frag%n_basis(ifrag_col, ispin), &
                       size(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val, 3), &
                       size(dg_frag%index_basis, 1))
            do io = 1, nrow
              row_gid = dg_frag%index_basis(io, ifrag_row, ispin)
              if (row_gid < 1 .or. row_gid > n_global) cycle
              row_pos = needed_pos(row_gid)
              if (row_pos <= 0) cycle
              do jo = 1, ncol
                col_gid = dg_frag%index_basis(jo, ifrag_col, ispin)
                if (col_gid < 1 .or. col_gid > n_global) cycle
                col_pos = needed_pos(col_gid)
                if (col_pos <= 0) cycle
                do idir = 1, 3
                  contrib = 0.0d0
                  do istate = 1, nbatch
                    occ = max(0.0d0, system%rocc(occ0 + istate - 1, 1, ispin))
                    if (occ <= 0.0d0) cycle
                    contrib = contrib - occ * real(conjg(coef_work(row_pos, istate)) * &
                      dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val(idir, io, jo, ispin) * &
                      coef_work(col_pos, istate), 8)
                  end do
                  pol_local(idir) = pol_local(idir) + contrib
                end do
              end do
            end do
          end do
        end if
      end do
    end do

    call comm_summation(pol_local, pol_sum, 3, dg_frag%icomm)
    if (pol_sum(1) /= pol_sum(1) .or. pol_sum(2) /= pol_sum(2) .or. pol_sum(3) /= pol_sum(3)) then
      stop "DG-Fragment RT: NaN in local Wannier polarization reduction"
    end if
    polarization_raw(:) = pol_sum(:)
  end subroutine calculate_local_wannier_polarization_dg

  subroutine calculate_global_wannier_center_polarization_dg(dg_frag, system, polarization_raw)
    use structures, only: s_dft_system
    use communication, only: comm_summation
    use salmon_global, only: yn_dg_mixed_z
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(out)   :: polarization_raw(3)

    integer, parameter :: pol_coef_cache_target_mb = 64
    integer :: ispin, ifrag, i_local, ib, iw, jw, idir
    integer :: istate, occ0, nbatch, state_work
    integer :: nocc_spin, occ_first, occ_last
    integer :: n_global, n_needed, gid, pos, nbasis
    integer :: n_wann, local_idx, local_col
    integer, allocatable :: needed_ids(:), needed_pos(:)
    logical, allocatable :: row_needed(:)
    complex(8), allocatable :: coef_work(:,:), cw_local(:,:), cw_sum(:,:)
    real(8) :: pol_local(3), pol_sum(3), occ, weight
    complex(8) :: pos_expect
    integer(8) :: target_bytes, bytes_per_state
    integer :: env_status, env_len
    logical :: accumulate_on_rank, force_full_position, use_mixed_z
    character(32) :: env_full_position, env_mixed_z

    polarization_raw(:) = 0.0d0
    pol_local(:) = 0.0d0
    pol_sum(:) = 0.0d0
    force_full_position = dg_frag%has_global_wannier_position .and. allocated(dg_frag%global_wannier_position)
    env_full_position = ''
    call get_environment_variable('SALMON_DG_POL_GLOBAL_WANNIER_FULL', env_full_position, &
      length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_full_position(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        force_full_position = dg_frag%has_global_wannier_position .and. allocated(dg_frag%global_wannier_position)
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        force_full_position = .false.
      end select
    end if
    env_mixed_z = ''
    call get_environment_variable('SALMON_DG_MIXED_Z', env_mixed_z, length=env_len, status=env_status)
    use_mixed_z = (yn_dg_mixed_z == 'y')
    if (env_status == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_mixed_z(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        use_mixed_z = dg_frag%has_mixed_wannier_bpw_position .and. allocated(dg_frag%mixed_wannier_bpw_z) .and. &
          allocated(dg_frag%mixed_wannier_bpw_pcoef)
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        use_mixed_z = .false.
      end select
    end if
    if (use_mixed_z) use_mixed_z = dg_frag%has_mixed_wannier_bpw_position .and. &
      allocated(dg_frag%mixed_wannier_bpw_z) .and. allocated(dg_frag%mixed_wannier_bpw_pcoef)
    if (use_mixed_z) then
      call calculate_global_mixed_wannier_bpw_polarization_dg(dg_frag, system, polarization_raw)
      return
    end if

    if (.not. force_full_position .and. dg_frag%has_formal_dg_wannier_basis .and. allocated(dg_frag%dg_wannier_basis_coef) .and. &
        allocated(dg_frag%dg_wannier_xi_local) .and. allocated(dg_frag%dg_wannier_ref_center) .and. &
        allocated(dg_frag%dg_wannier_nkeep)) then
      call calculate_formal_dg_wannier_polarization_dg(dg_frag, system, polarization_raw)
      return
    end if

    if (.not. dg_frag%has_global_wannier_basis) return
    if (.not. allocated(dg_frag%global_wannier_coef)) return
    if (.not. allocated(dg_frag%global_wannier_center)) return
    if (.not. force_full_position .and. dg_frag%has_global_wannier_local_basis .and. allocated(dg_frag%global_wannier_local_coef) .and. &
        allocated(dg_frag%global_wannier_local_center) .and. allocated(dg_frag%global_wannier_local_nkeep)) then
      call calculate_global_wannier_local_center_polarization_dg(dg_frag, system, polarization_raw)
      return
    end if
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. allocated(dg_frag%n_basis)) return

    n_wann = dg_frag%global_wannier_num_wann
    if (n_wann <= 0) return
    n_global = size(dg_frag%coef_owner, 1)
    allocate(row_needed(n_global), needed_pos(n_global))

    do ispin = 1, min(dg_frag%nspin, system%nspin)
      nocc_spin = 0
      if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
      if (dg_frag%coef_state_block_mode) then
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(system%rocc, 1))
      else
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(dg_frag%coef, 2), size(system%rocc, 1))
      end if
      if (nocc_spin <= 0) cycle

      occ_first = 1
      occ_last = nocc_spin
      if (dg_frag%coef_state_block_mode) then
        occ_first = max(1, dg_frag%coef_state_start)
        occ_last = min(nocc_spin, dg_frag%coef_state_end)
        if (occ_first > occ_last) cycle
      end if

      row_needed(:) = .false.
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = ifrag - dg_frag%ifrag_start + 1
        if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
        nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                     size(dg_frag%global_wannier_coef, 1))
        do ib = 1, nbasis
          gid = dg_frag%index_basis(ib, ifrag, ispin)
          if (gid >= 1 .and. gid <= n_global) row_needed(gid) = .true.
        end do
      end do

      n_needed = count(row_needed)
      if (n_needed <= 0) cycle
      if (allocated(needed_ids)) deallocate(needed_ids)
      allocate(needed_ids(n_needed))
      needed_pos(:) = 0
      n_needed = 0
      do gid = 1, n_global
        if (.not. row_needed(gid)) cycle
        n_needed = n_needed + 1
        needed_ids(n_needed) = gid
        needed_pos(gid) = n_needed
      end do

      target_bytes = int(pol_coef_cache_target_mb, kind=8) * 1024_8 * 1024_8
      bytes_per_state = 16_8 * int(max(1, n_needed + 2 * n_wann), kind=8)
      state_work = max(1, min(occ_last - occ_first + 1, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))

      if (allocated(coef_work)) deallocate(coef_work)
      if (allocated(cw_local)) deallocate(cw_local)
      if (allocated(cw_sum)) deallocate(cw_sum)
      allocate(coef_work(n_needed, state_work))
      allocate(cw_local(n_wann, state_work))
      allocate(cw_sum(n_wann, state_work))

      do occ0 = occ_first, occ_last, state_work
        nbatch = min(state_work, occ_last - occ0 + 1)
        coef_work(:, :) = (0.0d0, 0.0d0)
        if (dg_frag%coef_state_block_mode) then
          do pos = 1, n_needed
            local_idx = dg_frag%coef_global_to_local(needed_ids(pos), ispin)
            if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
            do istate = 1, nbatch
              local_col = occ0 + istate - 1 - dg_frag%coef_state_start + 1
              if (local_col < 1 .or. local_col > size(dg_frag%coef, 2)) cycle
              coef_work(pos, istate) = dg_frag%coef(local_idx, local_col, ispin)
            end do
          end do
        else
          call fetch_remote_coef_rows(dg_frag, ispin, needed_ids, coef_work(1:n_needed, 1:nbatch), &
                                      occ0, occ0 + nbatch - 1)
        end if

        cw_local(:, :) = (0.0d0, 0.0d0)
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = ifrag - dg_frag%ifrag_start + 1
          if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
          nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                       size(dg_frag%global_wannier_coef, 1))
          do ib = 1, nbasis
            gid = dg_frag%index_basis(ib, ifrag, ispin)
            if (gid < 1 .or. gid > n_global) cycle
            pos = needed_pos(gid)
            if (pos <= 0) cycle
            do iw = 1, n_wann
!$omp simd
              do istate = 1, nbatch
                cw_local(iw, istate) = cw_local(iw, istate) + &
                  conjg(dg_frag%global_wannier_coef(ib, iw, ispin, i_local)) * coef_work(pos, istate)
              end do
            end do
          end do
        end do

        if (dg_frag%isize > 1) then
          call comm_summation(cw_local, cw_sum, n_wann * nbatch, dg_frag%icomm)
        else
          cw_sum(:, 1:nbatch) = cw_local(:, 1:nbatch)
        end if

        accumulate_on_rank = dg_frag%coef_state_block_mode .or. (dg_frag%id == 0)
        if (accumulate_on_rank) then
          do istate = 1, nbatch
            occ = max(0.0d0, system%rocc(occ0 + istate - 1, 1, ispin))
            if (occ <= 0.0d0) cycle
            if (force_full_position) then
              do idir = 1, 3
                pos_expect = (0.0d0, 0.0d0)
                do iw = 1, n_wann
                  do jw = 1, n_wann
                    pos_expect = pos_expect + conjg(cw_sum(iw, istate)) * &
                      dg_frag%global_wannier_position(idir, iw, jw) * cw_sum(jw, istate)
                  end do
                end do
                pol_local(idir) = pol_local(idir) - occ * real(pos_expect, kind=8)
              end do
            else
              do iw = 1, n_wann
                weight = occ * real(conjg(cw_sum(iw, istate)) * cw_sum(iw, istate), 8)
                if (weight <= 0.0d0) cycle
                do idir = 1, 3
                  pol_local(idir) = pol_local(idir) - weight * dg_frag%global_wannier_center(idir, iw)
                end do
              end do
            end if
          end do
        end if
      end do
    end do

    call comm_summation(pol_local, pol_sum, 3, dg_frag%icomm)
    if (pol_sum(1) /= pol_sum(1) .or. pol_sum(2) /= pol_sum(2) .or. pol_sum(3) /= pol_sum(3)) then
      stop "DG-Fragment RT: NaN in global Wannier polarization reduction"
    end if
    polarization_raw(:) = pol_sum(:)
  end subroutine calculate_global_wannier_center_polarization_dg

  subroutine calculate_global_mixed_wannier_bpw_polarization_dg(dg_frag, system, polarization_raw)
    use structures, only: s_dft_system
    use communication, only: comm_summation
    use misc_routines, only: get_wtime
    use salmon_global, only: dg_mixed_z_polarization_branch
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(out)   :: polarization_raw(3)

    integer :: ispin, ifrag, i_local, ib, iw, jw, idir
    integer :: istate, nocc_spin, n_w, n_p, n_mix, nbasis
    integer :: global_row, local_row, n_wann_raw
    real(8) :: pol_local(3), pol_sum(3), occ
    real(8) :: pol_ww_local(3), pol_wp_local(3), pol_ww_sum(3), pol_wp_sum(3)
    real(8) :: pol_ww_diag_local(3), pol_ww_offdiag_local(3)
    real(8) :: pol_ww_same_owner_local(3), pol_ww_cross_owner_local(3)
    real(8) :: pol_ww_diag_sum(3), pol_ww_offdiag_sum(3)
    real(8) :: pol_ww_same_owner_sum(3), pol_ww_cross_owner_sum(3)
    real(8) :: ww_pair_count_total, ww_pair_count_diag, ww_pair_count_offdiag, ww_pair_count_cross_owner
    real(8) :: occ_local, occ_sum, w_occ_local, w_occ_sum
    real(8) :: zww_diag_sum, center_diag_sum, center_local_sum, zww_diag_val, center_val
    real(8) :: center_eig_weighted_sum
    real(8) :: zww_diag_min, zww_diag_max, center_min, center_max, diff_val
    real(8) :: diff_min, diff_max, diff_rms, rho_diag, weight_diag, weighted_zww, weighted_center
    real(8) :: weighted_diff, weighted_zww_sum, weighted_center_sum, weighted_diff_sum
    real(8) :: center_diag_weighted_local(3), center_diag_weighted_sum(3)
    real(8) :: center_eig_weighted_local(3), center_eig_weighted_sum_global(3)
    real(8) :: local_center_count, local_center_mismatch, cell_shift_min, cell_shift_max
    real(8), allocatable :: prod_weight_local(:), prod_weight_sum(:)
    real(8), allocatable :: prod_contrib_local(:), prod_contrib_sum(:)
    real(8), allocatable :: prod_rho_local(:), prod_rho_sum(:), prod_occ_local(:), prod_occ_sum_by_w(:)
    complex(8), allocatable :: cw_local(:,:), cw_sum(:,:), cw_eig(:,:), cmix(:,:), cmix_occ(:,:), rho_mix(:,:)
    complex(8), allocatable :: coef_block(:,:)
    complex(8), allocatable :: center_w_global(:,:), center_eig(:,:)
    complex(8) :: pos_expect, pos_ww, pos_wp, zterm
    character(32) :: env_trace
    integer :: env_status, env_len
    logical :: trace_mixed_pol, perf_count
    real(8) :: perf_t0

    polarization_raw(:) = 0.0d0
    pol_local(:) = 0.0d0
    pol_sum(:) = 0.0d0
    pol_ww_local(:) = 0.0d0
    pol_wp_local(:) = 0.0d0
    pol_ww_sum(:) = 0.0d0
    pol_wp_sum(:) = 0.0d0
    pol_ww_diag_local(:) = 0.0d0
    pol_ww_offdiag_local(:) = 0.0d0
    pol_ww_same_owner_local(:) = 0.0d0
    pol_ww_cross_owner_local(:) = 0.0d0
    pol_ww_diag_sum(:) = 0.0d0
    pol_ww_offdiag_sum(:) = 0.0d0
    pol_ww_same_owner_sum(:) = 0.0d0
    pol_ww_cross_owner_sum(:) = 0.0d0
    ww_pair_count_total = 0.0d0
    ww_pair_count_diag = 0.0d0
    ww_pair_count_offdiag = 0.0d0
    ww_pair_count_cross_owner = 0.0d0
    occ_local = 0.0d0
    occ_sum = 0.0d0
    w_occ_local = 0.0d0
    w_occ_sum = 0.0d0
    zww_diag_sum = 0.0d0
    center_diag_sum = 0.0d0
    center_local_sum = 0.0d0
    zww_diag_min = huge(1.0d0)
    zww_diag_max = -huge(1.0d0)
    center_min = huge(1.0d0)
    center_max = -huge(1.0d0)
    diff_min = huge(1.0d0)
    diff_max = -huge(1.0d0)
    diff_rms = 0.0d0
    weighted_zww = 0.0d0
    weighted_center = 0.0d0
    weighted_diff = 0.0d0
    center_eig_weighted_sum = 0.0d0
    center_diag_weighted_local(:) = 0.0d0
    center_diag_weighted_sum(:) = 0.0d0
    center_eig_weighted_local(:) = 0.0d0
    center_eig_weighted_sum_global(:) = 0.0d0
    weighted_zww_sum = 0.0d0
    weighted_center_sum = 0.0d0
    weighted_diff_sum = 0.0d0
    local_center_count = 0.0d0
    local_center_mismatch = 0.0d0
    cell_shift_min = huge(1.0d0)
    cell_shift_max = -huge(1.0d0)
    dg_frag%mixed_z_prod_pz_ww_raw = 0.0d0
    dg_frag%mixed_z_prod_pz_ww_diag_raw = 0.0d0
    dg_frag%mixed_z_prod_pz_ww_offdiag_raw = 0.0d0
    dg_frag%mixed_z_prod_pz_ww_same_owner_raw = 0.0d0
    dg_frag%mixed_z_prod_pz_ww_cross_owner_raw = 0.0d0
    dg_frag%mixed_z_prod_pz_ww_pair_count_total = 0.0d0
    dg_frag%mixed_z_prod_pz_ww_pair_count_diag = 0.0d0
    dg_frag%mixed_z_prod_pz_ww_pair_count_offdiag = 0.0d0
    dg_frag%mixed_z_prod_pz_ww_pair_count_cross_owner = 0.0d0
    dg_frag%mixed_z_prod_zww_diag_sum = 0.0d0
    dg_frag%mixed_z_prod_center_diag_sum = 0.0d0
    dg_frag%mixed_z_prod_center_diag_local_sum = 0.0d0
    dg_frag%mixed_z_prod_zww_diag_min = 0.0d0
    dg_frag%mixed_z_prod_zww_diag_max = 0.0d0
    dg_frag%mixed_z_prod_zww_diag_mean = 0.0d0
    dg_frag%mixed_z_prod_center_z_min = 0.0d0
    dg_frag%mixed_z_prod_center_z_max = 0.0d0
    dg_frag%mixed_z_prod_center_z_mean = 0.0d0
    dg_frag%mixed_z_prod_diag_minus_center_min = 0.0d0
    dg_frag%mixed_z_prod_diag_minus_center_max = 0.0d0
    dg_frag%mixed_z_prod_diag_minus_center_rms = 0.0d0
    dg_frag%mixed_z_prod_weighted_zww_diag_sum = 0.0d0
    dg_frag%mixed_z_prod_weighted_center_sum = 0.0d0
    dg_frag%mixed_z_prod_weighted_diff_sum = 0.0d0
    dg_frag%mixed_z_prod_wrap_shift_count = 0.0d0
    dg_frag%mixed_z_prod_cell_shift_min = 0.0d0
    dg_frag%mixed_z_prod_cell_shift_max = 0.0d0
    dg_frag%mixed_z_prod_owner_gid_mismatch_count = 0.0d0
    dg_frag%mixed_z_prod_center_source_mismatch_count = 0.0d0
    dg_frag%mixed_z_prod_pz_wp_raw = 0.0d0
    dg_frag%mixed_z_prod_pz_pp_raw = 0.0d0
    dg_frag%mixed_z_prod_pz_occ_sum = 0.0d0
    dg_frag%mixed_z_prod_pz_w_occ_weight = 0.0d0
    dg_frag%mixed_z_prod_pz_w_dim = 0.0d0
    dg_frag%mixed_z_prod_pz_decomp_ready = .false.
    if (allocated(dg_frag%mixed_z_prod_zww_diag_by_w)) deallocate(dg_frag%mixed_z_prod_zww_diag_by_w)
    if (allocated(dg_frag%mixed_z_prod_ww_diag_weight_by_w)) deallocate(dg_frag%mixed_z_prod_ww_diag_weight_by_w)
    if (allocated(dg_frag%mixed_z_prod_ww_diag_contrib_by_w)) deallocate(dg_frag%mixed_z_prod_ww_diag_contrib_by_w)
    if (allocated(dg_frag%mixed_z_prod_ww_diag_rho_by_w)) deallocate(dg_frag%mixed_z_prod_ww_diag_rho_by_w)
    if (allocated(dg_frag%mixed_z_prod_ww_diag_occ_by_w)) deallocate(dg_frag%mixed_z_prod_ww_diag_occ_by_w)
    env_trace = ''
    trace_mixed_pol = .false.
    perf_count = dg_frag%mixed_z_perf_count_enabled
    call get_environment_variable('SALMON_DG_MIXED_Z_TRACE', env_trace, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_trace(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        trace_mixed_pol = .true.
      end select
    end if
    if (.not. dg_frag%has_mixed_wannier_bpw_position) return
    if (.not. allocated(dg_frag%mixed_wannier_bpw_z) .or. .not. allocated(dg_frag%mixed_wannier_bpw_pcoef)) return
    if (.not. allocated(dg_frag%global_wannier_coef) .or. .not. allocated(dg_frag%coef_global_to_local)) return

    n_w = dg_frag%mixed_wannier_bpw_nw
    n_p = dg_frag%mixed_wannier_bpw_np
    n_mix = dg_frag%mixed_wannier_bpw_nmix
    if (n_w <= 0 .or. n_p < 0 .or. n_mix /= n_w + n_p) return
    n_wann_raw = n_w
    if (allocated(dg_frag%global_wannier_flux_evec)) n_wann_raw = size(dg_frag%global_wannier_flux_evec, 1)
    ww_pair_count_total = dble(n_w) * dble(n_w)
    ww_pair_count_diag = dble(n_w)
    ww_pair_count_offdiag = max(0.0d0, ww_pair_count_total - ww_pair_count_diag)
    if (allocated(dg_frag%global_wannier_owner_frag)) then
      do iw = 1, min(n_w, size(dg_frag%global_wannier_owner_frag))
        do jw = 1, min(n_w, size(dg_frag%global_wannier_owner_frag))
          if (dg_frag%global_wannier_owner_frag(iw) /= dg_frag%global_wannier_owner_frag(jw)) then
            ww_pair_count_cross_owner = ww_pair_count_cross_owner + 1.0d0
          end if
        end do
      end do
    end if
    do iw = 1, n_w
      zww_diag_val = real(dg_frag%mixed_wannier_bpw_z(3, iw, iw, 1), kind=8)
      zww_diag_sum = zww_diag_sum + zww_diag_val
      zww_diag_min = min(zww_diag_min, zww_diag_val)
      zww_diag_max = max(zww_diag_max, zww_diag_val)
      if (allocated(dg_frag%global_wannier_center) .and. iw <= size(dg_frag%global_wannier_center, 2)) then
        center_val = dg_frag%global_wannier_center(3, iw)
        center_diag_sum = center_diag_sum + center_val
        center_min = min(center_min, center_val)
        center_max = max(center_max, center_val)
        diff_val = zww_diag_val - center_val
        diff_min = min(diff_min, diff_val)
        diff_max = max(diff_max, diff_val)
        diff_rms = diff_rms + diff_val * diff_val
        cell_shift_min = min(cell_shift_min, diff_val)
        cell_shift_max = max(cell_shift_max, diff_val)
        if (abs(diff_val) > 1.0d-8) local_center_mismatch = local_center_mismatch + 1.0d0
      end if
    end do
    if (allocated(dg_frag%global_wannier_local_ids) .and. allocated(dg_frag%global_wannier_local_center)) then
      do i_local = 1, min(size(dg_frag%global_wannier_local_ids, 2), size(dg_frag%global_wannier_local_center, 3))
        do iw = 1, min(size(dg_frag%global_wannier_local_ids, 1), size(dg_frag%global_wannier_local_center, 2))
          if (dg_frag%global_wannier_local_ids(iw, i_local) < 1 .or. &
              dg_frag%global_wannier_local_ids(iw, i_local) > n_w) cycle
          center_local_sum = center_local_sum + dg_frag%global_wannier_local_center(3, iw, i_local)
          local_center_count = local_center_count + 1.0d0
        end do
      end do
    end if
    if (n_w > 0) then
      diff_rms = sqrt(max(diff_rms / dble(n_w), 0.0d0))
    end if
    if (zww_diag_min == huge(1.0d0)) zww_diag_min = 0.0d0
    if (zww_diag_max == -huge(1.0d0)) zww_diag_max = 0.0d0
    if (center_min == huge(1.0d0)) center_min = 0.0d0
    if (center_max == -huge(1.0d0)) center_max = 0.0d0
    if (diff_min == huge(1.0d0)) diff_min = 0.0d0
    if (diff_max == -huge(1.0d0)) diff_max = 0.0d0
    if (cell_shift_min == huge(1.0d0)) cell_shift_min = 0.0d0
    if (cell_shift_max == -huge(1.0d0)) cell_shift_max = 0.0d0

    allocate(dg_frag%mixed_z_prod_zww_diag_by_w(n_w), &
             dg_frag%mixed_z_prod_ww_diag_weight_by_w(n_w), &
             dg_frag%mixed_z_prod_ww_diag_contrib_by_w(n_w), &
             dg_frag%mixed_z_prod_ww_diag_rho_by_w(n_w), &
             dg_frag%mixed_z_prod_ww_diag_occ_by_w(n_w))
    allocate(prod_weight_local(n_w), prod_weight_sum(n_w), prod_contrib_local(n_w), prod_contrib_sum(n_w))
    allocate(prod_rho_local(n_w), prod_rho_sum(n_w), prod_occ_local(n_w), prod_occ_sum_by_w(n_w))
    dg_frag%mixed_z_prod_zww_diag_by_w(:) = 0.0d0
    dg_frag%mixed_z_prod_ww_diag_weight_by_w(:) = 0.0d0
    dg_frag%mixed_z_prod_ww_diag_contrib_by_w(:) = 0.0d0
    dg_frag%mixed_z_prod_ww_diag_rho_by_w(:) = 0.0d0
    dg_frag%mixed_z_prod_ww_diag_occ_by_w(:) = 0.0d0
    prod_weight_local(:) = 0.0d0
    prod_weight_sum(:) = 0.0d0
    prod_contrib_local(:) = 0.0d0
    prod_contrib_sum(:) = 0.0d0
    prod_rho_local(:) = 0.0d0
    prod_rho_sum(:) = 0.0d0
    prod_occ_local(:) = 0.0d0
    prod_occ_sum_by_w(:) = 0.0d0
    do iw = 1, n_w
      dg_frag%mixed_z_prod_zww_diag_by_w(iw) = real(dg_frag%mixed_wannier_bpw_z(3, iw, iw, 1), kind=8)
    end do

    allocate(cw_local(n_wann_raw,max(1,dg_frag%nstate_tot)), cw_sum(n_wann_raw,max(1,dg_frag%nstate_tot)))
    allocate(cw_eig(n_w,max(1,dg_frag%nstate_tot)))
    allocate(cmix(n_mix,max(1,dg_frag%nstate_tot)))
    allocate(cmix_occ(n_mix,max(1,dg_frag%nstate_tot)), rho_mix(n_mix,n_mix))
    allocate(coef_block(max(1,size(dg_frag%global_wannier_coef,1)), max(1,dg_frag%nstate_tot)))
    if (trim(dg_mixed_z_polarization_branch) == 'center_eig') then
      if (.not. allocated(dg_frag%global_wannier_center) .or. &
          .not. allocated(dg_frag%global_wannier_flux_evec) .or. &
          size(dg_frag%global_wannier_center, 1) < 3 .or. &
          size(dg_frag%global_wannier_center, 2) < n_w .or. &
          size(dg_frag%global_wannier_flux_evec, 1) < n_w .or. &
          size(dg_frag%global_wannier_flux_evec, 2) < n_w) then
        stop "DG-Fragment RT: center_eig polarization branch requires Wannier centers and Flux eigenvectors"
      end if
      allocate(center_w_global(n_w,n_w), center_eig(n_w,n_w))
    end if

    do ispin = 1, min(dg_frag%nspin, system%nspin, size(dg_frag%mixed_wannier_bpw_z, 4))
      nocc_spin = 0
      if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
      nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(dg_frag%coef, 2), size(system%rocc, 1))
      if (nocc_spin <= 0) cycle

      cw_local(:, :) = (0.0d0, 0.0d0)
      cw_eig(:, :) = (0.0d0, 0.0d0)
      if (perf_count) perf_t0 = get_wtime()
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = ifrag - dg_frag%ifrag_start + 1
        if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
        nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%global_wannier_coef, 1))
        if (nbasis <= 0) cycle
        coef_block(1:nbasis,1:nocc_spin) = (0.0d0, 0.0d0)
        do ib = 1, nbasis
          global_row = dg_frag%index_basis(ib, ifrag, ispin)
          if (global_row < 1 .or. global_row > dg_frag%n_mat_max) cycle
          local_row = dg_frag%coef_global_to_local(global_row, ispin)
          if (local_row < 1 .or. local_row > size(dg_frag%coef, 1)) cycle
          coef_block(ib,1:nocc_spin) = dg_frag%coef(local_row,1:nocc_spin,ispin)
        end do
        call zgemm('C', 'N', n_wann_raw, nocc_spin, nbasis, (1.0d0, 0.0d0), &
          dg_frag%global_wannier_coef(1:nbasis,1:n_wann_raw,ispin,i_local), max(1,nbasis), &
          coef_block, max(1,size(coef_block,1)), (1.0d0, 0.0d0), cw_local, max(1,n_wann_raw))
        if (perf_count) dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 1_8
      end do
      if (perf_count) dg_frag%mixed_z_perf_wall_pz_build_cw = &
        dg_frag%mixed_z_perf_wall_pz_build_cw + (get_wtime() - perf_t0)
      if (perf_count) perf_t0 = get_wtime()
      call comm_summation(cw_local, cw_sum, n_wann_raw * max(1,dg_frag%nstate_tot), dg_frag%icomm)
      if (perf_count) dg_frag%mixed_z_perf_wall_pz_reduce = &
        dg_frag%mixed_z_perf_wall_pz_reduce + (get_wtime() - perf_t0)

      if (dg_frag%id == 0) then
        if (perf_count) perf_t0 = get_wtime()
        if (allocated(dg_frag%global_wannier_flux_evec)) then
          cw_eig(1:n_w,1:nocc_spin) = matmul( &
            conjg(transpose(dg_frag%global_wannier_flux_evec(1:n_wann_raw,1:n_w))), &
            cw_sum(1:n_wann_raw,1:nocc_spin))
        else
          cw_eig(1:n_w,1:nocc_spin) = cw_sum(1:n_w,1:nocc_spin)
        end if
        cmix(:, :) = (0.0d0, 0.0d0)
        cmix_occ(:, :) = (0.0d0, 0.0d0)
        rho_mix(:, :) = (0.0d0, 0.0d0)
        cmix(1:n_w,1:nocc_spin) = cw_eig(1:n_w,1:nocc_spin)
        if (n_p > 0) cmix(n_w+1:n_w+n_p,1:nocc_spin) = dg_frag%mixed_wannier_bpw_pcoef(1:n_p,1:nocc_spin,ispin)
        do istate = 1, nocc_spin
          occ = max(0.0d0, system%rocc(istate, 1, ispin))
          if (occ <= 0.0d0) cycle
          occ_local = occ_local + occ
          if (allocated(dg_frag%global_wannier_center)) then
            do iw = 1, min(n_wann_raw, size(dg_frag%global_wannier_center, 2))
              weight_diag = occ * real(conjg(cw_sum(iw,istate)) * cw_sum(iw,istate), kind=8)
              weighted_center = weighted_center + weight_diag * dg_frag%global_wannier_center(3, iw)
              do idir = 1, 3
                center_diag_weighted_local(idir) = center_diag_weighted_local(idir) + &
                  weight_diag * dg_frag%global_wannier_center(idir, iw)
              end do
            end do
          end if
          cmix_occ(1:n_mix,istate) = cmix(1:n_mix,istate) * cmplx(occ, 0.0d0, kind=8)
          do iw = 1, n_w
            rho_diag = real(conjg(cmix(iw,istate)) * cmix(iw,istate), kind=8)
            prod_rho_local(iw) = prod_rho_local(iw) + rho_diag
            prod_occ_local(iw) = prod_occ_local(iw) + occ
          end do
        end do
        call zgemm('N', 'C', n_mix, n_mix, nocc_spin, (1.0d0, 0.0d0), cmix_occ, max(1,n_mix), &
          cmix, max(1,n_mix), (0.0d0, 0.0d0), rho_mix, max(1,n_mix))
        if (perf_count) dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 1_8
        do iw = 1, n_w
          weight_diag = real(rho_mix(iw,iw), kind=8)
          w_occ_local = w_occ_local + weight_diag
          weighted_zww = weighted_zww + weight_diag * real(dg_frag%mixed_wannier_bpw_z(3, iw, iw, ispin), kind=8)
            prod_weight_local(iw) = prod_weight_local(iw) + weight_diag
            prod_contrib_local(iw) = prod_contrib_local(iw) - weight_diag * &
              real(dg_frag%mixed_wannier_bpw_z(3, iw, iw, ispin), kind=8)
            if (allocated(dg_frag%global_wannier_center) .and. iw <= size(dg_frag%global_wannier_center, 2)) then
              weighted_diff = weighted_diff + weight_diag * &
                (real(dg_frag%mixed_wannier_bpw_z(3, iw, iw, ispin), kind=8) - &
                 dg_frag%global_wannier_center(3, iw))
            end if
        end do
        do idir = 1, 3
          pos_expect = (0.0d0, 0.0d0)
          pos_ww = (0.0d0, 0.0d0)
          pos_wp = (0.0d0, 0.0d0)
          do iw = 1, n_mix
            do jw = 1, n_mix
              zterm = dg_frag%mixed_wannier_bpw_z(idir,iw,jw,ispin) * rho_mix(jw,iw)
              pos_expect = pos_expect + zterm
              if (iw <= n_w .and. jw <= n_w) then
                pos_ww = pos_ww + zterm
                if (iw == jw) then
                  pol_ww_diag_local(idir) = pol_ww_diag_local(idir) - real(zterm, kind=8)
                else
                  pol_ww_offdiag_local(idir) = pol_ww_offdiag_local(idir) - real(zterm, kind=8)
                end if
                if (allocated(dg_frag%global_wannier_owner_frag)) then
                  if (iw <= size(dg_frag%global_wannier_owner_frag) .and. &
                      jw <= size(dg_frag%global_wannier_owner_frag)) then
                    if (dg_frag%global_wannier_owner_frag(iw) == dg_frag%global_wannier_owner_frag(jw)) then
                      pol_ww_same_owner_local(idir) = pol_ww_same_owner_local(idir) - real(zterm, kind=8)
                    else
                      pol_ww_cross_owner_local(idir) = pol_ww_cross_owner_local(idir) - real(zterm, kind=8)
                    end if
                  end if
                end if
              else if ((iw <= n_w .and. jw > n_w) .or. (iw > n_w .and. jw <= n_w)) then
                pos_wp = pos_wp + zterm
              end if
            end do
          end do
          pol_local(idir) = pol_local(idir) - real(pos_expect, kind=8)
          pol_ww_local(idir) = pol_ww_local(idir) - real(pos_ww, kind=8)
          pol_wp_local(idir) = pol_wp_local(idir) - real(pos_wp, kind=8)
        end do
        if (allocated(center_eig)) then
          do idir = 1, 3
            center_w_global(:, :) = (0.0d0, 0.0d0)
            do iw = 1, n_w
              center_w_global(iw,iw) = cmplx(dg_frag%global_wannier_center(idir,iw), 0.0d0, kind=8)
            end do
            center_eig(1:n_w,1:n_w) = matmul(conjg(transpose(dg_frag%global_wannier_flux_evec(1:n_w,1:n_w))), &
              matmul(center_w_global(1:n_w,1:n_w), dg_frag%global_wannier_flux_evec(1:n_w,1:n_w)))
            pos_ww = (0.0d0, 0.0d0)
            do iw = 1, n_w
              do jw = 1, n_w
                pos_ww = pos_ww + center_eig(iw,jw) * rho_mix(jw,iw)
              end do
            end do
            center_eig_weighted_local(idir) = center_eig_weighted_local(idir) + real(pos_ww, kind=8)
          end do
        end if
        if (perf_count) dg_frag%mixed_z_perf_wall_pz_contract_z = &
          dg_frag%mixed_z_perf_wall_pz_contract_z + (get_wtime() - perf_t0)
      end if
    end do

    if (perf_count) perf_t0 = get_wtime()
    call comm_summation(pol_local, pol_sum, 3, dg_frag%icomm)
    call comm_summation(pol_ww_local, pol_ww_sum, 3, dg_frag%icomm)
    call comm_summation(pol_wp_local, pol_wp_sum, 3, dg_frag%icomm)
    call comm_summation(pol_ww_diag_local, pol_ww_diag_sum, 3, dg_frag%icomm)
    call comm_summation(pol_ww_offdiag_local, pol_ww_offdiag_sum, 3, dg_frag%icomm)
    call comm_summation(pol_ww_same_owner_local, pol_ww_same_owner_sum, 3, dg_frag%icomm)
    call comm_summation(pol_ww_cross_owner_local, pol_ww_cross_owner_sum, 3, dg_frag%icomm)
    call comm_summation(occ_local, occ_sum, dg_frag%icomm)
    call comm_summation(w_occ_local, w_occ_sum, dg_frag%icomm)
    call comm_summation(weighted_zww, weighted_zww_sum, dg_frag%icomm)
    call comm_summation(weighted_center, weighted_center_sum, dg_frag%icomm)
    call comm_summation(weighted_diff, weighted_diff_sum, dg_frag%icomm)
    call comm_summation(center_diag_weighted_local, center_diag_weighted_sum, 3, dg_frag%icomm)
    call comm_summation(center_eig_weighted_local, center_eig_weighted_sum_global, 3, dg_frag%icomm)
    call comm_summation(prod_weight_local, prod_weight_sum, n_w, dg_frag%icomm)
    call comm_summation(prod_contrib_local, prod_contrib_sum, n_w, dg_frag%icomm)
    call comm_summation(prod_rho_local, prod_rho_sum, n_w, dg_frag%icomm)
    call comm_summation(prod_occ_local, prod_occ_sum_by_w, n_w, dg_frag%icomm)
    if (perf_count) dg_frag%mixed_z_perf_wall_pz_reduce = &
      dg_frag%mixed_z_perf_wall_pz_reduce + (get_wtime() - perf_t0)
    if (pol_sum(1) /= pol_sum(1) .or. pol_sum(2) /= pol_sum(2) .or. pol_sum(3) /= pol_sum(3)) then
      stop "DG-Fragment RT: NaN in mixed Wannier+BPW polarization reduction"
    end if
    if (trace_mixed_pol .and. dg_frag%id == 0) then
      write(*,'(1x,a,3(1x,1pe13.5),a,3(1x,1pe13.5),a,3(1x,1pe13.5))') &
        "[DG-MIXED-POL] total=", pol_sum(:), " WW=", pol_ww_sum(:), " WP+PW=", pol_wp_sum(:)
    end if
    dg_frag%mixed_z_prod_pz_ww_raw = pol_ww_sum(3)
    dg_frag%mixed_z_prod_pz_ww_diag_raw = pol_ww_diag_sum(3)
    dg_frag%mixed_z_prod_pz_ww_offdiag_raw = pol_ww_offdiag_sum(3)
    dg_frag%mixed_z_prod_pz_ww_same_owner_raw = pol_ww_same_owner_sum(3)
    dg_frag%mixed_z_prod_pz_ww_cross_owner_raw = pol_ww_cross_owner_sum(3)
    dg_frag%mixed_z_prod_pz_ww_pair_count_total = ww_pair_count_total
    dg_frag%mixed_z_prod_pz_ww_pair_count_diag = ww_pair_count_diag
    dg_frag%mixed_z_prod_pz_ww_pair_count_offdiag = ww_pair_count_offdiag
    dg_frag%mixed_z_prod_pz_ww_pair_count_cross_owner = ww_pair_count_cross_owner
    dg_frag%mixed_z_prod_zww_diag_sum = zww_diag_sum
    dg_frag%mixed_z_prod_center_diag_sum = center_diag_sum
    dg_frag%mixed_z_prod_center_diag_local_sum = center_local_sum
    dg_frag%mixed_z_prod_zww_diag_min = zww_diag_min
    dg_frag%mixed_z_prod_zww_diag_max = zww_diag_max
    dg_frag%mixed_z_prod_zww_diag_mean = zww_diag_sum / max(dble(n_w), 1.0d0)
    dg_frag%mixed_z_prod_center_z_min = center_min
    dg_frag%mixed_z_prod_center_z_max = center_max
    dg_frag%mixed_z_prod_center_z_mean = center_diag_sum / max(dble(n_w), 1.0d0)
    dg_frag%mixed_z_prod_diag_minus_center_min = diff_min
    dg_frag%mixed_z_prod_diag_minus_center_max = diff_max
    dg_frag%mixed_z_prod_diag_minus_center_rms = diff_rms
    dg_frag%mixed_z_prod_weighted_zww_diag_sum = weighted_zww_sum
    dg_frag%mixed_z_prod_ww_diag_weight_by_w(:) = prod_weight_sum(:)
    dg_frag%mixed_z_prod_ww_diag_contrib_by_w(:) = prod_contrib_sum(:)
    dg_frag%mixed_z_prod_ww_diag_rho_by_w(:) = prod_rho_sum(:)
    dg_frag%mixed_z_prod_ww_diag_occ_by_w(:) = prod_occ_sum_by_w(:)
    dg_frag%mixed_z_prod_weighted_center_sum = weighted_center_sum
    dg_frag%mixed_z_prod_weighted_diff_sum = weighted_diff_sum
    dg_frag%mixed_z_prod_wrap_shift_count = local_center_mismatch
    dg_frag%mixed_z_prod_cell_shift_min = cell_shift_min
    dg_frag%mixed_z_prod_cell_shift_max = cell_shift_max
    dg_frag%mixed_z_prod_owner_gid_mismatch_count = 0.0d0
    dg_frag%mixed_z_prod_center_source_mismatch_count = local_center_mismatch
    dg_frag%mixed_z_prod_pz_wp_raw = pol_wp_sum(3)
    dg_frag%mixed_z_prod_pz_pp_raw = pol_sum(3) - pol_ww_sum(3) - pol_wp_sum(3)
    dg_frag%mixed_z_prod_pz_occ_sum = occ_sum
    dg_frag%mixed_z_prod_pz_w_occ_weight = w_occ_sum
    dg_frag%mixed_z_prod_pz_w_dim = dble(n_w)
    dg_frag%mixed_z_prod_pz_decomp_ready = .true.
    polarization_raw(:) = pol_sum(:)
    if (trim(dg_mixed_z_polarization_branch) == 'center_diag') then
      if (.not. allocated(dg_frag%global_wannier_center) .or. &
          size(dg_frag%global_wannier_center, 1) < 3 .or. &
          size(dg_frag%global_wannier_center, 2) < n_w) then
        stop "DG-Fragment RT: center_diag polarization branch requires global Wannier centers"
      end if
      polarization_raw(:) = -center_diag_weighted_sum(:) + pol_ww_offdiag_sum(:) + pol_wp_sum(:) + &
        (pol_sum(:) - pol_ww_sum(:) - pol_wp_sum(:))
    else if (trim(dg_mixed_z_polarization_branch) == 'center_eig') then
      polarization_raw(:) = -center_eig_weighted_sum_global(:) + pol_wp_sum(:) + &
        (pol_sum(:) - pol_ww_sum(:) - pol_wp_sum(:))
    end if
    if (allocated(center_w_global)) deallocate(center_w_global, center_eig)
    deallocate(cw_local, cw_sum, cw_eig, cmix, cmix_occ, rho_mix, coef_block, prod_weight_local, prod_weight_sum, prod_contrib_local, prod_contrib_sum, &
               prod_rho_local, prod_rho_sum, prod_occ_local, prod_occ_sum_by_w)
  end subroutine calculate_global_mixed_wannier_bpw_polarization_dg

  subroutine calculate_formal_dg_wannier_polarization_dg(dg_frag, system, polarization_raw)
    use structures, only: s_dft_system
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(out)   :: polarization_raw(3)

    integer, parameter :: pol_coef_cache_target_mb = 64
    integer :: ispin, ifrag, i_local, ib, iw, jw, idir
    integer :: istate, occ0, nbatch, state_work
    integer :: nocc_spin, occ_first, occ_last
    integer :: n_global, n_needed, gid, pos, nbasis, nkeep
    integer :: local_idx, local_col
    integer, allocatable :: needed_ids(:), needed_pos(:)
    logical, allocatable :: row_needed(:)
    complex(8), allocatable :: coef_work(:,:), cw(:,:)
    complex(8) :: pos_mat, pos_expect
    real(8) :: pol_local(3), pol_sum(3), occ
    integer(8) :: target_bytes, bytes_per_state

    polarization_raw(:) = 0.0d0
    pol_local(:) = 0.0d0
    pol_sum(:) = 0.0d0

    if (.not. dg_frag%has_formal_dg_wannier_basis) return
    if (.not. allocated(dg_frag%dg_wannier_basis_coef)) return
    if (.not. allocated(dg_frag%dg_wannier_xi_local)) return
    if (.not. allocated(dg_frag%dg_wannier_ref_center)) return
    if (.not. allocated(dg_frag%dg_wannier_nkeep)) return
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. allocated(dg_frag%n_basis)) return

    n_global = size(dg_frag%coef_owner, 1)
    allocate(row_needed(n_global), needed_pos(n_global))

    do ispin = 1, min(dg_frag%nspin, system%nspin)
      nocc_spin = 0
      if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
      if (dg_frag%coef_state_block_mode) then
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(system%rocc, 1))
      else
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(dg_frag%coef, 2), size(system%rocc, 1))
      end if
      if (nocc_spin <= 0) cycle

      occ_first = 1
      occ_last = nocc_spin
      if (dg_frag%coef_state_block_mode) then
        occ_first = max(1, dg_frag%coef_state_start)
        occ_last = min(nocc_spin, dg_frag%coef_state_end)
        if (occ_first > occ_last) cycle
      end if

      row_needed(:) = .false.
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = ifrag - dg_frag%ifrag_start + 1
        if (i_local < 1 .or. i_local > size(dg_frag%dg_wannier_nkeep)) cycle
        nkeep = dg_frag%dg_wannier_nkeep(i_local)
        if (nkeep <= 0) cycle
        nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                     size(dg_frag%dg_wannier_basis_coef, 1))
        do ib = 1, nbasis
          gid = dg_frag%index_basis(ib, ifrag, ispin)
          if (gid >= 1 .and. gid <= n_global) row_needed(gid) = .true.
        end do
      end do
      n_needed = count(row_needed)
      if (n_needed <= 0) cycle
      if (allocated(needed_ids)) deallocate(needed_ids)
      allocate(needed_ids(n_needed))
      needed_pos(:) = 0
      n_needed = 0
      do gid = 1, n_global
        if (.not. row_needed(gid)) cycle
        n_needed = n_needed + 1
        needed_ids(n_needed) = gid
        needed_pos(gid) = n_needed
      end do

      target_bytes = int(pol_coef_cache_target_mb, kind=8) * 1024_8 * 1024_8
      bytes_per_state = 16_8 * int(max(1, n_needed + size(dg_frag%dg_wannier_basis_coef, 2)), kind=8)
      state_work = max(1, min(occ_last - occ_first + 1, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))

      if (allocated(coef_work)) deallocate(coef_work)
      if (allocated(cw)) deallocate(cw)
      allocate(coef_work(n_needed, state_work))
      allocate(cw(size(dg_frag%dg_wannier_basis_coef, 2), state_work))

      do occ0 = occ_first, occ_last, state_work
        nbatch = min(state_work, occ_last - occ0 + 1)
        coef_work(:, :) = (0.0d0, 0.0d0)
        if (dg_frag%coef_state_block_mode) then
          do pos = 1, n_needed
            local_idx = dg_frag%coef_global_to_local(needed_ids(pos), ispin)
            if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
            do istate = 1, nbatch
              local_col = occ0 + istate - 1 - dg_frag%coef_state_start + 1
              if (local_col < 1 .or. local_col > size(dg_frag%coef, 2)) cycle
              coef_work(pos, istate) = dg_frag%coef(local_idx, local_col, ispin)
            end do
          end do
        else
          call fetch_remote_coef_rows(dg_frag, ispin, needed_ids, coef_work(1:n_needed, 1:nbatch), &
                                      occ0, occ0 + nbatch - 1)
        end if

        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = ifrag - dg_frag%ifrag_start + 1
          if (i_local < 1 .or. i_local > size(dg_frag%dg_wannier_nkeep)) cycle
          nkeep = dg_frag%dg_wannier_nkeep(i_local)
          if (nkeep <= 0) cycle
          nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                       size(dg_frag%dg_wannier_basis_coef, 1))
          cw(1:nkeep, 1:nbatch) = (0.0d0, 0.0d0)
          do iw = 1, nkeep
            do ib = 1, nbasis
              gid = dg_frag%index_basis(ib, ifrag, ispin)
              if (gid < 1 .or. gid > n_global) cycle
              pos = needed_pos(gid)
              if (pos <= 0) cycle
!$omp simd
              do istate = 1, nbatch
                cw(iw, istate) = cw(iw, istate) + &
                  conjg(dg_frag%dg_wannier_basis_coef(ib, iw, ispin, i_local)) * coef_work(pos, istate)
              end do
            end do
          end do

          do istate = 1, nbatch
            occ = max(0.0d0, system%rocc(occ0 + istate - 1, 1, ispin))
            if (occ <= 0.0d0) cycle
            do idir = 1, 3
              pos_expect = (0.0d0, 0.0d0)
              do iw = 1, nkeep
                do jw = 1, nkeep
                  pos_mat = dg_frag%dg_wannier_xi_local(idir, iw, jw, ispin, i_local)
                  if (iw == jw) pos_mat = pos_mat + cmplx(dg_frag%dg_wannier_ref_center(idir, i_local), 0.0d0, kind=8)
                  pos_expect = pos_expect + conjg(cw(iw, istate)) * pos_mat * cw(jw, istate)
                end do
              end do
              pol_local(idir) = pol_local(idir) - occ * real(pos_expect, kind=8)
            end do
          end do
        end do
      end do
    end do

    call comm_summation(pol_local, pol_sum, 3, dg_frag%icomm)
    if (pol_sum(1) /= pol_sum(1) .or. pol_sum(2) /= pol_sum(2) .or. pol_sum(3) /= pol_sum(3)) then
      stop "DG-Fragment RT: NaN in formal DG-Wannier polarization reduction"
    end if
    polarization_raw(:) = pol_sum(:)
  end subroutine calculate_formal_dg_wannier_polarization_dg

  subroutine calculate_global_wannier_local_center_polarization_dg(dg_frag, system, polarization_raw)
    use structures, only: s_dft_system
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(out)   :: polarization_raw(3)

    integer, parameter :: pol_coef_cache_target_mb = 64
    integer :: ispin, ifrag, i_local, ib, iw, jw, idir
    integer :: istate, occ0, nbatch, state_work
    integer :: nocc_spin, occ_first, occ_last
    integer :: n_global, n_needed, gid, pos, nbasis, nkeep
    integer :: local_idx, local_col
    integer, allocatable :: needed_ids(:), needed_pos(:)
    logical, allocatable :: row_needed(:)
    complex(8), allocatable :: coef_work(:,:), cw(:,:)
    complex(8) :: pos_expect
    real(8) :: pol_local(3), pol_sum(3), occ, weight
    integer(8) :: target_bytes, bytes_per_state

    polarization_raw(:) = 0.0d0
    pol_local(:) = 0.0d0
    pol_sum(:) = 0.0d0

    if (.not. dg_frag%has_global_wannier_local_basis) return
    if (.not. allocated(dg_frag%global_wannier_local_coef)) return
    if (.not. allocated(dg_frag%global_wannier_local_center)) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. allocated(dg_frag%n_basis)) return

    n_global = size(dg_frag%coef_owner, 1)
    allocate(row_needed(n_global), needed_pos(n_global))

    do ispin = 1, min(dg_frag%nspin, system%nspin)
      nocc_spin = 0
      if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
      if (dg_frag%coef_state_block_mode) then
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(system%rocc, 1))
      else
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(dg_frag%coef, 2), size(system%rocc, 1))
      end if
      if (nocc_spin <= 0) cycle

      occ_first = 1
      occ_last = nocc_spin
      if (dg_frag%coef_state_block_mode) then
        occ_first = max(1, dg_frag%coef_state_start)
        occ_last = min(nocc_spin, dg_frag%coef_state_end)
        if (occ_first > occ_last) cycle
      end if

      row_needed(:) = .false.
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = ifrag - dg_frag%ifrag_start + 1
        if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
        nkeep = dg_frag%global_wannier_local_nkeep(i_local)
        if (nkeep <= 0) cycle
        nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                     size(dg_frag%global_wannier_local_coef, 1))
        do ib = 1, nbasis
          gid = dg_frag%index_basis(ib, ifrag, ispin)
          if (gid >= 1 .and. gid <= n_global) row_needed(gid) = .true.
        end do
      end do
      n_needed = count(row_needed)
      if (n_needed <= 0) cycle
      if (allocated(needed_ids)) deallocate(needed_ids)
      allocate(needed_ids(n_needed))
      needed_pos(:) = 0
      n_needed = 0
      do gid = 1, n_global
        if (.not. row_needed(gid)) cycle
        n_needed = n_needed + 1
        needed_ids(n_needed) = gid
        needed_pos(gid) = n_needed
      end do

      target_bytes = int(pol_coef_cache_target_mb, kind=8) * 1024_8 * 1024_8
      bytes_per_state = 16_8 * int(max(1, n_needed + size(dg_frag%global_wannier_local_coef, 2)), kind=8)
      state_work = max(1, min(occ_last - occ_first + 1, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))

      if (allocated(coef_work)) deallocate(coef_work)
      if (allocated(cw)) deallocate(cw)
      allocate(coef_work(n_needed, state_work))
      allocate(cw(size(dg_frag%global_wannier_local_coef, 2), state_work))

      do occ0 = occ_first, occ_last, state_work
        nbatch = min(state_work, occ_last - occ0 + 1)
        coef_work(:, :) = (0.0d0, 0.0d0)
        if (dg_frag%coef_state_block_mode) then
          do pos = 1, n_needed
            local_idx = dg_frag%coef_global_to_local(needed_ids(pos), ispin)
            if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
            do istate = 1, nbatch
              local_col = occ0 + istate - 1 - dg_frag%coef_state_start + 1
              if (local_col < 1 .or. local_col > size(dg_frag%coef, 2)) cycle
              coef_work(pos, istate) = dg_frag%coef(local_idx, local_col, ispin)
            end do
          end do
        else
          call fetch_remote_coef_rows(dg_frag, ispin, needed_ids, coef_work(1:n_needed, 1:nbatch), &
                                      occ0, occ0 + nbatch - 1)
        end if

        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = ifrag - dg_frag%ifrag_start + 1
          if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_local_nkeep)) cycle
          nkeep = dg_frag%global_wannier_local_nkeep(i_local)
          if (nkeep <= 0) cycle
          nbasis = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), &
                       size(dg_frag%global_wannier_local_coef, 1))
          cw(1:nkeep, 1:nbatch) = (0.0d0, 0.0d0)
          do iw = 1, nkeep
            do ib = 1, nbasis
              gid = dg_frag%index_basis(ib, ifrag, ispin)
              if (gid < 1 .or. gid > n_global) cycle
              pos = needed_pos(gid)
              if (pos <= 0) cycle
!$omp simd
              do istate = 1, nbatch
                cw(iw, istate) = cw(iw, istate) + &
                  conjg(dg_frag%global_wannier_local_coef(ib, iw, ispin, i_local)) * coef_work(pos, istate)
              end do
            end do
          end do

          do istate = 1, nbatch
            occ = max(0.0d0, system%rocc(occ0 + istate - 1, 1, ispin))
            if (occ <= 0.0d0) cycle
            if (dg_frag%has_global_wannier_position .and. allocated(dg_frag%global_wannier_local_position)) then
              do idir = 1, 3
                pos_expect = (0.0d0, 0.0d0)
                do iw = 1, nkeep
                  do jw = 1, nkeep
                    pos_expect = pos_expect + conjg(cw(iw, istate)) * &
                      dg_frag%global_wannier_local_position(idir, iw, jw, i_local) * cw(jw, istate)
                  end do
                end do
                pol_local(idir) = pol_local(idir) - occ * real(pos_expect, kind=8)
              end do
            else
              do iw = 1, nkeep
                weight = occ * real(conjg(cw(iw, istate)) * cw(iw, istate), 8)
                if (weight <= 0.0d0) cycle
                do idir = 1, 3
                  pol_local(idir) = pol_local(idir) - weight * &
                    dg_frag%global_wannier_local_center(idir, iw, i_local)
                end do
              end do
            end if
          end do
        end do
      end do
    end do

    call comm_summation(pol_local, pol_sum, 3, dg_frag%icomm)
    if (pol_sum(1) /= pol_sum(1) .or. pol_sum(2) /= pol_sum(2) .or. pol_sum(3) /= pol_sum(3)) then
      stop "DG-Fragment RT: NaN in local-selected global Wannier polarization reduction"
    end if
    polarization_raw(:) = pol_sum(:)
  end subroutine calculate_global_wannier_local_center_polarization_dg

  subroutine calculate_nonlocal_current_dg(dg_frag, system, mg, ppg, Ac_tot, current_raw)
    use structures
    use communication, only: comm_summation
    use misc_routines, only: get_wtime
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    real(8),                intent(in)    :: Ac_tot(3)
    real(8),                intent(out)   :: current_raw(3)

    integer :: ispin, ifrag, i_local, ib, ilma, idir, istate, occ0, nbatch
    integer :: nocc_spin, nbf, gid, n_global, n_needed, pos, state_work, max_nbf_work
    integer :: occ_first, occ_last, local_idx, local_col
    integer, allocatable :: needed_ids(:), needed_pos(:)
    logical, allocatable :: row_needed(:), row_output(:)
    complex(8), allocatable :: coef_work(:,:), uV_state(:,:), coef_block(:,:), r_weight(:,:), rV_state(:,:)
    real(8) :: curr_local(3), curr_sum(3), occ
    real(8) :: perf_t0
    integer(8) :: target_bytes, bytes_per_state
    integer, parameter :: current_coef_cache_target_mb = 64
    logical :: perf_count

    current_raw(:) = 0.0d0
    curr_local(:) = 0.0d0
    curr_sum(:) = 0.0d0
    if (ppg%Nlma <= 0) return
    if (.not. allocated(ppg%uV)) return
    if (.not. allocated(dg_frag%coef_owner)) return
    if (.not. allocated(dg_frag%coef_global_to_local)) return
    if (.not. allocated(dg_frag%index_basis)) return
    if (.not. allocated(dg_frag%n_basis)) return

    perf_count = dg_frag%mixed_z_perf_count_enabled
    if (perf_count) perf_t0 = get_wtime()
    call ensure_nonlocal_projector_overlap_cache(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot)
    if (perf_count) then
      dg_frag%mixed_z_perf_wall_current_nl_cache = &
        dg_frag%mixed_z_perf_wall_current_nl_cache + (get_wtime() - perf_t0)
    end if
    if (.not. allocated(dg_frag%nl_projector_overlap)) return
    if (.not. allocated(dg_frag%nl_projector_r_overlap)) return

    if (perf_count) perf_t0 = get_wtime()
    n_global = size(dg_frag%coef_owner, 1)
    allocate(row_needed(n_global), row_output(n_global), needed_pos(n_global))
    if (perf_count) then
      dg_frag%mixed_z_perf_wall_current_nl_setup = &
        dg_frag%mixed_z_perf_wall_current_nl_setup + (get_wtime() - perf_t0)
    end if

    do ispin = 1, min(dg_frag%nspin, system%nspin)
      nocc_spin = 0
      if (allocated(dg_frag%nocc_spin)) nocc_spin = dg_frag%nocc_spin(ispin)
      if (dg_frag%coef_state_block_mode) then
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(system%rocc, 1))
      else
        nocc_spin = min(nocc_spin, dg_frag%nstate_tot, size(dg_frag%coef, 2), size(system%rocc, 1))
      end if
      if (nocc_spin <= 0) cycle
      occ_first = 1
      occ_last = nocc_spin
      if (dg_frag%coef_state_block_mode) then
        occ_first = max(1, dg_frag%coef_state_start)
        occ_last = min(nocc_spin, dg_frag%coef_state_end)
        if (occ_first > occ_last) cycle
      end if

      if (perf_count) perf_t0 = get_wtime()
      row_needed(:) = .false.
      row_output(:) = .false.
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        nbf = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%nl_projector_overlap, 1))
        do ib = 1, nbf
          gid = dg_frag%index_basis(ib, ifrag, ispin)
          if (gid < 1 .or. gid > n_global) cycle
          row_needed(gid) = .true.
          if (dg_frag%coef_state_block_mode) then
            if (dg_frag%coef_global_to_local(gid, ispin) > 0) row_output(gid) = .true.
          else
            if (dg_frag%coef_owner(gid, ispin) == dg_frag%id) row_output(gid) = .true.
          end if
        end do
      end do

      n_needed = count(row_needed)
      if (n_needed <= 0) cycle
      if (allocated(needed_ids)) deallocate(needed_ids)
      allocate(needed_ids(n_needed))
      needed_pos(:) = 0
      n_needed = 0
      do gid = 1, n_global
        if (.not. row_needed(gid)) cycle
        n_needed = n_needed + 1
        needed_ids(n_needed) = gid
        needed_pos(gid) = n_needed
      end do

      target_bytes = int(current_coef_cache_target_mb, kind=8) * 1024_8 * 1024_8
      bytes_per_state = 16_8 * int(max(1, n_needed + 4 * ppg%Nlma), kind=8)
      state_work = max(1, min(occ_last - occ_first + 1, int(max(1_8, target_bytes / max(1_8, bytes_per_state)))))
      max_nbf_work = 1
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        max_nbf_work = max(max_nbf_work, &
          min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%nl_projector_overlap, 1)))
      end do
      if (allocated(coef_work)) deallocate(coef_work, uV_state, coef_block, r_weight, rV_state)
      allocate(coef_work(n_needed, state_work))
      allocate(uV_state(ppg%Nlma, state_work))
      allocate(coef_block(max_nbf_work, state_work))
      allocate(r_weight(max_nbf_work, ppg%Nlma))
      allocate(rV_state(max_nbf_work, state_work))
      if (perf_count) then
        dg_frag%mixed_z_perf_wall_current_nl_setup = &
          dg_frag%mixed_z_perf_wall_current_nl_setup + (get_wtime() - perf_t0)
      end if

      do occ0 = occ_first, occ_last, state_work
        nbatch = min(state_work, occ_last - occ0 + 1)
        coef_work(:, :) = (0.0d0, 0.0d0)
        uV_state(:, :) = (0.0d0, 0.0d0)
        if (dg_frag%coef_state_block_mode) then
          if (perf_count) perf_t0 = get_wtime()
          do pos = 1, n_needed
            local_idx = dg_frag%coef_global_to_local(needed_ids(pos), ispin)
            if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
            do istate = 1, nbatch
              local_col = occ0 + istate - 1 - dg_frag%coef_state_start + 1
              if (local_col < 1 .or. local_col > size(dg_frag%coef, 2)) cycle
              coef_work(pos, istate) = dg_frag%coef(local_idx, local_col, ispin)
            end do
          end do
          if (perf_count) then
            dg_frag%mixed_z_perf_wall_current_nl_fetch = &
              dg_frag%mixed_z_perf_wall_current_nl_fetch + (get_wtime() - perf_t0)
          end if
        else
          if (perf_count) perf_t0 = get_wtime()
          call fetch_remote_coef_rows(dg_frag, ispin, needed_ids, coef_work(1:n_needed, 1:nbatch), &
                                      occ0, occ0 + nbatch - 1)
          if (perf_count) then
            dg_frag%mixed_z_perf_wall_current_nl_fetch = &
              dg_frag%mixed_z_perf_wall_current_nl_fetch + (get_wtime() - perf_t0)
          end if
        end if

        if (perf_count) perf_t0 = get_wtime()
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = ifrag - dg_frag%ifrag_start + 1
          if (i_local < 1 .or. i_local > size(dg_frag%nl_projector_overlap, 4)) cycle
          nbf = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%nl_projector_overlap, 1))
          if (nbf <= 0) cycle
          coef_block(1:nbf,1:nbatch) = (0.0d0, 0.0d0)
          do ib = 1, nbf
            gid = dg_frag%index_basis(ib, ifrag, ispin)
            if (gid < 1 .or. gid > n_global) cycle
            pos = needed_pos(gid)
            if (pos <= 0) cycle
            coef_block(ib,1:nbatch) = coef_work(pos,1:nbatch)
          end do
          do ilma = 1, ppg%Nlma
            do ib = 1, nbf
              r_weight(ib, ilma) = dg_frag%nl_projector_overlap(ib, ilma, ispin, i_local)
            end do
          end do
          call zgemm('T', 'N', ppg%Nlma, nbatch, nbf, (1.0d0, 0.0d0), &
            r_weight, max(1,max_nbf_work), &
            coef_block, max(1,max_nbf_work), (1.0d0, 0.0d0), &
            uV_state, max(1,ppg%Nlma))
          if (perf_count) dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 1_8
        end do
        if (perf_count) then
          dg_frag%mixed_z_perf_wall_current_nl_project = &
            dg_frag%mixed_z_perf_wall_current_nl_project + (get_wtime() - perf_t0)
        end if

        if (perf_count) perf_t0 = get_wtime()
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = ifrag - dg_frag%ifrag_start + 1
          if (i_local < 1 .or. i_local > size(dg_frag%nl_projector_overlap, 4)) cycle
          nbf = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1), size(dg_frag%nl_projector_overlap, 1))
          if (nbf <= 0) cycle
          do idir = 1, 3
            do ilma = 1, ppg%Nlma
              do ib = 1, nbf
                r_weight(ib, ilma) = conjg(dg_frag%nl_projector_r_overlap(idir, ib, ilma, ispin, i_local)) * &
                  cmplx(ppg%rinv_uvu(ilma), 0.0d0, kind=8)
              end do
            end do
            call zgemm('N', 'N', nbf, nbatch, ppg%Nlma, (1.0d0, 0.0d0), &
              r_weight, max(1,max_nbf_work), uV_state, max(1,ppg%Nlma), &
              (0.0d0, 0.0d0), rV_state, max(1,max_nbf_work))
            if (perf_count) dg_frag%mixed_z_perf_zgemm_calls = dg_frag%mixed_z_perf_zgemm_calls + 1_8
            do ib = 1, nbf
              gid = dg_frag%index_basis(ib, ifrag, ispin)
              if (gid < 1 .or. gid > n_global) cycle
              if (.not. row_output(gid)) cycle
              pos = needed_pos(gid)
              if (pos <= 0) cycle
              do istate = 1, nbatch
                occ = system%rocc(occ0 + istate - 1, 1, ispin)
                if (occ <= 0.0d0) cycle
                curr_local(idir) = curr_local(idir) + occ * 2.0d0 * &
                  aimag(conjg(coef_work(pos, istate)) * rV_state(ib, istate))
              end do
            end do
          end do
        end do
        if (perf_count) then
          dg_frag%mixed_z_perf_wall_current_nl_contract = &
            dg_frag%mixed_z_perf_wall_current_nl_contract + (get_wtime() - perf_t0)
        end if
      end do
    end do

    if (perf_count) perf_t0 = get_wtime()
    call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
    if (perf_count) then
      dg_frag%mixed_z_perf_wall_current_nl_reduce = &
        dg_frag%mixed_z_perf_wall_current_nl_reduce + (get_wtime() - perf_t0)
    end if
    if (curr_sum(1) /= curr_sum(1) .or. curr_sum(2) /= curr_sum(2) .or. curr_sum(3) /= curr_sum(3)) then
      stop "DG-Fragment RT: NaN in nonlocal macroscopic current reduction"
    end if
    current_raw(:) = curr_sum(:)
    if (allocated(coef_work)) deallocate(coef_work, uV_state, coef_block, r_weight, rV_state)
    if (allocated(needed_ids)) deallocate(needed_ids)
    deallocate(row_needed, row_output, needed_pos)
  end subroutine calculate_nonlocal_current_dg


  subroutine calculate_macroscopic_current_dg(dg_frag, system, mg, stencil, current_raw)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    real(8),                intent(out)   :: current_raw(3)

    integer :: ifrag, i_local, ispin, istate_frag, jstate_frag
    integer :: ig_i, ig_j
    integer :: ix, iy, iz, igrid, valid_grid_count, point_idx
    integer :: nxyz(3), ifrag_count, ngrid_max, max_nbf, max_nocc
    integer :: nocc_spin
    integer(8) :: grid_checksum
    integer :: phi_lb1, phi_ub1, phi_lb2, phi_ub2, phi_lb3, phi_ub3
    integer :: grad_lb1, grad_ub1, grad_lb2, grad_ub2, grad_lb3, grad_ub3
    integer :: bad_nan
    real(8) :: occ_factor, current_tmp
    real(8) :: phi_i, grad_x, grad_y, grad_z
    real(8) :: curr_local(3), curr_sum(3)
    real(8) :: curr_x, curr_y, curr_z
	    complex(8), allocatable :: coef_occ_frag(:,:), coef_pair_mat(:,:)
	    integer, allocatable :: row_ids(:)
	    logical :: need_current_grid_cache

    current_raw(:) = 0.0d0
    curr_local(:) = 0.0d0
    curr_sum(:) = 0.0d0
    if (allocated(dg_frag%momentum_blocks)) then
      call calculate_macroscopic_current_from_momentum_blocks(dg_frag, system, current_raw)
      return
    end if
    write(*,'(1x,a)') '[FATAL] DG current requires block-sparse momentum blocks.'
    write(*,'(1x,a)') '[FATAL] The old real-space dense current fallback is disabled.'
    stop 'DG-Fragment RT: missing block-sparse momentum current operator'
    if (.not. allocated(dg_frag%phi_frag)) then
      call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
      current_raw(:) = curr_sum(:)
      return
    end if

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) then
      call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
      current_raw(:) = curr_sum(:)
      return
    end if

    call ensure_gradient_basis_cache(dg_frag, mg, stencil)
    if (.not. allocated(dg_frag%gradient_basis_cache)) then
      call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
      current_raw(:) = curr_sum(:)
      return
    end if

    occ_factor = 2.0d0 / real(system%nspin, 8)
    phi_lb1 = lbound(dg_frag%phi_frag, 1); phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2); phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3); phi_ub3 = ubound(dg_frag%phi_frag, 3)
    grad_lb1 = lbound(dg_frag%gradient_basis_cache, 1); grad_ub1 = ubound(dg_frag%gradient_basis_cache, 1)
    grad_lb2 = lbound(dg_frag%gradient_basis_cache, 2); grad_ub2 = ubound(dg_frag%gradient_basis_cache, 2)
    grad_lb3 = lbound(dg_frag%gradient_basis_cache, 3); grad_ub3 = ubound(dg_frag%gradient_basis_cache, 3)

    max_nbf = max(1, maxval(dg_frag%n_basis(dg_frag%ifrag_start:dg_frag%ifrag_end, 1:system%nspin)))
    max_nocc = max(1, maxval(dg_frag%nocc_spin(1:system%nspin)))
    allocate(coef_pair_mat(max_nbf, max_nbf))
    allocate(coef_occ_frag(max_nbf, max_nocc), row_ids(max_nbf))

    ngrid_max = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      ngrid_max = max(ngrid_max, product(dg_frag%nxyz_domain(1:3, ifrag)))
    end do
    if (ngrid_max <= 0) then
      if (allocated(coef_occ_frag)) deallocate(coef_occ_frag)
      if (allocated(row_ids)) deallocate(row_ids)
      deallocate(coef_pair_mat)
      call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
      current_raw(:) = curr_sum(:)
      return
    end if
    need_current_grid_cache = .false.
    if (.not. allocated(dg_frag%current_valid_grid_count)) then
      need_current_grid_cache = .true.
    else if (size(dg_frag%current_valid_grid_count) /= ifrag_count) then
      need_current_grid_cache = .true.
    else if (.not. allocated(dg_frag%current_valid_ix) .or. .not. allocated(dg_frag%current_valid_iy) .or. &
             .not. allocated(dg_frag%current_valid_iz) .or. .not. allocated(dg_frag%current_valid_ixg) .or. &
             .not. allocated(dg_frag%current_valid_iyg) .or. .not. allocated(dg_frag%current_valid_izg)) then
      need_current_grid_cache = .true.
    else if (.not. allocated(dg_frag%current_density_grid_point_count) .or. &
             .not. allocated(dg_frag%current_density_grid_checksum)) then
      need_current_grid_cache = .true.
    else if (size(dg_frag%current_valid_ix, 1) < ngrid_max .or. size(dg_frag%current_valid_ix, 2) /= ifrag_count) then
      need_current_grid_cache = .true.
    else if (size(dg_frag%current_density_grid_point_count) /= ifrag_count .or. &
             size(dg_frag%current_density_grid_checksum) /= ifrag_count) then
      need_current_grid_cache = .true.
    else
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        grid_checksum = int(dg_frag%density_grid_point_count(i_local), 8)
        do point_idx = 1, dg_frag%density_grid_point_count(i_local)
          grid_checksum = grid_checksum + int(point_idx, 8) * 1000003_8 + &
                          int(dg_frag%density_grid_points(point_idx, i_local)%ixg, 8) * 9176_8 + &
                          int(dg_frag%density_grid_points(point_idx, i_local)%iyg, 8) * 131_8 + &
                          int(dg_frag%density_grid_points(point_idx, i_local)%izg, 8)
        end do
        if (dg_frag%current_density_grid_point_count(i_local) /= dg_frag%density_grid_point_count(i_local) .or. &
            dg_frag%current_density_grid_checksum(i_local) /= grid_checksum) then
          need_current_grid_cache = .true.
          exit
        end if
      end do
    end if
    if (need_current_grid_cache) then
      if (allocated(dg_frag%current_valid_grid_count)) deallocate(dg_frag%current_valid_grid_count)
      if (allocated(dg_frag%current_density_grid_point_count)) deallocate(dg_frag%current_density_grid_point_count)
      if (allocated(dg_frag%current_density_grid_checksum)) deallocate(dg_frag%current_density_grid_checksum)
      if (allocated(dg_frag%current_valid_ix)) deallocate(dg_frag%current_valid_ix)
      if (allocated(dg_frag%current_valid_iy)) deallocate(dg_frag%current_valid_iy)
      if (allocated(dg_frag%current_valid_iz)) deallocate(dg_frag%current_valid_iz)
      if (allocated(dg_frag%current_valid_ixg)) deallocate(dg_frag%current_valid_ixg)
      if (allocated(dg_frag%current_valid_iyg)) deallocate(dg_frag%current_valid_iyg)
      if (allocated(dg_frag%current_valid_izg)) deallocate(dg_frag%current_valid_izg)
      allocate(dg_frag%current_valid_grid_count(ifrag_count))
      allocate(dg_frag%current_density_grid_point_count(ifrag_count))
      allocate(dg_frag%current_density_grid_checksum(ifrag_count))
      allocate(dg_frag%current_valid_ix(ngrid_max, ifrag_count), dg_frag%current_valid_iy(ngrid_max, ifrag_count), &
               dg_frag%current_valid_iz(ngrid_max, ifrag_count))
      allocate(dg_frag%current_valid_ixg(ngrid_max, ifrag_count), dg_frag%current_valid_iyg(ngrid_max, ifrag_count), &
               dg_frag%current_valid_izg(ngrid_max, ifrag_count))
      dg_frag%current_valid_grid_count = 0
      dg_frag%current_density_grid_point_count = 0
      dg_frag%current_density_grid_checksum = 0_8
      dg_frag%current_valid_ix = 0
      dg_frag%current_valid_iy = 0
      dg_frag%current_valid_iz = 0
      dg_frag%current_valid_ixg = 0
      dg_frag%current_valid_iyg = 0
      dg_frag%current_valid_izg = 0

      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        valid_grid_count = 0
        grid_checksum = int(dg_frag%density_grid_point_count(i_local), 8)
        do point_idx = 1, dg_frag%density_grid_point_count(i_local)
          grid_checksum = grid_checksum + int(point_idx, 8) * 1000003_8 + &
                          int(dg_frag%density_grid_points(point_idx, i_local)%ixg, 8) * 9176_8 + &
                          int(dg_frag%density_grid_points(point_idx, i_local)%iyg, 8) * 131_8 + &
                          int(dg_frag%density_grid_points(point_idx, i_local)%izg, 8)
          if (dg_frag%density_grid_points(point_idx, i_local)%ixg < mg%is(1) .or. &
              dg_frag%density_grid_points(point_idx, i_local)%ixg > mg%ie(1)) cycle
          if (dg_frag%density_grid_points(point_idx, i_local)%iyg < mg%is(2) .or. &
              dg_frag%density_grid_points(point_idx, i_local)%iyg > mg%ie(2)) cycle
          if (dg_frag%density_grid_points(point_idx, i_local)%izg < mg%is(3) .or. &
              dg_frag%density_grid_points(point_idx, i_local)%izg > mg%ie(3)) cycle
          ix = dg_frag%density_grid_points(point_idx, i_local)%ix
          iy = dg_frag%density_grid_points(point_idx, i_local)%iy
          iz = dg_frag%density_grid_points(point_idx, i_local)%iz
          if (ix < dg_frag%frag_core_lo(1, ifrag) .or. ix > dg_frag%frag_core_hi(1, ifrag)) cycle
          if (iy < dg_frag%frag_core_lo(2, ifrag) .or. iy > dg_frag%frag_core_hi(2, ifrag)) cycle
          if (iz < dg_frag%frag_core_lo(3, ifrag) .or. iz > dg_frag%frag_core_hi(3, ifrag)) cycle
          if (ix < phi_lb1 .or. ix > phi_ub1 .or. ix < grad_lb1 .or. ix > grad_ub1) cycle
          if (iy < phi_lb2 .or. iy > phi_ub2 .or. iy < grad_lb2 .or. iy > grad_ub2) cycle
          if (iz < phi_lb3 .or. iz > phi_ub3 .or. iz < grad_lb3 .or. iz > grad_ub3) cycle
          valid_grid_count = valid_grid_count + 1
          dg_frag%current_valid_ix(valid_grid_count, i_local) = ix
          dg_frag%current_valid_iy(valid_grid_count, i_local) = iy
          dg_frag%current_valid_iz(valid_grid_count, i_local) = iz
          dg_frag%current_valid_ixg(valid_grid_count, i_local) = dg_frag%density_grid_points(point_idx, i_local)%ixg
          dg_frag%current_valid_iyg(valid_grid_count, i_local) = dg_frag%density_grid_points(point_idx, i_local)%iyg
          dg_frag%current_valid_izg(valid_grid_count, i_local) = dg_frag%density_grid_points(point_idx, i_local)%izg
	        end do
	        dg_frag%current_valid_grid_count(i_local) = valid_grid_count
	        dg_frag%current_density_grid_point_count(i_local) = dg_frag%density_grid_point_count(i_local)
	        dg_frag%current_density_grid_checksum(i_local) = grid_checksum
	      end do
	    end if

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
      valid_grid_count = dg_frag%current_valid_grid_count(i_local)
      if (valid_grid_count <= 0) cycle

      do ispin = 1, system%nspin
        nocc_spin = dg_frag%nocc_spin(ispin)
        if (nocc_spin <= 0) cycle
        coef_pair_mat(:, :) = (0.0d0, 0.0d0)
        coef_occ_frag(:, :) = (0.0d0, 0.0d0)
        do jstate_frag = 1, dg_frag%n_basis(ifrag, ispin)
          row_ids(jstate_frag) = dg_frag%index_basis(jstate_frag, ifrag, ispin)
        end do
        call fetch_remote_coef_rows(dg_frag, ispin, row_ids(1:dg_frag%n_basis(ifrag, ispin)), &
                                    coef_occ_frag(1:dg_frag%n_basis(ifrag, ispin), 1:nocc_spin), &
                                    1, nocc_spin)
        call zgemm('N', 'C', dg_frag%n_basis(ifrag, ispin), dg_frag%n_basis(ifrag, ispin), nocc_spin, &
                   cmplx(occ_factor, 0.0d0, kind=8), coef_occ_frag, max_nbf, coef_occ_frag, max_nbf, &
                   cmplx(0.0d0, 0.0d0, kind=8), coef_pair_mat, max_nbf)

        curr_x = 0.0d0
        curr_y = 0.0d0
        curr_z = 0.0d0
        bad_nan = 0
!$omp parallel do collapse(2) private(jstate_frag, istate_frag, ig_j, ig_i, igrid, ix, iy, iz) &
!$omp& private(phi_i, current_tmp, grad_x, grad_y, grad_z) &
!$omp& reduction(+:curr_x,curr_y,curr_z) reduction(max:bad_nan) schedule(static)
        do jstate_frag = 1, dg_frag%n_basis(ifrag, ispin)
          do istate_frag = 1, dg_frag%n_basis(ifrag, ispin)
            ig_j = dg_frag%index_basis(jstate_frag, ifrag, ispin)
            ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            current_tmp = -aimag(coef_pair_mat(istate_frag, jstate_frag))
            if (current_tmp /= current_tmp) then
              bad_nan = 1
              cycle
            end if
            if (abs(current_tmp) < 1.0d-18) cycle
            do igrid = 1, valid_grid_count
              ix = dg_frag%current_valid_ix(igrid, i_local)
              iy = dg_frag%current_valid_iy(igrid, i_local)
              iz = dg_frag%current_valid_iz(igrid, i_local)
              phi_i = dg_frag%phi_frag(ix, iy, iz, istate_frag, i_local)
              grad_x = dg_frag%gradient_basis_cache(ix, iy, iz, 1, jstate_frag, i_local)
              grad_y = dg_frag%gradient_basis_cache(ix, iy, iz, 2, jstate_frag, i_local)
              grad_z = dg_frag%gradient_basis_cache(ix, iy, iz, 3, jstate_frag, i_local)
              if (phi_i /= phi_i .or. grad_x /= grad_x .or. grad_y /= grad_y .or. grad_z /= grad_z) then
                bad_nan = 1
                cycle
              end if
              curr_x = curr_x + current_tmp * phi_i * grad_x
              curr_y = curr_y + current_tmp * phi_i * grad_y
              curr_z = curr_z + current_tmp * phi_i * grad_z
            end do
          end do
        end do
!$omp end parallel do
        if (bad_nan /= 0) stop "DG-Fragment RT: NaN in macroscopic current inputs"
        curr_local(1) = curr_local(1) + curr_x
        curr_local(2) = curr_local(2) + curr_y
        curr_local(3) = curr_local(3) + curr_z
      end do
    end do

    curr_local(:) = curr_local(:) / dble(max(1, dg_frag%isize_frag))
    call comm_summation(curr_local, curr_sum, 3, dg_frag%icomm)
    if (curr_sum(1) /= curr_sum(1) .or. curr_sum(2) /= curr_sum(2) .or. curr_sum(3) /= curr_sum(3)) then
      stop "DG-Fragment RT: NaN in macroscopic current reduction"
    end if
    current_raw(:) = curr_sum(:)

    if (allocated(coef_occ_frag)) deallocate(coef_occ_frag)
    if (allocated(row_ids)) deallocate(row_ids)
    deallocate(coef_pair_mat)
  end subroutine calculate_macroscopic_current_dg

  subroutine calculate_microscopic_current_dg(dg_frag, system, mg, stencil, curr)
    ! NOTE: DG microscopic current here intentionally excludes non-local PP contribution.
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_vector),         intent(inout) :: curr

    integer :: ifrag, i_local, ispin, io, istate_frag, jstate_frag
    integer :: ig_i, ig_j
    integer :: ix, iy, iz, ixg, iyg, izg, igrid, valid_grid_count, point_idx
    integer :: grid_x_lo, grid_x_hi, grid_y_lo, grid_y_hi, grid_z_lo, grid_z_hi
    integer :: nxyz(3), ifrag_count, ngrid_max, max_nbf, max_nocc
    integer :: nocc_spin
    real(8) :: occ_factor, current_tmp
    real(8) :: phi_i
    real(8), allocatable :: curr_thread(:,:,:,:)
    real(8), allocatable :: curr_local(:,:,:,:), curr_sum(:,:,:,:)
    real(8), allocatable :: w_local(:,:,:), w_sum(:,:,:)
    complex(8), allocatable :: coef_occ_frag(:,:), coef_pair_mat(:,:)
    integer, allocatable :: row_ids(:)

    curr%v = 0.0d0
    if (.not. allocated(dg_frag%phi_frag)) return

    occ_factor = 2.0d0 / real(system%nspin, 8)

    allocate(curr_local(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(curr_sum(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_sum(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    curr_local = 0.0d0
    curr_sum = 0.0d0
    w_local = 0.0d0
    w_sum = 0.0d0
    call ensure_gradient_basis_cache(dg_frag, mg, stencil)
    max_nbf = max(1, maxval(dg_frag%n_basis(dg_frag%ifrag_start:dg_frag%ifrag_end, 1:system%nspin)))
    max_nocc = max(1, maxval(dg_frag%nocc_spin(1:system%nspin)))
    allocate(coef_occ_frag(max_nbf, max_nocc), coef_pair_mat(max_nbf, max_nbf))
    allocate(row_ids(max_nbf))
    grid_x_lo = mg%is(1)
    grid_x_hi = mg%ie(1)
    grid_y_lo = mg%is(2)
    grid_y_hi = mg%ie(2)
    grid_z_lo = mg%is(3)
    grid_z_hi = mg%ie(3)

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    ngrid_max = 0
    if (ifrag_count > 0) then
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        ngrid_max = max(ngrid_max, product(dg_frag%nxyz_domain(1:3, ifrag)))
      end do
    end if
    if (.not. allocated(dg_frag%current_valid_grid_count)) then
      if (ifrag_count > 0 .and. ngrid_max > 0) then
        allocate(dg_frag%current_valid_grid_count(ifrag_count))
        allocate(dg_frag%current_valid_ix(ngrid_max, ifrag_count), dg_frag%current_valid_iy(ngrid_max, ifrag_count), &
                 dg_frag%current_valid_iz(ngrid_max, ifrag_count))
        allocate(dg_frag%current_valid_ixg(ngrid_max, ifrag_count), dg_frag%current_valid_iyg(ngrid_max, ifrag_count), &
                 dg_frag%current_valid_izg(ngrid_max, ifrag_count))
        dg_frag%current_valid_grid_count = 0
        dg_frag%current_valid_ix = 0
        dg_frag%current_valid_iy = 0
        dg_frag%current_valid_iz = 0
        dg_frag%current_valid_ixg = 0
        dg_frag%current_valid_iyg = 0
        dg_frag%current_valid_izg = 0

        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          valid_grid_count = 0
          do point_idx = 1, dg_frag%density_grid_point_count(i_local)
            ixg = dg_frag%density_grid_points(point_idx, i_local)%ixg
            iyg = dg_frag%density_grid_points(point_idx, i_local)%iyg
            izg = dg_frag%density_grid_points(point_idx, i_local)%izg
            if (ixg < grid_x_lo .or. ixg > grid_x_hi) cycle
            if (iyg < grid_y_lo .or. iyg > grid_y_hi) cycle
            if (izg < grid_z_lo .or. izg > grid_z_hi) cycle
            if (dg_frag%density_grid_points(point_idx, i_local)%ix < dg_frag%frag_core_lo(1, ifrag) .or. &
                dg_frag%density_grid_points(point_idx, i_local)%ix > dg_frag%frag_core_hi(1, ifrag)) cycle
            if (dg_frag%density_grid_points(point_idx, i_local)%iy < dg_frag%frag_core_lo(2, ifrag) .or. &
                dg_frag%density_grid_points(point_idx, i_local)%iy > dg_frag%frag_core_hi(2, ifrag)) cycle
            if (dg_frag%density_grid_points(point_idx, i_local)%iz < dg_frag%frag_core_lo(3, ifrag) .or. &
                dg_frag%density_grid_points(point_idx, i_local)%iz > dg_frag%frag_core_hi(3, ifrag)) cycle
            valid_grid_count = valid_grid_count + 1
            dg_frag%current_valid_ix(valid_grid_count, i_local) = dg_frag%density_grid_points(point_idx, i_local)%ix
            dg_frag%current_valid_iy(valid_grid_count, i_local) = dg_frag%density_grid_points(point_idx, i_local)%iy
            dg_frag%current_valid_iz(valid_grid_count, i_local) = dg_frag%density_grid_points(point_idx, i_local)%iz
            dg_frag%current_valid_ixg(valid_grid_count, i_local) = ixg
            dg_frag%current_valid_iyg(valid_grid_count, i_local) = iyg
            dg_frag%current_valid_izg(valid_grid_count, i_local) = izg
          end do
          dg_frag%current_valid_grid_count(i_local) = valid_grid_count
        end do
      end if
    end if

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
      valid_grid_count = dg_frag%current_valid_grid_count(i_local)
      do igrid = 1, valid_grid_count
        ixg = dg_frag%current_valid_ixg(igrid, i_local)
        iyg = dg_frag%current_valid_iyg(igrid, i_local)
        izg = dg_frag%current_valid_izg(igrid, i_local)
        w_local(ixg, iyg, izg) = w_local(ixg, iyg, izg) + 1.0d0
      end do

      do ispin = 1, system%nspin
        nocc_spin = dg_frag%nocc_spin(ispin)
        if (nocc_spin <= 0) cycle
        coef_occ_frag(:, :) = (0.0d0, 0.0d0)
        coef_pair_mat(:, :) = (0.0d0, 0.0d0)
        do jstate_frag = 1, dg_frag%n_basis(ifrag, ispin)
          ig_j = dg_frag%index_basis(jstate_frag, ifrag, ispin)
          row_ids(jstate_frag) = ig_j
        end do
        call fetch_remote_coef_rows(dg_frag, ispin, row_ids(1:dg_frag%n_basis(ifrag, ispin)), &
                                    coef_occ_frag(1:dg_frag%n_basis(ifrag, ispin), 1:nocc_spin), &
                                    1, nocc_spin)
        call zgemm('N', 'C', dg_frag%n_basis(ifrag, ispin), dg_frag%n_basis(ifrag, ispin), nocc_spin, &
                   cmplx(occ_factor, 0.0d0, kind=8), coef_occ_frag, max_nbf, coef_occ_frag, max_nbf, &
                   cmplx(0.0d0, 0.0d0, kind=8), coef_pair_mat, max_nbf)
!$omp parallel private(jstate_frag, istate_frag, ig_j, ig_i, phi_i, iz, iy, ix, ixg, iyg, izg, curr_thread, current_tmp)
          allocate(curr_thread(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          curr_thread = 0.0d0

!$omp do schedule(static)
          do jstate_frag = 1, dg_frag%n_basis(ifrag, ispin)
            ig_j = dg_frag%index_basis(jstate_frag, ifrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

            do istate_frag = 1, dg_frag%n_basis(ifrag, ispin)
              ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

              current_tmp = aimag(coef_pair_mat(istate_frag, jstate_frag))
              if (abs(current_tmp) < 1.0d-18) cycle

              do igrid = 1, valid_grid_count
                ix = dg_frag%current_valid_ix(igrid, i_local)
                iy = dg_frag%current_valid_iy(igrid, i_local)
                iz = dg_frag%current_valid_iz(igrid, i_local)
                ixg = dg_frag%current_valid_ixg(igrid, i_local)
                iyg = dg_frag%current_valid_iyg(igrid, i_local)
                izg = dg_frag%current_valid_izg(igrid, i_local)

                phi_i = dg_frag%phi_frag(ix, iy, iz, istate_frag, i_local)
                curr_thread(1, ixg, iyg, izg) = curr_thread(1, ixg, iyg, izg) + &
                  current_tmp * phi_i * dg_frag%gradient_basis_cache(ix, iy, iz, 1, jstate_frag, i_local)
                curr_thread(2, ixg, iyg, izg) = curr_thread(2, ixg, iyg, izg) + &
                  current_tmp * phi_i * dg_frag%gradient_basis_cache(ix, iy, iz, 2, jstate_frag, i_local)
                curr_thread(3, ixg, iyg, izg) = curr_thread(3, ixg, iyg, izg) + &
                  current_tmp * phi_i * dg_frag%gradient_basis_cache(ix, iy, iz, 3, jstate_frag, i_local)
              end do
            end do
          end do
!$omp end do

!$omp critical
          curr_local(:, :, :, :) = curr_local(:, :, :, :) + curr_thread(:, :, :, :)
!$omp end critical

          deallocate(curr_thread)
!$omp end parallel
      end do
    end do

    call comm_summation(curr_local, curr_sum, size(curr_local), dg_frag%icomm)
    call comm_summation(w_local, w_sum, size(w_local), dg_frag%icomm)

    do izg = mg%is(3), mg%ie(3)
      do iyg = mg%is(2), mg%ie(2)
        do ixg = mg%is(1), mg%ie(1)
          if (w_sum(ixg, iyg, izg) > 1.0d-12) then
            curr%v(:, ixg, iyg, izg) = curr_sum(:, ixg, iyg, izg) / w_sum(ixg, iyg, izg)
          else
            curr%v(:, ixg, iyg, izg) = curr_sum(:, ixg, iyg, izg)
          end if
        end do
      end do
    end do

    deallocate(curr_local, curr_sum, w_local, w_sum)
    deallocate(coef_occ_frag, coef_pair_mat, row_ids)
  end subroutine calculate_microscopic_current_dg

end module rt_dg_fragment_ops
