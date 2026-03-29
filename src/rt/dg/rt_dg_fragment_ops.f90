module rt_dg_fragment_ops
  use communication, only: comm_bcast, comm_get_max, comm_is_root, comm_summation, COMM_GROUP_NULL
  use rt_dg_fragment_types, only: s_dg_fragment_rt, matrix_block_info, complex_matrix_block_info
  implicit none

  private
  public :: ensure_nonlocal_pp_matrix_A
  public :: ensure_overlap_prop_available
  public :: calculate_microscopic_current_dg
  public :: build_spatial_A_coupling_matrices
  public :: apply_gradient_to_basis
  public :: apply_momentum_blocks
  public :: apply_matrix_blocks
  public :: apply_matrix_blocks_batch
  public :: apply_complex_matrix_blocks_batch
  public :: apply_mixed_hamiltonian
  public :: rebuild_local_h_block_ids
  public :: copy_matrix_blocks_to_complex_dense
  public :: copy_matrix_blocks_metric_to_complex_dense
  public :: copy_matrix_blocks_metric_to_real_dense
  public :: copy_hamiltonian_metric_to_complex_dense
  public :: copy_momentum_blocks_to_complex_dense
  public :: symmetrize_real_matrix_blocks
  public :: mixed_fp_coupling_active
  public :: apply_overlap_operator
  public :: solve_overlap_operator_batch
  public :: solve_overlap_operator_batch_local
  public :: copy_overlap_operator_to_dense
  public :: pack_owned_coef
  public :: fetch_remote_coef_rows
  public :: pack_owned_coef_pw
  public :: fetch_remote_coef_pw_rows
  public :: refresh_pw_coef_cache
  public :: gather_full_coef_view
  public :: zero_nonowned_coefficients
  public :: zero_nonlocal_h_matrix_blocks
  public :: sync_raw_coef_from_mixed
  public :: sync_mixed_coef_from_raw

contains

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

  subroutine sync_raw_coef_from_mixed(dg_frag, ispin)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin

    integer :: n_frag, n_pw, n_tot, n_basis, nstate
    complex(8), allocatable :: raw_all(:,:)

    if (.not. dg_frag%mixed_basis_ready) return
    if (.not. allocated(dg_frag%mixed_transform)) return
    if (.not. allocated(dg_frag%mixed_basis_dim)) return
    if (.not. allocated(dg_frag%coef_mix)) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return

    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    nstate = dg_frag%nstate_tot
    n_basis = min(dg_frag%mixed_basis_dim(ispin), size(dg_frag%mixed_transform, 2), size(dg_frag%coef_mix, 1))
    if (n_basis <= 0 .or. n_tot <= 0 .or. nstate <= 0) return

    allocate(raw_all(n_tot, nstate))
    raw_all(:, :) = matmul(dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin), dg_frag%coef_mix(1:n_basis, 1:nstate, ispin))
    dg_frag%coef(1:n_frag, 1:nstate, ispin) = raw_all(1:n_frag, 1:nstate)
    if (n_pw > 0) dg_frag%coef_pw(1:n_pw, 1:nstate, ispin) = raw_all(n_frag+1:n_tot, 1:nstate)
    deallocate(raw_all)
  end subroutine sync_raw_coef_from_mixed

  subroutine sync_mixed_coef_from_raw(dg_frag, ispin)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin

    integer :: n_frag, n_pw, n_tot, n_basis, nstate
    complex(8), allocatable :: raw_all(:,:), overlap_all(:,:), mixed_all(:,:)

    if (.not. dg_frag%mixed_basis_ready) return
    if (.not. allocated(dg_frag%mixed_transform)) return
    if (.not. allocated(dg_frag%mixed_basis_dim)) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return

    n_frag = dg_frag%n_mat_max
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    n_tot = n_frag + n_pw
    nstate = dg_frag%nstate_tot
    n_basis = min(dg_frag%mixed_basis_dim(ispin), size(dg_frag%mixed_transform, 2))
    if (n_basis <= 0 .or. n_tot <= 0 .or. nstate <= 0) return

    if (.not. allocated(dg_frag%coef_mix)) then
      allocate(dg_frag%coef_mix(n_tot, nstate, dg_frag%nspin))
      dg_frag%coef_mix(:, :, :) = (0.0d0, 0.0d0)
    else if (size(dg_frag%coef_mix, 1) /= n_tot .or. size(dg_frag%coef_mix, 2) /= nstate .or. &
             size(dg_frag%coef_mix, 3) /= dg_frag%nspin) then
      deallocate(dg_frag%coef_mix)
      allocate(dg_frag%coef_mix(n_tot, nstate, dg_frag%nspin))
      dg_frag%coef_mix(:, :, :) = (0.0d0, 0.0d0)
    end if

    allocate(raw_all(n_tot, nstate), overlap_all(n_tot, nstate), mixed_all(n_basis, nstate))
    raw_all(:, :) = (0.0d0, 0.0d0)
    raw_all(1:n_frag, :) = dg_frag%coef(1:n_frag, 1:nstate, ispin)
    if (n_pw > 0) raw_all(n_frag+1:n_tot, :) = dg_frag%coef_pw(1:n_pw, 1:nstate, ispin)

    call apply_overlap_operator_batch(dg_frag, ispin, raw_all, overlap_all, .false.)
    mixed_all(:, :) = matmul(conjg(transpose(dg_frag%mixed_transform(1:n_tot, 1:n_basis, ispin))), overlap_all)
    dg_frag%coef_mix(:, :, ispin) = (0.0d0, 0.0d0)
    dg_frag%coef_mix(1:n_basis, 1:nstate, ispin) = mixed_all(:, :)

    deallocate(raw_all, overlap_all, mixed_all)
  end subroutine sync_mixed_coef_from_raw

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

    if (allocated(dg_frag%runtime_neighbor_pair_cache)) return
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

  !=======================================================================
  ! Build non-local pseudopotential matrix with vector potential A(t)
  ! V_NL(A) = e^{-i A·r} V_NL e^{i A·r}
  !
  ! Uses the SALMON approximation: A is nearly constant within PP cutoff
  !=======================================================================
  subroutine build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, nspin, hvol, Ac_tot, use_micro_A, Ac_micro, H_nl)
    use structures
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    integer,                intent(in) :: nspin
    real(8),                intent(in) :: hvol
    real(8),                intent(in) :: Ac_tot(3)
    logical,                intent(in) :: use_micro_A
    real(8),                intent(in), optional :: Ac_micro(:, :, :, :)
    complex(8),             intent(out) :: H_nl(dg_frag%n_mat_max, dg_frag%n_mat_max, nspin)

    integer :: ifrag, ispin, io, jo, i_local, ilma, ia, j, ix, iy, iz, ig_i, ig_j, nbf
    integer :: iorg(3), ndom(3), lx, ly, lz
    integer :: is(3), ie(3)
    real(8) :: x, y, z, phase
    real(8) :: A_local(3)
    complex(8), allocatable :: uVpsi(:,:,:)  ! (nstate_frag, Nlma, nspin)
    complex(8) :: overlap_i, overlap_j, nlpp_contrib

    if (ppg%Nlma == 0) then
      H_nl = (0.0d0, 0.0d0)
      return
    end if

    is = mg%is
    ie = mg%ie
    H_nl = (0.0d0, 0.0d0)

    allocate(uVpsi(dg_frag%nstate_frag, ppg%Nlma, nspin))

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
      uVpsi = (0.0d0, 0.0d0)

      do ispin = 1, nspin
        !$omp parallel do collapse(2) private(ilma, io, ia, j, ix, iy, iz, x, y, z, phase, overlap_i)
        do ilma = 1, ppg%Nlma
          do io = 1, dg_frag%nstate_frag
            ia = ppg%ia_tbl(ilma)
            overlap_i = (0.0d0, 0.0d0)
            do j = 1, ppg%mps(ia)
              ix = ppg%jxyz(1, j, ia)
              iy = ppg%jxyz(2, j, ia)
              iz = ppg%jxyz(3, j, ia)

              if (ix >= is(1) .and. ix <= ie(1) .and. &
                  iy >= is(2) .and. iy <= ie(2) .and. &
                  iz >= is(3) .and. iz <= ie(3)) then
                lx = mod(ix - iorg(1), dg_frag%lgnum_total(1)) + 1
                ly = mod(iy - iorg(2), dg_frag%lgnum_total(2)) + 1
                lz = mod(iz - iorg(3), dg_frag%lgnum_total(3)) + 1
                if (lx < 1 .or. lx > ndom(1)) cycle
                if (ly < 1 .or. ly > ndom(2)) cycle
                if (lz < 1 .or. lz > ndom(3)) cycle
                x = ppg%rxyz(1, j, ia)
                y = ppg%rxyz(2, j, ia)
                z = ppg%rxyz(3, j, ia)
                if (use_micro_A .and. present(Ac_micro)) then
                  A_local(1:3) = Ac_micro(1:3, ix, iy, iz)
                else
                  A_local(1:3) = Ac_tot(1:3)
                end if
                phase = A_local(1) * x + A_local(2) * y + A_local(3) * z
                overlap_i = overlap_i + dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                            ppg%uV(j, ilma) * exp(-zi * phase) * hvol
              end if
            end do
            uVpsi(io, ilma, ispin) = overlap_i
          end do
        end do
        !$omp end parallel do

        nbf = dg_frag%n_basis(ifrag, ispin)
        !$omp parallel do collapse(2) private(jo, io, ig_i, ig_j, ilma, nlpp_contrib, overlap_i, overlap_j)
        do jo = 1, nbf
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            nlpp_contrib = (0.0d0, 0.0d0)
            do ilma = 1, ppg%Nlma
              overlap_i = uVpsi(io, ilma, ispin)
              overlap_j = uVpsi(jo, ilma, ispin)
              nlpp_contrib = nlpp_contrib + overlap_i * overlap_j * ppg%rinv_uvu(ilma)
            end do
            H_nl(ig_i, ig_j, ispin) = H_nl(ig_i, ig_j, ispin) + nlpp_contrib
          end do
        end do
        !$omp end parallel do
      end do
    end do

    deallocate(uVpsi)

  end subroutine build_nonlocal_pp_matrix_A

  subroutine build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, nspin, hvol, Ac_tot, use_micro_A, Ac_micro, H_nl_blocks, H_nl_block_map)
    use structures
    use math_constants, only: zi
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid),          intent(in) :: mg
    type(s_pp_grid),        intent(in) :: ppg
    integer,                intent(in) :: nspin
    real(8),                intent(in) :: hvol
    real(8),                intent(in) :: Ac_tot(3)
    logical,                intent(in) :: use_micro_A
    real(8),                intent(in), optional :: Ac_micro(:, :, :, :)
    type(complex_matrix_block_info), intent(inout) :: H_nl_blocks(:)
    integer,                intent(in) :: H_nl_block_map(:, :)

    integer :: ifrag, ispin, io, jo, i_local, ilma, ia, j, ix, iy, iz, nbf, iblk
    integer :: iorg(3), ndom(3), lx, ly, lz
    integer :: is(3), ie(3)
    real(8) :: x, y, z, phase
    real(8) :: A_local(3)
    complex(8), allocatable :: uVpsi(:,:,:)  ! (nstate_frag, Nlma, nspin)
    complex(8) :: overlap_i, overlap_j, nlpp_contrib

    if (ppg%Nlma == 0) return

    is = mg%is
    ie = mg%ie
    do iblk = 1, size(H_nl_blocks)
      H_nl_blocks(iblk)%val(:, :, :) = (0.0d0, 0.0d0)
    end do

    allocate(uVpsi(dg_frag%nstate_frag, ppg%Nlma, nspin))

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
      uVpsi = (0.0d0, 0.0d0)

      do ispin = 1, nspin
        !$omp parallel do collapse(2) private(ilma, io, ia, j, ix, iy, iz, x, y, z, phase, overlap_i)
        do ilma = 1, ppg%Nlma
          do io = 1, dg_frag%nstate_frag
            ia = ppg%ia_tbl(ilma)
            overlap_i = (0.0d0, 0.0d0)
            do j = 1, ppg%mps(ia)
              ix = ppg%jxyz(1, j, ia)
              iy = ppg%jxyz(2, j, ia)
              iz = ppg%jxyz(3, j, ia)

              if (ix >= is(1) .and. ix <= ie(1) .and. &
                  iy >= is(2) .and. iy <= ie(2) .and. &
                  iz >= is(3) .and. iz <= ie(3)) then
                lx = mod(ix - iorg(1), dg_frag%lgnum_total(1)) + 1
                ly = mod(iy - iorg(2), dg_frag%lgnum_total(2)) + 1
                lz = mod(iz - iorg(3), dg_frag%lgnum_total(3)) + 1
                if (lx < 1 .or. lx > ndom(1)) cycle
                if (ly < 1 .or. ly > ndom(2)) cycle
                if (lz < 1 .or. lz > ndom(3)) cycle
                x = ppg%rxyz(1, j, ia)
                y = ppg%rxyz(2, j, ia)
                z = ppg%rxyz(3, j, ia)
                if (use_micro_A .and. present(Ac_micro)) then
                  A_local(1:3) = Ac_micro(1:3, ix, iy, iz)
                else
                  A_local(1:3) = Ac_tot(1:3)
                end if
                phase = A_local(1) * x + A_local(2) * y + A_local(3) * z
                overlap_i = overlap_i + dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                            ppg%uV(j, ilma) * exp(-zi * phase) * hvol
              end if
            end do
            uVpsi(io, ilma, ispin) = overlap_i
          end do
        end do
        !$omp end parallel do

        nbf = dg_frag%n_basis(ifrag, ispin)
        iblk = find_matrix_block_runtime(H_nl_block_map, ifrag, ifrag)
        if (iblk <= 0 .or. iblk > size(H_nl_blocks)) then
          write(*,'(a,i0,a,i0,a,i0)') "[FATAL] missing H_nl block: rank=", dg_frag%id, " ifrag=", ifrag, " ispin=", ispin
          stop "missing H_nl block for nonlocal PP"
        end if
        !$omp parallel do collapse(2) private(jo, io, ilma, nlpp_contrib, overlap_i, overlap_j)
        do jo = 1, nbf
          do io = 1, nbf
            nlpp_contrib = (0.0d0, 0.0d0)
            do ilma = 1, ppg%Nlma
              overlap_i = uVpsi(io, ilma, ispin)
              overlap_j = uVpsi(jo, ilma, ispin)
              nlpp_contrib = nlpp_contrib + overlap_i * overlap_j * ppg%rinv_uvu(ilma)
            end do
            H_nl_blocks(iblk)%val(io, jo, ispin) = H_nl_blocks(iblk)%val(io, jo, ispin) + nlpp_contrib
          end do
        end do
        !$omp end parallel do
      end do
    end do

    deallocate(uVpsi)
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
    logical :: reuse_allowed
    logical :: use_micro_A
    logical :: need_dense_cache
    logical :: want_dense_cache
    integer :: i
    logical, parameter :: enable_hmat_nl_progress = .true.

    if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "        hmat-nl trace: stage=entry"
      flush(6)
    end if

    if (ppg%Nlma == 0 .or. .not. allocated(ppg%uV)) then
      if (allocated(dg_frag%H_nl_cache)) deallocate(dg_frag%H_nl_cache)
      if (allocated(dg_frag%H_nl_blocks)) then
        do i = 1, size(dg_frag%H_nl_blocks)
          if (allocated(dg_frag%H_nl_blocks(i)%val)) deallocate(dg_frag%H_nl_blocks(i)%val)
        end do
        deallocate(dg_frag%H_nl_blocks)
      end if
      if (allocated(dg_frag%H_nl_block_map)) deallocate(dg_frag%H_nl_block_map)
      dg_frag%n_H_nl_blocks = 0
      dg_frag%has_nl_cache = .false.
      return
    end if

    use_micro_A = (trim(theory) == 'single_scale_maxwell_tddft' .and. allocated(system%Ac_micro%v))
    want_dense_cache = (.not. allocated(dg_frag%H_mat_blocks)) .or. allocated(dg_frag%H_mat_c)
    if (present(require_dense_cache)) want_dense_cache = require_dense_cache
    need_dense_cache = want_dense_cache

    reuse_allowed = (.not. use_micro_A) .and. &
                    (trim(ae_shape1) == 'impulse' .or. trim(ae_shape1) == 'none') .and. &
                    (trim(ae_shape2) == 'impulse' .or. trim(ae_shape2) == 'none')

    delta_A = maxval(abs(Ac_tot - dg_frag%Ac_nl_cache))
    if (.not. dg_frag%has_nl_cache .or. (.not. reuse_allowed) .or. delta_A > dg_frag%Ac_nl_cache_tol) then
      if (.not. allocated(dg_frag%H_nl_blocks) .or. .not. allocated(dg_frag%H_nl_block_map)) then
        call init_complex_matrix_blocks_runtime(dg_frag, dg_frag%H_nl_blocks, dg_frag%H_nl_block_map)
        dg_frag%n_H_nl_blocks = size(dg_frag%H_nl_blocks)
      end if
      if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        hmat-nl trace: stage=after-init-blocks"
        flush(6)
      end if
      if (use_micro_A) then
        call build_nonlocal_pp_matrix_A_blocks(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
             .true., system%Ac_micro%v, dg_frag%H_nl_blocks, dg_frag%H_nl_block_map)
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
      if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
        write(*,'(1x,a)') "        hmat-nl trace: stage=after-block-reduce"
        flush(6)
      end if
      if (need_dense_cache) then
        if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
          write(*,'(1x,a)') "        hmat-nl trace: stage=before-dense-cache-build"
          flush(6)
        end if
        if (.not. allocated(dg_frag%H_nl_cache)) allocate(dg_frag%H_nl_cache(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
        if (use_micro_A) then
          call build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
               .true., system%Ac_micro%v, dg_frag%H_nl_cache)
        else
          call build_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system%nspin, system%hvol, Ac_tot, &
               .false., H_nl=dg_frag%H_nl_cache)
        end if
        if (enable_hmat_nl_progress .and. dg_frag%id == 0) then
          write(*,'(1x,a)') "        hmat-nl trace: stage=after-dense-cache-build"
          flush(6)
        end if
      else if (allocated(dg_frag%H_nl_cache)) then
        deallocate(dg_frag%H_nl_cache)
      end if
      dg_frag%Ac_nl_cache = Ac_tot
      dg_frag%has_nl_cache = .true.
    end if

  end subroutine ensure_nonlocal_pp_matrix_A

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
      y(n_frag+1:n_tot, :) = y(n_frag+1:n_tot, :) + matmul(conjg(transpose(dg_frag%H_mat_frag_pw(1:n_frag, 1:n_pw, ispin))), x(1:n_frag, :))
    end if

    if (n_pw > 0) then
      if (allocated(dg_frag%H_mat_pw_diag)) then
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
    y(:) = (0.0d0, 0.0d0)

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
    y(:, :) = (0.0d0, 0.0d0)

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
      y(n_frag+1:n_tot, :) = x(n_frag+1:n_tot, :) + matmul(conjg(transpose(dg_frag%S_mat_frag_pw(1:n_frag, 1:n_pw, ispin))), x(1:n_frag, :))
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
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: rhs(:, :)
    complex(8), intent(out) :: sol(:, :)
    logical, intent(in) :: use_prop

    integer :: n_dim, n_rhs, icol, iter, max_iter
    real(8) :: rhs_norm, res_norm
    real(8), parameter :: diag_floor = 1.0d-10, tol_rel = 1.0d-10
    complex(8), allocatable :: r(:,:), z(:,:), p(:,:), ap(:,:)
    real(8), allocatable :: diag(:), tol_abs(:), rho(:), rho_new(:), denom(:)
    complex(8), allocatable :: alpha(:), beta(:)
    logical, allocatable :: active(:)

    n_dim = size(rhs, 1)
    n_rhs = size(rhs, 2)
    sol(:, :) = (0.0d0, 0.0d0)
    if (n_dim <= 0 .or. n_rhs <= 0) return

    allocate(diag(n_dim), r(n_dim, n_rhs), z(n_dim, n_rhs), p(n_dim, n_rhs), ap(n_dim, n_rhs))
    allocate(tol_abs(n_rhs), rho(n_rhs), rho_new(n_rhs), denom(n_rhs))
    allocate(alpha(n_rhs), beta(n_rhs), active(n_rhs))
    call build_overlap_operator_diagonal(dg_frag, ispin, use_prop, diag)
    where (diag < diag_floor) diag = diag_floor

    max_iter = max(20, min(6 * max(1, n_dim), 400))
    sol(:, :) = rhs(:, :)
    call apply_overlap_operator_batch(dg_frag, ispin, sol, ap, use_prop)
    r(:, :) = rhs(:, :) - ap(:, :)

    do icol = 1, n_rhs
      rhs_norm = sqrt(max(0.0d0, real(sum(conjg(rhs(:, icol)) * rhs(:, icol)), kind=8)))
      tol_abs(icol) = max(1.0d-12, tol_rel * max(rhs_norm, 1.0d0))
      res_norm = sqrt(max(0.0d0, real(sum(conjg(r(:, icol)) * r(:, icol)), kind=8)))
      active(icol) = (res_norm > tol_abs(icol))
    end do

    if (.not. any(active)) then
      deallocate(diag, r, z, p, ap, tol_abs, rho, rho_new, denom, alpha, beta, active)
      return
    end if

    do icol = 1, n_rhs
      if (.not. active(icol)) cycle
      z(:, icol) = r(:, icol) / diag(:)
      p(:, icol) = z(:, icol)
      rho(icol) = real(sum(conjg(r(:, icol)) * z(:, icol)), kind=8)
      if (rho(icol) <= 0.0d0) then
        sol(:, icol) = rhs(:, icol)
        active(icol) = .false.
      end if
    end do

    do iter = 1, max_iter
      if (.not. any(active)) exit
      do icol = 1, n_rhs
        if (.not. active(icol)) p(:, icol) = (0.0d0, 0.0d0)
      end do
      call apply_overlap_operator_batch(dg_frag, ispin, p, ap, use_prop)
      do icol = 1, n_rhs
        if (.not. active(icol)) cycle
        denom(icol) = real(sum(conjg(p(:, icol)) * ap(:, icol)), kind=8)
        if (abs(denom(icol)) <= 1.0d-30) then
          active(icol) = .false.
          cycle
        end if
        alpha(icol) = cmplx(rho(icol) / denom(icol), 0.0d0, kind=8)
        sol(:, icol) = sol(:, icol) + alpha(icol) * p(:, icol)
        r(:, icol) = r(:, icol) - alpha(icol) * ap(:, icol)
        res_norm = sqrt(max(0.0d0, real(sum(conjg(r(:, icol)) * r(:, icol)), kind=8)))
        if (res_norm <= tol_abs(icol)) then
          active(icol) = .false.
          cycle
        end if
        z(:, icol) = r(:, icol) / diag(:)
        rho_new(icol) = real(sum(conjg(r(:, icol)) * z(:, icol)), kind=8)
        if (rho(icol) <= 0.0d0) then
          active(icol) = .false.
          cycle
        end if
        beta(icol) = cmplx(rho_new(icol) / rho(icol), 0.0d0, kind=8)
        p(:, icol) = z(:, icol) + beta(icol) * p(:, icol)
        rho(icol) = rho_new(icol)
      end do
    end do

    deallocate(diag, r, z, p, ap, tol_abs, rho, rho_new, denom, alpha, beta, active)
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
    else
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
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin
    real(8), intent(in) :: scale_vec(3)
    complex(8), intent(inout) :: mat(:, :)

    integer :: iblk, idir, ib, jb, row_idx, col_idx, n_frag
    integer :: nrow, ncol, idx_ib, idx_jb, valid_row_count, valid_col_count
    integer :: scratch_size

    mat(:, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%momentum_blocks)) return
    if (.not. allocated(dg_frag%index_basis)) return
    n_frag = dg_frag%n_mat_max
    scratch_size = size(dg_frag%index_basis, 1)
    if (.not. allocated(dg_frag%momentum_dense_row_gid_cache)) then
      allocate(dg_frag%momentum_dense_row_gid_cache(scratch_size))
      allocate(dg_frag%momentum_dense_col_gid_cache(scratch_size))
      allocate(dg_frag%momentum_dense_valid_row_ids(scratch_size))
      allocate(dg_frag%momentum_dense_valid_col_ids(scratch_size))
    else if (size(dg_frag%momentum_dense_row_gid_cache) /= scratch_size) then
      deallocate(dg_frag%momentum_dense_row_gid_cache, dg_frag%momentum_dense_col_gid_cache, &
        dg_frag%momentum_dense_valid_row_ids, dg_frag%momentum_dense_valid_col_ids)
      allocate(dg_frag%momentum_dense_row_gid_cache(scratch_size))
      allocate(dg_frag%momentum_dense_col_gid_cache(scratch_size))
      allocate(dg_frag%momentum_dense_valid_row_ids(scratch_size))
      allocate(dg_frag%momentum_dense_valid_col_ids(scratch_size))
    end if
    do iblk = 1, dg_frag%n_momentum_blocks
      if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
      nrow = dg_frag%n_basis(dg_frag%momentum_blocks(iblk)%ifrag_row, ispin)
      ncol = dg_frag%n_basis(dg_frag%momentum_blocks(iblk)%ifrag_col, ispin)
      if (nrow <= 0 .or. ncol <= 0) cycle
      valid_row_count = 0
      do ib = 1, nrow
        dg_frag%momentum_dense_row_gid_cache(ib) = dg_frag%index_basis(ib, dg_frag%momentum_blocks(iblk)%ifrag_row, ispin)
        if (dg_frag%momentum_dense_row_gid_cache(ib) < 1 .or. dg_frag%momentum_dense_row_gid_cache(ib) > n_frag) cycle
        valid_row_count = valid_row_count + 1
        dg_frag%momentum_dense_valid_row_ids(valid_row_count) = ib
      end do
      valid_col_count = 0
      do jb = 1, ncol
        dg_frag%momentum_dense_col_gid_cache(jb) = dg_frag%index_basis(jb, dg_frag%momentum_blocks(iblk)%ifrag_col, ispin)
        if (dg_frag%momentum_dense_col_gid_cache(jb) < 1 .or. dg_frag%momentum_dense_col_gid_cache(jb) > n_frag) cycle
        valid_col_count = valid_col_count + 1
        dg_frag%momentum_dense_valid_col_ids(valid_col_count) = jb
      end do
      do idir = 1, 3
        if (abs(scale_vec(idir)) <= 0.0d0) cycle
        do idx_jb = 1, valid_col_count
          jb = dg_frag%momentum_dense_valid_col_ids(idx_jb)
          col_idx = dg_frag%momentum_dense_col_gid_cache(jb)
!$omp simd private(ib,row_idx)
          do idx_ib = 1, valid_row_count
            ib = dg_frag%momentum_dense_valid_row_ids(idx_ib)
            row_idx = dg_frag%momentum_dense_row_gid_cache(ib)
            mat(row_idx, col_idx) = mat(row_idx, col_idx) + &
              cmplx(scale_vec(idir) * dg_frag%momentum_blocks(iblk)%val(idir, ib, jb, ispin), 0.0d0, kind=8)
          end do
        end do
      end do
    end do
  end subroutine copy_momentum_blocks_to_complex_dense

  subroutine sync_dense_matrix_to_blocks_runtime(dg_frag, mat, blocks, block_map)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: mat(:, :, :)
    type(matrix_block_info), intent(inout) :: blocks(:)
    integer, intent(in) :: block_map(:, :)

    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

!$omp parallel do collapse(2) private(ifrag_col, ifrag_row, iblk, ispin, ii, jj, ig_i, ig_j, nrow, ncol) schedule(static)
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_matrix_block_runtime(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        blocks(iblk)%val(:, :, :) = 0.0d0
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          do jj = 1, ncol
            ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (ig_j < 1 .or. ig_j > size(mat, 2)) cycle
            do ii = 1, nrow
              ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
              if (ig_i < 1 .or. ig_i > size(mat, 1)) cycle
              blocks(iblk)%val(ii, jj, ispin) = mat(ig_i, ig_j, ispin)
            end do
          end do
        end do
      end do
    end do
!$omp end parallel do
  end subroutine sync_dense_matrix_to_blocks_runtime

  subroutine init_matrix_blocks_runtime(dg_frag, blocks, block_map)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)

    integer :: ifrag_row, ifrag_col, iblk, n_blocks, nrow_max, ncol_max

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)
    call ensure_runtime_neighbor_pair_cache(dg_frag)

    n_blocks = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (is_runtime_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) n_blocks = n_blocks + 1
      end do
    end do
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (.not. is_runtime_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
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

  subroutine sync_blocks_to_dense_matrix_runtime(dg_frag, blocks, block_map, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: block_map(:, :)
    real(8), intent(inout) :: mat(:, :, :)

    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer, allocatable :: row_gid(:), col_gid(:), valid_row_ids(:), valid_col_ids(:)

    mat(:, :, :) = 0.0d0
    allocate(row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1)))
    allocate(valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1)))
    do iblk = 1, size(blocks)
      ifrag_row = blocks(iblk)%ifrag_row
      ifrag_col = blocks(iblk)%ifrag_col
      if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
      if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
      do ispin = 1, dg_frag%nspin
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
        valid_col_count = 0
        do jj = 1, ncol
          col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
          if (col_gid(jj) < 1 .or. col_gid(jj) > size(mat, 2)) cycle
          valid_col_count = valid_col_count + 1
          valid_col_ids(valid_col_count) = jj
        end do
        do idx_jj = 1, valid_col_count
          jj = valid_col_ids(idx_jj)
          ig_j = col_gid(jj)
          do idx_ii = 1, valid_row_count
            ii = valid_row_ids(idx_ii)
            ig_i = row_gid(ii)
            mat(ig_i, ig_j, ispin) = blocks(iblk)%val(ii, jj, ispin)
          end do
        end do
      end do
    end do
    deallocate(row_gid, col_gid, valid_row_ids, valid_col_ids)
  end subroutine sync_blocks_to_dense_matrix_runtime

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
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat

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

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Matrix block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

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
    if (comm_is_root(dg_frag%id)) then
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
    integer :: ifrag_row, ifrag_col, iblk, n_blocks, nrow_max, ncol_max

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)
    call ensure_runtime_neighbor_pair_cache(dg_frag)

    n_blocks = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (.not. is_runtime_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
        n_blocks = n_blocks + 1
      end do
    end do
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (.not. is_runtime_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
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

  subroutine reduce_complex_matrix_blocks_runtime(dg_frag, blocks, label, icomm_reduce)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(complex_matrix_block_info), intent(inout) :: blocks(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat

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

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = real(blocks(iblk)%val(ii, jj, ispin), kind=8)
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
            blocks(iblk)%val(ii, jj, ispin) = cmplx(recv_block(offset_flat), aimag(blocks(iblk)%val(ii, jj, ispin)), kind=8)
            offset_flat = offset_flat + 1
          end do
        end do

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = aimag(blocks(iblk)%val(ii, jj, ispin))
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
            blocks(iblk)%val(ii, jj, ispin) = cmplx(real(blocks(iblk)%val(ii, jj, ispin), kind=8), recv_block(offset_flat), kind=8)
            offset_flat = offset_flat + 1
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (comm_is_root(dg_frag%id)) then
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

    integer :: irow, global_row

    packed(:, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_owner)) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return

!$omp parallel do private(irow, global_row) schedule(static)
    do irow = 1, min(size(row_ids), size(packed, 1))
      global_row = row_ids(irow)
      if (global_row < 1 .or. global_row > size(dg_frag%coef, 1)) cycle
      if (dg_frag%coef_owner(global_row, ispin) /= dg_frag%id) cycle
      packed(irow, 1:size(packed, 2)) = dg_frag%coef(global_row, 1:size(packed, 2), ispin)
    end do
!$omp end parallel do
  end subroutine pack_owned_coef

  subroutine fetch_remote_coef_rows(dg_frag, ispin, row_ids, fetched)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: fetched(:, :)

    integer :: irow, global_row, owner_rank
    complex(8), allocatable :: row_buf(:)

    fetched(:, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_owner)) return
    if (ispin < 1 .or. ispin > dg_frag%nspin) return

    allocate(row_buf(size(fetched, 2)))
    do irow = 1, min(size(row_ids), size(fetched, 1))
      global_row = row_ids(irow)
      row_buf(:) = (0.0d0, 0.0d0)
      if (global_row < 1 .or. global_row > size(dg_frag%coef, 1)) then
        fetched(irow, :) = row_buf(:)
        cycle
      end if
      owner_rank = dg_frag%coef_owner(global_row, ispin)
      if (owner_rank < 0) then
        fetched(irow, :) = row_buf(:)
        cycle
      end if
      if (owner_rank == dg_frag%id) row_buf(:) = dg_frag%coef(global_row, 1:size(fetched, 2), ispin)
      call comm_bcast(row_buf, dg_frag%icomm, owner_rank)
      fetched(irow, :) = row_buf(:)
    end do
    deallocate(row_buf)
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

  subroutine fetch_remote_coef_pw_rows(dg_frag, row_ids, fetched)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: row_ids(:)
    complex(8), intent(out) :: fetched(:, :, :)

    integer :: irow, pw_row, owner_rank
    complex(8), allocatable :: row_buf(:, :)

    fetched(:, :, :) = (0.0d0, 0.0d0)
    if (.not. allocated(dg_frag%coef_pw_owner)) return
    if (.not. allocated(dg_frag%coef_pw)) return

    allocate(row_buf(size(fetched, 2), size(fetched, 3)))
    do irow = 1, min(size(row_ids), size(fetched, 1))
      pw_row = row_ids(irow)
      row_buf(:, :) = (0.0d0, 0.0d0)
      if (pw_row < 1 .or. pw_row > size(dg_frag%coef_pw, 1)) then
        fetched(irow, :, :) = row_buf(:, :)
        cycle
      end if
      owner_rank = dg_frag%coef_pw_owner(pw_row)
      if (owner_rank < 0) then
        fetched(irow, :, :) = row_buf(:, :)
        cycle
      end if
      if (owner_rank == dg_frag%id) row_buf(:, :) = dg_frag%coef_pw(pw_row, 1:size(fetched, 2), 1:size(fetched, 3))
      call comm_bcast(row_buf, dg_frag%icomm, owner_rank)
      fetched(irow, :, :) = row_buf(:, :)
    end do
    deallocate(row_buf)
  end subroutine fetch_remote_coef_pw_rows

  subroutine refresh_pw_coef_cache(dg_frag, nstate_use)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in), optional :: nstate_use

    integer :: i, n_pw, nstate_cache
    integer, allocatable :: pw_row_ids(:)

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

    allocate(pw_row_ids(n_pw))
    do i = 1, n_pw
      pw_row_ids(i) = i
    end do
    dg_frag%coef_pw_full_cache(:, :, :) = (0.0d0, 0.0d0)
    call fetch_remote_coef_pw_rows(dg_frag, pw_row_ids, dg_frag%coef_pw_full_cache)
    dg_frag%coef_pw_full_cache_nstate = nstate_cache
    deallocate(pw_row_ids)
  end subroutine refresh_pw_coef_cache

  subroutine gather_full_coef_view(dg_frag, ispin, n_frag_rows, nstate_use, coef_frag, coef_pw)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ispin
    integer, intent(in) :: n_frag_rows
    integer, intent(in) :: nstate_use
    complex(8), allocatable, intent(out) :: coef_frag(:,:)
    complex(8), allocatable, intent(out) :: coef_pw(:,:)

    integer :: i, n_pw
    integer, allocatable :: frag_row_ids(:), pw_row_ids(:)
    complex(8), allocatable :: coef_pw_all(:,:,:)

    allocate(coef_frag(max(0, n_frag_rows), max(0, nstate_use)))
    coef_frag(:, :) = (0.0d0, 0.0d0)
    if (n_frag_rows > 0 .and. nstate_use > 0) then
      allocate(frag_row_ids(n_frag_rows))
      do i = 1, n_frag_rows
        frag_row_ids(i) = i
      end do
      call fetch_remote_coef_rows(dg_frag, ispin, frag_row_ids, coef_frag)
      deallocate(frag_row_ids)
    end if

    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) n_pw = dg_frag%n_plane_waves
    allocate(coef_pw(max(0, n_pw), max(0, nstate_use)))
    coef_pw(:, :) = (0.0d0, 0.0d0)
    if (n_pw > 0 .and. nstate_use > 0) then
      allocate(pw_row_ids(n_pw))
      do i = 1, n_pw
        pw_row_ids(i) = i
      end do
      allocate(coef_pw_all(n_pw, nstate_use, dg_frag%nspin))
      coef_pw_all(:, :, :) = (0.0d0, 0.0d0)
      call fetch_remote_coef_pw_rows(dg_frag, pw_row_ids, coef_pw_all)
      coef_pw(:, :) = coef_pw_all(:, :, min(max(ispin, 1), dg_frag%nspin))
      deallocate(coef_pw_all)
      deallocate(pw_row_ids)
    end if
  end subroutine gather_full_coef_view

  subroutine zero_nonowned_coefficients(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ispin, i

    if (allocated(dg_frag%coef_owner) .and. allocated(dg_frag%coef)) then
      do ispin = 1, min(size(dg_frag%coef_owner, 2), size(dg_frag%coef, 3))
        do i = 1, min(size(dg_frag%coef_owner, 1), size(dg_frag%coef, 1))
          if (dg_frag%coef_owner(i, ispin) == dg_frag%id) cycle
          dg_frag%coef(i, :, ispin) = (0.0d0, 0.0d0)
        end do
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
        if (abs(xj) == 0.0d0) cycle
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

    integer :: iblk, iblk_idx, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ist, ig_i, ig_j, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))
    complex(8) :: xj

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return

    if (present(block_ids)) then
      if (size(block_ids) <= 0) return
      do iblk_idx = 1, size(block_ids)
        iblk = block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(blocks)) cycle
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
          if (row_gid(ii) < 1 .or. row_gid(ii) > size(y, 1)) cycle
          valid_row_count = valid_row_count + 1
          valid_row_ids(valid_row_count) = ii
        end do
        if (valid_row_count <= 0) cycle
        valid_col_count = 0
        do jj = 1, ncol
          col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
          if (col_gid(jj) < 1 .or. col_gid(jj) > size(x, 1)) cycle
          valid_col_count = valid_col_count + 1
          valid_col_ids(valid_col_count) = jj
        end do
        if (valid_col_count <= 0) cycle
        do idx_jj = 1, valid_col_count
          jj = valid_col_ids(idx_jj)
          ig_j = col_gid(jj)
          do ist = 1, size(x, 2)
            xj = x(ig_j, ist)
            if (abs(xj) == 0.0d0) cycle
!$omp simd private(ii,ig_i)
            do idx_ii = 1, valid_row_count
              ii = valid_row_ids(idx_ii)
              ig_i = row_gid(ii)
              y(ig_i, ist) = y(ig_i, ist) + blocks(iblk)%val(ii, jj, ispin) * xj
            end do
          end do
        end do
      end do
      return
    end if

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
        if (row_gid(ii) < 1 .or. row_gid(ii) > size(y, 1)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ii
      end do
      if (valid_row_count <= 0) cycle
      valid_col_count = 0
      do jj = 1, ncol
        col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (col_gid(jj) < 1 .or. col_gid(jj) > size(x, 1)) cycle
        valid_col_count = valid_col_count + 1
        valid_col_ids(valid_col_count) = jj
      end do
      if (valid_col_count <= 0) cycle
      do idx_jj = 1, valid_col_count
        jj = valid_col_ids(idx_jj)
        ig_j = col_gid(jj)
        do ist = 1, size(x, 2)
          xj = x(ig_j, ist)
          if (abs(xj) == 0.0d0) cycle
!$omp simd private(ii,ig_i)
          do idx_ii = 1, valid_row_count
            ii = valid_row_ids(idx_ii)
            ig_i = row_gid(ii)
            y(ig_i, ist) = y(ig_i, ist) + blocks(iblk)%val(ii, jj, ispin) * xj
          end do
        end do
      end do
    end do
  end subroutine apply_matrix_blocks_batch

  subroutine apply_complex_matrix_blocks_batch(dg_frag, blocks, ispin, x, y, block_ids)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(complex_matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    complex(8), intent(in) :: x(:, :)
    complex(8), intent(inout) :: y(:, :)
    integer, intent(in), optional :: block_ids(:)

    integer :: iblk, iblk_idx, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ist, ig_i, ig_j, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))
    complex(8) :: xj

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (.not. allocated(dg_frag%index_basis)) return

    if (present(block_ids)) then
      if (size(block_ids) <= 0) return
      do iblk_idx = 1, size(block_ids)
        iblk = block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(blocks)) cycle
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
          if (row_gid(ii) < 1 .or. row_gid(ii) > size(y, 1)) cycle
          valid_row_count = valid_row_count + 1
          valid_row_ids(valid_row_count) = ii
        end do
        if (valid_row_count <= 0) cycle
        valid_col_count = 0
        do jj = 1, ncol
          col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
          if (col_gid(jj) < 1 .or. col_gid(jj) > size(x, 1)) cycle
          valid_col_count = valid_col_count + 1
          valid_col_ids(valid_col_count) = jj
        end do
        if (valid_col_count <= 0) cycle
        do idx_jj = 1, valid_col_count
          jj = valid_col_ids(idx_jj)
          ig_j = col_gid(jj)
          do ist = 1, size(x, 2)
            xj = x(ig_j, ist)
            if (abs(xj) == 0.0d0) cycle
!$omp simd private(ii,ig_i)
            do idx_ii = 1, valid_row_count
              ii = valid_row_ids(idx_ii)
              ig_i = row_gid(ii)
              y(ig_i, ist) = y(ig_i, ist) + blocks(iblk)%val(ii, jj, ispin) * xj
            end do
          end do
        end do
      end do
      return
    end if

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
        if (row_gid(ii) < 1 .or. row_gid(ii) > size(y, 1)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ii
      end do
      if (valid_row_count <= 0) cycle
      valid_col_count = 0
      do jj = 1, ncol
        col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (col_gid(jj) < 1 .or. col_gid(jj) > size(x, 1)) cycle
        valid_col_count = valid_col_count + 1
        valid_col_ids(valid_col_count) = jj
      end do
      if (valid_col_count <= 0) cycle
      do idx_jj = 1, valid_col_count
        jj = valid_col_ids(idx_jj)
        ig_j = col_gid(jj)
        do ist = 1, size(x, 2)
          xj = x(ig_j, ist)
          if (abs(xj) == 0.0d0) cycle
!$omp simd private(ii,ig_i)
          do idx_ii = 1, valid_row_count
            ii = valid_row_ids(idx_ii)
            ig_i = row_gid(ii)
            y(ig_i, ist) = y(ig_i, ist) + blocks(iblk)%val(ii, jj, ispin) * xj
          end do
        end do
      end do
    end do
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

  subroutine copy_matrix_blocks_metric_to_complex_dense(dg_frag, blocks, ispin, n_metric, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    integer, intent(in) :: n_metric
    complex(8), intent(inout) :: mat(:, :)

    integer :: iblk, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ig_i, ig_j, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (n_metric <= 0) return
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
        if (row_gid(ii) < 1 .or. row_gid(ii) > n_metric .or. row_gid(ii) > size(mat, 1)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ii
      end do
      if (valid_row_count <= 0) cycle

      valid_col_count = 0
      do jj = 1, ncol
        col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (col_gid(jj) < 1 .or. col_gid(jj) > n_metric .or. col_gid(jj) > size(mat, 2)) cycle
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
  end subroutine copy_matrix_blocks_metric_to_complex_dense

  subroutine copy_matrix_blocks_metric_to_real_dense(dg_frag, blocks, ispin, n_metric, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: ispin
    integer, intent(in) :: n_metric
    real(8), intent(inout) :: mat(:, :)

    integer :: iblk, ifrag_row, ifrag_col
    integer :: nrow, ncol, ii, jj, ig_i, ig_j, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    integer :: valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))

    if (ispin < 1 .or. ispin > dg_frag%nspin) return
    if (n_metric <= 0) return
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
        if (row_gid(ii) < 1 .or. row_gid(ii) > n_metric .or. row_gid(ii) > size(mat, 1)) cycle
        valid_row_count = valid_row_count + 1
        valid_row_ids(valid_row_count) = ii
      end do
      if (valid_row_count <= 0) cycle

      valid_col_count = 0
      do jj = 1, ncol
        col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
        if (col_gid(jj) < 1 .or. col_gid(jj) > n_metric .or. col_gid(jj) > size(mat, 2)) cycle
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
          mat(ig_i, ig_j) = blocks(iblk)%val(ii, jj, ispin)
        end do
      end do
    end do
  end subroutine copy_matrix_blocks_metric_to_real_dense

  subroutine copy_hamiltonian_metric_to_complex_dense(dg_frag, n_metric, mat)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: n_metric
    complex(8), intent(inout) :: mat(:, :, :)

    integer :: ispin, io, jo

    mat(:, :, :) = (0.0d0, 0.0d0)
    if (n_metric <= 0) return

    if (allocated(dg_frag%H_mat_c) .and. allocated(dg_frag%phi_frag_c)) then
      !$omp parallel do collapse(3) private(io,jo,ispin)
      do ispin = 1, dg_frag%nspin
        do jo = 1, n_metric
          do io = 1, n_metric
            mat(io, jo, ispin) = dg_frag%H_mat_c(io, jo, ispin)
          end do
        end do
      end do
      !$omp end parallel do
    else if (allocated(dg_frag%H_mat_blocks)) then
      do ispin = 1, dg_frag%nspin
        call copy_matrix_blocks_metric_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, n_metric, &
             mat(1:n_metric, 1:n_metric, ispin))
      end do
    else if (allocated(dg_frag%H_mat)) then
      !$omp parallel do collapse(3) private(io,jo,ispin)
      do ispin = 1, dg_frag%nspin
        do jo = 1, n_metric
          do io = 1, n_metric
            mat(io, jo, ispin) = cmplx(dg_frag%H_mat(io, jo, ispin), 0.0d0, kind=8)
          end do
        end do
      end do
      !$omp end parallel do
    end if
  end subroutine copy_hamiltonian_metric_to_complex_dense

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
    real(8) :: nabt(4,3)
    integer :: ndom(3)

    nabt = stencil%coef_nab
    ifrag = dg_frag%ifrag_start + i_local - 1
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)

    grad_phi = 0.0d0

    !$omp parallel do collapse(2) private(lx, ly, lz) schedule(static)
    do lz = 1, ndom(3)
      do ly = 1, ndom(2)
        !$omp simd
        do lx = 1, ndom(1)
          grad_phi(lx, ly, lz, 1) = &
              nabt(1,1) * (dg_frag%phi_frag(lx+1, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-1, ly, lz, jo, i_local)) + &
              nabt(2,1) * (dg_frag%phi_frag(lx+2, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-2, ly, lz, jo, i_local)) + &
              nabt(3,1) * (dg_frag%phi_frag(lx+3, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-3, ly, lz, jo, i_local)) + &
              nabt(4,1) * (dg_frag%phi_frag(lx+4, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-4, ly, lz, jo, i_local))

          grad_phi(lx, ly, lz, 2) = &
              nabt(1,2) * (dg_frag%phi_frag(lx, ly+1, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-1, lz, jo, i_local)) + &
              nabt(2,2) * (dg_frag%phi_frag(lx, ly+2, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-2, lz, jo, i_local)) + &
              nabt(3,2) * (dg_frag%phi_frag(lx, ly+3, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-3, lz, jo, i_local)) + &
              nabt(4,2) * (dg_frag%phi_frag(lx, ly+4, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-4, lz, jo, i_local))

          grad_phi(lx, ly, lz, 3) = &
              nabt(1,3) * (dg_frag%phi_frag(lx, ly, lz+1, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-1, jo, i_local)) + &
              nabt(2,3) * (dg_frag%phi_frag(lx, ly, lz+2, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-2, jo, i_local)) + &
              nabt(3,3) * (dg_frag%phi_frag(lx, ly, lz+3, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-3, jo, i_local)) + &
              nabt(4,3) * (dg_frag%phi_frag(lx, ly, lz+4, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-4, jo, i_local))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine apply_gradient_to_basis

  subroutine apply_momentum_blocks(dg_frag, ispin, scale_vec, x, y)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: ispin
    real(8),                intent(in) :: scale_vec(3)
    complex(8),             intent(in) :: x(:,:)
    complex(8),             intent(inout) :: y(:,:)

    integer :: iblk, idir, ib, jb, row_idx, col_idx, nstate, istate
    integer :: active_dir_count, valid_row_count, valid_col_count, idx_dir, idx_ib, idx_jb
    integer :: ifrag_row, ifrag_col, nrow, ncol
    integer :: active_dirs(3), valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1))
    integer :: row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1))
    real(8) :: scale

    if (.not. allocated(dg_frag%momentum_blocks)) return
    nstate = min(size(x, 2), size(y, 2))

!$omp parallel do private(istate,iblk,idir,scale,jb,col_idx,ib,row_idx,active_dir_count,valid_row_count,valid_col_count,idx_dir,idx_ib,idx_jb,ifrag_row,ifrag_col,nrow,ncol,active_dirs,valid_row_ids,valid_col_ids,row_gid,col_gid) schedule(static)
    do istate = 1, nstate
      do iblk = 1, dg_frag%n_momentum_blocks
        if (.not. allocated(dg_frag%momentum_blocks(iblk)%val)) cycle
        active_dir_count = 0
        do idir = 1, 3
          if (abs(scale_vec(idir)) < 1.0d-30) cycle
          active_dir_count = active_dir_count + 1
          active_dirs(active_dir_count) = idir
        end do
        if (active_dir_count <= 0) cycle

        ifrag_row = dg_frag%momentum_blocks(iblk)%ifrag_row
        ifrag_col = dg_frag%momentum_blocks(iblk)%ifrag_col
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

  subroutine calculate_microscopic_current_dg(dg_frag, system, mg, stencil, curr)
    ! NOTE: DG microscopic current here intentionally excludes non-local PP contribution.
    use structures
    use salmon_global, only: nelec
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_vector),         intent(inout) :: curr

    integer :: ifrag, i_local, ispin, io, istate_frag, jstate_frag
    integer :: ig_i, ig_j
    integer :: ix, iy, iz, ixg, iyg, izg, igrid, valid_grid_count
    integer :: nxyz(3), ixyz0(3), ifrag_count, ngrid_max
    integer :: nocc_per_spin
    real(8) :: occ_factor
    complex(8) :: coef_pair
    real(8) :: phi_i
    real(8), allocatable :: grad_phi(:,:,:,:), curr_thread(:,:,:,:)
    real(8), allocatable :: curr_local(:,:,:,:), curr_sum(:,:,:,:)
    real(8), allocatable :: w_local(:,:,:), w_sum(:,:,:)

    curr%v = 0.0d0
    if (.not. allocated(dg_frag%phi_frag)) return

    nocc_per_spin = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
    occ_factor = 2.0d0 / real(system%nspin, 8)

    allocate(curr_local(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(curr_sum(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(w_sum(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    curr_local = 0.0d0
    curr_sum = 0.0d0
    w_local = 0.0d0
    w_sum = 0.0d0

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
          nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
          ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)
          valid_grid_count = 0
          do iz = 1, nxyz(3)
            do iy = 1, nxyz(2)
              do ix = 1, nxyz(1)
                ixg = modulo(ixyz0(1) + ix - 2, mg%num(1)) + 1
                iyg = modulo(ixyz0(2) + iy - 2, mg%num(2)) + 1
                izg = modulo(ixyz0(3) + iz - 2, mg%num(3)) + 1
                if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
                if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
                if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle
                valid_grid_count = valid_grid_count + 1
                dg_frag%current_valid_ix(valid_grid_count, i_local) = ix
                dg_frag%current_valid_iy(valid_grid_count, i_local) = iy
                dg_frag%current_valid_iz(valid_grid_count, i_local) = iz
                dg_frag%current_valid_ixg(valid_grid_count, i_local) = ixg
                dg_frag%current_valid_iyg(valid_grid_count, i_local) = iyg
                dg_frag%current_valid_izg(valid_grid_count, i_local) = izg
              end do
            end do
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
!$omp parallel private(jstate_frag, istate_frag, ig_j, ig_i, io, coef_pair, phi_i, iz, iy, ix, ixg, iyg, izg, grad_phi, curr_thread)
          allocate(curr_thread(3, mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
          curr_thread = 0.0d0
          allocate(grad_phi(1:nxyz(1), 1:nxyz(2), 1:nxyz(3), 3))

!$omp do schedule(static)
          do jstate_frag = 1, dg_frag%n_basis(ifrag, ispin)
            ig_j = dg_frag%index_basis(jstate_frag, ifrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

            call apply_gradient_to_basis(dg_frag, i_local, jstate_frag, mg, stencil, grad_phi)

            do istate_frag = 1, dg_frag%n_basis(ifrag, ispin)
              ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

              coef_pair = (0.0d0, 0.0d0)
              do io = 1, nocc_per_spin
                coef_pair = coef_pair + occ_factor * conjg(dg_frag%coef(ig_i, io, ispin)) * dg_frag%coef(ig_j, io, ispin)
              end do
              if (abs(aimag(coef_pair)) < 1.0d-18) cycle

              do igrid = 1, valid_grid_count
                ix = dg_frag%current_valid_ix(igrid, i_local)
                iy = dg_frag%current_valid_iy(igrid, i_local)
                iz = dg_frag%current_valid_iz(igrid, i_local)
                ixg = dg_frag%current_valid_ixg(igrid, i_local)
                iyg = dg_frag%current_valid_iyg(igrid, i_local)
                izg = dg_frag%current_valid_izg(igrid, i_local)

                phi_i = dg_frag%phi_frag(ix, iy, iz, istate_frag, i_local)
                curr_thread(1, ixg, iyg, izg) = curr_thread(1, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ix, iy, iz, 1)
                curr_thread(2, ixg, iyg, izg) = curr_thread(2, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ix, iy, iz, 2)
                curr_thread(3, ixg, iyg, izg) = curr_thread(3, ixg, iyg, izg) + aimag(coef_pair) * phi_i * grad_phi(ix, iy, iz, 3)
              end do
            end do
          end do
!$omp end do

!$omp critical
          curr_local(:, :, :, :) = curr_local(:, :, :, :) + curr_thread(:, :, :, :)
!$omp end critical

          deallocate(grad_phi, curr_thread)
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

  end subroutine calculate_microscopic_current_dg

  !=======================================================================
  ! Build spatially resolved A·∇ and A^2/2 matrices from Ac_micro(r)
  !=======================================================================
  subroutine build_spatial_A_coupling_matrices(dg_frag, system, mg, stencil, ispin, Ap_mat, A2_mat)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    integer,                intent(in)    :: ispin
    real(8),                intent(out)   :: Ap_mat(dg_frag%n_mat_max, dg_frag%n_mat_max)
    real(8),                intent(out)   :: A2_mat(dg_frag%n_mat_max, dg_frag%n_mat_max)

    integer :: ifrag, i_local, io, jo
    integer :: ig_i, ig_j
    integer :: ix, iy, iz, ixg, iyg, izg
    integer :: nxyz(3), ixyz0(3)
    real(8) :: phi_i, A2val, Ap_int, A2_int
    real(8), allocatable :: grad_phi(:,:,:,:)
    real(8), allocatable :: Ap_spin(:,:,:), A2_spin(:,:,:)
    type(matrix_block_info), allocatable :: Ap_blocks(:), A2_blocks(:)
    integer, allocatable :: A_block_map(:,:)

    Ap_mat = 0.0d0
    A2_mat = 0.0d0
    if (.not. allocated(system%Ac_micro%v)) return

    allocate(Ap_spin(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    allocate(A2_spin(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    Ap_spin(:, :, :) = 0.0d0
    A2_spin(:, :, :) = 0.0d0

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
      ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)

      do jo = 1, dg_frag%n_basis(ifrag, ispin)
        ig_j = dg_frag%index_basis(jo, ifrag, ispin)
        if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

        allocate(grad_phi(1:nxyz(1), 1:nxyz(2), 1:nxyz(3), 3))
        call apply_gradient_to_basis(dg_frag, i_local, jo, mg, stencil, grad_phi)

        !$omp parallel do private(io, ig_i, Ap_int, A2_int, iz, iy, ix, ixg, iyg, izg, phi_i, A2val) schedule(static)
        do io = 1, dg_frag%n_basis(ifrag, ispin)
          ig_i = dg_frag%index_basis(io, ifrag, ispin)
          if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

          Ap_int = 0.0d0
          A2_int = 0.0d0
          do iz = 1, nxyz(3)
            do iy = 1, nxyz(2)
              do ix = 1, nxyz(1)
                ixg = modulo(ixyz0(1) + ix - 2, mg%num(1)) + 1
                iyg = modulo(ixyz0(2) + iy - 2, mg%num(2)) + 1
                izg = modulo(ixyz0(3) + iz - 2, mg%num(3)) + 1
                if (ixg < mg%is(1) .or. ixg > mg%ie(1)) cycle
                if (iyg < mg%is(2) .or. iyg > mg%ie(2)) cycle
                if (izg < mg%is(3) .or. izg > mg%ie(3)) cycle

                phi_i = dg_frag%phi_frag(ix, iy, iz, io, i_local)
                Ap_int = Ap_int + phi_i * ( &
                         system%Ac_micro%v(1, ixg, iyg, izg) * grad_phi(ix, iy, iz, 1) + &
                         system%Ac_micro%v(2, ixg, iyg, izg) * grad_phi(ix, iy, iz, 2) + &
                         system%Ac_micro%v(3, ixg, iyg, izg) * grad_phi(ix, iy, iz, 3) ) * system%hvol

                A2val = system%Ac_micro%v(1, ixg, iyg, izg)**2 + &
                        system%Ac_micro%v(2, ixg, iyg, izg)**2 + &
                        system%Ac_micro%v(3, ixg, iyg, izg)**2
                A2_int = A2_int + 0.5d0 * phi_i * A2val * dg_frag%phi_frag(ix, iy, iz, jo, i_local) * system%hvol
              end do
            end do
          end do

          Ap_spin(ig_i, ig_j, ispin) = Ap_spin(ig_i, ig_j, ispin) + Ap_int
          A2_spin(ig_i, ig_j, ispin) = A2_spin(ig_i, ig_j, ispin) + A2_int
        end do
        !$omp end parallel do

        deallocate(grad_phi)
      end do
    end do

    call init_matrix_blocks_runtime(dg_frag, Ap_blocks, A_block_map)
    call init_matrix_blocks_runtime(dg_frag, A2_blocks, A_block_map)
    call sync_dense_matrix_to_blocks_runtime(dg_frag, Ap_spin, Ap_blocks, A_block_map)
    call sync_dense_matrix_to_blocks_runtime(dg_frag, A2_spin, A2_blocks, A_block_map)
    call reduce_matrix_blocks_runtime(dg_frag, Ap_blocks, "spatial-Ap", dg_frag%icomm)
    call reduce_matrix_blocks_runtime(dg_frag, A2_blocks, "spatial-A2", dg_frag%icomm)
    call sync_blocks_to_dense_matrix_runtime(dg_frag, Ap_blocks, A_block_map, Ap_spin)
    call sync_blocks_to_dense_matrix_runtime(dg_frag, A2_blocks, A_block_map, A2_spin)

    Ap_mat(:, :) = Ap_spin(:, :, ispin)
    A2_mat(:, :) = A2_spin(:, :, ispin)

    if (allocated(Ap_blocks)) then
      do io = 1, size(Ap_blocks)
        if (allocated(Ap_blocks(io)%val)) deallocate(Ap_blocks(io)%val)
      end do
      deallocate(Ap_blocks)
    end if
    if (allocated(A2_blocks)) then
      do io = 1, size(A2_blocks)
        if (allocated(A2_blocks(io)%val)) deallocate(A2_blocks(io)%val)
      end do
      deallocate(A2_blocks)
    end if
    if (allocated(A_block_map)) deallocate(A_block_map)
    deallocate(Ap_spin, A2_spin)

  end subroutine build_spatial_A_coupling_matrices

end module rt_dg_fragment_ops
