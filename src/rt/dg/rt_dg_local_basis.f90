! Fragment-local eigen contraction for DG initial basis.
#include "config.h"
module rt_dg_local_basis
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A, apply_matrix_blocks_batch, &
                                apply_overlap_operator_batch
  use rt_dg_fragment_coefficients, only: rebuild_coef_owner_map
  use rt_dg_fragment_lifecycle, only: reset_basis_dependent_operator_storage
  use rt_dg_initial_state, only: solve_fragment_generalized_eigen, build_fd_occupations
  implicit none
  private
  public :: prepare_fragment_local_eigen_basis
  public :: solve_projected_lcfo_seed_coefficients

contains
  subroutine prepare_fragment_local_eigen_basis(dg_frag, system, mg, ppg, did_contract)
    use structures
    use communication, only: comm_is_root, comm_summation, comm_sync_all
    use salmon_global, only: nelec, temperature, dg_subspace_extra_states
    use phys_constants, only: kB_au
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    logical,                intent(out)   :: did_contract

    integer :: ifrag, i_local, ispin, io, jo, state_col, global_idx, local_idx
    integer :: iblk_h, iblk_s, iblk_nl, nfull, nkeep, nkeep_valid, nstate_old, nstate_new, nstate_seed_old
    integer :: nsolve, n_basis_keep, n_basis_rank, nocc_keep
    integer :: nocc_est, nocc_frag_est, n_basis_min, min_virtual_required, nvirt_est_min
    integer :: ifrag_count, max_keep, max_keep_basis, nstate_cap, nvalid
    integer :: c0, c1, nb, nstate_project, old_local_idx
    integer :: old_n_mat_max, old_coef_rows, new_n_mat_max
    integer :: lb1, ub1, lb2, ub2, lb3, ub3
    logical :: preserve_dense_seed, projection_ready, projection_nan
    integer :: nan_c0, nan_c1
    real(8) :: mu, electron_count, occ_eps, fd_temperature
    real(8) :: Ac_zero(3)
    real(8) :: core_res2_local, core_h2_local, core_rel
    real(8) :: core_sum_local(2), core_sum_global(2)
    real(8), allocatable :: H_local(:,:), S_local(:,:), H_full(:,:), S_full(:,:)
    real(8), allocatable :: HC(:,:), SC(:,:)
    real(8), allocatable :: eig(:), C_keep(:,:), transform(:,:,:)
    real(8), allocatable :: eig_local(:,:,:), eig_global(:,:,:), occ_global(:,:,:)
    real(8), allocatable :: phi_new(:,:,:,:,:), esp_seed(:,:)
    complex(8), allocatable :: sc_batch(:,:)
    complex(8), allocatable :: coef_seed_block_part(:,:,:), coef_seed_block(:,:,:)
    complex(8), allocatable :: coef_old(:,:,:)
    integer, allocatable :: nkeep_local(:,:), nkeep_global(:,:), nkeep_occ_local(:,:), nkeep_occ_global(:,:)
    integer, allocatable :: index_new(:,:,:), index_old(:,:,:)
    integer, allocatable :: n_basis_old(:,:), n_mat_old(:), old_global_to_local(:,:), old_coef_owner(:,:)
    integer, allocatable :: n_basis_new(:,:), n_mat_new(:), index_basis_new(:,:,:)
    integer, allocatable :: new_global_to_local(:,:), new_coef_owner(:,:)

    did_contract = .false.
    preserve_dense_seed = .not. dg_frag%identity_seed_coefficients
    if ((.not. preserve_dense_seed) .and. (.not. dg_frag%defer_fragment_cap_to_local_eigen)) return
    if (dg_frag%fragment_basis_contracted) then
      did_contract = preserve_dense_seed
      return
    end if

    if (dg_frag%nspin /= 1) then
      if (dg_frag%id == 0) then
        if (preserve_dense_seed) then
          write(*,'(1x,a)') "[FATAL] DGDFT/LCFO seed requires fragment-local core-S cleanup, but nspin /= 1 is unsupported."
        else
          write(*,'(1x,a)') "[WARN] fragment-local core-S cleanup currently supports non-SOI nspin=1 only; skipping."
        end if
      end if
      if (preserve_dense_seed) stop "DG-Fragment RT: required core-S cleanup unsupported for nspin"
      return
    end if
    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      if (dg_frag%id == 0) then
        if (preserve_dense_seed) then
          write(*,'(1x,a)') "[FATAL] DGDFT/LCFO seed requires fragment-local core-S cleanup, but mixed PW basis is unsupported."
        else
          write(*,'(1x,a)') "[WARN] fragment-local core-S cleanup is skipped for mixed PW basis in this rollout."
        end if
      end if
      if (preserve_dense_seed) stop "DG-Fragment RT: required core-S cleanup unsupported for mixed PW basis"
      return
    end if
    if (.not. dg_frag%parallel_mode_orbital) then
      if (dg_frag%id == 0) then
        if (preserve_dense_seed) then
          write(*,'(1x,a)') "[FATAL] DGDFT/LCFO seed requires fragment-local core-S cleanup in orbital fragment mode."
        else
          write(*,'(1x,a)') "[WARN] fragment-local core-S cleanup requires orbital fragment parallel mode; skipping."
        end if
      end if
      if (preserve_dense_seed) stop "DG-Fragment RT: required core-S cleanup needs orbital mode"
      return
    end if
    if (.not. allocated(dg_frag%H_mat_core_blocks) .or. .not. allocated(dg_frag%S_mat_blocks)) then
      if (preserve_dense_seed) then
        if (dg_frag%id == 0) write(*,'(1x,a)') "[FATAL] Missing H/S blocks for required DGDFT/LCFO core-S cleanup."
        stop "DG-Fragment RT: missing H/S blocks for core-S cleanup"
      end if
      return
    end if
    if (.not. allocated(dg_frag%H_block_map) .or. .not. allocated(dg_frag%S_block_map)) then
      if (preserve_dense_seed) then
        if (dg_frag%id == 0) write(*,'(1x,a)') "[FATAL] Missing H/S block maps for required DGDFT/LCFO core-S cleanup."
        stop "DG-Fragment RT: missing H/S block maps for core-S cleanup"
      end if
      return
    end if
    if (.not. allocated(dg_frag%phi_frag)) then
      if (preserve_dense_seed) then
        if (dg_frag%id == 0) write(*,'(1x,a)') "[FATAL] Missing real-space fragment basis for required DGDFT/LCFO core-S cleanup."
        stop "DG-Fragment RT: missing phi_frag for core-S cleanup"
      end if
      return
    end if

    nstate_old = dg_frag%nstate_frag
    nstate_seed_old = dg_frag%nstate_tot
    nstate_cap = min(max(1, dg_frag%requested_fragment_basis_cap), nstate_old)
    ifrag_count = max(0, dg_frag%ifrag_end - dg_frag%ifrag_start + 1)
    if (ifrag_count <= 0) then
      if (preserve_dense_seed) then
        if (dg_frag%id == 0) write(*,'(1x,a)') "[FATAL] No local fragment range for required DGDFT/LCFO core-S cleanup."
        stop "DG-Fragment RT: empty local fragment range for core-S cleanup"
      end if
      return
    end if
    max_keep = nstate_old
    occ_eps = 1.0d-12
    fd_temperature = temperature
    if (fd_temperature < 0.0d0) fd_temperature = 300.0d0 * kB_au
    Ac_zero(:) = 0.0d0

    ! The fragment-local core-S cleanup eigenproblem is for the core
    ! Hamiltonian on the buffered fragment basis.  Do not use H_mat_blocks here:
    ! those are the DG propagation blocks after surface-flux terms have been
    ! added.
    call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_zero)

    allocate(transform(nstate_old, max_keep, ifrag_count))
    allocate(nkeep_local(dg_frag%n_frag, dg_frag%nspin), nkeep_global(dg_frag%n_frag, dg_frag%nspin), &
             nkeep_occ_local(dg_frag%n_frag, dg_frag%nspin), nkeep_occ_global(dg_frag%n_frag, dg_frag%nspin))
    allocate(eig_local(max_keep, dg_frag%n_frag, dg_frag%nspin))
    allocate(eig_global(max_keep, dg_frag%n_frag, dg_frag%nspin))
    allocate(occ_global(max_keep, dg_frag%n_frag, dg_frag%nspin))
    transform(:, :, :) = 0.0d0
    nkeep_local(:, :) = 0
    nkeep_occ_local(:, :) = 0
    eig_local(:, :, :) = 0.0d0
    core_res2_local = 0.0d0
    core_h2_local = 0.0d0
    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      ispin = 1
      nfull = min(dg_frag%n_basis(ifrag, ispin), nstate_old)
      if (nfull <= 0) cycle
      ! First diagonalize the full core-cut fragment basis to remove
      ! overcomplete S-null directions.  DGDFT/LCFO seeds keep all S-positive
      ! directions; the legacy identity route may still request a smaller
      ! comparison basis after this cleanup.
      nsolve = nfull
      iblk_h = dg_frag%H_block_map(ifrag, ifrag)
      iblk_s = dg_frag%S_block_map(ifrag, ifrag)
      if (iblk_h <= 0 .or. iblk_s <= 0) then
        write(*,'(1x,a,i0,a,i0)') "[FATAL] missing local H/S block for fragment-local core-S cleanup: rank=", &
          dg_frag%id, " ifrag=", ifrag
        stop "DG-Fragment RT: missing local H/S block"
      end if

      allocate(H_local(nfull, nfull), S_local(nfull, nfull), H_full(nfull, nfull), S_full(nfull, nfull))
      H_local(:, :) = 0.0d0
      S_local(:, :) = 0.0d0
      H_local(1:nfull, 1:nfull) = dg_frag%H_mat_core_blocks(iblk_h)%val(1:nfull, 1:nfull, ispin)
      S_local(1:nfull, 1:nfull) = dg_frag%S_mat_blocks(iblk_s)%val(1:nfull, 1:nfull, ispin)
      if (dg_frag%parallel_mode_orbital) then
        ! Orbital fragment mode keeps the local H/S self block replicated on
        ! every subgroup rank.  Summing here would multiply H and S by the
        ! subgroup size; eigenvalues survive, but eigenvectors become
        ! normalized to isize_frag*S and the contracted basis gets S ~= I/N.
        H_full(:, :) = H_local(:, :)
        S_full(:, :) = S_local(:, :)
      else
        call comm_summation(H_local, H_full, nfull * nfull, dg_frag%icomm_frag)
        call comm_summation(S_local, S_full, nfull * nfull, dg_frag%icomm_frag)
      end if
      if (allocated(dg_frag%H_nl_blocks) .and. allocated(dg_frag%H_nl_block_map)) then
        iblk_nl = dg_frag%H_nl_block_map(ifrag, ifrag)
        if (iblk_nl > 0 .and. iblk_nl <= size(dg_frag%H_nl_blocks)) then
          H_full(1:nfull, 1:nfull) = H_full(1:nfull, 1:nfull) + &
            real(dg_frag%H_nl_blocks(iblk_nl)%val(1:nfull, 1:nfull, ispin), kind=8)
        end if
      end if

      do jo = 1, nfull
        do io = jo + 1, nfull
          H_full(io, jo) = 0.5d0 * (H_full(io, jo) + H_full(jo, io))
          H_full(jo, io) = H_full(io, jo)
          S_full(io, jo) = 0.5d0 * (S_full(io, jo) + S_full(jo, io))
          S_full(jo, io) = S_full(io, jo)
        end do
      end do

      allocate(eig(nsolve), C_keep(nfull, nsolve))
      call solve_fragment_generalized_eigen(H_full, S_full, nfull, nsolve, eig, C_keep, nkeep_valid)
      n_basis_rank = min(nkeep_valid, nsolve)
      if (n_basis_rank <= 0) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] local core-S cleanup removed all S-positive directions: rank=", &
          dg_frag%id, " ifrag=", ifrag, " nfull=", nfull
        stop "DG-Fragment RT: zero local eigen rank"
      end if
      n_basis_keep = min(max_keep, n_basis_rank)
      nocc_keep = n_basis_keep
      if (any(C_keep(1:nfull, 1:nocc_keep) /= C_keep(1:nfull, 1:nocc_keep))) then
        write(*,'(1x,a,i0,a,i0)') "[FATAL] NaN in local eigen transform: rank=", dg_frag%id, &
          " ifrag=", ifrag
        stop "DG-Fragment RT: NaN local eigen transform"
      end if
      allocate(HC(nfull, nocc_keep), SC(nfull, nocc_keep))
      HC(:, :) = matmul(H_full(:, :), C_keep(1:nfull, 1:nocc_keep))
      SC(:, :) = matmul(S_full(:, :), C_keep(1:nfull, 1:nocc_keep))
      core_h2_local = core_h2_local + sum(HC(:, :)**2)
      do jo = 1, nocc_keep
        HC(:, jo) = HC(:, jo) - eig(jo) * SC(:, jo)
      end do
      core_res2_local = core_res2_local + sum(HC(:, :)**2)
      deallocate(HC, SC)
      transform(1:nfull, 1:nocc_keep, i_local) = C_keep(1:nfull, 1:nocc_keep)
      if (dg_frag%is_frag_root) then
        nkeep_local(ifrag, ispin) = n_basis_keep
        nkeep_occ_local(ifrag, ispin) = nocc_keep
        eig_local(1:nocc_keep, ifrag, ispin) = eig(1:nocc_keep)
      end if

      deallocate(eig, C_keep)
      deallocate(H_local, S_local, H_full, S_full)
    end do

    core_sum_local(1) = core_res2_local
    core_sum_local(2) = core_h2_local
    call comm_summation(core_sum_local, core_sum_global, 2, dg_frag%icomm)
    if (comm_is_root(dg_frag%id) .and. core_sum_global(2) > 0.0d0) then
      core_rel = sqrt(core_sum_global(1) / core_sum_global(2))
      write(*,'(1x,a,es12.4,a,es12.4,a,es12.4)') &
        "[DG-CORE-EIG] rel||HC-eSC||/||HC||=", core_rel, &
        " res2=", core_sum_global(1), " Hc2=", core_sum_global(2)
    end if

    call comm_summation(nkeep_local, nkeep_global, dg_frag%n_frag * dg_frag%nspin, dg_frag%icomm)
    call comm_summation(nkeep_occ_local, nkeep_occ_global, dg_frag%n_frag * dg_frag%nspin, dg_frag%icomm)
    call comm_summation(eig_local, eig_global, max_keep * dg_frag%n_frag * dg_frag%nspin, dg_frag%icomm)
    max_keep_basis = max(1, maxval(nkeep_global(1:dg_frag%n_frag, 1:dg_frag%nspin)))
    occ_global(:, :, :) = 0.0d0
    if (preserve_dense_seed) then
      mu = 0.0d0
      electron_count = dble(nelec)
      nocc_est = max(0, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
      nocc_frag_est = (nocc_est + max(1, dg_frag%n_frag) - 1) / max(1, dg_frag%n_frag)
      n_basis_min = minval(nkeep_global(1:dg_frag%n_frag, 1:dg_frag%nspin))
      min_virtual_required = max(8, dg_subspace_extra_states)
      nvirt_est_min = n_basis_min - nocc_frag_est
      if (dg_frag%id == 0) then
        write(*,'(1x,a,4(a,i0))') "[INFO] DG LCFO core-S cleanup capacity:", &
          " nocc=", nocc_est, " nocc_per_frag_est=", nocc_frag_est, &
          " n_basis_min=", n_basis_min, " required_min=", min_virtual_required
      end if
      if (nvirt_est_min < min_virtual_required) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a)') "[FATAL] DGDFT/LCFO core-S cleanup left too few unoccupied fragment-basis states."
          write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] Estimated virtual states per fragment after cleanup: min=", &
            nvirt_est_min, " required_min=", min_virtual_required, " ispin=1"
        end if
        stop "DG-Fragment RT: insufficient unoccupied states after core-S cleanup"
      end if
    else
      call build_fd_occupations(eig_global, nkeep_occ_global, max_keep, dg_frag%n_frag, dg_frag%nspin, &
                                dble(nelec), fd_temperature, occ_global, mu, electron_count)
    end if
    max_keep = max_keep_basis

    if (preserve_dense_seed) then
      nstate_new = nstate_seed_old
    else
      nstate_new = 0
      do ispin = 1, dg_frag%nspin
        do ifrag = 1, dg_frag%n_frag
          do io = 1, nkeep_occ_global(ifrag, ispin)
            if (occ_global(io, ifrag, ispin) > occ_eps) nstate_new = nstate_new + 1
          end do
        end do
      end do
    end if
    if (nstate_new <= 0) then
      write(*,'(1x,a,i0)') "[FATAL] fragment-local core-S cleanup produced no occupied states: rank=", dg_frag%id
      stop "DG-Fragment RT: no occupied states after local core-S cleanup"
    end if
    if (allocated(system%rocc)) then
      if (nstate_new > size(system%rocc, 1)) then
        if (dg_frag%id == 0) then
          if (preserve_dense_seed) then
            write(*,'(1x,a,i0,a,i0)') &
              "[FATAL] RT input nstate is too small for DGDFT/LCFO seed states: nstate_seed=", &
              nstate_new, " system%rocc_dim=", size(system%rocc, 1)
          else
            write(*,'(1x,a,i0,a,i0)') &
              "[FATAL] RT input nstate is too small for FD occupied DG states: nstate_new=", &
              nstate_new, " system%rocc_dim=", size(system%rocc, 1)
          end if
        end if
        stop "DG-Fragment RT: nstate too small for local eigen seed"
      end if
    end if

    lb1 = lbound(dg_frag%phi_frag, 1); ub1 = ubound(dg_frag%phi_frag, 1)
    lb2 = lbound(dg_frag%phi_frag, 2); ub2 = ubound(dg_frag%phi_frag, 2)
    lb3 = lbound(dg_frag%phi_frag, 3); ub3 = ubound(dg_frag%phi_frag, 3)
    allocate(phi_new(lb1:ub1, lb2:ub2, lb3:ub3, max_keep, ifrag_count))
    phi_new(:, :, :, :, :) = 0.0d0
    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      nfull = min(dg_frag%n_basis(ifrag, 1), nstate_old)
      nkeep = min(max_keep, nkeep_global(ifrag, 1))
      do io = 1, nkeep
        do jo = 1, nfull
          if (abs(transform(jo, io, i_local)) <= 0.0d0) cycle
          phi_new(:, :, :, io, i_local) = phi_new(:, :, :, io, i_local) + &
            transform(jo, io, i_local) * dg_frag%phi_frag(:, :, :, jo, i_local)
        end do
      end do
    end do
    if (any(phi_new /= phi_new)) then
      write(*,'(1x,a,i0)') "[FATAL] NaN in contracted DG real-space basis: rank=", dg_frag%id
      stop "DG-Fragment RT: NaN contracted basis"
    end if
    deallocate(dg_frag%phi_frag)
    call move_alloc(phi_new, dg_frag%phi_frag)

    nstate_project = 0
    if (allocated(dg_frag%coef) .and. allocated(dg_frag%S_mat_blocks) .and. allocated(dg_frag%coef_global_to_local)) then
      nstate_project = min(nstate_new, dg_frag%nstate_tot, size(dg_frag%coef, 2))
    end if
    old_n_mat_max = dg_frag%n_mat_max
    old_coef_rows = 0
    if (allocated(dg_frag%coef)) old_coef_rows = size(dg_frag%coef, 1)
    allocate(n_basis_old(size(dg_frag%n_basis, 1), size(dg_frag%n_basis, 2)))
    n_basis_old(:, :) = dg_frag%n_basis(:, :)
    allocate(n_mat_old(size(dg_frag%n_mat)))
    n_mat_old(:) = dg_frag%n_mat(:)
    if (allocated(dg_frag%coef)) call move_alloc(dg_frag%coef, coef_old)
    if (allocated(dg_frag%coef_global_to_local)) then
      allocate(old_global_to_local(size(dg_frag%coef_global_to_local, 1), size(dg_frag%coef_global_to_local, 2)))
      old_global_to_local(:, :) = dg_frag%coef_global_to_local(:, :)
    end if
    if (allocated(dg_frag%coef_owner)) then
      allocate(old_coef_owner(size(dg_frag%coef_owner, 1), size(dg_frag%coef_owner, 2)))
      old_coef_owner(:, :) = dg_frag%coef_owner(:, :)
    end if
    allocate(index_old(size(dg_frag%index_basis, 1), size(dg_frag%index_basis, 2), size(dg_frag%index_basis, 3)))
    index_old(:, :, :) = dg_frag%index_basis(:, :, :)

    allocate(index_new(max_keep, dg_frag%n_frag, dg_frag%nspin))
    index_new(:, :, :) = 0
    do ispin = 1, dg_frag%nspin
      global_idx = 0
      do ifrag = 1, dg_frag%n_frag
        do io = 1, nkeep_global(ifrag, ispin)
          global_idx = global_idx + 1
          index_new(io, ifrag, ispin) = global_idx
        end do
      end do
      dg_frag%n_mat(ispin) = global_idx
    end do
    dg_frag%n_mat_max = max(1, maxval(dg_frag%n_mat(1:dg_frag%nspin)))
    dg_frag%nstate_frag = max_keep
    dg_frag%n_basis(:, :) = nkeep_global(:, :)
    deallocate(dg_frag%index_basis)
    call move_alloc(index_new, dg_frag%index_basis)

    if (allocated(dg_frag%coef_new)) deallocate(dg_frag%coef_new)
    if (allocated(dg_frag%coef_work)) deallocate(dg_frag%coef_work)
    if (preserve_dense_seed) then
      if (nstate_new /= nstate_seed_old) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a,i0,a,i0)') "[FATAL] LCFO seed state count changed during basis cap: before=", &
            nstate_seed_old, " after=", nstate_new
        end if
        stop "DG-Fragment RT: LCFO seed state count changed during basis cap"
      end if
      dg_frag%nstate_tot = nstate_seed_old
    else
      dg_frag%nstate_tot = nstate_new
    end if
    call rebuild_coef_owner_map(dg_frag, "fragment-local-eigen-cap")
    new_n_mat_max = dg_frag%n_mat_max
    allocate(n_basis_new(size(dg_frag%n_basis, 1), size(dg_frag%n_basis, 2)))
    n_basis_new(:, :) = dg_frag%n_basis(:, :)
    allocate(n_mat_new(size(dg_frag%n_mat)))
    n_mat_new(:) = dg_frag%n_mat(:)
    allocate(index_basis_new(size(dg_frag%index_basis, 1), size(dg_frag%index_basis, 2), size(dg_frag%index_basis, 3)))
    index_basis_new(:, :, :) = dg_frag%index_basis(:, :, :)
    if (allocated(dg_frag%coef_global_to_local)) then
      allocate(new_global_to_local(size(dg_frag%coef_global_to_local, 1), size(dg_frag%coef_global_to_local, 2)))
      new_global_to_local(:, :) = dg_frag%coef_global_to_local(:, :)
    end if
    if (allocated(dg_frag%coef_owner)) then
      allocate(new_coef_owner(size(dg_frag%coef_owner, 1), size(dg_frag%coef_owner, 2)))
      new_coef_owner(:, :) = dg_frag%coef_owner(:, :)
    end if
    nvalid = max(1, maxval(dg_frag%local_coef_count(1:dg_frag%nspin)))
    allocate(dg_frag%coef(nvalid, dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%coef(:, :, :) = (0.0d0, 0.0d0)

    if (.not. preserve_dense_seed) then
      if (allocated(system%rocc)) system%rocc(:, :, :) = 0.0d0
      if (allocated(dg_frag%nocc_spin)) dg_frag%nocc_spin(:) = 0
    end if
    if (preserve_dense_seed .and. allocated(dg_frag%esp)) then
      allocate(esp_seed(size(dg_frag%esp, 1), size(dg_frag%esp, 2)))
      esp_seed(:, :) = dg_frag%esp(:, :)
    end if
    if (allocated(dg_frag%esp)) deallocate(dg_frag%esp)
    allocate(dg_frag%esp(dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%esp(:, :) = 0.0d0
    if (allocated(esp_seed)) then
      dg_frag%esp(1:min(size(dg_frag%esp, 1), size(esp_seed, 1)), &
                  1:min(size(dg_frag%esp, 2), size(esp_seed, 2))) = &
        esp_seed(1:min(size(dg_frag%esp, 1), size(esp_seed, 1)), &
                 1:min(size(dg_frag%esp, 2), size(esp_seed, 2)))
      deallocate(esp_seed)
    end if

    ! apply_matrix_blocks_batch_orbital uses coef_owner to fetch old-basis
    ! rows from the rank that owns each sparse coefficient row.
    projection_ready = (nstate_project > 0 .and. allocated(coef_old) .and. &
                        allocated(old_global_to_local) .and. allocated(old_coef_owner) .and. &
                        allocated(new_global_to_local) .and. allocated(new_coef_owner))
    if (preserve_dense_seed .and. .not. projection_ready) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "[FATAL] Cannot project DGDFT/LCFO seed into contracted basis."
        write(*,'(1x,a,2(a,i0),5(a,l1))') "[FATAL] projection prerequisites:", &
          " nstate_project=", nstate_project, " old_coef_rows=", old_coef_rows, &
          " coef_old=", allocated(coef_old), " old_map=", allocated(old_global_to_local), &
          " old_owner=", allocated(old_coef_owner), " new_map=", allocated(new_global_to_local), &
          " new_owner=", allocated(new_coef_owner)
      end if
      stop "DG-Fragment RT: missing LCFO seed projection prerequisites"
    else if (.not. projection_ready) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "[INFO] DG identity/local seed projection skipped; rebuilding local seed."
      end if
      nstate_project = 0
    end if
    if (projection_ready) then
      ! S_mat_blocks still represent the pre-contraction basis.  Use the old
      ! row map while applying those old blocks, but write the projected result
      ! with the saved new map arrays.  Restore dg_frag's new map once at the end
      ! to avoid block-by-block allocatable component copies.
      projection_nan = .false.
      nan_c0 = 0
      nan_c1 = 0
      dg_frag%n_mat_max = old_n_mat_max
      dg_frag%n_basis = n_basis_old
      dg_frag%n_mat = n_mat_old
      dg_frag%index_basis = index_old
      dg_frag%coef_global_to_local = old_global_to_local
      dg_frag%coef_owner = old_coef_owner
      spin_projection_loop: do ispin = 1, dg_frag%nspin
        state_projection_loop: do c0 = 1, nstate_project, 64
          c1 = min(nstate_project, c0 + 63)
          nb = c1 - c0 + 1
          allocate(sc_batch(old_coef_rows, nb))
          allocate(coef_seed_block_part(max_keep, nb, ifrag_count))
          allocate(coef_seed_block(max_keep, nb, ifrag_count))
          sc_batch(:, :) = (0.0d0, 0.0d0)
          coef_seed_block_part(:, :, :) = (0.0d0, 0.0d0)
          coef_seed_block(:, :, :) = (0.0d0, 0.0d0)
          call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, coef_old(:, c0:c1, ispin), sc_batch)
          do i_local = 1, ifrag_count
            ifrag = dg_frag%ifrag_start + i_local - 1
            nfull = min(n_basis_old(ifrag, ispin), nstate_old)
            nkeep = min(max_keep, nkeep_global(ifrag, ispin))
            do io = 1, nkeep
              do jo = 1, nfull
                if (abs(transform(jo, io, i_local)) <= 0.0d0) cycle
                global_idx = index_old(jo, ifrag, ispin)
                old_local_idx = 0
                if (global_idx >= 1 .and. global_idx <= old_n_mat_max) &
                  old_local_idx = old_global_to_local(global_idx, ispin)
                if (old_local_idx > 0 .and. old_local_idx <= size(sc_batch, 1)) then
                  coef_seed_block_part(io, 1:nb, i_local) = &
                    coef_seed_block_part(io, 1:nb, i_local) + &
                    transform(jo, io, i_local) * sc_batch(old_local_idx, 1:nb)
                end if
              end do
            end do
          end do
          call comm_summation(coef_seed_block_part, coef_seed_block, size(coef_seed_block_part), dg_frag%icomm_frag)
          if (any(real(coef_seed_block) /= real(coef_seed_block)) .or. &
              any(aimag(coef_seed_block) /= aimag(coef_seed_block))) then
            projection_nan = .true.
            nan_c0 = c0
            nan_c1 = c1
            deallocate(sc_batch, coef_seed_block_part, coef_seed_block)
            exit state_projection_loop
          end if
          do i_local = 1, ifrag_count
            ifrag = dg_frag%ifrag_start + i_local - 1
            nkeep = min(max_keep, nkeep_global(ifrag, ispin))
            do io = 1, nkeep
              global_idx = index_basis_new(io, ifrag, ispin)
              local_idx = 0
              if (global_idx >= 1 .and. global_idx <= new_n_mat_max) &
                local_idx = new_global_to_local(global_idx, ispin)
              if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
                dg_frag%coef(local_idx, c0:c1, ispin) = coef_seed_block(io, 1:nb, i_local)
              end if
            end do
          end do
          deallocate(sc_batch, coef_seed_block_part, coef_seed_block)
        end do state_projection_loop
        if (projection_nan) exit spin_projection_loop
      end do spin_projection_loop
      dg_frag%n_mat_max = new_n_mat_max
      dg_frag%n_basis = n_basis_new
      dg_frag%n_mat = n_mat_new
      dg_frag%index_basis = index_basis_new
      dg_frag%coef_global_to_local = new_global_to_local
      dg_frag%coef_owner = new_coef_owner
      if (projection_nan) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] NaN in LCFO coefficient projection block: rank=", &
          dg_frag%id, " c0=", nan_c0, " c1=", nan_c1
        stop "DG-Fragment RT: NaN projected LCFO coefficient block"
      end if
      if (dg_frag%id == 0) then
        if (dg_frag%identity_seed_coefficients) then
          write(*,'(1x,a,i0)') "[INFO] DG projected identity/local seed RHS into contracted basis: nstate=", nstate_project
        else
          write(*,'(1x,a,i0)') "[INFO] DG projected LCFO coefficient RHS into contracted basis: nstate=", nstate_project
        end if
      end if
    end if
    if (allocated(coef_old)) deallocate(coef_old)
    deallocate(n_basis_old, n_mat_old)
    if (allocated(old_global_to_local)) deallocate(old_global_to_local)
    if (allocated(old_coef_owner)) deallocate(old_coef_owner)
    deallocate(n_basis_new, n_mat_new, index_basis_new)
    if (allocated(new_global_to_local)) deallocate(new_global_to_local)
    if (allocated(new_coef_owner)) deallocate(new_coef_owner)
    deallocate(index_old)

    if (.not. preserve_dense_seed) then
      state_col = 0
      do ispin = 1, dg_frag%nspin
        do ifrag = 1, dg_frag%n_frag
          do io = 1, nkeep_occ_global(ifrag, ispin)
            if (occ_global(io, ifrag, ispin) <= occ_eps) cycle
            state_col = state_col + 1
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            local_idx = 0
            if (global_idx >= 1 .and. global_idx <= dg_frag%n_mat_max) local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
            if (state_col > nstate_project .and. local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
              dg_frag%coef(local_idx, state_col, ispin) = (1.0d0, 0.0d0)
            end if
            if (allocated(system%rocc)) system%rocc(state_col, 1, ispin) = occ_global(io, ifrag, ispin)
            dg_frag%esp(state_col, ispin) = eig_global(io, ifrag, ispin)
            dg_frag%nocc_spin(ispin) = max(dg_frag%nocc_spin(ispin), state_col)
          end do
        end do
      end do
    end if
    allocate(dg_frag%coef_work(size(dg_frag%coef, 1), size(dg_frag%coef, 2), size(dg_frag%coef, 3)))
    dg_frag%coef_work(:, :, :) = dg_frag%coef(:, :, :)
    if (dg_frag%yn_adaptive_basis) then
      allocate(dg_frag%coef_new(size(dg_frag%coef, 1), size(dg_frag%coef, 2), size(dg_frag%coef, 3)))
      dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    end if

    call reset_basis_dependent_operator_storage(dg_frag)
    if (allocated(dg_frag%H_mat_old)) then
      deallocate(dg_frag%H_mat_old)
      allocate(dg_frag%H_mat_old(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin))
      dg_frag%H_mat_old(:, :, :) = (0.0d0, 0.0d0)
    end if
    if (allocated(dg_frag%eigenvalues)) then
      deallocate(dg_frag%eigenvalues)
      allocate(dg_frag%eigenvalues(dg_frag%nstate_frag, dg_frag%nspin))
      dg_frag%eigenvalues(:, :) = 0.0d0
    end if
    if (allocated(dg_frag%basis_overlap)) then
      deallocate(dg_frag%basis_overlap)
      allocate(dg_frag%basis_overlap(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin))
      dg_frag%basis_overlap(:, :, :) = 0.0d0
    end if

    dg_frag%fragment_basis_contracted = .true.
    dg_frag%defer_fragment_cap_to_local_eigen = .false.
    did_contract = .true.
    if (dg_frag%id == 0) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "[INFO] DG local core-S eigen basis contracted: final_keep=", max_keep, &
        " n_mat=", dg_frag%n_mat_max, " nstate_tot=", dg_frag%nstate_tot
      if (dg_frag%requested_fragment_basis_cap > 0 .and. .not. preserve_dense_seed) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[INFO] DG legacy identity core-S basis target: legacy_requested=", &
          nstate_cap, " current_keep=", max_keep, " file_basis=", nstate_old
      else if (preserve_dense_seed) then
        write(*,'(1x,a)') "[INFO] DG local core-S cleanup kept all S-positive directions for LCFO seed."
      end if
      if (preserve_dense_seed) then
        write(*,'(1x,a,1pe13.5)') &
          "[INFO] DG local eigen spectrum diagnostic only; LCFO occupations unchanged: Ne=", electron_count
      else
        write(*,'(1x,a,1pe13.5,a,1pe13.5,a,1pe13.5)') "[INFO] DG local eigen FD: mu=", mu, &
          " Ne=", electron_count, " kT=", fd_temperature
      end if
    end if

    call comm_sync_all(dg_frag%icomm)
    deallocate(transform, nkeep_local, nkeep_global, nkeep_occ_local, nkeep_occ_global, eig_local, eig_global, occ_global)
  end subroutine prepare_fragment_local_eigen_basis

  subroutine solve_projected_lcfo_seed_coefficients(dg_frag, did_solve)
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    logical, intent(out) :: did_solve

    integer, parameter :: seed_solve_block_size = 64
    integer :: ispin, c0, c1, nb, nrow
    complex(8), allocatable :: rhs(:, :), sol(:, :), ssol(:, :)
    real(8) :: res2_local, rhs2_local, res2_global, rhs2_global, rel
    real(8), parameter :: seed_solve_rel_tol = 1.0d-6

    did_solve = .false.
    if (dg_frag%identity_seed_coefficients) return
    if (.not. dg_frag%fragment_basis_contracted) return
    if (.not. allocated(dg_frag%coef)) return

    nrow = size(dg_frag%coef, 1)
    if (nrow <= 0 .or. dg_frag%nstate_tot <= 0) return
    if (.not. allocated(dg_frag%S_mat_blocks) .and. .not. allocated(dg_frag%S_mat_blocks_c)) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') '[FATAL] Cannot solve projected LCFO seed coefficients without rebuilt DG overlap blocks.'
        flush(6)
      end if
      stop 'DG-Fragment RT: missing overlap blocks for projected LCFO seed solve'
    end if

    res2_local = 0.0d0
    rhs2_local = 0.0d0
    do ispin = 1, dg_frag%nspin
      do c0 = 1, dg_frag%nstate_tot, seed_solve_block_size
        c1 = min(dg_frag%nstate_tot, c0 + seed_solve_block_size - 1)
        nb = c1 - c0 + 1
        allocate(rhs(nrow, nb), sol(nrow, nb), ssol(nrow, nb))
        rhs(:, :) = dg_frag%coef(:, c0:c1, ispin)
        sol(:, :) = (0.0d0, 0.0d0)
        call solve_projected_seed_overlap_full(dg_frag, ispin, rhs, sol)
        ssol(:, :) = (0.0d0, 0.0d0)
        call apply_overlap_operator_batch(dg_frag, ispin, sol, ssol, .false.)
        res2_local = res2_local + sum(abs(ssol(:, :) - rhs(:, :))**2)
        rhs2_local = rhs2_local + sum(abs(rhs(:, :))**2)
        dg_frag%coef(:, c0:c1, ispin) = sol(:, :)
        deallocate(rhs, sol, ssol)
      end do
    end do
    if (allocated(dg_frag%coef_work)) then
      if (size(dg_frag%coef_work, 1) /= size(dg_frag%coef, 1) .or. &
          size(dg_frag%coef_work, 2) /= size(dg_frag%coef, 2) .or. &
          size(dg_frag%coef_work, 3) /= size(dg_frag%coef, 3)) then
        deallocate(dg_frag%coef_work)
      end if
    end if
    if (.not. allocated(dg_frag%coef_work)) then
      allocate(dg_frag%coef_work(size(dg_frag%coef, 1), size(dg_frag%coef, 2), size(dg_frag%coef, 3)))
    end if
    dg_frag%coef_work(:, :, :) = dg_frag%coef(:, :, :)
    if (allocated(dg_frag%coef_new)) then
      if (size(dg_frag%coef_new, 1) /= size(dg_frag%coef, 1) .or. &
          size(dg_frag%coef_new, 2) /= size(dg_frag%coef, 2) .or. &
          size(dg_frag%coef_new, 3) /= size(dg_frag%coef, 3)) then
        deallocate(dg_frag%coef_new)
      end if
    end if
    if (dg_frag%yn_adaptive_basis) then
      if (.not. allocated(dg_frag%coef_new)) then
        allocate(dg_frag%coef_new(size(dg_frag%coef, 1), size(dg_frag%coef, 2), size(dg_frag%coef, 3)))
      end if
      dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    else if (allocated(dg_frag%coef_new)) then
      deallocate(dg_frag%coef_new)
    end if

    call comm_summation(res2_local, res2_global, dg_frag%icomm)
    call comm_summation(rhs2_local, rhs2_global, dg_frag%icomm)
    if (comm_is_root(dg_frag%id)) then
      rel = 0.0d0
      if (rhs2_global > 0.0d0) rel = sqrt(max(0.0d0, res2_global) / rhs2_global)
      write(*,'(1x,a,1pe13.5)') '[DG-DCDFT-SEED] solved projected LCFO coefficients with S_new residual=', rel
      flush(6)
    end if
    rel = 0.0d0
    if (rhs2_global > 0.0d0) rel = sqrt(max(0.0d0, res2_global) / rhs2_global)
    if (rel /= rel .or. rel > seed_solve_rel_tol) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,2(a,1pe13.5))') '[FATAL] projected LCFO coefficient solve did not converge:', &
          ' residual=', rel, ' tol=', seed_solve_rel_tol
        flush(6)
      end if
      stop 'DG-Fragment RT: projected LCFO seed solve residual too large'
    end if
    did_solve = .true.
  end subroutine solve_projected_lcfo_seed_coefficients

  subroutine solve_projected_seed_overlap_full(dg_frag, ispin, rhs, sol)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin
    complex(8), intent(in) :: rhs(:, :)
    complex(8), intent(out) :: sol(:, :)

    integer :: n_dim, n_rhs, icol, iter, max_iter
    real(8), parameter :: tol_rel = 1.0d-10
    complex(8), allocatable :: r(:, :), p(:, :), ap(:, :)
    real(8), allocatable :: rhs_norm_l(:), rhs_norm_g(:)
    real(8), allocatable :: rho_l(:), rho_g(:), rho_new_l(:), rho_new_g(:), denom_l(:), denom_g(:)
    real(8), allocatable :: tol_abs(:)
    logical, allocatable :: active(:)
    real(8) :: alpha, beta

    n_dim = size(rhs, 1)
    n_rhs = size(rhs, 2)
    sol(:, :) = (0.0d0, 0.0d0)
    if (n_dim <= 0 .or. n_rhs <= 0) return

    allocate(r(n_dim, n_rhs), p(n_dim, n_rhs), ap(n_dim, n_rhs))
    allocate(rhs_norm_l(n_rhs), rhs_norm_g(n_rhs), rho_l(n_rhs), rho_g(n_rhs))
    allocate(rho_new_l(n_rhs), rho_new_g(n_rhs), denom_l(n_rhs), denom_g(n_rhs), tol_abs(n_rhs), active(n_rhs))

    sol(:, :) = rhs(:, :)
    ap(:, :) = (0.0d0, 0.0d0)
    call apply_overlap_operator_batch(dg_frag, ispin, sol, ap, .false.)
    r(:, :) = rhs(:, :) - ap(:, :)
    p(:, :) = r(:, :)

    rhs_norm_l(:) = 0.0d0
    rho_l(:) = 0.0d0
    do icol = 1, n_rhs
      rhs_norm_l(icol) = real(sum(conjg(rhs(:, icol)) * rhs(:, icol)), kind=8)
      rho_l(icol) = real(sum(conjg(r(:, icol)) * r(:, icol)), kind=8)
    end do
    call comm_summation(rhs_norm_l, rhs_norm_g, n_rhs, dg_frag%icomm)
    call comm_summation(rho_l, rho_g, n_rhs, dg_frag%icomm)
    do icol = 1, n_rhs
      tol_abs(icol) = max(1.0d-20, (tol_rel**2) * max(rhs_norm_g(icol), 1.0d0))
      active(icol) = (rho_g(icol) > tol_abs(icol))
      if (.not. active(icol)) p(:, icol) = (0.0d0, 0.0d0)
    end do

    max_iter = max(20, min(6 * max(1, dg_frag%n_mat_max), 600))
    do iter = 1, max_iter
      if (.not. any(active)) exit
      ap(:, :) = (0.0d0, 0.0d0)
      call apply_overlap_operator_batch(dg_frag, ispin, p, ap, .false.)
      denom_l(:) = 0.0d0
      do icol = 1, n_rhs
        if (active(icol)) denom_l(icol) = real(sum(conjg(p(:, icol)) * ap(:, icol)), kind=8)
      end do
      call comm_summation(denom_l, denom_g, n_rhs, dg_frag%icomm)

      rho_new_l(:) = 0.0d0
      do icol = 1, n_rhs
        if (.not. active(icol)) cycle
        if (abs(denom_g(icol)) <= 1.0d-30) then
          active(icol) = .false.
          p(:, icol) = (0.0d0, 0.0d0)
          cycle
        end if
        alpha = rho_g(icol) / denom_g(icol)
        sol(:, icol) = sol(:, icol) + alpha * p(:, icol)
        r(:, icol) = r(:, icol) - alpha * ap(:, icol)
        rho_new_l(icol) = real(sum(conjg(r(:, icol)) * r(:, icol)), kind=8)
      end do
      call comm_summation(rho_new_l, rho_new_g, n_rhs, dg_frag%icomm)
      do icol = 1, n_rhs
        if (.not. active(icol)) cycle
        if (rho_new_g(icol) <= tol_abs(icol)) then
          active(icol) = .false.
          p(:, icol) = (0.0d0, 0.0d0)
          cycle
        end if
        if (rho_g(icol) <= 0.0d0) then
          active(icol) = .false.
          p(:, icol) = (0.0d0, 0.0d0)
          cycle
        end if
        beta = rho_new_g(icol) / rho_g(icol)
        p(:, icol) = r(:, icol) + beta * p(:, icol)
        rho_g(icol) = rho_new_g(icol)
      end do
    end do

    deallocate(r, p, ap, rhs_norm_l, rhs_norm_g, rho_l, rho_g, rho_new_l, rho_new_g)
    deallocate(denom_l, denom_g, tol_abs, active)
  end subroutine solve_projected_seed_overlap_full

end module rt_dg_local_basis
