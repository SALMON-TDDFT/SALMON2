! Fragment-local eigen contraction for DG initial basis.
#include "config.h"
module rt_dg_local_basis
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A, apply_matrix_blocks_batch
  use rt_dg_fragment_coefficients, only: rebuild_coef_owner_map
  use rt_dg_fragment_lifecycle, only: reset_basis_dependent_operator_storage
  use rt_dg_initial_state, only: solve_fragment_generalized_eigen, build_fd_occupations
  implicit none
  private
  public :: prepare_fragment_local_eigen_basis

contains
  subroutine prepare_fragment_local_eigen_basis(dg_frag, system, mg, ppg, did_contract)
    use structures
    use communication, only: comm_is_root, comm_summation, comm_sync_all, comm_get_min, comm_get_max
    use eigen_subdiag_sub, only: eigen_dsyev
    use salmon_global, only: nelec, temperature
    use phys_constants, only: kB_au
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_pp_grid),        intent(in)    :: ppg
    logical,                intent(out)   :: did_contract

    integer :: ifrag, i_local, ispin, io, jo, state_col, global_idx, local_idx
    integer :: iblk_h, iblk_s, iblk_nl, nfull, nkeep, nkeep_valid, nstate_old, nstate_new
    integer :: nsolve, n_basis_keep, nocc_keep
    integer :: ifrag_count, max_keep, nstate_cap, nvalid
    integer :: relax_cap, env_len, env_status, read_status
    integer :: c0, c1, nb, nstate_project, old_local_idx
    integer :: lb1, ub1, lb2, ub2, lb3, ub3
    real(8) :: mu, electron_count, occ_eps, fd_temperature
    real(8) :: Ac_zero(3)
    real(8) :: core_res2_local, core_h2_local, core_rel
    real(8) :: core_sum_local(2), core_sum_global(2)
    real(8), allocatable :: H_local(:,:), S_local(:,:), H_full(:,:), S_full(:,:)
    real(8), allocatable :: HC(:,:), SC(:,:)
    real(8), allocatable :: eig(:), C_keep(:,:), transform(:,:,:)
    real(8), allocatable :: eig_local(:,:,:), eig_global(:,:,:), occ_global(:,:,:)
    real(8), allocatable :: phi_new(:,:,:,:,:)
    complex(8), allocatable :: sc_batch(:,:), coef_seed_part(:,:,:,:), coef_seed(:,:,:,:)
    integer, allocatable :: nkeep_local(:,:), nkeep_global(:,:), nkeep_occ_local(:,:), nkeep_occ_global(:,:)
    integer, allocatable :: index_new(:,:,:)
    character(len=64) :: env_value

    did_contract = .false.
    if (.not. dg_frag%defer_fragment_cap_to_local_eigen) return
    if (dg_frag%fragment_basis_contracted) return

    if (dg_frag%nspin /= 1) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "[WARN] fragment-local eigen cap currently supports non-SOI nspin=1 only; skipping."
      end if
      return
    end if
    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "[WARN] fragment-local eigen cap is skipped for mixed PW basis in this rollout."
      end if
      return
    end if
    if (.not. dg_frag%parallel_mode_orbital) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "[WARN] fragment-local eigen cap requires orbital fragment parallel mode; skipping."
      end if
      return
    end if
    if (.not. allocated(dg_frag%H_mat_core_blocks) .or. .not. allocated(dg_frag%S_mat_blocks)) return
    if (.not. allocated(dg_frag%H_block_map) .or. .not. allocated(dg_frag%S_block_map)) return
    if (.not. allocated(dg_frag%phi_frag)) return

    nstate_old = dg_frag%nstate_frag
    nstate_cap = min(max(1, dg_frag%requested_fragment_basis_cap), nstate_old)
    ifrag_count = max(0, dg_frag%ifrag_end - dg_frag%ifrag_start + 1)
    if (ifrag_count <= 0) return
    relax_cap = nstate_old
    env_value = ''
    call get_environment_variable('SALMON_DG_FLUX_RELAX_CAP', env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_value(1:env_len), *, iostat=read_status) relax_cap
      if (read_status /= 0) relax_cap = nstate_old
    end if
    max_keep = min(max(relax_cap, nstate_cap), nstate_old)
    occ_eps = 1.0d-12
    fd_temperature = temperature
    if (fd_temperature < 0.0d0) fd_temperature = 300.0d0 * kB_au
    Ac_zero(:) = 0.0d0

    ! The fragment-local cap eigenproblem is for the core Hamiltonian on the
    ! buffered fragment basis.  Do not use H_mat_blocks here: those are the DG
    ! propagation blocks after surface-flux terms have been added.
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
      n_basis_keep = min(max_keep, nfull)
      nsolve = min(nstate_cap, n_basis_keep)
      iblk_h = dg_frag%H_block_map(ifrag, ifrag)
      iblk_s = dg_frag%S_block_map(ifrag, ifrag)
      if (iblk_h <= 0 .or. iblk_s <= 0) then
        write(*,'(1x,a,i0,a,i0)') "[FATAL] missing local H/S block for fragment-local eigen cap: rank=", &
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

      allocate(eig(nsolve), C_keep(nfull, nsolve))
      call solve_fragment_generalized_eigen(H_full, S_full, nfull, nsolve, eig, C_keep, nkeep_valid)
      nocc_keep = min(nkeep_valid, nsolve)
      if (nocc_keep <= 0) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] local eigen cap removed all S-positive directions: rank=", &
          dg_frag%id, " ifrag=", ifrag, " nfull=", nfull
        stop "DG-Fragment RT: zero local eigen rank"
      end if
      if (nocc_keep < nsolve) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] local eigen cap could not keep occupied basis: rank=", &
          dg_frag%id, " ifrag=", ifrag, " kept=", nocc_keep, " requested=", nsolve
        stop "DG-Fragment RT: insufficient local eigen occupied basis"
      end if
      if (any(C_keep(1:nfull, 1:nocc_keep) /= C_keep(1:nfull, 1:nocc_keep))) then
        write(*,'(1x,a,i0,a,i0)') "[FATAL] NaN in local eigen transform: rank=", dg_frag%id, &
          " ifrag=", ifrag
        stop "DG-Fragment RT: NaN local eigen transform"
      end if
      do jo = 1, nfull
        do io = jo + 1, nfull
          H_full(io, jo) = 0.5d0 * (H_full(io, jo) + H_full(jo, io))
          H_full(jo, io) = H_full(io, jo)
          S_full(io, jo) = 0.5d0 * (S_full(io, jo) + S_full(jo, io))
          S_full(jo, io) = S_full(io, jo)
        end do
      end do
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
      do io = nocc_keep + 1, n_basis_keep
        transform(io, io, i_local) = 1.0d0
      end do
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
    occ_global(:, :, :) = 0.0d0
    call build_fd_occupations(eig_global, nkeep_occ_global, max_keep, dg_frag%n_frag, dg_frag%nspin, &
                              dble(nelec), fd_temperature, occ_global, mu, electron_count)

    nstate_new = 0
    do ispin = 1, dg_frag%nspin
      do ifrag = 1, dg_frag%n_frag
        do io = 1, nkeep_occ_global(ifrag, ispin)
          if (occ_global(io, ifrag, ispin) > occ_eps) nstate_new = nstate_new + 1
        end do
      end do
    end do
    if (nstate_new <= 0) then
      write(*,'(1x,a,i0)') "[FATAL] fragment-local eigen cap produced no occupied states: rank=", dg_frag%id
      stop "DG-Fragment RT: no occupied states after local eigen cap"
    end if
    if (allocated(system%rocc)) then
      if (nstate_new > size(system%rocc, 1)) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a,i0,a,i0)') "[FATAL] RT input nstate is too small for FD occupied DG states: nstate_new=", &
            nstate_new, " system%rocc_dim=", size(system%rocc, 1)
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
    if (nstate_project > 0) then
      allocate(coef_seed_part(max_keep, nstate_new, dg_frag%nspin, ifrag_count))
      allocate(coef_seed(max_keep, nstate_new, dg_frag%nspin, ifrag_count))
      coef_seed_part(:, :, :, :) = (0.0d0, 0.0d0)
      coef_seed(:, :, :, :) = (0.0d0, 0.0d0)
      do ispin = 1, dg_frag%nspin
        do c0 = 1, nstate_project, 64
          c1 = min(nstate_project, c0 + 63)
          nb = c1 - c0 + 1
          allocate(sc_batch(size(dg_frag%coef, 1), nb))
          sc_batch(:, :) = (0.0d0, 0.0d0)
          call apply_matrix_blocks_batch(dg_frag, dg_frag%S_mat_blocks, ispin, &
                                         dg_frag%coef(:, c0:c1, ispin), sc_batch)
          do i_local = 1, ifrag_count
            ifrag = dg_frag%ifrag_start + i_local - 1
            nfull = min(dg_frag%n_basis(ifrag, ispin), nstate_old)
            nkeep = min(max_keep, nkeep_global(ifrag, ispin))
            do io = 1, nkeep
              do jo = 1, nfull
                if (abs(transform(jo, io, i_local)) <= 0.0d0) cycle
                global_idx = dg_frag%index_basis(jo, ifrag, ispin)
                old_local_idx = 0
                if (global_idx >= 1 .and. global_idx <= dg_frag%n_mat_max) &
                  old_local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
                if (old_local_idx > 0 .and. old_local_idx <= size(sc_batch, 1)) then
                  coef_seed_part(io, c0:c1, ispin, i_local) = &
                    coef_seed_part(io, c0:c1, ispin, i_local) + &
                    transform(jo, io, i_local) * sc_batch(old_local_idx, 1:nb)
                end if
              end do
            end do
          end do
          deallocate(sc_batch)
        end do
      end do
      call comm_summation(coef_seed_part, coef_seed, size(coef_seed_part), dg_frag%icomm)
      if (any(real(coef_seed) /= real(coef_seed)) .or. any(aimag(coef_seed) /= aimag(coef_seed))) then
        write(*,'(1x,a,i0)') "[FATAL] NaN in LCFO coefficient projection to contracted DG basis: rank=", dg_frag%id
        stop "DG-Fragment RT: NaN projected LCFO coefficients"
      end if
      deallocate(coef_seed_part)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,i0)') "[INFO] DG projected LCFO coefficients into contracted basis: nstate=", nstate_project
      end if
    end if

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

    if (allocated(dg_frag%coef)) deallocate(dg_frag%coef)
    if (allocated(dg_frag%coef_new)) deallocate(dg_frag%coef_new)
    if (allocated(dg_frag%coef_work)) deallocate(dg_frag%coef_work)
    dg_frag%nstate_tot = nstate_new
    call rebuild_coef_owner_map(dg_frag, "fragment-local-eigen-cap")
    nvalid = max(1, maxval(dg_frag%local_coef_count(1:dg_frag%nspin)))
    allocate(dg_frag%coef(nvalid, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(dg_frag%coef_work(nvalid, dg_frag%nstate_tot, dg_frag%nspin))
    if (dg_frag%yn_adaptive_basis) allocate(dg_frag%coef_new(nvalid, dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%coef(:, :, :) = (0.0d0, 0.0d0)
    dg_frag%coef_work(:, :, :) = (0.0d0, 0.0d0)
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = (0.0d0, 0.0d0)

    if (allocated(system%rocc)) system%rocc(:, :, :) = 0.0d0
    if (allocated(dg_frag%nocc_spin)) dg_frag%nocc_spin(:) = 0
    if (allocated(dg_frag%esp)) deallocate(dg_frag%esp)
    allocate(dg_frag%esp(dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%esp(:, :) = 0.0d0

    if (allocated(coef_seed)) then
      do ispin = 1, dg_frag%nspin
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = ifrag - dg_frag%ifrag_start + 1
          nkeep = min(max_keep, nkeep_global(ifrag, ispin))
          do io = 1, nkeep
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            local_idx = 0
            if (global_idx >= 1 .and. global_idx <= dg_frag%n_mat_max) local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
            if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
              dg_frag%coef(local_idx, 1:nstate_project, ispin) = &
                coef_seed(io, 1:nstate_project, ispin, i_local)
            end if
          end do
        end do
      end do
    end if

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
    if (allocated(dg_frag%coef_work)) dg_frag%coef_work(:, :, :) = dg_frag%coef(:, :, :)
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)

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
      write(*,'(1x,a,i0,a,i0,a,i0)') "[INFO] DG local eigen basis contracted: cap=", max_keep, &
        " n_mat=", dg_frag%n_mat_max, " nstate_tot=", dg_frag%nstate_tot
      if (max_keep /= nstate_cap) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[INFO] DG flux-relax basis kept above occupied cap: occ_cap=", &
          nstate_cap, " relax_keep=", max_keep, " file_basis=", nstate_old
      end if
      write(*,'(1x,a,1pe13.5,a,1pe13.5,a,1pe13.5)') "[INFO] DG local eigen FD: mu=", mu, &
        " Ne=", electron_count, " kT=", fd_temperature
    end if

    call comm_sync_all(dg_frag%icomm)
    deallocate(transform, nkeep_local, nkeep_global, nkeep_occ_local, nkeep_occ_global, eig_local, eig_global, occ_global)
    if (allocated(coef_seed)) deallocate(coef_seed)
  end subroutine prepare_fragment_local_eigen_basis

end module rt_dg_local_basis
