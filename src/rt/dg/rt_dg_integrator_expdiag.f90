  subroutine time_evolution_expdiag(dg_frag, system, info, rt, itt, dt, &
                                    lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                    rho, rho_s, Vh, Vxc, Vpsl, energy)
    use structures
    use salmon_global, only: yn_fix_func, yn_dg_length_gauge, ae_shape1, e_impulse, epdir_re1, yn_restart
    use sendrecv_grid, only: s_sendrecv_grid
    use salmon_xc, only: s_xc_functional
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A
    use eigen_subdiag_sub, only: eigen_dsyev
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(inout) :: system
    type(s_parallel_info),  intent(in)    :: info
    type(s_rt),             intent(inout) :: rt
    integer,                intent(in)    :: itt
    real(8),                intent(in)    :: dt
    type(s_rgrid),          intent(in)    :: lg, mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_xc_functional),  intent(in)    :: xc_func
    type(s_sendrecv_grid),  intent(inout) :: srg, srg_scalar
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson),        intent(inout) :: poisson
    type(s_pp_info),        intent(in)    :: pp
    type(s_pp_grid),        intent(in)    :: ppg
    type(s_pp_nlcc),        intent(in)    :: ppn
    type(s_scalar),         intent(inout) :: rho, Vh, Vpsl
    type(s_scalar),         intent(inout) :: rho_s(system%nspin), Vxc(system%nspin)
    type(s_dft_energy),     intent(inout) :: energy

    integer :: it0, it1, ispin, ifrag, i_local, iblk
    integer :: nw, nbf, nstate_prop, state_first, state_last
    integer :: io, jo, iw, jw, istate, global_idx, local_idx, local_col
    real(8) :: Ac_ham(3), E_mid(3), rdot_diag
    real(8) :: cphase, sphase
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    complex(8), allocatable :: c_w(:), tmp_w(:), next_w(:)
    real(8), allocatable :: h_eff(:,:), h_dg(:,:), hw_tmp(:,:), eval(:), evec(:,:)
    logical, save :: expdiag_warned = .false.
    logical, save :: xi_split_env_checked = .false.
    logical, save :: xi_split_enabled = .true.
    logical, save :: project_h_env_checked = .false.
    logical, save :: project_h_for_fixed_func = .false.
    character(32) :: xi_split_env
    character(32) :: project_h_env

    if (yn_dg_length_gauge /= 'y') stop "DG expdiag integrator currently requires length gauge"
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw)) &
      stop "DG expdiag integrator does not yet support BPW+PW mixed propagation"
    if (dg_frag%coef_state_block_mode) &
      stop "DG expdiag integrator does not yet support state-block coefficient ownership"
    if (.not. dg_frag%buffer_wannier_flux_seed_applied .or. &
        .not. dg_frag%has_buffer_periodic_wannier_basis .or. &
        .not. allocated(dg_frag%buffer_wannier_coef) .or. &
        .not. allocated(dg_frag%buffer_wannier_h_flux) .or. &
        .not. allocated(dg_frag%buffer_wannier_v)) then
      stop "DG expdiag integrator requires buffer-periodic Wannier seed data"
    end if
    if (.not. xi_split_env_checked) then
      xi_split_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_XI_SPLIT', xi_split_env)
      select case (trim(adjustl(xi_split_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        xi_split_enabled = .true.
      case ('0','n','N','no','NO','false','FALSE','off','OFF')
        xi_split_enabled = .false.
      case default
        xi_split_enabled = .true.
      end select
      xi_split_env_checked = .true.
    end if
    if (.not. project_h_env_checked) then
      project_h_env = ' '
      call get_environment_variable('SALMON_DG_EXPDIAG_PROJECT_H', project_h_env)
      select case (trim(adjustl(project_h_env)))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        project_h_for_fixed_func = .true.
      case default
        project_h_for_fixed_func = .false.
      end select
      project_h_env_checked = .true.
    end if

    it0 = max(lbound(rt%Ac_tot, 2), itt - 1)
    it1 = min(ubound(rt%Ac_tot, 2), itt)
    Ac_ham(:) = 0.0d0
    E_mid(:) = 0.5d0 * (rt%E_tot(:, it0) + rt%E_tot(:, it1))
    if (trim(ae_shape1) == 'impulse' .and. yn_restart == 'n' .and. itt == 1) then
      E_mid(:) = e_impulse * epdir_re1(:) / max(abs(dt), 1.0d-300)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,3(1x,1pe13.5),a,1pe13.5)') &
          '[DG-LG-IMPULSE] using first-step rectangular E field=', E_mid(:), ' dt=', dt
        flush(6)
      end if
    end if

    if (yn_fix_func == 'n') then
      call update_density_hamiltonian_stage(dg_frag, system, info, rt, itt, Ac_ham, &
                                            lg, mg, stencil, xc_func, srg, srg_scalar, fg, poisson, pp, ppg, ppn, &
                                            rho, rho_s, Vh, Vxc, Vpsl, energy, .false.)
    end if
    if (ppg%Nlma > 0 .and. allocated(ppg%uV)) then
      call ensure_nonlocal_pp_matrix_A(dg_frag, mg, ppg, system, Ac_ham, .false.)
    end if

    nstate_prop = dg_frag%nstate_tot
    if (allocated(dg_frag%nocc_spin)) then
      nstate_prop = min(dg_frag%nstate_tot, max(1, maxval(dg_frag%nocc_spin(1:dg_frag%nspin))))
    end if
    state_first = 1
    state_last = min(nstate_prop, size(dg_frag%coef, 2))

    if (xi_split_enabled) call apply_xi_flux_split_half(E_mid, 0.5d0 * dt, state_first, state_last)
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = ifrag - dg_frag%ifrag_start + 1
      if (i_local < 1 .or. i_local > size(dg_frag%buffer_wannier_nkeep)) cycle
      nw = dg_frag%buffer_wannier_nkeep(i_local)
      nbf = min(dg_frag%n_basis(ifrag, 1), dg_frag%nstate_frag, size(dg_frag%buffer_wannier_coef, 1))
      if (nw <= 0 .or. nbf <= 0) cycle

      allocate(h_eff(nw,nw), h_dg(nbf,nbf), hw_tmp(nbf,nw), eval(nw), evec(nw,nw), &
               c_w(nw), tmp_w(nw), next_w(nw))
      do ispin = 1, dg_frag%nspin
        h_eff(1:nw,1:nw) = dg_frag%buffer_wannier_h_flux(1:nw,1:nw,i_local)
        iblk = 0
        if (allocated(dg_frag%H_block_map)) iblk = find_matrix_block(dg_frag%H_block_map, ifrag, ifrag)
        if ((yn_fix_func == 'n' .or. project_h_for_fixed_func) .and. iblk > 0 .and. allocated(dg_frag%H_mat_blocks)) then
          if (iblk <= size(dg_frag%H_mat_blocks)) then
            h_dg(1:nbf,1:nbf) = dg_frag%H_mat_blocks(iblk)%val(1:nbf,1:nbf,ispin)
            hw_tmp(1:nbf,1:nw) = matmul(h_dg(1:nbf,1:nbf), &
              dg_frag%buffer_wannier_coef(1:nbf,1:nw,ispin,i_local))
            h_eff(1:nw,1:nw) = matmul(transpose(dg_frag%buffer_wannier_coef(1:nbf,1:nw,ispin,i_local)), &
              hw_tmp(1:nbf,1:nw))
          end if
        end if
        do iw = 1, nw
          do jw = 1, nw
            h_eff(jw,iw) = h_eff(jw,iw) &
              - E_mid(1) * dg_frag%buffer_wannier_v(1,jw,iw,i_local) &
              - E_mid(2) * dg_frag%buffer_wannier_v(2,jw,iw,i_local) &
              - E_mid(3) * dg_frag%buffer_wannier_v(3,jw,iw,i_local)
          end do
          rdot_diag = 0.0d0
          if (allocated(dg_frag%buffer_wannier_frag_center)) then
            if (i_local >= 1 .and. i_local <= size(dg_frag%buffer_wannier_frag_center, 2)) then
              rdot_diag = E_mid(1) * dg_frag%buffer_wannier_frag_center(1,i_local) &
                        + E_mid(2) * dg_frag%buffer_wannier_frag_center(2,i_local) &
                        + E_mid(3) * dg_frag%buffer_wannier_frag_center(3,i_local)
            end if
          end if
          h_eff(iw,iw) = h_eff(iw,iw) - rdot_diag
        end do
        evec(1:nw,1:nw) = h_eff(1:nw,1:nw)
        call eigen_dsyev(evec, eval, h_eff)

        do istate = state_first, state_last
          c_w(:) = (0.0d0, 0.0d0)
          do iw = 1, nw
            do io = 1, nbf
              global_idx = dg_frag%index_basis(io, ifrag, ispin)
              local_idx = 0
              if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) &
                local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
              if (local_idx <= 0 .or. local_idx > size(dg_frag%coef, 1)) cycle
              c_w(iw) = c_w(iw) + dg_frag%buffer_wannier_coef(io, iw, ispin, i_local) * &
                                    dg_frag%coef(local_idx, istate, ispin)
            end do
          end do

          tmp_w(:) = matmul(transpose(h_eff(1:nw,1:nw)), c_w(:))
          do iw = 1, nw
            cphase = cos(eval(iw) * dt)
            sphase = sin(eval(iw) * dt)
            tmp_w(iw) = cmplx(cphase, -sphase, kind=8) * tmp_w(iw)
          end do
          next_w(:) = matmul(h_eff(1:nw,1:nw), tmp_w(:))

          do io = 1, nbf
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            local_idx = 0
            if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) &
              local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
            if (local_idx <= 0 .or. local_idx > size(dg_frag%coef, 1)) cycle
            dg_frag%coef(local_idx, istate, ispin) = (0.0d0, 0.0d0)
            do iw = 1, nw
              dg_frag%coef(local_idx, istate, ispin) = dg_frag%coef(local_idx, istate, ispin) + &
                dg_frag%buffer_wannier_coef(io, iw, ispin, i_local) * next_w(iw)
            end do
          end do
        end do
      end do
      deallocate(h_eff, h_dg, hw_tmp, eval, evec, c_w, tmp_w, next_w)
    end do
    if (xi_split_enabled) call apply_xi_flux_split_half(E_mid, 0.5d0 * dt, state_first, state_last)

    if (.not. expdiag_warned .and. dg_frag%id == 0) then
      write(*,'(1x,a)') "[DG-EXPDIAG] experimental local BPW exponential integrator enabled."
      if (xi_split_enabled) then
        write(*,'(1x,a)') "[DG-EXPDIAG] dynamic neighbor xi_flux is applied by a Strang half-step split."
      else
        write(*,'(1x,a)') &
          "[DG-EXPDIAG] dynamic neighbor xi_flux split disabled; set SALMON_DG_EXPDIAG_XI_SPLIT=1 to test."
      end if
      if (yn_fix_func == 'y' .and. .not. project_h_for_fixed_func) then
        write(*,'(1x,a)') "[DG-EXPDIAG] fixed-function propagation uses the seed flux Hamiltonian."
      else
        write(*,'(1x,a)') "[DG-EXPDIAG] propagation projects the current DG Hamiltonian into the BPW basis."
      end if
      flush(6)
      expdiag_warned = .true.
    end if

  contains

    subroutine apply_xi_flux_split_half(E_use, dt_half, state_s, state_e)
      real(8), intent(in) :: E_use(3), dt_half
      integer, intent(in) :: state_s, state_e
      integer :: ispin_current, iblk_idx, iblk, ifrag_row, ifrag_col
      integer :: nstate_blk, state_col, io, jo, istate
      integer :: nrow, ncol, row_gid, col_gid, row_local, col_pos
      integer :: nfetched
      real(8) :: val
      complex(8), allocatable :: fetched(:,:), delta(:,:)
      integer, allocatable :: fetched_ids(:)

      if (sum(abs(E_use(1:3))) <= 1.0d-30) return
      if (abs(dt_half) <= 0.0d0) return
      nstate_blk = max(0, state_e - state_s + 1)
      if (nstate_blk <= 0) return

      if (.not. dg_frag%buffer_wannier_xi_flux_available .or. &
          .not. allocated(dg_frag%buffer_wannier_xi_flux_blocks) .or. &
          .not. allocated(dg_frag%buffer_wannier_xi_flux_local_block_ids)) then
        stop "DG expdiag requires neighbor xi_flux blocks"
      end if

      do ispin_current = 1, dg_frag%nspin
        call exchange_xi_flux_neighbor_coef(ispin_current, state_s, state_e, fetched_ids, fetched, nfetched)
        allocate(delta(size(dg_frag%coef,1),nstate_blk))
        delta = (0.0d0, 0.0d0)
        do iblk_idx = 1, size(dg_frag%buffer_wannier_xi_flux_local_block_ids)
          iblk = dg_frag%buffer_wannier_xi_flux_local_block_ids(iblk_idx)
          if (iblk < 1 .or. iblk > size(dg_frag%buffer_wannier_xi_flux_blocks)) cycle
          ifrag_row = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_row
          ifrag_col = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_col
          if (ifrag_row == ifrag_col) cycle
          if (.not. allocated(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val)) cycle
          nrow = min(dg_frag%n_basis(ifrag_row, ispin_current), &
                     size(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val, 2))
          ncol = min(dg_frag%n_basis(ifrag_col, ispin_current), &
                     size(dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val, 3))
          do io = 1, nrow
            row_gid = dg_frag%index_basis(io, ifrag_row, ispin_current)
            if (row_gid < 1 .or. row_gid > dg_frag%n_mat_max) cycle
            row_local = dg_frag%coef_global_to_local(row_gid, ispin_current)
            if (row_local <= 0 .or. row_local > size(dg_frag%coef, 1)) cycle
            do jo = 1, ncol
              col_gid = dg_frag%index_basis(jo, ifrag_col, ispin_current)
              col_pos = find_needed_row_pos(col_gid, fetched_ids, nfetched)
              if (col_pos <= 0) cycle
              val = E_use(1) * dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val(1,io,jo,ispin_current) &
                  + E_use(2) * dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val(2,io,jo,ispin_current) &
                  + E_use(3) * dg_frag%buffer_wannier_xi_flux_blocks(iblk)%val(3,io,jo,ispin_current)
              if (abs(val) <= 0.0d0) cycle
              do istate = 1, nstate_blk
                delta(row_local,istate) = delta(row_local,istate) - zi * dt_half * val * fetched(col_pos,istate)
              end do
            end do
          end do
        end do
        do istate = 1, nstate_blk
          state_col = state_s + istate - 1
          if (state_col < 1 .or. state_col > size(dg_frag%coef, 2)) cycle
          dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) = &
            dg_frag%coef(1:size(dg_frag%coef,1),state_col,ispin_current) + &
            delta(1:size(dg_frag%coef,1),istate)
        end do
        if (allocated(fetched_ids)) deallocate(fetched_ids)
        if (allocated(fetched)) deallocate(fetched)
        deallocate(delta)
      end do
    end subroutine apply_xi_flux_split_half

    subroutine exchange_xi_flux_neighbor_coef(ispin_use, state_s, state_e, fetched_ids, fetched_coef, nfetched)
      use mpi, only: MPI_Isend, MPI_Irecv, MPI_Waitall, MPI_INTEGER, MPI_DOUBLE_COMPLEX, &
                     MPI_SUCCESS, MPI_STATUSES_IGNORE
      integer, intent(in) :: ispin_use, state_s, state_e
      integer, allocatable, intent(out) :: fetched_ids(:)
      complex(8), allocatable, intent(out) :: fetched_coef(:,:)
      integer, intent(out) :: nfetched

      integer :: nstate_blk, peer, npeer, iblk, iblk_idx, ifrag_row, ifrag_col
      integer :: i, ist, local_idx, gid, ierr, nreq, total_recv, total_send
      integer :: data_pos, row_pos, peer_owner
      integer, allocatable :: peers(:), send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
      integer, allocatable :: requests(:), send_ids(:), recv_ids(:)
      complex(8), allocatable :: send_buf(:), recv_buf(:)

      nfetched = 0
      nstate_blk = max(0, state_e - state_s + 1)
      allocate(fetched_ids(1), fetched_coef(1,max(1,nstate_blk)))
      fetched_ids(:) = 0
      fetched_coef(:,:) = (0.0d0, 0.0d0)
      if (nstate_blk <= 0) return
      if (.not. allocated(dg_frag%coef)) return
      if (.not. allocated(dg_frag%coef_global_to_local)) return
      if (.not. allocated(dg_frag%coef_owner)) return
      if (.not. allocated(dg_frag%local_coef_count)) return
      if (.not. allocated(dg_frag%local_coef_global_ids)) return
      if (.not. allocated(dg_frag%id_array)) return

      allocate(peers(max(1, dg_frag%isize)))
      peers(:) = -1
      npeer = 0
      if (allocated(dg_frag%buffer_wannier_xi_flux_blocks) .and. &
          allocated(dg_frag%buffer_wannier_xi_flux_local_block_ids)) then
        do iblk_idx = 1, size(dg_frag%buffer_wannier_xi_flux_local_block_ids)
          iblk = dg_frag%buffer_wannier_xi_flux_local_block_ids(iblk_idx)
          if (iblk < 1 .or. iblk > size(dg_frag%buffer_wannier_xi_flux_blocks)) cycle
          ifrag_row = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_row
          ifrag_col = dg_frag%buffer_wannier_xi_flux_blocks(iblk)%ifrag_col
          if (ifrag_row < 1 .or. ifrag_row > dg_frag%n_frag) cycle
          if (ifrag_col < 1 .or. ifrag_col > dg_frag%n_frag) cycle
          if (ifrag_row == ifrag_col) cycle
          if (fragment_is_local(ifrag_row)) then
            peer_owner = dg_frag%id_array(ifrag_col)
            call add_exchange_peer(peer_owner, peers, npeer)
          end if
          if (fragment_is_local(ifrag_col)) then
            peer_owner = dg_frag%id_array(ifrag_row)
            call add_exchange_peer(peer_owner, peers, npeer)
          end if
        end do
      end if

      if (npeer <= 0) then
        deallocate(peers)
        return
      end if

      total_send = 0
      if (ispin_use >= 1 .and. ispin_use <= size(dg_frag%local_coef_count)) &
        total_send = max(0, dg_frag%local_coef_count(ispin_use))

      allocate(send_counts(npeer), recv_counts(npeer), send_displs(npeer), recv_displs(npeer))
      allocate(requests(max(1, 2*npeer)))
      send_counts(:) = total_send
      recv_counts(:) = 0
      nreq = 0
      do i = 1, npeer
        peer = peers(i)
        nreq = nreq + 1
        call MPI_Irecv(recv_counts(i), 1, MPI_INTEGER, peer, 8521, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux count recv failed"
      end do
      do i = 1, npeer
        peer = peers(i)
        nreq = nreq + 1
        call MPI_Isend(send_counts(i), 1, MPI_INTEGER, peer, 8521, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux count send failed"
      end do
      call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
      if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux count wait failed"

      send_displs(1) = 0
      recv_displs(1) = 0
      do i = 2, npeer
        send_displs(i) = send_displs(i-1) + send_counts(i-1)
        recv_displs(i) = recv_displs(i-1) + recv_counts(i-1)
      end do
      total_recv = recv_displs(npeer) + recv_counts(npeer)

      allocate(send_ids(max(1,total_send)), recv_ids(max(1,total_recv)))
      send_ids(:) = 0
      recv_ids(:) = 0
      if (total_send > 0) send_ids(1:total_send) = dg_frag%local_coef_global_ids(1:total_send, ispin_use)

      nreq = 0
      do i = 1, npeer
        peer = peers(i)
        if (recv_counts(i) <= 0) cycle
        nreq = nreq + 1
        call MPI_Irecv(recv_ids(recv_displs(i)+1), recv_counts(i), MPI_INTEGER, peer, 8522, dg_frag%icomm, &
                       requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux row recv failed"
      end do
      do i = 1, npeer
        peer = peers(i)
        if (send_counts(i) <= 0) cycle
        nreq = nreq + 1
        call MPI_Isend(send_ids(1), send_counts(i), MPI_INTEGER, peer, 8522, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux row send failed"
      end do
      if (nreq > 0) then
        call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux row wait failed"
      end if

      allocate(send_buf(max(1,total_send*nstate_blk)), recv_buf(max(1,total_recv*nstate_blk)))
      send_buf(:) = (0.0d0, 0.0d0)
      recv_buf(:) = (0.0d0, 0.0d0)
      do ist = 1, nstate_blk
        do i = 1, total_send
          gid = send_ids(i)
          if (gid < 1 .or. gid > size(dg_frag%coef_global_to_local, 1)) cycle
          local_idx = dg_frag%coef_global_to_local(gid, ispin_use)
          if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
          if (state_s + ist - 1 < 1 .or. state_s + ist - 1 > size(dg_frag%coef, 2)) cycle
          send_buf((ist-1)*total_send + i) = dg_frag%coef(local_idx, state_s + ist - 1, ispin_use)
        end do
      end do

      nreq = 0
      do i = 1, npeer
        peer = peers(i)
        if (recv_counts(i) <= 0) cycle
        nreq = nreq + 1
        call MPI_Irecv(recv_buf(recv_displs(i)*nstate_blk + 1), recv_counts(i)*nstate_blk, &
                       MPI_DOUBLE_COMPLEX, peer, 8523, dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux data recv failed"
      end do
      do i = 1, npeer
        peer = peers(i)
        if (send_counts(i) <= 0) cycle
        nreq = nreq + 1
        call MPI_Isend(send_buf(1), send_counts(i)*nstate_blk, MPI_DOUBLE_COMPLEX, peer, 8523, &
                       dg_frag%icomm, requests(nreq), ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux data send failed"
      end do
      if (nreq > 0) then
        call MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE, ierr)
        if (ierr /= MPI_SUCCESS) stop "DG expdiag xi_flux data wait failed"
      end if

      deallocate(fetched_ids, fetched_coef)
      nfetched = total_recv
      allocate(fetched_ids(max(1,nfetched)), fetched_coef(max(1,nfetched),nstate_blk))
      fetched_ids(:) = 0
      fetched_coef(:,:) = (0.0d0, 0.0d0)
      if (nfetched > 0) then
        fetched_ids(1:nfetched) = recv_ids(1:nfetched)
        do peer = 1, npeer
          do ist = 1, nstate_blk
            do i = 1, recv_counts(peer)
              row_pos = recv_displs(peer) + i
              data_pos = recv_displs(peer) * nstate_blk + (ist - 1) * recv_counts(peer) + i
              fetched_coef(row_pos,ist) = recv_buf(data_pos)
            end do
          end do
        end do
      end if

      deallocate(peers, send_counts, recv_counts, send_displs, recv_displs, requests)
      deallocate(send_ids, recv_ids, send_buf, recv_buf)
    end subroutine exchange_xi_flux_neighbor_coef

    logical function fragment_is_local(ifrag_use) result(is_local)
      integer, intent(in) :: ifrag_use
      is_local = (ifrag_use >= dg_frag%ifrag_start .and. ifrag_use <= dg_frag%ifrag_end)
    end function fragment_is_local

    subroutine add_exchange_peer(peer, peers, npeer)
      integer, intent(in) :: peer
      integer, intent(inout) :: peers(:)
      integer, intent(inout) :: npeer
      integer :: i
      if (peer < 0 .or. peer >= dg_frag%isize .or. peer == dg_frag%id) return
      do i = 1, npeer
        if (peers(i) == peer) return
      end do
      if (npeer >= size(peers)) return
      npeer = npeer + 1
      peers(npeer) = peer
    end subroutine add_exchange_peer

    integer function find_needed_row_pos(gid, row_ids, nneeded) result(pos)
      integer, intent(in) :: gid, nneeded
      integer, intent(in) :: row_ids(:)
      integer :: i

      pos = 0
      if (gid <= 0) return
      do i = 1, nneeded
        if (row_ids(i) == gid) then
          pos = i
          return
        end if
      end do
    end function find_needed_row_pos
  end subroutine time_evolution_expdiag
