  subroutine ensure_mixed_wannier_bpw_position(dg_frag)
    use communication, only: comm_is_root, comm_summation
    use eigen_subdiag_sub, only: eigen_zheev
    use rt_dg_plane_wave, only: compute_fragment_pw_overlap, compute_fragment_pw_position_overlap
    use salmon_global, only: dg_bpw_auto, dg_bpw_auto_accuracy, dg_bpw_auto_max_n, &
      dg_bpw_auto_min_n, dg_bpw_auto_report, n_plane_waves_dg, yn_dg_mixed_z_include_ww
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    complex(8), parameter :: zzero = (0.0d0, 0.0d0)
    real(8) :: sperp_tol
    integer :: nmat, npw, nwann, neig, nspin_file, nocc_ref
    integer :: ispin, idir, i_local, ifrag, ibasis, iw, jw, ipw, alpha, beta, eig
    integer :: imix, jmix, jdir, virt, total_dim
    integer :: nbf, global_row, local_row, nkeep, pslot, info_diag
    real(8) :: min_s_local, max_s_local, pw_energy_ref, gap, herm_max
    real(8) :: comm_rel_tol, comm_abs_tol, comm_min_fwp_floor
    real(8) :: bpw_ecut, shell_energy, selected_max_energy
    real(8) :: comm_old, comm_new, fwp_old, fwp_new, min_s_trial, max_s_trial
    real(8) :: best_comm, best_fwp
    real(8) :: fwp_axis_old(3), fwp_axis_new(3), fwp_axis_final(3)
    real(8) :: comm_pair_old(3), comm_pair_new(3), comm_pair_final(3)
    real(8) :: zwp_direct_norm(3)
    real(8) :: zwp_norm(3), zww_norm(3)
    real(8) :: zcov_ww(3,3), zcov_wp(3,3), zcov_tot(3,3)
    real(8) :: zcomm_ww(3,3), zcomm_wp(3,3), zcomm_tot(3,3)
    real(8) :: metric_sww_max, metric_sww_rms, metric_swz_max, metric_swz_rms
    real(8) :: metric_szz_max, metric_szz_rms, metric_stotal_max, metric_stotal_rms
    real(8) :: h_total_herm_max, h_total_herm_rms
    real(8) :: frag_center(3)
    real(8) :: sv_min_selected, sv_max_selected, sv_condition
    real(8) :: proj_min_selected, proj_max_rejected, proj_gap
    real(8) :: fsum_proxy_wp_avg, fsum_proxy_total_avg
    complex(8) :: coeff, amp, grad_corr
    character(len=32) :: env_sperp, env_direct_origin, env_h_evec
    character(len=32) :: env_legacy_h_vec
    character(len=32) :: env_select_mode, env_comm
    character(len=32) :: bpw_select_mode
    logical :: use_fragment_center_direct, use_comm_min_select
    logical :: use_shell_ecut_select, use_shell_nraw_select
    logical :: use_h_eigenvectors_for_final, use_legacy_h_vec_for_final
    logical :: bpw_ecut_is_set, comm_diag_enabled, fsum_diag_enabled
    integer :: env_len, env_stat
    integer :: target_nraw, max_shells, nshell, ishell, shell_id, shell_limit
    integer :: nraw_old, nraw_new, np_old, np_new, best_nraw, best_np
    integer :: selected_nshell, selected_g2_max, shell_count, accepted_by_ecut_int
    integer :: nraw_selected, report_unit, report_iostat
    logical :: auto_report_enabled
    complex(8), allocatable :: s_fp(:,:,:), r_fp(:,:,:,:)
    complex(8), allocatable :: b_local(:,:), b_global(:,:)
    complex(8), allocatable :: r_local(:,:), r_global(:,:), b_eig(:,:)
    complex(8), allocatable :: c_local(:,:), c_global(:,:), c_w_local(:,:), c_w_global(:,:), c_eig(:,:)
    complex(8), allocatable :: z_w(:,:), z_eig(:,:), o_eig(:,:)
    complex(8), allocatable :: s_perp(:,:), s_vec(:,:), s_work(:,:), qmat(:,:)
    complex(8), allocatable :: h_p(:,:), h_vec(:,:)
    complex(8), allocatable :: q_metric(:,:), q_left_metric(:,:), final_metric(:,:), transform_metric(:,:)
    complex(8), allocatable :: t_final(:,:), t_saved(:,:), metric_vec(:,:)
    complex(8), allocatable :: total_metric(:,:), h_total(:,:)
    complex(8), allocatable :: r_tilde(:,:), r_orth(:,:), r_final(:,:), z_wp_direct(:,:)
    real(8), allocatable :: s_eval(:), h_eval(:), pw_ekin(:), shell_ekin(:), diag_eval(:)
    logical, allocatable :: active_pw(:), trial_active_pw(:), best_active_pw(:)
    integer, allocatable :: shell_key(:), shell_values(:)

    if (dg_frag%has_mixed_wannier_bpw_position) return
    if (.not. dg_frag%use_plane_wave_basis .or. dg_frag%n_plane_waves <= 0) return
    if (.not. dg_frag%has_global_wannier_flux_eigen) return
    if (.not. dg_frag%has_global_wannier_position .or. .not. allocated(dg_frag%global_wannier_position)) return
    if (.not. allocated(dg_frag%global_wannier_coef) .or. .not. allocated(dg_frag%global_wannier_flux_evec) .or. &
        .not. allocated(dg_frag%global_wannier_flux_eval)) return

    nmat = dg_frag%n_mat_max
    npw = dg_frag%n_plane_waves
    nwann = size(dg_frag%global_wannier_flux_evec, 1)
    neig = min(size(dg_frag%global_wannier_flux_evec, 2), size(dg_frag%global_wannier_flux_eval, 1))
    nspin_file = min(dg_frag%nspin, size(dg_frag%global_wannier_flux_eval, 2))
    if (nmat <= 0 .or. npw <= 0 .or. nwann <= 0 .or. neig <= 0 .or. nspin_file <= 0) return
    auto_report_enabled = (dg_bpw_auto == 'y' .or. dg_bpw_auto == 'Y') .and. &
      (dg_bpw_auto_report /= 'n' .and. dg_bpw_auto_report /= 'N')

    sperp_tol = 1.0d-8
    env_sperp = ''
    call get_environment_variable('SALMON_DG_MIXED_SPERP_TOL', env_sperp, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_sperp(1:env_len), *, iostat=env_stat) sperp_tol
      if (env_stat /= 0 .or. sperp_tol < 0.0d0) sperp_tol = 1.0d-8
    end if

    use_fragment_center_direct = .false.
    env_direct_origin = ''
    call get_environment_variable('SALMON_DG_MIXED_Z_DIRECT_ORIGIN', env_direct_origin, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_direct_origin(1:env_len))))
      case ('fragment','frag','local','center','FRAGMENT','LOCAL')
        use_fragment_center_direct = .true.
      case default
        use_fragment_center_direct = .false.
      end select
    end if
    use_h_eigenvectors_for_final = .true.
    use_legacy_h_vec_for_final = .false.
    env_legacy_h_vec = ''
    call get_environment_variable('SALMON_DG_MIXED_Z_LEGACY_H_VEC', &
      env_legacy_h_vec, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_legacy_h_vec(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        use_legacy_h_vec_for_final = .true.
      end select
    end if
    env_h_evec = ''
    call get_environment_variable('SALMON_DG_MIXED_Z_USE_H_EIGENVEC', env_h_evec, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_h_evec(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        use_h_eigenvectors_for_final = .true.
      end select
    end if
    if (use_legacy_h_vec_for_final) use_h_eigenvectors_for_final = .false.

    bpw_select_mode = 'shell_ecut'
    use_comm_min_select = .false.
    use_shell_ecut_select = .true.
    use_shell_nraw_select = .false.
    env_select_mode = ''
    call get_environment_variable('SALMON_DG_BPW_SELECT_MODE', env_select_mode, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_select_mode(1:env_len))))
      case ('comm_min','COMM_MIN','comm-min','comm')
        bpw_select_mode = 'comm_min'
        use_comm_min_select = .true.
        use_shell_ecut_select = .false.
        use_shell_nraw_select = .false.
      case ('shell_nraw','SHELL_NRAW','shell-nraw','nraw')
        bpw_select_mode = 'shell_nraw'
        use_comm_min_select = .false.
        use_shell_ecut_select = .false.
        use_shell_nraw_select = .true.
      case ('shell_ecut','SHELL_ECUT','shell-ecut','ecut')
        bpw_select_mode = 'shell_ecut'
        use_comm_min_select = .false.
        use_shell_ecut_select = .true.
        use_shell_nraw_select = .false.
      case default
        bpw_select_mode = 'shell_ecut'
        use_comm_min_select = .false.
        use_shell_ecut_select = .true.
        use_shell_nraw_select = .false.
      end select
    end if
    bpw_ecut = huge(1.0d0)
    bpw_ecut_is_set = .false.
    env_comm = ''
    call get_environment_variable('SALMON_DG_BPW_ECUT', env_comm, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_comm(1:env_len), *, iostat=env_stat) bpw_ecut
      if (env_stat == 0 .and. bpw_ecut >= 0.0d0) then
        bpw_ecut_is_set = .true.
      else
        bpw_ecut = huge(1.0d0)
        bpw_ecut_is_set = .false.
      end if
    end if
    comm_diag_enabled = .false.
    env_comm = ''
    call get_environment_variable('SALMON_DG_BPW_COMM_DIAG', env_comm, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_comm(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        comm_diag_enabled = .true.
      case default
        comm_diag_enabled = .false.
      end select
    end if
    fsum_diag_enabled = .false.
    env_comm = ''
    call get_environment_variable('SALMON_DG_BPW_FSUM_DIAG', env_comm, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_comm(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        fsum_diag_enabled = .true.
      case default
        fsum_diag_enabled = .false.
      end select
    end if
    comm_rel_tol = 2.5d-1
    env_comm = ''
    call get_environment_variable('SALMON_DG_BPW_COMM_REL_TOL', env_comm, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_comm(1:env_len), *, iostat=env_stat) comm_rel_tol
      if (env_stat /= 0 .or. comm_rel_tol < 0.0d0) comm_rel_tol = 2.5d-1
    end if
    comm_abs_tol = 2.5d-1
    env_comm = ''
    call get_environment_variable('SALMON_DG_BPW_COMM_ABS_TOL', env_comm, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_comm(1:env_len), *, iostat=env_stat) comm_abs_tol
      if (env_stat /= 0 .or. comm_abs_tol < 0.0d0) comm_abs_tol = 2.5d-1
    end if
    comm_min_fwp_floor = 1.0d-10
    env_comm = ''
    call get_environment_variable('SALMON_DG_BPW_COMM_MIN_FWP_FLOOR', env_comm, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_comm(1:env_len), *, iostat=env_stat) comm_min_fwp_floor
      if (env_stat /= 0 .or. comm_min_fwp_floor < 0.0d0) comm_min_fwp_floor = 1.0d-10
    end if
    target_nraw = npw
    env_comm = ''
    call get_environment_variable('SALMON_DG_BPW_TARGET_NRAW', env_comm, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_comm(1:env_len), *, iostat=env_stat) target_nraw
      if (env_stat /= 0 .or. target_nraw <= 0) target_nraw = npw
    end if
    target_nraw = min(target_nraw, npw)
    max_shells = huge(1)
    env_comm = ''
    call get_environment_variable('SALMON_DG_BPW_MAX_SHELLS', env_comm, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_comm(1:env_len), *, iostat=env_stat) max_shells
      if (env_stat /= 0 .or. max_shells <= 0) max_shells = huge(1)
    end if

    allocate(s_fp(nmat,npw,dg_frag%nspin))
    allocate(r_fp(3,nmat,npw,dg_frag%nspin))
    call compute_fragment_pw_overlap(dg_frag, s_fp)
    call compute_fragment_pw_position_overlap(dg_frag, r_fp)

    allocate(b_local(nwann,npw), b_global(nwann,npw))
    allocate(r_local(nwann,npw), r_global(nwann,npw), b_eig(neig,npw))
    allocate(c_local(nwann,npw), c_global(nwann,npw), c_w_local(nwann,nwann), c_w_global(nwann,nwann), c_eig(neig,neig))
    allocate(z_w(nwann,nwann), z_eig(neig,neig), o_eig(neig,neig))
    allocate(s_perp(npw,npw), s_vec(npw,npw), s_work(npw,npw), s_eval(npw), pw_ekin(npw), shell_ekin(npw))
    allocate(qmat(npw,npw), h_p(npw,npw), h_vec(npw,npw), h_eval(npw))
    allocate(q_metric(npw,npw), q_left_metric(npw,npw), final_metric(npw,npw), transform_metric(npw,npw))
    allocate(t_final(npw,npw), t_saved(npw,npw), metric_vec(npw,npw), diag_eval(npw))
    allocate(total_metric(neig+npw,neig+npw), h_total(neig+npw,neig+npw))
    allocate(r_tilde(neig,npw), r_orth(neig,npw), r_final(neig,npw))
    allocate(z_wp_direct(neig,npw))
    allocate(active_pw(npw), trial_active_pw(npw), best_active_pw(npw), shell_key(npw), shell_values(npw))

    pw_ekin(:) = 0.0d0
    do ipw = 1, npw
      if (allocated(dg_frag%H_mat_pw_diag)) then
        pw_ekin(ipw) = dg_frag%H_mat_pw_diag(ipw, 1)
      else
        pw_ekin(ipw) = 0.5d0 * sum(dg_frag%k_pw(1:3,ipw)**2)
      end if
    end do

    call build_bpw_shells(nshell)
    shell_limit = min(nshell, max_shells)
    active_pw(:) = .false.
    selected_nshell = 0
    selected_g2_max = 0
    selected_max_energy = 0.0d0

    if (use_shell_ecut_select) then
      do ishell = 1, shell_limit
        if ((.not. bpw_ecut_is_set) .or. shell_ekin(ishell) <= bpw_ecut) then
          shell_id = shell_values(ishell)
          do ipw = 1, npw
            if (shell_key(ipw) == shell_id) active_pw(ipw) = .true.
          end do
          selected_nshell = ishell
          selected_g2_max = shell_id
          selected_max_energy = shell_ekin(ishell)
        end if
      end do
      if (comm_is_root(dg_frag%id)) then
        if (bpw_ecut_is_set) then
          write(*,'(1x,a,a,a,1pe13.5,a,i0,a,i0,a,i0,a,1pe13.5)') &
            "[DG-BPW-SELECT] mode=", trim(bpw_select_mode), " Ecut=", bpw_ecut, &
            " selected_shells=", selected_nshell, " selected_nraw=", count(active_pw), &
            " selected_G2_max=", selected_g2_max, " selected_Emax=", selected_max_energy
        else
          write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,1pe13.5)') &
            "[DG-BPW-SELECT] mode=", trim(bpw_select_mode), " Ecut=all_raw", &
            selected_nshell, " selected_nraw=", count(active_pw), &
            " selected_G2_max=", selected_g2_max, " selected_Emax=", selected_max_energy
        end if
      end if
    else if (use_shell_nraw_select) then
      nraw_old = 0
      do ishell = 1, shell_limit
        shell_id = shell_values(ishell)
        shell_count = count(shell_key(1:npw) == shell_id)
        nraw_new = nraw_old + shell_count
        if (nraw_new <= target_nraw) then
          do ipw = 1, npw
            if (shell_key(ipw) == shell_id) active_pw(ipw) = .true.
          end do
          selected_nshell = ishell
          selected_g2_max = shell_id
          selected_max_energy = shell_ekin(ishell)
          nraw_old = nraw_new
        else
          if (nraw_old <= 0 .or. abs(nraw_new - target_nraw) < abs(nraw_old - target_nraw)) then
            do ipw = 1, npw
              if (shell_key(ipw) == shell_id) active_pw(ipw) = .true.
            end do
            selected_nshell = ishell
            selected_g2_max = shell_id
            selected_max_energy = shell_ekin(ishell)
          end if
          exit
        end if
      end do
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,1pe13.5)') &
          "[DG-BPW-SELECT] mode=", trim(bpw_select_mode), " target_nraw=", target_nraw, &
          " selected_shells=", selected_nshell, " selected_nraw=", count(active_pw), &
          " selected_G2_max=", selected_g2_max, " selected_Emax=", selected_max_energy
      end if
    else if (use_comm_min_select) then
      active_pw(:) = .false.
      trial_active_pw(:) = .false.
      comm_old = 0.0d0
      fwp_old = 0.0d0
      fwp_axis_old(:) = 0.0d0
      comm_pair_old(:) = 0.0d0
      best_comm = huge(1.0d0)
      best_fwp = 0.0d0
      np_old = 0
      nraw_old = 0
      best_np = 0
      best_nraw = 0
      best_active_pw(:) = .false.
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,1pe13.5,a,1pe13.5,a,i0,a,i0)') &
          "[DG-MIXED-Z-COMM-MIN] start shells=", nshell, " rel_tol=", comm_rel_tol, &
          " abs_tol=", comm_abs_tol, " target_nraw=", target_nraw, " max_shells=", shell_limit
      end if
      do ishell = 1, shell_limit
        shell_id = shell_values(ishell)
        trial_active_pw(:) = active_pw(:)
        do ipw = 1, npw
          if (shell_key(ipw) == shell_id) trial_active_pw(ipw) = .true.
        end do
        nraw_new = count(trial_active_pw)
        if (nraw_new > target_nraw .and. nraw_old > 0 .and. fwp_old > comm_min_fwp_floor) cycle
        call evaluate_bpw_active_set(trial_active_pw, nraw_new, np_new, comm_new, fwp_new, min_s_trial, max_s_trial, &
          fwp_axis_new, comm_pair_new)
        if (np_new > 0 .and. fwp_new > comm_min_fwp_floor .and. comm_new < best_comm) then
          best_active_pw(:) = trial_active_pw(:)
          best_comm = comm_new
          best_fwp = fwp_new
          best_np = np_new
          best_nraw = nraw_new
        end if
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,7(a,1pe13.5),a)') &
            "[DG-MIXED-Z-COMM-MIN] shell=", ishell, " G2=", shell_id, &
            " nraw_if_accept=", nraw_new, " np_if_accept=", np_new, &
            " C_old=", comm_old, " C_new=", comm_new, " dC=", comm_new - comm_old, &
            " fWP_old=", fwp_old, " fWP_new=", fwp_new, &
            " min_Sperp=", min_s_trial, " max_Sperp=", max_s_trial, &
            " decision=", merge("accept", "reject", comm_new <= comm_old * (1.0d0 + comm_rel_tol) + comm_abs_tol)
        end if
        if (np_new > 0 .and. comm_new <= comm_old * (1.0d0 + comm_rel_tol) + comm_abs_tol) then
          active_pw(:) = trial_active_pw(:)
          nraw_old = nraw_new
          np_old = np_new
          comm_old = comm_new
          fwp_old = fwp_new
          fwp_axis_old(:) = fwp_axis_new(:)
          comm_pair_old(:) = comm_pair_new(:)
        end if
        if (count(active_pw) >= target_nraw) exit
      end do
      if (best_nraw > 0 .and. (fwp_old <= comm_min_fwp_floor .or. best_comm < comm_old)) then
        active_pw(:) = best_active_pw(:)
        nraw_old = best_nraw
        np_old = best_np
        comm_old = best_comm
        fwp_old = best_fwp
        if (comm_is_root(dg_frag%id)) &
          write(*,'(1x,a,i0,a,i0,a,1pe13.5,a,1pe13.5)') &
          "[DG-MIXED-Z-COMM-MIN] use best nonzero-coupling set: nraw=", nraw_old, &
          " np=", np_old, " C_comm=", comm_old, " fWP_avg_per_occ=", fwp_old
      end if
      if (count(active_pw) <= 0) then
        active_pw(:) = .true.
        if (comm_is_root(dg_frag%id)) &
          write(*,'(1x,a)') "[DG-MIXED-Z-COMM-MIN] no shell accepted; falling back to all active BPW."
      else if (fwp_old <= comm_min_fwp_floor) then
        active_pw(:) = .true.
        if (comm_is_root(dg_frag%id)) &
          write(*,'(1x,a,1pe13.5,a)') "[DG-MIXED-Z-COMM-MIN] selected set has negligible fWP=", &
          fwp_old, "; falling back to complete shell BPW."
      else if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,i0,a,1pe13.5,a,1pe13.5)') &
          "[DG-MIXED-Z-COMM-MIN] selected nraw=", count(active_pw), " np=", np_old, &
          " C_comm=", comm_old, " fWP_avg_per_occ=", fwp_old
      end if
    end if

    if (count(active_pw) <= 0) then
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a,a,a)') "[DG-BPW-SELECT] mode=", trim(bpw_select_mode), " selected no BPW raw states."
    end if

    if (comm_diag_enabled .or. fsum_diag_enabled) then
      trial_active_pw(:) = .false.
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') "BPW shell diag:"
        if (bpw_ecut_is_set) then
          write(*,'(3x,a,a,1x,a,1pe13.5)') "mode=", trim(bpw_select_mode), "Ecut=", bpw_ecut
        else
          write(*,'(3x,a,a,1x,a)') "mode=", trim(bpw_select_mode), "Ecut=all_raw"
        end if
      end if
      do ishell = 1, shell_limit
        shell_id = shell_values(ishell)
        do ipw = 1, npw
          if (shell_key(ipw) == shell_id) trial_active_pw(ipw) = .true.
        end do
        nraw_new = count(trial_active_pw)
        call evaluate_bpw_active_set(trial_active_pw, nraw_new, np_new, comm_new, fwp_new, min_s_trial, max_s_trial, &
          fwp_axis_new, comm_pair_new)
        if (use_shell_ecut_select) then
          accepted_by_ecut_int = merge(1, 0, (.not. bpw_ecut_is_set) .or. shell_ekin(ishell) <= bpw_ecut)
        else
          accepted_by_ecut_int = 1
          do ipw = 1, npw
            if (shell_key(ipw) == shell_id .and. .not. active_pw(ipw)) accepted_by_ecut_int = 0
          end do
        end if
        if (comm_is_root(dg_frag%id)) then
          write(*,'(3x,a,i0,1x,a,i0,1x,a,1pe13.5,1x,a,i0,1x,a,i0,1x,a,1pe13.5)') &
            "shell_id=", ishell, "G2=", shell_id, "E_shell=", shell_ekin(ishell), &
            "nraw=", nraw_new, "np=", np_new, "fWP_avg_occ=", fwp_new
          write(*,'(5x,a,3(1x,1pe13.5),1x,a,3(1x,1pe13.5),1x,a,i0)') &
            "fWP_xyz_occ=", fwp_axis_new(1:3) / max(1.0d0, dble(nocc_ref)), &
            "C_comm_xyz=", comm_pair_new(1:3), "accepted_by_ecut=", accepted_by_ecut_int
          write(*,'(5x,a,3(1x,1pe13.5),1x,a,1pe13.5,1x,a,1pe13.5)') &
            "fsum_xyz=", fwp_axis_new(1:3), "min_Sperp=", min_s_trial, "max_Sperp=", max_s_trial
        end if
      end do
    end if

    do ispin = 1, nspin_file
      b_local(:, :) = zzero
      do i_local = 1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1
        ifrag = dg_frag%ifrag_start + i_local - 1
        if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
        nbf = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%global_wannier_coef, 1))
        do ibasis = 1, nbf
          global_row = dg_frag%index_basis(ibasis, ifrag, ispin)
          if (global_row < 1 .or. global_row > nmat) cycle
          local_row = 0
          if (allocated(dg_frag%coef_global_to_local)) local_row = dg_frag%coef_global_to_local(global_row, ispin)
          if (local_row < 1 .or. local_row > nmat) cycle
          do iw = 1, nwann
            coeff = conjg(dg_frag%global_wannier_coef(ibasis, iw, ispin, i_local))
            b_local(iw,1:npw) = b_local(iw,1:npw) + coeff * s_fp(local_row,1:npw,ispin)
          end do
        end do
      end do
      call comm_summation(b_local, b_global, nwann*npw, dg_frag%icomm)
      do ipw = 1, npw
        if (.not. active_pw(ipw)) b_global(:,ipw) = zzero
      end do

      s_perp(:, :) = zzero
      do ipw = 1, npw
        if (active_pw(ipw)) s_perp(ipw,ipw) = (1.0d0, 0.0d0)
      end do
      s_perp(:, :) = s_perp(:, :) - matmul(conjg(transpose(b_global)), b_global)
      s_vec = s_perp
      call eigen_zheev(s_vec, s_eval, s_work)
      nkeep = count(s_eval > sperp_tol)
      if (nkeep <= 0) cycle

      if (allocated(dg_frag%mixed_wannier_bpw_eval)) deallocate(dg_frag%mixed_wannier_bpw_eval)
      if (allocated(dg_frag%mixed_wannier_bpw_z)) deallocate(dg_frag%mixed_wannier_bpw_z)
      if (allocated(dg_frag%mixed_wannier_bpw_pcoef)) deallocate(dg_frag%mixed_wannier_bpw_pcoef)
      if (allocated(dg_frag%mixed_wannier_bpw_p_transform)) deallocate(dg_frag%mixed_wannier_bpw_p_transform)
      if (allocated(dg_frag%mixed_wannier_bpw_p_metric)) deallocate(dg_frag%mixed_wannier_bpw_p_metric)
      allocate(dg_frag%mixed_wannier_bpw_eval(neig+nkeep,nspin_file))
      allocate(dg_frag%mixed_wannier_bpw_z(3,neig+nkeep,neig+nkeep,nspin_file))
      allocate(dg_frag%mixed_wannier_bpw_pcoef(nkeep,max(1,dg_frag%nstate_tot),nspin_file))
      allocate(dg_frag%mixed_wannier_bpw_p_transform(npw,nkeep,nspin_file))
      allocate(dg_frag%mixed_wannier_bpw_p_metric(nkeep,nkeep,nspin_file))
      dg_frag%mixed_wannier_bpw_eval(:, :) = 0.0d0
      dg_frag%mixed_wannier_bpw_z(:, :, :, :) = zzero
      dg_frag%mixed_wannier_bpw_pcoef(:, :, :) = zzero
      dg_frag%mixed_wannier_bpw_p_transform(:, :, :) = zzero
      dg_frag%mixed_wannier_bpw_p_metric(:, :, :) = zzero
      dg_frag%mixed_wannier_bpw_nw = neig
      dg_frag%mixed_wannier_bpw_np = nkeep
      dg_frag%mixed_wannier_bpw_nmix = neig + nkeep
      dg_frag%mixed_wannier_bpw_praw_dim = npw
      dg_frag%has_mixed_wannier_bpw_p_basis = .false.
      dg_frag%mixed_wannier_bpw_sraw_herm_diff = 0.0d0
      dg_frag%mixed_wannier_bpw_sperp_herm_diff = 0.0d0
      dg_frag%mixed_wannier_bpw_qmat_metric_herm_diff = 0.0d0
      dg_frag%mixed_wannier_bpw_qmat_metric_min_eval = huge(1.0d0)
      dg_frag%mixed_wannier_bpw_qmat_metric_max_eval = 0.0d0
      dg_frag%mixed_wannier_bpw_qmat_metric_cond = 0.0d0
      dg_frag%mixed_wannier_bpw_qmat_metric_diff_from_i = 0.0d0
      dg_frag%mixed_wannier_bpw_qleft_metric_diff_from_i = 0.0d0
      dg_frag%mixed_wannier_bpw_final_metric_herm_diff = 0.0d0
      dg_frag%mixed_wannier_bpw_final_metric_min_eval = huge(1.0d0)
      dg_frag%mixed_wannier_bpw_final_metric_max_eval = 0.0d0
      dg_frag%mixed_wannier_bpw_final_metric_cond = 0.0d0
      dg_frag%mixed_wannier_bpw_final_metric_diff_from_i = 0.0d0
      dg_frag%mixed_wannier_bpw_transform_metric_herm_diff = 0.0d0
      dg_frag%mixed_wannier_bpw_transform_metric_min_eval = huge(1.0d0)
      dg_frag%mixed_wannier_bpw_transform_metric_max_eval = 0.0d0
      dg_frag%mixed_wannier_bpw_transform_metric_cond = 0.0d0
      dg_frag%mixed_wannier_bpw_transform_metric_diff_from_i = 0.0d0
      dg_frag%mixed_wannier_bpw_transform_metric_diff_saved = 0.0d0
      dg_frag%mixed_wannier_bpw_h_input_herm_diff = 0.0d0
      dg_frag%mixed_wannier_bpw_h_evec_unitarity_diff = 0.0d0
      dg_frag%mixed_wannier_bpw_h_input_evec_diff = 0.0d0
      dg_frag%mixed_wannier_bpw_final_uses_h_evec = use_h_eigenvectors_for_final
      dg_frag%mixed_wannier_bpw_qmat_col_norm_min = huge(1.0d0)
      dg_frag%mixed_wannier_bpw_qmat_col_norm_max = 0.0d0
      dg_frag%mixed_wannier_bpw_qmat_row_norm_min = huge(1.0d0)
      dg_frag%mixed_wannier_bpw_qmat_row_norm_max = 0.0d0
      exit
    end do

    if (.not. allocated(dg_frag%mixed_wannier_bpw_eval)) then
      deallocate(s_fp, r_fp, b_local, b_global, r_local, r_global, b_eig)
      deallocate(c_local, c_global, c_w_local, c_w_global, c_eig)
      deallocate(z_w, z_eig, o_eig)
      deallocate(s_perp, s_vec, s_work, s_eval, pw_ekin, shell_ekin, qmat, h_p, h_vec, h_eval)
      deallocate(q_metric, q_left_metric, final_metric, transform_metric, t_final, t_saved, metric_vec, diag_eval)
      deallocate(total_metric, h_total)
      deallocate(r_tilde, r_orth, r_final, z_wp_direct)
      return
    end if

    min_s_local = huge(1.0d0)
    max_s_local = -huge(1.0d0)
    zwp_direct_norm(:) = 0.0d0
    do ispin = 1, nspin_file
      nocc_ref = min(neig - 1, max(1, dg_frag%nstate_tot))
      if (allocated(dg_frag%nocc_spin)) then
        if (ispin <= size(dg_frag%nocc_spin)) nocc_ref = min(neig - 1, max(1, dg_frag%nocc_spin(ispin)))
      end if

      b_local(:, :) = zzero
      do i_local = 1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1
        ifrag = dg_frag%ifrag_start + i_local - 1
        if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
        nbf = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%global_wannier_coef, 1))
        do ibasis = 1, nbf
          global_row = dg_frag%index_basis(ibasis, ifrag, ispin)
          if (global_row < 1 .or. global_row > nmat) cycle
          local_row = 0
          if (allocated(dg_frag%coef_global_to_local)) local_row = dg_frag%coef_global_to_local(global_row, ispin)
          if (local_row < 1 .or. local_row > nmat) cycle
          do iw = 1, nwann
            coeff = conjg(dg_frag%global_wannier_coef(ibasis, iw, ispin, i_local))
            b_local(iw,1:npw) = b_local(iw,1:npw) + coeff * s_fp(local_row,1:npw,ispin)
          end do
        end do
      end do
      call comm_summation(b_local, b_global, nwann*npw, dg_frag%icomm)
      do ipw = 1, npw
        if (.not. active_pw(ipw)) b_global(:,ipw) = zzero
      end do
      b_eig(:, :) = zzero
      do eig = 1, neig
        do ipw = 1, npw
          do iw = 1, nwann
            b_eig(eig,ipw) = b_eig(eig,ipw) + conjg(dg_frag%global_wannier_flux_evec(iw,eig)) * b_global(iw,ipw)
          end do
        end do
      end do

      s_perp(:, :) = zzero
      do ipw = 1, npw
        if (active_pw(ipw)) s_perp(ipw,ipw) = (1.0d0, 0.0d0)
      end do
      s_perp(:, :) = s_perp(:, :) - matmul(conjg(transpose(b_global)), b_global)
      s_vec = s_perp
      call eigen_zheev(s_vec, s_eval, s_work)
      qmat(:, :) = zzero
      pslot = 0
      do alpha = 1, npw
        if (s_eval(alpha) <= sperp_tol) cycle
        pslot = pslot + 1
        if (pslot > dg_frag%mixed_wannier_bpw_np) exit
        qmat(1:npw,pslot) = s_work(1:npw,alpha) / sqrt(s_eval(alpha))
      end do
      min_s_local = min(min_s_local, minval(s_eval))
      max_s_local = max(max_s_local, maxval(s_eval))

      h_p(:, :) = zzero
      do alpha = 1, dg_frag%mixed_wannier_bpw_np
        do beta = 1, dg_frag%mixed_wannier_bpw_np
          do ipw = 1, npw
            h_p(alpha,beta) = h_p(alpha,beta) + conjg(qmat(ipw,alpha)) * pw_ekin(ipw) * qmat(ipw,beta)
          end do
        end do
      end do
      h_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) = &
        h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np)
      dg_frag%mixed_wannier_bpw_h_input_herm_diff = max(dg_frag%mixed_wannier_bpw_h_input_herm_diff, &
        maxval(abs(h_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) - &
        conjg(transpose(h_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))))))
      call eigen_zheev(h_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np), &
        h_eval(1:dg_frag%mixed_wannier_bpw_np), h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
      dg_frag%mixed_wannier_bpw_h_input_evec_diff = max(dg_frag%mixed_wannier_bpw_h_input_evec_diff, &
        sqrt(sum(abs(h_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) - &
        h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))**2)))
      q_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) = matmul( &
        conjg(transpose(h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))), &
        h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
      do alpha = 1, dg_frag%mixed_wannier_bpw_np
        q_metric(alpha,alpha) = q_metric(alpha,alpha) - (1.0d0,0.0d0)
      end do
      dg_frag%mixed_wannier_bpw_h_evec_unitarity_diff = max(dg_frag%mixed_wannier_bpw_h_evec_unitarity_diff, &
        sqrt(sum(abs(q_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))**2)))
      if (allocated(dg_frag%mixed_wannier_bpw_p_transform) .and. allocated(dg_frag%mixed_wannier_bpw_p_metric)) then
        dg_frag%mixed_wannier_bpw_p_transform(1:npw,1:dg_frag%mixed_wannier_bpw_np,ispin) = &
          matmul(qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np), &
                 h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
        dg_frag%mixed_wannier_bpw_p_metric(1:dg_frag%mixed_wannier_bpw_np, &
          1:dg_frag%mixed_wannier_bpw_np,ispin) = matmul( &
          conjg(transpose(dg_frag%mixed_wannier_bpw_p_transform(1:npw,1:dg_frag%mixed_wannier_bpw_np,ispin))), &
          matmul(s_perp(1:npw,1:npw), &
            dg_frag%mixed_wannier_bpw_p_transform(1:npw,1:dg_frag%mixed_wannier_bpw_np,ispin)))
        dg_frag%has_mixed_wannier_bpw_p_basis = .true.
      end if
      dg_frag%mixed_wannier_bpw_sraw_herm_diff = max(dg_frag%mixed_wannier_bpw_sraw_herm_diff, &
        maxval(abs(s_vec(1:npw,1:npw) - conjg(transpose(s_vec(1:npw,1:npw))))))
      dg_frag%mixed_wannier_bpw_sperp_herm_diff = max(dg_frag%mixed_wannier_bpw_sperp_herm_diff, &
        maxval(abs(s_perp(1:npw,1:npw) - conjg(transpose(s_perp(1:npw,1:npw))))))
      q_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) = matmul( &
        conjg(transpose(qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np))), &
        matmul(s_perp(1:npw,1:npw), qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np)))
      q_left_metric(1:npw,1:npw) = matmul(qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np), &
        matmul(s_perp(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np), &
        conjg(transpose(qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np)))))
      if (use_h_eigenvectors_for_final) then
        t_final(1:npw,1:dg_frag%mixed_wannier_bpw_np) = matmul( &
          qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np), &
          h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
      else
        t_final(1:npw,1:dg_frag%mixed_wannier_bpw_np) = matmul( &
          qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np), &
          h_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
      end if
      t_saved(1:npw,1:dg_frag%mixed_wannier_bpw_np) = matmul( &
        qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np), &
        h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
      final_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) = matmul( &
        conjg(transpose(t_final(1:npw,1:dg_frag%mixed_wannier_bpw_np))), &
        matmul(s_perp(1:npw,1:npw), t_final(1:npw,1:dg_frag%mixed_wannier_bpw_np)))
      transform_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) = matmul( &
        conjg(transpose(t_saved(1:npw,1:dg_frag%mixed_wannier_bpw_np))), &
        matmul(s_perp(1:npw,1:npw), t_saved(1:npw,1:dg_frag%mixed_wannier_bpw_np)))
      total_dim = neig + dg_frag%mixed_wannier_bpw_np
      total_metric(1:total_dim,1:total_dim) = zzero
      total_metric(1:neig,1:neig) = matmul( &
        conjg(transpose(dg_frag%global_wannier_flux_evec(1:nwann,1:neig))), &
        dg_frag%global_wannier_flux_evec(1:nwann,1:neig))
      do alpha = 1, neig
        total_metric(alpha,alpha) = total_metric(alpha,alpha) - (1.0d0,0.0d0)
      end do
      total_metric(1:neig,neig+1:total_dim) = matmul( &
        b_eig(1:neig,1:npw), t_final(1:npw,1:dg_frag%mixed_wannier_bpw_np))
      total_metric(neig+1:total_dim,1:neig) = conjg(transpose( &
        total_metric(1:neig,neig+1:total_dim)))
      total_metric(neig+1:total_dim,neig+1:total_dim) = &
        final_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np)
      do alpha = 1, dg_frag%mixed_wannier_bpw_np
        total_metric(neig+alpha,neig+alpha) = total_metric(neig+alpha,neig+alpha) - (1.0d0,0.0d0)
      end do
      metric_sww_max = max(metric_sww_max, maxval(abs(total_metric(1:neig,1:neig))))
      metric_sww_rms = max(metric_sww_rms, sqrt(sum(abs(total_metric(1:neig,1:neig))**2) / &
        max(1.0d0, dble(neig*neig))))
      metric_swz_max = max(metric_swz_max, maxval(abs(total_metric(1:neig,neig+1:total_dim))))
      metric_swz_rms = max(metric_swz_rms, sqrt(sum(abs(total_metric(1:neig,neig+1:total_dim))**2) / &
        max(1.0d0, dble(neig*dg_frag%mixed_wannier_bpw_np))))
      metric_szz_max = max(metric_szz_max, maxval(abs(total_metric(neig+1:total_dim,neig+1:total_dim))))
      metric_szz_rms = max(metric_szz_rms, sqrt(sum(abs(total_metric(neig+1:total_dim,neig+1:total_dim))**2) / &
        max(1.0d0, dble(dg_frag%mixed_wannier_bpw_np*dg_frag%mixed_wannier_bpw_np))))
      metric_stotal_max = max(metric_stotal_max, maxval(abs(total_metric(1:total_dim,1:total_dim))))
      metric_stotal_rms = max(metric_stotal_rms, sqrt(sum(abs(total_metric(1:total_dim,1:total_dim))**2) / &
        max(1.0d0, dble(total_dim*total_dim))))
      h_total(1:total_dim,1:total_dim) = zzero
      do alpha = 1, neig
        h_total(alpha,alpha) = cmplx(dg_frag%global_wannier_flux_eval(alpha,ispin), 0.0d0, kind=8)
      end do
      do alpha = 1, dg_frag%mixed_wannier_bpw_np
        h_total(neig+alpha,neig+alpha) = cmplx(h_eval(alpha), 0.0d0, kind=8)
      end do
      h_total_herm_max = max(h_total_herm_max, maxval(abs(h_total(1:total_dim,1:total_dim) - &
        conjg(transpose(h_total(1:total_dim,1:total_dim))))))
      h_total_herm_rms = max(h_total_herm_rms, sqrt(sum(abs(h_total(1:total_dim,1:total_dim) - &
        conjg(transpose(h_total(1:total_dim,1:total_dim))))**2) / max(1.0d0, dble(total_dim*total_dim))))
      call eigen_zheev(q_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np), &
        diag_eval(1:dg_frag%mixed_wannier_bpw_np), metric_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
      dg_frag%mixed_wannier_bpw_qmat_metric_herm_diff = max(dg_frag%mixed_wannier_bpw_qmat_metric_herm_diff, &
        maxval(abs(q_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) - &
        conjg(transpose(q_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))))))
      dg_frag%mixed_wannier_bpw_qmat_metric_min_eval = min(dg_frag%mixed_wannier_bpw_qmat_metric_min_eval, &
        minval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)))
      dg_frag%mixed_wannier_bpw_qmat_metric_max_eval = max(dg_frag%mixed_wannier_bpw_qmat_metric_max_eval, &
        maxval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)))
      dg_frag%mixed_wannier_bpw_qmat_metric_cond = max(dg_frag%mixed_wannier_bpw_qmat_metric_cond, &
        maxval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)) / &
        max(minval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)), 1.0d-300))
      do alpha = 1, dg_frag%mixed_wannier_bpw_np
        q_metric(alpha,alpha) = q_metric(alpha,alpha) - (1.0d0,0.0d0)
      end do
      dg_frag%mixed_wannier_bpw_qmat_metric_diff_from_i = max(dg_frag%mixed_wannier_bpw_qmat_metric_diff_from_i, &
        sqrt(sum(abs(q_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))**2)))
      if (npw == dg_frag%mixed_wannier_bpw_np) then
        do alpha = 1, npw
          q_left_metric(alpha,alpha) = q_left_metric(alpha,alpha) - (1.0d0,0.0d0)
        end do
        dg_frag%mixed_wannier_bpw_qleft_metric_diff_from_i = max( &
          dg_frag%mixed_wannier_bpw_qleft_metric_diff_from_i, sqrt(sum(abs(q_left_metric(1:npw,1:npw))**2)))
      end if
      call eigen_zheev(final_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np), &
        diag_eval(1:dg_frag%mixed_wannier_bpw_np), metric_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
      dg_frag%mixed_wannier_bpw_final_metric_herm_diff = max(dg_frag%mixed_wannier_bpw_final_metric_herm_diff, &
        maxval(abs(final_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) - &
        conjg(transpose(final_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))))))
      dg_frag%mixed_wannier_bpw_final_metric_min_eval = min(dg_frag%mixed_wannier_bpw_final_metric_min_eval, &
        minval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)))
      dg_frag%mixed_wannier_bpw_final_metric_max_eval = max(dg_frag%mixed_wannier_bpw_final_metric_max_eval, &
        maxval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)))
      dg_frag%mixed_wannier_bpw_final_metric_cond = max(dg_frag%mixed_wannier_bpw_final_metric_cond, &
        maxval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)) / &
        max(minval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)), 1.0d-300))
      do alpha = 1, dg_frag%mixed_wannier_bpw_np
        final_metric(alpha,alpha) = final_metric(alpha,alpha) - (1.0d0,0.0d0)
      end do
      dg_frag%mixed_wannier_bpw_final_metric_diff_from_i = max(dg_frag%mixed_wannier_bpw_final_metric_diff_from_i, &
        sqrt(sum(abs(final_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))**2)))
      call eigen_zheev(transform_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np), &
        diag_eval(1:dg_frag%mixed_wannier_bpw_np), metric_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
      dg_frag%mixed_wannier_bpw_transform_metric_herm_diff = max( &
        dg_frag%mixed_wannier_bpw_transform_metric_herm_diff, &
        maxval(abs(transform_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) - &
        conjg(transpose(transform_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))))))
      dg_frag%mixed_wannier_bpw_transform_metric_min_eval = min(dg_frag%mixed_wannier_bpw_transform_metric_min_eval, &
        minval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)))
      dg_frag%mixed_wannier_bpw_transform_metric_max_eval = max(dg_frag%mixed_wannier_bpw_transform_metric_max_eval, &
        maxval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)))
      dg_frag%mixed_wannier_bpw_transform_metric_cond = max(dg_frag%mixed_wannier_bpw_transform_metric_cond, &
        maxval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)) / &
        max(minval(diag_eval(1:dg_frag%mixed_wannier_bpw_np)), 1.0d-300))
      if (allocated(dg_frag%mixed_wannier_bpw_p_metric)) then
        dg_frag%mixed_wannier_bpw_transform_metric_diff_saved = max( &
          dg_frag%mixed_wannier_bpw_transform_metric_diff_saved, sqrt(sum(abs( &
          transform_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np) - &
          dg_frag%mixed_wannier_bpw_p_metric(1:dg_frag%mixed_wannier_bpw_np, &
          1:dg_frag%mixed_wannier_bpw_np,ispin))**2)))
      end if
      do alpha = 1, dg_frag%mixed_wannier_bpw_np
        transform_metric(alpha,alpha) = transform_metric(alpha,alpha) - (1.0d0,0.0d0)
      end do
      dg_frag%mixed_wannier_bpw_transform_metric_diff_from_i = max( &
        dg_frag%mixed_wannier_bpw_transform_metric_diff_from_i, &
        sqrt(sum(abs(transform_metric(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))**2)))
      dg_frag%mixed_wannier_bpw_qmat_col_norm_min = min(dg_frag%mixed_wannier_bpw_qmat_col_norm_min, &
        minval(sum(abs(qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np))**2, dim=1)))
      dg_frag%mixed_wannier_bpw_qmat_col_norm_max = max(dg_frag%mixed_wannier_bpw_qmat_col_norm_max, &
        maxval(sum(abs(qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np))**2, dim=1)))
      dg_frag%mixed_wannier_bpw_qmat_row_norm_min = min(dg_frag%mixed_wannier_bpw_qmat_row_norm_min, &
        minval(sum(abs(qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np))**2, dim=2)))
      dg_frag%mixed_wannier_bpw_qmat_row_norm_max = max(dg_frag%mixed_wannier_bpw_qmat_row_norm_max, &
        maxval(sum(abs(qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np))**2, dim=2)))
      pw_energy_ref = dg_frag%global_wannier_flux_eval(nocc_ref+1,ispin)
      h_eval(1:dg_frag%mixed_wannier_bpw_np) = h_eval(1:dg_frag%mixed_wannier_bpw_np) + pw_energy_ref

      dg_frag%mixed_wannier_bpw_eval(1:neig,ispin) = dg_frag%global_wannier_flux_eval(1:neig,ispin)
      dg_frag%mixed_wannier_bpw_eval(neig+1:neig+dg_frag%mixed_wannier_bpw_np,ispin) = &
        h_eval(1:dg_frag%mixed_wannier_bpw_np)

      do idir = 1, 3
        z_w(1:nwann,1:nwann) = dg_frag%global_wannier_position(idir,1:nwann,1:nwann)
        z_eig(1:neig,1:neig) = matmul(conjg(transpose(dg_frag%global_wannier_flux_evec(1:nwann,1:neig))), &
          matmul(z_w, dg_frag%global_wannier_flux_evec(1:nwann,1:neig)))
        dg_frag%mixed_wannier_bpw_z(idir,1:neig,1:neig,ispin) = z_eig(1:neig,1:neig)
        if (yn_dg_mixed_z_include_ww /= 'y') &
          dg_frag%mixed_wannier_bpw_z(idir,1:neig,1:neig,ispin) = zzero

        r_local(:, :) = zzero
        c_local(:, :) = zzero
        c_w_local(:, :) = zzero
        do i_local = 1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1
          ifrag = dg_frag%ifrag_start + i_local - 1
          if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
          frag_center(1:3) = 0.0d0
          if (allocated(dg_frag%buffer_wannier_frag_center) .and. i_local <= size(dg_frag%buffer_wannier_frag_center, 2)) then
            frag_center(1:3) = dg_frag%buffer_wannier_frag_center(1:3,i_local)
          else
            do jdir = 1, 3
              frag_center(jdir) = dg_frag%hgs(jdir) * (dble(dg_frag%ixyz_frag(jdir,ifrag) - 1) + &
                0.5d0 * dble(dg_frag%nxyz_domain(jdir,ifrag) - 1))
            end do
          end if
          nbf = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%global_wannier_coef, 1))
          do ibasis = 1, nbf
            global_row = dg_frag%index_basis(ibasis, ifrag, ispin)
            if (global_row < 1 .or. global_row > nmat) cycle
            local_row = 0
            if (allocated(dg_frag%coef_global_to_local)) local_row = dg_frag%coef_global_to_local(global_row, ispin)
            if (local_row < 1 .or. local_row > nmat) cycle
            do iw = 1, nwann
              coeff = conjg(dg_frag%global_wannier_coef(ibasis, iw, ispin, i_local))
              r_local(iw,1:npw) = r_local(iw,1:npw) + coeff * r_fp(idir,local_row,1:npw,ispin)
              c_local(iw,1:npw) = c_local(iw,1:npw) + coeff * frag_center(idir) * s_fp(local_row,1:npw,ispin)
              do jw = 1, nwann
                c_w_local(iw,jw) = c_w_local(iw,jw) + coeff * frag_center(idir) * &
                  dg_frag%global_wannier_coef(ibasis, jw, ispin, i_local)
              end do
            end do
          end do
        end do
        call comm_summation(r_local, r_global, nwann*npw, dg_frag%icomm)
        call comm_summation(c_local, c_global, nwann*npw, dg_frag%icomm)
        call comm_summation(c_w_local, c_w_global, nwann*nwann, dg_frag%icomm)

        c_eig(1:neig,1:neig) = matmul(conjg(transpose(dg_frag%global_wannier_flux_evec(1:nwann,1:neig))), &
          matmul(c_w_global(1:nwann,1:nwann), dg_frag%global_wannier_flux_evec(1:nwann,1:neig)))
        if (use_fragment_center_direct) then
          o_eig(1:neig,1:neig) = z_eig(1:neig,1:neig) - c_eig(1:neig,1:neig)
        else
          o_eig(1:neig,1:neig) = z_eig(1:neig,1:neig)
        end if
        r_tilde(:, :) = zzero
        do eig = 1, neig
          do ipw = 1, npw
            amp = zzero
            do iw = 1, nwann
              amp = amp + conjg(dg_frag%global_wannier_flux_evec(iw,eig)) * r_global(iw,ipw)
            end do
            if (use_fragment_center_direct) then
              do iw = 1, nwann
                amp = amp - conjg(dg_frag%global_wannier_flux_evec(iw,eig)) * c_global(iw,ipw)
              end do
            end if
            grad_corr = zzero
            do jw = 1, neig
              grad_corr = grad_corr + o_eig(eig,jw) * b_eig(jw,ipw)
            end do
            r_tilde(eig,ipw) = amp - grad_corr
          end do
        end do
        r_orth(1:neig,1:dg_frag%mixed_wannier_bpw_np) = &
          matmul(r_tilde(1:neig,1:npw), qmat(1:npw,1:dg_frag%mixed_wannier_bpw_np))
        if (use_h_eigenvectors_for_final) then
          r_final(1:neig,1:dg_frag%mixed_wannier_bpw_np) = &
            matmul(r_orth(1:neig,1:dg_frag%mixed_wannier_bpw_np), &
                   h_p(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
        else
          r_final(1:neig,1:dg_frag%mixed_wannier_bpw_np) = &
            matmul(r_orth(1:neig,1:dg_frag%mixed_wannier_bpw_np), &
                   h_vec(1:dg_frag%mixed_wannier_bpw_np,1:dg_frag%mixed_wannier_bpw_np))
        end if
        z_wp_direct(:, :) = zzero
        do eig = 1, neig
          do alpha = 1, dg_frag%mixed_wannier_bpw_np
            z_wp_direct(eig,alpha) = r_final(eig,alpha)
            dg_frag%mixed_wannier_bpw_z(idir,eig,neig+alpha,ispin) = z_wp_direct(eig,alpha)
            dg_frag%mixed_wannier_bpw_z(idir,neig+alpha,eig,ispin) = &
              conjg(dg_frag%mixed_wannier_bpw_z(idir,eig,neig+alpha,ispin))
          end do
        end do
        zwp_direct_norm(idir) = zwp_direct_norm(idir) + sqrt(sum(abs(z_wp_direct(1:neig,1:dg_frag%mixed_wannier_bpw_np))**2))
      end do
    end do

    herm_max = 0.0d0
    zwp_norm(:) = 0.0d0
    zww_norm(:) = 0.0d0
    zcov_ww(:, :) = 0.0d0
    zcov_wp(:, :) = 0.0d0
    zcov_tot(:, :) = 0.0d0
    zcomm_ww(:, :) = 0.0d0
    zcomm_wp(:, :) = 0.0d0
    zcomm_tot(:, :) = 0.0d0
    metric_sww_max = 0.0d0
    metric_sww_rms = 0.0d0
    metric_swz_max = 0.0d0
    metric_swz_rms = 0.0d0
    metric_szz_max = 0.0d0
    metric_szz_rms = 0.0d0
    metric_stotal_max = 0.0d0
    metric_stotal_rms = 0.0d0
    h_total_herm_max = 0.0d0
    h_total_herm_rms = 0.0d0
    do ispin = 1, nspin_file
      nocc_ref = min(dg_frag%mixed_wannier_bpw_nw, max(0, dg_frag%nstate_tot))
      if (allocated(dg_frag%nocc_spin)) then
        if (ispin <= size(dg_frag%nocc_spin)) nocc_ref = min(dg_frag%mixed_wannier_bpw_nw, max(0, dg_frag%nocc_spin(ispin)))
      end if
      do idir = 1, 3
        zww_norm(idir) = zww_norm(idir) + sqrt(sum(abs( &
          dg_frag%mixed_wannier_bpw_z(idir,1:dg_frag%mixed_wannier_bpw_nw, &
          1:dg_frag%mixed_wannier_bpw_nw,ispin))**2))
        zwp_norm(idir) = zwp_norm(idir) + sqrt(sum(abs( &
          dg_frag%mixed_wannier_bpw_z(idir,1:dg_frag%mixed_wannier_bpw_nw, &
          dg_frag%mixed_wannier_bpw_nw+1:dg_frag%mixed_wannier_bpw_nmix,ispin))**2))
        do jmix = 1, dg_frag%mixed_wannier_bpw_nmix
          do imix = 1, dg_frag%mixed_wannier_bpw_nmix
            herm_max = max(herm_max, abs(dg_frag%mixed_wannier_bpw_z(idir,imix,jmix,ispin) - &
              conjg(dg_frag%mixed_wannier_bpw_z(idir,jmix,imix,ispin))))
          end do
        end do
      end do
      do idir = 1, 3
        do jdir = 1, 3
          do eig = 1, nocc_ref
            do virt = nocc_ref + 1, dg_frag%mixed_wannier_bpw_nw
              zcov_ww(idir,jdir) = zcov_ww(idir,jdir) + real( &
                conjg(dg_frag%mixed_wannier_bpw_z(idir,eig,virt,ispin)) * &
                      dg_frag%mixed_wannier_bpw_z(jdir,eig,virt,ispin), 8)
              zcomm_ww(idir,jdir) = zcomm_ww(idir,jdir) + aimag( &
                dg_frag%mixed_wannier_bpw_z(idir,eig,virt,ispin) * &
                dg_frag%mixed_wannier_bpw_z(jdir,virt,eig,ispin) - &
                dg_frag%mixed_wannier_bpw_z(jdir,eig,virt,ispin) * &
                dg_frag%mixed_wannier_bpw_z(idir,virt,eig,ispin))
            end do
            do virt = dg_frag%mixed_wannier_bpw_nw + 1, dg_frag%mixed_wannier_bpw_nmix
              zcov_wp(idir,jdir) = zcov_wp(idir,jdir) + real( &
                conjg(dg_frag%mixed_wannier_bpw_z(idir,eig,virt,ispin)) * &
                      dg_frag%mixed_wannier_bpw_z(jdir,eig,virt,ispin), 8)
              zcomm_wp(idir,jdir) = zcomm_wp(idir,jdir) + aimag( &
                dg_frag%mixed_wannier_bpw_z(idir,eig,virt,ispin) * &
                dg_frag%mixed_wannier_bpw_z(jdir,virt,eig,ispin) - &
                dg_frag%mixed_wannier_bpw_z(jdir,eig,virt,ispin) * &
                dg_frag%mixed_wannier_bpw_z(idir,virt,eig,ispin))
            end do
          end do
        end do
      end do
    end do
    zcov_tot(:, :) = zcov_ww(:, :) + zcov_wp(:, :)
    zcomm_tot(:, :) = zcomm_ww(:, :) + zcomm_wp(:, :)
    fsum_proxy_wp_avg = (zcov_wp(1,1) + zcov_wp(2,2) + zcov_wp(3,3)) / &
      max(1.0d0, 3.0d0 * dble(max(1, nocc_ref)))
    fsum_proxy_total_avg = (zcov_tot(1,1) + zcov_tot(2,2) + zcov_tot(3,3)) / &
      max(1.0d0, 3.0d0 * dble(max(1, nocc_ref)))

    dg_frag%mixed_wannier_bpw_sperp_min = min_s_local
    dg_frag%mixed_wannier_bpw_sperp_max = max_s_local
    dg_frag%has_mixed_wannier_bpw_position = .true.
    nraw_selected = count(active_pw)
    call summarize_sperp_spectrum()
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0,a,1pe13.5,a,1pe13.5)') &
        "[DG-MIXED-Z] built dense Wannier+BPW-perp position: nw=", neig, " np=", &
        dg_frag%mixed_wannier_bpw_np, " min_Sperp=", min_s_local, " max_Sperp=", max_s_local
      write(*,'(1x,a,1pe13.5)') "[DG-MIXED-Z] hermiticity max|Z-Z^H|=", herm_max
      write(*,'(1x,a,a)') "[DG-MIXED-Z] include WW position block=", yn_dg_mixed_z_include_ww
      write(*,'(1x,a,3(1x,1pe13.5))') "[DG-MIXED-Z] ||Z_WW||_F xyz=", zww_norm(1:3)
      write(*,'(1x,a,3(1x,1pe13.5))') "[DG-MIXED-Z] ||Z_WP||_F xyz=", zwp_norm(1:3)
      write(*,'(1x,a)') "[DG-MIXED-Z] Z_WP construction = direct <W|r|BPW_perp>"
      if (use_fragment_center_direct) then
        write(*,'(1x,a)') "[DG-MIXED-Z] direct Z_WP origin = fragment center"
      else
        write(*,'(1x,a)') "[DG-MIXED-Z] direct Z_WP origin = global coordinate"
      end if
      write(*,'(1x,a,3(1x,1pe13.5))') "[DG-MIXED-Z] ||Z_WP_direct||_F xyz=", zwp_direct_norm(1:3)
      write(*,'(1x,a,1pe13.5)') "[DG-MIXED-Z] Lowdin S_perp cutoff=", sperp_tol
      write(*,'(1x,a,2(1x,1pe13.5),2(a,1pe13.5))') &
        "[DG-MIXED-Z] fsum_proxy avg WP,total=", fsum_proxy_wp_avg, fsum_proxy_total_avg, &
        " C_comm_sum_WP=", sqrt(sum(zcomm_wp(:, :)**2)), &
        " C_comm_sum_total=", sqrt(sum(zcomm_tot(:, :)**2))
      write(*,'(1x,a,3(a,i0),8(a,1pe13.5))') &
        "[DG-MIXED-Z-METRIC]", &
        " dim_W=", dg_frag%mixed_wannier_bpw_nw, &
        " dim_Z=", dg_frag%mixed_wannier_bpw_np, &
        " dim_total=", dg_frag%mixed_wannier_bpw_nmix, &
        " S_WW_minus_I_max=", metric_sww_max, &
        " S_WW_minus_I_rms=", metric_sww_rms, &
        " S_WZ_max=", metric_swz_max, &
        " S_WZ_rms=", metric_swz_rms, &
        " S_ZZ_minus_I_max=", metric_szz_max, &
        " S_ZZ_minus_I_rms=", metric_szz_rms, &
        " S_total_minus_I_max=", metric_stotal_max, &
        " S_total_minus_I_rms=", metric_stotal_rms
      write(*,'(1x,a,3(a,i0),2(a,1pe13.5))') &
        "[DG-MIXED-Z-HMETRIC]", &
        " dim_W=", dg_frag%mixed_wannier_bpw_nw, &
        " dim_Z=", dg_frag%mixed_wannier_bpw_np, &
        " dim_total=", dg_frag%mixed_wannier_bpw_nmix, &
        " H_total_herm_max=", h_total_herm_max, &
        " H_total_herm_rms=", h_total_herm_rms
      do idir = 1, 3
        write(*,'(1x,a,i0,a,3(1x,1pe13.5))') "[DG-MIXED-Z] occ-virt metric WW row ", idir, " =", zcov_ww(idir,1:3)
      end do
      do idir = 1, 3
        write(*,'(1x,a,i0,a,3(1x,1pe13.5))') "[DG-MIXED-Z] occ-virt metric WP row ", idir, " =", zcov_wp(idir,1:3)
      end do
      do idir = 1, 3
        write(*,'(1x,a,i0,a,3(1x,1pe13.5))') "[DG-MIXED-Z] occ-virt metric total row ", idir, " =", zcov_tot(idir,1:3)
      end do
      do idir = 1, 3
        write(*,'(1x,a,i0,a,3(1x,1pe13.5))') "[DG-MIXED-Z] Im Tr_occ [Z_i,Z_j] WW row ", idir, " =", zcomm_ww(idir,1:3)
      end do
      do idir = 1, 3
        write(*,'(1x,a,i0,a,3(1x,1pe13.5))') "[DG-MIXED-Z] Im Tr_occ [Z_i,Z_j] WP row ", idir, " =", zcomm_wp(idir,1:3)
      end do
      do idir = 1, 3
        write(*,'(1x,a,i0,a,3(1x,1pe13.5))') "[DG-MIXED-Z] Im Tr_occ [Z_i,Z_j] total row ", idir, " =", zcomm_tot(idir,1:3)
      end do
      if (use_comm_min_select) then
        write(*,'(1x,a,i0,a,i0,a,1pe13.5)') "[DG-MIXED-Z-COMM-MIN] final nraw=", count(active_pw), &
          " np=", dg_frag%mixed_wannier_bpw_np, " C_comm=", &
          abs(zcomm_tot(1,2)) + abs(zcomm_tot(2,3)) + abs(zcomm_tot(3,1))
      end if
      if (auto_report_enabled) call write_bpw_auto_report()
    end if

    deallocate(s_fp, r_fp, b_local, b_global, r_local, r_global, b_eig)
    deallocate(c_local, c_global, c_w_local, c_w_global, c_eig)
    deallocate(z_w, z_eig, o_eig)
    deallocate(s_perp, s_vec, s_work, s_eval, pw_ekin, shell_ekin, qmat, h_p, h_vec, h_eval)
    deallocate(q_metric, q_left_metric, final_metric, transform_metric, t_final, t_saved, metric_vec, diag_eval)
    deallocate(total_metric, h_total)
    deallocate(r_tilde, r_orth, r_final, z_wp_direct)
    deallocate(active_pw, trial_active_pw, best_active_pw, shell_key, shell_values)

  contains

    subroutine summarize_sperp_spectrum()
      integer :: i
      real(8) :: sval

      proj_min_selected = huge(1.0d0)
      proj_max_rejected = -huge(1.0d0)
      proj_gap = 0.0d0
      sv_min_selected = huge(1.0d0)
      sv_max_selected = 0.0d0
      sv_condition = huge(1.0d0)
      do i = 1, npw
        if (s_eval(i) > sperp_tol) then
          proj_min_selected = min(proj_min_selected, s_eval(i))
          sv_min_selected = min(sv_min_selected, sqrt(max(0.0d0, s_eval(i))))
          sv_max_selected = max(sv_max_selected, sqrt(max(0.0d0, s_eval(i))))
        else
          proj_max_rejected = max(proj_max_rejected, s_eval(i))
        end if
      end do
      if (proj_min_selected == huge(1.0d0)) proj_min_selected = 0.0d0
      if (proj_max_rejected == -huge(1.0d0)) proj_max_rejected = 0.0d0
      proj_gap = proj_min_selected - proj_max_rejected
      if (sv_min_selected == huge(1.0d0)) sv_min_selected = 0.0d0
      if (sv_min_selected > 0.0d0) sv_condition = sv_max_selected / sv_min_selected
      sval = dg_bpw_auto_accuracy
      if (sval <= 0.0d0) sval = 1.0d-3
    end subroutine summarize_sperp_spectrum

    subroutine write_bpw_auto_report()
      character(len=256) :: warning_line
      character(len=32) :: recommendation
      integer :: warning_count
      integer :: severe_count
      real(8) :: c_comm_sum

      report_iostat = 0
      warning_count = 0
      severe_count = 0
      recommendation = 'acceptable'
      c_comm_sum = abs(zcomm_tot(1,2)) + abs(zcomm_tot(2,3)) + abs(zcomm_tot(3,1))
      open(newunit=report_unit, file='dg_bpw_auto_report.dat', status='replace', action='write', iostat=report_iostat)
      if (report_iostat /= 0) then
        write(*,'(1x,a,i0)') '[DG-BPW-AUTO][WARN] failed to open dg_bpw_auto_report.dat; iostat=', report_iostat
        return
      end if

      write(report_unit,'(a)') '# DG-BPW auto basis report'
      write(report_unit,'(a)') '# First-stage implementation: report only; propagation basis is unchanged.'
      write(report_unit,'(a,a)') 'dg_bpw_auto = ', dg_bpw_auto
      write(report_unit,'(a,1pe20.10)') 'dg_bpw_auto_accuracy = ', dg_bpw_auto_accuracy
      write(report_unit,'(a,i0)') 'dg_bpw_auto_min_n = ', dg_bpw_auto_min_n
      write(report_unit,'(a,i0)') 'dg_bpw_auto_max_n = ', dg_bpw_auto_max_n
      write(report_unit,'(a,i0)') 'manual_n_plane_waves_dg = ', n_plane_waves_dg
      write(report_unit,'(a,a)') 'select_mode = ', trim(bpw_select_mode)
      if (bpw_ecut_is_set) then
        write(report_unit,'(a,1pe20.10)') 'candidate_Ecut = ', bpw_ecut
      else
        write(report_unit,'(a)') 'candidate_Ecut = all_raw'
      end if
      write(report_unit,'(a,i0)') 'candidate_BPW_number = ', npw
      write(report_unit,'(a,i0)') 'complete_shells_available = ', nshell
      write(report_unit,'(a,i0)') 'complete_shells_selected = ', selected_nshell
      write(report_unit,'(a,i0)') 'selected_G2_max = ', selected_g2_max
      write(report_unit,'(a,1pe20.10)') 'selected_Emax = ', selected_max_energy
      write(report_unit,'(a,i0)') 'selected_raw_BPW_number = ', nraw_selected
      write(report_unit,'(a,i0)') 'selected_perp_BPW_number = ', dg_frag%mixed_wannier_bpw_np
      write(report_unit,'(a,1pe20.10)') 'Sperp_cutoff_used = ', sperp_tol
      write(report_unit,'(a,1pe20.10)') 'singular_value_min_selected = ', sv_min_selected
      write(report_unit,'(a,1pe20.10)') 'singular_value_max_selected = ', sv_max_selected
      write(report_unit,'(a,1pe20.10)') 'singular_value_condition = ', sv_condition
      write(report_unit,'(a,1pe20.10)') 'Sperp_eigenvalue_min_all = ', min_s_local
      write(report_unit,'(a,1pe20.10)') 'Sperp_eigenvalue_max_all = ', max_s_local
      write(report_unit,'(a,1pe20.10)') 'projectability_min_selected = ', proj_min_selected
      write(report_unit,'(a,1pe20.10)') 'projectability_max_rejected = ', proj_max_rejected
      write(report_unit,'(a,1pe20.10)') 'projectability_gap = ', proj_gap
      write(report_unit,'(a,1pe20.10)') 'hermiticity_max = ', herm_max
      write(report_unit,'(a,3(1x,1pe20.10))') 'Z_WP_norm_xyz = ', zwp_norm(1:3)
      write(report_unit,'(a,3(1x,1pe20.10))') 'C_comm_xyz = ', &
        abs(zcomm_tot(1,2)), abs(zcomm_tot(2,3)), abs(zcomm_tot(3,1))
      write(report_unit,'(a,1pe20.10)') 'C_comm_sum = ', c_comm_sum
      write(report_unit,'(a,1pe20.10)') 'fsum_proxy_wp_avg = ', fsum_proxy_wp_avg
      write(report_unit,'(a,1pe20.10)') 'fsum_proxy_total_avg = ', fsum_proxy_total_avg
      write(report_unit,'(a)') 'warnings:'
      if (proj_gap < max(1.0d-12, dg_bpw_auto_accuracy)) then
        warning_count = warning_count + 1
        write(warning_line,'(a,1pe12.4)') '  Projectability boundary is nearly degenerate; gap=', proj_gap
        write(report_unit,'(a)') trim(warning_line)
      end if
      if (sv_min_selected > 0.0d0 .and. sv_condition > 1.0d0 / max(1.0d-12, dg_bpw_auto_accuracy)) then
        warning_count = warning_count + 1
        write(warning_line,'(a,1pe12.4)') '  BPW-perp overlap condition number is large; cond=', sv_condition
        write(report_unit,'(a)') trim(warning_line)
      end if
      if (dg_bpw_auto_min_n > 0 .and. dg_frag%mixed_wannier_bpw_np < dg_bpw_auto_min_n) then
        warning_count = warning_count + 1
        write(warning_line,'(a,i0,a,i0)') '  Selected BPW count below dg_bpw_auto_min_n: ', &
          dg_frag%mixed_wannier_bpw_np, ' < ', dg_bpw_auto_min_n
        write(report_unit,'(a)') trim(warning_line)
      end if
      if (dg_bpw_auto_max_n > 0 .and. dg_frag%mixed_wannier_bpw_np > dg_bpw_auto_max_n) then
        warning_count = warning_count + 1
        write(warning_line,'(a,i0,a,i0)') '  Selected BPW count exceeds dg_bpw_auto_max_n: ', &
          dg_frag%mixed_wannier_bpw_np, ' > ', dg_bpw_auto_max_n
        write(report_unit,'(a)') trim(warning_line)
      end if
      if (dg_frag%mixed_wannier_bpw_np <= 0) then
        warning_count = warning_count + 1
        write(report_unit,'(a)') '  No BPW-perp states selected.'
      end if
      write(report_unit,'(a)') 'end_warnings'
      if (dg_frag%mixed_wannier_bpw_np <= 0) severe_count = severe_count + 1
      if (sv_condition > 1.0d8) severe_count = severe_count + 1
      if (sv_min_selected < 1.0d-12) severe_count = severe_count + 1
      if (dg_bpw_auto_min_n > 0 .and. dg_frag%mixed_wannier_bpw_np < dg_bpw_auto_min_n) severe_count = severe_count + 1
      if (dg_bpw_auto_max_n > 0 .and. dg_frag%mixed_wannier_bpw_np > dg_bpw_auto_max_n) severe_count = severe_count + 1
      if (severe_count > 0) then
        recommendation = 'rejected'
      else if (warning_count > 0 .or. sv_condition > 1.0d3 .or. &
               sv_min_selected < max(1.0d-12, dg_bpw_auto_accuracy) .or. &
               proj_gap < max(1.0d-12, dg_bpw_auto_accuracy) .or. &
               fsum_proxy_total_avg < 1.0d1 .or. c_comm_sum > 5.0d1) then
        recommendation = 'risky'
      else if (fsum_proxy_total_avg >= 2.0d1 .and. c_comm_sum <= 3.0d1) then
        recommendation = 'recommended'
      else
        recommendation = 'acceptable'
      end if
      write(report_unit,'(a,a)') 'recommendation = ', trim(recommendation)
      write(report_unit,'(a)') 'recommendation_reasons:'
      if (severe_count > 0) write(report_unit,'(a,i0)') '  severe numerical issues: ', severe_count
      if (warning_count > 0) write(report_unit,'(a,i0)') '  warnings: ', warning_count
      if (sv_condition > 1.0d3) write(report_unit,'(a,1pe12.4)') '  condition number is above risky threshold: ', sv_condition
      if (sv_min_selected < max(1.0d-12, dg_bpw_auto_accuracy)) &
        write(report_unit,'(a,1pe12.4)') '  sigma_min is below accuracy threshold: ', sv_min_selected
      if (proj_gap < max(1.0d-12, dg_bpw_auto_accuracy)) &
        write(report_unit,'(a,1pe12.4)') '  projectability gap is small: ', proj_gap
      if (fsum_proxy_total_avg < 1.0d1) write(report_unit,'(a,1pe12.4)') &
        '  f-sum proxy is low: ', fsum_proxy_total_avg
      if (c_comm_sum > 5.0d1) write(report_unit,'(a,1pe12.4)') &
        '  commutator diagnostic is large: ', c_comm_sum
      if (trim(recommendation) == 'recommended') write(report_unit,'(a)') &
        '  condition, sigma_min, f-sum proxy, and commutator diagnostics are within first-stage thresholds.'
      if (trim(recommendation) == 'acceptable') write(report_unit,'(a)') &
        '  numerically usable, but one or more diagnostics are not strong enough for recommended.'
      write(report_unit,'(a)') 'end_recommendation_reasons'
      write(report_unit,'(a,i0,a,i0,a,i0,a,1pe20.10,a,1pe20.10,a,1pe20.10,a,1pe20.10,a,1pe20.10,a,i0,a,a)') &
        'SUMMARY candidate_n=', npw, &
        ' selected_raw_n=', nraw_selected, &
        ' selected_perp_n=', dg_frag%mixed_wannier_bpw_np, &
        ' shell_ecut=', selected_max_energy, &
        ' sigma_min=', sv_min_selected, &
        ' sigma_max=', sv_max_selected, &
        ' cond=', sv_condition, &
        ' proj_gap=', proj_gap, &
        ' warnings=', warning_count, &
        ' recommendation=', trim(recommendation)
      close(report_unit)
      write(*,'(1x,a)') '[DG-BPW-AUTO] wrote dg_bpw_auto_report.dat'
    end subroutine write_bpw_auto_report

    subroutine build_bpw_shells(nshell_out)
      integer, intent(out) :: nshell_out
      integer :: i, j, key_tmp
      real(8), parameter :: twopi = 6.283185307179586476925286766559d0
      real(8) :: lbox(3)
      integer :: ik(3)

      do j = 1, 3
        lbox(j) = max(1.0d-12, dg_frag%hgs(j) * dble(max(1, dg_frag%nxyz_domain(j,1))))
      end do
      do i = 1, npw
        do j = 1, 3
          ik(j) = nint(dg_frag%k_pw(j,i) * lbox(j) / twopi)
        end do
        shell_key(i) = ik(1)*ik(1) + ik(2)*ik(2) + ik(3)*ik(3)
      end do
      shell_values(1:npw) = shell_key(1:npw)
      do i = 2, npw
        key_tmp = shell_values(i)
        j = i - 1
        do while (j >= 1 .and. shell_values(j) > key_tmp)
          shell_values(j+1) = shell_values(j)
          j = j - 1
        end do
        shell_values(j+1) = key_tmp
      end do
      nshell_out = 0
      do i = 1, npw
        if (i == 1 .or. shell_values(i) /= shell_values(i-1)) then
          nshell_out = nshell_out + 1
          shell_values(nshell_out) = shell_values(i)
        end if
      end do
      shell_ekin(:) = 0.0d0
      do j = 1, nshell_out
        do i = 1, npw
          if (shell_key(i) == shell_values(j)) shell_ekin(j) = max(shell_ekin(j), pw_ekin(i))
        end do
      end do
    end subroutine build_bpw_shells

    subroutine evaluate_bpw_active_set(mask, nraw_eval, np_eval, c_comm_eval, fwp_avg_eval, min_s_eval_out, &
      max_s_eval_out, fwp_axis_eval, comm_pair_eval)
      logical, intent(in) :: mask(:)
      integer, intent(in) :: nraw_eval
      integer, intent(out) :: np_eval
      real(8), intent(out) :: c_comm_eval, fwp_avg_eval, min_s_eval_out, max_s_eval_out
      real(8), intent(out) :: fwp_axis_eval(3), comm_pair_eval(3)
      integer :: ispin_eval
      integer :: occ, a, v
      real(8) :: fwp_axis(3), zcomm_eval(3,3), zcomm_ww_eval(3,3), zcomm_wp_eval(3,3)
      complex(8), allocatable :: zww_eval(:,:,:), zwp_eval(:,:,:)

      np_eval = 0
      c_comm_eval = huge(1.0d0)
      fwp_avg_eval = 0.0d0
      min_s_eval_out = 0.0d0
      max_s_eval_out = 0.0d0
      fwp_axis_eval(:) = 0.0d0
      comm_pair_eval(:) = huge(1.0d0)
      if (nraw_eval <= 0) return

      ispin_eval = 1
      nocc_ref = min(neig - 1, max(1, dg_frag%nstate_tot))
      if (allocated(dg_frag%nocc_spin)) then
        if (ispin_eval <= size(dg_frag%nocc_spin)) nocc_ref = min(neig - 1, max(1, dg_frag%nocc_spin(ispin_eval)))
      end if

      b_local(:, :) = zzero
      do i_local = 1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1
        ifrag = dg_frag%ifrag_start + i_local - 1
        if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
        nbf = min(dg_frag%n_basis(ifrag, ispin_eval), size(dg_frag%global_wannier_coef, 1))
        do ibasis = 1, nbf
          global_row = dg_frag%index_basis(ibasis, ifrag, ispin_eval)
          if (global_row < 1 .or. global_row > nmat) cycle
          local_row = 0
          if (allocated(dg_frag%coef_global_to_local)) local_row = dg_frag%coef_global_to_local(global_row, ispin_eval)
          if (local_row < 1 .or. local_row > nmat) cycle
          do iw = 1, nwann
            coeff = conjg(dg_frag%global_wannier_coef(ibasis, iw, ispin_eval, i_local))
            b_local(iw,1:npw) = b_local(iw,1:npw) + coeff * s_fp(local_row,1:npw,ispin_eval)
          end do
        end do
      end do
      call comm_summation(b_local, b_global, nwann*npw, dg_frag%icomm)
      do ipw = 1, npw
        if (.not. mask(ipw)) b_global(:,ipw) = zzero
      end do

      s_perp(:, :) = zzero
      do ipw = 1, npw
        if (mask(ipw)) s_perp(ipw,ipw) = (1.0d0, 0.0d0)
      end do
      s_perp(:, :) = s_perp(:, :) - matmul(conjg(transpose(b_global)), b_global)
      s_vec = s_perp
      call eigen_zheev(s_vec, s_eval, s_work)
      np_eval = count(s_eval > sperp_tol)
      min_s_eval_out = minval(s_eval)
      max_s_eval_out = maxval(s_eval)
      if (np_eval <= 0) return

      qmat(:, :) = zzero
      pslot = 0
      do alpha = 1, npw
        if (s_eval(alpha) <= sperp_tol) cycle
        pslot = pslot + 1
        if (pslot > np_eval) exit
        qmat(1:npw,pslot) = s_work(1:npw,alpha) / sqrt(s_eval(alpha))
      end do

      h_p(:, :) = zzero
      do alpha = 1, np_eval
        do beta = 1, np_eval
          do ipw = 1, npw
            h_p(alpha,beta) = h_p(alpha,beta) + conjg(qmat(ipw,alpha)) * pw_ekin(ipw) * qmat(ipw,beta)
          end do
        end do
      end do
      h_vec(1:np_eval,1:np_eval) = h_p(1:np_eval,1:np_eval)
      call eigen_zheev(h_vec(1:np_eval,1:np_eval), h_eval(1:np_eval), h_p(1:np_eval,1:np_eval))
      pw_energy_ref = dg_frag%global_wannier_flux_eval(nocc_ref+1,ispin_eval)
      h_eval(1:np_eval) = h_eval(1:np_eval) + pw_energy_ref

      allocate(zww_eval(3,neig,neig), zwp_eval(3,neig,np_eval))
      zww_eval(:, :, :) = zzero
      zwp_eval(:, :, :) = zzero

      b_eig(:, :) = zzero
      do eig = 1, neig
        do ipw = 1, npw
          do iw = 1, nwann
            b_eig(eig,ipw) = b_eig(eig,ipw) + conjg(dg_frag%global_wannier_flux_evec(iw,eig)) * b_global(iw,ipw)
          end do
        end do
      end do

      do idir = 1, 3
        z_w(1:nwann,1:nwann) = dg_frag%global_wannier_position(idir,1:nwann,1:nwann)
        z_eig(1:neig,1:neig) = matmul(conjg(transpose(dg_frag%global_wannier_flux_evec(1:nwann,1:neig))), &
          matmul(z_w, dg_frag%global_wannier_flux_evec(1:nwann,1:neig)))
        zww_eval(idir,1:neig,1:neig) = z_eig(1:neig,1:neig)

        r_local(:, :) = zzero
        c_local(:, :) = zzero
        c_w_local(:, :) = zzero
        do i_local = 1, dg_frag%ifrag_end - dg_frag%ifrag_start + 1
          ifrag = dg_frag%ifrag_start + i_local - 1
          if (i_local < 1 .or. i_local > size(dg_frag%global_wannier_coef, 4)) cycle
          frag_center(1:3) = 0.0d0
          if (allocated(dg_frag%buffer_wannier_frag_center) .and. i_local <= size(dg_frag%buffer_wannier_frag_center, 2)) then
            frag_center(1:3) = dg_frag%buffer_wannier_frag_center(1:3,i_local)
          else
            do jdir = 1, 3
              frag_center(jdir) = dg_frag%hgs(jdir) * (dble(dg_frag%ixyz_frag(jdir,ifrag) - 1) + &
                0.5d0 * dble(dg_frag%nxyz_domain(jdir,ifrag) - 1))
            end do
          end if
          nbf = min(dg_frag%n_basis(ifrag, ispin_eval), size(dg_frag%global_wannier_coef, 1))
          do ibasis = 1, nbf
            global_row = dg_frag%index_basis(ibasis, ifrag, ispin_eval)
            if (global_row < 1 .or. global_row > nmat) cycle
            local_row = 0
            if (allocated(dg_frag%coef_global_to_local)) local_row = dg_frag%coef_global_to_local(global_row, ispin_eval)
            if (local_row < 1 .or. local_row > nmat) cycle
            do iw = 1, nwann
              coeff = conjg(dg_frag%global_wannier_coef(ibasis, iw, ispin_eval, i_local))
              r_local(iw,1:npw) = r_local(iw,1:npw) + coeff * r_fp(idir,local_row,1:npw,ispin_eval)
              c_local(iw,1:npw) = c_local(iw,1:npw) + coeff * frag_center(idir) * s_fp(local_row,1:npw,ispin_eval)
              do jw = 1, nwann
                c_w_local(iw,jw) = c_w_local(iw,jw) + coeff * frag_center(idir) * &
                  dg_frag%global_wannier_coef(ibasis, jw, ispin_eval, i_local)
              end do
            end do
          end do
        end do
        call comm_summation(r_local, r_global, nwann*npw, dg_frag%icomm)
        call comm_summation(c_local, c_global, nwann*npw, dg_frag%icomm)
        call comm_summation(c_w_local, c_w_global, nwann*nwann, dg_frag%icomm)
        c_eig(1:neig,1:neig) = matmul(conjg(transpose(dg_frag%global_wannier_flux_evec(1:nwann,1:neig))), &
          matmul(c_w_global(1:nwann,1:nwann), dg_frag%global_wannier_flux_evec(1:nwann,1:neig)))
        if (use_fragment_center_direct) then
          o_eig(1:neig,1:neig) = z_eig(1:neig,1:neig) - c_eig(1:neig,1:neig)
        else
          o_eig(1:neig,1:neig) = z_eig(1:neig,1:neig)
        end if
        r_tilde(:, :) = zzero
        do eig = 1, neig
          do ipw = 1, npw
            amp = zzero
            do iw = 1, nwann
              amp = amp + conjg(dg_frag%global_wannier_flux_evec(iw,eig)) * r_global(iw,ipw)
            end do
            if (use_fragment_center_direct) then
              do iw = 1, nwann
                amp = amp - conjg(dg_frag%global_wannier_flux_evec(iw,eig)) * c_global(iw,ipw)
              end do
            end if
            grad_corr = zzero
            do jw = 1, neig
              grad_corr = grad_corr + o_eig(eig,jw) * b_eig(jw,ipw)
            end do
            r_tilde(eig,ipw) = amp - grad_corr
          end do
        end do
        r_orth(1:neig,1:np_eval) = matmul(r_tilde(1:neig,1:npw), qmat(1:npw,1:np_eval))
        if (use_h_eigenvectors_for_final) then
          r_final(1:neig,1:np_eval) = matmul(r_orth(1:neig,1:np_eval), h_p(1:np_eval,1:np_eval))
        else
          r_final(1:neig,1:np_eval) = matmul(r_orth(1:neig,1:np_eval), h_vec(1:np_eval,1:np_eval))
        end if
        zwp_eval(idir,1:neig,1:np_eval) = r_final(1:neig,1:np_eval)
      end do

      fwp_axis(:) = 0.0d0
      zcomm_ww_eval(:, :) = 0.0d0
      zcomm_wp_eval(:, :) = 0.0d0
      do idir = 1, 3
        do occ = 1, nocc_ref
          do a = 1, np_eval
            gap = h_eval(a) - dg_frag%global_wannier_flux_eval(occ,ispin_eval)
            if (gap > 1.0d-12) fwp_axis(idir) = fwp_axis(idir) + 2.0d0 * gap * abs(zwp_eval(idir,occ,a))**2
          end do
        end do
      end do
      do idir = 1, 3
        do jdir = 1, 3
          do occ = 1, nocc_ref
            do v = nocc_ref + 1, neig
              zcomm_ww_eval(idir,jdir) = zcomm_ww_eval(idir,jdir) + aimag( &
                zww_eval(idir,occ,v) * conjg(zww_eval(jdir,occ,v)) - &
                zww_eval(jdir,occ,v) * conjg(zww_eval(idir,occ,v)))
            end do
            do a = 1, np_eval
              zcomm_wp_eval(idir,jdir) = zcomm_wp_eval(idir,jdir) + aimag( &
                zwp_eval(idir,occ,a) * conjg(zwp_eval(jdir,occ,a)) - &
                zwp_eval(jdir,occ,a) * conjg(zwp_eval(idir,occ,a)))
            end do
          end do
        end do
      end do
      zcomm_eval(:, :) = zcomm_ww_eval(:, :) + zcomm_wp_eval(:, :)
      comm_pair_eval(1) = abs(zcomm_eval(2,3))
      comm_pair_eval(2) = abs(zcomm_eval(3,1))
      comm_pair_eval(3) = abs(zcomm_eval(1,2))
      c_comm_eval = sum(comm_pair_eval(1:3))
      fwp_axis_eval(:) = fwp_axis(:)
      fwp_avg_eval = sum(fwp_axis(1:3)) / (3.0d0 * max(1.0d0, dble(nocc_ref)))

      deallocate(zww_eval, zwp_eval)
    end subroutine evaluate_bpw_active_set
  end subroutine ensure_mixed_wannier_bpw_position

  subroutine maybe_build_mixed_wannier_bpw_position(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(len=32) :: env_value
    integer :: env_len, env_stat

    env_value = ''
    call get_environment_variable('SALMON_DG_MIXED_Z', env_value, length=env_len, status=env_stat)
    if (env_stat /= 0 .or. env_len <= 0) return
    select case (adjustl(trim(env_value(1:env_len))))
    case ('1','y','Y','yes','YES','true','TRUE','on','ON')
      call ensure_mixed_wannier_bpw_position(dg_frag)
    end select
  end subroutine maybe_build_mixed_wannier_bpw_position

  subroutine diagnose_mixed_wannier_pw_fsum(dg_frag)
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    real(8), parameter :: au_to_ev = 27.211386245988d0
    character(len=32) :: env_diag
    integer :: env_len, env_stat
    logical :: enabled
    integer :: nocc
    integer :: ispin, idir, occ, virt, alpha
    integer :: iwin, pairs_ww, pairs_wp, best_ww_occ, best_ww_virt, best_wp_occ, best_wp_virt
    real(8) :: gap, fsum_ww, fsum_wp, r2sum_ww, r2sum_wp
    real(8) :: max_ww, max_wp, max_ww_gap, max_wp_gap
    real(8) :: mean_gap_ww, mean_gap_wp, gap_weight_ww, gap_weight_wp
    real(8) :: win_lo(4), win_hi(4), win_r2_ww, win_r2_wp, win_f_ww, win_f_wp
    real(8) :: win_max_ww, win_max_wp, win_max_ww_gap, win_max_wp_gap
    complex(8) :: amp

    env_diag = ''
    call get_environment_variable('SALMON_DG_MIXED_FSUM', env_diag, length=env_len, status=env_stat)
    enabled = .false.
    if (env_stat == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_diag(1:env_len))))
      case ('1','y','Y','yes','YES','true','TRUE','on','ON')
        enabled = .true.
      end select
    end if
    if (.not. enabled) return

    win_lo = (/0.0d0, 3.0d0, 3.6d0, 0.0d0/)
    win_hi = (/1.0d0, 3.6d0, 5.0d0, 20.0d0/)

    if (.not. dg_frag%use_plane_wave_basis .or. dg_frag%n_plane_waves <= 0) then
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[DG-MIXED-FSUM] skipped: plane-wave basis is disabled."
      return
    end if
    if (.not. dg_frag%has_global_wannier_flux_eigen) then
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[DG-MIXED-FSUM] skipped: global Wannier Flux eigen seed is unavailable."
      return
    end if
    if (.not. dg_frag%has_global_wannier_position .or. .not. allocated(dg_frag%global_wannier_position)) then
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[DG-MIXED-FSUM] skipped: global Wannier AA_R position matrix is unavailable."
      return
    end if
    if (.not. allocated(dg_frag%global_wannier_coef) .or. .not. allocated(dg_frag%global_wannier_flux_evec) .or. &
        .not. allocated(dg_frag%global_wannier_flux_eval)) return

    call ensure_mixed_wannier_bpw_position(dg_frag)
    if (dg_frag%has_mixed_wannier_bpw_position .and. allocated(dg_frag%mixed_wannier_bpw_z) .and. &
        allocated(dg_frag%mixed_wannier_bpw_eval)) then
      do ispin = 1, min(dg_frag%nspin, size(dg_frag%mixed_wannier_bpw_eval, 2))
        nocc = min(dg_frag%mixed_wannier_bpw_nw, dg_frag%nstate_tot)
        if (allocated(dg_frag%nocc_spin)) then
          if (ispin <= size(dg_frag%nocc_spin)) nocc = min(dg_frag%mixed_wannier_bpw_nw, max(0, dg_frag%nocc_spin(ispin)))
        end if
        if (nocc <= 0) cycle
        do idir = 1, 3
          fsum_ww = 0.0d0
          fsum_wp = 0.0d0
          r2sum_ww = 0.0d0
          r2sum_wp = 0.0d0
          max_ww = 0.0d0
          max_wp = 0.0d0
          max_ww_gap = 0.0d0
          max_wp_gap = 0.0d0
          gap_weight_ww = 0.0d0
          gap_weight_wp = 0.0d0
          do virt = nocc + 1, dg_frag%mixed_wannier_bpw_nw
            do occ = 1, nocc
              gap = dg_frag%mixed_wannier_bpw_eval(virt,ispin) - dg_frag%mixed_wannier_bpw_eval(occ,ispin)
              if (gap <= 1.0d-12) cycle
              amp = dg_frag%mixed_wannier_bpw_z(idir,occ,virt,ispin)
              r2sum_ww = r2sum_ww + abs(amp)**2
              fsum_ww = fsum_ww + 2.0d0 * gap * abs(amp)**2
              gap_weight_ww = gap_weight_ww + gap * abs(amp)**2
              if (abs(amp)**2 > max_ww) then
                max_ww = abs(amp)**2
                max_ww_gap = gap
              end if
            end do
          end do
          do alpha = 1, dg_frag%mixed_wannier_bpw_np
            virt = dg_frag%mixed_wannier_bpw_nw + alpha
            do occ = 1, nocc
              gap = dg_frag%mixed_wannier_bpw_eval(virt,ispin) - dg_frag%mixed_wannier_bpw_eval(occ,ispin)
              if (gap <= 1.0d-12) cycle
              amp = dg_frag%mixed_wannier_bpw_z(idir,occ,virt,ispin)
              r2sum_wp = r2sum_wp + abs(amp)**2
              fsum_wp = fsum_wp + 2.0d0 * gap * abs(amp)**2
              gap_weight_wp = gap_weight_wp + gap * abs(amp)**2
              if (abs(amp)**2 > max_wp) then
                max_wp = abs(amp)**2
                max_wp_gap = gap
              end if
            end do
          end do
          mean_gap_ww = au_to_ev * gap_weight_ww / max(1.0d-300, r2sum_ww)
          mean_gap_wp = au_to_ev * gap_weight_wp / max(1.0d-300, r2sum_wp)
          if (comm_is_root(dg_frag%id)) then
            write(*,'(1x,a,i0,a,i0,a,i0,13(a,1pe13.5))') "[DG-MIXED-FSUM] ispin=", ispin, &
              " axis=", idir, " n_pperp=", dg_frag%mixed_wannier_bpw_np, &
              " f_WW_per_occ=", fsum_ww / max(1.0d0, dble(nocc)), &
              " f_WP_per_occ=", fsum_wp / max(1.0d0, dble(nocc)), &
              " f_tot_per_occ=", (fsum_ww + fsum_wp) / max(1.0d0, dble(nocc)), &
              " r2_WW=", r2sum_ww, &
              " r2_WP=", r2sum_wp, &
              " max_WW_r2=", max_ww, &
              " max_WW_gap_eV=", au_to_ev * max_ww_gap, &
              " max_WP_r2=", max_wp, &
              " max_WP_gap_eV=", au_to_ev * max_wp_gap, &
              " mean_WW_gap_eV=", mean_gap_ww, &
              " mean_WP_gap_eV=", mean_gap_wp, &
              " min_Sperp=", dg_frag%mixed_wannier_bpw_sperp_min, &
              " max_Sperp=", dg_frag%mixed_wannier_bpw_sperp_max
          end if
          do iwin = 1, size(win_lo)
            win_r2_ww = 0.0d0
            win_r2_wp = 0.0d0
            win_f_ww = 0.0d0
            win_f_wp = 0.0d0
            win_max_ww = 0.0d0
            win_max_wp = 0.0d0
            win_max_ww_gap = 0.0d0
            win_max_wp_gap = 0.0d0
            pairs_ww = 0
            pairs_wp = 0
            best_ww_occ = 0
            best_ww_virt = 0
            best_wp_occ = 0
            best_wp_virt = 0
            do virt = nocc + 1, dg_frag%mixed_wannier_bpw_nw
              do occ = 1, nocc
                gap = dg_frag%mixed_wannier_bpw_eval(virt,ispin) - dg_frag%mixed_wannier_bpw_eval(occ,ispin)
                if (gap <= 1.0d-12) cycle
                if (au_to_ev * gap < win_lo(iwin) .or. au_to_ev * gap > win_hi(iwin)) cycle
                amp = dg_frag%mixed_wannier_bpw_z(idir,occ,virt,ispin)
                pairs_ww = pairs_ww + 1
                win_r2_ww = win_r2_ww + abs(amp)**2
                win_f_ww = win_f_ww + 2.0d0 * gap * abs(amp)**2
                if (abs(amp)**2 > win_max_ww) then
                  win_max_ww = abs(amp)**2
                  win_max_ww_gap = gap
                  best_ww_occ = occ
                  best_ww_virt = virt
                end if
              end do
            end do
            do alpha = 1, dg_frag%mixed_wannier_bpw_np
              virt = dg_frag%mixed_wannier_bpw_nw + alpha
              do occ = 1, nocc
                gap = dg_frag%mixed_wannier_bpw_eval(virt,ispin) - dg_frag%mixed_wannier_bpw_eval(occ,ispin)
                if (gap <= 1.0d-12) cycle
                if (au_to_ev * gap < win_lo(iwin) .or. au_to_ev * gap > win_hi(iwin)) cycle
                amp = dg_frag%mixed_wannier_bpw_z(idir,occ,virt,ispin)
                pairs_wp = pairs_wp + 1
                win_r2_wp = win_r2_wp + abs(amp)**2
                win_f_wp = win_f_wp + 2.0d0 * gap * abs(amp)**2
                if (abs(amp)**2 > win_max_wp) then
                  win_max_wp = abs(amp)**2
                  win_max_wp_gap = gap
                  best_wp_occ = occ
                  best_wp_virt = alpha
                end if
              end do
            end do
            if (comm_is_root(dg_frag%id)) then
              write(*,*) '[DG-MIXED-FSUM-WINDOW]', &
                ' ispin=', ispin, &
                ' axis=', idir, &
                ' window_eV=[', win_lo(iwin), ',', win_hi(iwin), &
                '] pairs_WW=', pairs_ww, &
                ' pairs_WP=', pairs_wp, &
                ' r2_WW=', win_r2_ww, &
                ' r2_WP=', win_r2_wp, &
                ' f_WW=', win_f_ww, &
                ' f_WP=', win_f_wp, &
                ' max_WW_abs_r=', sqrt(win_max_ww), &
                ' max_WW_gap_eV=', au_to_ev * win_max_ww_gap, &
                ' WW_occ=', best_ww_occ, &
                ' WW_virt=', best_ww_virt, &
                ' max_WP_abs_r=', sqrt(win_max_wp), &
                ' max_WP_gap_eV=', au_to_ev * win_max_wp_gap, &
                ' WP_occ=', best_wp_occ, &
                ' WP_alpha=', best_wp_virt
            end if
          end do
        end do
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,1pe13.5,a,1pe13.5,a,i0)') "[DG-MIXED-FSUM-ORTH] min_Sperp=", &
            dg_frag%mixed_wannier_bpw_sperp_min, " max_Sperp=", dg_frag%mixed_wannier_bpw_sperp_max, &
            " nraw_pw=", dg_frag%n_plane_waves
        end if
      end do
      return
    end if

    if (comm_is_root(dg_frag%id)) &
      write(*,'(1x,a)') "[DG-MIXED-FSUM] skipped: mixed direct-Z matrix was not built."
    return
  end subroutine diagnose_mixed_wannier_pw_fsum
