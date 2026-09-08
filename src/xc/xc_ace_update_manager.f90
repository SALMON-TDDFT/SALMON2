module xc_ace_update_manager
  use structures, only: s_parallel_info, s_rgrid, s_reciprocal_grid, s_poisson
  use communication, only: comm_is_root, comm_bcast, comm_summation
  use xc_hse_grid_sr, only: ace_first_stage_update_trigger, build_exchange_W_truth, build_ACE_operator, &
                            apply_exchange_ACE, ace_stage2_residual_check
  implicit none

  type ace_update_state
    logical :: ace_enabled = .false.
    logical :: ace_built = .false.
    integer :: last_update_step = -huge(1)/2
    integer :: n_min = 3
    integer :: n_max = 30
    integer :: n_max_base = 30
    real(8) :: eps_d = 1.0d-3
    real(8) :: eps_R = 1.0d-3
    real(8) :: eps = 1.0d-14
    logical :: need_second_stage_check = .false.
    logical :: force_update = .false.
    logical :: rebuild_done = .false.
    integer :: log_level = 1
    integer :: nv = 6
    integer :: topm = 6
    integer :: nrand = 2
    integer :: seed = 12345
    integer :: ns_max = 16
    real(8) :: omega = 0.d0
    real(8) :: a = 0.d0
    real(8) :: d_max = 0.d0
    real(8) :: R_max = 0.d0
    integer :: i_max = 0
    logical :: log_ex = .false.
    integer :: ex_stride = 10
    real(8) :: diag_shift = 0.d0
    real(8) :: ksr_min = 1.0d-12
    real(8) :: stability_min_diag = 0.d0
    real(8) :: stability_cond_proxy = huge(1.0d0)
    complex(8), allocatable :: W(:,:,:,:,:)
    complex(8), allocatable :: Minv(:,:)
    real(8), allocatable :: d_values(:)
  end type ace_update_state

  interface ace_update_decision
    module procedure ace_update_decision_zwf7
    module procedure ace_update_decision_zwf5
  end interface ace_update_decision

contains

  subroutine ace_update_init_from_env(state)
    implicit none
    type(ace_update_state), intent(inout) :: state
    character(64) :: env
    integer :: ist, ios, iv
    real(8) :: rv

    state%ace_enabled = .false.
    state%ace_built = .false.
    state%last_update_step = -huge(1)/2
    state%n_min = 3
    state%n_max = 30
    state%n_max_base = 30
    state%eps_d = 1.0d-3
    state%eps_R = 1.0d-3
    state%eps = 1.0d-14
    state%need_second_stage_check = .false.
    state%force_update = .false.
    state%rebuild_done = .false.
    state%log_level = 1
    state%nv = 6
    state%topm = 6
    state%nrand = 2
    state%seed = 12345
    state%ns_max = 16
    state%omega = 0.d0
    state%a = 0.d0
    state%d_max = 0.d0
    state%R_max = 0.d0
    state%i_max = 0
    state%log_ex = .false.
    state%ex_stride = 10
    state%diag_shift = 0.d0
    state%ksr_min = 1.0d-12
    state%stability_min_diag = 0.d0
    state%stability_cond_proxy = huge(1.0d0)

    env = ''
    call get_environment_variable('SALMON_ACE_RT', env, status=ist)
    if (ist == 0) then
      select case(trim(adjustl(env)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        state%ace_enabled = .true.
      case default
        state%ace_enabled = .false.
      end select
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_EPSD', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) rv
      if (ios == 0) state%eps_d = max(0.d0, rv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_EPSR', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) rv
      if (ios == 0) state%eps_R = max(0.d0, rv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_NMIN', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) iv
      if (ios == 0) state%n_min = max(1, iv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_NMAX', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) iv
      if (ios == 0) state%n_max = max(1, iv)
    end if
    state%n_max_base = state%n_max

    env = ''
    call get_environment_variable('SALMON_ACE_NV', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) iv
      if (ios == 0) state%nv = max(0, iv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_TOPM', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) iv
      if (ios == 0) state%topm = max(0, iv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_NRAND', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) iv
      if (ios == 0) state%nrand = max(0, iv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_SEED', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) iv
      if (ios == 0) state%seed = max(1, iv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_LOG_EX', env, status=ist)
    if (ist == 0) then
      select case(trim(adjustl(env)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        state%log_ex = .true.
      case default
        state%log_ex = .false.
      end select
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_EX_STRIDE', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) iv
      if (ios == 0) state%ex_stride = max(1, iv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_DIAG_SHIFT', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) rv
      if (ios == 0) state%diag_shift = max(0.d0, rv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_KSR_MIN', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) rv
      if (ios == 0) state%ksr_min = max(0.d0, rv)
    end if

  end subroutine ace_update_init_from_env


  subroutine ace_update_decision_zwf7(step, psi_prev, psi_curr, state, info, hvol, comm)
    implicit none
    integer, intent(in) :: step
    complex(8), intent(in) :: psi_prev(:,:,:,:,:,:,:)
    complex(8), intent(in) :: psi_curr(:,:,:,:,:,:,:)
    type(ace_update_state), intent(inout) :: state
    type(s_parallel_info), intent(in) :: info
    real(8), intent(in), optional :: hvol
    integer, intent(in), optional :: comm

    call ace_update_decision_core(step, &
      psi_prev(:,:,:,:,info%io_s:info%io_e,info%ik_s,info%im_s), &
      psi_curr(:,:,:,:,info%io_s:info%io_e,info%ik_s,info%im_s), &
      state, info, hvol, comm)
  end subroutine ace_update_decision_zwf7


  subroutine ace_update_decision_zwf5(step, psi_prev, psi_curr, state, info, hvol, comm)
    implicit none
    integer, intent(in) :: step
    complex(8), intent(in) :: psi_prev(:,:,:,:,:)
    complex(8), intent(in) :: psi_curr(:,:,:,:,:)
    type(ace_update_state), intent(inout) :: state
    type(s_parallel_info), intent(in) :: info
    real(8), intent(in), optional :: hvol
    integer, intent(in), optional :: comm

    call ace_update_decision_core(step, psi_prev, psi_curr, state, info, hvol, comm)
  end subroutine ace_update_decision_zwf5


  subroutine ace_update_decision_core(step, phi_prev, phi_curr, state, info, hvol, comm)
    implicit none
    integer, intent(in) :: step
    complex(8), intent(in) :: phi_prev(:,:,:,:,:)
    complex(8), intent(in) :: phi_curr(:,:,:,:,:)
    type(ace_update_state), intent(inout) :: state
    type(s_parallel_info), intent(in) :: info
    real(8), intent(in), optional :: hvol
    integer, intent(in), optional :: comm

    integer :: gap, root_id
    logical :: trigger_true, allowed_check
    real(8) :: d_max_local, topk_vals(5)
    integer :: i_max_local, topk_idx(5)

    if (.not. state%ace_enabled) return

    d_max_local = 0.d0
    i_max_local = 0
    topk_vals = 0.d0
    topk_idx = 0
    state%force_update = .false.
    state%need_second_stage_check = .false.
    trigger_true = .false.

    if (present(comm)) continue
    root_id = info%id_rko

    if (state%last_update_step <= -huge(1)/4) then
      gap = huge(1)/8
    else
      gap = step - state%last_update_step
    end if

    allowed_check = (gap >= state%n_min)
    if (state%n_max > 0) state%force_update = (gap >= state%n_max)

    if (allowed_check) then
      call ace_first_stage_update_trigger(info, phi_prev, phi_curr, maxval(info%io_e_all), trigger_true, &
           d_max=d_max_local, i_max=i_max_local, d_topk=topk_vals, i_topk=topk_idx, &
           step=step, eps_d=state%eps_d, eps=state%eps, hvol=hvol, &
           print_summary=.false., print_topk=.false.)
    end if

    state%need_second_stage_check = state%force_update .or. trigger_true

    if (comm_is_root(root_id)) then
      write(*,'(A,I8,A,1PE11.3,A,I8,A,L1)') &
        'ACE manager step=', step, ' d_max=', d_max_local, ' i_max=', i_max_local, &
        ' need_second_stage=', state%need_second_stage_check
      if (state%need_second_stage_check) then
        write(*,'(A)', advance='no') 'ACE1 top-k:'
        write(*,'(A,I8,A,1PE12.4)', advance='no') ' (i=', topk_idx(1), ',d=', topk_vals(1), ')'
        write(*,'(A,I8,A,1PE12.4)', advance='no') ' (i=', topk_idx(2), ',d=', topk_vals(2), ')'
        write(*,'(A,I8,A,1PE12.4)', advance='no') ' (i=', topk_idx(3), ',d=', topk_vals(3), ')'
        write(*,'(A,I8,A,1PE12.4)', advance='no') ' (i=', topk_idx(4), ',d=', topk_vals(4), ')'
        write(*,'(A,I8,A,1PE12.4)')               ' (i=', topk_idx(5), ',d=', topk_vals(5), ')'
        write(*,'(A)') 'ACE second-stage/rebuild would run here (not implemented).'
      end if
    end if

  end subroutine ace_update_decision_core


  subroutine ace_update_step_rt(step, lg, mg, info, fg, poisson, phi_prev, phi_curr, state, &
                                omega, a, nocc, hvol, dt, comm)
    implicit none
    integer, intent(in) :: step
    type(s_rgrid), intent(in) :: lg, mg
    type(s_parallel_info), intent(in) :: info
    type(s_reciprocal_grid), intent(in) :: fg
    type(s_poisson), intent(inout) :: poisson
    complex(8), intent(in) :: phi_prev(:,:,:,:,:)
    complex(8), intent(inout) :: phi_curr(:,:,:,:,:)
    type(ace_update_state), intent(inout) :: state
    real(8), intent(in) :: omega, a
    integer, intent(in) :: nocc
    real(8), intent(in) :: hvol, dt
    integer, intent(in), optional :: comm

    integer :: nocc_eff, gap, root_id, n_topm, n_sample, nmax_orig
    logical :: trigger_true, allowed_check
    integer, allocatable :: sample_idx(:), topm_global(:)
    real(8), allocatable :: sample_res(:)
    logical :: ok_build
    real(8) :: ex_sr_log

    if (.not. state%ace_enabled) return
    if (present(comm)) continue

    nocc_eff = min(nocc, maxval(info%io_e_all))
    if (nocc_eff <= 0) return
    root_id = info%id_rko

    state%omega = omega
    state%a = a
    call ensure_ace_storage(state, phi_curr, info, nocc_eff)

    if (.not. state%ace_built) then
      call try_build_ace_with_safety(lg,mg,info,fg,poisson,phi_curr,state,omega,a,nocc_eff,hvol,ok_build)
      state%ace_built = ok_build
      if (ok_build) then
        state%last_update_step = step
        state%need_second_stage_check = .false.
        state%force_update = .false.
        state%rebuild_done = .true.
        state%d_max = 0.d0
        state%i_max = 0
        state%R_max = 0.d0
      else
        state%rebuild_done = .false.
      end if
    else
      if (state%last_update_step <= -huge(1)/4) then
        gap = huge(1)/8
      else
        gap = step - state%last_update_step
      end if
      allowed_check = (gap >= state%n_min)
      state%force_update = (gap >= state%n_max)
      trigger_true = .false.
      state%rebuild_done = .false.
      state%R_max = 0.d0

      if (allowed_check) then
        call ace_first_stage_update_trigger(info, phi_prev, phi_curr, nocc_eff, trigger_true, &
             d_max=state%d_max, i_max=state%i_max, step=step, eps_d=state%eps_d, eps=state%eps, hvol=hvol, &
             print_summary=.false., print_topk=.false., d_values=state%d_values)
      else
        state%d_max = 0.d0
        state%i_max = 0
        state%d_values = 0.d0
      end if

      state%need_second_stage_check = state%force_update .or. trigger_true
      if (state%need_second_stage_check) then
        allocate(topm_global(max(1,state%topm)))
        call select_global_topm_by_candidates(info, state%d_values, nocc_eff, state%topm, topm_global, n_topm)

        allocate(sample_idx(max(1,state%ns_max)), sample_res(max(1,state%ns_max)))
        call build_sample_set(nocc_eff, nocc_eff, state%nv, topm_global(1:max(1,n_topm)), n_topm, &
                              state%nrand, state%seed, state%ns_max, sample_idx, n_sample)
        if (n_sample > 0) then
          call ace_stage2_residual_check(lg,mg,info,fg,poisson,sample_idx,n_sample,phi_curr,state%W,state%Minv, &
                                         omega,a,nocc_eff,state%R_max,residuals=sample_res, &
                                         comm_orb=info%icomm_o, comm_space=info%icomm_rko, hvol=hvol, &
                                         ksr_min=state%ksr_min)
        else
          state%R_max = 0.d0
        end if

        if (state%force_update .or. state%R_max > state%eps_R) then
          nmax_orig = state%n_max
          call try_build_ace_with_safety(lg,mg,info,fg,poisson,phi_curr,state,omega,a,nocc_eff,hvol,ok_build)
          if (ok_build) then
            state%last_update_step = step
            state%rebuild_done = .true.
          else
            state%rebuild_done = .false.
            state%n_max = max(state%n_min, min(state%n_max, 5))
            if (comm_is_root(root_id)) then
              write(*,'(A,I8)') 'ACE warning: rebuild failed; tighten n_max to ', state%n_max
            end if
          end if
          if (ok_build) state%n_max = max(state%n_max_base, nmax_orig)
          state%ace_built = ok_build
        end if
        deallocate(topm_global, sample_idx, sample_res)
      end if
    end if

    if (state%log_ex .and. state%ace_built) then
      if (mod(max(step-1,0), max(state%ex_stride,1)) == 0) then
        call compute_sr_exchange_energy_from_W(info, phi_curr, state%W, nocc_eff, hvol, ex_sr_log)
        if (comm_is_root(root_id)) then
          write(*,'(A,I8,A,1PE14.6)') 'ACE Ex_sr step=', step, ' value=', ex_sr_log
        end if
      end if
    end if

    if (comm_is_root(root_id)) then
      write(*,'(A,I8,A,1PE11.3,A,I8,A,L1,A,L1,A,1PE11.3,A,L1,A,1PE11.3,A,1PE11.3)') &
        'ACE step=', step, ' d_max=', state%d_max, ' i_max=', state%i_max, &
        ' need2=', state%need_second_stage_check, ' force=', state%force_update, &
        ' R_max=', state%R_max, ' rebuild=', state%rebuild_done, &
        ' min|diag(M)|=', state%stability_min_diag, ' cond_proxy=', state%stability_cond_proxy
    end if

    if (state%ace_built) then
      call ace_apply_rt_exchange(info, state, phi_curr, nocc_eff, hvol, dt)
    end if

  end subroutine ace_update_step_rt


  subroutine ace_apply_rt_exchange(info, state, phi_occ, nocc, hvol, dt)
    implicit none
    type(s_parallel_info), intent(in) :: info
    type(ace_update_state), intent(inout) :: state
    complex(8), intent(inout) :: phi_occ(:,:,:,:,:)
    integer, intent(in) :: nocc
    real(8), intent(in) :: hvol, dt

    integer :: i, io, nspin, nloc_owned, nocc_eff
    complex(8), allocatable :: kpsi(:,:,:,:)
    complex(8) :: zcoef

    if (.not. state%ace_enabled .or. .not. state%ace_built) return

    nocc_eff = min(nocc, maxval(info%io_e_all))
    nspin = size(phi_occ,4)
    nloc_owned = min(size(phi_occ,5), max(0, info%numo))
    if (nocc_eff <= 0 .or. nloc_owned <= 0) return

    allocate(kpsi(lbound(phi_occ,1):ubound(phi_occ,1), lbound(phi_occ,2):ubound(phi_occ,2), &
                  lbound(phi_occ,3):ubound(phi_occ,3), 1:nspin))
    zcoef = cmplx(0.d0, -dt, kind=8)

    do i = 1, nloc_owned
      io = info%io_s + i - 1
      if (io > nocc_eff) cycle
      call apply_exchange_ACE(info, state%W, state%Minv, phi_occ(:,:,:,1:nspin,i), kpsi, nocc=nocc_eff, &
                              comm_orb=info%icomm_o, comm_space=info%icomm_rko, hvol=hvol)
      phi_occ(:,:,:,1:nspin,i) = phi_occ(:,:,:,1:nspin,i) + zcoef * kpsi(:,:,:,1:nspin)
    end do

    deallocate(kpsi)
  end subroutine ace_apply_rt_exchange


  subroutine ensure_ace_storage(state, phi_occ, info, nocc_eff)
    implicit none
    type(ace_update_state), intent(inout) :: state
    complex(8), intent(in) :: phi_occ(:,:,:,:,:)
    type(s_parallel_info), intent(in) :: info
    integer, intent(in) :: nocc_eff

    integer :: nspin, nloc

    nspin = size(phi_occ,4)
    nloc = size(phi_occ,5)

    if (.not. allocated(state%W)) then
      allocate(state%W(lbound(phi_occ,1):ubound(phi_occ,1), lbound(phi_occ,2):ubound(phi_occ,2), &
                       lbound(phi_occ,3):ubound(phi_occ,3), 1:nspin, 1:nloc))
      state%W = (0.d0, 0.d0)
    else if (size(state%W,1) /= size(phi_occ,1) .or. size(state%W,2) /= size(phi_occ,2) .or. &
             size(state%W,3) /= size(phi_occ,3) .or. size(state%W,4) /= nspin .or. size(state%W,5) /= nloc) then
      deallocate(state%W)
      allocate(state%W(lbound(phi_occ,1):ubound(phi_occ,1), lbound(phi_occ,2):ubound(phi_occ,2), &
                       lbound(phi_occ,3):ubound(phi_occ,3), 1:nspin, 1:nloc))
      state%W = (0.d0, 0.d0)
      state%ace_built = .false.
    end if

    if (.not. allocated(state%Minv)) then
      allocate(state%Minv(nocc_eff,nocc_eff))
      state%Minv = (0.d0, 0.d0)
    else if (size(state%Minv,1) /= nocc_eff .or. size(state%Minv,2) /= nocc_eff) then
      deallocate(state%Minv)
      allocate(state%Minv(nocc_eff,nocc_eff))
      state%Minv = (0.d0, 0.d0)
      state%ace_built = .false.
    end if

    if (.not. allocated(state%d_values)) then
      allocate(state%d_values(nocc_eff))
      state%d_values = 0.d0
    else if (size(state%d_values) /= nocc_eff) then
      deallocate(state%d_values)
      allocate(state%d_values(nocc_eff))
      state%d_values = 0.d0
    end if

    if (nocc_eff > maxval(info%io_e_all)) then
      stop "ensure_ace_storage: inconsistent nocc_eff"
    end if
  end subroutine ensure_ace_storage


  subroutine try_build_ace_with_safety(lg,mg,info,fg,poisson,phi_curr,state,omega,a,nocc_eff,hvol,ok_build)
    implicit none
    type(s_rgrid), intent(in) :: lg, mg
    type(s_parallel_info), intent(in) :: info
    type(s_reciprocal_grid), intent(in) :: fg
    type(s_poisson), intent(inout) :: poisson
    complex(8), intent(in) :: phi_curr(:,:,:,:,:)
    type(ace_update_state), intent(inout) :: state
    real(8), intent(in) :: omega, a, hvol
    integer, intent(in) :: nocc_eff
    logical, intent(out) :: ok_build

    logical :: ok_local
    real(8) :: dshift_try

    call build_exchange_W_truth(lg,mg,info,fg,poisson,phi_curr,state%W,omega,a,nocc_eff, &
                                comm_orb=info%icomm_o, comm_space=info%icomm_rko, hvol=hvol)

    dshift_try = state%diag_shift
    call build_ACE_operator(info, phi_curr, state%W, state%Minv, nocc=nocc_eff, &
                            comm_orb=info%icomm_o, comm_space=info%icomm_rko, hvol=hvol, &
                            success=ok_local, min_diag_abs=state%stability_min_diag, &
                            cond_proxy=state%stability_cond_proxy, diag_shift=dshift_try)

    if (.not. ok_local .and. dshift_try <= 0.d0) then
      dshift_try = 1.0d-10
      call build_ACE_operator(info, phi_curr, state%W, state%Minv, nocc=nocc_eff, &
                              comm_orb=info%icomm_o, comm_space=info%icomm_rko, hvol=hvol, &
                              success=ok_local, min_diag_abs=state%stability_min_diag, &
                              cond_proxy=state%stability_cond_proxy, diag_shift=dshift_try)
    end if

    ok_build = ok_local
    if (.not. ok_build .and. comm_is_root(info%id_rko)) then
      write(*,'(A,1PE11.3,A,1PE11.3)') 'ACE warning: build_ACE_operator failed. min|diag(M)|=', &
           state%stability_min_diag, ' cond_proxy=', state%stability_cond_proxy
    end if

  end subroutine try_build_ace_with_safety


  subroutine compute_sr_exchange_energy_from_W(info, phi_occ, W, nocc_eff, hvol, ex_sr)
    implicit none
    type(s_parallel_info), intent(in) :: info
    complex(8), intent(in) :: phi_occ(:,:,:,:,:)
    complex(8), intent(in) :: W(:,:,:,:,:)
    integer, intent(in) :: nocc_eff
    real(8), intent(in) :: hvol
    real(8), intent(out) :: ex_sr

    integer :: i, io, ispin, nspin, nloc_owned
    real(8) :: ex_local

    nspin = size(phi_occ,4)
    nloc_owned = min(size(phi_occ,5), max(0, info%numo))
    ex_local = 0.d0
    do i = 1, nloc_owned
      io = info%io_s + i - 1
      if (io > nocc_eff) cycle
      do ispin = 1, nspin
        ex_local = ex_local - 0.5d0 * real(sum(conjg(phi_occ(:,:,:,ispin,i))*W(:,:,:,ispin,i)),8) * hvol
      end do
    end do
    call comm_summation(ex_local, ex_sr, info%icomm_rko)
  end subroutine compute_sr_exchange_energy_from_W


  subroutine select_global_topm_by_candidates(info, dvals, nocc_eff, topm, topm_indices, nsel)
    implicit none
    type(s_parallel_info), intent(in) :: info
    real(8), intent(in) :: dvals(:)
    integer, intent(in) :: nocc_eff, topm
    integer, intent(out) :: topm_indices(:)
    integer, intent(out) :: nsel

    integer :: ncand, i, j, k, idx, nslot, ofs, n_take
    integer, allocatable :: idx_send(:), idx_all(:)
    real(8), allocatable :: val_send(:), val_all(:), dtmp(:)

    nsel = 0
    topm_indices = 0
    if (topm <= 0 .or. nocc_eff <= 0) return

    ncand = max(1, 2*topm)
    nslot = max(1, info%isize_o*ncand)
    allocate(idx_send(nslot), idx_all(nslot), val_send(nslot), val_all(nslot), dtmp(max(1,nocc_eff)))
    idx_send = 0
    idx_all = 0
    val_send = 0.d0
    val_all = 0.d0
    dtmp = -1.d0

    do i = info%io_s, min(info%io_e, nocc_eff)
      dtmp(i) = dvals(i)
    end do

    ofs = info%id_o*ncand
    n_take = 0
    do k = 1, ncand
      idx = maxloc(dtmp, dim=1)
      if (idx < 1 .or. idx > nocc_eff) exit
      if (dtmp(idx) < 0.d0) exit
      n_take = n_take + 1
      idx_send(ofs+n_take) = idx
      val_send(ofs+n_take) = dtmp(idx)
      dtmp(idx) = -1.d0
    end do

    call comm_summation(idx_send, idx_all, nslot, info%icomm_o)
    call comm_summation(val_send, val_all, nslot, info%icomm_o)

    if (comm_is_root(info%id_o)) then
      dtmp = -1.d0
      do i = 1, nslot
        idx = idx_all(i)
        if (idx < 1 .or. idx > nocc_eff) cycle
        if (val_all(i) > dtmp(idx)) dtmp(idx) = val_all(i)
      end do
      nsel = 0
      do k = 1, min(topm,size(topm_indices))
        idx = maxloc(dtmp, dim=1)
        if (idx < 1 .or. idx > nocc_eff) exit
        if (dtmp(idx) < 0.d0) exit
        nsel = nsel + 1
        topm_indices(nsel) = idx
        dtmp(idx) = -1.d0
      end do
    end if

    call comm_bcast(nsel, info%icomm_o, 0)
    call comm_bcast(topm_indices, info%icomm_o, 0)

    deallocate(idx_send, idx_all, val_send, val_all, dtmp)
  end subroutine select_global_topm_by_candidates


  subroutine build_sample_set(nocc, vbm, nv, topm_indices, n_topm, nrand, seed, ns_max, sample_idx, n_sample)
    implicit none
    integer, intent(in) :: nocc, vbm, nv, n_topm, nrand, ns_max
    integer, intent(in) :: topm_indices(:)
    integer, intent(inout) :: seed
    integer, intent(out) :: sample_idx(:)
    integer, intent(out) :: n_sample

    logical, allocatable :: chosen(:)
    integer :: i, idx, ibeg, iend, n_slots, n_remaining, nrand_target, pick, tmp
    integer, allocatable :: remaining(:)
    integer(kind=8) :: x

    n_sample = 0
    sample_idx = 0
    if (nocc <= 0) return
    n_slots = min(ns_max,size(sample_idx))

    allocate(chosen(nocc))
    chosen = .false.

    ibeg = max(1, vbm-max(0,nv)+1)
    iend = min(nocc, max(1,vbm))
    do i = ibeg, iend
      if (n_sample >= n_slots) exit
      n_sample = n_sample + 1
      sample_idx(n_sample) = i
      chosen(i) = .true.
    end do

    do i = 1, n_topm
      idx = topm_indices(i)
      if (idx < 1 .or. idx > nocc) cycle
      if (chosen(idx)) cycle
      if (n_sample >= n_slots) exit
      n_sample = n_sample + 1
      sample_idx(n_sample) = idx
      chosen(idx) = .true.
    end do

    x = int(max(1,seed),kind=8)
    n_remaining = count(.not. chosen)
    nrand_target = min(nrand, max(0, n_slots - n_sample), n_remaining)
    if (nrand_target > 0) then
      allocate(remaining(n_remaining))
      n_remaining = 0
      do i = 1, nocc
        if (.not. chosen(i)) then
          n_remaining = n_remaining + 1
          remaining(n_remaining) = i
        end if
      end do

      do i = 1, nrand_target
        x = mod(1103515245_8*x + 12345_8, 2147483647_8)
        pick = i + int(mod(x, int(n_remaining - i + 1, kind=8)))
        tmp = remaining(i)
        remaining(i) = remaining(pick)
        remaining(pick) = tmp

        idx = remaining(i)
      n_sample = n_sample + 1
      sample_idx(n_sample) = idx
      chosen(idx) = .true.
      end do

      deallocate(remaining)
    end if
    seed = int(x)

    call sort_indices(sample_idx, n_sample)
    deallocate(chosen)
  end subroutine build_sample_set


  subroutine sort_indices(a, n)
    implicit none
    integer, intent(inout) :: a(:)
    integer, intent(in) :: n
    integer :: i, j, t
    do i = 1, max(0,n-1)
      do j = i+1, n
        if (a(j) < a(i)) then
          t = a(i); a(i) = a(j); a(j) = t
        end if
      end do
    end do
  end subroutine sort_indices

end module xc_ace_update_manager
