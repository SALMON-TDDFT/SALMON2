!
!  Grid-based HSE short-range convolution utilities
!

module xc_hse_grid_sr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical, save :: sr_kernel_selftest_done = .false.
  logical, save :: ex_validation_initialized = .false.
  real(8), save :: ex_prev = 0.d0
  real(8), save :: omega_prev = 0.d0
  logical, save :: ace_w_runtime_initialized = .false.
  integer, save :: ace_w_batch_env = 4
  logical, save :: ace_profile_env = .false.

contains

  subroutine apply_sr_convolution(lg,mg,info,fg,poisson,rho,u,omega)
    use structures, only: s_rgrid, s_parallel_info, s_reciprocal_grid, s_poisson, s_scalar
    use inputoutput, only: iperiodic, yn_ffte
#ifdef USE_FFTW
    use inputoutput, only: yn_fftw
#endif
    use poisson_periodic, only: poisson_ft_hse_sr, poisson_ffte_hse_sr
#ifdef USE_FFTW
    use poisson_periodic, only: poisson_fftw_hse_sr
#endif
    implicit none
    type(s_rgrid)          ,intent(in)    :: lg
    type(s_rgrid)          ,intent(in)    :: mg
    type(s_parallel_info)  ,intent(in)    :: info
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson)        ,intent(inout) :: poisson
    type(s_scalar)         ,intent(in)    :: rho
    type(s_scalar)         ,intent(inout) :: u
    real(8)                ,intent(in)    :: omega ! screening parameter [bohr^-1]
    character(16) :: env_selftest
    integer :: env_status

    if (iperiodic /= 3) then
      stop "apply_sr_convolution: currently supports periodic systems only"
    end if

    if (.not. sr_kernel_selftest_done) then
      env_selftest = ''
      call get_environment_variable('SALMON_HSE_SR_KERNEL_SELFTEST', env_selftest, status=env_status)
      if (env_status == 0) then
        select case(trim(adjustl(env_selftest)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          call selftest_hse_sr_kernel_limits()
          sr_kernel_selftest_done = .true.
        end select
      end if
    end if

#ifdef USE_FFTW
    if (yn_fftw == 'y') then
      call poisson_fftw_hse_sr(lg,mg,info,fg,rho,u,poisson,omega)
      return
    end if
#endif
    if (yn_ffte == 'y') then
      call poisson_ffte_hse_sr(lg,mg,info,fg,rho,u,poisson,omega)
    else
      call poisson_ft_hse_sr(lg,mg,info,fg,rho,u,poisson,omega)
    end if

  end subroutine apply_sr_convolution


  subroutine selftest_hse_sr_kernel_limits()
    implicit none
    real(8), parameter :: g2_test = 1.0d0
    real(8), parameter :: tol_small = 1.0d-12
    real(8), parameter :: tol_close = 1.0d-6
    real(8) :: f_large_omega, f_small_omega

    ! SR factor: 1 - exp(-G^2/(4*omega^2))
    ! G is in bohr^-1 (fg%vec_G comes from reciprocal primitive vectors in a.u.).
    f_large_omega = sr_factor_from_g2(g2_test, 1.0d12)   ! omega -> infinity => factor -> 0
    f_small_omega = sr_factor_from_g2(g2_test, 1.0d-8)   ! omega -> 0        => factor -> 1

    if (abs(f_large_omega) > tol_small) then
      stop "HSE-SR kernel selftest failed: omega->inf limit is not 0"
    end if
    if (abs(f_small_omega - 1.0d0) > tol_close) then
      stop "HSE-SR kernel selftest failed: omega->0 limit is not Coulomb"
    end if
  end subroutine selftest_hse_sr_kernel_limits


  pure real(8) function sr_factor_from_g2(g2, omega) result(f)
    implicit none
    real(8), intent(in) :: g2, omega
    if (omega <= 0.d0) then
      f = 0.d0
    else
      f = 1.0d0 - exp(-g2/(4.0d0*omega*omega))
    end if
  end function sr_factor_from_g2


  subroutine init_ace_w_runtime_env()
    implicit none
    character(32) :: env
    integer :: ist, ios, iv

    if (ace_w_runtime_initialized) return

    ace_w_batch_env = 4
    ace_profile_env = .false.

    env = ''
    call get_environment_variable('SALMON_ACE_W_BATCH', env, status=ist)
    if (ist == 0) then
      read(env,*,iostat=ios) iv
      if (ios == 0) ace_w_batch_env = max(1, iv)
    end if

    env = ''
    call get_environment_variable('SALMON_ACE_PROFILE', env, status=ist)
    if (ist == 0) then
      select case(trim(adjustl(env)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        ace_profile_env = .true.
      case default
        ace_profile_env = .false.
      end select
    end if

    ace_w_runtime_initialized = .true.
  end subroutine init_ace_w_runtime_env


  subroutine apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets,target_global_idx, &
                                 Kpsi_targets,omega,a,nocc,comm_orb,comm_space,hvol,ex_sr)
    use structures, only: s_rgrid, s_parallel_info, s_reciprocal_grid, s_poisson, s_scalar
    use communication, only: comm_exchange, comm_summation, comm_is_root
    implicit none
    type(s_rgrid)          ,intent(in)    :: lg
    type(s_rgrid)          ,intent(in)    :: mg
    type(s_parallel_info)  ,intent(in)    :: info
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson)        ,intent(inout) :: poisson
    complex(8)             ,intent(in)    :: psi_occ(:,:,:,:,:)
    complex(8)             ,intent(in)    :: psi_targets(:,:,:,:,:)
    integer                ,intent(in)    :: target_global_idx(:)
    complex(8)             ,intent(out)   :: Kpsi_targets(:,:,:,:,:)
    real(8)                ,intent(in)    :: omega, a
    integer                ,intent(in)    :: nocc
    integer                ,intent(in), optional :: comm_orb, comm_space
    real(8)                ,intent(in), optional :: hvol
    real(8)                ,intent(out), optional :: ex_sr

    type(s_scalar) :: rho_re, rho_im, u_re, u_im
    complex(8), allocatable :: psi_work(:,:,:,:,:), psi_recv(:,:,:,:,:), u_c(:,:,:)
    real(8) :: vol, ex_local, ex_global, ex_term
    integer :: nspin, nloc_occ, nloc_owned, n_target, nocc_eff
    integer :: ix1, ix2, iy1, iy2, iz1, iz2
    integer :: it, ispin, step, owner, nowner, jb, io_global, jstart, jcount, b, nb_max, target_io
    integer :: ix, iy, iz
    integer :: next_o, prev_o, icomm_o_eff, icomm_s_eff
    integer, parameter :: tag_ring = 27413
    real(8) :: t_pack_sendrecv, t_rho_build, t_sr_conv, t_accum, t0, t1
    integer :: batch_size

    call init_ace_w_runtime_env()

    nspin = size(psi_occ,4)
    nloc_occ = size(psi_occ,5)
    n_target = size(psi_targets,5)
    nloc_owned = min(nloc_occ, max(0, info%numo))
    nocc_eff = min(nocc, maxval(info%io_e_all))
    ix1 = lbound(psi_occ,1); ix2 = ubound(psi_occ,1)
    iy1 = lbound(psi_occ,2); iy2 = ubound(psi_occ,2)
    iz1 = lbound(psi_occ,3); iz2 = ubound(psi_occ,3)

    if (size(Kpsi_targets,1) /= size(psi_targets,1) .or. size(Kpsi_targets,2) /= size(psi_targets,2) .or. &
        size(Kpsi_targets,3) /= size(psi_targets,3) .or. size(Kpsi_targets,4) /= size(psi_targets,4) .or. &
        size(Kpsi_targets,5) /= n_target) then
      stop "apply_sr_fock_exact: Kpsi_targets shape mismatch"
    end if
    if (size(target_global_idx) /= n_target) then
      stop "apply_sr_fock_exact: target_global_idx size mismatch"
    end if
    do it = 1, n_target
      if (target_global_idx(it) < 1 .or. target_global_idx(it) > maxval(info%io_e_all)) then
        stop "apply_sr_fock_exact: target_global_idx out of range"
      end if
    end do

    vol = 1.d0
    if (present(hvol)) vol = hvol
    icomm_o_eff = info%icomm_o
    if (present(comm_orb)) icomm_o_eff = comm_orb
    icomm_s_eff = info%icomm_rko
    if (present(comm_space)) icomm_s_eff = comm_space

    batch_size = ace_w_batch_env
    batch_size = max(1, batch_size)
    batch_size = min(batch_size, max(1, maxval(info%numo_all)))
    batch_size = max(1, batch_size)
    nb_max = 1
    do owner = 0, info%isize_o - 1
      nb_max = max(nb_max, (max(0, info%numo_all(owner)) + batch_size - 1) / batch_size)
    end do

    allocate(psi_work(ix1:ix2,iy1:iy2,iz1:iz2,1:nspin,1:batch_size))
    allocate(psi_recv(ix1:ix2,iy1:iy2,iz1:iz2,1:nspin,1:batch_size))
    allocate(rho_re%f(ix1:ix2,iy1:iy2,iz1:iz2))
    allocate(rho_im%f(ix1:ix2,iy1:iy2,iz1:iz2))
    allocate(u_re%f(ix1:ix2,iy1:iy2,iz1:iz2))
    allocate(u_im%f(ix1:ix2,iy1:iy2,iz1:iz2))
    allocate(u_c(ix1:ix2,iy1:iy2,iz1:iz2))

    next_o = modulo(info%id_o + 1, info%isize_o)
    prev_o = modulo(info%id_o - 1 + info%isize_o, info%isize_o)
    Kpsi_targets = (0.d0, 0.d0)

    t_pack_sendrecv = 0.d0
    t_rho_build = 0.d0
    t_sr_conv = 0.d0
    t_accum = 0.d0

    do it = 1, n_target
      target_io = target_global_idx(it)
      if (target_io > nocc_eff) cycle
      do b = 1, nb_max
        psi_work = (0.d0, 0.d0)
        jstart = (b - 1) * batch_size + 1
        if (jstart <= nloc_owned) then
          jcount = min(batch_size, nloc_owned - jstart + 1)
          psi_work(:,:,:,:,1:jcount) = psi_occ(:,:,:,:,jstart:jstart+jcount-1)
        else
          jcount = 0
        end if

        do step = 0, info%isize_o - 1
          owner = modulo(info%id_o - step + info%isize_o, info%isize_o)
          nowner = max(0, info%numo_all(owner))
          jstart = (b - 1) * batch_size + 1
          if (jstart <= nowner) then
            jcount = min(batch_size, nowner - jstart + 1)
          else
            jcount = 0
          end if

          do jb = 1, jcount
            io_global = info%io_s_all(owner) + jstart + jb - 2
            if (io_global > nocc_eff) cycle

            call cpu_time(t0)
!$OMP parallel do collapse(3) private(ispin) schedule(static)
            do iz = iz1, iz2
              do iy = iy1, iy2
                do ix = ix1, ix2
                  rho_re%f(ix,iy,iz) = 0.d0
                  rho_im%f(ix,iy,iz) = 0.d0
                  do ispin = 1, nspin
                    rho_re%f(ix,iy,iz) = rho_re%f(ix,iy,iz) + &
                      real(conjg(psi_work(ix,iy,iz,ispin,jb)) * psi_targets(ix,iy,iz,ispin,it),8)
                    rho_im%f(ix,iy,iz) = rho_im%f(ix,iy,iz) + &
                      aimag(conjg(psi_work(ix,iy,iz,ispin,jb)) * psi_targets(ix,iy,iz,ispin,it))
                  end do
                end do
              end do
            end do
            call cpu_time(t1)
            t_rho_build = t_rho_build + (t1 - t0)

            call cpu_time(t0)
            call apply_sr_convolution(lg,mg,info,fg,poisson,rho_re,u_re,omega)
            call apply_sr_convolution(lg,mg,info,fg,poisson,rho_im,u_im,omega)
            call cpu_time(t1)
            t_sr_conv = t_sr_conv + (t1 - t0)

            u_c = cmplx(u_re%f, u_im%f, kind=8)
            call cpu_time(t0)
!$OMP parallel do collapse(4) schedule(static)
            do ispin = 1, nspin
              do iz = iz1, iz2
                do iy = iy1, iy2
                  do ix = ix1, ix2
                    Kpsi_targets(ix,iy,iz,ispin,it) = Kpsi_targets(ix,iy,iz,ispin,it) + &
                      psi_work(ix,iy,iz,ispin,jb) * u_c(ix,iy,iz)
                  end do
                end do
              end do
            end do
            call cpu_time(t1)
            t_accum = t_accum + (t1 - t0)
          end do

          if (step < info%isize_o - 1 .and. info%isize_o > 1) then
            call cpu_time(t0)
            call comm_exchange(psi_work, next_o, psi_recv, prev_o, tag_ring, icomm_o_eff)
            psi_work = psi_recv
            call cpu_time(t1)
            t_pack_sendrecv = t_pack_sendrecv + (t1 - t0)
          end if
        end do
      end do
    end do

    Kpsi_targets = (-a) * Kpsi_targets

    if (present(ex_sr)) then
      ex_local = 0.d0
      do it = 1, n_target
        target_io = target_global_idx(it)
        if (target_io > nocc_eff) cycle
        ex_term = 0.d0
        do ispin = 1, nspin
          ex_term = ex_term + real(sum(conjg(psi_targets(:,:,:,ispin,it)) * Kpsi_targets(:,:,:,ispin,it)),8) * vol
        end do
        ex_local = ex_local + 0.5d0 * ex_term
      end do
      call comm_summation(ex_local, ex_global, icomm_s_eff)
      ex_sr = ex_global
    end if

    if (ace_profile_env .and. comm_is_root(info%id_rko)) then
      write(*,'(A,4(1X,1PE12.4),A,I0)') 'ACE-W profile[s]: comm rho conv accum =', &
        t_pack_sendrecv, t_rho_build, t_sr_conv, t_accum, ' batch=', batch_size
    end if

    deallocate(psi_work, psi_recv, u_c)
    deallocate(rho_re%f, rho_im%f, u_re%f, u_im%f)
  end subroutine apply_sr_fock_exact


  subroutine compute_stage2_residual(lg,mg,info,fg,poisson,sample_indices,psi_prev,psi_curr,omega,a,residuals, &
                                     nocc,comm_orb,comm_space,hvol,ex_sr)
    use structures, only: s_rgrid, s_parallel_info, s_reciprocal_grid, s_poisson
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_rgrid)          ,intent(in)    :: lg
    type(s_rgrid)          ,intent(in)    :: mg
    type(s_parallel_info)  ,intent(in)    :: info
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson)        ,intent(inout) :: poisson
    integer                ,intent(in)    :: sample_indices(:)
    complex(8)             ,intent(in)    :: psi_prev(:,:,:,:,:,:,:)
    complex(8)             ,intent(in)    :: psi_curr(:,:,:,:,:,:,:)
    real(8)                ,intent(in)    :: omega, a
    real(8)                ,intent(out)   :: residuals(:)
    integer                ,intent(in), optional :: nocc, comm_orb, comm_space
    real(8)                ,intent(in), optional :: hvol
    real(8)                ,intent(out), optional :: ex_sr

    complex(8), allocatable :: psi_targets(:,:,:,:,:), Kpsi(:,:,:,:,:)
    real(8), allocatable :: res_local(:), res_global(:)
    integer, allocatable :: local_samples(:)
    real(8) :: vol, ex_local
    integer :: nsample, nloc_samp, i, j, io, iloc, ispin, nocc_eff, icomm_s_eff
    logical :: env_on
    character(16) :: env
    integer :: est, pos

    nsample = size(sample_indices)
    if (size(residuals) /= nsample) stop "compute_stage2_residual: residual size mismatch"
    residuals = 0.d0
    if (nsample <= 0) then
      if (present(ex_sr)) ex_sr = 0.d0
      return
    end if

    vol = 1.d0
    if (present(hvol)) vol = hvol
    icomm_s_eff = info%icomm_rko
    if (present(comm_space)) icomm_s_eff = comm_space

    nocc_eff = maxval(info%io_e_all)
    if (present(nocc)) nocc_eff = min(nocc, nocc_eff)

    env = ''
    call get_environment_variable('SALMON_HSE_SR_FOCK_TEST', env, status=est)
    env_on = .false.
    if (est == 0) then
      select case(trim(adjustl(env)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        env_on = .true.
      end select
    end if

    nloc_samp = 0
    do i = 1, nsample
      io = sample_indices(i)
      if (io >= info%io_s .and. io <= info%io_e) nloc_samp = nloc_samp + 1
    end do
    allocate(local_samples(max(1,nloc_samp)))
    pos = 0
    do i = 1, nsample
      io = sample_indices(i)
      if (io >= info%io_s .and. io <= info%io_e) then
        pos = pos + 1
        local_samples(pos) = io
      end if
    end do

    allocate(psi_targets(lbound(psi_curr,1):ubound(psi_curr,1), &
                         lbound(psi_curr,2):ubound(psi_curr,2), &
                         lbound(psi_curr,3):ubound(psi_curr,3), &
                         1:size(psi_curr,4), 1:max(1,nloc_samp)))
    allocate(Kpsi(lbound(psi_curr,1):ubound(psi_curr,1), &
                  lbound(psi_curr,2):ubound(psi_curr,2), &
                  lbound(psi_curr,3):ubound(psi_curr,3), &
                  1:size(psi_curr,4), 1:max(1,nloc_samp)))
    psi_targets = (0.d0, 0.d0)

    do j = 1, nloc_samp
      io = local_samples(j)
      iloc = io - info%io_s + 1
      psi_targets(:,:,:,1:size(psi_curr,4),j) = psi_curr(:,:,:,1:size(psi_curr,4),iloc,info%ik_s,info%im_s)
    end do

    if (present(comm_orb) .and. present(comm_space) .and. present(hvol)) then
      call apply_sr_fock_exact(lg,mg,info,fg,poisson, &
           psi_curr(:,:,:,1:size(psi_curr,4),info%io_s:info%io_e,info%ik_s,info%im_s), &
           psi_targets(:,:,:,1:size(psi_curr,4),1:max(1,nloc_samp)), local_samples(1:max(1,nloc_samp)), &
           Kpsi(:,:,:,1:size(psi_curr,4),1:max(1,nloc_samp)), omega, a, nocc_eff, &
           comm_orb=comm_orb, comm_space=comm_space, hvol=hvol, ex_sr=ex_local)
    else if (present(comm_orb) .and. present(comm_space)) then
      call apply_sr_fock_exact(lg,mg,info,fg,poisson, &
           psi_curr(:,:,:,1:size(psi_curr,4),info%io_s:info%io_e,info%ik_s,info%im_s), &
           psi_targets(:,:,:,1:size(psi_curr,4),1:max(1,nloc_samp)), local_samples(1:max(1,nloc_samp)), &
           Kpsi(:,:,:,1:size(psi_curr,4),1:max(1,nloc_samp)), omega, a, nocc_eff, &
           comm_orb=comm_orb, comm_space=comm_space, ex_sr=ex_local)
    else if (present(hvol)) then
      call apply_sr_fock_exact(lg,mg,info,fg,poisson, &
           psi_curr(:,:,:,1:size(psi_curr,4),info%io_s:info%io_e,info%ik_s,info%im_s), &
           psi_targets(:,:,:,1:size(psi_curr,4),1:max(1,nloc_samp)), local_samples(1:max(1,nloc_samp)), &
           Kpsi(:,:,:,1:size(psi_curr,4),1:max(1,nloc_samp)), omega, a, nocc_eff, hvol=hvol, ex_sr=ex_local)
    else
      call apply_sr_fock_exact(lg,mg,info,fg,poisson, &
           psi_curr(:,:,:,1:size(psi_curr,4),info%io_s:info%io_e,info%ik_s,info%im_s), &
           psi_targets(:,:,:,1:size(psi_curr,4),1:max(1,nloc_samp)), local_samples(1:max(1,nloc_samp)), &
           Kpsi(:,:,:,1:size(psi_curr,4),1:max(1,nloc_samp)), omega, a, nocc_eff, ex_sr=ex_local)
    end if

    allocate(res_local(nsample), res_global(nsample))
    res_local = 0.d0
    do i = 1, nsample
      io = sample_indices(i)
      if (io < info%io_s .or. io > info%io_e) cycle
      do j = 1, nloc_samp
        if (local_samples(j) /= io) cycle
        do ispin = 1, size(Kpsi,4)
          res_local(i) = res_local(i) + sum(abs(Kpsi(:,:,:,ispin,j))**2) * vol
        end do
      end do
      res_local(i) = sqrt(max(res_local(i),0.d0))
    end do

    call comm_summation(res_local, res_global, nsample, icomm_s_eff)
    residuals = res_global
    if (present(ex_sr)) ex_sr = ex_local

    if (env_on .and. comm_is_root(info%id_rko)) then
      if (present(ex_sr)) then
        write(*,'(A,1PE16.8)') 'HSE-SR exact Ex diagnostic = ', ex_sr
      end if
      do i = 1, nsample
        write(*,'(A,I6,A,1PE16.8)') 'HSE-SR stage2 baseline ||Kphi|| sample=', sample_indices(i), ' : ', residuals(i)
      end do
    end if

    deallocate(psi_targets, Kpsi, local_samples, res_local, res_global)
  end subroutine compute_stage2_residual


  subroutine build_exchange_W(lg,mg,info,fg,poisson,phi_occ,W,omega,a,nocc,ex_estimate,hvol)
    use structures, only: s_rgrid, s_parallel_info, s_reciprocal_grid, s_poisson, s_scalar
    use communication, only: comm_exchange, comm_summation, comm_is_root
    implicit none
    type(s_rgrid)          ,intent(in)    :: lg
    type(s_rgrid)          ,intent(in)    :: mg
    type(s_parallel_info)  ,intent(in)    :: info
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson)        ,intent(inout) :: poisson
    complex(8)             ,intent(in)    :: phi_occ(:,:,:,:,:)
    complex(8)             ,intent(out)   :: W(:,:,:,:,:)
    real(8)                ,intent(in)    :: omega
    real(8)                ,intent(in)    :: a
    integer                ,intent(in)    :: nocc
    real(8)                ,intent(out),optional :: ex_estimate
    real(8)                ,intent(in), optional :: hvol

    type(s_scalar) :: rho_re, rho_im, u_re, u_im
    complex(8), allocatable :: phi_work(:,:,:,:,:), phi_recv(:,:,:,:,:)
    complex(8), allocatable :: u_c(:,:,:)
    real(8) :: ex_local, ex_global, vol
    real(8) :: overlap_re, overlap_im, overlap_local
    integer :: nx_l, ny_l, nz_l, nspin, nloc, nloc_owned, numo_max
    integer :: ix1, ix2, iy1, iy2, iz1, iz2
    integer :: i, jloc, ispin, step, owner, nowner, nocc_eff
    integer :: io_global, next_o, prev_o
    integer, parameter :: tag_ring = 17421
    logical :: has_nan
    character(16) :: env_dbg
    integer :: env_stat

    nx_l = size(phi_occ,1)
    ny_l = size(phi_occ,2)
    nz_l = size(phi_occ,3)
    nspin = size(phi_occ,4)
    nloc = size(phi_occ,5)

    if (size(W,1) /= nx_l .or. size(W,2) /= ny_l .or. size(W,3) /= nz_l .or. &
        size(W,4) /= nspin .or. size(W,5) /= nloc) then
      stop "build_exchange_W: W shape mismatch"
    end if
    if (nocc < 1) then
      W = (0.d0, 0.d0)
      if (present(ex_estimate)) ex_estimate = 0.d0
      return
    end if
    nocc_eff = min(nocc, maxval(info%io_e_all))

    ix1 = lbound(phi_occ,1); ix2 = ubound(phi_occ,1)
    iy1 = lbound(phi_occ,2); iy2 = ubound(phi_occ,2)
    iz1 = lbound(phi_occ,3); iz2 = ubound(phi_occ,3)
    numo_max = max(1, info%numo_max)

    nloc_owned = min(nloc, max(0, info%numo))

    allocate(phi_work(ix1:ix2,iy1:iy2,iz1:iz2,1:nspin,1:numo_max))
    allocate(phi_recv(ix1:ix2,iy1:iy2,iz1:iz2,1:nspin,1:numo_max))
    phi_work = (0.d0, 0.d0)
    if (nloc_owned > 0) phi_work(:,:,:,:,1:nloc_owned) = phi_occ(:,:,:,:,1:nloc_owned)

    allocate(rho_re%f(ix1:ix2,iy1:iy2,iz1:iz2))
    allocate(rho_im%f(ix1:ix2,iy1:iy2,iz1:iz2))
    allocate(u_re%f(ix1:ix2,iy1:iy2,iz1:iz2))
    allocate(u_im%f(ix1:ix2,iy1:iy2,iz1:iz2))
    allocate(u_c(ix1:ix2,iy1:iy2,iz1:iz2))

    if (present(hvol)) then
      vol = hvol
    else
      vol = 1.d0
    end if

    next_o = modulo(info%id_o + 1, info%isize_o)
    prev_o = modulo(info%id_o - 1 + info%isize_o, info%isize_o)
    W = (0.d0, 0.d0)
    ex_local = 0.d0

    do i = 1, nloc_owned
      phi_work(:,:,:,:,1:numo_max) = (0.d0, 0.d0)
      if (nloc_owned > 0) phi_work(:,:,:,:,1:nloc_owned) = phi_occ(:,:,:,:,1:nloc_owned)

      do step = 0, info%isize_o - 1
        owner = modulo(info%id_o - step + info%isize_o, info%isize_o)
        nowner = max(0, info%numo_all(owner))
        do jloc = 1, nowner
          io_global = info%io_s_all(owner) + jloc - 1
          if (io_global > nocc_eff) cycle

          rho_re%f = 0.d0
          rho_im%f = 0.d0
          do ispin = 1, nspin
            rho_re%f = rho_re%f + real(conjg(phi_work(:,:,:,ispin,jloc)) * phi_occ(:,:,:,ispin,i),8)
            rho_im%f = rho_im%f + aimag(conjg(phi_work(:,:,:,ispin,jloc)) * phi_occ(:,:,:,ispin,i))
          end do

          call apply_sr_convolution(lg,mg,info,fg,poisson,rho_re,u_re,omega)
          call apply_sr_convolution(lg,mg,info,fg,poisson,rho_im,u_im,omega)
          u_c = cmplx(u_re%f, u_im%f, kind=8)

          do ispin = 1, nspin
            W(:,:,:,ispin,i) = W(:,:,:,ispin,i) + phi_work(:,:,:,ispin,jloc) * u_c
          end do
        end do

        if (step < info%isize_o - 1 .and. info%isize_o > 1) then
          call comm_exchange(phi_work, next_o, phi_recv, prev_o, tag_ring, info%icomm_o)
          phi_work = phi_recv
        end if
      end do

      overlap_re = 0.d0
      overlap_im = 0.d0
      do ispin = 1, nspin
        overlap_re = overlap_re + real(sum(conjg(phi_occ(:,:,:,ispin,i)) * W(:,:,:,ispin,i)),8)
        overlap_im = overlap_im + aimag(sum(conjg(phi_occ(:,:,:,ispin,i)) * W(:,:,:,ispin,i)))
      end do
      overlap_local = overlap_re
      ex_local = ex_local - 0.5d0 * a * overlap_local * vol
    end do

    W = (-a) * W

    call comm_summation(ex_local, ex_global, info%icomm_rko)
    if (present(ex_estimate)) ex_estimate = ex_global

    has_nan = .false.
    if (.not. all(ieee_is_finite(real(W,8)))) has_nan = .true.
    if (.not. all(ieee_is_finite(aimag(W)))) has_nan = .true.
    if (has_nan) stop "build_exchange_W: NaN/Inf detected in W"

    call validate_exchange_estimate(ex_global, omega, info%id_rko)

    env_dbg = ''
    call get_environment_variable('SALMON_HSE_SR_EX_DEBUG', env_dbg, status=env_stat)
    if (env_stat == 0) then
      select case(trim(adjustl(env_dbg)))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        if (comm_is_root(info%id_rko)) then
          write(*,'(A,1PE20.12)') 'HSE-SR exact exchange estimate Ex = ', ex_global
        end if
      end select
    end if

    deallocate(phi_work, phi_recv, u_c)
    deallocate(rho_re%f, rho_im%f, u_re%f, u_im%f)
  end subroutine build_exchange_W


  subroutine build_exchange_W_truth(lg,mg,info,fg,poisson,psi_occ,W,omega,a,nocc, &
                                    comm_orb,comm_space,hvol,ex_estimate)
    use structures, only: s_rgrid, s_parallel_info, s_reciprocal_grid, s_poisson
    use communication, only: comm_summation
    implicit none
    type(s_rgrid)          ,intent(in)    :: lg
    type(s_rgrid)          ,intent(in)    :: mg
    type(s_parallel_info)  ,intent(in)    :: info
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson)        ,intent(inout) :: poisson
    complex(8)             ,intent(in)    :: psi_occ(:,:,:,:,:)
    complex(8)             ,intent(out)   :: W(:,:,:,:,:)
    real(8)                ,intent(in)    :: omega, a
    integer                ,intent(in)    :: nocc
    integer                ,intent(in), optional :: comm_orb, comm_space
    real(8)                ,intent(in), optional :: hvol
    real(8)                ,intent(out), optional :: ex_estimate

    integer :: nspin, nloc, nloc_owned, nocc_eff, io, i, nsel
    integer, allocatable :: target_idx(:)
    complex(8), allocatable :: psi_targets(:,:,:,:,:), Kpsi_targets(:,:,:,:,:)
    integer :: icomm_o_eff
    real(8) :: ex_local, ex_global

    nspin = size(psi_occ,4)
    nloc = size(psi_occ,5)
    nloc_owned = min(nloc, max(0, info%numo))
    nocc_eff = min(nocc, maxval(info%io_e_all))
    icomm_o_eff = info%icomm_o
    if (present(comm_orb)) icomm_o_eff = comm_orb

    if (size(W,1) /= size(psi_occ,1) .or. size(W,2) /= size(psi_occ,2) .or. size(W,3) /= size(psi_occ,3) .or. &
        size(W,4) /= nspin .or. size(W,5) /= nloc) then
      stop "build_exchange_W_truth: W shape mismatch"
    end if

    W = (0.d0, 0.d0)
    if (nocc_eff <= 0 .or. nloc_owned <= 0) then
      if (present(ex_estimate)) ex_estimate = 0.d0
      return
    end if

    nsel = 0
    do i = 1, nloc_owned
      io = info%io_s + i - 1
      if (io <= nocc_eff) nsel = nsel + 1
    end do
    if (nsel <= 0) then
      if (present(ex_estimate)) ex_estimate = 0.d0
      return
    end if

    allocate(target_idx(nsel))
    allocate(psi_targets(lbound(psi_occ,1):ubound(psi_occ,1), &
                         lbound(psi_occ,2):ubound(psi_occ,2), &
                         lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin, 1:nsel))
    allocate(Kpsi_targets(lbound(psi_occ,1):ubound(psi_occ,1), &
                          lbound(psi_occ,2):ubound(psi_occ,2), &
                          lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin, 1:nsel))

    nsel = 0
    do i = 1, nloc_owned
      io = info%io_s + i - 1
      if (io > nocc_eff) cycle
      nsel = nsel + 1
      target_idx(nsel) = io
      psi_targets(:,:,:,1:nspin,nsel) = psi_occ(:,:,:,1:nspin,i)
    end do

    if (present(comm_orb) .and. present(comm_space) .and. present(hvol)) then
      call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets,target_idx,Kpsi_targets,omega,a,nocc_eff, &
                               comm_orb=comm_orb, comm_space=comm_space, hvol=hvol, ex_sr=ex_local)
    else if (present(comm_orb) .and. present(comm_space)) then
      call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets,target_idx,Kpsi_targets,omega,a,nocc_eff, &
                               comm_orb=comm_orb, comm_space=comm_space, ex_sr=ex_local)
    else if (present(hvol)) then
      call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets,target_idx,Kpsi_targets,omega,a,nocc_eff, &
                               hvol=hvol, ex_sr=ex_local)
    else
      call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets,target_idx,Kpsi_targets,omega,a,nocc_eff, &
                               ex_sr=ex_local)
    end if

    nsel = 0
    do i = 1, nloc_owned
      io = info%io_s + i - 1
      if (io > nocc_eff) cycle
      nsel = nsel + 1
      W(:,:,:,1:nspin,i) = Kpsi_targets(:,:,:,1:nspin,nsel)
    end do

    if (present(ex_estimate)) then
      call comm_summation(ex_local, ex_global, icomm_o_eff)
      ex_estimate = ex_global
    end if

    deallocate(target_idx, psi_targets, Kpsi_targets)
  end subroutine build_exchange_W_truth


  subroutine build_ACE_operator(info, psi_occ, W, Minv, nocc, comm_orb, comm_space, hvol, M_mat, &
                                success, min_diag_abs, cond_proxy, diag_shift)
    use structures, only: s_parallel_info
    use communication, only: comm_exchange, comm_summation
    implicit none
    type(s_parallel_info), intent(in) :: info
    complex(8), intent(in) :: psi_occ(:,:,:,:,:)
    complex(8), intent(in) :: W(:,:,:,:,:)
    complex(8), intent(out) :: Minv(:,:)
    integer, intent(in), optional :: nocc, comm_orb, comm_space
    real(8), intent(in), optional :: hvol
    complex(8), intent(out), optional :: M_mat(:,:)
    logical, intent(out), optional :: success
    real(8), intent(out), optional :: min_diag_abs, cond_proxy
    real(8), intent(in), optional :: diag_shift

    integer :: nocc_eff, nspin, nloc, nloc_owned, numo_max
    integer :: ix1, ix2, iy1, iy2, iz1, iz2
    integer :: i, jloc, io, jo, step, owner, nowner
    integer :: next_o, prev_o, icomm_o_eff, icomm_s_eff
    integer, parameter :: tag_ring = 38421
    real(8) :: vol
    complex(8) :: ztmp
    complex(8), allocatable :: W_work(:,:,:,:,:), W_recv(:,:,:,:,:)
    complex(8), allocatable :: M_local(:,:), M_global(:,:), M_work(:,:)
    integer :: info_inv
    real(8) :: dmin, dmax, dshift

    nocc_eff = maxval(info%io_e_all)
    if (present(nocc)) nocc_eff = min(nocc, nocc_eff)
    nspin = size(psi_occ,4)
    nloc = size(psi_occ,5)
    nloc_owned = min(nloc, max(0, info%numo))
    numo_max = max(1, info%numo_max)
    ix1 = lbound(psi_occ,1); ix2 = ubound(psi_occ,1)
    iy1 = lbound(psi_occ,2); iy2 = ubound(psi_occ,2)
    iz1 = lbound(psi_occ,3); iz2 = ubound(psi_occ,3)
    vol = 1.d0
    if (present(hvol)) vol = hvol
    icomm_o_eff = info%icomm_o
    if (present(comm_orb)) icomm_o_eff = comm_orb
    icomm_s_eff = info%icomm_rko
    if (present(comm_space)) icomm_s_eff = comm_space

    Minv = (0.d0, 0.d0)
    if (present(success)) success = .false.
    if (present(min_diag_abs)) min_diag_abs = 0.d0
    if (present(cond_proxy)) cond_proxy = huge(1.0d0)
    if (nocc_eff <= 0) return
    if (size(Minv,1) < nocc_eff .or. size(Minv,2) < nocc_eff) then
      stop "build_ACE_operator: Minv size mismatch"
    end if
    if (size(W,1) /= size(psi_occ,1) .or. size(W,2) /= size(psi_occ,2) .or. size(W,3) /= size(psi_occ,3) .or. &
        size(W,4) /= nspin .or. size(W,5) /= nloc) then
      stop "build_ACE_operator: psi_occ/W shape mismatch"
    end if
    if (present(M_mat)) then
      if (size(M_mat,1) < nocc_eff .or. size(M_mat,2) < nocc_eff) then
        stop "build_ACE_operator: M_mat size mismatch"
      end if
    end if

    allocate(W_work(ix1:ix2,iy1:iy2,iz1:iz2,1:nspin,1:numo_max))
    allocate(W_recv(ix1:ix2,iy1:iy2,iz1:iz2,1:nspin,1:numo_max))
    allocate(M_local(nocc_eff,nocc_eff), M_global(nocc_eff,nocc_eff), M_work(nocc_eff,nocc_eff))
    W_work = (0.d0, 0.d0)
    if (nloc_owned > 0) W_work(:,:,:,:,1:nloc_owned) = W(:,:,:,:,1:nloc_owned)
    M_local = (0.d0, 0.d0)
    M_global = (0.d0, 0.d0)
    M_work = (0.d0, 0.d0)

    next_o = modulo(info%id_o + 1, info%isize_o)
    prev_o = modulo(info%id_o - 1 + info%isize_o, info%isize_o)

    do step = 0, info%isize_o - 1
      owner = modulo(info%id_o - step + info%isize_o, info%isize_o)
      nowner = max(0, info%numo_all(owner))
      do jloc = 1, nowner
        jo = info%io_s_all(owner) + jloc - 1
        if (jo > nocc_eff) cycle
        do i = 1, nloc_owned
          io = info%io_s + i - 1
          if (io > nocc_eff) cycle
          ztmp = sum(conjg(psi_occ(:,:,:,1:nspin,i)) * W_work(:,:,:,1:nspin,jloc))
          M_local(io,jo) = M_local(io,jo) + vol * ztmp
        end do
      end do
      if (step < info%isize_o - 1 .and. info%isize_o > 1) then
        call comm_exchange(W_work, next_o, W_recv, prev_o, tag_ring, icomm_o_eff)
        W_work = W_recv
      end if
    end do

    ! comm_space is expected to span all ranks contributing to matrix elements.
    call comm_summation(M_local, M_global, nocc_eff*nocc_eff, icomm_s_eff)

    M_work = 0.5d0 * (M_global + transpose(conjg(M_global)))
    dshift = 0.d0
    if (present(diag_shift)) dshift = max(0.d0, diag_shift)
    if (dshift > 0.d0) then
      do i = 1, nocc_eff
        M_work(i,i) = M_work(i,i) - dshift
      end do
    end if
    dmin = huge(1.0d0)
    dmax = 0.d0
    do i = 1, nocc_eff
      dmin = min(dmin, abs(M_work(i,i)))
      dmax = max(dmax, abs(M_work(i,i)))
    end do
    if (present(min_diag_abs)) min_diag_abs = dmin
    if (present(cond_proxy)) cond_proxy = dmax / max(dmin, 1.0d-30)

    Minv(1:nocc_eff,1:nocc_eff) = M_work
    call invert_complex_matrix_checked(Minv(1:nocc_eff,1:nocc_eff), nocc_eff, "ACE M", info_inv)
    if (info_inv /= 0) then
      Minv(1:nocc_eff,1:nocc_eff) = (0.d0, 0.d0)
      if (present(success)) success = .false.
      if (present(M_mat)) M_mat(1:nocc_eff,1:nocc_eff) = M_work
      deallocate(W_work, W_recv, M_local, M_global, M_work)
      return
    end if
    if (present(success)) success = .true.
    if (present(M_mat)) M_mat(1:nocc_eff,1:nocc_eff) = M_work

    deallocate(W_work, W_recv, M_local, M_global, M_work)
  end subroutine build_ACE_operator


  subroutine apply_exchange_ACE(info, W, Minv, psi, Kpsi, nocc, comm_orb, comm_space, hvol)
    use structures, only: s_parallel_info
    use communication, only: comm_exchange, comm_summation
    implicit none
    type(s_parallel_info), intent(in) :: info
    complex(8), intent(in) :: W(:,:,:,:,:)
    complex(8), intent(in) :: Minv(:,:)
    complex(8), intent(in) :: psi(:,:,:,:)
    complex(8), intent(out) :: Kpsi(:,:,:,:)
    integer, intent(in), optional :: nocc, comm_orb, comm_space
    real(8), intent(in), optional :: hvol

    integer :: nocc_eff, nspin, nloc, nloc_owned, numo_max
    integer :: ix1, ix2, iy1, iy2, iz1, iz2
    integer :: step, owner, nowner, jloc, jo, ispin
    integer :: next_o, prev_o, icomm_o_eff, icomm_s_eff
    integer, parameter :: tag_ring = 38422
    real(8) :: vol
    complex(8) :: ztmp
    complex(8), allocatable :: W_work(:,:,:,:,:), W_recv(:,:,:,:,:)
    complex(8), allocatable :: b_local(:), b_global(:), c_coef(:)

    nocc_eff = maxval(info%io_e_all)
    if (present(nocc)) nocc_eff = min(nocc, nocc_eff)
    nspin = size(psi,4)
    nloc = size(W,5)
    nloc_owned = min(nloc, max(0, info%numo))
    numo_max = max(1, info%numo_max)
    ix1 = lbound(W,1); ix2 = ubound(W,1)
    iy1 = lbound(W,2); iy2 = ubound(W,2)
    iz1 = lbound(W,3); iz2 = ubound(W,3)
    vol = 1.d0
    if (present(hvol)) vol = hvol
    icomm_o_eff = info%icomm_o
    if (present(comm_orb)) icomm_o_eff = comm_orb
    icomm_s_eff = info%icomm_rko
    if (present(comm_space)) icomm_s_eff = comm_space

    Kpsi = (0.d0, 0.d0)
    if (nocc_eff <= 0) return
    if (size(Minv,1) < nocc_eff .or. size(Minv,2) < nocc_eff) then
      stop "apply_exchange_ACE: Minv size mismatch"
    end if
    if (size(psi,1) /= size(W,1) .or. size(psi,2) /= size(W,2) .or. size(psi,3) /= size(W,3) .or. &
        size(psi,4) /= nspin .or. size(Kpsi,1) /= size(W,1) .or. size(Kpsi,2) /= size(W,2) .or. &
        size(Kpsi,3) /= size(W,3) .or. size(Kpsi,4) /= nspin) then
      stop "apply_exchange_ACE: psi/W/Kpsi shape mismatch"
    end if

    allocate(W_work(ix1:ix2,iy1:iy2,iz1:iz2,1:nspin,1:numo_max))
    allocate(W_recv(ix1:ix2,iy1:iy2,iz1:iz2,1:nspin,1:numo_max))
    allocate(b_local(nocc_eff), b_global(nocc_eff), c_coef(nocc_eff))
    b_local = (0.d0, 0.d0)
    b_global = (0.d0, 0.d0)
    c_coef = (0.d0, 0.d0)
    W_work = (0.d0, 0.d0)
    if (nloc_owned > 0) W_work(:,:,:,:,1:nloc_owned) = W(:,:,:,:,1:nloc_owned)

    next_o = modulo(info%id_o + 1, info%isize_o)
    prev_o = modulo(info%id_o - 1 + info%isize_o, info%isize_o)

    do step = 0, info%isize_o - 1
      owner = modulo(info%id_o - step + info%isize_o, info%isize_o)
      nowner = max(0, info%numo_all(owner))
      do jloc = 1, nowner
        jo = info%io_s_all(owner) + jloc - 1
        if (jo > nocc_eff) cycle
        ztmp = sum(conjg(W_work(:,:,:,1:nspin,jloc)) * psi(:,:,:,1:nspin))
        b_local(jo) = b_local(jo) + vol * ztmp
      end do
      if (step < info%isize_o - 1 .and. info%isize_o > 1) then
        call comm_exchange(W_work, next_o, W_recv, prev_o, tag_ring, icomm_o_eff)
        W_work = W_recv
      end if
    end do

    call comm_summation(b_local, b_global, nocc_eff, icomm_s_eff)
    c_coef = matmul(Minv(1:nocc_eff,1:nocc_eff), b_global)

    W_work = (0.d0, 0.d0)
    if (nloc_owned > 0) W_work(:,:,:,:,1:nloc_owned) = W(:,:,:,:,1:nloc_owned)
    do step = 0, info%isize_o - 1
      owner = modulo(info%id_o - step + info%isize_o, info%isize_o)
      nowner = max(0, info%numo_all(owner))
      do jloc = 1, nowner
        jo = info%io_s_all(owner) + jloc - 1
        if (jo > nocc_eff) cycle
        do ispin = 1, nspin
          Kpsi(:,:,:,ispin) = Kpsi(:,:,:,ispin) + W_work(:,:,:,ispin,jloc) * c_coef(jo)
        end do
      end do
      if (step < info%isize_o - 1 .and. info%isize_o > 1) then
        call comm_exchange(W_work, next_o, W_recv, prev_o, tag_ring, icomm_o_eff)
        W_work = W_recv
      end if
    end do

    deallocate(W_work, W_recv, b_local, b_global, c_coef)
  end subroutine apply_exchange_ACE


  subroutine ace_build_test_diagnostics(lg,mg,info,fg,poisson,psi_occ,sample_indices,omega,a,nocc, &
                                        residuals,comm_orb,comm_space,hvol,ex_sr)
    use structures, only: s_rgrid, s_parallel_info, s_reciprocal_grid, s_poisson
    use communication, only: comm_is_root, comm_summation
    implicit none
    type(s_rgrid)          ,intent(in)    :: lg
    type(s_rgrid)          ,intent(in)    :: mg
    type(s_parallel_info)  ,intent(in)    :: info
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson)        ,intent(inout) :: poisson
    complex(8)             ,intent(in)    :: psi_occ(:,:,:,:,:)
    integer                ,intent(in)    :: sample_indices(:)
    real(8)                ,intent(in)    :: omega, a
    integer                ,intent(in)    :: nocc
    real(8)                ,intent(out)   :: residuals(:)
    integer                ,intent(in), optional :: comm_orb, comm_space
    real(8)                ,intent(in), optional :: hvol
    real(8)                ,intent(out), optional :: ex_sr

    integer :: nocc_eff, nspin, nloc, nloc_owned, nsamp, nloc_samp
    integer :: i, j, io, iloc, ispin, pos, icomm_s_eff
    real(8) :: vol, ex_local
    integer, allocatable :: local_samples(:)
    complex(8), allocatable :: W(:,:,:,:,:), Minv(:,:), psi_targets(:,:,:,:,:), Ktrue(:,:,:,:,:)
    complex(8), allocatable :: Kace(:,:,:,:)
    real(8), allocatable :: rloc(:), rglob(:)
    real(8) :: num_v, den_v

    nsamp = size(sample_indices)
    if (size(residuals) /= nsamp) stop "ace_build_test_diagnostics: residual size mismatch"
    residuals = 0.d0
    if (nsamp <= 0) then
      if (present(ex_sr)) ex_sr = 0.d0
      return
    end if

    nocc_eff = min(nocc, maxval(info%io_e_all))
    nspin = size(psi_occ,4)
    nloc = size(psi_occ,5)
    nloc_owned = min(nloc, max(0, info%numo))
    vol = 1.d0
    if (present(hvol)) vol = hvol
    icomm_s_eff = info%icomm_rko
    if (present(comm_space)) icomm_s_eff = comm_space

    allocate(W(lbound(psi_occ,1):ubound(psi_occ,1), lbound(psi_occ,2):ubound(psi_occ,2), &
               lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin, 1:nloc))
    allocate(Minv(nocc_eff,nocc_eff))

    if (present(comm_orb) .and. present(comm_space) .and. present(hvol)) then
      call build_exchange_W_truth(lg,mg,info,fg,poisson,psi_occ,W,omega,a,nocc_eff, &
                                  comm_orb=comm_orb, comm_space=comm_space, hvol=hvol, ex_estimate=ex_local)
      call build_ACE_operator(info, psi_occ, W, Minv, nocc=nocc_eff, comm_orb=comm_orb, comm_space=comm_space, hvol=hvol)
    else if (present(comm_orb) .and. present(comm_space)) then
      call build_exchange_W_truth(lg,mg,info,fg,poisson,psi_occ,W,omega,a,nocc_eff, &
                                  comm_orb=comm_orb, comm_space=comm_space, ex_estimate=ex_local)
      call build_ACE_operator(info, psi_occ, W, Minv, nocc=nocc_eff, comm_orb=comm_orb, comm_space=comm_space)
    else if (present(hvol)) then
      call build_exchange_W_truth(lg,mg,info,fg,poisson,psi_occ,W,omega,a,nocc_eff, hvol=hvol, ex_estimate=ex_local)
      call build_ACE_operator(info, psi_occ, W, Minv, nocc=nocc_eff, hvol=hvol)
    else
      call build_exchange_W_truth(lg,mg,info,fg,poisson,psi_occ,W,omega,a,nocc_eff, ex_estimate=ex_local)
      call build_ACE_operator(info, psi_occ, W, Minv, nocc=nocc_eff)
    end if
    if (present(ex_sr)) ex_sr = ex_local

    nloc_samp = 0
    do i = 1, nsamp
      io = sample_indices(i)
      if (io >= info%io_s .and. io <= info%io_e .and. io <= nocc_eff) nloc_samp = nloc_samp + 1
    end do
    allocate(local_samples(max(1,nloc_samp)))
    allocate(psi_targets(lbound(psi_occ,1):ubound(psi_occ,1), lbound(psi_occ,2):ubound(psi_occ,2), &
                         lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin, 1:max(1,nloc_samp)))
    allocate(Ktrue(lbound(psi_occ,1):ubound(psi_occ,1), lbound(psi_occ,2):ubound(psi_occ,2), &
                   lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin, 1:max(1,nloc_samp)))
    psi_targets = (0.d0, 0.d0)
    pos = 0
    do i = 1, nsamp
      io = sample_indices(i)
      if (io < info%io_s .or. io > info%io_e .or. io > nocc_eff) cycle
      pos = pos + 1
      local_samples(pos) = io
      iloc = io - info%io_s + 1
      psi_targets(:,:,:,1:nspin,pos) = psi_occ(:,:,:,1:nspin,iloc)
    end do

    if (nloc_samp > 0) then
      if (present(comm_orb) .and. present(comm_space) .and. present(hvol)) then
        call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets(:,:,:,:,1:nloc_samp), &
                                 local_samples(1:nloc_samp), Ktrue(:,:,:,:,1:nloc_samp), omega,a,nocc_eff, &
                                 comm_orb=comm_orb, comm_space=comm_space, hvol=hvol)
      else if (present(comm_orb) .and. present(comm_space)) then
        call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets(:,:,:,:,1:nloc_samp), &
                                 local_samples(1:nloc_samp), Ktrue(:,:,:,:,1:nloc_samp), omega,a,nocc_eff, &
                                 comm_orb=comm_orb, comm_space=comm_space)
      else if (present(hvol)) then
        call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets(:,:,:,:,1:nloc_samp), &
                                 local_samples(1:nloc_samp), Ktrue(:,:,:,:,1:nloc_samp), omega,a,nocc_eff, hvol=hvol)
      else
        call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets(:,:,:,:,1:nloc_samp), &
                                 local_samples(1:nloc_samp), Ktrue(:,:,:,:,1:nloc_samp), omega,a,nocc_eff)
      end if
    end if

    allocate(rloc(nsamp), rglob(nsamp))
    rloc = 0.d0
    if (nloc_samp > 0) allocate(Kace(lbound(psi_occ,1):ubound(psi_occ,1), lbound(psi_occ,2):ubound(psi_occ,2), &
                                     lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin))

    do i = 1, nsamp
      io = sample_indices(i)
      if (io < info%io_s .or. io > info%io_e .or. io > nocc_eff) cycle
      iloc = io - info%io_s + 1
      if (present(comm_orb) .and. present(comm_space) .and. present(hvol)) then
        call apply_exchange_ACE(info, W, Minv, psi_occ(:,:,:,1:nspin,iloc), Kace, nocc=nocc_eff, &
                                comm_orb=comm_orb, comm_space=comm_space, hvol=hvol)
      else if (present(comm_orb) .and. present(comm_space)) then
        call apply_exchange_ACE(info, W, Minv, psi_occ(:,:,:,1:nspin,iloc), Kace, nocc=nocc_eff, &
                                comm_orb=comm_orb, comm_space=comm_space)
      else if (present(hvol)) then
        call apply_exchange_ACE(info, W, Minv, psi_occ(:,:,:,1:nspin,iloc), Kace, nocc=nocc_eff, hvol=hvol)
      else
        call apply_exchange_ACE(info, W, Minv, psi_occ(:,:,:,1:nspin,iloc), Kace, nocc=nocc_eff)
      end if

      num_v = 0.d0
      den_v = 0.d0
      do j = 1, nloc_samp
        if (local_samples(j) /= io) cycle
        do ispin = 1, nspin
          num_v = num_v + sum(abs(Kace(:,:,:,ispin) - Ktrue(:,:,:,ispin,j))**2) * vol
          den_v = den_v + sum(abs(Ktrue(:,:,:,ispin,j))**2) * vol
        end do
      end do
      rloc(i) = sqrt(max(num_v,0.d0)) / max(sqrt(max(den_v,0.d0)), 1.0d-14)
    end do
    call comm_summation(rloc, rglob, nsamp, icomm_s_eff)
    residuals = rglob

    if (comm_is_root(info%id_rko)) then
      do i = 1, nsamp
        write(*,'(A,I6,A,1PE12.4)') 'ACE build-test residual sample=', sample_indices(i), ' : ', residuals(i)
      end do
    end if

    if (allocated(Kace)) deallocate(Kace)
    deallocate(W, Minv, local_samples, psi_targets, Ktrue, rloc, rglob)
  end subroutine ace_build_test_diagnostics


  subroutine select_ace_sample_orbitals(nocc, vbm_index, d_i, nv, m, nr, seed, sample_indices, n_sample, ns_max)
    implicit none
    integer, intent(in) :: nocc, vbm_index, nv, m, nr
    real(8), intent(in) :: d_i(:)
    integer, intent(inout) :: seed
    integer, intent(out) :: sample_indices(:)
    integer, intent(out) :: n_sample
    integer, intent(in), optional :: ns_max

    integer :: cap, i, j, k, ibeg, iend, idx, nadd
    integer(kind=8) :: x
    logical, allocatable :: chosen(:)
    real(8), allocatable :: dwork(:)

    n_sample = 0
    sample_indices = 0
    if (nocc <= 0) return
    if (size(d_i) < nocc) stop "select_ace_sample_orbitals: d_i size mismatch"
    cap = size(sample_indices)
    if (present(ns_max)) cap = min(cap, max(1, ns_max))

    allocate(chosen(nocc))
    allocate(dwork(nocc))
    chosen = .false.
    dwork = d_i(1:nocc)

    ibeg = max(1, vbm_index - max(0,nv) + 1)
    iend = min(nocc, max(1, vbm_index))
    do i = ibeg, iend
      if (n_sample >= cap) exit
      if (.not. chosen(i)) then
        n_sample = n_sample + 1
        sample_indices(n_sample) = i
        chosen(i) = .true.
      end if
    end do

    do k = 1, max(0,m)
      if (n_sample >= cap) exit
      idx = maxloc(dwork, dim=1)
      if (idx < 1 .or. idx > nocc) exit
      if (dwork(idx) < 0.d0) exit
      if (.not. chosen(idx)) then
        n_sample = n_sample + 1
        sample_indices(n_sample) = idx
        chosen(idx) = .true.
      end if
      dwork(idx) = -1.d0
    end do

    x = int(max(1,seed),kind=8)
    nadd = 0
    do while (nadd < max(0,nr) .and. n_sample < cap)
      x = mod(1103515245_8 * x + 12345_8, 2147483647_8)
      idx = 1 + int(mod(x, int(nocc,kind=8)))
      if (.not. chosen(idx)) then
        n_sample = n_sample + 1
        sample_indices(n_sample) = idx
        chosen(idx) = .true.
        nadd = nadd + 1
      end if
    end do
    seed = int(x)

    do i = 1, max(0,n_sample-1)
      do j = i+1, n_sample
        if (sample_indices(j) < sample_indices(i)) then
          k = sample_indices(i)
          sample_indices(i) = sample_indices(j)
          sample_indices(j) = k
        end if
      end do
    end do

    deallocate(chosen, dwork)
  end subroutine select_ace_sample_orbitals


  subroutine ace_stage2_residual_check(lg,mg,info,fg,poisson,sample_indices,n_sample,psi_occ, &
                                       W,Minv,omega,a,nocc,R_max,residuals,comm_orb,comm_space,hvol,ksr_min)
    use structures, only: s_rgrid, s_parallel_info, s_reciprocal_grid, s_poisson
    use communication, only: comm_summation
    implicit none
    type(s_rgrid)          ,intent(in)    :: lg
    type(s_rgrid)          ,intent(in)    :: mg
    type(s_parallel_info)  ,intent(in)    :: info
    type(s_reciprocal_grid),intent(in)    :: fg
    type(s_poisson)        ,intent(inout) :: poisson
    integer                ,intent(in)    :: sample_indices(:)
    integer                ,intent(in)    :: n_sample
    complex(8)             ,intent(in)    :: psi_occ(:,:,:,:,:)
    complex(8)             ,intent(in)    :: W(:,:,:,:,:)
    complex(8)             ,intent(in)    :: Minv(:,:)
    real(8)                ,intent(in)    :: omega, a
    integer                ,intent(in)    :: nocc
    real(8)                ,intent(out)   :: R_max
    real(8)                ,intent(out), optional :: residuals(:)
    integer                ,intent(in), optional :: comm_orb, comm_space
    real(8)                ,intent(in), optional :: hvol
    real(8)                ,intent(in), optional :: ksr_min

    integer :: nspin, nocc_eff, nloc_samp, i, j, io, iloc, ispin, pos, icomm_s_eff
    real(8) :: vol, num_v, den_v, ksr_min_use, den_norm
    integer, allocatable :: local_samples(:)
    complex(8), allocatable :: psi_targets(:,:,:,:,:), Ktrue(:,:,:,:,:), Kace(:,:,:,:)
    real(8), allocatable :: rloc(:), rglob(:)

    R_max = 0.d0
    nocc_eff = min(nocc, maxval(info%io_e_all))
    nspin = size(psi_occ,4)
    vol = 1.d0
    if (present(hvol)) vol = hvol
    ksr_min_use = 1.0d-12
    if (present(ksr_min)) ksr_min_use = max(0.d0, ksr_min)
    icomm_s_eff = info%icomm_rko
    if (present(comm_space)) icomm_s_eff = comm_space

    if (n_sample <= 0 .or. nocc_eff <= 0) then
      if (present(residuals)) residuals = 0.d0
      return
    end if
    if (present(residuals)) then
      if (size(residuals) < n_sample) stop "ace_stage2_residual_check: residual size mismatch"
      residuals(1:n_sample) = 0.d0
    end if

    nloc_samp = 0
    do i = 1, n_sample
      io = sample_indices(i)
      if (io >= info%io_s .and. io <= info%io_e .and. io <= nocc_eff) nloc_samp = nloc_samp + 1
    end do
    allocate(local_samples(max(1,nloc_samp)))
    allocate(psi_targets(lbound(psi_occ,1):ubound(psi_occ,1), lbound(psi_occ,2):ubound(psi_occ,2), &
                         lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin, 1:max(1,nloc_samp)))
    allocate(Ktrue(lbound(psi_occ,1):ubound(psi_occ,1), lbound(psi_occ,2):ubound(psi_occ,2), &
                   lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin, 1:max(1,nloc_samp)))
    allocate(Kace(lbound(psi_occ,1):ubound(psi_occ,1), lbound(psi_occ,2):ubound(psi_occ,2), &
                  lbound(psi_occ,3):ubound(psi_occ,3), 1:nspin))
    allocate(rloc(n_sample), rglob(n_sample))
    psi_targets = (0.d0, 0.d0)
    Ktrue = (0.d0, 0.d0)
    Kace = (0.d0, 0.d0)
    rloc = 0.d0
    rglob = 0.d0

    pos = 0
    do i = 1, n_sample
      io = sample_indices(i)
      if (io < info%io_s .or. io > info%io_e .or. io > nocc_eff) cycle
      pos = pos + 1
      local_samples(pos) = io
      iloc = io - info%io_s + 1
      psi_targets(:,:,:,1:nspin,pos) = psi_occ(:,:,:,1:nspin,iloc)
    end do

    if (nloc_samp > 0) then
      if (present(comm_orb) .and. present(comm_space) .and. present(hvol)) then
        call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets(:,:,:,:,1:nloc_samp), &
                                 local_samples(1:nloc_samp), Ktrue(:,:,:,:,1:nloc_samp), omega,a,nocc_eff, &
                                 comm_orb=comm_orb, comm_space=comm_space, hvol=hvol)
      else if (present(comm_orb) .and. present(comm_space)) then
        call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets(:,:,:,:,1:nloc_samp), &
                                 local_samples(1:nloc_samp), Ktrue(:,:,:,:,1:nloc_samp), omega,a,nocc_eff, &
                                 comm_orb=comm_orb, comm_space=comm_space)
      else if (present(hvol)) then
        call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets(:,:,:,:,1:nloc_samp), &
                                 local_samples(1:nloc_samp), Ktrue(:,:,:,:,1:nloc_samp), omega,a,nocc_eff, hvol=hvol)
      else
        call apply_sr_fock_exact(lg,mg,info,fg,poisson,psi_occ,psi_targets(:,:,:,:,1:nloc_samp), &
                                 local_samples(1:nloc_samp), Ktrue(:,:,:,:,1:nloc_samp), omega,a,nocc_eff)
      end if
    end if

    do i = 1, n_sample
      io = sample_indices(i)
      if (io < info%io_s .or. io > info%io_e .or. io > nocc_eff) cycle
      iloc = io - info%io_s + 1
      if (present(comm_orb) .and. present(comm_space) .and. present(hvol)) then
        call apply_exchange_ACE(info, W, Minv, psi_occ(:,:,:,1:nspin,iloc), Kace, nocc=nocc_eff, &
                                comm_orb=comm_orb, comm_space=comm_space, hvol=hvol)
      else if (present(comm_orb) .and. present(comm_space)) then
        call apply_exchange_ACE(info, W, Minv, psi_occ(:,:,:,1:nspin,iloc), Kace, nocc=nocc_eff, &
                                comm_orb=comm_orb, comm_space=comm_space)
      else if (present(hvol)) then
        call apply_exchange_ACE(info, W, Minv, psi_occ(:,:,:,1:nspin,iloc), Kace, nocc=nocc_eff, hvol=hvol)
      else
        call apply_exchange_ACE(info, W, Minv, psi_occ(:,:,:,1:nspin,iloc), Kace, nocc=nocc_eff)
      end if

      num_v = 0.d0
      den_v = 0.d0
      do j = 1, nloc_samp
        if (local_samples(j) /= io) cycle
        do ispin = 1, nspin
          num_v = num_v + sum(abs(Kace(:,:,:,ispin) - Ktrue(:,:,:,ispin,j))**2) * vol
          den_v = den_v + sum(abs(Ktrue(:,:,:,ispin,j))**2) * vol
        end do
      end do
      den_norm = sqrt(max(den_v,0.d0))
      if (den_norm < ksr_min_use) then
        rloc(i) = 0.d0
      else
        rloc(i) = sqrt(max(num_v,0.d0)) / max(den_norm, 1.0d-14)
      end if
    end do

    call comm_summation(rloc, rglob, n_sample, icomm_s_eff)
    R_max = maxval(rglob(1:n_sample))
    if (present(residuals)) residuals(1:n_sample) = rglob(1:n_sample)

    deallocate(local_samples, psi_targets, Ktrue, Kace, rloc, rglob)
  end subroutine ace_stage2_residual_check


  subroutine invert_complex_matrix_checked(a, n, label, info_out)
    implicit none
    complex(8), intent(inout) :: a(:,:)
    integer, intent(in) :: n
    character(*), intent(in) :: label
    integer, intent(out), optional :: info_out
    integer :: info_lapack, lwork
    integer, allocatable :: ipiv(:)
    complex(8), allocatable :: work(:)

    if (n <= 0) then
      if (present(info_out)) info_out = 0
      return
    end if
    if (size(a,1) < n .or. size(a,2) < n) then
      stop "invert_complex_matrix_checked: matrix shape mismatch"
    end if

    lwork = max(1, n*max(n,64))
    allocate(ipiv(n), work(lwork))
    call zgetrf(n, n, a, n, ipiv, info_lapack)
    if (info_lapack /= 0) then
      write(*,'(A,1X,A,1X,A,I0)') 'invert_complex_matrix_checked:', trim(label), 'zgetrf info=', info_lapack
      if (present(info_out)) info_out = info_lapack
      deallocate(ipiv, work)
      return
    end if
    call zgetri(n, a, n, ipiv, work, lwork, info_lapack)
    if (info_lapack /= 0) then
      write(*,'(A,1X,A,1X,A,I0)') 'invert_complex_matrix_checked:', trim(label), 'zgetri info=', info_lapack
      if (present(info_out)) info_out = info_lapack
      deallocate(ipiv, work)
      return
    end if
    if (present(info_out)) info_out = 0
    deallocate(ipiv, work)
  end subroutine invert_complex_matrix_checked


  subroutine validate_exchange_estimate(ex_value, omega, id_rank)
    use communication, only: comm_is_root
    implicit none
    real(8), intent(in) :: ex_value, omega
    integer, intent(in) :: id_rank
    character(16) :: env_validate
    integer :: env_status
    real(8) :: rel_jump

    env_validate = ''
    call get_environment_variable('SALMON_HSE_SR_EX_VALIDATE', env_validate, status=env_status)
    if (env_status /= 0) return
    select case(trim(adjustl(env_validate)))
    case('1','y','Y','yes','YES','true','TRUE','on','ON')
      continue
    case default
      return
    end select

    if (.not. ieee_is_finite(ex_value)) then
      stop "HSE-SR Ex validation failed: Ex is NaN/Inf"
    end if
    if (ex_value > 1.0d-8) then
      stop "HSE-SR Ex validation failed: Ex is not negative"
    end if

    if (ex_validation_initialized) then
      if (abs(omega - omega_prev) > 1.0d-14) then
        rel_jump = abs(ex_value - ex_prev) / max(abs(ex_prev), 1.0d-12)
        if (rel_jump > 10.0d0) then
          stop "HSE-SR Ex validation failed: Ex jump too large vs previous omega"
        end if
      end if
    end if

    ex_prev = ex_value
    omega_prev = omega
    ex_validation_initialized = .true.

    if (comm_is_root(id_rank)) then
      write(*,'(A,1PE20.12,A,1PE20.12)') 'HSE-SR Ex validate: Ex=', ex_value, ' omega=', omega
    end if
  end subroutine validate_exchange_estimate


  subroutine ace_first_stage_update_trigger(info, phi_prev, phi_curr, nocc, need_second_stage_check, &
                                            d_max, i_max, d_topk, i_topk, step, time_fs, eps_d, eps, hvol, &
                                            print_summary, print_topk, d_values)
    use structures, only: s_parallel_info
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_parallel_info), intent(in) :: info
    complex(8), intent(in) :: phi_prev(:,:,:,:,:)
    complex(8), intent(in) :: phi_curr(:,:,:,:,:)
    integer, intent(in) :: nocc
    logical, intent(out) :: need_second_stage_check
    real(8), intent(out), optional :: d_max
    integer, intent(out), optional :: i_max
    real(8), intent(out), optional :: d_topk(:)
    integer, intent(out), optional :: i_topk(:)
    integer, intent(in), optional :: step
    real(8), intent(in), optional :: time_fs
    real(8), intent(in), optional :: eps_d
    real(8), intent(in), optional :: eps
    real(8), intent(in), optional :: hvol
    logical, intent(in), optional :: print_summary, print_topk
    real(8), intent(out), optional :: d_values(:)

    integer :: nocc_eff, nloc_owned, nspin, i, io, ispin, nlog
    integer :: idx_max, k, idx_sel
    real(8) :: eps_d_use, eps_use, vol, theta
    real(8) :: norm_curr, diff_norm, dmax_local
    complex(8) :: cval, phase_factor
    complex(8), allocatable :: c_local(:), c_global(:)
    real(8), allocatable :: n2_local(:), n2_global(:), d2_local(:), d2_global(:), dvals(:), dtmp(:)
    real(8), allocatable :: topv(:)
    integer, allocatable :: topi(:)
    logical :: do_print_summary, do_print_topk

    if (size(phi_prev,1) /= size(phi_curr,1) .or. size(phi_prev,2) /= size(phi_curr,2) .or. &
        size(phi_prev,3) /= size(phi_curr,3) .or. size(phi_prev,4) /= size(phi_curr,4) .or. &
        size(phi_prev,5) /= size(phi_curr,5)) then
      stop "ace_first_stage_update_trigger: phi_prev/phi_curr shape mismatch"
    end if
    if (nocc < 1) then
      need_second_stage_check = .false.
      if (present(d_max)) d_max = 0.d0
      if (present(i_max)) i_max = 0
      if (present(d_topk)) d_topk = 0.d0
      if (present(i_topk)) i_topk = 0
      return
    end if

    eps_d_use = 1.0d-3
    if (present(eps_d)) eps_d_use = eps_d
    eps_use = 1.0d-14
    if (present(eps)) eps_use = eps
    vol = 1.d0
    if (present(hvol)) vol = hvol
    do_print_summary = .true.
    if (present(print_summary)) do_print_summary = print_summary
    do_print_topk = .true.
    if (present(print_topk)) do_print_topk = print_topk

    nocc_eff = min(nocc, maxval(info%io_e_all))
    nspin = size(phi_curr,4)
    nloc_owned = min(size(phi_curr,5), max(0, info%numo))
    nlog = min(5, nocc_eff)

    allocate(c_local(nocc_eff), c_global(nocc_eff))
    allocate(n2_local(nocc_eff), n2_global(nocc_eff))
    allocate(d2_local(nocc_eff), d2_global(nocc_eff))
    allocate(dvals(nocc_eff), dtmp(nocc_eff))
    allocate(topv(max(1,nlog)), topi(max(1,nlog)))

    c_local = (0.d0, 0.d0)
    n2_local = 0.d0
    d2_local = 0.d0

    do i = 1, nloc_owned
      io = info%io_s + i - 1
      if (io > nocc_eff) cycle
      do ispin = 1, nspin
        c_local(io) = c_local(io) + sum(conjg(phi_prev(:,:,:,ispin,i)) * phi_curr(:,:,:,ispin,i)) * vol
        n2_local(io) = n2_local(io) + sum(abs(phi_curr(:,:,:,ispin,i))**2) * vol
      end do
    end do

    call comm_summation(c_local, c_global, nocc_eff, info%icomm_rko)
    call comm_summation(n2_local, n2_global, nocc_eff, info%icomm_rko)

    do i = 1, nloc_owned
      io = info%io_s + i - 1
      if (io > nocc_eff) cycle
      cval = c_global(io)
      if (abs(cval) > eps_use) then
        theta = atan2(aimag(cval), real(cval,8))
      else
        theta = 0.d0
      end if
      phase_factor = cmplx(cos(theta), sin(theta), kind=8)
      do ispin = 1, nspin
        d2_local(io) = d2_local(io) + sum(abs(phi_curr(:,:,:,ispin,i) - phase_factor*phi_prev(:,:,:,ispin,i))**2) * vol
      end do
    end do

    call comm_summation(d2_local, d2_global, nocc_eff, info%icomm_rko)

    dvals = 0.d0
    do io = 1, nocc_eff
      norm_curr = sqrt(max(n2_global(io), 0.d0))
      diff_norm = sqrt(max(d2_global(io), 0.d0))
      dvals(io) = diff_norm / (norm_curr + eps_use)
    end do

    dmax_local = maxval(dvals)
    idx_max = maxloc(dvals, dim=1)
    need_second_stage_check = (dmax_local > eps_d_use)

    if (present(d_max)) d_max = dmax_local
    if (present(i_max)) i_max = idx_max

    dtmp = dvals
    topv = 0.d0
    topi = 0
    do k = 1, nlog
      idx_sel = maxloc(dtmp, dim=1)
      topv(k) = dtmp(idx_sel)
      topi(k) = idx_sel
      dtmp(idx_sel) = -1.d0
    end do

    if (present(d_topk)) then
      d_topk = 0.d0
      d_topk(1:min(size(d_topk),nlog)) = topv(1:min(size(d_topk),nlog))
    end if
    if (present(i_topk)) then
      i_topk = 0
      i_topk(1:min(size(i_topk),nlog)) = topi(1:min(size(i_topk),nlog))
    end if
    if (present(d_values)) then
      d_values = 0.d0
      d_values(1:min(size(d_values),nocc_eff)) = dvals(1:min(size(d_values),nocc_eff))
    end if

    if (do_print_summary .and. comm_is_root(info%id_rko)) then
      if (present(step) .and. present(time_fs)) then
        write(*,'(A,I10,A,F14.6,A,1PE12.4,A,I8,A,L1)') &
          'ACE1 trigger step=', step, ' time_fs=', time_fs, ' d_max=', dmax_local, &
          ' i_max=', idx_max, ' need_second_stage=', need_second_stage_check
      else if (present(step)) then
        write(*,'(A,I10,A,1PE12.4,A,I8,A,L1)') &
          'ACE1 trigger step=', step, ' d_max=', dmax_local, ' i_max=', idx_max, &
          ' need_second_stage=', need_second_stage_check
      else
        write(*,'(A,1PE12.4,A,I8,A,L1)') &
          'ACE1 trigger d_max=', dmax_local, ' i_max=', idx_max, &
          ' need_second_stage=', need_second_stage_check
      end if
      if (do_print_topk) then
        write(*,'(A)', advance='no') 'ACE1 top-k:'
        do k = 1, nlog
          write(*,'(A,I8,A,1PE12.4)', advance='no') ' (i=', topi(k), ',d=', topv(k), ')'
        end do
        write(*,*)
      end if
    end if

    deallocate(c_local, c_global, n2_local, n2_global, d2_local, d2_global, dvals, dtmp, topv, topi)
  end subroutine ace_first_stage_update_trigger


end module xc_hse_grid_sr
