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
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  implicit none

  private
  public :: init_plane_wave_basis
  public :: compute_fragment_pw_overlap
  public :: compute_mean_potential
  public :: compute_fragment_pw_hamiltonian
  public :: build_mixed_hamiltonian
  public :: diagonalize_mixed_basis
  public :: assemble_mixed_hamiltonian_dense
  public :: prepare_mixed_basis_startup

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
    integer :: n_pw_candidate, n_selected

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

    ! Sort k-points by energy (|k|)
    block
      real(8), allocatable :: k_norms(:)
      integer, allocatable :: sort_indices(:)
      integer :: ii, jj, itemp(3)
      real(8) :: rtemp

      allocate(k_norms(n_selected))
      allocate(sort_indices(n_selected))

      do ipw = 1, n_selected
        ikx = k_indices(1, ipw)
        iky = k_indices(2, ipw)
        ikz = k_indices(3, ipw)
        k_vec(1) = 2.0d0*pi/Lbox(1) * dble(ikx)
        k_vec(2) = 2.0d0*pi/Lbox(2) * dble(iky)
        k_vec(3) = 2.0d0*pi/Lbox(3) * dble(ikz)
        k_norms(ipw) = sqrt(sum(k_vec**2))
        sort_indices(ipw) = ipw
      end do

      do ii = 2, n_selected
        do jj = ii, 2, -1
          if (k_norms(jj) < k_norms(jj-1)) then
            rtemp = k_norms(jj)
            k_norms(jj) = k_norms(jj-1)
            k_norms(jj-1) = rtemp
            itemp = k_indices(:, jj)
            k_indices(:, jj) = k_indices(:, jj-1)
            k_indices(:, jj-1) = itemp
          else
            exit
          end if
        end do
      end do

      deallocate(k_norms, sort_indices)
    end block

    dg_frag%n_plane_waves = min(n_selected, n_plane_waves_dg)

    if (comm_is_root(info%id_rko)) then
      write(*,'(1x,a,i0)') "k-points within cutoff: ", n_selected
      write(*,'(1x,a,i0)') "Using plane waves: ", dg_frag%n_plane_waves
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
  ! Prepare mixed-basis operator blocks at startup without dense EVP.
  ! Fragment initial states are kept as-is and PW coefficients start from zero.
  !=======================================================================
  subroutine prepare_mixed_basis_startup(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    real(8), intent(in) :: Ac_tot(3)

    integer :: n_frag, n_pw
    complex(8), allocatable :: S_frag_pw(:,:,:), H_frag_pw(:,:,:), P_frag_pw(:,:,:,:)

    n_frag = dg_frag%n_mat_max
    n_pw = dg_frag%n_plane_waves
    if (.not. dg_frag%use_plane_wave_basis .or. n_pw <= 0) return

    allocate(S_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(H_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(P_frag_pw(3, n_frag, n_pw, dg_frag%nspin))

    call compute_fragment_pw_overlap(dg_frag, S_frag_pw)
    call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_frag_pw)
    call compute_fragment_pw_gradient_from_overlap(dg_frag, S_frag_pw, P_frag_pw)

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

    if (allocated(dg_frag%coef_pw)) dg_frag%coef_pw(:, :, :) = (0.0d0, 0.0d0)

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "  [init] Prepared mixed FP/PP operators without dense diagonalization"
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
      call compute_fragment_pw_overlap_legacy(dg_frag, S_complex)
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

    integer :: ifrag, i_local, ispin, io, ig, ipw
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
    logical, save :: warned_missing_buffer = .false.
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
        if (plane_wave_trace_enabled() .and. use_buffered_domain .and. .not. warned_missing_buffer .and. dg_frag%id == 0) then
          write(*,'(1x,a)') '[FP-DOMAIN] buffered requested but frag_buf bounds are unavailable; fallback to core domain'
          warned_missing_buffer = .true.
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
          if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
          S_complex(ig, 1:dg_frag%n_plane_waves, ispin) = frag_block(io, 1:dg_frag%n_plane_waves, ispin)
        end do
      end do
    end do

    s_fp_norm = sqrt(sum(abs(S_complex(1:dg_frag%n_mat_max, 1:dg_frag%n_plane_waves, 1:dg_frag%nspin))**2))
    if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
      write(*,'(1x,a,1x,1pe14.6)') '[FP-DOMAIN] ||S_fp||_F =', s_fp_norm
    end if

    if (fft_info_ready) call finalize_fragment_parallel(info_fft)

    deallocate(fft_in, fft_out, frag_block, k_fft_slot, k_fft_slot_neg)
  end subroutine compute_fragment_pw_overlap_fft_gspace

  subroutine compute_fragment_pw_overlap_legacy(dg_frag, S_complex)
    use structures
    use communication, only: comm_bcast, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(out) :: S_complex(:,:,:)  ! (n_mat_max, n_pw, nspin)

    integer :: ipw, ifrag, i_local, io, ig, ispin, ix, iy, iz
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
    logical, save :: warned_missing_buffer = .false.
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
        if (plane_wave_trace_enabled() .and. use_buffered_domain .and. .not. warned_missing_buffer .and. dg_frag%id == 0) then
          write(*,'(1x,a)') '[FP-DOMAIN] buffered requested but frag_buf bounds are unavailable; fallback to core domain'
          warned_missing_buffer = .true.
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
            if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
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

            S_complex(ig, ipw, ispin) = S_complex(ig, ipw, ispin) + overlap_local
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
            if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
            frag_block(io, 1:dg_frag%n_plane_waves, ispin) = S_complex(ig, 1:dg_frag%n_plane_waves, ispin)
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
          if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
          S_complex(ig, 1:dg_frag%n_plane_waves, ispin) = frag_block(io, 1:dg_frag%n_plane_waves, ispin)
        end do
      end do
    end do
    deallocate(frag_block, frag_block_sum)

    s_fp_norm = sqrt(sum(abs(S_complex(1:dg_frag%n_mat_max, 1:dg_frag%n_plane_waves, 1:dg_frag%nspin))**2))
    if (plane_wave_trace_enabled() .and. dg_frag%id == 0) then
      write(*,'(1x,a,1x,1pe14.6)') '[FP-DOMAIN] ||S_fp||_F =', s_fp_norm
    end if

  end subroutine compute_fragment_pw_overlap_legacy

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
!$omp parallel do collapse(3) private(ix,iy,iz) reduction(+:integral_local) schedule(static)
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

    integer :: ipw, ifrag, i_local, io, ig, ispin, ix, iy, iz
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
    logical, save :: warned_missing_buffer = .false.
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
      if (plane_wave_trace_enabled() .and. use_buffered_domain .and. .not. warned_missing_buffer .and. dg_frag%id == 0) then
        write(*,'(1x,a)') '[FP-DOMAIN] buffered requested but frag_buf bounds are unavailable; fallback to core domain'
        warned_missing_buffer = .true.
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
            if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
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

            H_complex(ig, ipw, ispin) = H_complex(ig, ipw, ispin) + hamiltonian_local
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
            if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
            frag_block(io, 1:dg_frag%n_plane_waves, ispin) = H_complex(ig, 1:dg_frag%n_plane_waves, ispin)
          end do
        end do
        ! H_fp follows the same ownership as S_fp: legacy reduces spatial
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
          if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
          H_complex(ig, 1:dg_frag%n_plane_waves, ispin) = frag_block(io, 1:dg_frag%n_plane_waves, ispin)
        end do
      end do
    end do
    deallocate(frag_block, frag_block_sum)

    h_fp_norm = sqrt(sum(abs(H_complex(1:dg_frag%n_mat_max, 1:dg_frag%n_plane_waves, 1:dg_frag%nspin))**2))
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

    n_pw = dg_frag%n_plane_waves
    H_pw(:, :, :) = (0.0d0, 0.0d0)
    if (n_pw <= 0) return

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

      do ipw = 1, n_pw
        H_pw(ipw, ipw, ispin) = H_pw(ipw, ipw, ispin) + 0.5d0 * sum(dg_frag%k_pw(:, ipw)**2)
      end do

      do jpw = 1, n_pw
        do ipw = jpw + 1, n_pw
          H_pw(ipw, jpw, ispin) = 0.5d0 * (H_pw(ipw, jpw, ispin) + conjg(H_pw(jpw, ipw, ispin)))
          H_pw(jpw, ipw, ispin) = conjg(H_pw(ipw, jpw, ispin))
        end do
      end do
    end do

    deallocate(pw_phase, local_block, global_block)
  end subroutine compute_pw_hamiltonian_local_potential

  !=======================================================================
  ! Build mixed Hamiltonian matrix with plane wave basis
  !=======================================================================
  subroutine build_mixed_hamiltonian(dg_frag, lg, Vh, Vxc, Vpsl, Ac_tot, &
                                      S_frag_pw, H_frag_pw)
    use structures
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: lg
    type(s_scalar),         intent(in)    :: Vh, Vxc(:), Vpsl
    real(8),                intent(in)    :: Ac_tot(3)
    complex(8),             intent(in)    :: S_frag_pw(:,:,:)
    complex(8),             intent(inout) :: H_frag_pw(:,:,:)

    integer :: ispin, i, ipw, n_frag, n_pw
    integer :: env_len, env_stat
    real(8) :: k_vec(3), kinetic_energy
    real(8), allocatable :: V_mean(:)
    complex(8), allocatable :: H_pw(:,:,:)
    complex(8), allocatable :: H_nl_fp(:,:), H_nl_pp(:,:), H_nl_diag(:,:)
    character(len=64) :: env_hpp_mode
    character(len=64) :: hpp_mode
    logical :: use_local_nondiag, include_nl_mixed

    if (.not. dg_frag%use_plane_wave_basis) return

    n_frag = dg_frag%n_mat_max
    n_pw = dg_frag%n_plane_waves

    hpp_mode = 'local'
    env_hpp_mode = ''
    call get_environment_variable('SALMON_DG_HPP_MODE', env_hpp_mode, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) hpp_mode = adjustl(env_hpp_mode(1:env_len))
    use_local_nondiag = .true.
    select case (hpp_mode(1:1))
    case ('d','D','l','L','0')
      use_local_nondiag = .false.
    case default
      use_local_nondiag = .true.
    end select

    allocate(V_mean(dg_frag%nspin))
    call compute_mean_potential(dg_frag, Vh, Vxc, Vpsl, V_mean)

    include_nl_mixed = dg_frag%has_nl_cache .and. allocated(dg_frag%H_nl_cache) .and. n_frag > 0 .and. n_pw > 0
    if (include_nl_mixed) then
      allocate(H_nl_fp(n_frag, n_pw), H_nl_pp(n_pw, n_pw), H_nl_diag(n_pw, dg_frag%nspin))
      H_nl_fp(:, :) = (0.0d0, 0.0d0)
      H_nl_pp(:, :) = (0.0d0, 0.0d0)
      H_nl_diag(:, :) = (0.0d0, 0.0d0)
    end if

    if (.not. allocated(dg_frag%H_mat_pw_diag)) then
      allocate(dg_frag%H_mat_pw_diag(n_pw, dg_frag%nspin))
    else if (size(dg_frag%H_mat_pw_diag,1) /= n_pw .or. size(dg_frag%H_mat_pw_diag,2) /= dg_frag%nspin) then
      deallocate(dg_frag%H_mat_pw_diag)
      allocate(dg_frag%H_mat_pw_diag(n_pw, dg_frag%nspin))
    end if

    if (use_local_nondiag .and. n_pw > 0) then
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
          H_nl_fp(:, :) = matmul(dg_frag%H_nl_cache(1:n_frag, 1:n_frag, ispin), S_frag_pw(1:n_frag, 1:n_pw, ispin))
          H_frag_pw(1:n_frag, 1:n_pw, ispin) = H_frag_pw(1:n_frag, 1:n_pw, ispin) + H_nl_fp(:, :)
          H_nl_pp(:, :) = matmul(conjg(transpose(S_frag_pw(1:n_frag, 1:n_pw, ispin))), H_nl_fp(:, :))
          H_pw(1:n_pw, 1:n_pw, ispin) = H_pw(1:n_pw, 1:n_pw, ispin) + H_nl_pp(:, :)
          H_nl_diag(1:n_pw, ispin) = [(H_nl_pp(i, i), i=1,n_pw)]
        end do
      end if
      dg_frag%H_mat_pw(:, :, :) = H_pw(:, :, :)
    else
      if (allocated(dg_frag%H_mat_pw)) deallocate(dg_frag%H_mat_pw)
      if (include_nl_mixed) then
        do ispin = 1, dg_frag%nspin
          H_nl_fp(:, :) = matmul(dg_frag%H_nl_cache(1:n_frag, 1:n_frag, ispin), S_frag_pw(1:n_frag, 1:n_pw, ispin))
          H_frag_pw(1:n_frag, 1:n_pw, ispin) = H_frag_pw(1:n_frag, 1:n_pw, ispin) + H_nl_fp(:, :)
          H_nl_pp(:, :) = matmul(conjg(transpose(S_frag_pw(1:n_frag, 1:n_pw, ispin))), H_nl_fp(:, :))
          H_nl_diag(1:n_pw, ispin) = [(H_nl_pp(i, i), i=1,n_pw)]
        end do
      end if
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
        if (use_local_nondiag) then
          dg_frag%H_mat_pw_diag(ipw, ispin) = dg_frag%H_mat_pw(ipw, ipw, ispin)
        else
          k_vec(1:3) = dg_frag%k_pw(1:3, ipw)
          kinetic_energy = 0.5d0 * sum(k_vec**2)
          i = n_frag + ipw
          dg_frag%H_mat_pw_diag(ipw, ispin) = cmplx(kinetic_energy + V_mean(ispin), 0.0d0, kind=8)
          if (include_nl_mixed) dg_frag%H_mat_pw_diag(ipw, ispin) = dg_frag%H_mat_pw_diag(ipw, ispin) + H_nl_diag(ipw, ispin)
        end if
      end do
    end do

    if (allocated(H_pw)) deallocate(H_pw)
    if (allocated(H_nl_fp)) deallocate(H_nl_fp)
    if (allocated(H_nl_pp)) deallocate(H_nl_pp)
    if (allocated(H_nl_diag)) deallocate(H_nl_diag)

    deallocate(V_mean)

    if (plane_wave_trace_enabled() .and. comm_is_root(dg_frag%id) .and. n_pw > 0) then
      if (use_local_nondiag) then
        write(*,'(1x,a)') '[HPP-MODE] local non-diagonal PW-PW potential enabled'
      else
        write(*,'(1x,a)') '[HPP-MODE] diagonal PW-PW approximation (legacy)'
      end if
    end if

  end subroutine build_mixed_hamiltonian

  !=======================================================================
  ! Diagonalize mixed basis using generalized eigenproblem:
  !   H c = eps S c
  !=======================================================================
  subroutine diagonalize_mixed_basis(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use communication, only: comm_is_root
    use rt_dg_fragment_ops, only: copy_matrix_blocks_to_complex_dense, sync_raw_coef_from_mixed, sync_mixed_coef_from_raw
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    real(8), intent(in) :: Ac_tot(3)

    integer :: ispin, n_total, n_frag, n_pw, lda, lwork, info, i, j, k
    integer :: n_floor, n_keep_s, n_drop_s, n_below_abs, n_below_tau
    integer :: env_len, env_stat
    real(8) :: sij, s_eff_min, s_eff_max, s_eff_cond
    real(8) :: tau_s, trunc_ratio, lam_keep_min, lam_max
    real(8) :: s_mix_min, s_mix_max, s_mix_cond
    real(8) :: pw_dim_eff, frag_dim_eff, pw_weight, pw_weight_drop_max
    real(8) :: eps_s_rel_cfg
    character(len=64) :: env_eps_rel
    character(len=64) :: env_init_mode
    character(len=64) :: init_mode_norm
    character(len=64) :: env_pw_weight_protect_n, env_apply_pw_weight_keep, env_append_pw_protected
    character(len=64) :: env_fragment_qr_keep
    character(len=64) :: env_pw_perp_project
    real(8), parameter :: pw_dom_thresh = 0.5d0, pw_nontriv_thresh = 0.1d0
    integer :: n_pw_dom_modes, n_pw_nontriv_modes, n_pw_drop_dom
    logical :: use_raw_prop_s
    complex(8), allocatable :: H_work(:,:), S_work(:,:), work_c(:)
    complex(8), allocatable :: X(:,:), H_ortho(:,:), tmp_mat(:,:), eigvec(:,:)
    real(8), allocatable :: eigenvalues_tmp(:), eval_s(:), rwork(:), work(:)
    real(8), allocatable :: S_eff(:,:), S_eff_work(:,:), eval_eff(:), work_eff(:)
    real(8), allocatable :: A_qr(:,:), tau_qr(:), work_qr(:)
    complex(8), allocatable :: S_frag_pw(:,:,:)  ! Complex overlap matrix
    complex(8), allocatable :: H_frag_pw(:,:,:)  ! Hamiltonian coupling matrix
    complex(8), allocatable :: P_frag_pw(:,:,:,:) ! Gradient coupling matrix
    integer, allocatable :: jpvt(:), keep_idx(:), keep_s_idx(:)
    integer :: m_qr, n_qr, lwork_qr, info_qr, n_keep_pw, n_keep_pw_base, ndiag
    integer :: n_keep_pw_premerge, n_keep_pw_before_final
    integer :: n_drop_final_rank, n_drop_real_rank, n_drop_metric_rank
    integer :: n_protect_pw, n_keep_protect
    integer :: n_pw_weight_protect, n_keep_pw_pwprotect
    real(8) :: diag_max, tau_rr
    real(8) :: k_protect_thr
    real(8), parameter :: eps_s_abs = 1.0d-10
    real(8), parameter :: eps_s_rel_default = 1.0d-8
    real(8), parameter :: tau_pw_rank_rel = 1.0d-6
    real(8) :: eps_s
    logical :: init_identity, init_truncated, init_occupied_projection, apply_pw_weight_keep, append_pw_protected
    logical :: use_fragment_overlap_qr_keep
    logical :: apply_pw_perp_projection
    integer :: nocc_init
    complex(8), allocatable :: S_init_metric(:,:), U_keep(:,:), UH_scaled(:,:)
    logical, allocatable :: protect_low_g(:), selected_keep(:)
    real(8), allocatable :: pw_col_proxy(:)
    integer, allocatable :: keep_idx_jpvt(:), keep_idx_pw_protect(:)
    external :: dsyev, dgeqp3, zheev

    n_frag = dg_frag%n_mat_max
    n_pw = dg_frag%n_plane_waves
    n_total = n_frag + n_pw
    init_identity = .false.
    init_truncated = .false.
    init_occupied_projection = .false.
    eps_s_rel_cfg = eps_s_rel_default
    env_eps_rel = ''
    call get_environment_variable('SALMON_DG_MIXED_S_EPS_REL', env_eps_rel, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_eps_rel(1:env_len), *, iostat=info) eps_s_rel_cfg
      if (info /= 0 .or. eps_s_rel_cfg <= 0.0d0) eps_s_rel_cfg = eps_s_rel_default
    end if
    env_init_mode = ''
    call get_environment_variable('SALMON_DG_MIXED_INIT_MODE', env_init_mode, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      init_mode_norm = trim(adjustl(env_init_mode(1:env_len)))
      if (init_mode_norm == 'identity') then
        init_identity = .true.
      else if (init_mode_norm == 'truncated_projection' .or. init_mode_norm == 'truncated_metric_projection') then
        init_truncated = .true.
      else if (init_mode_norm == 'occupied_projection' .or. init_mode_norm == 'conservative_occupied') then
        init_occupied_projection = .true.
      end if
    end if
    n_pw_weight_protect = 0
    env_pw_weight_protect_n = ''
    call get_environment_variable('SALMON_DG_PW_WEIGHT_PROTECT_N', env_pw_weight_protect_n, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      read(env_pw_weight_protect_n(1:env_len), *, iostat=info) n_pw_weight_protect
      if (info /= 0 .or. n_pw_weight_protect < 0) n_pw_weight_protect = 0
    end if
    apply_pw_weight_keep = .false.
    env_apply_pw_weight_keep = ''
    call get_environment_variable('SALMON_DG_APPLY_PW_WEIGHT_KEEP', env_apply_pw_weight_keep, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (env_apply_pw_weight_keep(1:1))
      case ('1','y','Y','t','T')
        apply_pw_weight_keep = .true.
      case default
        apply_pw_weight_keep = .false.
      end select
    end if
    append_pw_protected = .false.
    env_append_pw_protected = ''
    call get_environment_variable('SALMON_DG_PW_APPEND_PROTECTED', env_append_pw_protected, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (env_append_pw_protected(1:1))
      case ('1','y','Y','t','T')
        append_pw_protected = .true.
      case default
        append_pw_protected = .false.
      end select
    end if
    use_fragment_overlap_qr_keep = .false.
    env_fragment_qr_keep = ''
    call get_environment_variable('SALMON_DG_PW_FRAGMENT_QR_KEEP', env_fragment_qr_keep, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (env_fragment_qr_keep(1:1))
      case ('1','y','Y','t','T')
        use_fragment_overlap_qr_keep = .true.
      case default
        use_fragment_overlap_qr_keep = .false.
      end select
    end if
    apply_pw_perp_projection = .false.
    env_pw_perp_project = ''
    call get_environment_variable('SALMON_DG_PW_PERP_PROJECT', env_pw_perp_project, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      select case (env_pw_perp_project(1:1))
      case ('1','y','Y','t','T')
        apply_pw_perp_projection = .true.
      case default
        apply_pw_perp_projection = .false.
      end select
    end if

    if (comm_is_root(dg_frag%id)) then
      write(*,*)
      write(*,*) "=== Diagonalizing mixed basis (Fragment + PW) ==="
      write(*,'(1x,a,i0)') "Fragment basis size: ", n_frag
      write(*,'(1x,a,i0)') "Plane wave basis size: ", n_pw
      write(*,'(1x,a,i0)') "Total mixed basis size: ", n_total
      write(*,'(1x,a,1pe11.3)') "Mixed-S eps_s_rel: ", eps_s_rel_cfg
      if (init_identity) then
        write(*,'(1x,a)') "Mixed-init mode: identity"
      else if (init_truncated) then
        write(*,'(1x,a)') "Mixed-init mode: truncated_projection"
      else if (init_occupied_projection) then
        write(*,'(1x,a)') "Mixed-init mode: occupied_projection"
      else
        write(*,'(1x,a)') "Mixed-init mode: raw_projection"
      end if
      write(*,'(1x,a,i0,a,l1)') 'PW weight protection: n=', n_pw_weight_protect, &
           ' apply=', apply_pw_weight_keep
      write(*,'(1x,a,l1)') 'PW append protected modes: ', append_pw_protected
      write(*,'(1x,a,l1)') 'PW fragment-overlap QR keep: ', use_fragment_overlap_qr_keep
      write(*,'(1x,a,l1)') 'PW_perp projection: ', apply_pw_perp_projection
    end if

    allocate(S_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(H_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(P_frag_pw(3, n_frag, n_pw, dg_frag%nspin))
    call compute_fragment_pw_overlap(dg_frag, S_frag_pw)
    call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_frag_pw)
    call compute_fragment_pw_gradient_from_overlap(dg_frag, S_frag_pw, P_frag_pw)
    call select_protected_low_g_modes(dg_frag, protect_low_g, n_protect_pw, k_protect_thr)
    if (comm_is_root(dg_frag%id) .and. n_pw > 0) then
      write(*,'(1x,a,i0,a,i0,a,1pe11.3)') 'PW low-|G| protection: ', n_protect_pw, ' / ', n_pw, ' (k_thr=', k_protect_thr, ')'
    end if
    if (apply_pw_perp_projection) then
      call project_pw_to_fragment_orthogonal_complement(dg_frag, eps_s_abs, eps_s_rel_cfg, S_frag_pw, H_frag_pw, P_frag_pw, protect_low_g)
    else if (comm_is_root(dg_frag%id) .and. n_pw > 0) then
      ! Keep the matrix basis in the same representation used by density/current
      ! reconstruction.  The mixed metric below handles any FP/PW overlap.
      write(*,'(1x,a)') 'PW_perp projection skipped: using raw PW basis with full mixed metric'
    end if
    if (n_pw > 0 .and. n_frag > 0) then
      if (.not. use_fragment_overlap_qr_keep) then
        allocate(keep_idx(n_pw))
        do i = 1, n_pw
          keep_idx(i) = i
        end do
        n_keep_pw = n_pw
        n_keep_pw_premerge = n_keep_pw
        n_keep_pw_before_final = n_keep_pw
        n_drop_metric_rank = 0
        call filter_pw_keep_by_mixed_metric(dg_frag, S_frag_pw, keep_idx, n_keep_pw, n_keep_pw, &
             eps_s_abs, eps_s_rel_cfg, n_drop_metric_rank)
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,i0,a,i0)') 'PW metric keep stats: total_kept=', n_keep_pw, ' / ', n_keep_pw_before_final, &
               ' metric_dropped=', n_drop_metric_rank
        end if
        if (n_keep_pw < n_pw) then
          call compact_plane_wave_basis(dg_frag, S_frag_pw, H_frag_pw, P_frag_pw, keep_idx, n_keep_pw)
          n_pw = dg_frag%n_plane_waves
          n_total = n_frag + n_pw
          if (comm_is_root(dg_frag%id)) then
            write(*,'(1x,a,i0,a,i0)') "PW metric rank selection: kept ", n_keep_pw, " / ", n_keep_pw_before_final
          end if
        end if
        deallocate(keep_idx)
      else
        m_qr = 2 * n_frag * dg_frag%nspin
        n_qr = n_pw
        allocate(A_qr(m_qr, n_qr), jpvt(n_qr), tau_qr(min(m_qr, n_qr)))
        allocate(pw_col_proxy(n_qr))
        call compute_pw_weight_proxy_from_overlap(S_frag_pw, pw_col_proxy)
        tau_rr = tau_pw_rank_rel
        A_qr(:, :) = 0.0d0
        do ispin = 1, dg_frag%nspin
          A_qr((ispin-1)*2*n_frag + 1:(ispin-1)*2*n_frag + n_frag, :) = real(S_frag_pw(:, :, ispin), kind=8)
          A_qr((ispin-1)*2*n_frag + n_frag + 1:ispin*2*n_frag, :) = aimag(S_frag_pw(:, :, ispin))
        end do
        jpvt(:) = 0
        lwork_qr = -1
        allocate(work_qr(1))
        call dgeqp3(m_qr, n_qr, A_qr, m_qr, jpvt, tau_qr, work_qr, lwork_qr, info_qr)
        lwork_qr = max(1, int(work_qr(1)))
        deallocate(work_qr)
        allocate(work_qr(lwork_qr))
        call dgeqp3(m_qr, n_qr, A_qr, m_qr, jpvt, tau_qr, work_qr, lwork_qr, info_qr)
        n_keep_pw = n_pw
        if (info_qr == 0) then
          ! Legacy experimental path: keep PW columns by their fragment-overlap
          ! rank. After PW_perp projection this can remove valid complement
          ! directions, so it is opt-in only.
          ndiag = min(m_qr, n_qr)
          diag_max = 0.0d0
          do i = 1, ndiag
            diag_max = max(diag_max, abs(A_qr(i, i)))
          end do
          tau_rr = tau_pw_rank_rel * max(diag_max, 1.0d0)
          n_keep_pw = 0
          do i = 1, ndiag
            if (abs(A_qr(i, i)) >= tau_rr) n_keep_pw = i
          end do
          if (n_keep_pw <= 0 .and. n_pw > 0) n_keep_pw = 1
        else if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0)') "WARN: dgeqp3 failed in PW rank selection, info=", info_qr
        end if

        n_keep_pw_base = n_keep_pw
        if (n_pw_weight_protect <= 0) n_pw_weight_protect = n_protect_pw
        ! Protection is applied after QR so physically important low-|G| or
        ! high-overlap PW modes are not discarded only because they are nearly
        ! dependent on the fragment subspace.
        call build_keep_sets_from_jpvt(n_qr, jpvt, n_keep_pw_base, pw_col_proxy, n_pw_weight_protect, &
             keep_idx_jpvt, keep_idx_pw_protect, n_keep_pw_pwprotect)

        if (allocated(keep_idx)) deallocate(keep_idx)
        allocate(keep_idx(n_qr))
        if (apply_pw_weight_keep .and. n_keep_pw_pwprotect > 0) then
          n_keep_pw = n_keep_pw_pwprotect
          keep_idx(1:n_keep_pw) = keep_idx_pw_protect(1:n_keep_pw)
        else
          n_keep_pw = n_keep_pw_base
          keep_idx(1:n_keep_pw) = keep_idx_jpvt(1:n_keep_pw)
        end if

        n_keep_pw_premerge = n_keep_pw
        if (n_protect_pw > 0 .and. append_pw_protected) then
          ! Merge explicit low-|G| protection with the selected QR/proxy keep set
          ! while preserving uniqueness and the compacted PW ordering.
          allocate(selected_keep(n_pw))
          selected_keep(:) = .false.
          n_keep_pw = 0
          do i = 1, n_keep_pw_premerge
            if (keep_idx(i) < 1 .or. keep_idx(i) > n_pw) cycle
            if (selected_keep(keep_idx(i))) cycle
            n_keep_pw = n_keep_pw + 1
            keep_idx(n_keep_pw) = keep_idx(i)
            selected_keep(keep_idx(i)) = .true.
            if (n_keep_pw >= n_qr) exit
          end do
          do i = 1, n_pw
            if (.not. protect_low_g(i)) cycle
            if (.not. selected_keep(i)) then
              n_keep_pw = n_keep_pw + 1
              keep_idx(n_keep_pw) = i
              selected_keep(i) = .true.
            end if
          end do
          deallocate(selected_keep)
        end if

        n_keep_pw_before_final = n_keep_pw
        call filter_pw_keep_by_real_rank(S_frag_pw, keep_idx, n_keep_pw, tau_pw_rank_rel, n_drop_real_rank)
        call filter_pw_keep_by_mixed_metric(dg_frag, S_frag_pw, keep_idx, n_keep_pw, &
             min(n_keep_pw_premerge, n_keep_pw), eps_s_abs, eps_s_rel_cfg, n_drop_metric_rank)
        n_drop_final_rank = n_drop_real_rank + n_drop_metric_rank
        if (comm_is_root(dg_frag%id) .and. n_drop_final_rank > 0) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') 'PW final rank filter: kept=', n_keep_pw, &
               ' / ', n_keep_pw_before_final, ' dropped=', n_drop_final_rank, &
               ' real=', n_drop_real_rank, ' metric=', n_drop_metric_rank
          write(*,'(1x,a)', advance='no') 'PW final keep indices:'
          do i = 1, n_keep_pw
            write(*,'(1x,i0)', advance='no') keep_idx(i)
          end do
          write(*,*)
        end if

        n_keep_protect = 0
        if (n_keep_pw > 0 .and. n_protect_pw > 0) then
          if (allocated(keep_idx)) then
            n_keep_protect = count(protect_low_g(keep_idx(1:n_keep_pw)))
          else
            n_keep_protect = count(protect_low_g(jpvt(1:n_keep_pw)))
          end if
        end if
        if (comm_is_root(dg_frag%id) .and. n_protect_pw > 0) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') 'PW QR keep stats: total_kept=', n_keep_pw, ' / ', n_pw, &
               ' protected_kept=', n_keep_protect, ' / ', n_protect_pw
        end if

        if (n_keep_pw < n_pw) then
          call compact_plane_wave_basis(dg_frag, S_frag_pw, H_frag_pw, P_frag_pw, keep_idx, n_keep_pw)
          deallocate(keep_idx)
          n_pw = dg_frag%n_plane_waves
          n_total = n_frag + n_pw
          if (comm_is_root(dg_frag%id)) then
            write(*,'(1x,a,i0,a,i0,a,1pe11.3)') "PW rank selection: kept ", n_keep_pw, " / ", n_qr, " (tau=", tau_rr, ")"
          end if
        end if
        if (allocated(keep_idx)) deallocate(keep_idx)
        if (allocated(keep_idx_jpvt)) deallocate(keep_idx_jpvt)
        if (allocated(keep_idx_pw_protect)) deallocate(keep_idx_pw_protect)
        deallocate(A_qr, jpvt, tau_qr, work_qr, pw_col_proxy)
      end if
    end if
    if (allocated(protect_low_g)) deallocate(protect_low_g)

    if (.not. allocated(dg_frag%mixed_basis_dim)) then
      allocate(dg_frag%mixed_basis_dim(dg_frag%nspin))
    end if
    dg_frag%mixed_basis_dim(:) = 0
    if (allocated(dg_frag%mixed_transform)) deallocate(dg_frag%mixed_transform)
    allocate(dg_frag%mixed_transform(n_total, n_total, dg_frag%nspin))
    dg_frag%mixed_transform(:, :, :) = (0.0d0, 0.0d0)
    if (allocated(dg_frag%coef_mix)) deallocate(dg_frag%coef_mix)
    allocate(dg_frag%coef_mix(n_total, dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%coef_mix(:, :, :) = (0.0d0, 0.0d0)
    dg_frag%mixed_basis_ready = .false.

    if (.not. allocated(dg_frag%S_mat_frag_pw)) then
      allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    else if (size(dg_frag%S_mat_frag_pw,1) /= n_frag .or. size(dg_frag%S_mat_frag_pw,2) /= n_pw .or. &
             size(dg_frag%S_mat_frag_pw,3) /= dg_frag%nspin) then
      deallocate(dg_frag%S_mat_frag_pw)
      allocate(dg_frag%S_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    end if
    dg_frag%S_mat_frag_pw(1:n_frag,1:n_pw,1:dg_frag%nspin) = S_frag_pw(1:n_frag,1:n_pw,1:dg_frag%nspin)
    if (.not. allocated(dg_frag%H_mat_frag_pw)) then
      allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    else if (size(dg_frag%H_mat_frag_pw,1) /= n_frag .or. size(dg_frag%H_mat_frag_pw,2) /= n_pw .or. &
             size(dg_frag%H_mat_frag_pw,3) /= dg_frag%nspin) then
      deallocate(dg_frag%H_mat_frag_pw)
      allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    end if
    dg_frag%H_mat_frag_pw(1:n_frag,1:n_pw,1:dg_frag%nspin) = H_frag_pw(1:n_frag,1:n_pw,1:dg_frag%nspin)
    if (.not. allocated(dg_frag%P_mat_frag_pw)) then
      allocate(dg_frag%P_mat_frag_pw(3, n_frag, n_pw, dg_frag%nspin))
    else if (size(dg_frag%P_mat_frag_pw,1) /= 3 .or. size(dg_frag%P_mat_frag_pw,2) /= n_frag .or. &
             size(dg_frag%P_mat_frag_pw,3) /= n_pw .or. size(dg_frag%P_mat_frag_pw,4) /= dg_frag%nspin) then
      deallocate(dg_frag%P_mat_frag_pw)
      allocate(dg_frag%P_mat_frag_pw(3, n_frag, n_pw, dg_frag%nspin))
    end if
    dg_frag%P_mat_frag_pw(:,1:n_frag,1:n_pw,1:dg_frag%nspin) = P_frag_pw(:,1:n_frag,1:n_pw,1:dg_frag%nspin)
    call build_mixed_hamiltonian(dg_frag, dg_frag%lg, Vh, Vxc, Vpsl, Ac_tot, S_frag_pw, H_frag_pw)

    do ispin = 1, dg_frag%nspin
      lda = n_total

      allocate(H_work(n_total, n_total), S_work(n_total, n_total))
      allocate(eigenvalues_tmp(n_total))
      allocate(eval_s(n_total))
      allocate(S_eff(n_frag, n_frag), S_eff_work(n_frag, n_frag), eval_eff(n_frag))

      call assemble_mixed_hamiltonian_dense(dg_frag, ispin, H_frag_pw, H_work)
      S_work(:, :) = (0.0d0, 0.0d0)
      if (allocated(dg_frag%S_mat_c)) then
        S_work(1:n_frag, 1:n_frag) = dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin)
      else if (allocated(dg_frag%S_mat_blocks)) then
        call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%S_mat_blocks, ispin, S_work(1:n_frag, 1:n_frag))
      else if (allocated(dg_frag%S_mat)) then
        S_work(1:n_frag, 1:n_frag) = cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
      else
        do i = 1, n_frag
          S_work(i, i) = (1.0d0, 0.0d0)
        end do
      end if
      do j = 1, n_pw
        do i = 1, n_frag
          S_work(i, n_frag + j) = S_frag_pw(i, j, ispin)
          S_work(n_frag + j, i) = conjg(S_work(i, n_frag + j))
          H_work(i, n_frag + j) = H_frag_pw(i, j, ispin)
          H_work(n_frag + j, i) = conjg(H_work(i, n_frag + j))
        end do
      end do
      do i = n_frag + 1, n_total
        if (real(S_work(i, i), kind=8) < 1.0d0) S_work(i, i) = (1.0d0, 0.0d0)
      end do
      do i = 1, n_total
        if (real(S_work(i, i), kind=8) < eps_s_abs) S_work(i, i) = cmplx(eps_s_abs, 0.0d0, kind=8)
      end do

      ! Build propagation overlap from the same S definition used in mixed
      ! diagonalization (fragment block of S_mixed), then regularize.
      if (.not. allocated(dg_frag%S_mat_prop)) then
        allocate(dg_frag%S_mat_prop(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
        dg_frag%S_mat_prop(:, :, :) = 0.0d0
      end if
      if (.not. allocated(dg_frag%S_mat_prop_c)) then
        allocate(dg_frag%S_mat_prop_c(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
        dg_frag%S_mat_prop_c(:, :, :) = (0.0d0, 0.0d0)
      end if
      S_eff(:, :) = real(S_work(1:n_frag, 1:n_frag), kind=8)
      do j = 1, n_frag
        do i = j + 1, n_frag
          sij = 0.5d0 * (S_eff(i, j) + S_eff(j, i))
          S_eff(i, j) = sij
          S_eff(j, i) = sij
        end do
      end do
      S_eff_work(:, :) = S_eff(:, :)
      lwork = max(1, 3 * n_frag)
      allocate(work_eff(lwork))
      call dsyev('V', 'U', n_frag, S_eff_work, n_frag, eval_eff, work_eff, lwork, info)
      deallocate(work_eff)
      use_raw_prop_s = .false.
      s_eff_min = 0.0d0
      s_eff_cond = huge(1.0d0)
      if (info == 0) then
        s_eff_min = eval_eff(1)
        s_eff_max = eval_eff(n_frag)
        eps_s = max(eps_s_abs, eps_s_rel_cfg * max(s_eff_max, 1.0d0))
        if (s_eff_min <= 0.0d0) then
          s_eff_cond = huge(1.0d0)
        else
          s_eff_cond = s_eff_max / s_eff_min
        end if
        if (s_eff_min < 1.0d-6 .or. s_eff_cond > 1.0d8) then
          use_raw_prop_s = .true.
        end if
      end if
      if (info == 0 .and. .not. use_raw_prop_s) then
        n_floor = 0
        do i = 1, n_frag
          if (eval_eff(i) < eps_s) then
            eval_eff(i) = eps_s
            n_floor = n_floor + 1
          end if
        end do
        if (n_floor > 0 .and. comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,1pe11.3)') "Prop-S regularization: floored ", n_floor, " eigenvalues below ", eps_s
        end if
        S_eff(:, :) = 0.0d0
        do j = 1, n_frag
          do i = 1, n_frag
            S_eff(i, :) = S_eff(i, :) + S_eff_work(i, j) * eval_eff(j) * S_eff_work(:, j)
          end do
        end do
      else
        use_raw_prop_s = .true.
        if (info /= 0 .and. comm_is_root(dg_frag%id)) then
          write(*,*) "WARN: Prop-S diagonalization failed, keeping unregularized S_eff, info=", info
        end if
      end if
      if (use_raw_prop_s) then
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,1pe11.3,a,1pe11.3)') "Prop-S fallback to raw S_ff (min/cond): ", s_eff_min, " / ", s_eff_cond
        end if
        if (allocated(dg_frag%S_mat_c)) then
          dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin) = dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin)
          dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin) = real(dg_frag%S_mat_c(1:n_frag, 1:n_frag, ispin), kind=8)
        else if (allocated(dg_frag%S_mat_prop_blocks)) then
          call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%S_mat_prop_blocks, ispin, dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin))
          dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin) = real(dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin), kind=8)
        else if (allocated(dg_frag%S_mat_blocks)) then
          call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%S_mat_blocks, ispin, dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin))
          dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin) = real(dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin), kind=8)
        else if (allocated(dg_frag%S_mat)) then
          dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin) = dg_frag%S_mat(1:n_frag, 1:n_frag, ispin)
          dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin) = cmplx(dg_frag%S_mat(1:n_frag, 1:n_frag, ispin), 0.0d0, kind=8)
        else
          dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin) = S_eff(:, :)
          dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin) = cmplx(S_eff(:, :), 0.0d0, kind=8)
        end if
      else
        dg_frag%S_mat_prop(1:n_frag, 1:n_frag, ispin) = S_eff(:, :)
        dg_frag%S_mat_prop_c(1:n_frag, 1:n_frag, ispin) = cmplx(S_eff(:, :), 0.0d0, kind=8)
      end if
      if (n_frag < dg_frag%n_mat_max) then
        dg_frag%S_mat_prop(n_frag+1:dg_frag%n_mat_max, :, ispin) = 0.0d0
        dg_frag%S_mat_prop(:, n_frag+1:dg_frag%n_mat_max, ispin) = 0.0d0
        dg_frag%S_mat_prop_c(n_frag+1:dg_frag%n_mat_max, :, ispin) = (0.0d0, 0.0d0)
        dg_frag%S_mat_prop_c(:, n_frag+1:dg_frag%n_mat_max, ispin) = (0.0d0, 0.0d0)
      end if

      do j = 1, n_total
        S_work(j, j) = cmplx(real(S_work(j, j), kind=8), 0.0d0, kind=8)
        do i = j + 1, n_total
          S_work(i, j) = conjg(S_work(j, i))
        end do
      end do

      ! Regularize generalized EVP by Lowdin transform with eigencut.
      ! S = U diag(s) U^H, keep s >= tau_s, X = U_keep diag(1/sqrt(s_keep)).
      lwork = -1
      allocate(work_c(1), rwork(max(1, 3*n_total-2)))
      call ZHEEV('V', 'U', n_total, S_work, lda, eval_s, work_c, lwork, rwork, info)
      lwork = int(real(work_c(1), kind=8)) + 1
      deallocate(work_c)
      allocate(work_c(lwork))
      call ZHEEV('V', 'U', n_total, S_work, lda, eval_s, work_c, lwork, rwork, info)
      if (info /= 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,*) "WARN: S diagonalization failed in mixed basis, info=", info
          write(*,*) "WARN: Skipping mixed-basis refresh for this spin channel"
        end if
        deallocate(H_work, S_work, eigenvalues_tmp, eval_s, work_c, rwork)
        deallocate(S_eff, S_eff_work, eval_eff)
        cycle
      end if
      eps_s = max(eps_s_abs, eps_s_rel_cfg * max(eval_s(n_total), 1.0d0))
      s_mix_min = eval_s(1)
      s_mix_max = eval_s(n_total)
      if (s_mix_min > 0.0d0) then
        s_mix_cond = s_mix_max / s_mix_min
      else
        s_mix_cond = huge(1.0d0)
      end if

      tau_s = eps_s
      n_keep_s = 0
      do i = 1, n_total
        if (eval_s(i) >= tau_s) n_keep_s = n_keep_s + 1
      end do
      if (n_keep_s <= 0) n_keep_s = 1
      n_drop_s = n_total - n_keep_s
      n_below_abs = count(eval_s(1:n_total) < eps_s_abs)
      n_below_tau = count(eval_s(1:n_total) < tau_s)
      allocate(keep_s_idx(n_keep_s))
      k = 0
      do i = 1, n_total
        if (eval_s(i) >= tau_s) then
          k = k + 1
          if (k <= n_keep_s) keep_s_idx(k) = i
        end if
      end do
      if (k < n_keep_s) then
        do i = k + 1, n_keep_s
          keep_s_idx(i) = n_total - n_keep_s + i
        end do
      end if
      lam_keep_min = eval_s(keep_s_idx(1))
      lam_max = eval_s(n_total)
      trunc_ratio = real(n_drop_s, kind=8) / max(1.0d0, real(n_total, kind=8))
      pw_dim_eff = 0.0d0
      n_pw_dom_modes = 0
      n_pw_nontriv_modes = 0
      if (n_pw > 0) then
        do j = 1, n_keep_s
          i = keep_s_idx(j)
          pw_weight = sum(abs(S_work(n_frag+1:n_total, i))**2)
          pw_dim_eff = pw_dim_eff + pw_weight
          if (pw_weight >= pw_dom_thresh) n_pw_dom_modes = n_pw_dom_modes + 1
          if (pw_weight >= pw_nontriv_thresh) n_pw_nontriv_modes = n_pw_nontriv_modes + 1
        end do
      end if
      frag_dim_eff = real(n_keep_s, kind=8) - pw_dim_eff
      n_pw_drop_dom = 0
      pw_weight_drop_max = 0.0d0
      if (n_pw > 0 .and. n_drop_s > 0) then
        do i = 1, n_total
          if (eval_s(i) >= tau_s) cycle
          pw_weight = sum(abs(S_work(n_frag+1:n_total, i))**2)
          pw_weight_drop_max = max(pw_weight_drop_max, pw_weight)
          if (pw_weight >= pw_dom_thresh) n_pw_drop_dom = n_pw_drop_dom + 1
        end do
      end if
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,1pe11.3,a,1pe11.3,a,1pe11.3,a,i0,a,i0)') "Mixed-S spectrum: spin=", ispin, &
             " lam_min=", s_mix_min, " lam_max=", s_mix_max, " cond=", s_mix_cond, &
             " below_abs=", n_below_abs, " below_tau=", n_below_tau
        write(*,'(1x,a,i0,a,i0,a,1pe11.3,a,1pe11.3,a,i0,a,i0,a,1pe11.3)') "Mixed-S effective rank: ", n_keep_s, " / ", n_total, &
             " tau=", tau_s, " lam_min_keep=", lam_keep_min, " drop=", n_drop_s, " / ", n_total, " trunc_ratio=", trunc_ratio
        write(*,'(1x,a,1pe11.3,a,1pe11.3,a,i0,a,i0,a,i0,a,i0)') "Mixed-S composition: pw_eff=", pw_dim_eff, &
             " frag_eff=", frag_dim_eff, " pw_dom_modes=", n_pw_dom_modes, " / ", n_keep_s, &
             " pw_nontriv_modes=", n_pw_nontriv_modes, " / ", n_keep_s
        write(*,'(1x,a,i0,a,i0,a,1pe11.3)') "Mixed-S dropped PW-like: ", n_pw_drop_dom, " / ", n_drop_s, &
             " max_pw_weight_dropped=", pw_weight_drop_max
      end if

      allocate(X(n_total, n_keep_s), tmp_mat(n_keep_s, n_total), H_ortho(n_keep_s, n_keep_s), eigvec(n_total, n_keep_s))
      do j = 1, n_keep_s
        i = keep_s_idx(j)
        X(:, j) = S_work(:, i) / sqrt(eval_s(i))
      end do
      tmp_mat = matmul(conjg(transpose(X)), H_work)
      H_ortho = matmul(tmp_mat, X)
      do j = 1, n_keep_s
        H_ortho(j, j) = cmplx(real(H_ortho(j, j), kind=8), 0.0d0, kind=8)
        do i = j + 1, n_keep_s
          H_ortho(i, j) = conjg(H_ortho(j, i))
        end do
      end do

      deallocate(work_c)
      lwork = -1
      allocate(work_c(1))
      call ZHEEV('V', 'U', n_keep_s, H_ortho, n_keep_s, eigenvalues_tmp, work_c, lwork, rwork, info)
      lwork = int(real(work_c(1), kind=8)) + 1
      deallocate(work_c)
      allocate(work_c(lwork))
      call ZHEEV('V', 'U', n_keep_s, H_ortho, n_keep_s, eigenvalues_tmp, work_c, lwork, rwork, info)
      if (info /= 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,*) "WARN: H_ortho diagonalization failed in mixed basis, info=", info
          write(*,*) "WARN: Skipping mixed-basis refresh for this spin channel"
        end if
        deallocate(H_work, S_work, eigenvalues_tmp, eval_s, X, H_ortho, tmp_mat, eigvec, work_c, rwork, keep_s_idx)
        deallocate(S_eff, S_eff_work, eval_eff)
        cycle
      end if

      eigvec = matmul(X, H_ortho)

      do i = 1, min(dg_frag%nstate_tot, n_keep_s)
        dg_frag%esp(i, ispin) = eigenvalues_tmp(i)
      end do
      ! The number of propagated states and the rank of the mixed basis are
      ! different quantities.  Keep the full S-effective mixed basis so a raw
      ! fragment/PW state can be represented without an artificial projection
      ! onto only the lowest nstate_tot mixed eigenmodes.
      dg_frag%mixed_basis_dim(ispin) = min(n_total, n_keep_s)
      dg_frag%mixed_transform(1:n_total, 1:dg_frag%mixed_basis_dim(ispin), ispin) = &
        eigvec(1:n_total, 1:dg_frag%mixed_basis_dim(ispin))
      dg_frag%coef_mix(:, :, ispin) = (0.0d0, 0.0d0)
      if (init_identity) then
        do i = 1, min(dg_frag%nstate_tot, dg_frag%mixed_basis_dim(ispin))
          dg_frag%coef_mix(i, i, ispin) = (1.0d0, 0.0d0)
        end do
      else if (init_truncated) then
        allocate(U_keep(n_total, n_keep_s), UH_scaled(n_keep_s, n_total), S_init_metric(n_total, n_total))
        do j = 1, n_keep_s
          k = keep_s_idx(j)
          U_keep(:, j) = S_work(:, k)
          UH_scaled(j, :) = eval_s(k) * conjg(S_work(:, k))
        end do
        S_init_metric(:, :) = matmul(U_keep, UH_scaled)
        call sync_mixed_coef_from_raw(dg_frag, ispin, overlap_metric=S_init_metric)
        deallocate(S_init_metric, UH_scaled, U_keep)
      else
        call sync_mixed_coef_from_raw(dg_frag, ispin)
      end if
      if (init_occupied_projection) then
        nocc_init = dg_frag%nstate_tot
        if (allocated(dg_frag%nocc_spin)) then
          if (ispin <= size(dg_frag%nocc_spin)) nocc_init = min(nocc_init, max(0, dg_frag%nocc_spin(ispin)))
        end if
        nocc_init = min(nocc_init, dg_frag%mixed_basis_dim(ispin))
        if (nocc_init < dg_frag%nstate_tot) then
          dg_frag%coef_mix(:, nocc_init+1:dg_frag%nstate_tot, ispin) = (0.0d0, 0.0d0)
        end if
        if (nocc_init < dg_frag%mixed_basis_dim(ispin)) then
          dg_frag%coef_mix(nocc_init+1:dg_frag%mixed_basis_dim(ispin), 1:dg_frag%nstate_tot, ispin) = (0.0d0, 0.0d0)
        end if
      end if
      call sync_raw_coef_from_mixed(dg_frag, ispin)

      deallocate(H_work, S_work, eigenvalues_tmp, eval_s, X, H_ortho, tmp_mat, eigvec, work_c, rwork, keep_s_idx)
      deallocate(S_eff, S_eff_work, eval_eff)
    end do

    dg_frag%mixed_basis_ready = .true.

    deallocate(S_frag_pw, H_frag_pw, P_frag_pw)

    if (comm_is_root(dg_frag%id)) then
      write(*,*) "Mixed basis diagonalization complete"
      write(*,'(1x,a,f12.6,a)') "Lowest eigenvalue: ", dg_frag%esp(1, 1), " a.u."
      if (dg_frag%nstate_tot > 1) then
        write(*,'(1x,a,f12.6,a)') "Highest occupied energy: ", &
             dg_frag%esp(min(system%no, dg_frag%nstate_tot), 1), " a.u."
      end if
      write(*,*)
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
!$omp parallel do private(iz,iy,ix,gz,gy,gx,vz,vy,vx) schedule(static)
    do iz = lbound(V_box, 3), ubound(V_box, 3)
      gz = gz0 + iz - 1
      vz = map_global_to_phi_box_coord_pw(gz, v_lb3, v_ub3, dg_frag%lgnum_total(3))
      do iy = lbound(V_box, 2), ubound(V_box, 2)
        gy = gy0 + iy - 1
        vy = map_global_to_phi_box_coord_pw(gy, v_lb2, v_ub2, dg_frag%lgnum_total(2))
!$omp simd private(ix,gx,vx)
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
  pure integer function map_global_to_phi_box_coord_pw(ig, phi_lo, phi_hi, lgnum) result(local_idx)
    implicit none
    integer, intent(in) :: ig, phi_lo, phi_hi, lgnum

    local_idx = modulo(ig - 1, lgnum) + 1
    if (local_idx < phi_lo) local_idx = local_idx + lgnum
    if (local_idx > phi_hi) local_idx = local_idx - lgnum
  end function map_global_to_phi_box_coord_pw

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
