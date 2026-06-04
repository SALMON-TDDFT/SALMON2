! DG fragment HSE/RI exchange support shared by non-SOI and SOI RT modules.
#include "config.h"
module rt_dg_hse_exchange
  use structures, only: s_dft_system, s_parallel_info
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use xc_hse_ri, only: hse_ri_data_t, init_hse_ri_fragment, &
                       calc_exact_exchange_hse_ri, deallocate_hse_ri_fragment
  implicit none
  private
  public :: init_hse_ri_data, add_exact_exchange_hse, finalize_hse_ri_data

  type(hse_ri_data_t), allocatable, save :: hse_ri_data_frag(:)
  logical, save :: dg_hse_ace_initialized = .false.
  real(8), allocatable, save :: hse_ace_vx_cache(:,:,:,:)
  complex(8), allocatable, save :: hse_ace_coef_snapshot(:,:,:,:)
  integer, allocatable, save :: hse_ace_last_rebuild(:,:)
  integer, allocatable, save :: hse_ace_call_count(:,:)
  logical, allocatable, save :: hse_ace_cache_valid(:,:)

contains

  subroutine init_hse_ri_data(dg_frag, system, info)
    use salmon_global, only: yn_hse, yn_hse_ri, hse_omega, hse_ri_ratio, &
                             yn_hse_cd_ri, hse_cd_ri_threshold, &
                             yn_dg_hse_ace, dg_hse_ace_max_age, dg_hse_ace_coef_thresh
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info

    integer :: ifrag, ifrag_local, n_frag_local
    integer :: natom_frag, icount, iatom
    real(8), allocatable :: atom_coords_frag(:,:)
    integer, allocatable :: atom_types_frag(:)

    dg_frag%use_hse_ri = (yn_hse == 'y' .and. yn_hse_ri == 'y')
    if (.not. dg_frag%use_hse_ri) return

    if (.not. dg_frag%has_real_space_basis) then
      if (comm_is_root(info%id_rko)) then
        write(*,*) "WARNING: HSE-RI requires real-space basis (phi_frag)"
        write(*,*) "         Disabling HSE-RI approximation"
      end if
      dg_frag%use_hse_ri = .false.
      return
    end if

    n_frag_local = dg_frag%ifrag_end - dg_frag%ifrag_start + 1

    if (comm_is_root(info%id_rko)) then
      write(*,*)
      write(*,*) "Initializing RI/DF approximation for HSE exchange..."
      write(*,'(1x,a,i0)') "  Local fragments: ", n_frag_local
      write(*,'(1x,a,f6.2)') "  N_aux/N_basis ratio: ", hse_ri_ratio
    end if

    allocate(hse_ri_data_frag(n_frag_local))
    call reset_hse_ace_storage()
    if (yn_dg_hse_ace == 'y') then
      allocate(hse_ace_vx_cache(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin, n_frag_local))
      allocate(hse_ace_coef_snapshot(dg_frag%nstate_frag, dg_frag%nstate_tot, dg_frag%nspin, n_frag_local))
      allocate(hse_ace_last_rebuild(n_frag_local, dg_frag%nspin))
      allocate(hse_ace_call_count(n_frag_local, dg_frag%nspin))
      allocate(hse_ace_cache_valid(n_frag_local, dg_frag%nspin))
      hse_ace_vx_cache = 0.0d0
      hse_ace_coef_snapshot = (0.0d0, 0.0d0)
      hse_ace_last_rebuild = 0
      hse_ace_call_count = 0
      hse_ace_cache_valid = .false.
      dg_hse_ace_initialized = .true.
    else
      dg_hse_ace_initialized = .false.
    end if

    do ifrag_local = 1, n_frag_local
      ifrag = dg_frag%ifrag_start + ifrag_local - 1
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a,i0,a,i0)') "  Initializing fragment ", ifrag, " (local ", ifrag_local, ")..."
      end if

      natom_frag = 0
      do iatom = 1, system%nion
        if (atom_belongs_to_fragment(dg_frag, system, ifrag, iatom)) natom_frag = natom_frag + 1
      end do

      allocate(atom_coords_frag(3, natom_frag))
      allocate(atom_types_frag(natom_frag))
      icount = 0
      do iatom = 1, system%nion
        if (.not. atom_belongs_to_fragment(dg_frag, system, ifrag, iatom)) cycle
        icount = icount + 1
        atom_coords_frag(:, icount) = system%Rion(:, iatom)
        atom_types_frag(icount) = system%kion(iatom)
      end do

      call init_hse_ri_fragment(hse_ri_data_frag(ifrag_local), &
                                dg_frag%phi_frag(:,:,:,:,ifrag_local), &
                                dg_frag%lg, dg_frag%mg, &
                                dg_frag%nstate_frag, system%hvol, &
                                natom_frag, atom_coords_frag, atom_types_frag, &
                                hse_omega, &
                                (yn_hse_cd_ri == 'y'), hse_cd_ri_threshold, &
                                max(dg_frag%mg%is(1:3), dg_frag%ixyz_frag(1:3, ifrag)), &
                                min(dg_frag%mg%ie(1:3), dg_frag%ixyz_frag(1:3, ifrag) + &
                                    dg_frag%nxyz_domain(1:3, ifrag) - 1))

      deallocate(atom_coords_frag, atom_types_frag)
    end do

    if (comm_is_root(info%id_rko)) then
      write(*,*) "RI/DF initialization complete!"
      if (yn_dg_hse_ace == 'y') then
        write(*,'(1x,a)') "DG-HSE-ACE cache: ENABLED"
        write(*,'(1x,a,i0)') "  max age (calls): ", dg_hse_ace_max_age
        write(*,'(1x,a,1pe12.4)') "  coef threshold : ", dg_hse_ace_coef_thresh
      else
        write(*,'(1x,a)') "DG-HSE-ACE cache: DISABLED"
      end if
      write(*,*)
    end if
  end subroutine init_hse_ri_data

  subroutine add_exact_exchange_hse(dg_frag, system, H_mat_spin, ifrag, ispin)
    use salmon_global, only: hse_alpha, hse_omega, nelec, nelec_spin, yn_hse_ri, &
                             yn_dg_hse_ace, dg_hse_ace_max_age, dg_hse_ace_coef_thresh, yn_spinorbit
    use xc_hse, only: calc_exact_exchange_hse_fragment
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    real(8), intent(inout) :: H_mat_spin(:, :)
    integer, intent(in) :: ifrag, ispin

    integer :: ifrag_local, n_base_frag, n_occ_frag
    integer :: is(3), ie(3)
    real(8) :: hvol, ace_metric
    real(8), allocatable :: density_matrix(:,:), v_x_tmp(:,:)
    logical :: rebuild_ace

    if (.not. dg_frag%has_real_space_basis) return

    hvol = system%hvol
    is = 1
    ie = dg_frag%nxyz_domain(1:3, ifrag)
    n_base_frag = dg_frag%nstate_frag
    ifrag_local = ifrag - dg_frag%ifrag_start + 1
    n_occ_frag = fragment_hse_nocc(dg_frag, ispin, n_base_frag, nelec, nelec_spin, yn_spinorbit)

    if (dg_frag%use_hse_ri .and. yn_hse_ri == 'y') then
      if (dg_hse_ace_initialized .and. yn_dg_hse_ace == 'y') then
        hse_ace_call_count(ifrag_local, ispin) = hse_ace_call_count(ifrag_local, ispin) + 1
        rebuild_ace = .not. hse_ace_cache_valid(ifrag_local, ispin)
        ace_metric = 0.0d0
        if (.not. rebuild_ace) then
          call compute_dg_hse_ace_metric(dg_frag, ifrag, ifrag_local, ispin, n_occ_frag, ace_metric)
          if (ace_metric > dg_hse_ace_coef_thresh) rebuild_ace = .true.
          if ((hse_ace_call_count(ifrag_local, ispin) - hse_ace_last_rebuild(ifrag_local, ispin)) >= &
              dg_hse_ace_max_age) rebuild_ace = .true.
        end if

        if (rebuild_ace) then
          allocate(density_matrix(n_base_frag, n_base_frag))
          allocate(v_x_tmp(n_base_frag, n_base_frag))
          call build_density_matrix(dg_frag, ifrag, ispin, n_occ_frag, density_matrix)
          v_x_tmp = 0.0d0
          call calc_exact_exchange_hse_ri(v_x_tmp, hse_ri_data_frag(ifrag_local), density_matrix, hse_alpha, n_occ_frag)
          H_mat_spin(1:n_base_frag, 1:n_base_frag) = H_mat_spin(1:n_base_frag, 1:n_base_frag) + v_x_tmp
          hse_ace_vx_cache(1:n_base_frag, 1:n_base_frag, ispin, ifrag_local) = v_x_tmp
          call update_dg_hse_ace_snapshot(dg_frag, ifrag, ifrag_local, ispin, n_occ_frag)
          hse_ace_cache_valid(ifrag_local, ispin) = .true.
          hse_ace_last_rebuild(ifrag_local, ispin) = hse_ace_call_count(ifrag_local, ispin)
          deallocate(v_x_tmp, density_matrix)
        else
          H_mat_spin(1:n_base_frag, 1:n_base_frag) = H_mat_spin(1:n_base_frag, 1:n_base_frag) + &
            hse_ace_vx_cache(1:n_base_frag, 1:n_base_frag, ispin, ifrag_local)
        end if
      else
        allocate(density_matrix(n_base_frag, n_base_frag))
        call build_density_matrix(dg_frag, ifrag, ispin, n_occ_frag, density_matrix)
        call calc_exact_exchange_hse_ri(H_mat_spin, hse_ri_data_frag(ifrag_local), &
                                        density_matrix, hse_alpha, n_occ_frag)
        deallocate(density_matrix)
      end if
    else
      call calc_exact_exchange_hse_fragment(H_mat_spin, dg_frag%phi_frag, ifrag_local, &
                                            dg_frag%hgs, hvol, hse_alpha, &
                                            is, ie, n_base_frag, n_occ_frag, nelec)
    end if
  end subroutine add_exact_exchange_hse

  subroutine finalize_hse_ri_data()
    implicit none
    integer :: ifrag_local

    if (allocated(hse_ri_data_frag)) then
      do ifrag_local = 1, size(hse_ri_data_frag)
        call deallocate_hse_ri_fragment(hse_ri_data_frag(ifrag_local))
      end do
      deallocate(hse_ri_data_frag)
    end if
    call reset_hse_ace_storage()
    dg_hse_ace_initialized = .false.
  end subroutine finalize_hse_ri_data

  integer function fragment_hse_nocc(dg_frag, ispin, n_base_frag, nelec, nelec_spin, yn_spinorbit) result(n_occ)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ispin, n_base_frag, nelec, nelec_spin(:)
    character(len=*), intent(in) :: yn_spinorbit

    if (yn_spinorbit == 'y') then
      n_occ = max(1, min(nelec, n_base_frag))
    else if (dg_frag%nspin == 1) then
      n_occ = max(1, min((nelec + 1) / 2, n_base_frag))
    else if (sum(nelec_spin(:)) > 0) then
      n_occ = max(1, min(nelec_spin(ispin), n_base_frag))
    else
      n_occ = max(1, min(int(nelec / 2.0d0 + 1.0d-12), n_base_frag))
    end if
  end function fragment_hse_nocc

  logical function atom_belongs_to_fragment(dg_frag, system, ifrag, iatom) result(is_owned)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    integer, intent(in) :: ifrag, iatom
    integer :: frag_idx(3), axis, g0, g1, ndiv, pos
    real(8) :: rmin, rmax, ratom

    is_owned = .true.
    frag_idx = fragment_cartesian_index(dg_frag, ifrag)
    do axis = 1, 3
      g0 = dg_frag%ixyz_frag(axis, ifrag)
      g1 = g0 + dg_frag%nxyz_domain(axis, ifrag) - 1
      rmin = dg_frag%mg%coordinate(g0, axis)
      rmax = dg_frag%mg%coordinate(g1, axis)
      ratom = system%Rion(axis, iatom)
      ndiv = max(1, dg_frag%num_fragment(axis))
      pos = frag_idx(axis)
      if (pos < ndiv) then
        if (ratom < rmin .or. ratom >= rmax) is_owned = .false.
      else
        if (ratom < rmin .or. ratom > rmax) is_owned = .false.
      end if
      if (.not. is_owned) return
    end do
  end function atom_belongs_to_fragment

  function fragment_cartesian_index(dg_frag, ifrag) result(idx)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    integer :: idx(3), rem

    rem = ifrag - 1
    idx(1) = mod(rem, dg_frag%num_fragment(1)) + 1
    rem = rem / dg_frag%num_fragment(1)
    idx(2) = mod(rem, dg_frag%num_fragment(2)) + 1
    rem = rem / dg_frag%num_fragment(2)
    idx(3) = mod(rem, dg_frag%num_fragment(3)) + 1
  end function fragment_cartesian_index

  subroutine compute_dg_hse_ace_metric(dg_frag, ifrag, ifrag_local, ispin, n_occ, metric)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, ifrag_local, ispin, n_occ
    real(8), intent(out) :: metric
    integer :: io, istate, ig, iloc, nb
    complex(8) :: c_now, c_prev
    real(8) :: amp_now, amp_prev, amp_num, amp_den
    real(8) :: phase_now, phase_prev, phase_diff, phase_num, phase_den, weight
    real(8), parameter :: pi = 3.14159265358979323846d0

    metric = 0.0d0
    amp_num = 0.0d0
    amp_den = 0.0d0
    phase_num = 0.0d0
    phase_den = 0.0d0
    nb = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)

    do io = 1, nb
      ig = dg_frag%index_basis(io, ifrag, ispin)
      if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
      if (.not. allocated(dg_frag%coef_global_to_local)) cycle
      iloc = dg_frag%coef_global_to_local(ig, ispin)
      if (iloc < 1 .or. iloc > size(dg_frag%coef, 1)) cycle
      do istate = 1, min(n_occ, dg_frag%nstate_tot)
        c_now = dg_frag%coef(iloc, istate, ispin)
        c_prev = hse_ace_coef_snapshot(io, istate, ispin, ifrag_local)
        amp_now = abs(c_now)
        amp_prev = abs(c_prev)
        amp_num = amp_num + (amp_now - amp_prev) * (amp_now - amp_prev)
        amp_den = amp_den + amp_prev * amp_prev

        if (amp_now > 1.0d-14 .and. amp_prev > 1.0d-14) then
          phase_now = atan2(aimag(c_now), real(c_now, kind=8))
          phase_prev = atan2(aimag(c_prev), real(c_prev, kind=8))
          phase_diff = phase_now - phase_prev
          if (phase_diff > pi) phase_diff = phase_diff - 2.0d0 * pi
          if (phase_diff < -pi) phase_diff = phase_diff + 2.0d0 * pi
          weight = 0.5d0 * (amp_now + amp_prev)
          phase_num = phase_num + (weight * phase_diff) * (weight * phase_diff)
          phase_den = phase_den + weight * weight
        end if
      end do
    end do

    if (amp_den > 1.0d-30) then
      metric = metric + (amp_num / amp_den)
    else
      metric = metric + amp_num
    end if
    if (phase_den > 1.0d-30) then
      metric = metric + (phase_num / phase_den)
    else
      metric = metric + phase_num
    end if
    metric = sqrt(max(metric, 0.0d0))
  end subroutine compute_dg_hse_ace_metric

  subroutine update_dg_hse_ace_snapshot(dg_frag, ifrag, ifrag_local, ispin, n_occ)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, ifrag_local, ispin, n_occ
    integer :: io, ig, iloc, nocc_eff

    nocc_eff = min(n_occ, dg_frag%nstate_tot)
    hse_ace_coef_snapshot(:, :, ispin, ifrag_local) = (0.0d0, 0.0d0)
    do io = 1, min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
      ig = dg_frag%index_basis(io, ifrag, ispin)
      if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
      if (.not. allocated(dg_frag%coef_global_to_local)) cycle
      iloc = dg_frag%coef_global_to_local(ig, ispin)
      if (iloc < 1 .or. iloc > size(dg_frag%coef, 1)) cycle
      hse_ace_coef_snapshot(io, 1:nocc_eff, ispin, ifrag_local) = dg_frag%coef(iloc, 1:nocc_eff, ispin)
    end do
  end subroutine update_dg_hse_ace_snapshot

  subroutine build_density_matrix(dg_frag, ifrag, ispin, n_occ, density_matrix)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, ispin, n_occ
    real(8), intent(out) :: density_matrix(:, :)

    integer :: i, j, iocc, istate, iidx, jidx, iloc, jloc, nb
    complex(8) :: coef_i, coef_j
    real(8) :: dval

    density_matrix = 0.0d0
    nb = min(size(density_matrix, 1), dg_frag%nstate_frag)
    do j = 1, nb
      jidx = dg_frag%index_basis(j, ifrag, ispin)
      if (jidx < 1 .or. jidx > dg_frag%n_mat_max) cycle
      if (.not. allocated(dg_frag%coef_global_to_local)) cycle
      jloc = dg_frag%coef_global_to_local(jidx, ispin)
      if (jloc < 1 .or. jloc > size(dg_frag%coef, 1)) cycle
      do i = 1, nb
        iidx = dg_frag%index_basis(i, ifrag, ispin)
        if (iidx < 1 .or. iidx > dg_frag%n_mat_max) cycle
        iloc = dg_frag%coef_global_to_local(iidx, ispin)
        if (iloc < 1 .or. iloc > size(dg_frag%coef, 1)) cycle
        do iocc = 1, n_occ
          istate = iocc
          if (istate > dg_frag%nstate_tot) cycle
          coef_i = dg_frag%coef(iloc, istate, ispin)
          coef_j = dg_frag%coef(jloc, istate, ispin)
          dval = real(conjg(coef_i) * coef_j, kind=8)
          density_matrix(i, j) = density_matrix(i, j) + dval
        end do
      end do
    end do

    do j = 1, nb
      do i = j + 1, nb
        dval = 0.5d0 * (density_matrix(i, j) + density_matrix(j, i))
        density_matrix(i, j) = dval
        density_matrix(j, i) = dval
      end do
    end do
  end subroutine build_density_matrix

  subroutine reset_hse_ace_storage()
    implicit none
    if (allocated(hse_ace_vx_cache)) deallocate(hse_ace_vx_cache)
    if (allocated(hse_ace_coef_snapshot)) deallocate(hse_ace_coef_snapshot)
    if (allocated(hse_ace_last_rebuild)) deallocate(hse_ace_last_rebuild)
    if (allocated(hse_ace_call_count)) deallocate(hse_ace_call_count)
    if (allocated(hse_ace_cache_valid)) deallocate(hse_ace_cache_valid)
  end subroutine reset_hse_ace_storage

end module rt_dg_hse_exchange
