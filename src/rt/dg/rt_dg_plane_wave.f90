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

    do ipw = 1, n_pw
      i = n_frag + ipw
      if (allocated(dg_frag%H_mat_pw_diag)) mat(i, i) = dg_frag%H_mat_pw_diag(ipw, ispin)
    end do
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
    complex(8), allocatable :: S_frag_pw(:,:,:), H_frag_pw(:,:,:)

    n_frag = dg_frag%n_mat_max
    n_pw = dg_frag%n_plane_waves
    if (.not. dg_frag%use_plane_wave_basis .or. n_pw <= 0) return

    allocate(S_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(H_frag_pw(n_frag, n_pw, dg_frag%nspin))

    call compute_fragment_pw_overlap(dg_frag, S_frag_pw)
    call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_frag_pw)

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

    call build_mixed_hamiltonian(dg_frag, dg_frag%lg, Vh, Vxc, Vpsl, Ac_tot, S_frag_pw, H_frag_pw)

    if (allocated(dg_frag%coef_pw)) dg_frag%coef_pw(:, :, :) = (0.0d0, 0.0d0)

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "  [init] Prepared mixed FP/PP operators without dense diagonalization"
      write(*,'(1x,a)') "  [init] PW coefficients start from zero and couple in during propagation"
    end if

    deallocate(S_frag_pw, H_frag_pw)
  end subroutine prepare_mixed_basis_startup

  !=======================================================================
  ! Compute overlap matrix between fragment basis and plane waves
  !=======================================================================
  subroutine compute_fragment_pw_overlap(dg_frag, S_complex)
    use structures
    use communication, only: comm_bcast
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), intent(out) :: S_complex(:,:,:)  ! (n_mat_max, n_pw, nspin)

    integer :: ipw, ifrag, i_local, io, ig, ispin, ix, iy, iz
    integer :: nx, ny, nz
    integer :: gx0, gy0, gz0
    real(8) :: k_vec(3), Lbox(3), sqrt_V, inv_sqrt_V
    real(8) :: vol_elem
    complex(8) :: pw_val, overlap_local, phase_yz
    complex(8) :: step_x, step_y, step_z
    complex(8) :: phase_x0, phase_y0, phase_z0
    complex(8), allocatable :: frag_block(:,:,:)
    complex(8), allocatable :: phase_x(:), phase_y(:), phase_z(:)

    if (.not. dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return

    vol_elem = product(dg_frag%hgs(1:3))
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))
    sqrt_V = sqrt(product(Lbox))
    inv_sqrt_V = 1.0d0 / sqrt_V

    S_complex = (0.0d0, 0.0d0)

    do ispin = 1, dg_frag%nspin
      do ipw = 1, dg_frag%n_plane_waves
        k_vec(1:3) = dg_frag%k_pw(1:3, ipw)

        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          nx = dg_frag%nxyz_domain(1, ifrag)
          ny = dg_frag%nxyz_domain(2, ifrag)
          nz = dg_frag%nxyz_domain(3, ifrag)
          gx0 = dg_frag%ixyz_frag(1, ifrag)
          gy0 = dg_frag%ixyz_frag(2, ifrag)
          gz0 = dg_frag%ixyz_frag(3, ifrag)

          if (.not. allocated(phase_x) .or. size(phase_x) < nx) then
            if (allocated(phase_x)) deallocate(phase_x)
            allocate(phase_x(nx))
          end if
          if (.not. allocated(phase_y) .or. size(phase_y) < ny) then
            if (allocated(phase_y)) deallocate(phase_y)
            allocate(phase_y(ny))
          end if
          if (.not. allocated(phase_z) .or. size(phase_z) < nz) then
            if (allocated(phase_z)) deallocate(phase_z)
            allocate(phase_z(nz))
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

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig = dg_frag%index_basis(io, ifrag, ispin)
            if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
            overlap_local = (0.0d0, 0.0d0)

            do iz = 1, nz
              do iy = 1, ny
                phase_yz = phase_y(iy) * phase_z(iz)
                do ix = 1, nx
                  pw_val = phase_x(ix) * phase_yz * inv_sqrt_V

                  if (allocated(dg_frag%phi_frag_c)) then
                    overlap_local = overlap_local + &
                         conjg(dg_frag%phi_frag_c(ix,iy,iz,io,i_local)) * pw_val * vol_elem
                  else
                    overlap_local = overlap_local + &
                         dg_frag%phi_frag(ix,iy,iz,io,i_local) * conjg(pw_val) * vol_elem
                  end if
                end do
              end do
            end do

            S_complex(ig, ipw, ispin) = S_complex(ig, ipw, ispin) + overlap_local
          end do
        end do
      end do
    end do

    if (allocated(phase_x)) deallocate(phase_x)
    if (allocated(phase_y)) deallocate(phase_y)
    if (allocated(phase_z)) deallocate(phase_z)

    allocate(frag_block(dg_frag%nstate_frag, dg_frag%n_plane_waves, dg_frag%nspin))
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
    deallocate(frag_block)

  end subroutine compute_fragment_pw_overlap

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
      do iz = dg_frag%lg%is(3), dg_frag%lg%ie(3)
        do iy = dg_frag%lg%is(2), dg_frag%lg%ie(2)
          do ix = dg_frag%lg%is(1), dg_frag%lg%ie(1)
            integral_local = integral_local + &
                 (Vpsl%f(ix,iy,iz) + Vh%f(ix,iy,iz) + Vxc(ispin)%f(ix,iy,iz)) * vol_elem
          end do
        end do
      end do
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
    use communication, only: comm_bcast
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_scalar),         intent(in)    :: Vh, Vxc(:), Vpsl
    complex(8),             intent(out)   :: H_complex(:,:,:)  ! (n_mat_max, n_pw, nspin)

    integer :: ipw, ifrag, i_local, io, ig, ispin, ix, iy, iz
    integer :: nx, ny, nz
    integer :: gx, gy, gz, gx0, gy0, gz0
    real(8) :: k_vec(3), Lbox(3), sqrt_V, inv_sqrt_V
    real(8) :: vol_elem, k_squared, V_total
    complex(8) :: pw_val, pw_laplacian, hamiltonian_local, phase_yz
    complex(8) :: step_x, step_y, step_z
    complex(8) :: phase_x0, phase_y0, phase_z0
    complex(8), allocatable :: frag_block(:,:,:)
    complex(8), allocatable :: phase_x(:), phase_y(:), phase_z(:)

    if (.not. dg_frag%use_plane_wave_basis) return
    if (.not. dg_frag%has_real_space_basis) return

    vol_elem = product(dg_frag%hgs(1:3))
    Lbox(1:3) = dg_frag%hgs(1:3) * dble(dg_frag%lgnum_total(1:3))
    sqrt_V = sqrt(product(Lbox))
    inv_sqrt_V = 1.0d0 / sqrt_V

    H_complex = (0.0d0, 0.0d0)

    do ispin = 1, dg_frag%nspin
      do ipw = 1, dg_frag%n_plane_waves
        k_vec(1:3) = dg_frag%k_pw(1:3, ipw)
        k_squared = sum(k_vec**2)

        i_local = 0
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          i_local = i_local + 1
          nx = dg_frag%nxyz_domain(1, ifrag)
          ny = dg_frag%nxyz_domain(2, ifrag)
          nz = dg_frag%nxyz_domain(3, ifrag)
          gx0 = dg_frag%ixyz_frag(1, ifrag)
          gy0 = dg_frag%ixyz_frag(2, ifrag)
          gz0 = dg_frag%ixyz_frag(3, ifrag)

          if (.not. allocated(phase_x) .or. size(phase_x) < nx) then
            if (allocated(phase_x)) deallocate(phase_x)
            allocate(phase_x(nx))
          end if
          if (.not. allocated(phase_y) .or. size(phase_y) < ny) then
            if (allocated(phase_y)) deallocate(phase_y)
            allocate(phase_y(ny))
          end if
          if (.not. allocated(phase_z) .or. size(phase_z) < nz) then
            if (allocated(phase_z)) deallocate(phase_z)
            allocate(phase_z(nz))
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

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig = dg_frag%index_basis(io, ifrag, ispin)
            if (ig < 1 .or. ig > dg_frag%n_mat_max) cycle
            hamiltonian_local = (0.0d0, 0.0d0)

            do iz = 1, nz
              gz = gz0 + iz - 1
              do iy = 1, ny
                gy = gy0 + iy - 1
                phase_yz = phase_y(iy) * phase_z(iz)
                do ix = 1, nx
                  gx = gx0 + ix - 1
                  pw_val = phase_x(ix) * phase_yz * inv_sqrt_V
                  pw_laplacian = (k_squared / 2.0d0) * pw_val

                  V_total = Vpsl%f(gx, gy, gz) + Vh%f(gx, gy, gz) + Vxc(ispin)%f(gx, gy, gz)

                  if (allocated(dg_frag%phi_frag_c)) then
                    hamiltonian_local = hamiltonian_local + &
                         conjg(dg_frag%phi_frag_c(ix,iy,iz,io,i_local)) * &
                         (pw_laplacian + V_total * pw_val) * vol_elem
                  else
                    hamiltonian_local = hamiltonian_local + &
                         dg_frag%phi_frag(ix,iy,iz,io,i_local) * &
                         conjg(pw_laplacian + V_total * pw_val) * vol_elem
                  end if
                end do
              end do
            end do

            H_complex(ig, ipw, ispin) = H_complex(ig, ipw, ispin) + hamiltonian_local
          end do
        end do
      end do
    end do

    if (allocated(phase_x)) deallocate(phase_x)
    if (allocated(phase_y)) deallocate(phase_y)
    if (allocated(phase_z)) deallocate(phase_z)

    allocate(frag_block(dg_frag%nstate_frag, dg_frag%n_plane_waves, dg_frag%nspin))
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
    deallocate(frag_block)

  end subroutine compute_fragment_pw_hamiltonian

  !=======================================================================
  ! Build mixed Hamiltonian matrix with plane wave basis
  !=======================================================================
  subroutine build_mixed_hamiltonian(dg_frag, lg, Vh, Vxc, Vpsl, Ac_tot, &
                                      S_frag_pw, H_frag_pw)
    use structures
    use rt_dg_fragment_ops, only: copy_matrix_blocks_to_complex_dense
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in)    :: lg
    type(s_scalar),         intent(in)    :: Vh, Vxc(:), Vpsl
    real(8),                intent(in)    :: Ac_tot(3)
    complex(8),             intent(in)    :: S_frag_pw(:,:,:)
    complex(8),             intent(in)    :: H_frag_pw(:,:,:)

    integer :: ispin, i, ipw, n_frag, n_pw
    real(8) :: k_vec(3), kinetic_energy
    real(8), allocatable :: V_mean(:)

    if (.not. dg_frag%use_plane_wave_basis) return

    n_frag = dg_frag%n_mat_max
    n_pw = dg_frag%n_plane_waves

    allocate(V_mean(dg_frag%nspin))
    call compute_mean_potential(dg_frag, Vh, Vxc, Vpsl, V_mean)
    if (.not. allocated(dg_frag%H_mat_pw_diag)) then
      allocate(dg_frag%H_mat_pw_diag(n_pw, dg_frag%nspin))
    else if (size(dg_frag%H_mat_pw_diag,1) /= n_pw .or. size(dg_frag%H_mat_pw_diag,2) /= dg_frag%nspin) then
      deallocate(dg_frag%H_mat_pw_diag)
      allocate(dg_frag%H_mat_pw_diag(n_pw, dg_frag%nspin))
    end if

    do ispin = 1, dg_frag%nspin
      if (n_pw <= 0) cycle

      do ipw = 1, n_pw
        k_vec(1:3) = dg_frag%k_pw(1:3, ipw)
        kinetic_energy = 0.5d0 * sum(k_vec**2)
        i = n_frag + ipw
        dg_frag%H_mat_pw_diag(ipw, ispin) = cmplx(kinetic_energy + V_mean(ispin), 0.0d0, kind=8)
      end do
    end do

    deallocate(V_mean)

  end subroutine build_mixed_hamiltonian

  !=======================================================================
  ! Diagonalize mixed basis using generalized eigenproblem:
  !   H c = eps S c
  !=======================================================================
  subroutine diagonalize_mixed_basis(dg_frag, system, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use communication, only: comm_is_root
    use rt_dg_fragment_ops, only: copy_matrix_blocks_to_complex_dense, sync_raw_coef_from_mixed
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_scalar), intent(in) :: Vh, Vxc(:), Vpsl
    real(8), intent(in) :: Ac_tot(3)

    integer :: ispin, n_total, n_frag, n_pw, lda, lwork, info, i, j
    integer :: n_floor
    real(8) :: sij, s_eff_min, s_eff_max, s_eff_cond
    logical :: use_raw_prop_s
    complex(8), allocatable :: H_work(:,:), S_work(:,:), work_c(:)
    complex(8), allocatable :: X(:,:), H_ortho(:,:), tmp_mat(:,:), eigvec(:,:), S_rebuild(:,:)
    real(8), allocatable :: eigenvalues_tmp(:), eval_s(:), rwork(:), work(:)
    real(8), allocatable :: S_eff(:,:), S_eff_work(:,:), eval_eff(:), work_eff(:)
    real(8), allocatable :: A_qr(:,:), tau_qr(:), work_qr(:)
    complex(8), allocatable :: coef_mixed(:,:,:)
    complex(8), allocatable :: S_frag_pw(:,:,:)  ! Complex overlap matrix
    complex(8), allocatable :: H_frag_pw(:,:,:)  ! Hamiltonian coupling matrix
    integer, allocatable :: jpvt(:), keep_idx(:)
    integer :: m_qr, n_qr, lwork_qr, info_qr, n_keep_pw, ndiag
    real(8) :: diag_max, tau_rr
    real(8), parameter :: eps_s_abs = 1.0d-10
    real(8), parameter :: eps_s_rel = 1.0d-8
    real(8), parameter :: tau_pw_rank_rel = 1.0d-6
    real(8) :: eps_s
    external :: dsyev, dgeqp3, zheev

    n_frag = dg_frag%n_mat_max
    n_pw = dg_frag%n_plane_waves
    n_total = n_frag + n_pw

    if (comm_is_root(dg_frag%id)) then
      write(*,*)
      write(*,*) "=== Diagonalizing mixed basis (Fragment + PW) ==="
      write(*,'(1x,a,i0)') "Fragment basis size: ", n_frag
      write(*,'(1x,a,i0)') "Plane wave basis size: ", n_pw
      write(*,'(1x,a,i0)') "Total mixed basis size: ", n_total
    end if

    allocate(S_frag_pw(n_frag, n_pw, dg_frag%nspin))
    allocate(H_frag_pw(n_frag, n_pw, dg_frag%nspin))
    call compute_fragment_pw_overlap(dg_frag, S_frag_pw)
    if (n_pw > 0 .and. n_frag > 0) then
      m_qr = 2 * n_frag * dg_frag%nspin
      n_qr = n_pw
      allocate(A_qr(m_qr, n_qr), jpvt(n_qr), tau_qr(min(m_qr, n_qr)))
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

      if (n_keep_pw < n_pw) then
        allocate(keep_idx(n_keep_pw))
        keep_idx(1:n_keep_pw) = jpvt(1:n_keep_pw)
        call compact_plane_wave_basis(dg_frag, S_frag_pw, H_frag_pw, keep_idx, n_keep_pw)
        deallocate(keep_idx)
        n_pw = dg_frag%n_plane_waves
        n_total = n_frag + n_pw
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a,i0,a,1pe11.3)') "PW rank selection: kept ", n_keep_pw, " / ", n_qr, " (tau=", tau_rr, ")"
        end if
      end if
      deallocate(A_qr, jpvt, tau_qr, work_qr)
    end if

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
    call compute_fragment_pw_hamiltonian(dg_frag, Vh, Vxc, Vpsl, H_frag_pw)
    if (.not. allocated(dg_frag%H_mat_frag_pw)) then
      allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    else if (size(dg_frag%H_mat_frag_pw,1) /= n_frag .or. size(dg_frag%H_mat_frag_pw,2) /= n_pw .or. &
             size(dg_frag%H_mat_frag_pw,3) /= dg_frag%nspin) then
      deallocate(dg_frag%H_mat_frag_pw)
      allocate(dg_frag%H_mat_frag_pw(n_frag, n_pw, dg_frag%nspin))
    end if
    dg_frag%H_mat_frag_pw(1:n_frag,1:n_pw,1:dg_frag%nspin) = H_frag_pw(1:n_frag,1:n_pw,1:dg_frag%nspin)
    call build_mixed_hamiltonian(dg_frag, dg_frag%lg, Vh, Vxc, Vpsl, Ac_tot, S_frag_pw, H_frag_pw)

    allocate(coef_mixed(n_total, dg_frag%nstate_tot, dg_frag%nspin))

    do ispin = 1, dg_frag%nspin
      lda = n_total

      allocate(H_work(n_total, n_total), S_work(n_total, n_total))
      allocate(eigenvalues_tmp(n_total))
      allocate(eval_s(n_total), X(n_total, n_total))
      allocate(H_ortho(n_total, n_total), tmp_mat(n_total, n_total), eigvec(n_total, n_total))
      allocate(S_rebuild(n_total, n_total))
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
        eps_s = max(eps_s_abs, eps_s_rel * max(s_eff_max, 1.0d0))
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

      ! Regularize generalized EVP by Lowdin transform with eigenvalue floor.
      ! S = U diag(s) U^T, X = U diag(1/sqrt(max(s, eps_s))), H_ortho = X^T H X.
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
        deallocate(H_work, S_work, eigenvalues_tmp, eval_s, X, H_ortho, tmp_mat, eigvec, S_rebuild, work_c, rwork)
        deallocate(S_eff, S_eff_work, eval_eff)
        cycle
      end if
      eps_s = max(eps_s_abs, eps_s_rel * max(eval_s(n_total), 1.0d0))

      n_floor = 0
      do i = 1, n_total
        if (eval_s(i) < eps_s) then
          eval_s(i) = eps_s
          n_floor = n_floor + 1
        end if
      end do
      if (n_floor > 0 .and. comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,1pe11.3)') "Mixed-S regularization: floored ", n_floor, " eigenvalues below ", eps_s
      end if

      X(:, :) = S_work(:, :)
      do j = 1, n_total
        X(:, j) = X(:, j) / sqrt(eval_s(j))
      end do
      tmp_mat = matmul(conjg(transpose(X)), H_work)
      H_ortho = matmul(tmp_mat, X)
      do j = 1, n_total
        H_ortho(j, j) = cmplx(real(H_ortho(j, j), kind=8), 0.0d0, kind=8)
        do i = j + 1, n_total
          H_ortho(i, j) = conjg(H_ortho(j, i))
        end do
      end do

      deallocate(work_c)
      lwork = -1
      allocate(work_c(1))
      call ZHEEV('V', 'U', n_total, H_ortho, lda, eigenvalues_tmp, work_c, lwork, rwork, info)
      lwork = int(real(work_c(1), kind=8)) + 1
      deallocate(work_c)
      allocate(work_c(lwork))
      call ZHEEV('V', 'U', n_total, H_ortho, lda, eigenvalues_tmp, work_c, lwork, rwork, info)
      if (info /= 0) then
        if (comm_is_root(dg_frag%id)) then
          write(*,*) "WARN: H_ortho diagonalization failed in mixed basis, info=", info
          write(*,*) "WARN: Skipping mixed-basis refresh for this spin channel"
        end if
        deallocate(H_work, S_work, eigenvalues_tmp, eval_s, X, H_ortho, tmp_mat, eigvec, S_rebuild, work_c, rwork)
        deallocate(S_eff, S_eff_work, eval_eff)
        cycle
      end if

      eigvec = matmul(X, H_ortho)

      do i = 1, min(dg_frag%nstate_tot, n_total)
        dg_frag%esp(i, ispin) = eigenvalues_tmp(i)
      end do
      dg_frag%mixed_basis_dim(ispin) = min(dg_frag%nstate_tot, n_total)
      dg_frag%mixed_transform(1:n_total, 1:dg_frag%mixed_basis_dim(ispin), ispin) = &
        eigvec(1:n_total, 1:dg_frag%mixed_basis_dim(ispin))
      dg_frag%coef_mix(:, :, ispin) = (0.0d0, 0.0d0)
      do i = 1, dg_frag%mixed_basis_dim(ispin)
        dg_frag%coef_mix(i, i, ispin) = (1.0d0, 0.0d0)
      end do
      call sync_raw_coef_from_mixed(dg_frag, ispin)

      coef_mixed(:, :, ispin) = (0.0d0, 0.0d0)
      coef_mixed(:, :, ispin) = dg_frag%mixed_transform(1:n_total, 1:n_total, ispin)

      deallocate(H_work, S_work, eigenvalues_tmp, eval_s, X, H_ortho, tmp_mat, eigvec, S_rebuild, work_c, rwork)
      deallocate(S_eff, S_eff_work, eval_eff)
    end do

    dg_frag%mixed_basis_ready = .true.

    deallocate(coef_mixed, S_frag_pw, H_frag_pw)

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
  ! Keep only rank-revealed independent PW columns.
  !=======================================================================
  subroutine compact_plane_wave_basis(dg_frag, S_frag_pw, H_frag_pw, keep_idx, n_keep)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    complex(8), allocatable, intent(inout) :: S_frag_pw(:,:,:)
    complex(8), allocatable, intent(inout) :: H_frag_pw(:,:,:)
    integer, intent(in) :: keep_idx(:)
    integer, intent(in) :: n_keep

    integer :: n_frag, n_pw_old, i
    real(8), allocatable :: k_new(:,:)
    complex(8), allocatable :: coef_new(:,:,:), Sfp_new(:,:,:), Hfp_new(:,:,:), Hpp_new(:,:)

    n_frag = dg_frag%n_mat_max
    n_pw_old = dg_frag%n_plane_waves
    if (n_keep <= 0 .or. n_keep >= n_pw_old) return

    allocate(k_new(3, n_keep), coef_new(n_keep, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(Sfp_new(n_frag, n_keep, dg_frag%nspin), Hfp_new(n_frag, n_keep, dg_frag%nspin), Hpp_new(n_keep, dg_frag%nspin))
    coef_new(:, :, :) = (0.0d0, 0.0d0)
    do i = 1, n_keep
      k_new(:, i) = dg_frag%k_pw(:, keep_idx(i))
      coef_new(i, :, :) = dg_frag%coef_pw(keep_idx(i), :, :)
      Sfp_new(:, i, :) = S_frag_pw(:, keep_idx(i), :)
      Hfp_new(:, i, :) = H_frag_pw(:, keep_idx(i), :)
      if (allocated(dg_frag%H_mat_pw_diag)) Hpp_new(i, :) = dg_frag%H_mat_pw_diag(keep_idx(i), :)
    end do

    deallocate(dg_frag%k_pw, dg_frag%coef_pw)
    if (allocated(dg_frag%coef_pw_owner)) deallocate(dg_frag%coef_pw_owner)
    if (allocated(dg_frag%coef_pw_full_cache)) deallocate(dg_frag%coef_pw_full_cache)
    dg_frag%coef_pw_full_cache_nstate = 0
    if (allocated(dg_frag%S_mat_frag_pw)) deallocate(dg_frag%S_mat_frag_pw)
    if (allocated(dg_frag%H_mat_frag_pw)) deallocate(dg_frag%H_mat_frag_pw)
    if (allocated(dg_frag%H_mat_pw_diag)) deallocate(dg_frag%H_mat_pw_diag)

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
    allocate(dg_frag%H_mat_frag_pw(n_frag, n_keep, dg_frag%nspin), dg_frag%H_mat_pw_diag(n_keep, dg_frag%nspin))
    dg_frag%S_mat_frag_pw(:, :, :) = Sfp_new(:, :, :)
    dg_frag%H_mat_frag_pw(:, :, :) = Hfp_new(:, :, :)
    dg_frag%H_mat_pw_diag(:, :) = Hpp_new(:, :)

    if (allocated(S_frag_pw)) deallocate(S_frag_pw)
    if (allocated(H_frag_pw)) deallocate(H_frag_pw)
    allocate(S_frag_pw(n_frag, n_keep, dg_frag%nspin), H_frag_pw(n_frag, n_keep, dg_frag%nspin))
    S_frag_pw(:, :, :) = Sfp_new(:, :, :)
    H_frag_pw(:, :, :) = Hfp_new(:, :, :)

    deallocate(k_new, coef_new, Sfp_new, Hfp_new, Hpp_new)
  end subroutine compact_plane_wave_basis

end module rt_dg_plane_wave
