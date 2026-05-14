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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

! DG-Fragment method for real-time TDDFT
! Time evolution in fragment basis coefficient space using Discontinuous Galerkin method
! Velocity gauge formulation: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
!
! ADAPTIVE BASIS UPDATES (NEW FEATURE):
! =====================================
! When strong external fields significantly modify the electronic structure,
! the initial fragment basis may become insufficient. This module implements
! adaptive basis updates to handle such situations:
!
! 1. MONITORING: Tracks ||H_new - H_old|| during RT propagation
! 2. TRIGGER: When change exceeds threshold, updates basis functions
! 3. UPDATE METHODS:
!    a) DC-LCFO CG solver (RECOMMENDED): Solves KS equation with current
!       potentials to get new eigenstates, expanding basis space
!    b) Diagonalization (FALLBACK): Rotates existing basis without expansion
! 4. PROJECTION: Projects current wave functions onto new basis: c'_j = Σ_i S_ji c_i
!
! USAGE:
! ------
! In input file, add:
!   yn_adaptive_basis = 'y'
!   basis_update_threshold = 0.1d0  ! in eV (converted to a.u. internally)
!
! IMPLEMENTATION STATUS:
! ---------------------
! ✓ Hamiltonian change monitoring (fragment-parallel)
! ✓ Overlap matrix calculation between old and new basis
! ✓ Wave function projection to new basis
! ✓ DC-LCFO CG integration framework
! ✓ ppg (pseudopotential grid) properly passed through processing chain
! ✓ Vpsl (local pseudopotential) included in all potential calculations
! ✓ DC-CG method enabled with full pseudopotential support
! ✓ Initial Hamiltonian matrix (H_mat) calculation with kinetic+potential
! ✓ Halo (ghost cell) exchange for phi_frag implemented
! ✓ Fragment boundary treatment via MPI communication
! ✓ System boundary: Periodic boundary conditions properly handled
! ✓ Non-local pseudopotential matrix elements <φ_i|V_NL|φ_j> implemented
! ✓ Input parameter yn_dc_cg_basis_update added for method selection
! ✓ OpenMP parallelization for performance optimization
! ✓ Memory efficiency improvements (reduced allocations)
!
! PERFORMANCE OPTIMIZATIONS:
! -------------------------
! ✓ OpenMP parallel regions:
!   - Non-local PP overlap calculation (collapse 2)
!   - Matrix element calculation (collapse 2)
!   - Time derivative calculation (collapse 3, dynamic scheduling)
!   - Coefficient updates (RK stages, final update)
!   - Basis projection (collapse 3)
!   - Hamiltonian change check (reduction)
! ✓ Memory optimizations:
!   - uVpsi array: Reduced from 4D to 3D (fragment dimension removed)
!   - Reuse uVpsi buffer for each fragment (no reallocation)
!   - Efficient temporary array management
!
! REMAINING TASKS:
! ----------------
! 1. Grid mapping (jxyz_tot) for non-cubic fragment geometries (minor issue)
! 2. Spin handling for matrix elements (currently averages, may need refinement)
! 3. Full runtime testing on HPC system with realistic systems

#include "config.h"
module rt_dg_fragment
  use rt_dg_fragment_types, only: s_dg_fragment_rt, halo_info
  use rt_dg_basis_projection, only: calculate_new_old_basis_overlap, &
                                   stabilize_basis_overlap, &
                                   project_wavefunction_to_new_basis, &
                                   diagonalize_and_update_basis
  use rt_dg_plane_wave, only: init_plane_wave_basis, compute_fragment_pw_overlap, compute_fragment_pw_hamiltonian, build_mixed_hamiltonian, &
                              diagonalize_mixed_basis_pw => diagonalize_mixed_basis
  use xc_hse_ri, only: hse_ri_data_t, init_hse_ri_fragment, &
                       calc_exact_exchange_hse_ri, deallocate_hse_ri_fragment
  use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A, ensure_overlap_prop_available, &
                                calculate_microscopic_current_dg, &
                                build_spatial_A_coupling_matrices, &
                                apply_gradient_to_basis_ops => apply_gradient_to_basis, &
                                rebuild_local_h_block_ids
  implicit none

  private
  public :: init_dg_fragment_rt, tddft_dg_fragment_iteration, finalize_dg_fragment_rt
  public :: calculate_hamiltonian_matrix
  public :: diagnose_density_from_fragments
  public :: get_dg_spin_occ_info
  public :: copy_periodic_global_scalar_to_rank_buffer
  public :: build_total_potential_grid_with_buffered_hartree
  public :: s_dg_fragment_rt, halo_info
  
  ! Types and data structures are defined in rt_dg_fragment_types
  
  ! RI/DF data for HSE exchange (one per local fragment)
  type(hse_ri_data_t), allocatable, save :: hse_ri_data_frag(:)  ! (n_frag_local)
  ! NOTE: dg_hse_ace_enabled, dg_hse_ace_max_age, dg_hse_ace_coef_thresh are now
  !       namelist variables in salmon_global (yn_dg_hse_ace, dg_hse_ace_max_age,
  !       dg_hse_ace_coef_thresh). Only runtime state is kept here.
  ! NOTE: This module is structurally parallel to rt_dg_fragment_soi (SOI variant).
  !       The HSE-ACE cache infrastructure is intentionally identical; any bug fix
  !       here must also be applied to rt_dg_fragment_soi and vice versa.
  logical, save :: dg_hse_ace_initialized = .false.
  real(8), allocatable, save :: hse_ace_vx_cache(:,:,:,:)         ! (nstate_frag,nstate_frag,nspin,n_frag_local)
  complex(8), allocatable, save :: hse_ace_coef_snapshot(:,:,:,:) ! (nstate_frag,nstate_tot,nspin,n_frag_local)
  integer, allocatable, save :: hse_ace_last_rebuild(:,:)         ! (n_frag_local,nspin)
  integer, allocatable, save :: hse_ace_call_count(:,:)           ! (n_frag_local,nspin)
  logical, allocatable, save :: hse_ace_cache_valid(:,:)          ! (n_frag_local,nspin)
  
contains

  subroutine get_dg_spin_occ_info(dg_frag, system, ispin, occ_weight, nocc, state_cap)
    use structures
    use salmon_global, only: nelec, nelec_spin
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    integer, intent(in) :: ispin
    real(8), intent(out) :: occ_weight(:)
    integer, intent(out) :: nocc
    integer, intent(in), optional :: state_cap

    integer :: cap, io, nocc_eff

    occ_weight(:) = 0.0d0
    cap = min(size(occ_weight), dg_frag%nstate_tot)
    if (present(state_cap)) cap = min(cap, max(0, state_cap))

    if (cap <= 0) then
      nocc = 0
      return
    end if

    if (allocated(system%rocc)) then
      if (ispin <= size(system%rocc, 3)) then
        occ_weight(1:min(cap, size(system%rocc, 1))) = &
          max(0.0d0, system%rocc(1:min(cap, size(system%rocc, 1)), 1, ispin))
      end if
    else
      if (dg_frag%nspin == 1) then
        nocc_eff = min(cap, int(nelec / 2.0d0 + 1.0d-12))
        if (nocc_eff > 0) occ_weight(1:nocc_eff) = 2.0d0
      else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
        nocc_eff = min(cap, nelec_spin(ispin))
        if (nocc_eff > 0) occ_weight(1:nocc_eff) = 1.0d0
      else
        nocc_eff = min(cap, int(nelec / 2.0d0 + 1.0d-12))
        if (nocc_eff > 0) occ_weight(1:nocc_eff) = 1.0d0
      end if
    end if

    nocc = 0
    do io = 1, cap
      if (occ_weight(io) > 0.0d0) nocc = io
    end do
  end subroutine get_dg_spin_occ_info

  !=======================================================================
  ! Initialize DG-Fragment RT calculation
  !=======================================================================
  subroutine init_dg_fragment_rt(dg_frag, system, rt, info, lg, mg, ppg)
    use structures
    use communication, only: comm_summation, comm_is_root, comm_create_group, COMM_GROUP_NULL
    use salmon_global, only: num_fragment, num_rgrid_buffer, nstate_frag, time_integrator_dg_fragment, &
                 yn_adaptive_basis, basis_update_threshold, yn_dg_fragment_from_dcdft, &
                 nproc_rgrid, nproc_ob, yn_dg_subspace_diag, dg_subspace_extra_states, yn_spinorbit, &
                 dg_nmat_cap_mode, dg_nmat_cap_fixed, &
                 dg_subspace_pw_vectors, dg_subspace_fallback_cond, nelec, nelec_spin, process_allocation
    use density_matrix_and_energy_plusU_sub, only: PLUS_U_ON
    use filesystem, only: get_filehandle
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rt),             intent(in)    :: rt
    type(s_parallel_info),  intent(in)    :: info
    type(s_rgrid), target,  intent(in)    :: lg, mg
    type(s_pp_grid),        intent(in)    :: ppg
    
    character(32), parameter :: bdir_frag='./data_dcdft/fragments/'
    integer :: i, io, ispin, ia, ip
    integer :: pp_buf(3), parent_rgrid(3)
    integer, parameter :: momentum_stencil_buf = 4
    real(8) :: abs_disp
    logical :: load_from_dcdft
    
    if (comm_is_root(info%id_rko)) then
      write(*,*) "=== Initializing DG-Fragment RT method ==="
    end if
    
    ! Check if we should load DC-LCFO data for DG-Fragment RT
    load_from_dcdft = (trim(yn_dg_fragment_from_dcdft) == 'y')
    
    if (load_from_dcdft .and. comm_is_root(info%id_rko)) then
      write(*,*)
      write(*,*) "  Loading DC-LCFO fragment basis data for DG-Fragment RT"
      write(*,*) "  (yn_dg_fragment_from_dcdft='y')"
      write(*,*)
    end if
    
    ! Set basic parameters
    dg_frag%n_frag = product(num_fragment)
    dg_frag%nspin = system%nspin
    dg_frag%nstate_frag = nstate_frag
    dg_frag%nstate_tot = system%no
    allocate(dg_frag%nocc_spin(dg_frag%nspin))
    dg_frag%nocc_spin(:) = 0
    if (allocated(system%rocc)) then
      do ispin = 1, dg_frag%nspin
        do io = 1, min(dg_frag%nstate_tot, size(system%rocc, 1))
          if (system%rocc(io, 1, ispin) > 0.0d0) dg_frag%nocc_spin(ispin) = io
        end do
      end do
    else
      do ispin = 1, dg_frag%nspin
        if (dg_frag%nspin == 1) then
          dg_frag%nocc_spin(ispin) = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
        else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
          dg_frag%nocc_spin(ispin) = min(dg_frag%nstate_tot, nelec_spin(ispin))
        else
          dg_frag%nocc_spin(ispin) = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
        end if
      end do
    end if

    if (dg_frag%nspin > 1 .and. yn_spinorbit /= 'y') then
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a,i0,a)') "ERROR: non-SOI DG-Fragment RT does not support nspin=", dg_frag%nspin, "."
        write(*,'(1x,a)') "       Use SOI/non-collinear DG-RT for multi-spin fragment propagation."
      end if
      stop "DG-Fragment RT requires SOI/non-collinear path when nspin > 1"
    end if

    dg_frag%icomm = info%icomm_rko
    dg_frag%id = info%id_rko
    dg_frag%isize = info%isize_rko
    dg_frag%icomm_frag = COMM_GROUP_NULL
    dg_frag%id_frag = 0
    dg_frag%isize_frag = 1
    dg_frag%ifrag_group = 0
    dg_frag%nproc_frag = 1
    dg_frag%is_frag_root = .true.
    ! DG fragment subgroups are used as orbital-parallel teams in this branch.
    dg_frag%parallel_mode = 'orbital'
    dg_frag%parallel_mode_orbital = .true.

    if (.not. dg_frag%parallel_mode_orbital .and. trim(dg_frag%parallel_mode) /= 'legacy_realspace') then
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a,a)') "ERROR: invalid dg_fragment_parallel_mode=", trim(dg_frag%parallel_mode)
      end if
      stop "DG-Fragment RT: dg_fragment_parallel_mode must be orbital or legacy_realspace"
    end if

    if (dg_frag%parallel_mode_orbital) then
      dg_frag%nproc_frag = max(1, nproc_ob)
    else
      dg_frag%nproc_frag = product(nproc_rgrid)
    end if
    if (dg_frag%nproc_frag < 1) then
      stop "DG-Fragment RT requires positive fragment subgroup size"
    end if
    ! initialization_rt may temporarily replace the namelist nproc_rgrid by
    ! nproc_rgrid_tot when building the parent DFT grid.  Use info%nprgrid
    ! here because it records the actual parent real-space communicator.
    parent_rgrid(1:3) = info%nprgrid(1:3)

    if (dg_frag%parallel_mode_orbital) then
      if (trim(process_allocation) /= 'orbital_sequential') then
        if (comm_is_root(info%id_rko)) then
          write(*,'(1x,a)') "ERROR: DG orbital fragment parallelism requires process_allocation='orbital_sequential'."
          write(*,'(1x,a,a)') "       Current process_allocation=", trim(process_allocation)
          write(*,'(1x,a)') "       Otherwise each fragment subgroup is laid out as real-space ranks, not orbital ranks."
        end if
        stop "DG-Fragment RT orbital mode requires orbital_sequential rank layout"
      end if
      if (any(parent_rgrid(1:3) /= num_fragment(1:3))) then
        if (comm_is_root(info%id_rko)) then
          write(*,'(1x,a)') "ERROR: DG orbital fragment parallelism requires one parent real-space rank per fragment."
          write(*,'(1x,a,3(i0,1x))') "       parent nproc_rgrid = ", parent_rgrid(1), parent_rgrid(2), parent_rgrid(3)
          write(*,'(1x,a,3(i0,1x))') "       num_fragment = ", num_fragment(1), num_fragment(2), num_fragment(3)
          write(*,'(1x,a,3(i0,1x))') "       namelist nproc_rgrid = ", nproc_rgrid(1), nproc_rgrid(2), nproc_rgrid(3)
          write(*,'(1x,a)') "       Use nproc_ob for fragment-internal orbital parallelism."
        end if
        stop "DG-Fragment RT orbital mode requires parent nproc_rgrid == num_fragment"
      end if
    end if

    if (dg_frag%isize < dg_frag%n_frag) then
      stop "DG-Fragment RT requires np >= n_frag"
    end if
    if (dg_frag%isize /= dg_frag%n_frag * dg_frag%nproc_frag) then
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a,i0,a,i0)') "ERROR: Invalid MPI setup for DG-Fragment RT: np=", dg_frag%isize, ", n_frag=", dg_frag%n_frag
        if (dg_frag%parallel_mode_orbital) then
          write(*,'(1x,a,i0)') "       nproc_ob = ", dg_frag%nproc_frag
        else
          write(*,'(1x,a,i0)') "       product(nproc_rgrid) = ", dg_frag%nproc_frag
        end if
        write(*,'(1x,a)') "       Stage-1 fragment-local MPI requires one subgroup per fragment."
        if (dg_frag%parallel_mode_orbital) then
          write(*,'(1x,a)') "       MPI process count must satisfy np = n_frag * nproc_ob."
        else
          write(*,'(1x,a)') "       MPI process count must satisfy np = n_frag * product(nproc_rgrid)."
        end if
        write(*,'(1x,a)') "       This is a current implementation restriction, not a general DG-RT requirement."
      end if
      stop "DG-Fragment RT stage-1 requires np = n_frag * nproc_frag"
    end if

    dg_frag%id_frag = mod(dg_frag%id, dg_frag%nproc_frag)
    if (dg_frag%parallel_mode_orbital) then
      dg_frag%ifrag_group = fragment_from_rank_address(info%iaddress(1:3), num_fragment)
    else
      dg_frag%ifrag_group = dg_frag%id / dg_frag%nproc_frag + 1
    end if
    dg_frag%isize_frag = dg_frag%nproc_frag
    dg_frag%is_frag_root = (dg_frag%id_frag == 0)
    dg_frag%icomm_frag = comm_create_group(dg_frag%icomm, dg_frag%ifrag_group - 1, dg_frag%id_frag)
    
    ! Check DFT+U status
    dg_frag%use_plusu = PLUS_U_ON
    if (dg_frag%use_plusu .and. comm_is_root(info%id_rko)) then
      write(*,*)
      write(*,*) "  *** DFT+U is ENABLED for DG-Fragment RT ***"
      write(*,*) "  +U update in DG coefficient-space loop is not fully wired yet"
      write(*,*) "  Current status: framework available, full RT coupling pending"
      write(*,*)
    end if
    
    ! Stage-1 two-level MPI layout: one fragment per subgroup.
    dg_frag%ifrag_start = dg_frag%ifrag_group
    dg_frag%ifrag_end = dg_frag%ifrag_group
    
    if (comm_is_root(info%id_rko)) then
      write(*,'(1x,a,i0,a,i0)') "  MPI parallelization: ", dg_frag%isize, " processes"
      write(*,'(1x,a,a)') "  DG fragment parallel mode: ", trim(dg_frag%parallel_mode)
      write(*,'(1x,a,i0)') "  MPI ranks per fragment subgroup: ", dg_frag%nproc_frag
      write(*,'(1x,a)') "  Fragment distribution across MPI ranks:"
    end if
    write(*,'(1x,a,i4,a,i4,a,i2,a,l1)') "    Rank", dg_frag%id, ": fragment ", &
                                    dg_frag%ifrag_start, ", subgroup rank ", dg_frag%id_frag, &
                                    ", root=", dg_frag%is_frag_root
    
    ! Set time integrator
    select case(trim(time_integrator_dg_fragment))
    case('ssprk3')
      dg_frag%time_integrator = 1
    case('aetrs')
      dg_frag%time_integrator = 2
    case('rk4')
      dg_frag%time_integrator = 3
    case default
      dg_frag%time_integrator = 3  ! default: RK4
    end select
    
    ! Store fragment geometry information for halo exchange
    dg_frag%num_fragment(1:3) = num_fragment(1:3)
    dg_frag%lgnum_total(1:3) = lg%num(1:3)
    pp_buf(1:3) = 0
    if (allocated(ppg%mps) .and. allocated(ppg%rxyz)) then
      do ia = 1, size(ppg%mps)
        do ip = 1, ppg%mps(ia)
          abs_disp = abs(ppg%rxyz(1, ip, ia))
          if (system%hgs(1) > 0d0) pp_buf(1) = max(pp_buf(1), ceiling(abs_disp / system%hgs(1)))
          abs_disp = abs(ppg%rxyz(2, ip, ia))
          if (system%hgs(2) > 0d0) pp_buf(2) = max(pp_buf(2), ceiling(abs_disp / system%hgs(2)))
          abs_disp = abs(ppg%rxyz(3, ip, ia))
          if (system%hgs(3) > 0d0) pp_buf(3) = max(pp_buf(3), ceiling(abs_disp / system%hgs(3)))
        end do
      end do
    end if
    dg_frag%nxyz_buffer(1:3) = max(num_rgrid_buffer(1:3), momentum_stencil_buf)
    dg_frag%nxyz_buffer(1:3) = max(dg_frag%nxyz_buffer(1:3), pp_buf(1:3))
    if (comm_is_root(info%id_rko)) then
      write(*,'(1x,a,3i5)') "  nxyz_buffer (runtime): ", dg_frag%nxyz_buffer(1:3)
      write(*,'(1x,a,3i5)') "  nxyz_buffer (pp-derived): ", pp_buf(1:3)
    end if
    
    ! Allocate fragment geometry arrays
    allocate(dg_frag%nxyz_domain(3, dg_frag%n_frag))
    allocate(dg_frag%ixyz_frag(3, dg_frag%n_frag))
    allocate(dg_frag%id_array(dg_frag%n_frag))
    do i = 1, dg_frag%n_frag
      dg_frag%id_array(i) = get_fragment_group_root_rank(i, dg_frag%nproc_frag)
    end do
    dg_frag%lg => lg
    dg_frag%mg => mg
    dg_frag%hgs = system%hgs
    
    ! Read fragment basis data from DC-LCFO calculation if requested
    if (load_from_dcdft) then
      call read_fragment_basis_data(dg_frag, bdir_frag)
      call read_dg_occupation_seed(dg_frag)
    else
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a)') "WARNING: yn_dg_fragment_from_dcdft='n'"
        write(*,'(1x,a)') "Fragment basis data will NOT be loaded from DC-LCFO"
        write(*,'(1x,a)') "You must provide initial data via yn_restart='y'"
      end if
      ! Set minimal required fields for non-DCDFT initialization
      dg_frag%has_real_space_basis = .false.
      if (.not. allocated(dg_frag%n_mat)) then
        allocate(dg_frag%n_mat(dg_frag%nspin))
        dg_frag%n_mat = dg_frag%nstate_frag
      end if
      dg_frag%n_mat_max = maxval(dg_frag%n_mat)
    end if
    
    ! Halo communication is explicitly disabled by design in this branch.
    dg_frag%n_halo = 0
    dg_frag%has_halo_exchange = .false.
    if (dg_frag%has_real_space_basis) then
      call rebuild_coef_owner_map(dg_frag, "post-halo-disabled")
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a)') "Fragment halo exchange disabled: local PBC-buffered stencil only"
      end if
    end if
    
    ! Note: has_real_space_basis flag is set in read_fragment_basis_data()
    ! after successfully loading basis functions from DC-LCFO data
    
    ! Coefficient arrays will be allocated by read_fragment_basis_data()
    ! with proper n_mat_max dimensions for global basis compatibility
    
    ! Allocate only ESP array (independent of basis function count)
    allocate(dg_frag%esp(dg_frag%nstate_tot, dg_frag%nspin))
    
    ! Initialize RK coefficients
    call init_rk_coefficients(dg_frag)
    
    ! Store grid pointers
    if (dg_frag%has_real_space_basis) call build_density_grid_owner_maps(dg_frag)
    
    ! Allocate density and potential arrays
    allocate(dg_frag%rho_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(dg_frag%rho_s_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), dg_frag%nspin))
    allocate(dg_frag%Vh_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(dg_frag%Vxc_frag(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), dg_frag%nspin))
    
    ! Initialize to zero
    dg_frag%rho_frag = 0.0d0
    dg_frag%rho_s_frag = 0.0d0
    dg_frag%Vh_frag = 0.0d0
    dg_frag%Vxc_frag = 0.0d0
    
    ! Initialize adaptive basis update functionality
    ! (Experimental feature for strong-field regimes where mean-field changes significantly)
    ! Initialize adaptive basis parameters from input file
    dg_frag%yn_adaptive_basis = (yn_adaptive_basis == 'y')
    dg_frag%basis_update_threshold = basis_update_threshold
    dg_frag%use_subspace_diag = (yn_dg_subspace_diag == 'y')
    dg_frag%subspace_extra_states = max(0, dg_subspace_extra_states)
    dg_frag%subspace_pw_vectors = max(0, dg_subspace_pw_vectors)
    dg_frag%subspace_fallback_cond = max(1.0d0, dg_subspace_fallback_cond)
    dg_frag%last_subspace_dim = 0
    dg_frag%subspace_fallback_count = 0
    dg_frag%has_nl_cache = .false.
    dg_frag%Ac_nl_cache = 0.0d0
    dg_frag%Ac_nl_cache_tol = 1.0d-12
    dg_frag%hamiltonian_change_norm = 0.0d0
    dg_frag%nbasis_update_count = 0
    dg_frag%last_basis_update_step = 0
    
    ! Allocate Hamiltonian tracking arrays (always allocated for monitoring)
    allocate(dg_frag%H_mat_old(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin))
    dg_frag%H_mat_old = (0.0d0, 0.0d0)
    
    ! Allocate overlap matrix only if adaptive basis is enabled
    if (dg_frag%yn_adaptive_basis) then
      allocate(dg_frag%eigenvalues(dg_frag%nstate_frag, dg_frag%nspin))
      dg_frag%eigenvalues = 0.0d0
      allocate(dg_frag%basis_overlap(dg_frag%nstate_frag, dg_frag%nstate_frag, dg_frag%nspin))
      dg_frag%basis_overlap = 0.0d0
      
      if (comm_is_root(info%id_rko)) then
        write(*,*)
        write(*,*) "=== Adaptive Basis Updates ENABLED ==="
        write(*,'(1x,a,f8.4,a)') "  Hamiltonian change threshold: ", &
                                 dg_frag%basis_update_threshold, " a.u."
        write(*,*) "  Basis will be recalculated when mean field changes significantly"
        write(*,*) "  DC-LCFO returns real basis (no gauge rotation needed)"
        write(*,*)
      end if
    end if
    if (comm_is_root(info%id_rko)) then
      write(*,'(1x,a,l1)') "  DG subspace diagonalization: ", dg_frag%use_subspace_diag
      write(*,'(1x,a,i0)') "  DG subspace extra states: ", dg_frag%subspace_extra_states
      write(*,'(1x,a,i0)') "  DG subspace PW vectors: ", dg_frag%subspace_pw_vectors
      write(*,'(1x,a,1pe11.3)') "  DG subspace fallback cond: ", dg_frag%subspace_fallback_cond
    end if
    
    ! NOTE: Dipole matrix NOT needed for velocity gauge formalism
    ! In velocity gauge, current is calculated from momentum operator <p> = <ψ|-i∇|ψ>
    
    ! Note: Momentum matrix will be calculated in calculate_hamiltonian_matrix
    !       when stencil and potentials are available
    
    ! Defer large global matrices until Hamiltonian construction to reduce init peak memory.
    if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
    
    ! H_mat and H_mat_kinetic are allocated in calculate_hamiltonian_matrix
    ! when they are first needed in RT propagation.
    
    if (comm_is_root(info%id_rko)) then
      write(*,*) "  Kinetic + ion Hamiltonian saved for incremental updates"
    end if
    
    ! Initialize +U if enabled (framework only, actual implementation requires ppg)
    ! Note: init_plusu_framework not implemented yet
    ! if (dg_frag%use_plusu) then
    !   call init_plusu_framework(dg_frag)
    ! end if
    
    ! Initialize plane wave basis if enabled
    call init_plane_wave_basis(dg_frag, system, lg, info)
    if (allocated(dg_frag%coef_pw_owner)) deallocate(dg_frag%coef_pw_owner)
    dg_frag%owned_coef_pw_start = 0
    dg_frag%owned_coef_pw_end = -1
    if (dg_frag%n_plane_waves > 0) then
      allocate(dg_frag%coef_pw_owner(dg_frag%n_plane_waves))
      do i = 1, dg_frag%n_plane_waves
        if (dg_frag%parallel_mode_orbital) then
          ! PW coefficients are fragment-local, just like fragment orbital
          ! coefficients.  Split their rows inside each fragment subgroup,
          ! not once over the whole MPI world.
          dg_frag%coef_pw_owner(i) = get_subgroup_block_owner_rank( &
            dg_frag%id_array(dg_frag%ifrag_group), dg_frag%isize_frag, i, dg_frag%n_plane_waves)
        else
          dg_frag%coef_pw_owner(i) = min(dg_frag%isize - 1, ((i - 1) * dg_frag%isize) / dg_frag%n_plane_waves)
        end if
      end do
      do i = 1, dg_frag%n_plane_waves
        if (dg_frag%coef_pw_owner(i) /= dg_frag%id) cycle
        if (dg_frag%owned_coef_pw_start == 0) dg_frag%owned_coef_pw_start = i
        dg_frag%owned_coef_pw_end = i
      end do
    end if
    
    ! Initialize RI/DF approximation for HSE if enabled
    call init_hse_ri_data(dg_frag, system, info)

    dg_frag%current(:) = 0.0d0
    dg_frag%current_para(:) = 0.0d0
    dg_frag%current_dia(:) = 0.0d0
    dg_frag%current_total(:) = 0.0d0
    dg_frag%elec_num_raw_t0 = 0.0d0
    dg_frag%elec_num_scaled_t0 = 0.0d0
    dg_frag%elec_num_baseline_ready = .false.
    dg_frag%rho_drift_indicator = 0.0d0
    dg_frag%rho_ff_elec = 0.0d0
    dg_frag%rho_fp_elec = 0.0d0
    dg_frag%rho_pp_elec = 0.0d0
    dg_frag%rho_owned_elec = 0.0d0
    dg_frag%rho_imported_elec = 0.0d0
    dg_frag%rho_local_all_elec = 0.0d0
    dg_frag%rho_global_raw_elec = 0.0d0
    dg_frag%rho_global_scaled_elec = 0.0d0
    dg_frag%rho_contract_residual_elec = 0.0d0
    dg_frag%coef_occ_norm0 = 0.0d0
    dg_frag%coef_occ_norm = 0.0d0
    dg_frag%coef_occ_norm_drift = 0.0d0
    dg_frag%csc_occ_identity_norm = 0.0d0
    dg_frag%csc_occ_identity_max = 0.0d0
    dg_frag%occvirt_leakage = 0.0d0
    dg_frag%occvirt_top_occ = 0
    dg_frag%occvirt_top_virt = 0
    dg_frag%occvirt_top_abs2 = 0.0d0
    dg_frag%jpara_top_occ_state = 0
    dg_frag%jpara_top_occ_value = 0.0d0
    dg_frag%selfexc_track_states(:) = (/ 99, 100, 129 /)
    dg_frag%selfexc_track_norm(:) = 0.0d0
    dg_frag%selfexc_csc_step_delta = 0.0d0
    dg_frag%selfexc_leak100_129_step_delta = 0.0d0
    dg_frag%selfexc_csc_stage_pre_mixed = 0.0d0
    dg_frag%selfexc_csc_stage_post_overlap = 0.0d0
    dg_frag%selfexc_csc_stage_post_raw = 0.0d0
    dg_frag%selfexc_leak100_129_pre_mixed = 0.0d0
    dg_frag%selfexc_leak100_129_post_overlap = 0.0d0
    dg_frag%selfexc_leak100_129_post_raw = 0.0d0
    dg_frag%selfexc_leak100_129_current = 0.0d0
    dg_frag%selfexc_csc_prev_step = 0.0d0
    dg_frag%selfexc_leak100_129_prev_step = 0.0d0
    dg_frag%selfexc_prev_step_itt = -1
    dg_frag%coef_ref_ready = .false.
    
    if (comm_is_root(info%id_rko)) then
      write(*,*) "DG-Fragment RT initialization complete"
      write(*,'(1x,a,i0)') "  Number of fragments: ", dg_frag%n_frag
      write(*,'(1x,a,i0)') "  States per fragment: ", dg_frag%nstate_frag
      write(*,'(1x,a,i0)') "  Total states: ", dg_frag%nstate_tot
      write(*,'(1x,a,a)')  "  Time integrator: ", trim(time_integrator_dg_fragment)
      if (dg_frag%use_plane_wave_basis) then
        write(*,'(1x,a,i0)') "  Plane waves added: ", dg_frag%n_plane_waves
        write(*,'(1x,a,f8.4,a)') "  PW cutoff: ", dg_frag%k_cutoff_pw, " a.u.^-1"
      end if
      if (dg_frag%use_hse_ri) then
        write(*,'(1x,a)') "  HSE RI/DF approximation: ENABLED"
      end if
    end if
    
  end subroutine init_dg_fragment_rt

  subroutine build_density_grid_owner_maps(dg_frag)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ifrag, i_local, ix, iy, iz, ixg, iyg, izg, owner_rank, source_rank, source_root_rank, i
    integer :: subgroup_target_rank
    integer :: ifrag_count
    integer :: nx_max, ny_max, nz_max
    integer, allocatable :: recv_count(:), recv_cursor(:)
    integer, allocatable :: id_tmp(:)

    if (allocated(dg_frag%density_owner_map)) deallocate(dg_frag%density_owner_map)
    if (allocated(dg_frag%density_primary_local_map)) deallocate(dg_frag%density_primary_local_map)
    if (allocated(dg_frag%density_ixg_map)) deallocate(dg_frag%density_ixg_map)
    if (allocated(dg_frag%density_iyg_map)) deallocate(dg_frag%density_iyg_map)
    if (allocated(dg_frag%density_izg_map)) deallocate(dg_frag%density_izg_map)
    if (allocated(dg_frag%density_send_count)) deallocate(dg_frag%density_send_count)
    if (allocated(dg_frag%density_send_slot_map)) deallocate(dg_frag%density_send_slot_map)
    if (allocated(dg_frag%density_subgroup_send_count)) deallocate(dg_frag%density_subgroup_send_count)
    if (allocated(dg_frag%density_subgroup_send_slot_map)) deallocate(dg_frag%density_subgroup_send_slot_map)
    if (allocated(dg_frag%density_grid_points)) deallocate(dg_frag%density_grid_points)
    if (allocated(dg_frag%density_grid_point_count)) deallocate(dg_frag%density_grid_point_count)
    if (allocated(dg_frag%density_grid_bx)) deallocate(dg_frag%density_grid_bx)
    if (allocated(dg_frag%density_grid_by)) deallocate(dg_frag%density_grid_by)
    if (allocated(dg_frag%density_grid_bz)) deallocate(dg_frag%density_grid_bz)
    dg_frag%density_rhobf_box_cache_valid = .false.
    if (allocated(dg_frag%density_subgroup_self_ixg)) deallocate(dg_frag%density_subgroup_self_ixg)
    if (allocated(dg_frag%density_subgroup_self_iyg)) deallocate(dg_frag%density_subgroup_self_iyg)
    if (allocated(dg_frag%density_subgroup_self_izg)) deallocate(dg_frag%density_subgroup_self_izg)
    if (allocated(dg_frag%density_block_nblocks)) deallocate(dg_frag%density_block_nblocks)
    if (allocated(dg_frag%density_block_first_offset)) deallocate(dg_frag%density_block_first_offset)
    if (allocated(dg_frag%density_block_step)) deallocate(dg_frag%density_block_step)
    if (allocated(dg_frag%current_valid_grid_count)) deallocate(dg_frag%current_valid_grid_count)
    if (allocated(dg_frag%current_valid_ix)) deallocate(dg_frag%current_valid_ix)
    if (allocated(dg_frag%current_valid_iy)) deallocate(dg_frag%current_valid_iy)
    if (allocated(dg_frag%current_valid_iz)) deallocate(dg_frag%current_valid_iz)
    if (allocated(dg_frag%current_valid_ixg)) deallocate(dg_frag%current_valid_ixg)
    if (allocated(dg_frag%current_valid_iyg)) deallocate(dg_frag%current_valid_iyg)
    if (allocated(dg_frag%current_valid_izg)) deallocate(dg_frag%current_valid_izg)
    if (allocated(dg_frag%runtime_neighbor_pair_cache)) deallocate(dg_frag%runtime_neighbor_pair_cache)
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) deallocate(dg_frag%momentum_neighbor_pair_cache)
    if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
    if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
    dg_frag%density_phi_block_size = 0
    dg_frag%density_phi_block_cache_valid = .false.
    if (allocated(dg_frag%density_phase_block_cache)) deallocate(dg_frag%density_phase_block_cache)
    dg_frag%density_phase_block_size = 0
    dg_frag%density_phase_block_npw = 0
    dg_frag%density_phase_block_cache_valid = .false.
    if (allocated(dg_frag%density_recv_map)) then
      do i = lbound(dg_frag%density_recv_map, 1), ubound(dg_frag%density_recv_map, 1)
        if (allocated(dg_frag%density_recv_map(i)%ixg)) deallocate(dg_frag%density_recv_map(i)%ixg)
        if (allocated(dg_frag%density_recv_map(i)%iyg)) deallocate(dg_frag%density_recv_map(i)%iyg)
        if (allocated(dg_frag%density_recv_map(i)%izg)) deallocate(dg_frag%density_recv_map(i)%izg)
        if (allocated(dg_frag%density_recv_map(i)%bx)) deallocate(dg_frag%density_recv_map(i)%bx)
        if (allocated(dg_frag%density_recv_map(i)%by)) deallocate(dg_frag%density_recv_map(i)%by)
        if (allocated(dg_frag%density_recv_map(i)%bz)) deallocate(dg_frag%density_recv_map(i)%bz)
      end do
      deallocate(dg_frag%density_recv_map)
    end if
    if (allocated(dg_frag%density_send_count)) deallocate(dg_frag%density_send_count)
    if (allocated(dg_frag%density_send_slot_map)) deallocate(dg_frag%density_send_slot_map)
    if (allocated(dg_frag%density_subgroup_send_count)) deallocate(dg_frag%density_subgroup_send_count)
    if (allocated(dg_frag%density_subgroup_send_slot_map)) deallocate(dg_frag%density_subgroup_send_slot_map)
    if (allocated(dg_frag%density_grid_points)) deallocate(dg_frag%density_grid_points)
    if (allocated(dg_frag%density_grid_point_count)) deallocate(dg_frag%density_grid_point_count)
    if (allocated(dg_frag%density_grid_bx)) deallocate(dg_frag%density_grid_bx)
    if (allocated(dg_frag%density_grid_by)) deallocate(dg_frag%density_grid_by)
    if (allocated(dg_frag%density_grid_bz)) deallocate(dg_frag%density_grid_bz)
    dg_frag%density_rhobf_box_cache_valid = .false.
    if (allocated(dg_frag%density_subgroup_self_ixg)) deallocate(dg_frag%density_subgroup_self_ixg)
    if (allocated(dg_frag%density_subgroup_self_iyg)) deallocate(dg_frag%density_subgroup_self_iyg)
    if (allocated(dg_frag%density_subgroup_self_izg)) deallocate(dg_frag%density_subgroup_self_izg)
    if (allocated(dg_frag%density_block_nblocks)) deallocate(dg_frag%density_block_nblocks)
    if (allocated(dg_frag%density_block_first_offset)) deallocate(dg_frag%density_block_first_offset)
    if (allocated(dg_frag%density_block_step)) deallocate(dg_frag%density_block_step)
    if (allocated(dg_frag%current_valid_grid_count)) deallocate(dg_frag%current_valid_grid_count)
    if (allocated(dg_frag%current_valid_ix)) deallocate(dg_frag%current_valid_ix)
    if (allocated(dg_frag%current_valid_iy)) deallocate(dg_frag%current_valid_iy)
    if (allocated(dg_frag%current_valid_iz)) deallocate(dg_frag%current_valid_iz)
    if (allocated(dg_frag%current_valid_ixg)) deallocate(dg_frag%current_valid_ixg)
    if (allocated(dg_frag%current_valid_iyg)) deallocate(dg_frag%current_valid_iyg)
    if (allocated(dg_frag%current_valid_izg)) deallocate(dg_frag%current_valid_izg)
    if (allocated(dg_frag%runtime_neighbor_pair_cache)) deallocate(dg_frag%runtime_neighbor_pair_cache)
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) deallocate(dg_frag%momentum_neighbor_pair_cache)
    if (allocated(dg_frag%density_recv_map)) then
      do ifrag = lbound(dg_frag%density_recv_map, 1), ubound(dg_frag%density_recv_map, 1)
        if (allocated(dg_frag%density_recv_map(ifrag)%ixg)) deallocate(dg_frag%density_recv_map(ifrag)%ixg)
        if (allocated(dg_frag%density_recv_map(ifrag)%iyg)) deallocate(dg_frag%density_recv_map(ifrag)%iyg)
        if (allocated(dg_frag%density_recv_map(ifrag)%izg)) deallocate(dg_frag%density_recv_map(ifrag)%izg)
        if (allocated(dg_frag%density_recv_map(ifrag)%bx)) deallocate(dg_frag%density_recv_map(ifrag)%bx)
        if (allocated(dg_frag%density_recv_map(ifrag)%by)) deallocate(dg_frag%density_recv_map(ifrag)%by)
        if (allocated(dg_frag%density_recv_map(ifrag)%bz)) deallocate(dg_frag%density_recv_map(ifrag)%bz)
      end do
      deallocate(dg_frag%density_recv_map)
    end if
    if (.not. associated(dg_frag%mg)) return
    if (.not. allocated(dg_frag%nxyz_domain)) return
    if (.not. allocated(dg_frag%ixyz_frag)) return
    if (dg_frag%ifrag_end < dg_frag%ifrag_start) return

    call build_fragment_global_boxes(dg_frag)

    ! Rebuild id_array (fragment -> root rank) via comm_summation so that all
    ! ranks have a consistent mapping regardless of halo-exchange init state.
    if (allocated(dg_frag%id_array)) then
      allocate(id_tmp(dg_frag%n_frag))
      id_tmp = 0
      if (dg_frag%is_frag_root) then
        do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
          id_tmp(ifrag) = dg_frag%id + 1  ! +1 to distinguish from unset (0)
        end do
      end if
      call comm_summation(id_tmp, dg_frag%id_array, dg_frag%n_frag, dg_frag%icomm)
      dg_frag%id_array = dg_frag%id_array - 1  ! Convert back to 0-based MPI rank
      deallocate(id_tmp)
    end if

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    nx_max = max(1, maxval(dg_frag%nxyz_domain(1, dg_frag%ifrag_start:dg_frag%ifrag_end)))
    ny_max = max(1, maxval(dg_frag%nxyz_domain(2, dg_frag%ifrag_start:dg_frag%ifrag_end)))
    nz_max = max(1, maxval(dg_frag%nxyz_domain(3, dg_frag%ifrag_start:dg_frag%ifrag_end)))

    ! The density maps are the contract between fragment-local coordinates and
    ! the parent real-space grid: every local fragment point gets a global grid
    ! index, an owner rank, and optional send/receive slots.
    allocate(dg_frag%density_owner_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_primary_local_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_ixg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_iyg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_izg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_send_count(0:dg_frag%isize-1))
    allocate(dg_frag%density_send_slot_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_subgroup_send_count(0:dg_frag%isize_frag-1))
    allocate(dg_frag%density_subgroup_send_slot_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_recv_map(0:dg_frag%isize-1))
    allocate(dg_frag%density_grid_points(nx_max * ny_max * nz_max, ifrag_count))
    allocate(dg_frag%density_grid_point_count(ifrag_count))
    dg_frag%density_owner_map = dg_frag%id
    dg_frag%density_primary_local_map = .false.
    dg_frag%density_ixg_map = 1
    dg_frag%density_iyg_map = 1
    dg_frag%density_izg_map = 1
    dg_frag%density_send_count = 0
    dg_frag%density_send_slot_map = 0
    dg_frag%density_subgroup_send_count = 0
    dg_frag%density_subgroup_send_slot_map = 0
    dg_frag%density_grid_point_count = 0

    ! First pass: build the maps for fragments local to this rank.  In orbital
    ! mode the fragment root owns the replicated density contribution, so the
    ! source rank and density owner must agree.
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        izg = wrap_global_grid_index(dg_frag%frag_core_lo(3, ifrag) + iz - 1, dg_frag%lgnum_total(3))
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          iyg = wrap_global_grid_index(dg_frag%frag_core_lo(2, ifrag) + iy - 1, dg_frag%lgnum_total(2))
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            ixg = wrap_global_grid_index(dg_frag%frag_core_lo(1, ifrag) + ix - 1, dg_frag%lgnum_total(1))
            dg_frag%density_ixg_map(ix, iy, iz, i_local) = ixg
            dg_frag%density_iyg_map(ix, iy, iz, i_local) = iyg
            dg_frag%density_izg_map(ix, iy, iz, i_local) = izg
            dg_frag%density_primary_local_map(ix, iy, iz, i_local) = .true.
            source_rank = dg_frag%id_array(ifrag)
            subgroup_target_rank = source_rank - dg_frag%id_array(ifrag)
            if (dg_frag%density_primary_local_map(ix, iy, iz, i_local)) then
              source_rank = get_fragment_grid_sender_rank(dg_frag%id_array(ifrag), dg_frag%nxyz_domain(:, ifrag), &
                                                           ix, iy, iz, dg_frag%parallel_mode_orbital)
              subgroup_target_rank = source_rank - dg_frag%id_array(ifrag)
              if (subgroup_target_rank < 0 .or. subgroup_target_rank > dg_frag%isize_frag - 1) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
                     "[FATAL] density subgroup target out of range: rank=", dg_frag%id, &
                     " ifrag=", ifrag, " ix=", ix, " iy=", iy, " iz=", iz
                flush(6)
                stop "DG-Fragment RT: density subgroup target out of range"
              end if
              if (subgroup_target_rank /= dg_frag%id_frag) then
                dg_frag%density_subgroup_send_count(subgroup_target_rank) = &
                  dg_frag%density_subgroup_send_count(subgroup_target_rank) + 1
                dg_frag%density_subgroup_send_slot_map(ix, iy, iz, i_local) = &
                  dg_frag%density_subgroup_send_count(subgroup_target_rank)
              end if
            end if
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg, dg_frag%id_array(ifrag))
            dg_frag%density_owner_map(ix, iy, iz, i_local) = owner_rank
            source_rank = dg_frag%id_array(ifrag)
            ! Density computation may be orbital-parallel, but the final rho
            ! field is still owned by the parent real-space grid.  The fragment
            ! root emits each primary point once, either locally or through the
            ! packed send map below.
            if (dg_frag%density_primary_local_map(ix, iy, iz, i_local) .and. owner_rank /= source_rank .and. &
                dg_frag%is_frag_root) then
              dg_frag%density_send_count(owner_rank) = dg_frag%density_send_count(owner_rank) + 1
              dg_frag%density_send_slot_map(ix, iy, iz, i_local) = dg_frag%density_send_count(owner_rank)
            end if
            dg_frag%density_grid_point_count(i_local) = dg_frag%density_grid_point_count(i_local) + 1
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%ix = ix
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%iy = iy
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%iz = iz
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%ixg = ixg
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%iyg = iyg
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%izg = izg
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%owner_rank = owner_rank
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%is_primary = &
              dg_frag%density_primary_local_map(ix, iy, iz, i_local)
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%send_slot = &
              dg_frag%density_send_slot_map(ix, iy, iz, i_local)
            dg_frag%density_grid_points(dg_frag%density_grid_point_count(i_local), i_local)%subgroup_send_slot = &
              dg_frag%density_subgroup_send_slot_map(ix, iy, iz, i_local)
          end do
        end do
      end do
    end do

    allocate(recv_count(0:dg_frag%isize-1), recv_cursor(0:dg_frag%isize-1))
    recv_count = 0
    ! Second pass: count remote density points that this rank receives from
    ! other fragment roots.  Counts are needed before allocating recv buffers.
    do ifrag = 1, dg_frag%n_frag
      source_root_rank = dg_frag%id_array(ifrag)
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        izg = wrap_global_grid_index(dg_frag%frag_core_lo(3, ifrag) + iz - 1, dg_frag%lgnum_total(3))
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          iyg = wrap_global_grid_index(dg_frag%frag_core_lo(2, ifrag) + iy - 1, dg_frag%lgnum_total(2))
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            ixg = wrap_global_grid_index(dg_frag%frag_core_lo(1, ifrag) + ix - 1, dg_frag%lgnum_total(1))
            source_rank = source_root_rank
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg, source_root_rank)
            if (source_rank == dg_frag%id) then
              cycle
            end if
            if (owner_rank == dg_frag%id) then
              recv_count(source_rank) = recv_count(source_rank) + 1
            end if
          end do
        end do
      end do
    end do

    do source_rank = 0, dg_frag%isize - 1
      dg_frag%density_recv_map(source_rank)%npts = recv_count(source_rank)
      if (recv_count(source_rank) <= 0) cycle
      allocate(dg_frag%density_recv_map(source_rank)%ixg(recv_count(source_rank)))
      allocate(dg_frag%density_recv_map(source_rank)%iyg(recv_count(source_rank)))
      allocate(dg_frag%density_recv_map(source_rank)%izg(recv_count(source_rank)))
    end do

    recv_cursor = 0
    ! Third pass: materialize the receive coordinate lists using the counts
    ! above.  The reconstruction path later uses these lists directly.
    do ifrag = 1, dg_frag%n_frag
      source_root_rank = dg_frag%id_array(ifrag)
      do iz = 1, dg_frag%nxyz_domain(3, ifrag)
        izg = wrap_global_grid_index(dg_frag%frag_core_lo(3, ifrag) + iz - 1, dg_frag%lgnum_total(3))
        do iy = 1, dg_frag%nxyz_domain(2, ifrag)
          iyg = wrap_global_grid_index(dg_frag%frag_core_lo(2, ifrag) + iy - 1, dg_frag%lgnum_total(2))
          do ix = 1, dg_frag%nxyz_domain(1, ifrag)
            ixg = wrap_global_grid_index(dg_frag%frag_core_lo(1, ifrag) + ix - 1, dg_frag%lgnum_total(1))
            source_rank = source_root_rank
            if (source_rank == dg_frag%id) cycle
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg, source_root_rank)
            if (owner_rank == dg_frag%id) then
              recv_cursor(source_rank) = recv_cursor(source_rank) + 1
              dg_frag%density_recv_map(source_rank)%ixg(recv_cursor(source_rank)) = ixg
              dg_frag%density_recv_map(source_rank)%iyg(recv_cursor(source_rank)) = iyg
              dg_frag%density_recv_map(source_rank)%izg(recv_cursor(source_rank)) = izg
            end if
          end do
        end do
      end do
    end do

    deallocate(recv_count, recv_cursor)

  end subroutine build_density_grid_owner_maps

  logical function is_density_primary_fragment(dg_frag, ifrag_ref, ixg, iyg, izg) result(is_primary)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_ref, ixg, iyg, izg
    is_primary = is_density_core_point(dg_frag, ifrag_ref, ixg, iyg, izg)
  end function is_density_primary_fragment

  integer function wrap_global_grid_index(ig_raw, ngrid) result(ig)
    implicit none
    integer, intent(in) :: ig_raw, ngrid

    ig = modulo(ig_raw - 1, ngrid) + 1
  end function wrap_global_grid_index

  integer function get_fragment_grid_sender_rank(root_rank, ndom, ix, iy, iz, parallel_mode_orbital) result(sender_rank)
    use salmon_global, only: nproc_rgrid
    implicit none
    integer, intent(in) :: root_rank, ndom(3), ix, iy, iz
    logical, intent(in) :: parallel_mode_orbital
    integer :: ipx, ipy, ipz, coords(3), nsize

    if (parallel_mode_orbital) then
      sender_rank = root_rank
      return
    end if

    ipx = max(1, nproc_rgrid(1))
    ipy = max(1, nproc_rgrid(2))
    ipz = max(1, nproc_rgrid(3))

    nsize = max(1, (ndom(1) + ipx - 1) / ipx)
    coords(1) = min(ipx - 1, max(0, (ix - 1) / nsize))
    nsize = max(1, (ndom(2) + ipy - 1) / ipy)
    coords(2) = min(ipy - 1, max(0, (iy - 1) / nsize))
    nsize = max(1, (ndom(3) + ipz - 1) / ipz)
    coords(3) = min(ipz - 1, max(0, (iz - 1) / nsize))

    sender_rank = root_rank + coords(1) + ipx * (coords(2) + ipy * coords(3))
  end function get_fragment_grid_sender_rank

  subroutine build_fragment_global_boxes(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag

    if (.not. allocated(dg_frag%nxyz_domain) .or. .not. allocated(dg_frag%ixyz_frag)) return

    if (allocated(dg_frag%frag_core_lo)) deallocate(dg_frag%frag_core_lo)
    if (allocated(dg_frag%frag_core_hi)) deallocate(dg_frag%frag_core_hi)
    if (allocated(dg_frag%frag_buf_lo)) deallocate(dg_frag%frag_buf_lo)
    if (allocated(dg_frag%frag_buf_hi)) deallocate(dg_frag%frag_buf_hi)
    allocate(dg_frag%frag_core_lo(3, dg_frag%n_frag), dg_frag%frag_core_hi(3, dg_frag%n_frag))
    allocate(dg_frag%frag_buf_lo(3, dg_frag%n_frag), dg_frag%frag_buf_hi(3, dg_frag%n_frag))

    do ifrag = 1, dg_frag%n_frag
      dg_frag%frag_core_lo(:, ifrag) = dg_frag%ixyz_frag(:, ifrag)
      dg_frag%frag_core_hi(:, ifrag) = dg_frag%ixyz_frag(:, ifrag) + dg_frag%nxyz_domain(:, ifrag) - 1
      dg_frag%frag_buf_lo(:, ifrag) = dg_frag%frag_core_lo(:, ifrag) - dg_frag%nxyz_buffer(:)
      dg_frag%frag_buf_hi(:, ifrag) = dg_frag%frag_core_hi(:, ifrag) + dg_frag%nxyz_buffer(:)
    end do

    if (dg_frag%parallel_mode_orbital .and. dg_frag%ifrag_group >= 1 .and. dg_frag%ifrag_group <= dg_frag%n_frag) then
      dg_frag%rank_core_lo(:) = dg_frag%frag_core_lo(:, dg_frag%ifrag_group)
      dg_frag%rank_core_hi(:) = dg_frag%frag_core_hi(:, dg_frag%ifrag_group)
      dg_frag%rank_buf_lo(:) = dg_frag%frag_buf_lo(:, dg_frag%ifrag_group)
      dg_frag%rank_buf_hi(:) = dg_frag%frag_buf_hi(:, dg_frag%ifrag_group)
    else
      dg_frag%rank_core_lo(:) = dg_frag%mg%is(:)
      dg_frag%rank_core_hi(:) = dg_frag%mg%ie(:)
      dg_frag%rank_buf_lo(:) = dg_frag%mg%is(:) - dg_frag%nxyz_buffer(:)
      dg_frag%rank_buf_hi(:) = dg_frag%mg%ie(:) + dg_frag%nxyz_buffer(:)
    end if
  end subroutine build_fragment_global_boxes

  logical function is_density_core_point(dg_frag, ifrag, ixg, iyg, izg) result(is_core)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, ixg, iyg, izg
    integer :: xloc, yloc, zloc
    integer :: nx, ny, nz

    nx = dg_frag%nxyz_domain(1, ifrag)
    ny = dg_frag%nxyz_domain(2, ifrag)
    nz = dg_frag%nxyz_domain(3, ifrag)

    xloc = modulo(ixg - dg_frag%frag_core_lo(1, ifrag), dg_frag%lgnum_total(1)) + 1
    yloc = modulo(iyg - dg_frag%frag_core_lo(2, ifrag), dg_frag%lgnum_total(2)) + 1
    zloc = modulo(izg - dg_frag%frag_core_lo(3, ifrag), dg_frag%lgnum_total(3)) + 1
    if (xloc < 1 .or. xloc > nx) then
      is_core = .false.
      return
    end if
    if (yloc < 1 .or. yloc > ny) then
      is_core = .false.
      return
    end if
    if (zloc < 1 .or. zloc > nz) then
      is_core = .false.
      return
    end if

    is_core = .true.
  end function is_density_core_point

  integer function wrap_fragment_cartesian_index(i, ndiv) result(iwrap)
    implicit none
    integer, intent(in) :: i, ndiv

    iwrap = modulo(i - 1, max(1, ndiv)) + 1
  end function wrap_fragment_cartesian_index

  integer function cartesian_index_to_fragment(dg_frag, idx) result(ifrag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: idx(3)

    ifrag = ((idx(1) - 1) * dg_frag%num_fragment(2) + (idx(2) - 1)) * &
            dg_frag%num_fragment(3) + idx(3)
  end function cartesian_index_to_fragment

  integer function find_density_grid_owner(dg_frag, ixg, iyg, izg, hint_rank) result(owner)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ixg, iyg, izg
    integer, intent(in), optional :: hint_rank
    integer :: jrank, first_match, nfrag_ranks, hint_group, rank_group
    integer :: xlo, xhi, ylo, yhi, zlo, zhi

    owner = -1
    first_match = -1
    nfrag_ranks = max(1, dg_frag%isize_frag)

    hint_group = -1
    if (present(hint_rank) .and. .not. dg_frag%parallel_mode_orbital) hint_group = max(0, hint_rank) / nfrag_ranks
    do jrank = 0, dg_frag%isize - 1
      xlo = dg_frag%mg%is_all(1, jrank)
      xhi = dg_frag%mg%ie_all(1, jrank)
      if (ixg < xlo .or. ixg > xhi) cycle
      ylo = dg_frag%mg%is_all(2, jrank)
      yhi = dg_frag%mg%ie_all(2, jrank)
      if (iyg < ylo .or. iyg > yhi) cycle
      zlo = dg_frag%mg%is_all(3, jrank)
      zhi = dg_frag%mg%ie_all(3, jrank)
      if (izg < zlo .or. izg > zhi) cycle
      if (first_match < 0) first_match = jrank
      if (hint_group >= 0) then
        rank_group = jrank / nfrag_ranks
        if (rank_group == hint_group) then
          owner = jrank
          return
        end if
      end if
    end do
    if (first_match >= 0) then
      owner = first_match
    else
      owner = dg_frag%id
    end if
  end function find_density_grid_owner

  !=======================================================================
  ! Initialize RI/DF data for HSE exchange (called once at startup)
  !=======================================================================
  subroutine init_hse_ri_data(dg_frag, system, info)
    use salmon_global, only: yn_hse, yn_hse_ri, hse_omega, hse_ri_ratio, &
                             yn_hse_cd_ri, hse_cd_ri_threshold, &
                             yn_dg_hse_ace, dg_hse_ace_max_age, dg_hse_ace_coef_thresh
    use structures
    use communication, only: comm_is_root, comm_get_groupinfo
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info

    integer :: ifrag, ifrag_local, n_frag_local
    integer :: natom_frag, icount
    integer :: iatom
    real(8), allocatable :: atom_coords_frag(:,:)
    integer, allocatable :: atom_types_frag(:)
    
    ! Check if HSE with RI is enabled
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
    
    ! Calculate number of local fragments
    n_frag_local = dg_frag%ifrag_end - dg_frag%ifrag_start + 1

    ! ACE cache settings are read from the input namelist (yn_dg_hse_ace,
    ! dg_hse_ace_max_age, dg_hse_ace_coef_thresh via salmon_global).
    
    if (comm_is_root(info%id_rko)) then
      write(*,*)
      write(*,*) "Initializing RI/DF approximation for HSE exchange..."
      write(*,'(1x,a,i0)') "  Local fragments: ", n_frag_local
      write(*,'(1x,a,f6.2)') "  N_aux/N_basis ratio: ", hse_ri_ratio
    end if
    
    ! Allocate RI data for local fragments
    allocate(hse_ri_data_frag(n_frag_local))
    if (allocated(hse_ace_vx_cache)) deallocate(hse_ace_vx_cache)
    if (allocated(hse_ace_coef_snapshot)) deallocate(hse_ace_coef_snapshot)
    if (allocated(hse_ace_last_rebuild)) deallocate(hse_ace_last_rebuild)
    if (allocated(hse_ace_call_count)) deallocate(hse_ace_call_count)
    if (allocated(hse_ace_cache_valid)) deallocate(hse_ace_cache_valid)
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
    
    ! Initialize RI data for each local fragment
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
      
      ! Initialize RI data for this fragment.
      ! Pass fragment-interior bounds (1:nxyz_domain) so that init_hse_ri_fragment
      ! uses the correct assumed-shape phi_frag without the global lg%is:lg%ie mismatch.
      call init_hse_ri_fragment(hse_ri_data_frag(ifrag_local), &
                                dg_frag%phi_frag(:,:,:,:,ifrag_local), &
                                dg_frag%lg, dg_frag%mg, &
                                dg_frag%nstate_frag, system%hvol, &
                                natom_frag, atom_coords_frag, atom_types_frag, &
                                hse_omega, &
                                (yn_hse_cd_ri == 'y'), hse_cd_ri_threshold, &
                                max(dg_frag%mg%is(1:3), dg_frag%ixyz_frag(1:3, ifrag)), &
                                min(dg_frag%mg%ie(1:3), dg_frag%ixyz_frag(1:3, ifrag) + dg_frag%nxyz_domain(1:3, ifrag) - 1))
      
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

  logical function atom_belongs_to_fragment(dg_frag, system, ifrag, iatom) result(is_owned)
    use structures, only: s_dft_system
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
        if (ratom < rmin .or. ratom >= rmax) then
          is_owned = .false.
          return
        end if
      else
        if (ratom < rmin .or. ratom > rmax) then
          is_owned = .false.
          return
        end if
      end if
    end do
  end function atom_belongs_to_fragment

  function fragment_cartesian_index(dg_frag, ifrag) result(idx)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag
    integer :: idx(3)
    integer :: rem

    rem = ifrag - 1
    idx(1) = mod(rem, dg_frag%num_fragment(1)) + 1
    rem = rem / dg_frag%num_fragment(1)
    idx(2) = mod(rem, dg_frag%num_fragment(2)) + 1
    rem = rem / dg_frag%num_fragment(2)
    idx(3) = mod(rem, dg_frag%num_fragment(3)) + 1
  end function fragment_cartesian_index

#include "rt_dg_fragment_halo.f90"

  subroutine read_dg_occupation_seed(dg_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_bcast
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: iunit, nspin_file, iostat_open
    integer, allocatable :: occ_file(:)
    logical :: has_file
    character(256) :: filename

    has_file = .false.
    nspin_file = dg_frag%nspin
    if (comm_is_root(dg_frag%id)) then
      filename = './data_dcdft/total/dg_occupation.bin'
      iunit = get_filehandle()
      open(iunit, file=filename, form='unformatted', access='stream', status='old', iostat=iostat_open)
      if (iostat_open == 0) then
        has_file = .true.
        read(iunit) nspin_file
        if (nspin_file >= 1) then
          allocate(occ_file(nspin_file))
          read(iunit) occ_file(1:nspin_file)
        end if
        close(iunit)
      end if
    end if

    call comm_bcast(has_file, dg_frag%icomm, 0)
    call comm_bcast(nspin_file, dg_frag%icomm, 0)
    if (.not. has_file) return
    if (nspin_file < 1) return

    if (.not. allocated(occ_file)) allocate(occ_file(nspin_file))
    call comm_bcast(occ_file, dg_frag%icomm, 0)

    if (nspin_file /= dg_frag%nspin) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,i0,a)') "[WARN] dg_occupation.bin nspin mismatch: file=", nspin_file, &
          " runtime=", dg_frag%nspin, " (ignoring seed occupancy)"
      end if
      deallocate(occ_file)
      return
    end if

    dg_frag%nocc_spin(1:dg_frag%nspin) = min(dg_frag%nstate_tot, max(0, occ_file(1:dg_frag%nspin)))
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,10(1x,i0))') "[INFO] DG occupancy seed loaded:", dg_frag%nocc_spin(1:dg_frag%nspin)
    end if
    deallocate(occ_file)
  end subroutine read_dg_occupation_seed


  !=======================================================================
  ! Read fragment basis data from DC-LCFO calculation (MPI-parallelized)
  !=======================================================================
  subroutine read_fragment_basis_data(dg_frag, bdir_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, comm_summation, comm_get_max
    use salmon_global, only: dg_nmat_cap_mode, dg_nmat_cap_fixed, dg_nmat_cap_multiple, nelec, nelec_spin, &
                             yn_adaptive_basis
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: bdir_frag
    
    character(32), parameter :: binfile_wf = "wavefunctions.bin"
    character(32), parameter :: binfile_bfb = "basis_functions_buffer.bin"
    character(32), parameter :: binfile_rg = "rgrid_index.bin"
    integer, parameter :: basis_buffer_magic = -22022212
    integer, parameter :: basis_buffer_version = 1
    character(256) :: filename
    integer :: iunit, ifrag, ispin, n_frag_file, nspin_file
    integer :: nstate_frag_file, nstate_tot_file
    integer, allocatable :: n_basis_tmp(:,:), index_basis_tmp(:,:,:)
    integer :: n_mat_tmp(2)   ! nspin is expected to be 1 or 2
    integer :: ifrag_count, i_local, io, global_idx
    integer :: local_coef_max, local_idx
    integer :: nxyz_domain(3), nxyz_alloc(3), lgnum_frag(3), lgnum_total(3)
    integer :: nxyz_buffer_file(3), nxyz_box(3)
    integer :: magic_file, version_file
    integer, allocatable :: n_basis_frag(:)
    integer, allocatable :: jxyz_tot(:,:)
    integer :: ix, iy, iz, n
    integer :: ixg_store, iyg_store, izg_store
    integer :: ix_src, iy_src, iz_src
    integer :: ix_box, iy_box, iz_box
    integer :: nb  ! halo width
    integer :: nbasis_iter
    integer :: n_mat_cap, n_mat_cap_env, ienv
    integer :: nocc_max, nocc_eff, ifrag_best, occ_min, occ_max, cap_min, cap_max
    integer :: env_status, env_len
    character(len=64) :: env_n_mat_cap
    logical :: warned_spin_discard, has_buffer_file, identity_seed_coefficients
    real(8) :: cap_avg, weight_best
    real(8), allocatable :: frag_weight_local(:,:,:), frag_weight_sum(:,:,:)
    integer, allocatable :: occ_count(:,:), cap_frag(:,:)
    
    ! Step 1: Root reads metadata from first fragment and broadcasts
    if (comm_is_root(dg_frag%id)) then
      ifrag = 1
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_wf
      
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
      
      if (n_frag_file /= dg_frag%n_frag .or. nspin_file /= dg_frag%nspin) then
        write(*,*) "Error: Fragment basis data mismatch"
        write(*,*) "  Expected n_frag=", dg_frag%n_frag, ", nspin=", dg_frag%nspin, &
                   ", nstate_frag=", dg_frag%nstate_frag
        write(*,*) "  Found    n_frag=", n_frag_file, ", nspin=", nspin_file, &
                   ", nstate_frag=", nstate_frag_file
        stop "DG-Fragment RT: Fragment basis data mismatch"
      end if
      
      close(iunit)
    end if
    
    ! Broadcast metadata to all ranks
    call comm_bcast(n_frag_file, dg_frag%icomm, 0)
    call comm_bcast(nspin_file, dg_frag%icomm, 0)
    call comm_bcast(nstate_frag_file, dg_frag%icomm, 0)
    call comm_bcast(nstate_tot_file, dg_frag%icomm, 0)

    identity_seed_coefficients = (nstate_tot_file < 0)
    if (identity_seed_coefficients) nstate_tot_file = -nstate_tot_file
    if (.not. identity_seed_coefficients) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "[FATAL] dense DG seed coefficients are no longer supported for DG-RT."
        write(*,'(1x,a)') "Regenerate the DC seed so wavefunctions.bin stores identity coefficients via index_basis."
      end if
      stop "DG-Fragment RT: dense seed coefficient block is unsupported"
    end if

    if (nstate_frag_file /= dg_frag%nstate_frag) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a,i0,a,i0,a)') "[INFO] nstate_frag differs: file=", nstate_frag_file, &
          " runtime=", dg_frag%nstate_frag, " (using fragment-state count from file)"
      end if
      dg_frag%nstate_frag = nstate_frag_file
    end if

    ! Use the full state count stored in fragment files (disable occupied-state subset mode).
    if (nstate_tot_file /= dg_frag%nstate_tot) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a,i0,a,i0,a)') "[INFO] nstate_tot differs: file=", nstate_tot_file, &
          " runtime=", dg_frag%nstate_tot, " (using full-state count from file)"
      end if
      dg_frag%nstate_tot = nstate_tot_file
    end if
    
    ! Allocate arrays
    allocate(dg_frag%n_basis(dg_frag%n_frag, dg_frag%nspin))
    if (.not. allocated(dg_frag%index_basis)) then
      allocate(dg_frag%index_basis(dg_frag%nstate_frag, dg_frag%n_frag, dg_frag%nspin))
    end if
    if (.not. allocated(dg_frag%n_mat)) then
      allocate(dg_frag%n_mat(dg_frag%nspin))
    end if
    allocate(n_basis_tmp(dg_frag%n_frag, dg_frag%nspin))
    allocate(index_basis_tmp(dg_frag%nstate_frag, dg_frag%n_frag, dg_frag%nspin))
    
    ! All ranks read global metadata from fragment 1.
    ! index_basis maps local->global indices across ALL fragments; every rank
    ! needs the full table or ranks > 0 will produce zero/NaN current.
    iunit = get_filehandle()
    write(filename, '(a, i6.6, a, a)') trim(bdir_frag), 1, '/', binfile_wf
    open(iunit, file=filename, form='unformatted', access='stream', status='old')
    read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
    if (nstate_tot_file < 0) nstate_tot_file = -nstate_tot_file
    read(iunit) n_mat_tmp(1:dg_frag%nspin)
    read(iunit) n_basis_tmp(1:dg_frag%n_frag, 1:dg_frag%nspin)
    read(iunit) index_basis_tmp(1:dg_frag%nstate_frag, 1:dg_frag%n_frag, 1:dg_frag%nspin)
    close(iunit)
    
    ! Step 3: Gather metadata (now consistent across all ranks)
    dg_frag%n_basis = n_basis_tmp
    dg_frag%index_basis = index_basis_tmp
    dg_frag%n_mat(1:dg_frag%nspin) = n_mat_tmp(1:dg_frag%nspin)
    dg_frag%n_mat_max = maxval(n_mat_tmp(1:dg_frag%nspin))

    ! Optional basis-size cap for fragment-comparison studies.
    ! Normal production runs keep the full fragment basis unless
    ! SALMON_DG_NMAT_CAP is explicitly set.
    n_mat_cap = 0
    if (trim(dg_nmat_cap_mode) == 'fixed' .and. dg_nmat_cap_fixed >= 1) then
      n_mat_cap = dg_nmat_cap_fixed
      if (dg_frag%id == 0) then
        write(*,'(1x,a,a,a,i0)') "[INFO] DG fragment cap mode='", trim(dg_nmat_cap_mode), "' fixed=", n_mat_cap
      end if
    end if
    env_n_mat_cap = ""
    env_status = 1
    env_len = 0
    call get_environment_variable("SALMON_DG_NMAT_CAP", env_n_mat_cap, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_n_mat_cap(1:env_len), *, iostat=ienv) n_mat_cap_env
      if (ienv == 0 .and. n_mat_cap_env >= 1) then
        n_mat_cap = n_mat_cap_env
        if (dg_frag%id == 0) then
          write(*,'(1x,a,i0)') "[INFO] SALMON_DG_NMAT_CAP override applied: ", n_mat_cap
        end if
      else
        if (dg_frag%id == 0) then
          write(*,'(1x,a,a,a)') "[WARN] Ignoring invalid SALMON_DG_NMAT_CAP='", &
                                trim(env_n_mat_cap(1:env_len)), "' (must be integer >= 1)"
        end if
      end if
    end if
    if (n_mat_cap >= 1) then
      dg_frag%n_mat(1:dg_frag%nspin) = min(dg_frag%n_mat(1:dg_frag%nspin), n_mat_cap)
      dg_frag%n_mat_max = maxval(dg_frag%n_mat(1:dg_frag%nspin))
    end if

    ! Enforce cap consistency in index_basis:
    ! indices beyond capped n_mat are masked out to prevent OOB accesses.
    block
      integer :: ispin_cap, ifrag_cap, io_cap, idx_cap, max_keep
      do ispin_cap = 1, dg_frag%nspin
        max_keep = 0
        do ifrag_cap = 1, dg_frag%n_frag
          nbasis_iter = min(dg_frag%n_basis(ifrag_cap, ispin_cap), size(dg_frag%index_basis, 1))
          do io_cap = 1, nbasis_iter
            idx_cap = dg_frag%index_basis(io_cap, ifrag_cap, ispin_cap)
            if (idx_cap < 1 .or. idx_cap > dg_frag%n_mat(ispin_cap)) then
              dg_frag%index_basis(io_cap, ifrag_cap, ispin_cap) = 0
            else
              max_keep = max(max_keep, idx_cap)
            end if
          end do
        end do
        dg_frag%n_mat(ispin_cap) = max_keep
      end do
    end block
    dg_frag%n_mat_max = max(1, maxval(dg_frag%n_mat(1:dg_frag%nspin)))

    ! Compress fragmented/global basis indices to a dense contiguous range.
    ! The DC-LCFO metadata may retain large holes between fragment-local basis blocks,
    ! which inflates n_mat_max and all O(n_mat_max^2) operator matrices.
    block
      integer :: ispin_cmp, ifrag_cmp, io_cmp, idx_cmp, n_old, n_new
      integer, allocatable :: remap(:)
      do ispin_cmp = 1, dg_frag%nspin
        n_old = max(1, dg_frag%n_mat(ispin_cmp))
        allocate(remap(n_old))
        remap = 0
        n_new = 0
        do ifrag_cmp = 1, dg_frag%n_frag
          nbasis_iter = min(dg_frag%n_basis(ifrag_cmp, ispin_cmp), size(dg_frag%index_basis, 1))
          do io_cmp = 1, nbasis_iter
            idx_cmp = dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp)
            if (idx_cmp <= 0) cycle
            if (idx_cmp > n_old) then
              dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp) = 0
              cycle
            end if
            if (remap(idx_cmp) == 0) then
              n_new = n_new + 1
              remap(idx_cmp) = n_new
            end if
            dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp) = remap(idx_cmp)
          end do
        end do
        if (dg_frag%id == 0 .and. n_new < n_old) then
          write(*,'(1x,a,i0,a,i0,a,i0)') "[INFO] Compressed DG basis indices for ispin=", ispin_cmp, &
            ": old n_mat=", n_old, " new n_mat=", n_new
        end if
        if (n_new <= 0) then
          write(*,'(1x,a,i0,a)') "[FATAL] DG basis compression produced zero active basis for ispin=", ispin_cmp, "."
          stop "DG-Fragment RT: zero active basis after index compression"
        end if
        dg_frag%n_mat(ispin_cmp) = n_new
        deallocate(remap)
      end do
    end block
    dg_frag%n_mat_max = max(1, maxval(dg_frag%n_mat(1:dg_frag%nspin)))

    
    ! Validate index_basis uniqueness and coverage (root only)
    if (comm_is_root(dg_frag%id)) then
      block
        integer :: ispin_chk, ifrag_chk, io_chk, global_idx
        integer :: dup_count, out_count, miss_count
        integer, allocatable :: seen(:)
        do ispin_chk = 1, dg_frag%nspin
          allocate(seen(max(1, dg_frag%n_mat_max)))
          seen = 0
          dup_count = 0
          out_count = 0
          do ifrag_chk = 1, dg_frag%n_frag
            nbasis_iter = min(dg_frag%n_basis(ifrag_chk, ispin_chk), size(dg_frag%index_basis, 1))
            do io_chk = 1, nbasis_iter
              global_idx = dg_frag%index_basis(io_chk, ifrag_chk, ispin_chk)
              if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) then
                out_count = out_count + 1
              else
                if (seen(global_idx) == 1) dup_count = dup_count + 1
                seen(global_idx) = 1
              end if
            end do
          end do
          miss_count = count(seen == 0)
          if (dup_count > 0 .or. out_count > 0 .or. miss_count > 0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[WARN] index_basis check (ispin=", ispin_chk, &
              "): dup=", dup_count, " out_of_range=", out_count, " missing=", miss_count
          end if
          deallocate(seen)
        end do
      end block
    end if

    dg_frag%owned_coef_start = 0
    dg_frag%owned_coef_end = -1

    ! Build the row-owner map before allocating coefficient storage.  In
    ! orbital mode this keeps only the subgroup-owned basis rows on each rank,
    ! avoiding a temporary full-fragment coefficient replica during seed load.
    call rebuild_coef_owner_map(dg_frag, "read-fragment-basis")
    local_coef_max = max(1, maxval(dg_frag%local_coef_count(1:dg_frag%nspin)))
    
    ! Step 4: nstate_tot was aligned to file metadata above (full-state mode).

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (allocated(dg_frag%phi_frag_has_seed_buffer)) deallocate(dg_frag%phi_frag_has_seed_buffer)
    allocate(dg_frag%phi_frag_has_seed_buffer(ifrag_count))
    dg_frag%phi_frag_has_seed_buffer(:) = .false.

    if (n_mat_cap < 1 .and. trim(dg_nmat_cap_mode) == 'occ_multiple' .and. dg_nmat_cap_multiple > 0.0d0) then
      if (dg_frag%nspin == 1) then
        nocc_max = max(1, min((nelec + 1) / 2, dg_frag%nstate_tot))
      else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
        nocc_max = max(1, min(maxval(nelec_spin(1:dg_frag%nspin)), dg_frag%nstate_tot))
      else
        nocc_max = max(1, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
      end if

      allocate(frag_weight_local(dg_frag%n_frag, nocc_max, dg_frag%nspin))
      allocate(frag_weight_sum(dg_frag%n_frag, nocc_max, dg_frag%nspin))
      allocate(occ_count(dg_frag%n_frag, dg_frag%nspin))
      allocate(cap_frag(dg_frag%n_frag, dg_frag%nspin))
      frag_weight_local(:, :, :) = 0.0d0
      frag_weight_sum(:, :, :) = 0.0d0
      occ_count(:, :) = 0
      cap_frag(:, :) = 0

      do i_local = 1, ifrag_count
        ifrag = dg_frag%ifrag_start + i_local - 1
        do ispin = 1, dg_frag%nspin
          if (dg_frag%nspin == 1) then
            nocc_eff = max(1, min((nelec + 1) / 2, dg_frag%nstate_tot))
          else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
            nocc_eff = max(1, min(nelec_spin(ispin), dg_frag%nstate_tot))
          else
            nocc_eff = max(1, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
          end if
          nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
          do io = 1, nbasis_iter
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            if (global_idx >= 1 .and. global_idx <= nocc_eff) then
              frag_weight_local(ifrag, global_idx, ispin) = 1.0d0
            end if
          end do
        end do
      end do
      call comm_summation(frag_weight_local, frag_weight_sum, dg_frag%n_frag * nocc_max * dg_frag%nspin, dg_frag%icomm)

      do ispin = 1, dg_frag%nspin
        if (dg_frag%nspin == 1) then
          nocc_eff = max(1, min((nelec + 1) / 2, dg_frag%nstate_tot))
        else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
          nocc_eff = max(1, min(nelec_spin(ispin), dg_frag%nstate_tot))
        else
          nocc_eff = max(1, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
        end if
        do io = 1, nocc_eff
          ifrag_best = 1
          weight_best = frag_weight_sum(1, io, ispin)
          do ifrag = 2, dg_frag%n_frag
            if (frag_weight_sum(ifrag, io, ispin) > weight_best) then
              ifrag_best = ifrag
              weight_best = frag_weight_sum(ifrag, io, ispin)
            end if
          end do
          occ_count(ifrag_best, ispin) = occ_count(ifrag_best, ispin) + 1
        end do
        do ifrag = 1, dg_frag%n_frag
          nbasis_iter = dg_frag%n_basis(ifrag, ispin)
          cap_frag(ifrag, ispin) = min(dg_frag%n_basis(ifrag, ispin), &
                                       int(floor(dg_nmat_cap_multiple * dble(occ_count(ifrag, ispin)))))
          cap_frag(ifrag, ispin) = max(1, cap_frag(ifrag, ispin))
          dg_frag%n_basis(ifrag, ispin) = cap_frag(ifrag, ispin)
          if (.false. .and. dg_frag%id == 0 .and. (cap_frag(ifrag, ispin) >= 60 .or. nbasis_iter >= 60)) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,f8.3,a,i0)') "        basis cap diag: rank=", dg_frag%id, &
              " ifrag=", ifrag, " ispin=", ispin, " before=", nbasis_iter, " occ_count=", occ_count(ifrag, ispin), &
              " multiple=", dg_nmat_cap_multiple, " after=", cap_frag(ifrag, ispin)
            flush(6)
          end if
          do io = cap_frag(ifrag, ispin) + 1, min(dg_frag%nstate_frag, size(dg_frag%index_basis, 1))
            dg_frag%index_basis(io, ifrag, ispin) = 0
          end do
        end do
      end do

      do ispin = 1, dg_frag%nspin
        do ifrag = 1, dg_frag%n_frag
          dg_frag%n_mat(ispin) = max(dg_frag%n_mat(ispin), 1)
        end do
      end do

      block
        integer :: ispin_cmp, ifrag_cmp, io_cmp, idx_cmp, n_old, n_new
        integer, allocatable :: remap(:)
        do ispin_cmp = 1, dg_frag%nspin
          n_old = max(1, dg_frag%n_mat(ispin_cmp))
          allocate(remap(n_old))
          remap = 0
          n_new = 0
          do ifrag_cmp = 1, dg_frag%n_frag
            nbasis_iter = min(dg_frag%n_basis(ifrag_cmp, ispin_cmp), size(dg_frag%index_basis, 1))
            do io_cmp = 1, nbasis_iter
              idx_cmp = dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp)
              if (idx_cmp <= 0) cycle
              if (idx_cmp > n_old) then
                dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp) = 0
                cycle
              end if
              if (remap(idx_cmp) == 0) then
                n_new = n_new + 1
                remap(idx_cmp) = n_new
              end if
              dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp) = remap(idx_cmp)
            end do
          end do
          dg_frag%n_mat(ispin_cmp) = max(1, n_new)
          deallocate(remap)
        end do
      end block
      dg_frag%n_mat_max = max(1, maxval(dg_frag%n_mat(1:dg_frag%nspin)))

      dg_frag%owned_coef_start = 0
      dg_frag%owned_coef_end = -1

      if (dg_frag%id == 0) then
        occ_min = minval(occ_count(:, 1:dg_frag%nspin))
        occ_max = maxval(occ_count(:, 1:dg_frag%nspin))
        cap_min = minval(cap_frag(:, 1:dg_frag%nspin))
        cap_max = maxval(cap_frag(:, 1:dg_frag%nspin))
        cap_avg = sum(dble(cap_frag(:, 1:dg_frag%nspin))) / dble(dg_frag%n_frag * dg_frag%nspin)
        write(*,'(1x,a,a,a,f8.3)') "[INFO] DG fragment cap mode='", trim(dg_nmat_cap_mode), &
          "' multiple=", dg_nmat_cap_multiple
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,f8.3,a,i0)') "[INFO] DG occ/cap summary: occ_min=", occ_min, &
          " occ_max=", occ_max, " cap_min=", cap_min, " cap_max=", cap_max, " cap_avg=", cap_avg, &
          " n_mat_max=", dg_frag%n_mat_max
      end if

      deallocate(frag_weight_local, frag_weight_sum, occ_count, cap_frag)
    end if
    
    ! Reallocate coefficient arrays with correct n_mat_max dimension
    if (allocated(dg_frag%coef)) deallocate(dg_frag%coef)
    if (allocated(dg_frag%coef_new)) deallocate(dg_frag%coef_new)
    if (allocated(dg_frag%coef_work)) deallocate(dg_frag%coef_work)
    allocate(dg_frag%coef(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
    if (yn_adaptive_basis == 'y') then
      allocate(dg_frag%coef_new(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
    end if
    allocate(dg_frag%coef_work(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%coef = 0.0d0
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new = 0.0d0
    dg_frag%coef_work = 0.0d0
    
    ! Step 5: Reorganize coefficient data from fragment-local to global basis order
    ! Now coef is allocated as (n_mat_max, nstate_tot, nspin) matching momentum_mat dimensions
    ! Map fragment-local indices to global indices using index_basis
    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      do ispin = 1, dg_frag%nspin
        ! For each local basis function io in fragment ifrag,
        ! map it to global basis index: global_idx = index_basis(io, ifrag, ispin)
        nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
        do io = 1, nbasis_iter
          global_idx = dg_frag%index_basis(io, ifrag, ispin)
          local_idx = 0
          if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
          if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
            ! DC-DG seed coefficients are the identity after index_basis
            ! compression: local basis row io represents global state global_idx.
            if (global_idx <= dg_frag%nstate_tot) dg_frag%coef(local_idx, global_idx, ispin) = (1.0d0, 0.0d0)
          end if
        end do
      end do
    end do

    ! Keep coefficients only on the owning fragment ranks.

    ! Step 6: Read real-space basis functions and grid mapping
    ! This enables density reconstruction: rho(r) = sum c*_i c_j phi_i(r) phi_j(r)
    ! Also extract fragment geometry information (ixyz_frag, nxyz_domain) for halo exchange
    nb = dg_frag%nxyz_buffer(1)  ! Assume uniform buffer width (4 for 4th-order stencil)
    nxyz_alloc = 0
    warned_spin_discard = .false.
    
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      
      ! Read grid index mapping
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_rg
      
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) lgnum_frag(1:3), lgnum_total(1:3)
      
      if (.not. allocated(jxyz_tot)) then
        allocate(jxyz_tot(maxval(lgnum_frag), 3))
      end if
      do n = 1, 3
        read(iunit) jxyz_tot(1:lgnum_frag(n), n)
      end do
      close(iunit)
      
      ! Extract ixyz_frag (fragment origin in global coordinates, 1-based)
      ! jxyz_tot(1,:) gives the global index of the first grid point in this fragment
      dg_frag%ixyz_frag(1:3, ifrag) = jxyz_tot(1, 1:3)
      
      ! DG-RT requires the DC-exported buffer-aware basis.  The core-only
      ! basis_functions.bin cannot provide fragment-boundary stencil data.
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bfb
      inquire(file=filename, exist=has_buffer_file)
      
      if (.not. has_buffer_file) then
        write(*,'(1x,a,i0,a,a)') "[FATAL] missing DG buffer basis at ifrag=", ifrag, &
          " file=", trim(filename)
        write(*,'(1x,a)') "Regenerate the DC seed so basis_functions_buffer.bin is exported."
        stop "DG-Fragment RT: missing basis buffer file"
      end if
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) magic_file, version_file
      if (magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "Error: invalid basis buffer header at ifrag=", ifrag, &
          " magic=", magic_file, " version=", version_file
        stop "DG-Fragment RT: invalid basis buffer file"
      end if
      read(iunit) nxyz_domain(1:3), nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
      if (nspin_file < 1) then
        write(*,'(1x,a,i0,a,i0)') "Error: invalid nspin_file in basis buffer header at ifrag=", ifrag, &
                                   " nspin_file=", nspin_file
        stop "DG-Fragment RT: invalid nspin_file"
      end if
      do n = 1, 3
        if (dg_frag%num_fragment(n) > 1 .and. nxyz_buffer_file(n) < nb) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] DG seed buffer too small at ifrag=", ifrag, &
            " axis=", n, " seed_buffer=", nxyz_buffer_file(n), " required=", nb
          stop "DG-Fragment RT: insufficient basis buffer"
        end if
      end do
      nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer_file(1:3)
      if (allocated(n_basis_frag)) deallocate(n_basis_frag)
      allocate(n_basis_frag(nspin_file))
      read(iunit) n_basis_frag(1:nspin_file)
      
      ! Store domain size for this fragment
      dg_frag%nxyz_domain(1:3, ifrag) = nxyz_domain(1:3)
      
      ! Orbital fragment parallelism must not partition the fragment basis in
      ! real space.  Each subgroup rank keeps the same full fragment box, then
      ! matrix construction is split over basis rows/columns.
      if (.not. allocated(dg_frag%phi_frag)) then
        if (dg_frag%parallel_mode_orbital) then
          allocate(dg_frag%phi_frag(dg_frag%ixyz_frag(1, ifrag)-nb:dg_frag%ixyz_frag(1, ifrag)+nxyz_domain(1)-1+nb, &
                                     dg_frag%ixyz_frag(2, ifrag)-nb:dg_frag%ixyz_frag(2, ifrag)+nxyz_domain(2)-1+nb, &
                                     dg_frag%ixyz_frag(3, ifrag)-nb:dg_frag%ixyz_frag(3, ifrag)+nxyz_domain(3)-1+nb, &
                                     dg_frag%nstate_frag, ifrag_count))
        else
          allocate(dg_frag%phi_frag(dg_frag%mg%is(1)-nb:dg_frag%mg%ie(1)+nb, &
                                     dg_frag%mg%is(2)-nb:dg_frag%mg%ie(2)+nb, &
                                     dg_frag%mg%is(3)-nb:dg_frag%mg%ie(3)+nb, &
                                     dg_frag%nstate_frag, ifrag_count))
        end if
        dg_frag%phi_frag = 0.0d0  ! Initialize (including halo) to zero
      end if
      
      ! Read basis functions: f_basis(ix,iy,iz,ispin,istate)
      ! phi_frag has no spin dimension: keep spin-1 basis and discard extra spin channels
      ! while still consuming all records to keep stream alignment.
      block
        real(8), allocatable :: phi_box(:,:,:)
        if (nspin_file < 1 .or. nstate_frag_file < 1) then
          write(*,'(1x,a,i0,a,i0,a,i0)') "Error: invalid basis buffer header at ifrag=", ifrag, &
                                         " nspin_file=", nspin_file, " nstate_frag_file=", nstate_frag_file
          stop "DG-Fragment RT: invalid basis buffer header"
        end if
        allocate(phi_box(nxyz_box(1), nxyz_box(2), nxyz_box(3)))
        
        do ispin = 1, nspin_file
          do n = 1, nstate_frag_file
            read(iunit) phi_box(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3))
            
            if (ispin == 1 .and. n <= dg_frag%nstate_frag) then
              ! The buffer file is stored as an unwrapped box around the core.
              ! Unsplitted axes may still wrap within the fragment itself; split
              ! axes must be covered by the seed buffer checked above.
              do izg_store = lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3)
                iz_box = izg_store - (dg_frag%ixyz_frag(3, ifrag) - nxyz_buffer_file(3)) + 1
                if (iz_box < 1 .or. iz_box > nxyz_box(3)) then
                  iz_src = nxyz_buffer_file(3) + modulo(izg_store - dg_frag%ixyz_frag(3, ifrag), nxyz_domain(3)) + 1
                else
                  iz_src = iz_box
                end if
                do iyg_store = lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2)
                  iy_box = iyg_store - (dg_frag%ixyz_frag(2, ifrag) - nxyz_buffer_file(2)) + 1
                  if (iy_box < 1 .or. iy_box > nxyz_box(2)) then
                    iy_src = nxyz_buffer_file(2) + modulo(iyg_store - dg_frag%ixyz_frag(2, ifrag), nxyz_domain(2)) + 1
                  else
                    iy_src = iy_box
                  end if
                  do ixg_store = lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1)
                    ix_box = ixg_store - (dg_frag%ixyz_frag(1, ifrag) - nxyz_buffer_file(1)) + 1
                    if (ix_box < 1 .or. ix_box > nxyz_box(1)) then
                      ix_src = nxyz_buffer_file(1) + modulo(ixg_store - dg_frag%ixyz_frag(1, ifrag), nxyz_domain(1)) + 1
                    else
                      ix_src = ix_box
                    end if
                    dg_frag%phi_frag(ixg_store, iyg_store, izg_store, n, i_local) = &
                      phi_box(ix_src, iy_src, iz_src)
                  end do
                end do
              end do
              dg_frag%phi_frag_has_seed_buffer(i_local) = .true.
            end if
          end do
        end do

        if (nspin_file > 1 .and. .not. warned_spin_discard .and. comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a)') "[WARN] basis_functions_buffer.bin has nspin_file=", nspin_file, &
                                 "; using spin-1 basis only in phi_frag"
          warned_spin_discard = .true.
        end if
        
        if (allocated(phi_box)) deallocate(phi_box)
      end block

      close(iunit)
    end do
    
    ! Clean up temporary arrays
    if (allocated(jxyz_tot)) deallocate(jxyz_tot)
    if (allocated(n_basis_frag)) deallocate(n_basis_frag)
    
    ! Share fragment geometry metadata across all ranks.  The fragment root
    ! ranks are fixed by the stage-1 subgroup layout, so use the same
    ! deterministic mapping that initializes id_array before owner-map builds.
    block
      integer :: owner_rank
      do ifrag = 1, dg_frag%n_frag
        owner_rank = get_fragment_group_root_rank(ifrag, dg_frag%nproc_frag)
        call comm_bcast(dg_frag%ixyz_frag(1:3, ifrag), dg_frag%icomm, owner_rank)
        call comm_bcast(dg_frag%nxyz_domain(1:3, ifrag), dg_frag%icomm, owner_rank)
      end do
    end block
    
    ! Set flag indicating real-space basis functions are available
    dg_frag%has_real_space_basis = .true.
    
    ! Synchronize all ranks before proceeding
    call comm_sync_all(dg_frag%icomm)
    
    ! Clean up
    deallocate(n_basis_tmp, index_basis_tmp)
    
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "Fragment basis data loaded (coefficients + real-space basis)"
      write(*,'(1x,a,i0,a,i0,a,i0)') "  Domain size: ", nxyz_domain(1), " x ", nxyz_domain(2), " x ", nxyz_domain(3)
      write(*,'(1x,a,i0)') "  Number of basis functions per fragment: ", dg_frag%nstate_frag
      write(*,'(1x,a,i0)') "  Number of fragments loaded: ", ifrag_count
    end if
    
  end subroutine read_fragment_basis_data

  !=======================================================================
  ! Initialize Runge-Kutta coefficients
  !=======================================================================
  subroutine init_rk_coefficients(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    
    select case(dg_frag%time_integrator)
    case(1)  ! SSPRK3 (Strong Stability Preserving RK3)
      dg_frag%rk_stages = 3
      allocate(dg_frag%rk_alpha(3), dg_frag%rk_beta(3), dg_frag%rk_gamma(3))
      dg_frag%rk_alpha = [1.0d0, 0.75d0, 1.0d0/3.0d0]
      dg_frag%rk_beta  = [0.0d0, 0.25d0, 2.0d0/3.0d0]
      dg_frag%rk_gamma = [1.0d0, 0.25d0, 2.0d0/3.0d0]
      
    case(3)  ! Classical RK4
      dg_frag%rk_stages = 4
      allocate(dg_frag%rk_alpha(4), dg_frag%rk_beta(4), dg_frag%rk_gamma(4))
      dg_frag%rk_alpha = [0.5d0, 0.5d0, 1.0d0, 0.0d0]
      dg_frag%rk_beta  = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
      dg_frag%rk_gamma = [1.0d0/6.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/6.0d0]
      
    case default  ! AETRS or others - handled separately
      dg_frag%rk_stages = 0
    end select
    
  end subroutine init_rk_coefficients

#include "rt_dg_fragment_hamiltonian.f90"


  !=======================================================================
  ! Update density and Hamiltonian matrix (self-consistent)
  !=======================================================================
#include "rt_dg_density_hamiltonian_update.f90"
  ! Calculate electron density from fragment basis coefficients
  !=======================================================================
#include "rt_dg_density_reconstruct.f90"

  !=======================================================================
  ! Add exact exchange using HSE hybrid functional (Plan A: density matrix)
  ! Fragment-local computation: Each fragment computes exchange independently
  ! Delegates to xc_hse module for actual computation
  !=======================================================================
  subroutine add_exact_exchange_hse(dg_frag, system, H_mat_spin, ifrag, ispin)
    use structures
    use salmon_global, only: hse_alpha, hse_omega, nelec, nelec_spin, yn_hse_ri, &
                             yn_dg_hse_ace, dg_hse_ace_max_age, dg_hse_ace_coef_thresh
    use xc_hse, only: calc_exact_exchange_hse_fragment
    implicit none
    type(s_dg_fragment_rt), intent(in)    :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    real(8),                intent(inout) :: H_mat_spin(:, :)  ! H_mat(:,:,ispin)
    integer,                intent(in)    :: ifrag, ispin
    
    integer :: ifrag_local, n_base_frag, n_occ_frag
    integer :: is(3), ie(3)
    real(8) :: hvol
    real(8), allocatable :: density_matrix(:,:), v_x_tmp(:,:)
    real(8) :: ace_metric
    logical :: rebuild_ace
    
    if (.not. dg_frag%has_real_space_basis) return
    
    hvol = system%hvol
    ! Use fragment-local interior indices (1:nxyz_domain).
    ! phi_frag is allocated with halo as (1-nb:nxyz+nb,...); passing lg%is/lg%ie
    ! (global grid bounds) caused shape mismatch in calc_exact_exchange_hse_fragment.
    is = 1
    ie = dg_frag%nxyz_domain(1:3, ifrag)
    n_base_frag = dg_frag%nstate_frag
    
    ! Fragment-local index
    ifrag_local = ifrag - dg_frag%ifrag_start + 1
    
    ! Count occupied states per spin channel.
    ! nspin==1: spin-paired; each orbital holds 2 electrons -> (nelec+1)/2 avoids
    !           floor-rounding error for odd nelec (e.g. nelec=1 -> 0 without +1).
    ! nspin==2: use spin-resolved counts when available; otherwise fallback to nelec/2.
    if (dg_frag%nspin == 1) then
      n_occ_frag = max(1, min((nelec + 1) / 2, n_base_frag))
    else
      if (sum(nelec_spin(:)) > 0) then
        n_occ_frag = max(1, min(nelec_spin(ispin), n_base_frag))
      else
        n_occ_frag = max(1, min(int(nelec / 2.0d0 + 1.0d-12), n_base_frag))
      end if
    end if
    
    ! Choose method: RI/DF (Plan C) or direct integration (Plan A)
    if (dg_frag%use_hse_ri .and. yn_hse_ri == 'y') then
      if (dg_hse_ace_initialized .and. yn_dg_hse_ace == 'y') then
        hse_ace_call_count(ifrag_local, ispin) = hse_ace_call_count(ifrag_local, ispin) + 1
        rebuild_ace = .not. hse_ace_cache_valid(ifrag_local, ispin)
        ace_metric = 0.0d0
        if (.not. rebuild_ace) then
          call compute_dg_hse_ace_metric(dg_frag, ifrag, ifrag_local, ispin, n_occ_frag, ace_metric)
          if (ace_metric > dg_hse_ace_coef_thresh) rebuild_ace = .true.
          if ((hse_ace_call_count(ifrag_local, ispin) - hse_ace_last_rebuild(ifrag_local, ispin)) >= dg_hse_ace_max_age) then
            rebuild_ace = .true.
          end if
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
        ! Plan C: RI/DF approximation (fast)
        ! Build density matrix from occupied orbitals
        allocate(density_matrix(n_base_frag, n_base_frag))
        call build_density_matrix(dg_frag, ifrag, ispin, n_occ_frag, density_matrix)
        
        ! Call RI version (O(N²N_aux N_occ) complexity)
        call calc_exact_exchange_hse_ri(H_mat_spin, hse_ri_data_frag(ifrag_local), &
                                        density_matrix, hse_alpha, n_occ_frag)
        
        deallocate(density_matrix)
      end if
    else
      ! Plan A: Direct integration (slow but exact)
      call calc_exact_exchange_hse_fragment(H_mat_spin, dg_frag%phi_frag, ifrag_local, &
                                            dg_frag%hgs, hvol, hse_alpha, &
                                            is, ie, n_base_frag, n_occ_frag, nelec)
    end if
    
  end subroutine add_exact_exchange_hse

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

  !=======================================================================
  ! Build density matrix from current wave function coefficients
  ! D_ij = Σ_(occ) c_i^* c_j for occupied orbitals
  !=======================================================================
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

    ! Build density matrix: D_ij = Σ_(occupied states) c_i^* c_j
    ! c_i are read through index_basis(local_i, ifrag, ispin) to preserve
    ! local-fragment basis mapping into global coefficient storage.
    ! RI exchange currently expects a real density matrix; use Hermitian real part.
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

    ! Enforce exact real-symmetric form expected by RI exchange kernel.
    do j = 1, nb
      do i = j + 1, nb
        dval = 0.5d0 * (density_matrix(i, j) + density_matrix(j, i))
        density_matrix(i, j) = dval
        density_matrix(j, i) = dval
      end do
    end do

  end subroutine build_density_matrix

  !=======================================================================
  ! Legacy implementation details (kept for reference)
  !=======================================================================
  ! Original 6-fold loop structure:
  ! do jo = 1, n_base_frag
  !   do io = 1, n_base_frag
  !     V_x_ij = 0.0d0
  !     do lo = 1, n_occ_frag
  !       do ko = 1, n_occ_frag
  !         coulomb_integral_ijkl = 0.0d0
  !         do iz = is(3), ie(3)
  !           do iy = is(2), ie(2)
  !             do ix = is(1), ie(1)
  !               do jz = is(3), ie(3)
  !                 do jy = is(2), ie(2)
  !                   do jx = is(1), ie(1)
  !                     Coulomb kernel: 1/|r1-r2|
  !                   end do
  !                 end do
  !               end do
  !             end do
  !           end do
  !         end do
  !       end do
  !     end do
  !   end do
  ! end do

  !=======================================================================
  ! Reconstruct Hamiltonian matrix with updated potentials
  !=======================================================================
#include "rt_dg_hmat_reconstruct.f90"

  !=======================================================================
  ! Main time evolution iteration for DG-Fragment method
  !=======================================================================
#include "rt_dg_iteration.f90"

  !=======================================================================
  ! Time evolution using Runge-Kutta method
  !=======================================================================
#include "rt_dg_integrator_rk.f90"

  !=======================================================================
  ! Stage-wise density/Hamiltonian update for RK4 (paper-aligned self-consistency)
  !=======================================================================
#include "rt_dg_integrator_stage_update.f90"

  !=======================================================================
  ! Stabilize coefficient matrix by modified Gram-Schmidt orthonormalization
  !=======================================================================
  !=======================================================================
  ! Stabilize coefficient matrix by S-metric Gram-Schmidt orthonormalization
  !=======================================================================
#include "rt_dg_integrator_unitarity.f90"

  !=======================================================================
  ! Calculate time derivative of coefficients (velocity gauge)
  !=======================================================================
#include "rt_dg_integrator_derivative.f90"

  !=======================================================================
  ! Time evolution using AETRS (Enforced Time-Reversal Symmetry)
  !=======================================================================
#include "rt_dg_integrator_aetrs.f90"

  !=======================================================================
  ! Calculate observables in fragment basis
  !=======================================================================
#include "rt_dg_observables.f90"

  !=======================================================================
  ! Finalize DG-Fragment RT calculation
  !=======================================================================
  subroutine finalize_dg_fragment_rt(dg_frag)
    use communication, only: comm_free_group, COMM_GROUP_NULL
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: i
    
    if (allocated(dg_frag%coef)) deallocate(dg_frag%coef)
    if (allocated(dg_frag%coef_new)) deallocate(dg_frag%coef_new)
    if (allocated(dg_frag%coef_work)) deallocate(dg_frag%coef_work)
    if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    if (allocated(dg_frag%H_mat_c)) deallocate(dg_frag%H_mat_c)
    if (allocated(dg_frag%H_mat_blocks)) then
      do i = 1, size(dg_frag%H_mat_blocks)
        if (allocated(dg_frag%H_mat_blocks(i)%val)) deallocate(dg_frag%H_mat_blocks(i)%val)
      end do
      deallocate(dg_frag%H_mat_blocks)
      dg_frag%n_H_blocks = 0
    end if
    if (allocated(dg_frag%H_mat_kinetic_blocks)) then
      do i = 1, size(dg_frag%H_mat_kinetic_blocks)
        if (allocated(dg_frag%H_mat_kinetic_blocks(i)%val)) deallocate(dg_frag%H_mat_kinetic_blocks(i)%val)
      end do
      deallocate(dg_frag%H_mat_kinetic_blocks)
    end if
    if (allocated(dg_frag%H_nl_blocks)) then
      do i = 1, size(dg_frag%H_nl_blocks)
        if (allocated(dg_frag%H_nl_blocks(i)%val)) deallocate(dg_frag%H_nl_blocks(i)%val)
      end do
      deallocate(dg_frag%H_nl_blocks)
      dg_frag%n_H_nl_blocks = 0
    end if
    if (allocated(dg_frag%H_block_map)) deallocate(dg_frag%H_block_map)
    if (allocated(dg_frag%H_nl_block_map)) deallocate(dg_frag%H_nl_block_map)
    if (allocated(dg_frag%H_local_block_ids)) deallocate(dg_frag%H_local_block_ids)
    if (allocated(dg_frag%H_nl_local_block_ids)) deallocate(dg_frag%H_nl_local_block_ids)
    if (allocated(dg_frag%H_local_rows)) deallocate(dg_frag%H_local_rows)
    if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
    if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
    if (allocated(dg_frag%S_mat_c)) deallocate(dg_frag%S_mat_c)
    if (allocated(dg_frag%S_mat_prop_c)) deallocate(dg_frag%S_mat_prop_c)
    if (allocated(dg_frag%S_mat_blocks)) then
      do i = 1, size(dg_frag%S_mat_blocks)
        if (allocated(dg_frag%S_mat_blocks(i)%val)) deallocate(dg_frag%S_mat_blocks(i)%val)
      end do
      deallocate(dg_frag%S_mat_blocks)
      dg_frag%n_S_blocks = 0
    end if
    if (allocated(dg_frag%S_mat_prop_blocks)) then
      do i = 1, size(dg_frag%S_mat_prop_blocks)
        if (allocated(dg_frag%S_mat_prop_blocks(i)%val)) deallocate(dg_frag%S_mat_prop_blocks(i)%val)
      end do
      deallocate(dg_frag%S_mat_prop_blocks)
    end if
    if (allocated(dg_frag%S_mat_blocks_c)) then
      do i = 1, size(dg_frag%S_mat_blocks_c)
        if (allocated(dg_frag%S_mat_blocks_c(i)%val)) deallocate(dg_frag%S_mat_blocks_c(i)%val)
      end do
      deallocate(dg_frag%S_mat_blocks_c)
    end if
    if (allocated(dg_frag%S_mat_prop_blocks_c)) then
      do i = 1, size(dg_frag%S_mat_prop_blocks_c)
        if (allocated(dg_frag%S_mat_prop_blocks_c(i)%val)) deallocate(dg_frag%S_mat_prop_blocks_c(i)%val)
      end do
      deallocate(dg_frag%S_mat_prop_blocks_c)
    end if
    if (allocated(dg_frag%S_block_map)) deallocate(dg_frag%S_block_map)
    if (allocated(dg_frag%n_basis)) deallocate(dg_frag%n_basis)
    if (allocated(dg_frag%index_basis)) deallocate(dg_frag%index_basis)
    if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
    if (allocated(dg_frag%momentum_mat_c)) deallocate(dg_frag%momentum_mat_c)
    if (allocated(dg_frag%gradient_basis_cache)) deallocate(dg_frag%gradient_basis_cache)
    dg_frag%gradient_basis_cache_valid = .false.
    if (allocated(dg_frag%momentum_blocks)) then
      do i = 1, size(dg_frag%momentum_blocks)
        if (allocated(dg_frag%momentum_blocks(i)%val)) deallocate(dg_frag%momentum_blocks(i)%val)
      end do
      deallocate(dg_frag%momentum_blocks)
      dg_frag%n_momentum_blocks = 0
    end if
    if (allocated(dg_frag%momentum_blocks_c)) then
      do i = 1, size(dg_frag%momentum_blocks_c)
        if (allocated(dg_frag%momentum_blocks_c(i)%val)) deallocate(dg_frag%momentum_blocks_c(i)%val)
      end do
      deallocate(dg_frag%momentum_blocks_c)
    end if
    if (allocated(dg_frag%momentum_block_map)) deallocate(dg_frag%momentum_block_map)
    if (allocated(dg_frag%dipole_mat)) deallocate(dg_frag%dipole_mat)
    if (allocated(dg_frag%dipole_blocks)) then
      do i = 1, size(dg_frag%dipole_blocks)
        if (allocated(dg_frag%dipole_blocks(i)%val)) deallocate(dg_frag%dipole_blocks(i)%val)
      end do
      deallocate(dg_frag%dipole_blocks)
      dg_frag%n_dipole_blocks = 0
    end if
    if (allocated(dg_frag%dipole_block_map)) deallocate(dg_frag%dipole_block_map)
    if (allocated(dg_frag%esp)) deallocate(dg_frag%esp)
    if (allocated(dg_frag%nxyz_domain)) deallocate(dg_frag%nxyz_domain)
    if (allocated(dg_frag%density_owner_map)) deallocate(dg_frag%density_owner_map)
    if (allocated(dg_frag%density_primary_local_map)) deallocate(dg_frag%density_primary_local_map)
    if (allocated(dg_frag%density_ixg_map)) deallocate(dg_frag%density_ixg_map)
    if (allocated(dg_frag%density_iyg_map)) deallocate(dg_frag%density_iyg_map)
    if (allocated(dg_frag%density_izg_map)) deallocate(dg_frag%density_izg_map)
    if (allocated(dg_frag%density_subgroup_send_count)) deallocate(dg_frag%density_subgroup_send_count)
    if (allocated(dg_frag%density_subgroup_send_slot_map)) deallocate(dg_frag%density_subgroup_send_slot_map)
    if (allocated(dg_frag%density_grid_points)) deallocate(dg_frag%density_grid_points)
    if (allocated(dg_frag%density_grid_point_count)) deallocate(dg_frag%density_grid_point_count)
    if (allocated(dg_frag%density_subgroup_self_ixg)) deallocate(dg_frag%density_subgroup_self_ixg)
    if (allocated(dg_frag%density_subgroup_self_iyg)) deallocate(dg_frag%density_subgroup_self_iyg)
    if (allocated(dg_frag%density_subgroup_self_izg)) deallocate(dg_frag%density_subgroup_self_izg)
    if (allocated(dg_frag%density_block_nblocks)) deallocate(dg_frag%density_block_nblocks)
    if (allocated(dg_frag%density_block_first_offset)) deallocate(dg_frag%density_block_first_offset)
    if (allocated(dg_frag%density_block_step)) deallocate(dg_frag%density_block_step)
    if (allocated(dg_frag%current_valid_grid_count)) deallocate(dg_frag%current_valid_grid_count)
    if (allocated(dg_frag%current_valid_ix)) deallocate(dg_frag%current_valid_ix)
    if (allocated(dg_frag%current_valid_iy)) deallocate(dg_frag%current_valid_iy)
    if (allocated(dg_frag%current_valid_iz)) deallocate(dg_frag%current_valid_iz)
    if (allocated(dg_frag%current_valid_ixg)) deallocate(dg_frag%current_valid_ixg)
    if (allocated(dg_frag%current_valid_iyg)) deallocate(dg_frag%current_valid_iyg)
    if (allocated(dg_frag%current_valid_izg)) deallocate(dg_frag%current_valid_izg)
    if (allocated(dg_frag%runtime_neighbor_pair_cache)) deallocate(dg_frag%runtime_neighbor_pair_cache)
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) deallocate(dg_frag%momentum_neighbor_pair_cache)
    if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
    if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
    dg_frag%density_phi_block_size = 0
    dg_frag%density_phi_block_cache_valid = .false.
    if (allocated(dg_frag%density_phase_block_cache)) deallocate(dg_frag%density_phase_block_cache)
    dg_frag%density_phase_block_size = 0
    dg_frag%density_phase_block_npw = 0
    dg_frag%density_phase_block_cache_valid = .false.
    if (allocated(dg_frag%density_matrix_frag)) deallocate(dg_frag%density_matrix_frag)
    if (allocated(dg_frag%density_matrix_frag_valid)) deallocate(dg_frag%density_matrix_frag_valid)
    if (allocated(dg_frag%jxyz_tot)) deallocate(dg_frag%jxyz_tot)
    if (allocated(dg_frag%phi_frag)) deallocate(dg_frag%phi_frag)
    if (allocated(dg_frag%phi_frag_has_seed_buffer)) deallocate(dg_frag%phi_frag_has_seed_buffer)
    if (allocated(dg_frag%rk_alpha)) deallocate(dg_frag%rk_alpha)
    if (allocated(dg_frag%rk_beta)) deallocate(dg_frag%rk_beta)
    if (allocated(dg_frag%rk_gamma)) deallocate(dg_frag%rk_gamma)
    if (allocated(dg_frag%H_mat_old)) deallocate(dg_frag%H_mat_old)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
    if (allocated(dg_frag%eigenvalues)) deallocate(dg_frag%eigenvalues)
    if (allocated(dg_frag%basis_overlap)) deallocate(dg_frag%basis_overlap)
    if (allocated(dg_frag%nl_pp_phi_self)) deallocate(dg_frag%nl_pp_phi_self)
    if (allocated(dg_frag%nl_pp_phi_halo)) deallocate(dg_frag%nl_pp_phi_halo)
    dg_frag%nl_pp_phi_cache_valid = .false.
    if (allocated(dg_frag%halo_ireq_send)) deallocate(dg_frag%halo_ireq_send)
    if (allocated(dg_frag%halo_ireq_recv)) deallocate(dg_frag%halo_ireq_recv)
    
    ! Plane wave basis deallocations
    if (allocated(dg_frag%k_pw)) deallocate(dg_frag%k_pw)
    if (allocated(dg_frag%coef_pw)) deallocate(dg_frag%coef_pw)
    if (allocated(dg_frag%coef_pw_full_cache)) deallocate(dg_frag%coef_pw_full_cache)
    dg_frag%coef_pw_full_cache_nstate = 0
    if (allocated(dg_frag%coef_mix)) deallocate(dg_frag%coef_mix)
    if (allocated(dg_frag%mixed_basis_dim)) deallocate(dg_frag%mixed_basis_dim)
    if (allocated(dg_frag%mixed_transform)) deallocate(dg_frag%mixed_transform)
    dg_frag%mixed_basis_ready = .false.
    if (allocated(dg_frag%S_mat_frag_pw)) deallocate(dg_frag%S_mat_frag_pw)
    if (allocated(dg_frag%H_mat_frag_pw)) deallocate(dg_frag%H_mat_frag_pw)
    if (allocated(dg_frag%P_mat_frag_pw)) deallocate(dg_frag%P_mat_frag_pw)
    if (allocated(dg_frag%H_mat_pw_diag)) deallocate(dg_frag%H_mat_pw_diag)
    if (allocated(dg_frag%H_mat_pw)) deallocate(dg_frag%H_mat_pw)
    if (dg_frag%icomm_frag /= COMM_GROUP_NULL) then
      call comm_free_group(dg_frag%icomm_frag)
      dg_frag%icomm_frag = COMM_GROUP_NULL
    end if
    
    ! RI/DF data deallocations
    if (dg_frag%use_hse_ri .and. allocated(hse_ri_data_frag)) then
      call finalize_hse_ri_data()
    end if
    
  end subroutine finalize_dg_fragment_rt

  integer function get_fragment_group_root_rank(ifrag, nproc_frag) result(owner_rank)
    use salmon_global, only: num_fragment
    implicit none
    integer, intent(in) :: ifrag, nproc_frag
    integer :: ix_frag, iy_frag, iz_frag, rem, group_index

    if (ifrag < 1 .or. nproc_frag <= 0) then
      owner_rank = 0
    else
      ! DC fragment files are numbered with z as the fastest index:
      !   ifrag = ((ix-1) * ny + (iy-1)) * nz + iz.
      ! The orbital_sequential MPI layout numbers real-space rank groups with
      ! x as the fastest index.  Convert explicitly so each fragment subgroup
      ! reads the DC file matching its parent real-space box.
      ix_frag = (ifrag - 1) / max(1, num_fragment(2) * num_fragment(3)) + 1
      rem = modulo(ifrag - 1, max(1, num_fragment(2) * num_fragment(3)))
      iy_frag = rem / max(1, num_fragment(3)) + 1
      iz_frag = modulo(rem, max(1, num_fragment(3))) + 1
      group_index = ((iz_frag - 1) * max(1, num_fragment(2)) + (iy_frag - 1)) * &
                    max(1, num_fragment(1)) + (ix_frag - 1)
      owner_rank = group_index * nproc_frag
    end if
  end function get_fragment_group_root_rank

  integer function fragment_from_rank_address(iaddr, nfrag_axis) result(ifrag)
    implicit none
    integer, intent(in) :: iaddr(3), nfrag_axis(3)
    integer :: ix_frag, iy_frag, iz_frag

    ix_frag = iaddr(1) + 1
    iy_frag = iaddr(2) + 1
    iz_frag = iaddr(3) + 1
    ifrag = ((ix_frag - 1) * max(1, nfrag_axis(2)) + (iy_frag - 1)) * &
            max(1, nfrag_axis(3)) + iz_frag
  end function fragment_from_rank_address

  subroutine validate_coef_owner_map(dg_frag, stage_label)
    use communication, only: comm_get_max
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    character(*), intent(in) :: stage_label

    integer :: ispin, global_idx, probe_row, owner_rank
    integer :: owner_min, owner_max, invalid_local, invalid_global
    integer :: probe_count

    if (.not. allocated(dg_frag%coef_owner)) return
    if (size(dg_frag%coef_owner, 1) <= 0 .or. size(dg_frag%coef_owner, 2) <= 0) return

    do ispin = 1, min(dg_frag%nspin, size(dg_frag%coef_owner, 2))
      invalid_local = 0
      do global_idx = 1, size(dg_frag%coef_owner, 1)
        owner_rank = dg_frag%coef_owner(global_idx, ispin)
        if (owner_rank < -1 .or. owner_rank >= dg_frag%isize) invalid_local = invalid_local + 1
      end do
      invalid_global = invalid_local
      call comm_get_max(invalid_global, dg_frag%icomm)
      if (invalid_global > 0) then
        write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "[FATAL] invalid coef owner values after build: stage=", &
          trim(stage_label), " rank=", dg_frag%id, " ispin=", ispin, " invalid_local=", invalid_local
        stop "DG-Fragment RT: invalid coef owner values after build"
      end if

      probe_count = min(8, size(dg_frag%coef_owner, 1))
      do probe_row = 1, probe_count
        owner_rank = dg_frag%coef_owner(probe_row, ispin)
        owner_min = -owner_rank
        call comm_get_max(owner_min, dg_frag%icomm)
        owner_min = -owner_min
        owner_max = owner_rank
        call comm_get_max(owner_max, dg_frag%icomm)
        if (owner_min /= owner_max) then
          write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] inconsistent coef owner after build: stage=", &
            trim(stage_label), " rank=", dg_frag%id, " ispin=", ispin, " row=", probe_row, &
            " owner_min=", owner_min, " owner_max=", owner_max
          stop "DG-Fragment RT: inconsistent coef owner after build"
        end if
      end do
    end do
  end subroutine validate_coef_owner_map

  subroutine rebuild_local_coef_storage(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ispin, global_idx, local_idx, old_local_idx
    integer :: local_coef_max, old_nstate, old_nspin
    integer, allocatable :: new_count(:)
    integer, allocatable :: new_global_ids(:,:)
    integer, allocatable :: new_global_to_local(:,:)
    complex(8), allocatable :: coef_new_rows(:,:,:)
    complex(8), allocatable :: coef_newbuf_rows(:,:,:)
    complex(8), allocatable :: coef_work_rows(:,:,:)
    logical :: has_coef, has_coef_new, has_coef_work

    if (.not. allocated(dg_frag%coef_owner)) return

    allocate(new_count(dg_frag%nspin))
    new_count(:) = 0
    do ispin = 1, dg_frag%nspin
      do global_idx = 1, dg_frag%n_mat_max
        if (dg_frag%coef_owner(global_idx, ispin) == dg_frag%id) then
          new_count(ispin) = new_count(ispin) + 1
        end if
      end do
    end do

    local_coef_max = max(1, maxval(new_count(1:dg_frag%nspin)))
    allocate(new_global_ids(local_coef_max, dg_frag%nspin))
    allocate(new_global_to_local(dg_frag%n_mat_max, dg_frag%nspin))
    new_global_ids(:, :) = 0
    new_global_to_local(:, :) = 0

    do ispin = 1, dg_frag%nspin
      local_idx = 0
      do global_idx = 1, dg_frag%n_mat_max
        if (dg_frag%coef_owner(global_idx, ispin) /= dg_frag%id) cycle
        local_idx = local_idx + 1
        new_global_ids(local_idx, ispin) = global_idx
        new_global_to_local(global_idx, ispin) = local_idx
      end do
    end do

    has_coef = allocated(dg_frag%coef)
    has_coef_new = allocated(dg_frag%coef_new)
    has_coef_work = allocated(dg_frag%coef_work)
    if (has_coef) then
      old_nstate = size(dg_frag%coef, 2)
      old_nspin = min(dg_frag%nspin, size(dg_frag%coef, 3))
      allocate(coef_new_rows(local_coef_max, old_nstate, size(dg_frag%coef, 3)))
      coef_new_rows(:, :, :) = (0.0d0, 0.0d0)
      do ispin = 1, old_nspin
        do local_idx = 1, new_count(ispin)
          global_idx = new_global_ids(local_idx, ispin)
          old_local_idx = 0
          if (allocated(dg_frag%coef_global_to_local)) then
            old_local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
          end if
          if (old_local_idx < 1 .or. old_local_idx > size(dg_frag%coef, 1)) cycle
          coef_new_rows(local_idx, 1:old_nstate, ispin) = &
            dg_frag%coef(old_local_idx, 1:old_nstate, ispin)
        end do
      end do
    end if

    if (has_coef_new) then
      old_nstate = size(dg_frag%coef_new, 2)
      old_nspin = min(dg_frag%nspin, size(dg_frag%coef_new, 3))
      allocate(coef_newbuf_rows(local_coef_max, old_nstate, size(dg_frag%coef_new, 3)))
      coef_newbuf_rows(:, :, :) = (0.0d0, 0.0d0)
      do ispin = 1, old_nspin
        do local_idx = 1, new_count(ispin)
          global_idx = new_global_ids(local_idx, ispin)
          old_local_idx = 0
          if (allocated(dg_frag%coef_global_to_local)) then
            old_local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
          end if
          if (old_local_idx < 1 .or. old_local_idx > size(dg_frag%coef_new, 1)) cycle
          coef_newbuf_rows(local_idx, 1:old_nstate, ispin) = &
            dg_frag%coef_new(old_local_idx, 1:old_nstate, ispin)
        end do
      end do
    end if

    if (has_coef_work) then
      old_nstate = size(dg_frag%coef_work, 2)
      old_nspin = min(dg_frag%nspin, size(dg_frag%coef_work, 3))
      allocate(coef_work_rows(local_coef_max, old_nstate, size(dg_frag%coef_work, 3)))
      coef_work_rows(:, :, :) = (0.0d0, 0.0d0)
      do ispin = 1, old_nspin
        do local_idx = 1, new_count(ispin)
          global_idx = new_global_ids(local_idx, ispin)
          old_local_idx = 0
          if (allocated(dg_frag%coef_global_to_local)) then
            old_local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
          end if
          if (old_local_idx < 1 .or. old_local_idx > size(dg_frag%coef_work, 1)) cycle
          coef_work_rows(local_idx, 1:old_nstate, ispin) = &
            dg_frag%coef_work(old_local_idx, 1:old_nstate, ispin)
        end do
      end do
    end if

    if (allocated(dg_frag%local_coef_count)) deallocate(dg_frag%local_coef_count)
    if (allocated(dg_frag%local_coef_global_ids)) deallocate(dg_frag%local_coef_global_ids)
    if (allocated(dg_frag%coef_global_to_local)) deallocate(dg_frag%coef_global_to_local)
    call move_alloc(new_count, dg_frag%local_coef_count)
    call move_alloc(new_global_ids, dg_frag%local_coef_global_ids)
    call move_alloc(new_global_to_local, dg_frag%coef_global_to_local)

    if (has_coef) then
      deallocate(dg_frag%coef)
      call move_alloc(coef_new_rows, dg_frag%coef)
    end if
    if (has_coef_new) then
      deallocate(dg_frag%coef_new)
      call move_alloc(coef_newbuf_rows, dg_frag%coef_new)
    end if
    if (has_coef_work) then
      deallocate(dg_frag%coef_work)
      call move_alloc(coef_work_rows, dg_frag%coef_work)
    end if
  end subroutine rebuild_local_coef_storage

  subroutine rebuild_coef_owner_map(dg_frag, stage_label)
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: stage_label

    integer :: ispin, ifrag, io, global_idx, nbasis_iter, owner_new
    integer :: env_status, env_len
    integer :: conflict_local
    character(len=64) :: env_owner_io_split
    logical :: owner_root_only
    logical, save :: printed_owner_mode = .false.

    owner_root_only = .true.
    env_owner_io_split = ""
    env_status = 1
    env_len = 0
    call get_environment_variable("SALMON_DG_COEF_OWNER_IO_SPLIT", env_owner_io_split, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      select case(trim(adjustl(env_owner_io_split(1:env_len))))
      case('y','Y','yes','YES','true','TRUE','1')
        owner_root_only = .false.
      end select
    end if

    if (dg_frag%parallel_mode_orbital) then
      ! Fragment-internal orbital parallelism distributes coefficient rows
      ! across the subgroup.  Matrix construction is column-split, while the
      ! propagated coefficient vector is row-owned by the same subgroup ranks.
      owner_root_only = .false.
    end if

    if (dg_frag%id == 0 .and. .not. printed_owner_mode) then
      if (dg_frag%parallel_mode_orbital) then
        write(*,'(1x,a)') "[INFO] DG coef owner mode: orbital row-split"
      else if (owner_root_only) then
        write(*,'(1x,a)') "[INFO] DG coef owner mode: root-only (io split disabled)"
      else
        write(*,'(1x,a)') "[INFO] DG coef owner mode: io-split (legacy)"
      end if
      flush(6)
      printed_owner_mode = .true.
    end if

    if (.not. allocated(dg_frag%id_array)) then
      write(*,'(1x,a,a,a,i0)') "[FATAL] id_array is not allocated before coef owner build: stage=", &
        trim(stage_label), " rank=", dg_frag%id
      stop "DG-Fragment RT: missing id_array before coef owner build"
    end if

    if (allocated(dg_frag%coef_owner)) deallocate(dg_frag%coef_owner)
    allocate(dg_frag%coef_owner(dg_frag%n_mat_max, dg_frag%nspin))
    dg_frag%coef_owner(:, :) = -1
    dg_frag%H_local_block_ids_valid = .false.
    conflict_local = 0
    do ispin = 1, dg_frag%nspin
      do ifrag = 1, dg_frag%n_frag
        nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
        do io = 1, nbasis_iter
          global_idx = dg_frag%index_basis(io, ifrag, ispin)
          if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) cycle
          if (owner_root_only) then
            owner_new = dg_frag%id_array(ifrag)
          else
            owner_new = get_subgroup_block_owner_rank(dg_frag%id_array(ifrag), dg_frag%isize_frag, io, nbasis_iter)
          end if
          if (dg_frag%coef_owner(global_idx, ispin) == -1) then
            dg_frag%coef_owner(global_idx, ispin) = owner_new
          else if (dg_frag%coef_owner(global_idx, ispin) /= owner_new) then
            conflict_local = conflict_local + 1
          end if
        end do
      end do
    end do

    if (dg_frag%id == 0 .and. conflict_local > 0) then
      write(*,'(1x,a,a,a,i0)') "[WARN] coef_owner conflicts detected (kept first owner): stage=", &
        trim(stage_label), " conflicts=", conflict_local
    end if

    dg_frag%owned_coef_start = 0
    dg_frag%owned_coef_end = -1
    do global_idx = 1, dg_frag%n_mat_max
      if (any(dg_frag%coef_owner(global_idx, 1:dg_frag%nspin) == dg_frag%id)) then
        if (dg_frag%owned_coef_start == 0) dg_frag%owned_coef_start = global_idx
        dg_frag%owned_coef_end = global_idx
      end if
    end do

    call validate_coef_owner_map(dg_frag, stage_label)
    call rebuild_local_coef_storage(dg_frag)
  end subroutine rebuild_coef_owner_map

  integer function get_subgroup_block_owner_rank(root_rank, subgroup_size, local_row, n_local_rows) result(owner_rank)
    implicit none
    integer, intent(in) :: root_rank, subgroup_size, local_row, n_local_rows
    integer :: subgroup_root_rank, owner_offset
    integer :: block_base, block_rem, cutoff

    if (subgroup_size <= 1 .or. n_local_rows <= 0) then
      owner_rank = max(0, root_rank)
      return
    end if

    subgroup_root_rank = root_rank - mod(max(0, root_rank), subgroup_size)
    block_base = n_local_rows / subgroup_size
    block_rem = mod(n_local_rows, subgroup_size)
    cutoff = (block_base + 1) * block_rem

    if (local_row <= 0) then
      owner_offset = 0
    else if (block_base <= 0) then
      owner_offset = min(local_row - 1, subgroup_size - 1)
    else if (local_row <= cutoff) then
      owner_offset = (local_row - 1) / (block_base + 1)
    else
      owner_offset = block_rem + (local_row - cutoff - 1) / block_base
    end if

    owner_rank = subgroup_root_rank + min(owner_offset, subgroup_size - 1)
  end function get_subgroup_block_owner_rank

  !=======================================================================
  ! Finalize RI/DF data
  !=======================================================================
  subroutine finalize_hse_ri_data()
    implicit none
    integer :: ifrag_local
    
    if (allocated(hse_ri_data_frag)) then
      do ifrag_local = 1, size(hse_ri_data_frag)
        call deallocate_hse_ri_fragment(hse_ri_data_frag(ifrag_local))
      end do
      deallocate(hse_ri_data_frag)
    end if
    if (allocated(hse_ace_vx_cache)) deallocate(hse_ace_vx_cache)
    if (allocated(hse_ace_coef_snapshot)) deallocate(hse_ace_coef_snapshot)
    if (allocated(hse_ace_last_rebuild)) deallocate(hse_ace_last_rebuild)
    if (allocated(hse_ace_call_count)) deallocate(hse_ace_call_count)
    if (allocated(hse_ace_cache_valid)) deallocate(hse_ace_cache_valid)
    dg_hse_ace_initialized = .false.
    ! dg_hse_ace_enabled is now yn_dg_hse_ace in salmon_global; no reset needed here.
    
  end subroutine finalize_hse_ri_data

#include "rt_dg_fragment_basis_update.f90"

  subroutine diagnose_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)

    call calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
  end subroutine diagnose_density_from_fragments


end module rt_dg_fragment
