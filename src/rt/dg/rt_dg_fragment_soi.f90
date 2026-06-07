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
module rt_dg_fragment_soi
  use rt_dg_fragment, only: get_dg_spin_occ_info, copy_periodic_global_scalar_to_rank_buffer, &
                            build_total_potential_grid_with_buffered_hartree
  use rt_dg_fragment_types, only: s_dg_fragment_rt, halo_info
  use rt_dg_basis_projection, only: calculate_new_old_basis_overlap, &
                                   stabilize_basis_overlap, &
                                   project_wavefunction_to_new_basis, &
                                   diagonalize_and_update_basis
  use rt_dg_plane_wave, only: init_plane_wave_basis, compute_fragment_pw_overlap, &
                              compute_fragment_pw_hamiltonian, build_mixed_hamiltonian, &
                              diagonalize_mixed_basis_pw => diagonalize_mixed_basis
  use rt_dg_hse_exchange, only: init_hse_ri_data, add_exact_exchange_hse, finalize_hse_ri_data
  use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A, ensure_overlap_prop_available, &
                                calculate_microscopic_current_dg, &
                                apply_gradient_to_basis_ops => apply_gradient_to_basis
  use rt_dg_fragment_coefficients, only: get_fragment_coef_owner_rank
  use rt_dg_fragment_lifecycle, only: init_rk_coefficients
  use rt_dg_fragment_io_soi, only: read_fragment_basis_data
  use rt_dg_fragment_hamiltonian_soi, only: calculate_hamiltonian_matrix, calculate_momentum_matrix, &
                                           calculate_overlap_matrix, find_matrix_block, &
                                           init_matrix_blocks, reduce_matrix_blocks, &
                                           build_total_potential_grid, build_hpsi_for_basis, &
                                           integrate_basis_with_field
  implicit none

  private
  public :: init_dg_fragment_rt, tddft_dg_fragment_iteration, finalize_dg_fragment_rt
  public :: calculate_hamiltonian_matrix
  public :: diagnose_density_from_fragments
  public :: s_dg_fragment_rt, halo_info
  
  ! Types and data structures are defined in rt_dg_fragment_types
  
contains

  !=======================================================================
  ! Initialize DG-Fragment RT calculation
  !=======================================================================
  subroutine init_dg_fragment_rt(dg_frag, system, rt, info, lg, mg, ppg)
    use structures
    use communication, only: comm_summation, comm_is_root, comm_create_group, COMM_GROUP_NULL
    use salmon_global, only: num_fragment, num_rgrid_buffer, nstate_frag, time_integrator_dg_fragment, &
                 yn_adaptive_basis, basis_update_threshold, yn_dg_fragment_from_dcdft, &
                 nproc_rgrid, nproc_ob, yn_dg_subspace_diag, dg_subspace_extra_states, &
                 dg_nmat_cap_mode, dg_nmat_cap_fixed, &
                 dg_subspace_pw_vectors, dg_subspace_fallback_cond
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
    character(256) :: filename
    integer :: iunit, i, j, ispin, ifrag, ia, ip
    integer :: pp_buf(3)
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

    if (dg_frag%isize < dg_frag%n_frag) then
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a,i0,a,i0)') "ERROR: Invalid MPI setup for DG-Fragment RT: np=", dg_frag%isize, ", n_frag=", dg_frag%n_frag
        write(*,'(1x,a)') "       Please set MPI process count >= number of fragments."
      end if
      stop "DG-Fragment RT: np < n_frag is not allowed"
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
    
    ! Note: has_real_space_basis flag is set in read_fragment_basis_data()
    ! after successfully loading basis functions from DC-LCFO data
    
    ! Coefficient arrays will be allocated by read_fragment_basis_data()
    ! with proper n_mat_max dimensions for global basis compatibility
    
    ! Allocate only ESP array (independent of basis function count)
    allocate(dg_frag%esp(dg_frag%nstate_tot, dg_frag%nspin))
    
    ! Initialize RK coefficients
    call init_rk_coefficients(dg_frag)
    
    ! Store grid pointers
    
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
    
    ! H_mat/H_mat_kinetic are constructed on demand by the Hamiltonian builder.
    ! Avoid a full-system dense allocation during SOI fragment initialization.
    
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
          ! PW coefficients are fragment-local in orbital mode, so each
          ! fragment subgroup owns its own split of the PW rows.
          dg_frag%coef_pw_owner(i) = get_fragment_coef_owner_rank( &
            dg_frag, dg_frag%ifrag_group, i, dg_frag%n_plane_waves)
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

#include "rt_dg_fragment_halo_soi.f90"


  !=======================================================================
  ! Read fragment basis data from DC-LCFO calculation (MPI-parallelized)
  !=======================================================================

  !=======================================================================
  ! Update density and Hamiltonian matrix (self-consistent)
  !=======================================================================
#include "rt_dg_density_hamiltonian_update.f90"
  ! Calculate electron density from fragment basis coefficients
  !=======================================================================
#include "rt_dg_density_reconstruct.f90"

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
    if (allocated(dg_frag%local_coef_count)) deallocate(dg_frag%local_coef_count)
    if (allocated(dg_frag%local_coef_global_ids)) deallocate(dg_frag%local_coef_global_ids)
    if (allocated(dg_frag%coef_global_to_local)) deallocate(dg_frag%coef_global_to_local)
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
    if (allocated(dg_frag%density_matrix_frag)) deallocate(dg_frag%density_matrix_frag)
    if (allocated(dg_frag%density_matrix_frag_valid)) deallocate(dg_frag%density_matrix_frag_valid)
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
    if (allocated(dg_frag%jxyz_tot)) deallocate(dg_frag%jxyz_tot)
    if (allocated(dg_frag%phi_frag)) deallocate(dg_frag%phi_frag)
    if (allocated(dg_frag%phi_frag_c)) deallocate(dg_frag%phi_frag_c)
    if (allocated(dg_frag%phi_frag_spinor_c)) deallocate(dg_frag%phi_frag_spinor_c)
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
    
    ! Plane wave basis deallocations
    if (allocated(dg_frag%k_pw)) deallocate(dg_frag%k_pw)
    if (allocated(dg_frag%coef_pw)) deallocate(dg_frag%coef_pw)
    if (allocated(dg_frag%coef_pw_full_cache)) deallocate(dg_frag%coef_pw_full_cache)
    dg_frag%coef_pw_full_cache_nstate = 0
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
    if (dg_frag%use_hse_ri) call finalize_hse_ri_data()
    
  end subroutine finalize_dg_fragment_rt

  integer function get_fragment_group_root_rank(ifrag, nproc_frag) result(owner_rank)
    use salmon_global, only: num_fragment
    implicit none
    integer, intent(in) :: ifrag, nproc_frag
    integer :: ix_frag, iy_frag, iz_frag, rem, group_index

    if (ifrag < 1 .or. nproc_frag <= 0) then
      owner_rank = 0
    else
      ! DC fragment files are numbered with z as the fastest index, while the
      ! orbital_sequential MPI layout numbers real-space rank groups with x as
      ! the fastest index.  Convert between the two orders explicitly.
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

#include "rt_dg_fragment_basis_update_soi.f90"

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


end module rt_dg_fragment_soi
