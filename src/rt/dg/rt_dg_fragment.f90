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

  !=======================================================================
  ! Initialize DG-Fragment RT calculation
  !=======================================================================
  subroutine init_dg_fragment_rt(dg_frag, system, rt, info, lg, mg)
    use structures
    use communication, only: comm_summation, comm_is_root, comm_create_group, COMM_GROUP_NULL
    use salmon_global, only: num_fragment, nstate_frag, time_integrator_dg_fragment, &
                 yn_adaptive_basis, basis_update_threshold, yn_dg_fragment_from_dcdft, &
                 nproc_rgrid, yn_dg_subspace_diag, dg_subspace_extra_states, &
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
    
    character(32), parameter :: bdir_frag='./data_dcdft/fragments/'
    character(256) :: filename
    integer :: iunit, i, j, ispin, ifrag
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

    dg_frag%nproc_frag = product(nproc_rgrid)
    if (dg_frag%nproc_frag < 1) then
      stop "DG-Fragment RT requires product(nproc_rgrid) >= 1"
    end if

    if (dg_frag%isize < dg_frag%n_frag) then
      stop "DG-Fragment RT requires np >= n_frag"
    end if
    if (dg_frag%isize /= dg_frag%n_frag * dg_frag%nproc_frag) then
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a,i0,a,i0)') "ERROR: Invalid MPI setup for DG-Fragment RT: np=", dg_frag%isize, ", n_frag=", dg_frag%n_frag
        write(*,'(1x,a,i0)') "       product(nproc_rgrid) = ", dg_frag%nproc_frag
        write(*,'(1x,a)') "       Stage-1 fragment-local MPI requires one subgroup per fragment."
        write(*,'(1x,a)') "       MPI process count must satisfy np = n_frag * product(nproc_rgrid)."
        write(*,'(1x,a)') "       This is a current implementation restriction, not a general DG-RT requirement."
      end if
      stop "DG-Fragment RT stage-1 requires np = n_frag * product(nproc_rgrid)"
    end if

    dg_frag%ifrag_group = dg_frag%id / dg_frag%nproc_frag + 1
    dg_frag%id_frag = mod(dg_frag%id, dg_frag%nproc_frag)
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
    dg_frag%nxyz_buffer(1:3) = 4  ! 4th-order stencil requires ±4 points
    if (comm_is_root(info%id_rko)) then
      write(*,'(1x,a,3i5)') "  nxyz_buffer (runtime): ", dg_frag%nxyz_buffer(1:3)
    end if
    
    ! Allocate fragment geometry arrays
    allocate(dg_frag%nxyz_domain(3, dg_frag%n_frag))
    allocate(dg_frag%ixyz_frag(3, dg_frag%n_frag))
    allocate(dg_frag%id_array(dg_frag%n_frag))
    
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
    
    ! Initialize halo communication structures
    if (dg_frag%has_real_space_basis) then
      call init_halo_communication(dg_frag, info)
    end if
    
    ! Perform initial halo exchange to fill ghost cells
    ! This ensures phi_frag halo regions contain correct neighbor data
    ! Skip if phi_frag was not loaded
    if (dg_frag%has_real_space_basis) then
      call exchange_phi_frag_halo(dg_frag)
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a)') "Initial halo exchange completed for phi_frag"
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
    dg_frag%lg => lg
    dg_frag%mg => mg
    dg_frag%hgs = system%hgs
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
        dg_frag%coef_pw_owner(i) = min(dg_frag%isize - 1, ((i - 1) * dg_frag%isize) / dg_frag%n_plane_waves)
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

  subroutine build_density_grid_owner_maps(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ifrag, i_local, ix, iy, iz, ixg, iyg, izg, owner_rank, source_rank, i
    integer :: nxyz(3), ixyz0(3), ifrag_count
    integer :: nx_max, ny_max, nz_max
    integer, allocatable :: recv_count(:), recv_cursor(:)

    if (allocated(dg_frag%density_owner_map)) deallocate(dg_frag%density_owner_map)
    if (allocated(dg_frag%density_ixg_map)) deallocate(dg_frag%density_ixg_map)
    if (allocated(dg_frag%density_iyg_map)) deallocate(dg_frag%density_iyg_map)
    if (allocated(dg_frag%density_izg_map)) deallocate(dg_frag%density_izg_map)
    if (allocated(dg_frag%density_send_count)) deallocate(dg_frag%density_send_count)
    if (allocated(dg_frag%density_send_slot_map)) deallocate(dg_frag%density_send_slot_map)
    if (allocated(dg_frag%density_recv_map)) then
      do i = lbound(dg_frag%density_recv_map, 1), ubound(dg_frag%density_recv_map, 1)
        if (allocated(dg_frag%density_recv_map(i)%ixg)) deallocate(dg_frag%density_recv_map(i)%ixg)
        if (allocated(dg_frag%density_recv_map(i)%iyg)) deallocate(dg_frag%density_recv_map(i)%iyg)
        if (allocated(dg_frag%density_recv_map(i)%izg)) deallocate(dg_frag%density_recv_map(i)%izg)
      end do
      deallocate(dg_frag%density_recv_map)
    end if
    if (allocated(dg_frag%density_send_count)) deallocate(dg_frag%density_send_count)
    if (allocated(dg_frag%density_send_slot_map)) deallocate(dg_frag%density_send_slot_map)
    if (allocated(dg_frag%density_recv_map)) then
      do ifrag = lbound(dg_frag%density_recv_map, 1), ubound(dg_frag%density_recv_map, 1)
        if (allocated(dg_frag%density_recv_map(ifrag)%ixg)) deallocate(dg_frag%density_recv_map(ifrag)%ixg)
        if (allocated(dg_frag%density_recv_map(ifrag)%iyg)) deallocate(dg_frag%density_recv_map(ifrag)%iyg)
        if (allocated(dg_frag%density_recv_map(ifrag)%izg)) deallocate(dg_frag%density_recv_map(ifrag)%izg)
      end do
      deallocate(dg_frag%density_recv_map)
    end if
    if (.not. associated(dg_frag%mg)) return
    if (.not. allocated(dg_frag%nxyz_domain)) return
    if (.not. allocated(dg_frag%ixyz_frag)) return
    if (dg_frag%ifrag_end < dg_frag%ifrag_start) return

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    nx_max = max(1, maxval(dg_frag%nxyz_domain(1, dg_frag%ifrag_start:dg_frag%ifrag_end)))
    ny_max = max(1, maxval(dg_frag%nxyz_domain(2, dg_frag%ifrag_start:dg_frag%ifrag_end)))
    nz_max = max(1, maxval(dg_frag%nxyz_domain(3, dg_frag%ifrag_start:dg_frag%ifrag_end)))

    allocate(dg_frag%density_owner_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_ixg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_iyg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_izg_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_send_count(0:dg_frag%isize-1))
    allocate(dg_frag%density_send_slot_map(nx_max, ny_max, nz_max, ifrag_count))
    allocate(dg_frag%density_recv_map(0:dg_frag%isize-1))
    dg_frag%density_owner_map = dg_frag%id
    dg_frag%density_ixg_map = 1
    dg_frag%density_iyg_map = 1
    dg_frag%density_izg_map = 1
    dg_frag%density_send_count = 0
    dg_frag%density_send_slot_map = 0

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      nxyz(:) = dg_frag%nxyz_domain(:, ifrag)
      ixyz0(:) = dg_frag%ixyz_frag(:, ifrag)
      do iz = 1, nxyz(3)
        do iy = 1, nxyz(2)
          do ix = 1, nxyz(1)
            ixg = modulo(ixyz0(1) + ix - 2, dg_frag%mg%num(1)) + 1
            iyg = modulo(ixyz0(2) + iy - 2, dg_frag%mg%num(2)) + 1
            izg = modulo(ixyz0(3) + iz - 2, dg_frag%mg%num(3)) + 1
            dg_frag%density_ixg_map(ix, iy, iz, i_local) = ixg
            dg_frag%density_iyg_map(ix, iy, iz, i_local) = iyg
            dg_frag%density_izg_map(ix, iy, iz, i_local) = izg
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg)
            dg_frag%density_owner_map(ix, iy, iz, i_local) = owner_rank
            if (owner_rank /= dg_frag%id) then
              dg_frag%density_send_count(owner_rank) = dg_frag%density_send_count(owner_rank) + 1
              dg_frag%density_send_slot_map(ix, iy, iz, i_local) = dg_frag%density_send_count(owner_rank)
            end if
          end do
        end do
      end do
    end do

    allocate(recv_count(0:dg_frag%isize-1), recv_cursor(0:dg_frag%isize-1))
    recv_count = 0
    do ifrag = 1, dg_frag%n_frag
      source_rank = dg_frag%id_array(ifrag)
      if (source_rank == dg_frag%id) cycle
      nxyz(:) = dg_frag%nxyz_domain(:, ifrag)
      ixyz0(:) = dg_frag%ixyz_frag(:, ifrag)
      do iz = 1, nxyz(3)
        do iy = 1, nxyz(2)
          do ix = 1, nxyz(1)
            ixg = modulo(ixyz0(1) + ix - 2, dg_frag%mg%num(1)) + 1
            iyg = modulo(ixyz0(2) + iy - 2, dg_frag%mg%num(2)) + 1
            izg = modulo(ixyz0(3) + iz - 2, dg_frag%mg%num(3)) + 1
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg)
            if (owner_rank == dg_frag%id) recv_count(source_rank) = recv_count(source_rank) + 1
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
    do ifrag = 1, dg_frag%n_frag
      source_rank = dg_frag%id_array(ifrag)
      if (source_rank == dg_frag%id) cycle
      nxyz(:) = dg_frag%nxyz_domain(:, ifrag)
      ixyz0(:) = dg_frag%ixyz_frag(:, ifrag)
      do iz = 1, nxyz(3)
        do iy = 1, nxyz(2)
          do ix = 1, nxyz(1)
            ixg = modulo(ixyz0(1) + ix - 2, dg_frag%mg%num(1)) + 1
            iyg = modulo(ixyz0(2) + iy - 2, dg_frag%mg%num(2)) + 1
            izg = modulo(ixyz0(3) + iz - 2, dg_frag%mg%num(3)) + 1
            owner_rank = find_density_grid_owner(dg_frag, ixg, iyg, izg)
            if (owner_rank /= dg_frag%id) cycle
            recv_cursor(source_rank) = recv_cursor(source_rank) + 1
            dg_frag%density_recv_map(source_rank)%ixg(recv_cursor(source_rank)) = ixg
            dg_frag%density_recv_map(source_rank)%iyg(recv_cursor(source_rank)) = iyg
            dg_frag%density_recv_map(source_rank)%izg(recv_cursor(source_rank)) = izg
          end do
        end do
      end do
    end do

    deallocate(recv_count, recv_cursor)
  end subroutine build_density_grid_owner_maps

  integer function find_density_grid_owner(dg_frag, ixg, iyg, izg) result(owner)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ixg, iyg, izg
    integer :: jrank

    owner = dg_frag%id
    do jrank = 0, dg_frag%isize - 1
      if (ixg < dg_frag%mg%is_all(1, jrank) .or. ixg > dg_frag%mg%ie_all(1, jrank)) cycle
      if (iyg < dg_frag%mg%is_all(2, jrank) .or. iyg > dg_frag%mg%ie_all(2, jrank)) cycle
      if (izg < dg_frag%mg%is_all(3, jrank) .or. izg > dg_frag%mg%ie_all(3, jrank)) cycle
      owner = jrank
      return
    end do
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
                                [1, 1, 1], dg_frag%nxyz_domain(1:3, ifrag))
      
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


  !=======================================================================
  ! Read fragment basis data from DC-LCFO calculation (MPI-parallelized)
  !=======================================================================
  subroutine read_fragment_basis_data(dg_frag, bdir_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, comm_summation
    use salmon_global, only: dg_nmat_cap_mode, dg_nmat_cap_fixed, dg_nmat_cap_multiple, nelec, nelec_spin
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: bdir_frag
    
    character(32), parameter :: binfile_wf = "wavefunctions.bin"
    character(32), parameter :: binfile_bf = "basis_functions.bin"
    character(32), parameter :: binfile_rg = "rgrid_index.bin"
    character(256) :: filename
    integer :: iunit, ifrag, ispin, n_frag_file, nspin_file
    integer :: nstate_frag_file, nstate_tot_file
    integer, allocatable :: n_basis_tmp(:,:), index_basis_tmp(:,:,:)
    real(8), allocatable :: coef_local(:,:,:,:), coef_tmp(:,:,:)  ! local coef buffers
    integer :: n_mat_tmp(2)   ! nspin is expected to be 1 or 2
    integer :: ifrag_count, i_local, io, global_idx
    integer :: nxyz_domain(3), nxyz_alloc(3), nxyz_new(3), lgnum_frag(3), lgnum_total(3)
    integer, allocatable :: n_basis_frag(:)
    integer, allocatable :: jxyz_tot(:,:)
    integer :: ix, iy, iz, n
    integer :: nb  ! halo width
    integer :: nbasis_iter
    integer :: n_mat_cap, n_mat_cap_env, ienv
    integer :: nocc_max, nocc_eff, ifrag_best, occ_min, occ_max, cap_min, cap_max
    integer :: env_status, env_len
    character(len=64) :: env_n_mat_cap
    logical :: warned_spin_discard
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

    if (allocated(dg_frag%coef_owner)) deallocate(dg_frag%coef_owner)
    allocate(dg_frag%coef_owner(dg_frag%n_mat_max, dg_frag%nspin))
    dg_frag%coef_owner(:, :) = -1
    do ispin = 1, dg_frag%nspin
      do ifrag = 1, dg_frag%n_frag
        nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
        do io = 1, nbasis_iter
          global_idx = dg_frag%index_basis(io, ifrag, ispin)
          if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) cycle
          dg_frag%coef_owner(global_idx, ispin) = dg_frag%id_array(ifrag)
        end do
      end do
    end do
    dg_frag%owned_coef_start = 0
    dg_frag%owned_coef_end = -1
    do global_idx = 1, dg_frag%n_mat_max
      if (any(dg_frag%coef_owner(global_idx, 1:dg_frag%nspin) == dg_frag%id)) then
        if (dg_frag%owned_coef_start == 0) dg_frag%owned_coef_start = global_idx
        dg_frag%owned_coef_end = global_idx
      end if
    end do
    
    ! Step 4: nstate_tot was aligned to file metadata above (full-state mode).

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    allocate(coef_local(dg_frag%nstate_frag, dg_frag%nstate_tot, dg_frag%nspin, ifrag_count))
    coef_local = 0.0d0
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_wf
      
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      
      ! Read header
      read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
      
      ! Read metadata blocks
      read(iunit) n_mat_tmp(1:dg_frag%nspin)
      read(iunit) n_basis_tmp(1:dg_frag%n_frag, 1:dg_frag%nspin)
      read(iunit) index_basis_tmp(1:dg_frag%nstate_frag, 1:dg_frag%n_frag, 1:dg_frag%nspin)

      ! Read coefficient block with file dimensions and map into runtime dimensions
      allocate(coef_tmp(dg_frag%nstate_frag, nstate_tot_file, dg_frag%nspin))
      coef_tmp = 0.0d0
      read(iunit) coef_tmp(1:dg_frag%nstate_frag, 1:nstate_tot_file, 1:dg_frag%nspin)
      coef_local(1:dg_frag%nstate_frag, 1:min(dg_frag%nstate_tot, nstate_tot_file), 1:dg_frag%nspin, i_local) = &
        coef_tmp(1:dg_frag%nstate_frag, 1:min(dg_frag%nstate_tot, nstate_tot_file), 1:dg_frag%nspin)
      deallocate(coef_tmp)
      
      close(iunit)
    end do

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
          do io = 1, nocc_eff
            frag_weight_local(ifrag, io, ispin) = sum(abs(coef_local(1:nbasis_iter, io, ispin, i_local))**2)
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

      if (allocated(dg_frag%coef_owner)) deallocate(dg_frag%coef_owner)
      allocate(dg_frag%coef_owner(dg_frag%n_mat_max, dg_frag%nspin))
      dg_frag%coef_owner(:, :) = -1
      do ispin = 1, dg_frag%nspin
        do ifrag = 1, dg_frag%n_frag
          nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
          do io = 1, nbasis_iter
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) cycle
            dg_frag%coef_owner(global_idx, ispin) = dg_frag%id_array(ifrag)
          end do
        end do
      end do
      dg_frag%owned_coef_start = 0
      dg_frag%owned_coef_end = -1
      do global_idx = 1, dg_frag%n_mat_max
        if (any(dg_frag%coef_owner(global_idx, 1:dg_frag%nspin) == dg_frag%id)) then
          if (dg_frag%owned_coef_start == 0) dg_frag%owned_coef_start = global_idx
          dg_frag%owned_coef_end = global_idx
        end if
      end do

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
    allocate(dg_frag%coef(dg_frag%n_mat_max, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(dg_frag%coef_new(dg_frag%n_mat_max, dg_frag%nstate_tot, dg_frag%nspin))
    allocate(dg_frag%coef_work(dg_frag%n_mat_max, dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%coef = 0.0d0
    dg_frag%coef_new = 0.0d0
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
          if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) then
            ! Copy coefficient for this global basis function
            dg_frag%coef(global_idx, 1:dg_frag%nstate_tot, ispin) = &
              coef_local(io, 1:dg_frag%nstate_tot, ispin, i_local)
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
      
      ! Read basis functions
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bf
      
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) nxyz_domain(1:3), nspin_file, nstate_frag_file
      if (nspin_file < 1) then
        write(*,'(1x,a,i0,a,i0)') "Error: invalid nspin_file in basis_functions header at ifrag=", ifrag, &
                                   " nspin_file=", nspin_file
        stop "DG-Fragment RT: invalid nspin_file"
      end if
      if (allocated(n_basis_frag)) deallocate(n_basis_frag)
      allocate(n_basis_frag(nspin_file))
      read(iunit) n_basis_frag(1:nspin_file)
      
      ! Store domain size for this fragment
      dg_frag%nxyz_domain(1:3, ifrag) = nxyz_domain(1:3)
      
      ! Allocate/expand phi_frag with halo regions.
      ! Use max local fragment domain so heterogeneous fragment sizes are safe.
      if (.not. allocated(dg_frag%phi_frag)) then
        nxyz_alloc = nxyz_domain
        allocate(dg_frag%phi_frag(1-nb:nxyz_alloc(1)+nb, &
                                   1-nb:nxyz_alloc(2)+nb, &
                                   1-nb:nxyz_alloc(3)+nb, &
                                   dg_frag%nstate_frag, ifrag_count))
        dg_frag%phi_frag = 0.0d0  ! Initialize (including halo) to zero
      else
        nxyz_new = max(nxyz_alloc, nxyz_domain)
        if (any(nxyz_new /= nxyz_alloc)) then
          block
            real(8), allocatable :: phi_frag_grow(:,:,:,:,:)
            allocate(phi_frag_grow(1-nb:nxyz_new(1)+nb, &
                                   1-nb:nxyz_new(2)+nb, &
                                   1-nb:nxyz_new(3)+nb, &
                                   dg_frag%nstate_frag, ifrag_count))
            phi_frag_grow = 0.0d0
            phi_frag_grow(1-nb:nxyz_alloc(1)+nb, 1-nb:nxyz_alloc(2)+nb, 1-nb:nxyz_alloc(3)+nb, :, :) = &
              dg_frag%phi_frag(1-nb:nxyz_alloc(1)+nb, 1-nb:nxyz_alloc(2)+nb, 1-nb:nxyz_alloc(3)+nb, :, :)
            call move_alloc(phi_frag_grow, dg_frag%phi_frag)
          end block
          nxyz_alloc = nxyz_new
        end if
      end if
      
      ! Read basis functions: f_basis(ix,iy,iz,ispin,istate)
      ! phi_frag has no spin dimension: keep spin-1 basis and discard extra spin channels
      ! while still consuming all records to keep stream alignment.
      block
        real(8), allocatable :: phi_tmp(:,:,:)
        if (nspin_file < 1 .or. nstate_frag_file < 1) then
          write(*,'(1x,a,i0,a,i0,a,i0)') "Error: invalid basis_functions header at ifrag=", ifrag, &
                                         " nspin_file=", nspin_file, " nstate_frag_file=", nstate_frag_file
          stop "DG-Fragment RT: invalid basis_functions header"
        end if
        allocate(phi_tmp(nxyz_domain(1), nxyz_domain(2), nxyz_domain(3)))
        
        do ispin = 1, nspin_file
          do n = 1, nstate_frag_file
            ! Read basis function for interior domain (1:nx, 1:ny, 1:nz)
            read(iunit) phi_tmp(1:nxyz_domain(1), 1:nxyz_domain(2), 1:nxyz_domain(3))
            
            if (ispin == 1 .and. n <= dg_frag%nstate_frag) then
              ! Store in phi_frag (interior only; halo will be filled by exchange)
              dg_frag%phi_frag(1:nxyz_domain(1), 1:nxyz_domain(2), 1:nxyz_domain(3), n, i_local) = &
                phi_tmp(1:nxyz_domain(1), 1:nxyz_domain(2), 1:nxyz_domain(3))
            end if
          end do
        end do

        if (nspin_file > 1 .and. .not. warned_spin_discard .and. comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a)') "[WARN] basis_functions.bin has nspin_file=", nspin_file, &
                                 "; using spin-1 basis only in phi_frag"
          warned_spin_discard = .true.
        end if
        
        deallocate(phi_tmp)
      end block
      
      close(iunit)
    end do
    
    ! Clean up temporary arrays
    if (allocated(jxyz_tot)) deallocate(jxyz_tot)
    if (allocated(n_basis_frag)) deallocate(n_basis_frag)
    
    ! CRITICAL: Share fragment geometry metadata across all ranks for Halo initialization
    ! id_array is not initialized yet here, so reconstruct owner rank using
    ! the same block distribution rule as distribute_fragments().
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
    deallocate(n_basis_tmp, index_basis_tmp, coef_local)
    
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
    integer :: io, istate, ig, nb
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
      if (ig < 1 .or. ig > size(dg_frag%coef, 1)) cycle
      do istate = 1, min(n_occ, dg_frag%nstate_tot)
        c_now = dg_frag%coef(ig, istate, ispin)
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
    integer :: io, ig, nocc_eff

    nocc_eff = min(n_occ, dg_frag%nstate_tot)
    hse_ace_coef_snapshot(:, :, ispin, ifrag_local) = (0.0d0, 0.0d0)
    do io = 1, min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
      ig = dg_frag%index_basis(io, ifrag, ispin)
      if (ig < 1 .or. ig > size(dg_frag%coef, 1)) cycle
      hse_ace_coef_snapshot(io, 1:nocc_eff, ispin, ifrag_local) = dg_frag%coef(ig, 1:nocc_eff, ispin)
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

    integer :: i, j, iocc, istate, iidx, jidx, nb
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
      if (jidx < 1 .or. jidx > size(dg_frag%coef, 1)) cycle
      do i = 1, nb
        iidx = dg_frag%index_basis(i, ifrag, ispin)
        if (iidx < 1 .or. iidx > size(dg_frag%coef, 1)) cycle
        do iocc = 1, n_occ
          istate = iocc
          if (istate > dg_frag%nstate_tot) cycle

          coef_i = dg_frag%coef(iidx, istate, ispin)
          coef_j = dg_frag%coef(jidx, istate, ispin)

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
    if (allocated(dg_frag%density_phi_cache)) deallocate(dg_frag%density_phi_cache)
    if (allocated(dg_frag%jxyz_tot)) deallocate(dg_frag%jxyz_tot)
    if (allocated(dg_frag%phi_frag)) deallocate(dg_frag%phi_frag)
    if (allocated(dg_frag%rk_alpha)) deallocate(dg_frag%rk_alpha)
    if (allocated(dg_frag%rk_beta)) deallocate(dg_frag%rk_beta)
    if (allocated(dg_frag%rk_gamma)) deallocate(dg_frag%rk_gamma)
    if (allocated(dg_frag%H_mat_old)) deallocate(dg_frag%H_mat_old)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
    if (allocated(dg_frag%eigenvalues)) deallocate(dg_frag%eigenvalues)
    if (allocated(dg_frag%basis_overlap)) deallocate(dg_frag%basis_overlap)
    if (allocated(dg_frag%H_nl_cache)) deallocate(dg_frag%H_nl_cache)
    
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
    if (allocated(dg_frag%H_mat_pw_diag)) deallocate(dg_frag%H_mat_pw_diag)
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
    implicit none
    integer, intent(in) :: ifrag, nproc_frag

    if (ifrag < 1 .or. nproc_frag <= 0) then
      owner_rank = 0
    else
      owner_rank = (ifrag - 1) * nproc_frag
    end if
  end function get_fragment_group_root_rank

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


end module rt_dg_fragment
