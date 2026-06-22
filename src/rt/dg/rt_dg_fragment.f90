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
! RT BASIS POLICY:
! ================
! Runtime adaptive basis updates have been removed from the production RT path.
! The RT basis must be generated in the DC-LCFO export path instead.
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
  use rt_dg_fragment_types, only: s_dg_fragment_rt, halo_info, invalidate_coef_exchange_cache
  use rt_dg_plane_wave, only: init_plane_wave_basis, initialize_wpw_windows
  use rt_dg_hse_exchange, only: init_hse_ri_data, add_exact_exchange_hse, finalize_hse_ri_data
  use rt_dg_fragment_layout, only: build_density_grid_owner_maps, &
                                  wrap_global_grid_index, get_fragment_grid_sender_rank, &
                                  wrap_fragment_cartesian_index, cartesian_index_to_fragment, &
                                  find_density_grid_owner, get_fragment_group_root_rank, &
                                  fragment_from_rank_address
  use rt_dg_fragment_lifecycle, only: init_rk_coefficients
  use rt_dg_fragment_io, only: read_fragment_basis_data, initialize_coefficients_from_buffer_wannier_flux
  use rt_dg_fragment_coefficients, only: rebuild_coef_owner_map, get_subgroup_block_owner_rank
  use rt_dg_fragment_ops, only: ensure_nonlocal_pp_matrix_A, ensure_overlap_prop_available, &
                                calculate_microscopic_current_dg, &
                                apply_gradient_to_basis_ops => apply_gradient_to_basis, &
                                rebuild_local_h_block_ids, &
                                apply_matrix_blocks_batch, apply_complex_matrix_blocks_batch, &
                                apply_overlap_operator_batch, fetch_remote_coef_rows, &
                                solve_overlap_operator_batch
  implicit none

  private
  public :: init_dg_fragment_rt, tddft_dg_fragment_iteration, finalize_dg_fragment_rt
  public :: calculate_hamiltonian_matrix
  public :: calculate_observables
  public :: update_density_and_hamiltonian
  public :: diagnose_density_from_fragments
  public :: diagnose_dcdft_lcfo_seed_stationarity
  public :: calibrate_dcdft_lcfo_static_hamiltonian
  public :: refresh_buffer_wannier_flux_seed_from_current_hamiltonian
  public :: get_dg_spin_occ_info
  public :: copy_periodic_global_scalar_to_rank_buffer
  public :: build_total_potential_grid_with_buffered_hartree
  public :: s_dg_fragment_rt, halo_info
  
  ! Types and data structures are defined in rt_dg_fragment_types
  
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
                 dg_subspace_pw_vectors, nelec, nelec_spin, process_allocation
    use density_matrix_and_energy_plusU_sub, only: PLUS_U_ON
    use filesystem, only: get_filehandle
    use unusedvar_mod, only: salmon_unusedvar
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
    
    call salmon_unusedvar(rt)

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
    dg_frag%coef_state_block_mode = .false.
    dg_frag%coef_state_start = 1
    dg_frag%coef_state_end = dg_frag%nstate_tot
    dg_frag%coef_nstate_local = dg_frag%nstate_tot
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
        write(*,'(1x,a,a)') "ERROR: invalid DG fragment parallel mode=", trim(dg_frag%parallel_mode)
      end if
      stop "DG-Fragment RT: parallel mode must be orbital or legacy_realspace"
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
    call initialize_dg_state_block_coef_mode(dg_frag)
    
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
    case('taylor4pc')
      dg_frag%time_integrator = 4
    case('expdiag')
      dg_frag%time_integrator = 5
    case default
      dg_frag%time_integrator = 5  ! default: expdiag
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
    
    ! read_fragment_basis_data may already have loaded DC-LCFO EigenExa
    ! eigenvalues for static-phase removal.  Keep them when present.
    if (.not. allocated(dg_frag%esp)) then
      allocate(dg_frag%esp(dg_frag%nstate_tot, dg_frag%nspin))
      dg_frag%esp(:, :) = 0.0d0
    end if
    
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
    if (dg_frag%yn_adaptive_basis) then
      if (comm_is_root(info%id_rko)) then
        write(*,'(1x,a)') "[FATAL] yn_adaptive_basis is deprecated for DG-Fragment RT."
        write(*,'(1x,a)') "[FATAL] Runtime adaptive basis updates were removed."
        write(*,'(1x,a)') "[FATAL] Generate the RT basis in the DC-LCFO export path instead."
      end if
      stop "DG-Fragment RT: runtime adaptive basis update removed"
    end if
    dg_frag%use_subspace_diag = (yn_dg_subspace_diag == 'y')
    dg_frag%subspace_extra_states = max(0, dg_subspace_extra_states)
    dg_frag%subspace_pw_vectors = max(0, dg_subspace_pw_vectors)
    dg_frag%last_subspace_dim = 0
    dg_frag%has_nl_cache = .false.
    dg_frag%has_nl_projector_cache = .false.
    dg_frag%Ac_nl_cache = 0.0d0
    dg_frag%Ac_nl_cache_tol = 1.0d-12
    dg_frag%hamiltonian_change_norm = 0.0d0
    dg_frag%nbasis_update_count = 0
    dg_frag%last_basis_update_step = 0
    
    if (comm_is_root(info%id_rko)) then
      write(*,'(1x,a,l1)') "  DG subspace diagonalization: ", dg_frag%use_subspace_diag
      write(*,'(1x,a,i0)') "  DG subspace extra states: ", dg_frag%subspace_extra_states
      write(*,'(1x,a,i0)') "  DG subspace PW vectors: ", dg_frag%subspace_pw_vectors
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
    call initialize_wpw_windows(dg_frag, info)
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
    call diagnose_mixed_wannier_pw_fsum(dg_frag)
    call maybe_build_mixed_wannier_bpw_position(dg_frag)
    
    ! Initialize RI/DF approximation for HSE if enabled
    call init_hse_ri_data(dg_frag, system, info)

    dg_frag%current(:) = 0.0d0
    dg_frag%current_para(:) = 0.0d0
    dg_frag%current_para_source_norm(:) = 0.0d0
    dg_frag%current_para_bound(:) = 0.0d0
    dg_frag%current_coef_norm = 0.0d0
    dg_frag%current_coef_imag_norm = 0.0d0
    dg_frag%current_nl(:) = 0.0d0
    dg_frag%current_dia(:) = 0.0d0
    dg_frag%current_total(:) = 0.0d0
    dg_frag%polarization_lg(:) = 0.0d0
    dg_frag%polarization_lg_ref(:) = 0.0d0
    dg_frag%polarization_lg_ref_ready = .false.
    dg_frag%dipole_lg_raw(:) = 0.0d0
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
      if (trim(time_integrator_dg_fragment) == 'taylor4pc') then
        write(*,'(1x,a)')  "  Time integrator: Taylor4 predictor-corrector explicitly selected."
      end if
      if (dg_frag%use_plane_wave_basis) then
        write(*,'(1x,a,i0)') "  Plane waves added: ", dg_frag%n_plane_waves
        write(*,'(1x,a,f8.4,a)') "  PW cutoff: ", dg_frag%k_cutoff_pw, " a.u.^-1"
      end if
      if (dg_frag%use_hse_ri) then
        write(*,'(1x,a)') "  HSE RI/DF approximation: ENABLED"
      end if
    end if
    
  end subroutine init_dg_fragment_rt

#include "rt_dg_fragment_halo.f90"


  !=======================================================================
  ! Rotate the DC buffer basis into a fragment-local eigenbasis, then keep
  ! only the requested cap.  This is intentionally run after the first H/S
  ! construction and before the production RT operator build.
  !=======================================================================

#include "rt_dg_fragment_hamiltonian.f90"


  !=======================================================================
  ! Update density and Hamiltonian matrix (self-consistent)
  !=======================================================================
#include "rt_dg_density_hamiltonian_update.f90"
  ! Calculate electron density from fragment basis coefficients
  !=======================================================================
#include "rt_dg_density_reconstruct.f90"

  !=======================================================================
  ! Diagnose the f-sum recovered by a Wannier occupied + orthogonalized PW
  ! virtual space before enabling mixed-basis RT propagation.
  !=======================================================================
#include "rt_dg_mixed_fsum_diagnose.f90"

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
  ! Time evolution using fourth-order Taylor predictor-corrector
  !=======================================================================
#include "rt_dg_integrator_taylor.f90"

  ! Time evolution by direct diagonal exponential in local BPW blocks
#include "rt_dg_integrator_expdiag.f90"

  !=======================================================================
  ! Stage-wise density/Hamiltonian update for RK4 (paper-aligned self-consistency)
  !=======================================================================
#include "rt_dg_integrator_stage_update.f90"

  !=======================================================================
  ! Calculate time derivative of coefficients (velocity gauge)
  !=======================================================================
#include "rt_dg_integrator_derivative.f90"

  subroutine calibrate_dcdft_lcfo_static_hamiltonian(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
    use structures
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_stencil), intent(in) :: stencil
    type(s_scalar), intent(in) :: Vh
    type(s_scalar), intent(in) :: Vxc(system%nspin)
    type(s_scalar), intent(in) :: Vpsl
    real(8), intent(in) :: Ac_tot(3)

    integer :: iblk
    real(8) :: diff_local, diff_max
    real(8), allocatable :: Vh_buffer(:,:,:), Vpsl_buffer(:,:,:), Vxc_buffer(:,:,:,:)
    logical :: use_rank_buffered_potential

    if (.not. dg_frag%dc_lcfo_seed_basis_cleaned) return
    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%H_mat_core_blocks)) return
    if (.not. allocated(dg_frag%H_mat_kinetic_blocks)) return

    use_rank_buffered_potential = all(dg_frag%rank_buf_hi(:) >= dg_frag%rank_buf_lo(:))
    if (use_rank_buffered_potential) then
      allocate(Vh_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                         dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                         dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      allocate(Vpsl_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                           dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                           dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3)))
      allocate(Vxc_buffer(dg_frag%rank_buf_lo(1):dg_frag%rank_buf_hi(1), &
                          dg_frag%rank_buf_lo(2):dg_frag%rank_buf_hi(2), &
                          dg_frag%rank_buf_lo(3):dg_frag%rank_buf_hi(3), system%nspin))
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, dg_frag%mg, Vh, Vh_buffer)
      call copy_periodic_global_scalar_to_rank_buffer(dg_frag, dg_frag%mg, Vpsl, Vpsl_buffer)
      do iblk = 1, system%nspin
        call copy_periodic_global_scalar_to_rank_buffer(dg_frag, dg_frag%mg, Vxc(iblk), Vxc_buffer(:, :, :, iblk))
      end do
    end if

    do iblk = 1, size(dg_frag%H_mat_blocks)
      if (allocated(dg_frag%H_mat_blocks(iblk)%val) .and. allocated(dg_frag%H_mat_core_blocks(iblk)%val)) then
        dg_frag%H_mat_core_blocks(iblk)%val(:, :, :) = dg_frag%H_mat_blocks(iblk)%val(:, :, :)
      end if
      if (allocated(dg_frag%H_mat_kinetic_blocks(iblk)%val)) then
        dg_frag%H_mat_kinetic_blocks(iblk)%val(:, :, :) = 0.0d0
      end if
    end do

    if (use_rank_buffered_potential) then
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, &
                                          Vh_buffer, Vxc_buffer, Vpsl_buffer)
    else
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
    end if

    do iblk = 1, size(dg_frag%H_mat_blocks)
      if (allocated(dg_frag%H_mat_blocks(iblk)%val) .and. &
          allocated(dg_frag%H_mat_core_blocks(iblk)%val) .and. &
          allocated(dg_frag%H_mat_kinetic_blocks(iblk)%val)) then
        dg_frag%H_mat_kinetic_blocks(iblk)%val(:, :, :) = &
          dg_frag%H_mat_core_blocks(iblk)%val(:, :, :) - dg_frag%H_mat_blocks(iblk)%val(:, :, :)
        dg_frag%H_mat_blocks(iblk)%val(:, :, :) = dg_frag%H_mat_core_blocks(iblk)%val(:, :, :)
      end if
    end do
    dg_frag%H_blocks_include_nonlocal = .true.
    call rebuild_local_h_block_ids(dg_frag)

    if (use_rank_buffered_potential) then
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot, &
                                          Vh_buffer, Vxc_buffer, Vpsl_buffer)
    else
      call reconstruct_hamiltonian_matrix(dg_frag, system, stencil, Vh, Vxc, Vpsl, Ac_tot)
    end if
    diff_max = 0.0d0
    do iblk = 1, size(dg_frag%H_mat_blocks)
      if (allocated(dg_frag%H_mat_blocks(iblk)%val) .and. allocated(dg_frag%H_mat_core_blocks(iblk)%val)) then
        diff_local = maxval(abs(dg_frag%H_mat_blocks(iblk)%val(:, :, :) - &
                                dg_frag%H_mat_core_blocks(iblk)%val(:, :, :)))
        diff_max = max(diff_max, diff_local)
        dg_frag%H_mat_blocks(iblk)%val(:, :, :) = dg_frag%H_mat_core_blocks(iblk)%val(:, :, :)
      end if
    end do
    dg_frag%H_blocks_include_nonlocal = .true.
    call rebuild_local_h_block_ids(dg_frag)

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe13.5)') &
        '[DG-DCDFT-SEED] calibrated static H reference; max|H_static+V0-H_export|=', diff_max
      flush(6)
    end if
    if (use_rank_buffered_potential) then
      if (allocated(dg_frag%H_ref_Vh_buffer)) deallocate(dg_frag%H_ref_Vh_buffer)
      if (allocated(dg_frag%H_ref_Vxc_buffer)) deallocate(dg_frag%H_ref_Vxc_buffer)
      allocate(dg_frag%H_ref_Vh_buffer(lbound(Vh_buffer, 1):ubound(Vh_buffer, 1), &
                                       lbound(Vh_buffer, 2):ubound(Vh_buffer, 2), &
                                       lbound(Vh_buffer, 3):ubound(Vh_buffer, 3)))
      allocate(dg_frag%H_ref_Vxc_buffer(lbound(Vxc_buffer, 1):ubound(Vxc_buffer, 1), &
                                        lbound(Vxc_buffer, 2):ubound(Vxc_buffer, 2), &
                                        lbound(Vxc_buffer, 3):ubound(Vxc_buffer, 3), &
                                        lbound(Vxc_buffer, 4):ubound(Vxc_buffer, 4)))
      dg_frag%H_ref_Vh_buffer(:, :, :) = Vh_buffer(:, :, :)
      dg_frag%H_ref_Vxc_buffer(:, :, :, :) = Vxc_buffer(:, :, :, :)
      dg_frag%H_delta_reference_valid = .true.
    else
      dg_frag%H_delta_reference_valid = .false.
    end if
    if (allocated(Vh_buffer)) deallocate(Vh_buffer)
    if (allocated(Vpsl_buffer)) deallocate(Vpsl_buffer)
    if (allocated(Vxc_buffer)) deallocate(Vxc_buffer)
  end subroutine calibrate_dcdft_lcfo_static_hamiltonian

  subroutine refresh_buffer_wannier_flux_seed_from_current_hamiltonian(dg_frag, label)
    use communication, only: comm_is_root, comm_summation, comm_get_max
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in), optional :: label

    integer :: ifrag, i_local, ispin, iblk, nw, nbf
    real(8), allocatable :: h_dg(:,:), hw_tmp(:,:)
    real(8) :: stat_local(2), stat_global(2), diff_local, diff_max
    real(8) :: diff_local_arr(1), diff_global_arr(1)
    character(128) :: label_use

    if (.not. dg_frag%has_buffer_periodic_wannier_basis) return
    if (.not. allocated(dg_frag%buffer_wannier_coef)) return
    if (.not. allocated(dg_frag%buffer_wannier_h_flux)) return
    if (.not. allocated(dg_frag%buffer_wannier_nkeep)) return
    if (.not. allocated(dg_frag%H_block_map)) return
    if (.not. allocated(dg_frag%H_mat_blocks)) return

    label_use = '[DG-BPW-SCF-SEED]'
    if (present(label)) label_use = trim(label)

    stat_local(:) = 0.0d0
    stat_global(:) = 0.0d0
    diff_max = 0.0d0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = ifrag - dg_frag%ifrag_start + 1
      if (i_local < 1 .or. i_local > size(dg_frag%buffer_wannier_nkeep)) cycle
      nw = dg_frag%buffer_wannier_nkeep(i_local)
      if (nw <= 0) cycle
      do ispin = 1, dg_frag%nspin
        nbf = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag, &
                  size(dg_frag%buffer_wannier_coef, 1))
        if (nbf <= 0) cycle
        iblk = find_matrix_block(dg_frag%H_block_map, ifrag, ifrag)
        if (iblk <= 0 .or. iblk > size(dg_frag%H_mat_blocks)) cycle
        if (.not. allocated(dg_frag%H_mat_blocks(iblk)%val)) cycle
        allocate(h_dg(nbf,nbf), hw_tmp(nbf,nw))
        h_dg(1:nbf,1:nbf) = dg_frag%H_mat_blocks(iblk)%val(1:nbf,1:nbf,ispin)
        hw_tmp(1:nbf,1:nw) = matmul(h_dg(1:nbf,1:nbf), &
          dg_frag%buffer_wannier_coef(1:nbf,1:nw,ispin,i_local))
        if (ispin == 1) then
          diff_local = maxval(abs(dg_frag%buffer_wannier_h_flux(1:nw,1:nw,i_local) - &
            matmul(transpose(dg_frag%buffer_wannier_coef(1:nbf,1:nw,ispin,i_local)), &
                   hw_tmp(1:nbf,1:nw))))
          diff_max = max(diff_max, diff_local)
          dg_frag%buffer_wannier_h_flux(1:nw,1:nw,i_local) = &
            matmul(transpose(dg_frag%buffer_wannier_coef(1:nbf,1:nw,ispin,i_local)), &
                   hw_tmp(1:nbf,1:nw))
        end if
        stat_local(1) = stat_local(1) + 1.0d0
        stat_local(2) = stat_local(2) + dble(nw)
        deallocate(h_dg, hw_tmp)
      end do
    end do

    call comm_summation(stat_local, stat_global, 2, dg_frag%icomm)
    diff_local_arr(1) = diff_max
    call comm_get_max(diff_local_arr, diff_global_arr, 1, dg_frag%icomm)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,3(a,1pe13.5))') trim(label_use), &
        ' projected_blocks=', stat_global(1), &
        ' total_wannier=', stat_global(2), &
        ' max|Hproj-Hseed|=', diff_global_arr(1)
      flush(6)
    end if

    call initialize_coefficients_from_buffer_wannier_flux(dg_frag)
  end subroutine refresh_buffer_wannier_flux_seed_from_current_hamiltonian

  subroutine diagnose_dcdft_lcfo_seed_stationarity(dg_frag, system, mg, ppg, Ac_tot, label)
    use structures
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: mg
    type(s_pp_grid), intent(in) :: ppg
    real(8), intent(in) :: Ac_tot(3)
    character(*), intent(in), optional :: label

    integer, parameter :: max_probe = 512
    integer :: ispin, nocc, nprobe, iprobe, istate, irow, local_idx
    integer :: iblk_idx, iblk, ifrag_row, ifrag_col, n_needed, k
    integer :: max_rel_state
    integer :: sample_state(max_probe)
    integer, allocatable :: needed_row_ids(:)
    logical, allocatable :: row_needed(:)
    real(8), allocatable :: occ_weight(:)
    complex(8), allocatable :: coef_probe(:,:), coef_needed(:,:), hcoef_probe(:,:)
    real(8) :: res2_local, res2_global, c2_local, c2_global
    real(8) :: rel_res, max_rel_res, eps_absmax, eps_sample
    real(8) :: h2_local, h2_global, max_h_norm
    real(8) :: rayleigh_local, rayleigh_global, rayleigh
    real(8) :: max_rayleigh_abs, max_rayleigh_eps_abs
    complex(8) :: res_val
    real(8) :: rsend(4), rrecv(4)
    character(len=64) :: out_label

    if (.not. allocated(dg_frag%coef)) return
    if (.not. allocated(dg_frag%H_mat_blocks)) return
    if (.not. allocated(dg_frag%esp)) return
    if (dg_frag%nstate_tot <= 0) return
    allocate(occ_weight(max(1, dg_frag%nstate_tot)))
    max_rel_res = 0.0d0
    eps_absmax = 0.0d0
    max_h_norm = 0.0d0
    max_rayleigh_abs = 0.0d0
    max_rayleigh_eps_abs = 0.0d0
    max_rel_state = 0
    call rebuild_local_h_block_ids(dg_frag)
    if (.not. allocated(dg_frag%H_local_block_ids)) then
      deallocate(occ_weight)
      return
    end if

    do ispin = 1, dg_frag%nspin
      call build_hres_needed_rows(ispin)
      if (n_needed <= 0) cycle
      call get_dg_spin_occ_info(dg_frag, system, ispin, occ_weight, nocc)
      if (nocc <= 0) cycle
      nprobe = 0
      if (nocc <= max_probe) then
        do k = 1, nocc
          call append_probe_state(k)
        end do
      else
        call append_probe_state(1)
        call append_probe_state(max(1, (nocc + 1) / 2))
        call append_probe_state(nocc)
      end if
      do iprobe = 1, nprobe
        istate = sample_state(iprobe)
        if (allocated(dg_frag%esp)) then
          if (istate <= size(dg_frag%esp, 1) .and. ispin <= size(dg_frag%esp, 2)) then
            eps_sample = dg_frag%esp(istate, ispin)
            eps_absmax = max(eps_absmax, abs(eps_sample))
          end if
        end if
        call ensure_hres_work_arrays(n_needed)
        coef_probe(:, :) = (0.0d0, 0.0d0)
        coef_needed(:, :) = (0.0d0, 0.0d0)
        call fetch_remote_coef_rows(dg_frag, ispin, needed_row_ids, coef_needed, istate, istate)
        do k = 1, n_needed
          if (needed_row_ids(k) < 1 .or. needed_row_ids(k) > size(coef_probe, 1)) cycle
          coef_probe(needed_row_ids(k), 1) = coef_needed(k, 1)
        end do
        if (.not. allocated(hcoef_probe)) then
          allocate(hcoef_probe(dg_frag%n_mat_max, 1))
        else if (size(hcoef_probe, 1) /= dg_frag%n_mat_max) then
          deallocate(hcoef_probe)
          allocate(hcoef_probe(dg_frag%n_mat_max, 1))
        end if
        hcoef_probe(:, :) = (0.0d0, 0.0d0)
        call apply_matrix_blocks_batch(dg_frag, dg_frag%H_mat_blocks, ispin, &
                                       coef_probe(:, 1:1), hcoef_probe(:, 1:1), &
                                       dg_frag%H_local_block_ids)
        res2_local = 0.0d0
        c2_local = 0.0d0
        h2_local = 0.0d0
        rayleigh_local = 0.0d0
        do irow = 1, dg_frag%n_mat_max
          if (allocated(dg_frag%coef_owner)) then
            if (irow > size(dg_frag%coef_owner, 1)) cycle
            if (dg_frag%coef_owner(irow, ispin) /= dg_frag%id) cycle
            local_idx = dg_frag%coef_global_to_local(irow, ispin)
            if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
          else
            local_idx = irow
            if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
          end if
          res_val = hcoef_probe(irow, 1) - dg_frag%esp(istate, ispin) * coef_probe(irow, 1)
          res2_local = res2_local + abs(res_val)**2
          h2_local = h2_local + abs(hcoef_probe(irow, 1))**2
          c2_local = c2_local + abs(dg_frag%coef(local_idx, istate, ispin))**2
          rayleigh_local = rayleigh_local + real(conjg(coef_probe(irow, 1)) * hcoef_probe(irow, 1), 8)
        end do
        rsend(1) = res2_local
        rsend(2) = c2_local
        rsend(3) = h2_local
        rsend(4) = rayleigh_local
        call comm_summation(rsend, rrecv, 4, dg_frag%icomm)
        res2_global = rrecv(1)
        c2_global = rrecv(2)
        h2_global = rrecv(3)
        rayleigh_global = rrecv(4)
        if (c2_global > 0.0d0) then
          rel_res = sqrt(max(0.0d0, res2_global) / c2_global)
          rayleigh = rayleigh_global / c2_global
        else
          rel_res = 0.0d0
          rayleigh = 0.0d0
        end if
        if (rel_res > max_rel_res) then
          max_rel_res = rel_res
          max_rel_state = istate
        end if
        max_h_norm = max(max_h_norm, sqrt(max(0.0d0, h2_global) / max(1.0d-300, c2_global)))
        max_rayleigh_abs = max(max_rayleigh_abs, abs(rayleigh))
        max_rayleigh_eps_abs = max(max_rayleigh_eps_abs, abs(rayleigh - eps_sample))
      end do
    end do

    if (comm_is_root(dg_frag%id)) then
      out_label = '[DG-DCDFT-SEED-HRES]'
      if (present(label)) then
        if (len_trim(label) > 0) out_label = trim(label)
      end if
      write(*,'(1x,a,1x,a,1pe13.5,a,i0,a,1pe13.5,a,1pe13.5,a,1pe13.5,a,1pe13.5,a,3(1x,1pe13.5))') &
        trim(out_label), 'sampled ||(H-e)C||/||C||=', max_rel_res, &
        ' state=', max_rel_state, ' eps_absmax=', eps_absmax, ' h_absmax=', max_h_norm, &
        ' rayleigh_absmax=', max_rayleigh_abs, ' rayleigh_minus_eps_absmax=', max_rayleigh_eps_abs, &
        ' Ac=', Ac_tot
      flush(6)
    end if

    if (allocated(coef_probe)) deallocate(coef_probe)
    if (allocated(coef_needed)) deallocate(coef_needed)
    if (allocated(hcoef_probe)) deallocate(hcoef_probe)
    if (allocated(needed_row_ids)) deallocate(needed_row_ids)
    if (allocated(row_needed)) deallocate(row_needed)
    deallocate(occ_weight)

  contains

    subroutine ensure_hres_work_arrays(nrow_needed)
      integer, intent(in) :: nrow_needed

      if (.not. allocated(coef_probe)) then
        allocate(coef_probe(dg_frag%n_mat_max, 1))
      else if (size(coef_probe, 1) /= dg_frag%n_mat_max) then
        deallocate(coef_probe)
        allocate(coef_probe(dg_frag%n_mat_max, 1))
      end if
      if (.not. allocated(coef_needed)) then
        allocate(coef_needed(nrow_needed, 1))
      else if (size(coef_needed, 1) /= nrow_needed) then
        deallocate(coef_needed)
        allocate(coef_needed(nrow_needed, 1))
      end if
    end subroutine ensure_hres_work_arrays

    subroutine build_hres_needed_rows(ispin_current)
      integer, intent(in) :: ispin_current

      if (.not. allocated(row_needed)) then
        allocate(row_needed(max(1, dg_frag%n_mat_max)))
      else if (size(row_needed) /= max(1, dg_frag%n_mat_max)) then
        deallocate(row_needed)
        allocate(row_needed(max(1, dg_frag%n_mat_max)))
      end if
      row_needed(:) = .false.
      do iblk_idx = 1, size(dg_frag%H_local_block_ids)
        iblk = dg_frag%H_local_block_ids(iblk_idx)
        if (iblk < 1 .or. iblk > size(dg_frag%H_mat_blocks)) cycle
        ifrag_row = dg_frag%H_mat_blocks(iblk)%ifrag_row
        ifrag_col = dg_frag%H_mat_blocks(iblk)%ifrag_col
        call mark_fragment_rows(ifrag_row, ispin_current)
        call mark_fragment_rows(ifrag_col, ispin_current)
      end do
      n_needed = count(row_needed)
      if (allocated(needed_row_ids)) deallocate(needed_row_ids)
      allocate(needed_row_ids(max(1, n_needed)))
      n_needed = 0
      do irow = 1, size(row_needed)
        if (.not. row_needed(irow)) cycle
        n_needed = n_needed + 1
        needed_row_ids(n_needed) = irow
      end do
    end subroutine build_hres_needed_rows

    subroutine mark_fragment_rows(ifrag, ispin_current)
      integer, intent(in) :: ifrag, ispin_current
      integer :: ib, gid

      if (ifrag < 1 .or. ifrag > dg_frag%n_frag) return
      if (ispin_current < 1 .or. ispin_current > dg_frag%nspin) return
      do ib = 1, min(dg_frag%n_basis(ifrag, ispin_current), size(dg_frag%index_basis, 1))
        gid = dg_frag%index_basis(ib, ifrag, ispin_current)
        if (gid < 1 .or. gid > size(row_needed)) cycle
        row_needed(gid) = .true.
      end do
    end subroutine mark_fragment_rows

    subroutine append_probe_state(candidate)
      integer, intent(in) :: candidate
      integer :: i
      logical :: duplicate

      if (candidate < 1 .or. candidate > nocc) return
      duplicate = .false.
      do i = 1, nprobe
        if (sample_state(i) == candidate) duplicate = .true.
      end do
      if (duplicate) return
      if (nprobe >= max_probe) return
      nprobe = nprobe + 1
      sample_state(nprobe) = candidate
    end subroutine append_probe_state

  end subroutine diagnose_dcdft_lcfo_seed_stationarity

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
    if (allocated(dg_frag%H_mat_core_blocks)) then
      do i = 1, size(dg_frag%H_mat_core_blocks)
        if (allocated(dg_frag%H_mat_core_blocks(i)%val)) deallocate(dg_frag%H_mat_core_blocks(i)%val)
      end do
      deallocate(dg_frag%H_mat_core_blocks)
    end if
    call invalidate_coef_exchange_cache(dg_frag)

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
    if (allocated(dg_frag%nl_projector_overlap)) deallocate(dg_frag%nl_projector_overlap)
    if (allocated(dg_frag%nl_projector_r_overlap)) deallocate(dg_frag%nl_projector_r_overlap)
    if (allocated(dg_frag%nl_projector_overlap_halo)) deallocate(dg_frag%nl_projector_overlap_halo)
    if (allocated(dg_frag%nl_projector_r_overlap_halo)) deallocate(dg_frag%nl_projector_r_overlap_halo)
    dg_frag%nl_projector_cache_nlma = 0
    dg_frag%has_nl_projector_cache = .false.
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
    if (allocated(dg_frag%nl_projector_overlap)) deallocate(dg_frag%nl_projector_overlap)
    if (allocated(dg_frag%nl_projector_r_overlap)) deallocate(dg_frag%nl_projector_r_overlap)
    if (allocated(dg_frag%nl_projector_overlap_halo)) deallocate(dg_frag%nl_projector_overlap_halo)
    if (allocated(dg_frag%nl_projector_r_overlap_halo)) deallocate(dg_frag%nl_projector_r_overlap_halo)
    dg_frag%has_nl_projector_cache = .false.
    dg_frag%gradient_basis_cache_valid = .false.
    if (allocated(dg_frag%momentum_blocks)) then
      do i = 1, size(dg_frag%momentum_blocks)
        if (allocated(dg_frag%momentum_blocks(i)%val)) deallocate(dg_frag%momentum_blocks(i)%val)
      end do
      deallocate(dg_frag%momentum_blocks)
      dg_frag%n_momentum_blocks = 0
      dg_frag%momentum_blocks_include_dg_flux = .false.
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
    if (allocated(dg_frag%buffer_wannier_xi_flux_blocks)) then
      do i = 1, size(dg_frag%buffer_wannier_xi_flux_blocks)
        if (allocated(dg_frag%buffer_wannier_xi_flux_blocks(i)%val)) &
          deallocate(dg_frag%buffer_wannier_xi_flux_blocks(i)%val)
      end do
      deallocate(dg_frag%buffer_wannier_xi_flux_blocks)
      dg_frag%n_buffer_wannier_xi_flux_blocks = 0
      dg_frag%buffer_wannier_xi_flux_available = .false.
    end if
    if (allocated(dg_frag%buffer_wannier_xi_flux_local_block_ids)) &
      deallocate(dg_frag%buffer_wannier_xi_flux_local_block_ids)
    if (allocated(dg_frag%local_wannier_nbasis)) deallocate(dg_frag%local_wannier_nbasis)
    if (allocated(dg_frag%local_wannier_nproj)) deallocate(dg_frag%local_wannier_nproj)
    if (allocated(dg_frag%local_wannier_nkeep)) deallocate(dg_frag%local_wannier_nkeep)
    if (allocated(dg_frag%local_wannier_coef)) deallocate(dg_frag%local_wannier_coef)
    if (allocated(dg_frag%local_wannier_r)) deallocate(dg_frag%local_wannier_r)
    if (allocated(dg_frag%local_wannier_center)) deallocate(dg_frag%local_wannier_center)
    if (allocated(dg_frag%local_wannier_owner_fragment)) deallocate(dg_frag%local_wannier_owner_fragment)
    if (allocated(dg_frag%local_wannier_owned)) deallocate(dg_frag%local_wannier_owned)
    if (allocated(dg_frag%buffer_wannier_nkeep)) deallocate(dg_frag%buffer_wannier_nkeep)
    if (allocated(dg_frag%buffer_wannier_coef)) deallocate(dg_frag%buffer_wannier_coef)
    if (allocated(dg_frag%buffer_wannier_spread)) deallocate(dg_frag%buffer_wannier_spread)
    if (allocated(dg_frag%buffer_wannier_tail)) deallocate(dg_frag%buffer_wannier_tail)
    if (allocated(dg_frag%buffer_wannier_h_flux)) deallocate(dg_frag%buffer_wannier_h_flux)
    if (allocated(dg_frag%buffer_wannier_v)) deallocate(dg_frag%buffer_wannier_v)
    if (allocated(dg_frag%buffer_wannier_frag_center)) deallocate(dg_frag%buffer_wannier_frag_center)
    dg_frag%has_buffer_periodic_wannier_basis = .false.
    dg_frag%buffer_wannier_flux_seed_applied = .false.
    if (allocated(dg_frag%dg_wannier_nkeep)) deallocate(dg_frag%dg_wannier_nkeep)
    if (allocated(dg_frag%dg_wannier_global_ids)) deallocate(dg_frag%dg_wannier_global_ids)
    if (allocated(dg_frag%dg_wannier_owner_frag)) deallocate(dg_frag%dg_wannier_owner_frag)
    if (allocated(dg_frag%dg_wannier_ref_center)) deallocate(dg_frag%dg_wannier_ref_center)
    if (allocated(dg_frag%dg_wannier_basis_coef)) deallocate(dg_frag%dg_wannier_basis_coef)
    if (allocated(dg_frag%dg_wannier_h0_local)) deallocate(dg_frag%dg_wannier_h0_local)
    if (allocated(dg_frag%dg_wannier_xi_local)) deallocate(dg_frag%dg_wannier_xi_local)
    if (allocated(dg_frag%dg_wannier_coef)) deallocate(dg_frag%dg_wannier_coef)
    if (allocated(dg_frag%dg_wannier_neighbor_blocks)) then
      do i = 1, size(dg_frag%dg_wannier_neighbor_blocks)
        if (allocated(dg_frag%dg_wannier_neighbor_blocks(i)%h_flux)) &
          deallocate(dg_frag%dg_wannier_neighbor_blocks(i)%h_flux)
        if (allocated(dg_frag%dg_wannier_neighbor_blocks(i)%xi_flux)) &
          deallocate(dg_frag%dg_wannier_neighbor_blocks(i)%xi_flux)
      end do
      deallocate(dg_frag%dg_wannier_neighbor_blocks)
    end if
    if (allocated(dg_frag%dg_wannier_local_neighbor_block_ids)) &
      deallocate(dg_frag%dg_wannier_local_neighbor_block_ids)
    dg_frag%has_formal_dg_wannier_basis = .false.
    dg_frag%n_dg_wannier_neighbor_blocks = 0
    dg_frag%dg_wannier_blocks_gauge_consistent = .false.
    dg_frag%dg_wannier_uses_full_global_position = .false.
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
    if (allocated(dg_frag%current_density_grid_point_count)) deallocate(dg_frag%current_density_grid_point_count)
    if (allocated(dg_frag%current_density_grid_checksum)) deallocate(dg_frag%current_density_grid_checksum)
    if (allocated(dg_frag%current_valid_ix)) deallocate(dg_frag%current_valid_ix)
    if (allocated(dg_frag%current_valid_iy)) deallocate(dg_frag%current_valid_iy)
    if (allocated(dg_frag%current_valid_iz)) deallocate(dg_frag%current_valid_iz)
    if (allocated(dg_frag%current_valid_ixg)) deallocate(dg_frag%current_valid_ixg)
    if (allocated(dg_frag%current_valid_iyg)) deallocate(dg_frag%current_valid_iyg)
    if (allocated(dg_frag%current_valid_izg)) deallocate(dg_frag%current_valid_izg)
    if (allocated(dg_frag%runtime_neighbor_pair_cache)) deallocate(dg_frag%runtime_neighbor_pair_cache)
    if (allocated(dg_frag%runtime_neighbor_frag_count)) deallocate(dg_frag%runtime_neighbor_frag_count)
    if (allocated(dg_frag%runtime_neighbor_frag_ids)) deallocate(dg_frag%runtime_neighbor_frag_ids)
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) deallocate(dg_frag%momentum_neighbor_pair_cache)
    if (allocated(dg_frag%momentum_neighbor_frag_count)) deallocate(dg_frag%momentum_neighbor_frag_count)
    if (allocated(dg_frag%momentum_neighbor_frag_ids)) deallocate(dg_frag%momentum_neighbor_frag_ids)
    if (allocated(dg_frag%flux_face_trace_cache)) deallocate(dg_frag%flux_face_trace_cache)
    if (allocated(dg_frag%density_phi_block_cache)) deallocate(dg_frag%density_phi_block_cache)
    if (allocated(dg_frag%density_phi_block_count)) deallocate(dg_frag%density_phi_block_count)
    dg_frag%density_phi_block_size = 0
    dg_frag%density_phi_block_cache_valid = .false.
    if (allocated(dg_frag%density_phase_block_cache)) deallocate(dg_frag%density_phase_block_cache)
    dg_frag%density_phase_block_size = 0
    dg_frag%density_phase_block_npw = 0
    dg_frag%density_phase_block_cache_valid = .false.
    if (allocated(dg_frag%stage_Vh_buffer)) deallocate(dg_frag%stage_Vh_buffer)
    if (allocated(dg_frag%stage_Vpsl_buffer)) deallocate(dg_frag%stage_Vpsl_buffer)
    if (allocated(dg_frag%stage_Vxc_buffer)) deallocate(dg_frag%stage_Vxc_buffer)
    if (allocated(dg_frag%stage_gx_map)) deallocate(dg_frag%stage_gx_map)
    if (allocated(dg_frag%stage_gy_map)) deallocate(dg_frag%stage_gy_map)
    if (allocated(dg_frag%stage_gz_map)) deallocate(dg_frag%stage_gz_map)
    dg_frag%stage_vpsl_buffer_valid = .false.
    dg_frag%stage_map_valid = .false.
    if (allocated(dg_frag%jxyz_tot)) deallocate(dg_frag%jxyz_tot)
    if (allocated(dg_frag%phi_frag)) deallocate(dg_frag%phi_frag)
    if (allocated(dg_frag%phi_frag_has_seed_buffer)) deallocate(dg_frag%phi_frag_has_seed_buffer)
    if (allocated(dg_frag%rk_alpha)) deallocate(dg_frag%rk_alpha)
    if (allocated(dg_frag%rk_beta)) deallocate(dg_frag%rk_beta)
    if (allocated(dg_frag%rk_gamma)) deallocate(dg_frag%rk_gamma)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
    if (allocated(dg_frag%eigenvalues)) deallocate(dg_frag%eigenvalues)
    if (allocated(dg_frag%basis_overlap)) deallocate(dg_frag%basis_overlap)
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
    dg_frag%mixed_basis_identity_raw = .false.
    if (allocated(dg_frag%wpw_reduced_dim)) deallocate(dg_frag%wpw_reduced_dim)
    if (allocated(dg_frag%wpw_reduced_nself)) deallocate(dg_frag%wpw_reduced_nself)
    if (allocated(dg_frag%wpw_reduced_nkeep)) deallocate(dg_frag%wpw_reduced_nkeep)
    if (allocated(dg_frag%wpw_reduced_ndrop)) deallocate(dg_frag%wpw_reduced_ndrop)
    if (allocated(dg_frag%wpw_reduced_nraw)) deallocate(dg_frag%wpw_reduced_nraw)
    if (allocated(dg_frag%wpw_reduced_H)) deallocate(dg_frag%wpw_reduced_H)
    if (allocated(dg_frag%wpw_reduced_S)) deallocate(dg_frag%wpw_reduced_S)
    if (allocated(dg_frag%wpw_reduced_transform)) deallocate(dg_frag%wpw_reduced_transform)
    if (allocated(dg_frag%wpw_reduced_eval)) deallocate(dg_frag%wpw_reduced_eval)
    if (allocated(dg_frag%wpw_reduced_evec)) deallocate(dg_frag%wpw_reduced_evec)
    if (allocated(dg_frag%coef_wpw_self)) deallocate(dg_frag%coef_wpw_self)
    if (allocated(dg_frag%coef_wpw_neighbor_reduced)) deallocate(dg_frag%coef_wpw_neighbor_reduced)
    if (allocated(dg_frag%wpw_reproject_prev_coef)) deallocate(dg_frag%wpw_reproject_prev_coef)
    dg_frag%wpw_reduced_ready = .false.
    dg_frag%wpw_reduced_max_dim = 0
    dg_frag%wpw_reduced_coef_initialized = .false.
    dg_frag%wpw_reproject_prev_valid = .false.
    if (allocated(dg_frag%S_mat_frag_pw)) deallocate(dg_frag%S_mat_frag_pw)
    if (allocated(dg_frag%H_mat_frag_pw)) deallocate(dg_frag%H_mat_frag_pw)
    if (allocated(dg_frag%P_mat_frag_pw)) deallocate(dg_frag%P_mat_frag_pw)
    if (allocated(dg_frag%fp_local_row_ids)) deallocate(dg_frag%fp_local_row_ids)
    if (allocated(dg_frag%fp_local_pw_ids)) deallocate(dg_frag%fp_local_pw_ids)
    if (allocated(dg_frag%S_mat_frag_pw_local)) deallocate(dg_frag%S_mat_frag_pw_local)
    if (allocated(dg_frag%H_mat_frag_pw_local)) deallocate(dg_frag%H_mat_frag_pw_local)
    if (allocated(dg_frag%P_mat_frag_pw_local)) deallocate(dg_frag%P_mat_frag_pw_local)
    if (allocated(dg_frag%H_mat_pw_diag)) deallocate(dg_frag%H_mat_pw_diag)
    if (allocated(dg_frag%H_mat_pw)) deallocate(dg_frag%H_mat_pw)
    if (allocated(dg_frag%wpw_window_box_lo)) deallocate(dg_frag%wpw_window_box_lo)
    if (allocated(dg_frag%wpw_window_box_hi)) deallocate(dg_frag%wpw_window_box_hi)
    if (allocated(dg_frag%wpw_chi)) deallocate(dg_frag%wpw_chi)
    if (allocated(dg_frag%wpw_grad_chi)) deallocate(dg_frag%wpw_grad_chi)
    if (allocated(dg_frag%wpw_S_pp_blocks)) deallocate(dg_frag%wpw_S_pp_blocks)
    if (allocated(dg_frag%wpw_T_pp_volume_blocks)) deallocate(dg_frag%wpw_T_pp_volume_blocks)
    if (allocated(dg_frag%wpw_T_pp_interface_blocks)) deallocate(dg_frag%wpw_T_pp_interface_blocks)
    dg_frag%has_wpw_window = .false.
    dg_frag%wpw_pp_blocks_ready = .false.
    if (dg_frag%icomm_frag /= COMM_GROUP_NULL) then
      call comm_free_group(dg_frag%icomm_frag)
      dg_frag%icomm_frag = COMM_GROUP_NULL
    end if
    
    ! RI/DF data deallocations
    if (dg_frag%use_hse_ri) call finalize_hse_ri_data()
    
  end subroutine finalize_dg_fragment_rt

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

  subroutine initialize_dg_state_block_coef_mode(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    character(16) :: env_value
    integer :: env_status, nworker, rank_in_frag, base, extra

    dg_frag%coef_state_block_mode = .false.
    dg_frag%coef_state_start = 1
    dg_frag%coef_state_end = dg_frag%nstate_tot
    dg_frag%coef_nstate_local = dg_frag%nstate_tot

    env_value = ''
    call get_environment_variable('SALMON_DG_STATE_BLOCK_COEF', env_value, status=env_status)
    if (env_status == 0) then
      select case(trim(adjustl(env_value)))
      case('0','n','N','no','NO','false','FALSE','off','OFF')
        dg_frag%coef_state_block_mode = .false.
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        dg_frag%coef_state_block_mode = .true.
      end select
    end if
    if (.not. dg_frag%coef_state_block_mode) return

    if (.not. dg_frag%parallel_mode_orbital) then
      stop 'DG-Fragment RT: state-block coefficients require orbital fragment mode'
    end if
    nworker = max(1, dg_frag%isize_frag)
    rank_in_frag = max(0, min(dg_frag%id_frag, nworker - 1))
    base = dg_frag%nstate_tot / nworker
    extra = mod(dg_frag%nstate_tot, nworker)
    if (rank_in_frag < extra) then
      dg_frag%coef_state_start = rank_in_frag * (base + 1) + 1
      dg_frag%coef_state_end = dg_frag%coef_state_start + base
    else
      dg_frag%coef_state_start = extra * (base + 1) + (rank_in_frag - extra) * base + 1
      dg_frag%coef_state_end = dg_frag%coef_state_start + base - 1
    end if
    if (dg_frag%coef_state_start > dg_frag%nstate_tot) then
      dg_frag%coef_state_start = 1
      dg_frag%coef_state_end = 0
    else
      dg_frag%coef_state_end = min(dg_frag%coef_state_end, dg_frag%nstate_tot)
    end if
    dg_frag%coef_nstate_local = max(0, dg_frag%coef_state_end - dg_frag%coef_state_start + 1)
    if (dg_frag%id == 0) then
      write(*,'(1x,a)') '[INFO] DG coefficient mode: state-block local columns'
      flush(6)
    end if
  end subroutine initialize_dg_state_block_coef_mode


end module rt_dg_fragment
