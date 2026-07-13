!
!  Copyright 2019-2024 SALMON developers
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

! DC-LCFO method [Phys. Rev. B 95, 045106 (2017).]

#include "config.h"
module lcfo_flux
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dc_fragment_geometry, only: get_fragment_domain
  use lcfo_wannier_sawf, only: t_sawf_projection_channel, t_sawf_symop, &
    build_sawf_spd_projection_map, build_sawf_wannier_representation, &
    sawf_real_harmonic_value, sawf_spd_projection_count, sawf_projection_shell_lmax, &
    write_sawf_projection_block, &
    load_sawf_symmetry_auto, load_sawf_symmetry_file, build_sawf_operation_product_table
  use lcfo_wannier_sawf_band, only: map_sawf_periodic_grid_point, &
    locate_sawf_fragment_point, validate_sawf_fragment_tiling, &
    validate_sawf_fragment_symmetry_map, build_sawf_fragment_buffer_point_map, &
    accumulate_sawf_dmn_band_blocks, validate_sawf_dmn_band, &
    validate_sawf_seed_header, validate_sawf_seed_basis_metadata, &
    diagnose_sawf_fragment_basis_closure, validate_sawf_operation_set_products
  use lcfo_wannier_sawf_basis, only: t_sawf_closed_basis,append_sawf_mapped_basis, &
    build_sawf_closed_core_buffer_basis,clear_sawf_closed_basis
  use lcfo_wannier_sawf_dmn, only: t_sawf_dmn_writer, begin_sawf_dmn, &
    append_sawf_dmn_operation, finish_sawf_dmn, abort_sawf_dmn
  use lcfo_wannier_sawf_collective, only: reduce_sawf_fragment_alignment_failure, &
    validate_sawf_density_contribution,assemble_sawf_density_unique
  use lcfo_wannier_sawf_templates, only: validate_sawf_actual_group_operation, &
    build_sawf_environment_orbits, build_sawf_supercell_fingerprint, &
    build_sawf_local_environment_fingerprint, select_sawf_environment_materialization
  use lcfo_wannier_sawf_templates, only: validate_sawf_structure_class
  use lcfo_wannier_sawf_templates, only: build_sawf_file_content_digest
  use lcfo_wannier_sawf_templates, only: measure_sawf_vacuum_occupancy
  use lcfo_wannier_sawf_templates, only: sawf_closest_periodic_cartesian
  use lcfo_wannier_sawf_orchestrator, only: t_sawf_environment_receipt,t_sawf_seed_bundle, &
    build_sawf_environment_execution_plan,build_sawf_seed_bundles,select_sawf_environment_stabilizer, &
    complete_sawf_seed_bundle,propagate_sawf_representative_receipts,validate_sawf_environment_receipts
  use lcfo_wannier_sawf_seed, only: select_sawf_local_complete_shells, &
    build_sawf_local_seed_matrices,write_sawf_local_eig_amn, &
    solve_sawf_local_generalized_eigensystem,read_sawf_nnkp_neighbors, &
    write_sawf_local_eig_amn_mmn,restrict_sawf_stabilizer_representation, &
    build_sawf_local_band_representation
  use lcfo_wannier_sawf_win, only: activate_sawf_win, deactivate_sawf_win, &
    t_atomic_win_writer, begin_atomic_win, finish_atomic_win, abort_atomic_win, &
    write_sawf_local_preprocess_win,write_sawf_atomic_text
  use lcfo_wannier_command, only: select_wannier90_command, execute_wannier90_command, &
    is_wannier90_reuse_command, is_wannier90_export_only_command, &
    is_wannier90_import_only_command, cache_resolved_wannier90_command, &
    get_cached_wannier90_command
  implicit none

  private
  public :: dc_lcfo_flux, dc_lcfo_wannier_import_only, dc_lcfo_wannier_import_only_requested

  character(48),parameter :: binfile_wf = "wavefunctions.bin", &
  &                          binfile_wf_wannier_seed = "wavefunctions_wannier_seed.bin", &
  &                          binfile_rg = "rgrid_index.bin", &
  &                          binfile_bf = "basis_functions.bin", &
  &                          binfile_bfb = "basis_functions_buffer.bin", &
  &                          binfile_hl = "hamiltonian_local.bin", &
  &                          binfile_hc = "hamiltonian_flux_components.bin", &
  &                          binfile_hcw = "hamiltonian_flux_weak_components.bin", &
  &                          binfile_vl = "velocity_local.bin", &
  &                          binfile_xfl = "position_flux_local.bin", &
  &                          binfile_lw = "local_wannier_basis.bin", &
  &                          binfile_bpw = "buffer_periodic_wannier_basis.bin", &
  &                          binfile_bpwt = "buffer_periodic_wannier_trace.bin", &
  &                          binfile_w90g = "wannier90_global_basis.bin", &
  &                          binfile_w90g_persistent = "wannier90_global_basis_persistent.bin", &
  &                          binfile_w90seed = "wannier_flux_eigen_seed.bin", &
  &                          binfile_w90seed_persistent = "wannier_flux_eigen_seed_persistent.bin", &
  &                          binfile_wcl = "wannier_cluster_partition.bin"
  integer, parameter :: basis_buffer_magic = -22022212
  integer, parameter :: basis_buffer_version = 2
  integer, parameter :: local_wannier_magic = -22022214
  integer, parameter :: local_wannier_version = 2
  integer, parameter :: buffer_periodic_wannier_magic = -22022215
  integer, parameter :: buffer_periodic_wannier_version = 3
  integer, parameter :: buffer_periodic_wannier_trace_magic = -22022218
  integer, parameter :: buffer_periodic_wannier_trace_version = 1
  integer, parameter :: wannier90_global_magic = -22022216
  integer, parameter :: wannier90_global_version = 2
  integer, parameter :: wannier_flux_eigen_seed_magic = -22022219
  integer, parameter :: wannier_flux_eigen_seed_version = 1
  integer, parameter :: wannier_cluster_magic = -22022217
  integer, parameter :: wannier_cluster_version = 1
  integer, parameter :: sawf_seed_provenance_magic = -22022220
  integer, parameter :: sawf_seed_provenance_version = 1
  integer, parameter :: sawf_closed_seed_magic = -22022221
  integer, parameter :: sawf_closed_seed_version = 1
  character(64), save :: current_sawf_seed_token = ''

  type :: t_wannier_symop
    integer :: owner_frag = 0
    integer :: rot(3,3) = 0
    real(8) :: origin_bohr(3) = 0.0d0
    real(8) :: tau_local(3) = 0.0d0
    real(8) :: atom_residual = 0.0d0
    character(32) :: label = ''
  end type t_wannier_symop

  integer, parameter :: sawf_target_cache_capacity=2

  type :: t_sawf_fragment_state_entry
    integer :: fragment=0
    integer :: shape(3)=0
    integer :: buffer_shape(3)=0,buffer_width(3)=0
    integer :: nspin=0,nstate_frag=0,nstate_tot=0,n_basis=0
    real(8), allocatable :: basis(:,:),buffer_basis(:,:)
    complex(8), allocatable :: states(:,:)
  end type t_sawf_fragment_state_entry

  type :: t_sawf_fragment_state_cache
    type(t_sawf_fragment_state_entry) :: source
    type(t_sawf_fragment_state_entry) :: target(sawf_target_cache_capacity)
    integer :: next_target_slot=1
    integer :: source_seed_reads=0,source_reconstructions=0
    integer :: target_seed_reads=0,target_reconstructions=0,target_cache_hits=0
  end type t_sawf_fragment_state_cache

contains

  pure real(8) function wrap_periodic_delta(delta, cell_length) result(delta_wrapped)
    implicit none
    real(8), intent(in) :: delta, cell_length

    if(cell_length > 0.0d0) then
      delta_wrapped = delta - dnint(delta / cell_length) * cell_length
    else
      delta_wrapped = delta
    end if
  end function wrap_periodic_delta

  pure real(8) function local_distance2(delta) result(distance2)
    implicit none
    real(8), intent(in) :: delta(3)

    distance2 = delta(1) * delta(1) + delta(2) * delta(2) + delta(3) * delta(3)
  end function local_distance2

  pure integer function integer_det3(rot) result(det)
    implicit none
    integer, intent(in) :: rot(3,3)

    det = rot(1,1) * (rot(2,2) * rot(3,3) - rot(2,3) * rot(3,2)) &
      - rot(1,2) * (rot(2,1) * rot(3,3) - rot(2,3) * rot(3,1)) &
      + rot(1,3) * (rot(2,1) * rot(3,2) - rot(2,2) * rot(3,1))
  end function integer_det3

  pure function matmul_int_real(rot, vec) result(out)
    implicit none
    integer, intent(in) :: rot(3,3)
    real(8), intent(in) :: vec(3)
    real(8) :: out(3)
    integer :: i

    do i=1,3
      out(i) = dble(rot(i,1)) * vec(1) + dble(rot(i,2)) * vec(2) + dble(rot(i,3)) * vec(3)
    end do
  end function matmul_int_real

  logical function dc_lcfo_wannier_import_only_requested() result(import_only)
    use communication, only: comm_bcast
    use parallelization, only: nproc_group_global, nproc_id_global
    use salmon_global, only: wannier90_command
    implicit none
    character(1024) :: environment_command, resolved_command, compiled_default
    integer :: env_status

    environment_command = ''; resolved_command = ''
#ifdef WANNIER90_EXECUTABLE_PATH
    compiled_default = WANNIER90_EXECUTABLE_PATH
#else
    compiled_default = 'wannier90.x'
#endif
    if(nproc_id_global == 0) then
      call get_environment_variable('SALMON_WANNIER90_COMMAND',environment_command,status=env_status)
      if(env_status /= 0) environment_command = ''
    end if
    call comm_bcast(environment_command,nproc_group_global,0)
    call select_wannier90_command(wannier90_command,environment_command,compiled_default,resolved_command)
    call cache_resolved_wannier90_command(resolved_command)
    import_only = is_wannier90_import_only_command(resolved_command)
  end function dc_lcfo_wannier_import_only_requested

  subroutine dc_lcfo_wannier_import_only(dc)
    use communication, only: comm_sync_all
    use filesystem, only: get_filehandle
    use inputoutput, only: au_length_aa
    use salmon_global, only: base_directory, sysname, wannier_num_wann, wannier_projection, dg_wannier_symmetry_gauge
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer :: num_bands_chk, num_wann_chk, nstate_tot_file, nspin_file
    integer :: iunit, iw, istate, nstate_seed, position_available, nsym
    integer, allocatable :: owner_frag(:), bond_owner_frag(:)
    type(t_wannier_symop), allocatable :: symops(:)
    real(8), allocatable :: center_aa(:,:), center_bohr(:,:), owner_center_bohr(:,:), spread_aa2(:)
    real(8), allocatable :: esp_file(:,:), eval_seed(:,:)
    complex(8), allocatable :: v_matrix(:,:), v_direct(:,:), v_direct_global(:,:)
    complex(8), allocatable :: aa_global(:,:,:), aa_direct(:,:,:), seed_wannier_to_eigen(:,:)
    logical :: ok_position, ok_direct, direct_gauge_requested, direct_gauge_applied
    character(256) :: filename

    call read_dc_lcfo_esp_from_wavefunctions_import(dc, nstate_tot_file, nspin_file, esp_file)

    if(dc%id_tot == 0) then
      direct_gauge_requested = (trim(dg_wannier_symmetry_gauge) == 'direct_amn_bond_block' .or. &
        trim(dg_wannier_symmetry_gauge) == 'direct_amn_bond_global' .or. &
        trim(dg_wannier_symmetry_gauge) == 'direct_amn_global')
      direct_gauge_applied = .false.
      if(direct_gauge_requested .and. trim(dg_wannier_symmetry_gauge) /= 'direct_amn_global' .and. &
          .not. is_bond_center_projection_import(trim(wannier_projection))) then
        write(*,'(1x,a,a,a)') "[FATAL] dg_wannier_symmetry_gauge=", &
          trim(dg_wannier_symmetry_gauge), " requires wannier_projection='bond_centers'."
        stop "DC-LCFO Wannier import-only: direct AMN bond gauge requires bond_centers projection"
      end if
      call read_wannier90_checkpoint_transform_import(dc, num_bands_chk, num_wann_chk, &
        v_matrix, center_aa, spread_aa2)
      if(num_wann_chk /= wannier_num_wann) then
        write(*,'(1x,a,2(a,i0))') "[DC-LCFO-W90-IMPORT] checkpoint dimension mismatch:", &
          " chk_wann=", num_wann_chk, " expected_wann=", wannier_num_wann
        stop "DC-LCFO Wannier import-only: checkpoint dimension mismatch"
      end if
      if(num_bands_chk > nstate_tot_file) then
        write(*,'(1x,a,2(a,i0))') "[DC-LCFO-W90-IMPORT] wavefunction seed has too few bands:", &
          " chk_bands=", num_bands_chk, " seed_states=", nstate_tot_file
        stop "DC-LCFO Wannier import-only: insufficient wavefunction seed"
      end if

      allocate(owner_frag(num_wann_chk), center_bohr(3,num_wann_chk))
      center_bohr(1:3,1:num_wann_chk) = center_aa(1:3,1:num_wann_chk) / au_length_aa
      call wrap_wannier_centers_to_total_cell_import(dc, center_bohr, num_wann_chk)
      do iw=1,num_wann_chk
        owner_frag(iw) = find_owner_fragment_from_center_import(dc, center_bohr(1:3,iw))
      end do
      call rebalance_wannier_owner_fragments_import(dc, center_bohr, owner_frag, num_wann_chk)
      write(*,'(1x,a)') "[DC-LCFO-W90-IMPORT] owner source=wannier_centers"
      if(is_bond_center_projection_import(trim(wannier_projection))) then
        call build_bond_center_projection_map_import(dc, num_wann_chk, owner_center_bohr)
        allocate(bond_owner_frag(num_wann_chk))
        do iw=1,num_wann_chk
          bond_owner_frag(iw) = find_owner_fragment_from_center_import(dc, owner_center_bohr(1:3,iw))
        end do
        call rebalance_wannier_owner_fragments_import(dc, owner_center_bohr, bond_owner_frag, num_wann_chk)
        write(*,'(1x,a)') "[DC-LCFO-W90-IMPORT] projection source=bond_centers"
      end if
      call detect_wannier_fragment_symops(dc, nsym, symops)
      if(allocated(owner_center_bohr)) then
        write(*,'(1x,a)') "[DC-LCFO-W90-SYM] center diagnostic source=bond_center_seeds"
        call diagnose_fragment_wannier_center_symmetry(dc, owner_center_bohr, bond_owner_frag, num_wann_chk, nsym, symops)
      end if
      call diagnose_fragment_wannier_center_symmetry(dc, center_bohr, owner_frag, num_wann_chk, nsym, symops)
      call read_wannier90_global_rmn_gamma_block_import(dc, num_wann_chk, aa_global, ok_position)
      if(.not. ok_position) then
        if(direct_gauge_requested) then
          write(*,'(1x,a,a)') "[FATAL] requested direct AMN gauge requires Wannier90 position data: mode=", &
            trim(dg_wannier_symmetry_gauge)
          stop "DC-LCFO Wannier import-only: requested direct AMN gauge requires Wannier90 position data"
        end if
        allocate(aa_global(3,num_wann_chk,num_wann_chk))
        aa_global = (0d0,0d0)
      else
        call set_wannier_centers_from_position_diagonal(num_wann_chk, center_bohr, aa_global)
        do iw=1,num_wann_chk
          owner_frag(iw) = find_owner_fragment_from_center_import(dc, center_bohr(1:3,iw))
        end do
        call rebalance_wannier_owner_fragments_import(dc, center_bohr, owner_frag, num_wann_chk)
        write(*,'(1x,a)') "[DC-LCFO-W90-IMPORT] center source=AA_R diagonal"
        call diagnose_fragment_wannier_symmetry_representation(dc, center_bohr, owner_frag, num_wann_chk, &
          num_bands_chk, v_matrix, esp_file, nsym, symops, aa_global, ok_position)
        call diagnose_global_wannier_pbc_operator_symmetry(dc, center_bohr, num_wann_chk, &
          num_bands_chk, v_matrix, esp_file, nsym, symops, aa_global, ok_position)
        if(is_bond_center_projection_import(trim(wannier_projection))) then
          call read_wannier90_amn_direct_transform_import(dc, num_bands_chk, num_wann_chk, bond_owner_frag, v_direct, ok_direct)
          if(ok_direct) then
            write(*,'(1x,a)') "[DC-LCFO-W90-SYM] diagnostic basis=direct_amn_bond_projectors_block"
            call diagnose_fragment_wannier_symmetry_representation(dc, owner_center_bohr, bond_owner_frag, num_wann_chk, &
              num_bands_chk, v_direct, esp_file, nsym, symops, aa_global, .false.)
            call diagnose_global_wannier_pbc_operator_symmetry(dc, owner_center_bohr, num_wann_chk, &
              num_bands_chk, v_direct, esp_file, nsym, symops, aa_global, .false.)
            if(trim(dg_wannier_symmetry_gauge) == 'direct_amn_bond_block') then
              call transform_wannier_position_gauge(num_bands_chk, num_wann_chk, v_matrix, v_direct, aa_global, aa_direct)
              aa_global(1:3,1:num_wann_chk,1:num_wann_chk) = aa_direct(1:3,1:num_wann_chk,1:num_wann_chk)
              v_matrix(1:num_bands_chk,1:num_wann_chk) = v_direct(1:num_bands_chk,1:num_wann_chk)
              center_bohr(1:3,1:num_wann_chk) = owner_center_bohr(1:3,1:num_wann_chk)
              owner_frag(1:num_wann_chk) = bond_owner_frag(1:num_wann_chk)
              call set_wannier_centers_from_position_diagonal(num_wann_chk, center_bohr, aa_global)
              do iw=1,num_wann_chk
                owner_frag(iw) = find_owner_fragment_from_center_import(dc, center_bohr(1:3,iw))
              end do
              call rebalance_wannier_owner_fragments_import(dc, center_bohr, owner_frag, num_wann_chk)
              write(*,'(1x,a)') "[DC-LCFO-W90-IMPORT] basis source=direct_amn_bond_projectors_block"
              deallocate(aa_direct)
              direct_gauge_applied = .true.
            end if
            deallocate(v_direct)
          else if(trim(dg_wannier_symmetry_gauge) == 'direct_amn_bond_block') then
            stop "DC-LCFO Wannier import-only: requested direct_amn_bond_block but AMN read failed"
          end if
          call read_wannier90_amn_direct_transform_import(dc, num_bands_chk, num_wann_chk, bond_owner_frag, &
            v_direct_global, ok_direct, .false.)
          if(ok_direct) then
            write(*,'(1x,a)') "[DC-LCFO-W90-SYM] diagnostic basis=direct_amn_bond_projectors_global"
            call diagnose_fragment_wannier_symmetry_representation(dc, owner_center_bohr, bond_owner_frag, num_wann_chk, &
              num_bands_chk, v_direct_global, esp_file, nsym, symops, aa_global, .false.)
            call diagnose_global_wannier_pbc_operator_symmetry(dc, owner_center_bohr, num_wann_chk, &
              num_bands_chk, v_direct_global, esp_file, nsym, symops, aa_global, .false.)
            if(trim(dg_wannier_symmetry_gauge) == 'direct_amn_bond_global') then
              call transform_wannier_position_gauge(num_bands_chk, num_wann_chk, v_matrix, v_direct_global, aa_global, aa_direct)
              aa_global(1:3,1:num_wann_chk,1:num_wann_chk) = aa_direct(1:3,1:num_wann_chk,1:num_wann_chk)
              v_matrix(1:num_bands_chk,1:num_wann_chk) = v_direct_global(1:num_bands_chk,1:num_wann_chk)
              center_bohr(1:3,1:num_wann_chk) = owner_center_bohr(1:3,1:num_wann_chk)
              owner_frag(1:num_wann_chk) = bond_owner_frag(1:num_wann_chk)
              call set_wannier_centers_from_position_diagonal(num_wann_chk, center_bohr, aa_global)
              do iw=1,num_wann_chk
                owner_frag(iw) = find_owner_fragment_from_center_import(dc, center_bohr(1:3,iw))
              end do
              call rebalance_wannier_owner_fragments_import(dc, center_bohr, owner_frag, num_wann_chk)
              write(*,'(1x,a)') "[DC-LCFO-W90-IMPORT] basis source=direct_amn_bond_projectors_global"
              deallocate(aa_direct)
              direct_gauge_applied = .true.
            end if
            deallocate(v_direct_global)
          else if(trim(dg_wannier_symmetry_gauge) == 'direct_amn_bond_global') then
            stop "DC-LCFO Wannier import-only: requested direct_amn_bond_global but AMN read failed"
          end if
        end if
        if(trim(dg_wannier_symmetry_gauge) == 'direct_amn_global') then
          call read_wannier90_amn_direct_transform_import(dc, num_bands_chk, num_wann_chk, owner_frag, &
            v_direct_global, ok_direct, .false.)
          if(ok_direct) then
            call diagnose_wannier_transform_subspace_overlap(num_bands_chk, num_wann_chk, &
              v_matrix, v_direct_global, "W90_vs_direct_amn_global")
            call build_direct_wannier_buffer_position_from_lcfo_import(dc, num_bands_chk, num_wann_chk, &
              owner_frag, center_bohr, v_direct_global, aa_direct, ok_position)
            if(.not. ok_position) stop "DC-LCFO Wannier import-only: failed to build direct AMN buffer-local position"
            aa_global(1:3,1:num_wann_chk,1:num_wann_chk) = aa_direct(1:3,1:num_wann_chk,1:num_wann_chk)
            v_matrix(1:num_bands_chk,1:num_wann_chk) = v_direct_global(1:num_bands_chk,1:num_wann_chk)
            write(*,'(1x,a)') "[DC-LCFO-W90-IMPORT] position source=direct_amn_buffer_local_r"
            write(*,'(1x,a)') "[DC-LCFO-W90-IMPORT] basis source=direct_amn_projectors_global"
            deallocate(aa_direct, v_direct_global)
            direct_gauge_applied = .true.
          else
            stop "DC-LCFO Wannier import-only: requested direct_amn_global but AMN read failed"
          end if
        end if
        if(direct_gauge_requested .and. .not. direct_gauge_applied) &
          stop "DC-LCFO Wannier import-only: requested direct AMN gauge was not applied"
        if(trim(dg_wannier_symmetry_gauge) == 'local_inversion_position') then
          call symmetrize_fragment_wannier_position_import(dc, center_bohr, owner_frag, num_wann_chk, &
            num_bands_chk, v_matrix, nsym, symops, aa_global)
          call diagnose_fragment_wannier_center_symmetry(dc, center_bohr, owner_frag, num_wann_chk, nsym, symops)
        else
          write(*,'(1x,a,a)') "[DC-LCFO-W90-SYM] position sym mode=", trim(dg_wannier_symmetry_gauge)
          if(direct_gauge_applied) then
            if(trim(dg_wannier_symmetry_gauge) == 'direct_amn_global') then
              write(*,'(1x,a)') "[DC-LCFO-W90-SYM] direct_amn_global operator diagnostic deferred: "// &
                "unitary symmetry representation required"
            else
              call diagnose_global_wannier_pbc_operator_symmetry(dc, center_bohr, num_wann_chk, &
                num_bands_chk, v_matrix, esp_file, nsym, symops, aa_global, ok_position)
            end if
          end if
        end if
      end if
      if(allocated(symops)) deallocate(symops)
      position_available = merge(1, 0, ok_position)

      filename = trim(import_run_root_dir())//'data_dcdft/total/'//binfile_w90g
      call write_wannier90_global_basis_stream_import(filename, num_bands_chk, num_wann_chk, dc%n_frag, &
        owner_frag, center_bohr, spread_aa2, v_matrix, position_available, aa_global)
      write(*,'(1x,a,i0,a,a)') "[DC-LCFO-W90-IMPORT] wrote ", num_wann_chk, &
        " Wannier functions to ", trim(filename)
      filename = trim(import_run_root_dir())//'data_dcdft/total/'//binfile_w90g_persistent
      call write_wannier90_global_basis_stream_import(filename, num_bands_chk, num_wann_chk, dc%n_frag, &
        owner_frag, center_bohr, spread_aa2, v_matrix, position_available, aa_global)

      call build_wannier_flux_eigen_seed_from_transform(num_bands_chk, num_wann_chk, nstate_tot_file, &
        nspin_file, v_matrix, esp_file, seed_wannier_to_eigen, eval_seed, nstate_seed)
      filename = trim(import_run_root_dir())//'data_dcdft/total/'//binfile_w90seed
      call write_wannier_flux_eigen_seed_stream_import(filename, num_bands_chk, num_wann_chk, nstate_seed, &
        nspin_file, dc%n_frag, eval_seed, seed_wannier_to_eigen)
      write(*,'(1x,a,i0,a,i0,a,a)') "[DC-LCFO-W90-IMPORT] wrote Flux-LCFO eigen seed: states=", &
        nstate_seed, " wann=", num_wann_chk, " file=", trim(filename)
      filename = trim(import_run_root_dir())//'data_dcdft/total/'//binfile_w90seed_persistent
      call write_wannier_flux_eigen_seed_stream_import(filename, num_bands_chk, num_wann_chk, nstate_seed, &
        nspin_file, dc%n_frag, eval_seed, seed_wannier_to_eigen)
      deallocate(seed_wannier_to_eigen, eval_seed)

      if(allocated(owner_center_bohr)) deallocate(owner_center_bohr)
      if(allocated(bond_owner_frag)) deallocate(bond_owner_frag)
      deallocate(owner_frag, center_bohr, center_aa, spread_aa2, v_matrix, aa_global)
    end if

    if(allocated(esp_file)) deallocate(esp_file)
    call comm_sync_all(dc%icomm_tot)
    if(dc%id_tot == 0) write(*,'(1x,a)') "[DC-LCFO-W90-IMPORT] import-only completed without SCF."
  end subroutine dc_lcfo_wannier_import_only

  subroutine transform_wannier_position_gauge(num_bands, num_wann, v_from, v_to, aa_from, aa_to)
    implicit none
    integer, intent(in) :: num_bands, num_wann
    complex(8), intent(in) :: v_from(num_bands,num_wann), v_to(num_bands,num_wann)
    complex(8), intent(in) :: aa_from(3,num_wann,num_wann)
    complex(8), allocatable, intent(out) :: aa_to(:,:,:)
    integer :: idir
    complex(8), allocatable :: tmp_bw(:,:), band_pos(:,:), tmp_to(:,:)

    allocate(aa_to(3,num_wann,num_wann))
    allocate(tmp_bw(num_bands,num_wann), band_pos(num_bands,num_bands), tmp_to(num_bands,num_wann))
    aa_to = (0.0d0, 0.0d0)
    do idir=1,3
      tmp_bw(1:num_bands,1:num_wann) = matmul(v_from(1:num_bands,1:num_wann), &
        aa_from(idir,1:num_wann,1:num_wann))
      band_pos(1:num_bands,1:num_bands) = matmul(tmp_bw(1:num_bands,1:num_wann), &
        conjg(transpose(v_from(1:num_bands,1:num_wann))))
      tmp_to(1:num_bands,1:num_wann) = matmul(band_pos(1:num_bands,1:num_bands), &
        v_to(1:num_bands,1:num_wann))
      aa_to(idir,1:num_wann,1:num_wann) = matmul(conjg(transpose(v_to(1:num_bands,1:num_wann))), &
        tmp_to(1:num_bands,1:num_wann))
    end do
    call hermitize_wannier_position_matrix(num_wann, aa_to)
    deallocate(tmp_bw, band_pos, tmp_to)
  end subroutine transform_wannier_position_gauge

  subroutine write_wannier90_global_basis_stream_import(filename, num_bands, num_wann, n_frag_file, &
      owner_frag, center_bohr, spread_aa2, v_matrix, position_available, aa_global)
    use filesystem, only: get_filehandle
    implicit none
    character(*), intent(in) :: filename
    integer, intent(in) :: num_bands, num_wann, n_frag_file, position_available
    integer, intent(in) :: owner_frag(num_wann)
    real(8), intent(in) :: center_bohr(3,num_wann), spread_aa2(num_wann)
    complex(8), intent(in) :: v_matrix(num_bands,num_wann), aa_global(3,num_wann,num_wann)
    integer :: iunit

    iunit = get_filehandle()
    open(iunit,file=filename,form='unformatted',access='stream',status='replace')
    write(iunit) wannier90_global_magic, wannier90_global_version
    write(iunit) num_bands, num_wann, n_frag_file
    write(iunit) owner_frag(1:num_wann)
    write(iunit) center_bohr(1:3,1:num_wann)
    write(iunit) spread_aa2(1:num_wann)
    write(iunit) v_matrix(1:num_bands,1:num_wann)
    write(iunit) position_available
    write(iunit) aa_global(1:3,1:num_wann,1:num_wann)
    close(iunit)
  end subroutine write_wannier90_global_basis_stream_import

  subroutine write_wannier_flux_eigen_seed_stream_import(filename, num_bands, num_wann, nstate_seed, &
      nspin_file, n_frag_file, eval_seed, seed_wannier_to_eigen)
    use filesystem, only: get_filehandle
    implicit none
    character(*), intent(in) :: filename
    integer, intent(in) :: num_bands, num_wann, nstate_seed, nspin_file, n_frag_file
    real(8), intent(in) :: eval_seed(nstate_seed,nspin_file)
    complex(8), intent(in) :: seed_wannier_to_eigen(num_wann,nstate_seed)
    integer :: iunit

    iunit = get_filehandle()
    open(iunit,file=filename,form='unformatted',access='stream',status='replace')
    write(iunit) wannier_flux_eigen_seed_magic, wannier_flux_eigen_seed_version
    write(iunit) num_bands, num_wann, nstate_seed, nspin_file, n_frag_file
    write(iunit) eval_seed(1:nstate_seed,1:nspin_file)
    write(iunit) seed_wannier_to_eigen(1:num_wann,1:nstate_seed)
    close(iunit)
  end subroutine write_wannier_flux_eigen_seed_stream_import

  subroutine build_direct_wannier_buffer_position_from_lcfo_import(dc, num_bands, num_wann, &
      owner_frag, center_bohr, v_matrix, aa_direct, ok)
    use salmon_global, only: num_fragment, lambda_cut
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_bands, num_wann
    integer, intent(in) :: owner_frag(num_wann)
    real(8), intent(in) :: center_bohr(3,num_wann)
    complex(8), intent(in) :: v_matrix(num_bands,num_wann)
    complex(8), allocatable, intent(out) :: aa_direct(:,:,:)
    logical, intent(out) :: ok
    integer :: ifrag, nxyz_domain(3), nxyz_buffer(3), nxyz_box(3)
    integer :: nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag
    integer :: npts, ibx, iby, ibz, p, axis
    integer :: idx0(3)
    real(8), allocatable :: phi_basis(:,:), coef_wf(:,:), psi_state(:,:)
    real(8), allocatable :: coord(:), weight(:)
    complex(8), allocatable :: wannier_frag(:,:), weighted(:,:), overlap_direct(:,:)
    real(8) :: hvol, hgrid(3), coord_val, global_origin(3)
    logical :: seed_ok

    ok = .false.
    if(num_wann > 0 .and. size(owner_frag) < num_wann) return
    if(num_wann > 0 .and. size(center_bohr, 2) < num_wann) return
    hvol = dc%system_tot%hvol
    hgrid(1:3) = dc%system_tot%hgs(1:3)
    allocate(aa_direct(3,num_wann,num_wann), overlap_direct(num_wann,num_wann))
    aa_direct = (0.0d0, 0.0d0)
    overlap_direct = (0.0d0, 0.0d0)

    do ifrag=1,dc%n_frag
      call read_fragment_lcfo_buffer_seed_for_wannier_import(dc, ifrag, num_bands, &
        nxyz_domain, nxyz_buffer, nxyz_box, nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag, &
        phi_basis, coef_wf, seed_ok)
      if(.not. seed_ok) then
        if(allocated(phi_basis)) deallocate(phi_basis)
        if(allocated(coef_wf)) deallocate(coef_wf)
        if(allocated(aa_direct)) deallocate(aa_direct)
        return
      end if

      npts = product(nxyz_box)
      idx0(1:3) = dc%ixyz_frag(1:3,ifrag)
      do axis=1,3
        global_origin(axis) = dc%lg_tot%coordinate(idx0(axis),axis)
      end do
      allocate(psi_state(npts,num_bands), wannier_frag(npts,num_wann), weighted(npts,num_wann), &
        coord(npts), weight(npts))
      psi_state = matmul(phi_basis(1:npts,1:n_basis_frag), coef_wf(1:n_basis_frag,1:num_bands))
      wannier_frag = matmul(cmplx(psi_state(1:npts,1:num_bands), 0.0d0, kind=8), &
        v_matrix(1:num_bands,1:num_wann))

      p = 0
      do ibz=1,nxyz_box(3)
        do iby=1,nxyz_box(2)
          do ibx=1,nxyz_box(1)
            p = p + 1
            weight(p) = buffer_partition_weight_import(nxyz_domain, nxyz_buffer, num_fragment, ibx, iby, ibz)
          end do
        end do
      end do
      do axis=1,3
        p = 0
        do ibz=1,nxyz_box(3)
          do iby=1,nxyz_box(2)
            do ibx=1,nxyz_box(1)
              p = p + 1
              select case(axis)
              case(1)
                coord_val = global_origin(axis) + dble(ibx - nxyz_buffer(axis) - 1) * hgrid(axis)
              case(2)
                coord_val = global_origin(axis) + dble(iby - nxyz_buffer(axis) - 1) * hgrid(axis)
              case default
                coord_val = global_origin(axis) + dble(ibz - nxyz_buffer(axis) - 1) * hgrid(axis)
              end select
              coord(p) = coord_val
            end do
          end do
        end do
        do p=1,npts
          weighted(p,1:num_wann) = hvol * weight(p) * coord(p) * wannier_frag(p,1:num_wann)
        end do
        aa_direct(axis,1:num_wann,1:num_wann) = aa_direct(axis,1:num_wann,1:num_wann) + &
          matmul(conjg(transpose(wannier_frag(1:npts,1:num_wann))), weighted(1:npts,1:num_wann))
      end do
      do p=1,npts
        weighted(p,1:num_wann) = hvol * weight(p) * wannier_frag(p,1:num_wann)
      end do
      overlap_direct(1:num_wann,1:num_wann) = overlap_direct(1:num_wann,1:num_wann) + &
        matmul(conjg(transpose(wannier_frag(1:npts,1:num_wann))), weighted(1:npts,1:num_wann))

      deallocate(phi_basis, coef_wf, psi_state, wannier_frag, weighted, coord, weight)
    end do

    call normalize_wannier_position_by_overlap_import(num_wann, overlap_direct, aa_direct, lambda_cut, ok)
    if(.not. ok) then
      if(allocated(aa_direct)) deallocate(aa_direct)
      if(allocated(overlap_direct)) deallocate(overlap_direct)
      return
    end if
    deallocate(overlap_direct)
    call hermitize_wannier_position_matrix(num_wann, aa_direct)
    ok = .true.
    write(*,'(1x,a,i0,a,i0)') "[DC-LCFO-W90-IMPORT] direct AMN buffer-local position built: bands=", &
      num_bands, " wann=", num_wann
  end subroutine build_direct_wannier_buffer_position_from_lcfo_import

  real(8) function buffer_partition_weight_import(nxyz_domain, nxyz_buffer, nfrag_axis, ibx, iby, ibz) result(weight)
    implicit none
    integer, intent(in) :: nxyz_domain(3), nxyz_buffer(3), nfrag_axis(3)
    integer, intent(in) :: ibx, iby, ibz
    integer :: mult(3)

    mult(1) = buffer_axis_multiplicity_import(nxyz_domain(1), nxyz_buffer(1), nfrag_axis(1), ibx)
    mult(2) = buffer_axis_multiplicity_import(nxyz_domain(2), nxyz_buffer(2), nfrag_axis(2), iby)
    mult(3) = buffer_axis_multiplicity_import(nxyz_domain(3), nxyz_buffer(3), nfrag_axis(3), ibz)
    weight = 1.0d0 / dble(max(1, mult(1) * mult(2) * mult(3)))
  end function buffer_partition_weight_import

  integer function buffer_axis_multiplicity_import(ndom, nbuf, nfrag_axis, ibox) result(mult)
    implicit none
    integer, intent(in) :: ndom, nbuf, nfrag_axis, ibox
    integer :: jcore

    mult = 1
    if(ndom <= 0 .or. nbuf <= 0 .or. nfrag_axis <= 1) return
    jcore = modulo(ibox - nbuf - 1, ndom) + 1
    if(jcore <= nbuf) mult = mult + 1
    if(jcore > ndom - nbuf) mult = mult + 1
  end function buffer_axis_multiplicity_import

  subroutine normalize_wannier_position_by_overlap_import(num_wann, overlap_mat, aa_direct, cutoff, ok)
    use eigen_subdiag_sub, only: eigen_zheev
    implicit none
    integer, intent(in) :: num_wann
    complex(8), intent(in) :: overlap_mat(num_wann,num_wann)
    complex(8), intent(inout) :: aa_direct(3,num_wann,num_wann)
    real(8), intent(in) :: cutoff
    logical, intent(out) :: ok
    integer :: axis, i, j, k
    real(8), allocatable :: eval(:)
    complex(8), allocatable :: smat(:,:), eigvec(:,:), sinv(:,:), tmp(:,:)
    real(8) :: min_eval, max_eval, threshold

    ok = .false.
    if(num_wann <= 0) return
    allocate(smat(num_wann,num_wann), eigvec(num_wann,num_wann), sinv(num_wann,num_wann), &
      tmp(num_wann,num_wann), eval(num_wann))
    smat = overlap_mat
    call hermitize_square_matrix_import(num_wann, smat)
    call eigen_zheev(smat, eval, eigvec)
    min_eval = minval(eval(1:num_wann))
    max_eval = maxval(eval(1:num_wann))
    threshold = max(cutoff, 1.0d-12 * max(max_eval, 1.0d0))
    if(min_eval <= threshold) then
      write(*,'(1x,a,3(a,es12.5))') "[DC-LCFO-W90-IMPORT] buffer position overlap is rank deficient:", &
        " min_eval=", min_eval, " max_eval=", max_eval, " threshold=", threshold
      deallocate(smat, eigvec, sinv, tmp, eval)
      return
    end if
    sinv = (0.0d0, 0.0d0)
    do i=1,num_wann
      do j=1,num_wann
        do k=1,num_wann
          sinv(i,j) = sinv(i,j) + eigvec(i,k) * (1.0d0 / sqrt(eval(k))) * conjg(eigvec(j,k))
        end do
      end do
    end do
    do axis=1,3
      tmp(1:num_wann,1:num_wann) = matmul(aa_direct(axis,1:num_wann,1:num_wann), sinv)
      aa_direct(axis,1:num_wann,1:num_wann) = matmul(conjg(transpose(sinv)), tmp)
    end do
    write(*,'(1x,a,2(a,es12.5))') "[DC-LCFO-W90-IMPORT] buffer position overlap normalized:", &
      " min_eval=", min_eval, " max_eval=", max_eval
    deallocate(smat, eigvec, sinv, tmp, eval)
    ok = .true.
  end subroutine normalize_wannier_position_by_overlap_import

  subroutine hermitize_square_matrix_import(n, amat)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(inout) :: amat(n,n)
    integer :: i, j
    complex(8) :: zij, zji

    do i=1,n
      amat(i,i) = cmplx(real(amat(i,i), kind=8), 0.0d0, kind=8)
      do j=i+1,n
        zij = amat(i,j)
        zji = amat(j,i)
        amat(i,j) = 0.5d0 * (zij + conjg(zji))
        amat(j,i) = conjg(amat(i,j))
      end do
    end do
  end subroutine hermitize_square_matrix_import

  subroutine build_direct_wannier_phase_position_from_lcfo_import(dc, num_bands, num_wann, &
      v_matrix, aa_direct, ok)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_bands, num_wann
    complex(8), intent(in) :: v_matrix(num_bands,num_wann)
    complex(8), allocatable, intent(out) :: aa_direct(:,:,:)
    logical, intent(out) :: ok
    integer :: ifrag, nxyz_domain(3), nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag
    integer :: npts, ix, iy, iz, p, axis, idx0
    real(8), allocatable :: phi_basis(:,:), coef_wf(:,:), psi_state(:,:), phase_arg(:)
    complex(8), allocatable :: wannier_frag(:,:), weighted(:,:), u_phase(:,:,:), x_axis(:,:)
    real(8) :: hvol, cell_length(3), hgrid(3), total_length(3), fragment_origin(3), coord_val, gvec
    logical :: seed_ok, axis_ok
    complex(8), parameter :: zi = (0.0d0, 1.0d0)

    ok = .false.
    hvol = dc%system_tot%hvol
    allocate(aa_direct(3,num_wann,num_wann), u_phase(3,num_wann,num_wann))
    aa_direct = (0.0d0, 0.0d0)
    u_phase = (0.0d0, 0.0d0)
    do axis=1,3
      total_length(axis) = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
        + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
      if(total_length(axis) <= 0.0d0) return
    end do

    do ifrag=1,dc%n_frag
      call read_fragment_lcfo_seed_for_wannier_import(dc, ifrag, num_bands, &
        nxyz_domain, nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag, &
        phi_basis, coef_wf, seed_ok)
      if(.not. seed_ok) then
        if(allocated(phi_basis)) deallocate(phi_basis)
        if(allocated(coef_wf)) deallocate(coef_wf)
        if(allocated(aa_direct)) deallocate(aa_direct)
        if(allocated(u_phase)) deallocate(u_phase)
        return
      end if

      npts = product(nxyz_domain)
      call fragment_cell_lengths_import(dc, ifrag, cell_length)
      hgrid(1:3) = cell_length(1:3) / dble(nxyz_domain(1:3))
      do axis=1,3
        idx0 = dc%ixyz_frag(axis,ifrag)
        fragment_origin(axis) = dc%lg_tot%coordinate(idx0,axis)
      end do

      allocate(psi_state(npts,num_bands), wannier_frag(npts,num_wann), weighted(npts,num_wann), &
        phase_arg(npts))
      psi_state = matmul(phi_basis(1:npts,1:n_basis_frag), coef_wf(1:n_basis_frag,1:num_bands))
      wannier_frag = matmul(cmplx(psi_state(1:npts,1:num_bands), 0.0d0, kind=8), &
        v_matrix(1:num_bands,1:num_wann))

      do axis=1,3
        gvec = 2.0d0 * acos(-1.0d0) / total_length(axis)
        p = 0
        do iz=1,nxyz_domain(3)
          do iy=1,nxyz_domain(2)
            do ix=1,nxyz_domain(1)
              p = p + 1
              select case(axis)
              case(1)
                coord_val = fragment_origin(axis) + dble(ix - 1) * hgrid(axis)
              case(2)
                coord_val = fragment_origin(axis) + dble(iy - 1) * hgrid(axis)
              case default
                coord_val = fragment_origin(axis) + dble(iz - 1) * hgrid(axis)
              end select
              phase_arg(p) = gvec * coord_val
            end do
          end do
        end do
        do p=1,npts
          weighted(p,1:num_wann) = hvol * exp(zi * phase_arg(p)) * wannier_frag(p,1:num_wann)
        end do
        u_phase(axis,1:num_wann,1:num_wann) = u_phase(axis,1:num_wann,1:num_wann) + &
          matmul(conjg(transpose(wannier_frag(1:npts,1:num_wann))), weighted(1:npts,1:num_wann))
      end do

      deallocate(phi_basis, coef_wf, psi_state, wannier_frag, weighted, phase_arg)
    end do

    allocate(x_axis(num_wann,num_wann))
    do axis=1,3
      call phase_overlap_to_position_matrix(num_wann, total_length(axis), &
        u_phase(axis,1:num_wann,1:num_wann), x_axis, axis_ok)
      if(.not. axis_ok) then
        deallocate(x_axis)
        if(allocated(aa_direct)) deallocate(aa_direct)
        if(allocated(u_phase)) deallocate(u_phase)
        return
      end if
      aa_direct(axis,1:num_wann,1:num_wann) = x_axis(1:num_wann,1:num_wann)
    end do
    deallocate(x_axis, u_phase)
    call hermitize_wannier_position_matrix(num_wann, aa_direct)
    ok = .true.
    write(*,'(1x,a,i0,a,i0)') "[DC-LCFO-W90-IMPORT] direct AMN phase-log position built: bands=", &
      num_bands, " wann=", num_wann
  end subroutine build_direct_wannier_phase_position_from_lcfo_import

  subroutine phase_overlap_to_position_matrix(num_wann, length_axis, u_phase, x_mat, ok)
    use eigen_subdiag_sub, only: eigen_zheev
    implicit none
    integer, intent(in) :: num_wann
    real(8), intent(in) :: length_axis
    complex(8), intent(in) :: u_phase(num_wann,num_wann)
    complex(8), intent(out) :: x_mat(num_wann,num_wann)
    logical, intent(out) :: ok
    integer :: iw, jw, k, info, lwork
    real(8), allocatable :: eval(:), rwork(:), phase(:)
    complex(8), allocatable :: metric(:,:), eigvec(:,:), invsqrt(:,:), u_polar(:,:), work(:), &
      wr(:), vl(:,:), vr(:,:), tmp(:,:)
    complex(8) :: zij, zji, query(1)
    real(8) :: theta
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    external :: zgeev

    ok = .false.
    if(length_axis <= 0.0d0) return
    allocate(metric(num_wann,num_wann), eigvec(num_wann,num_wann), invsqrt(num_wann,num_wann), &
      u_polar(num_wann,num_wann), eval(num_wann))
    metric = matmul(conjg(transpose(u_phase)), u_phase)
    do iw=1,num_wann
      metric(iw,iw) = cmplx(real(metric(iw,iw), kind=8), 0.0d0, kind=8)
      do jw=iw+1,num_wann
        zij = metric(iw,jw)
        zji = metric(jw,iw)
        metric(iw,jw) = 0.5d0 * (zij + conjg(zji))
        metric(jw,iw) = conjg(metric(iw,jw))
      end do
    end do
    call eigen_zheev(metric, eval, eigvec)
    invsqrt = (0.0d0, 0.0d0)
    do jw=1,num_wann
      do iw=1,num_wann
        do k=1,num_wann
          if(eval(k) > 1.0d-12) invsqrt(iw,jw) = invsqrt(iw,jw) + &
            eigvec(iw,k) * (1.0d0 / sqrt(eval(k))) * conjg(eigvec(jw,k))
        end do
      end do
    end do
    u_polar = matmul(u_phase, invsqrt)

    allocate(wr(num_wann), vl(1,1), vr(num_wann,num_wann), rwork(2*num_wann), phase(num_wann))
    lwork = -1
    call zgeev('N','V',num_wann,u_polar,num_wann,wr,vl,1,vr,num_wann,query,lwork,rwork,info)
    if(info /= 0) then
      deallocate(metric, eigvec, invsqrt, u_polar, eval, wr, vl, vr, rwork, phase)
      return
    end if
    lwork = max(1, int(real(query(1), kind=8)))
    allocate(work(lwork))
    call zgeev('N','V',num_wann,u_polar,num_wann,wr,vl,1,vr,num_wann,work,lwork,rwork,info)
    if(info /= 0) then
      deallocate(metric, eigvec, invsqrt, u_polar, eval, wr, vl, vr, rwork, phase, work)
      return
    end if

    do k=1,num_wann
      theta = atan2(aimag(wr(k)), real(wr(k), kind=8))
      phase(k) = theta * length_axis / (2.0d0 * acos(-1.0d0))
    end do
    allocate(tmp(num_wann,num_wann))
    tmp = (0.0d0, 0.0d0)
    do k=1,num_wann
      tmp(k,1:num_wann) = cmplx(phase(k), 0.0d0, kind=8) * conjg(vr(1:num_wann,k))
    end do
    x_mat = matmul(vr, tmp)
    do iw=1,num_wann
      x_mat(iw,iw) = cmplx(real(x_mat(iw,iw), kind=8), 0.0d0, kind=8)
      do jw=iw+1,num_wann
        zij = x_mat(iw,jw)
        zji = x_mat(jw,iw)
        x_mat(iw,jw) = 0.5d0 * (zij + conjg(zji))
        x_mat(jw,iw) = conjg(x_mat(iw,jw))
      end do
    end do
    ok = .true.
    deallocate(metric, eigvec, invsqrt, u_polar, eval, wr, vl, vr, rwork, phase, work, tmp)
  end subroutine phase_overlap_to_position_matrix

  subroutine diagnose_wannier_transform_subspace_overlap(num_bands, num_wann, v_from, v_to, label)
    use eigen_subdiag_sub, only: eigen_zheev
    implicit none
    integer, intent(in) :: num_bands, num_wann
    complex(8), intent(in) :: v_from(num_bands,num_wann), v_to(num_bands,num_wann)
    character(*), intent(in) :: label
    integer :: iw, jw, ib
    complex(8) :: zij, zji
    complex(8), allocatable :: overlap(:,:), metric(:,:), eigvec(:,:)
    real(8), allocatable :: eval(:)
    real(8) :: sval_min, sval_max, trace_loss

    allocate(overlap(num_wann,num_wann), metric(num_wann,num_wann), eigvec(num_wann,num_wann), eval(num_wann))
    overlap = (0.0d0, 0.0d0)
    do jw=1,num_wann
      do iw=1,num_wann
        do ib=1,num_bands
          overlap(iw,jw) = overlap(iw,jw) + conjg(v_from(ib,iw)) * v_to(ib,jw)
        end do
      end do
    end do
    metric = matmul(conjg(transpose(overlap)), overlap)
    do iw=1,num_wann
      metric(iw,iw) = cmplx(real(metric(iw,iw), kind=8), 0.0d0, kind=8)
      do jw=iw+1,num_wann
        zij = metric(iw,jw)
        zji = metric(jw,iw)
        metric(iw,jw) = 0.5d0 * (zij + conjg(zji))
        metric(jw,iw) = conjg(metric(iw,jw))
      end do
    end do
    call eigen_zheev(metric, eval, eigvec)
    sval_min = sqrt(max(0.0d0, minval(eval(1:num_wann))))
    sval_max = sqrt(max(0.0d0, maxval(eval(1:num_wann))))
    trace_loss = 1.0d0 - sum(eval(1:num_wann)) / dble(max(1,num_wann))
    write(*,'(1x,a,a,3(a,1pe13.5))') "[DC-LCFO-W90-SYM] subspace overlap ", trim(label), &
      " smin=", sval_min, " smax=", sval_max, " trace_loss=", trace_loss
    deallocate(overlap, metric, eigvec, eval)
  end subroutine diagnose_wannier_transform_subspace_overlap

  subroutine build_wannier_flux_eigen_seed_from_transform(num_bands, num_wann, nstate_available, &
    nspin_in, v_matrix, esp_in, seed_wannier_to_eigen, eval_seed, nstate_seed)
    use eigen_subdiag_sub, only: eigen_zheev
    implicit none
    integer, intent(in) :: num_bands, num_wann, nstate_available, nspin_in
    complex(8), intent(in) :: v_matrix(num_bands,num_wann)
    real(8), intent(in) :: esp_in(nstate_available,nspin_in)
    complex(8), allocatable, intent(out) :: seed_wannier_to_eigen(:,:)
    real(8), allocatable, intent(out) :: eval_seed(:,:)
    integer, intent(out) :: nstate_seed
    integer :: ib, ispin, istate, iw, jw
    complex(8), allocatable :: h_wann(:,:), eigvec(:,:)
    real(8), allocatable :: eval(:)

    if(num_bands > nstate_available) then
      write(*,'(1x,a,2(a,i0))') "[DC-LCFO-W90-SEED] wavefunction seed has too few bands:", &
        " bands=", num_bands, " seed_states=", nstate_available
      stop "DC-LCFO Wannier flux seed: insufficient wavefunction seed"
    end if

    nstate_seed = num_wann
    allocate(seed_wannier_to_eigen(num_wann,nstate_seed))
    allocate(eval_seed(nstate_seed,nspin_in))
    seed_wannier_to_eigen = (0.0d0, 0.0d0)
    eval_seed = 0.0d0

    if(num_bands == num_wann) then
      do istate=1,nstate_seed
        do iw=1,num_wann
          seed_wannier_to_eigen(iw,istate) = conjg(v_matrix(istate,iw))
        end do
      end do
      eval_seed(1:nstate_seed,1:nspin_in) = esp_in(1:nstate_seed,1:nspin_in)
    else
      allocate(h_wann(num_wann,num_wann), eigvec(num_wann,num_wann), eval(num_wann))
      do ispin=1,nspin_in
        h_wann = (0.0d0, 0.0d0)
        do jw=1,num_wann
          do iw=1,num_wann
            do ib=1,num_bands
              h_wann(iw,jw) = h_wann(iw,jw) + conjg(v_matrix(ib,iw)) * esp_in(ib,ispin) * v_matrix(ib,jw)
            end do
          end do
        end do
        call eigen_zheev(h_wann, eval, eigvec)
        if(ispin == 1) seed_wannier_to_eigen(1:num_wann,1:nstate_seed) = eigvec(1:num_wann,1:num_wann)
        eval_seed(1:nstate_seed,ispin) = eval(1:nstate_seed)
      end do
      deallocate(h_wann, eigvec, eval)
      write(*,'(1x,a,3(a,i0))') "[DC-LCFO-W90-SEED] rectangular transform projected and diagonalized:", &
        " bands=", num_bands, " wann=", num_wann, " states=", nstate_seed
    end if
  end subroutine build_wannier_flux_eigen_seed_from_transform

  subroutine read_dc_lcfo_esp_from_wavefunctions_import(dc, nstate_tot_file, nspin_file, esp_file)
    use communication, only: comm_bcast
    use filesystem, only: get_filehandle
    use salmon_global, only: base_directory
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(out) :: nstate_tot_file, nspin_file
    real(8), allocatable, intent(out) :: esp_file(:,:)
    integer :: iunit, io, n_frag_file, nstate_frag_file
    integer, allocatable :: n_mat_tmp(:), n_basis_tmp(:,:), index_basis_tmp(:,:,:)
    real(8), allocatable :: coef_tmp(:,:,:)
    character(256) :: filename

    nstate_tot_file = 0
    nspin_file = 0
    if(dc%id_tot == 0) then
      filename = wannier_wavefunction_seed_filename_import(1)
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='old',iostat=io)
      if(io /= 0) then
        write(*,'(1x,2a)') "[DC-LCFO-W90-IMPORT] failed to open wavefunction seed: ", trim(filename)
        stop "DC-LCFO Wannier import-only: missing wavefunction seed"
      end if
      read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
      allocate(n_mat_tmp(nspin_file), n_basis_tmp(n_frag_file,nspin_file), &
        index_basis_tmp(nstate_frag_file,n_frag_file,nspin_file), &
        coef_tmp(nstate_frag_file,nstate_tot_file,nspin_file), &
        esp_file(nstate_tot_file,nspin_file))
      read(iunit) n_mat_tmp(1:nspin_file)
      read(iunit) n_basis_tmp(1:n_frag_file,1:nspin_file)
      read(iunit) index_basis_tmp(1:nstate_frag_file,1:n_frag_file,1:nspin_file)
      read(iunit) coef_tmp(1:nstate_frag_file,1:nstate_tot_file,1:nspin_file)
      read(iunit) esp_file(1:nstate_tot_file,1:nspin_file)
      close(iunit)
      write(*,'(1x,a,i0,a,i0)') "[DC-LCFO-W90-IMPORT] read eigenvalue seed: states=", &
        nstate_tot_file, " nspin=", nspin_file
      deallocate(n_mat_tmp, n_basis_tmp, index_basis_tmp, coef_tmp)
    end if
    call comm_bcast(nstate_tot_file, dc%icomm_tot)
    call comm_bcast(nspin_file, dc%icomm_tot)
    if(dc%id_tot /= 0) allocate(esp_file(max(1,nstate_tot_file),max(1,nspin_file)))
    if(nstate_tot_file > 0 .and. nspin_file > 0) &
      call comm_bcast(esp_file, dc%icomm_tot)
  end subroutine read_dc_lcfo_esp_from_wavefunctions_import

  character(256) function wannier_wavefunction_seed_filename_import(ifrag) result(filename)
    implicit none
    integer, intent(in) :: ifrag
    character(256) :: candidate
    logical :: exists

    write(candidate, '(a,a,i6.6,a,a)') trim(import_run_root_dir()), 'data_dcdft/fragments/', ifrag, '/', &
      binfile_wf_wannier_seed
    inquire(file=candidate, exist=exists)
    if(exists) then
      filename = candidate
    else
      write(filename, '(a,a,i6.6,a,a)') trim(import_run_root_dir()), 'data_dcdft/fragments/', ifrag, '/', binfile_wf
    end if
  end function wannier_wavefunction_seed_filename_import

  subroutine read_wannier90_checkpoint_transform_import(dc, num_bands_chk, num_wann_chk, &
      v_matrix, center_aa, spread_aa2)
    use filesystem, only: get_filehandle
    use salmon_global, only: base_directory, sysname
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(out) :: num_bands_chk, num_wann_chk
    complex(8), allocatable, intent(out) :: v_matrix(:,:)
    real(8), allocatable, intent(out) :: center_aa(:,:), spread_aa2(:)
    integer :: iunit, io, i, j, k, nkp
    integer :: num_exclude_bands_chk, num_kpts_chk, nntot_chk
    integer :: mp_grid_chk(3)
    integer, allocatable :: exclude_bands_chk(:), ndimwin(:)
    real(8) :: real_lattice_chk(3,3), recip_lattice_chk(3,3), omega_invariant_chk
    real(8), allocatable :: kpt_latt_chk(:,:)
    logical :: have_disentangled_chk
    logical, allocatable :: lwindow(:,:)
    complex(8), allocatable :: u_matrix(:,:,:), u_matrix_opt(:,:,:), m_matrix(:,:,:,:)
    character(256) :: filename
    character(33) :: header_chk
    character(20) :: checkpoint_chk

    filename = trim(import_run_root_dir())//'data_dcdft/total/'//trim(sysname)//".chk"
    iunit = get_filehandle()
    open(iunit,file=filename,status='old',form='unformatted',iostat=io)
    if(io /= 0) then
      write(*,'(1x,2a)') "[DC-LCFO-W90-IMPORT] failed to open checkpoint: ", trim(filename)
      stop "DC-LCFO Wannier import-only: missing Wannier90 checkpoint"
    end if
    read(iunit) header_chk
    read(iunit) num_bands_chk
    read(iunit) num_exclude_bands_chk
    allocate(exclude_bands_chk(max(0,num_exclude_bands_chk)))
    if(num_exclude_bands_chk > 0) then
      read(iunit) (exclude_bands_chk(i), i=1,num_exclude_bands_chk)
    else
      read(iunit)
    end if
    read(iunit) ((real_lattice_chk(i,j), i=1,3), j=1,3)
    read(iunit) ((recip_lattice_chk(i,j), i=1,3), j=1,3)
    read(iunit) num_kpts_chk
    read(iunit) (mp_grid_chk(i), i=1,3)
    allocate(kpt_latt_chk(3,num_kpts_chk))
    read(iunit) ((kpt_latt_chk(i,nkp), i=1,3), nkp=1,num_kpts_chk)
    read(iunit) nntot_chk
    read(iunit) num_wann_chk
    read(iunit) checkpoint_chk
    read(iunit) have_disentangled_chk

    if(num_kpts_chk /= 1) stop "DC-LCFO Wannier import-only: only Gamma checkpoint is supported."
    if(have_disentangled_chk) then
      read(iunit) omega_invariant_chk
      allocate(lwindow(num_bands_chk,num_kpts_chk), ndimwin(num_kpts_chk))
      read(iunit) ((lwindow(i,nkp), i=1,num_bands_chk), nkp=1,num_kpts_chk)
      read(iunit) (ndimwin(nkp), nkp=1,num_kpts_chk)
      allocate(u_matrix_opt(num_bands_chk,num_wann_chk,num_kpts_chk))
      read(iunit) (((u_matrix_opt(i,j,nkp), i=1,num_bands_chk), j=1,num_wann_chk), nkp=1,num_kpts_chk)
    end if

    allocate(u_matrix(num_wann_chk,num_wann_chk,num_kpts_chk))
    read(iunit) (((u_matrix(i,j,k), i=1,num_wann_chk), j=1,num_wann_chk), k=1,num_kpts_chk)
    allocate(m_matrix(num_wann_chk,num_wann_chk,nntot_chk,num_kpts_chk))
    read(iunit) ((((m_matrix(i,j,k,nkp), i=1,num_wann_chk), j=1,num_wann_chk), k=1,nntot_chk), nkp=1,num_kpts_chk)
    allocate(center_aa(3,num_wann_chk), spread_aa2(num_wann_chk))
    read(iunit) ((center_aa(i,j), i=1,3), j=1,num_wann_chk)
    read(iunit) (spread_aa2(i), i=1,num_wann_chk)
    close(iunit)

    allocate(v_matrix(num_bands_chk,num_wann_chk))
    v_matrix = (0d0,0d0)
    if(have_disentangled_chk) then
      v_matrix(1:num_bands_chk,1:num_wann_chk) = matmul(u_matrix_opt(:,:,1), u_matrix(:,:,1))
    else
      if(num_bands_chk < num_wann_chk) &
        stop "DC-LCFO Wannier import-only: invalid checkpoint without disentanglement."
      v_matrix(1:num_wann_chk,1:num_wann_chk) = u_matrix(:,:,1)
    end if
    write(*,'(1x,a,i0,a,i0,a,l1)') "[DC-LCFO-W90-IMPORT] read checkpoint: bands=", &
      num_bands_chk, " wann=", num_wann_chk, " disentangled=", have_disentangled_chk

    deallocate(exclude_bands_chk, kpt_latt_chk, u_matrix, m_matrix)
    if(allocated(lwindow)) deallocate(lwindow)
    if(allocated(ndimwin)) deallocate(ndimwin)
    if(allocated(u_matrix_opt)) deallocate(u_matrix_opt)
  end subroutine read_wannier90_checkpoint_transform_import

  subroutine read_wannier90_amn_direct_transform_import(dc, num_bands_expected, num_wann_expected, owner_frag, &
      v_direct, ok, block_by_owner)
    use eigen_subdiag_sub, only: eigen_zheev
    use filesystem, only: get_filehandle
    use salmon_global, only: sysname, lambda_cut
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_bands_expected, num_wann_expected
    integer, intent(in) :: owner_frag(num_wann_expected)
    complex(8), allocatable, intent(out) :: v_direct(:,:)
    logical, intent(out) :: ok
    logical, intent(in), optional :: block_by_owner
    character(256) :: filename, header
    integer :: iunit, io, num_bands_file, num_kpts_file, num_wann_file
    integer :: irec, iband, iwann, ikpt, i, j, k, ifrag, nowned, iloc, nreject
    real(8) :: re_val, im_val, min_eval, max_eval, min_eval_all, max_eval_all
    integer, allocatable :: widx(:)
    real(8), allocatable :: eval(:)
    complex(8), allocatable :: amat(:,:), ablock(:,:), gram(:,:), eigvec(:,:), sinv(:,:)
    logical :: use_block_by_owner

    ok = .false.
    use_block_by_owner = .true.
    if(present(block_by_owner)) use_block_by_owner = block_by_owner
    if(allocated(v_direct)) deallocate(v_direct)
    if(num_bands_expected <= 0 .or. num_wann_expected <= 0) return

    if(dc%id_tot == 0) then
      filename = trim(import_run_root_dir())//'data_dcdft/total/'//trim(sysname)//".amn"
      iunit = get_filehandle()
      open(iunit,file=filename,status='old',action='read',iostat=io)
      if(io /= 0) then
        write(*,'(1x,2a)') "[DC-LCFO-W90-IMPORT] direct AMN diagnostic skipped; missing file: ", trim(filename)
      else
        read(iunit,'(a)',iostat=io) header
        if(io == 0) read(iunit,*,iostat=io) num_bands_file, num_kpts_file, num_wann_file
        if(io /= 0 .or. num_bands_file /= num_bands_expected .or. num_wann_file /= num_wann_expected .or. &
            num_kpts_file /= 1) then
          write(*,'(1x,a,6(a,i0))') "[DC-LCFO-W90-IMPORT] direct AMN diagnostic skipped; dimensions:", &
            " bands=", num_bands_file, " expected_bands=", num_bands_expected, &
            " wann=", num_wann_file, " expected_wann=", num_wann_expected, &
            " kpts=", num_kpts_file, " expected_kpts=", 1
        else
          allocate(amat(num_bands_expected,num_wann_expected))
          amat = (0.0d0,0.0d0)
          do irec=1,num_bands_expected*num_wann_expected
            read(iunit,*,iostat=io) iband, iwann, ikpt, re_val, im_val
            if(io /= 0) exit
            if(iband < 1 .or. iband > num_bands_expected) cycle
            if(iwann < 1 .or. iwann > num_wann_expected) cycle
            if(ikpt /= 1) cycle
            amat(iband,iwann) = cmplx(re_val, im_val, kind=8)
          end do
          if(io == 0) then
            allocate(v_direct(num_bands_expected,num_wann_expected))
            v_direct = (0.0d0,0.0d0)
            min_eval_all = huge(1.0d0)
            max_eval_all = -huge(1.0d0)
            nreject = 0
            if(use_block_by_owner) then
              do ifrag=1,dc%n_frag
                nowned = count(owner_frag(1:num_wann_expected) == ifrag)
                if(nowned <= 0) cycle
                allocate(widx(nowned), ablock(num_bands_expected,nowned), gram(nowned,nowned), &
                  eigvec(nowned,nowned), sinv(nowned,nowned), eval(nowned))
                iloc = 0
                do iwann=1,num_wann_expected
                  if(owner_frag(iwann) /= ifrag) cycle
                  iloc = iloc + 1
                  widx(iloc) = iwann
                  ablock(1:num_bands_expected,iloc) = amat(1:num_bands_expected,iwann)
                end do
                gram = matmul(conjg(transpose(ablock)), ablock)
                call eigen_zheev(gram, eval, eigvec)
                min_eval = minval(eval)
                max_eval = maxval(eval)
                min_eval_all = min(min_eval_all, min_eval)
                max_eval_all = max(max_eval_all, max_eval)
                nreject = nreject + count(eval(1:nowned) <= lambda_cut)
                sinv = (0.0d0,0.0d0)
                do i=1,nowned
                  do j=1,nowned
                    do k=1,nowned
                      if(eval(k) > lambda_cut) then
                        sinv(i,j) = sinv(i,j) + eigvec(i,k) * (1.0d0 / sqrt(eval(k))) * conjg(eigvec(j,k))
                      end if
                    end do
                  end do
                end do
                ablock = matmul(ablock, sinv)
                do iloc=1,nowned
                  v_direct(1:num_bands_expected,widx(iloc)) = ablock(1:num_bands_expected,iloc)
                end do
                deallocate(widx, ablock, gram, eigvec, sinv, eval)
              end do
            else
              nowned = num_wann_expected
              allocate(ablock(num_bands_expected,nowned), gram(nowned,nowned), &
                eigvec(nowned,nowned), sinv(nowned,nowned), eval(nowned))
              ablock(1:num_bands_expected,1:nowned) = amat(1:num_bands_expected,1:nowned)
              gram = matmul(conjg(transpose(ablock)), ablock)
              call eigen_zheev(gram, eval, eigvec)
              min_eval = minval(eval)
              max_eval = maxval(eval)
              min_eval_all = min(min_eval_all, min_eval)
              max_eval_all = max(max_eval_all, max_eval)
              nreject = nreject + count(eval(1:nowned) <= lambda_cut)
              sinv = (0.0d0,0.0d0)
              do i=1,nowned
                do j=1,nowned
                  do k=1,nowned
                    if(eval(k) > lambda_cut) then
                      sinv(i,j) = sinv(i,j) + eigvec(i,k) * (1.0d0 / sqrt(eval(k))) * conjg(eigvec(j,k))
                    end if
                  end do
                end do
              end do
              ablock = matmul(ablock, sinv)
              v_direct(1:num_bands_expected,1:nowned) = ablock(1:num_bands_expected,1:nowned)
              deallocate(ablock, gram, eigvec, sinv, eval)
            end if
            if(nreject > 0) then
              write(*,'(1x,a,i0,3(a,es12.5))') &
                "[DC-LCFO-W90-IMPORT] direct AMN projection matrix is rank deficient:", &
                nreject, " min_proj_s=", min_eval_all, " max_proj_s=", max_eval_all, &
                " lambda_cut=", lambda_cut
              stop "DC-LCFO Wannier import: direct AMN projection matrix is rank deficient"
            end if
            ok = .true.
            if(use_block_by_owner) then
              write(*,'(1x,a,2(a,es12.5))') "[DC-LCFO-W90-IMPORT] direct AMN block-orthonormalized:", &
                " min_proj_s=", min_eval_all, " max_proj_s=", max_eval_all
            else
              write(*,'(1x,a,2(a,es12.5))') "[DC-LCFO-W90-IMPORT] direct AMN global-orthonormalized:", &
                " min_proj_s=", min_eval_all, " max_proj_s=", max_eval_all
            end if
          end if
          deallocate(amat)
        end if
        close(iunit)
      end if
    end if
  end subroutine read_wannier90_amn_direct_transform_import

  subroutine hermitize_wannier_position_matrix(num_wann, aa_global)
    implicit none
    integer, intent(in) :: num_wann
    complex(8), intent(inout) :: aa_global(:,:,:)
    integer :: axis, iw, jw
    complex(8) :: zij, zji

    if(num_wann <= 0) return
    if(size(aa_global, 1) < 3 .or. size(aa_global, 2) < num_wann .or. &
       size(aa_global, 3) < num_wann) return
    do axis = 1, 3
      do iw = 1, num_wann
        aa_global(axis, iw, iw) = cmplx(real(aa_global(axis, iw, iw), kind=8), 0d0, kind=8)
        do jw = iw + 1, num_wann
          zij = aa_global(axis, iw, jw)
          zji = aa_global(axis, jw, iw)
          aa_global(axis, iw, jw) = 0.5d0 * (zij + conjg(zji))
          aa_global(axis, jw, iw) = conjg(aa_global(axis, iw, jw))
        end do
      end do
    end do
  end subroutine hermitize_wannier_position_matrix

  subroutine set_wannier_centers_from_position_diagonal(num_wann, center_bohr, aa_global)
    implicit none
    integer, intent(in) :: num_wann
    real(8), intent(inout) :: center_bohr(3,num_wann)
    complex(8), intent(in) :: aa_global(:,:,:)
    integer :: axis, iw

    if(num_wann <= 0) return
    if(size(aa_global, 1) < 3 .or. size(aa_global, 2) < num_wann .or. &
       size(aa_global, 3) < num_wann) return
    do iw=1,num_wann
      do axis=1,3
        center_bohr(axis,iw) = real(aa_global(axis,iw,iw), kind=8)
      end do
    end do
  end subroutine set_wannier_centers_from_position_diagonal

  subroutine wrap_wannier_centers_to_total_cell_import(dc, center_bohr, num_wann)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_wann
    real(8), intent(inout) :: center_bohr(3,num_wann)
    integer :: iw, axis
    real(8) :: total_length(3)

    call total_cell_lengths_import(dc, total_length)
    do iw=1,num_wann
      do axis=1,3
        if(total_length(axis) <= 0.0d0) cycle
        center_bohr(axis,iw) = center_bohr(axis,iw) - floor(center_bohr(axis,iw) / total_length(axis)) &
          * total_length(axis)
      end do
    end do
  end subroutine wrap_wannier_centers_to_total_cell_import

  subroutine read_wannier90_global_rmn_gamma_block_import(dc, num_wann_expected, aa_global, ok)
    use filesystem, only: get_filehandle
    use inputoutput, only: au_length_aa
    use salmon_global, only: base_directory, sysname
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_wann_expected
    complex(8), allocatable, intent(out) :: aa_global(:,:,:)
    logical, intent(out) :: ok
    integer :: iunit, io, num_wann_file, nrpts_file, ir
    integer :: rvec(3), n, m
    real(8) :: rx_re, rx_im, ry_re, ry_im, rz_re, rz_im
    character(256) :: filename, header
    logical :: exists

    ok = .false.
    filename = trim(import_run_root_dir())//'data_dcdft/total/'//trim(sysname)//"_r.dat"
    inquire(file=filename, exist=exists)
    if(.not. exists) return
    iunit = get_filehandle()
    open(iunit,file=filename,status='old',action='read')
    read(iunit,'(a)',iostat=io) header
    if(io == 0) read(iunit,*,iostat=io) num_wann_file
    if(io == 0) read(iunit,*,iostat=io) nrpts_file
    if(io /= 0 .or. num_wann_file /= num_wann_expected) then
      close(iunit)
      return
    end if
    allocate(aa_global(3,num_wann_file,num_wann_file))
    aa_global = (0d0,0d0)
    do ir=1,nrpts_file*num_wann_file*num_wann_file
      read(iunit,*,iostat=io) rvec(1:3), n, m, rx_re, rx_im, ry_re, ry_im, rz_re, rz_im
      if(io /= 0) exit
      if(any(rvec(1:3) /= 0)) cycle
      if(n < 1 .or. n > num_wann_file .or. m < 1 .or. m > num_wann_file) cycle
      aa_global(1,n,m) = cmplx(rx_re, rx_im, kind=8) / au_length_aa
      aa_global(2,n,m) = cmplx(ry_re, ry_im, kind=8) / au_length_aa
      aa_global(3,n,m) = cmplx(rz_re, rz_im, kind=8) / au_length_aa
    end do
    close(iunit)
    call hermitize_wannier_position_matrix(num_wann_file, aa_global)
    ok = .true.
  end subroutine read_wannier90_global_rmn_gamma_block_import

  logical function is_bond_center_projection_import(text) result(enabled)
    implicit none
    character(*), intent(in) :: text
    character(256) :: work

    work = adjustl(text)
    enabled = (trim(work) == 'bond_centers' .or. trim(work) == 'BOND_CENTERS')
  end function is_bond_center_projection_import

  subroutine build_bond_center_projection_map_import(dc, nproj, center_bohr)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: nproj
    real(8), allocatable, intent(out) :: center_bohr(:,:)
    integer :: ia, ja, axis, ip, ibond, nbond
    real(8) :: cutoff, min_dist, dist2, delta(3), center(3), length_axis(3)

    if(nproj <= 0) stop "DC-LCFO Wannier import: bond_centers requires positive num_wann."
    call find_bond_center_cutoff_import(dc, min_dist, cutoff)
    nbond = count_bond_centers_with_cutoff_import(dc, cutoff)
    if(nbond <= 0) stop "DC-LCFO Wannier import: no bond-center projection candidates were generated."
    if(nproj > nbond) stop "DC-LCFO Wannier import: target_wann exceeds unique bond-center candidates"
    call total_cell_lengths_import(dc, length_axis)

    allocate(center_bohr(3,nproj))
    ip = 0
    do while(ip < nproj)
      ibond = 0
      do ia=1,dc%system_tot%nion-1
        do ja=ia+1,dc%system_tot%nion
          dist2 = 0d0
          do axis=1,3
            delta(axis) = periodic_delta_import(dc%system_tot%rion(axis,ja) - dc%system_tot%rion(axis,ia), &
              length_axis(axis))
            dist2 = dist2 + delta(axis) * delta(axis)
          end do
          if(sqrt(dist2) > cutoff) cycle
          ibond = ibond + 1
          ip = ip + 1
          do axis=1,3
            center(axis) = dc%system_tot%rion(axis,ia) + 0.5d0 * delta(axis)
            if(length_axis(axis) > 0d0) center(axis) = center(axis) - floor(center(axis) / length_axis(axis)) &
              * length_axis(axis)
            center_bohr(axis,ip) = center(axis)
          end do
          if(ip >= nproj) exit
        end do
        if(ip >= nproj) exit
      end do
      if(ibond <= 0) exit
    end do
    if(ip < nproj) stop "DC-LCFO Wannier import: failed to complete bond-center projection map."
    if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0,a,es12.5,a,es12.5)') &
      "[DC-LCFO-WANNIER] bond-center candidates=", nbond, " target_wann=", nproj, &
      " nearest=", min_dist, " cutoff=", cutoff
  end subroutine build_bond_center_projection_map_import

  subroutine find_bond_center_cutoff_import(dc, min_dist, cutoff)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    real(8), intent(out) :: min_dist, cutoff
    integer :: ia, ja, axis
    real(8) :: dist2, dist, delta_axis, length_axis(3)

    call total_cell_lengths_import(dc, length_axis)
    min_dist = huge(1d0)
    do ia=1,dc%system_tot%nion-1
      do ja=ia+1,dc%system_tot%nion
        dist2 = 0d0
        do axis=1,3
          delta_axis = periodic_delta_import(dc%system_tot%rion(axis,ja) - dc%system_tot%rion(axis,ia), &
            length_axis(axis))
          dist2 = dist2 + delta_axis * delta_axis
        end do
        dist = sqrt(dist2)
        if(dist > 1d-8) min_dist = min(min_dist, dist)
      end do
    end do
    if(min_dist >= huge(1d0) * 0.5d0) &
      stop "DC-LCFO Wannier import: failed to find nearest-neighbor bond distance."
    cutoff = 1.20d0 * min_dist
  end subroutine find_bond_center_cutoff_import

  integer function count_bond_centers_with_cutoff_import(dc, cutoff) result(nbond)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    real(8), intent(in) :: cutoff
    integer :: ia, ja, axis
    real(8) :: dist2, delta_axis, length_axis(3)

    call total_cell_lengths_import(dc, length_axis)
    nbond = 0
    do ia=1,dc%system_tot%nion-1
      do ja=ia+1,dc%system_tot%nion
        dist2 = 0d0
        do axis=1,3
          delta_axis = periodic_delta_import(dc%system_tot%rion(axis,ja) - dc%system_tot%rion(axis,ia), &
            length_axis(axis))
          dist2 = dist2 + delta_axis * delta_axis
        end do
        if(sqrt(dist2) <= cutoff) nbond = nbond + 1
      end do
    end do
  end function count_bond_centers_with_cutoff_import

  integer function find_owner_fragment_from_center_import(dc, center_bohr) result(owner)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    real(8), intent(in) :: center_bohr(3)
    integer :: ifrag_try, axis, idx0, idx1, nxyz_domain(3)
    real(8) :: frag_center(3), dist2, best_dist2, delta_axis, length_axis

    owner = 1
    best_dist2 = huge(1d0)
    do ifrag_try=1,dc%n_frag
      call get_fragment_domain(dc, ifrag_try, nxyz_domain)
      do axis=1,3
        idx0 = dc%ixyz_frag(axis,ifrag_try)
        idx1 = idx0 + nxyz_domain(axis) - 1
        frag_center(axis) = 0.5d0 * (dc%lg_tot%coordinate(idx0,axis) + dc%lg_tot%coordinate(idx1,axis))
      end do
      dist2 = 0d0
      do axis=1,3
        length_axis = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
          + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
        delta_axis = periodic_delta_import(center_bohr(axis) - frag_center(axis), length_axis)
        dist2 = dist2 + delta_axis * delta_axis
      end do
      if(dist2 < best_dist2) then
        best_dist2 = dist2
        owner = ifrag_try
      end if
    end do
  end function find_owner_fragment_from_center_import

  subroutine rebalance_wannier_owner_fragments_import(dc, center_bohr, owner_frag, num_wann)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_wann
    real(8), intent(in) :: center_bohr(3,num_wann)
    integer, intent(inout) :: owner_frag(num_wann)
    integer :: target_base, target_rem, ifrag, iw, donor, receiver, moved, iw_best
    integer, allocatable :: target(:), count_frag(:)
    real(8) :: dist2, best_dist2

    if(dc%n_frag <= 0) return
    allocate(target(dc%n_frag), count_frag(dc%n_frag))
    target_base = num_wann / dc%n_frag
    target_rem = mod(num_wann, dc%n_frag)
    do ifrag=1,dc%n_frag
      target(ifrag) = target_base
      if(ifrag <= target_rem) target(ifrag) = target(ifrag) + 1
      count_frag(ifrag) = count(owner_frag(1:num_wann) == ifrag)
    end do
    moved = 0
    do
      receiver = 0
      donor = 0
      do ifrag=1,dc%n_frag
        if(count_frag(ifrag) < target(ifrag)) then
          receiver = ifrag
          exit
        end if
      end do
      if(receiver == 0) exit
      do ifrag=1,dc%n_frag
        if(count_frag(ifrag) > target(ifrag)) then
          donor = ifrag
          exit
        end if
      end do
      if(donor == 0) exit
      iw_best = 0
      best_dist2 = huge(1d0)
      do iw=1,num_wann
        if(owner_frag(iw) /= donor) cycle
        dist2 = distance_to_fragment_center_import(dc, center_bohr(1:3,iw), receiver)
        if(dist2 < best_dist2) then
          best_dist2 = dist2
          iw_best = iw
        end if
      end do
      if(iw_best == 0) exit
      owner_frag(iw_best) = receiver
      count_frag(donor) = count_frag(donor) - 1
      count_frag(receiver) = count_frag(receiver) + 1
      moved = moved + 1
    end do
    if(dc%id_tot == 0 .and. moved > 0) write(*,'(1x,a,i0,a,i0,a,i0)') &
      "[DC-LCFO-W90-IMPORT] rebalanced Wannier owners: moved=", moved, &
      " min_count=", minval(count_frag), " max_count=", maxval(count_frag)
    deallocate(target, count_frag)
  end subroutine rebalance_wannier_owner_fragments_import

  real(8) function distance_to_fragment_center_import(dc, center_bohr, ifrag_target) result(dist2)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    real(8), intent(in) :: center_bohr(3)
    integer, intent(in) :: ifrag_target
    integer :: axis, idx0, idx1, nxyz_domain(3)
    real(8) :: frag_center(3), delta_axis, length_axis

    call get_fragment_domain(dc, ifrag_target, nxyz_domain)
    do axis=1,3
      idx0 = dc%ixyz_frag(axis,ifrag_target)
      idx1 = idx0 + nxyz_domain(axis) - 1
      frag_center(axis) = 0.5d0 * (dc%lg_tot%coordinate(idx0,axis) + dc%lg_tot%coordinate(idx1,axis))
    end do
    dist2 = 0d0
    do axis=1,3
      length_axis = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
        + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
      delta_axis = periodic_delta_import(center_bohr(axis) - frag_center(axis), length_axis)
      dist2 = dist2 + delta_axis * delta_axis
    end do
  end function distance_to_fragment_center_import

  real(8) function periodic_delta_import(delta, length) result(dout)
    implicit none
    real(8), intent(in) :: delta, length
    dout = delta
    if(length <= 0d0) return
    dout = delta - dnint(delta / length) * length
  end function periodic_delta_import

  subroutine detect_wannier_fragment_symops(dc, nsym, symops)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(out) :: nsym
    type(t_wannier_symop), allocatable, intent(out) :: symops(:)
    integer :: ifrag, ia, ja, natom_frag, natom_buffer
    integer, allocatable :: atom_index(:), atom_index_buffer(:)
    real(8), allocatable :: atom_pos_local(:,:), atom_pos_buffer(:,:)
    real(8) :: cell_length(3), buffer_length(3), origin(3), residual, best_residual, best_origin(3)
    real(8) :: buffer_origin(3)
    logical :: accepted, found_inversion

    nsym = 0
    if(dc%n_frag <= 0) then
      allocate(symops(0))
      return
    end if

    allocate(symops(2 * dc%n_frag))

    do ifrag=1,dc%n_frag
      call collect_fragment_core_atoms_import(dc, ifrag, atom_index, natom_frag)
      call fragment_cell_lengths_import(dc, ifrag, cell_length)
      allocate(atom_pos_local(3,max(1,natom_frag)))
      do ia=1,natom_frag
        call atom_fragment_local_position_import(dc, ifrag, atom_index(ia), atom_pos_local(1:3,ia))
      end do

      nsym = nsym + 1
      symops(nsym)%owner_frag = ifrag
      symops(nsym)%rot = 0
      symops(nsym)%rot(1,1) = 1
      symops(nsym)%rot(2,2) = 1
      symops(nsym)%rot(3,3) = 1
      call fragment_center_bohr_import(dc, ifrag, symops(nsym)%origin_bohr)
      symops(nsym)%tau_local = 0.0d0
      symops(nsym)%atom_residual = 0.0d0
      symops(nsym)%label = 'identity'
      write(*,'(1x,a,i0,2a,es12.5,a,i0)') &
        "[DC-LCFO-W90-SYM] fragment=", ifrag, " detected symop label=", &
        trim(symops(nsym)%label), symops(nsym)%atom_residual, " natom=", natom_frag

      found_inversion = .false.
      best_residual = huge(1d0)
      best_origin = 0.0d0
      do ia=1,natom_frag
        do ja=ia,natom_frag
          call periodic_pair_midpoint_import(atom_pos_local(1:3,ia), atom_pos_local(1:3,ja), &
            cell_length, origin)
          call test_fragment_inversion_import(dc, atom_index, atom_pos_local, natom_frag, &
            cell_length, origin, accepted, residual)
          if(accepted .and. residual < best_residual) then
            found_inversion = .true.
            best_residual = residual
            best_origin = origin
          end if
        end do
      end do

      if(found_inversion) then
        nsym = nsym + 1
        symops(nsym)%owner_frag = ifrag
        symops(nsym)%rot = 0
        symops(nsym)%rot(1,1) = -1
        symops(nsym)%rot(2,2) = -1
        symops(nsym)%rot(3,3) = -1
        symops(nsym)%origin_bohr = best_origin
        symops(nsym)%tau_local = 0.0d0
        symops(nsym)%atom_residual = best_residual
        symops(nsym)%label = 'inversion'
        write(*,'(1x,a,i0,2a,es12.5,a,3(1x,es12.5))') &
          "[DC-LCFO-W90-SYM] fragment=", ifrag, " detected symop label=", &
          trim(symops(nsym)%label), symops(nsym)%atom_residual, " origin=", best_origin(1:3)
      end if

      call collect_fragment_buffer_atoms_import(dc, ifrag, atom_index_buffer, atom_pos_buffer, natom_buffer)
      call fragment_buffer_lengths_import(dc, ifrag, buffer_length)
      if(found_inversion) then
        call fragment_buffer_shift_import(dc, buffer_origin)
        buffer_origin(1:3) = best_origin(1:3) + buffer_origin(1:3)
        call test_fragment_inversion_nonperiodic_import(dc, atom_index_buffer, atom_pos_buffer, natom_buffer, &
          buffer_origin, accepted, residual)
        write(*,'(1x,a,i0,a,es12.5,a,l1,a,i0,a,3(1x,es12.5))') &
          "[DC-LCFO-W90-SYM] fragment=", ifrag, " buffer finite_box from_core label=inversion", &
          residual, " closed=", accepted, " natom_buffer=", natom_buffer, " origin=", buffer_origin(1:3)
        call test_fragment_inversion_global_pbc_import(dc, ifrag, atom_pos_buffer, natom_buffer, &
          buffer_origin, accepted, residual)
        write(*,'(1x,a,i0,a,es12.5,a,l1,a,i0,a,3(1x,es12.5))') &
          "[DC-LCFO-W90-SYM] fragment=", ifrag, " buffer global_pbc from_core label=inversion", &
          residual, " closed=", accepted, " natom_buffer=", natom_buffer, " origin=", buffer_origin(1:3)
      end if
      found_inversion = .false.
      best_residual = huge(1d0)
      best_origin = 0.0d0
      do ia=1,natom_buffer
        do ja=ia,natom_buffer
          origin(1:3) = 0.5d0 * (atom_pos_buffer(1:3,ia) + atom_pos_buffer(1:3,ja))
          call test_fragment_inversion_nonperiodic_import(dc, atom_index_buffer, atom_pos_buffer, natom_buffer, &
            origin, accepted, residual)
          if(accepted .and. residual < best_residual) then
            found_inversion = .true.
            best_residual = residual
            best_origin = origin
          end if
        end do
      end do
      if(found_inversion) then
        write(*,'(1x,a,i0,a,es12.5,a,i0,a,3(1x,es12.5),a,3(1x,es12.5))') &
          "[DC-LCFO-W90-SYM] fragment=", ifrag, " buffer finite_box detected label=inversion", &
          best_residual, " natom_buffer=", natom_buffer, " origin=", best_origin(1:3), &
          " box=", buffer_length(1:3)
      else
        write(*,'(1x,a,i0,a,i0,a,3(1x,es12.5))') &
          "[DC-LCFO-W90-SYM] fragment=", ifrag, &
          " buffer finite_box detected label=none natom_buffer=", natom_buffer, " box=", buffer_length(1:3)
      end if

      if(allocated(atom_index)) deallocate(atom_index)
      if(allocated(atom_index_buffer)) deallocate(atom_index_buffer)
      if(allocated(atom_pos_local)) deallocate(atom_pos_local)
      if(allocated(atom_pos_buffer)) deallocate(atom_pos_buffer)
    end do
  end subroutine detect_wannier_fragment_symops

  subroutine collect_fragment_core_atoms_import(dc, ifrag, atom_index, natom_frag)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag
    integer, allocatable, intent(out) :: atom_index(:)
    integer, intent(out) :: natom_frag
    integer :: ia

    natom_frag = 0
    allocate(atom_index(max(1,dc%system_tot%nion)))
    do ia=1,dc%system_tot%nion
      if(atom_in_fragment_core_import(dc, ifrag, ia)) then
        natom_frag = natom_frag + 1
        atom_index(natom_frag) = ia
      end if
    end do
  end subroutine collect_fragment_core_atoms_import

  subroutine collect_fragment_buffer_atoms_import(dc, ifrag, atom_index, atom_pos_local, natom_buffer)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag
    integer, allocatable, intent(out) :: atom_index(:)
    real(8), allocatable, intent(out) :: atom_pos_local(:,:)
    integer, intent(out) :: natom_buffer
    integer :: ia, sx, sy, sz, axis, nmax
    real(8) :: total_length(3), fragment_origin(3), buffer_shift(3), box_length(3)
    real(8) :: shifted(3), local_pos(3)

    nmax = max(1, 27 * dc%system_tot%nion)
    allocate(atom_index(nmax), atom_pos_local(3,nmax))
    natom_buffer = 0
    call total_cell_lengths_import(dc, total_length)
    call fragment_buffer_shift_import(dc, buffer_shift)
    call fragment_buffer_lengths_import(dc, ifrag, box_length)
    do axis=1,3
      fragment_origin(axis) = dc%lg_tot%coordinate(dc%ixyz_frag(axis,ifrag),axis) - buffer_shift(axis)
    end do

    do ia=1,dc%system_tot%nion
      do sx=-1,1
      do sy=-1,1
      do sz=-1,1
        shifted(1) = dc%system_tot%rion(1,ia) + dble(sx) * total_length(1)
        shifted(2) = dc%system_tot%rion(2,ia) + dble(sy) * total_length(2)
        shifted(3) = dc%system_tot%rion(3,ia) + dble(sz) * total_length(3)
        local_pos(1:3) = shifted(1:3) - fragment_origin(1:3)
        if(all(local_pos(1:3) >= -1.0d-8) .and. all(local_pos(1:3) <= box_length(1:3) + 1.0d-8)) then
          natom_buffer = natom_buffer + 1
          if(natom_buffer > nmax) stop "DC-LCFO-W90-SYM: buffer atom list overflow"
          atom_index(natom_buffer) = ia
          atom_pos_local(1:3,natom_buffer) = local_pos(1:3)
        end if
      end do
      end do
      end do
    end do
  end subroutine collect_fragment_buffer_atoms_import

  logical function atom_in_fragment_core_import(dc, ifrag, ia) result(in_core)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag, ia
    integer :: axis, i, ig, idx0, nxyz_domain(3), ibest
    real(8) :: r_atom, r_grid, dist, best_dist, length_axis, spacing_axis

    in_core = .true.
    call get_fragment_domain(dc, ifrag, nxyz_domain)
    do axis=1,3
      r_atom = dc%system_tot%rion(axis,ia)
      length_axis = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
        + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
      spacing_axis = length_axis / dble(dc%lg_tot%num(axis))
      idx0 = dc%ixyz_frag(axis,ifrag)
      best_dist = huge(1d0)
      ibest = 0
      do i=1,nxyz_domain(axis)
        ig = idx0 + i - 1
        r_grid = dc%lg_tot%coordinate(ig,axis)
        dist = abs(periodic_delta_import(r_grid - r_atom, length_axis))
        if(dist < best_dist) then
          best_dist = dist
          ibest = i
        end if
      end do
      if(ibest <= 0 .or. best_dist > 0.75d0 * spacing_axis) in_core = .false.
    end do
  end function atom_in_fragment_core_import

  subroutine test_fragment_inversion_import(dc, atom_index, atom_pos_local, natom_frag, cell_length, origin, &
      accepted, max_residual)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: natom_frag, atom_index(:)
    real(8), intent(in) :: atom_pos_local(3,*)
    real(8), intent(in) :: cell_length(3), origin(3)
    logical, intent(out) :: accepted
    real(8), intent(out) :: max_residual
    integer :: ia, ja, axis, iatom, jatom
    real(8) :: target(3), delta(3), dist2, best_dist2
    real(8), parameter :: atom_match_tol2 = 1.0d-6

    accepted = .false.
    max_residual = huge(1d0)
    if(natom_frag <= 0) return

    max_residual = 0.0d0
    do ia=1,natom_frag
      iatom = atom_index(ia)
      do axis=1,3
        target(axis) = origin(axis) - periodic_delta_import(atom_pos_local(axis,ia) - origin(axis), cell_length(axis))
      end do
      best_dist2 = huge(1d0)
      do ja=1,natom_frag
        jatom = atom_index(ja)
        if(dc%system_tot%kion(jatom) /= dc%system_tot%kion(iatom)) cycle
        do axis=1,3
          delta(axis) = periodic_delta_import(atom_pos_local(axis,ja) - target(axis), cell_length(axis))
        end do
        dist2 = local_distance2(delta)
        best_dist2 = min(best_dist2, dist2)
      end do
      max_residual = max(max_residual, sqrt(best_dist2))
      if(best_dist2 > atom_match_tol2) return
    end do

    accepted = .true.
  end subroutine test_fragment_inversion_import

  subroutine test_fragment_inversion_nonperiodic_import(dc, atom_index, atom_pos_local, natom_frag, origin, &
      accepted, max_residual)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: natom_frag, atom_index(:)
    real(8), intent(in) :: atom_pos_local(3,*), origin(3)
    logical, intent(out) :: accepted
    real(8), intent(out) :: max_residual
    integer :: ia, ja, iatom, jatom
    real(8) :: target(3), delta(3), dist2, best_dist2
    real(8), parameter :: atom_match_tol2 = 1.0d-6

    accepted = .false.
    max_residual = huge(1d0)
    if(natom_frag <= 0) return

    max_residual = 0.0d0
    do ia=1,natom_frag
      iatom = atom_index(ia)
      target(1:3) = 2.0d0 * origin(1:3) - atom_pos_local(1:3,ia)
      best_dist2 = huge(1d0)
      do ja=1,natom_frag
        jatom = atom_index(ja)
        if(dc%system_tot%kion(jatom) /= dc%system_tot%kion(iatom)) cycle
        delta(1:3) = atom_pos_local(1:3,ja) - target(1:3)
        dist2 = local_distance2(delta)
        best_dist2 = min(best_dist2, dist2)
      end do
      max_residual = max(max_residual, sqrt(best_dist2))
      if(best_dist2 > atom_match_tol2) return
    end do

    accepted = .true.
  end subroutine test_fragment_inversion_nonperiodic_import

  subroutine test_fragment_inversion_global_pbc_import(dc, ifrag, atom_pos_local, natom_frag, origin, &
      accepted, max_residual)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag, natom_frag
    real(8), intent(in) :: atom_pos_local(3,*), origin(3)
    logical, intent(out) :: accepted
    real(8), intent(out) :: max_residual
    integer :: ia, ja, axis
    real(8) :: total_length(3), fragment_origin(3), buffer_shift(3)
    real(8) :: target_global(3), delta(3), dist2, best_dist2
    real(8), parameter :: atom_match_tol2 = 1.0d-6

    accepted = .false.
    max_residual = huge(1d0)
    if(natom_frag <= 0) return

    call total_cell_lengths_import(dc, total_length)
    call fragment_buffer_shift_import(dc, buffer_shift)
    do axis=1,3
      fragment_origin(axis) = dc%lg_tot%coordinate(dc%ixyz_frag(axis,ifrag),axis) - buffer_shift(axis)
    end do

    max_residual = 0.0d0
    do ia=1,natom_frag
      target_global(1:3) = fragment_origin(1:3) + 2.0d0 * origin(1:3) - atom_pos_local(1:3,ia)
      do axis=1,3
        if(total_length(axis) > 0.0d0) then
          target_global(axis) = target_global(axis) - floor(target_global(axis) / total_length(axis)) * total_length(axis)
        end if
      end do
      best_dist2 = huge(1d0)
      do ja=1,dc%system_tot%nion
        do axis=1,3
          delta(axis) = periodic_delta_import(dc%system_tot%rion(axis,ja) - target_global(axis), total_length(axis))
        end do
        dist2 = local_distance2(delta)
        best_dist2 = min(best_dist2, dist2)
      end do
      max_residual = max(max_residual, sqrt(best_dist2))
      if(best_dist2 > atom_match_tol2) return
    end do

    accepted = .true.
  end subroutine test_fragment_inversion_global_pbc_import

  subroutine fragment_cell_lengths_import(dc, ifrag, cell_length)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag
    real(8), intent(out) :: cell_length(3)
    integer :: axis, nxyz_domain(3)
    real(8) :: total_length

    call get_fragment_domain(dc, ifrag, nxyz_domain)
    do axis=1,3
      total_length = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
        + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
      cell_length(axis) = total_length * dble(nxyz_domain(axis)) / dble(dc%lg_tot%num(axis))
    end do
  end subroutine fragment_cell_lengths_import

  subroutine fragment_buffer_lengths_import(dc, ifrag, buffer_length)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag
    real(8), intent(out) :: buffer_length(3)
    real(8) :: cell_length(3), hgrid(3)

    call fragment_cell_lengths_import(dc, ifrag, cell_length)
    call grid_spacings_import(dc, hgrid)
    buffer_length(1:3) = cell_length(1:3) + 2.0d0 * hgrid(1:3) * dble(dc%nxyz_buffer(1:3))
  end subroutine fragment_buffer_lengths_import

  subroutine fragment_buffer_shift_import(dc, buffer_shift)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    real(8), intent(out) :: buffer_shift(3)
    real(8) :: hgrid(3)

    call grid_spacings_import(dc, hgrid)
    buffer_shift(1:3) = hgrid(1:3) * dble(dc%nxyz_buffer(1:3))
  end subroutine fragment_buffer_shift_import

  subroutine grid_spacings_import(dc, hgrid)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    real(8), intent(out) :: hgrid(3)
    real(8) :: total_length(3)

    call total_cell_lengths_import(dc, total_length)
    hgrid(1:3) = total_length(1:3) / dble(dc%lg_tot%num(1:3))
  end subroutine grid_spacings_import

  subroutine total_cell_lengths_import(dc, total_length)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    real(8), intent(out) :: total_length(3)
    integer :: axis

    do axis=1,3
      total_length(axis) = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
        + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
    end do
  end subroutine total_cell_lengths_import

  subroutine atom_fragment_local_position_import(dc, ifrag, ia, local_pos)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag, ia
    real(8), intent(out) :: local_pos(3)
    integer :: axis, idx0
    real(8) :: total_length, fragment_origin

    do axis=1,3
      idx0 = dc%ixyz_frag(axis,ifrag)
      fragment_origin = dc%lg_tot%coordinate(idx0,axis)
      total_length = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
        + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
      local_pos(axis) = periodic_delta_import(dc%system_tot%rion(axis,ia) - fragment_origin, total_length)
      if(local_pos(axis) < 0.0d0) local_pos(axis) = local_pos(axis) + total_length
    end do
  end subroutine atom_fragment_local_position_import

  subroutine fragment_center_bohr_import(dc, ifrag, center_bohr)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag
    real(8), intent(out) :: center_bohr(3)
    integer :: axis, idx0, idx1, nxyz_domain(3)

    call get_fragment_domain(dc, ifrag, nxyz_domain)
    do axis=1,3
      idx0 = dc%ixyz_frag(axis,ifrag)
      idx1 = idx0 + nxyz_domain(axis) - 1
      center_bohr(axis) = 0.5d0 * (dc%lg_tot%coordinate(idx0,axis) + dc%lg_tot%coordinate(idx1,axis))
    end do
  end subroutine fragment_center_bohr_import

  subroutine periodic_pair_midpoint_import(ri, rj, cell_length, midpoint)
    implicit none
    real(8), intent(in) :: ri(3), rj(3), cell_length(3)
    real(8), intent(out) :: midpoint(3)
    integer :: axis

    do axis=1,3
      midpoint(axis) = ri(axis) + 0.5d0 * periodic_delta_import(rj(axis) - ri(axis), cell_length(axis))
      if(cell_length(axis) > 0.0d0) midpoint(axis) = midpoint(axis) - floor(midpoint(axis) / cell_length(axis)) * cell_length(axis)
    end do
  end subroutine periodic_pair_midpoint_import

  subroutine diagnose_fragment_wannier_center_symmetry(dc, center_bohr, owner_frag, num_wann, nsym, symops)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_wann, owner_frag(num_wann), nsym
    real(8), intent(in) :: center_bohr(3,num_wann)
    type(t_wannier_symop), intent(in) :: symops(:)
    integer :: isym, iw, jw, axis, ifrag, nowned, jbest
    real(8), allocatable :: center_local(:,:)
    real(8) :: cell_length(3), relative(3), mapped(3), delta(3)
    real(8) :: dist2, best_dist2, max_residual, rms_residual
    logical :: closed
    real(8), parameter :: center_match_tol = 1.0d-3

    if(nsym <= 0) return
    allocate(center_local(3,num_wann))

    do isym=1,nsym
      ifrag = symops(isym)%owner_frag
      if(ifrag <= 0) cycle
      call fragment_cell_lengths_import(dc, ifrag, cell_length)
      nowned = count(owner_frag(1:num_wann) == ifrag)
      if(nowned <= 0) cycle

      do iw=1,num_wann
        if(owner_frag(iw) /= ifrag) cycle
        call point_fragment_local_position_import(dc, ifrag, center_bohr(1:3,iw), center_local(1:3,iw))
      end do

      max_residual = 0.0d0
      rms_residual = 0.0d0
      closed = .true.
      do iw=1,num_wann
        if(owner_frag(iw) /= ifrag) cycle
        relative(1:3) = center_local(1:3,iw) - symops(isym)%origin_bohr(1:3)
        mapped(1:3) = symops(isym)%origin_bohr(1:3) + matmul_int_real(symops(isym)%rot, relative) &
          + symops(isym)%tau_local(1:3)
        call wrap_fragment_local_point(mapped, cell_length)

        best_dist2 = huge(1d0)
        jbest = 0
        do jw=1,num_wann
          if(owner_frag(jw) /= ifrag) cycle
          do axis=1,3
            delta(axis) = periodic_delta_import(center_local(axis,jw) - mapped(axis), cell_length(axis))
          end do
          dist2 = local_distance2(delta)
          if(dist2 < best_dist2) then
            best_dist2 = dist2
            jbest = jw
          end if
        end do
        if(jbest <= 0) then
          closed = .false.
        else
          max_residual = max(max_residual, sqrt(best_dist2))
          rms_residual = rms_residual + best_dist2
          if(sqrt(best_dist2) > center_match_tol) closed = .false.
        end if
      end do
      rms_residual = sqrt(rms_residual / dble(max(1,nowned)))
      write(*,'(1x,a,i0,2a,2(a,es12.5),a,l1,a,i0)') &
        "[DC-LCFO-W90-SYM] fragment=", ifrag, " center closure label=", &
        trim(symops(isym)%label), " max=", max_residual, " rms=", rms_residual, &
        " closed=", closed, " ncenter=", nowned
    end do

    deallocate(center_local)
  end subroutine diagnose_fragment_wannier_center_symmetry

  subroutine diagnose_local_bond_center_orbit_closure(dc, ncenter, center_bohr)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ncenter
    real(8), intent(in) :: center_bohr(3,*)

    call diagnose_local_center_orbit_closure(dc, ncenter, center_bohr, 'seed_bond')
  end subroutine diagnose_local_bond_center_orbit_closure

  subroutine diagnose_local_wannier_center_orbit_closure(dc, ncenter, center_bohr, label)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ncenter
    real(8), intent(in) :: center_bohr(3,*)
    character(*), intent(in) :: label

    call diagnose_local_center_orbit_closure(dc, ncenter, center_bohr, label)
  end subroutine diagnose_local_wannier_center_orbit_closure

  subroutine diagnose_local_center_orbit_closure(dc, ncenter, center_bohr, label)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ncenter
    real(8), intent(in) :: center_bohr(3,*)
    character(*), intent(in) :: label
    integer :: ia, ja, ic, jc, axis, ifrag
    integer :: natom_frag
    integer, allocatable :: atom_index(:)
    real(8), allocatable :: atom_pos_local(:,:), center_local(:,:)
    real(8) :: cell_length(3), origin(3), best_origin(3)
    real(8) :: residual, best_residual, target(3), delta(3)
    real(8) :: dist2, best_dist2, max_residual, rms_residual
    logical :: accepted, found_inversion, closed
    real(8), parameter :: center_match_tol = 1.0d-3

    ifrag = dc%i_frag
    if(ncenter <= 0 .or. ifrag <= 0) return

    call collect_fragment_core_atoms_import(dc, ifrag, atom_index, natom_frag)
    call fragment_cell_lengths_import(dc, ifrag, cell_length)
    allocate(atom_pos_local(3,max(1,natom_frag)))
    do ia=1,natom_frag
      call atom_fragment_local_position_import(dc, ifrag, atom_index(ia), atom_pos_local(1:3,ia))
    end do

    found_inversion = .false.
    best_residual = huge(1.0d0)
    best_origin = 0.0d0
    do ia=1,natom_frag
      do ja=ia,natom_frag
        call periodic_pair_midpoint_import(atom_pos_local(1:3,ia), atom_pos_local(1:3,ja), &
          cell_length, origin)
        call test_fragment_inversion_import(dc, atom_index, atom_pos_local, natom_frag, &
          cell_length, origin, accepted, residual)
        if(accepted .and. residual < best_residual) then
          found_inversion = .true.
          best_residual = residual
          best_origin = origin
        end if
      end do
    end do

    if(.not. found_inversion) then
      write(*,'(1x,a,i0,3a,i0)') "[DC-LCFO-LOCAL-WANNIER-SYM] fragment=", ifrag, &
        " label=", trim(label), " symop=none ncenter=", ncenter
      deallocate(atom_index, atom_pos_local)
      return
    end if

    allocate(center_local(3,ncenter))
    do ic=1,ncenter
      call point_fragment_local_position_import(dc, ifrag, center_bohr(1:3,ic), center_local(1:3,ic))
    end do

    max_residual = 0.0d0
    rms_residual = 0.0d0
    closed = .true.
    do ic=1,ncenter
      do axis=1,3
        target(axis) = best_origin(axis) - periodic_delta_import(center_local(axis,ic) - best_origin(axis), &
          cell_length(axis))
      end do
      best_dist2 = huge(1.0d0)
      do jc=1,ncenter
        do axis=1,3
          delta(axis) = periodic_delta_import(center_local(axis,jc) - target(axis), cell_length(axis))
        end do
        dist2 = local_distance2(delta)
        best_dist2 = min(best_dist2, dist2)
      end do
      max_residual = max(max_residual, sqrt(best_dist2))
      rms_residual = rms_residual + best_dist2
      if(sqrt(best_dist2) > center_match_tol) closed = .false.
    end do
    rms_residual = sqrt(rms_residual / dble(max(1,ncenter)))

    write(*,'(1x,a,i0,3a,2(a,es12.5),a,l1,a,i0,a,3(1x,es12.5))') &
      "[DC-LCFO-LOCAL-WANNIER-SYM] fragment=", ifrag, " label=", trim(label), &
      " symop=inversion", " max=", max_residual, " rms=", rms_residual, &
      " closed=", closed, " ncenter=", ncenter, " origin=", best_origin(1:3)

    deallocate(atom_index, atom_pos_local, center_local)
  end subroutine diagnose_local_center_orbit_closure

  subroutine point_fragment_local_position_import(dc, ifrag, point_bohr, local_pos)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag
    real(8), intent(in) :: point_bohr(3)
    real(8), intent(out) :: local_pos(3)
    integer :: axis, idx0
    real(8) :: total_length, fragment_origin, cell_length(3)

    call fragment_cell_lengths_import(dc, ifrag, cell_length)
    do axis=1,3
      idx0 = dc%ixyz_frag(axis,ifrag)
      fragment_origin = dc%lg_tot%coordinate(idx0,axis)
      total_length = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
        + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
      local_pos(axis) = periodic_delta_import(point_bohr(axis) - fragment_origin, total_length)
      if(local_pos(axis) < 0.0d0) local_pos(axis) = local_pos(axis) + total_length
      if(cell_length(axis) > 0.0d0) local_pos(axis) = local_pos(axis) - floor(local_pos(axis) / cell_length(axis)) * cell_length(axis)
    end do
  end subroutine point_fragment_local_position_import

  subroutine fragment_symmetry_origin_global_import(dc, ifrag, symop, origin_global)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag
    type(t_wannier_symop), intent(in) :: symop
    real(8), intent(out) :: origin_global(3)
    integer :: axis, idx0
    real(8) :: total_length, fragment_origin

    do axis=1,3
      idx0 = dc%ixyz_frag(axis,ifrag)
      fragment_origin = dc%lg_tot%coordinate(idx0,axis)
      total_length = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
        + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
      origin_global(axis) = fragment_origin + symop%origin_bohr(axis)
      if(total_length > 0.0d0) origin_global(axis) = origin_global(axis) &
        - floor(origin_global(axis) / total_length) * total_length
    end do
  end subroutine fragment_symmetry_origin_global_import

  subroutine wrap_fragment_local_point(point, cell_length)
    implicit none
    real(8), intent(inout) :: point(3)
    real(8), intent(in) :: cell_length(3)
    integer :: axis

    do axis=1,3
      if(cell_length(axis) <= 0.0d0) cycle
      point(axis) = point(axis) - floor(point(axis) / cell_length(axis)) * cell_length(axis)
    end do
  end subroutine wrap_fragment_local_point

  subroutine diagnose_fragment_wannier_symmetry_representation(dc, center_bohr, owner_frag, num_wann, &
      num_bands, v_matrix, esp_file, nsym, symops, aa_global, position_available)
    use eigen_subdiag_sub, only: eigen_zheev
    use filesystem, only: get_filehandle
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_wann, owner_frag(num_wann), num_bands, nsym
    real(8), intent(in) :: center_bohr(3,num_wann), esp_file(:,:)
    complex(8), intent(in) :: v_matrix(num_bands,num_wann)
    complex(8), intent(in) :: aa_global(3,num_wann,num_wann)
    logical, intent(in) :: position_available
    type(t_wannier_symop), intent(in) :: symops(:)
    integer :: ifrag, isym, nowned, iw, ib, npts, nspin_file, nstate_frag_file, nstate_tot_file
    integer :: nxyz_domain(3), n_basis_frag, p, ix, iy, iz, jb, jw, axis, nocc
    integer :: n_s2_gt_09, n_s2_gt_05
    integer, allocatable :: wann_index(:), pmap(:), center_perm(:)
    real(8), allocatable :: phi_basis(:,:), coef_wf(:,:), psi_state(:,:)
    complex(8), allocatable :: wannier_frag(:,:), srep(:,:), gram(:,:), eigvec(:,:), polar(:,:), invsqrt(:,:)
    complex(8), allocatable :: zloc(:,:), hloc(:,:), rholoc(:,:), work(:,:)
    real(8), allocatable :: eval(:)
    real(8) :: hvol, min_eval, max_eval, gram_residual, polar_residual, map_residual, eval_floor
    real(8) :: subspace_leakage, closure_tol
    real(8) :: center_perm_max, center_perm_rms, target_residual
    real(8) :: origin_global(3), z_odd_res2, z_norm2, h_even_res, rho_even_res
    logical :: map_ok, seed_ok, center_perm_ok

    if(nsym <= 0) return
    hvol = dc%system_tot%hvol
    eval_floor = 1.0d-10
    closure_tol = 1.0d-3
    do ifrag=1,dc%n_frag
      nowned = count(owner_frag(1:num_wann) == ifrag)
      if(nowned <= 0) cycle
      allocate(wann_index(nowned))
      ib = 0
      do iw=1,num_wann
        if(owner_frag(iw) /= ifrag) cycle
        ib = ib + 1
        wann_index(ib) = iw
      end do

      call read_fragment_lcfo_seed_for_wannier_import(dc, ifrag, num_bands, &
        nxyz_domain, nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag, &
        phi_basis, coef_wf, seed_ok)
      if(.not. seed_ok) then
        write(*,'(1x,a,i0,a)') "[DC-LCFO-W90-SYM] fragment=", ifrag, &
          " representation skipped: missing LCFO seed"
        deallocate(wann_index)
        cycle
      end if

      npts = product(nxyz_domain)
      allocate(psi_state(npts,num_bands), wannier_frag(npts,nowned))
      psi_state = matmul(phi_basis(1:npts,1:n_basis_frag), coef_wf(1:n_basis_frag,1:num_bands))
      wannier_frag = matmul(cmplx(psi_state(1:npts,1:num_bands), 0.0d0, kind=8), &
        v_matrix(1:num_bands,wann_index(1:nowned)))

      do isym=1,nsym
        if(symops(isym)%owner_frag /= ifrag) cycle
        allocate(pmap(npts))
        call build_fragment_symmetry_grid_map(dc, ifrag, nxyz_domain, symops(isym), &
          pmap, map_ok, map_residual)
        if(.not. map_ok) then
          write(*,'(1x,a,i0,2a,a,es12.5)') &
            "[DC-LCFO-W90-SYM] fragment=", ifrag, " representation label=", &
            trim(symops(isym)%label), " skipped: grid map residual=", map_residual
          deallocate(pmap)
          cycle
        end if

        allocate(center_perm(nowned))
        call build_fragment_center_permutation_import(dc, ifrag, center_bohr, owner_frag, num_wann, &
          symops(isym), wann_index, nowned, center_perm, center_perm_ok, center_perm_max, center_perm_rms)

        allocate(srep(nowned,nowned), gram(nowned,nowned), eigvec(nowned,nowned), &
          polar(nowned,nowned), invsqrt(nowned,nowned), eval(nowned))
        srep = (0.0d0,0.0d0)
        do iw=1,nowned
          do ib=1,nowned
            do p=1,npts
              srep(iw,ib) = srep(iw,ib) + conjg(wannier_frag(p,iw)) * wannier_frag(pmap(p),ib) * hvol
            end do
          end do
        end do

        gram = matmul(conjg(transpose(srep)), srep)
        call eigen_zheev(gram, eval, eigvec)
        min_eval = minval(eval)
        max_eval = maxval(eval)
        subspace_leakage = sqrt(max(0.0d0, 1.0d0 - min(1.0d0, min_eval)))
        n_s2_gt_09 = count(eval(1:nowned) > 0.9d0)
        n_s2_gt_05 = count(eval(1:nowned) > 0.5d0)
        gram_residual = hermitian_identity_residual(gram)
        polar_residual = -1.0d0
        target_residual = -1.0d0
        invsqrt = (0.0d0,0.0d0)
        polar = (0.0d0,0.0d0)
        do iw=1,nowned
          do ib=1,nowned
            invsqrt(iw,ib) = sum(eigvec(iw,1:nowned) * merge(1.0d0 / sqrt(eval(1:nowned)), 0.0d0, &
              eval(1:nowned) > eval_floor) * conjg(eigvec(ib,1:nowned)))
            polar(iw,ib) = sum(eigvec(iw,1:nowned) * merge(0.0d0, 1.0d0, &
              eval(1:nowned) > eval_floor) * conjg(eigvec(ib,1:nowned)))
          end do
        end do
        polar = matmul(srep, invsqrt) + polar
        polar_residual = hermitian_identity_residual(matmul(conjg(transpose(polar)), polar))
        if(center_perm_ok) target_residual = permutation_target_residual(polar, center_perm)
        write(*,'(1x,a,i0,2a,9(a,es12.5),2(a,l1),3(a,i0))') &
          "[DC-LCFO-W90-SYM] fragment=", ifrag, " representation label=", &
          trim(symops(isym)%label), " min_s2=", min_eval, " max_s2=", max_eval, &
          " subspace_leakage=", subspace_leakage, &
          " s_unit_res=", gram_residual, " polar_unit_res=", polar_residual, &
          " map_res=", map_residual, " center_perm_max=", center_perm_max, &
          " center_perm_rms=", center_perm_rms, " target_res=", target_residual, &
          " subspace_closed=", subspace_leakage <= closure_tol, &
          " center_perm_ok=", center_perm_ok, " n_s2_gt_09=", n_s2_gt_09, &
          " n_s2_gt_05=", n_s2_gt_05, " ncenter=", nowned

        if(trim(symops(isym)%label) == 'inversion') then
          call fragment_symmetry_origin_global_import(dc, ifrag, symops(isym), origin_global)
          allocate(zloc(nowned,nowned), hloc(nowned,nowned), rholoc(nowned,nowned), &
            work(nowned,nowned))
          hloc = (0.0d0,0.0d0)
          rholoc = (0.0d0,0.0d0)
          nocc = min(num_bands, max(1, num_wann / 2))
          do jb=1,nowned
            jw = wann_index(jb)
            do ib=1,nowned
              iw = wann_index(ib)
              hloc(ib,jb) = sum(conjg(v_matrix(1:num_bands,iw)) * &
                cmplx(esp_file(1:num_bands,1), 0.0d0, kind=8) * v_matrix(1:num_bands,jw))
              rholoc(ib,jb) = sum(conjg(v_matrix(1:nocc,iw)) * v_matrix(1:nocc,jw))
            end do
          end do
          work = matmul(conjg(transpose(polar)), matmul(hloc, polar)) - hloc
          h_even_res = sqrt(sum(abs(work)**2) / max(sum(abs(hloc)**2), 1.0d-300))
          work = matmul(conjg(transpose(polar)), matmul(rholoc, polar)) - rholoc
          rho_even_res = sqrt(sum(abs(work)**2) / max(sum(abs(rholoc)**2), 1.0d-300))
          z_odd_res2 = 0.0d0
          z_norm2 = 0.0d0
          if(position_available) then
            do axis=1,3
              do jb=1,nowned
                jw = wann_index(jb)
                do ib=1,nowned
                  iw = wann_index(ib)
                  zloc(ib,jb) = aa_global(axis,iw,jw)
                end do
                zloc(jb,jb) = zloc(jb,jb) - cmplx(origin_global(axis), 0.0d0, kind=8)
              end do
              work = matmul(conjg(transpose(polar)), matmul(zloc, polar)) + zloc
              z_odd_res2 = z_odd_res2 + sum(abs(work)**2)
              z_norm2 = z_norm2 + sum(abs(zloc)**2)
            end do
          end if
          write(*,'(1x,a,i0,a,5(a,es12.5),a,i0)') &
            "[DC-LCFO-W90-SYM-OP] fragment=", ifrag, " label=inversion", &
            " z_odd_res=", sqrt(z_odd_res2 / max(z_norm2, 1.0d-300)), &
            " h_even_res=", h_even_res, " rho_even_res=", rho_even_res, &
            " polar_unit_res=", polar_residual, " map_res=", map_residual, " nocc=", nocc
          deallocate(zloc, hloc, rholoc, work)
        end if

        deallocate(srep, gram, eigvec, polar, invsqrt, eval, pmap, center_perm)
      end do

      deallocate(wann_index, phi_basis, coef_wf, psi_state, wannier_frag)
    end do
  end subroutine diagnose_fragment_wannier_symmetry_representation

  subroutine diagnose_global_wannier_pbc_operator_symmetry(dc, center_bohr, num_wann, &
      num_bands, v_matrix, esp_file, nsym, symops, aa_global, position_available)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_wann, num_bands, nsym
    real(8), intent(in) :: center_bohr(3,num_wann), esp_file(:,:)
    complex(8), intent(in) :: v_matrix(num_bands,num_wann)
    complex(8), intent(in) :: aa_global(3,num_wann,num_wann)
    logical, intent(in) :: position_available
    type(t_wannier_symop), intent(in) :: symops(:)
    integer :: isym, iw, jw, ib, axis, nocc
    integer, allocatable :: perm(:)
    real(8) :: origin_global(3), total_length(3), perm_max, perm_rms
    real(8) :: h_norm2, h_res2, rho_norm2, rho_res2, z_norm2, z_res2
    complex(8), allocatable :: hmat(:,:), rhomat(:,:), zmat(:,:)
    complex(8) :: diff
    logical :: perm_ok

    if(nsym <= 0) return
    call total_cell_lengths_import(dc, total_length)
    allocate(perm(num_wann))
    allocate(hmat(num_wann,num_wann), rhomat(num_wann,num_wann), zmat(num_wann,num_wann))

    hmat = (0.0d0,0.0d0)
    rhomat = (0.0d0,0.0d0)
    nocc = min(num_bands, max(1, num_wann / 2))
    do jw=1,num_wann
      do iw=1,num_wann
        hmat(iw,jw) = sum(conjg(v_matrix(1:num_bands,iw)) * &
          cmplx(esp_file(1:num_bands,1), 0.0d0, kind=8) * v_matrix(1:num_bands,jw))
        rhomat(iw,jw) = sum(conjg(v_matrix(1:nocc,iw)) * v_matrix(1:nocc,jw))
      end do
    end do

    do isym=1,nsym
      if(trim(symops(isym)%label) /= 'inversion') cycle
      call fragment_symmetry_origin_global_import(dc, symops(isym)%owner_frag, symops(isym), origin_global)
      call build_global_center_pbc_permutation_import(center_bohr, num_wann, total_length, &
        origin_global, perm, perm_ok, perm_max, perm_rms)
      write(*,'(1x,a,i0,2a,2(a,es12.5),a,l1,a,i0)') &
        "[DC-LCFO-W90-SYM-GLOBAL] origin_frag=", symops(isym)%owner_frag, &
        " label=", trim(symops(isym)%label), " perm_max=", perm_max, &
        " perm_rms=", perm_rms, " perm_ok=", perm_ok, " ncenter=", num_wann
      if(.not. perm_ok) cycle

      h_norm2 = 0.0d0
      h_res2 = 0.0d0
      rho_norm2 = 0.0d0
      rho_res2 = 0.0d0
      do jw=1,num_wann
        do iw=1,num_wann
          diff = hmat(perm(iw),perm(jw)) - hmat(iw,jw)
          h_res2 = h_res2 + abs(diff)**2
          h_norm2 = h_norm2 + abs(hmat(iw,jw))**2
          diff = rhomat(perm(iw),perm(jw)) - rhomat(iw,jw)
          rho_res2 = rho_res2 + abs(diff)**2
          rho_norm2 = rho_norm2 + abs(rhomat(iw,jw))**2
        end do
      end do

      z_norm2 = 0.0d0
      z_res2 = 0.0d0
      if(position_available) then
        do axis=1,3
          zmat(1:num_wann,1:num_wann) = aa_global(axis,1:num_wann,1:num_wann)
          do ib=1,num_wann
            zmat(ib,ib) = zmat(ib,ib) - cmplx(origin_global(axis), 0.0d0, kind=8)
          end do
          do jw=1,num_wann
            do iw=1,num_wann
              diff = zmat(perm(iw),perm(jw)) + zmat(iw,jw)
              z_res2 = z_res2 + abs(diff)**2
              z_norm2 = z_norm2 + abs(zmat(iw,jw))**2
            end do
          end do
        end do
      end if

      if(position_available) then
        write(*,'(1x,a,i0,a,3(a,es12.5),a,i0)') &
          "[DC-LCFO-W90-SYM-GLOBAL-OP] origin_frag=", symops(isym)%owner_frag, &
          " label=inversion", " z_odd_res=", sqrt(z_res2 / max(z_norm2, 1.0d-300)), &
          " h_even_res=", sqrt(h_res2 / max(h_norm2, 1.0d-300)), &
          " rho_even_res=", sqrt(rho_res2 / max(rho_norm2, 1.0d-300)), " nocc=", nocc
      else
        write(*,'(1x,a,i0,2a,2(a,es12.5),a,i0)') &
          "[DC-LCFO-W90-SYM-GLOBAL-OP] origin_frag=", symops(isym)%owner_frag, &
          " label=inversion", " z_odd_available=F", &
          " h_even_res=", sqrt(h_res2 / max(h_norm2, 1.0d-300)), &
          " rho_even_res=", sqrt(rho_res2 / max(rho_norm2, 1.0d-300)), " nocc=", nocc
      end if
    end do

    deallocate(perm, hmat, rhomat, zmat)
  end subroutine diagnose_global_wannier_pbc_operator_symmetry

  subroutine build_global_center_pbc_permutation_import(center_bohr, num_wann, total_length, &
      origin_global, perm, ok, max_residual, rms_residual)
    implicit none
    integer, intent(in) :: num_wann
    real(8), intent(in) :: center_bohr(3,num_wann), total_length(3), origin_global(3)
    integer, intent(out) :: perm(num_wann)
    logical, intent(out) :: ok
    real(8), intent(out) :: max_residual, rms_residual
    integer :: iw, jw, jbest, axis
    integer, allocatable :: used(:)
    real(8) :: mapped(3), delta(3), dist2, best_dist2
    real(8), parameter :: center_match_tol = 1.0d-3

    allocate(used(num_wann))
    used = 0
    perm = 0
    ok = .true.
    max_residual = 0.0d0
    rms_residual = 0.0d0

    do iw=1,num_wann
      do axis=1,3
        mapped(axis) = origin_global(axis) - periodic_delta_import(center_bohr(axis,iw) - origin_global(axis), &
          total_length(axis))
        if(total_length(axis) > 0.0d0) mapped(axis) = mapped(axis) - floor(mapped(axis) / total_length(axis)) &
          * total_length(axis)
      end do

      best_dist2 = huge(1.0d0)
      jbest = 0
      do jw=1,num_wann
        if(used(jw) /= 0) cycle
        do axis=1,3
          delta(axis) = periodic_delta_import(center_bohr(axis,jw) - mapped(axis), total_length(axis))
        end do
        dist2 = local_distance2(delta)
        if(dist2 < best_dist2) then
          best_dist2 = dist2
          jbest = jw
        end if
      end do

      if(jbest <= 0) then
        ok = .false.
      else
        perm(iw) = jbest
        used(jbest) = 1
        max_residual = max(max_residual, sqrt(best_dist2))
        rms_residual = rms_residual + best_dist2
        if(sqrt(best_dist2) > center_match_tol) ok = .false.
      end if
    end do

    rms_residual = sqrt(rms_residual / dble(max(1,num_wann)))
    if(any(perm(1:num_wann) <= 0)) ok = .false.
    deallocate(used)
  end subroutine build_global_center_pbc_permutation_import

  subroutine symmetrize_fragment_wannier_position_import(dc, center_bohr, owner_frag, num_wann, &
      num_bands, v_matrix, nsym, symops, aa_global)
    use eigen_subdiag_sub, only: eigen_zheev
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: num_wann, owner_frag(num_wann), num_bands, nsym
    real(8), intent(inout) :: center_bohr(3,num_wann)
    complex(8), intent(in) :: v_matrix(num_bands,num_wann)
    type(t_wannier_symop), intent(in) :: symops(:)
    complex(8), intent(inout) :: aa_global(3,num_wann,num_wann)
    integer :: ifrag, isym, nowned, iw, jw, ib, jb, npts, nspin_file, nstate_frag_file, nstate_tot_file
    integer :: nxyz_domain(3), n_basis_frag, axis, p, n_s2_gt_09, n_s2_gt_05
    integer, allocatable :: wann_index(:), pmap(:)
    real(8), allocatable :: phi_basis(:,:), coef_wf(:,:), psi_state(:,:), eval(:)
    complex(8), allocatable :: wannier_frag(:,:), srep(:,:), gram(:,:), eigvec(:,:), drep(:,:), invsqrt(:,:)
    complex(8), allocatable :: ac(:,:), work(:,:), asym(:,:)
    real(8) :: hvol, min_eval, max_eval, map_residual, polar_residual
    real(8) :: origin_global(3), change_norm2, base_norm2, herm_residual, eval_floor
    real(8) :: subspace_leakage, closure_tol
    real(8) :: sym_residual2, sym_norm2
    logical :: seed_ok, map_ok

    if(nsym <= 0) return
    hvol = dc%system_tot%hvol
    eval_floor = 1.0d-10
    closure_tol = 1.0d-3
    do ifrag=1,dc%n_frag
      nowned = count(owner_frag(1:num_wann) == ifrag)
      if(nowned <= 0) cycle
      allocate(wann_index(nowned))
      ib = 0
      do iw=1,num_wann
        if(owner_frag(iw) /= ifrag) cycle
        ib = ib + 1
        wann_index(ib) = iw
      end do

      call read_fragment_lcfo_seed_for_wannier_import(dc, ifrag, num_bands, &
        nxyz_domain, nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag, &
        phi_basis, coef_wf, seed_ok)
      if(.not. seed_ok) then
        write(*,'(1x,a,i0,a)') "[DC-LCFO-W90-SYM] fragment=", ifrag, &
          " position sym skipped: missing LCFO seed"
        deallocate(wann_index)
        cycle
      end if

      npts = product(nxyz_domain)
      allocate(psi_state(npts,num_bands), wannier_frag(npts,nowned))
      psi_state = matmul(phi_basis(1:npts,1:n_basis_frag), coef_wf(1:n_basis_frag,1:num_bands))
      wannier_frag = matmul(cmplx(psi_state(1:npts,1:num_bands), 0.0d0, kind=8), &
        v_matrix(1:num_bands,wann_index(1:nowned)))

      do isym=1,nsym
        if(symops(isym)%owner_frag /= ifrag) cycle
        if(trim(symops(isym)%label) /= 'inversion') cycle
        call fragment_symmetry_origin_global_import(dc, ifrag, symops(isym), origin_global)
        allocate(pmap(npts))
        call build_fragment_symmetry_grid_map(dc, ifrag, nxyz_domain, symops(isym), &
          pmap, map_ok, map_residual)
        if(.not. map_ok) then
          write(*,'(1x,a,i0,a,es12.5)') &
            "[DC-LCFO-W90-SYM] fragment=", ifrag, " position sym skipped: grid map residual=", map_residual
          deallocate(pmap)
          cycle
        end if

        allocate(srep(nowned,nowned), gram(nowned,nowned), eigvec(nowned,nowned), &
          drep(nowned,nowned), invsqrt(nowned,nowned), eval(nowned))
        srep = (0.0d0,0.0d0)
        do iw=1,nowned
          do ib=1,nowned
            do p=1,npts
              srep(iw,ib) = srep(iw,ib) + conjg(wannier_frag(p,iw)) * wannier_frag(pmap(p),ib) * hvol
            end do
          end do
        end do
        gram = matmul(conjg(transpose(srep)), srep)
        call eigen_zheev(gram, eval, eigvec)
        min_eval = minval(eval)
        max_eval = maxval(eval)
        subspace_leakage = sqrt(max(0.0d0, 1.0d0 - min(1.0d0, min_eval)))
        n_s2_gt_09 = count(eval(1:nowned) > 0.9d0)
        n_s2_gt_05 = count(eval(1:nowned) > 0.5d0)
        invsqrt = (0.0d0,0.0d0)
        drep = (0.0d0,0.0d0)
        do iw=1,nowned
          do ib=1,nowned
            invsqrt(iw,ib) = sum(eigvec(iw,1:nowned) * merge(1.0d0 / sqrt(eval(1:nowned)), 0.0d0, &
              eval(1:nowned) > eval_floor) * conjg(eigvec(ib,1:nowned)))
            drep(iw,ib) = sum(eigvec(iw,1:nowned) * merge(0.0d0, 1.0d0, &
              eval(1:nowned) > eval_floor) * conjg(eigvec(ib,1:nowned)))
          end do
        end do
        drep = matmul(srep, invsqrt) + drep
        polar_residual = hermitian_identity_residual(matmul(conjg(transpose(drep)), drep))
        if(subspace_leakage > closure_tol) then
          write(*,'(1x,a,i0,3(a,es12.5),2(a,i0))') &
            "[DC-LCFO-W90-SYM] fragment=", ifrag, &
            " position sym skipped: subspace_leakage=", subspace_leakage, &
            " min_s2=", min_eval, " closure_tol=", closure_tol, &
            " n_s2_gt_09=", n_s2_gt_09, " ncenter=", nowned
          deallocate(srep, gram, eigvec, drep, invsqrt, eval, pmap)
          cycle
        end if

        allocate(ac(nowned,nowned), work(nowned,nowned), asym(nowned,nowned))
        change_norm2 = 0.0d0
        base_norm2 = 0.0d0
        herm_residual = 0.0d0
        sym_residual2 = 0.0d0
        sym_norm2 = 0.0d0
        do axis=1,3
          ac = (0.0d0,0.0d0)
          do jb=1,nowned
            jw = wann_index(jb)
            do ib=1,nowned
              iw = wann_index(ib)
              ac(ib,jb) = aa_global(axis,iw,jw)
            end do
            ac(jb,jb) = ac(jb,jb) - cmplx(origin_global(axis), 0.0d0, kind=8)
          end do
          work = matmul(conjg(transpose(drep)), matmul(ac, drep))
          asym = 0.5d0 * (ac - work)
          asym = 0.5d0 * (asym + conjg(transpose(asym)))
          work = asym + matmul(conjg(transpose(drep)), matmul(asym, drep))
          sym_residual2 = sym_residual2 + sum(abs(work)**2)
          sym_norm2 = sym_norm2 + sum(abs(asym)**2)
          do jb=1,nowned
            asym(jb,jb) = asym(jb,jb) + cmplx(origin_global(axis), 0.0d0, kind=8)
          end do
          change_norm2 = change_norm2 + sum(abs(asym - aa_global(axis,wann_index(1:nowned),wann_index(1:nowned)))**2)
          base_norm2 = base_norm2 + sum(abs(aa_global(axis,wann_index(1:nowned),wann_index(1:nowned)))**2)
          herm_residual = max(herm_residual, maxval(abs(asym - conjg(transpose(asym)))))
          do jb=1,nowned
            jw = wann_index(jb)
            do ib=1,nowned
              iw = wann_index(ib)
              aa_global(axis,iw,jw) = asym(ib,jb)
            end do
          end do
        end do
        do ib=1,nowned
          iw = wann_index(ib)
          center_bohr(1:3,iw) = real(aa_global(1:3,iw,iw), kind=8)
        end do
        write(*,'(1x,a,i0,a,7(a,es12.5),3(a,i0))') &
          "[DC-LCFO-W90-SYM] fragment=", ifrag, " position sym label=inversion", &
          " min_s2=", min_eval, " max_s2=", max_eval, " polar_unit_res=", polar_residual, &
          " rel_change=", sqrt(change_norm2 / max(base_norm2, 1.0d-300)), &
          " sym_res=", sqrt(sym_residual2 / max(sym_norm2, 1.0d-300)), &
          " herm_res=", herm_residual, " map_res=", map_residual, &
          " n_s2_gt_09=", n_s2_gt_09, " n_s2_gt_05=", n_s2_gt_05, " ncenter=", nowned

        deallocate(ac, work, asym, srep, gram, eigvec, drep, invsqrt, eval, pmap)
      end do

      deallocate(wann_index, phi_basis, coef_wf, psi_state, wannier_frag)
    end do
  end subroutine symmetrize_fragment_wannier_position_import

  subroutine read_fragment_lcfo_seed_for_wannier_import(dc, ifrag, num_bands, &
      nxyz_domain, nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag, &
      phi_basis, coef_wf, ok, strict_current_run, failure_message)
    use filesystem, only: get_filehandle
    use structures, only: s_dcdft
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag, num_bands
    integer, intent(out) :: nxyz_domain(3), nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag
    real(8), allocatable, intent(out) :: phi_basis(:,:), coef_wf(:,:)
    logical, intent(out) :: ok
    logical, intent(in), optional :: strict_current_run
    character(*), intent(out), optional :: failure_message
    integer :: iunit, io, allocation_status, n_frag_file, ispin, ibasis, ix, iy, iz, p
    integer :: provenance_magic,provenance_version,provenance_fragments,provenance_fragment
    integer :: provenance_spins,provenance_fragment_states,provenance_total_states
    integer, allocatable :: n_mat_tmp(:), n_basis_tmp(:,:), index_basis_tmp(:,:,:), n_basis_core(:)
    integer, allocatable :: provenance_basis(:)
    real(8), allocatable :: coef_all(:,:,:), f_basis(:,:,:,:,:), esp_tmp(:,:)
    integer(int64) :: point_count64,basis_count64,flattened_count64
    character(64) :: provenance_token
    character(256) :: filename,provenance_filename
    character(512) :: message
    logical :: strict_mode,metadata_ok

    ok = .false.
    strict_mode=.false.
    if(present(strict_current_run)) strict_mode=strict_current_run
    if(present(failure_message)) failure_message=''
    message=''
    nxyz_domain = 0
    nspin_file = 0
    nstate_frag_file = 0
    nstate_tot_file = 0
    n_basis_frag = 0

    if(strict_mode) then
      write(filename,'(a,a,i6.6,a,a)') trim(import_run_root_dir()), &
        'data_dcdft/fragments/',ifrag,'/',binfile_wf_wannier_seed
      write(provenance_filename,'(a,a,i6.6,a)') trim(import_run_root_dir()), &
        'data_dcdft/fragments/',ifrag,'/wavefunctions_wannier_seed.provenance'
      iunit=get_filehandle()
      open(iunit,file=provenance_filename,status='old',action='read',iostat=io)
      if(io /= 0) then
        message='missing current-run SAWF seed provenance sidecar'
        goto 900
      end if
      read(iunit,*,iostat=io) provenance_magic,provenance_version
      if(io == 0) read(iunit,'(a)',iostat=io) provenance_token
      if(io == 0) read(iunit,*,iostat=io) provenance_fragments,provenance_fragment, &
        provenance_spins,provenance_fragment_states,provenance_total_states
      if(io /= 0) then
        close(iunit)
        message='malformed current-run SAWF seed provenance sidecar'
        goto 900
      end if
      if(provenance_spins <= 0) then
        close(iunit)
        message='invalid spin count in SAWF seed provenance sidecar'
        goto 900
      end if
      allocate(provenance_basis(provenance_spins),stat=allocation_status)
      if(allocation_status /= 0) then
        close(iunit)
        message='SAWF provenance basis allocation failed'
        goto 900
      end if
      read(iunit,*,iostat=io) provenance_basis
      close(iunit)
      if(io /= 0 .or. provenance_magic /= sawf_seed_provenance_magic .or. &
          provenance_version /= sawf_seed_provenance_version .or. &
          trim(provenance_token) /= trim(current_sawf_seed_token) .or. &
          provenance_fragments /= dc%n_frag .or. provenance_fragment /= ifrag) then
        message='stale or mixed current-run SAWF seed provenance'
        goto 900
      end if
    else
      filename = wannier_wavefunction_seed_filename_import(ifrag)
    end if
    iunit = get_filehandle()
    open(iunit,file=filename,form='unformatted',access='stream',status='old',iostat=io)
    if(io /= 0) then
      message='failed to open LCFO wavefunction seed'
      goto 900
    end if
    read(iunit,iostat=io) n_frag_file,nspin_file,nstate_frag_file,nstate_tot_file
    if(io /= 0) then
      close(iunit)
      message='failed to read LCFO wavefunction seed header'
      goto 900
    end if
    call validate_sawf_seed_header(dc%n_frag,ifrag,num_bands,n_frag_file,nspin_file, &
      nstate_frag_file,nstate_tot_file,metadata_ok,message)
    if(.not.metadata_ok) then
      close(iunit)
      goto 900
    end if
    allocate(n_mat_tmp(nspin_file), n_basis_tmp(n_frag_file,nspin_file), &
      index_basis_tmp(nstate_frag_file,n_frag_file,nspin_file), &
      coef_all(nstate_frag_file,nstate_tot_file,nspin_file),esp_tmp(nstate_tot_file,nspin_file), &
      stat=allocation_status)
    if(allocation_status /= 0) then
      close(iunit)
      message='LCFO wavefunction seed metadata allocation failed'
      goto 900
    end if
    read(iunit,iostat=io) n_mat_tmp
    if(io == 0) read(iunit,iostat=io) n_basis_tmp
    if(io == 0) read(iunit,iostat=io) index_basis_tmp
    if(io == 0) read(iunit,iostat=io) coef_all
    if(io == 0) read(iunit,iostat=io) esp_tmp
    close(iunit)
    if(io /= 0) then
      message='truncated LCFO wavefunction seed metadata or coefficients'
      goto 900
    end if
    call validate_sawf_seed_basis_metadata(nstate_frag_file,nstate_tot_file,n_mat_tmp, &
      n_basis_tmp,metadata_ok,message,index_basis_tmp)
    if(.not.metadata_ok) goto 900
    n_basis_frag = n_basis_tmp(ifrag,1)
    if(n_basis_frag < 1) then
      message='LCFO wavefunction seed has no basis for the requested fragment'
      goto 900
    end if
    if(strict_mode) then
      if(provenance_spins /= nspin_file .or. provenance_fragment_states /= nstate_frag_file .or. &
          provenance_total_states /= nstate_tot_file .or. any(provenance_basis /= n_basis_tmp(ifrag,:))) then
        message='current-run SAWF provenance and wavefunction metadata disagree'
        goto 900
      end if
    end if
    allocate(coef_wf(nstate_frag_file,num_bands),stat=allocation_status)
    if(allocation_status /= 0) then
      message='LCFO wavefunction coefficient allocation failed'
      goto 900
    end if
    coef_wf(1:nstate_frag_file,1:num_bands) = coef_all(1:nstate_frag_file,1:num_bands,1)

    write(filename, '(a,a,i6.6,a,a)') trim(import_run_root_dir()), 'data_dcdft/fragments/', ifrag, '/', binfile_bf
    iunit = get_filehandle()
    open(iunit,file=filename,form='unformatted',access='stream',status='old',iostat=io)
    if(io /= 0) then
      message='failed to open LCFO fragment basis file'
      goto 900
    end if
    read(iunit,iostat=io) nxyz_domain(1:3),ispin,ibasis
    if(io /= 0 .or. any(nxyz_domain <= 0)) then
      close(iunit)
      message='invalid LCFO fragment basis header'
      goto 900
    end if
    if(ispin /= nspin_file .or. ibasis /= nstate_frag_file) then
      close(iunit)
      message='LCFO basis and wavefunction seed dimensions disagree'
      goto 900
    end if
    point_count64=1_int64
    do ix=1,3
      if(point_count64 > huge(0_int64)/int(nxyz_domain(ix),int64)) then
        close(iunit); message='LCFO basis grid size overflows int64'; goto 900
      end if
      point_count64=point_count64*int(nxyz_domain(ix),int64)
    end do
    if(point_count64 > int(huge(0),int64)) then
      close(iunit); message='LCFO basis grid exceeds default integer indexing'; goto 900
    end if
    if(point_count64 > huge(0_int64)/int(nspin_file,int64)) then
      close(iunit); message='LCFO basis grid-spin size overflows int64'; goto 900
    end if
    basis_count64=point_count64*int(nspin_file,int64)
    if(basis_count64 > huge(0_int64)/int(nstate_frag_file,int64)) then
      close(iunit); message='LCFO basis allocation size overflows int64'; goto 900
    end if
    basis_count64=basis_count64*int(nstate_frag_file,int64)
    if(point_count64 > huge(0_int64)/int(nstate_frag_file,int64)) then
      close(iunit); message='LCFO flattened basis size overflows int64'; goto 900
    end if
    flattened_count64=point_count64*int(nstate_frag_file,int64)
    if(basis_count64 > int(huge(0),int64) .or. flattened_count64 > int(huge(0),int64)) then
      close(iunit); message='LCFO basis allocation exceeds default integer indexing'; goto 900
    end if
    allocate(n_basis_core(nspin_file),stat=allocation_status)
    if(allocation_status == 0) allocate(f_basis(nxyz_domain(1),nxyz_domain(2),nxyz_domain(3), &
      nspin_file,nstate_frag_file),stat=allocation_status)
    if(allocation_status /= 0) then
      close(iunit); message='LCFO fragment basis allocation failed'; goto 900
    end if
    read(iunit,iostat=io) n_basis_core
    if(io == 0) read(iunit,iostat=io) f_basis
    close(iunit)
    if(io /= 0) then
      message='truncated LCFO fragment basis file'
      goto 900
    end if
    if(any(n_basis_core /= n_basis_tmp(ifrag,:))) then
      message='LCFO basis counts disagree with wavefunction seed metadata'
      goto 900
    end if

    allocate(phi_basis(int(point_count64),nstate_frag_file),stat=allocation_status)
    if(allocation_status /= 0) then
      message='LCFO flattened fragment basis allocation failed'
      goto 900
    end if
    p = 0
    do iz=1,nxyz_domain(3)
      do iy=1,nxyz_domain(2)
        do ix=1,nxyz_domain(1)
          p = p + 1
          phi_basis(p,1:nstate_frag_file) = f_basis(ix,iy,iz,1,1:nstate_frag_file)
        end do
      end do
    end do
    deallocate(n_mat_tmp,n_basis_tmp,index_basis_tmp,coef_all,esp_tmp,n_basis_core,f_basis)
    if(allocated(provenance_basis)) deallocate(provenance_basis)
    ok = .true.
    if(strict_mode) write(*,'(1x,a,i0,2a)') '[DC-LCFO-SAWF-DMN] strict current-run seed fragment=', &
      ifrag,' token=',trim(current_sawf_seed_token)
    return

900 continue
    if(allocated(n_mat_tmp)) deallocate(n_mat_tmp)
    if(allocated(n_basis_tmp)) deallocate(n_basis_tmp)
    if(allocated(index_basis_tmp)) deallocate(index_basis_tmp)
    if(allocated(coef_all)) deallocate(coef_all)
    if(allocated(esp_tmp)) deallocate(esp_tmp)
    if(allocated(n_basis_core)) deallocate(n_basis_core)
    if(allocated(f_basis)) deallocate(f_basis)
    if(allocated(provenance_basis)) deallocate(provenance_basis)
    if(allocated(phi_basis)) deallocate(phi_basis)
    if(allocated(coef_wf)) deallocate(coef_wf)
    if(present(failure_message)) failure_message=trim(message)
  end subroutine read_fragment_lcfo_seed_for_wannier_import

  subroutine read_fragment_lcfo_buffer_seed_for_wannier_import(dc, ifrag, num_bands, &
      nxyz_domain, nxyz_buffer, nxyz_box, nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag, &
      phi_basis, coef_wf, ok)
    use filesystem, only: get_filehandle
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag, num_bands
    integer, intent(out) :: nxyz_domain(3), nxyz_buffer(3), nxyz_box(3)
    integer, intent(out) :: nspin_file, nstate_frag_file, nstate_tot_file, n_basis_frag
    real(8), allocatable, intent(out) :: phi_basis(:,:), coef_wf(:,:)
    logical, intent(out) :: ok
    integer :: iunit, io, n_frag_file, ispin_read, ibasis_read, ix, iy, iz, p
    integer :: magic_file, version_file
    integer, allocatable :: n_mat_tmp(:), n_basis_tmp(:,:), index_basis_tmp(:,:,:), n_basis_file(:)
    real(8), allocatable :: coef_all(:,:,:), esp_tmp(:,:), phi_tmp(:,:,:)
    character(256) :: filename

    ok = .false.
    nxyz_domain = 0
    nxyz_buffer = 0
    nxyz_box = 0
    nspin_file = 0
    nstate_frag_file = 0
    nstate_tot_file = 0
    n_basis_frag = 0

    filename = wannier_wavefunction_seed_filename_import(ifrag)
    iunit = get_filehandle()
    open(iunit,file=filename,form='unformatted',access='stream',status='old',iostat=io)
    if(io /= 0) return
    read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
    if(nspin_file < 1 .or. nstate_tot_file < num_bands) then
      close(iunit)
      return
    end if
    allocate(n_mat_tmp(nspin_file), n_basis_tmp(n_frag_file,nspin_file), &
      index_basis_tmp(nstate_frag_file,n_frag_file,nspin_file), &
      coef_all(nstate_frag_file,nstate_tot_file,nspin_file), esp_tmp(nstate_tot_file,nspin_file))
    read(iunit) n_mat_tmp(1:nspin_file)
    read(iunit) n_basis_tmp(1:n_frag_file,1:nspin_file)
    read(iunit) index_basis_tmp(1:nstate_frag_file,1:n_frag_file,1:nspin_file)
    read(iunit) coef_all(1:nstate_frag_file,1:nstate_tot_file,1:nspin_file)
    read(iunit, iostat=io) esp_tmp(1:nstate_tot_file,1:nspin_file)
    close(iunit)
    if(io /= 0) then
      deallocate(n_mat_tmp, n_basis_tmp, index_basis_tmp, coef_all, esp_tmp)
      return
    end if
    n_basis_frag = n_basis_tmp(ifrag,1)
    if(n_basis_frag < 1) then
      deallocate(n_mat_tmp, n_basis_tmp, index_basis_tmp, coef_all, esp_tmp)
      return
    end if
    allocate(coef_wf(nstate_frag_file,num_bands))
    coef_wf(1:nstate_frag_file,1:num_bands) = coef_all(1:nstate_frag_file,1:num_bands,1)
    deallocate(n_mat_tmp, n_basis_tmp, index_basis_tmp, coef_all, esp_tmp)

    write(filename, '(a,a,i6.6,a,a)') trim(import_run_root_dir()), 'data_dcdft/fragments/', ifrag, '/', binfile_bfb
    iunit = get_filehandle()
    open(iunit,file=filename,form='unformatted',access='stream',status='old',iostat=io)
    if(io /= 0) then
      deallocate(coef_wf)
      return
    end if
    read(iunit) magic_file, version_file
    if(magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) then
      close(iunit)
      deallocate(coef_wf)
      return
    end if
    read(iunit) nxyz_domain(1:3), nxyz_buffer(1:3), nspin_file, nstate_frag_file
    if(any(nxyz_domain(1:3) <= 0) .or. any(nxyz_buffer(1:3) < 0)) then
      close(iunit)
      deallocate(coef_wf)
      return
    end if
    nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer(1:3)
    allocate(n_basis_file(nspin_file))
    read(iunit) n_basis_file(1:nspin_file)
    if(n_basis_file(1) /= n_basis_frag) then
      close(iunit)
      deallocate(n_basis_file, coef_wf)
      return
    end if

    allocate(phi_basis(product(nxyz_box),nstate_frag_file))
    allocate(phi_tmp(nxyz_box(1),nxyz_box(2),nxyz_box(3)))
    phi_basis = 0.0d0
    do ispin_read=1,nspin_file
      do ibasis_read=1,nstate_frag_file
        read(iunit, iostat=io) phi_tmp(1:nxyz_box(1),1:nxyz_box(2),1:nxyz_box(3))
        if(io /= 0) exit
        if(ispin_read == 1) then
          p = 0
          do iz=1,nxyz_box(3)
            do iy=1,nxyz_box(2)
              do ix=1,nxyz_box(1)
                p = p + 1
                phi_basis(p,ibasis_read) = phi_tmp(ix,iy,iz)
              end do
            end do
          end do
        end if
      end do
      if(io /= 0) exit
    end do
    close(iunit)
    if(io /= 0) then
      deallocate(n_basis_file, phi_tmp, phi_basis, coef_wf)
      return
    end if

    deallocate(n_basis_file, phi_tmp)
    ok = .true.
  end subroutine read_fragment_lcfo_buffer_seed_for_wannier_import

  subroutine build_fragment_symmetry_grid_map(dc, ifrag, nxyz_domain, symop, pmap, ok, max_residual)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag, nxyz_domain(3)
    type(t_wannier_symop), intent(in) :: symop
    integer, intent(out) :: pmap(product(nxyz_domain))
    logical, intent(out) :: ok
    real(8), intent(out) :: max_residual
    integer :: ix, iy, iz, axis, p, nearest(3), mapped_index(3), rot_inv(3,3)
    real(8) :: cell_length(3), hgrid(3), pos(3), mapped(3), relative(3), scaled, residual
    real(8), parameter :: grid_map_tol = 1.0d-7

    ok = .true.
    max_residual = 0.0d0
    call fragment_cell_lengths_import(dc, ifrag, cell_length)
    hgrid(1:3) = cell_length(1:3) / dble(nxyz_domain(1:3))
    rot_inv = transpose(symop%rot)
    p = 0
    do iz=1,nxyz_domain(3)
      pos(3) = dble(iz - 1) * hgrid(3)
      do iy=1,nxyz_domain(2)
        pos(2) = dble(iy - 1) * hgrid(2)
        do ix=1,nxyz_domain(1)
          pos(1) = dble(ix - 1) * hgrid(1)
          p = p + 1
          relative(1:3) = pos(1:3) - symop%origin_bohr(1:3) - symop%tau_local(1:3)
          mapped(1:3) = symop%origin_bohr(1:3) + matmul_int_real(rot_inv, relative)
          call wrap_fragment_local_point(mapped, cell_length)
          do axis=1,3
            scaled = mapped(axis) / hgrid(axis)
            nearest(axis) = nint(scaled)
            residual = abs(scaled - dble(nearest(axis))) * hgrid(axis)
            max_residual = max(max_residual, residual)
            if(residual > grid_map_tol) ok = .false.
            mapped_index(axis) = modulo(nearest(axis), nxyz_domain(axis)) + 1
          end do
          pmap(p) = mapped_index(1) + (mapped_index(2) - 1) * nxyz_domain(1) &
            + (mapped_index(3) - 1) * nxyz_domain(1) * nxyz_domain(2)
        end do
      end do
    end do
  end subroutine build_fragment_symmetry_grid_map

  subroutine build_fragment_center_permutation_import(dc, ifrag, center_bohr, owner_frag, num_wann, &
      symop, wann_index, nowned, center_perm, ok, max_residual, rms_residual)
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag, num_wann, owner_frag(num_wann), nowned, wann_index(nowned)
    real(8), intent(in) :: center_bohr(3,num_wann)
    type(t_wannier_symop), intent(in) :: symop
    integer, intent(out) :: center_perm(nowned)
    logical, intent(out) :: ok
    real(8), intent(out) :: max_residual, rms_residual
    integer :: i, j, iw, jw, axis, jbest
    integer, allocatable :: seen(:)
    real(8) :: cell_length(3), center_i(3), center_j(3), relative(3), mapped(3), delta(3)
    real(8) :: dist2, best_dist2
    real(8), parameter :: center_match_tol = 1.0d-3

    center_perm = 0
    ok = .true.
    max_residual = 0.0d0
    rms_residual = 0.0d0
    if(nowned <= 0) return
    allocate(seen(nowned))
    seen = 0
    call fragment_cell_lengths_import(dc, ifrag, cell_length)

    do i=1,nowned
      iw = wann_index(i)
      if(owner_frag(iw) /= ifrag) cycle
      call point_fragment_local_position_import(dc, ifrag, center_bohr(1:3,iw), center_i)
      relative(1:3) = center_i(1:3) - symop%origin_bohr(1:3)
      mapped(1:3) = symop%origin_bohr(1:3) + matmul_int_real(symop%rot, relative) &
        + symop%tau_local(1:3)
      call wrap_fragment_local_point(mapped, cell_length)

      best_dist2 = huge(1d0)
      jbest = 0
      do j=1,nowned
        jw = wann_index(j)
        call point_fragment_local_position_import(dc, ifrag, center_bohr(1:3,jw), center_j)
        do axis=1,3
          delta(axis) = periodic_delta_import(center_j(axis) - mapped(axis), cell_length(axis))
        end do
        dist2 = local_distance2(delta)
        if(dist2 < best_dist2) then
          best_dist2 = dist2
          jbest = j
        end if
      end do
      center_perm(i) = jbest
      if(jbest < 1) then
        ok = .false.
      else
        seen(jbest) = seen(jbest) + 1
        max_residual = max(max_residual, sqrt(best_dist2))
        rms_residual = rms_residual + best_dist2
        if(sqrt(best_dist2) > center_match_tol) ok = .false.
      end if
    end do
    if(any(seen(1:nowned) /= 1)) ok = .false.
    rms_residual = sqrt(rms_residual / dble(max(1,nowned)))
    deallocate(seen)
  end subroutine build_fragment_center_permutation_import

  real(8) function permutation_target_residual(mat, perm) result(residual)
    implicit none
    complex(8), intent(in) :: mat(:,:)
    integer, intent(in) :: perm(:)
    integer :: i, j, n
    complex(8) :: target

    n = size(mat,1)
    residual = 0.0d0
    do j=1,n
      do i=1,n
        target = (0.0d0,0.0d0)
        if(perm(j) == i) target = (1.0d0,0.0d0)
        residual = residual + abs(mat(i,j) - target)**2
      end do
    end do
    residual = sqrt(residual / dble(max(1,n*n)))
  end function permutation_target_residual

  real(8) function hermitian_identity_residual(mat) result(residual)
    implicit none
    complex(8), intent(in) :: mat(:,:)
    integer :: i, j, n
    complex(8) :: target

    n = size(mat,1)
    residual = 0.0d0
    do j=1,n
      do i=1,n
        target = (0.0d0,0.0d0)
        if(i == j) target = (1.0d0,0.0d0)
        residual = residual + abs(mat(i,j) - target)**2
      end do
    end do
    residual = sqrt(residual / dble(max(1,n*n)))
  end function hermitian_identity_residual

  function import_run_root_dir() result(root)
    use salmon_global, only: base_directory
    implicit none
    character(256) :: root
    integer :: ipos

    root = trim(base_directory)
    ipos = index(root, 'data_dcdft/fragments/')
    if(ipos <= 0) ipos = index(root, 'data_dcdft/total/')
    if(ipos > 0) then
      if(ipos == 1) then
        root = './'
      else
        root = root(1:ipos-1)
      end if
    end if
    if(len_trim(root) == 0) root = './'
  end function import_run_root_dir


  subroutine dc_lcfo_flux(lg,mg,system,info,stencil,ppg,energy,rho_s,v_local,spsi,shpsi,sttpsi,srg,dc)
    use communication, only: comm_summation
    use salmon_global, only: yn_dc_lcfo_diag
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),      intent(in) :: stencil
    type(s_pp_grid),      intent(in) :: ppg
    type(s_dft_energy),   intent(in) :: energy
    type(s_scalar),       intent(in) :: rho_s(system%nspin)
    type(s_scalar),       intent(in) :: V_local(system%nspin)
    type(s_orbital),      intent(in) :: spsi
    type(s_orbital)                  :: shpsi,sttpsi
    type(s_sendrecv_grid)            :: srg
    type(s_dcdft)                    :: dc
    !
    type halo_info
      integer :: id_src,id_dst,ifrag_src,dvec(3),length(3),dsp_send(3),dsp_recv(3),axis
      real(8),allocatable :: mat_H_local(:,:,:)
      real(8),allocatable :: mat_H_surface_cross(:,:,:)
      real(8),allocatable :: mat_V_local(:,:,:,:)
      real(8),allocatable :: mat_Xi_flux_local(:,:,:,:)
      real(8),allocatable :: trace_send(:,:,:),trace_recv(:,:,:)
    end type halo_info
    !
    type(halo_info) :: halo(26) ! 26 = 3^3-1
    integer :: nspin,n_halo
    integer :: stencil_radius(3)
    integer :: id_array(dc%n_frag)
    integer :: n_basis(dc%n_frag,system%nspin), n_mat(system%nspin)
    integer :: index_basis(dc%nstate_frag,dc%n_frag,system%nspin)
    real(8) :: hvol
    real(8),allocatable :: f_basis(:,:,:,:,:),hf(:,:,:,:,:),wrk_array(:,:,:,:,:) &
    & ,esp_tot(:,:),mat_H_local(:,:,:),mat_H_volume_local(:,:,:),mat_H_volume_weak_local(:,:,:) &
    & ,mat_H_weak_kinetic(:,:,:),mat_H_weak_potential(:,:,:),mat_H_weak_nonlocal(:,:,:) &
    & ,mat_H_surface_self(:,:,:) &
    & ,mat_V_local(:,:,:,:),coef_wf(:,:,:),basis_transform(:,:,:)
    real(8),allocatable :: sawf_explicit_buffer(:,:,:,:,:)
    logical :: sawf_explicit_basis_active
    !
    integer :: i,j,n,ix,iy,iz,io,jo,ispin,ifrag,jfrag,i_halo

    if(dc%id_tot==0) write(*,*) "start DC-LCFO-Flux"
    hvol = system%hvol
    nspin = system%nspin
    sawf_explicit_basis_active=.false.
    call init_lcfo
    call initialize_sawf_seed_provenance()
    call calc_basis
    call hpsi_basis
    if(dc%id_tot==0) write(*,*) "basis functions operation: done"

    call calc_hamiltonian_matrix
    if(dc%id_tot==0) write(*,*) "Hamiltonian matrix: done"

    if(yn_dc_lcfo_diag=='y') then
      allocate(esp_tot(maxval(n_mat),nspin))
      if(dc%id_frag==0) allocate(coef_wf(dc%nstate_frag,dc%nstate_tot,nspin))
      if(allocated(coef_wf)) coef_wf = 0d0
#ifdef USE_EIGENEXA
      call diag_eigenexa
#else
      stop "DC-LCFO-Flux requires EigenExa; dense LAPACK fallback is disabled."
#endif
      if(dc%id_tot==0) write(*,*) "diagonalization: done"
    end if

    call output

    if(allocated(coef_wf)) deallocate(coef_wf)
    if(allocated(f_basis)) deallocate(f_basis)
    if(allocated(basis_transform)) deallocate(basis_transform)
    if(allocated(sawf_explicit_buffer)) deallocate(sawf_explicit_buffer)
    if(allocated(esp_tot)) deallocate(esp_tot)
    if(allocated(mat_H_local)) deallocate(mat_H_local)
    if(allocated(mat_H_volume_local)) deallocate(mat_H_volume_local)
    if(allocated(mat_H_volume_weak_local)) deallocate(mat_H_volume_weak_local)
    if(allocated(mat_H_weak_kinetic)) deallocate(mat_H_weak_kinetic)
    if(allocated(mat_H_weak_potential)) deallocate(mat_H_weak_potential)
    if(allocated(mat_H_weak_nonlocal)) deallocate(mat_H_weak_nonlocal)
    if(allocated(mat_H_surface_self)) deallocate(mat_H_surface_self)
      if(allocated(mat_V_local)) deallocate(mat_V_local)
      do i=1,n_halo
        if(allocated(halo(i)%mat_H_local)) deallocate(halo(i)%mat_H_local)
        if(allocated(halo(i)%mat_H_surface_cross)) deallocate(halo(i)%mat_H_surface_cross)
        if(allocated(halo(i)%mat_V_local)) deallocate(halo(i)%mat_V_local)
        if(allocated(halo(i)%mat_Xi_flux_local)) deallocate(halo(i)%mat_Xi_flux_local)
      end do
    if(dc%id_tot==0) write(*,*) "end DC-LCFO-Flux"

  contains

    subroutine initialize_sawf_seed_provenance()
      use communication, only: comm_bcast
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      integer(int64) :: clock_count
      integer :: values(8)

      current_sawf_seed_token=''
      if(dc%id_tot == 0) then
        call system_clock(count=clock_count)
        call date_and_time(values=values)
        write(current_sawf_seed_token,'(7(i0,a),i0)') values(1),'-',values(2),'-',values(3),'-', &
          values(5),'-',values(6),'-',values(7),'-',values(8),'-',clock_count
      end if
      call comm_bcast(current_sawf_seed_token,dc%icomm_tot,0)
    end subroutine initialize_sawf_seed_provenance

    subroutine init_lcfo
      use salmon_global, only: num_fragment
      implicit none
      integer :: lx,ly,lz,nonzero_dirs
      integer :: nxyz_domain(3)
      integer,dimension(3) :: nh,ir1,ir2,d
      integer :: id_tmp(dc%n_frag)

      id_tmp = 0
      if(.not. stencil%if_orthogonal) stop "DC-LCFO-Flux: nonorthogonal lattice is not supported"
      if(dc%optimized_fragment_geometry) stop "DC-LCFO-Flux: optimized fragment geometry is not supported"
      if(mod(dc%isize_tot,dc%n_frag) /= 0) stop "DC-LCFO-Flux: MPI size must be divisible by number of fragments"
      if(dc%id_frag==0) id_tmp(dc%i_frag) = dc%id_tot + 1
      call comm_summation(id_tmp,id_array,dc%n_frag,dc%icomm_tot)
      id_array = id_array - 1
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      do n=1,3
        stencil_radius(n) = active_laplacian_radius(n)
      end do

      nh = 0
      do n=1,3 ! x,y,z
        if(dc%nxyz_buffer(n) > nxyz_domain(n)) stop "DC-LCFO: buffer > domain"
        if(num_fragment(n) > 1 .and. dc%nxyz_buffer(n) < stencil_radius(n)) &
        & stop "DC-LCFO-Flux: buffer is smaller than the active stencil radius"
        if(num_fragment(n) > 1) nh(n) = 1
      end do

      i = 0
      do lx=-nh(1),nh(1)
      do ly=-nh(2),nh(2)
      do lz=-nh(3),nh(3)
        if(lx==0 .and. ly==0 .and. lz==0) cycle
        nonzero_dirs = count([lx,ly,lz] /= 0)
        if(nonzero_dirs /= 1) cycle
        i = i + 1
        halo(i)%dvec(1:3) = [lx, ly, lz]
        halo(i)%axis = 0
        do n=1,3
          if(halo(i)%dvec(n) /= 0) halo(i)%axis = n
        end do
        halo(i)%id_dst = -1
        halo(i)%id_src = -1
        do ifrag=1,dc%n_frag
        ! dc%ixyz_frag: r-grid index of the fragment origin
          ir1(1:3) = dc%ixyz_frag(1:3,ifrag) ! position of fragment ifrag
        ! dst neighbor (+)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) + halo(i)%dvec(1:3)*nxyz_domain(1:3) ! neighbor fragment
          d(1:3) = mod( ir1(1:3) - ir2(1:3) , dc%lg_tot%num(1:3) )
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_dst < 0) then
            halo(i)%id_dst = id_array(ifrag) ! process ID of the communication destination
          end if
        ! src neighbor (-)
          ir2(1:3) = dc%ixyz_frag(1:3,dc%i_frag) - halo(i)%dvec(1:3)*nxyz_domain(1:3) ! neighbor fragment
          d(1:3) = mod( ir1(1:3) - ir2(1:3) , dc%lg_tot%num(1:3) )
          if(d(1)==0 .and. d(2)==0 .and. d(3)==0 .and. halo(i)%id_src < 0) then
            halo(i)%id_src = id_array(ifrag) ! process ID of the communication source
            halo(i)%ifrag_src = ifrag
          end if
        end do ! ifrag
        if(halo(i)%id_dst < 0 .or. halo(i)%id_src < 0) stop "DC-LCFO: dst, src"
      end do
      end do
      end do
      n_halo = i ! # of the halo regions (neighbor fragments)

      do i=1,n_halo
        do n=1,3 ! x,y,z
          select case (halo(i)%dvec(n))
          case(0)
            halo(i)%length(n) = nxyz_domain(n)
            halo(i)%dsp_send(n) = 0
            halo(i)%dsp_recv(n) = 0
          case(1)
            halo(i)%length(n) = dc%nxyz_buffer(n)
            halo(i)%dsp_send(n) = nxyz_domain(n) - dc%nxyz_buffer(n)
            halo(i)%dsp_recv(n) = nxyz_domain(n) + dc%nxyz_buffer(n)
          case(-1)
            halo(i)%length(n) = dc%nxyz_buffer(n)
            halo(i)%dsp_send(n) = 0
            halo(i)%dsp_recv(n) = nxyz_domain(n)
          end select
        end do
      end do

    end subroutine init_lcfo

    subroutine calc_basis
      use eigen_subdiag_sub, only: eigen_dsyev
      use salmon_global, only: energy_cut,lambda_cut
      implicit none
      integer, parameter :: n_spectrum_thr = 7
      real(8), parameter :: spectrum_thr(n_spectrum_thr) = &
        [1.0d-3, 1.0d-6, 1.0d-8, 1.0d-10, 1.0d-12, 1.0d-14, 1.0d-16]
      integer :: nb(nspin),itmp(dc%n_frag,nspin)
      integer :: nxyz_domain(3),ix1,ix2,iy1,iy2,iz1,iz2
      integer :: ithr, count_by_frag(n_spectrum_thr,dc%n_frag,nspin)
      integer :: count_all(n_spectrum_thr,dc%n_frag,nspin)
      real(8),dimension(dc%nstate_frag,dc%nstate_frag,system%nspin) :: mat_S,mat_U
      real(8),dimension(dc%nstate_frag,system%nspin) :: lambda
      real(8) :: alpha_gs, norm_basis
      logical :: active_state

      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      ix1 = max(mg%is(1),1); ix2 = min(mg%ie(1),nxyz_domain(1))
      iy1 = max(mg%is(2),1); iy2 = min(mg%ie(2),nxyz_domain(2))
      iz1 = max(mg%is(3),1); iz2 = min(mg%ie(3),nxyz_domain(3))
      if(dc%id_tot==0 .and. energy_cut <= 0d0) then
        write(*,'(1x,a)') &
          "[DC-LCFO-FLUX-BASIS] warning: energy_cut<=0 keeps only local states below mu; RT f-sum may be incomplete"
      end if

      allocate(f_basis  (nxyz_domain(1),nxyz_domain(2),nxyz_domain(3),nspin,dc%nstate_frag))
      allocate(wrk_array(nxyz_domain(1),nxyz_domain(2),nxyz_domain(3),nspin,dc%nstate_frag))
      if(.not. allocated(basis_transform)) allocate(basis_transform(dc%nstate_frag,dc%nstate_frag,nspin))
      basis_transform = 0d0

    ! f_basis <-- | \bar{\phi} > (projected fragment orbitals)
      wrk_array = 0d0
!$omp parallel do collapse(2) private(io,ispin,iz,iy,ix,active_state) schedule(static)
      do io=info%io_s,info%io_e
      do ispin=1,nspin
        active_state = energy%esp(io,1,ispin) - system%mu < energy_cut
        if(active_state) then
      do iz=iz1,iz2
      do iy=iy1,iy2
      do ix=ix1,ix2
          wrk_array(ix,iy,iz,ispin,io) = spsi%rwf(ix,iy,iz,ispin,io,1,1) ! | \phi > @ core domain
      end do
      end do
      end do
        end if
      end do
      end do
!$omp end parallel do
      call comm_summation(wrk_array,f_basis,product(nxyz_domain)*nspin*dc%nstate_frag,info%icomm_rko)

    ! mat_S <-- S_{ij} = < \bar{\phi}_i | \bar{\phi}_j > (overlap matrix)
!$omp parallel do collapse(3) private(ispin,io,jo) schedule(static)
      do ispin=1,nspin
      do io=1,dc%nstate_frag
      do jo=1,dc%nstate_frag
        mat_S(io,jo,ispin) = sum(f_basis(:,:,:,ispin,io)*f_basis(:,:,:,ispin,jo)) * hvol
      end do
      end do
      end do
!$omp end parallel do

    ! diagonalize mat_S
      do ispin=1,nspin
        call eigen_dsyev(mat_S(:,:,ispin),lambda(:,ispin),mat_U(:,:,ispin))
      end do
      count_by_frag = 0
      if(dc%id_frag==0) then
        do ispin=1,nspin
          do ithr=1,n_spectrum_thr
            count_by_frag(ithr,dc%i_frag,ispin) = count(lambda(:,ispin) > spectrum_thr(ithr))
          end do
        end do
      end if
      call comm_summation(count_by_frag,count_all,n_spectrum_thr*dc%n_frag*nspin,dc%icomm_tot)
      if(dc%id_tot==0) then
        do ispin=1,nspin
          write(*,'(1x,a,i0,7(a,1pe9.1,a,i0,a,i0))') &
            "[DC-LCFO-FLUX-BASIS-SPECTRUM] ispin=", ispin, &
            " thr=", spectrum_thr(1), " min/max=", minval(count_all(1,:,ispin)), "/", maxval(count_all(1,:,ispin)), &
            " thr=", spectrum_thr(2), " min/max=", minval(count_all(2,:,ispin)), "/", maxval(count_all(2,:,ispin)), &
            " thr=", spectrum_thr(3), " min/max=", minval(count_all(3,:,ispin)), "/", maxval(count_all(3,:,ispin)), &
            " thr=", spectrum_thr(4), " min/max=", minval(count_all(4,:,ispin)), "/", maxval(count_all(4,:,ispin)), &
            " thr=", spectrum_thr(5), " min/max=", minval(count_all(5,:,ispin)), "/", maxval(count_all(5,:,ispin)), &
            " thr=", spectrum_thr(6), " min/max=", minval(count_all(6,:,ispin)), "/", maxval(count_all(6,:,ispin)), &
            " thr=", spectrum_thr(7), " min/max=", minval(count_all(7,:,ispin)), "/", maxval(count_all(7,:,ispin))
        end do
      end if

    ! f_basis <-- | lambda > (basis functions)
      wrk_array = f_basis
      f_basis = 0d0
      do ispin=1,nspin
        i = 0 ! count # of basis functions
        do io=dc%nstate_frag,1,-1
          if( lambda(io,ispin) > lambda_cut ) then ! cutoff for the eigenvalues of the overlap matrix
            i = i + 1 ! count # of basis functions
            do jo=1,dc%nstate_frag
              f_basis(:,:,:,ispin,i) = f_basis(:,:,:,ispin,i) &
              & + wrk_array(:,:,:,ispin,jo) * mat_U(jo,io,ispin) / sqrt(lambda(io,ispin))
              basis_transform(jo,i,ispin) = mat_U(jo,io,ispin) / sqrt(lambda(io,ispin))
            end do
          end if
        end do ! io
        nb(ispin) = i ! # of basis functions
      end do ! ispin

    ! Gram–Schmidt orthonormalization
      wrk_array = f_basis
      do ispin=1,nspin
        do io=1,nb(ispin)
          do jo=1,io-1
            alpha_gs = (sum(f_basis(:,:,:,ispin,jo)*wrk_array(:,:,:,ispin,io)) * hvol) &
            & / (sum(f_basis(:,:,:,ispin,jo)*f_basis(:,:,:,ispin,jo)) * hvol)
            wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
            & - f_basis(:,:,:,ispin,jo) * alpha_gs
            basis_transform(:,io,ispin) = basis_transform(:,io,ispin) &
            & - basis_transform(:,jo,ispin) * alpha_gs
          end do
          norm_basis = sqrt( sum(wrk_array(:,:,:,ispin,io)*wrk_array(:,:,:,ispin,io)) * hvol )
          if(norm_basis <= 1.0d-12) stop "DC-LCFO-Flux: null basis after core-S cleanup"
          basis_transform(:,io,ispin) = basis_transform(:,io,ispin) &
          & / norm_basis
          wrk_array(:,:,:,ispin,io) = wrk_array(:,:,:,ispin,io) &
          & / norm_basis
        end do
      end do ! ispin
      f_basis = wrk_array

    ! sttpsi <-- f_basis == | lambda > (basis functions)
      sttpsi%rwf = 0d0
!$omp parallel do collapse(2) private(io,ispin,iz,iy,ix) schedule(static)
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=iz1,iz2
      do iy=iy1,iy2
      do ix=ix1,ix2
        sttpsi%rwf(ix,iy,iz,ispin,io,1,1) = f_basis(ix,iy,iz,ispin,io)
      end do
      end do
      end do
      end do
      end do
!$omp end parallel do

    ! n_basis: # of basis functions
      itmp = 0
      if(dc%id_frag==0) itmp(dc%i_frag,1:nspin) = nb(1:nspin)
      call comm_summation(itmp,n_basis,dc%n_frag*nspin,dc%icomm_tot)
      index_basis = 0
      do ispin=1,nspin
        i = 0
        do ifrag=1,dc%n_frag
          do io=1,n_basis(ifrag,ispin)
            i = i + 1
            index_basis(io,ifrag,ispin) = i ! index_basis: index for the total matrix
          end do
        end do
        n_mat(ispin) = i ! n_mat: dimension of the total matrix
      end do

      deallocate(wrk_array)

    end subroutine calc_basis

    subroutine hpsi_basis
      use hamiltonian, only: hpsi
      implicit none
      integer :: ibx,iby,ibz,sx,sy,sz,nxyz_box(3)

      allocate(hf       (lg%num(1),lg%num(2),lg%num(3),nspin,dc%nstate_frag))
      allocate(wrk_array(lg%num(1),lg%num(2),lg%num(3),nspin,dc%nstate_frag))

      if(sawf_explicit_basis_active) then
        sttpsi%rwf=0.0d0
        nxyz_box=dc%nxyz_domain+2*dc%nxyz_buffer
        do io=info%io_s,info%io_e
          if(io>n_basis(dc%i_frag,1)) cycle
          do ispin=1,nspin
            do ibz=1,nxyz_box(3)
              sz=dc_buffer_box_to_local_index(ibz,dc%nxyz_domain(3),dc%nxyz_buffer(3))
              if(sz<lbound(sttpsi%rwf,3) .or. sz>ubound(sttpsi%rwf,3)) cycle
              do iby=1,nxyz_box(2)
                sy=dc_buffer_box_to_local_index(iby,dc%nxyz_domain(2),dc%nxyz_buffer(2))
                if(sy<lbound(sttpsi%rwf,2) .or. sy>ubound(sttpsi%rwf,2)) cycle
                do ibx=1,nxyz_box(1)
                  sx=dc_buffer_box_to_local_index(ibx,dc%nxyz_domain(1),dc%nxyz_buffer(1))
                  if(sx<lbound(sttpsi%rwf,1) .or. sx>ubound(sttpsi%rwf,1)) cycle
                  sttpsi%rwf(sx,sy,sz,ispin,io,1,1)= &
                    sawf_explicit_buffer(ibx,iby,ibz,ispin,io)
                end do
              end do
            end do
          end do
        end do
      end if
    ! shpsi <-- H | lambda > (Hamiltonian operation)
      call hpsi(sttpsi,shpsi,info,mg,v_local,system,stencil,srg,ppg)

    ! hf <-- shpsi == H | lambda >
      wrk_array = 0d0
!$omp parallel do collapse(2) private(io,ispin,iz,iy,ix) schedule(static)
      do io=info%io_s,info%io_e
      do ispin=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        wrk_array(ix,iy,iz,ispin,io) = shpsi%rwf(ix,iy,iz,ispin,io,1,1)
      end do
      end do
      end do
      end do
      end do
!$omp end parallel do
      call comm_summation(wrk_array,hf,product(lg%num)*nspin*dc%nstate_frag,info%icomm_rko)

      deallocate(wrk_array)

    end subroutine hpsi_basis

    subroutine exchange_surface_trace_halo()
      use communication, only: comm_isend, comm_irecv, comm_wait_all
      implicit none
      integer :: axis, side, send_side, face_pt, npts
      integer :: ix, iy, iz, io, ispin
      integer :: itag_send, itag_recv, itag_dir
      integer, dimension(n_halo) :: ireq_send, ireq_recv

      if(dc%id_frag /= 0) return

      do i_halo=1,n_halo
        axis = halo(i_halo)%axis
        send_side = halo(i_halo)%dvec(axis)
        ! The receiver uses dnu_r with the receiver-side outward normal.
        ! Convert the local face derivative before sending so the surface
        ! formula is identical on both endpoints of the face.
        side = -send_side
        npts = face_point_count(axis)
        if(allocated(halo(i_halo)%trace_send)) deallocate(halo(i_halo)%trace_send)
        if(allocated(halo(i_halo)%trace_recv)) deallocate(halo(i_halo)%trace_recv)
        allocate(halo(i_halo)%trace_send(npts,nspin,2*dc%nstate_frag))
        allocate(halo(i_halo)%trace_recv(npts,nspin,2*dc%nstate_frag))
        halo(i_halo)%trace_send = 0d0
        halo(i_halo)%trace_recv = 0d0

        face_pt = 0
        do iz=1,dc%nxyz_domain(3)
        do iy=1,dc%nxyz_domain(2)
        do ix=1,dc%nxyz_domain(1)
          if(face_axis_index([ix,iy,iz], axis) /= face_coord(axis, send_side)) cycle
          face_pt = face_pt + 1
          do ispin=1,nspin
            do io=1,n_basis(dc%i_frag,ispin)
              halo(i_halo)%trace_send(face_pt,ispin,io) = &
              & local_basis_value(ix,iy,iz,ispin,io)
              halo(i_halo)%trace_send(face_pt,ispin,dc%nstate_frag+io) = &
              & local_basis_dn(ix,iy,iz,ispin,io,axis,side)
            end do
          end do
        end do
        end do
        end do

        itag_dir = halo_tag_offset(halo(i_halo)%dvec)
        itag_send = dc%i_frag + itag_dir*dc%n_frag
        ireq_send(i_halo) = comm_isend(halo(i_halo)%trace_send,halo(i_halo)%id_dst,itag_send,dc%icomm_tot)
        itag_recv = halo(i_halo)%ifrag_src + itag_dir*dc%n_frag
        ireq_recv(i_halo) = comm_irecv(halo(i_halo)%trace_recv,halo(i_halo)%id_src,itag_recv,dc%icomm_tot)
      end do
      call comm_wait_all(ireq_recv)
      call comm_wait_all(ireq_send)

    end subroutine exchange_surface_trace_halo

    subroutine release_surface_trace_halo()
      implicit none

      if(dc%id_frag /= 0) return
      do i_halo=1,n_halo
        if(allocated(halo(i_halo)%trace_send)) deallocate(halo(i_halo)%trace_send)
        if(allocated(halo(i_halo)%trace_recv)) deallocate(halo(i_halo)%trace_recv)
      end do

    end subroutine release_surface_trace_halo

    integer function halo_tag_offset(dvec) result(offset)
      implicit none
      integer,intent(in) :: dvec(3)

      offset = (dvec(1) + 1)*9 + (dvec(2) + 1)*3 + (dvec(3) + 1)
    end function halo_tag_offset

    subroutine calc_hamiltonian_matrix
      implicit none
      integer :: axis, side, face_pt, ix_face, iy_face, iz_face
      integer :: l(3), idir
      real(8), parameter :: surface_penalty_factor = 10.0d0
      real(8) :: area_weight, alpha, u_l, v_l, dnu_l, dnv_l, u_r, dnu_r
      real(8) :: term_sum, term_local, term_nonlocal, term_face, term_vsum, term_vface, pavg
      real(8) :: xi_rel(3), term_xisum(3), term_xiface
      real(8), allocatable :: trace_local(:,:,:)
      real(8), allocatable :: grad_work(:,:,:)
      real(8), allocatable :: basis_grad_all(:,:,:,:,:,:)
      real(8), allocatable :: basis_nonlocal_core(:,:,:,:,:)

    ! diagonal block < lambda_{ifrag,io} | H | lambda_{ifrag,jo} >
      allocate(mat_H_local(dc%nstate_frag,dc%nstate_frag,nspin))
      allocate(mat_H_volume_local(dc%nstate_frag,dc%nstate_frag,nspin))
      allocate(mat_H_volume_weak_local(dc%nstate_frag,dc%nstate_frag,nspin))
      allocate(mat_H_weak_kinetic(dc%nstate_frag,dc%nstate_frag,nspin))
      allocate(mat_H_weak_potential(dc%nstate_frag,dc%nstate_frag,nspin))
      allocate(mat_H_weak_nonlocal(dc%nstate_frag,dc%nstate_frag,nspin))
      allocate(mat_H_surface_self(dc%nstate_frag,dc%nstate_frag,nspin))
      allocate(mat_V_local(3,dc%nstate_frag,dc%nstate_frag,nspin))
      mat_H_local = 0d0
      mat_H_volume_local = 0d0
      mat_H_volume_weak_local = 0d0
      mat_H_weak_kinetic = 0d0
      mat_H_weak_potential = 0d0
      mat_H_weak_nonlocal = 0d0
      mat_H_surface_self = 0d0
      mat_V_local = 0d0
      l = dc%nxyz_domain
      allocate(basis_grad_all(l(1),l(2),l(3),nspin,dc%nstate_frag,3))
      allocate(basis_nonlocal_core(l(1),l(2),l(3),nspin,dc%nstate_frag))
      basis_grad_all = 0d0
      call build_fragment_nonlocal_basis_action(basis_nonlocal_core)
      do ispin=1,nspin
      do io=1,n_basis(dc%i_frag,ispin)
!$omp parallel do collapse(3) private(ix,iy,iz,idir) schedule(static)
        do iz=1,l(3)
        do iy=1,l(2)
        do ix=1,l(1)
          do idir=1,3
            basis_grad_all(ix,iy,iz,ispin,io,idir) = local_basis_grad(ix,iy,iz,ispin,io,idir)
          end do
        end do
        end do
        end do
!$omp end parallel do
      end do
      end do
      do ispin=1,nspin
!$omp parallel do private(io,jo,idir,term_sum,term_local,term_nonlocal) schedule(static)
      do io=1,n_basis(dc%i_frag,ispin)
      do jo=1,n_basis(dc%i_frag,ispin)
        mat_H_volume_local(io,jo,ispin) = &
        & + sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)*hf(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
        term_sum = 0d0
        do idir=1,3
          term_sum = term_sum &
          & + 0.5d0 * sum(basis_grad_all(1:l(1),1:l(2),1:l(3),ispin,io,idir) &
          &              * basis_grad_all(1:l(1),1:l(2),1:l(3),ispin,jo,idir)) * hvol
        end do
        term_local=sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io) &
        &       * V_local(ispin)%f(1:l(1),1:l(2),1:l(3)) &
        &       * f_basis(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
        term_nonlocal=sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io) &
        &       * basis_nonlocal_core(1:l(1),1:l(2),1:l(3),ispin,jo)) * hvol
        mat_H_weak_kinetic(io,jo,ispin)=term_sum
        mat_H_weak_potential(io,jo,ispin)=term_local
        mat_H_weak_nonlocal(io,jo,ispin)=term_nonlocal
        mat_H_volume_weak_local(io,jo,ispin)=term_sum+term_local+term_nonlocal
        if(use_weak_volume_hamiltonian_mode()) then
          mat_H_local(io,jo,ispin) = mat_H_volume_weak_local(io,jo,ispin)
        else
          mat_H_local(io,jo,ispin) = mat_H_volume_local(io,jo,ispin)
        end if
      end do
      end do
!$omp end parallel do
      end do
      allocate(grad_work(l(1),l(2),l(3)))
      do ispin=1,nspin
      do idir=1,3
      do jo=1,n_basis(dc%i_frag,ispin)
!$omp parallel do collapse(3) private(ix,iy,iz) schedule(static)
        do iz=1,l(3)
        do iy=1,l(2)
        do ix=1,l(1)
          grad_work(ix,iy,iz) = basis_grad_all(ix,iy,iz,ispin,jo,idir)
        end do
        end do
        end do
!$omp end parallel do
!$omp parallel do private(io) schedule(static)
        do io=1,n_basis(dc%i_frag,ispin)
          mat_V_local(idir,io,jo,ispin) = &
          & sum(f_basis(1:l(1),1:l(2),1:l(3),ispin,io)*grad_work(1:l(1),1:l(2),1:l(3))) * hvol
        end do
!$omp end parallel do
      end do
      end do
      end do
      deallocate(grad_work)

    ! DG surface jump/average/penalty terms on face-neighbor blocks.
      if(dc%id_frag==0) then
        call exchange_surface_trace_halo()
        do i_halo=1,n_halo
          axis = halo(i_halo)%axis
          side = -halo(i_halo)%dvec(axis)
          jfrag = halo(i_halo)%ifrag_src
          allocate(halo(i_halo)%mat_H_local(dc%nstate_frag,dc%nstate_frag,nspin))
          allocate(halo(i_halo)%mat_H_surface_cross(dc%nstate_frag,dc%nstate_frag,nspin))
          allocate(halo(i_halo)%mat_V_local(3,dc%nstate_frag,dc%nstate_frag,nspin))
          allocate(halo(i_halo)%mat_Xi_flux_local(3,dc%nstate_frag,dc%nstate_frag,nspin))
          halo(i_halo)%mat_H_local = 0d0
          halo(i_halo)%mat_H_surface_cross = 0d0
          halo(i_halo)%mat_V_local = 0d0
          halo(i_halo)%mat_Xi_flux_local = 0d0
          area_weight = system%hvol / system%hgs(axis)
          alpha = surface_penalty_factor / system%hgs(axis)
          allocate(trace_local(face_point_count(axis),nspin,2*dc%nstate_frag))
          trace_local = 0d0
          face_pt = 0
          do iz_face=1,dc%nxyz_domain(3)
          do iy_face=1,dc%nxyz_domain(2)
          do ix_face=1,dc%nxyz_domain(1)
            if(face_axis_index([ix_face,iy_face,iz_face], axis) /= face_coord(axis, side)) cycle
            face_pt = face_pt + 1
            do ispin=1,nspin
              do io=1,n_basis(dc%i_frag,ispin)
                trace_local(face_pt,ispin,io) = local_basis_value(ix_face,iy_face,iz_face,ispin,io)
                trace_local(face_pt,ispin,dc%nstate_frag+io) = &
                & local_basis_dn(ix_face,iy_face,iz_face,ispin,io,axis,side)
              end do
            end do
          end do
          end do
          end do
          do ispin=1,nspin
!$omp parallel do private(io,jo,face_pt,ix_face,iy_face,iz_face,u_l,v_l,dnu_l,dnv_l, &
!$omp& term_sum,term_face,term_vsum,term_vface) schedule(static)
          do io=1,n_basis(dc%i_frag,ispin)
          do jo=1,n_basis(dc%i_frag,ispin)
            term_sum = 0d0
            term_vsum = 0d0
            face_pt = 0
            do iz_face=1,dc%nxyz_domain(3)
            do iy_face=1,dc%nxyz_domain(2)
            do ix_face=1,dc%nxyz_domain(1)
              if(face_axis_index([ix_face,iy_face,iz_face], axis) /= face_coord(axis, side)) cycle
              face_pt = face_pt + 1
              v_l = trace_local(face_pt,ispin,io)
              dnv_l = trace_local(face_pt,ispin,dc%nstate_frag+io)
              u_l = trace_local(face_pt,ispin,jo)
              dnu_l = trace_local(face_pt,ispin,dc%nstate_frag+jo)
              term_face = (-0.25d0 * v_l * dnu_l - 0.25d0 * dnv_l * u_l + alpha * v_l * u_l) * area_weight
              term_vface = -0.5d0 * real(side,8) * v_l * u_l * area_weight
              term_sum = term_sum + term_face
              term_vsum = term_vsum + term_vface
            end do
            end do
            end do
            if(use_surface_self_hamiltonian_mode()) &
              mat_H_local(io,jo,ispin) = mat_H_local(io,jo,ispin) + term_sum
            mat_H_surface_self(io,jo,ispin) = mat_H_surface_self(io,jo,ispin) + term_sum
            mat_V_local(axis,io,jo,ispin) = mat_V_local(axis,io,jo,ispin) + term_vsum
          end do
          end do
!$omp end parallel do
          end do
          ! halo%mat_H_local(io_remote, jo_local) stores the transpose of
          ! the local-to-remote surface block, i.e. H(remote,local).  The
          ! term below is evaluated as H(local,remote) with receiver-normal
          ! traces and is placed transposed by construction for EigenExa.
          do ispin=1,nspin
!$omp parallel do private(io,jo,face_pt,ix_face,iy_face,iz_face,v_l,dnv_l,u_r,dnu_r, &
!$omp& term_sum,term_face,term_vsum,term_vface,term_xisum,term_xiface,xi_rel,idir) schedule(static)
          do io=1,n_basis(jfrag,ispin)
          do jo=1,n_basis(dc%i_frag,ispin)
            term_sum = 0d0
            term_vsum = 0d0
            term_xisum(1:3) = 0d0
            face_pt = 0
            do iz_face=1,dc%nxyz_domain(3)
            do iy_face=1,dc%nxyz_domain(2)
            do ix_face=1,dc%nxyz_domain(1)
              if(face_axis_index([ix_face,iy_face,iz_face], axis) /= face_coord(axis, side)) cycle
              face_pt = face_pt + 1
              v_l = trace_local(face_pt,ispin,jo)
              dnv_l = trace_local(face_pt,ispin,dc%nstate_frag+jo)
              u_r = halo(i_halo)%trace_recv(face_pt,ispin,io)
              dnu_r = halo(i_halo)%trace_recv(face_pt,ispin,dc%nstate_frag+io)
              term_face = (-0.25d0 * v_l * dnu_r + 0.25d0 * dnv_l * u_r - alpha * v_l * u_r) * area_weight
              term_vface = 0.5d0 * real(side,8) * v_l * u_r * area_weight
              xi_rel(1) = (dble(ix_face) - 0.5d0 * dble(dc%nxyz_domain(1) + 1)) * system%hgs(1)
              xi_rel(2) = (dble(iy_face) - 0.5d0 * dble(dc%nxyz_domain(2) + 1)) * system%hgs(2)
              xi_rel(3) = (dble(iz_face) - 0.5d0 * dble(dc%nxyz_domain(3) + 1)) * system%hgs(3)
              xi_rel(axis) = 0d0
              term_xiface = v_l * u_r * area_weight
              term_sum = term_sum + term_face
              term_vsum = term_vsum + term_vface
              do idir=1,3
                term_xisum(idir) = term_xisum(idir) + xi_rel(idir) * term_xiface
              end do
            end do
            end do
            end do
            halo(i_halo)%mat_H_local(io,jo,ispin) = halo(i_halo)%mat_H_local(io,jo,ispin) + term_sum
            halo(i_halo)%mat_H_surface_cross(io,jo,ispin) = &
              halo(i_halo)%mat_H_surface_cross(io,jo,ispin) + term_sum
            halo(i_halo)%mat_V_local(axis,io,jo,ispin) = halo(i_halo)%mat_V_local(axis,io,jo,ispin) + term_vsum
            halo(i_halo)%mat_Xi_flux_local(1:3,io,jo,ispin) = &
              halo(i_halo)%mat_Xi_flux_local(1:3,io,jo,ispin) + term_xisum(1:3)
          end do
          end do
!$omp end parallel do
          end do
          deallocate(trace_local)
          if(dc%id_tot==0) write(*,*) "Halo communication #",i_halo,": done"
        end do
        call release_surface_trace_halo()
      end if ! dc%id_frag==0
      do ispin=1,nspin
        do idir=1,3
          do io=1,n_basis(dc%i_frag,ispin)
            mat_V_local(idir,io,io,ispin) = 0d0
            do jo=io+1,n_basis(dc%i_frag,ispin)
              pavg = 0.5d0 * (mat_V_local(idir,io,jo,ispin) - mat_V_local(idir,jo,io,ispin))
              mat_V_local(idir,io,jo,ispin) = pavg
              mat_V_local(idir,jo,io,ispin) = -pavg
            end do
          end do
        end do
      end do
      deallocate(hf)
      deallocate(basis_grad_all)
      deallocate(basis_nonlocal_core)

    end subroutine calc_hamiltonian_matrix

    subroutine build_fragment_nonlocal_basis_action(nonlocal_action)
      use communication, only: comm_summation
      use nonlocal_potential, only: dpseudo
      use salmon_global, only: yn_jm
      implicit none
      real(8), intent(out) :: nonlocal_action(:,:,:,:,:)
      real(8), allocatable :: local_action(:,:,:,:,:)
      integer :: ix1,ix2,iy1,iy2,iz1,iz2

      nonlocal_action=0.0d0
      if(yn_jm/='n') return
      ix1=max(mg%is(1),1); ix2=min(mg%ie(1),dc%nxyz_domain(1))
      iy1=max(mg%is(2),1); iy2=min(mg%ie(2),dc%nxyz_domain(2))
      iz1=max(mg%is(3),1); iz2=min(mg%ie(3),dc%nxyz_domain(3))
      sttpsi%rwf=0.0d0
      shpsi%rwf=0.0d0
      do io=info%io_s,info%io_e
        if(io>n_basis(dc%i_frag,1)) cycle
        do ispin=1,nspin
          sttpsi%rwf(ix1:ix2,iy1:iy2,iz1:iz2,ispin,io,1,1)= &
            f_basis(ix1:ix2,iy1:iy2,iz1:iz2,ispin,io)
        end do
      end do
      call dpseudo(sttpsi,shpsi,info,nspin,ppg)
      allocate(local_action,source=nonlocal_action)
      local_action=0.0d0
      do io=info%io_s,info%io_e
        if(io>n_basis(dc%i_frag,1)) cycle
        do ispin=1,nspin
          local_action(ix1:ix2,iy1:iy2,iz1:iz2,ispin,io)= &
            shpsi%rwf(ix1:ix2,iy1:iy2,iz1:iz2,ispin,io,1,1)
        end do
      end do
      call comm_summation(local_action,nonlocal_action,size(nonlocal_action),info%icomm_rko)
      deallocate(local_action)
    end subroutine build_fragment_nonlocal_basis_action

    integer function face_point_count(axis) result(npts)
      implicit none
      integer,intent(in) :: axis

      select case(axis)
      case(1)
        npts = dc%nxyz_domain(2) * dc%nxyz_domain(3)
      case(2)
        npts = dc%nxyz_domain(1) * dc%nxyz_domain(3)
      case(3)
        npts = dc%nxyz_domain(1) * dc%nxyz_domain(2)
      case default
        npts = 0
      end select
    end function face_point_count

    integer function face_coord(axis, side) result(idx)
      implicit none
      integer,intent(in) :: axis, side

      if(side > 0) then
        idx = dc%nxyz_domain(axis)
      else
        idx = 1
      end if
    end function face_coord

    integer function face_axis_index(idx3, axis) result(idx)
      implicit none
      integer,intent(in) :: idx3(3), axis

      idx = idx3(axis)
    end function face_axis_index

    real(8) function local_basis_value(ix,iy,iz,ispin,ibasis) result(val)
      implicit none
      integer,intent(in) :: ix,iy,iz,ispin,ibasis
      integer :: raw_io, io_lb, io_ub, sx, sy, sz
      real(8) :: coef_val

      val = 0d0
      if(ibasis < 1 .or. ibasis > dc%nstate_frag) return
      if(ispin < 1 .or. ispin > nspin) return
      if(sawf_explicit_basis_active) then
        sx=ix+dc%nxyz_buffer(1)
        sy=iy+dc%nxyz_buffer(2)
        sz=iz+dc%nxyz_buffer(3)
        if(sx<1 .or. sx>size(sawf_explicit_buffer,1)) return
        if(sy<1 .or. sy>size(sawf_explicit_buffer,2)) return
        if(sz<1 .or. sz>size(sawf_explicit_buffer,3)) return
        val=sawf_explicit_buffer(sx,sy,sz,ispin,ibasis)
        return
      end if
      if(ix >= 1 .and. ix <= dc%nxyz_domain(1) .and. &
         iy >= 1 .and. iy <= dc%nxyz_domain(2) .and. &
         iz >= 1 .and. iz <= dc%nxyz_domain(3)) then
        val = f_basis(ix,iy,iz,ispin,ibasis)
        return
      end if
      sx = buffered_local_index(ix, 1)
      sy = buffered_local_index(iy, 2)
      sz = buffered_local_index(iz, 3)
      if(sx < lbound(spsi%rwf,1) .or. sx > ubound(spsi%rwf,1)) return
      if(sy < lbound(spsi%rwf,2) .or. sy > ubound(spsi%rwf,2)) return
      if(sz < lbound(spsi%rwf,3) .or. sz > ubound(spsi%rwf,3)) return
      io_lb = lbound(spsi%rwf,5)
      io_ub = ubound(spsi%rwf,5)
      do raw_io=max(1,io_lb),min(dc%nstate_frag,io_ub)
        coef_val = basis_transform(raw_io,ibasis,ispin)
        if(abs(coef_val) <= 0d0) cycle
        val = val + coef_val * spsi%rwf(sx,sy,sz,ispin,raw_io,1,1)
      end do
    end function local_basis_value

    integer function buffered_local_index(idx, axis) result(mapped)
      implicit none
      integer,intent(in) :: idx, axis

      if(idx < 1) then
        mapped = dc%nxyz_domain(axis) + dc%nxyz_buffer(axis) + (1 - idx)
      else
        mapped = idx
      end if
    end function buffered_local_index

    real(8) function local_basis_grad(ix,iy,iz,ispin,ibasis,axis) result(grad_axis)
      implicit none
      integer,intent(in) :: ix,iy,iz,ispin,ibasis,axis
      integer :: dist

      grad_axis = 0d0
      do dist=1,size(stencil%coef_nab,1)
        select case(axis)
        case(1)
          grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
          & (local_basis_value(ix+dist,iy,iz,ispin,ibasis) - local_basis_value(ix-dist,iy,iz,ispin,ibasis))
        case(2)
          grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
          & (local_basis_value(ix,iy+dist,iz,ispin,ibasis) - local_basis_value(ix,iy-dist,iz,ispin,ibasis))
        case(3)
          grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
          & (local_basis_value(ix,iy,iz+dist,ispin,ibasis) - local_basis_value(ix,iy,iz-dist,ispin,ibasis))
        end select
      end do
    end function local_basis_grad

    real(8) function local_basis_value_core(ix,iy,iz,ispin,ibasis) result(val)
      implicit none
      integer,intent(in) :: ix,iy,iz,ispin,ibasis

      val = 0d0
      if(ibasis < 1 .or. ibasis > dc%nstate_frag) return
      if(ispin < 1 .or. ispin > nspin) return
      if(ix < 1 .or. ix > dc%nxyz_domain(1)) return
      if(iy < 1 .or. iy > dc%nxyz_domain(2)) return
      if(iz < 1 .or. iz > dc%nxyz_domain(3)) return
      val = f_basis(ix,iy,iz,ispin,ibasis)
    end function local_basis_value_core

    real(8) function local_basis_kinetic_core(ix,iy,iz,ispin,ibasis) result(tval)
      implicit none
      integer,intent(in) :: ix,iy,iz,ispin,ibasis
      integer :: dist, axis
      real(8) :: v

      v = 0d0
      do axis=1,3
        do dist=1,size(stencil%coef_lap,1)
          select case(axis)
          case(1)
            v = v + stencil%coef_lap(dist,axis) * &
            & (local_basis_value_core(ix+dist,iy,iz,ispin,ibasis) + &
            &  local_basis_value_core(ix-dist,iy,iz,ispin,ibasis))
          case(2)
            v = v + stencil%coef_lap(dist,axis) * &
            & (local_basis_value_core(ix,iy+dist,iz,ispin,ibasis) + &
            &  local_basis_value_core(ix,iy-dist,iz,ispin,ibasis))
          case(3)
            v = v + stencil%coef_lap(dist,axis) * &
            & (local_basis_value_core(ix,iy,iz+dist,ispin,ibasis) + &
            &  local_basis_value_core(ix,iy,iz-dist,ispin,ibasis))
          end select
        end do
      end do
      tval = stencil%coef_lap0 * local_basis_value_core(ix,iy,iz,ispin,ibasis) - 0.5d0 * v
    end function local_basis_kinetic_core

    real(8) function strong_core_kinetic_integral(io,jo,ispin,l) result(val)
      implicit none
      integer,intent(in) :: io,jo,ispin,l(3)
      integer :: ix,iy,iz

      val = 0d0
      do iz=1,l(3)
      do iy=1,l(2)
      do ix=1,l(1)
        val = val + f_basis(ix,iy,iz,ispin,io) * local_basis_kinetic_core(ix,iy,iz,ispin,jo)
      end do
      end do
      end do
      val = val * hvol
    end function strong_core_kinetic_integral

    real(8) function weak_buffer_kinetic_integral(io,jo,ispin,l) result(val)
      implicit none
      integer,intent(in) :: io,jo,ispin,l(3)
      integer :: ix,iy,iz,axis

      val = 0d0
      do axis=1,3
      do iz=1,l(3)
      do iy=1,l(2)
      do ix=1,l(1)
        val = val + 0.5d0 * local_basis_grad(ix,iy,iz,ispin,io,axis) &
        &                 * local_basis_grad(ix,iy,iz,ispin,jo,axis)
      end do
      end do
      end do
      end do
      val = val * hvol
    end function weak_buffer_kinetic_integral

    real(8) function local_basis_dn(ix,iy,iz,ispin,ibasis,axis,side) result(dn)
      implicit none
      integer,intent(in) :: ix,iy,iz,ispin,ibasis,axis,side

      dn = real(side,8) * local_basis_grad(ix,iy,iz,ispin,ibasis,axis)
    end function local_basis_dn

    real(8) function volume_velocity_integral(io,jo,ispin,axis,l) result(val)
      implicit none
      integer,intent(in) :: io,jo,ispin,axis,l(3)
      integer :: ix,iy,iz

      val = 0d0
      do iz=1,l(3)
        do iy=1,l(2)
          do ix=1,l(1)
            val = val + f_basis(ix,iy,iz,ispin,io) * local_basis_grad(ix,iy,iz,ispin,jo,axis)
          end do
        end do
      end do
      val = val * hvol
    end function volume_velocity_integral

    integer function active_laplacian_radius(axis) result(radius)
      implicit none
      integer,intent(in) :: axis
      integer :: dist

      radius = 0
      do dist=1,size(stencil%coef_lap,1)
        if(abs(stencil%coef_lap(dist,axis)) > 0d0) radius = dist
      end do
    end function active_laplacian_radius

#ifdef USE_EIGENEXA
    subroutine diag_eigenexa
      use communication, only: comm_bcast, comm_summation
      use eigen_subdiag_sub, only: eigen_dsyev
      use eigen_libs_mod
      use salmon_global, only: nelec, nelec_spin
      implicit none
      integer, parameter :: coef_gather_target_elems = 2000000
      integer, parameter :: velocity_diag_max_dim = 1024
      integer :: n,nx,ny,ix_s,ix_e,iy_s,iy_e,ix_loc,iy_loc,ifrag_x,ifrag_y,io_x,io_y
      integer :: nnod,x_nnod,y_nnod,inod,x_inod,y_inod
      integer :: jfrag_halo(n_halo)
      integer, allocatable :: io_array(:),ifrag_array(:)
      integer :: c0, c1, ncol, nstate_chunk
      integer :: target_frag, nbasis_diag, level, max_level, state_col, n_entry
      integer :: pass, best, tmp_i
      integer, parameter :: nsample_max = 3
      integer :: nsample, isample, sample_state(nsample_max)
      integer :: nstate_use, nocc_diag, occ, virt, idir_diag
      integer :: nocc_nelec
      real(8) :: eps, hval, rel_col, rel_row, rel_col_max, rel_row_max
      real(8) :: occ_weight, gap, amp, strength_pair
      real(8) :: strength_total, strength_max, sum_gap_weighted, sum_inv_gap, occ_sum
      real(8) :: mean_gap, fsum_ratio
      real(8), allocatable :: h_div(:,:), h_ref_div(:,:), v_div(:,:), h(:,:,:)
      real(8), allocatable :: h_full_local(:,:), h_full(:,:)
      real(8), allocatable :: h_block(:,:,:), h_local_diag(:,:), evec_local(:,:), eval_local(:)
      real(8), allocatable :: eval_list(:)
      real(8), allocatable :: coef_state_norm(:), coef_state_norm_alt(:)
      real(8), allocatable :: v_tmp1(:,:), v_tmp2(:,:)
      integer, allocatable :: frag_list(:), level_list(:)
      real(8), allocatable :: v_col_local(:,:), v_col(:,:), hv_col_local(:,:), hv_col(:,:)
      real(8), allocatable :: v_row_local(:,:), v_row(:,:), hv_row_local(:,:), hv_row(:,:)
      real(8), allocatable :: eigvec_local(:,:), eigvec(:,:), vop(:,:), vwork(:)
      logical, allocatable :: repair_state(:)
      logical :: use_transpose_export, block_diag_h

      allocate(h(dc%nstate_frag,dc%nstate_frag,0:n_halo))
      block_diag_h = use_block_diag_hamiltonian_mode()
      if(dc%id_tot==0 .and. block_diag_h) then
        write(*,'(1x,a)') "[DC-LCFO-FLUX] block-diag H diagonalization: enabled"
      end if
      do ispin=1,nspin
        if(dc%id_tot==0) write(*,*) "eigenexa diag, #dim=",n_mat(ispin)
        n = n_mat(ispin)

        allocate(io_array(n),ifrag_array(n))
        do ifrag=1,dc%n_frag
          do io=1,n_basis(ifrag,ispin) ; i = index_basis(io,ifrag,ispin)
            io_array(i) = io
            ifrag_array(i) = ifrag
          end do
        end do

        call eigen_init(dc%icomm_tot)
        call eigen_get_matdims( n, nx, ny )
        call eigen_get_procs( nnod, x_nnod, y_nnod )
        call eigen_get_id   ( inod, x_inod, y_inod )
        allocate( h_div(nx,ny), h_ref_div(nx,ny), v_div(nx,ny) )
        ix_s = eigen_loop_start( 1, x_nnod, x_inod )
        ix_e = eigen_loop_end  ( n, x_nnod, x_inod )
        iy_s = eigen_loop_start( 1, y_nnod, y_inod )
        iy_e = eigen_loop_end  ( n, y_nnod, y_inod )
        nstate_chunk = max(1, min(dc%nstate_tot, &
          max(1, coef_gather_target_elems / max(1, dc%nstate_frag))))
        allocate(v_tmp1(dc%nstate_frag,nstate_chunk))
        allocate(v_tmp2(dc%nstate_frag,nstate_chunk))
        allocate(coef_state_norm(dc%nstate_tot), coef_state_norm_alt(dc%nstate_tot))
        allocate(repair_state(dc%nstate_tot))

        h_div = 0d0
        do ifrag=1,dc%n_frag
          if(ifrag==dc%i_frag .and. dc%id_frag==0) then
            h(:,:,0) = mat_H_local(:,:,ispin)
            do i_halo=1,n_halo
              jfrag_halo(i_halo) = halo(i_halo)%ifrag_src ! src fragment (recv)
              h(:,:,i_halo) = halo(i_halo)%mat_H_local(:,:,ispin)
            end do
            if(sawf_explicit_basis_active) write(*,'(1x,a,i0,a,6(i0,1x))') &
              '[DC-LCFO-SAWF-HALO] owner=',ifrag,' sources=',jfrag_halo
          end if
          call comm_bcast( h, dc%icomm_tot, id_array(ifrag) )
          call comm_bcast( jfrag_halo, dc%icomm_tot, id_array(ifrag) )
          do iy_loc=iy_s,iy_e
            iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
            if(iy > n) cycle
            ifrag_y = ifrag_array(iy)
            io_y = io_array(iy)
            do ix_loc=ix_s,ix_e
              ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
              if(ix > n) cycle
              ifrag_x = ifrag_array(ix)
              io_x = io_array(ix)
              if(block_diag_h) then
                if(ifrag_x == ifrag_y) then
                  if(ifrag_x == ifrag) then
                    h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) + 0.5d0 * (h(io_x,io_y,0) + h(io_y,io_x,0))
                    do i_halo=1,n_halo
                      h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                      & + 0.25d0 * (h(io_y,io_x,i_halo) + h(io_x,io_y,i_halo))
                    end do
                  end if
                  do i_halo=1,n_halo
                    if(ifrag_x == jfrag_halo(i_halo)) then
                      h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                      & + 0.25d0 * (h(io_x,io_y,i_halo) + h(io_y,io_x,i_halo))
                    end if
                  end do
                end if
              else
                if(ifrag_x == ifrag .and. ifrag_y == ifrag) then
                  h_div(ix_loc,iy_loc) = h(io_x,io_y,0)
                end if
                do i_halo=1,n_halo
                  if( ifrag_x == jfrag_halo(i_halo) .and. ifrag_y == ifrag ) then
                    h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                    & + 0.5d0 * h(io_x,io_y,i_halo) ! 0.5d0: avoid double counting face-neighbor pairs
                  else if( ifrag_x == ifrag .and. ifrag_y == jfrag_halo(i_halo) ) then
                    h_div(ix_loc,iy_loc) = h_div(ix_loc,iy_loc) &
                    & + 0.5d0 * h(io_y,io_x,i_halo) ! 0.5d0: avoid double counting face-neighbor pairs
                  end if
                end do ! i_halo
              end if
            end do ! ix_loc
          end do ! iy_loc
        end do ! ifrag
        if(dc%id_tot==0) write(*,*) "h_div: done"

        h_ref_div = h_div
        if(sawf_explicit_basis_active) then
          allocate(h_full_local(n,n),h_full(n,n))
          h_full_local=0.0d0
          do iy_loc=iy_s,iy_e
            iy=eigen_translate_l2g(iy_loc,y_nnod,y_inod)
            if(iy>n) cycle
            do ix_loc=ix_s,ix_e
              ix=eigen_translate_l2g(ix_loc,x_nnod,x_inod)
              if(ix>n) cycle
              h_full_local(ix,iy)=h_ref_div(ix_loc,iy_loc)
            end do
          end do
          call comm_summation(h_full_local,h_full,n*n,dc%icomm_tot)
          if(dc%id_tot==0) call diagnose_sawf_hamiltonian_hermiticity( &
            h_full,ifrag_array,io_array,mat_H_local(:,:,ispin),n_basis(dc%i_frag,ispin))
          deallocate(h_full_local,h_full)
        end if
        call eigen_sx(n, n, h_div, nx, esp_tot(1:n,ispin), v_div, nx)
        if(dc%id_tot==0) write(*,*) "eigen_sx: done"
        nocc_nelec = occupied_index_from_input(ispin)
        if(dc%id_tot==0) call print_flux_gap_diagnostic("full", ispin, nocc_nelec, n)

        nsample = 1
        sample_state(1) = 1
        i = min(n, max(1, dc%nstate_tot / 2))
        if(i /= sample_state(1)) then
          nsample = nsample + 1
          sample_state(nsample) = i
        end if
        i = min(n, dc%nstate_tot)
        if(all(sample_state(1:nsample) /= i)) then
          nsample = nsample + 1
          sample_state(nsample) = i
        end if

        allocate(v_col_local(n,nsample), v_col(n,nsample), hv_col_local(n,nsample), hv_col(n,nsample))
        allocate(v_row_local(n,nsample), v_row(n,nsample), hv_row_local(n,nsample), hv_row(n,nsample))
        v_col_local = 0d0
        v_row_local = 0d0
        do iy_loc=iy_s,iy_e
          iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
          if(iy > n) cycle
          do ix_loc=ix_s,ix_e
            ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
            if(ix > n) cycle
            do isample=1,nsample
              if(iy == sample_state(isample)) v_col_local(ix,isample) = v_div(ix_loc,iy_loc)
              if(ix == sample_state(isample)) v_row_local(iy,isample) = v_div(ix_loc,iy_loc)
            end do
          end do
        end do
        call comm_summation(v_col_local, v_col, n*nsample, dc%icomm_tot)
        call comm_summation(v_row_local, v_row, n*nsample, dc%icomm_tot)

        hv_col_local = 0d0
        hv_row_local = 0d0
        do iy_loc=iy_s,iy_e
          iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
          if(iy > n) cycle
          do ix_loc=ix_s,ix_e
            ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
            if(ix > n) cycle
            hval = h_ref_div(ix_loc,iy_loc)
            if(hval == 0d0) cycle
            do isample=1,nsample
              hv_col_local(ix,isample) = hv_col_local(ix,isample) + hval * v_col(iy,isample)
              hv_row_local(ix,isample) = hv_row_local(ix,isample) + hval * v_row(iy,isample)
            end do
          end do
        end do
        call comm_summation(hv_col_local, hv_col, n*nsample, dc%icomm_tot)
        call comm_summation(hv_row_local, hv_row, n*nsample, dc%icomm_tot)

        rel_col_max = 0d0
        rel_row_max = 0d0
        do isample=1,nsample
          eps = esp_tot(sample_state(isample),ispin)
          rel_col = sqrt(sum((hv_col(:,isample) - eps * v_col(:,isample))**2) &
            / max(1.0d-300, sum(hv_col(:,isample)**2)))
          rel_row = sqrt(sum((hv_row(:,isample) - eps * v_row(:,isample))**2) &
            / max(1.0d-300, sum(hv_row(:,isample)**2)))
          rel_col_max = max(rel_col_max, rel_col)
          rel_row_max = max(rel_row_max, rel_row)
        end do
        use_transpose_export = (rel_row_max < rel_col_max * 1.0d-3)
        if(dc%id_tot==0) then
          write(*,'(1x,a,i0,2(a,1pe12.5),a,l1)') &
            "[DC-LCFO-FLUX] EigenExa vector check: ispin=", ispin, &
            " rel_col=", rel_col_max, " rel_row=", rel_row_max, &
            " transpose_export=", use_transpose_export
        end if
        deallocate(v_col_local, v_col, hv_col_local, hv_col)
        deallocate(v_row_local, v_row, hv_row_local, hv_row)

        nstate_use = min(dc%nstate_tot, n)
        if(n <= velocity_diag_max_dim .and. nstate_use <= velocity_diag_max_dim) then
          idir_diag = 3
          allocate(eigvec_local(n,nstate_use), eigvec(n,nstate_use))
          allocate(vop(n,n), vwork(n))
          eigvec_local = 0d0
          do iy_loc=iy_s,iy_e
            iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
            if(iy < 1 .or. iy > nstate_use) cycle
            do ix_loc=ix_s,ix_e
              ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
              if(ix > n) cycle
              eigvec_local(ix,iy) = v_div(ix_loc,iy_loc)
            end do
          end do
          call comm_summation(eigvec_local, eigvec, n*nstate_use, dc%icomm_tot)

          vop = 0d0
          do ifrag=1,dc%n_frag
            if(ifrag==dc%i_frag .and. dc%id_frag==0) then
              h(:,:,0) = mat_V_local(idir_diag,:,:,ispin)
              do i_halo=1,n_halo
                jfrag_halo(i_halo) = halo(i_halo)%ifrag_src
                h(:,:,i_halo) = halo(i_halo)%mat_V_local(idir_diag,:,:,ispin)
              end do
            end if
            call comm_bcast( h, dc%icomm_tot, id_array(ifrag) )
            call comm_bcast( jfrag_halo, dc%icomm_tot, id_array(ifrag) )
            do iy=1,n
              ifrag_y = ifrag_array(iy)
              io_y = io_array(iy)
              do ix=1,n
                ifrag_x = ifrag_array(ix)
                io_x = io_array(ix)
                if(ifrag_x == ifrag .and. ifrag_y == ifrag) then
                  if(block_diag_h) then
                    vop(ix,iy) = vop(ix,iy) + h(io_x,io_y,0)
                  else
                    vop(ix,iy) = h(io_x,io_y,0)
                  end if
                end if
                if(block_diag_h) then
                  if(ifrag_x == ifrag_y) then
                    if(ifrag_x == ifrag) then
                      do i_halo=1,n_halo
                        vop(ix,iy) = vop(ix,iy) + 0.5d0 * h(io_y,io_x,i_halo)
                      end do
                    end if
                    do i_halo=1,n_halo
                      if(ifrag_x == jfrag_halo(i_halo)) then
                        vop(ix,iy) = vop(ix,iy) + 0.5d0 * h(io_x,io_y,i_halo)
                      end if
                    end do
                  end if
                else
                  do i_halo=1,n_halo
                    if(ifrag_x == jfrag_halo(i_halo) .and. ifrag_y == ifrag) then
                      vop(ix,iy) = vop(ix,iy) + 0.5d0 * h(io_x,io_y,i_halo)
                    else if(ifrag_x == ifrag .and. ifrag_y == jfrag_halo(i_halo)) then
                      vop(ix,iy) = vop(ix,iy) + 0.5d0 * h(io_y,io_x,i_halo)
                    end if
                  end do
                end if
              end do
            end do
          end do

          nocc_diag = 0
          occ_sum = 0d0
          do occ=1,nstate_use
            occ_weight = max(0d0, system%rocc(occ,1,ispin))
            if(occ_weight > 1d-12) nocc_diag = occ
            occ_sum = occ_sum + occ_weight
          end do
          strength_total = 0d0
          strength_max = 0d0
          sum_gap_weighted = 0d0
          sum_inv_gap = 0d0
          do virt=nocc_diag+1,nstate_use
            vwork(:) = matmul(vop, eigvec(:,virt))
            do occ=1,nocc_diag
              occ_weight = max(0d0, system%rocc(occ,1,ispin))
              if(occ_weight <= 1d-12) cycle
              gap = esp_tot(virt,ispin) - esp_tot(occ,ispin)
              amp = sum(eigvec(:,occ) * vwork(:))
              strength_pair = occ_weight * amp * amp
              strength_total = strength_total + strength_pair
              strength_max = max(strength_max, strength_pair)
              sum_gap_weighted = sum_gap_weighted + strength_pair * gap
              if(gap > 1d-12) sum_inv_gap = sum_inv_gap + strength_pair / gap
            end do
          end do
          if(strength_total > 0d0) then
            mean_gap = sum_gap_weighted / strength_total
          else
            mean_gap = 0d0
          end if
          if(occ_sum > 0d0) then
            fsum_ratio = 2d0 * sum_inv_gap / occ_sum
          else
            fsum_ratio = 0d0
          end if
          if(dc%id_tot==0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,4(a,1pe13.5))') &
              "[DC-LCFO-FLUX-TRANSITION] idir=", idir_diag, " ispin=", ispin, &
              " nocc=", nocc_diag, " nvirt=", max(0,nstate_use-nocc_diag), &
              " total=", strength_total, " max_pair=", strength_max, &
              " mean_gap_eV=", mean_gap * 27.211386245988d0, &
              " fsum_ratio=", fsum_ratio
          end if
          deallocate(eigvec_local, eigvec, vop, vwork)
        else if(dc%id_tot==0) then
          write(*,'(1x,a,i0,a,i0,a,i0)') &
            "[DC-LCFO-FLUX-TRANSITION] skipped: n=", n, &
            " nstate=", nstate_use, " max_auto=", velocity_diag_max_dim
        end if

        coef_state_norm = 0d0
        do ifrag=1,dc%n_frag
          do c0=1,dc%nstate_tot,nstate_chunk
            c1 = min(dc%nstate_tot, c0+nstate_chunk-1)
            ncol = c1-c0+1
            v_tmp1(:,1:ncol) = 0d0
            if(use_transpose_export) then
              do iy_loc=iy_s,iy_e
                iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
                if(iy > n) cycle
                ifrag_y = ifrag_array(iy)
                if(ifrag_y /= ifrag) cycle
                io_y = io_array(iy)
                do ix_loc=ix_s,ix_e
                  ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
                  if(ix < c0 .or. ix > c1 .or. ix > n) cycle
                  v_tmp1(io_y,ix-c0+1) = v_div(ix_loc,iy_loc)
                end do
              end do
            else
              do iy_loc=iy_s,iy_e
                iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
                if(iy < c0 .or. iy > c1 .or. iy > n) cycle
                do ix_loc=ix_s,ix_e
                  ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
                  if(ix > n) cycle
                  ifrag_x = ifrag_array(ix)
                  io_x = io_array(ix)
                  if(ifrag_x == ifrag) then
                    v_tmp1(io_x,iy-c0+1) = v_div(ix_loc,iy_loc)
                  end if
                end do
              end do
            end if
            call comm_summation(v_tmp1(:,1:ncol),v_tmp2(:,1:ncol), &
              dc%nstate_frag*ncol,dc%icomm_tot)
            do i=1,ncol
              coef_state_norm(c0+i-1) = coef_state_norm(c0+i-1) + sum(v_tmp2(:,i)**2)
            end do
            if(ifrag==dc%i_frag .and. dc%id_frag==0) then
              coef_wf(:,c0:c1,ispin) = v_tmp2(:,1:ncol)
            end if
          end do
        end do ! ifrag

        repair_state(:) = (coef_state_norm(:) <= 1.0d-12)
        if(any(repair_state)) then
          if(dc%id_tot==0) then
            write(*,'(1x,a,i0,a,1pe13.5)') &
              "[DC-LCFO-FLUX] repairing near-zero exported eigenvector columns: count=", &
              count(repair_state), " min_norm2=", minval(coef_state_norm)
          end if
          coef_state_norm_alt = 0d0
          do ifrag=1,dc%n_frag
            do c0=1,dc%nstate_tot,nstate_chunk
              c1 = min(dc%nstate_tot, c0+nstate_chunk-1)
              ncol = c1-c0+1
              if(.not. any(repair_state(c0:c1))) cycle
              v_tmp1(:,1:ncol) = 0d0
              if(use_transpose_export) then
                do iy_loc=iy_s,iy_e
                  iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
                  if(iy < c0 .or. iy > c1 .or. iy > n) cycle
                  do ix_loc=ix_s,ix_e
                    ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
                    if(ix > n) cycle
                    ifrag_x = ifrag_array(ix)
                    io_x = io_array(ix)
                    if(ifrag_x == ifrag) then
                      v_tmp1(io_x,iy-c0+1) = v_div(ix_loc,iy_loc)
                    end if
                  end do
                end do
              else
                do iy_loc=iy_s,iy_e
                  iy = eigen_translate_l2g(iy_loc, y_nnod, y_inod)
                  if(iy > n) cycle
                  ifrag_y = ifrag_array(iy)
                  if(ifrag_y /= ifrag) cycle
                  io_y = io_array(iy)
                  do ix_loc=ix_s,ix_e
                    ix = eigen_translate_l2g(ix_loc, x_nnod, x_inod)
                    if(ix < c0 .or. ix > c1 .or. ix > n) cycle
                    v_tmp1(io_y,ix-c0+1) = v_div(ix_loc,iy_loc)
                  end do
                end do
              end if
              call comm_summation(v_tmp1(:,1:ncol),v_tmp2(:,1:ncol), &
                dc%nstate_frag*ncol,dc%icomm_tot)
              do i=1,ncol
                if(.not. repair_state(c0+i-1)) cycle
                coef_state_norm_alt(c0+i-1) = coef_state_norm_alt(c0+i-1) + sum(v_tmp2(:,i)**2)
                if(ifrag==dc%i_frag .and. dc%id_frag==0) then
                  coef_wf(:,c0+i-1,ispin) = v_tmp2(:,i)
                end if
              end do
            end do
          end do
          if(dc%id_tot==0) then
            write(*,'(1x,a,1pe13.5)') &
              "[DC-LCFO-FLUX] repaired eigenvector column min_norm2=", &
              minval(coef_state_norm_alt, mask=repair_state)
          end if
        end if

        if(block_diag_h) then
          if(dc%id_tot==0) then
            write(*,'(1x,a)') &
              "[DC-LCFO-FLUX] exporting fragment-local block eigenvectors for block-diag H"
          end if
          allocate(h_block(dc%nstate_frag,dc%nstate_frag,dc%n_frag))
          allocate(eval_list(n), frag_list(n), level_list(n))
          h_block = 0d0
          eval_list = huge(1.0d0)
          frag_list = 0
          level_list = 0
          do ifrag=1,dc%n_frag
            if(ifrag==dc%i_frag .and. dc%id_frag==0) then
              h(:,:,0) = mat_H_local(:,:,ispin)
              do i_halo=1,n_halo
                jfrag_halo(i_halo) = halo(i_halo)%ifrag_src
                h(:,:,i_halo) = halo(i_halo)%mat_H_local(:,:,ispin)
              end do
            end if
            call comm_bcast( h, dc%icomm_tot, id_array(ifrag) )
            call comm_bcast( jfrag_halo, dc%icomm_tot, id_array(ifrag) )
            do target_frag=1,dc%n_frag
              if(target_frag == ifrag) then
                do io_x=1,n_basis(target_frag,ispin)
                do io_y=1,n_basis(target_frag,ispin)
                  h_block(io_x,io_y,target_frag) = h_block(io_x,io_y,target_frag) &
                  & + 0.5d0 * (h(io_x,io_y,0) + h(io_y,io_x,0))
                  do i_halo=1,n_halo
                    h_block(io_x,io_y,target_frag) = h_block(io_x,io_y,target_frag) &
                    & + 0.25d0 * (h(io_y,io_x,i_halo) + h(io_x,io_y,i_halo))
                  end do
                end do
                end do
              end if
              do i_halo=1,n_halo
                if(target_frag == jfrag_halo(i_halo)) then
                  do io_x=1,n_basis(target_frag,ispin)
                  do io_y=1,n_basis(target_frag,ispin)
                    h_block(io_x,io_y,target_frag) = h_block(io_x,io_y,target_frag) &
                    & + 0.25d0 * (h(io_x,io_y,i_halo) + h(io_y,io_x,i_halo))
                  end do
                  end do
                end if
              end do
            end do
          end do
          if(dc%id_frag==0) coef_wf(:,:,ispin) = 0d0
          n_entry = 0
          do ifrag=1,dc%n_frag
            nbasis_diag = n_basis(ifrag,ispin)
            if(nbasis_diag <= 0) cycle
            allocate(h_local_diag(nbasis_diag,nbasis_diag))
            allocate(evec_local(nbasis_diag,nbasis_diag))
            allocate(eval_local(nbasis_diag))
            h_local_diag(:,:) = h_block(1:nbasis_diag,1:nbasis_diag,ifrag)
            if (any(.not. ieee_is_finite(h_local_diag(1:nbasis_diag,1:nbasis_diag)))) then
              if (dc%id_tot == 0) write(*,'(1x,a,i0,a,i0)') &
                "[FATAL] non-finite block Flux Hamiltonian before local diagonalization: ifrag=", &
                ifrag, " ispin=", ispin
              stop "DC-LCFO block Flux H contains non-finite values"
            end if
            call eigen_dsyev(h_local_diag, eval_local, evec_local)
            if (any(.not. ieee_is_finite(evec_local(1:nbasis_diag,1:nbasis_diag))) .or. &
                any(.not. ieee_is_finite(eval_local(1:nbasis_diag)))) then
              if (dc%id_tot == 0) write(*,'(1x,a,i0,a,i0)') &
                "[FATAL] non-finite block Flux eigenpair: ifrag=", ifrag, " ispin=", ispin
              stop "DC-LCFO block Flux eigensolver returned non-finite values"
            end if
            max_level = min(nbasis_diag, dc%nstate_tot)
            do level=1,max_level
              if(n_entry >= n) exit
              n_entry = n_entry + 1
              eval_list(n_entry) = eval_local(level)
              frag_list(n_entry) = ifrag
              level_list(n_entry) = level
              h_block(1:nbasis_diag,level,ifrag) = evec_local(1:nbasis_diag,level)
            end do
            deallocate(h_local_diag,evec_local,eval_local)
          end do
          do pass=1,max(0,n_entry-1)
            best = pass
            do level=pass+1,n_entry
              if(eval_list(level) < eval_list(best)) best = level
            end do
            if(best /= pass) then
              eps = eval_list(pass)
              eval_list(pass) = eval_list(best)
              eval_list(best) = eps
              tmp_i = frag_list(pass)
              frag_list(pass) = frag_list(best)
              frag_list(best) = tmp_i
              tmp_i = level_list(pass)
              level_list(pass) = level_list(best)
              level_list(best) = tmp_i
            end if
          end do
          do state_col=1,min(n_entry,dc%nstate_tot)
            ifrag = frag_list(state_col)
            level = level_list(state_col)
            if(ifrag < 1 .or. ifrag > dc%n_frag) cycle
            nbasis_diag = n_basis(ifrag,ispin)
            if(level < 1 .or. level > nbasis_diag) cycle
            esp_tot(state_col,ispin) = eval_list(state_col)
            if(ifrag==dc%i_frag .and. dc%id_frag==0) then
              coef_wf(1:nbasis_diag,state_col,ispin) = h_block(1:nbasis_diag,level,ifrag)
            end if
          end do
          if(dc%id_tot==0) call print_flux_gap_diagnostic("block", ispin, nocc_nelec, min(n_entry,dc%nstate_tot))
          deallocate(h_block,eval_list,frag_list,level_list)
        end if

        deallocate(h_div,h_ref_div,v_div,v_tmp1,v_tmp2,io_array,ifrag_array)
        deallocate(coef_state_norm,coef_state_norm_alt,repair_state)
        call eigen_free()
      end do ! ispin

      deallocate(h)
    end subroutine diag_eigenexa

    subroutine diagnose_sawf_hamiltonian_hermiticity( &
        h,fragment_index,basis_index,local_block,local_basis_count)
      implicit none
      real(8), intent(in) :: h(:,:)
      integer, intent(in) :: fragment_index(:),basis_index(:)
      real(8), intent(in) :: local_block(:,:)
      integer, intent(in) :: local_basis_count
      integer :: location(2),block_location(2),block_start,block_end
      real(8) :: scale,residual,block_difference

      scale=max(1.0d-300,maxval(abs(h)))
      location=maxloc(abs(h-transpose(h)))
      residual=abs(h(location(1),location(2))-h(location(2),location(1)))/scale
      write(*,'(1x,a,es13.5,2(a,i0),4(a,i0))') &
        '[DC-LCFO-SAWF-HERMITICITY] relative_max=',residual, &
        ' row=',location(1),' col=',location(2), &
        ' row_fragment=',fragment_index(location(1)),' row_basis=',basis_index(location(1)), &
        ' col_fragment=',fragment_index(location(2)),' col_basis=',basis_index(location(2))
      block_start=findloc(fragment_index,dc%i_frag,dim=1)
      block_end=block_start+local_basis_count-1
      block_location=maxloc(abs(h(block_start:block_end,block_start:block_end)- &
        local_block(1:local_basis_count,1:local_basis_count)))
      block_difference=maxval(abs(h(block_start:block_end,block_start:block_end)- &
        local_block(1:local_basis_count,1:local_basis_count)))
      write(*,'(1x,a,es13.5,2(a,i0))') &
        '[DC-LCFO-SAWF-HERMITICITY] local_block_difference=',block_difference, &
        ' row_basis=',block_location(1),' col_basis=',block_location(2)
    end subroutine diagnose_sawf_hamiltonian_hermiticity

    integer function occupied_index_from_input(ispin_in) result(nocc_input)
      use salmon_global, only: nelec_spin
      implicit none
      integer, intent(in) :: ispin_in

      if(nspin == 1) then
        nocc_input = int(0.5d0 * dc%elec_num_tot + 1.0d-12)
        if(abs(2.0d0 * dble(nocc_input) - dc%elec_num_tot) > 1.0d-8) nocc_input = nocc_input + 1
      else if(sum(nelec_spin(1:min(2,size(nelec_spin)))) > 0) then
        nocc_input = nelec_spin(ispin_in)
      else
        nocc_input = int(dc%elec_num_tot + 1.0d-12)
      end if
      nocc_input = max(1, min(nocc_input, dc%nstate_tot - 1))
    end function occupied_index_from_input

    subroutine print_flux_gap_diagnostic(label, ispin_in, nocc_input, nstate_available)
      implicit none
      character(*), intent(in) :: label
      integer, intent(in) :: ispin_in, nocc_input, nstate_available
      real(8) :: gap_ev

      if(nocc_input < 1 .or. nocc_input + 1 > nstate_available) then
        write(*,'(1x,a,a,a,i0,a,i0,a,i0,a)') &
          "[DC-LCFO-FLUX-GAP] mode=", trim(label), " ispin=", ispin_in, &
          " nocc_input=", nocc_input, " nstate_available=", nstate_available, " unavailable"
        return
      end if
      gap_ev = (esp_tot(nocc_input+1,ispin_in) - esp_tot(nocc_input,ispin_in)) * 27.211386245988d0
      write(*,'(1x,a,a,a,i0,a,i0,3(a,1pe13.5))') &
        "[DC-LCFO-FLUX-GAP] mode=", trim(label), " ispin=", ispin_in, &
        " nocc_input=", nocc_input, &
        " homo_eV=", esp_tot(nocc_input,ispin_in) * 27.211386245988d0, &
        " lumo_eV=", esp_tot(nocc_input+1,ispin_in) * 27.211386245988d0, &
        " gap_eV=", gap_ev
    end subroutine print_flux_gap_diagnostic
#endif

    logical function use_block_diag_hamiltonian_mode() result(enabled)
      use salmon_global, only: yn_dc_lcfo_block_diag_h
      implicit none
      character(16) :: env_value
      integer :: env_status

      enabled = (yn_dc_lcfo_block_diag_h == 'y')
      env_value = ''
      call get_environment_variable('SALMON_DG_BLOCK_DIAG_H', env_value, status=env_status)
      if(env_status == 0) then
        select case(trim(adjustl(env_value)))
        case('1','y','Y','yes','YES','true','TRUE','on','ON')
          enabled = .true.
        case default
          enabled = .false.
        end select
      end if
    end function use_block_diag_hamiltonian_mode

    logical function use_weak_volume_hamiltonian_mode() result(enabled)
      use salmon_global, only: yn_dc_lcfo_flux_weak_volume
      implicit none
      logical, save :: initialized = .false.
      logical, save :: enabled_save = .false.
      character(16) :: env_value
      integer :: env_status

      enabled_save = (yn_dc_lcfo_flux_weak_volume == 'y')
      if(.not. initialized) then
        env_value = ''
        call get_environment_variable('SALMON_DG_FLUX_WEAK_VOLUME', env_value, status=env_status)
        if(env_status == 0) then
          select case(trim(adjustl(env_value)))
          case('1','y','Y','yes','YES','true','TRUE','on','ON')
            enabled_save = .true.
          end select
        end if
        initialized = .true.
      end if
      enabled = (yn_dc_lcfo_flux_weak_volume == 'y') .or. enabled_save
    end function use_weak_volume_hamiltonian_mode

    logical function use_surface_self_hamiltonian_mode() result(enabled)
      implicit none
      logical, save :: initialized = .false.
      logical, save :: enabled_save = .true.
      character(16) :: env_value
      integer :: env_status

      if(.not. initialized) then
        env_value = ''
        call get_environment_variable('SALMON_DG_FLUX_SELF_SURFACE', env_value, status=env_status)
        if(env_status == 0) then
          select case(trim(adjustl(env_value)))
          case('0','n','N','no','NO','false','FALSE','off','OFF')
            enabled_save = .false.
          case('1','y','Y','yes','YES','true','TRUE','on','ON')
            enabled_save = .true.
          end select
        end if
        initialized = .true.
      end if
      enabled = enabled_save
    end function use_surface_self_hamiltonian_mode

    subroutine output
      use salmon_global, only: base_directory, sysname, unit_energy, yn_dc_lcfo_wannier, &
        yn_dc_lcfo_local_wannier, yn_dc_lcfo_wannier_cluster, wannier_cluster_size
      use filesystem, only: get_filehandle
      use inputoutput, only: uenergy_from_au
      use communication, only: comm_sync_all
      implicit none
      integer :: iunit,i_halo
      integer :: nxyz_domain(3)
      integer :: nxyz_box(3), nxyz_buffer_seed(3)
      integer :: lb_rwf(3), ub_rwf(3), io_lb, io_ub
      integer :: ibasis, raw_io, ibx, iby, ibz, sx, sy, sz
      character(1024) :: filename
      real(8) :: coef_val
      real(8), allocatable :: phi_box_local(:,:,:), phi_box_sum(:,:,:)

    ! total system data
      if(dc%id_tot==0 .and. yn_dc_lcfo_diag=='y') then
      ! eigen.data
        iunit = get_filehandle()
        filename = trim(dc%base_directory)//trim(sysname)//"_eigen.data" ! @ ./data_dcdft/total/
        open(iunit,file=filename)
        write(iunit,'("#esp: single-particle energies (eigen energies) calculated by DC-LCFO method")')
        write(iunit,'("#io: orbital index")')
        select case(unit_energy)
        case('au','a.u.')
          write(iunit,'("# 1:io, 2:esp[a.u.]")')
        case('ev','eV')
          write(iunit,'("# 1:io, 2:esp[eV]")')
        end select
        do ispin=1,nspin
          write(iunit,'("# spin=",1x,i5)') ispin
          do i=1,dc%nstate_tot
            write(iunit,'(1x,i5,e26.16e3)') i,esp_tot(i,ispin)*uenergy_from_au
          end do
        end do
        close(iunit)
      end if
      if(yn_dc_lcfo_wannier_cluster == 'y' .and. yn_dc_lcfo_diag == 'y') call write_wannier_cluster_partition_file()

    ! fragment data
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      if(dc%id_frag==0) then
      ! r-grid index
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_rg ! base_directory==./data_dcdft/fragments/dc%i_frag/
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) lg%num(1:3), dc%lg_tot%num(1:3)
        do n=1,3 ! x,y,z
          write(iunit) dc%jxyz_tot(1:lg%num(n),n)
        end do
        close(iunit)
      ! basis functions | lambda >
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_bf ! base_directory==./data_dcdft/fragments/dc%i_frag/
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) nxyz_domain(1:3),nspin,dc%nstate_frag
        write(iunit) n_basis(dc%i_frag,1:nspin) ! # of basis functions
        write(iunit) f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3) &
        & ,1:nspin,1:dc%nstate_frag) ! basis functions | lambda >
        close(iunit)
      end if

      ! buffered basis functions for DG surface Flux/Flow terms
        nxyz_buffer_seed(1:3) = dc%nxyz_buffer(1:3)
        nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer_seed(1:3)
        lb_rwf(1) = lbound(spsi%rwf, 1)
        lb_rwf(2) = lbound(spsi%rwf, 2)
        lb_rwf(3) = lbound(spsi%rwf, 3)
        ub_rwf(1) = ubound(spsi%rwf, 1)
        ub_rwf(2) = ubound(spsi%rwf, 2)
        ub_rwf(3) = ubound(spsi%rwf, 3)
        io_lb = lbound(spsi%rwf, 5)
        io_ub = ubound(spsi%rwf, 5)
        allocate(phi_box_local(nxyz_box(1), nxyz_box(2), nxyz_box(3)))
        allocate(phi_box_sum(nxyz_box(1), nxyz_box(2), nxyz_box(3)))
        if(dc%id_frag==0) then
          iunit = get_filehandle()
          filename = trim(base_directory)//binfile_bfb
          open(iunit,file=filename,form='unformatted',access='stream')
          write(iunit) basis_buffer_magic, basis_buffer_version
          write(iunit) nxyz_domain(1:3), nxyz_buffer_seed(1:3), nspin, dc%nstate_frag
          write(iunit) n_basis(dc%i_frag,1:nspin)
        end if
        do ispin=1,nspin
          do ibasis=1,dc%nstate_frag
            phi_box_local(:,:,:) = 0d0
            if(ibasis <= n_basis(dc%i_frag,ispin)) then
              do raw_io=max(1,io_lb),min(dc%nstate_frag,io_ub)
                coef_val = basis_transform(raw_io,ibasis,ispin)
                if(abs(coef_val) <= 0d0) cycle
                do ibz=1,nxyz_box(3)
                  sz = dc_buffer_box_to_local_index(ibz,nxyz_domain(3),nxyz_buffer_seed(3))
                  if(sz < lb_rwf(3) .or. sz > ub_rwf(3)) cycle
                  do iby=1,nxyz_box(2)
                    sy = dc_buffer_box_to_local_index(iby,nxyz_domain(2),nxyz_buffer_seed(2))
                    if(sy < lb_rwf(2) .or. sy > ub_rwf(2)) cycle
                    do ibx=1,nxyz_box(1)
                      sx = dc_buffer_box_to_local_index(ibx,nxyz_domain(1),nxyz_buffer_seed(1))
                      if(sx < lb_rwf(1) .or. sx > ub_rwf(1)) cycle
                      phi_box_local(ibx,iby,ibz) = phi_box_local(ibx,iby,ibz) &
                      & + coef_val * spsi%rwf(sx,sy,sz,ispin,raw_io,1,1)
                    end do
                  end do
                end do
              end do
            end if
            call comm_summation(phi_box_local,phi_box_sum,product(nxyz_box),dc%icomm_frag)
            if(dc%id_frag==0) write(iunit) phi_box_sum(1:nxyz_box(1),1:nxyz_box(2),1:nxyz_box(3))
          end do
        end do
        if(dc%id_frag==0) close(iunit)
        deallocate(phi_box_local,phi_box_sum)
      if(dc%id_frag==0) then
      ! local hamiltonian matrix
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_hl
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) mat_H_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_H_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)
      ! diagnostic decomposition of the Flux Hamiltonian.
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_hc
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) mat_H_volume_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) mat_H_surface_self(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_H_surface_cross(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)
      ! diagnostic decomposition with the kinetic volume replaced by the
      ! weak-form fragment volume term.
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_hcw
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) mat_H_volume_weak_local(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) mat_H_surface_self(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_H_surface_cross(1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)
      ! local velocity matrix dH/dA at A=0, built from the same Flux-LCFO basis
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_vl
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) mat_V_local(1:3,1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_V_local(1:3,1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)
      ! local DG boundary position correction xi_flux for length-gauge RT
        iunit = get_filehandle()
        filename = trim(base_directory)//binfile_xfl
        open(iunit,file=filename,form='unformatted',access='stream')
        write(iunit) 0d0 * mat_V_local(1:3,1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        write(iunit) n_halo
        do i_halo=1,n_halo
          write(iunit) halo(i_halo)%mat_Xi_flux_local(1:3,1:dc%nstate_frag,1:dc%nstate_frag,1:nspin)
        end do
        close(iunit)
        if(yn_dc_lcfo_diag=='y') then
        ! coefficients of the wavefunctions
          iunit = get_filehandle()
          filename = trim(base_directory)//binfile_wf
          open(iunit,file=filename,form='unformatted',access='stream')
          write(iunit) dc%n_frag, nspin, dc%nstate_frag, dc%nstate_tot
          write(iunit) n_mat(1:nspin)
          write(iunit) n_basis(1:dc%n_frag,1:nspin)
          write(iunit) index_basis(1:dc%nstate_frag,1:dc%n_frag,1:nspin)
          write(iunit) coef_wf(1:dc%nstate_frag,1:dc%nstate_tot,1:nspin)
          write(iunit) esp_tot(1:dc%nstate_tot,1:nspin)
          close(iunit)
          iunit = get_filehandle()
          filename = trim(base_directory)//binfile_wf_wannier_seed
          open(iunit,file=filename,form='unformatted',access='stream',status='replace')
          write(iunit) dc%n_frag, nspin, dc%nstate_frag, dc%nstate_tot
          write(iunit) n_mat(1:nspin)
          write(iunit) n_basis(1:dc%n_frag,1:nspin)
          write(iunit) index_basis(1:dc%nstate_frag,1:dc%n_frag,1:nspin)
          write(iunit) coef_wf(1:dc%nstate_frag,1:dc%nstate_tot,1:nspin)
          write(iunit) esp_tot(1:dc%nstate_tot,1:nspin)
          close(iunit)
          filename = trim(base_directory)//'wavefunctions_wannier_seed.provenance'
          open(newunit=iunit,file=filename,status='replace',action='write')
          write(iunit,'(i0,1x,i0)') sawf_seed_provenance_magic,sawf_seed_provenance_version
          write(iunit,'(a)') trim(current_sawf_seed_token)
          write(iunit,'(5(i0,1x))') dc%n_frag,dc%i_frag,nspin,dc%nstate_frag,dc%nstate_tot
          write(iunit,'(*(i0,1x))') n_basis(dc%i_frag,1:nspin)
          close(iunit)
        end if
      end if

      call comm_sync_all(dc%icomm_tot)
      if(yn_dc_lcfo_wannier == 'y' .and. yn_dc_lcfo_diag == 'y') call write_wannier_seed_files()
      if(yn_dc_lcfo_local_wannier == 'y' .and. yn_dc_lcfo_diag == 'y') &
        call write_local_wannier_seed()
      if(yn_dc_lcfo_wannier == 'y' .and. yn_dc_lcfo_wannier_cluster == 'y' .and. &
         yn_dc_lcfo_local_wannier /= 'y' .and. yn_dc_lcfo_diag == 'y') then
        if(all(wannier_cluster_size(1:3) == 1)) then
          call write_wannier90_global_bpw_seed()
        else if(dc%id_tot == 0) then
          write(*,'(1x,a,3(i0,1x),a)') &
            "[DC-LCFO-W90-BPW] skip fragment BPW: cluster_size=", &
            wannier_cluster_size(1:3), " requires cluster/global BPW export."
        end if
      end if

    end subroutine output

    subroutine write_wannier_cluster_partition_file()
      use salmon_global, only: num_fragment, num_wannier_cluster, wannier_cluster_size
      use filesystem, only: get_filehandle
      implicit none
      integer :: iunit, ifrag, cid, ix_frag, iy_frag, iz_frag, rem
      integer :: icx, icy, icz, num_cluster(3), ncluster_tot
      integer, allocatable :: frag_cluster_id(:), cluster_frag_range(:,:)
      character(256) :: filename

      if(dc%id_tot /= 0) return
      num_cluster(1:3) = num_fragment(1:3) / wannier_cluster_size(1:3)
      ncluster_tot = product(num_cluster(1:3))
      allocate(frag_cluster_id(dc%n_frag), cluster_frag_range(6,ncluster_tot))
      cluster_frag_range = 0

      do icx=1,num_cluster(1)
      do icy=1,num_cluster(2)
      do icz=1,num_cluster(3)
        cid = ((icx - 1) * num_cluster(2) + (icy - 1)) * num_cluster(3) + icz
        cluster_frag_range(1,cid) = (icx - 1) * wannier_cluster_size(1) + 1
        cluster_frag_range(2,cid) = icx * wannier_cluster_size(1)
        cluster_frag_range(3,cid) = (icy - 1) * wannier_cluster_size(2) + 1
        cluster_frag_range(4,cid) = icy * wannier_cluster_size(2)
        cluster_frag_range(5,cid) = (icz - 1) * wannier_cluster_size(3) + 1
        cluster_frag_range(6,cid) = icz * wannier_cluster_size(3)
      end do
      end do
      end do

      do ifrag=1,dc%n_frag
        ix_frag = (ifrag - 1) / max(1, num_fragment(2) * num_fragment(3)) + 1
        rem = modulo(ifrag - 1, max(1, num_fragment(2) * num_fragment(3)))
        iy_frag = rem / max(1, num_fragment(3)) + 1
        iz_frag = modulo(rem, max(1, num_fragment(3))) + 1
        icx = (ix_frag - 1) / wannier_cluster_size(1) + 1
        icy = (iy_frag - 1) / wannier_cluster_size(2) + 1
        icz = (iz_frag - 1) / wannier_cluster_size(3) + 1
        frag_cluster_id(ifrag) = ((icx - 1) * num_cluster(2) + (icy - 1)) * num_cluster(3) + icz
      end do

      filename = trim(dc%base_directory)//binfile_wcl
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='replace')
      write(iunit) wannier_cluster_magic, wannier_cluster_version
      write(iunit) num_fragment(1:3), wannier_cluster_size(1:3), num_cluster(1:3)
      write(iunit) dc%n_frag, ncluster_tot
      write(iunit) frag_cluster_id(1:dc%n_frag)
      write(iunit) cluster_frag_range(1:6,1:ncluster_tot)
      close(iunit)
      write(*,'(1x,a,i0,a,a)') "[DC-LCFO-WANNIER-CLUSTER] wrote ", ncluster_tot, &
        " clusters to ", trim(filename)
      call write_wannier_cluster_manifest_file(num_cluster, ncluster_tot, frag_cluster_id, cluster_frag_range)
      deallocate(frag_cluster_id, cluster_frag_range)
    end subroutine write_wannier_cluster_partition_file

    subroutine write_wannier_cluster_manifest_file(num_cluster, ncluster_tot, frag_cluster_id, cluster_frag_range)
      use salmon_global, only: num_fragment, num_wannier_cluster, wannier_cluster_size
      use filesystem, only: get_filehandle
      implicit none
      integer, intent(in) :: num_cluster(3), ncluster_tot
      integer, intent(in) :: frag_cluster_id(dc%n_frag), cluster_frag_range(6,ncluster_tot)
      integer :: iunit, cid, ifrag, nmember
      character(256) :: filename

      filename = trim(dc%base_directory)//"wannier_cluster_partition.txt"
      iunit = get_filehandle()
      open(iunit,file=filename,status='replace')
      write(iunit,'(a)') "# SALMON DC-LCFO Wannier cluster partition"
      write(iunit,'(a,3(1x,i0))') "# num_fragment", num_fragment(1:3)
      write(iunit,'(a,3(1x,i0))') "# input_num_wannier_cluster", num_wannier_cluster(1:3)
      write(iunit,'(a,3(1x,i0))') "# wannier_cluster_size", wannier_cluster_size(1:3)
      write(iunit,'(a,3(1x,i0))') "# effective_num_wannier_cluster", num_cluster(1:3)
      write(iunit,'(a)') "# columns: cluster_id owner_fragment xlo xhi ylo yhi zlo zhi nmember members..."
      do cid=1,ncluster_tot
        nmember = count(frag_cluster_id(1:dc%n_frag) == cid)
        write(iunit,'(i8,1x,i8,6(1x,i8),1x,i8)', advance='no') cid, &
          first_fragment_in_cluster(cid, frag_cluster_id), cluster_frag_range(1:6,cid), nmember
        do ifrag=1,dc%n_frag
          if(frag_cluster_id(ifrag) == cid) write(iunit,'(1x,i8)', advance='no') ifrag
        end do
        write(iunit,*)
      end do
      close(iunit)
      write(*,'(1x,2a)') "[DC-LCFO-WANNIER-CLUSTER] manifest: ", trim(filename)
    end subroutine write_wannier_cluster_manifest_file

    integer function first_fragment_in_cluster(cid, frag_cluster_id) result(ifrag_first)
      implicit none
      integer, intent(in) :: cid, frag_cluster_id(dc%n_frag)
      integer :: ifrag

      ifrag_first = 1
      do ifrag=1,dc%n_frag
        if(frag_cluster_id(ifrag) == cid) then
          ifrag_first = ifrag
          return
        end if
      end do
    end function first_fragment_in_cluster

    subroutine write_local_wannier_seed()
      use eigen_subdiag_sub, only: eigen_dsyev
      use filesystem, only: get_filehandle
      use salmon_global, only: izatom, sysname, wannier_projection, wannier_projection_width, &
        lambda_cut, base_directory, yn_dc_lcfo_wannier, yn_dc_lcfo_wannier_pw, &
        yn_dc_lcfo_wannier_cluster, wannier_cluster_size, yn_dc_lcfo_wannier_symmetry_gauge
      implicit none
      integer :: nproj, nproj_seed, nkeep, nkeep_legacy, nbasis, iunit, ip, jp, iw, jw, io, jo, axis
      integer :: aa_wann_source
      integer :: ix, iy, iz, ixg, iyg, izg, ispin_local
      integer :: ibx, iby, ibz, sx, sy, sz, ispin_read, ibasis_read
      integer :: proj_l, proj_m
      integer :: nxyz_domain(3), nxyz_buffer_seed(3), nxyz_box(3)
      integer :: nxyz_domain_file(3), nxyz_buffer_file(3), nspin_file, nstate_frag_file
      integer :: magic_file, version_file
      integer, allocatable :: n_basis_file(:)
      integer, allocatable :: proj_atom(:), proj_hybrid(:), keep_index(:), keep_index_legacy(:)
      real(8), allocatable :: bproj(:,:), sw(:,:), uw(:,:), lambda_w(:)
      real(8), allocatable :: bond_center_bohr(:,:)
      real(8), allocatable :: wcoef(:,:), r_basis(:,:,:), r_wann(:,:,:), tmp(:,:), wcenter(:,:)
      real(8), allocatable :: sw_legacy(:,:), uw_legacy(:,:), lambda_legacy(:)
      real(8), allocatable :: wcoef_legacy(:,:), r_wann_legacy(:,:,:), tmp_legacy(:,:), wcenter_legacy(:,:)
      real(8), allocatable :: h_seed(:,:), v_seed(:,:,:)
      real(8), allocatable :: h_wann(:,:), v_wann(:,:,:), aa_wann(:,:,:), spread_est(:), tail_est(:)
      real(8), allocatable :: phi_box(:,:,:,:), phi_tmp(:,:,:)
      real(8), allocatable :: psi_w(:,:,:,:), rho_w_sum(:), cos_sum(:,:), sin_sum(:,:)
      real(8) :: x, y, z, gval, box_origin(3), cell_length(3)
      real(8) :: theta, pair_center, frag_center_axis, rel_coord, norm_w, pi_twice
      character(256) :: filename
      logical :: use_pseudo_projection, use_bond_projection

      if(nspin /= 1) stop "DC-LCFO local Wannier export: spin-polarized mode is not implemented."
      use_pseudo_projection = is_pseudo_channel_projection(trim(wannier_projection))
      use_bond_projection = is_bond_center_projection(trim(wannier_projection))
      if(.not. is_sp3_projection(trim(wannier_projection)) .and. .not. use_pseudo_projection .and. &
         .not. use_bond_projection) &
        stop "DC-LCFO local Wannier export: supported projections are C:sp3, Si:sp3, pseudo_channels, and bond_centers."
      if(yn_dc_lcfo_wannier_pw == 'y') &
        stop "DC-LCFO local Wannier PW augmentation is not connected yet."
      if(yn_dc_lcfo_wannier_cluster == 'y' .and. any(wannier_cluster_size(1:3) /= 1)) then
        if(dc%id_tot == 0) then
          write(*,'(1x,a)') "[DC-LCFO-LOCAL-WANNIER] cluster partition is enabled."
          write(*,'(1x,a)') "  Fragment-local BPW generation is disabled because it would ignore the requested larger Wannier space."
          write(*,'(1x,a)') "  Use the global Wannier90 export path for now, or generate cluster BPW seeds in a follow-up step."
        end if
        stop "DC-LCFO local Wannier cluster BPW generation is not implemented yet."
      end if
      if(dc%id_frag /= 0) return

      ispin_local = 1
      nbasis = n_basis(dc%i_frag,ispin_local)
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      nxyz_buffer_seed(1:3) = dc%nxyz_buffer(1:3)
      nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer_seed(1:3)
      do axis=1,3
        cell_length(axis) = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
          + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
        box_origin(axis) = dc%lg_tot%coordinate(dc%jxyz_tot(1,axis),axis) &
          - dble(nxyz_buffer_seed(axis)) * system%hgs(axis)
      end do
      pi_twice = 8d0 * atan(1d0)
      if(use_bond_projection) then
        call build_local_bond_center_projection_map(nproj_seed, bond_center_bohr)
        if(nproj_seed <= 0) stop "DC-LCFO local Wannier export: no local bond-center projections in fragment core."
        call diagnose_local_bond_center_orbit_closure(dc, nproj_seed, bond_center_bohr)
      else if(use_pseudo_projection) then
        nproj_seed = count_local_pseudo_channel_ao_candidates(nxyz_domain)
        if(nproj_seed <= 0) stop "DC-LCFO local Wannier export: no local pseudo-channel projections in fragment core."
      else
        nproj_seed = count_local_c_sp3_projections(nxyz_domain)
        if(nproj_seed <= 0) stop "DC-LCFO local Wannier export: no local sp3 projections in fragment core."
      end if
      nproj = nproj_seed + nbasis

      allocate(proj_atom(nproj), proj_hybrid(nproj))
      proj_atom = 0
      proj_hybrid = 0
      if(use_bond_projection) then
        proj_atom(1:nproj_seed) = 0
        proj_hybrid(1:nproj_seed) = 100
      else if(use_pseudo_projection) then
        call build_local_pseudo_channel_ao_candidate_map(nxyz_domain, proj_atom, proj_hybrid)
      else
        call build_local_c_sp3_projection_map(nxyz_domain, proj_atom, proj_hybrid)
      end if
      allocate(bproj(nbasis,nproj), sw(nproj,nproj), uw(nproj,nproj), lambda_w(nproj))
      allocate(wcoef(nbasis,nproj), keep_index(nproj))
      allocate(r_basis(3,nbasis,nbasis), r_wann(3,nproj,nproj), tmp(nbasis,nproj), wcenter(3,nproj))
      allocate(sw_legacy(nproj_seed,nproj_seed), uw_legacy(nproj_seed,nproj_seed), lambda_legacy(nproj_seed))
      allocate(wcoef_legacy(nbasis,nproj_seed), keep_index_legacy(nproj_seed))
      allocate(r_wann_legacy(3,nproj_seed,nproj_seed), tmp_legacy(nbasis,nproj_seed), wcenter_legacy(3,nproj_seed))
      allocate(h_seed(nbasis,nbasis), v_seed(3,nbasis,nbasis))
      allocate(h_wann(nproj,nproj), v_wann(3,nproj,nproj), aa_wann(3,nproj,nproj), spread_est(nproj), tail_est(nproj))
      bproj = 0d0
      sw = 0d0
      sw_legacy = 0d0
      r_basis = 0d0
      wcoef_legacy = 0d0
      r_wann_legacy = 0d0
      wcenter_legacy = 0d0
      h_wann = 0d0
      v_wann = 0d0
      aa_wann = 0d0
      aa_wann_source = 0
      h_seed = 0d0
      v_seed = 0d0
      spread_est = 0d0
      tail_est = 0d0
      allocate(phi_box(nxyz_box(1),nxyz_box(2),nxyz_box(3),nbasis))
      allocate(phi_tmp(nxyz_box(1),nxyz_box(2),nxyz_box(3)))
      allocate(n_basis_file(nspin))
      phi_box = 0d0

      filename = trim(base_directory)//binfile_bfb
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='old')
      read(iunit) magic_file, version_file
      if(magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) &
        stop "DC-LCFO buffer-periodic Wannier export: invalid buffered basis header."
      read(iunit) nxyz_domain_file(1:3), nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
      if(any(nxyz_domain_file(1:3) /= nxyz_domain(1:3)) .or. &
         any(nxyz_buffer_file(1:3) /= nxyz_buffer_seed(1:3)) .or. &
         nspin_file /= nspin .or. nstate_frag_file /= dc%nstate_frag) &
        stop "DC-LCFO buffer-periodic Wannier export: buffered basis metadata mismatch."
      read(iunit) n_basis_file(1:nspin)
      do ispin_read=1,nspin_file
        do ibasis_read=1,nstate_frag_file
          read(iunit) phi_tmp(1:nxyz_box(1),1:nxyz_box(2),1:nxyz_box(3))
          if(ispin_read == ispin_local .and. ibasis_read <= nbasis) &
            phi_box(1:nxyz_box(1),1:nxyz_box(2),1:nxyz_box(3),ibasis_read) = &
              phi_tmp(1:nxyz_box(1),1:nxyz_box(2),1:nxyz_box(3))
        end do
      end do
      close(iunit)

!$omp parallel do collapse(2) private(ip,io,ibz,iby,ibx,z,y,x,gval,proj_l,proj_m) schedule(static)
      do ip=1,nproj_seed
        do io=1,nbasis
          do ibz=1,nxyz_box(3)
            z = box_origin(3) + dble(ibz - 1) * system%hgs(3)
            do iby=1,nxyz_box(2)
              y = box_origin(2) + dble(iby - 1) * system%hgs(2)
              do ibx=1,nxyz_box(1)
                x = box_origin(1) + dble(ibx - 1) * system%hgs(1)
                if(use_bond_projection) then
                  gval = bond_center_projection_value_local_periodic(x, y, z, &
                    bond_center_bohr(1:3,ip), wannier_projection_width, cell_length)
                else if(use_pseudo_projection) then
                  proj_l = proj_hybrid(ip) / 10
                  proj_m = mod(proj_hybrid(ip), 10)
                  gval = pseudo_channel_projection_value_local_periodic(x, y, z, proj_atom(ip), &
                    proj_l, proj_m, wannier_projection_width, cell_length)
                else
                  gval = c_sp3_projection_value_local_periodic(x, y, z, proj_atom(ip), &
                    proj_hybrid(ip), wannier_projection_width, cell_length)
                end if
                if(abs(gval) <= 0d0) cycle
                bproj(io,ip) = bproj(io,ip) &
                  + phi_box(ibx,iby,ibz,io) * gval * hvol
              end do
            end do
          end do
        end do
      end do
!$omp end parallel do

      sw_legacy(1:nproj_seed,1:nproj_seed) = &
        matmul(transpose(bproj(1:nbasis,1:nproj_seed)), bproj(1:nbasis,1:nproj_seed))
      call eigen_dsyev(sw_legacy, lambda_legacy, uw_legacy)
      wcoef_legacy = 0d0
      nkeep_legacy = 0
      do ip=nproj_seed,1,-1
        if(lambda_legacy(ip) <= lambda_cut) cycle
        nkeep_legacy = nkeep_legacy + 1
        keep_index_legacy(nkeep_legacy) = ip
        do io=1,nbasis
          do jp=1,nproj_seed
            wcoef_legacy(io,nkeep_legacy) = wcoef_legacy(io,nkeep_legacy) &
              + bproj(io,jp) * uw_legacy(jp,ip) / sqrt(lambda_legacy(ip))
          end do
        end do
      end do
      if(nkeep_legacy <= 0) stop "DC-LCFO local Wannier export: all local projection directions were removed."
      if(use_bond_projection .and. yn_dc_lcfo_wannier_symmetry_gauge == 'y') then
        call apply_local_bond_center_symmetry_gauge(nbasis, nproj_seed, lambda_legacy, uw_legacy, &
          bproj(1:nbasis,1:nproj_seed), wcoef_legacy, keep_index_legacy, nkeep_legacy)
      end if

      do io=1,nbasis
        bproj(io,nproj_seed+io) = 1d0
      end do

      sw(1:nproj,1:nproj) = matmul(transpose(bproj(1:nbasis,1:nproj)), bproj(1:nbasis,1:nproj))
      call eigen_dsyev(sw, lambda_w, uw)
      wcoef = 0d0
      nkeep = 0
      do ip=nproj,1,-1
        if(lambda_w(ip) <= lambda_cut) cycle
        nkeep = nkeep + 1
        keep_index(nkeep) = ip
        do io=1,nbasis
          do jp=1,nproj
            wcoef(io,nkeep) = wcoef(io,nkeep) &
              + bproj(io,jp) * uw(jp,ip) / sqrt(lambda_w(ip))
          end do
        end do
      end do
      if(nkeep <= 0) stop "DC-LCFO local Wannier export: all local projection directions were removed."

!$omp parallel do collapse(2) private(jo,io,ibz,iby,ibx,z,y,x) schedule(static)
      do jo=1,nbasis
        do io=1,nbasis
          do ibz=1,nxyz_box(3)
            z = box_origin(3) + dble(ibz - 1) * system%hgs(3)
            do iby=1,nxyz_box(2)
              y = box_origin(2) + dble(iby - 1) * system%hgs(2)
              do ibx=1,nxyz_box(1)
                x = box_origin(1) + dble(ibx - 1) * system%hgs(1)
                r_basis(1,io,jo) = r_basis(1,io,jo) &
                  + phi_box(ibx,iby,ibz,io) * x &
                  * phi_box(ibx,iby,ibz,jo) * hvol
                r_basis(2,io,jo) = r_basis(2,io,jo) &
                  + phi_box(ibx,iby,ibz,io) * y &
                  * phi_box(ibx,iby,ibz,jo) * hvol
                r_basis(3,io,jo) = r_basis(3,io,jo) &
                  + phi_box(ibx,iby,ibz,io) * z &
                  * phi_box(ibx,iby,ibz,jo) * hvol
              end do
            end do
          end do
        end do
      end do
!$omp end parallel do

      allocate(psi_w(nxyz_box(1),nxyz_box(2),nxyz_box(3),nkeep))
      allocate(rho_w_sum(nkeep), cos_sum(3,nkeep), sin_sum(3,nkeep))
      psi_w = 0d0
      rho_w_sum = 0d0
      cos_sum = 0d0
      sin_sum = 0d0
!$omp parallel do collapse(4) private(iw,io,ibz,iby,ibx) schedule(static)
      do iw=1,nkeep
        do ibz=1,nxyz_box(3)
          do iby=1,nxyz_box(2)
            do ibx=1,nxyz_box(1)
              do io=1,nbasis
                psi_w(ibx,iby,ibz,iw) = psi_w(ibx,iby,ibz,iw) &
                  + phi_box(ibx,iby,ibz,io) * wcoef(io,iw)
              end do
            end do
          end do
        end do
      end do
!$omp end parallel do

      do iw=1,nkeep
        do ibz=1,nxyz_box(3)
          z = box_origin(3) + dble(ibz - 1) * system%hgs(3)
          do iby=1,nxyz_box(2)
            y = box_origin(2) + dble(iby - 1) * system%hgs(2)
            do ibx=1,nxyz_box(1)
              x = box_origin(1) + dble(ibx - 1) * system%hgs(1)
              gval = psi_w(ibx,iby,ibz,iw) * psi_w(ibx,iby,ibz,iw) * hvol
              rho_w_sum(iw) = rho_w_sum(iw) + gval
              theta = pi_twice * modulo(x, cell_length(1)) / cell_length(1)
              cos_sum(1,iw) = cos_sum(1,iw) + cos(theta) * gval
              sin_sum(1,iw) = sin_sum(1,iw) + sin(theta) * gval
              theta = pi_twice * modulo(y, cell_length(2)) / cell_length(2)
              cos_sum(2,iw) = cos_sum(2,iw) + cos(theta) * gval
              sin_sum(2,iw) = sin_sum(2,iw) + sin(theta) * gval
              theta = pi_twice * modulo(z, cell_length(3)) / cell_length(3)
              cos_sum(3,iw) = cos_sum(3,iw) + cos(theta) * gval
              sin_sum(3,iw) = sin_sum(3,iw) + sin(theta) * gval
            end do
          end do
        end do
        do axis=1,3
          if (rho_w_sum(iw) > 0d0) then
            wcenter(axis,iw) = modulo(atan2(sin_sum(axis,iw), cos_sum(axis,iw)) / pi_twice, 1d0) &
              * cell_length(axis)
          else
            wcenter(axis,iw) = 0d0
          end if
        end do
      end do

      r_wann = 0d0
      do axis=1,3
        frag_center_axis = 0.5d0 * (dc%lg_tot%coordinate(dc%jxyz_tot(1,axis),axis) + &
          dc%lg_tot%coordinate(dc%jxyz_tot(1,axis) + nxyz_domain(axis) - 1,axis))
        do jw=1,nkeep
          do iw=1,nkeep
            pair_center = 0.5d0 * (wcenter(axis,iw) + wcenter(axis,jw))
            if (abs(wcenter(axis,iw) - wcenter(axis,jw)) > 0.5d0 * cell_length(axis)) &
              pair_center = modulo(pair_center + 0.5d0 * cell_length(axis), cell_length(axis))
            do ibz=1,nxyz_box(3)
              z = box_origin(3) + dble(ibz - 1) * system%hgs(3)
              do iby=1,nxyz_box(2)
                y = box_origin(2) + dble(iby - 1) * system%hgs(2)
                do ibx=1,nxyz_box(1)
                  x = box_origin(1) + dble(ibx - 1) * system%hgs(1)
                  select case(axis)
                  case(1)
                    rel_coord = periodic_delta_import(x - pair_center, cell_length(axis))
                  case(2)
                    rel_coord = periodic_delta_import(y - pair_center, cell_length(axis))
                  case default
                    rel_coord = periodic_delta_import(z - pair_center, cell_length(axis))
                  end select
                  r_wann(axis,iw,jw) = r_wann(axis,iw,jw) + &
                    psi_w(ibx,iby,ibz,iw) * rel_coord * psi_w(ibx,iby,ibz,jw) * hvol
                end do
              end do
            end do
            norm_w = sqrt(max(rho_w_sum(iw) * rho_w_sum(jw), 1d-60))
            r_wann(axis,iw,jw) = r_wann(axis,iw,jw) / norm_w
            if (iw == jw) then
              r_wann(axis,iw,iw) = r_wann(axis,iw,iw) + &
                periodic_delta_import(wcenter(axis,iw) - frag_center_axis, cell_length(axis))
            end if
          end do
        end do
      end do

      if(yn_dc_lcfo_wannier == 'y') then
        call project_wannier90_aa_to_bpw(nbasis, nkeep, wcoef, aa_wann, aa_wann_source)
      end if

      do axis=1,3
        tmp_legacy(1:nbasis,1:nkeep_legacy) = &
          matmul(r_basis(axis,1:nbasis,1:nbasis), wcoef_legacy(1:nbasis,1:nkeep_legacy))
        r_wann_legacy(axis,1:nkeep_legacy,1:nkeep_legacy) = &
          matmul(transpose(wcoef_legacy(1:nbasis,1:nkeep_legacy)), tmp_legacy(1:nbasis,1:nkeep_legacy))
      end do
      wcenter_legacy(1:3,1:nkeep_legacy) = 0d0
      do iw=1,nkeep_legacy
        wcenter_legacy(1:3,iw) = r_wann_legacy(1:3,iw,iw)
      end do
      if(use_bond_projection .and. yn_dc_lcfo_wannier_symmetry_gauge == 'y' .and. &
         allocated(bond_center_bohr)) then
        do iw=1,nkeep_legacy
          if(keep_index_legacy(iw) < 1 .or. keep_index_legacy(iw) > nproj_seed) cycle
          wcenter_legacy(1:3,iw) = bond_center_bohr(1:3,keep_index_legacy(iw))
          r_wann_legacy(1:3,iw,iw) = wcenter_legacy(1:3,iw)
        end do
        write(*,'(1x,a,i0,a,i0)') &
          "[DC-LCFO-LOCAL-WANNIER-SYM-GAUGE] center override=bond_center fragment=", &
          dc%i_frag, " keep=", nkeep_legacy
      end if
      call diagnose_local_wannier_center_orbit_closure(dc, nkeep_legacy, wcenter_legacy, 'local_wcenter')
      call diagnose_local_wannier_center_orbit_closure(dc, nkeep, wcenter, 'bpw_wcenter')

      h_seed(1:nbasis,1:nbasis) = 0d0
      v_seed(1:3,1:nbasis,1:nbasis) = 0d0
      do io=1,nbasis
        do jo=1,nbasis
          h_seed(io,jo) = 0.5d0 * (mat_H_local(io,jo,ispin_local) + mat_H_local(jo,io,ispin_local))
          v_seed(1:3,io,jo) = mat_V_local(1:3,io,jo,ispin_local)
          do i_halo=1,n_halo
            h_seed(io,jo) = h_seed(io,jo) + 0.25d0 * &
              (halo(i_halo)%mat_H_local(jo,io,ispin_local) + halo(i_halo)%mat_H_local(io,jo,ispin_local))
            v_seed(1:3,io,jo) = v_seed(1:3,io,jo) + 0.5d0 * &
              (halo(i_halo)%mat_V_local(1:3,jo,io,ispin_local) + &
               halo(i_halo)%mat_V_local(1:3,io,jo,ispin_local))
          end do
        end do
      end do
      do axis=1,3
        do io=1,nbasis
          v_seed(axis,io,io) = 0d0
          do jo=io+1,nbasis
            x = 0.5d0 * (v_seed(axis,io,jo) - v_seed(axis,jo,io))
            v_seed(axis,io,jo) = x
            v_seed(axis,jo,io) = -x
          end do
        end do
      end do

      tmp(1:nbasis,1:nkeep) = matmul(h_seed(1:nbasis,1:nbasis), wcoef(1:nbasis,1:nkeep))
      h_wann(1:nkeep,1:nkeep) = matmul(transpose(wcoef(1:nbasis,1:nkeep)), tmp(1:nbasis,1:nkeep))
      do axis=1,3
        tmp(1:nbasis,1:nkeep) = matmul(v_seed(axis,1:nbasis,1:nbasis), wcoef(1:nbasis,1:nkeep))
        v_wann(axis,1:nkeep,1:nkeep) = matmul(transpose(wcoef(1:nbasis,1:nkeep)), tmp(1:nbasis,1:nkeep))
      end do

      filename = trim(base_directory)//binfile_lw
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream')
      write(iunit) local_wannier_magic, local_wannier_version
      write(iunit) nxyz_domain(1:3), nspin, nbasis, nproj_seed, nkeep_legacy
      write(iunit) proj_atom(1:nproj_seed), proj_hybrid(1:nproj_seed)
      write(iunit) lambda_legacy(1:nproj_seed), keep_index_legacy(1:nkeep_legacy)
      write(iunit) wcoef_legacy(1:nbasis,1:nkeep_legacy)
      write(iunit) r_wann_legacy(1:3,1:nkeep_legacy,1:nkeep_legacy)
      write(iunit) wcenter_legacy(1:3,1:nkeep_legacy)
      close(iunit)
      if(use_bond_projection) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[DC-LCFO-LOCAL-WANNIER] fragment=", &
          dc%i_frag, " bond_centers=", nproj_seed, " keep=", nkeep_legacy
      else if(use_pseudo_projection) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[DC-LCFO-LOCAL-WANNIER] fragment=", &
          dc%i_frag, " pseudo_channels=", nproj_seed, " candidates=", nproj_seed, " keep=", nkeep_legacy
      else
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[DC-LCFO-LOCAL-WANNIER] fragment=", &
          dc%i_frag, " sp3=", nproj_seed, " candidates=", nproj_seed, " keep=", nkeep_legacy
      end if

      filename = trim(base_directory)//binfile_bpw
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream')
      write(iunit) buffer_periodic_wannier_magic, buffer_periodic_wannier_version
      write(iunit) dc%i_frag, nxyz_domain(1:3), nxyz_buffer_seed(1:3), nxyz_box(1:3)
      write(iunit) nspin, nbasis, nproj, nkeep
      write(iunit) proj_atom(1:nproj), proj_hybrid(1:nproj)
      write(iunit) lambda_w(1:nproj), keep_index(1:nkeep)
      write(iunit) spread_est(1:nkeep), tail_est(1:nkeep)
      write(iunit) wcoef(1:nbasis,1:nkeep)
      write(iunit) r_wann(1:3,1:nkeep,1:nkeep)
      write(iunit) wcenter(1:3,1:nkeep)
      write(iunit) h_wann(1:nkeep,1:nkeep)
      write(iunit) v_wann(1:3,1:nkeep,1:nkeep)
      write(iunit) aa_wann(1:3,1:nkeep,1:nkeep)
      write(iunit) aa_wann_source
      close(iunit)
      write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') "[DC-LCFO-BUFFER-WANNIER] fragment=", &
        dc%i_frag, " candidates=", nproj, " keep=", nkeep, &
        " lambda_min_keep=", minval(lambda_w(keep_index(1:nkeep)))

      deallocate(proj_atom, proj_hybrid, bproj, sw, uw, lambda_w)
      deallocate(wcoef, keep_index, r_basis, r_wann, tmp, wcenter)
      deallocate(h_seed, v_seed, h_wann, v_wann, aa_wann, spread_est, tail_est)
      deallocate(sw_legacy, uw_legacy, lambda_legacy, wcoef_legacy, keep_index_legacy)
      deallocate(r_wann_legacy, tmp_legacy, wcenter_legacy)
      deallocate(phi_box, phi_tmp, psi_w, rho_w_sum, cos_sum, sin_sum)
      deallocate(n_basis_file)
      if(allocated(bond_center_bohr)) deallocate(bond_center_bohr)
    end subroutine write_local_wannier_seed

    subroutine apply_local_bond_center_symmetry_gauge(nbasis, nproj_seed, lambda_w, eigvec_w, &
        bproj_seed, wcoef, keep_index, nkeep)
      use salmon_global, only: lambda_cut
      implicit none
      integer, intent(in) :: nbasis, nproj_seed
      integer, intent(inout) :: keep_index(:), nkeep
      real(8), intent(in) :: lambda_w(:), eigvec_w(:,:), bproj_seed(:,:)
      real(8), intent(inout) :: wcoef(:,:)
      integer :: ip, jp, kp, io
      real(8), allocatable :: sinv_half(:,:), work(:,:)
      real(8) :: min_lambda_keep

      if(nbasis <= 0 .or. nproj_seed <= 0) return
      if(nkeep /= nproj_seed) then
        write(*,'(1x,a,i0,a,i0)') &
          "[DC-LCFO-LOCAL-WANNIER-SYM-GAUGE] skipped: keep=", nkeep, " nproj_seed=", nproj_seed
        return
      end if
      min_lambda_keep = minval(lambda_w(1:nproj_seed))
      if(min_lambda_keep <= lambda_cut) then
        write(*,'(1x,a,2(a,1pe12.4))') &
          "[DC-LCFO-LOCAL-WANNIER-SYM-GAUGE] skipped:", &
          " min_lambda=", min_lambda_keep, " lambda_cut=", lambda_cut
        return
      end if

      allocate(sinv_half(nproj_seed,nproj_seed), work(nbasis,nproj_seed))
      sinv_half = 0.0d0
      do ip=1,nproj_seed
        do jp=1,nproj_seed
          do kp=1,nproj_seed
            sinv_half(ip,jp) = sinv_half(ip,jp) + eigvec_w(ip,kp) * &
              (1.0d0 / sqrt(max(lambda_w(kp), 1.0d-300))) * eigvec_w(jp,kp)
          end do
        end do
      end do

      work(1:nbasis,1:nproj_seed) = 0.0d0
      do io=1,nbasis
        do jp=1,nproj_seed
          do ip=1,nproj_seed
            work(io,jp) = work(io,jp) + bproj_seed(io,ip) * sinv_half(ip,jp)
          end do
        end do
      end do
      wcoef(1:nbasis,1:nproj_seed) = work(1:nbasis,1:nproj_seed)
      do ip=1,nproj_seed
        keep_index(ip) = ip
      end do
      write(*,'(1x,a,i0,a,1pe12.4)') &
        "[DC-LCFO-LOCAL-WANNIER-SYM-GAUGE] applied lowdin_bond_center keep=", &
        nkeep, " min_lambda=", min_lambda_keep
      deallocate(sinv_half, work)
    end subroutine apply_local_bond_center_symmetry_gauge

    subroutine project_wannier90_aa_to_bpw(nbasis, nkeep, wcoef, aa_wann, aa_wann_source)
      use inputoutput, only: au_length_aa
      implicit none
      integer, intent(in) :: nbasis, nkeep
      real(8), intent(in) :: wcoef(nbasis,nkeep)
      real(8), intent(inout) :: aa_wann(3,nkeep,nkeep)
      integer, intent(out) :: aa_wann_source
      integer :: num_bands_file, num_wann_file
      integer :: axis, n, m
      logical :: ok_transform, ok_position
      complex(8), allocatable :: v_matrix(:,:), aa_global(:,:,:)
      complex(8), allocatable :: coef_local(:,:), wcoef_local(:,:), psi_wann(:,:)
      complex(8), allocatable :: tmat(:,:), tmp_c(:,:), aa_proj(:,:)

      aa_wann_source = 0
      call read_wannier90_global_transform_file(num_bands_file, num_wann_file, v_matrix, ok_transform)
      if(.not. ok_transform) return
      call read_wannier90_global_rmn_gamma_block(num_wann_file, aa_global, ok_position)
      if(.not. ok_position) then
        deallocate(v_matrix)
        return
      end if

      if(num_bands_file /= dc%nstate_tot) then
        write(*,'(1x,a,i0,a,i0,a)') "[DC-LCFO-BUFFER-WANNIER] global Wannier AA_R band space is not full BPW: bands=", &
          num_bands_file, " dc_nstate_tot=", dc%nstate_tot, "; BPW AA_R block is left zero."
        write(*,'(1x,a)') &
          "[DC-LCFO-BUFFER-WANNIER] Regenerate Wannier90 with num_bands=num_wann=the full BPW target subspace."
        deallocate(v_matrix, aa_global)
        return
      end if

      if(num_wann_file /= dc%nstate_tot) then
        write(*,'(1x,a,i0,a,i0,a)') "[DC-LCFO-BUFFER-WANNIER] global Wannier AA_R dimension is not full BPW: wann=", &
          num_wann_file, " dc_nstate_tot=", dc%nstate_tot, "; BPW AA_R block is left zero."
        write(*,'(1x,a)') &
          "[DC-LCFO-BUFFER-WANNIER] Regenerate Wannier90 with num_bands=num_wann=the full BPW target subspace."
        deallocate(v_matrix, aa_global)
        return
      end if

      allocate(coef_local(nbasis,num_bands_file), wcoef_local(nbasis,nkeep))
      allocate(psi_wann(nbasis,num_wann_file), tmat(nkeep,num_wann_file))
      allocate(tmp_c(nkeep,num_wann_file), aa_proj(nkeep,nkeep))
      coef_local(1:nbasis,1:num_bands_file) = &
        cmplx(coef_wf(1:nbasis,1:num_bands_file,1), 0d0, kind=8)
      wcoef_local(1:nbasis,1:nkeep) = cmplx(wcoef(1:nbasis,1:nkeep), 0d0, kind=8)

      psi_wann(1:nbasis,1:num_wann_file) = &
        matmul(coef_local(1:nbasis,1:num_bands_file), &
               v_matrix(1:num_bands_file,1:num_wann_file))
      tmat(1:nkeep,1:num_wann_file) = &
        matmul(transpose(conjg(wcoef_local(1:nbasis,1:nkeep))), &
               psi_wann(1:nbasis,1:num_wann_file))

      aa_wann(1:3,1:nkeep,1:nkeep) = 0d0
      do axis=1,3
        tmp_c(1:nkeep,1:num_wann_file) = &
          matmul(tmat(1:nkeep,1:num_wann_file), &
                 aa_global(axis,1:num_wann_file,1:num_wann_file))
        aa_proj(1:nkeep,1:nkeep) = &
          matmul(tmp_c(1:nkeep,1:num_wann_file), &
                 transpose(conjg(tmat(1:nkeep,1:num_wann_file))))
        do n=1,nkeep
          aa_wann(axis,n,n) = real(aa_proj(n,n), kind=8)
          do m=n+1,nkeep
            aa_wann(axis,n,m) = 0.5d0 * real(aa_proj(n,m) + conjg(aa_proj(m,n)), kind=8)
            aa_wann(axis,m,n) = aa_wann(axis,n,m)
          end do
        end do
      end do
      aa_wann_source = 1
      write(*,'(1x,a,i0,a,i0,a,i0,a,1pe12.4)') &
        "[DC-LCFO-BUFFER-WANNIER] projected Wannier90 AA_R to BPW: fragment=", dc%i_frag, &
        " keep=", nkeep, " global_wann=", num_wann_file, " max_abs=", &
        maxval(abs(aa_wann(1:3,1:nkeep,1:nkeep)))

      deallocate(v_matrix, aa_global, coef_local, wcoef_local, psi_wann, tmat, tmp_c, aa_proj)
    end subroutine project_wannier90_aa_to_bpw

    subroutine read_wannier90_global_transform_file(num_bands_file, num_wann_file, v_matrix, ok)
      use filesystem, only: get_filehandle
      implicit none
      integer, intent(out) :: num_bands_file, num_wann_file
      complex(8), allocatable, intent(out) :: v_matrix(:,:)
      logical, intent(out) :: ok
      integer :: iunit, io, magic_file, version_file, nfrag_file
      integer, allocatable :: owner_frag(:)
      real(8), allocatable :: center_bohr(:,:), spread_aa2(:)
      character(256) :: filename
      logical :: exists

      ok = .false.
      num_bands_file = 0
      num_wann_file = 0
      filename = trim(dc%base_directory)//binfile_w90g
      inquire(file=filename, exist=exists)
      if(.not. exists) then
        write(*,'(1x,3a)') "[DC-LCFO-BUFFER-WANNIER] global Wannier basis not found: ", &
          trim(filename), "; BPW AA_R block is left zero."
        return
      end if

      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='old',iostat=io)
      if(io /= 0) then
        write(*,'(1x,2a)') "[DC-LCFO-BUFFER-WANNIER] failed to open global Wannier basis: ", trim(filename)
        return
      end if
      read(iunit,iostat=io) magic_file, version_file
      if(io /= 0 .or. magic_file /= wannier90_global_magic .or. &
          version_file < 1 .or. version_file > wannier90_global_version) then
        close(iunit)
        write(*,'(1x,2a)') "[DC-LCFO-BUFFER-WANNIER] invalid global Wannier basis header: ", trim(filename)
        return
      end if
      read(iunit,iostat=io) num_bands_file, num_wann_file, nfrag_file
      if(io /= 0 .or. num_bands_file <= 0 .or. num_wann_file <= 0) then
        close(iunit)
        write(*,'(1x,2a)') "[DC-LCFO-BUFFER-WANNIER] invalid global Wannier basis dimensions: ", trim(filename)
        return
      end if
      allocate(owner_frag(num_wann_file), center_bohr(3,num_wann_file), spread_aa2(num_wann_file))
      allocate(v_matrix(num_bands_file,num_wann_file))
      read(iunit,iostat=io) owner_frag(1:num_wann_file)
      if(io == 0) read(iunit,iostat=io) center_bohr(1:3,1:num_wann_file)
      if(io == 0) read(iunit,iostat=io) spread_aa2(1:num_wann_file)
      if(io == 0) read(iunit,iostat=io) v_matrix(1:num_bands_file,1:num_wann_file)
      close(iunit)
      deallocate(owner_frag, center_bohr, spread_aa2)
      if(io /= 0) then
        write(*,'(1x,2a)') "[DC-LCFO-BUFFER-WANNIER] failed to read global Wannier transform: ", trim(filename)
        if(allocated(v_matrix)) deallocate(v_matrix)
        num_bands_file = 0
        num_wann_file = 0
        return
      end if
      ok = .true.
    end subroutine read_wannier90_global_transform_file

    subroutine read_wannier90_global_rmn_gamma_block(num_wann_expected, aa_global, ok)
      use filesystem, only: get_filehandle
      use inputoutput, only: au_length_aa
      use salmon_global, only: sysname
      implicit none
      integer, intent(in) :: num_wann_expected
      complex(8), allocatable, intent(out) :: aa_global(:,:,:)
      logical, intent(out) :: ok
      integer :: iunit, io, num_wann_file, nrpts_file, ir
      integer :: rvec(3), n, m
      real(8) :: rx_re, rx_im, ry_re, ry_im, rz_re, rz_im
      character(256) :: filename
      character(256) :: header
      logical :: exists

      ok = .false.
      filename = trim(dc%base_directory)//trim(sysname)//"_r.dat"
      inquire(file=filename, exist=exists)
      if(.not. exists) then
        write(*,'(1x,3a)') "[DC-LCFO-BUFFER-WANNIER] Wannier90 r-matrix not found: ", &
          trim(filename), "; BPW AA_R block is left zero."
        return
      end if

      iunit = get_filehandle()
      open(iunit,file=filename,status='old',action='read')
      read(iunit,'(a)',iostat=io) header
      if(io /= 0) then
        close(iunit)
        write(*,'(1x,2a)') "[DC-LCFO-BUFFER-WANNIER] failed to read Wannier90 r header: ", trim(filename)
        return
      end if
      read(iunit,*,iostat=io) num_wann_file
      if(io /= 0) then
        close(iunit)
        write(*,'(1x,2a)') "[DC-LCFO-BUFFER-WANNIER] failed to read Wannier90 r num_wann: ", trim(filename)
        return
      end if
      read(iunit,*,iostat=io) nrpts_file
      if(io /= 0) then
        close(iunit)
        write(*,'(1x,2a)') "[DC-LCFO-BUFFER-WANNIER] failed to read Wannier90 r nrpts: ", trim(filename)
        return
      end if

      if(num_wann_file /= num_wann_expected) then
        close(iunit)
        write(*,'(1x,a,i0,a,i0,a)') "[DC-LCFO-BUFFER-WANNIER] Wannier90 r-matrix size mismatch: file=", &
          num_wann_file, " expected=", num_wann_expected, "; BPW AA_R block is left zero."
        return
      end if

      allocate(aa_global(3,num_wann_file,num_wann_file))
      aa_global = (0d0,0d0)
      do ir=1,nrpts_file*num_wann_file*num_wann_file
        read(iunit,*,iostat=io) rvec(1:3), n, m, rx_re, rx_im, ry_re, ry_im, rz_re, rz_im
        if(io /= 0) exit
        if(any(rvec(1:3) /= 0)) cycle
        if(n < 1 .or. n > num_wann_file .or. m < 1 .or. m > num_wann_file) cycle
        aa_global(1,n,m) = cmplx(rx_re, rx_im, kind=8) / au_length_aa
        aa_global(2,n,m) = cmplx(ry_re, ry_im, kind=8) / au_length_aa
        aa_global(3,n,m) = cmplx(rz_re, rz_im, kind=8) / au_length_aa
      end do
      close(iunit)
      call hermitize_wannier_position_matrix(num_wann_file, aa_global)
      ok = .true.
    end subroutine read_wannier90_global_rmn_gamma_block

    integer function count_local_c_sp3_projections(nxyz_domain) result(nproj)
      use salmon_global, only: izatom, wannier_projection
      implicit none
      integer, intent(in) :: nxyz_domain(3)
      integer :: ia, ix_local(3), target_iz

      nproj = 0
      target_iz = sp3_projection_iz(trim(wannier_projection))
      do ia=1,dc%system_tot%nion
        if(izatom(dc%system_tot%kion(ia)) /= target_iz) cycle
        call find_atom_core_grid(nxyz_domain, ia, ix_local)
        if(all(ix_local(1:3) > 0)) nproj = nproj + 4
      end do
    end function count_local_c_sp3_projections

    subroutine build_local_c_sp3_projection_map(nxyz_domain, proj_atom, proj_hybrid)
      use salmon_global, only: izatom, wannier_projection
      implicit none
      integer, intent(in) :: nxyz_domain(3)
      integer, intent(out) :: proj_atom(:), proj_hybrid(:)
      integer :: ia, ih, ip, ix_local(3), target_iz

      ip = 0
      target_iz = sp3_projection_iz(trim(wannier_projection))
      do ia=1,dc%system_tot%nion
        if(izatom(dc%system_tot%kion(ia)) /= target_iz) cycle
        call find_atom_core_grid(nxyz_domain, ia, ix_local)
        if(any(ix_local(1:3) <= 0)) cycle
        do ih=1,4
          ip = ip + 1
          proj_atom(ip) = ia
          proj_hybrid(ip) = ih
        end do
      end do
    end subroutine build_local_c_sp3_projection_map

    integer function count_local_pseudo_channel_ao_candidates(nxyz_domain) result(nproj)
      use salmon_global, only: izatom
      implicit none
      integer, intent(in) :: nxyz_domain(3)
      integer :: ia, iz, ix_local(3), lmax_ao

      nproj = 0
      do ia=1,dc%system_tot%nion
        call find_atom_core_grid(nxyz_domain, ia, ix_local)
        if(any(ix_local(1:3) <= 0)) cycle
        iz = izatom(dc%system_tot%kion(ia))
        lmax_ao = pseudo_channel_ao_lmax_for_species(iz)
        nproj = nproj + count_real_ao_for_lmax(lmax_ao)
      end do
    end function count_local_pseudo_channel_ao_candidates

    subroutine build_local_pseudo_channel_ao_candidate_map(nxyz_domain, proj_atom, proj_code)
      use salmon_global, only: izatom
      implicit none
      integer, intent(in) :: nxyz_domain(3)
      integer, intent(out) :: proj_atom(:), proj_code(:)
      integer :: ia, iz, lmax_ao, ip, im, ix_local(3)

      ip = 0
      do ia=1,dc%system_tot%nion
        call find_atom_core_grid(nxyz_domain, ia, ix_local)
        if(any(ix_local(1:3) <= 0)) cycle
        iz = izatom(dc%system_tot%kion(ia))
        lmax_ao = pseudo_channel_ao_lmax_for_species(iz)
        if(lmax_ao >= 0) then
          ip = ip + 1
          proj_atom(ip) = ia
          proj_code(ip) = 1
        end if
        if(lmax_ao >= 1) then
          do im=1,3
            ip = ip + 1
            proj_atom(ip) = ia
            proj_code(ip) = 10 + im
          end do
        end if
        if(lmax_ao >= 2) then
          do im=1,5
            ip = ip + 1
            proj_atom(ip) = ia
            proj_code(ip) = 20 + im
          end do
        end if
      end do
    end subroutine build_local_pseudo_channel_ao_candidate_map

    subroutine build_local_bond_center_projection_map(nproj, center_bohr)
      implicit none
      integer, intent(out) :: nproj
      real(8), allocatable, intent(out) :: center_bohr(:,:)
      real(8), allocatable :: all_center_bohr(:,:)
      integer :: iw, ip, nall

      nall = max(1, count_bond_center_candidates())
      call build_bond_center_projection_map(nall, all_center_bohr)
      nproj = 0
      do iw=1,nall
        if(find_owner_fragment_from_center(all_center_bohr(1:3,iw)) == dc%i_frag) nproj = nproj + 1
      end do
      if(nproj <= 0) then
        allocate(center_bohr(3,1))
        center_bohr = 0.0d0
      else
        allocate(center_bohr(3,nproj))
        ip = 0
        do iw=1,nall
          if(find_owner_fragment_from_center(all_center_bohr(1:3,iw)) /= dc%i_frag) cycle
          ip = ip + 1
          center_bohr(1:3,ip) = all_center_bohr(1:3,iw)
        end do
      end if
      if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0)') &
        "[DC-LCFO-LOCAL-WANNIER] fragment=", dc%i_frag, " local bond-center seeds=", nproj
      deallocate(all_center_bohr)
    end subroutine build_local_bond_center_projection_map

    subroutine find_atom_core_grid(nxyz_domain, ia, ix_local)
      implicit none
      integer, intent(in) :: nxyz_domain(3), ia
      integer, intent(out) :: ix_local(3)
      integer :: axis, i, ig, ibest
      real(8) :: r_atom, r_grid, dist, best_dist, length_axis, spacing_axis

      ix_local = 0
      do axis=1,3
        r_atom = dc%system_tot%rion(axis,ia)
        length_axis = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
          + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
        spacing_axis = length_axis / dble(dc%lg_tot%num(axis))
        best_dist = huge(1d0)
        ibest = 0
        do i=1,nxyz_domain(axis)
          ig = dc%jxyz_tot(i,axis)
          r_grid = dc%lg_tot%coordinate(ig,axis)
          dist = abs(periodic_delta(r_grid - r_atom, length_axis))
          if(dist < best_dist) then
            best_dist = dist
            ibest = i
          end if
        end do
        if(best_dist <= 0.75d0 * spacing_axis) ix_local(axis) = ibest
      end do
    end subroutine find_atom_core_grid

    subroutine prepare_sawf_closed_seed_eigensystem(nband_wann)
      use communication, only: comm_get_max,comm_bcast,comm_summation
      use salmon_global, only: izatom,wannier_site_symmetry,wannier_symmetry_file, &
        wannier_symmetry_tolerance
      use salmon_math, only: matrix_inverse
      implicit none
      integer, intent(in) :: nband_wann
      integer :: mesh(3),ifrag,isym,ia,failure,local_failure,allocation_status
      integer :: max_targets,nbasis_closed,npoint_core,npoint_buffer
      integer :: itmp(dc%n_frag,nspin)
      integer, allocatable :: species(:),fragment_origin(:,:),fragment_shape(:,:)
      integer, allocatable :: fragment_map(:),fragment_maps(:,:)
      real(8) :: a1(3),a2(3),a3(3),lattice(3,3),lattice_inverse(3,3)
      real(8) :: max_grid_residual,center_grid(3)
      real(8), allocatable :: fractional_positions(:,:)
      type(t_sawf_symop), allocatable :: operations(:)
      type(t_sawf_fragment_state_cache) :: cache
      type(t_sawf_closed_basis) :: closed
      logical :: local_ok,grid_ok,fragment_ok,center_available,split_fragment_global_mode
      character(512) :: message
      character(256) :: symmetry_filename
      integer :: nxyz_box(3),ibasis

      if(dc%n_frag==1) then
        if(dc%id_tot==0) write(*,'(1x,a)') &
          '[DC-LCFO-SAWF-GLOBAL-SINGLE] single fragment already spans the global LCFO eigensystem'
        return
      end if

      mesh=dc%lg_tot%num
      allocate(species(dc%system_tot%nion),fractional_positions(3,dc%system_tot%nion), &
        fragment_origin(3,dc%n_frag),fragment_shape(3,dc%n_frag),stat=allocation_status)
      failure=merge(0,1,allocation_status==0)
      call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) call lcfo_sawf_fatal('SAWF closed-seed context allocation failed')
      fragment_origin=dc%ixyz_frag
      split_fragment_global_mode=(dc%n_frag==1)
      do ifrag=1,dc%n_frag
        call get_fragment_domain(dc,ifrag,fragment_shape(:,ifrag))
      end do
      call get_lattice_vectors(a1,a2,a3)
      lattice(:,1)=a1; lattice(:,2)=a2; lattice(:,3)=a3
      lattice_inverse=lattice; call matrix_inverse(lattice_inverse)
      do ia=1,dc%system_tot%nion
        species(ia)=izatom(dc%system_tot%kion(ia))
        fractional_positions(:,ia)=modulo(matmul(lattice_inverse,dc%system_tot%rion(:,ia)),1.0d0)
      end do
      if(trim(wannier_site_symmetry)=='auto') then
        call load_sawf_symmetry_auto(lattice,fractional_positions,species, &
          wannier_symmetry_tolerance,operations,local_ok,message)
      else
        if(wannier_symmetry_file(1:1)=='/') then
          symmetry_filename=trim(wannier_symmetry_file)
        else
          symmetry_filename=trim(import_run_root_dir())//trim(wannier_symmetry_file)
        end if
        call load_sawf_symmetry_file(symmetry_filename,lattice,fractional_positions,species, &
          wannier_symmetry_tolerance,operations,local_ok,message)
      end if
      failure=merge(0,1,local_ok); call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) call lcfo_sawf_fatal('SAWF closed-seed symmetry loading failed')
      call put_sawf_identity_first(operations,wannier_symmetry_tolerance,local_ok,message)
      failure=merge(0,1,local_ok); call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) call lcfo_sawf_fatal('SAWF closed-seed operation normalization failed')
      if(size(operations)>=2) call diagnose_sawf_vlocal_symmetry(operations(2),mesh, &
        fragment_origin,fragment_shape,wannier_symmetry_tolerance)
      allocate(fragment_maps(dc%n_frag,size(operations)),stat=allocation_status)
      failure=merge(0,1,allocation_status==0); call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) call lcfo_sawf_fatal('SAWF closed-seed fragment-map allocation failed')
      do isym=1,size(operations)
        call validate_sawf_fragment_symmetry_map(operations(isym),mesh,fragment_origin, &
          fragment_shape,dc%nxyz_buffer,wannier_symmetry_tolerance,grid_ok,fragment_ok, &
          max_targets,fragment_map,max_grid_residual,center_available,center_grid,message)
        local_failure=merge(0,1,grid_ok)
        call reduce_sawf_fragment_alignment_failure(local_failure,dc%icomm_tot,dc%id_tot,isym, &
          grid_ok,fragment_ok,max_targets,max_grid_residual,message,failure)
        if(failure/=0) call lcfo_sawf_fatal('SAWF closed-seed grid alignment failed')
        local_failure=merge(0,1,fragment_ok)
        call comm_get_max(local_failure,dc%icomm_tot)
        if(local_failure/=0) then
          split_fragment_global_mode=.true.
          if(allocated(fragment_map)) deallocate(fragment_map)
          exit
        end if
        fragment_maps(:,isym)=fragment_map
        deallocate(fragment_map)
      end do
      if(split_fragment_global_mode) then
        if(dc%id_tot==0) write(*,'(1x,a)') &
          '[DC-LCFO-SAWF-GLOBAL-SPLIT] symmetry splits fragment cores; use global LCFO eigensystem'
        call clear_sawf_fragment_state_cache(cache)
        deallocate(species,fractional_positions,fragment_origin,fragment_shape,fragment_maps,operations)
        return
      end if
      call prepare_sawf_fragment_state_cache(nband_wann,fragment_shape,cache,local_ok,message)
      failure=merge(0,1,local_ok); call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) call lcfo_sawf_fatal('SAWF closed-seed source cache failed')
      call build_sawf_closed_fragment_seed_basis(nband_wann,mesh,fragment_origin,fragment_shape, &
        operations,fragment_maps,cache,closed,local_ok,message)
      failure=merge(0,1,local_ok); call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) call lcfo_sawf_fatal('SAWF closed-seed basis construction failed')

      nbasis_closed=closed%nbasis; npoint_core=closed%npoint_core
      npoint_buffer=closed%npoint_buffer
      call comm_bcast(nbasis_closed,dc%icomm_frag,0)
      call comm_bcast(npoint_core,dc%icomm_frag,0)
      call comm_bcast(npoint_buffer,dc%icomm_frag,0)
      if(nbasis_closed>dc%nstate_frag) call lcfo_sawf_fatal( &
        'SAWF closed basis exceeds nstate_frag; dynamic batching is required')
      if(dc%id_frag/=0) then
        allocate(closed%core(npoint_core,nbasis_closed),closed%buffer(npoint_buffer,nbasis_closed))
        closed%nbasis=nbasis_closed; closed%npoint_core=npoint_core; closed%npoint_buffer=npoint_buffer
      end if
      call comm_bcast(closed%core,dc%icomm_frag,0)
      call comm_bcast(closed%buffer,dc%icomm_frag,0)

      f_basis=0.0d0
      do ibasis=1,nbasis_closed
        f_basis(:,:,:,1,ibasis)=reshape(closed%core(:,ibasis),dc%nxyz_domain)
      end do
      nxyz_box=dc%nxyz_domain+2*dc%nxyz_buffer
      if(allocated(sawf_explicit_buffer)) deallocate(sawf_explicit_buffer)
      allocate(sawf_explicit_buffer(nxyz_box(1),nxyz_box(2),nxyz_box(3),nspin,dc%nstate_frag))
      sawf_explicit_buffer=0.0d0
      do ibasis=1,nbasis_closed
        sawf_explicit_buffer(:,:,:,1,ibasis)=reshape(closed%buffer(:,ibasis),nxyz_box)
      end do
      sawf_explicit_basis_active=.true.
      itmp=0
      if(dc%id_frag==0) itmp(dc%i_frag,1)=nbasis_closed
      call comm_summation(itmp,n_basis,dc%n_frag*nspin,dc%icomm_tot)
      index_basis=0
      do isym=1,nspin
        ia=0
        do ifrag=1,dc%n_frag
          do ibasis=1,n_basis(ifrag,isym)
            ia=ia+1; index_basis(ibasis,ifrag,isym)=ia
          end do
        end do
        n_mat(isym)=ia
      end do

      if(allocated(hf)) deallocate(hf)
      if(allocated(mat_H_local)) deallocate(mat_H_local)
      if(allocated(mat_H_volume_local)) deallocate(mat_H_volume_local)
      if(allocated(mat_H_volume_weak_local)) deallocate(mat_H_volume_weak_local)
      if(allocated(mat_H_weak_kinetic)) deallocate(mat_H_weak_kinetic)
      if(allocated(mat_H_weak_potential)) deallocate(mat_H_weak_potential)
      if(allocated(mat_H_weak_nonlocal)) deallocate(mat_H_weak_nonlocal)
      if(allocated(mat_H_surface_self)) deallocate(mat_H_surface_self)
      if(allocated(mat_V_local)) deallocate(mat_V_local)
      call release_surface_trace_halo()
      do isym=1,n_halo
        if(allocated(halo(isym)%mat_H_local)) deallocate(halo(isym)%mat_H_local)
        if(allocated(halo(isym)%mat_H_surface_cross)) deallocate(halo(isym)%mat_H_surface_cross)
        if(allocated(halo(isym)%mat_V_local)) deallocate(halo(isym)%mat_V_local)
        if(allocated(halo(isym)%mat_Xi_flux_local)) deallocate(halo(isym)%mat_Xi_flux_local)
      end do
      call hpsi_basis
      call calc_hamiltonian_matrix
      if(allocated(esp_tot)) deallocate(esp_tot)
      if(allocated(coef_wf)) deallocate(coef_wf)
      allocate(esp_tot(maxval(n_mat),nspin))
      if(dc%id_frag==0) allocate(coef_wf(dc%nstate_frag,dc%nstate_tot,nspin))
      if(allocated(coef_wf)) coef_wf=0.0d0
#ifdef USE_EIGENEXA
      call diag_eigenexa
#else
      call lcfo_sawf_fatal('SAWF closed-seed diagonalization requires EigenExa')
#endif
      call write_sawf_closed_seed_file()
      if(dc%id_tot==0) write(*,'(1x,a,i0,a,i0)') &
        '[DC-LCFO-SAWF-CLOSED] physical Flux eigensystem basis=',maxval(n_basis),' bands=',nband_wann
      call clear_sawf_closed_basis(closed)
      call clear_sawf_fragment_state_cache(cache)
      deallocate(species,fractional_positions,fragment_origin,fragment_shape,fragment_maps,operations)
    end subroutine prepare_sawf_closed_seed_eigensystem

    subroutine diagnose_sawf_vlocal_symmetry(operation,mesh,fragment_origin,fragment_shape,tolerance)
      use communication, only: comm_summation
      implicit none
      type(t_sawf_symop), intent(in) :: operation
      integer, intent(in) :: mesh(3),fragment_origin(:,:),fragment_shape(:,:)
      real(8), intent(in) :: tolerance
      real(8), allocatable :: value_local(:),value_global(:),difference(:)
      integer, allocatable :: count_local(:),count_global(:)
      integer :: ix,iy,iz,p,q,source(3),target(3),location(1),npoint
      integer :: source_location(3),target_location(3),failure
      real(8) :: scale,relative_max,relative_rms
      logical :: map_ok
      character(512) :: message

      npoint=product(mesh)
      allocate(value_local(npoint),value_global(npoint),difference(npoint), &
        count_local(npoint),count_global(npoint))
      value_local=0.0d0; count_local=0
      if(dc%id_frag==0) then
        do iz=1,fragment_shape(3,dc%i_frag)
          do iy=1,fragment_shape(2,dc%i_frag)
            do ix=1,fragment_shape(1,dc%i_frag)
              source=1+modulo(fragment_origin(:,dc%i_frag)+[ix-1,iy-1,iz-1],mesh)
              p=source(1)+(source(2)-1)*mesh(1)+(source(3)-1)*mesh(1)*mesh(2)
              value_local(p)=V_local(1)%f(ix,iy,iz)
              count_local(p)=1
            end do
          end do
        end do
      end if
      call comm_summation(value_local,value_global,npoint,dc%icomm_tot)
      call comm_summation(count_local,count_global,npoint,dc%icomm_tot)
      failure=merge(0,1,all(count_global==1))
      if(dc%id_tot==0 .and. failure/=0) write(*,'(1x,a,i0,a,i0)') &
        '[DC-LCFO-SAWF-VLOCAL-SYMMETRY] ownership_min=',minval(count_global), &
        ' ownership_max=',maxval(count_global)
      if(failure/=0) then
        deallocate(value_local,value_global,difference,count_local,count_global)
        return
      end if
      if(dc%id_tot==0) then
        difference=0.0d0
        do iz=1,mesh(3)
          do iy=1,mesh(2)
            do ix=1,mesh(1)
              source=[ix,iy,iz]
              call map_sawf_periodic_grid_point(operation,mesh,tolerance,source,target,map_ok,message)
              if(.not.map_ok) call lcfo_sawf_fatal('SAWF V_local symmetry grid map failed: '//trim(message))
              p=ix+(iy-1)*mesh(1)+(iz-1)*mesh(1)*mesh(2)
              q=target(1)+(target(2)-1)*mesh(1)+(target(3)-1)*mesh(1)*mesh(2)
              difference(p)=value_global(q)-value_global(p)
            end do
          end do
        end do
        scale=max(1.0d-300,maxval(abs(value_global)))
        relative_max=maxval(abs(difference))/scale
        relative_rms=sqrt(sum(difference*difference)/sum(value_global*value_global))
        location=maxloc(abs(difference)); p=location(1)
        source_location(1)=1+modulo(p-1,mesh(1))
        source_location(2)=1+modulo((p-1)/mesh(1),mesh(2))
        source_location(3)=1+(p-1)/(mesh(1)*mesh(2))
        call map_sawf_periodic_grid_point(operation,mesh,tolerance,source_location, &
          target_location,map_ok,message)
        write(*,'(1x,a,2(a,es13.5),2(a,3(i0,1x)))') &
          '[DC-LCFO-SAWF-VLOCAL-SYMMETRY]',' relative_max=',relative_max, &
          ' relative_rms=',relative_rms,' source=',source_location,' target=',target_location
      end if
      deallocate(value_local,value_global,difference,count_local,count_global)
    end subroutine diagnose_sawf_vlocal_symmetry

    subroutine write_sawf_closed_seed_file()
      use filesystem, only: get_filehandle
      use salmon_global, only: base_directory
      implicit none
      integer :: iunit
      character(256) :: filename
      if(dc%id_frag/=0) return
      filename=trim(base_directory)//'sawf_closed_seed.bin'
      iunit=get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='replace')
      write(iunit) sawf_closed_seed_magic,sawf_closed_seed_version
      write(iunit) current_sawf_seed_token
      write(iunit) dc%nxyz_domain,dc%nxyz_buffer,nspin,dc%nstate_frag,dc%nstate_tot
      write(iunit) n_basis(dc%i_frag,1)
      write(iunit) f_basis(:,:,:,1,1:n_basis(dc%i_frag,1))
      write(iunit) sawf_explicit_buffer(:,:,:,1,1:n_basis(dc%i_frag,1))
      write(iunit) coef_wf(1:n_basis(dc%i_frag,1),1:dc%nstate_tot,1)
      write(iunit) esp_tot(1:dc%nstate_tot,1)
      close(iunit)
    end subroutine write_sawf_closed_seed_file

    subroutine write_wannier_seed_files()
      use salmon_global, only: izatom, sysname, &
        wannier_projection, wannier_num_wann, wannier_num_iter, &
        wannier_projection_width, wannier_dis_froz_max, wannier_dis_win_max, &
        yn_dc_lcfo_wannier_cluster, yn_dc_lcfo_local_wannier, wannier_site_symmetry
      use filesystem, only: get_filehandle
      use communication, only: comm_bcast
      use communication, only: comm_sync_all
      implicit none
      real(8), parameter :: hartree_to_ev = 27.211386245988d0
      integer :: iunit, iband, ikpt, nband_wann, nproj_csp3
      integer :: projection_failure
      character(256) :: winfile, eigfile, projection_failure_message
      character(1024) :: resolved_wannier_command
      logical :: projection_ok
      real(8) :: a1(3), a2(3), a3(3)

      winfile = trim(dc%base_directory)//trim(sysname)//".win"
      if(trim(wannier_site_symmetry) /= 'off') call deactivate_sawf_win_collective(winfile)
      call get_cached_wannier90_command(resolved_wannier_command, projection_ok)
      if(.not. projection_ok) call lcfo_sawf_fatal( &
        'Wannier90 command was not resolved during initialization')
      if(trim(wannier_site_symmetry) /= 'off') &
        call reject_sawf_reuse_collective(resolved_wannier_command)
      if(nspin /= 1) stop "DC-LCFO Wannier export: spin-polarized Wannier seed files are not implemented."

      call validate_sawf_projection_configuration()
      projection_failure = 0
      projection_failure_message = ''

      nband_wann = determine_wannier_num_bands()
      nproj_csp3 = count_c_sp3_projections()
      if(trim(wannier_site_symmetry)/='off') call prepare_sawf_closed_seed_eigensystem(nband_wann)
      eigfile = trim(dc%base_directory)//trim(sysname)//".eig"
      call get_lattice_vectors(a1, a2, a3)
      if(yn_dc_lcfo_wannier_cluster == 'y') call log_wannier_cluster_partition()

      if(dc%id_tot == 0) then
        write(*,'(1x,a,i0,a,i0)') "[DC-LCFO-WANNIER] export bands=", nband_wann, &
          " wann=", wannier_num_wann
        call write_wannier_base_win_atomic(winfile, nband_wann, nproj_csp3, &
          a1, a2, a3, projection_ok, projection_failure_message)
        if(.not. projection_ok) projection_failure = 1

        if(projection_failure == 0) iunit = get_filehandle()
        if(projection_failure == 0) then
        open(iunit,file=eigfile,status='replace')
        do ikpt=1,1
          do iband=1,nband_wann
            write(iunit,'(2i8,1x,es23.15)') iband, ikpt, esp_tot(iband,1) * hartree_to_ev
          end do
        end do
        close(iunit)
        end if
      end if

      call comm_bcast(projection_failure, dc%icomm_tot, 0)
      call comm_bcast(projection_failure_message, dc%icomm_tot, 0)
      if(projection_failure /= 0) call lcfo_sawf_fatal(projection_failure_message)

      call diagnose_wannier_coef_rank(nband_wann)
      call write_wannier_amn_file(nband_wann)
      call diagnose_wannier_amn_conditioning(nband_wann, wannier_num_wann)
      call write_wannier_mmn_file(nband_wann)
      call comm_sync_all(dc%icomm_tot)
      if(trim(wannier_site_symmetry) /= 'off') then
        call generate_sawf_dmn(nband_wann,resolved_wannier_command)
        call activate_sawf_win_collective(winfile)
      end if

      if(dc%id_tot == 0) &
        write(*,'(1x,a,2a)') "[DC-LCFO-WANNIER] wrote seed files: ", trim(winfile), " and .eig/.amn/.mmn"
      if(is_wannier90_export_only_command(resolved_wannier_command)) then
        if(dc%id_tot == 0) then
          write(*,'(1x,a)') "[DC-LCFO-WANNIER] export-only mode: external Wannier90 and checkpoint import are skipped."
          if(trim(wannier_site_symmetry) /= 'off') then
            write(*,'(1x,a)') &
              "  Run Wannier90 separately in the seed directory, then rerun with wannier90_command='import_only'."
          else
            write(*,'(1x,a)') &
              "  Run Wannier90 separately in the seed directory, then rerun/import with SALMON_WANNIER90_COMMAND=reuse."
          end if
        end if
        return
      end if

      call run_wannier90_seed_files(resolved_wannier_command)
      call write_wannier90_global_basis_file(nband_wann)
      call write_wannier90_flux_eigen_seed_file(nband_wann)
    end subroutine write_wannier_seed_files

    subroutine write_wannier_base_win_atomic(winfile, nband_wann, nproj_csp3, &
        a1, a2, a3, ok, message)
      use salmon_global, only: izatom, wannier_projection, wannier_num_wann, wannier_num_iter, &
        wannier_dis_froz_max, wannier_dis_win_max
      implicit none
      character(*), intent(in) :: winfile
      integer, intent(in) :: nband_wann, nproj_csp3
      real(8), intent(in) :: a1(3), a2(3), a3(3)
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      real(8), parameter :: hartree_to_ev = 27.211386245988d0
      type(t_atomic_win_writer) :: writer
      real(8), allocatable :: bond_center_bohr(:,:)
      integer :: iunit, io_status, ia, ip
      character(512) :: io_message
      logical :: projection_ok

      call begin_atomic_win(writer, winfile, iunit, ok, message)
      if(.not. ok) return
      io_status = 0
      write(iunit,'(a)',iostat=io_status,iomsg=io_message) &
        'num_bands = '//trim(adjustl(int_to_string(nband_wann)))
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) &
        'num_wann = '//trim(adjustl(int_to_string(wannier_num_wann)))
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) &
        'num_iter = '//trim(adjustl(int_to_string(wannier_num_iter)))
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'mp_grid = 1 1 1'
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'gamma_only = true'
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'write_hr = true'
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'write_rmn = true'
      if(io_status == 0 .and. wannier_dis_froz_max > 0d0) write(iunit,'(a,1x,es23.15)', &
        iostat=io_status,iomsg=io_message) 'dis_froz_max =', wannier_dis_froz_max*hartree_to_ev
      if(io_status == 0 .and. wannier_dis_win_max > 0d0) write(iunit,'(a,1x,es23.15)', &
        iostat=io_status,iomsg=io_message) 'dis_win_max =', wannier_dis_win_max*hartree_to_ev
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) ''
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'begin unit_cell_cart'
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'bohr'
      if(io_status == 0) write(iunit,'(3es23.15)',iostat=io_status,iomsg=io_message) a1
      if(io_status == 0) write(iunit,'(3es23.15)',iostat=io_status,iomsg=io_message) a2
      if(io_status == 0) write(iunit,'(3es23.15)',iostat=io_status,iomsg=io_message) a3
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'end unit_cell_cart'
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) ''
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'begin atoms_cart'
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'bohr'
      do ia=1,dc%system_tot%nion
        if(io_status /= 0) exit
        write(iunit,'(a,3(1x,es23.15))',iostat=io_status,iomsg=io_message) &
          trim(element_symbol(izatom(dc%system_tot%kion(ia)))), dc%system_tot%rion(:,ia)
      end do
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'end atoms_cart'
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) ''
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'begin kpoints'
      if(io_status == 0) write(iunit,'(3f12.6)',iostat=io_status,iomsg=io_message) 0d0,0d0,0d0
      if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'end kpoints'
      if(io_status == 0 .and. len_trim(wannier_projection) > 0) then
        write(iunit,'(a)',iostat=io_status,iomsg=io_message) ''
        if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'begin projections'
        if(io_status == 0 .and. is_bond_center_projection(trim(wannier_projection))) then
          call build_bond_center_projection_map(wannier_num_wann, bond_center_bohr)
          do ip=1,wannier_num_wann
            write(iunit,'(a,f18.12,a,f18.12,a,f18.12,a)',iostat=io_status,iomsg=io_message) &
              'f=', bond_center_bohr(1,ip)/cell_length(a1), ',', &
              bond_center_bohr(2,ip)/cell_length(a2), ',', &
              bond_center_bohr(3,ip)/cell_length(a3), ':s'
            if(io_status /= 0) exit
          end do
          deallocate(bond_center_bohr)
        else if(io_status == 0 .and. is_pseudo_channel_projection(trim(wannier_projection))) then
          call write_pseudo_channel_projection_block(iunit, projection_ok, message)
          if(.not. projection_ok) io_status = 1
        else if(io_status == 0) then
          if(is_sp3_projection(trim(wannier_projection)) .and. nproj_csp3 < wannier_num_wann) &
            write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'random'
          if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) trim(wannier_projection)
        end if
        if(io_status == 0) write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'end projections'
      end if
      if(io_status /= 0) then
        if(len_trim(message) == 0) message = 'atomic base .win write failed: '//trim(io_message)
        call abort_atomic_win(writer); ok = .false.; return
      end if
      call finish_atomic_win(writer, ok, message)
    end subroutine write_wannier_base_win_atomic

    subroutine deactivate_sawf_win_collective(winfile)
      use communication, only: comm_bcast
      implicit none
      character(*), intent(in) :: winfile
      integer :: failure
      logical :: local_ok
      character(512) :: message
      failure = 0; message = ''
      if(dc%id_tot == 0) then
        call deactivate_sawf_win(winfile, local_ok, message)
        failure = merge(0,1,local_ok)
      end if
      call comm_bcast(failure,dc%icomm_tot,0)
      call comm_bcast(message,dc%icomm_tot,0)
      if(failure /= 0) call lcfo_sawf_fatal('SAWF stale .win deactivation failed: '//trim(message))
    end subroutine deactivate_sawf_win_collective

    subroutine reject_sawf_reuse_collective(resolved_command)
      use communication, only: comm_bcast
      implicit none
      character(*), intent(in) :: resolved_command
      integer :: failure
      failure = 0
      if(dc%id_tot == 0 .and. (is_wannier90_reuse_command(resolved_command).or. &
          is_wannier90_import_only_command(resolved_command))) failure = 1
      call comm_bcast(failure,dc%icomm_tot,0)
      if(failure /= 0) call lcfo_sawf_fatal( &
        'SAWF forbids reuse/import of an existing Wannier90 checkpoint without a seed/DMN fingerprint binding')
    end subroutine reject_sawf_reuse_collective

    subroutine activate_sawf_win_collective(winfile)
      use communication, only: comm_bcast
      use salmon_global, only: wannier_symmetry_tolerance
      implicit none
      character(*), intent(in) :: winfile
      integer :: activation_failure
      logical :: activation_ok
      character(512) :: activation_message

      activation_failure = 0
      activation_message = ''
      if(dc%id_tot == 0) then
        call activate_sawf_win(winfile, wannier_symmetry_tolerance, activation_ok, &
          activation_message)
        activation_failure = merge(0, 1, activation_ok)
      end if
      call comm_bcast(activation_failure, dc%icomm_tot, 0)
      call comm_bcast(activation_message, dc%icomm_tot, 0)
      if(activation_failure /= 0) call lcfo_sawf_fatal( &
        'SAWF .win activation failed after .dmn publication: '//trim(activation_message))
      if(dc%id_tot == 0) write(*,'(1x,a,1x,a)') &
        '[DC-LCFO-SAWF] activated site_symmetry in', trim(winfile)
    end subroutine activate_sawf_win_collective

    subroutine generate_sawf_dmn(nband_wann,resolved_wannier_command)
      use communication, only: comm_get_max,comm_bcast,comm_summation
      use filesystem, only: atomic_create_directory
      use salmon_global, only: izatom, sysname, wannier_num_wann, &
        wannier_site_symmetry, wannier_symmetry_file, wannier_symmetry_tolerance, &
        wannier_sawf_generation,wannier_sawf_structure_class
      use salmon_math, only: matrix_inverse
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      integer, intent(in) :: nband_wann
      character(*),intent(in)::resolved_wannier_command
      integer :: allocation_failure, allocation_status, failure, local_failure, ia, ifrag, iop, isym,ibundle
      integer :: projection_lmax
      integer :: mesh(3)
      integer, allocatable :: species(:), fragment_origin(:,:), fragment_shape(:,:)
      integer, allocatable :: symmetry_fragment_map(:),symmetry_fragment_maps(:,:)
      integer, allocatable :: product_left(:),product_right(:),product_result(:)
      integer, allocatable :: sawf_environment_orbit(:)
      integer, allocatable :: sawf_representative_fragment(:),sawf_materialize_operation(:)
      integer,allocatable :: sawf_expected_channels(:),sawf_selected_channels(:)
      integer,allocatable :: sawf_neighbor_gvec(:,:)
      integer,allocatable :: sawf_local_stabilizer(:)
      integer,allocatable :: sawf_local_point_map(:,:),sawf_point_map_column(:)
      integer,allocatable :: sawf_local_product_left(:),sawf_local_product_right(:), &
        sawf_local_product_result(:)
      integer,allocatable :: sawf_receipt_local(:),sawf_receipt_global(:), &
        sawf_bands_local(:),sawf_bands_global(:),sawf_wann_local(:),sawf_wann_global(:)
      character(256), allocatable :: sawf_environment_key(:)
      real(8) :: a1(3), a2(3), a3(3), lattice(3,3), lattice_inverse(3,3)
      real(8) :: singular_min, singular_max, closure_residual, closure_tolerance
      real(8) :: max_grid_residual,center_grid(3)
      real(8), allocatable :: fractional_positions(:,:),d_wann_real(:,:,:)
      real(8), allocatable :: sawf_vacuum_fraction(:)
      complex(8), allocatable :: d_band_local(:,:),d_band_sum(:,:),d_wann(:,:),amn(:,:)
      complex(8), allocatable :: d_band_set(:,:,:),d_wann_set(:,:,:)
      complex(8),allocatable :: sawf_local_d_wann(:,:,:)
      complex(8),allocatable :: sawf_local_d_band(:,:,:),sawf_local_states(:,:)
      complex(8),allocatable :: sawf_local_amn(:,:)
      complex(8),allocatable :: sawf_local_v_matrix(:,:),sawf_representative_orbitals(:,:)
      real(8),allocatable :: sawf_local_energy(:)
      real(8),allocatable :: sawf_local_centers(:,:),sawf_local_spreads(:)
      real(8) :: representation_residual
      type(t_sawf_projection_channel), allocatable :: channels(:)
      type(t_sawf_symop), allocatable :: symmetry_operations(:)
      type(t_sawf_dmn_writer) :: writer
      type(t_sawf_dmn_writer)::sawf_local_writer
      type(t_sawf_fragment_state_cache) :: state_cache
      type(t_sawf_closed_basis) :: closed_basis
      type(t_sawf_environment_receipt),allocatable :: sawf_environment_receipts(:)
      type(t_sawf_seed_bundle),allocatable :: sawf_seed_bundles(:)
      logical :: local_ok,grid_map_ok,fragment_map_ok,center_available,split_fragment_global_mode
      logical, allocatable :: sawf_environment_equivalent(:,:),sawf_defect_intersects(:), &
        sawf_regenerate_environment(:),sawf_generate_independently(:),sawf_inside_atom(:)
      integer :: max_targets_per_source,local_left,local_right,local_relation,global_relation, &
        num_bands_chk,num_wann_chk
      character(512) :: message
      character(256) :: symmetry_filename,allocation_message,dmn_filename,amn_filename, &
        sawf_supercell_fingerprint
      character(512)::local_chk_filename

      if(trim(wannier_site_symmetry) == 'off') return
      ! The scalable SAWF route is admitted by representation, provenance,
      ! gauge and operator gates; Wannier90 site_symmetry alone is insufficient.
      split_fragment_global_mode=(dc%n_frag==1)
      call sawf_projection_shell_lmax(dc%system_tot%nion, wannier_num_wann, projection_lmax, local_ok, message)
      if(.not. local_ok) call lcfo_sawf_fatal(message)
      mesh = dc%lg_tot%num(1:3)
      if(any(mesh <= 0) .or. nband_wann <= 0) &
        call lcfo_sawf_fatal('SAWF D_band requires positive global mesh and band dimensions')
      if(int(nband_wann,int64) > huge(0_int64)/int(nband_wann,int64)) &
        call lcfo_sawf_fatal('SAWF D_band matrix size overflows int64')

      allocation_failure = 0
      allocate(species(dc%system_tot%nion), stat=allocation_status, errmsg=allocation_message)
      if(allocation_status /= 0) allocation_failure = 1
      if(allocation_failure == 0) then
        allocate(fractional_positions(3,dc%system_tot%nion), stat=allocation_status, &
          errmsg=allocation_message)
        if(allocation_status /= 0) allocation_failure = 1
      end if
      if(allocation_failure == 0) then
        allocate(fragment_origin(3,dc%n_frag),fragment_shape(3,dc%n_frag), &
          d_band_local(nband_wann,nband_wann),d_band_sum(nband_wann,nband_wann), &
          stat=allocation_status,errmsg=allocation_message)
        if(allocation_status /= 0) allocation_failure = 1
      end if
      call comm_get_max(allocation_failure,dc%icomm_tot)
      if(allocation_failure /= 0) then
        if(allocated(species)) deallocate(species)
        if(allocated(fractional_positions)) deallocate(fractional_positions)
        if(allocated(fragment_origin)) deallocate(fragment_origin)
        if(allocated(fragment_shape)) deallocate(fragment_shape)
        if(allocated(d_band_local)) deallocate(d_band_local)
        if(allocated(d_band_sum)) deallocate(d_band_sum)
        call lcfo_sawf_fatal('SAWF D_band allocation failed on one or more ranks')
      end if

      fragment_origin=dc%ixyz_frag
      do ifrag=1,dc%n_frag
        call get_fragment_domain(dc,ifrag,fragment_shape(:,ifrag))
      end do
      call validate_sawf_fragment_tiling(mesh,fragment_origin,fragment_shape,local_ok,message)
      failure=merge(0,1,local_ok)
      if(failure /= 0) write(*,'(1x,a,i0,2a)') &
        '[DC-LCFO-SAWF-DMN] rank=',dc%id_tot,' fragment tiling failed: ',trim(message)
      call comm_get_max(failure,dc%icomm_tot)
      if(failure /= 0) call lcfo_sawf_fatal( &
        'SAWF fragment core domains do not tile the global mesh exactly once')

      call get_lattice_vectors(a1,a2,a3)
      lattice(:,1)=a1
      lattice(:,2)=a2
      lattice(:,3)=a3
      lattice_inverse=lattice
      call matrix_inverse(lattice_inverse)
      do ia=1,dc%system_tot%nion
        species(ia)=izatom(dc%system_tot%kion(ia))
        fractional_positions(:,ia)=modulo(matmul(lattice_inverse,dc%system_tot%rion(:,ia)),1.0d0)
      end do

      if(trim(wannier_site_symmetry) == 'auto') then
        call load_sawf_symmetry_auto(lattice,fractional_positions,species, &
          wannier_symmetry_tolerance,symmetry_operations,local_ok,message)
      else
        if(wannier_symmetry_file(1:1) == '/') then
          symmetry_filename=trim(wannier_symmetry_file)
        else
          symmetry_filename=trim(import_run_root_dir())//trim(wannier_symmetry_file)
        end if
        call load_sawf_symmetry_file(symmetry_filename,lattice,fractional_positions,species, &
          wannier_symmetry_tolerance,symmetry_operations,local_ok,message)
      end if
      failure=merge(0,1,local_ok)
      if(failure /= 0) write(*,'(1x,a,i0,2a)') &
        '[DC-LCFO-SAWF-DMN] rank=',dc%id_tot,' symmetry load failed: ',trim(message)
      call comm_get_max(failure,dc%icomm_tot)
      if(failure /= 0) call lcfo_sawf_fatal( &
        'SAWF D_band symmetry loading failed on one or more ranks')

      call put_sawf_identity_first(symmetry_operations,wannier_symmetry_tolerance,local_ok,message)
      failure=merge(0,1,local_ok)
      call comm_get_max(failure,dc%icomm_tot)
      if(failure /= 0) call lcfo_sawf_fatal('SAWF normalized operation set has no unique identity')
      allocate(symmetry_fragment_maps(dc%n_frag,size(symmetry_operations)),stat=allocation_status)
      failure=merge(0,1,allocation_status==0)
      call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) call lcfo_sawf_fatal('SAWF fragment-operation map allocation failed')
      symmetry_fragment_maps=0

      do isym=1,size(symmetry_operations)
        call validate_sawf_actual_group_operation(.true.,.false.,local_ok,message)
        if(.not.local_ok) then
          write(*,'(1x,a,i0,a,i0,2a)') '[DC-LCFO-SAWF-ACTUAL-GROUP] rank=',dc%id_tot, &
            ' operation=',isym,' ',trim(message)
          call lcfo_sawf_fatal('SAWF operation is not in the actual supercell group')
        end if
        call validate_sawf_fragment_symmetry_map(symmetry_operations(isym),mesh, &
          fragment_origin,fragment_shape,dc%nxyz_buffer,wannier_symmetry_tolerance, &
          grid_map_ok,fragment_map_ok,max_targets_per_source,symmetry_fragment_map, &
          max_grid_residual,center_available,center_grid,message)
        local_failure=merge(0,1,grid_map_ok)
        call reduce_sawf_fragment_alignment_failure(local_failure,dc%icomm_tot,dc%id_tot,isym, &
          grid_map_ok,fragment_map_ok,max_targets_per_source,max_grid_residual,message,failure)
        if(failure/=0) then
          if(allocated(symmetry_fragment_map)) deallocate(symmetry_fragment_map)
          call lcfo_sawf_fatal( &
            'SAWF symmetry operation is incompatible with the periodic grid')
        end if
        if(fragment_map_ok) then
          symmetry_fragment_maps(:,isym)=symmetry_fragment_map
        else
          split_fragment_global_mode=.true.
          symmetry_fragment_maps(:,isym)=0
        end if
        if(dc%id_tot==0) then
          if(center_available) then
            write(*,'(1x,a,i0,2(a,l1),a,i0,a,es13.5,a,3f12.5)') &
              '[DC-LCFO-SAWF-ALIGN] operation=',isym, &
              ' grid_map_ok=',grid_map_ok,' fragment_map_ok=',fragment_map_ok, &
              ' max_targets_per_source=',max_targets_per_source, &
              ' max_grid_residual=',max_grid_residual,' center_grid=',center_grid
          else
            write(*,'(1x,a,i0,2(a,l1),a,i0,a,es13.5)') &
              '[DC-LCFO-SAWF-ALIGN] operation=',isym, &
              ' grid_map_ok=',grid_map_ok,' fragment_map_ok=',fragment_map_ok, &
              ' max_targets_per_source=',max_targets_per_source, &
              ' max_grid_residual=',max_grid_residual
          end if
        end if
        if(allocated(symmetry_fragment_map)) deallocate(symmetry_fragment_map)
      end do

      call build_sawf_spd_projection_map(dc%system_tot%nion,channels,local_ok,message,projection_lmax)
      if(local_ok.and.size(channels)/=wannier_num_wann)then
        local_ok=.false.;message='SAWF D_wann channel count differs from num_wann'
      end if
      failure=merge(0,1,local_ok)
      if(failure/=0)write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-CHANNEL] rank=',dc%id_tot,' ',trim(message)
      call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0)call lcfo_sawf_fatal('SAWF projection-channel construction failed')

      if(trim(wannier_sawf_generation)=='hierarchical') then
        allocate(sawf_environment_equivalent(dc%n_frag,dc%n_frag), &
          sawf_defect_intersects(dc%n_frag),sawf_environment_orbit(dc%n_frag), &
          sawf_regenerate_environment(dc%n_frag),sawf_environment_key(dc%n_frag), &
          sawf_representative_fragment(dc%n_frag),sawf_materialize_operation(dc%n_frag), &
          sawf_generate_independently(dc%n_frag),sawf_vacuum_fraction(dc%n_frag),stat=allocation_status)
        if(allocation_status/=0) call lcfo_sawf_fatal( &
          'SAWF hierarchical environment-orbit allocation failed on this rank')
        call build_sawf_fragment_environment_fingerprints(lattice,fractional_positions,species,mesh, &
          fragment_origin,fragment_shape,sawf_environment_key,sawf_vacuum_fraction, &
          sawf_supercell_fingerprint,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-ENV] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF local-environment fingerprint construction failed')
        end if
        sawf_environment_equivalent=.false.; sawf_defect_intersects=.false.
        do ifrag=1,dc%n_frag
          sawf_environment_equivalent(ifrag,ifrag)=.true.
          do isym=1,size(symmetry_operations)
            if(symmetry_fragment_maps(ifrag,isym)>0) then
              sawf_environment_equivalent(ifrag,symmetry_fragment_maps(ifrag,isym))=.true.
              sawf_environment_equivalent(symmetry_fragment_maps(ifrag,isym),ifrag)=.true.
            end if
          end do
        end do
        call build_sawf_environment_orbits(sawf_environment_equivalent,sawf_defect_intersects, &
          sawf_environment_orbit,sawf_regenerate_environment,local_ok,message)
        if(.not.local_ok) then
          write(*,'(1x,a,i0,2a)') '[DC-LCFO-SAWF-ORBIT] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF hierarchical environment-orbit construction failed')
        end if
        call validate_sawf_structure_class(wannier_sawf_structure_class,sawf_environment_key, &
          sawf_vacuum_fraction,sawf_environment_orbit,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-CLASS] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF declared structure class disagrees with measured geometry')
        end if
        call select_sawf_environment_materialization(sawf_environment_orbit, &
          symmetry_fragment_maps,sawf_representative_fragment,sawf_materialize_operation, &
          sawf_generate_independently,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-MATERIALIZE] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF materialization provenance construction failed')
        end if
        call build_sawf_environment_execution_plan(sawf_representative_fragment, &
          sawf_materialize_operation,sawf_generate_independently,sawf_supercell_fingerprint, &
          sawf_environment_receipts,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-PLAN] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF hierarchical execution-plan construction failed')
        end if
        call build_sawf_seed_bundles(sawf_environment_receipts, &
          trim(dc%base_directory)//'sawf-hierarchical',trim(sysname),sawf_seed_bundles,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-SEEDS] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF representative seed-bundle construction failed')
        end if
        call atomic_create_directory(trim(dc%base_directory)//'sawf-hierarchical', &
          dc%icomm_tot,dc%id_tot)
        do ibundle=1,size(sawf_seed_bundles)
          call atomic_create_directory(sawf_seed_bundles(ibundle)%directory,dc%icomm_tot,dc%id_tot)
        end do
        if(sawf_environment_receipts(dc%i_frag)%requires_execution.and.dc%id_frag==0)then
          allocate(sawf_expected_channels(dc%system_tot%nion),sawf_inside_atom(dc%system_tot%nion))
          sawf_expected_channels=0;sawf_inside_atom=.false.
          do ia=1,size(channels)
            sawf_expected_channels(channels(ia)%atom)=sawf_expected_channels(channels(ia)%atom)+1
          end do
          do ia=1,dc%system_tot%nion
            sawf_inside_atom(ia)=sawf_atom_inside_fragment_buffer(fractional_positions(:,ia),mesh, &
              fragment_origin(:,dc%i_frag),fragment_shape(:,dc%i_frag),dc%nxyz_buffer)
          end do
          call select_sawf_local_complete_shells(channels%atom,sawf_expected_channels, &
            sawf_inside_atom,sawf_selected_channels,local_ok,message)
          if(.not.local_ok)then
            write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-LOCAL-SHELL] rank=',dc%id_tot,' ',trim(message)
            call lcfo_sawf_fatal('SAWF local complete-shell selection failed before seed writing')
          end if
        end if
      end if
      call build_sawf_operation_product_table(symmetry_operations,lattice,lattice_inverse, &
        wannier_symmetry_tolerance,product_left,product_right,product_result,local_ok,message)
      if(.not.local_ok) then
        write(*,'(1x,a,i0,2a)') '[DC-LCFO-SAWF-GROUP] rank=',dc%id_tot,' ',trim(message)
        call lcfo_sawf_fatal('SAWF actual-supercell group product table construction failed')
      end if
      if(trim(wannier_sawf_generation)=='hierarchical'.and.dc%id_frag==0.and. &
          sawf_environment_receipts(dc%i_frag)%requires_execution)then
        call select_sawf_environment_stabilizer(dc%i_frag,symmetry_fragment_maps,product_left, &
          product_right,product_result,sawf_local_stabilizer,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-STABILIZER] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF local actual-group stabilizer validation failed')
        end if
      end if
      local_failure=merge(0,1,.not.split_fragment_global_mode)
      call comm_get_max(local_failure,dc%icomm_tot)
      split_fragment_global_mode=(local_failure/=0)
      if(split_fragment_global_mode .and. dc%id_tot==0) write(*,'(1x,a)') &
        '[DC-LCFO-SAWF-GLOBAL-SPLIT] build D_band from split fragment blocks; skip local closed-basis gate'

      failure=0; message=''
      if(dc%id_tot == 0) then
        if(local_ok) then
          amn_filename=trim(dc%base_directory)//trim(sysname)//'.amn'
          call read_sawf_amn_matrix(amn_filename,nband_wann,wannier_num_wann,amn,local_ok,message)
        end if
        if(local_ok) then
          dmn_filename=trim(dc%base_directory)//trim(sysname)//'.dmn'
          call begin_sawf_dmn(writer,dmn_filename,nband_wann,wannier_num_wann, &
            size(symmetry_operations),max(1.0d-10,wannier_symmetry_tolerance),local_ok,message)
        end if
        failure=merge(0,1,local_ok)
      end if
      call comm_bcast(failure,dc%icomm_tot,0)
      call comm_bcast(message,dc%icomm_tot,0)
      if(failure /= 0) then
        if(dc%id_tot == 0) call abort_sawf_dmn(writer)
        call lcfo_sawf_fatal('SAWF DMN initialization failed: '//trim(message))
      end if

      closure_tolerance=max(1.0d-10,wannier_symmetry_tolerance)
      allocate(d_wann_set(wannier_num_wann,wannier_num_wann,size(symmetry_operations)), &
        stat=allocation_status)
      if(allocation_status/=0)call lcfo_sawf_fatal('SAWF D_wann representation-set allocation failed')
      if(dc%id_tot==0) then
        allocate(d_band_set(nband_wann,nband_wann,size(symmetry_operations)),stat=allocation_status)
        if(allocation_status/=0) call lcfo_sawf_fatal('SAWF representation-set allocation failed')
      end if
      call prepare_sawf_fragment_state_cache(nband_wann,fragment_shape,state_cache,local_ok,message)
      failure=merge(0,1,local_ok)
      if(failure/=0) write(*,'(1x,a,i0,2a)') '[DC-LCFO-SAWF-DMN] rank=',dc%id_tot, &
        ' source cache preparation failed: ',trim(message)
      call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) then
        if(dc%id_tot==0) call abort_sawf_dmn(writer)
        call clear_sawf_fragment_state_cache(state_cache)
        call lcfo_sawf_fatal('SAWF source fragment state cache preparation failed')
      end if
      if(trim(wannier_sawf_generation)=='hierarchical'.and.dc%id_frag==0.and. &
          sawf_environment_receipts(dc%i_frag)%requires_execution)then
        call write_sawf_representative_local_seed(state_cache%source,channels, &
          sawf_selected_channels,sawf_seed_bundles,local_ok,message,local_states_out=sawf_local_states, &
          local_energy_out=sawf_local_energy,local_amn_out=sawf_local_amn)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-LOCAL-SEED] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF representative local eigensystem/seed construction failed')
        end if
      end if
      if(trim(wannier_sawf_generation)=='hierarchical')then
        call run_sawf_local_preprocessing(resolved_wannier_command,sawf_seed_bundles)
      end if
      if(trim(wannier_sawf_generation)=='hierarchical'.and. &
          .not.is_wannier90_export_only_command(resolved_wannier_command).and.dc%id_frag==0.and. &
          sawf_environment_receipts(dc%i_frag)%requires_execution)then
        ibundle=findloc(sawf_seed_bundles%environment,dc%i_frag,dim=1)
        if(ibundle<=0)call lcfo_sawf_fatal('SAWF representative bundle lookup failed after preprocessing')
        call read_sawf_nnkp_neighbors(trim(sawf_seed_bundles(ibundle)%directory)//'/'// &
          trim(sawf_seed_bundles(ibundle)%seedname)//'.nnkp',sawf_neighbor_gvec,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-NNKP] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF representative preprocessing neighbor contract failed')
        end if
        call write_sawf_representative_local_seed(state_cache%source,channels, &
          sawf_selected_channels,sawf_seed_bundles,local_ok,message,sawf_neighbor_gvec)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-MMN] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF representative MMN construction failed')
        end if
      end if
      local_ok=.true.; message=''
      if(.not.split_fragment_global_mode) call build_sawf_closed_fragment_seed_basis( &
        nband_wann,mesh,fragment_origin,fragment_shape,symmetry_operations, &
        symmetry_fragment_maps,state_cache,closed_basis,local_ok,message)
      failure=merge(0,1,local_ok)
      if(failure/=0) write(*,'(1x,a,i0,2a)') '[DC-LCFO-SAWF-CLOSED] rank=',dc%id_tot, &
        ' construction failed: ',trim(message)
      call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) then
        if(dc%id_tot==0) call abort_sawf_dmn(writer)
        call clear_sawf_closed_basis(closed_basis)
        call clear_sawf_fragment_state_cache(state_cache)
        call lcfo_sawf_fatal('SAWF symmetry-closed seed basis construction failed')
      end if
      do iop=1,size(symmetry_operations)
        call build_sawf_dmn_operation_fragment_local(nband_wann,mesh,fragment_origin, &
          fragment_shape,iop,symmetry_operations(iop),wannier_symmetry_tolerance, &
          state_cache,d_band_local,local_ok,message)
        failure=merge(0,1,local_ok)
        if(failure /= 0) write(*,'(1x,a,i0,a,i0,2a)') &
          '[DC-LCFO-SAWF-DMN] rank=',dc%id_tot,' operation=',iop, &
          ' fragment-local build failed: ',trim(message)
        call comm_get_max(failure,dc%icomm_tot)
        if(failure /= 0) then
          if(dc%id_tot == 0) call abort_sawf_dmn(writer)
          call clear_sawf_fragment_state_cache(state_cache)
          call lcfo_sawf_fatal( &
            'SAWF D_band fragment-local construction failed on one or more ranks')
        end if
        call reduce_sawf_band_matrix(d_band_local,d_band_sum,nband_wann)
        if(iop==2 .and. dc%id_tot==0) call diagnose_sawf_hamiltonian_covariance_blocks( &
          d_band_sum,esp_tot(1:nband_wann,1),min(128,nband_wann), &
          min(wannier_num_wann,nband_wann))
        if(iop==2) call diagnose_sawf_hamiltonian_component_covariance(d_band_sum,nband_wann)

        failure=0; message=''
        singular_min=0.0d0; singular_max=0.0d0; closure_residual=0.0d0
        if(dc%id_tot == 0) then
          call build_sawf_wannier_representation(symmetry_operations(iop:iop),channels, &
            d_wann_real,local_ok,message)
          if(local_ok) then
            allocate(d_wann(wannier_num_wann,wannier_num_wann),stat=allocation_status)
            if(allocation_status/=0) then
              local_ok=.false.; message='SAWF D_wann complex conversion allocation failed'
            else
              d_wann=cmplx(d_wann_real(:,:,1),0.0d0,kind=8)
              d_band_set(:,:,iop)=d_band_sum
              d_wann_set(:,:,iop)=d_wann
              call append_sawf_dmn_operation(writer,iop,d_wann,d_band_sum, &
                esp_tot(1:nband_wann,1),amn,iop==1,local_ok,message, &
                singular_min,singular_max,closure_residual)
            end if
          end if
          if(allocated(d_wann_real)) deallocate(d_wann_real)
          if(allocated(d_wann)) deallocate(d_wann)
          failure=merge(0,1,local_ok)
        end if
        call comm_bcast(failure,dc%icomm_tot,0)
        call comm_bcast(message,dc%icomm_tot,0)
        call comm_bcast(singular_min,dc%icomm_tot,0)
        call comm_bcast(singular_max,dc%icomm_tot,0)
        call comm_bcast(closure_residual,dc%icomm_tot,0)
        if(dc%id_tot == 0) write(*,'(1x,a,i0,4(a,es13.5))') &
          '[DC-LCFO-SAWF-DMN] operation=',iop,' singular_min=',singular_min, &
          ' singular_max=',singular_max,' closure_residual=',closure_residual, &
          ' tolerance=',closure_tolerance
        if(failure /= 0) then
          if(dc%id_tot == 0) call abort_sawf_dmn(writer)
          call clear_sawf_fragment_state_cache(state_cache)
          call lcfo_sawf_fatal('SAWF DMN operation validation/write failed: '//trim(message))
        end if
        call comm_bcast(d_wann_set(:,:,iop),dc%icomm_tot,0)
      end do

      failure=0; message=''; representation_residual=0d0
      if(dc%id_tot==0) then
        call validate_sawf_operation_set_products(d_band_set,product_left,product_right, &
          product_result,closure_tolerance,representation_residual,local_ok,message)
        if(local_ok) call validate_sawf_operation_set_products(d_wann_set,product_left,product_right, &
          product_result,closure_tolerance,representation_residual,local_ok,message)
        failure=merge(0,1,local_ok)
      end if
      call comm_bcast(failure,dc%icomm_tot,0); call comm_bcast(message,dc%icomm_tot,0)
      call comm_bcast(representation_residual,dc%icomm_tot,0)
      if(failure/=0) then
        if(dc%id_tot==0) call abort_sawf_dmn(writer)
        call lcfo_sawf_fatal('SAWF D_band/D_wann group representation validation failed: '//trim(message))
      end if
      if(dc%id_tot==0) write(*,'(1x,a,es13.5)') &
        '[DC-LCFO-SAWF-GROUP] max_representation_residual=',representation_residual
      if(trim(wannier_sawf_generation)=='hierarchical'.and.dc%id_frag==0.and. &
          sawf_environment_receipts(dc%i_frag)%requires_execution)then
        allocate(sawf_local_d_wann(size(sawf_selected_channels),size(sawf_selected_channels), &
          size(sawf_local_stabilizer)),stat=allocation_status)
        if(allocation_status/=0)call lcfo_sawf_fatal('SAWF local D_wann allocation failed')
        call restrict_sawf_stabilizer_representation(d_wann_set,sawf_selected_channels, &
          sawf_local_stabilizer,closure_tolerance,sawf_local_d_wann,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-LOCAL-DWANN] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF local D_wann restriction failed')
        end if
        allocate(sawf_local_point_map(size(sawf_local_states,1),size(sawf_local_stabilizer)), &
          sawf_local_d_band(size(sawf_local_states,2),size(sawf_local_states,2), &
          size(sawf_local_stabilizer)),stat=allocation_status)
        if(allocation_status/=0)call lcfo_sawf_fatal('SAWF local D_band workspace allocation failed')
        do iop=1,size(sawf_local_stabilizer)
          call build_sawf_fragment_buffer_point_map(symmetry_operations(sawf_local_stabilizer(iop)), &
            mesh,fragment_origin(:,dc%i_frag),fragment_shape(:,dc%i_frag), &
            fragment_origin(:,dc%i_frag),fragment_shape(:,dc%i_frag),dc%nxyz_buffer, &
            wannier_symmetry_tolerance,sawf_point_map_column,local_ok,message)
          if(.not.local_ok.or.size(sawf_point_map_column)/=size(sawf_local_states,1))then
            write(*,'(1x,a,i0,a,i0,2a)')'[DC-LCFO-SAWF-LOCAL-DBAND] rank=',dc%id_tot, &
              ' operation=',sawf_local_stabilizer(iop),' ',trim(message)
            call lcfo_sawf_fatal('SAWF local stabilizer buffer-grid map failed')
          end if
          sawf_local_point_map(:,iop)=sawf_point_map_column;deallocate(sawf_point_map_column)
        end do
        call build_sawf_local_band_representation(sawf_local_states,sawf_local_point_map,hvol, &
          closure_tolerance,sawf_local_d_band,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-LOCAL-DBAND] rank=',dc%id_tot,' ',trim(message)
          call lcfo_sawf_fatal('SAWF local D_band construction failed')
        end if
        allocate(sawf_local_product_left(size(sawf_local_stabilizer)**2), &
          sawf_local_product_right(size(sawf_local_stabilizer)**2), &
          sawf_local_product_result(size(sawf_local_stabilizer)**2))
        local_relation=0
        do local_left=1,size(sawf_local_stabilizer);do local_right=1,size(sawf_local_stabilizer)
          local_relation=local_relation+1;global_relation=0
          do iop=1,size(product_left)
            if(product_left(iop)==sawf_local_stabilizer(local_left).and. &
                product_right(iop)==sawf_local_stabilizer(local_right))then
              global_relation=iop;exit
            end if
          end do
          if(global_relation==0)call lcfo_sawf_fatal('SAWF local stabilizer product lookup failed')
          sawf_local_product_left(local_relation)=local_left
          sawf_local_product_right(local_relation)=local_right
          sawf_local_product_result(local_relation)=findloc(sawf_local_stabilizer== &
            product_result(global_relation),.true.,dim=1)
          if(sawf_local_product_result(local_relation)<=0) &
            call lcfo_sawf_fatal('SAWF local stabilizer product escapes subgroup')
        end do;end do
        call validate_sawf_operation_set_products(sawf_local_d_band,sawf_local_product_left, &
          sawf_local_product_right,sawf_local_product_result,closure_tolerance, &
          representation_residual,local_ok,message)
        if(local_ok)call validate_sawf_operation_set_products(sawf_local_d_wann, &
          sawf_local_product_left,sawf_local_product_right,sawf_local_product_result, &
          closure_tolerance,representation_residual,local_ok,message)
        if(.not.local_ok)then
          write(*,'(1x,a,i0,a,es13.5,2a)')'[DC-LCFO-SAWF-LOCAL-GROUP] rank=',dc%id_tot, &
            ' residual=',representation_residual,' ',trim(message)
          call lcfo_sawf_fatal('SAWF local D_band/D_wann representation closure failed')
        end if
        ibundle=findloc(sawf_seed_bundles%environment,dc%i_frag,dim=1)
        if(ibundle<=0)call lcfo_sawf_fatal('SAWF local DMN bundle lookup failed')
        dmn_filename=trim(sawf_seed_bundles(ibundle)%directory)//'/'// &
          trim(sawf_seed_bundles(ibundle)%seedname)//'.dmn'
        call begin_sawf_dmn(sawf_local_writer,dmn_filename,size(sawf_local_states,2), &
          size(sawf_selected_channels),size(sawf_local_stabilizer),closure_tolerance,local_ok,message)
        if(.not.local_ok)call lcfo_sawf_fatal('SAWF local DMN initialization failed: '//trim(message))
        do iop=1,size(sawf_local_stabilizer)
          call append_sawf_dmn_operation(sawf_local_writer,iop,sawf_local_d_wann(:,:,iop), &
            sawf_local_d_band(:,:,iop),sawf_local_energy,sawf_local_amn, &
            sawf_local_stabilizer(iop)==1,local_ok,message,singular_min,singular_max,closure_residual)
          if(.not.local_ok)then
            call abort_sawf_dmn(sawf_local_writer)
            call lcfo_sawf_fatal('SAWF local DMN operation failed: '//trim(message))
          end if
        end do
        call finish_sawf_dmn(sawf_local_writer,symmetry_operations(sawf_local_stabilizer), &
          local_ok,message)
        if(.not.local_ok)then
          call abort_sawf_dmn(sawf_local_writer)
          call lcfo_sawf_fatal('SAWF local DMN publication failed: '//trim(message))
        end if
        call activate_sawf_win(trim(sawf_seed_bundles(ibundle)%directory)//'/'// &
          trim(sawf_seed_bundles(ibundle)%seedname)//'.win',closure_tolerance,local_ok,message,dc%id_tot)
        if(.not.local_ok)call lcfo_sawf_fatal('SAWF local WIN activation failed: '//trim(message))
      end if
      if(trim(wannier_sawf_generation)=='hierarchical'.and. &
          .not.is_wannier90_export_only_command(resolved_wannier_command)) &
        call run_sawf_local_wannier(resolved_wannier_command,sawf_seed_bundles)
      if(trim(wannier_sawf_generation)=='hierarchical'.and. &
          .not.is_wannier90_export_only_command(resolved_wannier_command))then
        allocate(sawf_receipt_local(dc%n_frag),sawf_receipt_global(dc%n_frag), &
          sawf_bands_local(dc%n_frag),sawf_bands_global(dc%n_frag), &
          sawf_wann_local(dc%n_frag),sawf_wann_global(dc%n_frag))
        sawf_receipt_local=0;sawf_bands_local=0;sawf_wann_local=0
        if(dc%id_frag==0.and.sawf_environment_receipts(dc%i_frag)%requires_execution)then
          ibundle=findloc(sawf_seed_bundles%environment,dc%i_frag,dim=1)
          if(ibundle<=0)call lcfo_sawf_fatal('SAWF receipt bundle lookup failed')
          call write_sawf_atomic_text(trim(sawf_seed_bundles(ibundle)%directory)//'/'// &
            trim(sawf_seed_bundles(ibundle)%seedname)//'.sawf-fingerprint', &
            sawf_seed_bundles(ibundle)%same_supercell_fingerprint,local_ok,message)
          if(.not.local_ok)then
            write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-RECEIPT] rank=',dc%id_tot,' ',trim(message)
            call lcfo_sawf_fatal('SAWF fingerprint publication failed before receipt collective')
          end if
          call complete_sawf_seed_bundle(sawf_seed_bundles(ibundle), &
            sawf_environment_receipts(dc%i_frag),local_ok,message)
          if(.not.local_ok)then
            write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-RECEIPT] rank=',dc%id_tot,' ',trim(message)
            call lcfo_sawf_fatal('SAWF artifact validation failed before receipt collective')
          end if
          sawf_receipt_local(dc%i_frag)=1
          sawf_bands_local(dc%i_frag)=sawf_environment_receipts(dc%i_frag)%num_bands
          sawf_wann_local(dc%i_frag)=sawf_environment_receipts(dc%i_frag)%num_wann
        end if
        call comm_summation(sawf_receipt_local,sawf_receipt_global,dc%n_frag,dc%icomm_tot)
        call comm_summation(sawf_bands_local,sawf_bands_global,dc%n_frag,dc%icomm_tot)
        call comm_summation(sawf_wann_local,sawf_wann_global,dc%n_frag,dc%icomm_tot)
        do ifrag=1,dc%n_frag
          if(.not.sawf_environment_receipts(ifrag)%requires_execution)cycle
          if(sawf_receipt_global(ifrag)/=1.or.sawf_bands_global(ifrag)<=0.or. &
              sawf_wann_global(ifrag)<=0)call lcfo_sawf_fatal( &
            'SAWF representative receipt collective is missing or duplicated')
          sawf_environment_receipts(ifrag)%completed=.true.
          sawf_environment_receipts(ifrag)%num_bands=sawf_bands_global(ifrag)
          sawf_environment_receipts(ifrag)%num_wann=sawf_wann_global(ifrag)
        end do
        call propagate_sawf_representative_receipts(sawf_environment_receipts,local_ok,message)
        if(local_ok)call validate_sawf_environment_receipts(sawf_environment_receipts, &
          sawf_supercell_fingerprint,local_ok,message)
        if(.not.local_ok)call lcfo_sawf_fatal('SAWF collective receipt validation failed: '//trim(message))
        if(dc%id_frag==0.and.sawf_environment_receipts(dc%i_frag)%requires_execution)then
          ibundle=findloc(sawf_seed_bundles%environment,dc%i_frag,dim=1)
          local_chk_filename=trim(sawf_seed_bundles(ibundle)%directory)//'/'// &
            trim(sawf_seed_bundles(ibundle)%seedname)//'.chk'
          call read_wannier90_checkpoint_transform(num_bands_chk,num_wann_chk,sawf_local_v_matrix, &
            sawf_local_centers,sawf_local_spreads,local_chk_filename)
          if(num_bands_chk/=sawf_environment_receipts(dc%i_frag)%num_bands.or. &
              num_wann_chk/=sawf_environment_receipts(dc%i_frag)%num_wann.or. &
              size(sawf_local_states,2)/=num_bands_chk)call lcfo_sawf_fatal( &
            'SAWF local checkpoint dimensions disagree with validated receipt/eigensystem')
          allocate(sawf_representative_orbitals(size(sawf_local_states,1),num_wann_chk))
          sawf_representative_orbitals=matmul(sawf_local_states,sawf_local_v_matrix)
        end if
      end if

      failure=0
      if(dc%id_frag==0) then
        write(*,'(1x,a,i0,5(a,i0))') '[DC-LCFO-SAWF-DMN-CACHE] rank=',dc%id_tot, &
          ' source_seed_reads=',state_cache%source_seed_reads, &
          ' source_reconstructions=',state_cache%source_reconstructions, &
          ' target_seed_reads=',state_cache%target_seed_reads, &
          ' target_reconstructions=',state_cache%target_reconstructions, &
          ' target_cache_hits=',state_cache%target_cache_hits
        if(state_cache%source_seed_reads/=1 .or. state_cache%source_reconstructions/=1) failure=1
      end if
      call comm_get_max(failure,dc%icomm_tot)
      if(failure/=0) then
        if(dc%id_tot==0) call abort_sawf_dmn(writer)
        call clear_sawf_fragment_state_cache(state_cache)
        call lcfo_sawf_fatal('SAWF source fragment was not reconstructed exactly once')
      end if

      failure=0; message=''
      if(dc%id_tot == 0) then
        call finish_sawf_dmn(writer,symmetry_operations,local_ok,message)
        failure=merge(0,1,local_ok)
      end if
      call comm_bcast(failure,dc%icomm_tot,0)
      call comm_bcast(message,dc%icomm_tot,0)
      if(failure /= 0) then
        if(dc%id_tot == 0) call abort_sawf_dmn(writer)
        call clear_sawf_fragment_state_cache(state_cache)
        call lcfo_sawf_fatal('SAWF DMN group validation/publication failed: '//trim(message))
      end if
      if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0,5(a,es13.5))') &
        '[DC-LCFO-SAWF-DMN] published operations=',size(symmetry_operations), &
        ' bands=',nband_wann,' unitarity_max=',writer%max_unitarity, &
        ' hamiltonian_max=',writer%max_hamiltonian,' amn_max=',writer%max_amn, &
        ' group_wann_max=',writer%max_group_wann,' group_band_max=',writer%max_group_band
      call clear_sawf_fragment_state_cache(state_cache)
      call clear_sawf_closed_basis(closed_basis)
      deallocate(species,fractional_positions,fragment_origin,fragment_shape, &
        d_band_local,d_band_sum,symmetry_operations,symmetry_fragment_maps)
      deallocate(channels)
      if(dc%id_tot == 0) deallocate(amn)
    end subroutine generate_sawf_dmn

    subroutine build_sawf_fragment_environment_fingerprints(lattice,fractional_positions,species, &
        mesh,fragment_origin,fragment_shape,environment_key,vacuum_by_environment, &
        supercell_fingerprint,ok,message)
      use salmon_global, only: file_pseudo,xc,wannier_projection,wannier_num_bands,wannier_num_wann, &
        wannier_symmetry_tolerance,wannier_sawf_vacuum_density_threshold
      real(8),intent(in)::lattice(3,3),fractional_positions(:,:)
      integer,intent(in)::species(:),mesh(3),fragment_origin(:,:),fragment_shape(:,:)
      character(256),intent(out)::environment_key(:)
      real(8),intent(out)::vacuum_by_environment(:)
      character(*),intent(out)::supercell_fingerprint
      logical,intent(out)::ok;character(*),intent(out)::message
      integer::ifrag,ia,axis,count,k,idx(3),start(3),extent(3),delta
      integer,allocatable::local_species(:),local_count(:)
      real(8),allocatable::relative(:,:),rho_local(:,:,:),rho_global(:,:,:),density_buffer(:),density_collective(:)
      real(8)::center_fractional(3),vacuum_fraction
      integer::ix,iy,iz,ispin,low_density_count,total_point_count,gx,gy,gz
      logical::inside,pbc(3)
      character(256)::supercell_key
      character(1024)::pseudo_signature
      character(256)::pseudo_digest,pseudo_filename
      character(64)::band_window,shell
      ok=.false.;message='';pbc=.true.;pseudo_signature=''
      if(size(environment_key)/=size(fragment_origin,2).or.size(fragment_shape,2)/=size(environment_key).or. &
          size(vacuum_by_environment)/=size(environment_key))then
        message='SAWF fragment fingerprint dimensions are inconsistent';return
      end if
      do ia=1,min(size(file_pseudo),maxval(species))
        if(len_trim(file_pseudo(ia))>0.and.trim(file_pseudo(ia))/='none')then
          if(file_pseudo(ia)(1:1)=='/')then
            pseudo_filename=trim(file_pseudo(ia))
          else
            pseudo_filename=trim(import_run_root_dir())//trim(file_pseudo(ia))
          end if
          call build_sawf_file_content_digest(pseudo_filename,pseudo_digest,ok,message)
          if(.not.ok)return
          pseudo_signature=trim(pseudo_signature)//'|'//trim(pseudo_digest)
        end if
      end do
      write(band_window,'(a,i0,a,i0)')'bands=',wannier_num_bands,':wann=',wannier_num_wann
      shell=trim(wannier_projection)
      call build_sawf_supercell_fingerprint(lattice,pbc,species,fractional_positions,mesh, &
        dc%nxyz_buffer,trim(pseudo_signature),trim(band_window),trim(shell),trim(xc), &
        'SALMON-SAWF-schema-2',supercell_key,ok,message)
      if(.not.ok)return
      supercell_fingerprint=supercell_key
      allocate(local_species(size(species)),relative(3,size(species)),local_count(size(environment_key)))
      allocate(rho_local(mesh(1),mesh(2),mesh(3)),rho_global(mesh(1),mesh(2),mesh(3)))
      allocate(density_collective(2*product(mesh)))
      rho_local=0d0
      if(info%id_ko==0)then
        do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
          do ispin=1,system%nspin
            rho_local(ix,iy,iz)=rho_local(ix,iy,iz)+rho_s(ispin)%f(ix,iy,iz)
          end do
        end do;end do;end do
      end if
      density_collective(1:product(mesh))=reshape(rho_local,[product(mesh)])
      call validate_sawf_density_contribution(density_collective(1:product(mesh)), &
        info%id_ko==0,dc%id_tot,ok,message)
      if(.not.ok)call lcfo_sawf_fatal(trim(message))
      call assemble_sawf_density_unique(density_collective(1:product(mesh)),info%id_ko==0, &
        dc%icomm_tot,density_collective(product(mesh)+1:2*product(mesh)))
      rho_global=reshape(density_collective(product(mesh)+1:2*product(mesh)),shape(rho_global))
      local_count=0
      do ifrag=1,size(environment_key)
        start=fragment_origin(:,ifrag)-dc%nxyz_buffer;extent=fragment_shape(:,ifrag)+2*dc%nxyz_buffer
        do ia=1,size(species)
          idx=1+modulo(floor(fractional_positions(:,ia)*real(mesh,8)),mesh);inside=.true.
          do axis=1,3
            delta=modulo((idx(axis)-1)-start(axis),mesh(axis));if(delta>=extent(axis))inside=.false.
          end do
          if(inside)local_count(ifrag)=local_count(ifrag)+1
        end do
      end do
      do ifrag=1,size(environment_key)
        count=0;start=fragment_origin(:,ifrag)-dc%nxyz_buffer
        extent=fragment_shape(:,ifrag)+2*dc%nxyz_buffer
        center_fractional=(real(fragment_origin(:,ifrag),8)+0.5d0*real(fragment_shape(:,ifrag),8))/real(mesh,8)
        do ia=1,size(species)
          idx=1+modulo(floor(fractional_positions(:,ia)*real(mesh,8)),mesh)
          inside=.true.
          do axis=1,3
            delta=modulo((idx(axis)-1)-start(axis),mesh(axis))
            if(delta>=extent(axis))inside=.false.
          end do
          if(.not.inside)cycle
          count=count+1;local_species(count)=species(ia)
          relative(:,count)=fractional_positions(:,ia)-center_fractional
          relative(:,count)=relative(:,count)-dnint(relative(:,count))
        end do
        low_density_count=0;total_point_count=product(extent);allocate(density_buffer(total_point_count));k=0
        do iz=0,extent(3)-1;do iy=0,extent(2)-1;do ix=0,extent(1)-1
          gx=1+modulo(start(1)+ix,mesh(1));gy=1+modulo(start(2)+iy,mesh(2))
          gz=1+modulo(start(3)+iz,mesh(3))
          k=k+1;density_buffer(k)=rho_global(gx,gy,gz)
        end do;end do;end do
        call measure_sawf_vacuum_occupancy(density_buffer,wannier_sawf_vacuum_density_threshold, &
          vacuum_fraction,ok,message);deallocate(density_buffer)
        if(.not.ok)return
        vacuum_by_environment(ifrag)=vacuum_fraction
        call build_sawf_local_environment_fingerprint(supercell_key,lattice, &
          wannier_symmetry_tolerance,local_species(:count),relative(:,:count),vacuum_fraction, &
          environment_key(ifrag),ok,message)
        if(.not.ok)return
      end do
      ok=.true.
    end subroutine

    logical function sawf_atom_inside_fragment_buffer(fractional_position,mesh,origin,fragment_shape,buffer) &
        result(inside)
      real(8),intent(in)::fractional_position(3)
      integer,intent(in)::mesh(3),origin(3),fragment_shape(3),buffer(3)
      integer::axis,index(3),delta,extent(3),start(3)
      index=1+modulo(floor(fractional_position*real(mesh,8)),mesh)
      start=origin-buffer;extent=fragment_shape+2*buffer;inside=.true.
      do axis=1,3
        delta=modulo((index(axis)-1)-start(axis),mesh(axis))
        if(delta>=extent(axis))inside=.false.
      end do
    end function sawf_atom_inside_fragment_buffer

    subroutine put_sawf_identity_first(operations,tolerance,ok,message)
      type(t_sawf_symop), intent(inout) :: operations(:)
      real(8), intent(in) :: tolerance
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      type(t_sawf_symop) :: temporary
      integer :: iop,identity_count,identity_index,i
      real(8) :: wrapped(3)
      identity_count=0; identity_index=0
      do iop=1,size(operations)
        wrapped=operations(iop)%tau-anint(operations(iop)%tau)
        if(all(operations(iop)%W==reshape([1,0,0,0,1,0,0,0,1],[3,3])) .and. &
            maxval(abs(wrapped))<=tolerance) then
          identity_count=identity_count+1; identity_index=iop
        end if
      end do
      if(identity_count/=1) then
        ok=.false.; write(message,'(a,i0)') 'SAWF identity operation count=',identity_count; return
      end if
      if(identity_index/=1) then
        temporary=operations(1); operations(1)=operations(identity_index); operations(identity_index)=temporary
      end if
      do i=1,size(operations(1)%atom_map)
        if(operations(1)%atom_map(i)/=i) then
          ok=.false.; message='SAWF identity operation has a non-identity atom map'; return
        end if
      end do
      ok=.true.; message=''
    end subroutine put_sawf_identity_first

    subroutine read_sawf_amn_matrix(filename,expected_bands,expected_wann,amn,ok,message)
      use filesystem, only: get_filehandle
      implicit none
      character(*), intent(in) :: filename
      integer, intent(in) :: expected_bands,expected_wann
      complex(8), allocatable, intent(out) :: amn(:,:)
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      integer :: iunit,io_status,nb,nk,nw,entry,ib,iw,ik,allocation_status
      real(8) :: re,im
      logical, allocatable :: seen(:,:)
      character(512) :: header,io_message
      ok=.false.; message=''
      iunit=get_filehandle()
      open(iunit,file=filename,status='old',action='read',iostat=io_status,iomsg=io_message)
      if(io_status/=0) then; message='SAWF AMN open failed: '//trim(io_message); return; endif
      read(iunit,'(a)',iostat=io_status,iomsg=io_message) header
      if(io_status==0) read(iunit,*,iostat=io_status,iomsg=io_message) nb,nk,nw
      if(io_status/=0 .or. nb/=expected_bands .or. nk/=1 .or. nw/=expected_wann) then
        message='SAWF AMN header does not match current-run Gamma seed'; close(iunit); return
      end if
      allocate(amn(nb,nw),seen(nb,nw),stat=allocation_status)
      if(allocation_status/=0) then; message='SAWF AMN allocation failed'; close(iunit); return; endif
      amn=(0.0d0,0.0d0); seen=.false.
      do entry=1,nb*nw
        read(iunit,*,iostat=io_status,iomsg=io_message) ib,iw,ik,re,im
        if(io_status/=0 .or. ib<1 .or. ib>nb .or. iw<1 .or. iw>nw .or. ik/=1) then
          message='SAWF AMN entry is malformed or out of range'; close(iunit); deallocate(amn,seen); return
        end if
        if(seen(ib,iw) .or. .not.ieee_is_finite(re) .or. .not.ieee_is_finite(im)) then
          message='SAWF AMN contains a duplicate or non-finite entry'; close(iunit); deallocate(amn,seen); return
        end if
        seen(ib,iw)=.true.; amn(ib,iw)=cmplx(re,im,kind=8)
      end do
      read(iunit,*,iostat=io_status)
      if(io_status==0 .or. .not.all(seen)) then
        message='SAWF AMN has trailing or missing entries'; close(iunit); deallocate(amn,seen); return
      end if
      close(iunit); deallocate(seen); ok=.true.
    end subroutine read_sawf_amn_matrix

    subroutine prepare_sawf_fragment_state_cache(nband_wann,fragment_shape,cache,ok,message)
      implicit none
      integer, intent(in) :: nband_wann,fragment_shape(:,:)
      type(t_sawf_fragment_state_cache), intent(inout) :: cache
      logical, intent(out) :: ok
      character(*), intent(out) :: message

      call clear_sawf_fragment_state_cache(cache)
      ok=.true.; message=''
      if(dc%id_frag/=0 .and. dc%n_frag/=1) return
      cache%source_seed_reads=1
      call load_sawf_fragment_state_entry(dc%i_frag,nband_wann,fragment_shape(:,dc%i_frag), &
        cache%source,ok,message)
      if(ok) cache%source_reconstructions=1
    end subroutine prepare_sawf_fragment_state_cache

    subroutine get_sawf_target_cache_slot(target_frag,nband_wann,fragment_shape,cache, &
        use_source,target_slot,ok,message)
      implicit none
      integer, intent(in) :: target_frag,nband_wann,fragment_shape(:,:)
      type(t_sawf_fragment_state_cache), intent(inout) :: cache
      logical, intent(out) :: use_source,ok
      integer, intent(out) :: target_slot
      character(*), intent(out) :: message
      integer :: slot

      ok=.false.; message=''; use_source=.false.; target_slot=0
      if(target_frag==cache%source%fragment) then
        use_source=.true.; ok=.true.; cache%target_cache_hits=cache%target_cache_hits+1
        return
      end if
      do slot=1,sawf_target_cache_capacity
        if(cache%target(slot)%fragment==target_frag .and. allocated(cache%target(slot)%states)) then
          target_slot=slot; ok=.true.; cache%target_cache_hits=cache%target_cache_hits+1
          return
        end if
      end do
      target_slot=cache%next_target_slot
      cache%next_target_slot=1+modulo(cache%next_target_slot,sawf_target_cache_capacity)
      call clear_sawf_fragment_state_entry(cache%target(target_slot))
      cache%target_seed_reads=cache%target_seed_reads+1
      call load_sawf_fragment_state_entry(target_frag,nband_wann,fragment_shape(:,target_frag), &
        cache%target(target_slot),ok,message)
      if(.not.ok) return
      cache%target_reconstructions=cache%target_reconstructions+1
      if(cache%target(target_slot)%nspin/=cache%source%nspin .or. &
          cache%target(target_slot)%nstate_frag/=cache%source%nstate_frag .or. &
          cache%target(target_slot)%nstate_tot/=cache%source%nstate_tot) then
        message='SAWF cached target seed metadata differs from the source fragment'
        call clear_sawf_fragment_state_entry(cache%target(target_slot)); ok=.false.; return
      end if
      ok=.true.
    end subroutine get_sawf_target_cache_slot

    subroutine load_sawf_fragment_state_entry(ifrag,nband_wann,expected_shape,entry,ok,message)
      implicit none
      integer, intent(in) :: ifrag,nband_wann,expected_shape(3)
      type(t_sawf_fragment_state_entry), intent(inout) :: entry
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      integer :: shape_read(3),nspin_file,nstate_frag_file,nstate_tot_file,n_basis_frag
      integer :: shape_buffer_read(3),nxyz_buffer_read(3),nxyz_box_read(3)
      integer :: nspin_buffer,nstate_frag_buffer,nstate_tot_buffer,n_basis_buffer
      integer :: npoints,allocation_status
      real(8), allocatable :: basis(:,:),coef(:,:),buffer_basis(:,:),coef_buffer(:,:)
      logical :: read_ok
      character(512) :: local_message

      call clear_sawf_fragment_state_entry(entry)
      if(sawf_explicit_basis_active) then
        call read_sawf_closed_fragment_state_entry(ifrag,nband_wann,expected_shape,entry,ok,message)
        return
      end if
      call read_fragment_lcfo_seed_for_wannier_import(dc,ifrag,nband_wann,shape_read, &
        nspin_file,nstate_frag_file,nstate_tot_file,n_basis_frag,basis,coef,read_ok,.true.,local_message)
      if(.not.read_ok .or. any(shape_read/=expected_shape) .or. nspin_file/=1 .or. &
          nstate_tot_file<nband_wann .or. n_basis_frag<1 .or. n_basis_frag>nstate_frag_file) then
        write(message,'(a,i0,2a)') 'failed to read consistent current-run LCFO fragment ', &
          ifrag,': ',trim(local_message)
        if(allocated(basis)) deallocate(basis)
        if(allocated(coef)) deallocate(coef)
        ok=.false.; return
      end if
      call checked_sawf_fragment_point_count(shape_read,npoints,ok,local_message)
      if(.not.ok) then
        message=trim(local_message); deallocate(basis,coef); return
      end if
      call read_fragment_lcfo_buffer_seed_for_wannier_import(dc,ifrag,nband_wann, &
        shape_buffer_read,nxyz_buffer_read,nxyz_box_read,nspin_buffer,nstate_frag_buffer, &
        nstate_tot_buffer,n_basis_buffer,buffer_basis,coef_buffer,read_ok)
      if(.not.read_ok .or. any(shape_buffer_read/=shape_read) .or. &
          any(nxyz_buffer_read/=dc%nxyz_buffer) .or. nspin_buffer/=nspin_file .or. &
          nstate_frag_buffer/=nstate_frag_file .or. nstate_tot_buffer/=nstate_tot_file .or. &
          n_basis_buffer/=n_basis_frag) then
        message='SAWF buffered fragment seed metadata differs from the current core seed'
        if(allocated(buffer_basis)) deallocate(buffer_basis)
        if(allocated(coef_buffer)) deallocate(coef_buffer)
        deallocate(basis,coef); ok=.false.; return
      end if
      allocate(entry%basis(npoints,n_basis_frag),entry%states(npoints,nband_wann), &
        entry%buffer_basis(product(nxyz_box_read),n_basis_frag), &
        stat=allocation_status)
      if(allocation_status/=0) then
        call clear_sawf_fragment_state_entry(entry)
        message='SAWF cached fragment-state allocation failed'
        deallocate(basis,coef,buffer_basis,coef_buffer); ok=.false.; return
      end if
      entry%basis=basis(:,1:n_basis_frag)
      entry%buffer_basis=buffer_basis(:,1:n_basis_frag)
      entry%states=matmul(cmplx(basis(:,1:n_basis_frag),0.0d0,kind=8), &
        cmplx(coef(1:n_basis_frag,1:nband_wann),0.0d0,kind=8))
      entry%fragment=ifrag; entry%shape=shape_read; entry%buffer_shape=nxyz_box_read
      entry%buffer_width=nxyz_buffer_read; entry%nspin=nspin_file
      entry%nstate_frag=nstate_frag_file; entry%nstate_tot=nstate_tot_file; entry%n_basis=n_basis_frag
      deallocate(basis,coef,buffer_basis,coef_buffer); ok=.true.; message=''
    end subroutine load_sawf_fragment_state_entry

    subroutine read_sawf_closed_fragment_state_entry(ifrag,nband_wann,expected_shape,entry,ok,message)
      use filesystem, only: get_filehandle
      implicit none
      integer, intent(in) :: ifrag,nband_wann,expected_shape(3)
      type(t_sawf_fragment_state_entry), intent(inout) :: entry
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      integer :: iunit,io,magic_file,version_file,shape_file(3),buffer_file(3)
      integer :: nspin_file,nstate_frag_file,nstate_tot_file,n_basis_file
      integer :: npoint_core,npoint_buffer,allocation_status
      character(64) :: token_file
      character(256) :: filename
      real(8), allocatable :: coef(:,:),esp(:)

      ok=.false.; message=''
      write(filename,'(a,a,i6.6,a,a)') trim(import_run_root_dir()), &
        'data_dcdft/fragments/',ifrag,'/','sawf_closed_seed.bin'
      iunit=get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='old',iostat=io)
      if(io/=0) then; message='SAWF closed seed file is missing'; return; end if
      read(iunit,iostat=io) magic_file,version_file
      if(io==0) read(iunit,iostat=io) token_file
      if(io==0) read(iunit,iostat=io) shape_file,buffer_file,nspin_file,nstate_frag_file,nstate_tot_file
      if(io==0) read(iunit,iostat=io) n_basis_file
      if(io/=0 .or. magic_file/=sawf_closed_seed_magic .or. &
          version_file/=sawf_closed_seed_version .or. trim(token_file)/=trim(current_sawf_seed_token)) then
        close(iunit); message='SAWF closed seed header or current-run token is invalid'; return
      end if
      if(any(shape_file/=expected_shape) .or. any(buffer_file/=dc%nxyz_buffer) .or. &
          nspin_file/=1 .or. nstate_frag_file<1 .or. nstate_tot_file<nband_wann .or. &
          n_basis_file<1 .or. n_basis_file>nstate_frag_file) then
        close(iunit); message='SAWF closed seed dimensions are inconsistent'; return
      end if
      npoint_core=product(shape_file); npoint_buffer=product(shape_file+2*buffer_file)
      allocate(entry%basis(npoint_core,n_basis_file), &
        entry%buffer_basis(npoint_buffer,n_basis_file),entry%states(npoint_core,nband_wann), &
        coef(n_basis_file,nstate_tot_file),esp(nstate_tot_file),stat=allocation_status)
      if(allocation_status/=0) then
        close(iunit); call clear_sawf_fragment_state_entry(entry)
        message='SAWF closed seed cache allocation failed'; return
      end if
      read(iunit,iostat=io) entry%basis
      if(io==0) read(iunit,iostat=io) entry%buffer_basis
      if(io==0) read(iunit,iostat=io) coef
      if(io==0) read(iunit,iostat=io) esp
      close(iunit)
      if(io/=0 .or. .not.all(ieee_is_finite(entry%basis)) .or. &
          .not.all(ieee_is_finite(entry%buffer_basis)) .or. .not.all(ieee_is_finite(coef))) then
        call clear_sawf_fragment_state_entry(entry); deallocate(coef,esp)
        message='SAWF closed seed payload is truncated or non-finite'; return
      end if
      entry%states=matmul(cmplx(entry%basis,0.0d0,kind=8), &
        cmplx(coef(:,1:nband_wann),0.0d0,kind=8))
      entry%fragment=ifrag; entry%shape=shape_file; entry%buffer_width=buffer_file
      entry%buffer_shape=shape_file+2*buffer_file; entry%nspin=nspin_file
      entry%nstate_frag=nstate_frag_file; entry%nstate_tot=nstate_tot_file; entry%n_basis=n_basis_file
      deallocate(coef,esp); ok=.true.; message=''
    end subroutine read_sawf_closed_fragment_state_entry

    subroutine clear_sawf_fragment_state_entry(entry)
      implicit none
      type(t_sawf_fragment_state_entry), intent(inout) :: entry
      if(allocated(entry%basis)) deallocate(entry%basis)
      if(allocated(entry%buffer_basis)) deallocate(entry%buffer_basis)
      if(allocated(entry%states)) deallocate(entry%states)
      entry%fragment=0; entry%shape=0; entry%buffer_shape=0; entry%buffer_width=0
      entry%nspin=0; entry%nstate_frag=0
      entry%nstate_tot=0; entry%n_basis=0
    end subroutine clear_sawf_fragment_state_entry

    subroutine clear_sawf_fragment_state_cache(cache)
      implicit none
      type(t_sawf_fragment_state_cache), intent(inout) :: cache
      integer :: slot
      call clear_sawf_fragment_state_entry(cache%source)
      do slot=1,sawf_target_cache_capacity
        call clear_sawf_fragment_state_entry(cache%target(slot))
      end do
      cache%next_target_slot=1; cache%source_seed_reads=0; cache%source_reconstructions=0
      cache%target_seed_reads=0; cache%target_reconstructions=0; cache%target_cache_hits=0
    end subroutine clear_sawf_fragment_state_cache

    subroutine build_sawf_closed_fragment_seed_basis(nband_wann,mesh,fragment_origin, &
        fragment_shape,operations,fragment_maps,cache,closed,ok,message)
      use salmon_global, only: wannier_symmetry_tolerance
      implicit none
      integer, intent(in) :: nband_wann,mesh(3),fragment_origin(:,:),fragment_shape(:,:)
      type(t_sawf_symop), intent(in) :: operations(:)
      integer, intent(in) :: fragment_maps(:,:)
      type(t_sawf_fragment_state_cache), intent(inout) :: cache
      type(t_sawf_closed_basis), intent(inout) :: closed
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      real(8), allocatable :: core_candidate(:,:),buffer_candidate(:,:)
      integer, allocatable :: core_map(:),buffer_map(:)
      integer :: iop,source_frag,target_frag,target_slot,ncandidate,column_first
      integer :: npoint_core,npoint_buffer,allocation_status
      logical :: use_source,local_ok
      character(512) :: detail

      call clear_sawf_closed_basis(closed)
      ok=.true.; message=''
      if(dc%id_frag/=0) return
      target_frag=dc%i_frag
      if(size(fragment_maps,1)/=dc%n_frag .or. size(fragment_maps,2)/=size(operations)) then
        message='SAWF closed-basis fragment-map dimensions are inconsistent'
        ok=.false.; return
      end if
      ncandidate=0
      do iop=1,size(operations)
        if(count(fragment_maps(:,iop)==target_frag)/=1) then
          message='SAWF closed-basis operation has no unique inverse source fragment'
          ok=.false.; return
        end if
        source_frag=findloc(fragment_maps(:,iop)==dc%i_frag,.true.,dim=1)
        call get_sawf_target_cache_slot(source_frag,nband_wann,fragment_shape,cache, &
          use_source,target_slot,local_ok,detail)
        if(.not.local_ok) then; message=trim(detail); ok=.false.; return; end if
        if(use_source) then
          ncandidate=ncandidate+cache%source%n_basis
        else
          ncandidate=ncandidate+cache%target(target_slot)%n_basis
        end if
      end do
      npoint_core=product(fragment_shape(:,target_frag))
      npoint_buffer=product(fragment_shape(:,target_frag)+2*dc%nxyz_buffer)
      allocate(core_candidate(npoint_core,ncandidate), &
        buffer_candidate(npoint_buffer,ncandidate),stat=allocation_status)
      if(allocation_status/=0) then
        message='SAWF closed-basis candidate allocation failed'
        ok=.false.; return
      end if
      core_candidate=0.0d0; buffer_candidate=0.0d0; column_first=1
      do iop=1,size(operations)
        source_frag=findloc(fragment_maps(:,iop)==dc%i_frag,.true.,dim=1)
        call get_sawf_target_cache_slot(source_frag,nband_wann,fragment_shape,cache, &
          use_source,target_slot,local_ok,detail)
        if(.not.local_ok) then
          message=trim(detail); deallocate(core_candidate,buffer_candidate); ok=.false.; return
        end if
        call build_sawf_fragment_buffer_point_map(operations(iop),mesh, &
          fragment_origin(:,source_frag),fragment_shape(:,source_frag), &
          fragment_origin(:,target_frag),fragment_shape(:,target_frag),[0,0,0], &
          wannier_symmetry_tolerance,core_map,local_ok,detail)
        if(.not.local_ok) then
          message=trim(detail); deallocate(core_candidate,buffer_candidate); ok=.false.; return
        end if
        call build_sawf_fragment_buffer_point_map(operations(iop),mesh, &
          fragment_origin(:,source_frag),fragment_shape(:,source_frag), &
          fragment_origin(:,target_frag),fragment_shape(:,target_frag),dc%nxyz_buffer, &
          wannier_symmetry_tolerance,buffer_map,local_ok,detail)
        if(.not.local_ok) then
          message=trim(detail); deallocate(core_map,core_candidate,buffer_candidate); ok=.false.; return
        end if
        if(use_source) then
          associate(source_entry=>cache%source)
            call append_sawf_mapped_basis(source_entry%basis,core_map,core_candidate, &
              column_first,local_ok,detail)
            if(local_ok) call append_sawf_mapped_basis(source_entry%buffer_basis,buffer_map, &
              buffer_candidate,column_first,local_ok,detail)
            column_first=column_first+source_entry%n_basis
          end associate
        else
          associate(source_entry=>cache%target(target_slot))
            call append_sawf_mapped_basis(source_entry%basis,core_map,core_candidate, &
              column_first,local_ok,detail)
            if(local_ok) call append_sawf_mapped_basis(source_entry%buffer_basis,buffer_map, &
              buffer_candidate,column_first,local_ok,detail)
            column_first=column_first+source_entry%n_basis
          end associate
        end if
        deallocate(core_map,buffer_map)
        if(.not.local_ok) then
          message=trim(detail); deallocate(core_candidate,buffer_candidate); ok=.false.; return
        end if
      end do
      call build_sawf_closed_core_buffer_basis(core_candidate,buffer_candidate,hvol, &
        max(1.0d-10,wannier_symmetry_tolerance),npoint_core,closed,local_ok,detail)
      deallocate(core_candidate,buffer_candidate)
      if(.not.local_ok) then; message=trim(detail); ok=.false.; return; end if
      write(*,'(1x,a,i0,4(a,i0),2(a,es13.5))') '[DC-LCFO-SAWF-CLOSED] rank=',dc%id_tot, &
        ' fragment=',target_frag,' candidates=',closed%ncandidate,' basis=',closed%nbasis, &
        ' buffer_points=',closed%npoint_buffer,' singular_max=',closed%singular_values(1), &
        ' singular_min=',closed%singular_values(closed%nbasis)
      ok=.true.; message=''
    end subroutine build_sawf_closed_fragment_seed_basis

    subroutine build_sawf_dmn_operation_fragment_local(nband_wann,mesh,fragment_origin, &
        fragment_shape,operation_index,operation,grid_tolerance,cache,d_band_local,ok,message)
      implicit none
      integer, intent(in) :: nband_wann,mesh(3),fragment_origin(:,:),fragment_shape(:,:),operation_index
      type(t_sawf_symop), intent(in) :: operation
      real(8), intent(in) :: grid_tolerance
      type(t_sawf_fragment_state_cache), intent(inout) :: cache
      complex(8), intent(out) :: d_band_local(nband_wann,nband_wann)
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      integer :: source_frag,target_frag,nsource,nmapped,npair,target_slot,expected_local_points
      integer, parameter :: histogram_entry_cap=32
      integer :: source_shape(3),source_global(3),target_global(3)
      integer :: target_owner,target_local(3),ix,iy,iz,p,allocation_status
      integer :: listed_entries,nonzero_entries,truncated_entries
      integer, allocatable :: source_points(:),target_points(:),target_owner_list(:), &
        target_local_flat(:),target_histogram(:)
      complex(8), allocatable :: block_contribution(:,:)
      logical :: map_ok,use_source
      character(512) :: local_message
      character(1024) :: histogram_text
      character(64) :: histogram_entry
      real(8) :: block_leakage,max_block_leakage,aggregate_leakage
      real(8) :: block_residual_norm2,block_transformed_norm2
      real(8) :: residual_norm2_sum,transformed_norm2_sum

      ok=.false.; message=''; d_band_local=(0.0d0,0.0d0)
      if(dc%id_frag/=0 .and. dc%n_frag/=1) then
        ok=.true.
        return
      end if
      source_frag=cache%source%fragment
      if(source_frag/=dc%i_frag .or. .not.allocated(cache%source%states)) then
        message='SAWF source fragment state cache is not prepared'
        return
      end if
      source_shape=cache%source%shape; nsource=size(cache%source%states,1)
      allocate(source_points(nsource),target_points(nsource), &
        target_owner_list(nsource),target_local_flat(nsource), &
        target_histogram(dc%n_frag),block_contribution(nband_wann,nband_wann), &
        stat=allocation_status)
      if(allocation_status /= 0) then
        message='source fragment D_band allocation failed'
        return
      end if
      p=0
      do iz=1,source_shape(3)
        do iy=1,source_shape(2)
          do ix=1,source_shape(1)
            p=p+1
            source_global=1+modulo(fragment_origin(:,source_frag)+[ix-1,iy-1,iz-1],mesh)
            call map_sawf_periodic_grid_point(operation,mesh,grid_tolerance, &
              source_global,target_global,map_ok,local_message)
            if(.not.map_ok) then
              message=trim(local_message)
              return
            end if
            call locate_sawf_fragment_point(target_global,mesh,fragment_origin,fragment_shape, &
              target_owner,target_local,map_ok,local_message)
            if(.not.map_ok) then
              message=trim(local_message)
              return
            end if
            target_owner_list(p)=target_owner
            target_local_flat(p)=target_local(1)+(target_local(2)-1)*fragment_shape(1,target_owner) + &
              (target_local(3)-1)*fragment_shape(1,target_owner)*fragment_shape(2,target_owner)
          end do
        end do
      end do
      nmapped=0
      target_histogram=0
      do p=1,nsource
        if(dc%n_frag==1 .and. modulo(p-1,dc%isize_frag)/=dc%id_frag) cycle
        target_histogram(target_owner_list(p))=target_histogram(target_owner_list(p))+1
        nmapped=nmapped+1
      end do
      max_block_leakage=0.0d0
      residual_norm2_sum=0.0d0
      transformed_norm2_sum=0.0d0
      do target_frag=1,dc%n_frag
        npair=target_histogram(target_frag)
        if(npair == 0) cycle
        call get_sawf_target_cache_slot(target_frag,nband_wann,fragment_shape,cache, &
          use_source,target_slot,map_ok,local_message)
        if(.not.map_ok) then; message=trim(local_message); return; end if
        npair=0
        do p=1,nsource
          if(dc%n_frag==1 .and. modulo(p-1,dc%isize_frag)/=dc%id_frag) cycle
          if(target_owner_list(p) /= target_frag) cycle
          npair=npair+1
          source_points(npair)=p
          target_points(npair)=target_local_flat(p)
        end do
        if(use_source) then
          call diagnose_sawf_fragment_basis_closure(cache%source%basis,cache%source%basis, &
            source_points(1:npair),target_points(1:npair),dc%system_tot%hvol, &
            block_leakage,block_residual_norm2,block_transformed_norm2,map_ok,local_message)
          if(.not.map_ok) then; message=trim(local_message); return; end if
          call accumulate_sawf_dmn_band_blocks(cache%source%states,cache%source%states, &
            source_points(1:npair),target_points(1:npair),dc%system_tot%hvol, &
            block_contribution,map_ok,local_message)
        else
          call diagnose_sawf_fragment_basis_closure(cache%source%basis, &
            cache%target(target_slot)%basis,source_points(1:npair),target_points(1:npair), &
            dc%system_tot%hvol,block_leakage,block_residual_norm2,block_transformed_norm2, &
            map_ok,local_message)
          if(.not.map_ok) then; message=trim(local_message); return; end if
          call accumulate_sawf_dmn_band_blocks(cache%source%states,cache%target(target_slot)%states, &
            source_points(1:npair),target_points(1:npair),dc%system_tot%hvol, &
            block_contribution,map_ok,local_message)
        end if
        if(.not.map_ok) then
          message=trim(local_message)
          return
        end if
        max_block_leakage=max(max_block_leakage,block_leakage)
        residual_norm2_sum=residual_norm2_sum+block_residual_norm2
        transformed_norm2_sum=transformed_norm2_sum+block_transformed_norm2
        d_band_local=d_band_local+block_contribution
      end do
      expected_local_points=nsource
      if(dc%n_frag==1) expected_local_points=(nsource+dc%isize_frag-1-dc%id_frag)/dc%isize_frag
      if(nmapped /= expected_local_points) then
        write(message,'(a,i0,a,i0)') 'fragment symmetry map count=',nmapped, &
          ' expected=',expected_local_points
        return
      end if
      aggregate_leakage=residual_norm2_sum/transformed_norm2_sum
      histogram_text=''; listed_entries=0
      nonzero_entries=count(target_histogram>0)
      do target_frag=1,dc%n_frag
        if(target_histogram(target_frag)>0 .and. listed_entries<histogram_entry_cap) then
          write(histogram_entry,'(i0,a,i0)') target_frag,':',target_histogram(target_frag)
          if(listed_entries==0) then
            histogram_text=trim(histogram_entry)
          else
            histogram_text=trim(histogram_text)//','//trim(histogram_entry)
          end if
          listed_entries=listed_entries+1
        end if
      end do
      truncated_entries=nonzero_entries-listed_entries
      if(truncated_entries>0) then
        write(histogram_entry,'(a,i0)') 'truncated:',truncated_entries
        histogram_text=trim(histogram_text)//','//trim(histogram_entry)
      end if
      write(*,'(1x,a,i0,a,i0,a,i0,3a,es13.5,3(a,es13.5))') &
        '[DC-LCFO-SAWF-CLOSURE-LOCAL] rank=',dc%id_tot,' operation=',operation_index, &
        ' source_fragment=',source_frag,' histogram=',trim(histogram_text), &
        ' aggregate_leakage=',aggregate_leakage,' max_block_leakage=',max_block_leakage, &
        ' residual_norm2=',residual_norm2_sum,' transformed_norm2=',transformed_norm2_sum
      deallocate(source_points,target_points, &
        target_owner_list,target_local_flat,target_histogram,block_contribution)
      ok=.true.
    end subroutine build_sawf_dmn_operation_fragment_local

    subroutine diagnose_sawf_hamiltonian_covariance_blocks(d_band,eigenvalues,noccupied,nwann)
      implicit none
      complex(8), intent(in) :: d_band(:,:)
      real(8), intent(in) :: eigenvalues(:)
      integer, intent(in) :: noccupied,nwann
      complex(8), allocatable :: weighted(:,:),covariance(:,:)
      real(8), allocatable :: centered(:)
      real(8) :: center,scale,residual,leakage
      integer :: i,nblock

      allocate(weighted(size(d_band,1),size(d_band,2)), &
        covariance(size(d_band,1),size(d_band,2)),centered(size(eigenvalues)))
      center=minval(eigenvalues)+0.5d0*(maxval(eigenvalues)-minval(eigenvalues))
      centered=eigenvalues-center
      weighted=d_band
      weighted=spread(centered,2,size(d_band,2))*weighted
      covariance=matmul(conjg(transpose(d_band)),weighted)
      do i=1,size(centered); covariance(i,i)=covariance(i,i)-centered(i); end do
      scale=max(1.0d-300,maxval(abs(centered)))
      nblock=noccupied
      residual=maxval(abs(covariance(1:nblock,1:nblock)))/scale
      leakage=sqrt(sum(abs(d_band(nblock+1:,1:nblock))**2)/dble(nblock))
      write(*,'(1x,a,a,a,i0,2(a,es13.5))') &
        '[DC-LCFO-SAWF-H-COVARIANCE-BLOCK] label=','occupied',' bands=',nblock, &
        ' relative=',residual,' outside_leakage=',leakage
      nblock=nwann
      residual=maxval(abs(covariance(1:nblock,1:nblock)))/scale
      leakage=sqrt(sum(abs(d_band(nblock+1:,1:nblock))**2)/dble(nblock))
      write(*,'(1x,a,a,a,i0,2(a,es13.5))') &
        '[DC-LCFO-SAWF-H-COVARIANCE-BLOCK] label=','wannier',' bands=',nblock, &
        ' relative=',residual,' outside_leakage=',leakage
      nblock=size(d_band,1)
      residual=maxval(abs(covariance))/scale
      write(*,'(1x,a,a,a,i0,2(a,es13.5))') &
        '[DC-LCFO-SAWF-H-COVARIANCE-BLOCK] label=','all',' bands=',nblock, &
        ' relative=',residual,' outside_leakage=',0.0d0
      deallocate(weighted,covariance,centered)
    end subroutine diagnose_sawf_hamiltonian_covariance_blocks

    subroutine diagnose_sawf_hamiltonian_component_covariance(d_band,nband_wann)
      use communication, only: comm_summation
      implicit none
      integer, intent(in) :: nband_wann
      complex(8), intent(in) :: d_band(nband_wann,nband_wann)
      integer :: n,ifrag,io,jo,i_halo,ig,jg,allocation_status
      real(8), allocatable :: c_local(:,:),c_global(:,:),h_local(:,:),h_global(:,:)

      n=n_mat(1)
      allocate(c_local(n,nband_wann),c_global(n,nband_wann), &
        h_local(n,n),h_global(n,n),stat=allocation_status)
      if(allocation_status/=0) call lcfo_sawf_fatal('SAWF H covariance diagnostic allocation failed')
      c_local=0.0d0
      if(dc%id_frag==0 .and. allocated(coef_wf)) then
        ifrag=dc%i_frag
        do io=1,n_basis(ifrag,1)
          ig=index_basis(io,ifrag,1)
          c_local(ig,:)=coef_wf(io,1:nband_wann,1)
        end do
      end if
      call comm_summation(c_local,c_global,n*nband_wann,dc%icomm_tot)

      call diagnose_sawf_diagonal_h_component(mat_H_weak_kinetic,h_local,h_global,c_global, &
        d_band,nband_wann,label='volume_kinetic')
      call diagnose_sawf_diagonal_h_component(mat_H_weak_potential,h_local,h_global,c_global, &
        d_band,nband_wann,label='volume_local')
      call diagnose_sawf_diagonal_h_component(mat_H_weak_nonlocal,h_local,h_global,c_global, &
        d_band,nband_wann,label='volume_nonlocal')

      h_local=0.0d0
      if(dc%id_frag==0) then
        ifrag=dc%i_frag
        do jo=1,n_basis(ifrag,1)
          jg=index_basis(jo,ifrag,1)
          do io=1,n_basis(ifrag,1)
            ig=index_basis(io,ifrag,1)
            if(use_weak_volume_hamiltonian_mode()) then
              h_local(ig,jg)=mat_H_volume_weak_local(io,jo,1)
            else
              h_local(ig,jg)=mat_H_volume_local(io,jo,1)
            end if
          end do
        end do
      end if
      call diagnose_sawf_h_component(h_local,h_global,c_global,d_band,nband_wann,label='volume')

      h_local=0.0d0
      if(dc%id_frag==0 .and. use_surface_self_hamiltonian_mode()) then
        ifrag=dc%i_frag
        do jo=1,n_basis(ifrag,1)
          jg=index_basis(jo,ifrag,1)
          do io=1,n_basis(ifrag,1)
            ig=index_basis(io,ifrag,1)
            h_local(ig,jg)=mat_H_surface_self(io,jo,1)
          end do
        end do
      end if
      call diagnose_sawf_h_component(h_local,h_global,c_global,d_band,nband_wann,label='surface_self')

      h_local=0.0d0
      if(dc%id_frag==0) then
        ifrag=dc%i_frag
        do i_halo=1,n_halo
          do jo=1,n_basis(ifrag,1)
            jg=index_basis(jo,ifrag,1)
            do io=1,n_basis(halo(i_halo)%ifrag_src,1)
              ig=index_basis(io,halo(i_halo)%ifrag_src,1)
              h_local(ig,jg)=h_local(ig,jg)+0.5d0*halo(i_halo)%mat_H_surface_cross(io,jo,1)
              h_local(jg,ig)=h_local(jg,ig)+0.5d0*halo(i_halo)%mat_H_surface_cross(io,jo,1)
            end do
          end do
        end do
      end if
      call diagnose_sawf_h_component(h_local,h_global,c_global,d_band,nband_wann,label='surface_cross')
      deallocate(c_local,c_global,h_local,h_global)
    end subroutine diagnose_sawf_hamiltonian_component_covariance

    subroutine diagnose_sawf_diagonal_h_component(component,local_matrix,global_matrix, &
        eigenvectors,representation,nband_wann,label)
      implicit none
      integer, intent(in) :: nband_wann
      real(8), intent(in) :: component(:,:,:),eigenvectors(:,:)
      real(8), intent(inout) :: local_matrix(:,:),global_matrix(:,:)
      complex(8), intent(in) :: representation(:,:)
      character(*), intent(in) :: label
      integer :: io,jo,ig,jg,ifrag

      local_matrix=0.0d0
      if(dc%id_frag==0) then
        ifrag=dc%i_frag
        do jo=1,n_basis(ifrag,1)
          jg=index_basis(jo,ifrag,1)
          do io=1,n_basis(ifrag,1)
            ig=index_basis(io,ifrag,1)
            local_matrix(ig,jg)=component(io,jo,1)
          end do
        end do
      end if
      call diagnose_sawf_h_component(local_matrix,global_matrix,eigenvectors, &
        representation,nband_wann,label)
    end subroutine diagnose_sawf_diagonal_h_component

    subroutine diagnose_sawf_h_component(local_matrix,global_matrix,eigenvectors, &
        representation,nband_wann,label)
      use communication, only: comm_summation
      implicit none
      integer, intent(in) :: nband_wann
      real(8), intent(in) :: local_matrix(:,:),eigenvectors(:,:)
      real(8), intent(out) :: global_matrix(:,:)
      complex(8), intent(in) :: representation(:,:)
      character(*), intent(in) :: label
      real(8), allocatable :: work(:,:),band_matrix(:,:)
      complex(8), allocatable :: transformed(:,:)
      real(8) :: scale,residual

      call comm_summation(local_matrix,global_matrix,size(local_matrix),dc%icomm_tot)
      if(dc%id_tot/=0) return
      allocate(work(size(global_matrix,1),nband_wann),band_matrix(nband_wann,nband_wann), &
        transformed(nband_wann,nband_wann))
      work=matmul(global_matrix,eigenvectors)
      band_matrix=matmul(transpose(eigenvectors),work)
      transformed=matmul(conjg(transpose(representation)), &
        matmul(cmplx(band_matrix,0.0d0,kind=8),representation))-cmplx(band_matrix,0.0d0,kind=8)
      scale=max(1.0d-300,maxval(abs(band_matrix)))
      residual=maxval(abs(transformed))/scale
      write(*,'(1x,3a,es13.5,a,es13.5)') &
        '[DC-LCFO-SAWF-H-COVARIANCE] label=',trim(label),' relative=',residual,' scale=',scale
      deallocate(work,band_matrix,transformed)
    end subroutine diagnose_sawf_h_component

    subroutine checked_sawf_fragment_point_count(fragment_shape,npoints,ok,message)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      integer, intent(in) :: fragment_shape(3)
      integer, intent(out) :: npoints
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      integer(int64) :: product64
      integer :: axis

      ok=.false.; message=''; npoints=0
      if(any(fragment_shape <= 0)) then
        message='SAWF fragment shape contains a nonpositive dimension'
        return
      end if
      product64=1_int64
      do axis=1,3
        if(product64 > int(huge(0),int64)/int(fragment_shape(axis),int64)) then
          message='SAWF fragment point count overflows default integer'
          return
        end if
        product64=product64*int(fragment_shape(axis),int64)
      end do
      npoints=int(product64)
      ok=.true.
    end subroutine checked_sawf_fragment_point_count

    subroutine write_sawf_representative_local_seed(entry,projection_channels,selected_channels, &
        bundles,ok,message,neighbor_gvec,local_states_out,local_energy_out,local_amn_out)
      use salmon_global, only: wannier_projection_width,wannier_num_iter
      use salmon_math, only: matrix_inverse
      type(t_sawf_fragment_state_entry),intent(in)::entry
      type(t_sawf_projection_channel),intent(in)::projection_channels(:)
      integer,intent(in)::selected_channels(:)
      type(t_sawf_seed_bundle),intent(in)::bundles(:)
      logical,intent(out)::ok
      character(*),intent(out)::message
      integer,intent(in),optional::neighbor_gvec(:,:)
      complex(8),allocatable,intent(out),optional::local_states_out(:,:)
      real(8),allocatable,intent(out),optional::local_energy_out(:)
      complex(8),allocatable,intent(out),optional::local_amn_out(:,:)
      real(8),parameter::hartree_to_ev=27.211386245988d0
      real(8),allocatable::h_basis(:,:),energy_hartree(:),energy_ev(:)
      complex(8),allocatable::states(:,:),projection(:,:),phase_factor(:,:),amn_local(:,:),mmn_dummy(:,:,:)
      real(8),allocatable::atom_fractional(:,:)
      real(8),allocatable::gcart(:,:)
      real(8)::a1(3),a2(3),a3(3),lattice(3,3),local_lattice(3,3),lattice_inverse(3,3),x,y,z
      integer::bundle_index,channel,ix,iy,iz,gx,gy,gz,start(3),k,npoint,io,jo,nbasis,natom,atom,last_atom
      integer::nneighbor

      ok=.false.;message='';bundle_index=0;nbasis=entry%n_basis
      do k=1,size(bundles)
        if(bundles(k)%environment==entry%fragment)then;bundle_index=k;exit;end if
      end do
      if(bundle_index==0.or..not.allocated(entry%basis).or.nbasis<=0.or. &
          size(entry%basis,2)/=nbasis.or.size(selected_channels)<=0)then
        message='SAWF representative local seed payload dimensions are inconsistent';return
      end if
      allocate(h_basis(nbasis,nbasis));h_basis=0d0
      do jo=1,nbasis;do io=1,nbasis
        h_basis(io,jo)=0.5d0*(mat_H_local(io,jo,1)+mat_H_local(jo,io,1))
      end do;end do
      call solve_sawf_local_generalized_eigensystem(entry%basis,hvol,h_basis, &
        1d-10,states,energy_hartree,ok,message)
      if(.not.ok)return
      if(present(local_states_out))then
        allocate(local_states_out(size(states,1),size(states,2)));local_states_out=states
      end if
      if(size(selected_channels)>size(states,2))then
        message='SAWF local complete projection shell exceeds generalized eigensystem rank';ok=.false.;return
      end if
      npoint=size(states,1)
      nneighbor=1
      if(present(neighbor_gvec))then
        if(size(neighbor_gvec,1)/=3.or.size(neighbor_gvec,2)<=0)then
          message='SAWF NNKP neighbor payload dimensions are invalid';ok=.false.;return
        end if
        nneighbor=size(neighbor_gvec,2)
      end if
      allocate(projection(npoint,size(selected_channels)),phase_factor(npoint,nneighbor), &
        amn_local(size(states,2),size(selected_channels)),mmn_dummy(size(states,2),size(states,2),nneighbor), &
        energy_ev(size(energy_hartree)))
      call get_lattice_vectors(a1,a2,a3);lattice(:,1)=a1;lattice(:,2)=a2;lattice(:,3)=a3
      lattice_inverse=lattice;call matrix_inverse(lattice_inverse)
      local_lattice=lattice
      do channel=1,3;local_lattice(:,channel)=lattice(:,channel)* &
        real(entry%shape(channel),8)/real(dc%lg_tot%num(channel),8);end do
      allocate(gcart(3,nneighbor));gcart=0d0
      if(present(neighbor_gvec))then
        do channel=1,nneighbor
          call reciprocal_vector_from_index(neighbor_gvec(:,channel),local_lattice(:,1), &
            local_lattice(:,2),local_lattice(:,3),gcart(1,channel),gcart(2,channel),gcart(3,channel))
        end do
      end if
      start=dc%ixyz_frag(:,entry%fragment);k=0
      do iz=1,entry%shape(3);do iy=1,entry%shape(2);do ix=1,entry%shape(1)
        gx=1+modulo(start(1)+ix-1,dc%lg_tot%num(1));gy=1+modulo(start(2)+iy-1,dc%lg_tot%num(2))
        gz=1+modulo(start(3)+iz-1,dc%lg_tot%num(3));k=k+1
        x=dc%lg_tot%coordinate(gx,1);y=dc%lg_tot%coordinate(gy,2);z=dc%lg_tot%coordinate(gz,3)
        do channel=1,size(selected_channels)
          projection(k,channel)=cmplx(sawf_local_projection_value([x,y,z], &
            projection_channels(selected_channels(channel)),lattice,lattice_inverse, &
            wannier_projection_width),0d0,kind=8)
        end do
        do channel=1,nneighbor
          phase_factor(k,channel)=exp(cmplx(0d0,-dot_product(gcart(:,channel),[x,y,z]),kind=8))
        end do
      end do;end do;end do
      call build_sawf_local_seed_matrices(states,projection,phase_factor,hvol,amn_local,mmn_dummy,ok,message)
      if(.not.ok)return
      if(present(local_energy_out))then
        allocate(local_energy_out(size(energy_hartree)));local_energy_out=energy_hartree
      end if
      if(present(local_amn_out))then
        allocate(local_amn_out(size(amn_local,1),size(amn_local,2)));local_amn_out=amn_local
      end if
      energy_ev=energy_hartree*hartree_to_ev
      natom=1;do channel=2,size(selected_channels)
        if(projection_channels(selected_channels(channel))%atom/= &
            projection_channels(selected_channels(channel-1))%atom)natom=natom+1
      end do
      allocate(atom_fractional(3,natom));natom=0;last_atom=0
      do channel=1,size(selected_channels)
        atom=projection_channels(selected_channels(channel))%atom
        if(atom==last_atom)cycle
        natom=natom+1;last_atom=atom
        atom_fractional(:,natom)=modulo((matmul(lattice_inverse,dc%system_tot%rion(:,atom))* &
          real(dc%lg_tot%num,8)-real(start,8))/real(entry%shape,8),1d0)
      end do
      if(present(neighbor_gvec))then
        call write_sawf_local_eig_amn_mmn(bundles(bundle_index)%directory, &
          bundles(bundle_index)%seedname,energy_ev,amn_local,mmn_dummy,neighbor_gvec,ok,message)
      else
        call write_sawf_local_preprocess_win(trim(bundles(bundle_index)%directory)//'/'// &
          trim(bundles(bundle_index)%seedname)//'.win',size(states,2),size(selected_channels),wannier_num_iter, &
          local_lattice,atom_fractional,ok,message)
        if(.not.ok)return
        call write_sawf_local_eig_amn(bundles(bundle_index)%directory,bundles(bundle_index)%seedname, &
          energy_ev,amn_local,ok,message)
      end if
    end subroutine write_sawf_representative_local_seed

    real(8) function sawf_local_projection_value(position,channel,lattice,lattice_inverse,sigma) result(value)
      real(8),intent(in)::position(3),lattice(3,3),lattice_inverse(3,3),sigma
      type(t_sawf_projection_channel),intent(in)::channel
      real(8)::fractional_delta(3),cartesian_delta(3),scaled(3),r2
      logical::image_ok
      fractional_delta=matmul(lattice_inverse,position-dc%system_tot%rion(:,channel%atom))
      call sawf_closest_periodic_cartesian(lattice,fractional_delta,cartesian_delta,image_ok)
      if(.not.image_ok.or.sigma<=0d0)then;value=0d0;return;end if
      r2=sum(cartesian_delta**2)
      if(r2>(8d0*sigma)**2)then;value=0d0;return;end if
      scaled=cartesian_delta/sigma
      value=sawf_real_harmonic_value(channel%l,channel%m,scaled)*exp(-0.5d0*r2/(sigma*sigma))
    end function sawf_local_projection_value

    subroutine reduce_sawf_band_matrix(local_matrix,sum_matrix,nband_wann)
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      integer, intent(in) :: nband_wann
      complex(8), intent(in) :: local_matrix(nband_wann,nband_wann)
      complex(8), intent(out) :: sum_matrix(nband_wann,nband_wann)
      integer, parameter :: mpi_chunk_limit=1000000
      integer :: first_column,last_column,column_chunk,count

      if(nband_wann <= 0 .or. int(nband_wann,int64) > int(huge(0),int64)) &
        call lcfo_sawf_fatal('SAWF D_band MPI row count exceeds default integer')
      column_chunk=max(1,mpi_chunk_limit/nband_wann)
      do first_column=1,nband_wann,column_chunk
        last_column=min(nband_wann,first_column+column_chunk-1)
        if(int(nband_wann,int64)*int(last_column-first_column+1,int64) > int(huge(0),int64)) &
          call lcfo_sawf_fatal('SAWF D_band MPI chunk count overflows default integer')
        count=nband_wann*(last_column-first_column+1)
        call comm_summation(local_matrix(:,first_column:last_column), &
          sum_matrix(:,first_column:last_column),count,dc%icomm_tot)
      end do
    end subroutine reduce_sawf_band_matrix

    subroutine validate_sawf_projection_configuration()
      use salmon_global, only: wannier_site_symmetry, wannier_projection, wannier_num_wann, &
        wannier_dis_froz_max
      implicit none
      integer :: channel_count, projection_lmax
      logical :: ok
      character(256) :: message

      if(trim(wannier_site_symmetry) == 'off') return
      if(wannier_dis_froz_max > 0.0d0) then
        call lcfo_sawf_fatal('SAWF does not support a frozen window; set wannier_dis_froz_max=0')
      end if
      if(.not. is_pseudo_channel_projection(trim(wannier_projection))) then
        call lcfo_sawf_fatal("SAWF requires wannier_projection='pseudo_channels'")
      end if
      call sawf_projection_shell_lmax(dc%system_tot%nion, wannier_num_wann, projection_lmax, ok, message)
      if(.not. ok) call lcfo_sawf_fatal(message)
      call sawf_spd_projection_count(dc%system_tot%nion, channel_count, ok, message, projection_lmax)
      if(.not. ok .or. wannier_num_wann /= channel_count) call lcfo_sawf_fatal(message)
      if(dc%id_tot == 0) write(*,'(1x,a,i0)') &
        '[DC-LCFO-SAWF] SAWF pseudo_channels ordering: complete atom-major shell lmax=', projection_lmax
    end subroutine validate_sawf_projection_configuration

    subroutine lcfo_sawf_fatal(message)
      use parallelization, only: end_parallel
      implicit none
      character(*), intent(in) :: message

      if(dc%id_tot == 0) write(*,'(1x,a,a)') '[FATAL] ', trim(message)
      call end_parallel
      stop 1
    end subroutine lcfo_sawf_fatal

    subroutine write_pseudo_channel_projection_block(iunit, ok, message)
      use salmon_global, only: wannier_site_symmetry, wannier_num_wann
      use salmon_math, only: matrix_inverse
      implicit none
      integer, intent(in) :: iunit
      logical, intent(out) :: ok
      character(*), intent(out) :: message
      integer :: allocation_status, ia, projection_lmax
      character(256) :: allocation_message
      real(8) :: lattice(3,3), lattice_inverse(3,3)
      real(8) :: a1(3), a2(3), a3(3)
      real(8), allocatable :: fractional_positions(:,:)

      ok = .false.
      message = ''
      projection_lmax = 2
      if(trim(wannier_site_symmetry) /= 'off') then
        call sawf_projection_shell_lmax(dc%system_tot%nion, wannier_num_wann, projection_lmax, ok, message)
        if(.not. ok) return
      end if
      if(trim(wannier_site_symmetry) == 'off') then
        allocate(fractional_positions(3,0), stat=allocation_status, errmsg=allocation_message)
        if(allocation_status /= 0) then
          message = 'SAWF off projection allocation failed: '//trim(allocation_message)
          return
        end if
        call write_sawf_projection_block(iunit, wannier_site_symmetry, &
          fractional_positions, ok, message, projection_lmax)
      else
        allocate(fractional_positions(3,dc%system_tot%nion), &
          stat=allocation_status, errmsg=allocation_message)
        if(allocation_status /= 0) then
          message = 'SAWF atom projection allocation failed: '//trim(allocation_message)
          return
        end if
        call get_lattice_vectors(a1, a2, a3)
        lattice(:,1) = a1
        lattice(:,2) = a2
        lattice(:,3) = a3
        lattice_inverse = lattice
        call matrix_inverse(lattice_inverse)
        do ia=1,dc%system_tot%nion
          fractional_positions(:,ia) = matmul(lattice_inverse, dc%system_tot%rion(:,ia))
          fractional_positions(:,ia) = modulo(fractional_positions(:,ia), 1.0d0)
        end do
        call write_sawf_projection_block(iunit, wannier_site_symmetry, &
          fractional_positions, ok, message, projection_lmax)
      end if
      if(allocated(fractional_positions)) deallocate(fractional_positions)
    end subroutine write_pseudo_channel_projection_block

    integer function pseudo_channel_win_projection_lmax(target_wann) result(lmax_write)
      implicit none
      integer, intent(in) :: target_wann
      integer :: ltry

      lmax_write = -1
      do ltry=0,2
        if(count_pseudo_channel_ao_shells_upto(ltry) == target_wann) then
          lmax_write = ltry
          return
        end if
      end do
    end function pseudo_channel_win_projection_lmax

    integer function count_pseudo_channel_ao_shells_upto(lmax_limit) result(nproj)
      use salmon_global, only: izatom
      implicit none
      integer, intent(in) :: lmax_limit
      integer :: ia, iz, lmax_ao

      nproj = 0
      do ia=1,dc%system_tot%nion
        iz = izatom(dc%system_tot%kion(ia))
        lmax_ao = min(pseudo_channel_ao_lmax_for_species(iz), lmax_limit)
        nproj = nproj + count_real_ao_for_lmax(lmax_ao)
      end do
    end function count_pseudo_channel_ao_shells_upto

    subroutine diagnose_wannier_coef_rank(nband_wann)
      use communication, only: comm_summation
      use eigen_subdiag_sub, only: eigen_dsyev
      implicit none
      integer, intent(in) :: nband_wann
      integer :: i, near_null
      real(8) :: diag_min, diag_max, offdiag_max, gram_min, gram_max
      real(8), allocatable :: gram_local(:,:), gram_sum(:,:), eigvec(:,:), eval(:)

      if(nband_wann <= 0) return
      allocate(gram_local(nband_wann,nband_wann), gram_sum(nband_wann,nband_wann))
      allocate(eigvec(nband_wann,nband_wann), eval(nband_wann))
      gram_local = 0d0
      if(dc%id_frag == 0) then
        gram_local(1:nband_wann,1:nband_wann) = &
          matmul(transpose(coef_wf(1:dc%nstate_frag,1:nband_wann,1)), &
          coef_wf(1:dc%nstate_frag,1:nband_wann,1))
      end if
      call comm_summation(gram_local, gram_sum, nband_wann*nband_wann, dc%icomm_tot)

      diag_min = huge(1d0)
      diag_max = -huge(1d0)
      offdiag_max = 0d0
      do i=1,nband_wann
        diag_min = min(diag_min, gram_sum(i,i))
        diag_max = max(diag_max, gram_sum(i,i))
        gram_sum(i,i) = gram_sum(i,i) - 1d0
      end do
      offdiag_max = maxval(abs(gram_sum(1:nband_wann,1:nband_wann)))
      do i=1,nband_wann
        gram_sum(i,i) = gram_sum(i,i) + 1d0
      end do

      call eigen_dsyev(gram_sum, eval, eigvec)
      gram_min = minval(eval(1:nband_wann))
      gram_max = maxval(eval(1:nband_wann))
      near_null = count(eval(1:nband_wann) < 1d-8 * max(gram_max, 1d-300))
      if(dc%id_tot == 0) then
        write(*,'(1x,a,i0,5(a,es12.5),a,i0)') &
          "[DC-LCFO-WANNIER-COEF] bands=", nband_wann, &
          " diag_min=", diag_min, " diag_max=", diag_max, &
          " offdiag_max=", offdiag_max, " gram_min=", gram_min, &
          " gram_max=", gram_max, " near_null=", near_null
      end if
      deallocate(gram_local, gram_sum, eigvec, eval)
    end subroutine diagnose_wannier_coef_rank

    subroutine diagnose_wannier_amn_conditioning(nband_wann, num_wann)
      use communication, only: comm_bcast
      use eigen_subdiag_sub, only: eigen_zheev
      use filesystem, only: get_filehandle
      use salmon_global, only: sysname, wannier_amn_svd_tol, wannier_amn_reject_tol
      implicit none
      integer, intent(in) :: nband_wann, num_wann
      integer :: iunit, io, iband, iwann, ikpt, irec
      integer :: num_bands_file, num_kpts_file, num_wann_file
      integer :: near_null, reject_flag
      real(8) :: re_val, im_val, smax, smin, min_sv_rel, sv_tol, reject_tol
      real(8), allocatable :: eval(:)
      complex(8), allocatable :: amat(:,:), gram(:,:), eigvec(:,:)
      character(256) :: amnfile, header

      reject_flag = 0
      if(dc%id_tot == 0) then
        amnfile = trim(dc%base_directory)//trim(sysname)//".amn"
        iunit = get_filehandle()
        open(iunit,file=amnfile,status='old',action='read',iostat=io)
        if(io /= 0) then
          write(*,'(1x,2a)') "[DC-LCFO-WANNIER-AMN] failed to open AMN file: ", trim(amnfile)
          reject_flag = 1
        else
          read(iunit,'(a)',iostat=io) header
          if(io == 0) read(iunit,*,iostat=io) num_bands_file, num_kpts_file, num_wann_file
          if(io /= 0 .or. num_bands_file /= nband_wann .or. num_wann_file /= num_wann .or. &
              num_kpts_file /= 1) then
            write(*,'(1x,a,6(a,i0))') "[DC-LCFO-WANNIER-AMN] dimension mismatch:", &
              " bands=", num_bands_file, " expected_bands=", nband_wann, &
              " wann=", num_wann_file, " expected_wann=", num_wann, &
              " kpts=", num_kpts_file, " expected_kpts=", 1
            reject_flag = 1
          else
            allocate(amat(nband_wann,num_wann), gram(num_wann,num_wann), eigvec(num_wann,num_wann), eval(num_wann))
            amat = (0.0d0,0.0d0)
            do irec=1,nband_wann*num_wann
              read(iunit,*,iostat=io) iband, iwann, ikpt, re_val, im_val
              if(io /= 0) exit
              if(iband < 1 .or. iband > nband_wann) cycle
              if(iwann < 1 .or. iwann > num_wann) cycle
              if(ikpt /= 1) cycle
              amat(iband,iwann) = cmplx(re_val, im_val, kind=8)
            end do
            if(io /= 0) then
              write(*,'(1x,a)') "[DC-LCFO-WANNIER-AMN] failed to read complete AMN matrix."
              reject_flag = 1
            else
              gram = matmul(conjg(transpose(amat)), amat)
              call eigen_zheev(gram, eval, eigvec)
              eval(1:num_wann) = max(eval(1:num_wann), 0.0d0)
              smax = sqrt(maxval(eval(1:num_wann)))
              smin = sqrt(minval(eval(1:num_wann)))
              if(smax > 0.0d0) then
                min_sv_rel = smin / smax
              else
                min_sv_rel = 0.0d0
              end if
              sv_tol = max(wannier_amn_svd_tol, 0.0d0)
              reject_tol = max(wannier_amn_reject_tol, 0.0d0)
              if(smax > 0.0d0 .and. sv_tol > 0.0d0) then
                near_null = count(sqrt(eval(1:num_wann)) < sv_tol * smax)
              else
                near_null = 0
              end if
              write(*,'(1x,a,i0,a,i0,4(a,es12.5),a,i0)') &
                "[DC-LCFO-WANNIER-AMN] bands=", nband_wann, " wann=", num_wann, &
                " smax=", smax, " smin=", smin, " min_sv_rel=", min_sv_rel, &
                " sv_tol=", sv_tol, " near_null=", near_null
              if(reject_tol > 0.0d0 .and. min_sv_rel < reject_tol) then
                write(*,'(1x,a,2(a,es12.5))') &
                  "[DC-LCFO-WANNIER-AMN] reject: min_sv_rel below wannier_amn_reject_tol", &
                  " min_sv_rel=", min_sv_rel, " wannier_amn_reject_tol=", reject_tol
                reject_flag = 1
              end if
            end if
            deallocate(amat, gram, eigvec, eval)
          end if
          close(iunit)
        end if
      end if

      call comm_bcast(reject_flag, dc%icomm_tot, 0)
      if(reject_flag /= 0) &
        stop "DC-LCFO Wannier export: AMN projection matrix is rank deficient; adjust projections or wannier_amn_reject_tol."
    end subroutine diagnose_wannier_amn_conditioning

    subroutine log_wannier_cluster_partition()
      use salmon_global, only: num_fragment, num_wannier_cluster, wannier_cluster_size
      implicit none
      integer :: num_cluster(3), icx, icy, icz, cid, frag_lo(3), frag_hi(3)
      integer :: ncluster_tot

      if(dc%id_tot /= 0) return
      num_cluster(1:3) = num_fragment(1:3) / wannier_cluster_size(1:3)
      ncluster_tot = num_cluster(1) * num_cluster(2) * num_cluster(3)
      write(*,'(1x,a,3(i0,1x),a,3(i0,1x),a,3(i0,1x),a,i0)') &
        "[DC-LCFO-WANNIER-CLUSTER] fragment_grid=", num_fragment(1:3), &
        " input_num_wannier_cluster=", num_wannier_cluster(1:3), &
        " cluster_size=", wannier_cluster_size(1:3), " ncluster=", ncluster_tot
      do icx=1,num_cluster(1)
      do icy=1,num_cluster(2)
      do icz=1,num_cluster(3)
        cid = ((icx - 1) * num_cluster(2) + (icy - 1)) * num_cluster(3) + icz
        if(ncluster_tot > 64 .and. cid > 8) cycle
        frag_lo(1) = (icx - 1) * wannier_cluster_size(1) + 1
        frag_hi(1) = icx * wannier_cluster_size(1)
        frag_lo(2) = (icy - 1) * wannier_cluster_size(2) + 1
        frag_hi(2) = icy * wannier_cluster_size(2)
        frag_lo(3) = (icz - 1) * wannier_cluster_size(3) + 1
        frag_hi(3) = icz * wannier_cluster_size(3)
        write(*,'(1x,a,i0,a,3(i0,a,i0,1x))') &
          "[DC-LCFO-WANNIER-CLUSTER] cluster=", cid, " frag_ranges=", &
          frag_lo(1), ":", frag_hi(1), frag_lo(2), ":", frag_hi(2), frag_lo(3), ":", frag_hi(3)
      end do
      end do
      end do
      if(ncluster_tot > 64) write(*,'(1x,a)') &
        "[DC-LCFO-WANNIER-CLUSTER] cluster list truncated after the first 8 clusters."
    end subroutine log_wannier_cluster_partition

    subroutine write_wannier90_global_basis_file(nband_wann)
      use communication, only: comm_sync_all
      use filesystem, only: get_filehandle
      use inputoutput, only: au_length_aa
      use salmon_global, only: sysname, wannier_num_wann, wannier_projection
      implicit none
      integer, intent(in) :: nband_wann
      integer :: iunit, iw
      integer :: num_bands_chk, num_wann_chk
      integer :: position_available
      integer, allocatable :: owner_frag(:)
      real(8), allocatable :: center_aa(:,:), center_bohr(:,:), spread_aa2(:)
      complex(8), allocatable :: v_matrix(:,:), aa_global(:,:,:)
      logical :: ok_position
      character(256) :: filename

      if(dc%id_tot == 0) then
        call read_wannier90_checkpoint_transform(num_bands_chk, num_wann_chk, v_matrix, center_aa, spread_aa2)
        if(num_bands_chk /= nband_wann .or. num_wann_chk /= wannier_num_wann) then
          write(*,'(1x,a,4(a,i0))') "[DC-LCFO-W90-GLOBAL] checkpoint dimension mismatch:", &
            " chk_bands=", num_bands_chk, " expected_bands=", nband_wann, &
            " chk_wann=", num_wann_chk, " expected_wann=", wannier_num_wann
          stop "DC-LCFO Wannier export: checkpoint dimension mismatch"
        end if
        allocate(owner_frag(num_wann_chk), center_bohr(3,num_wann_chk))
        center_bohr(1:3,1:num_wann_chk) = center_aa(1:3,1:num_wann_chk) / au_length_aa
        call wrap_wannier_centers_to_total_cell_import(dc, center_bohr, num_wann_chk)
        do iw=1,num_wann_chk
          owner_frag(iw) = find_owner_fragment_from_center(center_bohr(1:3,iw))
        end do
        call rebalance_wannier_owner_fragments(center_bohr, owner_frag, num_wann_chk)
        write(*,'(1x,a)') "[DC-LCFO-W90-GLOBAL] owner source=wannier_centers"
        if(is_bond_center_projection(trim(wannier_projection))) then
          write(*,'(1x,a)') "[DC-LCFO-W90-GLOBAL] projection source=bond_centers"
        end if
        call read_wannier90_global_rmn_gamma_block(num_wann_chk, aa_global, ok_position)
        if(.not. ok_position) then
          allocate(aa_global(3,num_wann_chk,num_wann_chk))
          aa_global = (0d0,0d0)
        else
          call set_wannier_centers_from_position_diagonal(num_wann_chk, center_bohr, aa_global)
          do iw=1,num_wann_chk
            owner_frag(iw) = find_owner_fragment_from_center(center_bohr(1:3,iw))
          end do
          call rebalance_wannier_owner_fragments(center_bohr, owner_frag, num_wann_chk)
          write(*,'(1x,a)') "[DC-LCFO-W90-GLOBAL] center source=AA_R diagonal"
        end if
        position_available = merge(1, 0, ok_position)

        filename = trim(dc%base_directory)//binfile_w90g
        call write_wannier90_global_basis_stream_import(filename, num_bands_chk, num_wann_chk, dc%n_frag, &
          owner_frag, center_bohr, spread_aa2, v_matrix, position_available, aa_global)
        write(*,'(1x,a,i0,a,a)') "[DC-LCFO-W90-GLOBAL] wrote ", num_wann_chk, &
          " Wannier functions to ", trim(filename)
        filename = trim(dc%base_directory)//binfile_w90g_persistent
        call write_wannier90_global_basis_stream_import(filename, num_bands_chk, num_wann_chk, dc%n_frag, &
          owner_frag, center_bohr, spread_aa2, v_matrix, position_available, aa_global)
        if(ok_position) then
          write(*,'(1x,a)') "[DC-LCFO-W90-GLOBAL] stored Wannier90 AA_R position matrix."
        else
          write(*,'(1x,a)') "[DC-LCFO-W90-GLOBAL] AA_R position matrix unavailable; stored zero matrix."
        end if
        deallocate(owner_frag, center_bohr, center_aa, spread_aa2, v_matrix, aa_global)
      end if
      call comm_sync_all(dc%icomm_tot)
    end subroutine write_wannier90_global_basis_file

    subroutine write_wannier90_flux_eigen_seed_file(nband_wann)
      use communication, only: comm_sync_all
      use filesystem, only: get_filehandle
      use salmon_global, only: wannier_num_wann
      implicit none
      integer, intent(in) :: nband_wann
      integer :: iunit, iw, istate
      integer :: num_bands_chk, num_wann_chk, nstate_seed
      real(8), allocatable :: center_aa(:,:), spread_aa2(:), eval_seed(:,:)
      complex(8), allocatable :: v_matrix(:,:), seed_wannier_to_eigen(:,:)
      character(256) :: filename

      if(dc%id_tot == 0) then
        call read_wannier90_checkpoint_transform(num_bands_chk, num_wann_chk, v_matrix, center_aa, spread_aa2)
        if(num_bands_chk /= nband_wann .or. num_wann_chk /= wannier_num_wann) then
          write(*,'(1x,a,4(a,i0))') "[DC-LCFO-W90-SEED] checkpoint dimension mismatch:", &
            " chk_bands=", num_bands_chk, " expected_bands=", nband_wann, &
            " chk_wann=", num_wann_chk, " expected_wann=", wannier_num_wann
          stop "DC-LCFO Wannier flux seed: checkpoint dimension mismatch"
        end if
        call build_wannier_flux_eigen_seed_from_transform(num_bands_chk, num_wann_chk, dc%nstate_tot, &
          nspin, v_matrix, esp_tot, seed_wannier_to_eigen, eval_seed, nstate_seed)

        filename = trim(dc%base_directory)//binfile_w90seed
        call write_wannier_flux_eigen_seed_stream_import(filename, num_bands_chk, num_wann_chk, nstate_seed, &
          nspin, dc%n_frag, eval_seed, seed_wannier_to_eigen)
        write(*,'(1x,a,i0,a,i0,a,a)') "[DC-LCFO-W90-SEED] wrote Flux-LCFO eigen seed in Wannier basis: states=", &
          nstate_seed, " wann=", num_wann_chk, " file=", trim(filename)
        filename = trim(dc%base_directory)//binfile_w90seed_persistent
        call write_wannier_flux_eigen_seed_stream_import(filename, num_bands_chk, num_wann_chk, nstate_seed, &
          nspin, dc%n_frag, eval_seed, seed_wannier_to_eigen)
        deallocate(seed_wannier_to_eigen, eval_seed, center_aa, spread_aa2, v_matrix)
      end if
      call comm_sync_all(dc%icomm_tot)
    end subroutine write_wannier90_flux_eigen_seed_file

    subroutine write_wannier90_global_bpw_seed()
      use eigen_subdiag_sub, only: eigen_dsyev
      use communication, only: comm_sync_all
      use filesystem, only: get_filehandle
      use salmon_global, only: base_directory, lambda_cut
      implicit none
      integer :: iunit, nbasis, num_bands_file, num_wann_file, nfrag_file
      integer :: iw, jw, iband, axis, nkeep, io, jo, i_full, j_full, nsel, imode, i_halo
      integer :: nxyz_domain(3), nxyz_buffer_seed(3), nxyz_box(3)
      integer :: magic_file, version_file, io_stat
      integer, allocatable :: owner_frag(:), proj_atom(:), proj_hybrid(:), keep_index(:)
      real(8), allocatable :: center_bohr(:,:), spread_aa2(:), lambda_w(:), spread_est(:), tail_est(:)
      real(8), allocatable :: wcoef(:,:), r_wann(:,:,:), wcenter(:,:), h_wann(:,:), v_wann(:,:,:), aa_wann(:,:,:)
      real(8), allocatable :: h_seed(:,:), v_seed(:,:,:), tmp(:,:)
      real(8), allocatable :: p_frag(:,:), p_eval(:), p_evec(:,:), tmat(:,:)
      complex(8), allocatable :: v_matrix(:,:), aa_global(:,:,:)
      complex(8), allocatable :: psi_wann_c(:,:)
      real(8) :: h_imag_max, w_imag_max, bpw_lambda_cut, x
      logical :: ok_position
      character(256) :: filename

      if(dc%id_frag /= 0) return
      if(nspin /= 1) stop "DC-LCFO global Wannier BPW seed: spin-polarized mode is not implemented."

      call read_wannier90_global_basis_file_full(num_bands_file, num_wann_file, nfrag_file, &
        owner_frag, center_bohr, spread_aa2, v_matrix)
      if(num_bands_file <= 0 .or. num_wann_file <= 0) return
      if(num_bands_file /= dc%nstate_tot) then
        if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0)') &
          "[DC-LCFO-W90-BPW] skip: global bands=", num_bands_file, " dc_nstate_tot=", dc%nstate_tot
        deallocate(owner_frag, center_bohr, spread_aa2, v_matrix)
        return
      end if
      if(nfrag_file /= dc%n_frag) then
        if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0)') &
          "[DC-LCFO-W90-BPW] skip: global nfrag=", nfrag_file, " dc_n_frag=", dc%n_frag
        deallocate(owner_frag, center_bohr, spread_aa2, v_matrix)
        return
      end if

      nkeep = num_wann_file / max(1, dc%n_frag)
      if(mod(num_wann_file, max(1, dc%n_frag)) /= 0) then
        if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0)') &
          "[FATAL] global Wannier count is not divisible by fragments: wann=", &
          num_wann_file, " nfrag=", dc%n_frag
        stop "DC-LCFO global Wannier BPW seed: nonuniform target count"
      end if

      call read_wannier90_global_rmn_gamma_block(num_wann_file, aa_global, ok_position)
      if(.not. ok_position) allocate(aa_global(3,num_wann_file,num_wann_file))
      if(.not. ok_position) aa_global = (0d0,0d0)

      nbasis = n_basis(dc%i_frag,1)
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      nxyz_buffer_seed(1:3) = dc%nxyz_buffer(1:3)
      nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer_seed(1:3)

      allocate(proj_atom(nkeep), proj_hybrid(nkeep), keep_index(nkeep))
      allocate(lambda_w(nkeep), spread_est(nkeep), tail_est(nkeep))
      allocate(wcoef(nbasis,nkeep), psi_wann_c(nbasis,num_wann_file))
      allocate(r_wann(3,nkeep,nkeep), wcenter(3,nkeep), h_wann(nkeep,nkeep), &
        v_wann(3,nkeep,nkeep), aa_wann(3,nkeep,nkeep))
      allocate(h_seed(nbasis,nbasis), v_seed(3,nbasis,nbasis), tmp(nbasis,nkeep))
      allocate(p_frag(nbasis,nbasis), p_eval(nbasis), p_evec(nbasis,nbasis), &
        tmat(nkeep,num_wann_file))

      proj_atom = 0
      proj_hybrid = 0
      keep_index = 0

      psi_wann_c = (0d0,0d0)
      do iw=1,num_wann_file
        do iband=1,num_bands_file
          psi_wann_c(1:nbasis,iw) = psi_wann_c(1:nbasis,iw) + &
            cmplx(coef_wf(1:nbasis,iband,1), 0d0, kind=8) * &
            v_matrix(iband,iw)
        end do
      end do
      w_imag_max = maxval(abs(aimag(psi_wann_c(1:nbasis,1:num_wann_file))))

      p_frag = 0d0
      do iw=1,num_wann_file
        do io=1,nbasis
          do jo=1,nbasis
            p_frag(io,jo) = p_frag(io,jo) + &
              real(psi_wann_c(io,iw),kind=8) * real(psi_wann_c(jo,iw),kind=8) + &
              aimag(psi_wann_c(io,iw)) * aimag(psi_wann_c(jo,iw))
          end do
        end do
      end do
      call eigen_dsyev(p_frag, p_eval, p_evec)
      bpw_lambda_cut = max(1d-14, min(lambda_cut, maxval(p_eval(1:nbasis)) * 1d-12))
      wcoef = 0d0
      nsel = 0
      do imode=nbasis,1,-1
        if(p_eval(imode) <= bpw_lambda_cut) cycle
        nsel = nsel + 1
        if(nsel > nkeep) exit
        do io=1,nbasis
          wcoef(io,nsel) = p_evec(io,imode)
        end do
        lambda_w(nsel) = p_eval(imode)
        keep_index(nsel) = imode
      end do
      if(nsel < nkeep) then
        write(*,'(1x,a,i0,a,i0,a,i0,3(a,1pe12.4))') &
          "[WARN] local Wannier projection rank is small; completing with local LCFO vectors: fragment=", dc%i_frag, &
          " keep=", nsel, " expected=", nkeep, " bpw_cut=", bpw_lambda_cut, &
          " lambda_cut=", lambda_cut, " max_eval=", maxval(p_eval(1:nbasis))
      end if
      do imode=1,nbasis
        if(nsel >= nkeep) exit
        p_evec(1:nbasis,1) = 0d0
        p_evec(imode,1) = 1d0
        do iw=1,nsel
          p_evec(1:nbasis,1) = p_evec(1:nbasis,1) - &
            sum(p_evec(1:nbasis,1) * wcoef(1:nbasis,iw)) * wcoef(1:nbasis,iw)
        end do
        if(sum(p_evec(1:nbasis,1) * p_evec(1:nbasis,1)) <= 1d-20) cycle
        nsel = nsel + 1
        wcoef(1:nbasis,nsel) = p_evec(1:nbasis,1) / &
          sqrt(sum(p_evec(1:nbasis,1) * p_evec(1:nbasis,1)))
        lambda_w(nsel) = 0d0
        keep_index(nsel) = -imode
      end do
      if(nsel < nkeep) then
        write(*,'(1x,a,i0,a,i0,a,i0)') &
          "[FATAL] failed to complete BPW local basis: fragment=", dc%i_frag, &
          " keep=", nsel, " expected=", nkeep
        stop "DC-LCFO global Wannier BPW seed: insufficient local completion rank"
      end if
      do iw=1,nkeep
        if(keep_index(iw) > 0 .and. lambda_w(iw) <= bpw_lambda_cut) then
          if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0,a,1pe12.4)') &
            "[FATAL] owned Wannier projection is singular: fragment=", dc%i_frag, &
            " mode=", iw, " lambda=", lambda_w(iw)
          stop "DC-LCFO global Wannier BPW seed: singular owned Wannier projection"
        end if
      end do
      tail_est(1:nkeep) = 0d0
      tmat(1:nkeep,1:num_wann_file) = matmul(transpose(wcoef(1:nbasis,1:nkeep)), &
        real(psi_wann_c(1:nbasis,1:num_wann_file), kind=8))
      spread_est(1:nkeep) = 0d0
      do iw=1,nkeep
        do j_full=1,num_wann_file
          spread_est(iw) = spread_est(iw) + tmat(iw,j_full) * tmat(iw,j_full) * &
            spread_aa2(j_full)
        end do
      end do

      h_imag_max = 0d0
      h_seed = 0d0
      v_seed = 0d0
      do io=1,nbasis
        do jo=1,nbasis
          h_seed(io,jo) = 0.5d0 * (mat_H_local(io,jo,1) + mat_H_local(jo,io,1))
          v_seed(1:3,io,jo) = mat_V_local(1:3,io,jo,1)
          do i_halo=1,n_halo
            h_seed(io,jo) = h_seed(io,jo) + 0.25d0 * &
              (halo(i_halo)%mat_H_local(jo,io,1) + halo(i_halo)%mat_H_local(io,jo,1))
            v_seed(1:3,io,jo) = v_seed(1:3,io,jo) + 0.5d0 * &
              (halo(i_halo)%mat_V_local(1:3,jo,io,1) + &
               halo(i_halo)%mat_V_local(1:3,io,jo,1))
          end do
        end do
      end do
      do axis=1,3
        do io=1,nbasis
          v_seed(axis,io,io) = 0d0
          do jo=io+1,nbasis
            x = 0.5d0 * (v_seed(axis,io,jo) - v_seed(axis,jo,io))
            v_seed(axis,io,jo) = x
            v_seed(axis,jo,io) = -x
          end do
        end do
      end do

      tmp(1:nbasis,1:nkeep) = matmul(h_seed(1:nbasis,1:nbasis), wcoef(1:nbasis,1:nkeep))
      h_wann(1:nkeep,1:nkeep) = matmul(transpose(wcoef(1:nbasis,1:nkeep)), tmp(1:nbasis,1:nkeep))
      do axis=1,3
        tmp(1:nbasis,1:nkeep) = matmul(v_seed(axis,1:nbasis,1:nbasis), wcoef(1:nbasis,1:nkeep))
        v_wann(axis,1:nkeep,1:nkeep) = matmul(transpose(wcoef(1:nbasis,1:nkeep)), tmp(1:nbasis,1:nkeep))
      end do

      r_wann = 0d0
      aa_wann = 0d0
      do axis=1,3
        do iw=1,nkeep
          wcenter(axis,iw) = 0d0
          do j_full=1,num_wann_file
            wcenter(axis,iw) = wcenter(axis,iw) + tmat(iw,j_full) * tmat(iw,j_full) * &
              center_bohr(axis,j_full)
          end do
          r_wann(axis,iw,iw) = wcenter(axis,iw)
        end do
      end do
      do axis=1,3
        do iw=1,nkeep
          do jw=1,nkeep
            aa_wann(axis,iw,jw) = 0d0
            do i_full=1,num_wann_file
              do j_full=1,num_wann_file
                aa_wann(axis,iw,jw) = aa_wann(axis,iw,jw) + tmat(iw,i_full) * &
                  real(aa_global(axis,i_full,j_full), kind=8) * tmat(jw,j_full)
              end do
            end do
          end do
        end do
      end do

      filename = trim(base_directory)//binfile_bpw
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='replace',iostat=io_stat)
      if(io_stat /= 0) stop "DC-LCFO global Wannier BPW seed: failed to open output file"
      write(iunit) buffer_periodic_wannier_magic, buffer_periodic_wannier_version
      write(iunit) dc%i_frag, nxyz_domain(1:3), nxyz_buffer_seed(1:3), nxyz_box(1:3)
      write(iunit) nspin, nbasis, nkeep, nkeep
      write(iunit) proj_atom(1:nkeep), proj_hybrid(1:nkeep)
      write(iunit) lambda_w(1:nkeep), keep_index(1:nkeep)
      write(iunit) spread_est(1:nkeep), tail_est(1:nkeep)
      write(iunit) wcoef(1:nbasis,1:nkeep)
      write(iunit) r_wann(1:3,1:nkeep,1:nkeep)
      write(iunit) wcenter(1:3,1:nkeep)
      write(iunit) h_wann(1:nkeep,1:nkeep)
      write(iunit) v_wann(1:3,1:nkeep,1:nkeep)
      write(iunit) aa_wann(1:3,1:nkeep,1:nkeep)
      write(iunit) merge(1, 0, ok_position)
      close(iunit)

      call write_wannier90_global_trace_file(nbasis, nkeep, wcoef, &
        nxyz_domain, nxyz_buffer_seed, nxyz_box)

      write(*,'(1x,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "[DC-LCFO-W90-BPW] fragment=", dc%i_frag, " keep=", nkeep, &
        " min_proj_eval=", minval(lambda_w(1:nkeep)), &
        " max_coef_imag=", w_imag_max, " max_h_imag_term=", h_imag_max

      deallocate(owner_frag, center_bohr, spread_aa2, v_matrix, aa_global)
      deallocate(proj_atom, proj_hybrid, keep_index, lambda_w, spread_est, tail_est)
      deallocate(wcoef, psi_wann_c, r_wann, wcenter, h_wann, v_wann, aa_wann)
      deallocate(h_seed, v_seed, tmp)
      deallocate(p_frag, p_eval, p_evec, tmat)
      call comm_sync_all(dc%icomm_tot)
    end subroutine write_wannier90_global_bpw_seed

    subroutine write_wannier90_global_trace_file(nbasis, nkeep, wcoef, &
        nxyz_domain, nxyz_buffer_seed, nxyz_box)
      use filesystem, only: get_filehandle
      use salmon_global, only: base_directory
      implicit none
      integer, intent(in) :: nbasis, nkeep
      integer, intent(in) :: nxyz_domain(3), nxyz_buffer_seed(3), nxyz_box(3)
      real(8), intent(in) :: wcoef(nbasis,nkeep)
      integer :: iunit, ibasis_read, ispin_read, iw
      integer :: magic_file, version_file, nspin_file, nstate_frag_file
      integer :: nxyz_domain_file(3), nxyz_buffer_file(3)
      integer :: axis, side, face, npts, face_pt
      integer :: ix_face, iy_face, iz_face, ibx, iby, ibz, dist
      integer, allocatable :: n_basis_file(:)
      real(8), allocatable :: phi_tmp(:,:,:), phi_wann(:,:,:,:)
      real(8), allocatable :: face_u(:,:), face_dn(:,:)
      real(8) :: grad_axis, area_weight, alpha
      character(256) :: filename
      real(8), parameter :: surface_penalty_factor = 10.0d0

      if(dc%id_frag /= 0) return
      if(nspin /= 1) return

      allocate(phi_tmp(nxyz_box(1),nxyz_box(2),nxyz_box(3)))
      allocate(phi_wann(nxyz_box(1),nxyz_box(2),nxyz_box(3),nkeep))
      allocate(n_basis_file(nspin))
      phi_wann = 0d0

      filename = trim(base_directory)//binfile_bfb
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='old')
      read(iunit) magic_file, version_file
      if(magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) &
        stop "DC-LCFO global Wannier trace export: invalid buffered basis header."
      read(iunit) nxyz_domain_file(1:3), nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
      if(any(nxyz_domain_file(1:3) /= nxyz_domain(1:3)) .or. &
         any(nxyz_buffer_file(1:3) /= nxyz_buffer_seed(1:3)) .or. &
         nspin_file /= nspin .or. nstate_frag_file /= dc%nstate_frag) &
        stop "DC-LCFO global Wannier trace export: buffered basis metadata mismatch."
      read(iunit) n_basis_file(1:nspin)
      if(n_basis_file(1) < nbasis) stop "DC-LCFO global Wannier trace export: basis count mismatch."
      do ispin_read=1,nspin_file
        do ibasis_read=1,nstate_frag_file
          read(iunit) phi_tmp(1:nxyz_box(1),1:nxyz_box(2),1:nxyz_box(3))
          if(ispin_read /= 1 .or. ibasis_read > nbasis) cycle
          do iw=1,nkeep
            if(abs(wcoef(ibasis_read,iw)) <= 0d0) cycle
            phi_wann(:,:,:,iw) = phi_wann(:,:,:,iw) + &
              wcoef(ibasis_read,iw) * phi_tmp(:,:,:)
          end do
        end do
      end do
      close(iunit)

      filename = trim(base_directory)//binfile_bpwt
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='replace')
      write(iunit) buffer_periodic_wannier_trace_magic, buffer_periodic_wannier_trace_version
      write(iunit) dc%i_frag, nxyz_domain(1:3), nxyz_buffer_seed(1:3), nxyz_box(1:3)
      write(iunit) system%hgs(1:3), hvol, nkeep
      do face=1,6
        axis = (face + 1) / 2
        if(mod(face,2) == 1) then
          side = -1
        else
          side = 1
        end if
        npts = face_point_count(axis)
        area_weight = hvol / system%hgs(axis)
        alpha = surface_penalty_factor / system%hgs(axis)
        allocate(face_u(npts,nkeep), face_dn(npts,nkeep))
        face_u = 0d0
        face_dn = 0d0
        face_pt = 0
        do iz_face=1,nxyz_domain(3)
        do iy_face=1,nxyz_domain(2)
        do ix_face=1,nxyz_domain(1)
          if(face_axis_index([ix_face,iy_face,iz_face], axis) /= face_coord(axis, side)) cycle
          face_pt = face_pt + 1
          ibx = ix_face + nxyz_buffer_seed(1)
          iby = iy_face + nxyz_buffer_seed(2)
          ibz = iz_face + nxyz_buffer_seed(3)
          do iw=1,nkeep
            face_u(face_pt,iw) = phi_wann(ibx,iby,ibz,iw)
            grad_axis = 0d0
            do dist=1,size(stencil%coef_nab,1)
              select case(axis)
              case(1)
                grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
                  (phi_wann(ibx+dist,iby,ibz,iw) - phi_wann(ibx-dist,iby,ibz,iw))
              case(2)
                grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
                  (phi_wann(ibx,iby+dist,ibz,iw) - phi_wann(ibx,iby-dist,ibz,iw))
              case(3)
                grad_axis = grad_axis + stencil%coef_nab(dist,axis) * &
                  (phi_wann(ibx,iby,ibz+dist,iw) - phi_wann(ibx,iby,ibz-dist,iw))
              end select
            end do
            face_dn(face_pt,iw) = real(side,8) * grad_axis
          end do
        end do
        end do
        end do
        write(iunit) axis, side, npts, area_weight, alpha
        write(iunit) face_u(1:npts,1:nkeep)
        write(iunit) face_dn(1:npts,1:nkeep)
        deallocate(face_u, face_dn)
      end do
      close(iunit)
      write(*,'(1x,a,i0,a,a)') "[DC-LCFO-W90-BPW-TRACE] fragment=", dc%i_frag, &
        " wrote ", trim(filename)

      deallocate(phi_tmp, phi_wann, n_basis_file)
    end subroutine write_wannier90_global_trace_file

    subroutine read_wannier90_global_basis_file_full(num_bands_file, num_wann_file, nfrag_file, &
        owner_frag, center_bohr, spread_aa2, v_matrix)
      use filesystem, only: get_filehandle
      implicit none
      integer, intent(out) :: num_bands_file, num_wann_file, nfrag_file
      integer, allocatable, intent(out) :: owner_frag(:)
      real(8), allocatable, intent(out) :: center_bohr(:,:), spread_aa2(:)
      complex(8), allocatable, intent(out) :: v_matrix(:,:)
      integer :: iunit_local, io_local, magic_file, version_file
      character(256) :: filename_local
      logical :: exists_local

      num_bands_file = 0
      num_wann_file = 0
      nfrag_file = 0
      filename_local = trim(dc%base_directory)//binfile_w90g
      inquire(file=filename_local, exist=exists_local)
      if(.not. exists_local) then
        write(*,'(1x,2a)') "[DC-LCFO-W90-BPW] global Wannier basis not found: ", trim(filename_local)
        return
      end if
      iunit_local = get_filehandle()
      open(iunit_local,file=filename_local,form='unformatted',access='stream',status='old',iostat=io_local)
      if(io_local /= 0) return
      read(iunit_local,iostat=io_local) magic_file, version_file
      if(io_local /= 0 .or. magic_file /= wannier90_global_magic .or. &
          version_file < 1 .or. version_file > wannier90_global_version) then
        close(iunit_local)
        return
      end if
      read(iunit_local,iostat=io_local) num_bands_file, num_wann_file, nfrag_file
      if(io_local /= 0 .or. num_bands_file <= 0 .or. num_wann_file <= 0) then
        close(iunit_local)
        num_bands_file = 0
        num_wann_file = 0
        return
      end if
      allocate(owner_frag(num_wann_file), center_bohr(3,num_wann_file), spread_aa2(num_wann_file))
      allocate(v_matrix(num_bands_file,num_wann_file))
      read(iunit_local,iostat=io_local) owner_frag(1:num_wann_file)
      if(io_local == 0) read(iunit_local,iostat=io_local) center_bohr(1:3,1:num_wann_file)
      if(io_local == 0) read(iunit_local,iostat=io_local) spread_aa2(1:num_wann_file)
      if(io_local == 0) read(iunit_local,iostat=io_local) v_matrix(1:num_bands_file,1:num_wann_file)
      close(iunit_local)
      if(io_local /= 0) then
        if(allocated(owner_frag)) deallocate(owner_frag)
        if(allocated(center_bohr)) deallocate(center_bohr)
        if(allocated(spread_aa2)) deallocate(spread_aa2)
        if(allocated(v_matrix)) deallocate(v_matrix)
        num_bands_file = 0
        num_wann_file = 0
        nfrag_file = 0
      end if
    end subroutine read_wannier90_global_basis_file_full

    subroutine read_wannier90_checkpoint_transform(num_bands_chk,num_wann_chk,v_matrix,center_aa, &
        spread_aa2,filename_in)
      use filesystem, only: get_filehandle
      use salmon_global, only: sysname
      implicit none
      integer, intent(out) :: num_bands_chk, num_wann_chk
      complex(8), allocatable, intent(out) :: v_matrix(:,:)
      real(8), allocatable, intent(out) :: center_aa(:,:), spread_aa2(:)
      character(*),intent(in),optional::filename_in
      integer :: iunit, io, i, j, k, nkp
      integer :: num_exclude_bands_chk, num_kpts_chk, nntot_chk
      integer :: mp_grid_chk(3)
      integer, allocatable :: exclude_bands_chk(:), ndimwin(:)
      real(8) :: real_lattice_chk(3,3), recip_lattice_chk(3,3), omega_invariant_chk
      real(8), allocatable :: kpt_latt_chk(:,:)
      logical :: have_disentangled_chk
      logical, allocatable :: lwindow(:,:)
      complex(8), allocatable :: u_matrix(:,:,:), u_matrix_opt(:,:,:), m_matrix(:,:,:,:)
      character(1024) :: filename
      character(33) :: header_chk
      character(20) :: checkpoint_chk

      filename = trim(dc%base_directory)//trim(sysname)//".chk"
      if(present(filename_in))filename=trim(filename_in)
      iunit = get_filehandle()
      open(iunit,file=filename,status='old',form='unformatted',iostat=io)
      if(io /= 0) then
        write(*,'(1x,2a)') "[DC-LCFO-W90-GLOBAL] failed to open checkpoint: ", trim(filename)
        stop "DC-LCFO Wannier export: missing Wannier90 checkpoint"
      end if
      read(iunit) header_chk
      read(iunit) num_bands_chk
      read(iunit) num_exclude_bands_chk
      allocate(exclude_bands_chk(max(0,num_exclude_bands_chk)))
      if(num_exclude_bands_chk > 0) then
        read(iunit) (exclude_bands_chk(i), i=1,num_exclude_bands_chk)
      else
        read(iunit)
      end if
      read(iunit) ((real_lattice_chk(i,j), i=1,3), j=1,3)
      read(iunit) ((recip_lattice_chk(i,j), i=1,3), j=1,3)
      read(iunit) num_kpts_chk
      read(iunit) (mp_grid_chk(i), i=1,3)
      allocate(kpt_latt_chk(3,num_kpts_chk))
      read(iunit) ((kpt_latt_chk(i,nkp), i=1,3), nkp=1,num_kpts_chk)
      read(iunit) nntot_chk
      read(iunit) num_wann_chk
      read(iunit) checkpoint_chk
      read(iunit) have_disentangled_chk

      if(num_kpts_chk /= 1) stop "DC-LCFO Wannier export: only Gamma checkpoint is supported."
      if(have_disentangled_chk) then
        read(iunit) omega_invariant_chk
        allocate(lwindow(num_bands_chk,num_kpts_chk), ndimwin(num_kpts_chk))
        read(iunit) ((lwindow(i,nkp), i=1,num_bands_chk), nkp=1,num_kpts_chk)
        read(iunit) (ndimwin(nkp), nkp=1,num_kpts_chk)
        allocate(u_matrix_opt(num_bands_chk,num_wann_chk,num_kpts_chk))
        read(iunit) (((u_matrix_opt(i,j,nkp), i=1,num_bands_chk), j=1,num_wann_chk), nkp=1,num_kpts_chk)
      end if

      allocate(u_matrix(num_wann_chk,num_wann_chk,num_kpts_chk))
      read(iunit) (((u_matrix(i,j,k), i=1,num_wann_chk), j=1,num_wann_chk), k=1,num_kpts_chk)
      allocate(m_matrix(num_wann_chk,num_wann_chk,nntot_chk,num_kpts_chk))
      read(iunit) ((((m_matrix(i,j,k,nkp), i=1,num_wann_chk), j=1,num_wann_chk), k=1,nntot_chk), nkp=1,num_kpts_chk)
      allocate(center_aa(3,num_wann_chk), spread_aa2(num_wann_chk))
      read(iunit) ((center_aa(i,j), i=1,3), j=1,num_wann_chk)
      read(iunit) (spread_aa2(i), i=1,num_wann_chk)
      close(iunit)

      allocate(v_matrix(num_bands_chk,num_wann_chk))
      v_matrix = (0d0,0d0)
      if(have_disentangled_chk) then
        v_matrix(1:num_bands_chk,1:num_wann_chk) = matmul(u_matrix_opt(:,:,1), u_matrix(:,:,1))
      else
        if(num_bands_chk < num_wann_chk) &
          stop "DC-LCFO Wannier export: invalid checkpoint without disentanglement."
        v_matrix(1:num_wann_chk,1:num_wann_chk) = u_matrix(:,:,1)
      end if
      write(*,'(1x,a,i0,a,i0,a,l1)') "[DC-LCFO-W90-GLOBAL] read checkpoint: bands=", &
        num_bands_chk, " wann=", num_wann_chk, " disentangled=", have_disentangled_chk

      deallocate(exclude_bands_chk, kpt_latt_chk, u_matrix, m_matrix)
      if(allocated(lwindow)) deallocate(lwindow)
      if(allocated(ndimwin)) deallocate(ndimwin)
      if(allocated(u_matrix_opt)) deallocate(u_matrix_opt)
    end subroutine read_wannier90_checkpoint_transform

    integer function find_owner_fragment_from_center(center_bohr) result(owner)
      implicit none
      real(8), intent(in) :: center_bohr(3)
      integer :: ifrag_try, axis, idx0, idx1
      real(8) :: frag_center(3), dist2, best_dist2, delta_axis, length_axis
      integer :: nxyz_domain(3)

      owner = 1
      best_dist2 = huge(1d0)
      do ifrag_try=1,dc%n_frag
        call get_fragment_domain(dc, ifrag_try, nxyz_domain)
        do axis=1,3
          idx0 = dc%ixyz_frag(axis,ifrag_try)
          idx1 = idx0 + nxyz_domain(axis) - 1
          frag_center(axis) = 0.5d0 * (dc%lg_tot%coordinate(idx0,axis) + dc%lg_tot%coordinate(idx1,axis))
        end do
        dist2 = 0d0
        do axis=1,3
          length_axis = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
            + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
          delta_axis = periodic_delta(center_bohr(axis) - frag_center(axis), length_axis)
          dist2 = dist2 + delta_axis * delta_axis
        end do
        if(dist2 < best_dist2) then
          best_dist2 = dist2
          owner = ifrag_try
        end if
      end do
    end function find_owner_fragment_from_center

    subroutine rebalance_wannier_owner_fragments(center_bohr, owner_frag, num_wann)
      implicit none
      integer, intent(in) :: num_wann
      real(8), intent(in) :: center_bohr(3,num_wann)
      integer, intent(inout) :: owner_frag(num_wann)
      integer :: target_base, target_rem, ifrag, iw, donor, receiver, moved, iw_best
      integer, allocatable :: target(:), count_frag(:)
      real(8) :: dist2, best_dist2

      if(dc%n_frag <= 0) return
      allocate(target(dc%n_frag), count_frag(dc%n_frag))
      target_base = num_wann / dc%n_frag
      target_rem = mod(num_wann, dc%n_frag)
      do ifrag=1,dc%n_frag
        target(ifrag) = target_base
        if(ifrag <= target_rem) target(ifrag) = target(ifrag) + 1
        count_frag(ifrag) = count(owner_frag(1:num_wann) == ifrag)
      end do

      moved = 0
      do
        receiver = 0
        donor = 0
        do ifrag=1,dc%n_frag
          if(count_frag(ifrag) < target(ifrag)) then
            receiver = ifrag
            exit
          end if
        end do
        if(receiver == 0) exit
        do ifrag=1,dc%n_frag
          if(count_frag(ifrag) > target(ifrag)) then
            donor = ifrag
            exit
          end if
        end do
        if(donor == 0) exit

        iw_best = 0
        best_dist2 = huge(1d0)
        do iw=1,num_wann
          if(owner_frag(iw) /= donor) cycle
          dist2 = distance_to_fragment_center(center_bohr(1:3,iw), receiver)
          if(dist2 < best_dist2) then
            best_dist2 = dist2
            iw_best = iw
          end if
        end do
        if(iw_best == 0) exit
        owner_frag(iw_best) = receiver
        count_frag(donor) = count_frag(donor) - 1
        count_frag(receiver) = count_frag(receiver) + 1
        moved = moved + 1
      end do

      if(dc%id_tot == 0 .and. moved > 0) then
        write(*,'(1x,a,i0,a,i0,a,i0)') &
          "[DC-LCFO-W90-GLOBAL] rebalanced Wannier owners: moved=", moved, &
          " min_count=", minval(count_frag), " max_count=", maxval(count_frag)
      end if
      deallocate(target, count_frag)
    end subroutine rebalance_wannier_owner_fragments

    real(8) function distance_to_fragment_center(center_bohr, ifrag_target) result(dist2)
      implicit none
      real(8), intent(in) :: center_bohr(3)
      integer, intent(in) :: ifrag_target
      integer :: axis, idx0, idx1
      integer :: nxyz_domain(3)
      real(8) :: frag_center(3), delta_axis, length_axis

      call get_fragment_domain(dc, ifrag_target, nxyz_domain)
      do axis=1,3
        idx0 = dc%ixyz_frag(axis,ifrag_target)
        idx1 = idx0 + nxyz_domain(axis) - 1
        frag_center(axis) = 0.5d0 * (dc%lg_tot%coordinate(idx0,axis) + dc%lg_tot%coordinate(idx1,axis))
      end do
      dist2 = 0d0
      do axis=1,3
        length_axis = dc%lg_tot%coordinate(dc%lg_tot%num(axis),axis) &
          + (dc%lg_tot%coordinate(2,axis) - dc%lg_tot%coordinate(1,axis))
        delta_axis = periodic_delta(center_bohr(axis) - frag_center(axis), length_axis)
        dist2 = dist2 + delta_axis * delta_axis
      end do
    end function distance_to_fragment_center

    subroutine run_sawf_local_preprocessing(resolved_command,bundles)
      use communication,only:comm_sync_all
      character(*),intent(in)::resolved_command
      type(t_sawf_seed_bundle),intent(in)::bundles(:)
      integer::bundle
      call comm_sync_all(dc%icomm_tot)
      if(is_wannier90_export_only_command(resolved_command))then
        if(dc%id_tot==0)write(*,'(1x,a)') &
          '[DC-LCFO-SAWF-PREPROCESS] export-only: local WIN/EIG/AMN staged; NNKP/MMN deferred'
        return
      end if
      do bundle=1,size(bundles)
        call run_wannier90_seed_files(resolved_command,bundles(bundle)%seedname, &
          bundles(bundle)%directory,preprocess_only=.true.)
      end do
    end subroutine run_sawf_local_preprocessing

    subroutine run_sawf_local_wannier(resolved_command,bundles)
      character(*),intent(in)::resolved_command
      type(t_sawf_seed_bundle),intent(in)::bundles(:)
      integer::bundle
      do bundle=1,size(bundles)
        call run_wannier90_seed_files(resolved_command,bundles(bundle)%seedname, &
          bundles(bundle)%directory)
      end do
    end subroutine run_sawf_local_wannier

    subroutine run_wannier90_seed_files(resolved_command,seedname_in,directory_in,preprocess_only)
      use communication, only: comm_sync_all, comm_bcast
      use salmon_global, only: sysname
      implicit none
      character(*), intent(in) :: resolved_command
      character(*), intent(in), optional :: seedname_in,directory_in
      logical,intent(in),optional::preprocess_only
      character(1024) :: seedname
      character(4096) :: command_line, change_dir_command
      character(512) :: command_message
      integer :: command_failure
      logical :: command_ok
      logical :: run_preprocess

      seedname = trim(sysname)
      run_preprocess=.false.;if(present(preprocess_only))run_preprocess=preprocess_only
      if(present(seedname_in))seedname=trim(seedname_in)
      call comm_sync_all(dc%icomm_tot)
      if(is_wannier90_reuse_command(resolved_command)) then
        if(dc%id_tot == 0) write(*,'(1x,a)') &
          '[DC-LCFO-WANNIER] skip external Wannier90; reuse existing checkpoint.'
        call comm_sync_all(dc%icomm_tot)
        return
      end if
      command_failure = 0; command_message = ''
      if(dc%id_tot == 0) then
        change_dir_command = 'cd '//trim(shell_quote(dc%base_directory))//' && '
        if(present(directory_in))change_dir_command='cd '//trim(shell_quote(directory_in))//' && '
        command_line = trim(change_dir_command)//trim(mpi_clean_env_prefix())//' '//trim(resolved_command)
        if(run_preprocess)command_line=trim(command_line)//' -pp'
        command_line=trim(command_line)//' '//trim(shell_quote(seedname))
        write(*,'(1x,a,1x,a)') "[DC-LCFO-WANNIER] run:", trim(command_line)
        call execute_wannier90_command(command_line,command_ok,command_message)
        command_failure = merge(0,1,command_ok)
      end if
      call comm_bcast(command_failure,dc%icomm_tot,0)
      call comm_bcast(command_message,dc%icomm_tot,0)
      if(command_failure /= 0) call lcfo_sawf_fatal( &
        'DC-LCFO Wannier export failed collectively: '//trim(command_message))
      if(dc%id_tot == 0) write(*,'(1x,a)') '[DC-LCFO-WANNIER] Wannier90 completed.'
      call comm_sync_all(dc%icomm_tot)
    end subroutine run_wannier90_seed_files

    function mpi_clean_env_prefix() result(prefix)
      implicit none
      character(2048) :: prefix

      prefix = 'env' &
        //' -u OMPI_COMM_WORLD_SIZE' &
        //' -u OMPI_COMM_WORLD_RANK' &
        //' -u OMPI_COMM_WORLD_LOCAL_SIZE' &
        //' -u OMPI_COMM_WORLD_LOCAL_RANK' &
        //' -u OMPI_UNIVERSE_SIZE' &
        //' -u PMIX_NAMESPACE' &
        //' -u PMIX_RANK' &
        //' -u PMIX_SERVER_URI' &
        //' -u PMIX_SERVER_URI2' &
        //' -u PMIX_SERVER_URI21' &
        //' -u PMIX_SECURITY_MODE' &
        //' -u PMIX_GDS_MODULE' &
        //' '
    end function mpi_clean_env_prefix

    function shell_quote(text) result(quoted)
      implicit none
      character(*), intent(in) :: text
      character(2048) :: quoted

      quoted = "'"//trim(text)//"'"
    end function shell_quote

    integer function determine_wannier_num_bands() result(nband)
      use salmon_global, only: wannier_num_bands, wannier_num_wann, wannier_dis_win_max
      implicit none
      integer :: iband

      if(wannier_num_bands > 0) then
        if(wannier_num_bands > dc%nstate_tot) &
          stop "DC-LCFO Wannier export: wannier_num_bands exceeds DC-LCFO state count."
        nband = wannier_num_bands
      else if(wannier_dis_win_max > 0d0) then
        nband = 0
        do iband=1,dc%nstate_tot
          if(esp_tot(iband,1) <= wannier_dis_win_max) nband = iband
        end do
        nband = max(nband, wannier_num_wann)
      else
        nband = dc%nstate_tot
      end if
      if(nband < wannier_num_wann) &
        stop "DC-LCFO Wannier export: selected band count is smaller than wannier_num_wann."
    end function determine_wannier_num_bands

    subroutine write_wannier_mmn_file(nband_wann)
      use salmon_global, only: sysname
      use filesystem, only: get_filehandle
      implicit none
      integer, intent(in) :: nband_wann
      integer, parameter :: nntot_gamma = 3
      integer, parameter :: mmn_target_chunk_elems = 1000000
      integer :: gvec(3,nntot_gamma)
      integer :: nxyz_domain(3), iunit, inn, ibasis, jbasis
      integer :: ix, iy, iz, ixg, iyg, izg, iband, j0, j1, nchunk, jband
      integer :: state_chunk_size, ispin_local
      real(8) :: x, y, z, phase_arg, cphase, sphase
      real(8) :: gx, gy, gz
      real(8) :: a1(3), a2(3), a3(3)
      real(8), allocatable :: s_re(:,:), s_im(:,:), tmp_re(:,:), tmp_im(:,:)
      real(8), allocatable :: m_re_local(:,:), m_im_local(:,:), m_re_sum(:,:), m_im_sum(:,:)
      character(256) :: mmnfile

      gvec(:,1) = (/ 1, 0, 0 /)
      gvec(:,2) = (/ 0, 1, 0 /)
      gvec(:,3) = (/ 0, 0, 1 /)
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      call get_lattice_vectors(a1, a2, a3)

      mmnfile = trim(dc%base_directory)//trim(sysname)//".mmn"
      if(dc%id_tot == 0) then
        iunit = get_filehandle()
        open(iunit,file=mmnfile,status='replace')
        write(iunit,'(a)') "SALMON DC-LCFO generated overlaps"
        write(iunit,'(3i10)') nband_wann, 1, nntot_gamma
      end if

      state_chunk_size = max(1, min(nband_wann, &
        max(1, mmn_target_chunk_elems / max(1, nband_wann))))
      allocate(s_re(dc%nstate_frag,dc%nstate_frag), s_im(dc%nstate_frag,dc%nstate_frag))
      allocate(tmp_re(dc%nstate_frag,state_chunk_size), tmp_im(dc%nstate_frag,state_chunk_size))
      allocate(m_re_local(nband_wann,state_chunk_size), m_im_local(nband_wann,state_chunk_size))
      allocate(m_re_sum(nband_wann,state_chunk_size), m_im_sum(nband_wann,state_chunk_size))

      ispin_local = 1
      do inn=1,nntot_gamma
        if(dc%id_tot == 0) write(iunit,'(5i8)') 1, 1, gvec(1,inn), gvec(2,inn), gvec(3,inn)

        s_re(:,:) = 0d0
        s_im(:,:) = 0d0
        if(dc%id_frag == 0) then
          call reciprocal_vector_from_index(gvec(:,inn), a1, a2, a3, gx, gy, gz)
          do iz=1,nxyz_domain(3)
            izg = dc%jxyz_tot(iz,3)
            z = dc%lg_tot%coordinate(izg,3)
            do iy=1,nxyz_domain(2)
              iyg = dc%jxyz_tot(iy,2)
              y = dc%lg_tot%coordinate(iyg,2)
              do ix=1,nxyz_domain(1)
                ixg = dc%jxyz_tot(ix,1)
                x = dc%lg_tot%coordinate(ixg,1)
                phase_arg = gx*x + gy*y + gz*z
                cphase = cos(phase_arg)
                sphase = -sin(phase_arg)
                do jbasis=1,n_basis(dc%i_frag,ispin_local)
                  do ibasis=1,n_basis(dc%i_frag,ispin_local)
                    s_re(ibasis,jbasis) = s_re(ibasis,jbasis) &
                      + f_basis(ix,iy,iz,ispin_local,ibasis) &
                      * f_basis(ix,iy,iz,ispin_local,jbasis) * cphase * hvol
                    s_im(ibasis,jbasis) = s_im(ibasis,jbasis) &
                      + f_basis(ix,iy,iz,ispin_local,ibasis) &
                      * f_basis(ix,iy,iz,ispin_local,jbasis) * sphase * hvol
                  end do
                end do
              end do
            end do
          end do
        end if

        do j0=1,nband_wann,state_chunk_size
          j1 = min(nband_wann, j0 + state_chunk_size - 1)
          nchunk = j1 - j0 + 1
          m_re_local(:,:) = 0d0
          m_im_local(:,:) = 0d0
          if(dc%id_frag == 0) then
            tmp_re(1:dc%nstate_frag,1:nchunk) = &
              matmul(s_re(1:dc%nstate_frag,1:dc%nstate_frag), &
              coef_wf(1:dc%nstate_frag,j0:j1,1))
            tmp_im(1:dc%nstate_frag,1:nchunk) = &
              matmul(s_im(1:dc%nstate_frag,1:dc%nstate_frag), &
              coef_wf(1:dc%nstate_frag,j0:j1,1))
            m_re_local(1:nband_wann,1:nchunk) = &
              matmul(transpose(coef_wf(1:dc%nstate_frag,1:nband_wann,1)), &
              tmp_re(1:dc%nstate_frag,1:nchunk))
            m_im_local(1:nband_wann,1:nchunk) = &
              matmul(transpose(coef_wf(1:dc%nstate_frag,1:nband_wann,1)), &
              tmp_im(1:dc%nstate_frag,1:nchunk))
          end if
          call comm_summation(m_re_local,m_re_sum,nband_wann*state_chunk_size,dc%icomm_tot)
          call comm_summation(m_im_local,m_im_sum,nband_wann*state_chunk_size,dc%icomm_tot)

          if(dc%id_tot == 0) then
            do jband=1,nchunk
              do iband=1,nband_wann
                write(iunit,'(2(1x,es23.15))') m_re_sum(iband,jband), m_im_sum(iband,jband)
              end do
            end do
          end if
        end do
      end do

      if(dc%id_tot == 0) close(iunit)
      deallocate(s_re, s_im, tmp_re, tmp_im, m_re_local, m_im_local, m_re_sum, m_im_sum)
    end subroutine write_wannier_mmn_file

    subroutine write_wannier_amn_file(nband_wann)
      use salmon_global, only: izatom, sysname, &
        wannier_projection, wannier_num_wann, wannier_projection_width
      use filesystem, only: get_filehandle
      implicit none
      integer, intent(in) :: nband_wann
      integer :: nproj, nproj_csp3, nproj_pseudo, chunk_size, p0, p1, nchunk, iunit
      integer :: iband, ip, ibasis, ix, iy, iz, ispin_local
      integer :: ip_global, ip_base, ip_shell
      integer :: ixg, iyg, izg
      integer :: nxyz_domain(3)
      integer, parameter :: amn_target_chunk_elems = 1000000
      integer, allocatable :: proj_atom(:), proj_hybrid(:)
      integer, allocatable :: proj_atom_pseudo(:), proj_l_pseudo(:), proj_m_pseudo(:)
      real(8) :: x, y, z, gval, projection_width_eff
      real(8), allocatable :: bproj(:,:), a_local(:,:), a_sum(:,:)
      real(8), allocatable :: norm_local(:), norm_sum(:)
      character(256) :: amnfile

      if(.not. is_sp3_projection(trim(wannier_projection)) .and. &
         .not. is_pseudo_channel_projection(trim(wannier_projection)) .and. &
         .not. is_bond_center_projection(trim(wannier_projection))) &
        stop "DC-LCFO Wannier export: supported projections are C:sp3, Si:sp3, pseudo_channels, and bond_centers."
      if(is_pseudo_channel_projection(trim(wannier_projection))) then
        call write_wannier_amn_file_pseudo_channels(nband_wann)
        return
      end if
      if(is_bond_center_projection(trim(wannier_projection))) then
        call write_wannier_amn_file_bond_centers(nband_wann)
        return
      end if

      nproj_csp3 = count_c_sp3_projections()
      if(nproj_csp3 <= 0) stop "DC-LCFO Wannier export: no atoms found for sp3 projections."
      if(nproj_csp3 > wannier_num_wann) &
        stop "DC-LCFO Wannier export: sp3 projection count exceeds wannier_num_wann."
      if(wannier_num_wann > nband_wann) &
        stop "DC-LCFO Wannier export: wannier_num_wann must not exceed exported band count."
      nproj = wannier_num_wann

      allocate(proj_atom(nproj_csp3), proj_hybrid(nproj_csp3))
      call build_c_sp3_projection_map(proj_atom, proj_hybrid)
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)

      amnfile = trim(dc%base_directory)//trim(sysname)//".amn"
      if(dc%id_tot == 0) then
        iunit = get_filehandle()
        open(iunit,file=amnfile,status='replace')
        write(iunit,'(a)') "SALMON DC-LCFO generated projections"
        write(iunit,'(3i10)') nband_wann, 1, nproj
      end if

      chunk_size = max(1, min(nproj, max(1, amn_target_chunk_elems / max(1, nband_wann))))
      allocate(bproj(dc%nstate_frag,chunk_size))
      allocate(a_local(nband_wann,chunk_size))
      allocate(a_sum(nband_wann,chunk_size))
      allocate(norm_local(chunk_size), norm_sum(chunk_size))

      do p0=1,nproj,chunk_size
        p1 = min(nproj, p0 + chunk_size - 1)
        nchunk = p1 - p0 + 1
        bproj(:,:) = 0d0
        norm_local(:) = 0d0

        ispin_local = 1
        if(dc%id_frag == 0) then
          do ip=1,nchunk
            ip_global = p0 + ip - 1
            ip_base = modulo(ip_global - 1, nproj_csp3) + 1
            ip_shell = (ip_global - 1) / nproj_csp3
            projection_width_eff = wannier_projection_width / sqrt(1d0 + 0.75d0 * dble(ip_shell))
            do iz=1,nxyz_domain(3)
              izg = dc%jxyz_tot(iz,3)
              z = dc%lg_tot%coordinate(izg,3)
              do iy=1,nxyz_domain(2)
                iyg = dc%jxyz_tot(iy,2)
                y = dc%lg_tot%coordinate(iyg,2)
                do ix=1,nxyz_domain(1)
                  ixg = dc%jxyz_tot(ix,1)
                  x = dc%lg_tot%coordinate(ixg,1)
                  gval = c_sp3_projection_value(x, y, z, proj_atom(ip_base), &
                    proj_hybrid(ip_base), projection_width_eff)
                  if(abs(gval) <= 0d0) cycle
                  norm_local(ip) = norm_local(ip) + gval * gval * hvol
                  do ibasis=1,n_basis(dc%i_frag,ispin_local)
                    bproj(ibasis,ip) = bproj(ibasis,ip) &
                      + f_basis(ix,iy,iz,ispin_local,ibasis) * gval * hvol
                  end do
                end do
              end do
            end do
          end do
        end if

        a_local(:,:) = 0d0
        if(dc%id_frag == 0) then
          a_local(1:nband_wann,1:nchunk) = &
            matmul(transpose(coef_wf(1:dc%nstate_frag,1:nband_wann,1)), &
            bproj(1:dc%nstate_frag,1:nchunk))
        end if
        call comm_summation(a_local,a_sum,nband_wann*chunk_size,dc%icomm_tot)
        call comm_summation(norm_local,norm_sum,chunk_size,dc%icomm_tot)

        if(dc%id_tot == 0) then
          do ip=1,nchunk
            ip_global = p0 + ip - 1
            if(norm_sum(ip) <= 1d-300) &
              stop "DC-LCFO Wannier export: zero norm projection in .amn generation."
            a_sum(1:nband_wann,ip) = a_sum(1:nband_wann,ip) / sqrt(norm_sum(ip))
            do iband=1,nband_wann
              write(iunit,'(3i8,2(1x,es23.15))') iband, ip_global, 1, a_sum(iband,ip), 0d0
            end do
          end do
        end if
      end do

      if(dc%id_tot == 0) close(iunit)
      deallocate(proj_atom, proj_hybrid)
      deallocate(bproj, a_local, a_sum, norm_local, norm_sum)
    end subroutine write_wannier_amn_file

    subroutine write_wannier_amn_file_pseudo_channels(nband_wann)
      use salmon_global, only: sysname, wannier_num_wann, wannier_site_symmetry
      use filesystem, only: get_filehandle
      use communication, only: comm_get_max
      implicit none
      integer, intent(in) :: nband_wann
      integer :: allocation_failure, allocation_status, zero_norm_failure
      integer :: nproj_raw, nproj, chunk_size, p0, p1, nchunk, iunit
      integer :: ip, ip_raw, iband, projection_lmax, selected_count
      integer, parameter :: amn_target_chunk_elems = 1000000
      integer, allocatable :: proj_atom_raw(:), proj_l_raw(:), proj_m_raw(:)
      integer, allocatable :: raw_index(:), selected_raw(:)
      real(8), allocatable :: projectability(:), a_sum(:,:), norm_sum(:)
      logical :: complete_shell
      character(512) :: shell_message
      character(256) :: allocation_message, amnfile

      if(wannier_num_wann > nband_wann) &
        call lcfo_sawf_fatal( &
          'DC-LCFO Wannier export: wannier_num_wann must not exceed exported band count.')

      nproj_raw = count_pseudo_channel_ao_candidates()
      if(nproj_raw <= 0) call lcfo_sawf_fatal( &
        'DC-LCFO Wannier export: no pseudo-channel AO candidates were generated.')
      if(nproj_raw < wannier_num_wann) &
        call lcfo_sawf_fatal( &
          'DC-LCFO Wannier export: pseudo-channel candidate count is smaller than wannier_num_wann.')
      nproj = wannier_num_wann
      call sawf_projection_shell_lmax(dc%system_tot%nion, nproj, projection_lmax, &
        complete_shell, shell_message)

      allocation_failure = 0
      allocate(proj_atom_raw(nproj_raw), stat=allocation_status, errmsg=allocation_message)
      if(allocation_status /= 0) allocation_failure = 1
      if(allocation_failure == 0) then
        allocate(proj_l_raw(nproj_raw), stat=allocation_status, errmsg=allocation_message)
        if(allocation_status /= 0) allocation_failure = 1
      end if
      if(allocation_failure == 0) then
        allocate(proj_m_raw(nproj_raw), stat=allocation_status, errmsg=allocation_message)
        if(allocation_status /= 0) allocation_failure = 1
      end if
      if(allocation_failure == 0) then
        allocate(projectability(nproj_raw), stat=allocation_status, errmsg=allocation_message)
        if(allocation_status /= 0) allocation_failure = 1
      end if
      if(allocation_failure == 0) then
        allocate(selected_raw(nproj), stat=allocation_status, errmsg=allocation_message)
        if(allocation_status /= 0) allocation_failure = 1
      end if
      call comm_get_max(allocation_failure, dc%icomm_tot)
      if(allocation_failure /= 0) then
        call deallocate_pseudo_channel_amn_arrays(proj_atom_raw, proj_l_raw, proj_m_raw, &
          projectability, selected_raw, raw_index, a_sum, norm_sum)
        call lcfo_sawf_fatal( &
          'DC-LCFO Wannier export: pseudo-channel metadata allocation failed on one or more ranks.')
      end if
      call build_pseudo_channel_ao_candidate_map(proj_atom_raw, proj_l_raw, proj_m_raw)

      if(dc%id_tot == 0) then
        write(*,'(1x,a,i0,a,i0)') "[DC-LCFO-WANNIER] pseudo_channels candidates=", &
          nproj_raw, " target_wann=", nproj
        write(*,'(1x,a)') "[DC-LCFO-WANNIER] pseudo_channels selection: projectability top candidates"
      end if

      chunk_size = max(1, min(nproj_raw, max(1, amn_target_chunk_elems / max(1, nband_wann))))
      allocation_failure = 0
      allocate(raw_index(chunk_size), stat=allocation_status, errmsg=allocation_message)
      if(allocation_status /= 0) allocation_failure = 1
      if(allocation_failure == 0) then
        allocate(a_sum(nband_wann,chunk_size), stat=allocation_status, errmsg=allocation_message)
        if(allocation_status /= 0) allocation_failure = 1
      end if
      if(allocation_failure == 0) then
        allocate(norm_sum(chunk_size), stat=allocation_status, errmsg=allocation_message)
        if(allocation_status /= 0) allocation_failure = 1
      end if
      call comm_get_max(allocation_failure, dc%icomm_tot)
      if(allocation_failure /= 0) then
        call deallocate_pseudo_channel_amn_arrays(proj_atom_raw, proj_l_raw, proj_m_raw, &
          projectability, selected_raw, raw_index, a_sum, norm_sum)
        call lcfo_sawf_fatal( &
          'DC-LCFO Wannier export: pseudo-channel work allocation failed on one or more ranks.')
      end if
      projectability(:) = 0d0

      do p0=1,nproj_raw,chunk_size
        p1 = min(nproj_raw, p0 + chunk_size - 1)
        nchunk = p1 - p0 + 1
        do ip=1,nchunk
          raw_index(ip) = p0 + ip - 1
        end do
        call compute_pseudo_channel_projection_chunk(nband_wann, nchunk, raw_index, &
          proj_atom_raw, proj_l_raw, proj_m_raw, a_sum, norm_sum)
        zero_norm_failure = 0
        if(any(norm_sum(1:nchunk) <= 1d-300)) zero_norm_failure = 1
        call comm_get_max(zero_norm_failure, dc%icomm_tot)
        if(zero_norm_failure /= 0) call lcfo_sawf_fatal( &
          'DC-LCFO Wannier export: zero norm pseudo-channel projection.')
        do ip=1,nchunk
          ip_raw = raw_index(ip)
          a_sum(1:nband_wann,ip) = a_sum(1:nband_wann,ip) / sqrt(norm_sum(ip))
          projectability(ip_raw) = sum(a_sum(1:nband_wann,ip)**2)
        end do
      end do

      if(complete_shell) then
        selected_count = 0
        do ip_raw=1,nproj_raw
          if(proj_l_raw(ip_raw) > projection_lmax) cycle
          selected_count = selected_count + 1
          if(selected_count <= nproj) selected_raw(selected_count) = ip_raw
        end do
        if(selected_count /= nproj) call lcfo_sawf_fatal( &
          'DC-LCFO Wannier export: complete pseudo-channel shell count is inconsistent.')
        if(dc%id_tot == 0) write(*,'(1x,a,i0)') &
          '[DC-LCFO-WANNIER] pseudo_channels selection: complete atom-major shell lmax=', &
          projection_lmax
      else
        call select_top_projectors(projectability, nproj_raw, nproj, selected_raw)
      end if
      call write_pseudo_channel_projection_diagnostics(nproj_raw, nproj, proj_l_raw, proj_m_raw, &
        projectability, selected_raw)

      if(dc%id_tot == 0) then
        do ip=1,nproj
          ip_raw = selected_raw(ip)
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,es14.6)') &
            "[DC-LCFO-WANNIER] pseudo selected ip=", ip, &
            " raw=", ip_raw, " atom=", proj_atom_raw(ip_raw), &
            " l=", proj_l_raw(ip_raw), " m=", proj_m_raw(ip_raw), &
            " projectability=", projectability(ip_raw)
        end do
      end if

      amnfile = trim(dc%base_directory)//trim(sysname)//".amn"
      if(dc%id_tot == 0) then
        iunit = get_filehandle()
        open(iunit,file=amnfile,status='replace')
        write(iunit,'(a)') "SALMON DC-LCFO generated pseudo-channel projections"
        write(iunit,'(3i10)') nband_wann, 1, nproj
      end if

      do p0=1,nproj,chunk_size
        p1 = min(nproj, p0 + chunk_size - 1)
        nchunk = p1 - p0 + 1
        do ip=1,nchunk
          raw_index(ip) = selected_raw(p0 + ip - 1)
        end do
        call compute_pseudo_channel_projection_chunk(nband_wann, nchunk, raw_index, &
          proj_atom_raw, proj_l_raw, proj_m_raw, a_sum, norm_sum)
        zero_norm_failure = 0
        if(any(norm_sum(1:nchunk) <= 1d-300)) zero_norm_failure = 1
        call comm_get_max(zero_norm_failure, dc%icomm_tot)
        if(zero_norm_failure /= 0) call lcfo_sawf_fatal( &
          'DC-LCFO Wannier export: zero norm selected pseudo-channel projection.')
        if(dc%id_tot == 0) then
          do ip=1,nchunk
            a_sum(1:nband_wann,ip) = a_sum(1:nband_wann,ip) / sqrt(norm_sum(ip))
            do iband=1,nband_wann
              write(iunit,'(3i8,2(1x,es23.15))') iband, p0 + ip - 1, 1, a_sum(iband,ip), 0d0
            end do
          end do
        end if
      end do

      if(dc%id_tot == 0) close(iunit)
      call deallocate_pseudo_channel_amn_arrays(proj_atom_raw, proj_l_raw, proj_m_raw, &
        projectability, selected_raw, raw_index, a_sum, norm_sum)
    end subroutine write_wannier_amn_file_pseudo_channels

    subroutine deallocate_pseudo_channel_amn_arrays(proj_atom, proj_l, proj_m, &
        projectability, selected, raw_index, a_sum, norm_sum)
      implicit none
      integer, allocatable, intent(inout) :: proj_atom(:), proj_l(:), proj_m(:)
      integer, allocatable, intent(inout) :: selected(:), raw_index(:)
      real(8), allocatable, intent(inout) :: projectability(:), a_sum(:,:), norm_sum(:)

      if(allocated(proj_atom)) deallocate(proj_atom)
      if(allocated(proj_l)) deallocate(proj_l)
      if(allocated(proj_m)) deallocate(proj_m)
      if(allocated(projectability)) deallocate(projectability)
      if(allocated(selected)) deallocate(selected)
      if(allocated(raw_index)) deallocate(raw_index)
      if(allocated(a_sum)) deallocate(a_sum)
      if(allocated(norm_sum)) deallocate(norm_sum)
    end subroutine deallocate_pseudo_channel_amn_arrays

    subroutine write_wannier_amn_file_bond_centers(nband_wann)
      use salmon_global, only: sysname, wannier_num_wann, wannier_projection_width
      use filesystem, only: get_filehandle
      implicit none
      integer, intent(in) :: nband_wann
      integer :: nproj, chunk_size, p0, p1, nchunk, iunit
      integer :: iband, ip, ibasis, ix, iy, iz, ispin_local
      integer :: ixg, iyg, izg, nxyz_domain(3)
      integer, parameter :: amn_target_chunk_elems = 1000000
      real(8) :: x, y, z, gval, projection_width_eff
      real(8), allocatable :: bond_center_bohr(:,:), bproj(:,:), a_local(:,:), a_sum(:,:)
      real(8), allocatable :: norm_local(:), norm_sum(:)
      character(256) :: amnfile

      if(wannier_num_wann > nband_wann) &
        stop "DC-LCFO Wannier export: wannier_num_wann must not exceed exported band count."

      nproj = wannier_num_wann
      call build_bond_center_projection_map(nproj, bond_center_bohr)
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)

      amnfile = trim(dc%base_directory)//trim(sysname)//".amn"
      if(dc%id_tot == 0) then
        write(*,'(1x,a,i0)') "[DC-LCFO-WANNIER] bond-center projections target_wann=", nproj
        iunit = get_filehandle()
        open(iunit,file=amnfile,status='replace')
        write(iunit,'(a)') "SALMON DC-LCFO generated bond-center projections"
        write(iunit,'(3i10)') nband_wann, 1, nproj
      end if

      chunk_size = max(1, min(nproj, max(1, amn_target_chunk_elems / max(1, nband_wann))))
      allocate(bproj(dc%nstate_frag,chunk_size))
      allocate(a_local(nband_wann,chunk_size))
      allocate(a_sum(nband_wann,chunk_size))
      allocate(norm_local(chunk_size), norm_sum(chunk_size))

      do p0=1,nproj,chunk_size
        p1 = min(nproj, p0 + chunk_size - 1)
        nchunk = p1 - p0 + 1
        bproj(:,:) = 0d0
        norm_local(:) = 0d0

        ispin_local = 1
        if(dc%id_frag == 0) then
          do ip=1,nchunk
            projection_width_eff = wannier_projection_width / &
              sqrt(1d0 + 0.75d0 * dble((p0 + ip - 2) / max(1, count_bond_center_candidates())))
            do iz=1,nxyz_domain(3)
              izg = dc%jxyz_tot(iz,3)
              z = dc%lg_tot%coordinate(izg,3)
              do iy=1,nxyz_domain(2)
                iyg = dc%jxyz_tot(iy,2)
                y = dc%lg_tot%coordinate(iyg,2)
                do ix=1,nxyz_domain(1)
                  ixg = dc%jxyz_tot(ix,1)
                  x = dc%lg_tot%coordinate(ixg,1)
                  gval = bond_center_projection_value(x, y, z, &
                    bond_center_bohr(1:3,p0 + ip - 1), projection_width_eff)
                  if(abs(gval) <= 0d0) cycle
                  norm_local(ip) = norm_local(ip) + gval * gval * hvol
                  do ibasis=1,n_basis(dc%i_frag,ispin_local)
                    bproj(ibasis,ip) = bproj(ibasis,ip) &
                      + f_basis(ix,iy,iz,ispin_local,ibasis) * gval * hvol
                  end do
                end do
              end do
            end do
          end do
        end if

        a_local(:,:) = 0d0
        if(dc%id_frag == 0) then
          a_local(1:nband_wann,1:nchunk) = &
            matmul(transpose(coef_wf(1:dc%nstate_frag,1:nband_wann,1)), &
            bproj(1:dc%nstate_frag,1:nchunk))
        end if
        call comm_summation(a_local,a_sum,nband_wann*chunk_size,dc%icomm_tot)
        call comm_summation(norm_local,norm_sum,chunk_size,dc%icomm_tot)

        if(dc%id_tot == 0) then
          do ip=1,nchunk
            if(norm_sum(ip) <= 1d-300) &
              stop "DC-LCFO Wannier export: zero norm bond-center projection."
            a_sum(1:nband_wann,ip) = a_sum(1:nband_wann,ip) / sqrt(norm_sum(ip))
            do iband=1,nband_wann
              write(iunit,'(3i8,2(1x,es23.15))') iband, p0 + ip - 1, 1, a_sum(iband,ip), 0d0
            end do
          end do
        end if
      end do

      if(dc%id_tot == 0) close(iunit)
      deallocate(bond_center_bohr, bproj, a_local, a_sum, norm_local, norm_sum)
    end subroutine write_wannier_amn_file_bond_centers

    subroutine compute_pseudo_channel_projection_chunk(nband_wann, nchunk, raw_index, &
        proj_atom, proj_l, proj_m, a_sum, norm_sum)
      use salmon_global, only: wannier_projection_width
      implicit none
      integer, intent(in) :: nband_wann, nchunk
      integer, intent(in) :: raw_index(:), proj_atom(:), proj_l(:), proj_m(:)
      real(8), intent(out) :: a_sum(:,:), norm_sum(:)
      integer :: ip, ip_raw, ibasis, ix, iy, iz, ixg, iyg, izg, ispin_local
      integer :: nxyz_domain(3)
      real(8) :: x, y, z, gval
      real(8), allocatable :: bproj(:,:), a_local(:,:), norm_local(:)

      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
      allocate(bproj(dc%nstate_frag,nchunk))
      allocate(a_local(nband_wann,nchunk))
      allocate(norm_local(nchunk))
      bproj(:,:) = 0d0
      a_local(:,:) = 0d0
      norm_local(:) = 0d0

      ispin_local = 1
      if(dc%id_frag == 0) then
        do ip=1,nchunk
          ip_raw = raw_index(ip)
          do iz=1,nxyz_domain(3)
            izg = dc%jxyz_tot(iz,3)
            z = dc%lg_tot%coordinate(izg,3)
            do iy=1,nxyz_domain(2)
              iyg = dc%jxyz_tot(iy,2)
              y = dc%lg_tot%coordinate(iyg,2)
              do ix=1,nxyz_domain(1)
                ixg = dc%jxyz_tot(ix,1)
                x = dc%lg_tot%coordinate(ixg,1)
                gval = pseudo_channel_projection_value(x, y, z, proj_atom(ip_raw), &
                  proj_l(ip_raw), proj_m(ip_raw), wannier_projection_width)
                if(abs(gval) <= 0d0) cycle
                norm_local(ip) = norm_local(ip) + gval * gval * hvol
                do ibasis=1,n_basis(dc%i_frag,ispin_local)
                  bproj(ibasis,ip) = bproj(ibasis,ip) &
                    + f_basis(ix,iy,iz,ispin_local,ibasis) * gval * hvol
                end do
              end do
            end do
          end do
        end do
        a_local(1:nband_wann,1:nchunk) = &
          matmul(transpose(coef_wf(1:dc%nstate_frag,1:nband_wann,1)), &
          bproj(1:dc%nstate_frag,1:nchunk))
      end if

      call comm_summation(a_local,a_sum,nband_wann*nchunk,dc%icomm_tot)
      call comm_summation(norm_local,norm_sum,nchunk,dc%icomm_tot)

      deallocate(bproj, a_local, norm_local)
    end subroutine compute_pseudo_channel_projection_chunk

    subroutine select_top_projectors(projectability, nraw, nselect, selected)
      implicit none
      integer, intent(in) :: nraw, nselect
      real(8), intent(in) :: projectability(nraw)
      integer, intent(out) :: selected(nselect)
      logical, allocatable :: used(:)
      integer :: isel, ip, best_ip
      real(8) :: best_val

      allocate(used(nraw))
      used(:) = .false.
      selected(:) = 0
      do isel=1,nselect
        best_ip = 0
        best_val = -1d0
        do ip=1,nraw
          if(used(ip)) cycle
          if(projectability(ip) > best_val) then
            best_val = projectability(ip)
            best_ip = ip
          end if
        end do
        if(best_ip <= 0) stop "DC-LCFO Wannier export: failed to select pseudo-channel projector."
        selected(isel) = best_ip
        used(best_ip) = .true.
      end do
      deallocate(used)
    end subroutine select_top_projectors

    subroutine write_pseudo_channel_projection_diagnostics(nraw, nselect, proj_l, proj_m, &
        projectability, selected)
      use filesystem, only: get_filehandle
      implicit none
      integer, intent(in) :: nraw, nselect
      integer, intent(in) :: proj_l(nraw), proj_m(nraw), selected(nselect)
      real(8), intent(in) :: projectability(nraw)
      integer :: iunit, ip, l, ip_raw
      integer :: selected_l(0:2), raw_l(0:2), count_l(0:2)
      logical, allocatable :: is_selected(:)
      real(8) :: sum_l(0:2), avg_l(0:2)
      real(8) :: min_selected, max_rejected
      character(256) :: diagfile

      if(dc%id_tot /= 0) return

      allocate(is_selected(nraw))
      is_selected(:) = .false.
      do ip=1,nselect
        if(selected(ip) >= 1 .and. selected(ip) <= nraw) is_selected(selected(ip)) = .true.
      end do

      selected_l(:) = 0
      raw_l(:) = 0
      count_l(:) = 0
      sum_l(:) = 0d0
      avg_l(:) = 0d0
      min_selected = huge(1d0)
      max_rejected = -huge(1d0)

      do ip=1,nraw
        l = proj_l(ip)
        if(l >= 0 .and. l <= 2) then
          raw_l(l) = raw_l(l) + 1
          count_l(l) = count_l(l) + 1
          sum_l(l) = sum_l(l) + projectability(ip)
        end if
        if(is_selected(ip)) then
          min_selected = min(min_selected, projectability(ip))
          if(l >= 0 .and. l <= 2) selected_l(l) = selected_l(l) + 1
        else
          max_rejected = max(max_rejected, projectability(ip))
        end if
      end do

      do l=0,2
        if(count_l(l) > 0) avg_l(l) = sum_l(l) / dble(count_l(l))
      end do
      if(nselect <= 0) min_selected = 0d0
      if(nraw <= nselect) max_rejected = 0d0

      diagfile = trim(dc%base_directory)//"pseudo_channel_projection_diag.csv"
      iunit = get_filehandle()
      open(iunit,file=diagfile,status='replace')
      write(iunit,'(a)') "nproj_raw,num_wann,raw_l0,raw_l1,raw_l2,selected_l0,selected_l1,selected_l2,"// &
        "avg_projectability_l0,avg_projectability_l1,avg_projectability_l2,"// &
        "min_selected_projectability,max_rejected_projectability"
      write(iunit,'(8(i0,","),5(es23.15,:,","))') nraw, nselect, raw_l(0), raw_l(1), raw_l(2), &
        selected_l(0), selected_l(1), selected_l(2), avg_l(0), avg_l(1), avg_l(2), &
        min_selected, max_rejected
      write(iunit,'(a)') "selected_index,raw_index,l,m,projectability"
      do ip=1,nselect
        ip_raw = selected(ip)
        write(iunit,'(4(i0,","),es23.15)') ip, ip_raw, proj_l(ip_raw), &
          proj_m(ip_raw), projectability(ip_raw)
      end do
      close(iunit)

      write(*,'(1x,a,a)') "[DC-LCFO-WANNIER] pseudo projection diagnostics: ", trim(diagfile)
      write(*,'(1x,a,3(i0,1x),a,3(i0,1x),a,es14.6,a,es14.6)') &
        "[DC-LCFO-WANNIER] pseudo l-count raw=", raw_l(0), raw_l(1), raw_l(2), &
        " selected=", selected_l(0), selected_l(1), selected_l(2), &
        " min_selected=", min_selected, " max_rejected=", max_rejected
      deallocate(is_selected)
    end subroutine write_pseudo_channel_projection_diagnostics

    logical function is_sp3_projection(text) result(enabled)
      implicit none
      character(*), intent(in) :: text

      enabled = sp3_projection_iz(text) > 0
    end function is_sp3_projection

    integer function sp3_projection_iz(text) result(target_iz)
      implicit none
      character(*), intent(in) :: text
      character(256) :: work

      work = adjustl(text)
      select case(trim(work))
      case('C:sp3', 'c:sp3')
        target_iz = 6
      case('Si:sp3', 'si:sp3', 'SI:sp3')
        target_iz = 14
      case default
        target_iz = 0
      end select
    end function sp3_projection_iz

    logical function is_pseudo_channel_projection(text) result(enabled)
      implicit none
      character(*), intent(in) :: text
      character(256) :: work

      work = adjustl(text)
      enabled = (trim(work) == 'pseudo_channels' .or. trim(work) == 'PSEUDO_CHANNELS')
    end function is_pseudo_channel_projection

    integer function count_c_sp3_projections() result(nproj)
      use salmon_global, only: izatom, wannier_projection
      implicit none
      integer :: ia, target_iz

      nproj = 0
      target_iz = sp3_projection_iz(trim(wannier_projection))
      do ia=1,dc%system_tot%nion
        if(izatom(dc%system_tot%kion(ia)) == target_iz) nproj = nproj + 4
      end do
    end function count_c_sp3_projections

    logical function is_bond_center_projection(text) result(enabled)
      implicit none
      character(*), intent(in) :: text
      character(256) :: work

      work = adjustl(text)
      enabled = (trim(work) == 'bond_centers' .or. trim(work) == 'BOND_CENTERS')
    end function is_bond_center_projection

    integer function count_bond_center_candidates() result(nbond)
      implicit none
      real(8) :: cutoff, min_dist

      call find_bond_center_cutoff(min_dist, cutoff)
      nbond = count_bond_centers_with_cutoff(cutoff)
    end function count_bond_center_candidates

    subroutine build_bond_center_projection_map(nproj, center_bohr)
      implicit none
      integer, intent(in) :: nproj
      real(8), allocatable, intent(out) :: center_bohr(:,:)
      integer :: ia, ja, axis, ip, ibond, nbond, shell
      real(8) :: cutoff, min_dist, dist2, delta(3), center(3), length_axis(3)
      real(8) :: a1(3), a2(3), a3(3)

      if(nproj <= 0) stop "DC-LCFO Wannier export: bond_centers requires positive num_wann."
      call get_lattice_vectors(a1, a2, a3)
      length_axis(1) = cell_length(a1)
      length_axis(2) = cell_length(a2)
      length_axis(3) = cell_length(a3)
      call find_bond_center_cutoff(min_dist, cutoff)
      nbond = count_bond_centers_with_cutoff(cutoff)
      if(nbond <= 0) stop "DC-LCFO Wannier export: no bond-center projection candidates were generated."

      allocate(center_bohr(3,nproj))
      ip = 0
      shell = 0
      do while(ip < nproj)
        ibond = 0
        do ia=1,dc%system_tot%nion-1
          do ja=ia+1,dc%system_tot%nion
            dist2 = 0d0
            do axis=1,3
              delta(axis) = periodic_delta(dc%system_tot%rion(axis,ja) - dc%system_tot%rion(axis,ia), &
                length_axis(axis))
              dist2 = dist2 + delta(axis) * delta(axis)
            end do
            if(sqrt(dist2) > cutoff) cycle
            ibond = ibond + 1
            ip = ip + 1
            do axis=1,3
              center(axis) = dc%system_tot%rion(axis,ia) + 0.5d0 * delta(axis)
              if(length_axis(axis) > 0d0) center(axis) = center(axis) - floor(center(axis) / length_axis(axis)) &
                * length_axis(axis)
              center_bohr(axis,ip) = center(axis)
            end do
            if(ip >= nproj) exit
          end do
          if(ip >= nproj) exit
        end do
        shell = shell + 1
        if(ibond <= 0) exit
      end do
      if(ip < nproj) stop "DC-LCFO Wannier export: failed to complete bond-center projection map."
      if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0,a,es12.5,a,es12.5)') &
        "[DC-LCFO-WANNIER] bond-center candidates=", nbond, " target_wann=", nproj, &
        " nearest=", min_dist, " cutoff=", cutoff
    end subroutine build_bond_center_projection_map

    subroutine find_bond_center_cutoff(min_dist, cutoff)
      implicit none
      real(8), intent(out) :: min_dist, cutoff
      integer :: ia, ja, axis
      real(8) :: dist2, dist, delta_axis, length_axis(3)
      real(8) :: a1(3), a2(3), a3(3)

      call get_lattice_vectors(a1, a2, a3)
      length_axis(1) = cell_length(a1)
      length_axis(2) = cell_length(a2)
      length_axis(3) = cell_length(a3)
      min_dist = huge(1d0)
      do ia=1,dc%system_tot%nion-1
        do ja=ia+1,dc%system_tot%nion
          dist2 = 0d0
          do axis=1,3
            delta_axis = periodic_delta(dc%system_tot%rion(axis,ja) - dc%system_tot%rion(axis,ia), &
              length_axis(axis))
            dist2 = dist2 + delta_axis * delta_axis
          end do
          dist = sqrt(dist2)
          if(dist > 1d-8) min_dist = min(min_dist, dist)
        end do
      end do
      if(min_dist >= huge(1d0) * 0.5d0) &
        stop "DC-LCFO Wannier export: failed to find nearest-neighbor bond distance."
      cutoff = 1.20d0 * min_dist
    end subroutine find_bond_center_cutoff

    integer function count_bond_centers_with_cutoff(cutoff) result(nbond)
      implicit none
      real(8), intent(in) :: cutoff
      integer :: ia, ja, axis
      real(8) :: dist2, delta_axis, length_axis(3)
      real(8) :: a1(3), a2(3), a3(3)

      call get_lattice_vectors(a1, a2, a3)
      length_axis(1) = cell_length(a1)
      length_axis(2) = cell_length(a2)
      length_axis(3) = cell_length(a3)
      nbond = 0
      do ia=1,dc%system_tot%nion-1
        do ja=ia+1,dc%system_tot%nion
          dist2 = 0d0
          do axis=1,3
            delta_axis = periodic_delta(dc%system_tot%rion(axis,ja) - dc%system_tot%rion(axis,ia), &
              length_axis(axis))
            dist2 = dist2 + delta_axis * delta_axis
          end do
          if(sqrt(dist2) <= cutoff) nbond = nbond + 1
        end do
      end do
    end function count_bond_centers_with_cutoff

    integer function count_pseudo_channel_ao_candidates() result(nproj)
      use salmon_global, only: izatom, wannier_site_symmetry, wannier_num_wann
      implicit none
      integer :: ia, iz, lmax_ao, projection_lmax
      logical :: ok
      character(256) :: message

      if(trim(wannier_site_symmetry) /= 'off') then
        call sawf_projection_shell_lmax(dc%system_tot%nion, wannier_num_wann, projection_lmax, ok, message)
        if(ok) call sawf_spd_projection_count(dc%system_tot%nion, nproj, ok, message, projection_lmax)
        if(.not. ok) call lcfo_sawf_fatal(message)
        return
      end if
      nproj = 0
      do ia=1,dc%system_tot%nion
        iz = izatom(dc%system_tot%kion(ia))
        lmax_ao = pseudo_channel_ao_lmax_for_species(iz)
        nproj = nproj + count_real_ao_for_lmax(lmax_ao)
      end do
    end function count_pseudo_channel_ao_candidates

    integer function pseudo_channel_ao_lmax_for_species(iz) result(lmax_ao)
      implicit none
      integer, intent(in) :: iz
      integer :: lmax_ps

      ! Initial conservative heuristic until pseudopotential-channel metadata is wired in:
      ! H gets s/p, while heavier elements get s/p/d polarization candidates.
      if(iz == 1) then
        lmax_ps = 0
      else
        lmax_ps = 1
      end if
      lmax_ao = min(lmax_ps + 1, 2)
    end function pseudo_channel_ao_lmax_for_species

    integer function count_real_ao_for_lmax(lmax_ao) result(norb)
      implicit none
      integer, intent(in) :: lmax_ao

      norb = 0
      if(lmax_ao >= 0) norb = norb + 1
      if(lmax_ao >= 1) norb = norb + 3
      if(lmax_ao >= 2) norb = norb + 5
    end function count_real_ao_for_lmax

    subroutine build_pseudo_channel_ao_candidate_map(proj_atom, proj_l, proj_m)
      use salmon_global, only: izatom, wannier_site_symmetry
      use communication, only: comm_get_max
      implicit none
      integer, intent(out) :: proj_atom(:), proj_l(:), proj_m(:)
      type(t_sawf_projection_channel), allocatable :: channels(:)
      integer :: failure_flag, ia, iz, lmax_ao, ip, im, projection_lmax
      logical :: ok
      character(256) :: message

      if(trim(wannier_site_symmetry) /= 'off') then
        call sawf_projection_shell_lmax(dc%system_tot%nion, size(proj_atom), projection_lmax, ok, message)
        if(ok) call build_sawf_spd_projection_map(dc%system_tot%nion, channels, ok, message, projection_lmax)
        failure_flag = 0
        if(.not. ok) then
          failure_flag = 1
        else if(size(proj_atom) /= size(channels) .or. size(proj_l) /= size(channels) .or. &
            size(proj_m) /= size(channels)) then
          failure_flag = 1
        end if
        call comm_get_max(failure_flag, dc%icomm_tot)
        if(failure_flag /= 0) call lcfo_sawf_fatal( &
          'SAWF pseudo-channel map allocation or size validation failed on one or more ranks')
        do ip=1,size(channels)
          proj_atom(ip) = channels(ip)%atom
          proj_l(ip) = channels(ip)%l
          proj_m(ip) = channels(ip)%m
        end do
        return
      end if

      ip = 0
      do ia=1,dc%system_tot%nion
        iz = izatom(dc%system_tot%kion(ia))
        lmax_ao = pseudo_channel_ao_lmax_for_species(iz)
        if(lmax_ao >= 0) then
          ip = ip + 1
          proj_atom(ip) = ia
          proj_l(ip) = 0
          proj_m(ip) = 1
        end if
        if(lmax_ao >= 1) then
          do im=1,3
            ip = ip + 1
            proj_atom(ip) = ia
            proj_l(ip) = 1
            proj_m(ip) = im
          end do
        end if
        if(lmax_ao >= 2) then
          do im=1,5
            ip = ip + 1
            proj_atom(ip) = ia
            proj_l(ip) = 2
            proj_m(ip) = im
          end do
        end if
      end do
    end subroutine build_pseudo_channel_ao_candidate_map

    subroutine build_c_sp3_projection_map(proj_atom, proj_hybrid)
      use salmon_global, only: izatom, wannier_projection
      implicit none
      integer, intent(out) :: proj_atom(:), proj_hybrid(:)
      integer :: ia, ih, ip, target_iz

      ip = 0
      target_iz = sp3_projection_iz(trim(wannier_projection))
      do ia=1,dc%system_tot%nion
        if(izatom(dc%system_tot%kion(ia)) /= target_iz) cycle
        do ih=1,4
          ip = ip + 1
          proj_atom(ip) = ia
          proj_hybrid(ip) = ih
        end do
      end do
    end subroutine build_c_sp3_projection_map

    real(8) function pseudo_channel_projection_value(x, y, z, ia, l, m, sigma) result(val)
      use salmon_global, only: wannier_site_symmetry
      implicit none
      real(8), intent(in) :: x, y, z, sigma
      integer, intent(in) :: ia, l, m
      real(8) :: dx, dy, dz, r2, gaussian, inv_sigma, scaled_r(3)
      real(8) :: a1(3), a2(3), a3(3)

      call get_lattice_vectors(a1, a2, a3)
      dx = periodic_delta(x - dc%system_tot%rion(1,ia), cell_length(a1))
      dy = periodic_delta(y - dc%system_tot%rion(2,ia), cell_length(a2))
      dz = periodic_delta(z - dc%system_tot%rion(3,ia), cell_length(a3))
      r2 = dx*dx + dy*dy + dz*dz
      if(r2 > (8d0*sigma)**2) then
        val = 0d0
        return
      end if

      gaussian = exp(-0.5d0 * r2 / (sigma*sigma))
      inv_sigma = 1d0 / sigma
      scaled_r = [dx,dy,dz] * inv_sigma
      if(trim(wannier_site_symmetry) /= 'off') then
        val = sawf_real_harmonic_value(l, m, scaled_r) * gaussian
        return
      end if

      ! Preserve the historical non-SAWF pseudo-channel ordering exactly.
      select case(l)
      case(0)
        val = gaussian
      case(1)
        if(m >= 1 .and. m <= 3) then
          val = scaled_r(m) * gaussian
        else
          val = 0.0d0
        end if
      case(2)
        select case(m)
        case(1)
          val = scaled_r(1)*scaled_r(2)*gaussian
        case(2)
          val = scaled_r(2)*scaled_r(3)*gaussian
        case(3)
          val = scaled_r(3)*scaled_r(1)*gaussian
        case(4)
          val = (scaled_r(1)**2-scaled_r(2)**2)*gaussian
        case(5)
          val = (2.0d0*scaled_r(3)**2-scaled_r(1)**2-scaled_r(2)**2)*gaussian
        case default
          val = 0.0d0
        end select
      case default
        val = 0.0d0
      end select
    end function pseudo_channel_projection_value

    real(8) function bond_center_projection_value(x, y, z, center_bohr, sigma) result(val)
      implicit none
      real(8), intent(in) :: x, y, z, center_bohr(3), sigma
      real(8) :: dx, dy, dz, r2
      real(8) :: a1(3), a2(3), a3(3)

      call get_lattice_vectors(a1, a2, a3)
      dx = periodic_delta(x - center_bohr(1), cell_length(a1))
      dy = periodic_delta(y - center_bohr(2), cell_length(a2))
      dz = periodic_delta(z - center_bohr(3), cell_length(a3))
      r2 = dx*dx + dy*dy + dz*dz
      if(r2 > (8d0*sigma)**2) then
        val = 0d0
        return
      end if
      val = exp(-0.5d0 * r2 / (sigma*sigma))
    end function bond_center_projection_value

    real(8) function pseudo_channel_projection_value_local_periodic(x, y, z, ia, l, m, sigma, box_length) result(val)
      implicit none
      real(8), intent(in) :: x, y, z, sigma, box_length(3)
      integer, intent(in) :: ia, l, m
      real(8) :: dx, dy, dz, r2, gaussian, inv_sigma, inv_sigma2

      dx = periodic_delta(x - dc%system_tot%rion(1,ia), box_length(1))
      dy = periodic_delta(y - dc%system_tot%rion(2,ia), box_length(2))
      dz = periodic_delta(z - dc%system_tot%rion(3,ia), box_length(3))
      r2 = dx*dx + dy*dy + dz*dz
      if(r2 > (8d0*sigma)**2) then
        val = 0d0
        return
      end if

      gaussian = exp(-0.5d0 * r2 / (sigma*sigma))
      inv_sigma = 1d0 / sigma
      inv_sigma2 = inv_sigma * inv_sigma
      select case(l)
      case(0)
        val = gaussian
      case(1)
        select case(m)
        case(1)
          val = dx * inv_sigma * gaussian
        case(2)
          val = dy * inv_sigma * gaussian
        case(3)
          val = dz * inv_sigma * gaussian
        case default
          val = 0d0
        end select
      case(2)
        select case(m)
        case(1)
          val = dx * dy * inv_sigma2 * gaussian
        case(2)
          val = dy * dz * inv_sigma2 * gaussian
        case(3)
          val = dz * dx * inv_sigma2 * gaussian
        case(4)
          val = (dx*dx - dy*dy) * inv_sigma2 * gaussian
        case(5)
          val = (2d0*dz*dz - dx*dx - dy*dy) * inv_sigma2 * gaussian
        case default
          val = 0d0
        end select
      case default
        val = 0d0
      end select
    end function pseudo_channel_projection_value_local_periodic

    real(8) function bond_center_projection_value_local_periodic(x, y, z, center_bohr, sigma, box_length) result(val)
      implicit none
      real(8), intent(in) :: x, y, z, center_bohr(3), sigma, box_length(3)
      real(8) :: dx, dy, dz, r2

      dx = periodic_delta(x - center_bohr(1), box_length(1))
      dy = periodic_delta(y - center_bohr(2), box_length(2))
      dz = periodic_delta(z - center_bohr(3), box_length(3))
      r2 = dx*dx + dy*dy + dz*dz
      if(r2 > (8d0*sigma)**2) then
        val = 0d0
        return
      end if
      val = exp(-0.5d0 * r2 / (sigma*sigma))
    end function bond_center_projection_value_local_periodic

    real(8) function c_sp3_projection_value(x, y, z, ia, ih, sigma) result(val)
      implicit none
      real(8), intent(in) :: x, y, z, sigma
      integer, intent(in) :: ia, ih
      real(8) :: dx, dy, dz, r2, gaussian, s_part, px_part, py_part, pz_part
      real(8) :: ns, np, pi_const, sx, sy, sz
      real(8) :: a1(3), a2(3), a3(3)

      call c_sp3_signs(ih, sx, sy, sz)
      call get_lattice_vectors(a1, a2, a3)
      dx = periodic_delta(x - dc%system_tot%rion(1,ia), cell_length(a1))
      dy = periodic_delta(y - dc%system_tot%rion(2,ia), cell_length(a2))
      dz = periodic_delta(z - dc%system_tot%rion(3,ia), cell_length(a3))
      r2 = dx*dx + dy*dy + dz*dz
      if(r2 > (8d0*sigma)**2) then
        val = 0d0
        return
      end if

      pi_const = acos(-1d0)
      gaussian = exp(-0.5d0 * r2 / (sigma*sigma))
      ns = 1d0 / (pi_const**0.75d0 * sigma**1.5d0)
      np = sqrt(2d0) / (pi_const**0.75d0 * sigma**2.5d0)
      s_part = ns * gaussian
      px_part = np * dx * gaussian
      py_part = np * dy * gaussian
      pz_part = np * dz * gaussian
      val = 0.5d0 * (s_part + sx*px_part + sy*py_part + sz*pz_part)
    end function c_sp3_projection_value

    real(8) function c_sp3_projection_value_local_periodic(x, y, z, ia, ih, sigma, box_length) result(val)
      implicit none
      real(8), intent(in) :: x, y, z, sigma, box_length(3)
      integer, intent(in) :: ia, ih
      real(8) :: dx, dy, dz, r2, gaussian, s_part, px_part, py_part, pz_part
      real(8) :: ns, np, pi_const, sx, sy, sz

      call c_sp3_signs(ih, sx, sy, sz)
      dx = periodic_delta(x - dc%system_tot%rion(1,ia), box_length(1))
      dy = periodic_delta(y - dc%system_tot%rion(2,ia), box_length(2))
      dz = periodic_delta(z - dc%system_tot%rion(3,ia), box_length(3))
      r2 = dx*dx + dy*dy + dz*dz
      if(r2 > (8d0*sigma)**2) then
        val = 0d0
        return
      end if

      pi_const = acos(-1d0)
      gaussian = exp(-0.5d0 * r2 / (sigma*sigma))
      ns = 1d0 / (pi_const**0.75d0 * sigma**1.5d0)
      np = sqrt(2d0) / (pi_const**0.75d0 * sigma**2.5d0)
      s_part = ns * gaussian
      px_part = np * dx * gaussian
      py_part = np * dy * gaussian
      pz_part = np * dz * gaussian
      val = 0.5d0 * (s_part + sx*px_part + sy*py_part + sz*pz_part)
    end function c_sp3_projection_value_local_periodic

    subroutine get_lattice_vectors(a1, a2, a3)
      implicit none
      real(8), intent(out) :: a1(3), a2(3), a3(3)

      a1(1:3) = dc%system_tot%primitive_a(1:3,1)
      a2(1:3) = dc%system_tot%primitive_a(1:3,2)
      a3(1:3) = dc%system_tot%primitive_a(1:3,3)
    end subroutine get_lattice_vectors

    subroutine reciprocal_vector_from_index(idx, a1, a2, a3, gx, gy, gz)
      implicit none
      integer, intent(in) :: idx(3)
      real(8), intent(in) :: a1(3), a2(3), a3(3)
      real(8), intent(out) :: gx, gy, gz
      real(8) :: b1(3), b2(3), b3(3), volume, twopi

      twopi = 2d0 * acos(-1d0)
      volume = dot_product(a1, cross_product3(a2, a3))
      if(abs(volume) <= 1d-300) stop "DC-LCFO Wannier export: invalid lattice volume."
      b1(1:3) = twopi * cross_product3(a2, a3) / volume
      b2(1:3) = twopi * cross_product3(a3, a1) / volume
      b3(1:3) = twopi * cross_product3(a1, a2) / volume
      gx = dble(idx(1)) * b1(1) + dble(idx(2)) * b2(1) + dble(idx(3)) * b3(1)
      gy = dble(idx(1)) * b1(2) + dble(idx(2)) * b2(2) + dble(idx(3)) * b3(2)
      gz = dble(idx(1)) * b1(3) + dble(idx(2)) * b2(3) + dble(idx(3)) * b3(3)
    end subroutine reciprocal_vector_from_index

    function cross_product3(a, b) result(c)
      implicit none
      real(8), intent(in) :: a(3), b(3)
      real(8) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product3

    subroutine c_sp3_signs(ih, sx, sy, sz)
      implicit none
      integer, intent(in) :: ih
      real(8), intent(out) :: sx, sy, sz

      select case(ih)
      case(1)
        sx = 1d0; sy = 1d0; sz = 1d0
      case(2)
        sx = 1d0; sy = -1d0; sz = -1d0
      case(3)
        sx = -1d0; sy = 1d0; sz = -1d0
      case default
        sx = -1d0; sy = -1d0; sz = 1d0
      end select
    end subroutine c_sp3_signs

    real(8) function periodic_delta(delta, length) result(dout)
      implicit none
      real(8), intent(in) :: delta, length

      if(length > 0d0) then
        dout = delta - dnint(delta / length) * length
      else
        dout = delta
      end if
    end function periodic_delta

    real(8) function cell_length(vec) result(length)
      implicit none
      real(8), intent(in) :: vec(3)

      length = sqrt(sum(vec(1:3)**2))
    end function cell_length

    function int_to_string(value) result(text)
      implicit none
      integer, intent(in) :: value
      character(16) :: text
      write(text,'(i0)') value
    end function int_to_string

    function element_symbol(z) result(sym)
      implicit none
      integer, intent(in) :: z
      character(2) :: sym

      select case(z)
      case(1);  sym = 'H '
      case(2);  sym = 'He'
      case(3);  sym = 'Li'
      case(4);  sym = 'Be'
      case(5);  sym = 'B '
      case(6);  sym = 'C '
      case(7);  sym = 'N '
      case(8);  sym = 'O '
      case(9);  sym = 'F '
      case(10); sym = 'Ne'
      case(11); sym = 'Na'
      case(12); sym = 'Mg'
      case(13); sym = 'Al'
      case(14); sym = 'Si'
      case(15); sym = 'P '
      case(16); sym = 'S '
      case(17); sym = 'Cl'
      case(18); sym = 'Ar'
      case default
        sym = 'X '
      end select
    end function element_symbol

  end subroutine dc_lcfo_flux

  integer function dc_buffer_box_to_local_index(ibox, ndom, nbuf) result(iloc)
    implicit none
    integer, intent(in) :: ibox, ndom, nbuf

    if (nbuf <= 0) then
      iloc = ibox
    else if (ibox <= nbuf) then
      iloc = ndom + nbuf + ibox
    else
      iloc = ibox - nbuf
    end if
  end function dc_buffer_box_to_local_index

end module lcfo_flux
