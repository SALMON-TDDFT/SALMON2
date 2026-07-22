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
  use lcfo_wannier_sawf_templates, only: t_sawf_template_checkpoint,t_sawf_template_fingerprint, &
    write_sawf_template_checkpoint,read_sawf_template_checkpoint,t_sawf_ragged_local_basis, &
    materialize_sawf_ragged_local_basis,write_sawf_materialized_basis_checkpoint, &
    read_sawf_materialized_basis_checkpoint,build_sawf_shared_buffer_point_maps, &
    stitch_sawf_materialized_neighbor_pair,build_sawf_fragment_gauge_tree
  use lcfo_wannier_sawf_templates, only: validate_sawf_materialized_neighbor_closure, &
    sawf_fragments_share_face
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


  subroutine dc_lcfo_flux(lg,mg,system,info,stencil,xc_func,pp,ppn,ppg,energy,rho_s,v_local,&
      spsi,shpsi,sttpsi,srg,dc)
    use communication, only: comm_summation
    use mpi, only: MPI_Comm_split,MPI_Comm_free,MPI_Allreduce,MPI_UNDEFINED,MPI_COMM_NULL,MPI_SUCCESS,&
      MPI_INTEGER,MPI_INTEGER8,MPI_MAX,MPI_MIN,MPI_Abort,MPI_Bcast,MPI_Reduce,MPI_Gather,&
      MPI_DOUBLE_COMPLEX,MPI_DOUBLE_PRECISION,MPI_LOGICAL,MPI_SUM
    use salmon_global, only: yn_dc_lcfo_diag,yn_dg_wpw_production,yn_dg_wpw_fixed_h_relaxation,&
      n_plane_waves_dg,k_cutoff_plane_wave,&
      dg_wpw_window_buffer,dg_wpw_window_width,dg_wpw_extra_states,dg_wpw_scf_max_iter,&
      dg_wpw_gap_threshold,dg_wpw_metric_cutoff,dg_wpw_scf_mix,dg_wpw_scf_residual_tolerance,num_fragment
    use dg_wpw_g_modes, only: select_dg_wpw_g_modes
    use dg_wpw_fragment_support, only: s_dg_wpw_support_record,build_dg_wpw_fragment_support,&
      dg_wpw_support_overlap_box
    use dg_wpw_canonical_face_schedule,only:s_dg_wpw_canonical_face_record,&
      build_dg_wpw_canonical_face_schedule
    use dg_wpw_w_row_layout,only:build_dg_wpw_w_row_layout_from_owned_ids
    use dg_wpw_production_context, only: s_dg_wpw_production_context,&
      s_dg_wpw_production_context_snapshot,initialize_dg_wpw_fragment_root_context,&
      snapshot_dg_wpw_production_context,validate_dg_wpw_production_context_snapshot,&
      release_dg_wpw_production_context_snapshot
    use rt_dg_wpw_volume_halo_provider,only:s_dg_wpw_volume_halo_send,s_dg_wpw_volume_halo_state,&
      pack_dg_wpw_volume_halo_send,exchange_dg_wpw_volume_halo_schedule
    use rt_dg_wpw_volume_halo_provider,only:broadcast_dg_wpw_volume_halos
    use rt_dg_wpw_face_side_halo,only:s_dg_wpw_face_side_send,s_dg_wpw_face_side_state,&
      reduce_dg_wpw_face_side_parts,exchange_dg_wpw_face_side_schedule
    use rt_dg_wpw_trace_halo_provider,only:s_dg_wpw_trace_halo_set
    use rt_dg_wpw_face_trace_provider,only:s_wpw_face_trace_provider
    use rt_dg_wpw_face_trace_binding,only:prepare_dg_wpw_trace_halo=>bind_dg_wpw_canonical_face_sides
    use dg_wpw_core_w_provider,only:evaluate_dg_wpw_core_w_support,reconstruct_dg_wpw_core_w_support
    use dg_wpw_core_point_builder,only:evaluate_dg_wpw_core_point
    use dg_wpw_rank_local_quadrature,only:s_dg_wpw_volume_accumulator,s_dg_wpw_core_p_accumulator,&
      initialize_dg_wpw_volume_accumulator,add_dg_wpw_core_point,build_dg_wpw_rank_local_quadrature,&
      accumulate_dg_wpw_core_volume,initialize_dg_wpw_core_p_accumulator,add_dg_wpw_core_p_point,&
      s_dg_wpw_metric_accumulator,initialize_dg_wpw_metric_accumulator,add_dg_wpw_metric_point,&
      build_dg_wpw_metric_gram,build_dg_wpw_hamiltonian_volume,build_dg_wpw_ww_volume_components
    use rt_dg_wpw_candidate_halo,only:s_dg_wpw_staged_candidate,s_dg_wpw_owned_candidates,&
      route_dg_wpw_candidate_halo,wpw_candidate_kind_ww,wpw_candidate_kind_wp,wpw_candidate_kind_pp
    use dg_wpw_lcfo_operator_adapter,only:s_dg_wpw_lcfo_ww_components,&
      import_dg_wpw_lcfo_ww_components,publish_dg_wpw_lcfo_ww_components,&
      validate_dg_wpw_surface_self_policy
    use dg_wpw_bounded_operator,only:release_dg_wpw_bounded_operator,fetch_dg_wpw_support_coefficients,&
      apply_h_dg_wpw_bounded,apply_s_dg_wpw_bounded,global_gram_dg_wpw_bounded,&
      get_dg_wpw_owned_metric_diagonal,reduce_dg_wpw_metric_rhs_partials,set_dg_wpw_interface_lambda,&
      initialize_dg_wpw_fragment_block_preconditioner,apply_dg_wpw_fragment_block_preconditioner
    use dg_wpw_nonlocal_projector,only:s_dg_wpw_projector_overlap,&
      build_dg_wpw_projector_overlap_partials,reduce_dg_wpw_p_projector_partials,&
      exchange_dg_wpw_projector_overlaps,assemble_dg_wpw_nonlocal_blocks,&
      reduce_dg_wpw_projector_records,validate_dg_wpw_projector_provenance,&
      canonicalize_dg_wpw_projector_records
    use dg_wpw_projector_nonlocal,only:validate_dg_wpw_projector_channels_collective
    use dg_wpw_scf_driver,only:s_dg_wpw_scf_iteration_state,initialize_dg_wpw_scf_iteration_state,&
      advance_dg_wpw_scf_iteration,verify_dg_wpw_scf_fixed_point
    use dg_wpw_matrix_free_scf,only:run_dg_wpw_matrix_free_algebra_step,&
      initialize_dg_wpw_deterministic_seed,initialize_dg_wpw_projected_occupied
    use dg_wpw_matrix_free_scf,only:complete_dg_wpw_projected_subspace,&
      solve_dg_wpw_metric_projection
    use dg_wpw_checkpoint,only:s_dg_wpw_checkpoint_state,write_dg_wpw_checkpoint,&
      read_dg_wpw_checkpoint,write_dg_wpw_checkpoint_manifest,dg_wpw_checkpoint_checksum
    use dg_wpw_lcfo_potential_map,only:s_dg_wpw_grid_route,build_dg_wpw_core_density,&
      prepare_dg_wpw_core_density_route,return_dg_wpw_core_potential,&
      evaluate_and_run_dg_wpw_lcfo_potential_map
    use dg_wpw_lda_hartree,only:update_wpw_owned_lda_hartree,hartree_energy_global
    use lcfo_wannier_sawf_seed,only:build_sawf_projected_wannier_from_overlap,&
      apply_sawf_projected_wannier_transform,build_sawf_projected_buffer_gradients,&
      apply_sawf_projected_wannier_gradient_transform,&
      build_sawf_wannier_density,transform_sawf_wannier_occupation,&
      qualify_sawf_wannier_density_projection,canonicalize_sawf_wannier_center,&
      canonicalize_sawf_bond_identity,assemble_sawf_diagonal_periodic_links,&
      diagnose_sawf_discrete_wannier_spread,normalize_sawf_projected_wannier_columns,&
      validate_sawf_projected_wannier_columns
    use dg_wpw_wannier_tail_halo,only:exchange_sawf_wannier_tail_values,&
      exchange_sawf_discovered_wannier_tails,locate_sawf_wannier_tail_core,&
      locate_sawf_wannier_tail_rank,classify_sawf_wannier_buffer_tail,is_sawf_outer_buffer_shell
    use dg_wpw_occupied_w_basis,only:s_dg_wpw_occupied_w_basis,t_dg_wpw_periodic_image_mismatch,&
      gather_dg_wpw_occupied_w_payload,initialize_dg_wpw_occupied_w_basis_collective,&
      broadcast_dg_wpw_occupied_w_basis,evaluate_dg_wpw_occupied_w_point,&
      dg_wpw_unwrapped_to_storage_index,reorder_dg_wpw_fragment_buffer,&
      periodize_dg_wpw_fragment_buffer,assemble_dg_wpw_canonical_buffer_norm
    use dg_wpw_occupied_w_basis,only:extract_dg_wpw_canonical_cell
    use rt_dg_wpw_production_builder,only:route_dg_wpw_staged_candidates,&
      scan_dg_wpw_canonical_faces=>scan_and_route_dg_wpw_canonical_faces,&
      build_dg_wpw_production_operator,bind_dg_wpw_hs_callbacks,replace_dg_wpw_potential_volume,&
      route_and_replace_dg_wpw_potential_volume,install_dg_wpw_projector_nonlocal,&
      publish_dg_wpw_realspace_metric
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),      intent(in) :: stencil
    type(s_xc_functional),intent(in) :: xc_func
    type(s_pp_info),      intent(in) :: pp
    type(s_pp_nlcc),      intent(in) :: ppn
    type(s_pp_grid),      intent(in) :: ppg
    type(s_dft_energy),   intent(in) :: energy
    type(s_scalar),       intent(in) :: rho_s(system%nspin)
    type(s_scalar),       intent(in) :: V_local(system%nspin)
    type(s_orbital),      intent(inout) :: spsi
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
    integer,allocatable :: wpw_g_indices(:,:)
    integer,allocatable :: wpw_owned_w_ids(:),wpw_support_w_ids(:),wpw_support_w_owner(:)
    integer,allocatable :: wpw_support_fragment_ids(:),wpw_support_core_lo(:,:),wpw_support_core_hi(:,:)
    integer,allocatable :: wpw_support_overlap_lo(:,:),wpw_support_overlap_hi(:,:)
    type(s_dg_wpw_support_record),allocatable::wpw_support_records(:)
    type(s_dg_wpw_canonical_face_record),allocatable::wpw_canonical_faces(:)
    type(s_dg_wpw_volume_halo_send),allocatable::wpw_volume_send(:)
    type(s_dg_wpw_volume_halo_state),allocatable::wpw_volume_halos(:)
    type(s_dg_wpw_face_side_send),allocatable::wpw_face_sides(:)
    type(s_dg_wpw_face_side_state),allocatable::wpw_remote_face_sides(:)
    type(s_dg_wpw_trace_halo_set),target::wpw_trace_set
    type(s_wpw_face_trace_provider)::wpw_trace_provider
    real(8),allocatable :: wpw_g_vectors(:,:)
    complex(8),allocatable::wpw_volume_wp_h(:,:),wpw_volume_wp_s(:,:),wpw_volume_pp_h(:,:),wpw_volume_pp_s(:,:)
    complex(8),allocatable::wpw_metric_ww_s(:,:),wpw_metric_wp_s(:,:),wpw_metric_pp_s(:,:)
    complex(8),allocatable::wpw_metric_ww_h(:,:),wpw_metric_wp_h(:,:),wpw_metric_pp_h(:,:)
    complex(8),allocatable::wpw_metric_ww_kinetic(:,:),wpw_metric_ww_potential(:,:)
    integer :: wpw_mode_info
    integer :: wpw_production_comm,wpw_mpi_info,wpw_color
    integer :: wpw_failure_local,wpw_failure_global
    real(8) :: wpw_box_length(3)
    integer::wpw_window_buffer(3),wpw_window_width(3)
    type(s_dg_wpw_production_context) :: wpw_context
    type(s_dg_wpw_volume_accumulator)::wpw_volume_accumulator
    type(s_dg_wpw_core_p_accumulator)::wpw_core_p_accumulator
    type(s_dg_wpw_metric_accumulator)::wpw_metric_accumulator
    type(s_dg_wpw_occupied_w_basis)::wpw_occupied_w_basis
    type(s_dg_wpw_lcfo_ww_components)::wpw_ww_components
    type(s_dg_wpw_lcfo_ww_components)::wpw_frozen_ww_components
    type(s_dg_wpw_production_context_snapshot)::wpw_frozen_production_context
    type(s_dg_wpw_scf_iteration_state)::wpw_scf_state
    type(s_dg_wpw_grid_route)::wpw_density_route
    complex(8),allocatable::wpw_qw(:,:),wpw_qp(:,:),wpw_q_old_occ(:,:)
    complex(8),allocatable::wpw_bootstrap_source_values(:,:)
    integer,allocatable::wpw_bootstrap_source_keys(:,:)
    real(8),allocatable::wpw_eigenvalues(:),wpw_occupations(:)
    real(8)::wpw_bootstrap_source_condition
    integer::wpw_bootstrap_source_count,wpw_bootstrap_source_info
    complex(8),allocatable::wpw_saved_wp_volume(:),wpw_saved_pp_volume(:),wpw_saved_ww_potential(:)
    complex(8),allocatable::wpw_frozen_wp_volume(:),wpw_frozen_wp_nonlocal(:),wpw_frozen_wp_face(:),&
      wpw_frozen_pp_volume(:),wpw_frozen_pp_nonlocal(:)
    complex(8),allocatable::wpw_frozen_ww_projector_nonlocal(:,:),&
      wpw_frozen_ww_projector_cross_value(:)
    integer,allocatable::wpw_frozen_ww_projector_cross_row_id(:),&
      wpw_frozen_ww_projector_cross_col_id(:)
    complex(8),allocatable::wpw_frozen_ww_h0_dense(:,:),wpw_frozen_ww_interface_dense(:,:),&
      wpw_frozen_wp_h0_dense(:,:),wpw_frozen_wp_interface_dense(:,:),&
      wpw_frozen_pp_h0_dense(:,:),wpw_frozen_pp_interface_dense(:,:)
    integer,allocatable::wpw_frozen_owned_w_ids(:),wpw_frozen_required_w_ids(:),&
      wpw_frozen_owned_p_ids(:),wpw_frozen_required_p_ids(:)
    real(8),allocatable::wpw_saved_accumulator_potential(:),wpw_saved_accumulator_density(:)
    real(8),allocatable::wpw_saved_rho_tot(:,:,:),wpw_saved_rho_tot_s(:,:,:),&
      wpw_saved_vh_tot(:,:,:),wpw_saved_vxc_tot(:,:,:)
    real(8),allocatable::wpw_frozen_rho_tot(:,:,:),wpw_frozen_vh_tot(:,:,:),&
      wpw_frozen_vxc_tot(:,:,:),wpw_frozen_vloc_tot(:,:,:)
    integer,allocatable::wpw_total_owned_grid_ids(:),wpw_core_destinations(:)
    integer::wpw_nocc,wpw_nretain
    integer::wpw_frozen_halo_epoch,wpw_frozen_scan_epoch,wpw_frozen_operator_epoch
    logical::wpw_frozen_callbacks_bound,wpw_frozen_operator_valid,&
      wpw_frozen_ww_projector_nonlocal_valid
    integer::wpw_projection_rank
    real(8)::wpw_projection_residual,wpw_projection_captured_norm,&
      wpw_projection_orthogonality,wpw_projection_charge
    real(8)::wpw_fixed_h_final_residual,wpw_fixed_h_final_orthogonality,&
      wpw_fixed_h_final_projector,wpw_fixed_h_density_baseline,wpw_fixed_h_charge_error
    logical::wpw_fixed_h_qualified,wpw_metric_diagnostic_only
    logical :: sawf_explicit_basis_active
    logical::wpw_nonowned_candidates_pending
    logical::wpw_fixed_point_mode
    !
    integer :: i,j,n,ix,iy,iz,io,jo,ispin,ifrag,jfrag,i_halo

    if(dc%id_tot==0) write(*,*) "start DC-LCFO-Flux"
    hvol = system%hvol
    nspin = system%nspin
    sawf_explicit_basis_active=.false.
    wpw_nonowned_candidates_pending=.false.
    wpw_fixed_point_mode=.false.
    wpw_metric_diagnostic_only=.false.
    wpw_window_buffer=0
    wpw_window_width=merge(dg_wpw_window_width,0,num_fragment>1)
    wpw_production_comm=MPI_COMM_NULL
    call init_lcfo
    wpw_window_buffer=merge(dc%nxyz_buffer-stencil_radius,0,num_fragment>1)
    call initialize_sawf_seed_provenance()
    call calc_basis
    if(yn_dg_wpw_production=='y') then
      wpw_box_length=system%hgs*dble(dc%lg_tot%num)
      call select_dg_wpw_g_modes(wpw_box_length,k_cutoff_plane_wave,n_plane_waves_dg,&
        wpw_g_indices,wpw_g_vectors,wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0.and.size(wpw_g_indices,2)>0,&
        'DG WPW production G-mode selection failed')
      call build_wpw_support_geometry(wpw_support_records,wpw_support_fragment_ids,wpw_support_core_lo,&
        wpw_support_core_hi,wpw_support_overlap_lo,wpw_support_overlap_hi,wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW bounded support geometry initialization failed')
      call build_wpw_canonical_face_geometry(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW canonical-face schedule initialization failed')
      wpw_color=merge(1,MPI_UNDEFINED,dc%id_frag==0)
      call MPI_Comm_split(dc%icomm_tot,wpw_color,dc%i_frag,wpw_production_comm,wpw_mpi_info)
      call wpw_collective_require(wpw_mpi_info==MPI_SUCCESS,'DG WPW production communicator split failed collectively')
      call assemble_wpw_core_p_bootstrap(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW W-independent core/P bootstrap failed collectively')
      call build_core_owned_projected_wannier_density_seed(wpw_bootstrap_source_values,&
        wpw_bootstrap_source_count,wpw_bootstrap_source_condition,wpw_occupied_w_basis,&
        wpw_bootstrap_source_info)
      call wpw_collective_require(wpw_bootstrap_source_info==0,&
        'DG WPW occupied-W descriptor bootstrap failed collectively')
      call build_wpw_w_row_layout(wpw_owned_w_ids,wpw_support_w_ids,wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW W-row layout initialization failed')
      call broadcast_wpw_w_row_layout(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW fragment W-row layout broadcast failed collectively')
      wpw_mode_info=0
      if(dc%id_frag==0)then
        call initialize_dg_wpw_fragment_root_context(wpw_context,wpw_production_comm,dc%n_frag,&
          size(wpw_g_indices,2),wpw_occupied_w_basis%global_count,dc%i_frag,&
          wpw_support_fragment_ids,wpw_owned_w_ids,&
          wpw_support_w_ids,wpw_mode_info)
      endif
      call wpw_collective_require(wpw_mode_info==0,'DG WPW fragment-root context initialization failed collectively')
      wpw_mode_info=0
      if(dc%id_frag==0)call prepare_wpw_volume_halo(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW support-W halo preparation failed collectively')
      if(.not.allocated(wpw_volume_halos))allocate(wpw_volume_halos(0))
      call broadcast_dg_wpw_volume_halos(dc%icomm_frag,0,wpw_volume_halos,wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW fragment-rank support-W halo broadcast failed collectively')
      call assemble_wpw_core_volume_quadrature(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW core volume quadrature failed collectively')
      wpw_mode_info=0
      if(dc%id_frag==0)call publish_wpw_core_volume_candidates(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW core volume candidate publication failed collectively')
      call prepare_wpw_canonical_face_trace_provider(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW canonical-face trace preparation failed collectively')
      wpw_mode_info=0
      if(dc%id_frag==0)call scan_wpw_canonical_faces(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW canonical-face scan and routing failed collectively')
    endif
    call check_lcfo_basis_potential_inputs_finite
    call hpsi_basis
    if(dc%id_tot==0) write(*,*) "basis functions operation: done"

    call calc_hamiltonian_matrix
    if(dc%id_tot==0) write(*,*) "Hamiltonian matrix: done"
    if(yn_dg_wpw_production=='y')then
      call assemble_wpw_projector_nonlocal(wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW projector nonlocal assembly failed collectively')
      wpw_mode_info=0
      if(dc%id_frag==0)call build_dg_wpw_production_operator(wpw_context,wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW sparse production operator build failed collectively')
      wpw_mode_info=0
      if(dc%id_frag==0)call bind_dg_wpw_hs_callbacks(wpw_context,wpw_mode_info)
      call wpw_collective_require(wpw_mode_info==0,'DG WPW H/S callback binding failed collectively')
      if(yn_dg_wpw_fixed_h_relaxation=='y')then
        call snapshot_wpw_frozen_h_state(wpw_mode_info)
        call wpw_collective_require(wpw_mode_info==0,'DG WPW frozen-H snapshot failed collectively')
        call run_wpw_fixed_h_relaxation(wpw_mode_info)
        call wpw_collective_require(wpw_mode_info==0,'DG WPW fixed-H relaxation failed collectively')
      else
        call run_wpw_production_scf(wpw_mode_info)
        call wpw_collective_require(wpw_mode_info==0,'DG WPW production SCF failed collectively')
      endif
    endif

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
    if(allocated(wpw_g_indices)) deallocate(wpw_g_indices)
    if(allocated(wpw_g_vectors)) deallocate(wpw_g_vectors)
    if(allocated(wpw_owned_w_ids)) deallocate(wpw_owned_w_ids)
    if(allocated(wpw_support_w_ids)) deallocate(wpw_support_w_ids)
    if(allocated(wpw_support_w_owner)) deallocate(wpw_support_w_owner)
    if(allocated(wpw_support_fragment_ids)) deallocate(wpw_support_fragment_ids)
    if(allocated(wpw_support_core_lo)) deallocate(wpw_support_core_lo)
    if(allocated(wpw_support_core_hi)) deallocate(wpw_support_core_hi)
    if(allocated(wpw_support_overlap_lo)) deallocate(wpw_support_overlap_lo)
    if(allocated(wpw_support_overlap_hi)) deallocate(wpw_support_overlap_hi)
    if(allocated(wpw_support_records)) deallocate(wpw_support_records)
    if(allocated(wpw_canonical_faces))deallocate(wpw_canonical_faces)
    if(allocated(wpw_volume_send)) deallocate(wpw_volume_send)
    if(allocated(wpw_volume_halos)) deallocate(wpw_volume_halos)
    if(allocated(wpw_face_sides))deallocate(wpw_face_sides)
    if(allocated(wpw_remote_face_sides))deallocate(wpw_remote_face_sides)
    if(allocated(wpw_volume_wp_h))deallocate(wpw_volume_wp_h)
    if(allocated(wpw_volume_wp_s))deallocate(wpw_volume_wp_s)
    if(allocated(wpw_volume_pp_h))deallocate(wpw_volume_pp_h)
    if(allocated(wpw_volume_pp_s))deallocate(wpw_volume_pp_s)
    if(wpw_production_comm/=MPI_COMM_NULL.and.dc%id_frag==0) &
      call release_dg_wpw_bounded_operator(wpw_context%bounded_operator)
    if(wpw_production_comm/=MPI_COMM_NULL)call MPI_Comm_free(wpw_production_comm,wpw_mpi_info)
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

    subroutine assemble_wpw_projector_nonlocal(nonlocal_info)
      integer,intent(out)::nonlocal_info
      type(s_dg_wpw_projector_overlap),allocatable::local_records(:),w_records(:),p_partial(:),&
        p_owned(:),owned_records(:),support_records(:),canonical_records(:)
      integer::nw,np,nlma,nsample,ilma,ia,j,sample,i,ib,ip,nwrec,nprec,ilocal,n_w_global,&
        ierr_nl,root_info,local_info,npeer,ncross,k,production_bad,astat,ik,ll,lproj,m,l0,&
        channel_scan,angular_code,projector_ordinal
      integer,allocatable::p_ids(:),sample_grid_ids(:),sample_projector_ids(:),peers(:),projector_ids(:),&
        ww_rows(:),ww_cols(:),cross_rows(:),cross_cols(:),projector_channel_key(:,:),projector_owner(:)
      complex(8),allocatable::sample_values(:),ww_values(:),&
        cross_values(:),wp_values(:),pp_values(:)
      real(8)::scale,projector_stage_clock
      real(8),allocatable::projector_weights(:)
      integer::projector_fallback(1),atom_fallback(1)
      real(8)::weight_fallback(1)
      integer(8)::np8,nsample8,nentry8

      nonlocal_info=1;root_info=0;nsample=0
      call cpu_time(projector_stage_clock)
      nw=size(wpw_owned_w_ids);n_w_global=wpw_occupied_w_basis%global_count;np=0
      np8=int(size(wpw_support_fragment_ids),8)*int(size(wpw_g_indices,2),8)
      if(np8<=0_8.or.np8>int(huge(0),8))root_info=1
      if(root_info==0)np=int(np8)
      nlma=ppg%Nlma
      if(nw<=0.or.n_w_global<=0.or.np<=0.or.nlma<=0.or..not.allocated(ppg%uV).or.&
         .not.allocated(ppg%rinv_uvu).or..not.allocated(ppg%ia_tbl).or.&
         .not.allocated(ppg%mps).or..not.allocated(ppg%jxyz).or.&
         .not.allocated(ppg%ilocal_nlma2ilma).or.&
         .not.allocated(ppg%ilocal_nlma2ia))root_info=1
      if(root_info==0)then
      if(ppg%ilocal_nlma<0.or.size(ppg%ilocal_nlma2ilma)/=ppg%ilocal_nlma.or.&
         size(ppg%ilocal_nlma2ia)/=ppg%ilocal_nlma.or.size(ppg%ia_tbl)<nlma.or.&
         size(ppg%rinv_uvu)<nlma.or.size(ppg%uV,2)<nlma)root_info=1
      if(root_info==0)then
        do ilocal=1,ppg%ilocal_nlma
          ilma=ppg%ilocal_nlma2ilma(ilocal);ia=ppg%ilocal_nlma2ia(ilocal)
          if(ilma<1.or.ilma>nlma.or.ia<1.or.ia>size(ppg%mps))root_info=1
        enddo
      endif
      if(root_info==0)then
        if(.not.allocated(ppg%lma_tbl).or..not.allocated(system%kion).or.&
          .not.allocated(pp%mlps).or..not.allocated(pp%nproj).or..not.allocated(pp%inorm))root_info=1
      endif
      if(root_info==0)then
        if(any(ppg%ia_tbl(1:nlma)<1))root_info=1
      endif
      if(root_info==0)then
        if(size(system%kion)<maxval(ppg%ia_tbl(1:nlma)).or.&
          size(ppg%lma_tbl,2)<maxval(ppg%ia_tbl(1:nlma)))root_info=1
      endif
      if(root_info==0)then
        do ia=1,maxval(ppg%ia_tbl(1:nlma))
          if(count(ppg%ia_tbl(1:nlma)==ia)>size(ppg%lma_tbl,1))then;root_info=1;exit;endif
        enddo
      endif
      if(root_info==0)then
        if(any(system%kion(ppg%ia_tbl(1:nlma))<1).or.&
          any(system%kion(ppg%ia_tbl(1:nlma))>size(pp%mlps)))root_info=1
      endif
      if(root_info==0)then
        if(size(pp%nproj,2)<maxval(system%kion(ppg%ia_tbl(1:nlma))).or.&
          size(pp%inorm,2)<maxval(system%kion(ppg%ia_tbl(1:nlma))))root_info=1
      endif
      if(root_info==0)then
        do i=1,nlma
          ik=system%kion(ppg%ia_tbl(i))
          if(lbound(pp%nproj,1)>0.or.ubound(pp%nproj,1)<pp%mlps(ik).or.&
            lbound(pp%inorm,1)>0)then;root_info=1;exit;endif
          if(ubound(pp%inorm,1)<sum(pp%nproj(0:pp%mlps(ik),ik))-1)then;root_info=1;exit;endif
        enddo
      endif
      if(root_info==0)then
        do ilocal=1,ppg%ilocal_nlma
          ia=ppg%ilocal_nlma2ia(ilocal)
          if(ppg%mps(ia)<0.or.size(ppg%jxyz,1)<3.or.size(ppg%jxyz,2)<ppg%mps(ia).or.&
             size(ppg%jxyz,3)<ia.or.size(ppg%uV,1)<ppg%mps(ia))root_info=1
        enddo
      endif
      endif
      if(root_info==0)then
      allocate(p_ids(np),stat=astat);if(astat/=0)root_info=1
      if(root_info==0)then;do i=1,np
        nentry8=(int(wpw_support_fragment_ids((i-1)/size(wpw_g_indices,2)+1),8)-1_8)*&
          int(size(wpw_g_indices,2),8)+int(modulo(i-1,size(wpw_g_indices,2))+1,8)
        if(nentry8<=0_8.or.nentry8>int(huge(0),8))then;root_info=1;exit;endif
        p_ids(i)=int(nentry8)
      enddo;endif
      nsample8=0_8
      do ilocal=1,ppg%ilocal_nlma
        ia=ppg%ilocal_nlma2ia(ilocal);nsample8=nsample8+int(ppg%mps(ia),8)
      enddo
      if(nsample8>int(huge(0),8))root_info=1
      if(root_info==0)nsample=int(nsample8)
      endif
      if(root_info==0)then
      allocate(sample_grid_ids(nsample),sample_projector_ids(nsample),sample_values(nsample),stat=astat)
      if(astat/=0)root_info=1;sample=0
      if(root_info==0)then
      do ilocal=1,ppg%ilocal_nlma
        ilma=ppg%ilocal_nlma2ilma(ilocal);ia=ppg%ilocal_nlma2ia(ilocal)
        do j=1,ppg%mps(ia)
          sample=sample+1;sample_projector_ids(sample)=ilma
          sample_grid_ids(sample)=wpw_core_global_grid_id(dc%ixyz_frag(:,dc%i_frag)+ppg%jxyz(:,j,ia)-1)
          sample_values(sample)=cmplx(ppg%uV(j,ilma),0d0,8)
        enddo
      enddo
      endif
      endif
      allocate(local_records(0));local_info=root_info
      if(local_info==0.and.wpw_volume_accumulator%npoint>0.and.nsample>0)then
        deallocate(local_records)
        call build_dg_wpw_projector_overlap_partials(&
          wpw_volume_accumulator%grid_ids(1:wpw_volume_accumulator%npoint),&
          wpw_metric_accumulator%weights(1:wpw_metric_accumulator%npoint),wpw_support_w_ids,&
          wpw_metric_accumulator%w_points(:,1:wpw_metric_accumulator%npoint),p_ids,&
          wpw_metric_accumulator%p_points(:,1:wpw_metric_accumulator%npoint),sample_grid_ids,&
          sample_projector_ids,sample_values,n_w_global,local_records,local_info)
        if(local_info/=0)write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] projector_quadrature rank=',&
          dc%id_tot,' info=',local_info
      endif
      call reduce_dg_wpw_projector_records(dc%icomm_frag,0,local_records,support_records,&
        local_info,root_info)
      if(dc%id_frag==0)then
        if(root_info/=0)write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] projector_fragment_reduce fragment=',&
          dc%i_frag,' info=',root_info
        nwrec=0;nprec=0
        if(root_info==0)then
          nwrec=count(support_records%basis_id<=n_w_global)
          nprec=size(support_records)-nwrec
        endif
        allocate(w_records(nwrec),p_partial(nprec),stat=astat)
        if(astat/=0)root_info=1
        call synchronize_wpw_projector_root_info(root_info)
        nwrec=0;nprec=0
        if(root_info==0)then
          do i=1,size(support_records)
            if(support_records(i)%basis_id<=n_w_global)then
              nwrec=nwrec+1;w_records(nwrec)=support_records(i)
            else
              nprec=nprec+1;p_partial(nprec)=support_records(i)
            endif
          enddo
        endif
        allocate(projector_ids(max(1,nlma)),projector_weights(max(1,nlma)),&
          projector_channel_key(6,max(1,nlma)),projector_owner(max(1,nlma)),stat=astat)
        if(astat/=0)root_info=1
        call synchronize_wpw_projector_root_info(root_info)
        if(root_info==0)then
          projector_ids=0;projector_weights=0d0;projector_channel_key=0;projector_owner=0
          projector_weights=ppg%rinv_uvu(1:nlma)/hvol
          if(.not.all(ieee_is_finite(projector_weights)))root_info=1
        endif
        if(root_info==0)then
          do i=1,nlma
            ia=ppg%ia_tbl(i)
            ik=system%kion(ia);channel_scan=0;angular_code=-1;projector_ordinal=0;l0=0
            channel_lookup: do ll=0,pp%mlps(ik)
              do lproj=l0,l0+pp%nproj(ll,ik)-1
                if(pp%inorm(lproj,ik)==0)cycle
                do m=-ll,ll
                  channel_scan=channel_scan+1
                  if(ppg%lma_tbl(channel_scan,ia)==i)then
                    angular_code=ll*(ll+1)+m;projector_ordinal=lproj-l0+1
                    exit channel_lookup
                  endif
                enddo
              enddo
              l0=lproj
            enddo channel_lookup
            projector_channel_key(:,i)=[ia,0,0,0,angular_code,projector_ordinal]
            projector_owner(i)=modulo(ia-1,dc%n_frag)
            projector_ids(i)=i
          enddo
        endif
        call validate_dg_wpw_projector_channels_collective(wpw_production_comm,&
          projector_channel_key(:,1:nlma),projector_owner(1:nlma),root_info)
        if(root_info==0)then
          call validate_dg_wpw_projector_provenance(wpw_production_comm,projector_ids,&
            ppg%ia_tbl(1:nlma),ppg%rinv_uvu(1:nlma),local_info,root_info)
        else
          projector_fallback=0;atom_fallback=0;weight_fallback=0d0
          call validate_dg_wpw_projector_provenance(wpw_production_comm,projector_fallback,&
            atom_fallback,weight_fallback,local_info,root_info)
        endif
        if(local_info/=0)root_info=1
        call MPI_Allreduce(root_info,production_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_nl)
        if(ierr_nl/=MPI_SUCCESS.or.production_bad/=0)root_info=1
        if(root_info==0)call canonicalize_dg_wpw_projector_records(wpw_production_comm,&
          projector_ids(1:nlma),projector_owner(1:nlma),support_records,canonical_records,root_info)
        call synchronize_wpw_projector_root_info(root_info)
        if(root_info==0)call move_alloc(canonical_records,support_records)
        call synchronize_wpw_projector_root_info(root_info)
        call trace_wpw_projector_stage('canonical_projector_owner',projector_stage_clock,size(support_records),&
          size(wpw_context%pp_r))
        nentry8=int(nw,8)*int(nw,8)
        if(nentry8>int(huge(0),8))root_info=1
        call synchronize_wpw_projector_root_info(root_info)
        if(root_info==0)then
          allocate(ww_rows(nw*nw),ww_cols(nw*nw),ww_values(nw*nw),stat=astat)
          if(astat/=0)root_info=1
        endif
        call synchronize_wpw_projector_root_info(root_info)
        if(root_info==0)then;k=0
          do j=1,nw;do i=1,nw;k=k+1;ww_rows(k)=wpw_owned_w_ids(i);ww_cols(k)=wpw_owned_w_ids(j);enddo;enddo
          call assemble_dg_wpw_nonlocal_blocks(support_records,projector_ids,projector_weights,&
            ww_rows,ww_cols,ww_values,root_info)
          if(root_info/=0)write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] projector_ww_local fragment=',&
            dc%i_frag,' info=',root_info
        endif
        call synchronize_wpw_projector_root_info(root_info)
        if(allocated(ww_values))call trace_wpw_projector_stage('ww_local',projector_stage_clock,&
          size(support_records),size(ww_values))
        if(root_info==0)then
          ncross=0
          do j=1,size(wpw_context%support_w_ids)
            if(.not.any(wpw_owned_w_ids==wpw_context%support_w_ids(j)))then
              nentry8=int(ncross,8)+int(nw,8)
              if(nentry8>int(huge(0),8))then;root_info=1;exit;endif
              ncross=int(nentry8)
            endif
          enddo
          if(root_info/=0)then
            ncross=0
          endif
        endif
        call synchronize_wpw_projector_root_info(root_info)
        if(root_info==0)then
          allocate(cross_rows(ncross),cross_cols(ncross),cross_values(ncross),stat=astat)
          if(astat/=0)root_info=1
        endif
        call synchronize_wpw_projector_root_info(root_info)
        if(root_info==0)then;k=0
          do j=1,size(wpw_context%support_w_ids)
            if(any(wpw_owned_w_ids==wpw_context%support_w_ids(j)))cycle
            do i=1,nw;k=k+1;cross_rows(k)=wpw_owned_w_ids(i);cross_cols(k)=wpw_context%support_w_ids(j);enddo
          enddo
          call assemble_dg_wpw_nonlocal_blocks(support_records,projector_ids,projector_weights,&
            cross_rows,cross_cols,cross_values,root_info)
          if(root_info/=0)write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] projector_ww_cross fragment=',&
            dc%i_frag,' info=',root_info
        endif
        call synchronize_wpw_projector_root_info(root_info)
        if(allocated(cross_values))call trace_wpw_projector_stage('ww_cross',projector_stage_clock,&
          size(support_records),size(cross_values))
        if(root_info==0)then
          if(any(int(wpw_context%wp_p,8)>int(huge(0),8)-int(n_w_global,8)).or.&
             any(int(wpw_context%pp_r,8)>int(huge(0),8)-int(n_w_global,8)).or.&
             any(int(wpw_context%pp_c,8)>int(huge(0),8)-int(n_w_global,8)))root_info=1
        endif
        call synchronize_wpw_projector_root_info(root_info)
        if(root_info==0)then
          allocate(wp_values(size(wpw_context%wp_w)),pp_values(size(wpw_context%pp_r)),stat=astat)
          if(astat/=0)root_info=1
        endif
        call synchronize_wpw_projector_root_info(root_info)
        if(root_info==0)then
          call assemble_dg_wpw_nonlocal_blocks(support_records,projector_ids,projector_weights,&
            wpw_context%wp_w,n_w_global+wpw_context%wp_p,wp_values,root_info)
          if(root_info/=0)write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] projector_wp fragment=',&
            dc%i_frag,' info=',root_info
          if(allocated(wp_values))call trace_wpw_projector_stage('wp',projector_stage_clock,&
            size(support_records),size(wp_values))
          if(root_info==0)call assemble_dg_wpw_nonlocal_blocks(support_records,projector_ids,&
            projector_weights,n_w_global+wpw_context%pp_r,n_w_global+wpw_context%pp_c,pp_values,root_info)
          if(root_info/=0)write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] projector_pp fragment=',&
            dc%i_frag,' info=',root_info
        endif
        call synchronize_wpw_projector_root_info(root_info)
        if(allocated(pp_values))call trace_wpw_projector_stage('pp',projector_stage_clock,&
          size(support_records),size(pp_values))
        if(root_info==0)then
          call install_dg_wpw_projector_nonlocal(wpw_context,reshape(ww_values,[nw,nw]),&
            cross_rows,cross_cols,cross_values,wp_values,pp_values,root_info)
          if(root_info/=0)write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] projector_install fragment=',&
            dc%i_frag,' info=',root_info
        endif
        call trace_wpw_projector_stage('install',projector_stage_clock,size(support_records),&
          size(wpw_context%pp_r))
      endif
      call MPI_Bcast(root_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr_nl)
      if(ierr_nl/=MPI_SUCCESS.or.root_info/=0)return
      nonlocal_info=0
    end subroutine assemble_wpw_projector_nonlocal

    subroutine trace_wpw_projector_stage(stage,stage_clock,nrecords,nvalues)
      character(*),intent(in)::stage
      real(8),intent(inout)::stage_clock
      integer,intent(in)::nrecords,nvalues
      real(8)::now
      call cpu_time(now)
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,f10.3)')'[DG-WPW-PROJECTOR-STAGE] stage=',trim(stage),&
        ' fragment=',dc%i_frag,' records=',nrecords,' values=',nvalues,' cpu_seconds=',now-stage_clock
      flush(6)
      stage_clock=now
    end subroutine trace_wpw_projector_stage

    subroutine synchronize_wpw_projector_root_info(stage_info)
      integer,intent(inout)::stage_info
      integer::global_info,ierr_sync
      call MPI_Allreduce(stage_info,global_info,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_sync)
      if(ierr_sync/=MPI_SUCCESS)then;stage_info=1
      else;stage_info=global_info;endif
    end subroutine synchronize_wpw_projector_root_info

    subroutine wpw_collective_require(local_ok,message)
      logical,intent(in)::local_ok
      character(*),intent(in)::message
      integer::local_failure,global_failure,ierr_require,ierr_abort
      local_failure=merge(0,1,local_ok)
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr_require)
      if(ierr_require==MPI_SUCCESS.and.global_failure==0)return
      if(dc%id_tot==0)write(*,'(1x,2a)')'[FATAL] ',trim(message)
      call MPI_Abort(dc%icomm_tot,1,ierr_abort)
      error stop 'DG WPW communicator abort returned unexpectedly'
    end subroutine wpw_collective_require

    subroutine snapshot_wpw_frozen_h_state(snapshot_info)
      integer,intent(out)::snapshot_info
      integer::astat,local_bad

      snapshot_info=1;local_bad=0
      call release_wpw_frozen_h_state()
      allocate(wpw_frozen_rho_tot,source=dc%rho_tot%f,stat=astat);if(astat/=0)local_bad=1
      if(local_bad==0)then;allocate(wpw_frozen_vh_tot,source=dc%vh_tot%f,stat=astat);if(astat/=0)local_bad=1;endif
      if(local_bad==0)then;allocate(wpw_frozen_vxc_tot,source=dc%vxc_tot(1)%f,stat=astat);if(astat/=0)local_bad=1;endif
      if(local_bad==0)then;allocate(wpw_frozen_vloc_tot,source=dc%vloc_tot(1)%f,stat=astat);if(astat/=0)local_bad=1;endif
      if(dc%id_frag==0.and.local_bad==0)then
        if(.not.wpw_context%realspace_ww_metric_valid.or..not.allocated(wpw_context%wp_h_volume).or.&
           .not.allocated(wpw_context%wp_h_nonlocal).or..not.allocated(wpw_context%wp_h_face).or.&
           .not.allocated(wpw_context%pp_h_volume).or..not.allocated(wpw_context%pp_h_nonlocal).or.&
           .not.wpw_context%ww_projector_nonlocal_valid.or.&
           .not.allocated(wpw_context%ww_projector_nonlocal).or.&
           .not.allocated(wpw_context%ww_projector_cross_value).or.&
           .not.allocated(wpw_context%ww_projector_cross_row_id).or.&
           .not.allocated(wpw_context%ww_projector_cross_col_id).or.&
           .not.wpw_context%bounded_operator%valid)local_bad=1
      endif
      if(dc%id_frag==0.and.local_bad==0)then
        allocate(wpw_frozen_wp_volume,source=wpw_context%wp_h_volume,stat=astat);if(astat/=0)local_bad=1
        if(local_bad==0)then;allocate(wpw_frozen_wp_nonlocal,source=wpw_context%wp_h_nonlocal,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_wp_face,source=wpw_context%wp_h_face,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_pp_volume,source=wpw_context%pp_h_volume,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_pp_nonlocal,source=wpw_context%pp_h_nonlocal,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then
          allocate(wpw_frozen_ww_projector_nonlocal,&
            source=wpw_context%ww_projector_nonlocal,stat=astat);if(astat/=0)local_bad=1
        endif
        if(local_bad==0)then
          allocate(wpw_frozen_ww_projector_cross_value,&
            source=wpw_context%ww_projector_cross_value,stat=astat);if(astat/=0)local_bad=1
        endif
        if(local_bad==0)then
          allocate(wpw_frozen_ww_projector_cross_row_id,&
            source=wpw_context%ww_projector_cross_row_id,stat=astat);if(astat/=0)local_bad=1
        endif
        if(local_bad==0)then
          allocate(wpw_frozen_ww_projector_cross_col_id,&
            source=wpw_context%ww_projector_cross_col_id,stat=astat);if(astat/=0)local_bad=1
        endif
        if(local_bad==0)then;allocate(wpw_frozen_ww_h0_dense,source=wpw_context%bounded_operator%ww_h0_dense,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_ww_interface_dense,source=wpw_context%bounded_operator%ww_interface_dense,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_wp_h0_dense,source=wpw_context%bounded_operator%wp_h0_dense,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_wp_interface_dense,source=wpw_context%bounded_operator%wp_interface_dense,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_pp_h0_dense,source=wpw_context%bounded_operator%pp_h0_dense,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_pp_interface_dense,source=wpw_context%bounded_operator%pp_interface_dense,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_owned_w_ids,source=wpw_context%bounded_operator%owned_w_ids,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_required_w_ids,source=wpw_context%bounded_operator%required_w_ids,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_owned_p_ids,source=wpw_context%bounded_operator%owned_p_ids,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then;allocate(wpw_frozen_required_p_ids,source=wpw_context%bounded_operator%required_p_ids,stat=astat);if(astat/=0)local_bad=1;endif
        if(local_bad==0)then
          wpw_frozen_halo_epoch=wpw_context%halo_epoch
          wpw_frozen_scan_epoch=wpw_context%scan_epoch
          wpw_frozen_operator_epoch=wpw_context%operator_epoch
          wpw_frozen_callbacks_bound=wpw_context%callbacks_bound
          wpw_frozen_operator_valid=wpw_context%operator_valid
          wpw_frozen_ww_projector_nonlocal_valid=wpw_context%ww_projector_nonlocal_valid
        endif
      endif
      if(dc%id_frag==0)then
        call synchronize_wpw_projector_root_info(local_bad)
        if(local_bad==0)call snapshot_dg_wpw_production_context(&
          wpw_context,wpw_frozen_production_context,local_bad)
      endif
      snapshot_info=local_bad
      if(.not.wpw_potential_stage_ok(snapshot_info))snapshot_info=1
      if(snapshot_info/=0)then
        call release_wpw_frozen_h_state()
      endif
    end subroutine snapshot_wpw_frozen_h_state

    subroutine snapshot_wpw_frozen_ww_components(info)
      integer,intent(out)::info
      integer::astat
      info=1
      allocate(wpw_frozen_ww_components%owned_w_ids,source=wpw_ww_components%owned_w_ids,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%kinetic,source=wpw_ww_components%kinetic,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%potential,source=wpw_ww_components%potential,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%nonlocal,source=wpw_ww_components%nonlocal,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%face_self,source=wpw_ww_components%face_self,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%cross_face_id,source=wpw_ww_components%cross_face_id,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%cross_row_id,source=wpw_ww_components%cross_row_id,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%cross_col_id,source=wpw_ww_components%cross_col_id,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%cross_axis,source=wpw_ww_components%cross_axis,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%cross_side,source=wpw_ww_components%cross_side,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%cross_image,source=wpw_ww_components%cross_image,stat=astat);if(astat/=0)return
      allocate(wpw_frozen_ww_components%cross_value,source=wpw_ww_components%cross_value,stat=astat);if(astat/=0)return
      wpw_frozen_ww_components%metric_convention=wpw_ww_components%metric_convention
      wpw_frozen_ww_components%provenance_fingerprint=wpw_ww_components%provenance_fingerprint
      wpw_frozen_ww_components%valid=wpw_ww_components%valid
      info=0
    end subroutine snapshot_wpw_frozen_ww_components

    subroutine release_wpw_frozen_h_state()
      if(allocated(wpw_frozen_rho_tot))deallocate(wpw_frozen_rho_tot)
      if(allocated(wpw_frozen_vh_tot))deallocate(wpw_frozen_vh_tot)
      if(allocated(wpw_frozen_vxc_tot))deallocate(wpw_frozen_vxc_tot)
      if(allocated(wpw_frozen_vloc_tot))deallocate(wpw_frozen_vloc_tot)
      if(allocated(wpw_frozen_wp_volume))deallocate(wpw_frozen_wp_volume)
      if(allocated(wpw_frozen_wp_nonlocal))deallocate(wpw_frozen_wp_nonlocal)
      if(allocated(wpw_frozen_wp_face))deallocate(wpw_frozen_wp_face)
      if(allocated(wpw_frozen_pp_volume))deallocate(wpw_frozen_pp_volume)
      if(allocated(wpw_frozen_pp_nonlocal))deallocate(wpw_frozen_pp_nonlocal)
      if(allocated(wpw_frozen_ww_projector_nonlocal))deallocate(wpw_frozen_ww_projector_nonlocal)
      if(allocated(wpw_frozen_ww_projector_cross_value))deallocate(wpw_frozen_ww_projector_cross_value)
      if(allocated(wpw_frozen_ww_projector_cross_row_id))deallocate(wpw_frozen_ww_projector_cross_row_id)
      if(allocated(wpw_frozen_ww_projector_cross_col_id))deallocate(wpw_frozen_ww_projector_cross_col_id)
      if(allocated(wpw_frozen_ww_h0_dense))deallocate(wpw_frozen_ww_h0_dense)
      if(allocated(wpw_frozen_ww_interface_dense))deallocate(wpw_frozen_ww_interface_dense)
      if(allocated(wpw_frozen_wp_h0_dense))deallocate(wpw_frozen_wp_h0_dense)
      if(allocated(wpw_frozen_wp_interface_dense))deallocate(wpw_frozen_wp_interface_dense)
      if(allocated(wpw_frozen_pp_h0_dense))deallocate(wpw_frozen_pp_h0_dense)
      if(allocated(wpw_frozen_pp_interface_dense))deallocate(wpw_frozen_pp_interface_dense)
      if(allocated(wpw_frozen_owned_w_ids))deallocate(wpw_frozen_owned_w_ids)
      if(allocated(wpw_frozen_required_w_ids))deallocate(wpw_frozen_required_w_ids)
      if(allocated(wpw_frozen_owned_p_ids))deallocate(wpw_frozen_owned_p_ids)
      if(allocated(wpw_frozen_required_p_ids))deallocate(wpw_frozen_required_p_ids)
      if(allocated(wpw_frozen_ww_components%owned_w_ids))deallocate(wpw_frozen_ww_components%owned_w_ids)
      if(allocated(wpw_frozen_ww_components%kinetic))deallocate(wpw_frozen_ww_components%kinetic)
      if(allocated(wpw_frozen_ww_components%potential))deallocate(wpw_frozen_ww_components%potential)
      if(allocated(wpw_frozen_ww_components%nonlocal))deallocate(wpw_frozen_ww_components%nonlocal)
      if(allocated(wpw_frozen_ww_components%face_self))deallocate(wpw_frozen_ww_components%face_self)
      if(allocated(wpw_frozen_ww_components%cross_face_id))deallocate(wpw_frozen_ww_components%cross_face_id)
      if(allocated(wpw_frozen_ww_components%cross_row_id))deallocate(wpw_frozen_ww_components%cross_row_id)
      if(allocated(wpw_frozen_ww_components%cross_col_id))deallocate(wpw_frozen_ww_components%cross_col_id)
      if(allocated(wpw_frozen_ww_components%cross_axis))deallocate(wpw_frozen_ww_components%cross_axis)
      if(allocated(wpw_frozen_ww_components%cross_side))deallocate(wpw_frozen_ww_components%cross_side)
      if(allocated(wpw_frozen_ww_components%cross_image))deallocate(wpw_frozen_ww_components%cross_image)
      if(allocated(wpw_frozen_ww_components%cross_value))deallocate(wpw_frozen_ww_components%cross_value)
      call release_dg_wpw_production_context_snapshot(wpw_frozen_production_context)
      wpw_frozen_ww_components%valid=.false.
    end subroutine release_wpw_frozen_h_state

    logical function wpw_ww_components_fully_allocated()result(ok)
      ok=wpw_ww_components%valid.and.allocated(wpw_ww_components%owned_w_ids).and.&
        allocated(wpw_ww_components%kinetic).and.allocated(wpw_ww_components%potential).and.&
        allocated(wpw_ww_components%nonlocal).and.allocated(wpw_ww_components%face_self).and.&
        allocated(wpw_ww_components%cross_face_id).and.allocated(wpw_ww_components%cross_row_id).and.&
        allocated(wpw_ww_components%cross_col_id).and.allocated(wpw_ww_components%cross_axis).and.&
        allocated(wpw_ww_components%cross_side).and.allocated(wpw_ww_components%cross_image).and.&
        allocated(wpw_ww_components%cross_value)
    end function wpw_ww_components_fully_allocated

    subroutine validate_wpw_frozen_h_state(validation_info)
      integer,intent(out)::validation_info
      integer::local_bad,operator_info

      operator_info=0
      if(dc%id_frag==0)call validate_dg_wpw_production_context_snapshot(&
        wpw_context,wpw_frozen_production_context,operator_info,allow_interface_lambda_change=.true.)
      local_bad=merge(0,1,wpw_frozen_state_local_matches())
      if(operator_info/=0)local_bad=1
      validation_info=local_bad
      if(.not.wpw_potential_stage_ok(validation_info))validation_info=1
    end subroutine validate_wpw_frozen_h_state

    logical function wpw_frozen_state_local_matches()result(ok)
      ok=.false.
      if(.not.same_real3(wpw_frozen_rho_tot,dc%rho_tot%f))return
      if(.not.same_real3(wpw_frozen_vh_tot,dc%vh_tot%f))return
      if(.not.same_real3(wpw_frozen_vxc_tot,dc%vxc_tot(1)%f))return
      if(.not.same_real3(wpw_frozen_vloc_tot,dc%vloc_tot(1)%f))return
      if(dc%id_frag/=0)then;ok=.true.;return;endif
      if(.not.same_complex1(wpw_frozen_wp_volume,wpw_context%wp_h_volume))return
      if(.not.same_complex1(wpw_frozen_wp_nonlocal,wpw_context%wp_h_nonlocal))return
      if(.not.same_complex1(wpw_frozen_wp_face,wpw_context%wp_h_face))return
      if(.not.same_complex1(wpw_frozen_pp_volume,wpw_context%pp_h_volume))return
      if(.not.same_complex1(wpw_frozen_pp_nonlocal,wpw_context%pp_h_nonlocal))return
      if(.not.same_complex2(wpw_frozen_ww_projector_nonlocal,&
        wpw_context%ww_projector_nonlocal))return
      if(.not.same_complex1(wpw_frozen_ww_projector_cross_value,&
        wpw_context%ww_projector_cross_value))return
      if(.not.same_integer1(wpw_frozen_ww_projector_cross_row_id,&
        wpw_context%ww_projector_cross_row_id))return
      if(.not.same_integer1(wpw_frozen_ww_projector_cross_col_id,&
        wpw_context%ww_projector_cross_col_id))return
      if(wpw_frozen_halo_epoch/=wpw_context%halo_epoch.or.&
         wpw_frozen_scan_epoch/=wpw_context%scan_epoch.or.&
         wpw_frozen_operator_epoch/=wpw_context%operator_epoch.or.&
         wpw_frozen_callbacks_bound.neqv.wpw_context%callbacks_bound.or.&
         wpw_frozen_operator_valid.neqv.wpw_context%operator_valid.or.&
         wpw_frozen_ww_projector_nonlocal_valid.neqv.&
           wpw_context%ww_projector_nonlocal_valid)return
      if(.not.same_complex2(wpw_frozen_ww_h0_dense,wpw_context%bounded_operator%ww_h0_dense))return
      if(.not.same_complex2(wpw_frozen_ww_interface_dense,wpw_context%bounded_operator%ww_interface_dense))return
      if(.not.same_complex2(wpw_frozen_wp_h0_dense,wpw_context%bounded_operator%wp_h0_dense))return
      if(.not.same_complex2(wpw_frozen_wp_interface_dense,wpw_context%bounded_operator%wp_interface_dense))return
      if(.not.same_complex2(wpw_frozen_pp_h0_dense,wpw_context%bounded_operator%pp_h0_dense))return
      if(.not.same_complex2(wpw_frozen_pp_interface_dense,wpw_context%bounded_operator%pp_interface_dense))return
      if(.not.same_integer1(wpw_frozen_owned_w_ids,wpw_context%bounded_operator%owned_w_ids))return
      if(.not.same_integer1(wpw_frozen_required_w_ids,wpw_context%bounded_operator%required_w_ids))return
      if(.not.same_integer1(wpw_frozen_owned_p_ids,wpw_context%bounded_operator%owned_p_ids))return
      if(.not.same_integer1(wpw_frozen_required_p_ids,wpw_context%bounded_operator%required_p_ids))return
      ok=.true.
    end function wpw_frozen_state_local_matches

    logical function same_real3(saved,current)result(ok)
      real(8),allocatable,intent(in)::saved(:,:,:)
      real(8),intent(in)::current(:,:,:)
      ok=allocated(saved);if(.not.ok)return
      ok=all(shape(saved)==shape(current));if(.not.ok)return
      ok=all(saved==current)
    end function same_real3

    logical function same_complex1(saved,current)result(ok)
      complex(8),allocatable,intent(in)::saved(:)
      complex(8),allocatable,intent(in)::current(:)
      ok=allocated(saved).and.allocated(current);if(.not.ok)return
      ok=size(saved)==size(current);if(.not.ok)return
      ok=all(saved==current)
    end function same_complex1

    logical function same_complex2(saved,current)result(ok)
      complex(8),allocatable,intent(in)::saved(:,:)
      complex(8),allocatable,intent(in)::current(:,:)
      ok=allocated(saved).and.allocated(current);if(.not.ok)return
      ok=all(shape(saved)==shape(current));if(.not.ok)return
      ok=all(saved==current)
    end function same_complex2

    logical function same_integer1(saved,current)result(ok)
      integer,allocatable,intent(in)::saved(:)
      integer,allocatable,intent(in)::current(:)
      ok=allocated(saved).and.allocated(current);if(.not.ok)return
      ok=size(saved)==size(current);if(.not.ok)return
      ok=all(saved==current)
    end function same_integer1

    logical function same_real2_alloc(saved,current)result(ok)
      real(8),allocatable,intent(in)::saved(:,:),current(:,:)
      ok=allocated(saved).and.allocated(current);if(.not.ok)return
      ok=all(shape(saved)==shape(current));if(.not.ok)return
      ok=all(saved==current)
    end function same_real2_alloc

    logical function same_wpw_ww_components(saved,current)result(ok)
      type(s_dg_wpw_lcfo_ww_components),intent(in)::saved,current
      ok=saved%valid.and.current%valid;if(.not.ok)return
      ok=saved%provenance_fingerprint==current%provenance_fingerprint;if(.not.ok)return
      ok=trim(saved%metric_convention)==trim(current%metric_convention);if(.not.ok)return
      if(.not.same_integer1(saved%owned_w_ids,current%owned_w_ids))then;ok=.false.;return;endif
      if(.not.same_real2_alloc(saved%kinetic,current%kinetic))then;ok=.false.;return;endif
      if(.not.same_real2_alloc(saved%potential,current%potential))then;ok=.false.;return;endif
      if(.not.same_real2_alloc(saved%nonlocal,current%nonlocal))then;ok=.false.;return;endif
      if(.not.same_real2_alloc(saved%face_self,current%face_self))then;ok=.false.;return;endif
      if(.not.same_integer1(saved%cross_face_id,current%cross_face_id))then;ok=.false.;return;endif
      if(.not.same_integer1(saved%cross_row_id,current%cross_row_id))then;ok=.false.;return;endif
      if(.not.same_integer1(saved%cross_col_id,current%cross_col_id))then;ok=.false.;return;endif
      if(.not.same_integer1(saved%cross_axis,current%cross_axis))then;ok=.false.;return;endif
      if(.not.same_integer1(saved%cross_side,current%cross_side))then;ok=.false.;return;endif
      ok=allocated(saved%cross_image).and.allocated(current%cross_image);if(.not.ok)return
      ok=all(shape(saved%cross_image)==shape(current%cross_image));if(.not.ok)return
      ok=all(saved%cross_image==current%cross_image);if(.not.ok)return
      ok=allocated(saved%cross_value).and.allocated(current%cross_value);if(.not.ok)return
      ok=size(saved%cross_value)==size(current%cross_value);if(.not.ok)return
      ok=all(saved%cross_value==current%cross_value)
    end function same_wpw_ww_components

    subroutine run_wpw_fixed_h_relaxation(fixed_info)
      integer,intent(out)::fixed_info
      integer,parameter::fixed_h_max_iter=8
      integer::iter,root_info,ierr
      real(8)::gap,residual,orth,projector,occupied_trace,previous_trace,trace_change
      real(8)::density_baseline,charge_error
      logical::accepted

      fixed_info=1;root_info=0;previous_trace=huge(1d0);accepted=.false.;wpw_fixed_h_qualified=.false.
      wpw_fixed_h_final_residual=huge(1d0);wpw_fixed_h_final_orthogonality=huge(1d0)
      wpw_fixed_h_final_projector=huge(1d0);wpw_fixed_h_density_baseline=huge(1d0)
      wpw_fixed_h_charge_error=huge(1d0)
      call validate_wpw_frozen_h_state(fixed_info)
      if(fixed_info/=0)then
        if(dc%id_tot==0)write(*,'(1x,a,i0)')'[DG-WPW-LOCAL-FAIL] fixed_h_stage=initial_frozen info=',fixed_info
        return
      endif
      if(dc%id_frag==0)call set_dg_wpw_interface_lambda(wpw_context%bounded_operator,0d0,root_info)
      call MPI_Bcast(root_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS.or.root_info/=0)return
      call validate_wpw_frozen_h_state(fixed_info)
      if(fixed_info/=0)return
      call build_wpw_density_carrying_fragment_seed(fixed_info)
      if(fixed_info/=0)then
        if(dc%id_tot==0)write(*,'(1x,a,i0)')'[DG-WPW-LOCAL-FAIL] fixed_h_stage=density_seed info=',fixed_info
        return
      endif
      call validate_wpw_frozen_h_state(fixed_info)
      if(fixed_info/=0)then
        if(dc%id_tot==0)write(*,'(1x,a,i0)')'[DG-WPW-LOCAL-FAIL] fixed_h_stage=post_seed_frozen info=',fixed_info
        return
      endif
      do iter=1,fixed_h_max_iter
        root_info=0;gap=0d0;residual=huge(1d0);orth=huge(1d0);projector=huge(1d0)
        call validate_wpw_frozen_h_state(fixed_info)
        if(fixed_info/=0)return
        if(dc%id_frag==0)then
          call run_dg_wpw_matrix_free_algebra_step(wpw_context,wpw_production_comm,wpw_apply_h,&
            wpw_apply_s,wpw_global_gram,size(wpw_qw,1),size(wpw_qp,1),wpw_nocc,wpw_nretain,&
            iter+1,dg_wpw_metric_cutoff,dg_wpw_scf_residual_tolerance,wpw_qw,wpw_qp,&
            wpw_q_old_occ,wpw_eigenvalues,gap,residual,orth,projector,root_info,wpw_precondition)
          if(root_info==0)then
            occupied_trace=sum(wpw_occupations*wpw_eigenvalues(1:wpw_nocc))
            trace_change=merge(abs(occupied_trace-previous_trace),huge(1d0),iter>1)
            accepted=iter>1.and.gap>=dg_wpw_gap_threshold.and.&
              max(residual,max(orth,projector))<dg_wpw_scf_residual_tolerance.and.&
              trace_change<dg_wpw_scf_residual_tolerance.and.ieee_is_finite(occupied_trace)
            previous_trace=occupied_trace
            wpw_q_old_occ(1:size(wpw_qw,1),:)=wpw_qw(:,1:wpw_nocc)
            wpw_q_old_occ(size(wpw_qw,1)+1:,:)=wpw_qp(:,1:wpw_nocc)
          endif
        endif
        call MPI_Bcast(root_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
        if(ierr/=MPI_SUCCESS.or.root_info/=0)then
          if(dc%id_tot==0)write(*,'(1x,a,i0,a,i0,a,i0)')&
            '[DG-WPW-LOCAL-FAIL] fixed_h_stage=algebra iter=',iter,' info=',root_info,' mpi=',ierr
          return
        endif
        call MPI_Bcast(accepted,1,MPI_LOGICAL,0,dc%icomm_frag,ierr)
        if(ierr==MPI_SUCCESS)call MPI_Bcast(wpw_qw,size(wpw_qw),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
        if(ierr==MPI_SUCCESS)call MPI_Bcast(wpw_qp,size(wpw_qp),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
        if(ierr==MPI_SUCCESS)call MPI_Bcast(wpw_eigenvalues,size(wpw_eigenvalues),MPI_DOUBLE_PRECISION,0,&
          dc%icomm_frag,ierr)
        if(ierr/=MPI_SUCCESS)return
        call validate_wpw_frozen_h_state(fixed_info)
        if(fixed_info/=0)return
        if(dc%id_tot==0)write(*,'(1x,a,i0,6(a,es12.4),a,l1)')'[DG-WPW-FIXED-H] iter=',iter,&
          ' lambda=',0d0,' trace=',occupied_trace,' trace_change=',trace_change,' residual=',residual,&
          ' orth=',orth,' projector=',projector,' accepted=',accepted
        if(accepted)then
          call diagnose_wpw_fixed_h_density(density_baseline,charge_error,fixed_info)
          if(fixed_info/=0)return
          call validate_wpw_frozen_h_state(fixed_info)
          if(fixed_info/=0)return
          if(dc%id_tot==0)write(*,'(1x,a,2(a,es12.4))')'[DG-WPW-FIXED-H-DENSITY]',&
            ' representation_baseline=',density_baseline,' charge_error=',charge_error
          call continue_wpw_fixed_h_interface(fixed_info)
          if(fixed_info/=0)return
          call validate_wpw_frozen_h_state(fixed_info)
          if(fixed_info/=0)return
          call diagnose_wpw_fixed_h_density(wpw_fixed_h_density_baseline,wpw_fixed_h_charge_error,fixed_info)
          if(fixed_info/=0.or..not.ieee_is_finite(wpw_fixed_h_density_baseline).or.&
            wpw_fixed_h_density_baseline>1d0.or.wpw_fixed_h_charge_error>dg_wpw_scf_residual_tolerance)return
          wpw_fixed_h_qualified=.true.
          if(wpw_metric_diagnostic_only)then
            if(dc%id_tot==0)write(*,'(1x,a,5(a,es12.4))')'[DG-WPW-DIAGNOSTIC-COMPLETE]',&
              ' metric_residual=',wpw_projection_residual,&
              ' captured_norm=',wpw_projection_captured_norm,&
              ' fixed_h_residual=',wpw_fixed_h_final_residual,&
              ' density_baseline=',wpw_fixed_h_density_baseline,&
              ' charge_error=',wpw_fixed_h_charge_error
            fixed_info=0
          else
            call publish_wpw_production_checkpoint(fixed_info)
          endif
          return
        endif
      enddo
      if(dc%id_tot==0)write(*,'(1x,a,i0)')'[DG-WPW-LOCAL-FAIL] fixed_h_max_iter=',fixed_h_max_iter
      fixed_info=2
    end subroutine run_wpw_fixed_h_relaxation

    subroutine continue_wpw_fixed_h_interface(continuation_info)
      integer,intent(out)::continuation_info
      integer,parameter::max_continuation_trials=16
      real(8),parameter::minimum_lambda_step=1d0/64d0
      complex(8),allocatable::accepted_qw(:,:),accepted_qp(:,:),accepted_old_occ(:,:)
      real(8),allocatable::accepted_eigenvalues(:)
      real(8)::accepted_lambda,trial_lambda,lambda_step,accepted_merit,trial_merit
      real(8)::gap,residual,orth,projector
      integer::trial,root_info,global_info,ierr
      logical::trial_accepted

      continuation_info=1;accepted_lambda=0d0;lambda_step=0.5d0
      accepted_merit=dg_wpw_scf_residual_tolerance
      allocate(accepted_qw,source=wpw_qw)
      allocate(accepted_qp,source=wpw_qp)
      allocate(accepted_old_occ,source=wpw_q_old_occ)
      allocate(accepted_eigenvalues,source=wpw_eigenvalues)
      do trial=1,max_continuation_trials
        if(accepted_lambda>=1d0-10d0*epsilon(1d0))then;continuation_info=0;return;endif
        trial_lambda=min(1d0,accepted_lambda+lambda_step);root_info=0;trial_accepted=.false.
        call validate_wpw_frozen_h_state(continuation_info);if(continuation_info/=0)return
        if(dc%id_frag==0)call set_dg_wpw_interface_lambda(wpw_context%bounded_operator,trial_lambda,root_info)
        call MPI_Bcast(root_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
        if(ierr/=MPI_SUCCESS.or.root_info/=0)return
        call validate_wpw_frozen_h_state(continuation_info);if(continuation_info/=0)return
        if(dc%id_frag==0)then
          call run_dg_wpw_matrix_free_algebra_step(wpw_context,wpw_production_comm,&
            wpw_apply_h,wpw_apply_s,wpw_global_gram,size(wpw_qw,1),size(wpw_qp,1),wpw_nocc,&
            wpw_nretain,trial+1,dg_wpw_metric_cutoff,dg_wpw_scf_residual_tolerance,wpw_qw,wpw_qp,&
            accepted_old_occ,wpw_eigenvalues,gap,residual,orth,projector,root_info,wpw_precondition)
          if(root_info==0)then
            trial_merit=max(residual,orth)
            trial_accepted=ieee_is_finite(trial_merit).and.ieee_is_finite(projector).and.&
              ieee_is_finite(gap).and.gap>=dg_wpw_gap_threshold.and.&
              trial_merit<=max(10d0*accepted_merit,dg_wpw_scf_residual_tolerance)
            root_info=merge(0,1,trial_accepted)
          endif
          call MPI_Allreduce(root_info,global_info,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr)
          if(ierr/=MPI_SUCCESS)global_info=1
          root_info=global_info;trial_accepted=root_info==0
        endif
        call MPI_Bcast(root_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
        if(ierr/=MPI_SUCCESS)return
        call MPI_Bcast(trial_accepted,1,MPI_LOGICAL,0,dc%icomm_frag,ierr)
        if(ierr/=MPI_SUCCESS)return
        if(trial_accepted)then
          accepted_lambda=trial_lambda;accepted_merit=trial_merit
          accepted_qw=wpw_qw;accepted_qp=wpw_qp;accepted_eigenvalues=wpw_eigenvalues
          accepted_old_occ(1:size(wpw_qw,1),:)=wpw_qw(:,1:wpw_nocc)
          accepted_old_occ(size(wpw_qw,1)+1:,:)=wpw_qp(:,1:wpw_nocc)
          wpw_q_old_occ=accepted_old_occ
          lambda_step=min(1d0-accepted_lambda,1.5d0*lambda_step)
          if(accepted_lambda>=1d0-10d0*epsilon(1d0))then
            wpw_fixed_h_final_residual=residual
            wpw_fixed_h_final_orthogonality=orth
            wpw_fixed_h_final_projector=projector
          endif
        else
          wpw_qw=accepted_qw;wpw_qp=accepted_qp;wpw_eigenvalues=accepted_eigenvalues
          wpw_q_old_occ=accepted_old_occ
          if(dc%id_frag==0)call set_dg_wpw_interface_lambda(wpw_context%bounded_operator,&
            accepted_lambda,root_info)
          call MPI_Bcast(root_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
          if(ierr/=MPI_SUCCESS.or.root_info/=0)return
          lambda_step=0.5d0*lambda_step
          if(lambda_step<minimum_lambda_step)return
        endif
        call MPI_Bcast(wpw_qw,size(wpw_qw),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
        if(ierr==MPI_SUCCESS)call MPI_Bcast(wpw_qp,size(wpw_qp),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
        if(ierr==MPI_SUCCESS)call MPI_Bcast(wpw_eigenvalues,size(wpw_eigenvalues),MPI_DOUBLE_PRECISION,0,&
          dc%icomm_frag,ierr)
        if(ierr/=MPI_SUCCESS)return
        call validate_wpw_frozen_h_state(continuation_info);if(continuation_info/=0)return
        if(dc%id_tot==0)write(*,'(1x,a,i0,3(a,es12.4),a,l1)')'[DG-WPW-CONTINUATION] trial=',trial,&
          ' lambda=',trial_lambda,' accepted_lambda=',accepted_lambda,' merit=',accepted_merit,&
          ' accepted=',trial_accepted
      enddo
    end subroutine continue_wpw_fixed_h_interface

    subroutine diagnose_wpw_fixed_h_density(representation_baseline,charge_error,diagnostic_info)
      real(8),intent(out)::representation_baseline,charge_error
      integer,intent(out)::diagnostic_info
      complex(8),allocatable::rw_support(:,:),rw(:,:),rp(:,:)
      real(8),allocatable::rho_raw(:)
      real(8)::charge_local,totals_local(3),totals_global(3)
      integer::iw,wpos,stage_info,ierr
      diagnostic_info=1;representation_baseline=huge(1d0);charge_error=huge(1d0)
      allocate(rw_support(size(wpw_support_w_ids),wpw_nocc),rw(size(wpw_owned_w_ids),wpw_nocc),&
        rp(size(wpw_support_fragment_ids)*size(wpw_g_indices,2),wpw_nocc),&
        rho_raw(wpw_volume_accumulator%npoint))
      rw_support=0;rw=0;rp=0;stage_info=0
      if(dc%id_frag==0)call fetch_dg_wpw_support_coefficients(wpw_context%bounded_operator,&
        wpw_qw(:,1:wpw_nocc),wpw_qp(:,1:wpw_nocc),rw_support,rp,stage_info)
      call MPI_Bcast(stage_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS.or.stage_info/=0)return
      call MPI_Bcast(rw_support,size(rw_support),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(rp,size(rp),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)return
      do iw=1,size(wpw_owned_w_ids)
        wpos=find_integer_id(wpw_support_w_ids,wpw_owned_w_ids(iw));if(wpos==0)return
        rw(iw,:)=rw_support(wpos,:)
      enddo
      call build_dg_wpw_core_density(wpw_volume_accumulator%w_points(:,1:size(rho_raw)),&
        wpw_volume_accumulator%p_points(:,1:size(rho_raw)),rw,rp,wpw_occupations,&
        wpw_volume_accumulator%weights(1:size(rho_raw)),rho_raw,charge_local,stage_info)
      if(stage_info/=0)return
      totals_local=[sum((rho_raw-wpw_volume_accumulator%densities(1:size(rho_raw)))**2),&
        sum(rho_raw**2),charge_local]
      call MPI_Allreduce(totals_local,totals_global,3,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or..not.all(ieee_is_finite(totals_global)))return
      representation_baseline=sqrt(totals_global(1)/max(1d-30,totals_global(2)))
      charge_error=abs(totals_global(3)-wpw_projection_charge)
      if(.not.all(ieee_is_finite([representation_baseline,charge_error])))return
      diagnostic_info=0
    end subroutine diagnose_wpw_fixed_h_density

    logical function wpw_seed_collective_stage_ok(stage,local_info)result(ok)
      character(*),intent(in)::stage
      integer,intent(in)::local_info
      integer::local_bad,global_bad,ierr
      local_bad=merge(0,1,local_info==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
      ok=ierr==MPI_SUCCESS.and.global_bad==0
      if(.not.ok.and.dc%id_tot==0)write(*,'(1x,a,a,a,i0,a,i0)')&
        '[DG-WPW-SEED-FAIL] stage=',trim(stage),' global_bad=',global_bad,' mpi=',ierr
    end function wpw_seed_collective_stage_ok

    subroutine build_core_owned_projected_wannier_density_seed(source_values,source_count,&
        source_condition,occupied_w_basis,source_info)
      use communication,only:comm_summation
      use salmon_global,only:wannier_projection_width
      use inputoutput,only:au_length_aa
      complex(8),allocatable,intent(out)::source_values(:,:)
      integer,intent(out)::source_count,source_info
      real(8),intent(out)::source_condition
      type(s_dg_wpw_occupied_w_basis),intent(inout)::occupied_w_basis
      real(8),allocatable::bond_center(:,:)
      integer,allocatable::bond_atoms(:,:),bond_images(:,:)
      complex(8),allocatable::overlap_local(:,:),overlap(:,:),occupied_values(:,:),&
        source_partial(:,:),source_core(:,:),polar_transform(:,:),occupied_buffer(:,:),&
        buffer_partial(:,:),buffer_values(:,:),occupied_stencil_partial(:,:,:,:),&
        occupied_storage(:,:,:,:),occupied_stencil(:,:,:,:),&
        occupied_w_p(:,:,:,:),canonical_occupied_w(:,:,:,:),&
        buffer_gradient_partial(:,:,:),buffer_gradient_values(:,:,:),&
        compact_gradient_partial(:,:,:),compact_gradient_values(:,:,:),&
        descriptor_buffer_values(:,:),descriptor_buffer_gradients(:,:,:),send_value(:),received_value(:)
      integer::source_counts(dc%n_frag),source_counts_all(dc%n_frag),source_offset,global_source_count
      integer,allocatable::source_key_local(:,:),source_key_partial(:,:),source_key_global(:,:)
      integer,allocatable::local_buffer_point_ids(:),gathered_core_grid_ids(:)
      complex(8),allocatable::gathered_core_values(:,:),gathered_buffer_values(:,:),&
        gathered_buffer_gradients(:,:,:)
      integer,allocatable::send_source(:,:),send_destination(:),send_point(:),send_image(:,:),&
        received_source(:,:),received_point(:),received_image(:,:),all_rank_meta(:,:),&
        fragment_start(:,:),fragment_shape(:,:)
      real(8)::x,y,z,gval,wrapped_center(3),fractional_center(3),core_lower(3),core_upper(3)
      real(8)::tail_norm_local(2),tail_norm_global(2),outer_shell_ratio,omitted_tail_ratio
      real(8),allocatable::spread_coordinates(:,:),spread_bvec(:,:),spread_weight(:),&
        spread_norm(:),spread_norm_partial(:),spread_norm_global(:),spread_center(:,:),spread_omega(:),&
        spread_omega_a2(:),spread_width_a(:),spread_width_sorted(:)
      complex(8),allocatable::spread_link(:,:),spread_link_partial(:,:),spread_link_global(:,:)
      logical,allocatable::spread_center_valid(:)
      integer,allocatable::spread_order(:),spread_stable_id(:),spread_owner(:)
      integer::ixp,iyp,izp,ixb,iyb,izb,sx,sy,sz,io_state,isource,jsource,point,local_info,center_image(3),&
        canonical_atoms(2),canonical_image(3),buffer_point,buffer_count,record,max_record,&
        destination_fragment,destination_local(3),destination_rank,global_point(3),global_grid_id,&
        source_index,core_point,core_descriptor_point,rank_meta(8),ierr,stencil_radius,stencil_extent(3),stencil_lo(3),&
        descriptor_buffer_lo(3),descriptor_buffer_hi(3),descriptor_gradient_lo(3),&
        descriptor_gradient_hi(3),descriptor_grid(3),gradient_point,gradient_extent(3),&
        unwrapped_point(3),periodic_image(3),local_grid(3),axis,&
        full_cell_coverage_local,full_cell_coverage_global
      integer::iw,jw,kw,widest_w,count_above_1p2a,valid_spread_count,comparison_axis,owner_offset,&
        worst_norm_w
      real(8)::spread_swap,median_a,p90_a
      logical::spread_key_less
      logical::center_owned,full_cell_coverage,tail_warning
      type(t_dg_wpw_periodic_image_mismatch)::periodic_mismatch

      source_info=1;source_count=0;source_condition=huge(1d0)
      call build_local_bond_center_projection_map(source_count,bond_center,bond_atoms,bond_images)
      local_info=merge(0,1,source_count>0.and.wpw_core_p_accumulator%npoint>0)
      core_lower=dc%rxyz_frag(:,dc%i_frag)/wpw_box_length
      core_upper=core_lower+system%hgs*dble(dc%nxyz_domain)/wpw_box_length
      do isource=1,merge(source_count,0,local_info==0)
        fractional_center=bond_center(:,isource)/wpw_box_length
        call canonicalize_sawf_wannier_center(fractional_center,[1d0,1d0,1d0],core_lower,&
          core_upper,1d-10,wrapped_center,center_image,center_owned,local_info)
        if(local_info/=0.or..not.center_owned)then;local_info=1;exit;endif
        call canonicalize_sawf_bond_identity(bond_atoms(1,isource),bond_atoms(2,isource),&
          bond_images(:,isource),canonical_atoms,canonical_image,local_info)
        if(local_info/=0)exit
        bond_atoms(:,isource)=canonical_atoms;bond_images(:,isource)=canonical_image
        do jsource=1,isource-1
          if(all(bond_atoms(:,jsource)==bond_atoms(:,isource)).and.&
            all(bond_images(:,jsource)==bond_images(:,isource)))local_info=1
        enddo
        if(local_info/=0)exit
      enddo
      if(.not.wpw_seed_collective_stage_ok('canonical_sources',local_info))return
      source_counts=0
      if(dc%id_frag==0)source_counts(dc%i_frag)=source_count
      call comm_summation(source_counts,source_counts_all,dc%n_frag,dc%icomm_tot)
      global_source_count=sum(source_counts_all);source_offset=sum(source_counts_all(1:dc%i_frag-1))
      if(global_source_count/=occupied_index_from_input(1))return
      allocate(source_key_local(5,source_count),source_key_partial(5,global_source_count),&
        source_key_global(5,global_source_count))
      source_key_local(1:2,:)=bond_atoms;source_key_local(3:5,:)=bond_images
      if(allocated(wpw_bootstrap_source_keys))deallocate(wpw_bootstrap_source_keys)
      allocate(wpw_bootstrap_source_keys,source=source_key_local)
      source_key_partial=0
      if(dc%id_frag==0)source_key_partial(:,source_offset+1:source_offset+source_count)=source_key_local
      call comm_summation(source_key_partial,source_key_global,size(source_key_global),dc%icomm_tot)
      do isource=1,global_source_count;do jsource=isource+1,global_source_count
        if(all(source_key_global(:,isource)==source_key_global(:,jsource)))return
      enddo;enddo
      allocate(overlap_local(system%no,source_count),overlap(system%no,source_count))
      overlap_local=(0d0,0d0);overlap=(0d0,0d0)
      do izp=mg%is(3),mg%ie(3)
        z=dc%lg_tot%coordinate(dc%jxyz_tot(izp,3),3)
        do iyp=mg%is(2),mg%ie(2)
          y=dc%lg_tot%coordinate(dc%jxyz_tot(iyp,2),2)
          do ixp=mg%is(1),mg%ie(1)
            x=dc%lg_tot%coordinate(dc%jxyz_tot(ixp,1),1)
            do isource=1,source_count
              gval=bond_center_projection_value_local_periodic(x,y,z,bond_center(:,isource),&
                wannier_projection_width,wpw_box_length)
              if(gval==0d0)cycle
              do io_state=info%io_s,info%io_e
                if(system%rocc(io_state,1,1)<=1d-12)cycle
                overlap_local(io_state,isource)=overlap_local(io_state,isource)+&
                  cmplx(spsi%rwf(ixp,iyp,izp,1,io_state,1,1)*gval*system%hvol,0d0,8)
              enddo
            enddo
          enddo
        enddo
      enddo
      call comm_summation(overlap_local,overlap,size(overlap),info%icomm_rko)
      allocate(occupied_values(wpw_core_p_accumulator%npoint,system%no),&
        source_partial(wpw_core_p_accumulator%npoint,source_count),&
        source_core(wpw_core_p_accumulator%npoint,source_count),&
        source_values(wpw_core_p_accumulator%npoint,global_source_count),&
        polar_transform(source_count,source_count))
      occupied_values=(0d0,0d0);point=0
      do izp=max(mg%is(3),1),min(mg%ie(3),size(f_basis,3))
        do iyp=max(mg%is(2),1),min(mg%ie(2),size(f_basis,2))
          do ixp=max(mg%is(1),1),min(mg%ie(1),size(f_basis,1))
            point=point+1
            do io_state=info%io_s,info%io_e
              if(system%rocc(io_state,1,1)<=1d-12)cycle
              occupied_values(point,io_state)=cmplx(spsi%rwf(ixp,iyp,izp,1,io_state,1,1),0d0,8)
            enddo
          enddo
        enddo
      enddo
      local_info=merge(0,1,point==wpw_core_p_accumulator%npoint)
      if(local_info==0)call build_sawf_projected_wannier_from_overlap(occupied_values,overlap,system%hvol,&
        dg_wpw_metric_cutoff,source_partial,polar_transform,source_condition,local_info)
      if(.not.wpw_seed_collective_stage_ok('core_projection',local_info))return
      call comm_summation(source_partial,source_core,size(source_core),info%icomm_o)
      source_values=(0d0,0d0)
      source_values(:,source_offset+1:source_offset+source_count)=source_core

      buffer_count=product(dc%nxyz_domain+2*dc%nxyz_buffer)
      allocate(occupied_buffer(buffer_count,system%no),buffer_partial(buffer_count,source_count),&
        buffer_values(buffer_count,source_count),buffer_gradient_partial(3,buffer_count,source_count),&
        buffer_gradient_values(3,buffer_count,source_count))
      stencil_radius=size(stencil%coef_nab,1)
      local_info=merge(0,1,all(dc%nxyz_buffer>stencil_radius))
      stencil_lo=1;stencil_extent=dc%nxyz_domain+2*dc%nxyz_buffer
      gradient_extent=stencil_extent-2*stencil_radius
      if(any(gradient_extent<=0))local_info=1
      if(.not.wpw_seed_collective_stage_ok('buffer_gradient_window',local_info))return
      allocate(occupied_stencil_partial(stencil_extent(1),stencil_extent(2),stencil_extent(3),system%no),&
        occupied_storage(stencil_extent(1),stencil_extent(2),stencil_extent(3),system%no),&
        occupied_stencil(stencil_extent(1),stencil_extent(2),stencil_extent(3),system%no))
      occupied_stencil_partial=(0d0,0d0)
      do izp=max(1,mg%is(3)),min(stencil_extent(3),mg%ie(3))
        do iyp=max(1,mg%is(2)),min(stencil_extent(2),mg%ie(2))
          do ixp=max(1,mg%is(1)),min(stencil_extent(1),mg%ie(1))
            do io_state=info%io_s,info%io_e
              if(system%rocc(io_state,1,1)<=1d-12)cycle
              occupied_stencil_partial(ixp,iyp,izp,io_state)=&
                cmplx(spsi%rwf(ixp,iyp,izp,1,io_state,1,1),0d0,8)
            enddo
          enddo
        enddo
      enddo
      call comm_summation(occupied_stencil_partial,occupied_storage,size(occupied_storage),info%icomm_r)
      call reorder_dg_wpw_fragment_buffer(occupied_storage,dc%nxyz_domain,dc%nxyz_buffer,&
        occupied_stencil,local_info)
      if(.not.wpw_seed_collective_stage_ok('unwrapped_buffer_order',local_info))return
      call periodize_dg_wpw_fragment_buffer(occupied_stencil,dc%nxyz_domain,dc%nxyz_buffer,&
        dc%lg_tot%num,dc%ixyz_frag(:,dc%i_frag),local_info)
      if(.not.wpw_seed_collective_stage_ok('physical_periodic_p',local_info))return
      occupied_buffer=(0d0,0d0);buffer_point=0
      do izb=1,stencil_extent(3)
        do iyb=1,stencil_extent(2)
          do ixb=1,stencil_extent(1)
        buffer_point=buffer_point+1
        do io_state=info%io_s,info%io_e
          if(system%rocc(io_state,1,1)<=1d-12)cycle
          occupied_buffer(buffer_point,io_state)=occupied_stencil(ixb,iyb,izb,io_state)
        enddo
      enddo;enddo;enddo
      call apply_sawf_projected_wannier_transform(occupied_buffer,overlap,polar_transform,&
        buffer_partial,local_info)
      if(.not.wpw_seed_collective_stage_ok('buffer_projection',local_info))return
      call comm_summation(buffer_partial,buffer_values,size(buffer_values),info%icomm_o)
      local_info=merge(0,1,global_source_count==128.and.dc%n_frag==8.and.source_count==16)
      if(.not.wpw_seed_collective_stage_ok('occupied_w_link_shape',local_info))return
      allocate(occupied_w_p(stencil_extent(1),stencil_extent(2),stencil_extent(3),source_count),&
        canonical_occupied_w(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),source_count),&
        spread_coordinates(3,product(dc%lg_tot%num)),spread_bvec(3,6),spread_weight(6),&
        spread_norm(source_count),spread_link(source_count,6),spread_center(3,global_source_count),&
        spread_omega(global_source_count),spread_center_valid(global_source_count),&
        spread_omega_a2(global_source_count),spread_width_a(global_source_count),&
        spread_width_sorted(global_source_count),&
        spread_order(global_source_count),spread_stable_id(global_source_count),&
        spread_owner(global_source_count),&
        spread_norm_partial(global_source_count),spread_norm_global(global_source_count),&
        spread_link_partial(global_source_count,6),spread_link_global(global_source_count,6))
      occupied_w_p=reshape(buffer_values,shape(occupied_w_p));canonical_occupied_w=(0d0,0d0)
      spread_coordinates=0d0;spread_bvec=0d0;spread_weight=0d0
      spread_norm=0d0;spread_link=(0d0,0d0);spread_center=0d0
      spread_norm_partial=0d0;spread_norm_global=0d0
      spread_link_partial=(0d0,0d0);spread_link_global=(0d0,0d0)
      spread_omega=huge(1d0);spread_center_valid=.false.;local_info=0
      do axis=1,3
        spread_bvec(axis,2*axis-1)=2d0*acos(-1d0)/wpw_box_length(axis)
        spread_bvec(axis,2*axis)=-spread_bvec(axis,2*axis-1)
        spread_weight(2*axis-1:2*axis)=0.5d0/spread_bvec(axis,2*axis-1)**2
      enddo
      if(dc%id_frag==0)call assemble_dg_wpw_canonical_buffer_norm(occupied_w_p,&
        dc%nxyz_domain,dc%nxyz_buffer,dc%lg_tot%num,dc%ixyz_frag(:,dc%i_frag),&
        system%hvol,spread_norm,local_info)
      if(.not.wpw_seed_collective_stage_ok('occupied_w_normalization_norm_assembly',local_info))return
      if(dc%id_frag==0)spread_norm_partial(source_offset+1:source_offset+source_count)=spread_norm
      call comm_summation(spread_norm_partial,spread_norm_global,size(spread_norm_global),dc%icomm_tot)
      local_info=merge(0,1,all(spread_norm_global>0d0).and.all(ieee_is_finite(spread_norm_global)))
      spread_norm=spread_norm_global(source_offset+1:source_offset+source_count)
      if(local_info==0)call validate_sawf_projected_wannier_columns(spread_norm,polar_transform,&
        source_core,buffer_values,local_info)
      if(.not.wpw_seed_collective_stage_ok('occupied_w_normalization_precheck',local_info))return
      call normalize_sawf_projected_wannier_columns(spread_norm,polar_transform,source_core,&
        buffer_values,local_info)
      if(.not.wpw_seed_collective_stage_ok('occupied_w_physical_normalization',local_info))return
      source_values(:,source_offset+1:source_offset+source_count)=source_core
      occupied_w_p=reshape(buffer_values,shape(occupied_w_p));canonical_occupied_w=(0d0,0d0)
      spread_norm=0d0;spread_link=(0d0,0d0)
      spread_norm_partial=0d0;spread_norm_global=0d0
      spread_link_partial=(0d0,0d0);spread_link_global=(0d0,0d0)
      if(dc%id_frag==0)call assemble_dg_wpw_canonical_buffer_norm(occupied_w_p,&
        dc%nxyz_domain,dc%nxyz_buffer,dc%lg_tot%num,dc%ixyz_frag(:,dc%i_frag),&
        system%hvol,spread_norm,local_info)
      if(.not.wpw_seed_collective_stage_ok('occupied_w_post_normalization_norm_assembly',local_info))return
      if(dc%id_frag==0)spread_norm_partial(source_offset+1:source_offset+source_count)=spread_norm
      call comm_summation(spread_norm_partial,spread_norm_global,size(spread_norm_global),dc%icomm_tot)
      worst_norm_w=maxloc(abs(spread_norm_global-1d0),dim=1)
      if(dc%id_tot==0)write(*,'(1x,a,3(a,es13.5),a,i0)')'[DG-WPW-WANNIER-NORM]',&
        ' min=',minval(spread_norm_global),' max=',maxval(spread_norm_global),&
        ' max_abs_unit_error=',maxval(abs(spread_norm_global-1d0)),' row=',worst_norm_w
      local_info=merge(0,1,all(spread_norm_global>0d0).and.&
        all(ieee_is_finite(spread_norm_global)).and.maxval(abs(spread_norm_global-1d0))<=1d-8)
      if(.not.wpw_seed_collective_stage_ok('occupied_w_post_normalization_norm',local_info))return
      full_cell_coverage_local=merge(1,0,all(stencil_extent>=dc%lg_tot%num))
      call MPI_Allreduce(full_cell_coverage_local,full_cell_coverage_global,1,MPI_INTEGER,MPI_MIN,&
        dc%icomm_tot,ierr)
      local_info=merge(0,1,ierr==MPI_SUCCESS)
      if(.not.wpw_seed_collective_stage_ok('full_cell_coverage_globalization',local_info))return
      full_cell_coverage=full_cell_coverage_global==1
      if(full_cell_coverage.and.dc%id_frag==0)then
        call extract_dg_wpw_canonical_cell(occupied_w_p,dc%nxyz_domain,dc%nxyz_buffer,&
          dc%lg_tot%num,dc%ixyz_frag(:,dc%i_frag),canonical_occupied_w,local_info,periodic_mismatch)
        if(periodic_mismatch%valid)write(*,'(1x,a,i0,4(a,3(i0,1x)),a,i0,a,es13.5)')&
          '[DG-WPW-CANONICAL-IMAGE-MISMATCH] fragment=',dc%i_frag,&
          ' canonical=',periodic_mismatch%canonical,' first_p=',periodic_mismatch%first_p,&
          ' second_p=',periodic_mismatch%second_p,' domain=',dc%nxyz_domain,&
          ' w_row=',periodic_mismatch%w_row,' max_abs_diff=',periodic_mismatch%abs_diff
      endif
      if(full_cell_coverage)then
        if(.not.wpw_seed_collective_stage_ok('occupied_w_canonical_cell_extraction',local_info))return
        local_info=0;spread_norm=0d0;spread_link=(0d0,0d0)
        spread_norm_partial=0d0;spread_norm_global=0d0
        spread_link_partial=(0d0,0d0);spread_link_global=(0d0,0d0)
      if(dc%id_frag==0)then
        point=0
        do izp=1,dc%lg_tot%num(3);do iyp=1,dc%lg_tot%num(2);do ixp=1,dc%lg_tot%num(1)
          point=point+1
          spread_coordinates(:,point)=system%hgs*dble([ixp-1,iyp-1,izp-1])
        enddo;enddo;enddo
        call assemble_sawf_diagonal_periodic_links(&
          reshape(canonical_occupied_w,[product(dc%lg_tot%num),source_count]),&
          spread_coordinates,spread_bvec,system%hvol,spread_norm,spread_link,local_info)
      endif
      if(.not.wpw_seed_collective_stage_ok('occupied_w_link_assembly',local_info))return
      if(dc%id_frag==0)then
        spread_norm_partial(source_offset+1:source_offset+source_count)=spread_norm
        spread_link_partial(source_offset+1:source_offset+source_count,:)=spread_link
      endif
      call comm_summation(spread_norm_partial,spread_norm_global,size(spread_norm_global),dc%icomm_tot)
      call comm_summation(spread_link_partial,spread_link_global,size(spread_link_global),dc%icomm_tot)
      local_info=merge(0,1,all(spread_norm_global>0d0).and.&
        all(ieee_is_finite(spread_norm_global)).and.&
        all(ieee_is_finite(real(spread_link_global))).and.&
        all(ieee_is_finite(aimag(spread_link_global))))
      if(.not.wpw_seed_collective_stage_ok('occupied_w_link_collection',local_info))return
      call diagnose_sawf_discrete_wannier_spread(spread_link_global,spread_norm_global,&
        spread_bvec,spread_weight,spread_center,spread_omega,spread_center_valid,local_info,&
        require_unit_norm=.true.)
      if(local_info==0)then
        spread_omega_a2=spread_omega*au_length_aa**2
        spread_width_a=sqrt(max(0d0,spread_omega_a2))
        spread_order=[(iw,iw=1,global_source_count)]
        do iw=1,global_source_count-1
          kw=iw
          do jw=iw+1,global_source_count
            spread_key_less=.false.
            do comparison_axis=1,5
              if(source_key_global(comparison_axis,spread_order(jw))<&
                  source_key_global(comparison_axis,spread_order(kw)))then
                spread_key_less=.true.;exit
              elseif(source_key_global(comparison_axis,spread_order(jw))>&
                  source_key_global(comparison_axis,spread_order(kw)))then
                exit
              endif
            enddo
            if(spread_key_less)kw=jw
          enddo
          if(kw/=iw)then
            jw=spread_order(iw);spread_order(iw)=spread_order(kw);spread_order(kw)=jw
          endif
        enddo
        spread_stable_id=0;spread_owner=0;owner_offset=0
        do iw=1,global_source_count
          spread_stable_id(spread_order(iw))=iw
        enddo
        do iw=1,dc%n_frag
          spread_owner(owner_offset+1:owner_offset+source_counts_all(iw))=iw
          owner_offset=owner_offset+source_counts_all(iw)
        enddo
        valid_spread_count=count(spread_center_valid);spread_width_sorted=0d0;jw=0
        do iw=1,global_source_count
          if(.not.spread_center_valid(iw))cycle
          jw=jw+1;spread_width_sorted(jw)=spread_width_a(iw)
        enddo
        do iw=2,valid_spread_count
          spread_swap=spread_width_sorted(iw);jw=iw-1
          do while(jw>=1)
            if(spread_width_sorted(jw)<=spread_swap)exit
            spread_width_sorted(jw+1)=spread_width_sorted(jw);jw=jw-1
          enddo
          spread_width_sorted(jw+1)=spread_swap
        enddo
        if(valid_spread_count>0)then
          if(mod(valid_spread_count,2)==0)then
            median_a=0.5d0*(spread_width_sorted(valid_spread_count/2)+&
              spread_width_sorted(valid_spread_count/2+1))
          else
            median_a=spread_width_sorted((valid_spread_count+1)/2)
          endif
          p90_a=spread_width_sorted(ceiling(0.9d0*dble(valid_spread_count)))
          count_above_1p2a=count(spread_width_a>1.2d0.and.spread_center_valid)
          widest_w=0
          do iw=1,global_source_count
            if(.not.spread_center_valid(iw))cycle
            if(widest_w==0.or.spread_width_a(iw)>spread_width_a(widest_w))widest_w=iw
          enddo
        else
          median_a=huge(1d0);p90_a=huge(1d0);count_above_1p2a=0;widest_w=0
        endif
        if(dc%id_tot==0)then
          do iw=1,global_source_count
            if(spread_center_valid(iw))then
              write(*,'(1x,a,i0,a,i0,a,l1,5(a,es13.5))')'[DG-WPW-WANNIER-SPREAD] id=',&
                spread_stable_id(iw),' fragment=',spread_owner(iw),' center_valid=',.true.,&
                ' center_x_A=',spread_center(1,iw)*au_length_aa,&
                ' center_y_A=',spread_center(2,iw)*au_length_aa,&
                ' center_z_A=',spread_center(3,iw)*au_length_aa,&
                ' omega_A2=',spread_omega_a2(iw),' width_A=',spread_width_a(iw)
            else
              write(*,'(1x,a,i0,a,i0,a,l1)')'[DG-WPW-WANNIER-SPREAD] id=',&
                spread_stable_id(iw),' fragment=',spread_owner(iw),' center_valid=',.false.
            endif
          enddo
          if(valid_spread_count>0)write(*,'(1x,a,5(a,es13.5),4(a,i0))')&
            '[DG-WPW-WANNIER-SPREAD-SUMMARY]',&
            ' min_A=',spread_width_sorted(1),&
            ' mean_A=',sum(spread_width_a,mask=spread_center_valid)/dble(valid_spread_count),&
            ' median_A=',median_a,' p90_A=',p90_a,' max_A=',spread_width_sorted(valid_spread_count),&
            ' count_above_1p2A=',count_above_1p2a,' valid_count=',valid_spread_count,&
            ' invalid_count=',global_source_count-valid_spread_count,&
            ' widest_id=',merge(spread_stable_id(widest_w),0,widest_w>0)
        endif
        if(.not.all(spread_center_valid))local_info=1
      endif
      if(.not.wpw_seed_collective_stage_ok('occupied_w_spread_diagnostic',local_info))return
      elseif(dc%id_tot==0)then
        write(*,'(1x,a,3(i0,1x),a,3(i0,1x))')'[DG-WPW-WANNIER-SPREAD-SKIPPED] p_extent=',&
          stencil_extent,' physical_extent=',dc%lg_tot%num
      endif
      allocate(compact_gradient_partial(3,product(gradient_extent),source_count),&
        compact_gradient_values(3,product(gradient_extent),source_count))
      call build_sawf_projected_buffer_gradients(occupied_stencil,stencil%coef_nab,overlap,&
        polar_transform,compact_gradient_partial,local_info)
      if(.not.wpw_seed_collective_stage_ok('buffer_gradient_projection',local_info))return
      call comm_summation(compact_gradient_partial,compact_gradient_values,&
        size(compact_gradient_values),info%icomm_o)
      buffer_gradient_values=(0d0,0d0);gradient_point=0
      do izp=1+stencil_radius,stencil_extent(3)-stencil_radius
        do iyp=1+stencil_radius,stencil_extent(2)-stencil_radius
          do ixp=1+stencil_radius,stencil_extent(1)-stencil_radius
            gradient_point=gradient_point+1
            buffer_point=ixp+(iyp-1)*stencil_extent(1)+(izp-1)*stencil_extent(1)*stencil_extent(2)
            buffer_gradient_values(:,buffer_point,:)=compact_gradient_values(:,gradient_point,:)
          enddo
        enddo
      enddo
      tail_norm_local=0d0;buffer_point=0
      if(info%id_o==0)then
        do izb=1,stencil_extent(3);do iyb=1,stencil_extent(2);do ixb=1,stencil_extent(1)
          buffer_point=buffer_point+1
          tail_norm_local(1)=tail_norm_local(1)+sum(abs(buffer_values(buffer_point,:))**2)*system%hvol
          if(is_sawf_outer_buffer_shell([ixb,iyb,izb],stencil_extent,dc%lg_tot%num))&
            tail_norm_local(2)=tail_norm_local(2)+&
              sum(abs(buffer_values(buffer_point,:))**2)*system%hvol
        enddo;enddo;enddo
      endif
      call MPI_Allreduce(tail_norm_local,tail_norm_global,2,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
      local_info=merge(0,1,ierr==MPI_SUCCESS)
      if(local_info==0.and.dc%id_tot==0)write(*,'(1x,a,4(a,es12.4))')&
        '[DG-WPW-SEED-BUFFER]',' total_norm=',tail_norm_global(1),&
        ' outer_shell_norm=',tail_norm_global(2),' outer_ratio=',&
        tail_norm_global(2)/max(1d-300,tail_norm_global(1)),&
        ' tolerance=',dg_wpw_scf_residual_tolerance
      tail_warning=.false.
      if(local_info==0)call classify_sawf_wannier_buffer_tail(tail_norm_global(1),&
        tail_norm_global(2),dg_wpw_scf_residual_tolerance,outer_shell_ratio,tail_warning,local_info)
      if(.not.wpw_seed_collective_stage_ok('buffer_tail_validity',local_info))return
      if(tail_warning.and.dc%id_tot==0)write(*,'(1x,a,2(a,es12.4))')&
        '[DG-WPW-SEED-BUFFER-WARNING] status=warning',' outer_ratio=',outer_shell_ratio,&
        ' tolerance=',dg_wpw_scf_residual_tolerance

      ! P already contains its periodic buffer in one unwrapped-contiguous
      ! ordering.  Do not canonicalize and add aliases back into core values.
      local_info=merge(0,1,all(ieee_is_finite(real(source_values))).and.&
        all(ieee_is_finite(aimag(source_values))))
      if(.not.wpw_seed_collective_stage_ok('unwrapped_source_values',local_info))return
      descriptor_buffer_lo=1-dc%nxyz_buffer+stencil_radius
      descriptor_buffer_hi=dc%nxyz_domain+dc%nxyz_buffer-stencil_radius
      descriptor_gradient_lo=descriptor_buffer_lo
      descriptor_gradient_hi=descriptor_buffer_hi
      allocate(local_buffer_point_ids(product(gradient_extent)),&
        descriptor_buffer_values(product(gradient_extent),source_count),&
        descriptor_buffer_gradients(3,product(gradient_extent),source_count))
      gradient_point=0
      do izb=1+stencil_radius,stencil_extent(3)-stencil_radius
        do iyb=1+stencil_radius,stencil_extent(2)-stencil_radius
          do ixb=1+stencil_radius,stencil_extent(1)-stencil_radius
        gradient_point=gradient_point+1
        buffer_point=ixb+(iyb-1)*stencil_extent(1)+(izb-1)*stencil_extent(1)*stencil_extent(2)
        descriptor_grid=[ixb,iyb,izb]-dc%nxyz_buffer
        local_buffer_point_ids(gradient_point)=descriptor_grid(1)-descriptor_buffer_lo(1)+1+&
          (descriptor_grid(2)-descriptor_buffer_lo(2))*(descriptor_buffer_hi(1)-descriptor_buffer_lo(1)+1)+&
          (descriptor_grid(3)-descriptor_buffer_lo(3))*(descriptor_buffer_hi(1)-descriptor_buffer_lo(1)+1)*&
            (descriptor_buffer_hi(2)-descriptor_buffer_lo(2)+1)
        descriptor_buffer_values(gradient_point,:)=buffer_values(buffer_point,:)
        descriptor_buffer_gradients(:,gradient_point,:)=compact_gradient_values(:,gradient_point,:)
      enddo;enddo;enddo
      point=0
      do izp=max(mg%is(3),1),min(mg%ie(3),dc%nxyz_domain(3))
        do iyp=max(mg%is(2),1),min(mg%ie(2),dc%nxyz_domain(2))
          do ixp=max(mg%is(1),1),min(mg%ie(1),dc%nxyz_domain(1))
            point=point+1;descriptor_grid=[ixp,iyp,izp]
            core_descriptor_point=descriptor_grid(1)-descriptor_buffer_lo(1)+1+&
              (descriptor_grid(2)-descriptor_buffer_lo(2))*gradient_extent(1)+&
              (descriptor_grid(3)-descriptor_buffer_lo(3))*gradient_extent(1)*gradient_extent(2)
            if(core_descriptor_point<1.or.core_descriptor_point>size(descriptor_buffer_values,1).or.&
              maxval(abs(source_core(point,:)-descriptor_buffer_values(core_descriptor_point,:)))>&
                1d-12*max(1d0,maxval(abs(source_core(point,:)))))local_info=1
          enddo
        enddo
      enddo
      if(point/=size(source_core,1))local_info=1
      if(.not.wpw_seed_collective_stage_ok('descriptor_core_d_invariant',local_info))return
      write(*,'(1x,a,3(a,i0),4(a,3(i0,1x)),2(a,l1))')'[DG-WPW-GRADIENT-WINDOW]',&
        ' rank=',dc%id_tot,' fragment=',dc%i_frag,' orbital_rank=',info%id_o,&
        ' mg_lo=',mg%is,' mg_hi=',mg%ie,' dg_lo=',descriptor_buffer_lo,&
        ' dg_hi=',descriptor_buffer_hi,' values_finite=',&
        all(ieee_is_finite(real(descriptor_buffer_values))).and.&
        all(ieee_is_finite(aimag(descriptor_buffer_values))),&
        ' gradients_finite=',all(ieee_is_finite(real(descriptor_buffer_gradients))).and.&
        all(ieee_is_finite(aimag(descriptor_buffer_gradients)))
      call gather_dg_wpw_occupied_w_payload(dc%icomm_frag,info%id_o==0,&
        wpw_core_p_accumulator%grid_ids(1:wpw_core_p_accumulator%npoint),&
        source_values(:,source_offset+1:source_offset+source_count),&
        local_buffer_point_ids,descriptor_buffer_values,descriptor_buffer_gradients,&
        product(dc%nxyz_domain),&
        product(descriptor_buffer_hi-descriptor_buffer_lo+1),gathered_core_grid_ids,&
        gathered_core_values,gathered_buffer_values,gathered_buffer_gradients,local_info)
      if(.not.wpw_seed_collective_stage_ok('occupied_w_payload_gather',local_info))return
      local_info=0
      if(dc%id_frag==0)then
        if(global_source_count==128.and.dc%n_frag==8.and.source_count/=16)local_info=1
        if(local_info==0)call initialize_dg_wpw_occupied_w_basis_collective(occupied_w_basis,&
          wpw_production_comm,dc%i_frag,source_key_local,gathered_core_grid_ids,gathered_core_values,&
          descriptor_buffer_lo,descriptor_buffer_hi,gathered_buffer_values,gathered_buffer_gradients,&
          source_condition,1,global_source_count,local_info,&
          descriptor_gradient_lo,descriptor_gradient_hi)
      endif
      if(.not.wpw_seed_collective_stage_ok('occupied_w_descriptor_commit',local_info))return
      call broadcast_dg_wpw_occupied_w_basis(dc%icomm_frag,0,occupied_w_basis,local_info)
      if(.not.wpw_seed_collective_stage_ok('occupied_w_descriptor_broadcast',local_info))return
      source_info=0
    end subroutine build_core_owned_projected_wannier_density_seed

    subroutine build_wpw_density_carrying_fragment_seed(seed_info)
      integer,intent(out)::seed_info
      integer::local_nsource,isource,point,ixp,iyp,izp,iw,offset,ierr,root_info,point_info,source_info
      integer::grid_point(3)
      integer::occ_counts(dc%n_frag),occ_counts_all(dc%n_frag)
      integer,allocatable::source_w_ids_local(:),source_w_ids(:)
      complex(8),allocatable::local_overlap_w(:,:),root_overlap_w(:,:),support_w_values(:),&
        support_w_gradients(:,:),owned_w_values(:),owned_w_gradients(:,:),&
        local_overlap_p(:,:),root_overlap_p(:,:),metric_rhs_w(:,:),metric_rhs_p(:,:),&
        metric_partial_w(:,:),metric_partial_p(:,:),rhs_all(:,:),q_all(:,:),capture(:,:),source_values(:,:),&
        local_w_norm(:,:),root_w_norm(:,:),zero_p_norm(:,:),owned_w_norm(:,:),owned_p_norm(:,:)
      complex(8),allocatable::raw_qw(:,:),raw_qp(:,:),normalization_transform(:,:),f_source(:,:),f_q(:,:),&
        density_rw_support(:,:),density_rp(:,:),density_values(:,:),metric_s_w(:,:),metric_s_p(:,:),&
        metric_zero_w(:,:),metric_zero_p(:,:)
      complex(8),allocatable::density_w_values(:,:),density_p_values(:,:)
      real(8),allocatable::local_source_norm(:),root_source_norm(:),metric_diagonal_w(:),&
        metric_diagonal_p(:),metric_rhs_residuals(:),occupations_local(:),source_norm_local_global(:),&
        source_norm_global(:)
      real(8),allocatable::metric_rhs_residual_history(:,:)
      real(8),allocatable::source_density(:),projected_density(:),normalized_density(:)
      real(8)::projection_norms_local(2),projection_norms_global(2),metric_diagonal_spread
      real(8)::projection_orth,source_condition,source_charge,projected_charge,normalized_charge,&
        density_projection_residual,density_normalization_residual,density_charge_error,&
        density_dc_local(2),density_dc_global(2),direct_cross_local,direct_cross_global,&
        direct_captured_norm,direct_w_local,direct_w_global,direct_p_local,direct_p_global,&
        routed_w_local,routed_w_global,routed_p_local,routed_p_global,capture_denominator,&
        metric_s_norm_local,metric_s_norm_global,gram_split_local(3),gram_split_global(3),&
        assembled_split_local(3),assembled_split_global(3),w_norm_stats_local(3),w_norm_stats_global(3)
      integer::projection_rank,metric_iterations,metric_history_iter,w_norm_count_local,w_norm_count_global,&
        source_w_position,local_source_position
      integer::w_norm_max_location(1)

      seed_info=1;root_info=0
      local_nsource=wpw_bootstrap_source_count;source_condition=wpw_bootstrap_source_condition
      source_info=wpw_bootstrap_source_info
      if(allocated(wpw_bootstrap_source_values))then
        allocate(source_values,source=wpw_bootstrap_source_values)
      else
        source_info=1
      endif
      if(.not.wpw_seed_collective_stage_ok('source_build',source_info))return
      occ_counts=0
      if(dc%id_frag==0)occ_counts(dc%i_frag)=local_nsource
      call comm_summation(occ_counts,occ_counts_all,dc%n_frag,dc%icomm_tot)
      wpw_nocc=occupied_index_from_input(1);wpw_nretain=wpw_nocc+dg_wpw_extra_states
      if(sum(occ_counts_all)/=wpw_nocc.or.local_nsource<1.or.&
         wpw_nretain>n_mat(1)+dc%n_frag*size(wpw_g_indices,2))root_info=1
      allocate(source_w_ids_local(wpw_nocc),source_w_ids(wpw_nocc));source_w_ids_local=0
      offset=sum(occ_counts_all(1:dc%i_frag-1))
      if(dc%id_frag==0)then
        do isource=1,local_nsource
          local_source_position=0
          do iw=1,wpw_occupied_w_basis%local_count
            if(all(wpw_occupied_w_basis%stable_keys(:,iw)==wpw_bootstrap_source_keys(:,isource)))then
              local_source_position=iw;exit
            endif
          enddo
          if(local_source_position==0)then
            root_info=1
          else
            source_w_ids_local(offset+isource)=wpw_occupied_w_basis%owned_ids(local_source_position)
          endif
        enddo
      endif
      call comm_summation(source_w_ids_local,source_w_ids,wpw_nocc,dc%icomm_tot)
      if(any(source_w_ids<=0))root_info=1
      do isource=1,wpw_nocc;do iw=isource+1,wpw_nocc
        if(source_w_ids(isource)==source_w_ids(iw))root_info=1
      enddo;enddo
      allocate(local_overlap_w(size(wpw_support_w_ids),wpw_nocc),&
        root_overlap_w(size(wpw_support_w_ids),wpw_nocc),&
        local_overlap_p(size(wpw_support_fragment_ids)*size(wpw_g_indices,2),wpw_nocc),&
        root_overlap_p(size(wpw_support_fragment_ids)*size(wpw_g_indices,2),wpw_nocc),&
        local_source_norm(wpw_nocc),root_source_norm(wpw_nocc),&
        local_w_norm(size(wpw_support_w_ids),3),root_w_norm(size(wpw_support_w_ids),3),&
        support_w_values(size(wpw_support_w_ids)),&
        support_w_gradients(3,size(wpw_support_w_ids)),owned_w_values(size(wpw_owned_w_ids)),&
        owned_w_gradients(3,size(wpw_owned_w_ids)))
      local_overlap_w=0;root_overlap_w=0;local_overlap_p=0;root_overlap_p=0;local_w_norm=0;root_w_norm=0
      local_source_norm=0;root_source_norm=0;point=0
      do izp=max(mg%is(3),1),min(mg%ie(3),size(f_basis,3))
        do iyp=max(mg%is(2),1),min(mg%ie(2),size(f_basis,2))
          do ixp=max(mg%is(1),1),min(mg%ie(1),size(f_basis,1))
            point=point+1
            owned_w_values=wpw_volume_accumulator%w_points(:,point)
            owned_w_gradients=wpw_volume_accumulator%grad_w_points(:,:,point)
            grid_point=dc%ixyz_frag(:,dc%i_frag)+[ixp,iyp,izp]-1
            call evaluate_dg_wpw_core_w_support(wpw_owned_w_ids,owned_w_values,owned_w_gradients,&
              wpw_support_w_ids,wpw_volume_halos,grid_point,1,support_w_values,support_w_gradients,&
              point_info,zero_outside_halo=.true.)
            if(point_info/=0)then
              if(root_info==0)write(*,'(1x,a,i0,a,3(i0,1x),a,i0)')&
                '[DG-WPW-LOCAL-FAIL] seed_w_support rank=',dc%id_tot,' grid=',grid_point,&
                ' info=',point_info
              root_info=1;cycle
            endif
            local_w_norm(:,1)=local_w_norm(:,1)+wpw_volume_accumulator%weights(point)*&
              abs(support_w_values)**2
            do iw=1,size(wpw_support_w_ids)
              if(find_integer_id(wpw_owned_w_ids,wpw_support_w_ids(iw))>0)then
                local_w_norm(iw,2)=local_w_norm(iw,2)+wpw_volume_accumulator%weights(point)*&
                  abs(support_w_values(iw))**2
              else
                local_w_norm(iw,3)=local_w_norm(iw,3)+wpw_volume_accumulator%weights(point)*&
                  abs(support_w_values(iw))**2
              endif
            enddo
            do isource=1,wpw_nocc
              if(info%id_o/=0)cycle
              source_w_position=find_integer_id(wpw_support_w_ids,source_w_ids(isource))
              if(source_w_position>0)then
                source_values(point,isource)=support_w_values(source_w_position)
              else
                source_values(point,isource)=(0d0,0d0)
              endif
              local_source_norm(isource)=local_source_norm(isource)+&
                abs(source_values(point,isource))**2*&
                wpw_volume_accumulator%weights(point)
              do iw=1,size(wpw_support_w_ids)
                local_overlap_w(iw,isource)=local_overlap_w(iw,isource)+conjg(support_w_values(iw))*&
                  source_values(point,isource)*&
                  wpw_volume_accumulator%weights(point)
              enddo
              do iw=1,size(local_overlap_p,1)
                local_overlap_p(iw,isource)=local_overlap_p(iw,isource)+conjg(&
                  wpw_volume_accumulator%p_points(iw,point))*&
                  source_values(point,isource)*&
                  wpw_volume_accumulator%weights(point)
              enddo
            enddo
          enddo
        enddo
      enddo
      if(point/=wpw_volume_accumulator%npoint)root_info=1
      if(.not.wpw_seed_collective_stage_ok('source_overlap_assembly',root_info))return
      call MPI_Reduce(local_overlap_w,root_overlap_w,size(local_overlap_w),MPI_DOUBLE_COMPLEX,MPI_SUM,0,&
        dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)root_info=1
      call MPI_Reduce(local_overlap_p,root_overlap_p,size(local_overlap_p),MPI_DOUBLE_COMPLEX,MPI_SUM,0,&
        dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)root_info=1
      call MPI_Reduce(local_source_norm,root_source_norm,size(local_source_norm),MPI_DOUBLE_PRECISION,&
        MPI_SUM,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)root_info=1
      call MPI_Reduce(local_w_norm,root_w_norm,size(local_w_norm),MPI_DOUBLE_COMPLEX,MPI_SUM,0,&
        dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)root_info=1
      allocate(wpw_qw(size(wpw_owned_w_ids),wpw_nretain),&
        wpw_qp(size(wpw_g_indices,2),wpw_nretain),wpw_q_old_occ(&
        size(wpw_owned_w_ids)+size(wpw_g_indices,2),wpw_nocc),wpw_eigenvalues(wpw_nretain),&
        wpw_occupations(wpw_nocc),raw_qw(size(wpw_owned_w_ids),wpw_nocc),&
        raw_qp(size(wpw_g_indices,2),wpw_nocc),normalization_transform(wpw_nocc,wpw_nocc),&
        f_source(wpw_nocc,wpw_nocc),f_q(wpw_nocc,wpw_nocc))
      wpw_qw=0;wpw_qp=0;wpw_q_old_occ=0;wpw_eigenvalues=0;wpw_occupations=0
      raw_qw=0;raw_qp=0;normalization_transform=0;f_source=0;f_q=0
      if(dc%id_frag==0)then
        call MPI_Allreduce(root_info,projection_rank,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr)
        if(ierr/=MPI_SUCCESS)then;root_info=1;else;root_info=projection_rank;endif
      endif
      if(dc%id_frag==0.and.root_info==0)then
        call initialize_dg_wpw_deterministic_seed(wpw_owned_w_ids,&
          wpw_context%bounded_operator%owned_p_ids,wpw_qw,wpw_qp,root_info)
        offset=sum(occ_counts_all(1:dc%i_frag-1))
        allocate(metric_rhs_w(size(wpw_qw,1),wpw_nocc),metric_rhs_p(size(wpw_qp,1),wpw_nocc),&
          metric_partial_w(size(root_overlap_w,1),wpw_nocc),&
          metric_partial_p(size(root_overlap_p,1),wpw_nocc),&
          metric_diagonal_w(size(wpw_qw,1)),metric_diagonal_p(size(wpw_qp,1)),&
          metric_rhs_residuals(wpw_nocc),metric_rhs_residual_history(256,wpw_nocc),&
          occupations_local(wpw_nocc),source_norm_local_global(wpw_nocc),source_norm_global(wpw_nocc))
        allocate(zero_p_norm(size(wpw_context%bounded_operator%required_p_ids),3),&
          owned_w_norm(size(wpw_qw,1),3),owned_p_norm(size(wpw_qp,1),3))
        metric_rhs_w=0;metric_rhs_p=0;metric_partial_w=0;metric_partial_p=0
        metric_iterations=0;wpw_projection_residual=huge(1d0)
        occupations_local=0;source_norm_local_global=0;source_norm_global=0
        metric_partial_w=root_overlap_w
        metric_partial_p=root_overlap_p
        zero_p_norm=0;owned_w_norm=0;owned_p_norm=0
        call reduce_dg_wpw_metric_rhs_partials(wpw_context%bounded_operator,root_w_norm,zero_p_norm,&
          owned_w_norm,owned_p_norm,root_info)
        w_norm_stats_local=[minval(real(owned_w_norm(:,1),8)),maxval(real(owned_w_norm(:,1),8)),&
          sum(real(owned_w_norm(:,1),8))]
        w_norm_count_local=size(owned_w_norm,1)
        call MPI_Allreduce(w_norm_stats_local(1),w_norm_stats_global(1),1,MPI_DOUBLE_PRECISION,MPI_MIN,&
          wpw_production_comm,ierr)
        call MPI_Allreduce(w_norm_stats_local(2),w_norm_stats_global(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,&
          wpw_production_comm,ierr)
        call MPI_Allreduce(w_norm_stats_local(3),w_norm_stats_global(3),1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          wpw_production_comm,ierr)
        call MPI_Allreduce(w_norm_count_local,w_norm_count_global,1,MPI_INTEGER,MPI_SUM,&
          wpw_production_comm,ierr)
        if(ierr/=MPI_SUCCESS.or.root_info/=0.or.w_norm_count_global<=0.or.&
          .not.all(ieee_is_finite(w_norm_stats_global)).or.w_norm_stats_global(1)<=0d0)root_info=1
        if(maxval(abs(real(owned_w_norm(:,1)-owned_w_norm(:,2)-owned_w_norm(:,3),8)))>&
          1d-12*max(1d0,maxval(real(owned_w_norm(:,1),8))))root_info=1
        if(root_info==0.and.wpw_context%rank_id==0)write(*,'(1x,a,3(a,es12.4),a,i0)')&
          '[DG-WPW-W-REALSPACE-NORM]',' min=',w_norm_stats_global(1),&
          ' max=',w_norm_stats_global(2),' mean=',w_norm_stats_global(3)/dble(w_norm_count_global),&
          ' count=',w_norm_count_global
        if(root_info==0)then
          w_norm_max_location=maxloc(real(owned_w_norm(:,1),8))
          write(*,'(1x,a,3(a,i0),3(a,es12.4))')'[DG-WPW-W-NORM-LOCAL-MAX]',&
            ' rank=',wpw_context%rank_id,' fragment=',dc%i_frag,&
            ' w_id=',wpw_owned_w_ids(w_norm_max_location(1)),&
            ' total=',real(owned_w_norm(w_norm_max_location(1),1),8),&
            ' owner_core=',real(owned_w_norm(w_norm_max_location(1),2),8),&
            ' halo_tail=',real(owned_w_norm(w_norm_max_location(1),3),8)
        endif
        occupations_local(offset+1:offset+local_nsource)=2d0
        source_norm_local_global=root_source_norm
        call reduce_dg_wpw_metric_rhs_partials(wpw_context%bounded_operator,metric_partial_w,&
          metric_partial_p,metric_rhs_w,metric_rhs_p,root_info)
        if(root_info/=0)write(*,'(1x,a,i0)')'[DG-WPW-LOCAL-FAIL] seed_metric_rhs rank=',dc%id_tot
        projection_norms_local=[sum(abs(root_overlap_w)**2)+sum(abs(root_overlap_p)**2),&
          sum(root_source_norm)]
        call MPI_Allreduce(projection_norms_local,projection_norms_global,2,MPI_DOUBLE_PRECISION,MPI_SUM,&
          wpw_production_comm,ierr)
        if(ierr/=MPI_SUCCESS.or..not.all(ieee_is_finite(projection_norms_global)).or.&
          projection_norms_global(2)<=0d0)then
          write(*,'(1x,a,i0,a,i0,2(a,es12.4))')'[DG-WPW-LOCAL-FAIL] seed_projection_norm rank=',&
            dc%id_tot,' mpi=',ierr,' rhs_norm=',projection_norms_global(1),&
            ' source_norm=',projection_norms_global(2)
          root_info=1
        endif
        call MPI_Allreduce(occupations_local,wpw_occupations,wpw_nocc,MPI_DOUBLE_PRECISION,MPI_SUM,&
          wpw_production_comm,ierr)
        if(ierr/=MPI_SUCCESS.or.any(wpw_occupations<=0d0))then
          write(*,'(1x,a,i0,a,i0,a,2(es12.4,1x))')'[DG-WPW-LOCAL-FAIL] seed_occupations rank=',&
            dc%id_tot,' mpi=',ierr,' minmax=',minval(wpw_occupations),maxval(wpw_occupations)
          root_info=1
        endif
        call MPI_Allreduce(source_norm_local_global,source_norm_global,wpw_nocc,MPI_DOUBLE_PRECISION,&
          MPI_SUM,wpw_production_comm,ierr)
        if(ierr/=MPI_SUCCESS.or.any(source_norm_global<=0d0))then
          write(*,'(1x,a,i0,a,i0,a,2(es12.4,1x))')'[DG-WPW-LOCAL-FAIL] seed_source_norm rank=',&
            dc%id_tot,' mpi=',ierr,' minmax=',minval(source_norm_global),maxval(source_norm_global)
          root_info=1
        endif
        if(root_info==0)call get_dg_wpw_owned_metric_diagonal(wpw_context%bounded_operator,&
          metric_diagonal_w,metric_diagonal_p,root_info)
        if(root_info/=0)write(*,'(1x,a,i0)')'[DG-WPW-LOCAL-FAIL] seed_metric_diagonal rank=',dc%id_tot
        if(root_info==0)call initialize_dg_wpw_fragment_block_preconditioner(&
          wpw_context%bounded_operator,dg_wpw_metric_cutoff,&
          wpw_context%metric_block_preconditioner,root_info)
        if(root_info==0)write(*,'(1x,a,i0,4(a,es12.4))')'[DG-WPW-METRIC-BLOCK] dimension=',&
          wpw_context%metric_block_preconditioner%dimension,&
          ' hermitian_defect=',wpw_context%metric_block_preconditioner%hermitian_defect,&
          ' min=',wpw_context%metric_block_preconditioner%minimum_eigenvalue,&
          ' max=',wpw_context%metric_block_preconditioner%maximum_eigenvalue,&
          ' condition=',wpw_context%metric_block_preconditioner%condition_number
        if(root_info/=0)write(*,'(1x,a,i0)')'[DG-WPW-LOCAL-FAIL] seed_metric_block rank=',dc%id_tot
        if(root_info==0)call solve_dg_wpw_metric_projection(wpw_context,&
          wpw_production_comm,wpw_apply_s,wpw_global_gram,size(wpw_qw,1),size(wpw_qp,1),&
          wpw_nocc,dg_wpw_metric_cutoff,256,metric_diagonal_w,metric_diagonal_p,metric_rhs_w,metric_rhs_p,&
          wpw_qw(:,1:wpw_nocc),wpw_qp(:,1:wpw_nocc),wpw_projection_residual,&
          metric_rhs_residuals,metric_rhs_residual_history,metric_iterations,metric_diagonal_spread,&
          root_info,stagnation_limit=257,diagnose_recurrence=.true.,&
          apply_metric_preconditioner=wpw_apply_metric_block_preconditioner,&
          allow_diagnostic_continuation=.true.,diagnostic_continuation=wpw_metric_diagnostic_only)
        if(root_info==0.and.wpw_metric_diagnostic_only.and.wpw_context%rank_id==0)then
          write(*,'(1x,a,3(a,es12.4))')'[DG-WPW-METRIC-DIAGNOSTIC-CONTINUATION]',&
            ' strict_target=',dg_wpw_metric_cutoff,' aggregate_residual=',wpw_projection_residual,&
            ' worst_rhs=',maxval(metric_rhs_residuals)
        endif
        if(root_info/=0)write(*,'(1x,a,i0,a,i0,a,es12.4)')&
          '[DG-WPW-LOCAL-FAIL] seed_metric_solve rank=',dc%id_tot,' iterations=',metric_iterations,&
          ' residual=',wpw_projection_residual
        if(root_info/=0.and.wpw_context%rank_id==0)then
          write(*,'(1x,a,es12.4)')'[DG-WPW-METRIC-DIAG] diagonal_spread=',metric_diagonal_spread
          do metric_history_iter=1,metric_iterations
            write(*,'(1x,a,i0,3(a,es12.4),a,i0)')'[DG-WPW-METRIC-HISTORY] iter=',metric_history_iter,&
              ' min=',minval(metric_rhs_residual_history(metric_history_iter,:)),&
              ' max=',maxval(metric_rhs_residual_history(metric_history_iter,:)),&
              ' mean=',sum(metric_rhs_residual_history(metric_history_iter,:))/dble(wpw_nocc),&
              ' worst_rhs=',maxloc(metric_rhs_residual_history(metric_history_iter,:),dim=1)
          enddo
        endif
        if(root_info==0)then
          raw_qw=wpw_qw(:,1:wpw_nocc);raw_qp=wpw_qp(:,1:wpw_nocc)
          allocate(metric_s_w(size(raw_qw,1),wpw_nocc),metric_s_p(size(raw_qp,1),wpw_nocc),&
            metric_zero_w(size(raw_qw,1),wpw_nocc),metric_zero_p(size(raw_qp,1),wpw_nocc))
          metric_s_w=0;metric_s_p=0;metric_zero_w=0;metric_zero_p=0
          call wpw_apply_s(wpw_context,raw_qw,raw_qp,metric_s_w,metric_s_p,source_info)
          call MPI_Allreduce(source_info,projection_rank,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr)
          if(ierr/=MPI_SUCCESS.or.projection_rank/=0)root_info=1
          metric_s_norm_local=0d0
          if(root_info==0)then
            do iw=1,wpw_nocc
              metric_s_norm_local=metric_s_norm_local+wpw_occupations(iw)*real(&
                sum(conjg(raw_qw(:,iw))*metric_s_w(:,iw))+&
                sum(conjg(raw_qp(:,iw))*metric_s_p(:,iw)),8)
            enddo
          endif
          call MPI_Allreduce(metric_s_norm_local,metric_s_norm_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
            wpw_production_comm,ierr)
          if(ierr/=MPI_SUCCESS.or..not.ieee_is_finite(metric_s_norm_global))root_info=1
          assembled_split_local=0d0
          if(root_info==0)call wpw_apply_s(wpw_context,raw_qw,metric_zero_p,metric_s_w,metric_s_p,source_info)
          call MPI_Allreduce(source_info,projection_rank,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr)
          if(ierr/=MPI_SUCCESS.or.projection_rank/=0)root_info=1
          if(root_info==0)then
            do iw=1,wpw_nocc
              assembled_split_local(1)=assembled_split_local(1)+wpw_occupations(iw)*&
                real(sum(conjg(raw_qw(:,iw))*metric_s_w(:,iw)),8)
              assembled_split_local(2)=assembled_split_local(2)+wpw_occupations(iw)*&
                real(sum(conjg(raw_qp(:,iw))*metric_s_p(:,iw)),8)
            enddo
          endif
          if(root_info==0)call wpw_apply_s(wpw_context,metric_zero_w,raw_qp,metric_s_w,metric_s_p,source_info)
          call MPI_Allreduce(source_info,projection_rank,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr)
          if(ierr/=MPI_SUCCESS.or.projection_rank/=0)root_info=1
          if(root_info==0)then
            do iw=1,wpw_nocc
              assembled_split_local(2)=assembled_split_local(2)+wpw_occupations(iw)*&
                real(sum(conjg(raw_qw(:,iw))*metric_s_w(:,iw)),8)
              assembled_split_local(3)=assembled_split_local(3)+wpw_occupations(iw)*&
                real(sum(conjg(raw_qp(:,iw))*metric_s_p(:,iw)),8)
            enddo
          endif
          call MPI_Allreduce(assembled_split_local,assembled_split_global,3,MPI_DOUBLE_PRECISION,MPI_SUM,&
            wpw_production_comm,ierr)
          if(ierr/=MPI_SUCCESS.or..not.all(ieee_is_finite(assembled_split_global)))root_info=1
          do iw=1,wpw_nocc;f_source(iw,iw)=wpw_occupations(iw);enddo
          allocate(rhs_all(size(wpw_qw,1)+size(wpw_qp,1),wpw_nocc),&
            q_all(size(wpw_qw,1)+size(wpw_qp,1),wpw_nocc),capture(wpw_nocc,wpw_nocc))
          rhs_all(1:size(wpw_qw,1),:)=metric_rhs_w
          rhs_all(size(wpw_qw,1)+1:,:)=metric_rhs_p
          q_all(1:size(wpw_qw,1),:)=wpw_qw(:,1:wpw_nocc)
          q_all(size(wpw_qw,1)+1:,:)=wpw_qp(:,1:wpw_nocc)
          call wpw_global_gram(rhs_all,q_all,size(rhs_all,1),wpw_nocc,wpw_nocc,capture,root_info)
          wpw_projection_captured_norm=sum([(wpw_occupations(iw)*real(capture(iw,iw),8),iw=1,wpw_nocc)])/&
            sum(wpw_occupations*source_norm_global)
          routed_w_local=0d0;routed_p_local=0d0
          do iw=1,wpw_nocc
            routed_w_local=routed_w_local+wpw_occupations(iw)*real(sum(&
              conjg(metric_rhs_w(:,iw))*raw_qw(:,iw)),8)
            routed_p_local=routed_p_local+wpw_occupations(iw)*real(sum(&
              conjg(metric_rhs_p(:,iw))*raw_qp(:,iw)),8)
          enddo
          call MPI_Allreduce(routed_w_local,routed_w_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
            wpw_production_comm,ierr)
          call MPI_Allreduce(routed_p_local,routed_p_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
            wpw_production_comm,ierr)
          if(.not.ieee_is_finite(wpw_projection_captured_norm).or.&
             wpw_projection_captured_norm<=0d0)root_info=1
        endif
        if(root_info==0)call initialize_dg_wpw_projected_occupied(wpw_context,wpw_production_comm,&
          wpw_apply_s,wpw_global_gram,size(wpw_qw,1),size(wpw_qp,1),wpw_nocc,dg_wpw_metric_cutoff,&
          wpw_qw(:,1:wpw_nocc),wpw_qp(:,1:wpw_nocc),projection_rank,projection_orth,root_info,&
          normalization_transform)
        if(root_info==0)call transform_sawf_wannier_occupation(normalization_transform,f_source,f_q,root_info)
        if(root_info==0.and.(projection_rank/=wpw_nocc.or.&
          projection_orth>100d0*dg_wpw_metric_cutoff))root_info=1
        if(root_info==0)call complete_dg_wpw_projected_subspace(wpw_context,wpw_production_comm,&
          wpw_apply_s,wpw_global_gram,size(wpw_qw,1),size(wpw_qp,1),wpw_nocc,wpw_nretain,&
          dg_wpw_metric_cutoff,wpw_qw,wpw_qp,projection_rank,projection_orth,root_info)
        if(root_info==0.and.(projection_rank/=wpw_nretain.or.&
          projection_orth>100d0*dg_wpw_metric_cutoff))root_info=1
        wpw_projection_rank=projection_rank;wpw_projection_orthogonality=projection_orth
        wpw_projection_charge=sum(wpw_occupations)
        if(wpw_context%rank_id==0.and.root_info==0)write(*,'(1x,a,i0,4(a,es12.4))')&
          '[DG-WPW-PROJECTION] source=density_carrying_fragment_seed rank=',wpw_projection_rank,&
          ' metric_residual=',wpw_projection_residual,' captured_norm=',wpw_projection_captured_norm,&
          ' s_orth=',wpw_projection_orthogonality,' charge=',wpw_projection_charge
      endif
      call MPI_Bcast(root_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS.or.root_info/=0)return
      call MPI_Bcast(wpw_metric_diagnostic_only,1,MPI_LOGICAL,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(wpw_qw,size(wpw_qw),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(wpw_qp,size(wpw_qp),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(wpw_occupations,size(wpw_occupations),MPI_DOUBLE_PRECISION,0,&
        dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(raw_qw,size(raw_qw),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(raw_qp,size(raw_qp),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(f_source,size(f_source),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(f_q,size(f_q),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)return

      allocate(density_rw_support(size(wpw_support_w_ids),wpw_nocc),&
        density_rp(size(wpw_support_fragment_ids)*size(wpw_g_indices,2),wpw_nocc),&
        density_values(wpw_volume_accumulator%npoint,wpw_nocc),&
        density_w_values(wpw_volume_accumulator%npoint,wpw_nocc),&
        density_p_values(wpw_volume_accumulator%npoint,wpw_nocc),&
        source_density(wpw_volume_accumulator%npoint),projected_density(wpw_volume_accumulator%npoint),&
        normalized_density(wpw_volume_accumulator%npoint))
      density_rw_support=0;density_rp=0;source_info=0
      if(dc%id_frag==0)call fetch_dg_wpw_support_coefficients(wpw_context%bounded_operator,&
        raw_qw,raw_qp,density_rw_support,density_rp,source_info)
      call MPI_Bcast(source_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS.or.source_info/=0)return
      call MPI_Bcast(density_rw_support,size(density_rw_support),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(density_rp,size(density_rp),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)return
      density_w_values=0;point=0
      do izp=max(mg%is(3),1),min(mg%ie(3),size(f_basis,3))
        do iyp=max(mg%is(2),1),min(mg%ie(2),size(f_basis,2))
          do ixp=max(mg%is(1),1),min(mg%ie(1),size(f_basis,1))
            point=point+1
            owned_w_values=wpw_volume_accumulator%w_points(:,point)
            owned_w_gradients=wpw_volume_accumulator%grad_w_points(:,:,point)
            grid_point=dc%ixyz_frag(:,dc%i_frag)+[ixp,iyp,izp]-1
            call reconstruct_dg_wpw_core_w_support(wpw_owned_w_ids,owned_w_values,owned_w_gradients,&
              wpw_support_w_ids,wpw_volume_halos,grid_point,1,density_rw_support,support_w_values,&
              support_w_gradients,density_w_values(point,:),point_info,zero_outside_halo=.true.)
            if(point_info/=0)source_info=1
          enddo
        enddo
      enddo
      if(point/=wpw_volume_accumulator%npoint)source_info=1
      if(.not.wpw_seed_collective_stage_ok('raw_support_w_reconstruction',source_info))return
      density_p_values=matmul(transpose(wpw_volume_accumulator%p_points),density_rp)
      density_values=density_w_values+density_p_values
      direct_cross_local=0d0;direct_w_local=0d0;direct_p_local=0d0
      gram_split_local=0d0
      if(info%id_o==0)then
        do iw=1,wpw_nocc
          direct_w_local=direct_w_local+wpw_occupations(iw)*real(sum(&
            wpw_volume_accumulator%weights*conjg(source_values(:,iw))*density_w_values(:,iw)),8)
          direct_p_local=direct_p_local+wpw_occupations(iw)*real(sum(&
            wpw_volume_accumulator%weights*conjg(source_values(:,iw))*density_p_values(:,iw)),8)
          gram_split_local(1)=gram_split_local(1)+wpw_occupations(iw)*sum(&
            wpw_volume_accumulator%weights*abs(density_w_values(:,iw))**2)
          gram_split_local(2)=gram_split_local(2)+wpw_occupations(iw)*sum(&
            wpw_volume_accumulator%weights*2d0*real(conjg(density_w_values(:,iw))*density_p_values(:,iw),8))
          gram_split_local(3)=gram_split_local(3)+wpw_occupations(iw)*sum(&
            wpw_volume_accumulator%weights*abs(density_p_values(:,iw))**2)
        enddo
        direct_cross_local=direct_w_local+direct_p_local
      endif
      call MPI_Allreduce(direct_cross_local,direct_cross_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        dc%icomm_tot,ierr)
      call MPI_Allreduce(direct_w_local,direct_w_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
      call MPI_Allreduce(direct_p_local,direct_p_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
      call MPI_Allreduce(gram_split_local,gram_split_global,3,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(dc%id_tot==0)then
        capture_denominator=sum(wpw_occupations*source_norm_global)
        direct_captured_norm=direct_cross_global/capture_denominator
        write(*,'(1x,a,7(a,es12.4))')'[DG-WPW-RHS-REALSPACE-IDENTITY]',&
          ' routed_B_capture=',wpw_projection_captured_norm,&
          ' direct_A_capture=',direct_captured_norm,&
          ' routed_W=',routed_w_global/capture_denominator,&
          ' direct_W=',direct_w_global/capture_denominator,&
          ' routed_P=',routed_p_global/capture_denominator,&
          ' direct_P=',direct_p_global/capture_denominator,&
          ' relative_defect=',abs(direct_captured_norm-wpw_projection_captured_norm)/&
            max(1d-300,abs(wpw_projection_captured_norm))
        write(*,'(1x,a,5(a,es12.4))')'[DG-WPW-METRIC-REALSPACE-GRAM]',&
          ' assembled_S=',metric_s_norm_global/capture_denominator,&
          ' realspace_total=',sum(gram_split_global)/capture_denominator,&
          ' WW=',gram_split_global(1)/capture_denominator,&
          ' WP=',gram_split_global(2)/capture_denominator,&
          ' PP=',gram_split_global(3)/capture_denominator
        write(*,'(1x,a,4(a,es12.4))')'[DG-WPW-ASSEMBLED-METRIC-SPLIT]',&
          ' total=',sum(assembled_split_global)/capture_denominator,&
          ' WW=',assembled_split_global(1)/capture_denominator,&
          ' WP=',assembled_split_global(2)/capture_denominator,&
          ' PP=',assembled_split_global(3)/capture_denominator
      endif
      source_density=0;projected_density=0;normalized_density=0
      if(info%id_o==0)then
        call build_sawf_wannier_density(source_values,f_source,wpw_volume_accumulator%weights,&
          source_density,source_charge,source_info)
        if(source_info==0)call build_sawf_wannier_density(density_values,f_source,&
          wpw_volume_accumulator%weights,projected_density,projected_charge,source_info)
      endif
      call MPI_Allreduce(source_info,root_info,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or.root_info/=0)return

      density_rw_support=0;density_rp=0;source_info=0
      if(dc%id_frag==0)call fetch_dg_wpw_support_coefficients(wpw_context%bounded_operator,&
        wpw_qw(:,1:wpw_nocc),wpw_qp(:,1:wpw_nocc),density_rw_support,density_rp,source_info)
      call MPI_Bcast(source_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS.or.source_info/=0)return
      call MPI_Bcast(density_rw_support,size(density_rw_support),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(density_rp,size(density_rp),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)return
      density_w_values=0;point=0
      do izp=max(mg%is(3),1),min(mg%ie(3),size(f_basis,3))
        do iyp=max(mg%is(2),1),min(mg%ie(2),size(f_basis,2))
          do ixp=max(mg%is(1),1),min(mg%ie(1),size(f_basis,1))
            point=point+1
            owned_w_values=wpw_volume_accumulator%w_points(:,point)
            owned_w_gradients=wpw_volume_accumulator%grad_w_points(:,:,point)
            grid_point=dc%ixyz_frag(:,dc%i_frag)+[ixp,iyp,izp]-1
            call reconstruct_dg_wpw_core_w_support(wpw_owned_w_ids,owned_w_values,owned_w_gradients,&
              wpw_support_w_ids,wpw_volume_halos,grid_point,1,density_rw_support,support_w_values,&
              support_w_gradients,density_w_values(point,:),point_info,zero_outside_halo=.true.)
            if(point_info/=0)source_info=1
          enddo
        enddo
      enddo
      if(point/=wpw_volume_accumulator%npoint)source_info=1
      if(.not.wpw_seed_collective_stage_ok('normalized_support_w_reconstruction',source_info))return
      density_values=density_w_values+matmul(transpose(wpw_volume_accumulator%p_points),density_rp)
      source_info=0
      if(info%id_o==0)then
        call build_sawf_wannier_density(density_values,f_q,wpw_volume_accumulator%weights,&
          normalized_density,normalized_charge,source_info)
        if(source_info==0)call qualify_sawf_wannier_density_projection(source_density,projected_density,&
          normalized_density,wpw_volume_accumulator%weights,source_charge,wpw_projection_captured_norm,&
          dg_wpw_scf_residual_tolerance,density_projection_residual,density_normalization_residual,&
          density_charge_error,source_info)
        density_dc_local=[sum(wpw_volume_accumulator%weights*(source_density-&
          wpw_volume_accumulator%densities(1:size(source_density)))**2),&
          sum(wpw_volume_accumulator%weights*wpw_volume_accumulator%densities(1:size(source_density))**2)]
      else
        density_dc_local=0d0
      endif
      call MPI_Allreduce(density_dc_local,density_dc_global,2,MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
      if(dc%id_tot==0)write(*,'(1x,a,8(a,es12.4),a,i0)')'[DG-WPW-PHYSICAL-DIAGNOSTIC]',&
        ' captured_norm=',wpw_projection_captured_norm,&
        ' projection_residual=',density_projection_residual,&
        ' normalization_residual=',density_normalization_residual,&
        ' projected_charge_error=',density_charge_error,&
        ' source_charge=',source_charge,' projected_charge=',projected_charge,&
        ' normalized_charge=',normalized_charge,' dc_density_residual=',&
        sqrt(max(0d0,density_dc_global(1))/max(1d-300,density_dc_global(2))),&
        ' local_info=',source_info
      if(ierr/=MPI_SUCCESS.or.density_dc_global(2)<=0d0.or.&
        sqrt(density_dc_global(1)/density_dc_global(2))>dg_wpw_scf_residual_tolerance)source_info=1
      call MPI_Allreduce(source_info,root_info,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or.root_info/=0)return
      wpw_q_old_occ(1:size(wpw_qw,1),:)=wpw_qw(:,1:wpw_nocc)
      wpw_q_old_occ(size(wpw_qw,1)+1:,:)=wpw_qp(:,1:wpw_nocc)
      seed_info=0
    end subroutine build_wpw_density_carrying_fragment_seed

    subroutine run_wpw_production_scf(scf_info)
      integer,intent(out)::scf_info
      integer::i,iter,state_info,nw_local,np_local
      real(8),allocatable::density_initial(:)
      real(8)::remaining_electrons,scf_stage_clock

      scf_info=1;call cpu_time(scf_stage_clock)
      if(nspin/=1.or.wpw_volume_accumulator%npoint<=0)return
      wpw_nocc=occupied_index_from_input(1)
      wpw_nretain=wpw_nocc+dg_wpw_extra_states
      if(wpw_nocc<1.or.wpw_nretain>n_mat(1)+dc%n_frag*size(wpw_g_indices,2))return
      nw_local=size(wpw_owned_w_ids);np_local=size(wpw_g_indices,2)
      call trace_wpw_scf_stage('dimensions',scf_stage_clock,0,nw_local,np_local,wpw_nretain)
      allocate(wpw_occupations(wpw_nocc),wpw_eigenvalues(wpw_nretain))
      remaining_electrons=dc%elec_num_tot
      do i=1,wpw_nocc
        wpw_occupations(i)=min(2d0,remaining_electrons)
        remaining_electrons=max(0d0,remaining_electrons-wpw_occupations(i))
      enddo
      if(remaining_electrons>1d-10.or.any(wpw_occupations<=0d0))return
      wpw_eigenvalues=0d0
      allocate(wpw_qw(nw_local,wpw_nretain),wpw_qp(np_local,wpw_nretain),&
        wpw_q_old_occ(nw_local+np_local,wpw_nocc))
      wpw_qw=(0d0,0d0);wpw_qp=(0d0,0d0);wpw_q_old_occ=(0d0,0d0)
      state_info=0
      if(dc%id_frag==0)then
        call initialize_dg_wpw_deterministic_seed(wpw_owned_w_ids,&
          wpw_context%bounded_operator%owned_p_ids,wpw_qw,wpw_qp,state_info)
      endif
      call MPI_Bcast(state_info,1,MPI_INTEGER,0,dc%icomm_frag,i)
      if(i/=MPI_SUCCESS.or.state_info/=0)return
      call trace_wpw_scf_stage('seed',scf_stage_clock,0,nw_local,np_local,wpw_nretain)
      allocate(density_initial(wpw_volume_accumulator%npoint))
      density_initial=wpw_volume_accumulator%densities(1:wpw_volume_accumulator%npoint)
      call initialize_dg_wpw_scf_iteration_state(wpw_scf_state,dc%icomm_frag,0,nw_local,np_local,&
        wpw_nocc,density_initial,state_info)
      if(state_info/=0)then
        write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] scf_initialize rank=',dc%id_tot,' info=',state_info
        return
      endif
      call trace_wpw_scf_stage('state',scf_stage_clock,0,nw_local,np_local,wpw_nretain)
      do iter=1,dg_wpw_scf_max_iter
        call trace_wpw_scf_stage('algebra_begin',scf_stage_clock,iter,nw_local,np_local,wpw_nretain)
        call advance_dg_wpw_scf_iteration(wpw_scf_state,wpw_algebra_step,wpw_potential_step,state_info)
        call trace_wpw_scf_stage('algebra_end',scf_stage_clock,iter,nw_local,np_local,wpw_nretain)
        if(state_info/=0)then
          write(*,'(1x,a,i0,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] scf_advance rank=',dc%id_tot,&
            ' iter=',iter,' info=',state_info
          return
        endif
        if(.not.wpw_scf_state%algebra_ready)cycle
        call snapshot_wpw_fixed_point_state(state_info)
        if(state_info/=0)then
          write(*,'(1x,a,i0,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] scf_snapshot rank=',dc%id_tot,&
            ' iter=',iter,' info=',state_info
          return
        endif
        wpw_fixed_point_mode=.true.
        call verify_dg_wpw_scf_fixed_point(wpw_scf_state,dg_wpw_gap_threshold,&
          dg_wpw_scf_residual_tolerance,wpw_potential_step,state_info)
        wpw_fixed_point_mode=.false.
        if(state_info==0.and.wpw_scf_state%converged)then
          call publish_wpw_production_checkpoint(state_info)
          if(state_info/=0)return
          scf_info=0;return
        endif
        call restore_wpw_fixed_point_state(state_info)
        if(state_info/=0)then
          write(*,'(1x,a,i0,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] scf_restore rank=',dc%id_tot,&
            ' iter=',iter,' info=',state_info
          return
        endif
      enddo
      write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] scf_max_iter rank=',dc%id_tot,&
        ' iterations=',dg_wpw_scf_max_iter
    end subroutine run_wpw_production_scf

    subroutine trace_wpw_scf_stage(stage,stage_clock,iteration,nw,np,nretain)
      character(*),intent(in)::stage
      real(8),intent(inout)::stage_clock
      integer,intent(in)::iteration,nw,np,nretain
      real(8)::now
      call cpu_time(now)
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,f10.3)')'[DG-WPW-SCF-STAGE] stage=',trim(stage),&
        ' rank=',dc%id_tot,' iter=',iteration,' nw=',nw,' np=',np,' nretain=',nretain,&
        ' cpu_seconds=',now-stage_clock
      flush(6);stage_clock=now
    end subroutine trace_wpw_scf_stage

    subroutine publish_wpw_production_checkpoint(checkpoint_info)
      use salmon_global,only:sysname
      type(s_dg_wpw_checkpoint_state)::checkpoint,verified_checkpoint
      real(8),allocatable::retained_occupations(:)
      complex(8),allocatable::ww_z_sparse(:),wp_z_sparse(:),pp_z_dense(:,:)
      character(1024)::checkpoint_path
      integer,intent(out)::checkpoint_info
      integer::ierr,local_bad,global_bad,manifest_unit,manifest_ios
      logical::manifest_exists,publication_allowed
      integer(8)::local_checksum
      integer(8),allocatable::checksums(:)

      checkpoint_info=1;local_bad=0
      publication_allowed=.not.wpw_metric_diagnostic_only
      if(.not.publication_allowed)then
        if(dc%id_tot==0)write(*,'(1x,a)')'[DG-WPW-CHECKPOINT-BLOCKED] diagnostic metric continuation'
        return
      endif
      if(yn_dg_wpw_fixed_h_relaxation=='y')then
        call validate_wpw_frozen_h_state(checkpoint_info)
        if(checkpoint_info/=0)return
        local_bad=0
        if(dc%id_frag==0)then
          if(wpw_context%bounded_operator%interface_lambda/=1d0.or.&
            wpw_projection_rank<wpw_nocc.or.wpw_projection_residual<0d0.or.&
            .not.ieee_is_finite(wpw_projection_residual).or.&
            .not.ieee_is_finite(wpw_projection_captured_norm).or.wpw_projection_captured_norm<=0d0.or.&
            .not.wpw_fixed_h_qualified.or.&
            max(wpw_fixed_h_final_residual,wpw_fixed_h_final_orthogonality)>&
              dg_wpw_scf_residual_tolerance.or..not.ieee_is_finite(wpw_fixed_h_final_projector).or.&
            .not.ieee_is_finite(wpw_fixed_h_density_baseline).or.wpw_fixed_h_density_baseline>1d0.or.&
            wpw_fixed_h_charge_error>dg_wpw_scf_residual_tolerance) local_bad=1
        endif
        checkpoint_info=local_bad
        if(.not.wpw_potential_stage_ok(checkpoint_info))return
      endif
      local_bad=0
      if(dc%id_frag==0)then
        if(wpw_context%rank_id==0)then
          checkpoint_path=trim(import_run_root_dir())//trim(sysname)//'.dg_wpw.manifest'
          inquire(file=trim(checkpoint_path),exist=manifest_exists)
          if(manifest_exists)then
            open(newunit=manifest_unit,file=trim(checkpoint_path),status='old',iostat=manifest_ios)
            if(manifest_ios==0)close(manifest_unit,status='delete',iostat=manifest_ios)
            if(manifest_ios/=0)local_bad=1
          endif
        endif
      endif
      checkpoint_info=local_bad
      if(.not.wpw_potential_stage_ok(checkpoint_info))return
      call build_wpw_checkpoint_position_volume(ww_z_sparse,wp_z_sparse,pp_z_dense,checkpoint_info)
      if(checkpoint_info/=0)return
      if(dc%id_frag==0)then
        allocate(retained_occupations(wpw_nretain));retained_occupations=0d0
        retained_occupations(1:wpw_nocc)=wpw_occupations
        checkpoint%schema_version=3
        checkpoint%operator_epoch=wpw_context%bounded_operator%operator_epoch
        checkpoint%layout_fingerprint=wpw_context%bounded_operator%layout_fingerprint
        checkpoint%ownership_kind=wpw_context%bounded_operator%ownership_kind
        checkpoint%metric_convention=wpw_context%bounded_operator%metric_convention
        checkpoint%operator_convention=wpw_context%bounded_operator%operator_convention
        checkpoint%n_occ=wpw_nocc
        if(yn_dg_wpw_fixed_h_relaxation=='y')then
          checkpoint%fixed_h_mode=1
          checkpoint%seed_provenance='density_carrying_fragment_seed'
          checkpoint%metric_residual=wpw_projection_residual
          checkpoint%captured_norm=wpw_projection_captured_norm
          checkpoint%projection_rank=wpw_projection_rank
          checkpoint%projection_charge=wpw_projection_charge
          checkpoint%final_interface_lambda=wpw_context%bounded_operator%interface_lambda
          checkpoint%tolerance_profile='metric_and_scf_residual_tolerance'
          checkpoint%frozen_layout_fingerprint=wpw_frozen_production_context%bounded_operator%layout_fingerprint
          checkpoint%frozen_ww_provenance_fingerprint=&
            wpw_frozen_production_context%ww_provenance_fingerprint
        endif
        checkpoint%peer_ranks=pack(wpw_context%support_fragment_ids-1,&
          wpw_context%support_fragment_ids/=wpw_context%owned_fragment_id)
        checkpoint%owned_w_ids=wpw_context%bounded_operator%owned_w_ids
        checkpoint%owned_p_ids=wpw_context%bounded_operator%owned_p_ids
        checkpoint%support_w_ids=wpw_context%bounded_operator%required_w_ids
        checkpoint%support_p_ids=wpw_context%bounded_operator%required_p_ids
        checkpoint%eigenvalues=wpw_eigenvalues
        checkpoint%occupations=retained_occupations
        checkpoint%coeff_w=wpw_qw
        checkpoint%coeff_p=wpw_qp
        checkpoint%potential=wpw_volume_accumulator%potentials(1:wpw_volume_accumulator%npoint)
        checkpoint%ww_r=wpw_context%bounded_operator%ww_r
        checkpoint%ww_c=wpw_context%bounded_operator%ww_c
        checkpoint%wp_w=wpw_context%bounded_operator%wp_w
        checkpoint%wp_p=wpw_context%bounded_operator%wp_p
        checkpoint%pp_r=wpw_context%bounded_operator%pp_r
        checkpoint%pp_c=wpw_context%bounded_operator%pp_c
        checkpoint%ww_h=wpw_context%bounded_operator%ww_h
        checkpoint%ww_s=wpw_context%bounded_operator%ww_s
        checkpoint%wp_h=wpw_context%bounded_operator%wp_h
        checkpoint%wp_s=wpw_context%bounded_operator%wp_s
        checkpoint%pp_h=wpw_context%bounded_operator%pp_h
        checkpoint%pp_s=wpw_context%bounded_operator%pp_s
        allocate(checkpoint%ww_z(size(checkpoint%ww_r)),checkpoint%wp_z(size(checkpoint%wp_w)),&
          checkpoint%pp_z(size(checkpoint%pp_r)))
        call map_wpw_checkpoint_position_sparse(checkpoint,ww_z_sparse,wp_z_sparse,pp_z_dense,checkpoint_info)
        if(checkpoint_info/=0)local_bad=1
        local_checksum=dg_wpw_checkpoint_checksum(checkpoint)
        write(checkpoint_path,'(a,a,a,i6.6,a)')trim(import_run_root_dir()),trim(sysname),&
          '.dg_wpw.rank_',wpw_context%rank_id,'.chk'
        call write_dg_wpw_checkpoint(trim(checkpoint_path),checkpoint,checkpoint_info)
        if(checkpoint_info==0)call read_dg_wpw_checkpoint(trim(checkpoint_path),verified_checkpoint,&
          checkpoint_info,expected_fingerprint=checkpoint%layout_fingerprint)
        if(checkpoint_info==0)then
          if(dg_wpw_checkpoint_checksum(verified_checkpoint)/=local_checksum)checkpoint_info=1
        endif
        local_bad=merge(0,1,checkpoint_info==0)
        call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr)
        if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
          checkpoint_info=1
        else
          allocate(checksums(wpw_context%nrank))
          call MPI_Gather(local_checksum,1,MPI_INTEGER8,checksums,1,MPI_INTEGER8,0,wpw_production_comm,ierr)
          if(ierr/=MPI_SUCCESS)checkpoint_info=1
          if(wpw_context%rank_id==0.and.checkpoint_info==0)then
            checkpoint_path=trim(import_run_root_dir())//trim(sysname)//'.dg_wpw.manifest'
            call write_dg_wpw_checkpoint_manifest(trim(checkpoint_path),checksums,checkpoint_info)
          endif
          call MPI_Bcast(checkpoint_info,1,MPI_INTEGER,0,wpw_production_comm,ierr)
          if(ierr/=MPI_SUCCESS)checkpoint_info=1
        endif
      endif
      call MPI_Bcast(checkpoint_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)checkpoint_info=1
    end subroutine publish_wpw_production_checkpoint

    subroutine build_wpw_checkpoint_position_volume(ww_z,wp_z,pp_z,position_info)
      complex(8),allocatable,intent(out)::ww_z(:),wp_z(:),pp_z(:,:)
      integer,intent(out)::position_info
      complex(8),allocatable::local_ww(:,:),root_ww(:,:),local_wp(:,:),root_wp(:,:),local_pp(:,:)
      type(s_dg_wpw_staged_candidate),allocatable::staged(:)
      type(s_dg_wpw_owned_candidates)::owned
      integer::point,iw,jw,ip,jp,gid,iz,ierr,row,entry,match,nmatch,route_info,ierr_ww,ierr_wp,ierr_pp
      real(8)::z
      allocate(local_ww(size(wpw_support_w_ids),size(wpw_support_w_ids)),&
        root_ww(size(wpw_support_w_ids),size(wpw_support_w_ids)),&
        local_wp(size(wpw_support_w_ids),size(wpw_context%support_column_ids)),&
        root_wp(size(wpw_support_w_ids),size(wpw_context%support_column_ids)),&
        local_pp(n_plane_waves_dg,size(wpw_context%support_column_ids)))
      allocate(pp_z,source=local_pp)
      local_ww=(0d0,0d0);root_ww=(0d0,0d0);local_wp=(0d0,0d0);root_wp=(0d0,0d0);local_pp=(0d0,0d0)
      do point=1,wpw_volume_accumulator%npoint
        gid=wpw_volume_accumulator%grid_ids(point)
        iz=(gid-1)/(dc%lg_tot%num(1)*dc%lg_tot%num(2))+1
        z=modulo(dble(iz-1)*system%hgs(3)+0.5d0*wpw_box_length(3),wpw_box_length(3))-0.5d0*wpw_box_length(3)
        do jw=1,size(wpw_support_w_ids);do iw=1,size(wpw_support_w_ids)
          local_ww(iw,jw)=local_ww(iw,jw)+conjg(wpw_volume_accumulator%w_points(iw,point))*z*&
            wpw_volume_accumulator%w_points(jw,point)*wpw_volume_accumulator%weights(point)
        enddo;enddo
        do jp=1,size(wpw_context%support_column_ids)
          do iw=1,size(wpw_support_w_ids)
            local_wp(iw,jp)=local_wp(iw,jp)+conjg(wpw_volume_accumulator%w_points(iw,point))*z*&
              wpw_volume_accumulator%p_points(jp,point)*wpw_volume_accumulator%weights(point)
          enddo
          do ip=1,n_plane_waves_dg
            row=find_integer_id(wpw_context%support_column_ids,(dc%i_frag-1)*n_plane_waves_dg+ip)
            if(row<=0)cycle
            local_pp(ip,jp)=local_pp(ip,jp)+conjg(wpw_volume_accumulator%p_points(row,point))*z*&
              wpw_volume_accumulator%p_points(jp,point)*wpw_volume_accumulator%weights(point)
          enddo
        enddo
      enddo
      call MPI_Reduce(local_ww,root_ww,size(root_ww),MPI_DOUBLE_COMPLEX,MPI_SUM,0,dc%icomm_frag,ierr_ww)
      call MPI_Reduce(local_wp,root_wp,size(root_wp),MPI_DOUBLE_COMPLEX,MPI_SUM,0,dc%icomm_frag,ierr_wp)
      call MPI_Reduce(local_pp,pp_z,size(pp_z),MPI_DOUBLE_COMPLEX,MPI_SUM,0,dc%icomm_frag,ierr_pp)
      position_info=merge(0,1,ierr_ww==MPI_SUCCESS.and.ierr_wp==MPI_SUCCESS.and.ierr_pp==MPI_SUCCESS)
      if(dc%id_frag==0.and.position_info==0)then
        allocate(staged(size(root_ww)));entry=0
        do jw=1,size(wpw_support_w_ids);do iw=1,size(wpw_support_w_ids)
          entry=entry+1;staged(entry)%kind=wpw_candidate_kind_ww;staged(entry)%image_id=dc%i_frag
          staged(entry)%target_fragment=wpw_support_w_owner(iw)
          staged(entry)%ww_r=wpw_support_w_ids(iw);staged(entry)%ww_c=wpw_support_w_ids(jw)
          staged(entry)%ww_h=root_ww(iw,jw)
        enddo;enddo
        call route_dg_wpw_candidate_halo(wpw_production_comm,wpw_context%halo_epoch,dc%n_frag,&
          size(wpw_g_indices,2),wpw_support_fragment_ids,staged,owned,route_info,size(staged))
        position_info=route_info
        if(position_info==0)then
          allocate(ww_z(size(wpw_context%bounded_operator%ww_r)));ww_z=(0d0,0d0)
          do iw=1,size(owned%ww_r)
            match=0;nmatch=0
            do entry=1,size(wpw_context%bounded_operator%ww_r)
              if(wpw_context%bounded_operator%ww_r(entry)==owned%ww_r(iw).and.&
                wpw_context%bounded_operator%ww_c(entry)==owned%ww_c(iw))then
                match=entry;nmatch=nmatch+1
              endif
            enddo
            if(nmatch/=1)then;position_info=1;exit;endif
            ww_z(match)=owned%ww_h(iw)
          enddo
          if(size(owned%ww_r)/=size(ww_z))position_info=1
        endif
        if(position_info==0)then
          deallocate(staged);allocate(staged(size(root_wp)));entry=0
          do jp=1,size(wpw_context%support_column_ids);do iw=1,size(wpw_support_w_ids)
            entry=entry+1;staged(entry)%kind=wpw_candidate_kind_wp;staged(entry)%image_id=dc%i_frag
            staged(entry)%wp_w=wpw_support_w_ids(iw)
            staged(entry)%wp_p=wpw_context%support_column_ids(jp);staged(entry)%wp_h=root_wp(iw,jp)
          enddo;enddo
          call route_dg_wpw_candidate_halo(wpw_production_comm,wpw_context%halo_epoch,dc%n_frag,&
            size(wpw_g_indices,2),wpw_support_fragment_ids,staged,owned,route_info,size(staged))
          position_info=route_info
        endif
        if(position_info==0)then
          allocate(wp_z(size(wpw_context%bounded_operator%wp_w)));wp_z=(0d0,0d0)
          do iw=1,size(owned%wp_w)
            match=0;nmatch=0
            do entry=1,size(wpw_context%bounded_operator%wp_w)
              if(wpw_context%bounded_operator%wp_w(entry)==owned%wp_w(iw).and.&
                wpw_context%bounded_operator%wp_p(entry)==owned%wp_p(iw))then
                match=entry;nmatch=nmatch+1
              endif
            enddo
            if(nmatch/=1)then;position_info=1;exit;endif
            wp_z(match)=owned%wp_h(iw)
          enddo
          if(size(owned%wp_w)/=size(wp_z))position_info=1
        endif
      else
        allocate(ww_z(0),wp_z(0))
      endif
      call MPI_Bcast(position_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)position_info=1
    end subroutine build_wpw_checkpoint_position_volume

    subroutine map_wpw_checkpoint_position_sparse(checkpoint,ww_z,wp_z,pp_z,position_info)
      type(s_dg_wpw_checkpoint_state),intent(inout)::checkpoint
      complex(8),intent(in)::ww_z(:),wp_z(:),pp_z(:,:)
      integer,intent(out)::position_info
      integer::i,row,column
      checkpoint%ww_z=(0d0,0d0);checkpoint%wp_z=(0d0,0d0);checkpoint%pp_z=(0d0,0d0)
      if(size(ww_z)/=size(checkpoint%ww_z))then;position_info=1;return;endif
      checkpoint%ww_z=ww_z
      if(size(wp_z)/=size(checkpoint%wp_z))then;position_info=1;return;endif
      checkpoint%wp_z=wp_z
      do i=1,size(checkpoint%pp_z)
        row=checkpoint%pp_r(i)-(dc%i_frag-1)*n_plane_waves_dg
        column=find_integer_id(wpw_context%support_column_ids,checkpoint%pp_c(i))
        if(row>0.and.column>0)checkpoint%pp_z(i)=pp_z(row,column)
      enddo
      position_info=0
    end subroutine map_wpw_checkpoint_position_sparse

    subroutine snapshot_wpw_fixed_point_state(snapshot_info)
      integer,intent(out)::snapshot_info
      snapshot_info=0
      if(dc%id_frag==0)then
        wpw_saved_wp_volume=wpw_context%wp_h_volume
        wpw_saved_pp_volume=wpw_context%pp_h_volume
        wpw_saved_ww_potential=wpw_context%realspace_ww_potential
      endif
      wpw_saved_accumulator_potential=&
        wpw_volume_accumulator%potentials(1:wpw_volume_accumulator%npoint)
      wpw_saved_accumulator_density=&
        wpw_volume_accumulator%densities(1:wpw_volume_accumulator%npoint)
      wpw_saved_rho_tot=dc%rho_tot%f;wpw_saved_rho_tot_s=dc%rho_tot_s(1)%f
      wpw_saved_vh_tot=dc%vh_tot%f
      wpw_saved_vxc_tot=dc%vxc_tot(1)%f
      if(.not.wpw_potential_stage_ok(snapshot_info))snapshot_info=1
    end subroutine snapshot_wpw_fixed_point_state

    subroutine restore_wpw_fixed_point_state(restore_info)
      integer,intent(out)::restore_info
      integer::ierr
      restore_info=0
      wpw_volume_accumulator%potentials(1:wpw_volume_accumulator%npoint)=&
        wpw_saved_accumulator_potential
      wpw_volume_accumulator%densities(1:wpw_volume_accumulator%npoint)=&
        wpw_saved_accumulator_density
      dc%rho_tot%f=wpw_saved_rho_tot;dc%rho_tot_s(1)%f=wpw_saved_rho_tot_s
      dc%vh_tot%f=wpw_saved_vh_tot
      dc%vxc_tot(1)%f=wpw_saved_vxc_tot
      if(dc%id_frag==0)then
        call replace_dg_wpw_potential_volume(wpw_context,wpw_saved_wp_volume,&
          wpw_saved_pp_volume,restore_info,ww_potential=wpw_saved_ww_potential)
        if(restore_info==0)call build_dg_wpw_production_operator(wpw_context,restore_info)
        if(restore_info==0)call bind_dg_wpw_hs_callbacks(wpw_context,restore_info)
      endif
      call MPI_Bcast(restore_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)restore_info=1
      if(.not.wpw_potential_stage_ok(restore_info))restore_info=1
    end subroutine restore_wpw_fixed_point_state

    subroutine wpw_algebra_step(iteration,cw,cp,gap,residual,orth,projector,converged,algebra_info)
      integer,intent(in)::iteration
      complex(8),intent(out)::cw(:,:),cp(:,:)
      real(8),intent(out)::gap,residual,orth,projector
      logical,intent(out)::converged
      integer,intent(out)::algebra_info
      call run_dg_wpw_matrix_free_algebra_step(wpw_context,wpw_production_comm,wpw_apply_h,&
        wpw_apply_s,wpw_global_gram,size(wpw_qw,1),size(wpw_qp,1),wpw_nocc,wpw_nretain,&
        iteration,dg_wpw_metric_cutoff,dg_wpw_scf_residual_tolerance,wpw_qw,wpw_qp,&
        wpw_q_old_occ,wpw_eigenvalues,gap,residual,orth,projector,algebra_info)
      if(algebra_info/=0)write(*,'(1x,a,i0,a,i0,a,i0,4(a,es12.4))')'[DG-WPW-LOCAL-FAIL] algebra rank=',&
        dc%id_tot,' iter=',iteration,' info=',algebra_info,' residual=',residual,' orth=',orth,&
        ' eval_min=',minval(wpw_eigenvalues),' eval_max=',maxval(wpw_eigenvalues)
      converged=algebra_info==0.and.gap>=dg_wpw_gap_threshold.and.&
        max(residual,max(orth,projector))<dg_wpw_scf_residual_tolerance
      if(algebra_info==0)then
        cw=wpw_qw(:,1:wpw_nocc);cp=wpw_qp(:,1:wpw_nocc)
        wpw_q_old_occ(1:size(wpw_qw,1),:)=cw
        wpw_q_old_occ(size(wpw_qw,1)+1:,:)=cp
      endif
    end subroutine wpw_algebra_step

    subroutine wpw_precondition(context,eigenvalues,rw,rp,zw,zp,precondition_info)
      class(*),intent(inout)::context
      real(8),intent(in)::eigenvalues(:)
      complex(8),intent(in)::rw(:,:),rp(:,:)
      complex(8),intent(out)::zw(:,:),zp(:,:)
      integer,intent(out)::precondition_info
      integer::iw,ip,state,wpos,ppos,local_bad,global_bad,ierr
      real(8)::hdiag,sdiag,denominator,diagonal_scale,preconditioner_floor
      zw=0;zp=0;local_bad=0
      select type(context)
      type is(s_dg_wpw_production_context)
        if(size(eigenvalues)/=size(rw,2).or.size(rp,2)/=size(rw,2).or.&
           any(shape(zw)/=shape(rw)).or.any(shape(zp)/=shape(rp)).or.&
           size(rw,1)/=size(context%bounded_operator%owned_w_ids).or.&
           size(rp,1)/=size(context%bounded_operator%owned_p_ids))then
          local_bad=1
        else
          do state=1,size(eigenvalues)
            do iw=1,size(rw,1)
              wpos=find_integer_id(context%bounded_operator%required_w_ids,&
                context%bounded_operator%owned_w_ids(iw))
              if(wpos<=0)then;local_bad=1;cycle;endif
              hdiag=real(context%bounded_operator%ww_h_dense(wpos,wpos),8)
              sdiag=real(context%bounded_operator%ww_s_dense(wpos,wpos),8)
              denominator=hdiag-eigenvalues(state)*sdiag
              diagonal_scale=max(1d0,abs(hdiag),abs(eigenvalues(state)*sdiag))
              preconditioner_floor=sqrt(epsilon(1d0))*diagonal_scale
              if(abs(denominator)<preconditioner_floor)&
                denominator=merge(preconditioner_floor,-preconditioner_floor,denominator>=0d0)
              zw(iw,state)=-rw(iw,state)/denominator
            enddo
            do ip=1,size(rp,1)
              ppos=find_integer_id(context%bounded_operator%required_p_ids,&
                context%bounded_operator%owned_p_ids(ip))
              if(ppos<=0)then;local_bad=1;cycle;endif
              hdiag=real(context%bounded_operator%pp_h_dense(ip,ppos),8)
              sdiag=real(context%bounded_operator%pp_s_dense(ip,ppos),8)
              denominator=hdiag-eigenvalues(state)*sdiag
              diagonal_scale=max(1d0,abs(hdiag),abs(eigenvalues(state)*sdiag))
              preconditioner_floor=sqrt(epsilon(1d0))*diagonal_scale
              if(abs(denominator)<preconditioner_floor)&
                denominator=merge(preconditioner_floor,-preconditioner_floor,denominator>=0d0)
              zp(ip,state)=-rp(ip,state)/denominator
            enddo
          enddo
        endif
      class default
        local_bad=1
      end select
      if(.not.all(ieee_is_finite(real(zw,8))).or..not.all(ieee_is_finite(aimag(zw))).or.&
         .not.all(ieee_is_finite(real(zp,8))).or..not.all(ieee_is_finite(aimag(zp))))local_bad=1
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr)
      precondition_info=merge(0,1,ierr==MPI_SUCCESS.and.global_bad==0)
      if(precondition_info/=0)then;zw=0;zp=0;endif
    end subroutine wpw_precondition

    subroutine wpw_apply_h(context,xw,xp,yw,yp,apply_info)
      class(*),intent(inout)::context
      complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
      integer,intent(out)::apply_info
      select type(context)
      type is(s_dg_wpw_production_context)
        call apply_h_dg_wpw_bounded(context%bounded_operator,context%bounded_operator%operator_epoch,&
          context%bounded_operator%layout_fingerprint,xw,xp,yw,yp,apply_info)
      class default
        yw=0;yp=0;apply_info=1
      end select
    end subroutine wpw_apply_h

    subroutine wpw_apply_s(context,xw,xp,yw,yp,apply_info)
      class(*),intent(inout)::context
      complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
      integer,intent(out)::apply_info
      select type(context)
      type is(s_dg_wpw_production_context)
        call apply_s_dg_wpw_bounded(context%bounded_operator,context%bounded_operator%operator_epoch,&
          context%bounded_operator%layout_fingerprint,xw,xp,yw,yp,apply_info)
      class default
        yw=0;yp=0;apply_info=1
      end select
    end subroutine wpw_apply_s

    subroutine wpw_apply_metric_block_preconditioner(context,rw,rp,zw,zp,apply_info)
      class(*),intent(inout)::context
      complex(8),intent(in)::rw(:,:),rp(:,:)
      complex(8),intent(out)::zw(:,:),zp(:,:)
      integer,intent(out)::apply_info
      select type(c=>context)
      type is(s_dg_wpw_production_context)
        call apply_dg_wpw_fragment_block_preconditioner(c%bounded_operator,&
          c%metric_block_preconditioner,rw,rp,zw,zp,apply_info)
      class default
        zw=0;zp=0;apply_info=1
      end select
    end subroutine wpw_apply_metric_block_preconditioner

    subroutine wpw_global_gram(x,y,nrow,nx,ny,g,gram_info)
      integer,intent(in)::nrow,nx,ny
      complex(8),intent(in)::x(nrow,nx),y(nrow,ny);complex(8),intent(out)::g(nx,ny)
      integer,intent(out)::gram_info
      integer::nw
      nw=size(wpw_qw,1)
      call global_gram_dg_wpw_bounded(wpw_context%bounded_operator,x(1:nw,:),x(nw+1:,:),&
        y(1:nw,:),y(nw+1:,:),g,gram_info)
    end subroutine wpw_global_gram

    subroutine wpw_potential_step(iteration,cw,cp,density_in,density_out,potential_residual,&
        total_energy,charge_error,potential_info)
      integer,intent(in)::iteration
      complex(8),intent(in)::cw(:,:),cp(:,:)
      real(8),intent(in)::density_in(:)
      real(8),intent(out)::density_out(:),potential_residual,total_energy,charge_error
      integer,intent(out)::potential_info
      complex(8),allocatable::rw(:,:),rw_support(:,:),rp(:,:),local_wp_h(:,:),local_wp_s(:,:),local_pp_h(:,:),local_pp_s(:,:)
      complex(8),allocatable::root_wp_h(:,:),root_wp_s(:,:),root_pp_h(:,:),root_pp_s(:,:)
      complex(8),allocatable::local_ww_potential(:,:),root_ww_potential(:,:)
      real(8),allocatable::rho_raw(:),rho_candidate(:),rho_owned(:,:),potential_owned(:),potential_core(:)
      integer,allocatable::owned_ids(:),destinations(:)
      real(8)::charge,e_xc,e_xc_owned,n_vxc,e_h,band_energy,wpw_mix,energy_min_total,energy_max_total
      real(8)::potential_norms_local(2),potential_norms_global(2)
      integer::stage_info,ierr,iw,jw,ip,jp,point,nw_support,np_support,np_owned,wpos_i,wpos_j

      potential_info=1;density_out=0d0;potential_residual=huge(1d0);total_energy=huge(1d0);charge_error=huge(1d0)
      nw_support=size(wpw_support_w_ids);np_owned=size(wpw_g_indices,2)
      np_support=size(wpw_support_fragment_ids)*np_owned
      allocate(rw_support(nw_support,wpw_nocc),rw(size(wpw_owned_w_ids),wpw_nocc),rp(np_support,wpw_nocc))
      rw_support=0;rw=0;rp=0
      stage_info=0
      if(dc%id_frag==0)call fetch_dg_wpw_support_coefficients(wpw_context%bounded_operator,cw,cp,rw_support,rp,stage_info)
      call MPI_Bcast(stage_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr);if(ierr/=MPI_SUCCESS.or.stage_info/=0)return
      call MPI_Bcast(rw_support,size(rw_support),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Bcast(rp,size(rp),MPI_DOUBLE_COMPLEX,0,dc%icomm_frag,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Bcast(wpw_eigenvalues,size(wpw_eigenvalues),MPI_DOUBLE_PRECISION,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)return
      do iw=1,size(wpw_owned_w_ids)
        wpos_i=find_integer_id(wpw_support_w_ids,wpw_owned_w_ids(iw));if(wpos_i==0)return
        rw(iw,:)=rw_support(wpos_i,:)
      enddo
      allocate(rho_raw(size(density_in)),rho_candidate(size(density_in)))
      call build_dg_wpw_core_density(wpw_volume_accumulator%w_points(:,1:size(density_in)),&
        wpw_volume_accumulator%p_points(:,1:size(density_in)),rw,rp,wpw_occupations,&
        wpw_volume_accumulator%weights(1:size(density_in)),rho_raw,charge,stage_info)
      if(.not.wpw_potential_stage_ok(stage_info))then
        write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] potential_density rank=',dc%id_tot,' info=',stage_info;return
      endif
      wpw_mix=merge(1d0,dg_wpw_scf_mix,wpw_fixed_point_mode)
      rho_candidate=(1d0-wpw_mix)*density_in+wpw_mix*rho_raw
      call build_wpw_total_grid_route(rho_candidate,owned_ids,destinations,rho_owned,charge,stage_info)
      if(.not.wpw_potential_stage_ok(stage_info))then
        write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] potential_grid_route rank=',dc%id_tot,' info=',stage_info;return
      endif
      call update_wpw_owned_lda_hartree(dc%lg_tot,dc%mg_tot,dc%system_tot,dc%info_tot,stencil,xc_func,pp,ppn,&
        spsi,dc%srg_tot,dc%srg_scalar_tot,dc%poisson_tot,dc%fg_tot,dc%rho_tot,dc%rho_tot_s,dc%vh_tot,&
        dc%vxc_tot,rho_owned,size(owned_ids),e_xc,e_xc_owned,n_vxc,stage_info)
      if(.not.wpw_potential_stage_ok(stage_info))then
        write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] potential_lda_hartree rank=',dc%id_tot,' info=',stage_info;return
      endif
      allocate(potential_owned(size(owned_ids)),potential_core(size(density_in)))
      potential_owned=reshape(dc%vloc_tot(1)%f+dc%vh_tot%f+dc%vxc_tot(1)%f,[size(owned_ids)])
      call return_dg_wpw_core_potential(wpw_density_route,owned_ids,potential_owned,potential_core,stage_info)
      if(.not.wpw_potential_stage_ok(stage_info))then
        write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] potential_return rank=',dc%id_tot,' info=',stage_info;return
      endif
      allocate(local_wp_h(size(wpw_owned_w_ids),np_support),local_wp_s(size(wpw_owned_w_ids),np_support),&
        local_pp_h(np_support,np_support),local_pp_s(np_support,np_support))
      call accumulate_dg_wpw_core_volume(wpw_volume_accumulator%w_points(:,1:size(density_in)),&
        wpw_volume_accumulator%grad_w_points(:,:,1:size(density_in)),&
        wpw_volume_accumulator%p_points(:,1:size(density_in)),&
        wpw_volume_accumulator%grad_p_points(:,:,1:size(density_in)),&
        wpw_volume_accumulator%p_points(:,1:size(density_in)),&
        wpw_volume_accumulator%grad_p_points(:,:,1:size(density_in)),potential_core,&
        wpw_volume_accumulator%weights(1:size(density_in)),local_wp_h,local_wp_s,local_pp_h,local_pp_s,stage_info)
      if(.not.wpw_potential_stage_ok(stage_info))then
        write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] potential_volume rank=',dc%id_tot,' info=',stage_info;return
      endif
      allocate(root_wp_h,source=local_wp_h);allocate(root_wp_s,source=local_wp_s)
      allocate(root_pp_h,source=local_pp_h);allocate(root_pp_s,source=local_pp_s)
      allocate(local_ww_potential(nw_support,nw_support),root_ww_potential(nw_support,nw_support))
      local_ww_potential=(0d0,0d0);root_ww_potential=(0d0,0d0)
      do jw=1,nw_support
        do iw=1,nw_support
          do point=1,size(density_in)
            local_ww_potential(iw,jw)=local_ww_potential(iw,jw)+conjg(&
              wpw_volume_accumulator%w_points(iw,point))*potential_core(point)*&
              wpw_volume_accumulator%w_points(jw,point)*wpw_volume_accumulator%weights(point)
          enddo
        enddo
      enddo
      call MPI_Reduce(local_wp_h,root_wp_h,size(local_wp_h),MPI_DOUBLE_COMPLEX,MPI_SUM,0,dc%icomm_frag,ierr)
      call MPI_Reduce(local_wp_s,root_wp_s,size(local_wp_s),MPI_DOUBLE_COMPLEX,MPI_SUM,0,dc%icomm_frag,ierr)
      call MPI_Reduce(local_pp_h,root_pp_h,size(local_pp_h),MPI_DOUBLE_COMPLEX,MPI_SUM,0,dc%icomm_frag,ierr)
      call MPI_Reduce(local_pp_s,root_pp_s,size(local_pp_s),MPI_DOUBLE_COMPLEX,MPI_SUM,0,dc%icomm_frag,ierr)
      call MPI_Reduce(local_ww_potential,root_ww_potential,size(local_ww_potential),MPI_DOUBLE_COMPLEX,&
        MPI_SUM,0,dc%icomm_frag,ierr)
      if(ierr/=MPI_SUCCESS)return
      density_out=rho_candidate
      potential_norms_local=[sum((potential_core-wpw_volume_accumulator%potentials(1:size(density_in)))**2),&
        sum(potential_core**2)]
      call MPI_Allreduce(potential_norms_local,potential_norms_global,2,MPI_DOUBLE_PRECISION,MPI_SUM,&
        dc%icomm_tot,ierr);if(ierr/=MPI_SUCCESS)return
      potential_residual=sqrt(max(0d0,potential_norms_global(1)))/&
        max(1d-30,sqrt(max(0d0,potential_norms_global(2))))
      call hartree_energy_global(reshape(rho_owned(:,1),[size(owned_ids)]),reshape(dc%vh_tot%f,[size(owned_ids)]),&
        [(system%hvol,iw=1,size(owned_ids))],size(owned_ids),dc%icomm_tot,e_h,stage_info)
      if(.not.wpw_potential_stage_ok(stage_info))then
        write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] potential_hartree_energy rank=',dc%id_tot,' info=',stage_info;return
      endif
      band_energy=sum(wpw_occupations*wpw_eigenvalues(1:wpw_nocc))
      total_energy=band_energy-e_h-n_vxc+e_xc+energy%E_ion_ion
      call MPI_Allreduce(total_energy,energy_min_total,1,MPI_DOUBLE_PRECISION,MPI_MIN,dc%icomm_tot,ierr)
      call MPI_Allreduce(total_energy,energy_max_total,1,MPI_DOUBLE_PRECISION,MPI_MAX,dc%icomm_tot,ierr)
      if(ierr/=MPI_SUCCESS.or.abs(energy_max_total-energy_min_total)>&
        100d0*epsilon(1d0)*max(1d0,abs(energy_max_total)))return
      total_energy=energy_max_total
      charge_error=charge-sum(wpw_occupations)
      stage_info=merge(0,1,abs(charge_error)<huge(1d0).and.ieee_is_finite(total_energy))
      if(.not.wpw_potential_stage_ok(stage_info))then
        write(*,'(1x,a,i0,a,i0)')'[DG-WPW-LOCAL-FAIL] potential_energy_charge rank=',dc%id_tot,' info=',stage_info;return
      endif
      if(dc%id_frag==0)call publish_wpw_iterated_operator(root_wp_h,root_pp_h,root_ww_potential,stage_info)
      call MPI_Bcast(stage_info,1,MPI_INTEGER,0,dc%icomm_frag,ierr);if(ierr/=MPI_SUCCESS.or.stage_info/=0)return
      wpw_volume_accumulator%potentials(1:size(density_in))=potential_core
      wpw_volume_accumulator%densities(1:size(density_in))=rho_candidate
      potential_info=0
    end subroutine wpw_potential_step

    logical function wpw_potential_stage_ok(stage_info)result(ok)
      integer,intent(in)::stage_info
      integer::local_bad,global_bad,ierr
      local_bad=merge(0,1,stage_info==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
      ok=ierr==MPI_SUCCESS.and.global_bad==0
    end function wpw_potential_stage_ok

    subroutine build_wpw_total_grid_route(rho_candidate,owned_ids,destinations,rho_owned,charge,route_info)
      real(8),intent(in)::rho_candidate(:)
      integer,allocatable,intent(out)::owned_ids(:),destinations(:)
      real(8),allocatable,intent(out)::rho_owned(:,:)
      real(8),intent(out)::charge
      integer,intent(out)::route_info
      integer::ixg,iyg,izg,p,entry,gid,owner,nowned
      nowned=product(dc%mg_tot%num);allocate(owned_ids(nowned),rho_owned(nowned,nspin));entry=0
      do izg=dc%mg_tot%is(3),dc%mg_tot%ie(3);do iyg=dc%mg_tot%is(2),dc%mg_tot%ie(2)
        do ixg=dc%mg_tot%is(1),dc%mg_tot%ie(1)
          entry=entry+1;owned_ids(entry)=wpw_core_global_grid_id([ixg,iyg,izg])
        enddo
      enddo;enddo
      allocate(destinations(size(rho_candidate)))
      do entry=1,size(rho_candidate)
        gid=wpw_volume_accumulator%grid_ids(entry)
        ixg=modulo(gid-1,dc%lg_tot%num(1))+1
        iyg=modulo((gid-1)/dc%lg_tot%num(1),dc%lg_tot%num(2))+1
        izg=(gid-1)/(dc%lg_tot%num(1)*dc%lg_tot%num(2))+1;owner=-1
        do p=0,dc%isize_tot-1
          if(all([ixg,iyg,izg]>=dc%mg_tot%is_all(:,p)).and.&
             all([ixg,iyg,izg]<=dc%mg_tot%ie_all(:,p)))then;owner=p;exit;endif
        enddo
        destinations(entry)=owner
      enddo
      call prepare_dg_wpw_core_density_route(dc%icomm_tot,wpw_volume_accumulator%grid_ids(1:size(rho_candidate)),&
        rho_candidate,wpw_volume_accumulator%weights(1:size(rho_candidate)),product(dc%lg_tot%num),owned_ids,&
        wpw_density_route,rho_owned(:,1),charge,route_info,destinations)
    end subroutine build_wpw_total_grid_route

    subroutine publish_wpw_iterated_operator(root_wp_h,root_pp_h,root_ww_potential,publish_info)
      complex(8),intent(in)::root_wp_h(:,:),root_pp_h(:,:)
      complex(8),intent(in)::root_ww_potential(:,:)
      integer,intent(out)::publish_info
      type(s_dg_wpw_staged_candidate),allocatable::ww_staged(:)
      type(s_dg_wpw_owned_candidates)::ww_owned
      complex(8),allocatable::wp_volume(:),pp_volume(:),old_wp_volume(:),old_pp_volume(:),&
        old_ww_potential(:),candidate_ww_potential(:)
      integer,allocatable::wp_w(:),wp_p(:),pp_r(:),pp_c(:)
      integer::entry,iw,jw,ip,rollback_info,nmatch,match,route_info
      logical::old_operator_valid,old_callbacks_bound
      publish_info=1
      route_info=merge(0,1,all(ieee_is_finite(real(root_ww_potential,8))).and.&
        all(ieee_is_finite(aimag(root_ww_potential))).and.&
        maxval(abs(root_ww_potential-conjg(transpose(root_ww_potential))))<=1d-12)
      call MPI_Allreduce(route_info,publish_info,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,rollback_info)
      if(rollback_info/=MPI_SUCCESS.or.publish_info/=0)return
      allocate(old_wp_volume,source=wpw_context%wp_h_volume)
      allocate(old_pp_volume,source=wpw_context%pp_h_volume)
      allocate(old_ww_potential,source=wpw_context%realspace_ww_potential)
      old_operator_valid=wpw_context%operator_valid;old_callbacks_bound=wpw_context%callbacks_bound
      allocate(wp_w(size(root_wp_h)),wp_p(size(root_wp_h)),wp_volume(size(root_wp_h)))
      entry=0
      do ip=1,size(wpw_context%support_column_ids);do iw=1,size(wpw_owned_w_ids)
        entry=entry+1;wp_w(entry)=wpw_owned_w_ids(iw);wp_p(entry)=wpw_context%support_column_ids(ip)
        wp_volume(entry)=root_wp_h(iw,ip)
      enddo;enddo
      allocate(pp_r(size(root_pp_h)),pp_c(size(root_pp_h)),pp_volume(size(root_pp_h)))
      entry=0
      do ip=1,size(wpw_context%support_column_ids);do iw=1,size(wpw_context%support_column_ids)
        entry=entry+1;pp_r(entry)=wpw_context%support_column_ids(iw)
        pp_c(entry)=wpw_context%support_column_ids(ip);pp_volume(entry)=root_pp_h(iw,ip)
      enddo;enddo
      call route_and_replace_dg_wpw_potential_volume(wpw_context,wpw_context%halo_epoch,&
        wpw_support_fragment_ids,wp_w,wp_p,wp_volume,pp_r,pp_c,pp_volume,&
        size(wp_w)+size(pp_r),publish_info)
      if(publish_info/=0)return
      allocate(ww_staged(size(root_ww_potential)));entry=0
      do jw=1,size(wpw_support_w_ids);do iw=1,size(wpw_support_w_ids)
        entry=entry+1;ww_staged(entry)%kind=wpw_candidate_kind_ww
        ww_staged(entry)%image_id=dc%i_frag;ww_staged(entry)%target_fragment=wpw_support_w_owner(iw)
        ww_staged(entry)%ww_r=wpw_support_w_ids(iw);ww_staged(entry)%ww_c=wpw_support_w_ids(jw)
        ww_staged(entry)%ww_potential=root_ww_potential(iw,jw)
      enddo;enddo
      call route_dg_wpw_candidate_halo(wpw_production_comm,wpw_context%halo_epoch,dc%n_frag,&
        size(wpw_g_indices,2),wpw_support_fragment_ids,ww_staged,ww_owned,route_info,size(ww_staged))
      if(route_info/=0)then
        call replace_dg_wpw_potential_volume(wpw_context,old_wp_volume,old_pp_volume,rollback_info,&
          ww_potential=old_ww_potential);return
      endif
      allocate(candidate_ww_potential(size(wpw_context%realspace_ww_potential)))
      candidate_ww_potential=(0d0,0d0)
      do iw=1,size(ww_owned%ww_r)
        match=0;nmatch=0
        do entry=1,size(wpw_context%realspace_ww_r)
          if(wpw_context%realspace_ww_r(entry)==ww_owned%ww_r(iw).and.&
            wpw_context%realspace_ww_c(entry)==ww_owned%ww_c(iw))then
            match=entry;nmatch=nmatch+1
          endif
        enddo
        if(nmatch/=1)then;route_info=1;exit;endif
        candidate_ww_potential(match)=ww_owned%ww_potential(iw)
      enddo
      if(size(ww_owned%ww_r)/=size(candidate_ww_potential))route_info=1
      call MPI_Allreduce(route_info,publish_info,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,rollback_info)
      if(rollback_info/=MPI_SUCCESS.or.publish_info/=0)then
        call replace_dg_wpw_potential_volume(wpw_context,old_wp_volume,old_pp_volume,rollback_info,&
          ww_potential=old_ww_potential);return
      endif
      call replace_dg_wpw_potential_volume(wpw_context,wpw_context%wp_h_volume,&
        wpw_context%pp_h_volume,publish_info,ww_potential=candidate_ww_potential)
      if(publish_info/=0)then
        call replace_dg_wpw_potential_volume(wpw_context,old_wp_volume,old_pp_volume,rollback_info,&
          ww_potential=old_ww_potential);return
      endif
      call build_dg_wpw_production_operator(wpw_context,publish_info)
      if(publish_info/=0)then
        call replace_dg_wpw_potential_volume(wpw_context,old_wp_volume,old_pp_volume,rollback_info,&
          ww_potential=old_ww_potential)
        wpw_context%pending_ww_valid=.false.;wpw_context%operator_valid=old_operator_valid
        wpw_context%callbacks_bound=old_callbacks_bound
        return
      endif
      call bind_dg_wpw_hs_callbacks(wpw_context,publish_info)
      if(publish_info==0)then
      else
        call replace_dg_wpw_potential_volume(wpw_context,old_wp_volume,old_pp_volume,rollback_info,&
          ww_potential=old_ww_potential)
        wpw_context%pending_ww_valid=.false.;wpw_context%operator_valid=old_operator_valid
        wpw_context%callbacks_bound=old_callbacks_bound
      endif
    end subroutine publish_wpw_iterated_operator

    integer function find_support_fragment(fragment)result(position)
      integer,intent(in)::fragment
      position=find_integer_id(wpw_support_fragment_ids,fragment)
    end function find_support_fragment

    integer function find_integer_id(ids,target)result(position)
      integer,intent(in)::ids(:),target
      integer::ii
      position=0
      do ii=1,size(ids);if(ids(ii)==target)then;position=ii;return;endif;enddo
    end function find_integer_id

    subroutine prepare_wpw_volume_halo(halo_info)
      use salmon_global,only:num_fragment
      integer,intent(out)::halo_info
      integer::domain(3),local_lo(3),local_hi(3),buffer(3),buffer_lo(3),buffer_hi(3),buffer_extent(3),&
        descriptor_extent(3)
      integer::target_lo(3),target_hi(3),send_lo(3),send_hi(3),cell_extent(3)
      integer::k,ixb,iyb,izb,point,local_grid(3),pack_info,local_failure,global_failure,ierr_halo,&
        descriptor_point
      integer::max_location(2),best_location(2),best_send,best_grid(3),best_unwrapped(3),source_point
      complex(8),allocatable::buffer_values(:,:),buffer_gradients(:,:,:)
      complex(8)::packed_max_value,prepack_max_value
      real(8)::packed_max_abs
      character(32)::buffer_path

      halo_info=1;local_failure=0
      if(dc%id_frag/=0.or.wpw_production_comm==MPI_COMM_NULL)local_failure=1
      if(.not.wpw_occupied_w_basis%valid.or.&
        wpw_occupied_w_basis%local_fragment/=dc%i_frag.or.&
        wpw_occupied_w_basis%local_count/=size(wpw_owned_w_ids))local_failure=1
      if(local_failure==0)then
        if(any(wpw_occupied_w_basis%owned_ids/=wpw_owned_w_ids))local_failure=1
      endif
      buffer=wpw_window_buffer
      if(any(buffer<0).or.any(buffer>dc%nxyz_buffer))local_failure=1
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_halo)
      if(ierr_halo/=MPI_SUCCESS.or.global_failure/=0)return
      call get_fragment_domain(dc,dc%i_frag,domain)
      local_lo=dc%ixyz_frag(:,dc%i_frag);local_hi=local_lo+domain-1
      buffer_lo=local_lo-buffer;buffer_hi=local_hi+buffer;buffer_extent=buffer_hi-buffer_lo+1
      descriptor_extent=wpw_occupied_w_basis%buffer_hi-wpw_occupied_w_basis%buffer_lo+1
      cell_extent=dc%lg_tot%num
      allocate(buffer_values(size(wpw_owned_w_ids),product(buffer_extent)),&
        buffer_gradients(3,size(wpw_owned_w_ids),product(buffer_extent)))
      buffer_values=(0d0,0d0);buffer_gradients=(0d0,0d0)
      point=0
      do izb=buffer_lo(3),buffer_hi(3)
        do iyb=buffer_lo(2),buffer_hi(2)
          do ixb=buffer_lo(1),buffer_hi(1)
            point=point+1;local_grid=[ixb,iyb,izb]-local_lo+1
            descriptor_point=local_grid(1)-wpw_occupied_w_basis%buffer_lo(1)+1+&
              (local_grid(2)-wpw_occupied_w_basis%buffer_lo(2))*descriptor_extent(1)+&
              (local_grid(3)-wpw_occupied_w_basis%buffer_lo(3))*descriptor_extent(1)*descriptor_extent(2)
            if(descriptor_point<1.or.descriptor_point>size(wpw_occupied_w_basis%buffer_values,1))then
              local_failure=1;cycle
            endif
            if(any(local_grid<wpw_occupied_w_basis%gradient_lo).or.&
              any(local_grid>wpw_occupied_w_basis%gradient_hi))then
              local_failure=1;cycle
            endif
            buffer_values(:,point)=wpw_occupied_w_basis%buffer_values(descriptor_point,:)
            buffer_gradients(:,:,point)=wpw_occupied_w_basis%buffer_gradients(:,descriptor_point,:)
          enddo
        enddo
      enddo
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_halo)
      if(ierr_halo/=MPI_SUCCESS.or.global_failure/=0)return
      if(size(wpw_support_records)==0)then
        allocate(wpw_volume_send(0),wpw_volume_halos(0));halo_info=0;return
      endif
      allocate(wpw_volume_send(size(wpw_support_records)))
      packed_max_abs=-1d0;best_send=0;best_location=0
      do k=1,size(wpw_support_records)
        target_lo=local_lo+wpw_support_records(k)%image_shift*domain
        target_hi=target_lo+domain-1
        call dg_wpw_support_overlap_box(target_lo,target_hi,-wpw_support_records(k)%image_shift,&
          buffer,send_lo,send_hi,pack_info)
        if(pack_info/=0)then;local_failure=1;exit;endif
        call pack_dg_wpw_volume_halo_send(wpw_support_records(k)%fragment_id-1,&
          wpw_support_records(k)%image_shift,wpw_support_records(k)%periodic_shift,send_lo,send_hi,&
          cell_extent,buffer_lo,buffer_hi,wpw_owned_w_ids,buffer_values,buffer_gradients,&
          wpw_volume_send(k),pack_info)
        if(pack_info/=0)then;local_failure=1;exit;endif
        max_location=maxloc(abs(wpw_volume_send(k)%values))
        if(abs(wpw_volume_send(k)%values(max_location(1),max_location(2)))>packed_max_abs)then
          packed_max_abs=abs(wpw_volume_send(k)%values(max_location(1),max_location(2)))
          best_send=k;best_location=max_location
        endif
      enddo
      if(local_failure==0.and.best_send>0)then
        local_grid=wpw_volume_send(best_send)%box_hi-wpw_volume_send(best_send)%box_lo+1
        point=best_location(2)-1
        best_grid(1)=wpw_volume_send(best_send)%box_lo(1)+modulo(point,local_grid(1))
        best_grid(2)=wpw_volume_send(best_send)%box_lo(2)+modulo(point/local_grid(1),local_grid(2))
        best_grid(3)=wpw_volume_send(best_send)%box_lo(3)+point/(local_grid(1)*local_grid(2))
        best_unwrapped=best_grid+wpw_support_records(best_send)%periodic_shift*cell_extent
        source_point=best_unwrapped(1)-buffer_lo(1)+1+&
          (best_unwrapped(2)-buffer_lo(2))*buffer_extent(1)+&
          (best_unwrapped(3)-buffer_lo(3))*buffer_extent(1)*buffer_extent(2)
        packed_max_value=wpw_volume_send(best_send)%values(best_location(1),best_location(2))
        prepack_max_value=buffer_values(best_location(1),source_point)
        buffer_path='occupied_w_descriptor'
        write(*,'(1x,a,5(a,i0),4(a,3(i0,1x)),3(a,es12.4),2a)')'[DG-WPW-W-HALO-PACK-MAX]',&
          ' rank=',wpw_context%rank_id,' source_fragment=',dc%i_frag,&
          ' target_rank=',wpw_volume_send(best_send)%peer,&
          ' target_fragment=',wpw_support_records(best_send)%fragment_id,&
          ' w_id=',wpw_volume_send(best_send)%w_ids(best_location(1)),&
          ' canonical_grid=',best_grid,' unwrapped_grid=',best_unwrapped,&
          ' route_image=',wpw_support_records(best_send)%image_shift,&
          ' periodic_image=',wpw_support_records(best_send)%periodic_shift,&
          ' abs_value=',packed_max_abs,' prepack_abs=',abs(prepack_max_value),&
          ' pack_defect=',abs(packed_max_value-prepack_max_value),' path=',trim(buffer_path)
      endif
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_halo)
      if(ierr_halo/=MPI_SUCCESS.or.global_failure/=0)return
      call exchange_dg_wpw_volume_halo_schedule(wpw_production_comm,1,wpw_volume_send,&
        wpw_volume_halos,halo_info)
    end subroutine prepare_wpw_volume_halo

    subroutine build_wpw_canonical_face_geometry(face_info)
      use salmon_global,only:num_fragment
      integer,intent(out)::face_info
      call build_dg_wpw_canonical_face_schedule(dc%i_frag,num_fragment,wpw_canonical_faces,face_info)
    end subroutine build_wpw_canonical_face_geometry

    subroutine broadcast_wpw_w_row_layout(layout_info)
      use mpi,only:MPI_Bcast
      integer,intent(out)::layout_info
      integer::counts(2),ierr_layout
      layout_info=1;counts=0
      if(dc%id_frag==0)counts=[size(wpw_owned_w_ids),size(wpw_support_w_ids)]
      call MPI_Bcast(counts,2,MPI_INTEGER,0,dc%icomm_frag,ierr_layout)
      if(ierr_layout/=MPI_SUCCESS.or.any(counts<=0))return
      if(dc%id_frag/=0)then
        if(allocated(wpw_owned_w_ids))deallocate(wpw_owned_w_ids)
        if(allocated(wpw_support_w_ids))deallocate(wpw_support_w_ids)
        allocate(wpw_owned_w_ids(counts(1)),wpw_support_w_ids(counts(2)))
      endif
      call MPI_Bcast(wpw_owned_w_ids,counts(1),MPI_INTEGER,0,dc%icomm_frag,ierr_layout)
      if(ierr_layout/=MPI_SUCCESS)return
      call MPI_Bcast(wpw_support_w_ids,counts(2),MPI_INTEGER,0,dc%icomm_frag,ierr_layout)
      if(ierr_layout/=MPI_SUCCESS)return
      if(any(wpw_owned_w_ids<=0).or.any(wpw_support_w_ids<=0))return
      do i=2,size(wpw_owned_w_ids);if(wpw_owned_w_ids(i)<=wpw_owned_w_ids(i-1))return;enddo
      do i=2,size(wpw_support_w_ids);if(wpw_support_w_ids(i)<=wpw_support_w_ids(i-1))return;enddo
      layout_info=0
    end subroutine broadcast_wpw_w_row_layout

    subroutine assemble_wpw_core_volume_quadrature(volume_info)
      integer,intent(out)::volume_info
      integer::domain(3),local_lo(3),ix1,ix2,iy1,iy2,iz1,iz2,ixq,iyq,izq,ibq,ng,nps,owned_k(1)
      integer::grid(3),local_grid(3),point_info,owned_p_lo,owned_p_hi
      integer::local_failure,global_failure,ierr_quadrature,wpw_local_core_point_count,&
        core_descriptor_point,buffer_descriptor_point,descriptor_extent(3),astat
      integer,allocatable::column_ids(:)
      real(8),allocatable::chi(:),grad_chi(:,:)
      complex(8),allocatable::owned_w(:),owned_grad_w(:,:),support_w(:),support_grad_w(:,:)
      complex(8),allocatable::support_p(:),support_grad_p(:,:)
      real(8)::omega_cell

      volume_info=1;local_failure=0;ng=size(wpw_g_vectors,2);nps=size(wpw_support_fragment_ids)*ng
      if(ng<=0.or.nps<=0.or.size(wpw_owned_w_ids)<=0.or.size(wpw_support_w_ids)<=0)local_failure=1
      if(.not.wpw_occupied_w_basis%valid.or.&
        wpw_occupied_w_basis%local_count/=size(wpw_owned_w_ids))local_failure=1
      if(local_failure==0)then
        if(any(wpw_occupied_w_basis%owned_ids/=wpw_owned_w_ids))local_failure=1
      endif
      owned_k=findloc(wpw_support_fragment_ids,dc%i_frag)
      if(owned_k(1)<=0)local_failure=1
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      owned_p_lo=(owned_k(1)-1)*ng+1;owned_p_hi=owned_p_lo+ng-1
      astat=0;allocate(column_ids(nps),chi(size(wpw_support_fragment_ids)),&
        grad_chi(3,size(wpw_support_fragment_ids)),owned_w(size(wpw_owned_w_ids)),&
        owned_grad_w(3,size(wpw_owned_w_ids)),support_w(size(wpw_support_w_ids)),&
        support_grad_w(3,size(wpw_support_w_ids)),support_p(nps),support_grad_p(3,nps),stat=astat)
      local_failure=merge(0,1,astat==0)
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      call get_fragment_domain(dc,dc%i_frag,domain);local_lo=dc%ixyz_frag(:,dc%i_frag)
      ix1=max(mg%is(1),1);ix2=min(mg%ie(1),domain(1))
      iy1=max(mg%is(2),1);iy2=min(mg%ie(2),domain(2))
      iz1=max(mg%is(3),1);iz2=min(mg%ie(3),domain(3))
      wpw_local_core_point_count=max(0,ix2-ix1+1)*max(0,iy2-iy1+1)*max(0,iz2-iz1+1)
      call initialize_dg_wpw_volume_accumulator(wpw_volume_accumulator,size(owned_w),nps,nps,point_info,&
        point_capacity=wpw_local_core_point_count)
      if(point_info/=0)local_failure=1
      call initialize_dg_wpw_metric_accumulator(wpw_metric_accumulator,size(support_w),nps,point_info,&
        point_capacity=wpw_local_core_point_count)
      if(point_info/=0)local_failure=1
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      omega_cell=product(wpw_box_length)
      descriptor_extent=wpw_occupied_w_basis%buffer_hi-wpw_occupied_w_basis%buffer_lo+1
      if(info%id_o==0)then
      quadrature_z: do izq=iz1,iz2
        do iyq=iy1,iy2
          do ixq=ix1,ix2
            local_grid=[ixq,iyq,izq];grid=local_lo+local_grid-1
            buffer_descriptor_point=local_grid(1)-wpw_occupied_w_basis%buffer_lo(1)+1+&
              (local_grid(2)-wpw_occupied_w_basis%buffer_lo(2))*descriptor_extent(1)+&
              (local_grid(3)-wpw_occupied_w_basis%buffer_lo(3))*descriptor_extent(1)*descriptor_extent(2)
            if(buffer_descriptor_point<1.or.&
              buffer_descriptor_point>size(wpw_occupied_w_basis%buffer_values,1))then
              local_failure=1;exit quadrature_z
            endif
            if(any(local_grid<wpw_occupied_w_basis%gradient_lo).or.&
              any(local_grid>wpw_occupied_w_basis%gradient_hi))then
              local_failure=1;exit quadrature_z
            endif
            call evaluate_dg_wpw_occupied_w_point(wpw_occupied_w_basis,local_grid,&
              owned_w,owned_grad_w,point_info)
            if(point_info/=0)then;local_failure=1;exit quadrature_z;endif
            call evaluate_dg_wpw_core_w_support(wpw_owned_w_ids,owned_w,owned_grad_w,wpw_support_w_ids,&
              wpw_volume_halos,grid,1,support_w,support_grad_w,point_info,zero_outside_halo=.true.)
            if(point_info/=0)then;local_failure=1;exit quadrature_z;endif
            call evaluate_dg_wpw_core_point(wpw_support_core_lo,wpw_support_core_hi,&
              wpw_support_fragment_ids,dc%lg_tot%num,system%hgs,&
              wpw_window_buffer,wpw_window_width,wpw_g_vectors,grid,&
              omega_cell,column_ids,chi,grad_chi,support_p,support_grad_p,point_info)
            if(point_info/=0)then
              write(*,'(1x,a,i0,a,3(i0,1x),a,i0)')'[DG-WPW-LOCAL-FAIL] volume_core_point rank=',&
                dc%id_tot,' grid=',grid,' info=',point_info
              local_failure=1;exit quadrature_z
            endif
            if(any(column_ids(owned_p_lo:owned_p_hi)/=&
              [( (dc%i_frag-1)*ng+ibq,ibq=1,ng )]))then
              write(*,'(1x,a,i0,a,3(i0,1x))')'[DG-WPW-LOCAL-FAIL] volume_column_ids rank=',&
                dc%id_tot,' grid=',grid
              local_failure=1;exit quadrature_z
            endif
            call add_dg_wpw_core_point(wpw_volume_accumulator,owned_w,owned_grad_w,support_p,&
              support_grad_p,support_p,support_grad_p,V_local(1)%f(ixq,iyq,izq),hvol,point_info,&
              grid_id=wpw_core_global_grid_id(grid),density=rho_s(1)%f(ixq,iyq,izq))
            if(point_info/=0)then
              write(*,'(1x,a,i0,a,3(i0,1x),a,i0)')'[DG-WPW-LOCAL-FAIL] volume_accumulate rank=',&
                dc%id_tot,' grid=',grid,' info=',point_info
              local_failure=1;exit quadrature_z
            endif
            call add_dg_wpw_metric_point(wpw_metric_accumulator,support_w,support_p,hvol,point_info,&
              grid_id=wpw_core_global_grid_id(grid),w_grad=support_grad_w,p_grad=support_grad_p,&
              potential=V_local(1)%f(ixq,iyq,izq))
            if(point_info/=0)then;local_failure=1;exit quadrature_z;endif
          enddo
        enddo
      enddo quadrature_z
      endif
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      if(info%id_o==0)then
        if(wpw_core_p_accumulator%npoint/=wpw_volume_accumulator%npoint.or.&
           any(wpw_core_p_accumulator%grid_ids/=wpw_volume_accumulator%grid_ids).or.&
           maxval(abs(wpw_core_p_accumulator%pp_h-wpw_volume_accumulator%pp_h))>1d-12.or.&
           maxval(abs(wpw_core_p_accumulator%pp_s-wpw_volume_accumulator%pp_s))>1d-12)then
          local_failure=1
        endif
      endif
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      astat=0;allocate(wpw_volume_wp_h(size(owned_w),nps),wpw_volume_wp_s(size(owned_w),nps),&
        wpw_volume_pp_h(nps,nps),wpw_volume_pp_s(nps,nps),stat=astat)
      local_failure=merge(0,1,astat==0)
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      call build_dg_wpw_rank_local_quadrature(dc%icomm_frag,0,wpw_volume_accumulator,&
        wpw_volume_wp_h,wpw_volume_wp_s,wpw_volume_pp_h,wpw_volume_pp_s,volume_info)
      if(volume_info/=0)return
      astat=0;allocate(wpw_metric_ww_s(size(support_w),size(support_w)),&
        wpw_metric_wp_s(size(support_w),nps),wpw_metric_pp_s(nps,nps),stat=astat)
      local_failure=merge(0,1,astat==0)
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      call build_dg_wpw_metric_gram(dc%icomm_frag,0,wpw_metric_accumulator,&
        wpw_metric_ww_s,wpw_metric_wp_s,wpw_metric_pp_s,volume_info)
      astat=0;allocate(wpw_metric_ww_h(size(support_w),size(support_w)),&
        wpw_metric_wp_h(size(support_w),nps),wpw_metric_pp_h(nps,nps),&
        wpw_metric_ww_kinetic(size(support_w),size(support_w)),&
        wpw_metric_ww_potential(size(support_w),size(support_w)),stat=astat)
      local_failure=merge(0,1,astat==0)
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      if(volume_info==0)call build_dg_wpw_hamiltonian_volume(dc%icomm_frag,0,wpw_metric_accumulator,&
        wpw_metric_ww_h,wpw_metric_wp_h,wpw_metric_pp_h,volume_info)
      if(volume_info==0)call build_dg_wpw_ww_volume_components(dc%icomm_frag,0,wpw_metric_accumulator,&
        wpw_metric_ww_kinetic,wpw_metric_ww_potential,volume_info)
      if(volume_info==0.and.dc%id_frag==0)then
        if(maxval(abs(wpw_metric_pp_s-wpw_volume_pp_s))>1d-12.or.&
          maxval(abs(wpw_metric_pp_h-wpw_volume_pp_h))>1d-12)volume_info=1
      endif
      if(volume_info==0)wpw_nonowned_candidates_pending=.true.
    end subroutine assemble_wpw_core_volume_quadrature

    subroutine assemble_wpw_core_p_bootstrap(bootstrap_info)
      integer,intent(out)::bootstrap_info
      integer::domain(3),local_lo(3),ix1,ix2,iy1,iy2,iz1,iz2,ixq,iyq,izq,ng,nps,owned_k(1),ibq
      integer::grid(3),local_grid(3),point_info,owned_p_lo,owned_p_hi,local_failure,global_failure,&
        ierr_quadrature,local_core_point_count
      integer,allocatable::column_ids(:)
      real(8),allocatable::chi(:),grad_chi(:,:)
      complex(8),allocatable::support_p(:),support_grad_p(:,:)
      real(8)::omega_cell

      bootstrap_info=1;local_failure=0;ng=size(wpw_g_vectors,2)
      nps=size(wpw_support_fragment_ids)*ng
      if(ng<=0.or.nps<=0)local_failure=1
      owned_k=findloc(wpw_support_fragment_ids,dc%i_frag)
      if(owned_k(1)<=0)local_failure=1
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      owned_p_lo=(owned_k(1)-1)*ng+1;owned_p_hi=owned_p_lo+ng-1
      allocate(column_ids(nps),chi(size(wpw_support_fragment_ids)),&
        grad_chi(3,size(wpw_support_fragment_ids)),support_p(nps),support_grad_p(3,nps))
      call get_fragment_domain(dc,dc%i_frag,domain);local_lo=dc%ixyz_frag(:,dc%i_frag)
      ix1=max(mg%is(1),1);ix2=min(mg%ie(1),domain(1))
      iy1=max(mg%is(2),1);iy2=min(mg%ie(2),domain(2))
      iz1=max(mg%is(3),1);iz2=min(mg%ie(3),domain(3))
      local_core_point_count=max(0,ix2-ix1+1)*max(0,iy2-iy1+1)*max(0,iz2-iz1+1)
      call initialize_dg_wpw_core_p_accumulator(wpw_core_p_accumulator,nps,nps,point_info,&
        point_capacity=local_core_point_count)
      if(point_info/=0)local_failure=1
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      omega_cell=product(wpw_box_length)
      bootstrap_z: do izq=iz1,iz2;do iyq=iy1,iy2;do ixq=ix1,ix2
        local_grid=[ixq,iyq,izq];grid=local_lo+local_grid-1
        call evaluate_dg_wpw_core_point(wpw_support_core_lo,wpw_support_core_hi,&
          wpw_support_fragment_ids,dc%lg_tot%num,system%hgs,wpw_window_buffer,wpw_window_width,&
          wpw_g_vectors,grid,omega_cell,column_ids,chi,grad_chi,support_p,support_grad_p,point_info)
        if(point_info/=0.or.any(column_ids(owned_p_lo:owned_p_hi)/=&
          [((dc%i_frag-1)*ng+ibq,ibq=1,ng)]))then
          local_failure=1;exit bootstrap_z
        endif
        call add_dg_wpw_core_p_point(wpw_core_p_accumulator,support_p,support_grad_p,support_p,&
          support_grad_p,V_local(1)%f(ixq,iyq,izq),hvol,point_info,&
          grid_id=wpw_core_global_grid_id(grid),density=rho_s(1)%f(ixq,iyq,izq))
        if(point_info/=0)then;local_failure=1;exit bootstrap_z;endif
      enddo;enddo;enddo bootstrap_z
      call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,dc%icomm_frag,ierr_quadrature)
      if(ierr_quadrature/=MPI_SUCCESS.or.global_failure/=0)return
      bootstrap_info=0
    end subroutine assemble_wpw_core_p_bootstrap

    integer function wpw_core_global_grid_id(grid)result(grid_id)
      integer,intent(in)::grid(3)
      integer::wrapped(3)
      wrapped=modulo(grid-1,dc%lg_tot%num)+1
      grid_id=wrapped(1)+(wrapped(2)-1)*dc%lg_tot%num(1)+&
        (wrapped(3)-1)*dc%lg_tot%num(1)*dc%lg_tot%num(2)
    end function wpw_core_global_grid_id

    subroutine publish_wpw_core_volume_candidates(candidate_info)
      integer,intent(out)::candidate_info
      integer::nw,np_support,iwp,ipp,iww,iw,jw,ip,max_candidates,entry,owner_fragment,&
        local_bad,global_bad,ierr_publish,local_info,astat
      integer,allocatable::wp_w(:),wp_p(:),pp_r(:),pp_c(:)
      complex(8),allocatable::wp_h(:),wp_s(:),pp_h(:),pp_s(:)
      type(s_dg_wpw_staged_candidate),allocatable::staged(:),staged_conjugate(:)
      type(s_dg_wpw_owned_candidates)::owned,owned_conjugate
      candidate_info=1;local_bad=0
      if(dc%id_frag/=0.or..not.allocated(wpw_volume_wp_h).or..not.allocated(wpw_volume_pp_h).or.&
        .not.allocated(wpw_metric_ww_s).or..not.allocated(wpw_metric_wp_s).or.&
        .not.allocated(wpw_metric_pp_s))local_bad=1
      if(local_bad==0)then
        nw=size(wpw_support_w_ids);np_support=size(wpw_context%support_column_ids)
        if(any(shape(wpw_metric_ww_s)/=[nw,nw]).or.any(shape(wpw_metric_wp_s)/=[nw,np_support]).or.&
          any(shape(wpw_metric_pp_s)/=[np_support,np_support]).or.&
          any(shape(wpw_volume_wp_h)/=[size(wpw_owned_w_ids),np_support]).or.&
          any(shape(wpw_volume_pp_h)/=[np_support,np_support]))local_bad=1
      endif
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_publish)
      if(ierr_publish/=MPI_SUCCESS.or.global_bad/=0)return
      astat=0;allocate(wpw_support_w_owner(nw),stat=astat)
      local_bad=merge(0,1,astat==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_publish)
      if(ierr_publish/=MPI_SUCCESS.or.global_bad/=0)return
      wpw_support_w_owner=0
      do iw=1,nw
        if(any(wpw_owned_w_ids==wpw_support_w_ids(iw)))wpw_support_w_owner(iw)=dc%i_frag
        do ip=1,size(wpw_volume_halos)
          if(.not.wpw_volume_halos(ip)%valid.or..not.allocated(wpw_volume_halos(ip)%w_ids))cycle
          if(.not.any(wpw_volume_halos(ip)%w_ids==wpw_support_w_ids(iw)))cycle
          owner_fragment=wpw_volume_halos(ip)%source_rank+1
          if(wpw_support_w_owner(iw)==0)then;wpw_support_w_owner(iw)=owner_fragment
          else if(wpw_support_w_owner(iw)/=owner_fragment)then;wpw_support_w_owner(iw)=-1;endif
        enddo
      enddo
      local_bad=merge(0,1,all(wpw_support_w_owner>0))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_publish)
      if(ierr_publish/=MPI_SUCCESS.or.global_bad/=0)return
      astat=0
      allocate(wp_w(nw*np_support),wp_p(nw*np_support),wp_h(nw*np_support),wp_s(nw*np_support),stat=astat)
      local_bad=merge(0,1,astat==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_publish)
      if(ierr_publish/=MPI_SUCCESS.or.global_bad/=0)return
      iwp=0
      do ip=1,np_support
        do iw=1,nw
          iwp=iwp+1;wp_w(iwp)=wpw_support_w_ids(iw);wp_p(iwp)=wpw_context%support_column_ids(ip)
          wp_h(iwp)=wpw_metric_wp_h(iw,ip);wp_s(iwp)=wpw_metric_wp_s(iw,ip)
        enddo
      enddo
      astat=0
      allocate(pp_r(np_support*np_support),pp_c(np_support*np_support),&
        pp_h(np_support*np_support),pp_s(np_support*np_support),stat=astat)
      local_bad=merge(0,1,astat==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_publish)
      if(ierr_publish/=MPI_SUCCESS.or.global_bad/=0)return
      ipp=0
      do ip=1,np_support
        do iw=1,np_support
          ipp=ipp+1;pp_r(ipp)=wpw_context%support_column_ids(ip);pp_c(ipp)=wpw_context%support_column_ids(iw)
          pp_h(ipp)=wpw_metric_pp_h(ip,iw);pp_s(ipp)=wpw_metric_pp_s(ip,iw)
        enddo
      enddo
      astat=0;allocate(staged(nw*nw+size(wp_w)+size(pp_r)),stat=astat)
      local_bad=merge(0,1,astat==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_publish)
      if(ierr_publish/=MPI_SUCCESS.or.global_bad/=0)return
      entry=0;iww=0
      do jw=1,nw;do iw=1,nw
        entry=entry+1;iww=iww+1;staged(entry)%kind=wpw_candidate_kind_ww
        staged(entry)%image_id=iww;staged(entry)%target_fragment=wpw_support_w_owner(iw)
        staged(entry)%ww_r=wpw_support_w_ids(iw);staged(entry)%ww_c=wpw_support_w_ids(jw)
        staged(entry)%ww_s=wpw_metric_ww_s(iw,jw)
        staged(entry)%ww_h=wpw_metric_ww_h(iw,jw)
        staged(entry)%ww_kinetic=wpw_metric_ww_kinetic(iw,jw)
        staged(entry)%ww_potential=wpw_metric_ww_potential(iw,jw)
      enddo;enddo
      do iwp=1,size(wp_w)
        entry=entry+1;staged(entry)%kind=wpw_candidate_kind_wp;staged(entry)%image_id=entry
        staged(entry)%wp_w=wp_w(iwp);staged(entry)%wp_p=wp_p(iwp)
        staged(entry)%wp_h=wp_h(iwp);staged(entry)%wp_s=wp_s(iwp)
      enddo
      do ipp=1,size(pp_r)
        entry=entry+1;staged(entry)%kind=wpw_candidate_kind_pp;staged(entry)%image_id=entry
        staged(entry)%pp_r=pp_r(ipp);staged(entry)%pp_c=pp_c(ipp)
        staged(entry)%pp_h=pp_h(ipp);staged(entry)%pp_s=pp_s(ipp)
      enddo
      max_candidates=size(staged)
      call route_dg_wpw_candidate_halo(wpw_production_comm,1,dc%n_frag,size(wpw_g_indices,2),&
        wpw_support_fragment_ids,staged,owned,candidate_info,max_candidates)
      if(candidate_info/=0)return
      astat=0;allocate(staged_conjugate(nw*nw),stat=astat)
      local_bad=merge(0,1,astat==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_publish)
      if(ierr_publish/=MPI_SUCCESS.or.global_bad/=0)return
      entry=0
      do jw=1,nw;do iw=1,nw
        entry=entry+1;staged_conjugate(entry)%kind=wpw_candidate_kind_ww
        staged_conjugate(entry)%image_id=entry;staged_conjugate(entry)%target_fragment=wpw_support_w_owner(jw)
        staged_conjugate(entry)%ww_r=wpw_support_w_ids(jw)
        staged_conjugate(entry)%ww_c=wpw_support_w_ids(iw)
        staged_conjugate(entry)%ww_s=conjg(wpw_metric_ww_s(iw,jw))
        staged_conjugate(entry)%ww_h=conjg(wpw_metric_ww_h(iw,jw))
        staged_conjugate(entry)%ww_kinetic=conjg(wpw_metric_ww_kinetic(iw,jw))
        staged_conjugate(entry)%ww_potential=conjg(wpw_metric_ww_potential(iw,jw))
      enddo;enddo
      call route_dg_wpw_candidate_halo(wpw_production_comm,1,dc%n_frag,size(wpw_g_indices,2),&
        wpw_support_fragment_ids,staged_conjugate,owned_conjugate,local_info,max_candidates)
      candidate_info=max(candidate_info,local_info)
      if(candidate_info/=0)return
      local_bad=0
      if(size(owned%ww_r)/=size(owned_conjugate%ww_r))then;local_bad=1
      else if(any(owned%ww_r/=owned_conjugate%ww_r).or.any(owned%ww_c/=owned_conjugate%ww_c))then;local_bad=1;endif
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,wpw_production_comm,ierr_publish)
      if(ierr_publish/=MPI_SUCCESS.or.global_bad/=0)return
      call publish_dg_wpw_realspace_metric(wpw_context,owned%wp_w,owned%wp_p,owned%wp_h,owned%wp_s,&
        owned%pp_r,owned%pp_c,owned%pp_h,owned%pp_s,owned%ww_r,owned%ww_c,owned%ww_s,&
        owned_conjugate%ww_s,owned%ww_kinetic,owned_conjugate%ww_kinetic,&
        owned%ww_potential,owned_conjugate%ww_potential,'tail_carrying_realspace_gram',candidate_info)
      if(candidate_info==0)then
        wpw_context%ww_provenance_fingerprint=wpw_occupied_w_basis%fingerprint
        wpw_nonowned_candidates_pending=.false.
      endif
    end subroutine publish_wpw_core_volume_candidates

    subroutine prepare_wpw_canonical_face_trace_provider(trace_info)
      integer,intent(out)::trace_info
      integer::iface,axis,tangent(2),domain(3),local_lo(3),local_grid(3),grid_point(3)
      integer::it,jt,point,npoint,point_info,ng,np,ib,epoch,owned_count,w_fail,p_fail
      integer,allocatable::grid(:,:),coverage(:),column_ids(:),p_ids(:)
      real(8),allocatable::chi(:),grad_chi(:,:)
      complex(8),allocatable::owned_w(:),owned_grad_w(:,:),side_w(:,:),side_grad_w(:,:,:)
      complex(8),allocatable::side_p(:,:),side_grad_p(:,:,:),point_p(:),point_grad_p(:,:)
      real(8)::omega_cell

      trace_info=1;epoch=1;ng=size(wpw_g_vectors,2);np=size(wpw_support_fragment_ids)*ng
      if(ng<=0.or.np<=0.or..not.allocated(wpw_canonical_faces))return
      allocate(p_ids(np),column_ids(np),chi(size(wpw_support_fragment_ids)),&
        grad_chi(3,size(wpw_support_fragment_ids)))
      do point=1,np
        p_ids(point)=(wpw_support_fragment_ids((point-1)/ng+1)-1)*ng+modulo(point-1,ng)+1
      enddo
      allocate(wpw_face_sides(size(wpw_canonical_faces)))
      call get_fragment_domain(dc,dc%i_frag,domain);local_lo=dc%ixyz_frag(:,dc%i_frag)
      omega_cell=product(wpw_box_length)
      do iface=1,size(wpw_canonical_faces)
        axis=wpw_canonical_faces(iface)%axis
        select case(axis);case(1);tangent=[2,3];case(2);tangent=[1,3];case(3);tangent=[1,2];case default;return;end select
        npoint=domain(tangent(1))*domain(tangent(2))
        allocate(grid(3,npoint),coverage(npoint),side_w(size(wpw_support_w_ids),npoint),&
          side_grad_w(3,size(wpw_support_w_ids),npoint),side_p(np,npoint),side_grad_p(3,np,npoint))
        allocate(owned_w(size(wpw_owned_w_ids)),owned_grad_w(3,size(wpw_owned_w_ids)),&
          point_p(np),point_grad_p(3,np))
        coverage=0;side_w=0;side_grad_w=0;side_p=0;side_grad_p=0;point=0
        owned_count=0;w_fail=0;p_fail=0
        do jt=1,domain(tangent(2));do it=1,domain(tangent(1))
          point=point+1;local_grid=1
          local_grid(axis)=wpw_trace_face_coord(axis,wpw_canonical_faces(iface)%side_from_local)
          local_grid(tangent(1))=it;local_grid(tangent(2))=jt;grid_point=local_lo+local_grid-1
          grid(:,point)=grid_point
          if(any(local_grid<mg%is).or.any(local_grid>mg%ie))cycle
          owned_count=owned_count+1;coverage(point)=1
          call evaluate_dg_wpw_occupied_w_point(wpw_occupied_w_basis,local_grid,&
            owned_w,owned_grad_w,point_info)
          if(point_info/=0)then;w_fail=w_fail+1;coverage(point)=0;cycle;endif
          call evaluate_dg_wpw_core_w_support(wpw_owned_w_ids,owned_w,owned_grad_w,wpw_support_w_ids,&
            wpw_volume_halos,grid_point,epoch,side_w(:,point),side_grad_w(:,:,point),point_info,&
            zero_outside_halo=.true.)
          if(point_info/=0)then;w_fail=w_fail+1;coverage(point)=0;cycle;endif
          call evaluate_dg_wpw_core_point(wpw_support_core_lo,wpw_support_core_hi,&
            wpw_support_fragment_ids,dc%lg_tot%num,system%hgs,&
            wpw_window_buffer,wpw_window_width,wpw_g_vectors,grid_point,&
            omega_cell,column_ids,chi,grad_chi,point_p,point_grad_p,point_info)
          if(point_info/=0.or.any(column_ids/=p_ids))then;p_fail=p_fail+1;coverage(point)=0;cycle;endif
          side_p(:,point)=point_p;side_grad_p(:,:,point)=point_grad_p
        enddo;enddo
        call reduce_dg_wpw_face_side_parts(dc%icomm_frag,0,wpw_canonical_faces(iface)%neighbor_fragment-1,&
          epoch,dc%i_frag,wpw_canonical_faces(iface)%neighbor_fragment,axis,&
          wpw_canonical_faces(iface)%side_from_local,wpw_support_w_ids,p_ids,grid,coverage,&
          side_w,side_grad_w,side_p,side_grad_p,wpw_face_sides(iface),point_info)
        if(point_info/=0)then
          write(*,'(1x,a,i0,8(a,i0),a,3(i0,1x),a,3(i0,1x))')&
            '[DG-WPW-LOCAL-FAIL] face_reduce rank=',dc%id_tot,&
            ' iface=',iface,' axis=',axis,' side=',wpw_canonical_faces(iface)%side_from_local,&
            ' covered=',count(coverage==1),' npoint=',npoint,' owned=',owned_count,&
            ' w_fail=',w_fail,' p_fail=',p_fail,' mg_is=',mg%is,' mg_ie=',mg%ie
          return
        endif
        deallocate(grid,coverage,side_w,side_grad_w,side_p,side_grad_p,owned_w,owned_grad_w,point_p,point_grad_p)
      enddo
      if(dc%id_frag==0)then
        if(size(wpw_face_sides)>0)then
          call exchange_dg_wpw_face_side_schedule(wpw_production_comm,epoch,wpw_face_sides,&
            wpw_remote_face_sides,point_info)
          if(point_info/=0)then
            write(*,'(1x,a,i0,2(a,i0))')'[DG-WPW-LOCAL-FAIL] face_exchange rank=',dc%id_tot,&
              ' info=',point_info,' nface=',size(wpw_face_sides)
            return
          endif
          wpw_context%halo_epoch=epoch
          do iface=1,size(wpw_face_sides)
            if(.not.wpw_canonical_faces(iface)%canonical_owner)cycle
            call prepare_dg_wpw_trace_halo(wpw_context,wpw_trace_set,wpw_trace_provider,&
              wpw_face_sides(iface),wpw_remote_face_sides(iface),point_info)
            if(point_info/=0)then
              write(*,'(1x,a,i0,4(a,i0))')'[DG-WPW-LOCAL-FAIL] face_bind rank=',dc%id_tot,&
                ' iface=',iface,' neighbor=',wpw_canonical_faces(iface)%neighbor_fragment,&
                ' axis=',wpw_canonical_faces(iface)%axis,&
                ' side=',wpw_canonical_faces(iface)%side_from_local
              return
            endif
          enddo
        else
          allocate(wpw_remote_face_sides(0));wpw_context%halo_epoch=epoch
        endif
      endif
      trace_info=0
    end subroutine prepare_wpw_canonical_face_trace_provider

    subroutine scan_wpw_canonical_faces(scan_info)
      integer,intent(out)::scan_info
      integer::nface,iface,jface,domain(3),local_lo(3),local_hi(3),max_candidates
      integer,allocatable::kminus(:),kplus(:),axes(:),sides(:)
      integer,allocatable::minus_lo(:,:),minus_hi(:,:),plus_lo(:,:),plus_hi(:,:)
      scan_info=1
      if(dc%id_frag/=0.or.wpw_production_comm==MPI_COMM_NULL)return
      nface=count(wpw_canonical_faces%canonical_owner)
      allocate(kminus(nface),kplus(nface),axes(nface),sides(nface),minus_lo(3,nface),&
        minus_hi(3,nface),plus_lo(3,nface),plus_hi(3,nface))
      call get_fragment_domain(dc,dc%i_frag,domain);local_lo=dc%ixyz_frag(:,dc%i_frag)
      local_hi=local_lo+domain-1;jface=0
      do iface=1,size(wpw_canonical_faces)
        if(.not.wpw_canonical_faces(iface)%canonical_owner)cycle
        jface=jface+1;kminus(jface)=dc%i_frag
        kplus(jface)=wpw_canonical_faces(iface)%neighbor_fragment
        axes(jface)=wpw_canonical_faces(iface)%axis
        sides(jface)=wpw_canonical_faces(iface)%side_from_local
        minus_lo(:,jface)=local_lo;minus_hi(:,jface)=local_hi
        plus_lo(:,jface)=local_lo+wpw_canonical_faces(iface)%neighbor_displacement*domain
        plus_hi(:,jface)=plus_lo(:,jface)+domain-1
        if(sides(jface)>0)minus_hi(axes(jface),jface)=minus_hi(axes(jface),jface)+1
      enddo
      max_candidates=max(1,nface*size(wpw_support_w_ids)*(&
        size(wpw_context%support_column_ids)+size(wpw_support_w_ids)))
      call scan_dg_wpw_canonical_faces(wpw_context,1,wpw_support_fragment_ids,wpw_trace_provider,&
        kminus,kplus,axes,sides,minus_lo,minus_hi,plus_lo,plus_hi,system%hgs,&
        wpw_support_w_ids,wpw_support_w_owner,wpw_context%support_column_ids,max_candidates,scan_info)
      if(scan_info/=0)write(*,'(1x,a,i0,3(a,i0))')'[DG-WPW-LOCAL-FAIL] face_scan rank=',dc%id_tot,&
        ' info=',scan_info,' nface=',nface,' max_candidates=',max_candidates
    end subroutine scan_wpw_canonical_faces

    subroutine build_wpw_support_geometry(records,fragment_ids,core_lo,core_hi,overlap_lo,overlap_hi,support_info)
      use salmon_global,only:num_fragment
      type(s_dg_wpw_support_record),allocatable,intent(out)::records(:)
      integer,allocatable,intent(out)::fragment_ids(:),core_lo(:,:),core_hi(:,:)
      integer,allocatable,intent(out)::overlap_lo(:,:),overlap_hi(:,:)
      integer,intent(out)::support_info
      integer::extent(3),domain(3),k,local_lo(3),local_hi(3),buffer(3)
      extent=merge(1,0,num_fragment>1 .and. wpw_window_buffer>0)
      call build_dg_wpw_fragment_support(dc%i_frag,num_fragment,extent,records,fragment_ids,support_info)
      if(support_info/=0)return
      allocate(core_lo(3,size(fragment_ids)),core_hi(3,size(fragment_ids)))
      do k=1,size(fragment_ids)
        call get_fragment_domain(dc,fragment_ids(k),domain)
        core_lo(:,k)=dc%ixyz_frag(:,fragment_ids(k))
        core_hi(:,k)=core_lo(:,k)+domain-1
      enddo
      call get_fragment_domain(dc,dc%i_frag,domain)
      local_lo=dc%ixyz_frag(:,dc%i_frag);local_hi=local_lo+domain-1
      buffer=wpw_window_buffer
      allocate(overlap_lo(3,size(records)),overlap_hi(3,size(records)))
      do k=1,size(records)
        call dg_wpw_support_overlap_box(local_lo,local_hi,records(k)%image_shift,buffer,&
          overlap_lo(:,k),overlap_hi(:,k),support_info)
        if(support_info/=0)return
      enddo
    end subroutine build_wpw_support_geometry

    subroutine build_wpw_w_row_layout(owned_ids,support_ids,layout_info)
      integer,allocatable,intent(out)::owned_ids(:),support_ids(:)
      integer,intent(out)::layout_info
      if(dc%id_frag/=0)then
        allocate(owned_ids(0),support_ids(0));layout_info=0;return
      endif
      if(.not.wpw_occupied_w_basis%valid.or.wpw_occupied_w_basis%local_fragment/=dc%i_frag)then
        layout_info=1;return
      endif
      call build_dg_wpw_w_row_layout_from_owned_ids(wpw_production_comm,dc%i_frag,&
        wpw_support_fragment_ids,wpw_occupied_w_basis%owned_ids,owned_ids,support_ids,layout_info)
    end subroutine

    subroutine check_lcfo_basis_potential_inputs_finite
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      integer :: ix0,iy0,iz0,is0,io0

      do io0=1,size(f_basis,5); do is0=1,size(f_basis,4); do iz0=1,size(f_basis,3)
        do iy0=1,size(f_basis,2); do ix0=1,size(f_basis,1)
          if(.not.ieee_is_finite(f_basis(ix0,iy0,iz0,is0,io0))) then
            write(*,'(1x,a,i0,5(a,i0))') '[DC-LCFO-INPUT-LOCAL] label=basis rank=',dc%id_tot, &
              ' ix=',ix0,' iy=',iy0,' iz=',iz0,' spin=',is0,' state=',io0
            stop 'DC-LCFO rank-local basis is non-finite'
          end if
        end do; end do
      end do; end do; end do
      do is0=1,nspin
        do iz0=lbound(V_local(is0)%f,3),ubound(V_local(is0)%f,3)
          do iy0=lbound(V_local(is0)%f,2),ubound(V_local(is0)%f,2)
            do ix0=lbound(V_local(is0)%f,1),ubound(V_local(is0)%f,1)
              if(.not.ieee_is_finite(V_local(is0)%f(ix0,iy0,iz0))) then
                write(*,'(1x,a,i0,4(a,i0))') '[DC-LCFO-INPUT-LOCAL] label=vlocal rank=',dc%id_tot, &
                  ' ix=',ix0,' iy=',iy0,' iz=',iz0,' spin=',is0
                stop 'DC-LCFO rank-local potential is non-finite'
              end if
            end do
          end do
        end do
      end do
    end subroutine check_lcfo_basis_potential_inputs_finite

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

    subroutine import_wpw_lcfo_ww_components(import_info)
      integer,intent(out)::import_info
      integer::iface,ih,io_local,io_remote,entry,nentry,nlocal,side,face_key
      integer,allocatable::face_id(:),row_id(:),col_id(:),axis_id(:),side_id(:),image(:,:)
      real(8),allocatable::value(:)
      logical::found
      call validate_dg_wpw_surface_self_policy(use_surface_self_hamiltonian_mode(),import_info)
      if(import_info/=0)return
      import_info=1;nlocal=n_basis(dc%i_frag,1);nentry=0
      do iface=1,size(wpw_canonical_faces)
        if(.not.wpw_canonical_faces(iface)%canonical_owner)cycle
        nentry=nentry+nlocal*n_basis(wpw_canonical_faces(iface)%neighbor_fragment,1)
      enddo
      allocate(face_id(nentry),row_id(nentry),col_id(nentry),axis_id(nentry),side_id(nentry),&
        image(3,nentry),value(nentry));entry=0
      do iface=1,size(wpw_canonical_faces)
        if(.not.wpw_canonical_faces(iface)%canonical_owner)cycle
        side=wpw_canonical_faces(iface)%side_from_local;found=.false.
        do ih=1,n_halo
          if(halo(ih)%ifrag_src/=wpw_canonical_faces(iface)%neighbor_fragment.or.&
             halo(ih)%axis/=wpw_canonical_faces(iface)%axis.or.-halo(ih)%dvec(halo(ih)%axis)/=side)cycle
          found=.true.
          face_key=(dc%i_frag-1)*6+(wpw_canonical_faces(iface)%axis-1)*2+merge(1,2,side<0)
          do io_local=1,nlocal
            do io_remote=1,n_basis(halo(ih)%ifrag_src,1)
              entry=entry+1;face_id(entry)=face_key
              row_id(entry)=index_basis(io_remote,halo(ih)%ifrag_src,1)
              col_id(entry)=index_basis(io_local,dc%i_frag,1)
              axis_id(entry)=halo(ih)%axis;side_id(entry)=side
              image(:,entry)=wpw_canonical_faces(iface)%periodic_shift
              value(entry)=halo(ih)%mat_H_surface_cross(io_remote,io_local,1)
            enddo
          enddo
          exit
        enddo
        if(.not.found)return
      enddo
      if(entry/=nentry)return
      call import_dg_wpw_lcfo_ww_components(dc%i_frag,n_basis(:,1),index_basis(:,:,1),&
        mat_H_weak_kinetic(1:nlocal,1:nlocal,1),mat_H_weak_potential(1:nlocal,1:nlocal,1),&
        mat_H_weak_nonlocal(1:nlocal,1:nlocal,1),mat_H_surface_self(1:nlocal,1:nlocal,1),&
        face_id,row_id,col_id,axis_id,side_id,image,value,'orthonormal_ww',wpw_ww_components,import_info)
      if(import_info/=0)return
      call publish_dg_wpw_lcfo_ww_components(wpw_context,wpw_ww_components,import_info)
    end subroutine import_wpw_lcfo_ww_components

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
      call check_lcfo_h_component_finite(mat_H_weak_kinetic,'volume_kinetic')
      call check_lcfo_h_component_finite(mat_H_weak_potential,'volume_local_potential')
      call check_lcfo_h_component_finite(mat_H_weak_nonlocal,'volume_nonlocal')
      call check_lcfo_h_component_finite(mat_H_volume_weak_local,'volume_weak_total')
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
          call check_lcfo_h_component_finite(halo(i_halo)%mat_H_surface_cross, &
            'surface_cross')
          if(dc%id_tot==0) write(*,*) "Halo communication #",i_halo,": done"
        end do
        call check_lcfo_h_component_finite(mat_H_surface_self,'surface_self')
        call check_lcfo_h_component_finite(mat_H_local,'assembled_local')
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

    subroutine check_lcfo_h_component_finite(matrix,label)
      implicit none
      real(8),intent(in)::matrix(:,:,:)
      character(*),intent(in)::label
      integer::i,j,spin

      if(all(ieee_is_finite(matrix)))return
      do spin=1,size(matrix,3)
        do j=1,size(matrix,2)
          do i=1,size(matrix,1)
            if(ieee_is_finite(matrix(i,j,spin)))cycle
            write(*,'(1x,3a,i0,3(a,i0),a,es24.16)') &
              '[DC-LCFO-H-COMPONENT-LOCAL] label=',trim(label),' rank=',dc%id_tot, &
              ' i=',i,' j=',j,' spin=',spin,' value=',matrix(i,j,spin)
            stop 'DC-LCFO rank-local Hamiltonian component is non-finite'
          end do
        end do
      end do
    end subroutine check_lcfo_h_component_finite

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

    integer function wpw_trace_face_coord(axis, side) result(idx)
      implicit none
      integer,intent(in) :: axis, side

      ! A canonical trace is sampled at one shared unwrapped grid point.  The
      ! plus trace therefore comes from the first buffered point beyond the
      ! local core, which is the neighbor's first core point.  Sampling the
      ! local core endpoint here would compare adjacent points as if they were
      ! the two limits of one face.
      if(side > 0) then
        idx = dc%nxyz_domain(axis) + 1
      else
        idx = 1
      end if
    end function wpw_trace_face_coord

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

        if(any(.not.ieee_is_finite(h_div))) then
          do iy_loc=1,size(h_div,2)
            do ix_loc=1,size(h_div,1)
              if(ieee_is_finite(h_div(ix_loc,iy_loc))) cycle
              write(*,'(1x,a,i0,2(a,i0),a,es24.16)') &
                '[DC-LCFO-EIGENEXA-LOCAL-H] rank=',dc%id_tot, &
                ' local_i=',ix_loc,' local_j=',iy_loc,' value=',h_div(ix_loc,iy_loc)
              stop 'DC-LCFO rank-local Hamiltonian is non-finite before EigenExa'
            end do
          end do
        end if

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
        call write_wannier90_global_bpw_seed()
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

    subroutine build_local_bond_center_projection_map(nproj, center_bohr,atom_ids,bond_images)
      implicit none
      integer, intent(out) :: nproj
      real(8), allocatable, intent(out) :: center_bohr(:,:)
      integer,allocatable,optional,intent(out)::atom_ids(:,:),bond_images(:,:)
      real(8), allocatable :: all_center_bohr(:,:)
      integer,allocatable::all_atom_ids(:,:),all_bond_images(:,:)
      integer :: iw, ip, nall

      nall = max(1, count_bond_center_candidates())
      call build_bond_center_projection_map(nall,all_center_bohr,all_atom_ids,all_bond_images)
      nproj = 0
      do iw=1,nall
        if(find_owner_fragment_from_center(all_center_bohr(1:3,iw)) == dc%i_frag) nproj = nproj + 1
      end do
      if(nproj <= 0) then
        allocate(center_bohr(3,1))
        center_bohr = 0.0d0
        if(present(atom_ids))allocate(atom_ids(2,1),source=0)
        if(present(bond_images))allocate(bond_images(3,1),source=0)
      else
        allocate(center_bohr(3,nproj))
        if(present(atom_ids))allocate(atom_ids(2,nproj))
        if(present(bond_images))allocate(bond_images(3,nproj))
        ip = 0
        do iw=1,nall
          if(find_owner_fragment_from_center(all_center_bohr(1:3,iw)) /= dc%i_frag) cycle
          ip = ip + 1
          center_bohr(1:3,ip) = all_center_bohr(1:3,iw)
          if(present(atom_ids))atom_ids(:,ip)=all_atom_ids(:,iw)
          if(present(bond_images))bond_images(:,ip)=all_bond_images(:,iw)
        end do
      end if
      if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0)') &
        "[DC-LCFO-LOCAL-WANNIER] fragment=", dc%i_frag, " local bond-center seeds=", nproj
      deallocate(all_center_bohr,all_atom_ids,all_bond_images)
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
      symmetry_filename='<auto>'
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
      if(.not.local_ok) write(*,'(1x,a,i0,3(a))') &
        '[FATAL] SAWF closed-seed symmetry load rank=',dc%id_tot, &
        ' source=',trim(symmetry_filename),' reason='//trim(message)
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
      if(.not.local_ok) write(*,'(1x,a,i0,a,i0,2a)') &
        '[FATAL] SAWF closed-seed basis rank=',dc%id_tot,' fragment=',dc%i_frag, &
        ' reason=',trim(message)
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
      use communication, only: comm_get_max,comm_bcast,comm_summation,comm_sync_all
      use filesystem, only: atomic_create_directory
      use salmon_global, only: izatom, sysname, wannier_num_wann, &
        wannier_site_symmetry, wannier_symmetry_file, wannier_symmetry_tolerance, &
        wannier_sawf_generation,wannier_sawf_structure_class,wannier_sawf_gauge_tolerance, &
        wannier_sawf_hamiltonian_tolerance
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
      integer,allocatable::sawf_source_channels(:),sawf_core_materialize_map(:), &
        sawf_buffer_materialize_map(:),sawf_gauge_parent(:),sawf_parent_shared_map(:), &
        sawf_child_shared_map(:)
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
      real(8) :: hamiltonian_tolerance
      real(8) :: max_grid_residual,center_grid(3)
      real(8), allocatable :: fractional_positions(:,:),d_wann_real(:,:,:)
      real(8), allocatable :: sawf_vacuum_fraction(:)
      complex(8), allocatable :: d_band_local(:,:),d_band_sum(:,:),d_wann(:,:),amn(:,:)
      complex(8), allocatable :: d_band_set(:,:,:),d_wann_set(:,:,:)
      complex(8),allocatable :: sawf_local_d_wann(:,:,:)
      complex(8),allocatable :: sawf_local_d_band(:,:,:),sawf_local_states(:,:)
      complex(8),allocatable :: sawf_local_amn(:,:)
      complex(8),allocatable :: sawf_local_v_matrix(:,:),sawf_representative_orbitals(:,:)
      complex(8),allocatable :: sawf_representative_buffer_orbitals(:,:)
      complex(8),allocatable::sawf_materialize_d_wann(:,:)
      real(8),allocatable :: sawf_local_energy(:)
      real(8),allocatable::sawf_local_coefficients(:,:)
      real(8),allocatable :: sawf_local_centers(:,:),sawf_local_spreads(:)
      real(8) :: representation_residual
      real(8)::sawf_representation_max,sawf_local_representation_local(1), &
        sawf_local_representation_global(1)
      real(8)::sawf_gauge_closure_residual,sawf_gauge_alignment_residual
      real(8)::sawf_gauge_residual_local(2),sawf_gauge_residual_global(2)
      type(t_sawf_projection_channel), allocatable :: channels(:)
      type(t_sawf_symop), allocatable :: symmetry_operations(:)
      type(t_sawf_dmn_writer) :: writer
      type(t_sawf_dmn_writer)::sawf_local_writer
      type(t_sawf_template_checkpoint)::sawf_template
      type(t_sawf_template_fingerprint)::sawf_template_fingerprint
      type(t_sawf_ragged_local_basis)::sawf_materialized_basis
      type(t_sawf_ragged_local_basis)::sawf_parent_basis
      type(t_sawf_fragment_state_cache) :: state_cache
      type(t_sawf_closed_basis) :: closed_basis
      type(t_sawf_environment_receipt),allocatable :: sawf_environment_receipts(:)
      type(t_sawf_seed_bundle),allocatable :: sawf_seed_bundles(:)
      logical :: local_ok,grid_map_ok,fragment_map_ok,center_available,split_fragment_global_mode
      logical::sawf_template_reuse
      logical, allocatable :: sawf_environment_equivalent(:,:),sawf_defect_intersects(:), &
        sawf_regenerate_environment(:),sawf_generate_independently(:),sawf_inside_atom(:)
      integer :: max_targets_per_source,local_left,local_right,local_relation,global_relation, &
        num_bands_chk,num_wann_chk,representative,materialize_operation,parent_fragment,parent_representative
      character(512) :: message
      character(256) :: symmetry_filename,allocation_message,dmn_filename,amn_filename, &
        sawf_supercell_fingerprint
      character(512)::local_chk_filename,local_basis_filename,parent_basis_filename

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
        if(dc%id_frag==0)then
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
      if(split_fragment_global_mode) call lcfo_sawf_fatal( &
        'SAWF split-fragment symmetry requires a constructed global symmetry-closed basis')
      closure_tolerance=max(1.0d-10,wannier_symmetry_tolerance)
      hamiltonian_tolerance=closure_tolerance
      if(wannier_sawf_hamiltonian_tolerance>0d0) &
        hamiltonian_tolerance=max(closure_tolerance,wannier_sawf_hamiltonian_tolerance)
      failure=0; message=''
      if(dc%id_tot == 0) then
        if(local_ok) then
          amn_filename=trim(dc%base_directory)//trim(sysname)//'.amn'
          call read_sawf_amn_matrix(amn_filename,nband_wann,wannier_num_wann,amn,local_ok,message)
        end if
        if(local_ok) then
          dmn_filename=trim(dc%base_directory)//trim(sysname)//'.dmn'
          call begin_sawf_dmn(writer,dmn_filename,nband_wann,wannier_num_wann, &
            size(symmetry_operations),closure_tolerance,local_ok,message, &
            hamiltonian_tolerance=hamiltonian_tolerance)
        end if
        failure=merge(0,1,local_ok)
      end if
      call comm_bcast(failure,dc%icomm_tot,0)
      call comm_bcast(message,dc%icomm_tot,0)
      if(failure /= 0) then
        if(dc%id_tot == 0) call abort_sawf_dmn(writer)
        call lcfo_sawf_fatal('SAWF DMN initialization failed: '//trim(message))
      end if

      allocate(d_wann_set(wannier_num_wann,wannier_num_wann,size(symmetry_operations)), &
        stat=allocation_status)
      if(allocation_status/=0)call lcfo_sawf_fatal('SAWF D_wann representation-set allocation failed')
      if(dc%id_tot==0) then
        allocate(d_band_set(nband_wann,nband_wann,size(symmetry_operations)),stat=allocation_status)
        if(allocation_status/=0) call lcfo_sawf_fatal('SAWF representation-set allocation failed')
      end if
      local_ok=.true.; message=''
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
          local_energy_out=sawf_local_energy,local_amn_out=sawf_local_amn, &
          local_coeff_out=sawf_local_coefficients)
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

      failure=0; message=''; representation_residual=0d0;sawf_representation_max=0d0
      if(dc%id_tot==0) then
        call validate_sawf_operation_set_products(d_band_set,product_left,product_right, &
          product_result,closure_tolerance,representation_residual,local_ok,message)
        if(local_ok)sawf_representation_max=max(sawf_representation_max,representation_residual)
        if(local_ok)then
          call validate_sawf_operation_set_products(d_wann_set,product_left,product_right, &
            product_result,closure_tolerance,representation_residual,local_ok,message)
          if(local_ok)sawf_representation_max=max(sawf_representation_max,representation_residual)
        end if
        failure=merge(0,1,local_ok)
      end if
      call comm_bcast(failure,dc%icomm_tot,0); call comm_bcast(message,dc%icomm_tot,0)
      call comm_bcast(sawf_representation_max,dc%icomm_tot,0)
      if(failure/=0) then
        if(dc%id_tot==0) call abort_sawf_dmn(writer)
        call lcfo_sawf_fatal('SAWF D_band/D_wann group representation validation failed: '//trim(message))
      end if
      if(dc%id_tot==0) write(*,'(1x,a,es13.5)') &
        '[DC-LCFO-SAWF-GROUP] max_representation_residual=',sawf_representation_max
      sawf_local_representation_local=0d0;sawf_local_representation_global=0d0
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
        if(local_ok)sawf_local_representation_local(1)=max( &
          sawf_local_representation_local(1),representation_residual)
        if(local_ok)then
          call validate_sawf_operation_set_products(sawf_local_d_wann, &
            sawf_local_product_left,sawf_local_product_right,sawf_local_product_result, &
            closure_tolerance,representation_residual,local_ok,message)
          if(local_ok)sawf_local_representation_local(1)=max( &
            sawf_local_representation_local(1),representation_residual)
        end if
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
          size(sawf_selected_channels),size(sawf_local_stabilizer),closure_tolerance,local_ok,message, &
          hamiltonian_tolerance=hamiltonian_tolerance)
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
      if(trim(wannier_sawf_generation)=='hierarchical')then
        call comm_get_max(sawf_local_representation_local,sawf_local_representation_global,1,dc%icomm_tot)
        if(dc%id_tot==0)write(*,'(1x,a,a,es13.5)')'[DC-LCFO-SAWF-SYMMETRY-SUMMARY]', &
          ' D_band_D_wann_closure_max=',sawf_local_representation_global(1)
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
          allocate(sawf_representative_buffer_orbitals(size(state_cache%source%buffer_basis,1),num_wann_chk))
          sawf_representative_buffer_orbitals=matmul(cmplx(state_cache%source%buffer_basis,0d0,kind=8), &
            matmul(cmplx(sawf_local_coefficients,0d0,kind=8),sawf_local_v_matrix))
          sawf_template_fingerprint%geometry=sawf_supercell_fingerprint
          sawf_template_fingerprint%pseudopotential=sawf_supercell_fingerprint
          sawf_template_fingerprint%grid=sawf_environment_key(dc%i_frag)
          write(sawf_template_fingerprint%band_window,'(a,i0,a,i0)')'bands=',num_bands_chk,':wann=',num_wann_chk
          write(sawf_template_fingerprint%complete_projection_shell,'(a,i0)') &
            'selected_channels=',size(sawf_selected_channels)
          write(sawf_template_fingerprint%symmetry,'(a,i0)')'actual_stabilizer=',size(sawf_local_stabilizer)
          write(sawf_template_fingerprint%buffer,'(3(i0,1x))')state_cache%source%buffer_width
          sawf_template_fingerprint%generator='SALMON hierarchical SAWF template schema 2'
          sawf_template%fingerprint=sawf_template_fingerprint
          allocate(sawf_template%centers(3,num_wann_chk),sawf_template%spreads(num_wann_chk), &
            sawf_template%basis(size(state_cache%source%basis,1),size(state_cache%source%basis,2)), &
            sawf_template%buffer_basis(size(state_cache%source%buffer_basis,1), &
              size(state_cache%source%buffer_basis,2)), &
            sawf_template%orbitals(size(sawf_representative_orbitals,1),num_wann_chk), &
            sawf_template%buffer_orbitals(size(sawf_representative_buffer_orbitals,1),num_wann_chk), &
            sawf_template%band_to_wannier(num_bands_chk,num_wann_chk), &
            sawf_template%d_band(size(sawf_local_d_band,1),size(sawf_local_d_band,2), &
              size(sawf_local_d_band,3)), &
            sawf_template%d_wann(size(sawf_local_d_wann,1),size(sawf_local_d_wann,2), &
              size(sawf_local_d_wann,3)),sawf_template%gauge_unitary(num_wann_chk,num_wann_chk))
          sawf_template%centers=sawf_local_centers;sawf_template%spreads=sawf_local_spreads
          sawf_template%basis=state_cache%source%basis
          sawf_template%buffer_basis=state_cache%source%buffer_basis
          sawf_template%orbitals=sawf_representative_orbitals
          sawf_template%buffer_orbitals=sawf_representative_buffer_orbitals
          sawf_template%band_to_wannier=sawf_local_v_matrix
          sawf_template%d_band=sawf_local_d_band;sawf_template%d_wann=sawf_local_d_wann
          sawf_template%gauge_unitary=(0d0,0d0)
          do iop=1,num_wann_chk;sawf_template%gauge_unitary(iop,iop)=(1d0,0d0);end do
          sawf_template%gauge_residual=0d0
          call write_sawf_template_checkpoint(trim(sawf_seed_bundles(ibundle)%directory)//'/'// &
            trim(sawf_seed_bundles(ibundle)%seedname)//'.sawf-template',sawf_template,local_ok,message)
          if(.not.local_ok)call lcfo_sawf_fatal('SAWF representative template publication failed: '//trim(message))
        end if
        call comm_sync_all(dc%icomm_tot)
        if(dc%id_frag==0)then
          representative=sawf_environment_receipts(dc%i_frag)%representative_fragment
          materialize_operation=sawf_environment_receipts(dc%i_frag)%operation_index
          ibundle=findloc(sawf_seed_bundles%environment,representative,dim=1)
          if(ibundle<=0)call lcfo_sawf_fatal('SAWF representative template bundle lookup failed')
          sawf_inside_atom=.false.
          do ia=1,dc%system_tot%nion
            sawf_inside_atom(ia)=sawf_atom_inside_fragment_buffer(fractional_positions(:,ia),mesh, &
              fragment_origin(:,representative),fragment_shape(:,representative),dc%nxyz_buffer)
          end do
          call select_sawf_local_complete_shells(channels%atom,sawf_expected_channels, &
            sawf_inside_atom,sawf_source_channels,local_ok,message)
          if(.not.local_ok)call lcfo_sawf_fatal('SAWF representative source-shell selection failed')
          call select_sawf_environment_stabilizer(representative,symmetry_fragment_maps,product_left, &
            product_right,product_result,sawf_local_stabilizer,local_ok,message)
          if(.not.local_ok)call lcfo_sawf_fatal('SAWF representative template stabilizer rebuild failed')
          sawf_template_fingerprint%geometry=sawf_supercell_fingerprint
          sawf_template_fingerprint%pseudopotential=sawf_supercell_fingerprint
          sawf_template_fingerprint%grid=sawf_environment_key(representative)
          write(sawf_template_fingerprint%band_window,'(a,i0,a,i0)')'bands=', &
            sawf_environment_receipts(representative)%num_bands,':wann=', &
            sawf_environment_receipts(representative)%num_wann
          write(sawf_template_fingerprint%complete_projection_shell,'(a,i0)') &
            'selected_channels=',size(sawf_source_channels)
          write(sawf_template_fingerprint%symmetry,'(a,i0)')'actual_stabilizer=',size(sawf_local_stabilizer)
          write(sawf_template_fingerprint%buffer,'(3(i0,1x))')dc%nxyz_buffer
          sawf_template_fingerprint%generator='SALMON hierarchical SAWF template schema 2'
          call read_sawf_template_checkpoint(trim(sawf_seed_bundles(ibundle)%directory)//'/'// &
            trim(sawf_seed_bundles(ibundle)%seedname)//'.sawf-template',sawf_template_fingerprint, &
            sawf_template,sawf_template_reuse,local_ok,message)
          if(.not.local_ok.or..not.sawf_template_reuse)call lcfo_sawf_fatal( &
            'SAWF representative template read/fingerprint validation failed: '//trim(message))
          call build_sawf_fragment_buffer_point_map(symmetry_operations(materialize_operation),mesh, &
            fragment_origin(:,representative),fragment_shape(:,representative), &
            fragment_origin(:,dc%i_frag),fragment_shape(:,dc%i_frag),[0,0,0], &
            wannier_symmetry_tolerance,sawf_core_materialize_map,local_ok,message)
          if(local_ok)call build_sawf_fragment_buffer_point_map(symmetry_operations(materialize_operation), &
            mesh,fragment_origin(:,representative),fragment_shape(:,representative), &
            fragment_origin(:,dc%i_frag),fragment_shape(:,dc%i_frag),dc%nxyz_buffer, &
            wannier_symmetry_tolerance,sawf_buffer_materialize_map,local_ok,message)
          if(.not.local_ok)call lcfo_sawf_fatal('SAWF representative-to-target grid map failed: '//trim(message))
          if(size(sawf_source_channels)/=size(sawf_selected_channels))call lcfo_sawf_fatal( &
            'SAWF representative and target ragged channel counts disagree')
          allocate(sawf_materialize_d_wann(size(sawf_source_channels),size(sawf_selected_channels)))
          sawf_materialize_d_wann=d_wann_set(sawf_source_channels,sawf_selected_channels,materialize_operation)
          call materialize_sawf_ragged_local_basis(sawf_template%orbitals,sawf_template%buffer_orbitals, &
            sawf_core_materialize_map,sawf_buffer_materialize_map,sawf_materialize_d_wann, &
            representative,materialize_operation,sawf_generate_independently(dc%i_frag), &
            sawf_materialized_basis,local_ok,message)
          if(.not.local_ok)call lcfo_sawf_fatal('SAWF ragged production materialization failed: '//trim(message))
          write(local_basis_filename,'(a,"/fragment-",i0,".sawf-local-basis")') &
            trim(sawf_seed_bundles(ibundle)%directory),dc%i_frag
          call write_sawf_materialized_basis_checkpoint(trim(local_basis_filename), &
            sawf_supercell_fingerprint,dc%i_frag,sawf_materialized_basis,local_ok,message)
          if(.not.local_ok)then
            write(*,'(1x,a,i0,2a)')'[DC-LCFO-SAWF-LOCAL-BASIS] rank=',dc%id_tot,' ',trim(message)
            call lcfo_sawf_fatal('SAWF local basis publication failed')
          end if
        end if
        call comm_sync_all(dc%icomm_tot)
        allocate(sawf_gauge_parent(dc%n_frag))
        call build_sawf_fragment_gauge_tree(mesh,fragment_origin,fragment_shape,sawf_gauge_parent, &
          local_ok,message)
        if(.not.local_ok)call lcfo_sawf_fatal('SAWF gauge spanning tree construction failed: '//trim(message))
        do ifrag=2,dc%n_frag
          if(dc%id_frag==0.and.dc%i_frag==ifrag)then
            parent_fragment=sawf_gauge_parent(ifrag)
            parent_representative=sawf_environment_receipts(parent_fragment)%representative_fragment
            ibundle=findloc(sawf_seed_bundles%environment,parent_representative,dim=1)
            if(ibundle<=0)call lcfo_sawf_fatal('SAWF parent gauge checkpoint bundle lookup failed')
            write(parent_basis_filename,'(a,"/fragment-",i0,".sawf-local-basis")') &
              trim(sawf_seed_bundles(ibundle)%directory),parent_fragment
            call read_sawf_materialized_basis_checkpoint(trim(parent_basis_filename), &
              sawf_supercell_fingerprint,parent_fragment,sawf_parent_basis,sawf_template_reuse,local_ok,message)
            if(.not.local_ok.or..not.sawf_template_reuse)call lcfo_sawf_fatal( &
              'SAWF parent gauge checkpoint validation failed: '//trim(message))
            call build_sawf_shared_buffer_point_maps(mesh,fragment_origin(:,parent_fragment), &
              fragment_shape(:,parent_fragment),fragment_origin(:,ifrag),fragment_shape(:,ifrag), &
              dc%nxyz_buffer,sawf_parent_shared_map,sawf_child_shared_map,local_ok,message)
            if(.not.local_ok)call lcfo_sawf_fatal('SAWF neighbor shared-buffer map failed: '//trim(message))
            call stitch_sawf_materialized_neighbor_pair(sawf_parent_basis,sawf_materialized_basis, &
              sawf_parent_shared_map,sawf_child_shared_map,hvol,1d-12,wannier_sawf_gauge_tolerance, &
              local_ok,message)
            if(.not.local_ok)call lcfo_sawf_fatal('SAWF neighbor gauge stitching failed: '//trim(message))
            representative=sawf_environment_receipts(ifrag)%representative_fragment
            ibundle=findloc(sawf_seed_bundles%environment,representative,dim=1)
            if(ibundle<=0)call lcfo_sawf_fatal('SAWF child gauge checkpoint bundle lookup failed')
            write(local_basis_filename,'(a,"/fragment-",i0,".sawf-local-basis")') &
              trim(sawf_seed_bundles(ibundle)%directory),ifrag
            call write_sawf_materialized_basis_checkpoint(trim(local_basis_filename), &
              sawf_supercell_fingerprint,ifrag,sawf_materialized_basis,local_ok,message)
            if(.not.local_ok)call lcfo_sawf_fatal('SAWF stitched local basis publication failed: '//trim(message))
          end if
          call comm_sync_all(dc%icomm_tot)
        end do
        sawf_gauge_residual_local=0d0;sawf_gauge_residual_global=0d0
        do ifrag=2,dc%n_frag
          if(dc%id_frag==0.and.dc%i_frag==ifrag)then
            representative=sawf_environment_receipts(ifrag)%representative_fragment
            ibundle=findloc(sawf_seed_bundles%environment,representative,dim=1)
            if(ibundle<=0)call lcfo_sawf_fatal('SAWF gauge-closure child bundle lookup failed')
            write(local_basis_filename,'(a,"/fragment-",i0,".sawf-local-basis")') &
              trim(sawf_seed_bundles(ibundle)%directory),ifrag
            call read_sawf_materialized_basis_checkpoint(trim(local_basis_filename), &
              sawf_supercell_fingerprint,ifrag,sawf_materialized_basis,sawf_template_reuse,local_ok,message)
            if(.not.local_ok.or..not.sawf_template_reuse)call lcfo_sawf_fatal( &
              'SAWF gauge-closure child checkpoint validation failed: '//trim(message))
            do parent_fragment=1,ifrag-1
              if(.not.sawf_fragments_share_face(mesh,fragment_origin(:,parent_fragment), &
                  fragment_shape(:,parent_fragment),fragment_origin(:,ifrag),fragment_shape(:,ifrag)))cycle
              parent_representative=sawf_environment_receipts(parent_fragment)%representative_fragment
              ibundle=findloc(sawf_seed_bundles%environment,parent_representative,dim=1)
              if(ibundle<=0)call lcfo_sawf_fatal('SAWF gauge-closure parent bundle lookup failed')
              write(parent_basis_filename,'(a,"/fragment-",i0,".sawf-local-basis")') &
                trim(sawf_seed_bundles(ibundle)%directory),parent_fragment
              call read_sawf_materialized_basis_checkpoint(trim(parent_basis_filename), &
                sawf_supercell_fingerprint,parent_fragment,sawf_parent_basis,sawf_template_reuse,local_ok,message)
              if(.not.local_ok.or..not.sawf_template_reuse)call lcfo_sawf_fatal( &
                'SAWF gauge-closure parent checkpoint validation failed: '//trim(message))
              call build_sawf_shared_buffer_point_maps(mesh,fragment_origin(:,parent_fragment), &
                fragment_shape(:,parent_fragment),fragment_origin(:,ifrag),fragment_shape(:,ifrag), &
                dc%nxyz_buffer,sawf_parent_shared_map,sawf_child_shared_map,local_ok,message)
              if(.not.local_ok)call lcfo_sawf_fatal('SAWF gauge-closure shared map failed: '//trim(message))
              call validate_sawf_materialized_neighbor_closure(sawf_parent_basis,sawf_materialized_basis, &
                sawf_parent_shared_map,sawf_child_shared_map,hvol,1d-12,wannier_sawf_gauge_tolerance, &
                sawf_gauge_closure_residual,sawf_gauge_alignment_residual,local_ok,message)
              if(.not.local_ok)then
                write(*,'(1x,a,i0,a,i0,2(a,es13.5))')'[DC-LCFO-SAWF-GAUGE-CLOSURE] left=', &
                  parent_fragment,' right=',ifrag,' closure=',sawf_gauge_closure_residual, &
                  ' alignment=',sawf_gauge_alignment_residual
                call lcfo_sawf_fatal('SAWF all-neighbor gauge closure failed: '//trim(message))
              end if
              sawf_gauge_residual_local(1)=max(sawf_gauge_residual_local(1),sawf_gauge_closure_residual)
              sawf_gauge_residual_local(2)=max(sawf_gauge_residual_local(2),sawf_gauge_alignment_residual)
            end do
          end if
        end do
        call comm_get_max(sawf_gauge_residual_local,sawf_gauge_residual_global,2,dc%icomm_tot)
        if(dc%id_tot==0)write(*,'(1x,a,2(a,es13.5))')'[DC-LCFO-SAWF-GAUGE-SUMMARY]', &
          ' closure_max=',sawf_gauge_residual_global(1), &
          ' alignment_max=',sawf_gauge_residual_global(2)
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
      if(dc%id_tot == 0) write(*,'(1x,a,i0,a,i0,6(a,es13.5))') &
        '[DC-LCFO-SAWF-DMN] published operations=',size(symmetry_operations), &
        ' bands=',nband_wann,' unitarity_max=',writer%max_unitarity, &
        ' hamiltonian_max=',writer%max_hamiltonian, &
        ' hamiltonian_tolerance=',writer%hamiltonian_tolerance,' amn_max=',writer%max_amn, &
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
        bundles,ok,message,neighbor_gvec,local_states_out,local_energy_out,local_amn_out,local_coeff_out)
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
      real(8),allocatable,intent(out),optional::local_coeff_out(:,:)
      real(8),parameter::hartree_to_ev=27.211386245988d0
      real(8),allocatable::h_basis(:,:),energy_hartree(:),energy_ev(:),local_coefficients(:,:)
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
        1d-10,states,energy_hartree,ok,message,local_coefficients)
      if(.not.ok)return
      if(present(local_states_out))then
        allocate(local_states_out(size(states,1),size(states,2)));local_states_out=states
      end if
      if(present(local_coeff_out))then
        allocate(local_coeff_out(size(local_coefficients,1),size(local_coefficients,2)))
        local_coeff_out=local_coefficients
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

    subroutine lcfo_sawf_rank_local_fatal(message)
      implicit none
      character(*), intent(in) :: message

      write(*,'(1x,a,i0,a,i0,2a)') '[FATAL] SAWF rank-local rank=',dc%id_tot, &
        ' fragment=',dc%i_frag,' reason=',trim(message)
      stop 1
    end subroutine lcfo_sawf_rank_local_fatal

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
      integer :: iunit, ibasis_read, ispin_read, iw, io_status
      integer :: magic_file, version_file, nspin_file, nstate_frag_file
      integer :: nxyz_domain_file(3), nxyz_buffer_file(3)
      integer :: axis, side, face, npts, face_pt
      integer :: ix_face, iy_face, iz_face, ibx, iby, ibz, dist
      integer, allocatable :: n_basis_file(:)
      real(8), allocatable :: phi_tmp(:,:,:), phi_wann(:,:,:,:)
      real(8), allocatable :: face_u(:,:), face_dn(:,:)
      real(8) :: grad_axis, area_weight, alpha
      character(256) :: filename
      character(512) :: io_message
      real(8), parameter :: surface_penalty_factor = 10.0d0

      if(dc%id_frag /= 0) return
      if(nspin /= 1) call lcfo_sawf_rank_local_fatal( &
        'global Wannier trace export supports nspin=1 only')
      if(any(nxyz_buffer_seed(1:3)<size(stencil%coef_nab,1))) &
        call lcfo_sawf_rank_local_fatal( &
          'global Wannier trace export buffer is smaller than finite-difference stencil radius')

      allocate(phi_tmp(nxyz_box(1),nxyz_box(2),nxyz_box(3)))
      allocate(phi_wann(nxyz_box(1),nxyz_box(2),nxyz_box(3),nkeep))
      allocate(n_basis_file(nspin))
      phi_wann = 0d0

      if(sawf_explicit_basis_active) then
        if(.not.allocated(sawf_explicit_buffer) .or. &
           any(shape(sawf_explicit_buffer(:,:,:,1,1))/=nxyz_box) .or. &
           size(sawf_explicit_buffer,5)<nbasis) &
          call lcfo_sawf_rank_local_fatal( &
            'global Wannier trace export closed SAWF buffer basis count mismatch')
        do ibasis_read=1,nbasis
          do iw=1,nkeep
            if(abs(wcoef(ibasis_read,iw)) <= 0d0) cycle
            phi_wann(:,:,:,iw) = phi_wann(:,:,:,iw) + &
              wcoef(ibasis_read,iw) * sawf_explicit_buffer(:,:,:,1,ibasis_read)
          end do
        end do
      else
        filename = trim(base_directory)//binfile_bfb
        iunit = get_filehandle()
        open(iunit,file=filename,form='unformatted',access='stream',status='old', &
          iostat=io_status,iomsg=io_message)
        if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
          'global Wannier trace export cannot open buffered basis: '//trim(io_message))
        read(iunit,iostat=io_status,iomsg=io_message) magic_file, version_file
        if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
          'global Wannier trace export cannot read buffered basis header: '//trim(io_message))
        if(magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) &
          call lcfo_sawf_rank_local_fatal('global Wannier trace export invalid buffered basis header')
        read(iunit,iostat=io_status,iomsg=io_message) nxyz_domain_file(1:3), &
          nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
        if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
          'global Wannier trace export cannot read buffered basis metadata: '//trim(io_message))
        if(any(nxyz_domain_file(1:3) /= nxyz_domain(1:3)) .or. &
           any(nxyz_buffer_file(1:3) /= nxyz_buffer_seed(1:3)) .or. &
           nspin_file /= nspin .or. nstate_frag_file /= dc%nstate_frag) &
          call lcfo_sawf_rank_local_fatal('global Wannier trace export buffered basis metadata mismatch')
        read(iunit,iostat=io_status,iomsg=io_message) n_basis_file(1:nspin)
        if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
          'global Wannier trace export cannot read buffered basis count: '//trim(io_message))
        if(n_basis_file(1) < nbasis) call lcfo_sawf_rank_local_fatal( &
          'global Wannier trace export basis count mismatch')
        do ispin_read=1,nspin_file
          do ibasis_read=1,nstate_frag_file
            read(iunit,iostat=io_status,iomsg=io_message) &
              phi_tmp(1:nxyz_box(1),1:nxyz_box(2),1:nxyz_box(3))
            if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
              'global Wannier trace export cannot read buffered basis payload: '//trim(io_message))
            if(ispin_read /= 1 .or. ibasis_read > nbasis) cycle
            do iw=1,nkeep
              if(abs(wcoef(ibasis_read,iw)) <= 0d0) cycle
              phi_wann(:,:,:,iw) = phi_wann(:,:,:,iw) + wcoef(ibasis_read,iw) * phi_tmp(:,:,:)
            end do
          end do
        end do
        close(iunit)
      end if

      filename = trim(base_directory)//binfile_bpwt
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='replace', &
        iostat=io_status,iomsg=io_message)
      if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
        'global Wannier trace export cannot open output: '//trim(io_message))
      write(iunit,iostat=io_status,iomsg=io_message) &
        buffer_periodic_wannier_trace_magic, buffer_periodic_wannier_trace_version
      if(io_status==0) write(iunit,iostat=io_status,iomsg=io_message) &
        dc%i_frag, nxyz_domain(1:3), nxyz_buffer_seed(1:3), nxyz_box(1:3)
      if(io_status==0) write(iunit,iostat=io_status,iomsg=io_message) system%hgs(1:3), hvol, nkeep
      if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
        'global Wannier trace export cannot write output header: '//trim(io_message))
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
        write(iunit,iostat=io_status,iomsg=io_message) axis, side, npts, area_weight, alpha
        if(io_status==0) write(iunit,iostat=io_status,iomsg=io_message) face_u(1:npts,1:nkeep)
        if(io_status==0) write(iunit,iostat=io_status,iomsg=io_message) face_dn(1:npts,1:nkeep)
        if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
          'global Wannier trace export cannot write face payload: '//trim(io_message))
        deallocate(face_u, face_dn)
      end do
      close(iunit,iostat=io_status,iomsg=io_message)
      if(io_status/=0) call lcfo_sawf_rank_local_fatal( &
        'global Wannier trace export cannot close output: '//trim(io_message))
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

    subroutine build_bond_center_projection_map(nproj,center_bohr,atom_ids,bond_images)
      implicit none
      integer, intent(in) :: nproj
      real(8), allocatable, intent(out) :: center_bohr(:,:)
      integer,allocatable,optional,intent(out)::atom_ids(:,:),bond_images(:,:)
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
      if(present(atom_ids))allocate(atom_ids(2,nproj))
      if(present(bond_images))allocate(bond_images(3,nproj))
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
            if(present(atom_ids))atom_ids(:,ip)=[ia,ja]
            do axis=1,3
              center(axis) = dc%system_tot%rion(axis,ia) + 0.5d0 * delta(axis)
              if(length_axis(axis) > 0d0) center(axis) = center(axis) - floor(center(axis) / length_axis(axis)) &
                * length_axis(axis)
              center_bohr(axis,ip) = center(axis)
              if(present(bond_images))bond_images(axis,ip)=nint((delta(axis)-&
                (dc%system_tot%rion(axis,ja)-dc%system_tot%rion(axis,ia)))/length_axis(axis))
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
