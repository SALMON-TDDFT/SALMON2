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
! Types for DG-Fragment RT-TDDFT
!=======================================================================

module rt_dg_fragment_types
  use structures, only: s_rgrid
  implicit none

  private
  public :: halo_info, matrix_block_info, complex_matrix_block_info, vector_block_info, complex_vector_block_info, &
            momentum_block_info, flux_face_trace_info, density_recv_map_info, density_grid_point_info, &
            real_buffer_info, dg_wannier_neighbor_block_info, s_dg_fragment_rt, invalidate_coef_exchange_cache

  ! Halo communication structure (for phi_frag exchange between fragments)
  type :: halo_info
    integer :: id_src, id_dst                ! MPI ranks for communication
    integer :: ifrag_src                     ! Source fragment index
    integer :: ifrag_dst                     ! Destination (local) fragment index
    integer :: dvec(3)                       ! Direction vector to neighbor (-1,0,+1)
    integer :: length(3)                     ! Size of halo region in each direction
    integer :: dsp_send(3), dsp_recv(3)      ! Displacement for send/recv buffers
    integer :: send_lo(3), send_hi(3)        ! Local phi_frag send block bounds
    integer :: recv_lo(3), recv_hi(3)        ! Local phi_frag recv block bounds
    real(8), allocatable :: buf_send(:,:,:,:,:)  ! (lx,ly,lz,nstate,ifrag_count)
    real(8), allocatable :: buf_recv(:,:,:,:,:)  ! (lx,ly,lz,nstate,ifrag_count)
    complex(8), allocatable :: buf_send_c(:,:,:,:,:)  ! complex halo buffer for SOI path
    complex(8), allocatable :: buf_recv_c(:,:,:,:,:)  ! complex halo buffer for SOI path
  end type halo_info

  type :: matrix_block_info
    integer :: ifrag_row = 0
    integer :: ifrag_col = 0
    integer :: nrow_max = 0
    integer :: ncol_max = 0
    real(8), allocatable :: val(:,:,:) ! (nrow_max, ncol_max, nspin)
  end type matrix_block_info

  type :: complex_matrix_block_info
    integer :: ifrag_row = 0
    integer :: ifrag_col = 0
    integer :: nrow_max = 0
    integer :: ncol_max = 0
    complex(8), allocatable :: val(:,:,:) ! (nrow_max, ncol_max, nspin)
  end type complex_matrix_block_info

  type :: vector_block_info
    integer :: ifrag_row = 0
    integer :: ifrag_col = 0
    integer :: nrow_max = 0
    integer :: ncol_max = 0
    real(8), allocatable :: val(:,:,:,:) ! (ndir, nrow_max, ncol_max, nspin)
  end type vector_block_info

  type :: complex_vector_block_info
    integer :: ifrag_row = 0
    integer :: ifrag_col = 0
    integer :: nrow_max = 0
    integer :: ncol_max = 0
    complex(8), allocatable :: val(:,:,:,:) ! (ndir, nrow_max, ncol_max, nspin)
  end type complex_vector_block_info

  type :: momentum_block_info
    integer :: ifrag_row = 0
    integer :: ifrag_col = 0
    integer :: nrow_max = 0
    integer :: ncol_max = 0
    real(8), allocatable :: val(:,:,:,:) ! (3, nrow_max, ncol_max, nspin)
  end type momentum_block_info

  type :: flux_face_trace_info
    integer :: ifrag = 0
    integer :: jfrag = 0
    integer :: axis = 0
    integer :: side = 0
    integer :: npts = 0
    integer :: ncol_max = 0
    integer :: nspin = 0
    logical :: initialized = .false.
    real(8), allocatable :: u(:,:,:)   ! (face point, ket basis, spin)
    real(8), allocatable :: dn(:,:,:)  ! (face point, ket basis, spin)
  end type flux_face_trace_info

  type :: density_recv_map_info
    integer :: npts = 0
    integer, allocatable :: ixg(:)
    integer, allocatable :: iyg(:)
    integer, allocatable :: izg(:)
    integer, allocatable :: bx(:)
    integer, allocatable :: by(:)
    integer, allocatable :: bz(:)
  end type density_recv_map_info

  type :: density_grid_point_info
    integer :: ix = 1
    integer :: iy = 1
    integer :: iz = 1
    integer :: ixg = 1
    integer :: iyg = 1
    integer :: izg = 1
    integer :: owner_rank = -1
    integer :: send_slot = 0
    integer :: subgroup_send_slot = 0
    logical :: is_primary = .false.
  end type density_grid_point_info

  type :: real_buffer_info
    real(8), allocatable :: val(:)
  end type real_buffer_info

  type :: dg_wannier_neighbor_block_info
    integer :: ifrag_row = 0
    integer :: ifrag_col = 0
    integer :: nrow_max = 0
    integer :: ncol_max = 0
    integer :: peer_rank = -1
    complex(8), allocatable :: h_flux(:,:,:)     ! (nrow_max,ncol_max,nspin)
    complex(8), allocatable :: xi_flux(:,:,:,:)  ! (3,nrow_max,ncol_max,nspin)
  end type dg_wannier_neighbor_block_info

  ! Fragment basis data structure
  !
  ! Spin convention in DG-RT:
  ! - phi_frag / phi_frag_c store the fragment real-space basis itself.
  !   In the current non-SOI DG path this basis is not duplicated per spin.
  ! - n_basis, index_basis, n_mat, coef, H_mat, S_mat, momentum_mat, etc.
  !   carry an nspin axis because the occupied subspace bookkeeping and the
  !   operator representation on that basis are resolved by spin channel.
  ! - Therefore "spin-dependent basis" in this module means spin-resolved
  !   basis indexing / coefficients / projected operators, not necessarily
  !   distinct real-space basis functions for each spin channel.
  type :: s_dg_fragment_rt
    integer :: n_frag                          ! number of fragments
    integer :: nspin                           ! number of spin
    integer :: nstate_frag                     ! number of states per fragment
    integer :: nstate_tot                      ! total number of states
    integer, allocatable :: nocc_spin(:)       ! (nspin), occupied-state count per spin channel
    integer :: time_integrator                 ! 1: SSPRK3, 2: AETRS, 3: RK4

    ! Spin-resolved coefficients on the fragment basis (MUST be complex).
    ! The first index is the global fragment-basis index after index_basis mapping,
    ! not a separate real-space basis copy for each spin channel.
    complex(8), allocatable :: coef(:,:,:)        ! (nstate_frag, nstate_tot, nspin)
    complex(8), allocatable :: coef_new(:,:,:)    ! for time propagation
    complex(8), allocatable :: coef_work(:,:,:)   ! work array
    logical :: coef_state_block_mode = .false.     ! store only local state columns, with fragment-local rows
    integer :: coef_state_start = 1                ! first global state column stored on this rank
    integer :: coef_state_end = 0                  ! last global state column stored on this rank
    integer :: coef_nstate_local = 0               ! local state-column count
    ! Orbital parallelism stores only the coefficient rows owned by this rank.
    integer, allocatable :: local_coef_count(:)        ! (nspin), local coefficient-row count
    integer, allocatable :: local_coef_global_ids(:,:) ! (local_coef_max,nspin), local row -> global basis id
    integer, allocatable :: coef_global_to_local(:,:)  ! (n_mat_max,nspin), global basis id -> local row, 0 if absent
    integer, allocatable :: coef_owner(:,:)       ! (n_mat_max, nspin), owning rank of each fragment-basis row
    integer, allocatable :: coef_row_fragment(:,:) ! (n_mat_max,nspin), global basis id -> fragment id
    integer, allocatable :: coef_exchange_peer_ranks(:,:) ! cached neighbor coefficient exchange peers
    integer, allocatable :: coef_exchange_peer_count(:)   ! (nspin), number of cached exchange peers
    logical, allocatable :: coef_allowed_request_frag(:)  ! fragments reachable by local row-owner stencil
    logical :: coef_exchange_peer_cache_valid = .false.
    integer :: owned_coef_start = 0               ! first owned fragment-basis row (contiguous hint)
    integer :: owned_coef_end = -1                ! last owned fragment-basis row (contiguous hint)

    ! Spin-channel-resolved operator matrices on the fragment basis.
    ! These arrays are spin-resolved projected representations, even when the
    ! underlying real-space fragment basis functions are spin independent.
    real(8), allocatable :: H_mat(:,:,:)       ! (nmat, nmat, nspin)
    complex(8), allocatable :: H_mat_c(:,:,:)      ! complex Hamiltonian (SOI/mixed propagation path)
    type(matrix_block_info), allocatable :: H_mat_blocks(:)
    type(matrix_block_info), allocatable :: H_mat_core_blocks(:)
    type(matrix_block_info), allocatable :: H_mat_kinetic_blocks(:)
    type(complex_matrix_block_info), allocatable :: H_nl_blocks(:)
    logical :: H_blocks_include_nonlocal = .false. ! true when H_mat_blocks are read from DC exported full H
    integer, allocatable :: H_block_map(:,:)
    integer :: n_H_blocks = 0
    integer, allocatable :: H_nl_block_map(:,:)
    integer :: n_H_nl_blocks = 0
    integer, allocatable :: H_local_block_ids(:) ! row-owner-local H block ids for RT apply
    integer, allocatable :: H_nl_local_block_ids(:) ! row-owner-local nonlocal H block ids for RT apply
    integer, allocatable :: H_local_rows(:)      ! fragment rows owned by this rank
    logical :: H_local_block_ids_valid = .false.
    type(flux_face_trace_info), allocatable :: flux_face_trace_cache(:)
    logical :: flux_face_trace_mix_enabled = .false.
    logical, allocatable :: runtime_neighbor_pair_cache(:,:)    ! static fragment-pair runtime adjacency
    integer, allocatable :: runtime_neighbor_frag_count(:)       ! runtime-neighbor list length for each fragment
    integer, allocatable :: runtime_neighbor_frag_ids(:,:)       ! compact runtime-neighbor fragment ids
    logical, allocatable :: momentum_neighbor_pair_cache(:,:)   ! static fragment-pair momentum adjacency
    integer, allocatable :: momentum_neighbor_frag_count(:)      ! face-neighbor list length for each fragment
    integer, allocatable :: momentum_neighbor_frag_ids(:,:)      ! compact face-neighbor fragment ids, including self
    real(8), allocatable :: S_mat(:,:,:)       ! raw fragment overlap matrix
    real(8), allocatable :: S_mat_prop(:,:,:)  ! overlap matrix used in propagation/unitarity
    complex(8), allocatable :: S_mat_c(:,:,:)      ! complex raw fragment overlap matrix (SOI path)
    complex(8), allocatable :: S_mat_prop_c(:,:,:) ! complex overlap used in propagation/unitarity
    type(matrix_block_info), allocatable :: S_mat_blocks(:)
    type(matrix_block_info), allocatable :: S_mat_prop_blocks(:)
    ! Complex overlap blocks are used by SOI and mixed/PW propagation paths.
    type(complex_matrix_block_info), allocatable :: S_mat_blocks_c(:)
    type(complex_matrix_block_info), allocatable :: S_mat_prop_blocks_c(:)
    integer, allocatable :: S_block_map(:,:)
    integer :: n_S_blocks = 0
    logical :: has_global_overlap_copy = .true.
    logical :: overlap_prop_root_authoritative = .false.
    integer, allocatable :: n_basis(:,:)       ! (n_frag, nspin), spin-resolved active basis count per fragment
    integer, allocatable :: index_basis(:,:,:) ! (nstate_frag, n_frag, nspin), spin-resolved local->global basis map
    integer, allocatable :: n_mat(:)           ! (nspin), spin-resolved projected-matrix dimension
    integer :: n_mat_max                       ! max projected-matrix dimension over spin channels
    logical :: identity_seed_coefficients = .false.      ! DC seed used fragment-local identity coefficients
    logical :: fragment_basis_contracted = .false.       ! basis is RT-ready after DC cleanup or legacy contraction
    logical :: dc_lcfo_seed_basis_cleaned = .false.      ! DC export says basis is core-S-cleaned and RT-ready

    ! Time-dependent external field coupling (velocity gauge: H = H_0 - i*A·∇ + A^2/2)
    real(8), allocatable :: momentum_mat(:,:,:,:) ! momentum matrix elements p_ij = <phi_i|p|phi_j> (x,y,z)
    complex(8), allocatable :: momentum_mat_c(:,:,:,:) ! complex momentum matrix for SOI/mixed propagation
    type(vector_block_info), allocatable :: momentum_blocks(:)
    type(complex_vector_block_info), allocatable :: momentum_blocks_c(:)
    integer, allocatable :: momentum_block_map(:,:)
    integer :: n_momentum_blocks = 0
    logical :: momentum_blocks_include_dg_flux = .false. ! true after covariant DG surface-Flux A operator is added
    integer, allocatable :: momentum_dense_row_gid_cache(:)    ! reusable scratch for momentum dense materialization
    integer, allocatable :: momentum_dense_col_gid_cache(:)    ! reusable scratch for momentum dense materialization
    integer, allocatable :: momentum_dense_valid_row_ids(:)    ! reusable scratch for momentum dense materialization
    integer, allocatable :: momentum_dense_valid_col_ids(:)    ! reusable scratch for momentum dense materialization
    real(8), allocatable :: gradient_basis_cache(:,:,:,:,:,:)  ! (nx,ny,nz,3,nstate_frag,ifrag_local)
    logical :: gradient_basis_cache_valid = .false.
    real(8), allocatable :: dipole_mat(:,:,:,:)   ! dipole matrix elements for observables (x,y,z)
    type(vector_block_info), allocatable :: dipole_blocks(:)
    integer, allocatable :: dipole_block_map(:,:)
    integer :: n_dipole_blocks = 0
    complex(8), allocatable :: full_h_seed_evec(:,:,:) ! (DG basis, Full-H eigenstate, spin)
    real(8), allocatable :: full_h_seed_eval(:,:)       ! (Full-H eigenstate, spin)
    complex(8), allocatable :: full_h_seed_xi(:,:,:,:)  ! (x/y/z, Full-H eigenstate, Full-H eigenstate, spin)
    integer :: full_h_seed_nstate = 0
    logical :: has_full_h_seed_eigen = .false.
    logical :: has_full_h_seed_xi = .false.
    logical :: has_local_wannier_basis = .false.
    integer, allocatable :: local_wannier_nbasis(:) ! (local fragment index)
    integer, allocatable :: local_wannier_nproj(:)  ! (local fragment index)
    integer, allocatable :: local_wannier_nkeep(:)  ! (local fragment index)
    real(8), allocatable :: local_wannier_coef(:,:,:,:) ! (DG basis, W basis, nspin, local fragment index)
    real(8), allocatable :: local_wannier_r(:,:,:,:,:) ! (3, nkeep_max, nkeep_max, nspin, local fragment index)
    real(8), allocatable :: local_wannier_center(:,:,:,:) ! (3, W basis, nspin, local fragment index)
    integer, allocatable :: local_wannier_owner_fragment(:,:,:) ! (W basis, nspin, local fragment index)
    logical, allocatable :: local_wannier_owned(:,:,:) ! true when center belongs to this rank's fragment
    logical :: has_buffer_periodic_wannier_basis = .false.
    logical :: buffer_wannier_flux_seed_applied = .false.
    integer, allocatable :: buffer_wannier_nkeep(:) ! (local fragment index)
    real(8), allocatable :: buffer_wannier_coef(:,:,:,:) ! (DG basis, W basis, nspin, local fragment index)
    real(8), allocatable :: buffer_wannier_spread(:,:) ! (W basis, local fragment index)
    real(8), allocatable :: buffer_wannier_tail(:,:) ! (W basis, local fragment index)
    real(8), allocatable :: buffer_wannier_h_flux(:,:,:) ! (W basis, W basis, local fragment index)
    real(8), allocatable :: buffer_wannier_v(:,:,:,:) ! xi: (3, W basis, W basis, local fragment index)
    real(8), allocatable :: buffer_wannier_center(:,:,:) ! (3, W basis, local fragment index)
    logical, allocatable :: buffer_wannier_position_has_relative_center(:) ! local fragment index
    real(8), allocatable :: buffer_wannier_frag_center(:,:) ! R_I: (3, local fragment index)
    type(vector_block_info), allocatable :: buffer_wannier_xi_flux_blocks(:) ! neighbor xi_flux blocks
    integer, allocatable :: buffer_wannier_xi_flux_local_block_ids(:)
    integer :: n_buffer_wannier_xi_flux_blocks = 0
    logical :: buffer_wannier_xi_flux_available = .false.
    logical :: buffer_wannier_xi_flux_warned = .false.
    logical :: has_global_wannier_basis = .false.
    integer :: global_wannier_num_bands = 0
    integer :: global_wannier_num_wann = 0
    integer :: global_wannier_n_frag = 0
    integer, allocatable :: global_wannier_owner_frag(:) ! (num_wann)
    real(8), allocatable :: global_wannier_center(:,:) ! (3,num_wann), bohr
    real(8), allocatable :: global_wannier_spread(:) ! (num_wann), Angstrom^2 from Wannier90
    complex(8), allocatable :: global_wannier_transform(:,:) ! (num_bands,num_wann)
    complex(8), allocatable :: global_wannier_position(:,:,:) ! (3,num_wann,num_wann), bohr
    logical :: has_global_wannier_position = .false.
    complex(8), allocatable :: global_wannier_coef(:,:,:,:) ! (DG basis,W basis,nspin,local fragment)
    complex(8), allocatable :: global_wannier_flux_evec(:,:) ! (W basis, Flux eigenstate)
    real(8), allocatable :: global_wannier_flux_eval(:,:) ! (Flux eigenstate,nspin)
    logical :: has_global_wannier_flux_eigen = .false.
    logical :: has_global_wannier_local_basis = .false.
    integer, allocatable :: global_wannier_local_nkeep(:) ! (local fragment)
    integer, allocatable :: global_wannier_local_ids(:,:) ! (local W, local fragment) -> global W id
    integer, allocatable :: global_wannier_local_owner_frag(:,:) ! (local W, local fragment)
    real(8), allocatable :: global_wannier_local_center(:,:,:) ! (3, local W, local fragment)
    complex(8), allocatable :: global_wannier_local_coef(:,:,:,:) ! (DG basis, local W, nspin, local fragment)
    complex(8), allocatable :: global_wannier_local_position(:,:,:,:) ! (3, local W, local W, local fragment)
    logical :: has_mixed_wannier_bpw_position = .false.
    integer :: mixed_wannier_bpw_nw = 0
    integer :: mixed_wannier_bpw_np = 0
    integer :: mixed_wannier_bpw_nmix = 0
    real(8) :: mixed_wannier_bpw_sperp_min = 0.0d0
    real(8) :: mixed_wannier_bpw_sperp_max = 0.0d0
    real(8), allocatable :: mixed_wannier_bpw_eval(:,:) ! (W+P state,nspin)
    complex(8), allocatable :: mixed_wannier_bpw_z(:,:,:,:) ! (3,W+P,W+P,nspin)
    complex(8), allocatable :: mixed_wannier_bpw_pcoef(:,:,:) ! (P state,propagated state,nspin)
    real(8), allocatable :: mixed_z_prod_zww_diag_by_w(:)
    real(8), allocatable :: mixed_z_prod_ww_diag_weight_by_w(:)
    real(8), allocatable :: mixed_z_prod_ww_diag_contrib_by_w(:)
    real(8), allocatable :: mixed_z_prod_ww_diag_rho_by_w(:)
    real(8), allocatable :: mixed_z_prod_ww_diag_occ_by_w(:)
    complex(8), allocatable :: mixed_wannier_bpw_p_transform(:,:,:) ! (raw PW,P state,nspin)
    complex(8), allocatable :: mixed_wannier_bpw_p_metric(:,:,:) ! (P state,P state,nspin)
    integer :: mixed_wannier_bpw_praw_dim = 0
    logical :: has_mixed_wannier_bpw_p_basis = .false.
    real(8) :: mixed_wannier_bpw_sraw_herm_diff = 0.0d0
    real(8) :: mixed_wannier_bpw_sperp_herm_diff = 0.0d0
    real(8) :: mixed_wannier_bpw_qmat_metric_herm_diff = 0.0d0
    real(8) :: mixed_wannier_bpw_qmat_metric_min_eval = 0.0d0
    real(8) :: mixed_wannier_bpw_qmat_metric_max_eval = 0.0d0
    real(8) :: mixed_wannier_bpw_qmat_metric_cond = 0.0d0
    real(8) :: mixed_wannier_bpw_qmat_metric_diff_from_i = 0.0d0
    real(8) :: mixed_wannier_bpw_qleft_metric_diff_from_i = 0.0d0
    real(8) :: mixed_wannier_bpw_final_metric_herm_diff = 0.0d0
    real(8) :: mixed_wannier_bpw_final_metric_min_eval = 0.0d0
    real(8) :: mixed_wannier_bpw_final_metric_max_eval = 0.0d0
    real(8) :: mixed_wannier_bpw_final_metric_cond = 0.0d0
    real(8) :: mixed_wannier_bpw_final_metric_diff_from_i = 0.0d0
    real(8) :: mixed_wannier_bpw_transform_metric_herm_diff = 0.0d0
    real(8) :: mixed_wannier_bpw_transform_metric_min_eval = 0.0d0
    real(8) :: mixed_wannier_bpw_transform_metric_max_eval = 0.0d0
    real(8) :: mixed_wannier_bpw_transform_metric_cond = 0.0d0
    real(8) :: mixed_wannier_bpw_transform_metric_diff_from_i = 0.0d0
    real(8) :: mixed_wannier_bpw_transform_metric_diff_saved = 0.0d0
    real(8) :: mixed_wannier_bpw_h_input_herm_diff = 0.0d0
    real(8) :: mixed_wannier_bpw_h_evec_unitarity_diff = 0.0d0
    real(8) :: mixed_wannier_bpw_h_input_evec_diff = 0.0d0
    logical :: mixed_wannier_bpw_final_uses_h_evec = .false.
    real(8) :: mixed_wannier_bpw_qmat_col_norm_min = 0.0d0
    real(8) :: mixed_wannier_bpw_qmat_col_norm_max = 0.0d0
    real(8) :: mixed_wannier_bpw_qmat_row_norm_min = 0.0d0
    real(8) :: mixed_wannier_bpw_qmat_row_norm_max = 0.0d0
    ! Fragment-local WPW reduced-neighbor propagation candidate.
    ! This is intentionally separate from the global Wannier+BPW-perp mixed-Z path.
    logical :: wpw_reduced_ready = .false.
    integer :: wpw_reduced_shell = 0
    integer :: wpw_reduced_max_dim = 0
    integer, allocatable :: wpw_reduced_dim(:)       ! (local fragment)
    integer, allocatable :: wpw_reduced_nself(:)     ! (local fragment)
    integer, allocatable :: wpw_reduced_nkeep(:)     ! kept S-orthogonal neighbor modes
    integer, allocatable :: wpw_reduced_ndrop(:)     ! dropped near-null neighbor modes
    complex(8), allocatable :: wpw_reduced_H(:,:,:,:) ! (max_dim,max_dim,nspin,local fragment)
    complex(8), allocatable :: wpw_reduced_S(:,:,:,:) ! (max_dim,max_dim,nspin,local fragment)
    complex(8), allocatable :: wpw_reduced_transform(:,:,:) ! raw extended basis <- reduced basis, (max_dim,max_dim,local fragment)
    complex(8), allocatable :: wpw_reduced_Hraw_build(:,:,:) ! raw H used to build reduced transform, (max_dim,max_dim,local fragment)
    complex(8), allocatable :: wpw_reduced_Sraw_build(:,:,:) ! raw S used to build reduced transform, (max_dim,max_dim,local fragment)
    integer, allocatable :: wpw_reduced_nraw(:)       ! raw extended basis size before S-orthogonal reduction
    real(8), allocatable :: wpw_reduced_eval(:,:,:) ! (max_dim,nspin,local fragment)
    complex(8), allocatable :: wpw_reduced_evec(:,:,:,:) ! (max_dim,max_dim,nspin,local fragment)
    complex(8), allocatable :: coef_wpw_self(:,:,:,:) ! (self basis,state,nspin,local fragment)
    complex(8), allocatable :: coef_wpw_neighbor_reduced(:,:,:,:) ! (reduced neighbor,state,nspin,local fragment)
    complex(8), allocatable :: wpw_reproject_prev_coef(:,:,:,:) ! previous reprojected reduced coef for diagnostics
    logical :: wpw_reproject_prev_valid = .false.
    logical :: wpw_reduced_coef_initialized = .false.
    logical :: has_formal_dg_wannier_basis = .false.
    integer, allocatable :: dg_wannier_nkeep(:) ! (local fragment)
    integer, allocatable :: dg_wannier_global_ids(:,:) ! (local W, local fragment)
    integer, allocatable :: dg_wannier_owner_frag(:,:) ! (local W, local fragment)
    real(8), allocatable :: dg_wannier_ref_center(:,:) ! R_I: (3, local fragment)
    complex(8), allocatable :: dg_wannier_basis_coef(:,:,:,:) ! (DG basis, local W, nspin, local fragment)
    complex(8), allocatable :: dg_wannier_h0_local(:,:,:,:) ! (local W, local W, nspin, local fragment)
    complex(8), allocatable :: dg_wannier_xi_local(:,:,:,:,:) ! (3, local W, local W, nspin, local fragment)
    complex(8), allocatable :: dg_wannier_coef(:,:,:,:) ! C_aI(t): (local W, state, nspin, local fragment)
    type(dg_wannier_neighbor_block_info), allocatable :: dg_wannier_neighbor_blocks(:)
    integer, allocatable :: dg_wannier_local_neighbor_block_ids(:)
    integer :: n_dg_wannier_neighbor_blocks = 0
    logical :: dg_wannier_blocks_gauge_consistent = .false.
    logical :: dg_wannier_uses_full_global_position = .false.
    real(8) :: Ac_nl_cache(3)                      ! cached vector potential for H_nl_blocks
    real(8) :: Ac_nl_cache_tol                     ! tolerance for cache reuse
    logical :: has_nl_cache                        ! flag: cached H_nl available
    complex(8), allocatable :: nl_projector_overlap(:,:,:,:) ! (nstate_frag,Nlma,nspin,ifrag_local)
    complex(8), allocatable :: nl_projector_overlap_halo(:,:,:,:) ! (nstate_frag,Nlma,nspin,halo)
    complex(8), allocatable :: nl_projector_r_overlap(:,:,:,:,:) ! (3,nstate_frag,Nlma,nspin,ifrag_local)
    complex(8), allocatable :: nl_projector_r_overlap_halo(:,:,:,:,:) ! (3,nstate_frag,Nlma,nspin,halo)
    real(8) :: Ac_nl_projector_cache(3) = 0.0d0
    integer :: nl_projector_cache_nlma = 0
    logical :: has_nl_projector_cache = .false.

    ! Observables
    real(8), allocatable :: esp(:,:)           ! eigenvalues (nstate_tot, nspin)
    real(8), allocatable :: eigenvalues(:,:)   ! eigenvalues per fragment (nstate_frag, nspin)
    real(8) :: dipole_moment(3)                ! total dipole moment
    real(8) :: current(3)                      ! current density
    real(8) :: current_para(3)                 ! paramagnetic current density
    real(8) :: current_para_source_norm(3)     ! density-scaled ||grad * C_occ|| diagnostic
    real(8) :: current_para_bound(3)           ! Cauchy bound for |<C|grad|C>|/volume
    real(8) :: current_coef_norm               ! density-scaled ||C_occ|| diagnostic
    real(8) :: current_coef_imag_norm          ! density-scaled ||Im(C_occ)|| diagnostic
    real(8) :: current_nl(3)                   ! nonlocal pseudopotential current density
    real(8) :: current_dia(3)                  ! diamagnetic current density proxy
    real(8) :: current_total(3)                ! total current density = para + nonlocal + dia
    real(8) :: current_momentum_self_raw(3) = 0.0d0  ! production momentum-block self contribution before /Ngrid
    real(8) :: current_momentum_cross_raw(3) = 0.0d0 ! production momentum-block cross contribution before /Ngrid
    logical :: current_momentum_decomp_ready = .false.
    real(8) :: polarization_lg(3)              ! length-gauge electronic polarization density
    real(8) :: polarization_lg_ref(3)          ! initial length-gauge polarization reference
    logical :: polarization_lg_ref_ready       ! true after initial reference is captured
    real(8) :: dipole_lg_raw(3)                ! unnormalized length-gauge electronic dipole
    logical :: mixed_z_local_prop_rho_ready = .false.
    logical :: mixed_z_local_prop_rho_bad = .true.
    integer :: mixed_z_local_prop_rho_step = -1
    real(8) :: mixed_z_local_prop_rho_prod_int = 0.0d0
    real(8) :: mixed_z_local_prop_rho_candidate_int = 0.0d0
    real(8) :: mixed_z_local_prop_rho_diff_int = 0.0d0
    real(8) :: mixed_z_local_prop_rho_max_abs_diff = 0.0d0
    real(8) :: mixed_z_local_prop_rho_rms_diff = 0.0d0
    logical :: mixed_z_local_prop_payload_ready = .false.
    logical :: mixed_z_local_prop_payload_bad = .true.
    integer :: mixed_z_local_prop_payload_step = -1
    real(8) :: mixed_z_local_prop_payload_coef_norm = 0.0d0
    real(8) :: mixed_z_local_prop_payload_prod_coef_norm = 0.0d0
    real(8) :: mixed_z_local_prop_payload_coef_diff_snorm = 0.0d0
    real(8) :: mixed_z_local_prop_payload_rel_coef_diff_snorm = 0.0d0
    real(8) :: mixed_z_local_prop_payload_dim = 0.0d0
    real(8) :: mixed_z_local_prop_payload_occ_trace = 0.0d0
    complex(8), allocatable :: mixed_z_local_prop_payload_wcoef(:,:,:) ! (global W,state,spin)
    complex(8), allocatable :: mixed_z_frag_local_wcoef(:,:,:) ! (owner-local W slot,state,spin)
    complex(8), allocatable :: mixed_z_frag_local_pself_coef(:,:,:) ! (owner-local P_self slot,state,spin)
    complex(8), allocatable :: mixed_z_frag_local_pneighbor_coef(:,:,:) ! (owner-local P_neighbor slot,state,spin)
    integer, allocatable :: mixed_z_frag_local_w_gid(:)
    integer, allocatable :: mixed_z_frag_local_pself_gid(:)
    integer, allocatable :: mixed_z_frag_local_pneighbor_gid(:)
    integer, allocatable :: mixed_z_frag_local_w_mix_gid(:)
    integer, allocatable :: mixed_z_frag_local_pself_mix_gid(:)
    integer, allocatable :: mixed_z_frag_local_pneighbor_mix_gid(:)
    integer :: mixed_z_frag_local_w_slots = 0
    integer :: mixed_z_frag_local_pself_slots = 0
    integer :: mixed_z_frag_local_pneighbor_slots = 0
    integer :: mixed_z_frag_local_nstate = 0
    integer :: mixed_z_frag_local_nspin = 0
    logical :: mixed_z_frag_local_storage_ready = .false.
    character(len=32) :: mixed_z_local_prop_payload_source = 'missing_payload'
    character(len=32) :: mixed_z_local_prop_payload_basis_kind = 'none'
    character(len=32) :: mixed_z_local_prop_payload_build_route = 'not_requested'
    character(len=32) :: mixed_z_local_prop_payload_block_reason = 'not_requested'
    integer(8) :: mixed_z_perf_nstep = 0_8
    integer(8) :: mixed_z_perf_prop_writeback_calls = 0_8
    integer(8) :: mixed_z_perf_rho_writeback_calls = 0_8
    integer(8) :: mixed_z_perf_pz_writeback_calls = 0_8
    integer(8) :: mixed_z_perf_current_writeback_calls = 0_8
    integer(8) :: mixed_z_perf_series_validation_calls = 0_8
    integer(8) :: mixed_z_perf_payload_current_calls = 0_8
    integer(8) :: mixed_z_perf_zgemm_calls = 0_8
    integer(8) :: mixed_z_perf_eigensolve_calls = 0_8
    integer(8) :: mixed_z_perf_expdiag_calls = 0_8
    integer(8) :: mixed_z_perf_mpi_reduce_calls = 0_8
    integer(8) :: mixed_z_perf_obs_mpi_reduce_calls = 0_8
    integer(8) :: mixed_z_perf_prop_pack_calls = 0_8
    integer(8) :: mixed_z_perf_prop_unpack_calls = 0_8
    integer(8) :: mixed_z_perf_prop_pack_allocs = 0_8
    integer(8) :: mixed_z_perf_prop_unpack_allocs = 0_8
    integer(8) :: mixed_z_perf_prop_pack_zero_bytes = 0_8
    integer(8) :: mixed_z_perf_prop_unpack_zero_bytes = 0_8
    integer(8) :: mixed_z_perf_prop_pack_w_copy_bytes = 0_8
    integer(8) :: mixed_z_perf_prop_pack_p_copy_bytes = 0_8
    integer(8) :: mixed_z_perf_prop_unpack_w_copy_bytes = 0_8
    integer(8) :: mixed_z_perf_prop_unpack_p_copy_bytes = 0_8
    integer :: mixed_z_perf_final_itt = -1
    logical :: mixed_z_perf_count_enabled = .false.
    logical :: mixed_z_frag_local_stability_baseline_ready = .false.
    character(len=16) :: mixed_z_frag_local_field_block_kind = 'all'
    real(8) :: mixed_z_frag_local_field_abs_max = 0.0d0
    real(8) :: mixed_z_frag_local_rho_baseline = 0.0d0
    real(8) :: mixed_z_frag_local_pz_baseline = 0.0d0
    real(8) :: mixed_z_frag_local_current_baseline = 0.0d0
    real(8) :: mixed_z_frag_local_norm_baseline = 0.0d0
    real(8) :: mixed_z_frag_local_energy_baseline = 0.0d0
    real(8) :: mixed_z_frag_local_rho_drift_max = 0.0d0
    real(8) :: mixed_z_frag_local_pz_drift_max = 0.0d0
    real(8) :: mixed_z_frag_local_current_drift_max = 0.0d0
    real(8) :: mixed_z_frag_local_norm_drift_max = 0.0d0
    real(8) :: mixed_z_frag_local_energy_drift_max = 0.0d0
    real(8) :: mixed_z_perf_wall_prop = 0.0d0
    real(8) :: mixed_z_perf_wall_prop_pack = 0.0d0
    real(8) :: mixed_z_perf_wall_prop_comm = 0.0d0
    real(8) :: mixed_z_perf_wall_prop_field_exp = 0.0d0
    real(8) :: mixed_z_perf_wall_prop_phase = 0.0d0
    real(8) :: mixed_z_perf_wall_prop_unpack = 0.0d0
    real(8) :: mixed_z_perf_wall_obs = 0.0d0
    real(8) :: mixed_z_perf_wall_series = 0.0d0
    real(8) :: mixed_z_perf_wall_obs_prepare = 0.0d0
    real(8) :: mixed_z_perf_wall_rho_writeback = 0.0d0
    real(8) :: mixed_z_perf_wall_pz_writeback = 0.0d0
    real(8) :: mixed_z_perf_wall_pz_prod = 0.0d0
    real(8) :: mixed_z_perf_wall_pz_reduced = 0.0d0
    real(8) :: mixed_z_perf_wall_pz_candidate = 0.0d0
    real(8) :: mixed_z_perf_wall_pz_build_cw = 0.0d0
    real(8) :: mixed_z_perf_wall_pz_contract_z = 0.0d0
    real(8) :: mixed_z_perf_wall_pz_reduce = 0.0d0
    real(8) :: mixed_z_perf_wall_current_writeback = 0.0d0
    real(8) :: mixed_z_perf_wall_payload_current = 0.0d0
    real(8) :: mixed_z_perf_wall_current_para = 0.0d0
    real(8) :: mixed_z_perf_wall_current_nl = 0.0d0
    real(8) :: mixed_z_perf_wall_current_nl_cache = 0.0d0
    real(8) :: mixed_z_perf_wall_current_nl_setup = 0.0d0
    real(8) :: mixed_z_perf_wall_current_nl_fetch = 0.0d0
    real(8) :: mixed_z_perf_wall_current_nl_project = 0.0d0
    real(8) :: mixed_z_perf_wall_current_nl_contract = 0.0d0
    real(8) :: mixed_z_perf_wall_current_nl_reduce = 0.0d0
    real(8) :: mixed_z_perf_wall_current_dia = 0.0d0
    real(8) :: mixed_z_perf_wall_obs_mpi_reduce = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_raw_z = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_ref_raw_z = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_diff = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_rel_diff = huge(1.0d0)
    real(8) :: mixed_z_local_pz_wcenter_weight = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_ref_weight = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_slot_count = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_block_count = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_owner_unique_raw_z = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_owner_unique_ref_raw_z = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_owner_unique_weight = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_owner_unique_ref_weight = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_owner_unique_count = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_missing_owner_count = 0.0d0
    real(8) :: mixed_z_local_pz_wcenter_duplicate_owner_count = 0.0d0
    real(8) :: mixed_z_local_pz_owner_unique_global_center_raw_z = 0.0d0
    real(8) :: mixed_z_local_pz_owner_unique_global_zww_diag_raw_z = 0.0d0
    real(8) :: mixed_z_local_pz_owner_unique_weighted_diff_global_center = 0.0d0
    real(8) :: mixed_z_local_pz_owner_unique_weighted_diff_global_zww_diag = 0.0d0
    real(8) :: mixed_z_local_pz_owner_unique_center_source_match_count = 0.0d0
    real(8) :: mixed_z_local_pz_owner_unique_wrap_match_count = 0.0d0
    real(8) :: mixed_z_local_pz_owner_unique_cell_shift_count = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_sum_z_prod = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_sum_z_lsp = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_sum_weight_prod = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_sum_weight_lsp = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_sum_contrib_prod = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_sum_contrib_lsp = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_max_abs_diff_z = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_max_abs_diff_weight = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_max_abs_diff_contrib = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_rms_diff_contrib = 0.0d0
    real(8) :: mixed_z_local_pz_zww_weight_owner_gid_mismatch_count = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_sum_occ_prod = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_sum_occ_lsp = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_sum_rho_prod = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_sum_rho_lsp = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_max_abs_diff_occ = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_max_abs_diff_rho = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_max_abs_diff_factor = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_max_abs_diff_weight = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_weighted_zww_prod = 0.0d0
    real(8) :: mixed_z_local_pz_weight_source_weighted_zww_lsp = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_prod_observable_abs2 = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_prod_expdiag_abs2 = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_source_abs2 = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_source_smetric = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_ref_after_abs2 = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_ref_after_smetric = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_lsp_after_abs2 = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_lsp_after_smetric = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_repacked_global_abs2 = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_repacked_global_smetric = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_max_abs_diff_abs2 = 0.0d0
    real(8) :: mixed_z_local_pz_rhodiag_max_abs_diff_smetric = 0.0d0
    logical :: mixed_z_local_pz_rhodiag_repacked_global_available = .false.
    real(8) :: mixed_z_local_pz_repacked_bridge_prod = 0.0d0
    real(8) :: mixed_z_local_pz_repacked_bridge_local_cref = 0.0d0
    real(8) :: mixed_z_local_pz_repacked_bridge_repacked_global = 0.0d0
    real(8) :: mixed_z_local_pz_repacked_bridge_rho_local_cref = 0.0d0
    real(8) :: mixed_z_local_pz_repacked_bridge_rho_repacked_global = 0.0d0
    real(8) :: mixed_z_local_pz_repacked_bridge_weight_repacked_global = 0.0d0
    real(8) :: mixed_z_local_pz_repacked_bridge_diff_prod_local = 0.0d0
    real(8) :: mixed_z_local_pz_repacked_bridge_diff_prod_repacked = 0.0d0
    logical :: mixed_z_local_pz_repacked_bridge_available = .false.
    logical :: mixed_z_local_pz_wcenter_ready = .false.
    logical :: mixed_z_local_pz_wcenter_bad = .true.
    real(8) :: mixed_z_prod_pz_ww_raw = 0.0d0
    real(8) :: mixed_z_prod_pz_ww_diag_raw = 0.0d0
    real(8) :: mixed_z_prod_pz_ww_offdiag_raw = 0.0d0
    real(8) :: mixed_z_prod_pz_ww_same_owner_raw = 0.0d0
    real(8) :: mixed_z_prod_pz_ww_cross_owner_raw = 0.0d0
    real(8) :: mixed_z_prod_pz_ww_pair_count_total = 0.0d0
    real(8) :: mixed_z_prod_pz_ww_pair_count_diag = 0.0d0
    real(8) :: mixed_z_prod_pz_ww_pair_count_offdiag = 0.0d0
    real(8) :: mixed_z_prod_pz_ww_pair_count_cross_owner = 0.0d0
    real(8) :: mixed_z_prod_zww_diag_sum = 0.0d0
    real(8) :: mixed_z_prod_center_diag_sum = 0.0d0
    real(8) :: mixed_z_prod_center_diag_local_sum = 0.0d0
    real(8) :: mixed_z_prod_zww_diag_min = 0.0d0
    real(8) :: mixed_z_prod_zww_diag_max = 0.0d0
    real(8) :: mixed_z_prod_zww_diag_mean = 0.0d0
    real(8) :: mixed_z_prod_center_z_min = 0.0d0
    real(8) :: mixed_z_prod_center_z_max = 0.0d0
    real(8) :: mixed_z_prod_center_z_mean = 0.0d0
    real(8) :: mixed_z_prod_diag_minus_center_min = 0.0d0
    real(8) :: mixed_z_prod_diag_minus_center_max = 0.0d0
    real(8) :: mixed_z_prod_diag_minus_center_rms = 0.0d0
    real(8) :: mixed_z_prod_weighted_zww_diag_sum = 0.0d0
    real(8) :: mixed_z_prod_weighted_center_sum = 0.0d0
    real(8) :: mixed_z_prod_weighted_diff_sum = 0.0d0
    real(8) :: mixed_z_prod_wrap_shift_count = 0.0d0
    real(8) :: mixed_z_prod_cell_shift_min = 0.0d0
    real(8) :: mixed_z_prod_cell_shift_max = 0.0d0
    real(8) :: mixed_z_prod_owner_gid_mismatch_count = 0.0d0
    real(8) :: mixed_z_prod_center_source_mismatch_count = 0.0d0
    real(8) :: mixed_z_prod_pz_wp_raw = 0.0d0
    real(8) :: mixed_z_prod_pz_pp_raw = 0.0d0
    real(8) :: mixed_z_prod_pz_occ_sum = 0.0d0
    real(8) :: mixed_z_prod_pz_w_occ_weight = 0.0d0
    real(8) :: mixed_z_prod_pz_w_dim = 0.0d0
    logical :: mixed_z_prod_pz_decomp_ready = .false.
    real(8) :: total_energy                    ! total energy
    real(8) :: energy_kinetic                  ! occupied expectation of kinetic block
    real(8) :: energy_nonlocal                 ! occupied expectation of nonlocal PP block
    real(8) :: energy_Ap                       ! occupied expectation of -i A.p term
    real(8) :: energy_A2                       ! occupied expectation of A^2/2 term
    real(8) :: elec_num_scaled                 ! total electrons from final real-space rho
    real(8) :: elec_num_raw                    ! raw reconstructed real-space electron count
    real(8) :: rho_scale_factor                ! kept at 1; RT density is not renormalized
    real(8) :: elec_num_raw_t0                 ! baseline raw electron count for drift metric
    real(8) :: elec_num_scaled_t0              ! baseline final electron count (diagnostic)
    logical :: elec_num_baseline_ready         ! baseline cache status
    real(8) :: rho_drift_indicator             ! time drift proxy using same convention: raw(t)-raw(t0)
    real(8) :: rho_ff_elec                     ! electron count from frag-frag block (c_f^dag S_ff c_f)
    real(8) :: rho_fp_elec                     ! electron count from frag-pw cross block
    real(8) :: rho_pp_elec                     ! electron count from pw-pw block (c_p^dag S_pp c_p)
    real(8) :: rho_owned_elec                  ! real-space integral from direct owner-local apply path
    real(8) :: rho_imported_elec               ! real-space integral from imported/unpacked contributions
    real(8) :: rho_local_all_elec              ! real-space integral over local rho before global reduction
    real(8) :: rho_global_raw_elec             ! global real-space integral (same as elec_num_raw)
    real(8) :: rho_global_scaled_elec          ! final global integral (same as elec_num_scaled)
    real(8) :: rho_contract_residual_elec      ! global_raw - (owned + imported)
    real(8) :: coef_occ_norm0                  ! occupied-coefficient norm baseline (itt=1)
    real(8) :: coef_occ_norm                   ! occupied-coefficient norm at current step
    real(8) :: coef_occ_norm_drift             ! occupied-coefficient norm drift from baseline
    real(8) :: csc_occ_identity_norm           ! Frobenius norm of (C^dagger S C - I) in occupied space
    real(8) :: csc_occ_identity_max            ! max absolute element of (C^dagger S C - I)
    real(8) :: occvirt_leakage                 ! occupied->virtual leakage wrt reference basis
    integer :: occvirt_top_occ                 ! top occupied state index (global) for leakage
    integer :: occvirt_top_virt                ! top virtual state index (global) for leakage pair
    real(8) :: occvirt_top_abs2                ! |U_vo|^2 of top leakage pair
    integer :: jpara_top_occ_state             ! top occupied state index contributing to J_para,z
    real(8) :: jpara_top_occ_value             ! signed contribution of top occupied state to J_para,z
    integer :: selfexc_track_states(3)         ! tracked state IDs for self-excitation diagnostics
    real(8) :: selfexc_track_norm(3)           ! tracked state coefficient norms at step sample
    real(8) :: selfexc_csc_step_delta          ! per-step increment of C^dagger S C - I norm
    real(8) :: selfexc_leak100_129_step_delta  ! per-step increment of 100->129 leakage amplitude
    real(8) :: selfexc_csc_stage_pre_mixed     ! C^dagger S C - I norm before mixed sync
    real(8) :: selfexc_csc_stage_post_overlap  ! C^dagger S C - I norm after overlap solve
    real(8) :: selfexc_csc_stage_post_raw      ! C^dagger S C - I norm after raw restore
    real(8) :: selfexc_leak100_129_pre_mixed   ! |U_129,100|^2 before mixed sync
    real(8) :: selfexc_leak100_129_post_overlap! |U_129,100|^2 after overlap solve
    real(8) :: selfexc_leak100_129_post_raw    ! |U_129,100|^2 after raw restore
    real(8) :: selfexc_leak100_129_current     ! current-step |U_129,100|^2 sample
    real(8) :: selfexc_csc_prev_step           ! previous-step C^dagger S C - I norm sample
    real(8) :: selfexc_leak100_129_prev_step   ! previous-step |U_129,100|^2 sample
    integer :: selfexc_prev_step_itt           ! step index cached for per-step delta
    real(8) :: pw_weight_raw                   ! simple diagnostic: sum |coef_pw|^2 over occupied states
    complex(8), allocatable :: coef_ref_all(:,:,:) ! local-row occ/virt diagnostic reference coefficients
    logical :: coef_ref_ready = .false.        ! reference coefficient cache status

    ! Fragment geometry information
    integer :: num_fragment(3)                 ! Fragment division in each direction (from salmon_global)
    logical :: has_wannier_cluster_partition = .false. ! DC provided cluster-level Wannier partition
    integer :: wannier_cluster_size(3) = 1     ! # of DG fragments grouped into one Wannier cluster
    integer :: num_wannier_cluster(3) = 0      ! Wannier cluster grid
    integer :: n_wannier_cluster = 0           ! total number of Wannier clusters
    integer, allocatable :: wannier_cluster_id(:) ! (n_frag) fragment -> Wannier cluster id
    integer, allocatable :: wannier_cluster_range(:,:) ! (6,n_wannier_cluster): xlo,xhi,ylo,yhi,zlo,zhi in fragment coords
    integer, allocatable :: nxyz_domain(:,:)   ! (3, n_frag) domain size of each fragment
    integer, allocatable :: ixyz_frag(:,:)     ! (3, n_frag) r-grid index of fragment origin
    integer, allocatable :: frag_core_lo(:,:)  ! (3, n_frag) unwrapped fragment-core lower bounds
    integer, allocatable :: frag_core_hi(:,:)  ! (3, n_frag) unwrapped fragment-core upper bounds
    integer, allocatable :: frag_buf_lo(:,:)   ! (3, n_frag) unwrapped fragment buffer-box lower bounds
    integer, allocatable :: frag_buf_hi(:,:)   ! (3, n_frag) unwrapped fragment buffer-box upper bounds
    integer :: rank_core_lo(3) = 0             ! unwrapped owner box lower bounds on this MPI rank
    integer :: rank_core_hi(3) = -1            ! unwrapped owner box upper bounds on this MPI rank
    integer :: rank_buf_lo(3) = 0              ! unwrapped local phi/rho buffer-box lower bounds on this MPI rank
    integer :: rank_buf_hi(3) = -1             ! unwrapped local phi/rho buffer-box upper bounds on this MPI rank
    integer, allocatable :: jxyz_tot(:,:)      ! r-grid mapping
    integer :: nxyz_buffer(3)                  ! # of halo points (4 for 4th-order stencil)
    integer, allocatable :: id_array(:)        ! (n_frag) MPI rank owning each fragment
    integer, allocatable :: density_owner_map(:,:,:,:) ! local-fragment interior grid -> owner rank
    logical, allocatable :: density_primary_local_map(:,:,:,:) ! local-fragment interior grid -> primary fragment selection
    integer, allocatable :: density_ixg_map(:,:,:,:)   ! local-fragment interior grid -> wrapped global x index
    integer, allocatable :: density_iyg_map(:,:,:,:)   ! local-fragment interior grid -> wrapped global y index
    integer, allocatable :: density_izg_map(:,:,:,:)   ! local-fragment interior grid -> wrapped global z index
    integer, allocatable :: density_send_count(:)      ! packed send length per target rank
    integer, allocatable :: density_send_slot_map(:,:,:,:) ! local-fragment grid -> packed send slot (0 if local)
    integer, allocatable :: density_subgroup_send_count(:)      ! packed subgroup-reduce length per target rank
    integer, allocatable :: density_subgroup_send_slot_map(:,:,:,:) ! local-fragment grid -> packed subgroup-reduce slot
    type(density_grid_point_info), allocatable :: density_grid_points(:,:) ! (max_points_per_local_fragment, ifrag_local)
    integer, allocatable :: density_grid_point_count(:)         ! valid point count per local fragment
    integer, allocatable :: density_grid_bx(:,:)                ! cached rho_bf x index for density_grid_points
    integer, allocatable :: density_grid_by(:,:)                ! cached rho_bf y index for density_grid_points
    integer, allocatable :: density_grid_bz(:,:)                ! cached rho_bf z index for density_grid_points
    integer :: density_rhobf_box_lo(3) = 0
    integer :: density_rhobf_box_hi(3) = -1
    logical :: density_rhobf_box_cache_valid = .false.
    integer, allocatable :: density_subgroup_self_ixg(:)        ! root-owned subgroup-reduce unpack x index
    integer, allocatable :: density_subgroup_self_iyg(:)        ! root-owned subgroup-reduce unpack y index
    integer, allocatable :: density_subgroup_self_izg(:)        ! root-owned subgroup-reduce unpack z index
    integer, allocatable :: density_block_nblocks(:)            ! subgroup density block count per local fragment
    integer, allocatable :: density_block_first_offset(:)       ! first subgroup-owned block offset per local fragment
    integer, allocatable :: density_block_step(:)               ! subgroup density block stride per local fragment
    integer, allocatable :: current_valid_grid_count(:)         ! valid local-current grid count per local fragment
    integer, allocatable :: current_density_grid_point_count(:) ! source density-grid count used to build current cache
    integer(8), allocatable :: current_density_grid_checksum(:) ! source density-grid signature used to build current cache
    integer, allocatable :: current_valid_ix(:,:)               ! valid fragment-local x indices for microscopic current
    integer, allocatable :: current_valid_iy(:,:)               ! valid fragment-local y indices for microscopic current
    integer, allocatable :: current_valid_iz(:,:)               ! valid fragment-local z indices for microscopic current
    integer, allocatable :: current_valid_ixg(:,:)              ! valid wrapped global x indices for microscopic current
    integer, allocatable :: current_valid_iyg(:,:)              ! valid wrapped global y indices for microscopic current
    integer, allocatable :: current_valid_izg(:,:)              ! valid wrapped global z indices for microscopic current
    type(density_recv_map_info), allocatable :: density_recv_map(:) ! packed recv unpack map per source rank
    real(8), allocatable :: density_phi_block_cache(:,:,:,:) ! block-packed basis (block_size,nstate_frag,nblock_max,ifrag_local)
    integer, allocatable :: density_phi_block_count(:)       ! block count per local fragment
    integer :: density_phi_block_size = 0
    logical :: density_phi_block_cache_valid = .false.
    complex(8), allocatable :: density_phase_block_cache(:,:,:,:) ! block-packed PW phases (block_size,n_pw,nblock_max,ifrag_local)
    integer :: density_phase_block_size = 0
    integer :: density_phase_block_npw = 0
    logical :: density_phase_block_cache_valid = .false.
    real(8), allocatable :: mixed_density_pp_cache(:,:,:,:) ! cached PP density (block_size,nblock_max,ifrag_local,nspin)
    integer :: mixed_density_pp_cache_mode = 0
    integer :: mixed_density_pp_cache_interval = 1
    logical :: mixed_density_pp_cache_valid = .false.
    integer :: lgnum_total(3)                  ! Total grid size (lg_tot%num)
    real(8) :: hgs(3)                           ! Grid spacing (a.u.)
    integer :: icomm                           ! MPI communicator for fragment RT
    integer :: id                              ! MPI rank in icomm
    integer :: isize                           ! MPI size in icomm
    integer :: icomm_frag                      ! MPI communicator for one fragment group
    integer :: id_frag                         ! MPI rank inside fragment group
    integer :: isize_frag                      ! MPI size inside fragment group
    integer :: ifrag_group                     ! global fragment index owned by this subgroup
    integer :: nproc_frag                      ! # of MPI ranks assigned to one fragment
    logical :: is_frag_root                    ! subgroup root responsible for fragment-global ownership
    character(32) :: parallel_mode = 'orbital' ! 'orbital' or 'legacy_realspace'
    logical :: parallel_mode_orbital = .true.

    type(halo_info), allocatable :: halo(:)    ! Halo regions (max 26 = 3^3-1 neighbors)
    integer :: n_halo                          ! Number of active halo regions
    integer, allocatable :: halo_ireq_send(:)  ! persistent send requests for halo exchange
    integer, allocatable :: halo_ireq_recv(:)  ! persistent recv requests for halo exchange

    ! Auxiliary arrays for self-consistent calculation
    real(8), allocatable :: phi_frag(:,:,:,:,:)    ! fragment basis functions in real space
                                                   ! With halo: (1-nb:nx+nb, 1-nb:ny+nb, 1-nb:nz+nb, nstate, ifrag)
                                                   ! where nb = nxyz_buffer = 4 for 4th-order stencil.
                                                   ! This is the shared real-space fragment basis itself.
                                                   ! Spin dependence enters through bookkeeping/operators, not here.
    complex(8), allocatable :: phi_frag_c(:,:,:,:,:) ! complex fragment basis for SOI/noncollinear path
    complex(8), allocatable :: phi_frag_spinor_c(:,:,:,:,:,:) ! (x,y,z,nspin,basis,ifrag_local) SOI spinor basis
    logical, allocatable :: phi_frag_has_seed_buffer(:) ! true when halo values came from DC buffer seed
    logical :: has_real_space_basis                ! flag: real-space basis available
    logical :: has_halo_exchange                   ! flag: halo exchange implemented

    ! Grid and potential storage for self-consistent updates
    type(s_rgrid), pointer :: lg => null()         ! large grid
    type(s_rgrid), pointer :: mg => null()         ! medium grid

    ! Density and potentials (allocated internally)
    real(8), allocatable :: rho_frag(:,:,:)        ! electron density on grid
    real(8), allocatable :: rho_s_frag(:,:,:,:)    ! spin-resolved density (ix,iy,iz,ispin)
    real(8), allocatable :: Vh_frag(:,:,:)         ! Hartree potential
    real(8), allocatable :: Vxc_frag(:,:,:,:)      ! XC potential (ix,iy,iz,ispin)
    real(8), allocatable :: stage_Vh_buffer(:,:,:) ! RK-stage scratch: fragment-buffered Hartree potential
    real(8), allocatable :: stage_Vpsl_buffer(:,:,:) ! RK-stage scratch: static local ionic potential
    real(8), allocatable :: stage_Vxc_buffer(:,:,:,:) ! RK-stage scratch: fragment-buffered XC potential
    real(8), allocatable :: H_ref_Vh_buffer(:,:,:) ! DC-LCFO seed Hartree reference for H_export + delta V
    real(8), allocatable :: H_ref_Vxc_buffer(:,:,:,:) ! DC-LCFO seed XC reference for H_export + delta V
    integer, allocatable :: stage_gx_map(:), stage_gy_map(:), stage_gz_map(:) ! buffer index -> parent-grid index
    logical :: stage_vpsl_buffer_valid = .false.  ! Vpsl is static during RT unless the scratch bounds change
    logical :: stage_map_valid = .false.
    logical :: H_delta_reference_valid = .false.

    ! Self-consistent basis update (adaptive basis)
    real(8), allocatable :: H_mat_kinetic(:,:,:)   ! Kinetic part only (constant, nstate,nstate,nspin)
    real(8) :: hamiltonian_change_norm             ! ||H_new - H_old|| Frobenius norm
    real(8) :: basis_update_threshold              ! Threshold for basis update (default: 0.1 a.u.)
    integer :: nbasis_update_count                 ! Number of times basis was updated
    integer :: last_basis_update_step              ! RT step index when last basis update was triggered
    logical :: yn_adaptive_basis                   ! Flag: enable adaptive basis updates

    ! Overlap matrix for basis projection (DC-LCFO returns real basis, no gauge rotation needed)
    real(8), allocatable :: basis_overlap(:,:,:)   ! <phi_new|phi_old> overlap matrix (real, no SVD)
    integer :: ifrag_start, ifrag_end     ! fragment distribution for this MPI rank

    ! DFT+U support
    logical :: use_plusu                  ! flag: DFT+U enabled
    complex(8), allocatable :: dm_plusu(:,:,:,:)  ! +U density matrix (m1,m2,ispin,ilma)
    real(8), allocatable :: U_mat(:,:,:)  ! +U potential matrix elements

    ! RK coefficients
    real(8), allocatable :: rk_alpha(:), rk_beta(:), rk_gamma(:)
    integer :: rk_stages

    ! ===================================================================
    ! Plane wave basis mixing (for metallic/delocalized states)
    ! ===================================================================
    logical :: use_plane_wave_basis           ! Enable plane wave mixing
    integer :: n_plane_waves                  ! Number of plane waves to add
    real(8) :: k_cutoff_pw                    ! Cutoff wave number for PWs (a.u.^-1)
    real(8), allocatable :: k_pw(:,:)         ! Wave vectors (3, n_plane_waves)
    complex(8), allocatable :: coef_pw(:,:,:) ! PW coefficients (n_plane_waves, nstate_tot, nspin)
    integer, allocatable :: coef_pw_owner(:)  ! (n_plane_waves), owning rank of each PW row
    complex(8), allocatable :: coef_pw_full_cache(:,:,:) ! subgroup-replicated PW coefficients for local-density evaluation
    integer :: coef_pw_full_cache_nstate = 0
    integer :: owned_coef_pw_start = 0         ! first owned PW row (contiguous hint)
    integer :: owned_coef_pw_end = -1          ! last owned PW row (contiguous hint)
    logical :: use_subspace_diag              ! Use compact DG subspace diagonalization path
    integer :: subspace_extra_states          ! Extra fragment states appended to occupied trial space
    integer :: subspace_pw_vectors            ! Number of PW helper vectors appended to trial space
    integer :: last_subspace_dim              ! Last accepted DG trial subspace dimension

    ! Mixed basis operators stored in FF/FP/PP form
    complex(8), allocatable :: S_mat_frag_pw(:,:,:) ! Overlap <fragment_i | PW_j>
    complex(8), allocatable :: S_mat_frag_pw_g(:,:,:) ! G-space overlap <fragment_i | PW_j>
    complex(8), allocatable :: H_mat_frag_pw(:,:,:) ! Hamiltonian <fragment_i | PW_j>
    complex(8), allocatable :: P_mat_frag_pw(:,:,:,:) ! Gradient <fragment_i | nabla_idir | PW_j>
                             ! (3, n_mat_max, n_plane_waves, nspin)
    integer, allocatable :: fp_local_row_ids(:) ! Local fragment-row ids for row-local FP operators
    integer, allocatable :: fp_local_pw_ids(:)  ! Requested PW row ids for row-local FP operators
    complex(8), allocatable :: S_mat_frag_pw_local(:,:,:) ! Row-local overlap <fragment_i | PW_j>
                             ! (n_local_fragment_rows, n_requested_pw_rows, nspin)
    complex(8), allocatable :: H_mat_frag_pw_local(:,:,:) ! Row-local Hamiltonian <fragment_i | PW_j>
                             ! (n_local_fragment_rows, n_requested_pw_rows, nspin)
    complex(8), allocatable :: P_mat_frag_pw_local(:,:,:,:) ! Row-local gradient <fragment_i | nabla_idir | PW_j>
                             ! (3, n_local_fragment_rows, n_requested_pw_rows, nspin)
    complex(8), allocatable :: P_mat_frag_pw_g(:,:,:,:) ! G-space gradient <fragment_i | nabla_idir | PW_j>
                  ! (3, n_mat_max, n_plane_waves, nspin)
    complex(8), allocatable :: H_mat_pw_diag(:,:)   ! PW-PW diagonal Hamiltonian block
    complex(8), allocatable :: H_mat_pw(:,:,:)      ! PW-PW Hamiltonian block (non-diagonal)
                             ! (n_plane_waves, n_plane_waves, nspin)
    logical :: has_wpw_window = .false.              ! Windowed-PW chi/grad(chi) are prepared
    integer, allocatable :: wpw_window_box_lo(:,:)   ! (3, n_local_frag), unwrapped storage lower bounds
    integer, allocatable :: wpw_window_box_hi(:,:)   ! (3, n_local_frag), unwrapped storage upper bounds
    real(8), allocatable :: wpw_chi(:,:,:,:)         ! (nx, ny, nz, n_local_frag), normalized chi_A
    real(8), allocatable :: wpw_grad_chi(:,:,:,:,:)  ! (3, nx, ny, nz, n_local_frag)
    real(8) :: wpw_partition_sum_chi2_min = 0.0d0
    real(8) :: wpw_partition_sum_chi2_max = 0.0d0
    real(8) :: wpw_partition_sum_chi2_maxdev = 0.0d0
    logical :: wpw_pp_blocks_ready = .false.          ! Fragment-local WPW PP block prototype is prepared
    type(complex_matrix_block_info), allocatable :: wpw_S_pp_blocks(:)
    type(complex_matrix_block_info), allocatable :: wpw_T_pp_volume_blocks(:)
    type(complex_matrix_block_info), allocatable :: wpw_T_pp_interface_blocks(:)
    real(8), allocatable :: pw_orthogonalized(:,:,:) ! [UNUSED] Orthogonalized PWs in real space
    logical :: mixed_basis_ready = .false.          ! startup/basis-update mixed diagonalization prepared
    logical :: mixed_basis_identity_raw = .false.   ! mixed_transform is raw FP+PW identity view
    integer, allocatable :: mixed_basis_dim(:)      ! retained mixed dimension per spin
    complex(8), allocatable :: mixed_transform(:,:,:) ! raw (F+P) basis <- orthonormal mixed basis
    complex(8), allocatable :: coef_mix(:,:,:)      ! canonical coefficients in orthonormal mixed basis

    ! ===================================================================
    ! RI/DF approximation for HSE exchange (Plan C)
    ! ===================================================================
    logical :: use_hse_ri                     ! Enable RI approximation for HSE
    ! Note: hse_ri_data type is defined in xc_hse_ri module
    ! Store as allocatable array indexed by fragment: hse_ri_data(ifrag_local)
    ! Actual type allocation handled in RT module after xc_hse_ri module is loaded

  end type s_dg_fragment_rt

contains

  subroutine invalidate_coef_exchange_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    if (allocated(dg_frag%coef_row_fragment)) deallocate(dg_frag%coef_row_fragment)
    if (allocated(dg_frag%coef_exchange_peer_ranks)) deallocate(dg_frag%coef_exchange_peer_ranks)
    if (allocated(dg_frag%coef_exchange_peer_count)) deallocate(dg_frag%coef_exchange_peer_count)
    if (allocated(dg_frag%coef_allowed_request_frag)) deallocate(dg_frag%coef_allowed_request_frag)
    dg_frag%coef_exchange_peer_cache_valid = .false.
  end subroutine invalidate_coef_exchange_cache

end module rt_dg_fragment_types
