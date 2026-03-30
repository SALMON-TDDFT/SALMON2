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
  public :: halo_info, matrix_block_info, complex_matrix_block_info, vector_block_info, momentum_block_info, &
            density_recv_map_info, density_grid_point_info, real_buffer_info, s_dg_fragment_rt

  ! Halo communication structure (for phi_frag exchange between fragments)
  type :: halo_info
    integer :: id_src, id_dst                ! MPI ranks for communication
    integer :: ifrag_src                     ! Source fragment index
    integer :: ifrag_dst                     ! Destination (local) fragment index
    integer :: dvec(3)                       ! Direction vector to neighbor (-1,0,+1)
    integer :: length(3)                     ! Size of halo region in each direction
    integer :: dsp_send(3), dsp_recv(3)      ! Displacement for send/recv buffers
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

  type :: momentum_block_info
    integer :: ifrag_row = 0
    integer :: ifrag_col = 0
    integer :: nrow_max = 0
    integer :: ncol_max = 0
    real(8), allocatable :: val(:,:,:,:) ! (3, nrow_max, ncol_max, nspin)
  end type momentum_block_info

  type :: density_recv_map_info
    integer :: npts = 0
    integer, allocatable :: ixg(:)
    integer, allocatable :: iyg(:)
    integer, allocatable :: izg(:)
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
    integer :: time_integrator                 ! 1: SSPRK3, 2: AETRS, 3: RK4

    ! Spin-resolved coefficients on the fragment basis (MUST be complex).
    ! The first index is the global fragment-basis index after index_basis mapping,
    ! not a separate real-space basis copy for each spin channel.
    complex(8), allocatable :: coef(:,:,:)        ! (nstate_frag, nstate_tot, nspin)
    complex(8), allocatable :: coef_new(:,:,:)    ! for time propagation
    complex(8), allocatable :: coef_work(:,:,:)   ! work array
    integer, allocatable :: coef_owner(:,:)       ! (n_mat_max, nspin), owning rank of each fragment-basis row
    integer :: owned_coef_start = 0               ! first owned fragment-basis row (contiguous hint)
    integer :: owned_coef_end = -1                ! last owned fragment-basis row (contiguous hint)

    ! Spin-channel-resolved operator matrices on the fragment basis.
    ! These arrays are spin-resolved projected representations, even when the
    ! underlying real-space fragment basis functions are spin independent.
    real(8), allocatable :: H_mat(:,:,:)       ! (nmat, nmat, nspin)
    complex(8), allocatable :: H_mat_c(:,:,:)      ! complex Hamiltonian (SOI/mixed propagation path)
    type(matrix_block_info), allocatable :: H_mat_blocks(:)
    type(matrix_block_info), allocatable :: H_mat_kinetic_blocks(:)
    type(complex_matrix_block_info), allocatable :: H_nl_blocks(:)
    integer, allocatable :: H_block_map(:,:)
    integer :: n_H_blocks = 0
    integer, allocatable :: H_nl_block_map(:,:)
    integer :: n_H_nl_blocks = 0
    integer, allocatable :: H_local_block_ids(:) ! row-owner-local H block ids for RT apply
    integer, allocatable :: H_nl_local_block_ids(:) ! row-owner-local nonlocal H block ids for RT apply
    integer, allocatable :: H_local_rows(:)      ! fragment rows owned by this rank
    logical :: H_local_block_ids_valid = .false.
    logical, allocatable :: runtime_neighbor_pair_cache(:,:)    ! static fragment-pair runtime adjacency
    logical, allocatable :: momentum_neighbor_pair_cache(:,:)   ! static fragment-pair momentum adjacency
    real(8), allocatable :: S_mat(:,:,:)       ! raw fragment overlap matrix
    real(8), allocatable :: S_mat_prop(:,:,:)  ! overlap matrix used in propagation/unitarity
    complex(8), allocatable :: S_mat_c(:,:,:)      ! complex raw fragment overlap matrix (SOI path)
    complex(8), allocatable :: S_mat_prop_c(:,:,:) ! complex overlap used in propagation/unitarity
    type(matrix_block_info), allocatable :: S_mat_blocks(:)
    type(matrix_block_info), allocatable :: S_mat_prop_blocks(:)
    integer, allocatable :: S_block_map(:,:)
    integer :: n_S_blocks = 0
    logical :: has_global_overlap_copy = .true.
    logical :: overlap_prop_root_authoritative = .false.
    integer, allocatable :: n_basis(:,:)       ! (n_frag, nspin), spin-resolved active basis count per fragment
    integer, allocatable :: index_basis(:,:,:) ! (nstate_frag, n_frag, nspin), spin-resolved local->global basis map
    integer, allocatable :: n_mat(:)           ! (nspin), spin-resolved projected-matrix dimension
    integer :: n_mat_max                       ! max projected-matrix dimension over spin channels

    ! Time-dependent external field coupling (velocity gauge: H = H_0 - i*A·∇ + A^2/2)
    real(8), allocatable :: momentum_mat(:,:,:,:) ! momentum matrix elements p_ij = <phi_i|p|phi_j> (x,y,z)
    complex(8), allocatable :: momentum_mat_c(:,:,:,:) ! complex momentum matrix for SOI/mixed propagation
    type(vector_block_info), allocatable :: momentum_blocks(:)
    integer, allocatable :: momentum_block_map(:,:)
    integer :: n_momentum_blocks = 0
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
    complex(8), allocatable :: H_nl_cache(:,:,:)   ! cached non-local PP matrix (A-dependent)
    real(8) :: Ac_nl_cache(3)                      ! cached vector potential for H_nl_cache
    real(8) :: Ac_nl_cache_tol                     ! tolerance for cache reuse
    logical :: has_nl_cache                        ! flag: cached H_nl available
    complex(8), allocatable :: nl_pp_phi_self(:,:,:,:) ! (nps,natom,nstate_frag,ifrag_local)
    complex(8), allocatable :: nl_pp_phi_halo(:,:,:,:) ! (nps,natom,nstate_frag,n_halo)
    logical :: nl_pp_phi_cache_valid = .false.

    ! Observables
    real(8), allocatable :: esp(:,:)           ! eigenvalues (nstate_tot, nspin)
    real(8), allocatable :: eigenvalues(:,:)   ! eigenvalues per fragment (nstate_frag, nspin)
    real(8) :: dipole_moment(3)                ! total dipole moment
    real(8) :: current(3)                      ! current density
    real(8) :: total_energy                    ! total energy
    real(8) :: elec_num_scaled                 ! total electrons after rho normalization
    real(8) :: elec_num_raw                    ! total electrons before rho normalization
    real(8) :: rho_scale_factor                ! density renormalization factor
    real(8) :: pw_weight_raw                   ! simple diagnostic: sum |coef_pw|^2 over occupied states

    ! Fragment geometry information
    integer :: num_fragment(3)                 ! Fragment division in each direction (from salmon_global)
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
    integer, allocatable :: density_subgroup_self_ixg(:)        ! root-owned subgroup-reduce unpack x index
    integer, allocatable :: density_subgroup_self_iyg(:)        ! root-owned subgroup-reduce unpack y index
    integer, allocatable :: density_subgroup_self_izg(:)        ! root-owned subgroup-reduce unpack z index
    integer, allocatable :: density_block_nblocks(:)            ! subgroup density block count per local fragment
    integer, allocatable :: density_block_first_offset(:)       ! first subgroup-owned block offset per local fragment
    integer, allocatable :: density_block_step(:)               ! subgroup density block stride per local fragment
    integer, allocatable :: current_valid_grid_count(:)         ! valid local-current grid count per local fragment
    integer, allocatable :: current_valid_ix(:,:)               ! valid fragment-local x indices for microscopic current
    integer, allocatable :: current_valid_iy(:,:)               ! valid fragment-local y indices for microscopic current
    integer, allocatable :: current_valid_iz(:,:)               ! valid fragment-local z indices for microscopic current
    integer, allocatable :: current_valid_ixg(:,:)              ! valid wrapped global x indices for microscopic current
    integer, allocatable :: current_valid_iyg(:,:)              ! valid wrapped global y indices for microscopic current
    integer, allocatable :: current_valid_izg(:,:)              ! valid wrapped global z indices for microscopic current
    type(density_recv_map_info), allocatable :: density_recv_map(:) ! packed recv unpack map per source rank
    real(8), allocatable :: density_weight_local(:,:,:) ! static overlap count on local grid
    real(8), allocatable :: density_inv_weight_local(:,:,:) ! inverse of static overlap count on local grid
    real(8), allocatable :: density_phi_block_cache(:,:,:,:) ! block-packed basis (block_size,nstate_frag,nblock_max,ifrag_local)
    integer, allocatable :: density_phi_block_count(:)       ! block count per local fragment
    integer :: density_phi_block_size = 0
    logical :: density_phi_block_cache_valid = .false.
    complex(8), allocatable :: density_matrix_frag(:,:,:,:) ! raw fragment density matrix cache (nbf,nbf,nspin,ifrag_local)
    logical, allocatable :: density_matrix_frag_valid(:,:)  ! raw fragment density matrix validity (nspin,ifrag_local)
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

    ! Self-consistent basis update (adaptive basis)
    complex(8), allocatable :: H_mat_old(:,:,:)    ! Previous Hamiltonian matrix (nstate,nstate,nspin)
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
    real(8) :: subspace_fallback_cond         ! Threshold that triggers full dense fallback
    integer :: last_subspace_dim              ! Last accepted DG trial subspace dimension
    integer :: subspace_fallback_count        ! Number of times full dense fallback was used

    ! Mixed basis operators stored in FF/FP/PP form
    complex(8), allocatable :: S_mat_frag_pw(:,:,:) ! Overlap <fragment_i | PW_j>
    complex(8), allocatable :: H_mat_frag_pw(:,:,:) ! Hamiltonian <fragment_i | PW_j>
    complex(8), allocatable :: H_mat_pw_diag(:,:)   ! PW-PW diagonal Hamiltonian block
    real(8), allocatable :: pw_orthogonalized(:,:,:) ! [UNUSED] Orthogonalized PWs in real space
    logical :: mixed_basis_ready = .false.          ! startup/basis-update mixed diagonalization prepared
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

end module rt_dg_fragment_types
