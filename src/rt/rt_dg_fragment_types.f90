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
  public :: halo_info, s_dg_fragment_rt

  ! Halo communication structure (for phi_frag exchange between fragments)
  type, public :: halo_info
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

  ! Fragment basis data structure
  type, public :: s_dg_fragment_rt
    integer :: n_frag                          ! number of fragments
    integer :: nspin                           ! number of spin
    integer :: nstate_frag                     ! number of states per fragment
    integer :: nstate_tot                      ! total number of states
    integer :: time_integrator                 ! 1: SSPRK3, 2: AETRS, 3: RK4

    ! Fragment basis coefficients and their time derivatives (MUST be complex)
    complex(8), allocatable :: coef(:,:,:)        ! (nstate_frag, nstate_tot, nspin)
    complex(8), allocatable :: coef_new(:,:,:)    ! for time propagation
    complex(8), allocatable :: coef_work(:,:,:)   ! work array

    ! Hamiltonian matrix in fragment basis
    real(8), allocatable :: H_mat(:,:,:)       ! (nmat, nmat, nspin)
    complex(8), allocatable :: H_mat_c(:,:,:)      ! complex Hamiltonian (SOI/mixed propagation path)
    real(8), allocatable :: S_mat(:,:,:)       ! raw fragment overlap matrix
    real(8), allocatable :: S_mat_prop(:,:,:)  ! overlap matrix used in propagation/unitarity
    complex(8), allocatable :: S_mat_c(:,:,:)      ! complex raw fragment overlap matrix (SOI path)
    complex(8), allocatable :: S_mat_prop_c(:,:,:) ! complex overlap used in propagation/unitarity
    integer, allocatable :: n_basis(:,:)       ! (n_frag, nspin)
    integer, allocatable :: index_basis(:,:,:) ! (nstate_frag, n_frag, nspin)
    integer, allocatable :: n_mat(:)           ! (nspin)
    integer :: n_mat_max                       ! max dimension of H_mat

    ! Time-dependent external field coupling (velocity gauge: H = H_0 - i*A·∇ + A^2/2)
    real(8), allocatable :: momentum_mat(:,:,:,:) ! momentum matrix elements p_ij = <phi_i|p|phi_j> (x,y,z)
    complex(8), allocatable :: momentum_mat_c(:,:,:,:) ! complex momentum matrix for SOI/mixed propagation
    real(8), allocatable :: dipole_mat(:,:,:,:)   ! dipole matrix elements for observables (x,y,z)
    complex(8), allocatable :: H_nl_cache(:,:,:)   ! cached non-local PP matrix (A-dependent)
    real(8) :: Ac_nl_cache(3)                      ! cached vector potential for H_nl_cache
    real(8) :: Ac_nl_cache_tol                     ! tolerance for cache reuse
    logical :: has_nl_cache                        ! flag: cached H_nl available

    ! Observables
    real(8), allocatable :: esp(:,:)           ! eigenvalues (nstate_tot, nspin)
    real(8), allocatable :: eigenvalues(:,:)   ! eigenvalues per fragment (nstate_frag, nspin)
    real(8) :: dipole_moment(3)                ! total dipole moment
    real(8) :: current(3)                      ! current density
    real(8) :: total_energy                    ! total energy

    ! Fragment geometry information
    integer :: num_fragment(3)                 ! Fragment division in each direction (from salmon_global)
    integer, allocatable :: nxyz_domain(:,:)   ! (3, n_frag) domain size of each fragment
    integer, allocatable :: ixyz_frag(:,:)     ! (3, n_frag) r-grid index of fragment origin
    integer, allocatable :: jxyz_tot(:,:)      ! r-grid mapping
    integer :: nxyz_buffer(3)                  ! # of halo points (4 for 4th-order stencil)
    integer, allocatable :: id_array(:)        ! (n_frag) MPI rank owning each fragment
    integer :: lgnum_total(3)                  ! Total grid size (lg_tot%num)
    real(8) :: hgs(3)                           ! Grid spacing (a.u.)
    integer :: icomm                           ! MPI communicator for fragment RT
    integer :: id                              ! MPI rank in icomm
    integer :: isize                           ! MPI size in icomm

    type(halo_info), allocatable :: halo(:)    ! Halo regions (max 26 = 3^3-1 neighbors)
    integer :: n_halo                          ! Number of active halo regions

    ! Auxiliary arrays for self-consistent calculation
    real(8), allocatable :: phi_frag(:,:,:,:,:)    ! fragment basis functions in real space
                                                   ! With halo: (1-nb:nx+nb, 1-nb:ny+nb, 1-nb:nz+nb, nstate, ifrag)
                                                   ! where nb = nxyz_buffer = 4 for 4th-order stencil
                                                   ! REAL from DC-LCFO (no complex phase)
    complex(8), allocatable :: phi_frag_c(:,:,:,:,:) ! complex fragment basis (SOI path)
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

    ! Mixed basis Hamiltonian and overlap matrices
    complex(8), allocatable :: H_mat_mixed(:,:,:)     ! Full Hamiltonian (fragment + PW)
    complex(8), allocatable :: S_mat_mixed_prop(:,:,:) ! propagation overlap in mixed basis
    complex(8), allocatable :: S_mat_frag_pw(:,:,:) ! Overlap <fragment_i | PW_j>
    real(8), allocatable :: pw_orthogonalized(:,:,:) ! [UNUSED] Orthogonalized PWs in real space

    ! ===================================================================
    ! RI/DF approximation for HSE exchange (Plan C)
    ! ===================================================================
    logical :: use_hse_ri                     ! Enable RI approximation for HSE
    ! Note: hse_ri_data type is defined in xc_hse_ri module
    ! Store as allocatable array indexed by fragment: hse_ri_data(ifrag_local)
    ! Actual type allocation handled in RT module after xc_hse_ri module is loaded

  end type s_dg_fragment_rt

end module rt_dg_fragment_types
