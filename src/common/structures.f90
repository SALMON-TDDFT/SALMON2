!
!  Copyright 2019 SALMON developers
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

#include "config.h"

module structures
#ifdef USE_LIBXC
  use xc_f90_types_m
  use xc_f90_lib_m
#endif
  implicit none

! scalar field
  type s_scalar
    real(8),allocatable :: f(:,:,:) ! f(x,y,z)
  end type s_scalar

! vector field
  type s_vector
    real(8),allocatable :: v(:,:,:,:) ! v(1:3,x,y,z)
  end type s_vector

! density matrix rho(r,r')
  type s_dmatrix
    complex(8),allocatable :: zrho_mat(:,:,:,:,:,:,:) ! (ii,dir,x,y,z,ispin,im), ii=1~Nd, dir=1~6(xx,yy,zz,yz,zx,xy)
  end type s_dmatrix

  type s_dft_system
    integer :: iperiodic              ! iperiodic==0 --> isolated system, iperiodic==3 --> 3D periodic system
    integer :: ngrid,nspin,no,nk,nion ! # of r-grid points, spin indices, orbitals, k points, and ions
    real(8) :: hvol,hgs(3),primitive_a(3,3),det_a,primitive_b(3,3)
    real(8) :: rmatrix_a(3,3),rmatrix_b(3,3)
    real(8),allocatable :: vec_k(:,:)    ! (1:3,1:nk), k-vector
    real(8),allocatable :: wtk(:)        ! (1:nk), weight of k points
    real(8),allocatable :: rocc(:,:,:)   ! (1:no,1:nk,1:nspin), occupation rate
  ! atomic ...
    real(8),allocatable :: Mass(:)       ! (1:nelem), Atomic weight
    real(8),allocatable :: Rion(:,:)     ! (1:3,1:nion), atom position
    real(8),allocatable :: Velocity(:,:) ! (1:3,1:nion), atomic velocity
    real(8),allocatable :: Force(:,:)    ! (1:3,1:nion), force on atom
  ! external field
    real(8) :: vec_Ac(3) ! A/c (spatially averaged), A: vector potential, c: speed of light
    type(s_vector) :: Ac_micro ! A/c (microscopic)      ! for single-scale Maxwell-TDDFT
    type(s_scalar) :: div_Ac   ! divergence of Ac_micro ! for single-scale Maxwell-TDDFT
  end type s_dft_system

  type s_dft_energy
    real(8),allocatable :: esp(:,:,:) ! (1:no,1:nk,1:nspin), single-particle energy
    real(8) :: E_tot,E_kin,E_h,E_xc,E_ion_ion,E_ion_loc,E_ion_nloc
    real(8) :: E_U
  end type s_dft_energy

  type s_rgrid
    integer              :: ndir,Nd                 ! ndir=3 --> dir=xx,yy,zz, ndir=6 --> dir=xx,yy,zz,yz,zx,xy
    integer,dimension(3) :: is,ie,num &             ! num=ie-is+1
                           ,is_overlap,ie_overlap & ! is_overlap=is-Nd, ie_overlap=ie+Nd
                           ,is_array,ie_array       ! allocate( array(is_array(1):ie_array(1), ...) )
    integer ,allocatable :: idx(:),idy(:),idz(:)    ! idx(is_overlap(1):ie_overlap(1))=is_array(1)~ie_array(1), ...
    integer ,allocatable :: is_all(:,:),ie_all(:,:) ! (1:3,0:nproc-1), is & ie for all MPI processes
    real(8) ,allocatable :: coordinate(:,:)         ! (minval(is_overlap):maxval(ie_overlap),1:3), coordinate of grids 
  end type s_rgrid

! for persistent communication
  type s_pcomm_cache
    real(8), allocatable :: dbuf(:, :, :, :)
    complex(8), allocatable :: zbuf(:, :, :, :)
#ifdef FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION
!dir$ attributes align : 2097152 :: dbuf, zbuf
#endif
  end type s_pcomm_cache

! for update_overlap
  type s_sendrecv_grid
    ! Number of orbitals (4-th dimension of grid)
    integer :: nb
    ! Communicator
    integer :: icomm
    ! Neightboring MPI id (1:upside,2:downside, 1:x,2:y,3:z):
    integer :: neig(1:2, 1:3) 
    ! Communication requests (1:send,2:recv, 1:upside,2:downside, 1:x,2:y,3:z):
    integer :: ireq_real8(1:2, 1:2, 1:3)
    integer :: ireq_complex8(1:2, 1:2, 1:3)
    ! PComm cache (1:src/2:dst, 1:upside,2:downside, 1:x,2:y,3:z)
    type(s_pcomm_cache) :: cache(1:2, 1:2, 1:3)
    ! Range (axis=1...3, 1:src/2:dst, dir=1:upside,2:downside, dim=1:x,2:y,3:z)
    integer :: is_block(1:3, 1:2, 1:2, 1:3)
    integer :: ie_block(1:3, 1:2, 1:2, 1:3)
    ! Initialization flags
    logical :: if_pcomm_real8_initialized
    logical :: if_pcomm_complex8_initialized
  end type s_sendrecv_grid

  type s_orbital_parallel
    logical :: if_divide_rspace
    logical :: if_divide_orbit
    integer :: icomm_r,   id_r,   isize_r   ! communicator, process ID, & # of processes for r-space
    integer :: icomm_k,   id_k,   isize_k   ! communicator, process ID, & # of processes for k-space
    integer :: icomm_o,   id_o,   isize_o   ! communicator, process ID, & # of processes for orbital
    integer :: icomm_ro,  id_ro,  isize_ro  ! communicator, process ID, & # of processes for r-space & orbital
    integer :: icomm_ko,  id_ko,  isize_ko  ! communicator, process ID, & # of processes for k-space & orbital
    integer :: icomm_rko, id_rko, isize_rko ! communicator, process ID, & # of processes for r-space, k-space & orbital
    integer :: im_s,im_e,numm ! im=im_s,...,im_e, numm=im_e-im_s+1
    integer :: ik_s,ik_e,numk ! ik=ik_s,...,ik_e, numk=ik_e-ik_s+1
    integer :: io_s,io_e,numo ! io=io_s,...,io_e, numo=io_e-io_s+1
                              ! For calc_mode='RT' and temperature<0, these values are calculated from nelec.
                              ! In other cases, these are calculated from nstate.
    integer,allocatable :: irank_io(:) ! MPI rank of the orbital index #io
    integer :: imr(3) ! for sendrecv
  end type s_orbital_parallel

  type s_field_parallel
    integer :: icomm_all,id_all,isize_all ! communicator, process ID, & # of processes
    integer :: icomm(3)  ! 1: x-direction, 2: y-direction, 3: z-direction
    integer :: id(3), isize(3)
    integer :: icomm_ffte(3) ! 1: x-direction, 2: y-direction, 3: z-direction
                             ! Inside core FFTE routine, x-direction is redundant and
                             ! yz-direction is parallel.
    integer :: id_ffte(3), isize_ffte(3)
    integer :: imr(3),imrs(3) ! for sendrecv
    integer :: icomm_v,ngo(3),ngo_xyz,nproc_o ! for allgatherv_vlocal
  end type s_field_parallel

  type s_orbital
  ! ispin=1~nspin, io=io_s~io_e, ik=ik_s~ik_e, im=im_s~im_e (cf. s_orbital_parallel)
    real(8)   ,allocatable :: rwf(:,:,:,:,:,:,:) ! (ix,iy,iz,ispin,io,ik,im)
    complex(8),allocatable :: zwf(:,:,:,:,:,:,:) ! (ix,iy,iz,ispin,io,ik,im)
    complex(8),allocatable :: ztmp(:,:,:,:)
  end type s_orbital

  type s_stencil
    logical :: if_orthogonal
    real(8) :: coef_lap0,coef_lap(4,3),coef_nab(4,3) ! (4,3) --> (Nd,3) (future work)
    real(8) :: coef_f(6) ! for non-orthogonal lattice
  end type s_stencil

! pseudopotential
  type s_pp_info
    real(8) :: zion
    integer :: lmax,lmax0
    integer :: nrmax,nrmax0
    logical :: flag_nlcc
    character(2),allocatable :: atom_symbol(:)
    real(8),allocatable :: rmass(:)
    integer,allocatable :: mr(:)
    integer,allocatable :: lref(:)
    integer,allocatable :: nrps(:)
    integer,allocatable :: mlps(:)
    integer,allocatable :: nproj(:,:)
    integer,allocatable :: zps(:)
    integer,allocatable :: nrloc(:)
    real(8),allocatable :: rloc(:)
    real(8),allocatable :: rps(:)
    real(8),allocatable :: anorm(:,:)
    integer,allocatable :: inorm(:,:)
    real(8),allocatable :: anorm_so(:,:) ! '*_so' means what is used in 
    integer,allocatable :: inorm_so(:,:) !   spin-orbit calculation
    real(8),allocatable :: rad(:,:)
    real(8),allocatable :: radnl(:,:)
    real(8),allocatable :: vloctbl(:,:)
    real(8),allocatable :: dvloctbl(:,:)
    real(8),allocatable :: udvtbl(:,:,:)
    real(8),allocatable :: dudvtbl(:,:,:)
    real(8),allocatable :: rho_nlcc_tbl(:,:)
    real(8),allocatable :: tau_nlcc_tbl(:,:)
    real(8),allocatable :: upp_f(:,:,:)
    real(8),allocatable :: vpp_f(:,:,:)
    real(8),allocatable :: vpp_f_so(:,:,:)
    real(8),allocatable :: upp(:,:)
    real(8),allocatable :: dupp(:,:)
    real(8),allocatable :: vpp(:,:)
    real(8),allocatable :: dvpp(:,:)
    real(8),allocatable :: vpp_so(:,:)
    real(8),allocatable :: dvpp_so(:,:)
    real(8),allocatable :: udvtbl_so(:,:,:)
    real(8),allocatable :: dudvtbl_so(:,:,:)
    real(8),allocatable :: rps_ao(:)
    integer,allocatable :: nrps_ao(:)
    real(8),allocatable :: upptbl_ao(:,:,:)
    real(8),allocatable :: dupptbl_ao(:,:,:)
  end type s_pp_info

! pseudopotential on r-space grid
  type s_pp_grid
    integer :: nps
    integer,allocatable :: mps(:)
    integer,allocatable :: jxyz(:,:,:)
    integer,allocatable :: jxx(:,:)
    integer,allocatable :: jyy(:,:)
    integer,allocatable :: jzz(:,:)
    real(8),allocatable :: rxyz(:,:,:)
    real(8),allocatable :: uv(:,:)
    real(8),allocatable :: duv(:,:,:)
    integer :: nlma
    integer,allocatable :: lma_tbl(:,:)
    integer,allocatable :: ia_tbl(:)
    real(8),allocatable :: rinv_uvu(:)
    complex(8),allocatable :: zekr_uv(:,:,:) ! (j,ilma,ik), j=1~Mps(ia), ilma=1~Nlma, zekr_uV = exp(-i(k+A/c)r)*uv
    real(8),allocatable :: Vpsl_atom(:,:,:,:)
    !
    integer,allocatable :: ia_tbl_so(:)
    complex(8),allocatable :: uv_so(:,:,:,:)
    complex(8),allocatable :: duv_so(:,:,:,:,:)
    complex(8),allocatable :: zekr_uv_so(:,:,:,:,:)
    !
    integer,allocatable :: proj_pairs_ao(:,:)
    integer,allocatable :: proj_pairs_info_ao(:,:)
    integer,allocatable :: ia_tbl_ao(:)
    real(8),allocatable :: phi_ao(:,:)
    real(8),allocatable :: dphi_ao(:,:,:)
    complex(8),allocatable :: zekr_phi_ao(:,:,:)
    integer :: nps_ao
    integer,allocatable :: mps_ao(:)
    integer,allocatable :: jxyz_ao(:,:,:)
    integer,allocatable :: jxx_ao(:,:)
    integer,allocatable :: jyy_ao(:,:)
    integer,allocatable :: jzz_ao(:,:)
    real(8),allocatable :: rxyz_ao(:,:,:)
    real(8) :: Hvol
    ! for localized communication when calculating non-local pseudo-pt.
    integer,allocatable :: irange_atom(:,:)  ! uVpsi range for atom: n = (1,ia), m = (2,ia)
    logical,allocatable :: ireferred_atom(:) ! uVpsi(n:m) is referred in this process
    integer,allocatable :: icomm_atom(:)     ! communicator for uVpsi(n:m)
  end type s_pp_grid

  type s_pp_nlcc
    real(8), allocatable :: rho_nlcc(:,:,:)
    real(8), allocatable :: tau_nlcc(:,:,:)
  end type s_pp_nlcc

! exchange-correlation functional
  type s_xc_functional
    integer :: xctype(3)
    integer :: ispin
    real(8) :: cval
    logical :: use_gradient
    logical :: use_laplacian
    logical :: use_kinetic_energy
    logical :: use_current
#ifdef USE_LIBXC
    type(xc_f90_pointer_t) :: func(3)
    type(xc_f90_pointer_t) :: info(3)
#endif
  end type

  type s_reciprocal_grid
    integer :: icomm_G
    integer :: ng,iG_s,iG_e,iGzero
    real(8),allocatable :: Gx(:),Gy(:),Gz(:)
    complex(8),allocatable :: zrhoG_ion(:),zrhoG_ele(:),zdVG_ion(:,:)
    complex(8),allocatable :: zrhoG_ion_tmp(:),zrhoG_ele_tmp(:),zdVG_ion_tmp(:,:) ! work arrays
  end type s_reciprocal_grid

  type s_poisson
  ! for poisson_cg (conjugate-gradient method)
    integer :: iterVh                              ! iteration number for poisson_cg
    integer :: npole_partial                       ! number of multipoles calculated in each node
    integer :: npole_total                         ! total number of multipoles
    integer,allocatable :: ipole_tbl(:)            ! table for multipoles
    integer,allocatable :: ig_num(:)               ! number of grids for domains to which each multipole belongs
    integer,allocatable :: ig(:,:,:)               ! grid table for domains to which each multipole belongs
    integer,allocatable :: ig_bound(:,:,:)         ! grid table for boundaries
    real(8),allocatable :: wkbound(:), wkbound2(:) ! values on boundary represented in one-dimentional grid
  ! for discrete Fourier transform (general)
    complex(8),allocatable :: ff1(:,:,:),ff1x(:,:,:),ff1y(:,:,:),ff1z(:,:,:) &
                           & ,ff2(:,:,:),ff2x(:,:,:),ff2y(:,:,:),ff2z(:,:,:)
    real(8),allocatable    :: trho2z(:,:,:),trho3z(:,:,:)
    complex(8),allocatable :: egx(:,:),egxc(:,:),egy(:,:),egyc(:,:),egz(:,:),egzc(:,:)
  ! for FFTE
    real(8),allocatable :: coef(:,:,:)             ! coefficient of Poisson equation
    complex(8),allocatable :: a_ffte(:,:,:)        ! input matrix for Fourier transformation
    complex(8),allocatable :: a_ffte_tmp(:,:,:)    ! work array to make input matrix
    complex(8),allocatable :: b_ffte(:,:,:)        ! output matrix for Fourier transformation
  end type s_poisson

  type s_fdtd_system
    type(s_rgrid)         :: lg, mg, ng   ! Structure for send and receive in fdtd
    type(s_sendrecv_grid) :: srg_ng       ! Structure for send and receive in fdtd
    real(8) :: rlsize(3)                  ! Size of Cell
    real(8) :: hgs(3)                     ! Grid Spacing
    real(8) :: origin(3)                  ! Coordinate of Origin Point (TBA)
    character(8)  :: a_bc(3,2)            ! Boundary Condition for 1:x, 2:y, 3:z and 1:bottom and 2:top
    integer, allocatable :: imedia(:,:,:) ! Material information
  end type s_fdtd_system

  type s_fdtd_field
    type(s_scalar) :: phi, rho_em
    type(s_vector) :: vec_e, vec_h, vec_a, vec_j_em
    ! Experimental implementation
    type(s_vector) :: vec_Ac, vec_Ac_old
  end type s_fdtd_field

  type s_md
     real(8) :: Tene, Temperature, E_work, xi_nh
     real(8),allocatable :: Rion_last(:,:), Force_last(:,:)
  end type s_md

! output files
  type s_ofile
     integer :: fh_rt, fh_rt_energy
     character(256) :: file_rt_data, file_rt_energy_data
     character(256) :: dir_out_restart, dir_out_checkpoint
  end type s_ofile

! for DFT ground state calculations

  type s_cg
    type(s_orbital) :: xk,hxk,gk,pk,pko,hwf
  end type s_cg
  
  type s_mixing
    integer :: num_rho_stock
    type(s_scalar),allocatable :: srho_in(:), srho_out(:), srho_s_in(:,:), srho_s_out(:,:)
  end type s_mixing

!===================================================================================================================================

contains

  subroutine allocate_scalar(rg,field)
    implicit none
    type(s_rgrid),intent(in) :: rg
    type(s_scalar)           :: field
    integer :: ix,iy,iz
    allocate(field%f(rg%is(1):rg%ie(1),rg%is(2):rg%ie(2),rg%is(3):rg%ie(3)))
!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=rg%is(3),rg%ie(3)
    do iy=rg%is(2),rg%ie(2)
    do ix=rg%is(1),rg%ie(1)
      field%f(ix,iy,iz) = 0d0
    end do
    end do
    end do
  end subroutine allocate_scalar

  subroutine allocate_scalar_array(rg,field)
    implicit none
    type(s_rgrid),intent(in) :: rg
    type(s_scalar)           :: field
    integer :: ix,iy,iz
    allocate(field%f(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3)))
!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=rg%is_array(3),rg%ie_array(3)
    do iy=rg%is_array(2),rg%ie_array(2)
    do ix=rg%is_array(1),rg%ie_array(1)
      field%f(ix,iy,iz) = 0d0
    end do
    end do
    end do
  end subroutine allocate_scalar_array

  subroutine allocate_vector(rg,field)
    implicit none
    type(s_rgrid),intent(in) :: rg
    type(s_vector)           :: field
    integer :: ix,iy,iz
    allocate(field%v(3,rg%is(1):rg%ie(1),rg%is(2):rg%ie(2),rg%is(3):rg%ie(3)))
!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=rg%is(3),rg%ie(3)
    do iy=rg%is(2),rg%ie(2)
    do ix=rg%is(1),rg%ie(1)
      field%v = 0d0
    end do
    end do
    end do
  end subroutine allocate_vector

  subroutine allocate_dmatrix(nspin,mg,info,dmat)
    implicit none
    integer                 ,intent(in) :: nspin
    type(s_rgrid)           ,intent(in) :: mg
    type(s_orbital_parallel),intent(in) :: info
    type(s_dmatrix)                     :: dmat
    integer :: im,is,ix,iy,iz
    allocate(dmat%zrho_mat(mg%Nd,mg%ndir,mg%is(1)-mg%Nd:mg%ie(1),mg%is(2)-mg%Nd:mg%ie(2),mg%is(3)-mg%Nd:mg%ie(3), &
    & nspin,info%im_s:info%im_e))
!$omp parallel do collapse(4) private(im,is,iz,iy,ix)
    do im=info%im_s,info%im_e
    do is=1,nspin
    do iz=mg%is(3)-mg%Nd,mg%ie(3)
    do iy=mg%is(2)-mg%Nd,mg%ie(2)
    do ix=mg%is(1)-mg%Nd,mg%ie(1)
      dmat%zrho_mat(:,:,ix,iy,iz,is,im) = 0d0
    end do
    end do
    end do
    end do
    end do
  end subroutine allocate_dmatrix

  subroutine allocate_orbital_real(nspin,mg,info,psi)
    implicit none
    integer                 ,intent(in) :: nspin
    type(s_rgrid)           ,intent(in) :: mg
    type(s_orbital_parallel),intent(in) :: info
    type(s_orbital)                     :: psi
    integer :: im,ik,io,is,iz,iy,ix
    allocate(psi%rwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),  &
                     nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
!$omp parallel do collapse(6) private(im,ik,io,is,iz,iy,ix)
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do is=1,nspin
    do iz=mg%is_array(3),mg%ie_array(3)
    do iy=mg%is_array(2),mg%ie_array(2)
    do ix=mg%is_array(1),mg%ie_array(1)
      psi%rwf(ix,iy,iz,is,io,ik,im) = 0d0
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end subroutine allocate_orbital_real

  subroutine allocate_orbital_complex(nspin,mg,info,psi)
    implicit none
    integer                 ,intent(in) :: nspin
    type(s_rgrid)           ,intent(in) :: mg
    type(s_orbital_parallel),intent(in) :: info
    type(s_orbital)                     :: psi
    integer :: im,ik,io,is,iz,iy,ix
    allocate(psi%zwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),  &
                     nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
!$omp parallel do collapse(6) private(im,ik,io,is,iz,iy,ix)
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do is=1,nspin
    do iz=mg%is_array(3),mg%ie_array(3)
    do iy=mg%is_array(2),mg%ie_array(2)
    do ix=mg%is_array(1),mg%ie_array(1)
      psi%zwf(ix,iy,iz,is,io,ik,im) = 0d0
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end subroutine allocate_orbital_complex

!===================================================================================================================================

# define DEAL(x) if(allocated(x)) deallocate(x)

  subroutine deallocate_dft_system(system)
    type(s_dft_system) :: system
    DEAL(system%rocc)
    DEAL(system%wtk)
    DEAL(system%Rion)
    DEAL(system%Velocity)
    DEAL(system%Force)
  end subroutine deallocate_dft_system

  subroutine deallocate_dft_energy(energy)
    type(s_dft_energy) :: energy
    DEAL(energy%esp)
  end subroutine deallocate_dft_energy

  subroutine deallocate_rgrid(rg)
    type(s_rgrid) :: rg
    DEAL(rg%idx)
    DEAL(rg%idy)
    DEAL(rg%idz)
  end subroutine deallocate_rgrid

  subroutine deallocate_orbital_parallel(info)
    type(s_orbital_parallel) :: info
    DEAL(info%irank_io)
  end subroutine deallocate_orbital_parallel

  subroutine deallocate_orbital(psi)
    type(s_orbital) :: psi
    DEAL(psi%rwf)
    DEAL(psi%zwf)
    DEAL(psi%ztmp)
  end subroutine deallocate_orbital

  subroutine deallocate_pp_info(pp)
    type(s_pp_info) :: pp
    DEAL(pp%atom_symbol)
    DEAL(pp%rmass)
    DEAL(pp%mr)
    DEAL(pp%lref)
    DEAL(pp%nrps)
    DEAL(pp%mlps)
    DEAL(pp%zps)
    DEAL(pp%nrloc)
    DEAL(pp%rloc)
    DEAL(pp%rps)
    DEAL(pp%anorm)
    DEAL(pp%inorm)
    DEAL(pp%rad)
    DEAL(pp%radnl)
    DEAL(pp%vloctbl)
    DEAL(pp%dvloctbl)
    DEAL(pp%udvtbl)
    DEAL(pp%dudvtbl)
    DEAL(pp%rho_nlcc_tbl)
    DEAL(pp%tau_nlcc_tbl)
    DEAL(pp%upp_f)
    DEAL(pp%vpp_f)
    DEAL(pp%upp)
    DEAL(pp%dupp)
    DEAL(pp%vpp)
    DEAL(pp%dvpp)
  end subroutine deallocate_pp_info

  subroutine deallocate_pp_grid(ppg)
    type(s_pp_grid) :: ppg
    DEAL(ppg%mps)
    DEAL(ppg%jxyz)
    DEAL(ppg%jxx)
    DEAL(ppg%jyy)
    DEAL(ppg%jzz)
    DEAL(ppg%uv)
    DEAL(ppg%duv)
    DEAL(ppg%lma_tbl)
    DEAL(ppg%ia_tbl)
    DEAL(ppg%rinv_uvu)
    DEAL(ppg%zekr_uV)
    DEAL(ppg%Vpsl_atom)
    DEAL(ppg%uv_so)
    DEAL(ppg%duv_so)
  end subroutine deallocate_pp_grid

  subroutine deallocate_scalar(x)
    type(s_scalar) :: x
    DEAL(x%f)
  end subroutine deallocate_scalar

  subroutine deallocate_vector(x)
    type(s_vector) :: x
    DEAL(x%v)
  end subroutine deallocate_vector

  subroutine deallocate_dmatrix(dm)
    type(s_dmatrix) :: dm
    DEAL(dm%zrho_mat)
  end subroutine deallocate_dmatrix

  subroutine deallocate_reciprocal_grid(fg)
    type(s_reciprocal_grid) :: fg
    DEAL(fg%Gx)
    DEAL(fg%Gy)
    DEAL(fg%Gz)
    DEAL(fg%zrhoG_ion)
    DEAL(fg%zrhoG_ele)
    DEAL(fg%zdVG_ion)
  end subroutine deallocate_reciprocal_grid

end module structures
