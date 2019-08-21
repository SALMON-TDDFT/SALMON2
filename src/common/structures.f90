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
module structures
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
    real(8),allocatable :: vec_k(:,:)    ! (1:3,1:nk), k-vector
    real(8),allocatable :: wtk(:)        ! (1:nk), weight of k points
    real(8),allocatable :: rocc(:,:,:)   ! (1:no,1:nk,1:nspin), occupation rate
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
  end type s_dft_energy

  type s_rgrid
    integer              :: ndir,Nd                 ! ndir=3 --> dir=xx,yy,zz, ndir=6 --> dir=xx,yy,zz,yz,zx,xy
    integer,dimension(3) :: is,ie,num &             ! num=ie-is+1
                           ,is_overlap,ie_overlap & ! is_overlap=is-Nd, ie_overlap=ie+Nd
                           ,is_array,ie_array       ! allocate( array(is_array(1):ie_array(1), ...) )
    integer ,allocatable :: idx(:),idy(:),idz(:)    ! idx(is_overlap(1):ie_overlap(1))=is_array(1)~ie_array(1), ...
  end type s_rgrid

  type s_pcomm_cache
    real(8), allocatable :: dbuf(:, :, :, :)
    complex(8), allocatable :: zbuf(:, :, :, :)

#ifdef SALMON_ENABLE_2MB_ALIGNED_ALLOCATE
!dir$ attributes align : 2097152 :: dbuf, zbuf
#endif
  end type s_pcomm_cache

  type s_sendrecv_grid
    ! Number of orbitals (4-th dimension of grid)
    integer :: nb
    ! Communicator
    integer :: icomm
    ! Neightboring MPI id (1:x,2:y,3:z, 1:upside,2:downside):
    integer :: neig(1:3, 1:2) 
    ! Communication requests (1:x,2:y,3:z, 1:upside,2:downside, 1:send,2:recv):
    integer :: ireq_real8(1:3, 1:2, 1:2)
    integer :: ireq_complex8(1:3, 1:2, 1:2)
    ! PComm cache (1:x,2:y,3:z, 1:upside,2:downside, 1:src/2:dst)
    type(s_pcomm_cache) :: cache(1:3, 1:2, 1:2)
    ! Range (dim=1:x,2:y,3:z, dir=1:upside,2:downside, 1:src/2:dst, axis=1...3)
    integer :: is_block(1:3, 1:2, 1:2, 1:3)
    integer :: ie_block(1:3, 1:2, 1:2, 1:3)
    ! Initialization flags
    logical :: if_pcomm_real8_initialized
    logical :: if_pcomm_complex8_initialized
  end type s_sendrecv_grid

  type s_orbital_parallel
    logical :: if_divide_rspace
    logical :: if_divide_orbit
    integer :: icomm_r,   id_r,   isize_r   ! communicator for r-space
    integer :: icomm_k,   id_k,   isize_k   ! communicator for k-space
    integer :: icomm_o,   id_o,   isize_o   ! communicator for orbital
    integer :: icomm_ro,  id_ro,  isize_ro  ! communicator for r-space & orbital
    integer :: icomm_ko,  id_ko,  isize_ko  ! communicator for k-space & orbital
    integer :: icomm_rko ! communicator for r-space, k-space & orbital
    integer :: im_s,im_e,numm ! im=im_s,...,im_e, numm=im_e-im_s+1
    integer :: ik_s,ik_e,numk ! ik=ik_s,...,ik_e, numk=ik_e-ik_s+1
    integer :: io_s,io_e,numo ! io=io_s,...,io_e, numo=io_e-io_s+1
    integer,allocatable :: io_tbl(:)  ! jo=io_tbl(io), io=io_s~io_e, jo=1~no
    integer,allocatable :: jo_tbl(:)  ! io=io_tbl(jo), jo=1~no, io=io_s~io_e
    integer,allocatable :: irank_jo(:) ! MPI rank of the orbital index #jo
    real(8),allocatable :: occ(:,:,:,:) ! (io,ik,ispin,im), occ = rocc*wk, occupation numbers
  end type s_orbital_parallel

  type s_field_parallel
    integer :: icomm(3)  ! 1: x-direction, 2: y-direction, 3: z-direction
    integer :: id(3), isize(3)
  end type s_field_parallel

  type s_orbital
  ! ispin=1~nspin, io=io_s~io_e, ik=ik_s~ik_e, im=im_s~im_e (cf. s_orbital_parallel)
    real(8)   ,allocatable :: rwf(:,:,:,:,:,:,:) ! (ix,iy,iz,ispin,io,ik,im)
    complex(8),allocatable :: zwf(:,:,:,:,:,:,:) ! (ix,iy,iz,ispin,io,ik,im)
    complex(8),allocatable :: ztmp(:,:,:,:)
  end type s_orbital

  type s_stencil
    real(8) :: coef_lap0,coef_lap(4,3),coef_nab(4,3) ! (4,3) --> (Nd,3) (future work)
    real(8),allocatable :: vec_kAc(:,:) ! (1:3,ik_s:ik_e)
    logical :: if_orthogonal
    real(8) :: rmatrix_a(3,3),rmatrix_b(3,3),coef_f(6) ! for non-orthogonal lattice
    integer,allocatable :: isign(:,:) ! sign(3,4:ndir) (for ndir=4~6)

  ! Experimental implementation of srg
    type(s_sendrecv_grid) :: srg
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
    real(8),allocatable :: upp(:,:)
    real(8),allocatable :: dupp(:,:)
    real(8),allocatable :: vpp(:,:)
    real(8),allocatable :: dvpp(:,:)
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

#ifdef SALMON_ENABLE_MPI3
    integer,allocatable :: irange_atom(:,:)  ! uVpsi range for atom: n = (1,ia), m = (2,ia)
    logical,allocatable :: ireferred_atom(:) ! uVpsi(n:m) is referred in this process
    integer,allocatable :: icomm_atom(:)     ! communicator for uVpsi(n:m)
#endif
  end type s_pp_grid

  type s_pp_nlcc
    real(8), allocatable :: rho_nlcc(:,:,:)
    real(8), allocatable :: tau_nlcc(:,:,:)
  end type s_pp_nlcc

  type s_reciprocal_grid
    integer :: icomm_G
    integer :: ng,iG_s,iG_e,iGzero
    real(8),allocatable :: Gx(:),Gy(:),Gz(:)
    complex(8),allocatable :: zrhoG_ion(:),zrhoG_ele(:),zdVG_ion(:,:)
  end type s_reciprocal_grid

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
  end type s_fdtd_field

  type s_md
     real(8) :: Tene, Temperature, E_work, xi_nh
     real(8),allocatable :: Rion_last(:,:), Force_last(:,:)
  end type s_md

  type s_ofile
     integer :: fh_rt, fh_rt_energy
     character(256) :: file_rt_data, file_rt_energy_data
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
    allocate(field%f(rg%is(1):rg%ie(1),rg%is(2):rg%ie(2),rg%is(3):rg%ie(3)))
    field%f = 0d0
  end subroutine allocate_scalar

  subroutine allocate_vector(rg,field)
    implicit none
    type(s_rgrid),intent(in) :: rg
    type(s_vector)           :: field
    allocate(field%v(3,rg%is(1):rg%ie(1),rg%is(2):rg%ie(2),rg%is(3):rg%ie(3)))
    field%v = 0d0
  end subroutine allocate_vector

  subroutine allocate_dmatrix(nspin,mg,info,dmat)
    implicit none
    integer                 ,intent(in) :: nspin
    type(s_rgrid)           ,intent(in) :: mg
    type(s_orbital_parallel),intent(in) :: info
    type(s_dmatrix)                     :: dmat
    allocate(dmat%zrho_mat(mg%Nd,mg%ndir,mg%is(1)-mg%Nd:mg%ie(1),mg%is(2)-mg%Nd:mg%ie(2),mg%is(3)-mg%Nd:mg%ie(3), &
    & nspin,info%im_s:info%im_e))
    dmat%zrho_mat = 0d0
  end subroutine allocate_dmatrix

  subroutine allocate_orbital_real(nspin,mg,info,psi)
    implicit none
    integer                 ,intent(in) :: nspin
    type(s_rgrid)           ,intent(in) :: mg
    type(s_orbital_parallel),intent(in) :: info
    type(s_orbital)                     :: psi
    allocate(psi%rwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),  &
                     nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
    psi%rwf = 0d0
  end subroutine allocate_orbital_real

  subroutine allocate_orbital_complex(nspin,mg,info,psi)
    implicit none
    integer                 ,intent(in) :: nspin
    type(s_rgrid)           ,intent(in) :: mg
    type(s_orbital_parallel),intent(in) :: info
    type(s_orbital)                     :: psi
    allocate(psi%zwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),  &
                     nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
    psi%zwf = 0d0
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
    DEAL(info%io_tbl)
    DEAL(info%jo_tbl)
    DEAL(info%occ)
    DEAL(info%irank_jo)
  end subroutine deallocate_orbital_parallel

  subroutine deallocate_orbital(psi)
    type(s_orbital) :: psi
    DEAL(psi%rwf)
    DEAL(psi%zwf)
    DEAL(psi%ztmp)
  end subroutine deallocate_orbital

  subroutine deallocate_stencil(stencil)
    type(s_stencil) :: stencil
    DEAL(stencil%vec_kAc)
    DEAL(stencil%isign)
  end subroutine deallocate_stencil

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
