!
!  Copyright 2017 SALMON developers
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

  type s_system
    integer :: iperiodic              ! iperiodic==0 --> isolated system, iperiodic==3 --> 3D periodic system
    integer :: ngrid,nspin,no,nk,nion ! # of r-grid points, spin indices, orbitals, k points, and ions
    real(8) :: Hvol,Hgs(3),al(3,3),det_al,brl(3,3)
    real(8),allocatable :: wtk(:) ! (1:nk), weight of k points
    real(8),allocatable :: Rion(:,:) ! (1:3,1:nion), atom position
    real(8),allocatable :: rocc(:,:,:) ! (1:no,1:nk,1:nspin), occupation rate
  end type s_system

  type s_energy
    real(8),allocatable :: esp(:,:,:) ! (1:no,1:nk,1:nspin), single-particle energy
    real(8) :: E_tot,E_kin,E_h,E_xc,E_ion_ion,E_ion_loc,E_ion_nloc
  end type s_energy

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
  end type s_pcomm_cache

  type s_sendrecv_grid
    ! Number of orbitals (4-th dimension of grid)
    integer :: nb
    ! Communicator
    integer :: icomm, myrank
    ! Neightboring MPI id (1:x,2:y,3:z, 1:upside,2:downside):
    integer :: neig(1:3, 1:2) 
    ! Communication requests (1:x,2:y,3:z, 1:upside,2:downside, 1:send,2:recv):
    integer :: ireq(1:3, 1:2, 1:2)
    ! PComm cache (1:x,2:y,3:z, 1:upside,2:downside, 1:src/2:dst)
    type(s_pcomm_cache) :: cache(1:3, 1:2, 1:2)
    ! Range (dim=1:x,2:y,3:z, dir=1:upside,2:downside, 1:src/2:dst, axis=1...3)
    integer :: is_block(1:3, 1:2, 1:2, 1:3)
    integer :: ie_block(1:3, 1:2, 1:2, 1:3)
    logical :: pcomm_initialized
  end type s_sendrecv_grid

  type s_wf_info
    logical :: if_divide_rspace
    integer :: irank_r(6)
    integer :: icomm_r   ! communicator for r-space
    integer :: icomm_ko  ! communicator for k-space & orbital
    integer :: icomm_rko ! communicator for r-space, k-space & orbital
    integer :: im_s,im_e,numm ! im=im_s,...,im_e, numm=im_e-im_s+1
    integer :: ik_s,ik_e,numk ! ik=ik_s,...,ik_e, numk=ik_e-ik_s+1
    integer :: io_s,io_e,numo ! io=io_s,...,io_e, numo=io_e-io_s+1
    real(8),allocatable :: occ(:,:,:) ! occ(io_s:io_e,ik_s:ik_e,1:nspin) = rocc*wk, occupation numbers
    integer,allocatable :: io_tbl(:)  ! jo=io_tbl(io), io=io_s~io_e, jo=1~no
    integer,allocatable :: irank_io(:) ! MPI rank of the orbital index #io
  end type s_wf_info

  type s_wavefunction
    real(8)   ,allocatable :: rwf(:,:,:,:,:,:,:) ! rwf(x,y,z,ispin,io,ik,im)
    complex(8),allocatable :: zwf(:,:,:,:,:,:,:) ! zwf(x,y,z,ispin,io,ik,im)
    complex(8),allocatable :: wrk(:,:,:,:)
  end type s_wavefunction

  type s_stencil
    real(8) :: lap0,lapt(4,3),nabt(4,3) !????? (4,3) --> (Nd,3)
    real(8),allocatable :: kAc(:,:) ! kAc(Nk,3)

  ! for non-orthogonal lattice
    logical :: if_orthogonal
    integer,allocatable :: sign(:,:) ! sign(3,4:ndir) (for ndir=4~6) 
    real(8),allocatable :: coef_lap(:,:),coef_nab(:,:) !?????? --> lapt,nabt (future work)
    real(8) :: matrix_A(3,3),matrix_B(3,3),coef_F(6)

  ! Experimental implementation of srg
    type(s_sendrecv_grid) :: srg
  end type s_stencil

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
    complex(8),allocatable :: zproj(:,:,:) ! zproj(j,ilma,ik) = exp(-i(k+A/c)r)*uv ! j=1~Mps(ia), ilma=1~Nlma
  end type s_pp_grid

! rho%f, V_local(1:nspin)%f, tau%f, V_H%f, V_xc%f, current(3)%f?
  type s_scalar
    real(8),allocatable :: f(:,:,:) ! f(x,y,z)
  end type s_scalar

! current%v, vec_A%v
  type s_vector
    real(8),allocatable :: v(:,:,:,:) ! v(1:3,x,y,z)
  end type s_vector

  type s_force
    real(8),allocatable :: force(:,:) ! force(1:3,1:NI)
  end type s_force

  type s_dmatrix
    complex(8),allocatable :: rho(:,:,:,:,:,:,:) ! rho(ii,dir,x,y,z,ispin,im), ii=1~Nd, dir=1~6(xx,yy,zz,yz,zx,xy)
  end type s_dmatrix

  type s_pp_nlcc
    real(8), allocatable :: rho_nlcc(:,:,:)
    real(8), allocatable :: tau_nlcc(:,:,:)
  end type s_pp_nlcc

  type s_fourier_grid
    integer :: icomm_fourier
    integer :: ng,iG_s,iG_e,iGzero
    real(8),allocatable :: Gx(:),Gy(:),Gz(:)
    complex(8),allocatable :: rhoG_ion(:),rhoG_elec(:),dVG_ion(:,:)
  end type s_fourier_grid

!===================================================================================================================================

contains

# define DEAL(x) if(allocated(x)) deallocate(x)

  subroutine deallocate_system(system)
    type(s_system) :: system
    DEAL(system%rocc)
    DEAL(system%wtk)
    DEAL(system%Rion)
  end subroutine deallocate_system

  subroutine deallocate_energy(energy)
    type(s_energy) :: energy
    DEAL(energy%esp)
  end subroutine deallocate_energy

  subroutine deallocate_rgrid(rg)
    type(s_rgrid) :: rg
    DEAL(rg%idx)
    DEAL(rg%idy)
    DEAL(rg%idz)
  end subroutine deallocate_rgrid

  subroutine deallocate_wf_info(info)
    type(s_wf_info) :: info
    DEAL(info%io_tbl)
    DEAL(info%occ)
  end subroutine deallocate_wf_info

  subroutine deallocate_wavefunction(psi)
    type(s_wavefunction) :: psi
    DEAL(psi%rwf)
    DEAL(psi%zwf)
    DEAL(psi%wrk)
  end subroutine deallocate_wavefunction

  subroutine deallocate_stencil(stencil)
    type(s_stencil) :: stencil
    DEAL(stencil%kAc)
    DEAL(stencil%sign)
    DEAL(stencil%coef_lap)
    DEAL(stencil%coef_nab)
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
    DEAL(ppg%zproj)
  end subroutine deallocate_pp_grid

  subroutine deallocate_scalar(x)
    type(s_scalar) :: x
    DEAL(x%f)
  end subroutine deallocate_scalar

  subroutine deallocate_vector(x)
    type(s_vector) :: x
    DEAL(x%v)
  end subroutine deallocate_vector

  subroutine deallocate_force(x)
    type(s_force) :: x
    DEAL(x%force)
  end subroutine deallocate_force

  subroutine deallocate_dmatrix(dm)
    type(s_dmatrix) :: dm
    DEAL(dm%rho)
  end subroutine deallocate_dmatrix

  subroutine deallocate_fourier_grid(fg)
    type(s_fourier_grid) :: fg
    DEAL(fg%Gx)
    DEAL(fg%Gy)
    DEAL(fg%Gz)
    DEAL(fg%rhoG_ion)
    DEAL(fg%rhoG_elec)
    DEAL(fg%dVG_ion)
  end subroutine deallocate_fourier_grid

end module structures
