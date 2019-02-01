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
module structures
  implicit none

  type s_system
    integer :: norb,nspin,nk
    real(8) :: Hvol,Hgs(3)
    real(8),allocatable :: occ(:,:),wk(:),esp(:,:) ! occ(NB,NK),wk(NK),esp(NB,NK)
    real(8),allocatable :: Rion(:,:) ! atom position
  end type s_system

  type s_rgrid
    integer,dimension(3) :: is,ie,num &             ! num=ie-is+1
                           ,is_overlap,ie_overlap & ! is_overlap=is-4, ie_overlap=ie+4
                           ,is_array,ie_array       ! allocate( array(is_array(1):ie_array(1), ...) )
    integer ,allocatable :: idx(:),idy(:),idz(:)    ! idx(is_overlap(1):ie_overlap(1))=is_array(1)~ie_array(1), ...
  end type s_rgrid

  type s_wf_info
    logical :: if_divide_rspace
    integer :: irank_overlap(6),icomm_overlap,icomm_pseudo
    integer :: im_s,im_e,numm ! im=im_s,...,im_e, numm=im_e-im_s+1
    integer :: ik_s,ik_e,numk ! ik=ik_s,...,ik_e, numk=ik_e-ik_s+1
    integer :: io_s,io_e,numo ! io=io_s,...,io_e, numo=io_e-io_s+1
  end type s_wf_info

  type s_wavefunction
    real(8)   ,allocatable :: rwf(:,:,:,:,:,:,:) ! rwf(x,y,z,ispin,io,ik,im)
    complex(8),allocatable :: zwf(:,:,:,:,:,:,:) ! zwf(x,y,z,ispin,io,ik,im)
  end type s_wavefunction

  type s_stencil
    real(8) :: lap0,lapt(4,3),nabt(4,3)
    real(8),allocatable :: kAc(:,:) ! kAc(Nk,3)
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

contains

# define DEAL(x) if(allocated(x)) deallocate(x)

  subroutine deallocate_system(system)
    type(s_system) :: system
    DEAL(system%occ)
    DEAL(system%wk)
    DEAL(system%esp)
    DEAL(system%Rion)
  end subroutine deallocate_system

  subroutine deallocate_rgrid(rg)
    type(s_rgrid) :: rg
    DEAL(rg%idx)
    DEAL(rg%idy)
    DEAL(rg%idz)
  end subroutine deallocate_rgrid

  subroutine deallocate_wavefunction(psi)
    type(s_wavefunction) :: psi
    DEAL(psi%rwf)
    DEAL(psi%zwf)
  end subroutine deallocate_wavefunction

  subroutine deallocate_stencil(stencil)
    type(s_stencil) :: stencil
    DEAL(stencil%kAc)
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

end module structures
