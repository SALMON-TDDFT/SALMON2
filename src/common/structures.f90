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

  type s_wavefunction
    integer :: i1_s,i1_e,num1 ! i1=i1_s,...,i1_e, num1=i1_e-i1_s+1
    integer :: ik_s,ik_e,numk ! ik=ik_s,...,ik_e
    integer :: io_s,io_e,numo ! io=io_s,...,io_e
    real(8)   ,allocatable :: rwf(:,:,:,:,:,:,:) ! rwf(x,y,z,ispin,io,ik,i1)
    complex(8),allocatable :: zwf(:,:,:,:,:,:,:) ! zwf(x,y,z,ispin,io,ik,i1)
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

! memo: structures, pp_grid, hpsi, GCEED/scf(total_energy, force)
! mn2007/SALMON (branch: develop-2.0.0)

end module structures
