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

  type s_rgrid
    integer,dimension(3) :: is,ie,num &             ! num=ie-is+1
                           ,is_overlap,ie_overlap & ! is_overlap=is-4, ie_overlap=ie+4
                           ,is_array,ie_array       ! allocate( array(is_array(1):ie_array(1), ...) )
    integer ,allocatable :: idx(:),idy(:),idz(:)    ! idx(is_overlap(1):ie_overlap(1))=is_array(1)~ie_array(1), ...
  end type s_rgrid

  type s_wavefunction
    integer :: i1_s,i1_e,n1,i2_s,i2_e,n2,i3_s,i3_e,n3 ! i1=i1_s,...,i1_e, n1=i1_e-i1_s+1
    real(8)   ,allocatable :: rwf(:,:,:,:,:,:,:) ! rwf(x,y,z,ispin,i3,i2,i1)
    complex(8),allocatable :: zwf(:,:,:,:,:,:,:) ! zwf(x,y,z,ispin,i3,i2,i1)
!    real(8)   ,allocatable :: rwf(:,:,:,:)   ! rwf(x,y,z,iorb)
!    complex(8),allocatable :: zwf(:,:,:,:,:) ! zwf(x,y,z,isigma,iorb)
  end type s_wavefunction

  type s_stencil
    real(8) :: lap0,lapt(4,3),nabt(4,3)
    real(8),allocatable :: kAc(:,:) ! kAc(Nk,3)
  end type s_stencil

  type struct_pseudopotential ! --> pp_grid
     integer :: NI,Nps,Nlma
     integer,allocatable :: index_atom(:),Mps(:),Jxyz(:,:,:)
     real(8),allocatable :: uV(:,:),uVu(:)
! exp_ikr=ekr --> zproj
     complex(8),allocatable :: zproj(:,:,:) ! zproj(j,ilma,ik) = exp(-i(k+A/c)r)*uV ! j=1~Mps(ia), ilma=1~Nlma
  end type struct_pseudopotential

! rho%f, V_local(1:nspin)%f, tau%f, V_H%f, V_xc%f, current(3)%f?
  type s_scalar
    real(8),allocatable :: f(:,:,:) ! f(x,y,z)
  end type s_scalar

! current%v, vec_A%v
  type s_vector
    real(8),allocatable :: v(:,:,:,:) ! v(1:3,x,y,z)
  end type s_vector

  type s_system
    integer :: norb,nspin,nk
    real(8) :: Hvol,Hgs(3)
    real(8),allocatable :: occ(:,:),wk(:),esp(:,:) ! occ(NB,NK),wk(NK),esp(NB,NK)
    real(8),allocatable :: Rion(:,:) ! atom position
  end type s_system

! real(8),allocatable :: force(:,:),Floc(:,:),Fnl(:,:),Fion(:,:),FionAc(:,:)
! Ftot%force, Floc%force, Fnl%force, Fion%force, FionAc%force ?
  type s_force
    real(8),allocatable :: force(:,:) ! force(1:3,1:NI)
  end type s_force

! memo: structures, pp_grid, hpsi, GCEED/scf(total_energy, force)
! mn2007/SALMON (branch: develop-2.0.0)

end module structures
