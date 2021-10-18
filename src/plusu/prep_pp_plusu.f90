!
!  Copyright 2019-2020 SALMON developers
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
module prep_pp_plusU_sub

  use plusU_global, only: PLUS_U_ON, read_Hubbard_parameters

  implicit none

  private
  public :: calc_uv_plusU
  public :: PLUS_U_ON

  integer :: Nlma_ao
  integer,allocatable :: lma_tbl_ao(:,:)

contains

  subroutine set_nlma_ao( pp, ppg )
    use salmon_global,only : natom,kion
    use structures,only : s_pp_info,s_pp_grid
    implicit none 
    type(s_pp_info) :: pp
    type(s_pp_grid) :: ppg
    integer :: lma, lm, lm_max
    integer :: a,ik,m,l,ll,l0
    integer :: lma2,lma1,lma2_0,m1,m2,nproj

    lm_max=0
    lma=0
    do a=1,natom
      ik=kion(a)
      lm=0
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        do m=-ll,ll
          lm=lm+1
          lma=lma+1
        end do
      end do
      l0=l
      end do
      lm_max=max(lm_max,lm)
    end do

    Nlma_ao = lma

!    write(*,*) "Nlma_ao (ppg%nlma)=", Nlma_ao, " (",ppg%nlma,")"

    allocate( lma_tbl_ao(lm_max,natom) ); lma_tbl_ao=0
    allocate( ppg%ia_tbl_ao(Nlma_ao) ); ppg%ia_tbl_ao=0

    call set_lma_tbl( lma_tbl_ao, ppg%ia_tbl_ao, pp )

    allocate( ppg%phi_ao(ppg%nps_ao,Nlma_ao)    ); ppg%phi_ao=0.0d0
    allocate( ppg%dphi_ao(ppg%nps_ao,Nlma_ao,3) ); ppg%dphi_ao=0.0d0

    nproj=0
    lma1=0
    lma2_0=0
    do a=1,natom
      ik=kion(a)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        do m1=-ll,ll
          lma1=lma1+1
          lma2=lma2_0
        do m2=-ll,ll
          lma2=lma2+1
          nproj=nproj+1
          !write(*,'(1x,8i6)') nproj,a,ll,l,m1,m2,lma1,lma2
        end do
        end do
        lma2_0=lma2
      end do
      l0=l
      end do
    end do

    allocate( ppg%proj_pairs_ao(2,nproj) ); ppg%proj_pairs_ao=0
    allocate( ppg%proj_pairs_info_ao(5,nproj) ); ppg%proj_pairs_info_ao=0

    nproj=0
    lma1=0
    lma2_0=0
    do a=1,natom
      ik=kion(a)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        do m1=-ll,ll
          lma1=lma1+1
          lma2=lma2_0
        do m2=-ll,ll
          lma2=lma2+1
          nproj=nproj+1
          ppg%proj_pairs_ao(1,nproj)=lma1
          ppg%proj_pairs_ao(2,nproj)=lma2
          ppg%proj_pairs_info_ao(1,nproj)=a
          ppg%proj_pairs_info_ao(2,nproj)=ll
          ppg%proj_pairs_info_ao(3,nproj)=l-l0+1
          ppg%proj_pairs_info_ao(4,nproj)=m1
          ppg%proj_pairs_info_ao(5,nproj)=m2
        end do
        end do
        lma2_0=lma2
      end do
      l0=l
      end do
    end do

!    write(*,*) "----- set_nlma_ao(end)"

  end subroutine set_nlma_ao

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  subroutine set_lma_tbl( lma_tbl, ia_tbl, pp )
    use salmon_global,only : natom, kion
    use structures,only : s_pp_info
    implicit none
    integer,intent(out) :: lma_tbl(:,:)
    integer,intent(out) :: ia_tbl(:)
    type(s_pp_info) :: pp
    integer :: lm,lma
    integer :: a,ik,m,l,l0,ll
    lma_tbl=0
    ia_tbl=0
    lma=0
    do a=1,natom
      ik=kion(a)
      lm=0
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        do m=-ll,ll
          lm=lm+1
          lma=lma+1
          lma_tbl(lm,a)=lma
          ia_tbl(lma)=a
        end do
      end do
      l0=l
      end do
    end do !a
  end subroutine set_lma_tbl

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

  subroutine calc_uv_plusU(pp,ppg,property)
    use salmon_global, only : natom, kion
    use salmon_math,   only : ylm, dylm, spline
    use structures,    only : s_pp_info, s_pp_grid
    implicit none
    type(s_pp_info),intent(in) :: pp
    type(s_pp_grid),intent(inout) :: ppg
    character(17),intent(in) :: property
    real(8) :: save_upptbl_a(pp%nrmax,0:2*pp%lmax+1,natom)
    real(8) :: save_upptbl_b(pp%nrmax,0:2*pp%lmax+1,natom)
    real(8) :: save_upptbl_c(pp%nrmax,0:2*pp%lmax+1,natom)
    real(8) :: save_upptbl_d(pp%nrmax,0:2*pp%lmax+1,natom)
    integer :: a,ik,j,l,lm,m,ll,l0
    integer :: ilma,intr,ir
    real(8),allocatable :: xn(:),yn(:),an(:),bn(:),cn(:),dn(:)  
    real(8) :: upp(0:2*pp%lmax+1)
    real(8) :: r,x,y,z, xx
    
    ppg%nps_ao = ppg%nps
    ppg%jxyz_ao = ppg%jxyz
    ppg%rxyz_ao = ppg%rxyz
    ppg%mps_ao = ppg%mps

    call set_nlma_ao( pp, ppg )

    if ( property /= 'update_wo_realloc' ) then

      do a=1,natom
        ik=kion(a)
        allocate( xn(0:pp%nrps_ao(ik)-1),yn(0:pp%nrps_ao(ik)-1),an(0:pp%nrps_ao(ik)-2) &
                 ,bn(0:pp%nrps_ao(ik)-2),cn(0:pp%nrps_ao(ik)-2),dn(0:pp%nrps_ao(ik)-2) )
        xn(0:pp%nrps_ao(ik)-1) = pp%rad(1:pp%nrps_ao(ik),ik)
        l0=0
        do ll=0,pp%mlps(ik)
        do l=l0,l0+pp%nproj(ll,ik)-1
          yn(0:pp%nrps_ao(ik)-1) = pp%upptbl_ao(1:pp%nrps_ao(ik),l,ik)
          call spline(pp%nrps_ao(ik),xn,yn,an,bn,cn,dn)
          save_upptbl_a(1:pp%nrps_ao(ik)-1,l,a) = an(0:pp%nrps_ao(ik)-2)
          save_upptbl_b(1:pp%nrps_ao(ik)-1,l,a) = bn(0:pp%nrps_ao(ik)-2)
          save_upptbl_c(1:pp%nrps_ao(ik)-1,l,a) = cn(0:pp%nrps_ao(ik)-2)
          save_upptbl_d(1:pp%nrps_ao(ik)-1,l,a) = dn(0:pp%nrps_ao(ik)-2)
        end do
        l0=l
        end do
        deallocate(xn,yn,an,bn,cn,dn)
      end do !a

    end if

    do a=1,natom

      ik=kion(a)

      do j=1,ppg%mps_ao(a)

        x=ppg%rxyz_ao(1,j,a)
        y=ppg%rxyz_ao(2,j,a)
        z=ppg%rxyz_ao(3,j,a)
        r=sqrt(x*x+y*y+z*z)+1d-50
        do ir=1,pp%nrps_ao(ik)
          if ( pp%rad(ir,ik) > r ) exit
        end do
        intr=ir-1
        if (intr < 0 .or. intr >= pp%nrps_ao(ik) ) stop 'bad intr at prep_ps'
        xx = r - pp%rad(intr,ik)

        l0=0
        do ll=0,pp%mlps(ik)
        do l=l0,l0+pp%nproj(ll,ik)-1
          upp(l) = save_upptbl_a(intr,l,a)*xx**3 + save_upptbl_b(intr,l,a)*xx**2 &
                 + save_upptbl_c(intr,l,a)*xx    + save_upptbl_d(intr,l,a)
        end do
        l0=l
        end do
 
        lm=0
        l0=0
        do ll=0,pp%mlps(ik)
        do l=l0,l0+pp%nproj(ll,ik)-1
          do m=-ll,ll
            lm=lm+1
            ilma=lma_tbl_ao(lm,a)
            ppg%phi_ao(j,ilma) = upp(l)*ylm(x,y,z,ll,m)
          end do !m
        end do
        l0=l
        end do

      end do !j

    end do !a

!    rewind 10
!    do j=1,ppg%mps_ao(1)
!       x=ppg%rxyz(1,j,1)
!       y=ppg%rxyz(2,j,1)
!       z=ppg%rxyz(3,j,1)
!       r=sqrt(x*x+y*y+z*z)
!       write(10,*) r, ppg%phi_ao(j,1)
!    end do

!write(*,*) "---------------- calc_uv_plusU(end)"
!stop "stop@calc_uv_plusU"

  end subroutine calc_uv_plusU

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  subroutine bisection(xx,inode,iak,nr,rad_psl)
  use salmon_global,only : nelem
  implicit none
  integer,intent(out) :: inode
  integer,intent(in)  :: iak
  integer,intent(in)  :: nr
  real(8),intent(in)  :: rad_psl(nr,nelem)
  real(8),intent(in)  :: xx
  integer :: imin,imax
  
  imin=1
  imax=nr
  do while (imax-imin>1)
    inode=(imin+imax)/2
    if(xx>rad_psl(inode,iak))then
      imin=inode
    else
      imax=inode
    end if
  end do
  inode=imin

  end subroutine bisection

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  subroutine init_uvpsi_summation(ppg,icomm_r)
  use structures,    only: s_pp_grid
  use salmon_global, only: natom
  use communication, only: comm_get_groupinfo &
                                 ,comm_allgather &
                                 ,comm_create_group_byid
  implicit none
  type(s_pp_grid),intent(inout) :: ppg
  integer,intent(in) :: icomm_r

  integer :: ilma,ia
  integer :: irank_r,isize_r,n,i
  integer :: nlma
  logical,allocatable :: ireferred_atom_comm_r(:,:)
  integer,allocatable :: iranklist(:)
stop "xxxx"
#if 0
  call comm_get_groupinfo(icomm_r, irank_r, isize_r)

  nlma = ppg%Nlma

  allocate(iranklist(isize_r))
  allocate(ireferred_atom_comm_r(natom,isize_r))
  allocate(ppg%ireferred_atom(natom))
  allocate(ppg%icomm_atom(natom))
  allocate(ppg%irange_atom(2,natom))

  ppg%ireferred_atom = .false.
  do ilma=1,nlma
    ia = ppg%ia_tbl(ilma)
    ppg%ireferred_atom(ia) = ppg%ireferred_atom(ia) .or. (ppg%mps(ia) > 0)
  end do
  call comm_allgather(ppg%ireferred_atom, ireferred_atom_comm_r, icomm_r)

  ppg%irange_atom(1,:) = 1
  ppg%irange_atom(2,:) = 0
  do ia=1,natom
    ! forward search
    do ilma=1,nlma
      if (ppg%ia_tbl(ilma) == ia) then
        ppg%irange_atom(1,ia) = ilma
        exit
      end if
    end do

    ! backward search
    do ilma=nlma,1,-1
      if (ppg%ia_tbl(ilma) == ia) then
        ppg%irange_atom(2,ia) = ilma
        exit
      end if
    end do
  end do

  do ia=1,natom
    n = 0
    do i=1,isize_r
      if (ireferred_atom_comm_r(ia,i)) then
        n = n + 1
        iranklist(n) = i - 1
      end if
    end do

    ppg%icomm_atom(ia) = comm_create_group_byid(icomm_r, iranklist(1:n))
  end do
#endif
  end subroutine init_uvpsi_summation

  subroutine finalize_uvpsi_summation(ppg)
  use structures,    only: s_pp_grid
  use communication, only: comm_free_group
  implicit none
  type(s_pp_grid),intent(inout) :: ppg
  integer :: ia

  if (allocated(ppg%irange_atom))    deallocate(ppg%irange_atom)
  if (allocated(ppg%ireferred_atom)) deallocate(ppg%ireferred_atom)
  if (allocated(ppg%icomm_atom)) then
    do ia=1,size(ppg%icomm_atom)
      call comm_free_group(ppg%icomm_atom(ia))
    end do
    deallocate(ppg%icomm_atom)
  end if
  end subroutine finalize_uvpsi_summation

end module prep_pp_plusU_sub
