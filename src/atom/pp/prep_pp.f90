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
module prep_pp_sub
  implicit none

contains

subroutine init_ps(lg,mg,ng,system,info,info_field,fg,poisson,pp,ppg,sVpsl)
  use structures
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use salmon_global, only: iperiodic,yn_ffte
  implicit none
  type(s_rgrid)           ,intent(in) :: lg,mg,ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_field_parallel)  ,intent(in) :: info_field
  type(s_reciprocal_grid)             :: fg
  type(s_poisson)                     :: poisson
  type(s_pp_info)         ,intent(in) :: pp
  type(s_pp_grid)                     :: ppg
  type(s_scalar)                      :: sVpsl
  !
  integer :: ix,iy,iz,i,nl
  integer :: mmx(mg%num(1)*mg%num(2)*mg%num(3))
  integer :: mmy(mg%num(1)*mg%num(2)*mg%num(3))
  integer :: mmz(mg%num(1)*mg%num(2)*mg%num(3))
  integer :: lx(lg%num(1)*lg%num(2)*lg%num(3))
  integer :: ly(lg%num(1)*lg%num(2)*lg%num(3))
  integer :: lz(lg%num(1)*lg%num(2)*lg%num(3))
  real(8) :: alx,aly,alz
  real(8) :: hx,hy,hz
  character(17) :: property
  real(8),allocatable :: save_udVtbl_a(:,:,:)
  real(8),allocatable :: save_udVtbl_b(:,:,:)
  real(8),allocatable :: save_udVtbl_c(:,:,:)
  real(8),allocatable :: save_udVtbl_d(:,:,:)
  logical :: flag_use_grad_wf_on_force

  if(comm_is_root(nproc_id_global))then
    write(*,*) ''
    write(*,*) '============init_ps=============='
  endif

  hx = system%Hgs(1)
  hy = system%Hgs(2)
  hz = system%Hgs(3)
  alx = system%Hgs(1)*dble(lg%num(1))
  aly = system%Hgs(2)*dble(lg%num(2))
  alz = system%Hgs(3)*dble(lg%num(3))
  nl = lg%num(1)*lg%num(2)*lg%num(3)

  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    i=(iz-mg%is(3))*mg%num(1)*mg%num(2)+(iy-mg%is(2))*mg%num(1)+ix-mg%is(1)+1
    mmx(i)=ix
    mmy(i)=iy
    mmz(i)=iz
  end do
  end do
  end do

  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    i=(iz-lg%is(3))*lg%num(1)*lg%num(2)+(iy-lg%is(2))*lg%num(1)+ix-lg%is(1)+1
    lx(i)=ix
    ly(i)=iy
    lz(i)=iz
  end do
  end do
  end do

  call calc_nps(pp,ppg,alx,aly,alz,lx,ly,lz,lg%num(1)*lg%num(2)*lg%num(3),   &
                                   mmx,mmy,mmz,mg%num(1)*mg%num(2)*mg%num(3),   &
                                   hx,hy,hz,system%primitive_a,system%rmatrix_A)

  call init_jxyz(ppg)

  call calc_jxyz(pp,ppg,alx,aly,alz,lx,ly,lz,lg%num(1)*lg%num(2)*lg%num(3),   &
                                    mmx,mmy,mmz,mg%num(1)*mg%num(2)*mg%num(3),   &
                                    hx,hy,hz,system%primitive_a,system%rmatrix_A)

  call set_nlma(pp,ppg)
  call init_lma_tbl(pp,ppg)
  call init_uv(pp,ppg)
  call set_lma_tbl(pp,ppg)

  allocate( save_udVtbl_a(pp%nrmax,0:2*pp%lmax+1,system%nion) )
  allocate( save_udVtbl_b(pp%nrmax,0:2*pp%lmax+1,system%nion) )
  allocate( save_udVtbl_c(pp%nrmax,0:2*pp%lmax+1,system%nion) )
  allocate( save_udVtbl_d(pp%nrmax,0:2*pp%lmax+1,system%nion) )

  property='initial'
  flag_use_grad_wf_on_force=.false.
  call calc_uv(pp,ppg,save_udvtbl_a,save_udvtbl_b,save_udvtbl_c,save_udvtbl_d, &
               lx,ly,lz,nl,hx,hy,hz,  &
               flag_use_grad_wf_on_force,property,system%Hvol)

  select case(iperiodic)
  case(0)
    call calc_Vpsl_isolated(mg,lg,system,pp,sVpsl,ppg)
  case(3)
    select case(yn_ffte)
    case('n')
      call calc_vpsl_periodic(mg,system,pp,fg,sVpsl,ppg)
    case('y')
      if(.not.allocated(sVpsl%f)) call allocate_scalar(mg,sVpsl)
      call calcVpsl_periodic_FFTE(lg,ng,system,fg,info_field,poisson,pp,ppg,sVpsl)
    end select
  end select

  call init_uvpsi_summation(ppg,info%icomm_r)

end subroutine init_ps

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

SUBROUTINE dealloc_init_ps(ppg)
  use structures, only: s_pp_grid
  implicit none
  type(s_pp_grid) :: ppg

  deallocate(ppg%jxyz, ppg%jxx, ppg%jyy, ppg%jzz, ppg%rxyz)
  deallocate(ppg%lma_tbl, ppg%ia_tbl)
  deallocate(ppg%rinv_uvu,ppg%uv,ppg%duv)

  if(allocated(ppg%Vpsl_atom)) deallocate(ppg%Vpsl_atom)
  if(allocated(ppg%zekr_uV)) deallocate(ppg%zekr_uV)

  call finalize_uvpsi_summation(ppg)
END SUBROUTINE dealloc_init_ps

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

! for ARTED
! future work: remove this subroutine
subroutine calc_vloc(pp,dvloc_g,gx,gy,gz,ng,ng_s,ng_e,ngzero)
  use salmon_global,only : nelem
  use math_constants,only : pi
  use structures,only : s_pp_info
  implicit none
  type(s_pp_info) :: pp
  integer,intent(in) :: ng,ng_s,ng_e
  integer,intent(in) :: ngzero
  real(8),intent(in) :: gx(ng),gy(ng),gz(ng)
  complex(8),intent(out) :: dvloc_g(ng_s:ng_e,nelem)
  integer :: i,ik,n
  real(8) :: dr
  real(8) :: g2sq
  real(8) :: vloc_av
  real(8) :: s,r

!$omp parallel
!$omp do private(ik,n,g2sq,s,r,dr,i,vloc_av) collapse(2)
  do ik=1,nelem
    do n=ng_s,ng_e
      g2sq=sqrt(gx(n)**2+gy(n)**2+gz(n)**2)
      s=0.d0
      if (n == ngzero) then
        do i=2,pp%nrloc(ik)
          r=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
          dr=pp%rad(i,ik)-pp%rad(i-1,ik)
          vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
          s=s+4d0*pi*(r**2*vloc_av+r*pp%zps(ik))*dr
        enddo
      else
        do i=2,pp%nrloc(ik)
          r=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
          dr=pp%rad(i,ik)-pp%rad(i-1,ik)
          vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
          s=s+4d0*pi*sin(g2sq*r)/g2sq*(r*vloc_av+pp%zps(ik))*dr !Vloc - coulomb
        enddo
      endif
      dvloc_g(n,ik)=s
    enddo
  enddo
!$omp end do
!$omp end parallel

end subroutine calc_vloc

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

! for ARTED
! future work: remove this subroutine
subroutine calc_vpsl(pp,rhoion_g,vpsl_ia,vpsl,dvloc_g,  &
                     ngzero,gx,gy,gz,ng,ng_s,ng_e,nl,alxyz,lx,ly,lz,hx,hy,hz,matrix_A0)
  use salmon_global,only : natom, nelem, kion, rion
  use salmon_parallel,only : nproc_group_tdks
  use salmon_communication, only: comm_summation
  use math_constants,only : pi
  use structures,only : s_pp_info
  implicit none
  type(s_pp_info) :: pp
  integer,intent(in) :: ngzero
  integer,intent(in) :: ng,ng_s,ng_e
  complex(8),intent(in) :: dvloc_g(ng_s:ng_e,nelem)
  real(8),intent(in) :: gx(ng),gy(ng),gz(ng)
  integer,intent(in) :: nl
  real(8),intent(in) :: alxyz
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  real(8),intent(in) :: hx,hy,hz
  complex(8),intent(out) :: rhoion_g(ng_s:ng_e)
  real(8),intent(out) :: vpsl_ia(nl,natom)
  real(8),intent(out) :: vpsl(nl)
  real(8),intent(in),optional :: matrix_A0(3,3)
  integer :: a,i,n,ik
  real(8) :: gd,r(3),matrix_A(3,3),u,v,w
  real(8) :: g2
  complex(8) :: vion_g_ia(ng_s:ng_e,natom),tmp_exp !, Vion_G(NG_s:NG_e)
  real(8) :: vpsl_ia_l(nl,natom)
  real(8) :: gr
  complex(8),parameter :: zi=(0.d0,1.d0)

  matrix_A = 0d0
  matrix_A(1,1) = 1d0
  matrix_A(2,2) = 1d0
  matrix_A(3,3) = 1d0
  if(present(matrix_A0)) matrix_A = matrix_A0

  !(Local pseudopotential: Vlocal in G-space(=Vion_G))
  vion_g_ia=0.d0
 !Vion_G   =0.d0
  rhoion_g =0.d0
!$omp parallel private(a,ik)
  do a=1,natom
    ik=kion(a)
!$omp do private(n,g2,gd,tmp_exp)
    do n=ng_s,ng_e
      gd=gx(n)*rion(1,a)+gy(n)*rion(2,a)+gz(n)*rion(3,a)
      tmp_exp = exp(-zi*gd)/alxyz
     !Vion_G(n)     = Vion_G(n)      + dvloc_g(n,ik)*tmp_exp
      vion_g_ia(n,a)= vion_g_ia(n,a) + dvloc_g(n,ik)*tmp_exp
      rhoion_g(n)   = rhoion_g(n) + pp%zps(ik)*tmp_exp
      if(n == ngzero) cycle
      !(add coulomb as dvloc_g is given by Vloc - coulomb)
      g2=gx(n)**2+gy(n)**2+gz(n)**2
     !Vion_G(n)     = Vion_G(n)      -4d0*pi/g2*pp%zps(ik)*tmp_exp
      vion_g_ia(n,a)= vion_g_ia(n,a) -4d0*pi/g2*pp%zps(ik)*tmp_exp
    enddo
!$omp end do
  enddo
!$omp end parallel

  !(Local pseudopotential: Vlocal(=Vpsl) in real-space)
  vpsl_ia_l=0.d0
 !Vpsl_l   =0.d0
!$omp parallel private(n)
  do n=ng_s,ng_e
!$omp do private(i,gr,a,tmp_exp,u,v,w,r)
    do i=1,NL
      u = lx(i)*hx
      v = ly(i)*hy
      w = lz(i)*hz
      r(1) = u*matrix_A(1,1) + v*matrix_A(1,2) + w*matrix_A(1,3)
      r(2) = u*matrix_A(2,1) + v*matrix_A(2,2) + w*matrix_A(2,3)
      r(3) = u*matrix_A(3,1) + v*matrix_A(3,2) + w*matrix_A(3,3)
      gr = gx(n)*r(1) + gy(n)*r(2) + gz(n)*r(3)
!      gr = gx(n)*lx(i)*hx+gy(n)*ly(i)*hy+gz(n)*lz(i)*hz
      tmp_exp = exp(zi*gr)
      !vpsl_l(i) = vpsl_l(i) + vion_G(n)*tmp_exp
      do a=1,natom
        vpsl_ia_l(i,a)= vpsl_ia_l(i,a)+ vion_g_ia(n,a)*tmp_exp
      enddo
    enddo
!$omp end do
  enddo
!$omp end parallel

 !call comm_summation(vpsl_l,vpsl,nl,nproc_group_tdks)
  call comm_summation(vpsl_ia_l,vpsl_ia,nl*natom,nproc_group_tdks)

!$omp parallel
!$omp do private(i)
  do i=1,nl
    vpsl(i) = sum(vpsl_ia(i,:))
  enddo
!$omp end do
!$omp end parallel

end subroutine calc_vpsl

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

SUBROUTINE calc_Vpsl_isolated(mg,lg,system,pp,vpsl,ppg)
  use structures
  use salmon_global,only : natom, kion
  use salmon_parallel, only: nproc_id_global
  implicit none
  type(s_rgrid)     ,intent(in) :: mg,lg
  type(s_dft_system),intent(in) :: system
  type(s_pp_info),intent(in) :: pp
  type(s_scalar)             :: vpsl
  type(s_pp_grid)            :: ppg
  !
  integer :: ix,iy,iz,ak
  integer :: j,a,intr
  real(8) :: ratio1,ratio2,r

  Vpsl%f=0.d0

  if(.not.allocated(vpsl%f)) call allocate_scalar(mg,Vpsl)
  allocate(ppg%Vpsl_atom(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),natom))
  do a=1,natom
    ak=Kion(a)
    do j=1,3
      if(abs(system%Rion(j,a))<lg%num(j)*system%Hgs(j))then
        continue
      else
        write(*,*) "Rion error",nproc_id_global,a,j,system%Rion(j,a)
      end if
    end do
    do ix=mg%is(1),mg%ie(1)
    do iy=mg%is(2),mg%ie(2)
    do iz=mg%is(3),mg%ie(3)
      r=sqrt( (lg%coordinate(ix,1)-system%Rion(1,a))**2      &
             +(lg%coordinate(iy,2)-system%Rion(2,a))**2      &
             +(lg%coordinate(iz,3)-system%Rion(3,a))**2 )+1.d-50
      call bisection(r,intr,ak,pp%nrmax,pp%rad)
      ratio1=(r-pp%rad(intr,ak))/(pp%rad(intr+1,ak)-pp%rad(intr,ak)) ; ratio2=1.d0-ratio1
      if(intr>0.and.intr<=pp%nrmax)then
        continue
      else
        write(*,*) "intr error",nproc_id_global,intr,r
      end if

      Vpsl%f(ix,iy,iz)=Vpsl%f(ix,iy,iz)      &
                  +ratio1*pp%vpp_f(intr,pp%Lref(ak),ak)      &
                  +ratio2*pp%vpp_f(intr-1,pp%Lref(ak),ak)  !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate

      ppg%Vpsl_atom(ix,iy,iz,a) = ratio1*pp%vpp_f(intr,pp%Lref(ak),ak) + ratio2*pp%vpp_f(intr-1,pp%Lref(ak),ak)
    end do
    end do
    end do
  end do

  allocate(ppg%zekr_uV(ppg%nps,ppg%nlma,1))
  ppg%zekr_uV(:,:,1) = cmplx(ppg%uV)

  return
END SUBROUTINE calc_Vpsl_isolated

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

subroutine calc_vpsl_periodic(mg,system,pp,fg,vpsl,ppg)
  use salmon_global,only : natom, nelem, kion
  use salmon_communication, only: comm_summation
  use math_constants,only : pi,zi
  use structures
  implicit none
  type(s_rgrid)     ,intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_pp_info),intent(in) :: pp
  type(s_reciprocal_grid)    :: fg
  type(s_scalar)             :: vpsl
  type(s_pp_grid)            :: ppg
  !
  integer :: a,i,n,ik,ix,iy,iz
  real(8) :: gd,r(3),matrix_A(3,3),u,v,w
  real(8) :: g2,gr,s,g2sq,r1,dr,vloc_av
  complex(8) :: vion_g_ia(fg%ng,natom),tmp_exp

  matrix_A = system%rmatrix_A
  if(.not.allocated(vpsl%f)) call allocate_scalar(mg,Vpsl)
  allocate(ppg%Vpsl_atom(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),natom))

!$omp parallel
!$omp do private(ik,n,g2sq,s,r1,dr,i,vloc_av) collapse(2)
  do ik=1,nelem
!    do n=fg%ig_s,fg%ig_e
      do n=1,fg%ng
      g2sq=sqrt(fg%gx(n)**2+fg%gy(n)**2+fg%gz(n)**2)
      s=0.d0
      if (n == fg%iGzero) then
        do i=2,pp%nrloc(ik)
          r1=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
          dr=pp%rad(i,ik)-pp%rad(i-1,ik)
          vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
          s=s+4d0*pi*(r1**2*vloc_av+r1*pp%zps(ik))*dr
        end do
      else
        do i=2,pp%nrloc(ik)
          r1=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
          dr=pp%rad(i,ik)-pp%rad(i-1,ik)
          vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
          s=s+4d0*pi*sin(g2sq*r1)/g2sq*(r1*vloc_av+pp%zps(ik))*dr !Vloc - coulomb
        end do
      end if
      fg%zdVG_ion(n,ik) = s
    end do
  end do
!$omp end do
!$omp end parallel

  !(Local pseudopotential: Vlocal in G-space(=Vion_G))
  vion_g_ia=0.d0
  fg%zrhoG_ion =0.d0
!$omp parallel private(a,ik)
  do a=1,natom
    ik=kion(a)
!$omp do private(n,g2,gd,tmp_exp)
    do n=1,fg%ng
!    do n=fg%ig_s,fg%ig_e
      gd = fg%gx(n)*system%rion(1,a) + fg%gy(n)*system%rion(2,a) + fg%gz(n)*system%rion(3,a)
      tmp_exp = exp(-zi*gd)/system%det_A
      vion_g_ia(n,a)  = vion_g_ia(n,a) + fg%zdVG_ion(n,ik)*tmp_exp
      fg%zrhoG_ion(n) = fg%zrhoG_ion(n) + pp%zps(ik)*tmp_exp
      if(n == fg%iGzero) cycle
      !(add coulomb as dvloc_g is given by Vloc - coulomb)
      g2 = fg%gx(n)**2+fg%gy(n)**2+fg%gz(n)**2
      vion_g_ia(n,a) = vion_g_ia(n,a) -4d0*pi/g2*pp%zps(ik)*tmp_exp
    end do
!$omp end do
  end do
!$omp end parallel

!(Local pseudopotential: Vlocal(=Vpsl) in real-space)
  ppg%Vpsl_atom = 0d0
!$omp parallel private(n)
  do n=1,fg%ng
!$omp do private(ix,iy,iz,gr,a,tmp_exp,u,v,w,r)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      u = dble(ix-1)*system%hgs(1)
      v = dble(iy-1)*system%hgs(2)
      w = dble(iz-1)*system%hgs(3)
      r(1) = u*matrix_A(1,1) + v*matrix_A(1,2) + w*matrix_A(1,3)
      r(2) = u*matrix_A(2,1) + v*matrix_A(2,2) + w*matrix_A(2,3)
      r(3) = u*matrix_A(3,1) + v*matrix_A(3,2) + w*matrix_A(3,3)
      gr = fg%gx(n)*r(1) + fg%gy(n)*r(2) + fg%gz(n)*r(3)
      tmp_exp = exp(zi*gr)
      do a=1,natom
        ppg%Vpsl_atom(ix,iy,iz,a) = ppg%Vpsl_atom(ix,iy,iz,a) + vion_g_ia(n,a)*tmp_exp
      end do
    end do
    end do
    end do
!$omp end do
  end do
!$omp end parallel

!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    Vpsl%f(ix,iy,iz) = sum(ppg%Vpsl_atom(ix,iy,iz,:))
  end do
  end do
  end do
!$omp end do
!$omp end parallel

end subroutine calc_vpsl_periodic

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_mps(ppg)
  use salmon_global,only : natom
  use structures,only : s_pp_grid
  implicit none 
  type(s_pp_grid) :: ppg

  allocate(ppg%mps(natom))

end subroutine init_mps
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_jxyz(ppg)
  use salmon_global,only : natom
  use structures,only : s_pp_grid
  implicit none 
  type(s_pp_grid) :: ppg

  allocate(ppg%jxyz(3,ppg%nps,natom))
  allocate(ppg%jxx( ppg%nps,natom))
  allocate(ppg%jyy( ppg%nps,natom))
  allocate(ppg%jzz( ppg%nps,natom))
  allocate(ppg%rxyz(3,ppg%nps,natom))

end subroutine init_jxyz
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine finalize_jxyz(ppg)
  use structures,only : s_pp_grid
  implicit none 
  type(s_pp_grid) :: ppg

  deallocate(ppg%jxyz)
  deallocate(ppg%jxx,ppg%jyy,ppg%jzz)
  deallocate(ppg%rxyz)

end subroutine finalize_jxyz
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

subroutine calc_nps(pp,ppg,alx,aly,alz,lx,ly,lz,nl,mx,my,mz,ml,hx,hy,hz,al0,matrix_A0)
  use salmon_global,only : natom,kion,rion,iperiodic,yn_domain_parallel
  use structures,only : s_pp_info,s_pp_grid
  implicit none
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  real(8),intent(in) :: alx,aly,alz
  integer,intent(in) :: nl,ml
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  integer,intent(in) :: mx(ml),my(ml),mz(ml)
  real(8),intent(in) :: hx,hy,hz
  real(8),intent(in),optional :: al0(3,3),matrix_A0(3,3)
  integer :: a,i,ik,ix,iy,iz,j
  integer :: nc
  real(8) :: tmpx,tmpy,tmpz
  real(8) :: x,y,z,r,u,v,w
  real(8) :: rshift(3),matrix_a(3,3),rr(3),al(3,3)
  integer :: mps_tmp(natom)

  matrix_a = 0d0
  matrix_a(1,1) = 1d0
  matrix_a(2,2) = 1d0
  matrix_a(3,3) = 1d0
  if(present(matrix_A0)) matrix_a = matrix_A0

  al = 0d0
  al(1,1) = alx
  al(2,2) = aly
  al(3,3) = alz
  if(present(al0)) al = al0

  if(iperiodic==0)then
    nc=0
  else if(iperiodic==3)then
    nc=2
  end if

  if(iperiodic==0)then
    if(mod(lx(nl)-lx(1)+1,2)==1)then
      rshift(1)=0.d0
    else
      rshift(1)=-0.5d0*Hx
    end if
    if(mod(ly(nl)-ly(1)+1,2)==1)then
      rshift(2)=0.d0
    else
      rshift(2)=-0.5d0*Hy
    end if
    if(mod(lz(nl)-lz(1)+1,2)==1)then
      rshift(3)=0.d0
    else
      rshift(3)=-0.5d0*Hz
    end if
  else if(iperiodic==3)then
    if(yn_domain_parallel=='y')then
      rshift(1)=-Hx
      rshift(2)=-Hy
      rshift(3)=-Hz
    else
      rshift(:)=0.d0
    end if
  end if


!$omp parallel
!$omp do private(a,ik,j,i,ix,iy,iz,tmpx,tmpy,tmpz,x,y,z,r,rr,u,v,w)
  do a=1,natom
    ik=kion(a)
    j=0
    do ix=-nc,nc
    do iy=-nc,nc
    do iz=-nc,nc
      rr(1) = ix*al(1,1) + iy*al(1,2) + iz*al(1,3)
      rr(2) = ix*al(2,1) + iy*al(2,2) + iz*al(2,3)
      rr(3) = ix*al(3,1) + iy*al(3,2) + iz*al(3,3)
      tmpx = rion(1,a)+ rr(1)
      tmpy = rion(2,a)+ rr(2)
      tmpz = rion(3,a)+ rr(3)
!      tmpx = rion(1,a)+ix*alx
!      tmpy = rion(2,a)+iy*aly
!      tmpz = rion(3,a)+iz*alz
      do i=1,ml
        u = mx(i)*hx + rshift(1)
        v = my(i)*hy + rshift(2)
        w = mz(i)*hz + rshift(3)
        rr(1) = u*matrix_a(1,1) + v*matrix_a(1,2) + w*matrix_a(1,3)
        rr(2) = u*matrix_a(2,1) + v*matrix_a(2,2) + w*matrix_a(2,3)
        rr(3) = u*matrix_a(3,1) + v*matrix_a(3,2) + w*matrix_a(3,3)
        x = rr(1) - tmpx
        y = rr(2) - tmpy
        z = rr(3) - tmpz
!        x=mx(i)*Hx+rshift(1)-tmpx
!        y=my(i)*Hy+rshift(2)-tmpy
!        z=mz(i)*Hz+rshift(3)-tmpz
        r=sqrt(x*x+y*y+z*z)
        if (r<pp%rps(ik)+1.d-12) j=j+1
      enddo
    enddo
    enddo
    enddo
    mps_tmp(a)=j
  end do
!$omp end do
!$omp end parallel

  ppg%nps=maxval(mps_tmp(:))

end subroutine calc_nps

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine calc_jxyz(pp,ppg,alx,aly,alz,lx,ly,lz,nl,mx,my,mz,ml,hx,hy,hz,al0,matrix_A0)
  use salmon_global,only : natom,kion,rion,iperiodic,yn_domain_parallel
  use structures,only : s_pp_info,s_pp_grid
  implicit none
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  real(8),intent(in) :: alx,aly,alz
  integer,intent(in) :: nl,ml
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  integer,intent(in) :: mx(ml),my(ml),mz(ml)
  real(8),intent(in) :: hx,hy,hz
  real(8),intent(in),optional :: al0(3,3),matrix_A0(3,3)
  integer :: a,i,ik,ix,iy,iz,j
  integer :: nc
  real(8) :: tmpx,tmpy,tmpz
  real(8) :: r,x,y,z,u,v,w
  real(8) :: rshift(3),matrix_a(3,3),rr(3),al(3,3)

  matrix_a = 0d0
  matrix_a(1,1) = 1d0
  matrix_a(2,2) = 1d0
  matrix_a(3,3) = 1d0
  if(present(matrix_A0)) matrix_a = matrix_A0

  al = 0d0
  al(1,1) = alx
  al(2,2) = aly
  al(3,3) = alz
  if(present(al0)) al = al0


  if(iperiodic==0)then
    nc=0
  else if(iperiodic==3)then
    nc=2
  end if

  if(iperiodic==0)then
    if(mod(lx(nl)-lx(1)+1,2)==1)then
      rshift(1)=0.d0
    else
      rshift(1)=-0.5d0*Hx
    end if
    if(mod(ly(nl)-ly(1)+1,2)==1)then
      rshift(2)=0.d0
    else
      rshift(2)=-0.5d0*Hy
    end if
    if(mod(lz(nl)-lz(1)+1,2)==1)then
      rshift(3)=0.d0
    else
      rshift(3)=-0.5d0*Hz
    end if
  else if(iperiodic==3)then 
    if(yn_domain_parallel=='y')then
      rshift(1)=-Hx
      rshift(2)=-Hy
      rshift(3)=-Hz
    else
      rshift(:)=0.d0
    end if
  end if

!$omp parallel
!$omp do private(a,ik,j,i,ix,iy,iz,tmpx,tmpy,tmpz,x,y,z,r,rr,u,v,w)
  do a=1,natom
    ik=kion(a)
    j=0
    do ix=-nc,nc
    do iy=-nc,nc
    do iz=-nc,nc
      rr(1) = ix*al(1,1) + iy*al(1,2) + iz*al(1,3)
      rr(2) = ix*al(2,1) + iy*al(2,2) + iz*al(2,3)
      rr(3) = ix*al(3,1) + iy*al(3,2) + iz*al(3,3)
      tmpx = rion(1,a)+ rr(1)
      tmpy = rion(2,a)+ rr(2)
      tmpz = rion(3,a)+ rr(3)
!      tmpx = rion(1,a)+ix*aLx
!      tmpy = rion(2,a)+iy*aLy
!      tmpz = rion(3,a)+iz*aLz
      do i=1,ml
        u = mx(i)*hx + rshift(1)
        v = my(i)*hy + rshift(2)
        w = mz(i)*hz + rshift(3)
        rr(1) = u*matrix_a(1,1) + v*matrix_a(1,2) + w*matrix_a(1,3)
        rr(2) = u*matrix_a(2,1) + v*matrix_a(2,2) + w*matrix_a(2,3)
        rr(3) = u*matrix_a(3,1) + v*matrix_a(3,2) + w*matrix_a(3,3)
        x = rr(1) - tmpx
        y = rr(2) - tmpy
        z = rr(3) - tmpz
!        x=mx(i)*Hx+rshift(1)-tmpx
!        y=my(i)*Hy+rshift(2)-tmpy
!        z=mz(i)*Hz+rshift(3)-tmpz
        r=sqrt(x*x+y*y+z*z)
        if (r<pp%rps(ik)) then
          j=j+1
          if (j<=ppg%nps) then
            ppg%jxyz(1,j,a)=mx(i)
            ppg%jxyz(2,j,a)=my(i)
            ppg%jxyz(3,j,a)=mz(i)
            ppg%jxx( j,a)=ix
            ppg%jyy( j,a)=iy
            ppg%jzz( j,a)=iz
            ppg%rxyz(1,j,a) = x
            ppg%rxyz(2,j,a) = y
            ppg%rxyz(3,j,a) = z
!            ppg%rxyz(1,j,a)=dble(mx(i))*hx+rshift(1)-(rion(1,a)+dble(ix)*alx)
!            ppg%rxyz(2,j,a)=dble(my(i))*hy+rshift(2)-(rion(2,a)+dble(iy)*aly)
!            ppg%rxyz(3,j,a)=dble(mz(i))*hz+rshift(3)-(rion(3,a)+dble(iz)*alz)
          endif
        endif
      enddo
    enddo
    enddo
    enddo
    ppg%mps(a)=j
  end do
!$omp end do
!$omp end parallel

end subroutine calc_jxyz
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_lma_tbl(pp,ppg)
  use salmon_global,only : natom
  use structures,only : s_pp_info,s_pp_grid
  implicit none 
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  integer :: n

  n=maxval(pp%nproj)*(pp%lmax+1)**2
  allocate(ppg%lma_tbl(n,natom)); ppg%lma_tbl=0

end subroutine init_lma_tbl
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine finalize_lma_tbl(ppg)
  use structures,only : s_pp_grid
  implicit none 
  type(s_pp_grid) :: ppg

  deallocate(ppg%lma_tbl)

end subroutine finalize_lma_tbl
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_uv(pp,ppg)
  use salmon_global,only : natom
  use structures,only : s_pp_info,s_pp_grid
  implicit none 
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  integer :: n

  n=maxval(pp%nproj)*(pp%lmax+1)**2
  allocate(ppg%ia_tbl(n*natom)); ppg%ia_tbl=0
  allocate(ppg%rinv_uvu(n*natom)); ppg%rinv_uvu=0.0d0
  allocate(ppg%uv(ppg%nps,ppg%nlma),ppg%duv(ppg%nps,ppg%nlma,3))

end subroutine init_uv
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine finalize_uv(ppg)
  use structures,only : s_pp_grid
  implicit none 
  type(s_pp_grid) :: ppg

  deallocate(ppg%ia_tbl,ppg%rinv_uvu)
  deallocate(ppg%uv,ppg%duv)

end subroutine finalize_uv
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine set_nlma(pp,ppg)
  use salmon_global,only : natom,kion
  use structures,only : s_pp_info,s_pp_grid
  implicit none 
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  integer :: lma
  integer :: a,ik,m,l,ll,l0

  lma=0
  do a=1,natom
    ik=kion(a)
    l0=0
    do ll=0,pp%mlps(ik)
    do l=l0,l0+pp%nproj(ll,ik)-1
      if(pp%inorm(l,ik)==0) cycle
      do m=-ll,ll
        lma=lma+1
      enddo
    enddo
    l0=l
    enddo
  enddo
  ppg%nlma=lma

end subroutine set_nlma
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine set_lma_tbl(pp,ppg)
  use salmon_global,only : natom,kion
  use structures,only : s_pp_info,s_pp_grid
  implicit none 
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  integer :: lm,lma
  integer :: a,ik,m,l,l0,ll

  lma=0
  do a=1,natom
    ik=kion(a)
    lm=0
    l0=0
    do ll=0,pp%mlps(ik)
    do l=l0,l0+pp%nproj(ll,ik)-1
      if(pp%inorm(l,ik)==0) cycle
      do m=-ll,ll
        lm=lm+1
        lma=lma+1
        ppg%lma_tbl(lm,a)=lma
        ppg%ia_tbl(lma)=a
      enddo
    enddo
    l0=l
    enddo
  enddo

end subroutine set_lma_tbl
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine calc_uv(pp,ppg,save_udvtbl_a,save_udvtbl_b,save_udvtbl_c,save_udvtbl_d, &
                   lx,ly,lz,nl,hx,hy,hz,  &
                   flag_use_grad_wf_on_force,property,hvol0)
  use salmon_global,only : natom,kion,iperiodic,yn_domain_parallel
  use math_constants,only : pi
  use structures,only : s_pp_info,s_pp_grid
  implicit none
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  integer,intent(in) :: nl
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  real(8),intent(in) :: hx,hy,hz
  logical,intent(in) :: flag_use_grad_wf_on_force
  character(17),intent(in) :: property
  real(8),intent(out) :: save_udvtbl_a(pp%nrmax,0:2*pp%lmax+1,natom)
  real(8),intent(out) :: save_udvtbl_b(pp%nrmax,0:2*pp%lmax+1,natom)
  real(8),intent(out) :: save_udvtbl_c(pp%nrmax,0:2*pp%lmax+1,natom)
  real(8),intent(out) :: save_udvtbl_d(pp%nrmax,0:2*pp%lmax+1,natom)
  real(8),intent(in),optional :: hvol0
  integer :: a,ik,j,l,lm,m,ll,l0
  integer :: ilma,intr,ir,lma
  real(8),allocatable :: xn(:),yn(:),an(:),bn(:),cn(:),dn(:)  
  real(8) :: dudvtbl_a(pp%nrmax,0:2*pp%lmax+1), dudvtbl_b(pp%nrmax,0:2*pp%lmax+1)
  real(8) :: dudvtbl_c(pp%nrmax,0:2*pp%lmax+1), dudvtbl_d(pp%nrmax,0:2*pp%lmax+1)
  real(8) :: uvr(0:2*pp%lmax+1),duvr(0:2*pp%lmax+1)
  real(8) :: r,x,y,z
  real(8) :: xx
  real(8) :: ylm,dylm
  real(8) :: rshift(3)
  real(8) :: hvol

  hvol=hx*hy*hz
  if(present(hvol0)) hvol = hvol0

  if(iperiodic==0)then
    if(mod(lx(nl)-lx(1)+1,2)==1)then
      rshift(1)=0.d0
    else
      rshift(1)=-0.5d0*Hx
    end if
    if(mod(ly(nl)-ly(1)+1,2)==1)then
      rshift(2)=0.d0
    else
      rshift(2)=-0.5d0*Hy
    end if
    if(mod(lz(nl)-lz(1)+1,2)==1)then
      rshift(3)=0.d0
    else
      rshift(3)=-0.5d0*Hz
    end if
  else if(iperiodic==3)then 
    if(yn_domain_parallel=='y')then
      rshift(1)=-Hx
      rshift(2)=-Hy
      rshift(3)=-Hz
    else
      rshift(:)=0.d0
    end if
  end if

  if(property /= 'update_wo_realloc') then

    do a=1,natom
      ik=kion(a)
      allocate(xn(0:pp%nrps(ik)-1),yn(0:pp%nrps(ik)-1),an(0:pp%nrps(ik)-2) &
              ,bn(0:pp%nrps(ik)-2),cn(0:pp%nrps(ik)-2),dn(0:pp%nrps(ik)-2))
  
      xn(0:pp%nrps(ik)-1) = pp%radnl(1:pp%nrps(ik),ik)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        yn(0:pp%nrps(ik)-1) = pp%udvtbl(1:pp%nrps(ik),l,ik)
        call spline(pp%nrps(ik),xn,yn,an,bn,cn,dn)
        save_udvtbl_a(1:pp%nrps(ik)-1,l,a) = an(0:pp%nrps(ik)-2)
        save_udvtbl_b(1:pp%nrps(ik)-1,l,a) = bn(0:pp%nrps(ik)-2)
        save_udvtbl_c(1:pp%nrps(ik)-1,l,a) = cn(0:pp%nrps(ik)-2)
        save_udvtbl_d(1:pp%nrps(ik)-1,l,a) = dn(0:pp%nrps(ik)-2)
      end do
      l0=l
      end do
      deallocate(xn,yn,an,bn,cn,dn)
    enddo
  end if
  
  do a=1,natom
    ik=kion(a)

    if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
      allocate(xn(0:pp%nrps(ik)-1),yn(0:pp%nrps(ik)-1),an(0:pp%nrps(ik)-2) &
              ,bn(0:pp%nrps(ik)-2),cn(0:pp%nrps(ik)-2),dn(0:pp%nrps(ik)-2))
      xn(0:pp%nrps(ik)-1) = pp%radnl(1:pp%nrps(ik),ik)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        yn(0:pp%nrps(ik)-1) = pp%dudvtbl(1:pp%nrps(ik),l,ik)
        call spline(pp%nrps(ik),xn,yn,an,bn,cn,dn)
        dudvtbl_a(1:pp%nrps(ik)-1,l) = an(0:pp%nrps(ik)-2)
        dudvtbl_b(1:pp%nrps(ik)-1,l) = bn(0:pp%nrps(ik)-2)
        dudvtbl_c(1:pp%nrps(ik)-1,l) = cn(0:pp%nrps(ik)-2)
        dudvtbl_d(1:pp%nrps(ik)-1,l) = dn(0:pp%nrps(ik)-2)
      end do
      l0=l
      end do
      deallocate(xn,yn,an,bn,cn,dn)
    endif
 
  !!$omp parallel
  !!$omp do private(j,x,y,z,r,ir,intr,xx,l,lm,m,uvr,duvr,ilma,l0,ll)
     do j=1,ppg%mps(a)
       x=ppg%rxyz(1,j,a)
       y=ppg%rxyz(2,j,a)
       z=ppg%rxyz(3,j,a)
       r=sqrt(x*x+y*y+z*z)+1d-50
       do ir=1,pp%nrps(ik)
         if(pp%radnl(ir,ik).gt.r) exit
       enddo
       intr=ir-1
       if (intr.lt.0.or.intr.ge.pp%nrps(ik))stop 'bad intr at prep_ps'
       xx = r - pp%radnl(intr,ik)
       l0=0
       do ll=0,pp%mlps(ik)
       do l=l0,l0+pp%nproj(ll,ik)-1
         uvr(l) = save_udvtbl_a(intr,l,a)*xx**3 + save_udvtbl_b(intr,l,a)*xx**2 &
                + save_udvtbl_c(intr,l,a)*xx    + save_udvtbl_d(intr,l,a)
         if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
           duvr(l)=dudvtbl_a(intr,l)*xx**3 +dudvtbl_b(intr,l)*xx**2 &
                  +dudvtbl_c(intr,l)*xx    +dudvtbl_d(intr,l)
         endif
       enddo
       l0=l
       enddo
 
       lm=0
       l0=0
       do ll=0,pp%mlps(ik)
       do l=l0,l0+pp%nproj(ll,ik)-1
         if(pp%inorm(l,ik)==0) cycle
         do m=-ll,ll
           lm=lm+1
           ilma=ppg%lma_tbl(lm,a)
           ppg%uv(j,ilma)=uvr(l)*ylm(x,y,z,ll,m)
           if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
             if(r>1d-6)then
               ppg%duv(j,ilma,1) = duvr(l)*(x/r)*ylm(x,y,z,ll,m)+uvr(l)*dylm(x,y,z,ll,m,1)
               ppg%duv(j,ilma,2) = duvr(l)*(y/r)*ylm(x,y,z,ll,m)+uvr(l)*dylm(x,y,z,ll,m,2)
               ppg%duv(j,ilma,3) = duvr(l)*(z/r)*ylm(x,y,z,ll,m)+uvr(l)*dylm(x,y,z,ll,m,3)
             else
               ppg%duv(j,ilma,1) = uvr(l)*dylm(x,y,z,ll,m,1)
               ppg%duv(j,ilma,2) = uvr(l)*dylm(x,y,z,ll,m,2)
               ppg%duv(j,ilma,3) = uvr(l)*dylm(x,y,z,ll,m,3)
             end if
           end if
         enddo
       enddo
       l0=l
       enddo
 
     enddo
 !!$omp end do
 !!$omp end parallel
 
  enddo

  lma=0
  do a=1,natom
    ik=kion(a)
    if(property /= 'update_wo_realloc') then
    l0=0
    do ll=0,pp%mlps(ik)
    do l=l0,l0+pp%nproj(ll,ik)-1
      if(pp%inorm(l,ik)==0) cycle
      do m=-ll,ll
        lma=lma+1
        ppg%rinv_uvu(lma)=dble(pp%inorm(l,ik))*hvol
      enddo
    enddo
    l0=l
    enddo
    endif
  enddo
 
end subroutine calc_uv


subroutine spline(Np,xn,yn,an,bn,cn,dn)
  integer,intent(in) :: Np
  real(8),intent(in) :: xn(0:Np-1),yn(0:Np-1)
  real(8),intent(out) :: an(0:Np-2),bn(0:Np-2),cn(0:Np-2),dn(0:Np-2)
  integer :: i,Npm2,info
  real(8) :: dxn(0:Np-1),dyn(0:Np-1),u(1:Np-2),v(1:Np-2),Amat(1:Np-2,1:Np-2)
  real(8) :: Amat_t(1:Np-2,1:Np-2)
! for lapack
  integer :: LWORK
  integer, allocatable :: IPIV(:) ! dimension N
  real(8), allocatable :: WORK(:) ! dimension LWORK
! for check inverse matrix problem
!  integer :: j,k
!  real(8) :: Amat_chk(1:Np-2,1:Np-2)
!  real(8) :: ss

  Npm2 = Np-2
  LWORK = Npm2*Npm2*6
  allocate(IPIV(Npm2),WORK(LWORK))


  do i = 0,Np-2
    dxn(i) = xn(i+1) - xn(i)
    dyn(i) = yn(i+1) - yn(i)
  end do

  do i = 1,Npm2
    v(i) = 6d0*(dyn(i)/dxn(i) - dyn(i-1)/dxn(i-1))
  end do

  Amat = 0d0
  Amat(1,1) = 2d0*(dxn(1) + dxn(0))
  Amat(1,2) = dxn(1)
  do i = 2,Npm2-1
    Amat(i,i+1) = dxn(i)
    Amat(i,i  ) = 2d0*(dxn(i)+dxn(i-1))
    Amat(i,i-1) = dxn(i-1)
  end do
  Amat(Npm2,Npm2  ) = 2d0*(dxn(Npm2)+dxn(Npm2-1))
  Amat(Npm2,Npm2-1) = dxn(Npm2-1)

! inverse matrix problem
  Amat_t = Amat


  call DGETRF(Npm2, Npm2, Amat_t, Npm2, IPIV, info)  ! factorize
  call DGETRI(Npm2, Amat_t, Npm2, IPIV, WORK, LWORK, info)  ! inverse

!  check inverse matrix problem
!  do i = 1,Npm2
!    do j = 1,Npm2
!      ss = 0d0
!      do k = 1,Npm2
!        ss = ss + Amat(i,k)*Amat_t(k,j)
!      end do
!      Amat_chk(i,j) = ss
!    end do
!  end do
!
!  do i = 1,Npm2
!    write(*,'(999e16.6e3)')(Amat_chk(i,j),j=1,Npm2)
!  end do
!
!  stop


  do i = 1,Npm2
    u(i) = sum(Amat_t(i,:)*v(:))
  end do

! for b
  bn(0) = 0d0
  bn(1:Np-2) = 0.5d0*u(1:Np-2)
! for a
  do i = 0,Npm2-1
    an(i) = (u(i+1) -2d0*bn(i))/(6d0*dxn(i))
  end do
  an(Npm2) = (0d0 -2d0*bn(Npm2))/(6d0*dxn(Npm2))
! for d
  dn(0:Npm2) = yn(0:Npm2)
! for c
  i=0
  cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*0.d0)/6d0
  do i = 1,Npm2-1
     cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*u(i))/6d0
  end do
  cn(Npm2) = dyn(Npm2)/dxn(Npm2) - dxn(Npm2)*(0d0+2d0*u(Npm2))/6d0

  return
end subroutine spline

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
  use salmon_communication, only: comm_get_groupinfo &
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
end subroutine init_uvpsi_summation

subroutine finalize_uvpsi_summation(ppg)
  use structures,    only: s_pp_grid
  use salmon_communication, only: comm_free_group
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
end subroutine

end module prep_pp_sub
