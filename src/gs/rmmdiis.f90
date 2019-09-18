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
module rmmdiis_sub
  implicit none

  public  :: rmmdiis
  private :: axpyv,axpyzv,scalev,copyv,dotv,normv,diffv,setv
  private :: hpsi_test_diis

  interface inner_product

    module procedure r_inner_product, c_inner_product

  end interface

  INTERFACE eigenval
    MODULE PROCEDURE R_eigenval, C_eigenval
  END INTERFACE

  INTERFACE gen_eigen
    MODULE PROCEDURE R_gen_eigen, C_gen_eigen
  END INTERFACE

contains

!=======================================================================
!============================================================== RMM-DIIS
! This routine is RMM-DIIS
! J. Soc. Mat. Sci., Japan, vol.52 (3), p.260-265. (in Japanese)

subroutine rmmdiis(mg,system,info,stencil,srg_ob_1,spsi,energy,itotmst  &
                  ,mst,iflag_diisjump,norm_diff_psi_stock  &
                  ,info_ob,ppg,vlocal)
  use inputoutput, only: ncg,lambda1_diis,lambda2_diis
  use structures, only: s_rgrid,s_dft_system,s_orbital_parallel,s_orbital,   &
                        s_dft_energy,s_stencil,s_scalar,s_pp_grid
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use calc_allob_sub
  use sendrecv_grid, only: s_sendrecv_grid
  implicit none
  
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel) :: info
  type(s_orbital) :: spsi
  type(s_dft_energy) :: energy
  type(s_stencil) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg_ob_1
  type(s_pp_grid) :: ppg
  integer,intent(in)    :: itotmst
  integer,intent(in)    :: mst(2)
  integer,intent(out)   :: iflag_diisjump
  real(8),intent(out)   :: norm_diff_psi_stock(itotmst,1)
  type(s_orbital_parallel)       :: info_ob
  type(s_scalar),intent(in) :: vlocal(system%nspin)
  integer,parameter :: nd=4
  integer :: iob,iob_allob,iter,ix,iy,iz
  integer :: nspin_1
  type(s_orbital)  :: stpsi
  type(s_orbital)  :: shtpsi
  type(s_scalar),allocatable :: v(:)
  integer,allocatable :: iflagdiis(:)
  integer,allocatable :: iobcheck(:,:)
  real(8),allocatable :: phi(:,:,:,:)
  real(8),allocatable :: htphi(:,:,:)
  real(8),allocatable :: R1(:,:,:,:)
  real(8),allocatable :: phibar(:,:,:,:),Rbar(:,:,:,:)
  real(8),allocatable :: phibox(:,:,:),Rbox(:,:,:)
  real(8),allocatable :: psi_stock(:,:,:,:,:,:,:)
  real(8) :: rbox1
  real(8),allocatable :: epsdiis(:,:),Rnorm(:,:)
  real(8) :: rnorm_diff_psi(itotmst,1)
  integer :: is

  integer :: mg_xs,mg_xe,mg_ys,mg_ye,mg_zs,mg_ze
  integer :: numo
  type(s_dft_system) :: system_spin1 ! temporary

  mg_xs = mg%is(1) ; mg_xe = mg%ie(1)
  mg_ys = mg%is(2) ; mg_ye = mg%ie(2)
  mg_zs = mg%is(3) ; mg_ze = mg%ie(3)
  numo = info%numo
  
  allocate(stpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  nspin_1=1
  system_spin1%nspin = 1
  allocate(v(nspin_1))
  allocate(v(nspin_1)%f(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))

  allocate (htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))
  allocate (phibox(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))
  allocate (Rbox(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))
  allocate (phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0:ncg))
  allocate (R1(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0:ncg))
  allocate (phibar(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0:ncg))
  allocate (Rbar(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0:ncg))
  allocate (psi_stock(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1:system%nspin,  &
                      1:system%nspin*numo,info%ik_s:info%ik_e,1))
  
  allocate (iobcheck(1:itotmst,0:ncg))
  iobcheck=0
  
!$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    stpsi%rwf(ix,iy,iz,1,1,1,1)=0.d0
  end do
  end do
  end do
  
  if(system%nspin*numo>=1)then
    allocate (iflagdiis(1:system%nspin*numo))
    allocate (epsdiis(1:system%nspin*numo,0:ncg))
    allocate (Rnorm(1:system%nspin*numo,0:ncg))
  end if 
  
  ! Flag for convergence
  if(system%nspin*numo >= 1) iflagdiis=1
  
  if(system%nspin*numo >= 1) then
    phi=0.d0
    psi_stock(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1:system%nspin,  &
                      1:numo,info%ik_s:info%ik_e,1)=   &
      spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1:system%nspin,  &
                      1:numo,info%ik_s:info%ik_e,1)
  end if
  
  iflag_diisjump=0
  
  do iob=1,system%nspin*numo
    call calc_allob(iob,info,iob_allob,itotmst,mst,system%nspin*numo)
    if(iob>numo)then
      is=2
    else
      is=1
    end if
  
    call setv(mg,vlocal,v,iob_allob,mst)

    Iteration : do iter=1,ncg
  
    if(iter == 1) then
  ! Obtain residual vector R_0
      call copyv(mg,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1) &
                   ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0))

      call copyv(mg,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0) &
                   ,stpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1))

      call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)

      call copyv(mg,shtpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1) &
                   ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))

      call inner_product3(mg,info,phi(mg_xs,mg_ys,mg_zs,0),htphi(mg_xs,mg_ys,mg_zs),rbox1)
  
      call axpyzv(mg,-rbox1*system%hvol &
                    ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0) &
                    ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze) &
                    ,R1(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0))

      epsdiis(iob,0)=rbox1*system%hvol
      call inner_product3(mg,info,R1(mg_xs,mg_ys,mg_zs,0),R1(mg_xs,mg_ys,mg_zs,0),rbox1)
      Rnorm(iob,0)=rbox1*system%hvol
  
    else
  ! Solve by Lagrange's method of undetermined multipliers, and obtain 
  ! Rbar from previous combinations of phi and R.
      if(iflagdiis(iob) == 1)then
        call diis_core(mg,info,itotmst,system%hvol,phi,R1,phibar,Rbar,iob,iter,iobcheck)
      end if
    end if
  
    if(iflagdiis(iob) == 1)then
  
      if(iter == 1) then
        call axpyv(mg,-lambda1_diis &
                     ,R1(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0) &
                     ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0))
      else
        call axpyzv(mg,-lambda2_diis &
                      ,Rbar(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter-1) &
                      ,phibar(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter-1) &
                      ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter))
      end if
  
  ! normalization
      call inner_product3(mg,info,phi(mg_xs,mg_ys,mg_zs,iter),phi(mg_xs,mg_ys,mg_zs,iter),rbox1)

      call scalev(mg,1.0d0/sqrt(rbox1*system%hvol) &
                    ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter))
  
      call copyv(mg,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter) &
                   ,stpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1))

      call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)

      call copyv(mg,shtpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1) &
                   ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))
  
      call inner_product3(mg,info,phi(mg_xs,mg_ys,mg_zs,iter),htphi(mg_xs,mg_ys,mg_zs),rbox1)

      call axpyzv(mg,-rbox1*system%hvol &
                    ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter) &
                    ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze) &
                    ,R1(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter))

      call inner_product3(mg,info,phi(mg_xs,mg_ys,mg_zs,iter),htphi(mg_xs,mg_ys,mg_zs),rbox1)
      epsdiis(iob,iter)=rbox1*system%hvol
  
      call inner_product3(mg,info,R1(mg_xs,mg_ys,mg_zs,iter),R1(mg_xs,mg_ys,mg_zs,iter),rbox1)
      Rnorm(iob,iter)=rbox1*system%hvol
  
  ! judgement for closing loop.
  ! The ratio of Rnorm is set to 0.3 as well as Kresse-Furthmuller.
      if(iter >= 2) then
        if(iter >= 3 .and. epsdiis(iob,iter) > epsdiis(iob,iter-1)) then
          call copyv(mg,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter-1) &
                       ,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1))
          iflagdiis(iob)=0
        else if(-(epsdiis(iob,iter)-epsdiis(iob,iter-1)) <= 1.0d-8 .or.      &
                Rnorm(iob,iter)/Rnorm(iob,0) <= 0.3d0) then
          call copyv(mg,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter) &
                       ,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1))
          iflagdiis(iob)=0
        end if
      end if
  
      if(iter == ncg) then
        call copyv(mg,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1) &
                     ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,ncg))
      end if

      if(iter == 1 .and. iflag_diisjump == 1) then
        call copyv(mg,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1) &
                     ,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1))
      end if
 
      call copyv(mg,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,iter) &
                   ,stpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1))

      call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)

      call copyv(mg,shtpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1) &
                   ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))

      call inner_product3(mg,info,phi(mg_xs,mg_ys,mg_zs,iter),htphi(mg_xs,mg_ys,mg_zs),rbox1)
      
      end if
  
    end do Iteration
  
  end do        ! loop for iob
  
  iflag_diisjump=0
  do iob=1,system%nspin*numo
    call calc_allob(iob,info,iob_allob,itotmst,mst,system%nspin*numo)
    if(iob>numo)then
      is=2
    else
      is=1
    end if

    call setv(mg,vlocal,v,iob_allob,mst)

    call copyv(mg,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1) &
                 ,stpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1))

    call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)

    call copyv(mg,shtpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1) &
                 ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))

    rbox1=system%hvol* &
          dotv(mg,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1) &
                 ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))

    if(rbox1-energy%esp(iob,1,1)>5.d0) iflag_diisjump=1
  end do
  
  if(iflag_diisjump==0)then
    continue
  else if(iflag_diisjump==1)then
    spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1:system%nspin,  &
                      1:numo,info%ik_s:info%ik_e,1)=   &
      psi_stock(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1:system%nspin,  &
                      1:numo,info%ik_s:info%ik_e,1)
    do iob=1,system%nspin*numo
      call calc_allob(iob,info,iob_allob,itotmst,mst,system%nspin*numo)
      if(iob>numo)then
        is=2
      else
        is=1
      end if
    
      call setv(mg,vlocal,v,iob_allob,mst)

      call copyv(mg,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1) &
                   ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0))

      call copyv(mg,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0) &
                   ,stpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1))

      call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)

      call copyv(mg,shtpsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,1,1,1,1) &
                   ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze))
  
      call inner_product3(mg,info,phi(mg_xs,mg_ys,mg_zs,0),htphi(mg_xs,mg_ys,mg_zs),rbox1)

      call axpyzv(mg,-rbox1*system%hvol &
                    ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0) &
                    ,htphi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze) &
                    ,R1(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0))

      call axpyzv(mg,lambda2_diis &
                    ,R1(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0) &
                    ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0) &
                    ,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1))
    end do
  
    rnorm_diff_psi=0.d0
    do iob=1,system%nspin*numo
      if(iob>numo)then
        is=2
      else
        is=1
      end if
      call diffv(mg,spsi%rwf(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1) &
                   ,psi_stock(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,is,iob-(is-1)*numo,1,1) &
                   ,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0))
      rbox1=normv(mg,phi(mg_xs:mg_xe,mg_ys:mg_ye,mg_zs:mg_ze,0))*system%hvol
      rnorm_diff_psi(iob,1)=rbox1
    end do
    call comm_summation(rnorm_diff_psi,norm_diff_psi_stock,itotmst,nproc_group_global)
  end if
  
  
  deallocate(htphi)
  deallocate(phibox,Rbox,phi,R1,phibar,Rbar)
  
  if(system%nspin*numo>=1)then
    deallocate (iflagdiis,epsdiis,Rnorm)
  end if 
  deallocate(iobcheck) 
  
  deallocate(stpsi%rwf,shtpsi%rwf)
  deallocate(v(nspin_1)%f)
  deallocate(v)

  return
  
end subroutine rmmdiis

function dotv(mg,x,y) result(ret)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in)       :: x(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(in)       :: y(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: ix,iy,iz
  integer :: mg_xs,mg_xe,mg_ys,mg_ye,mg_zs,mg_ze
  real(8) :: ret

  mg_xs = mg%is(1) ; mg_xe = mg%ie(1)
  mg_ys = mg%is(2) ; mg_ye = mg%ie(2)
  mg_zs = mg%is(3) ; mg_ze = mg%ie(3)

  ret = 0.d0
!$OMP parallel do private(iz,iy,ix) reduction(+:ret)
  do iz=mg_zs,mg_ze
  do iy=mg_ys,mg_ye
  do ix=mg_xs,mg_xe
    ret = ret + x(ix,iy,iz) * y(ix,iy,iz)
  end do
  end do
  end do
end function

function normv(mg,x) result(ret)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in)       :: x(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: ret
  ret = dotv(mg,x,x)
end function

subroutine scalev(mg,a,x)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in)       :: a
  real(8),intent(inout)    :: x(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: ix,iy,iz
  integer :: mg_xs,mg_xe,mg_ys,mg_ye,mg_zs,mg_ze

  mg_xs = mg%is(1) ; mg_xe = mg%ie(1)
  mg_ys = mg%is(2) ; mg_ye = mg%ie(2)
  mg_zs = mg%is(3) ; mg_ze = mg%ie(3)

!$OMP parallel do private(iz,iy,ix) firstprivate(a)
  do iz=mg_zs,mg_ze
  do iy=mg_ys,mg_ye
  do ix=mg_xs,mg_xe
    x(ix,iy,iz) = a * x(ix,iy,iz)
  end do
  end do
  end do
end subroutine

subroutine diffv(mg,x,y,z)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in)       :: x(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(in)       :: y(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out)      :: z(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: ix,iy,iz
  integer :: mg_xs,mg_xe,mg_ys,mg_ye,mg_zs,mg_ze

  mg_xs = mg%is(1) ; mg_xe = mg%ie(1)
  mg_ys = mg%is(2) ; mg_ye = mg%ie(2)
  mg_zs = mg%is(3) ; mg_ze = mg%ie(3)

!$OMP parallel do private(iz,iy,ix)
  do iz=mg_zs,mg_ze
  do iy=mg_ys,mg_ye
  do ix=mg_xs,mg_xe
    z(ix,iy,iz) = abs(x(ix,iy,iz) - y(ix,iy,iz))
  end do
  end do
  end do
end subroutine

subroutine axpyv(mg,a,x,y)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in)       :: a
  real(8),intent(in)       :: x(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(inout)    :: y(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: ix,iy,iz
  integer :: mg_xs,mg_xe,mg_ys,mg_ye,mg_zs,mg_ze

  mg_xs = mg%is(1) ; mg_xe = mg%ie(1)
  mg_ys = mg%is(2) ; mg_ye = mg%ie(2)
  mg_zs = mg%is(3) ; mg_ze = mg%ie(3)

!$OMP parallel do private(iz,iy,ix) firstprivate(a)
  do iz=mg_zs,mg_ze
  do iy=mg_ys,mg_ye
  do ix=mg_xs,mg_xe
    y(ix,iy,iz) = a * x(ix,iy,iz) + y(iz,iy,iz)
  end do
  end do
  end do
end subroutine

subroutine axpyzv(mg,a,x,y,z)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in)       :: a
  real(8),intent(in)       :: x(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(in)       :: y(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out)      :: z(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: ix,iy,iz
  integer :: mg_xs,mg_xe,mg_ys,mg_ye,mg_zs,mg_ze

  mg_xs = mg%is(1) ; mg_xe = mg%ie(1)
  mg_ys = mg%is(2) ; mg_ye = mg%ie(2)
  mg_zs = mg%is(3) ; mg_ze = mg%ie(3)

!$OMP parallel do private(iz,iy,ix) firstprivate(a)
  do iz=mg_zs,mg_ze
  do iy=mg_ys,mg_ye
  do ix=mg_xs,mg_xe
    z(ix,iy,iz) = a * x(ix,iy,iz) + y(iz,iy,iz)
  end do
  end do
  end do
end subroutine

subroutine copyv(mg,src,dst)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: mg
  real(8),intent(in)       :: src(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out)      :: dst(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: ix,iy,iz
  integer :: mg_xs,mg_xe,mg_ys,mg_ye,mg_zs,mg_ze

  mg_xs = mg%is(1) ; mg_xe = mg%ie(1)
  mg_ys = mg%is(2) ; mg_ye = mg%ie(2)
  mg_zs = mg%is(3) ; mg_ze = mg%ie(3)

!$OMP parallel do private(iz,iy,ix)
  do iz=mg_zs,mg_ze
  do iy=mg_ys,mg_ye
  do ix=mg_xs,mg_xe
    dst(ix,iy,iz) = src(ix,iy,iz)
  end do
  end do
  end do
end subroutine

subroutine setv(mg,vlocal,v,iob_allob,mst)
  use inputoutput, only: ispin
  use structures, only: s_rgrid,s_orbital_parallel,s_orbital,s_stencil,s_scalar,s_pp_grid
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_scalar)        :: v(1)
  type(s_scalar),intent(in) :: vlocal(ispin+1)
  integer,intent(in)    :: iob_allob
  integer,intent(in)    :: mst(2)
  integer :: ix,iy,iz
  integer :: mg_xs,mg_xe,mg_ys,mg_ye,mg_zs,mg_ze

  mg_xs = mg%is(1) ; mg_xe = mg%ie(1)
  mg_ys = mg%is(2) ; mg_ye = mg%ie(2)
  mg_zs = mg%is(3) ; mg_ze = mg%ie(3)

  if(iob_allob<=mst(1))then
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg_zs,mg_ze
    do iy=mg_ys,mg_ye
    do ix=mg_xs,mg_xe
      v(1)%f(ix,iy,iz) = vlocal(1)%f(ix,iy,iz)
    end do
    end do
    end do
  else
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg_zs,mg_ze
    do iy=mg_ys,mg_ye
    do ix=mg_xs,mg_xe
      v(1)%f(ix,iy,iz) = vlocal(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if

end subroutine

subroutine hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,srg_ob_1,ppg)
  use structures
  use hpsi_sub, only: hpsi
  use sendrecv_grid, only: s_sendrecv_grid
  implicit none
  type(s_orbital)  :: stpsi
  type(s_orbital)  :: shtpsi
  type(s_orbital_parallel)       :: info_ob
  type(s_rgrid),intent(in) :: mg
  type(s_scalar)        :: v(1)
  integer :: nspin_1
  type(s_stencil) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg_ob_1
  type(s_pp_grid) :: ppg
  type(s_dft_system) :: system_spin1 ! temporary
  system_spin1%nspin = nspin_1

  call hpsi(stpsi,shtpsi,info_ob,mg,v,system_spin1,stencil,srg_ob_1,ppg)

end subroutine hpsi_test_diis

subroutine r_inner_product(mg,info,matbox1,matbox2,rbox2)
  use structures, only: s_rgrid, s_orbital_parallel
  use salmon_communication, only: comm_summation
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_orbital_parallel),intent(in) :: info
  real(8),intent(in) :: matbox1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(in) :: matbox2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out) :: rbox2
  integer :: ix,iy,iz
  real(8) :: rbox

  rbox=0.d0
  !$omp parallel do reduction(+ : rbox) private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rbox=rbox+matbox1(ix,iy,iz)*matbox2(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(rbox,rbox2,info%icomm_r)

end subroutine r_inner_product

!=======================================================================
subroutine c_inner_product(mg,info,matbox1,matbox2,cbox2)
  use structures, only: s_rgrid, s_orbital_parallel
  use salmon_communication, only: comm_summation
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_orbital_parallel),intent(in) :: info
  complex(8),intent(in)  :: matbox1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  complex(8),intent(in)  :: matbox2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  complex(8),intent(out) :: cbox2
  integer :: ix,iy,iz
  complex(8) :: cbox

  cbox=0.d0
  !$omp parallel do reduction(+ : cbox) private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    cbox=cbox+conjg(matbox1(ix,iy,iz))*matbox2(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(cbox,cbox2,info%icomm_r)

end subroutine c_inner_product

subroutine inner_product3(mg,info,rmatbox1,rmatbox2,rbox2)
  use structures, only: s_rgrid, s_orbital_parallel
  use salmon_communication, only: comm_summation
  use timer
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_orbital_parallel),intent(in) :: info
  real(8),intent(in) :: rmatbox1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(in) :: rmatbox2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out) :: rbox2
  integer :: ix,iy,iz
  real(8) :: rbox

  rbox=0.d0
  !$omp parallel do reduction(+ : rbox) private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rbox=rbox+rmatbox1(ix,iy,iz)*rmatbox2(ix,iy,iz)
  end do
  end do
  end do

  call timer_begin(LOG_ALLREDUCE_INNER_PRODUCT3)
  call comm_summation(rbox,rbox2,info%icomm_r)
  call timer_end(LOG_ALLREDUCE_INNER_PRODUCT3)

end subroutine inner_product3

subroutine diis_core(mg,info,itotmst,hvol,phi,R1,phibar,Rbar,iob,iter,pcheck)
  use inputoutput, only: ncg
  use structures, only: s_rgrid, s_orbital_parallel
  !$ use omp_lib
  implicit none

  type(s_rgrid),intent(in) :: mg
  type(s_orbital_parallel),intent(in) :: info
  integer,intent(in) :: itotmst
  real(8),intent(in) :: hvol
  integer :: ii,jj,iob,iter,ix,iy,iz,ier2
  integer :: ibox,icount
  real(8),allocatable :: Rmat(:,:),Smat(:,:)
  real(8),allocatable :: betav(:)
  real(8),allocatable :: alpha(:),eval(:),evec(:,:)
  real(8) :: evalbox
  real(8) :: rnorm
  real(8) :: rbox
  integer :: pcheck(1:itotmst,0:ncg)

  real(8) :: phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),      &
                           mg%is(3):mg%ie(3),0:ncg)
  real(8) :: R1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),      &
                           mg%is(3):mg%ie(3),0:ncg)
  real(8) :: phibar(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),      &
                           mg%is(3):mg%ie(3),0:ncg)
  real(8) :: Rbar(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),      &
                           mg%is(3):mg%ie(3),0:ncg)

  allocate(Rmat(iter,iter))
  allocate(Smat(iter,iter))
  allocate(alpha(iter),betav(iter),eval(iter))
  allocate(evec(iter,iter))

  do ii=0,iter-1
    do jj=0,iter-1
      call inner_product(mg,info,R1(:,:,:,ii),R1(:,:,:,jj),rbox)
      Rmat(ii+1,jj+1)=rbox*hvol

      call inner_product(mg,info,phi(:,:,:,ii),phi(:,:,:,jj),rbox)
      Smat(ii+1,jj+1)=rbox*hvol
    end do
  end do

  call eigenval(Smat,eval,iter)

  do ii=1,iter-1
    evalbox=eval(ii)
    ibox=ii
    do jj=ii+1,iter
      if(eval(jj) > evalbox)then
        evalbox=eval(jj)
        ibox=jj
      end if
    end do
    if(ibox /= ii)then
      eval(ibox) = eval(ii)
      eval(ii) = evalbox
    end if
  end do

  icount=0
  do ii=1,iter
    if(ii == 1)then
      if(abs(eval(1)-dble(iter)) < 1.0d-13)      &
        icount = icount + 1
    else
      if(eval(ii) < 1.0d-13) icount = icount + 1
    end if
  end do

  if(icount == iter)then
  ! if phibar is estimated to be equal mixture of previous phis,
  ! update phibar by newest phi
  !$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      phibar(ix,iy,iz,iter-1)=phi(ix,iy,iz,iter-1)
    end do
    end do
    end do

    call inner_product(mg,info,phibar(:,:,:,iter-1),phibar(:,:,:,iter-1),rbox)
    rnorm=sqrt(rbox*hvol)
  !$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      phibar(ix,iy,iz,iter-1)=phibar(ix,iy,iz,iter-1)/rnorm
      Rbar(ix,iy,iz,iter-1)=R1(ix,iy,iz,iter-1)
    end do
    end do
    end do

  else

    call gen_eigen(Rmat,Smat,alpha,betav,evec,iter,ier2)

    if(ier2 .ne. 0) then
  ! if Smat is not positive-definite,
  ! update phibar by newest phi
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        phibar(ix,iy,iz,iter-1)=phi(ix,iy,iz,iter-1)
      end do
      end do
      end do
      call inner_product(mg,info,phibar(:,:,:,iter-1),phibar(:,:,:,iter-1),rbox)

      rnorm=sqrt(rbox*hvol)
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        phibar(ix,iy,iz,iter-1)=phibar(ix,iy,iz,iter-1)/rnorm
      end do
      end do
      end do

      pcheck(iob,iter)=ier2

    else
      eval=alpha/betav

      evalbox=eval(1)
      ibox=1
      do ii=2,iter
        if(eval(ii) <= evalbox)then
          evalbox=eval(ii)
          ibox=ii
        end if
      end do

  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        phibar(ix,iy,iz,iter-1)=0d0
      end do
      end do
      end do
      do ii=1,iter
  !$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          phibar(ix,iy,iz,iter-1)=phibar(ix,iy,iz,iter-1)      &
                        +dble(evec(ii,ibox))*phi(ix,iy,iz,ii-1)
        end do
        end do
        end do
      end do

      call inner_product(mg,info,phibar(:,:,:,iter-1),phibar(:,:,:,iter-1),rbox)
      rnorm=sqrt(rbox*hvol)
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        phibar(ix,iy,iz,iter-1)=phibar(ix,iy,iz,iter-1)/rnorm
      end do
      end do
      end do

  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        Rbar(ix,iy,iz,iter-1)=0d0
      end do
      end do
      end do
      do ii=1,iter
  !$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          Rbar(ix,iy,iz,iter-1)=Rbar(ix,iy,iz,iter-1)      &
                       +dble(evec(ii,ibox))*R1(ix,iy,iz,ii-1)/rnorm
        end do
        end do
        end do
      end do

    end if
  end if

  deallocate(Rmat, Smat)
  deallocate(alpha, betav, eval, evec)

  return

end subroutine diis_core

subroutine R_eigenval(Smat, eval, iter)

implicit none
integer :: ii,iter
real(8) :: Smat(iter,iter),eval(iter)
character :: JOBVL*1, JOBVR*1
integer :: INFO, LDA, LDB, LDVL, LDVR, LWORK, N
real(8),allocatable :: A( :, : ), ALPHAI( : ), ALPHAR( : )
real(8),allocatable :: B( :, : ), BETA( : ), VL( :, : )
real(8),allocatable :: VR( :, : ), WORK( : )

JOBVL='N'
JOBVR='V'
N=iter
LDA=iter
LDB=iter

LDVL=iter
LDVR=iter

LWORK=8*iter

allocate ( A( LDA, N ), ALPHAI( N ), ALPHAR( N ) )
allocate ( B( LDB, N ), BETA( N ), VL( LDVL, N ) )
allocate ( VR( LDVR, N ), WORK( LWORK ) )

A=Smat

B=0.d0
do ii=1,iter
  B(ii,ii)=1.d0
end do

call DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,      &
       BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

! check whether eigenvalues are real.

do ii=1,iter
  if(ALPHAI(ii) > 0.05d0)then
    write(*,*) "========== DIIS error =========="
    write(*,*) "One of eigenvalue is imaginary."
    write(*,*) "================================"
    stop
  end if
end do

eval=ALPHAR

deallocate (A,ALPHAI,ALPHAR,B,BETA,VL,VR,WORK)

return

end subroutine R_eigenval

! This is a routine to solve generalised eigenvalue problem.

SUBROUTINE R_gen_eigen(Rmat,Smat,alpha,betav,evec,iter,ier2)

implicit none

integer :: iter,ii,ier2
real(8) :: Rmat(iter,iter)
real(8) :: Smat(iter,iter)
real(8) :: alpha(iter),betav(iter)
real(8) :: evec(iter,iter)
character :: JOBVL*1, JOBVR*1
integer :: INFO, LDA, LDB, LDVL, LDVR, LWORK, N
real(8),allocatable :: A( :, : ), ALPHAI( : ), ALPHAR( : )
real(8),allocatable :: B( :, : ), BETA( : ), VL( :, : )
real(8),allocatable :: VR( :, : ), WORK( : )

ier2=0

JOBVL='N'
JOBVR='V'
N=iter
LDA=iter
LDB=iter

LDVL=iter
LDVR=iter

LWORK=8*iter

allocate ( A( LDA, N ), ALPHAI( N ), ALPHAR( N ) )
allocate ( B( LDB, N ), BETA( N ), VL( LDVL, N ) )
allocate ( VR( LDVR, N ), WORK( LWORK ) )

A=Rmat
B=Smat

call DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,      &
       BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

! check whether eigenvalues and eigenvectors take real value.
do ii=1,iter
  if(ALPHAI(ii) > 1.d-8)then
    write(*,*) "========== DIIS error =========="
    write(*,*) "One of eigenvector is imaginary."
    write(*,*) "================================"
    stop
  end if
end do

alpha=ALPHAR
betav=BETA
evec=VR

deallocate (A,ALPHAI,ALPHAR,B,BETA,VL,VR,WORK)

return

end subroutine R_gen_eigen

!=======================================================================
!======================================================= RMM-DIIS_lapack
!==================================================== eigenvalue problem
! This is a routine to solve eigenvalue problem.
subroutine C_eigenval(Smat, eval, iter)

implicit none
integer :: ii,iter
complex(8) :: Smat(iter,iter)
real(8) :: eval(iter)
character :: JOBVL*1, JOBVR*1
integer :: INFO, LDA, LDB, LDVL, LDVR, LWORK, N
complex(8),allocatable :: A( :, : ), ALPHA( : )
complex(8),allocatable :: B( :, : ), BETA( : ), VL( :, : )
complex(8),allocatable :: VR( :, : ), WORK( : )

JOBVL='N'
JOBVR='V'
N=iter
LDA=iter
LDB=iter

LDVL=iter
LDVR=iter

LWORK=8*iter

allocate ( A( LDA, N ), ALPHA( N ) )
allocate ( B( LDB, N ), BETA( N ), VL( LDVL, N ) )
allocate ( VR( LDVR, N ), WORK( LWORK ) )

A=Smat

B=0.d0
do ii=1,iter
  B(ii,ii)=1.d0
end do

call ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA,      &
       BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

eval=real(ALPHA)

deallocate (A,ALPHA,B,BETA,VL,VR,WORK)

return

end subroutine C_eigenval

! This is a routine to solve generalised eigenvalue problem.

SUBROUTINE C_gen_eigen(Rmat,Smat,alpha,betav,evec,iter,ier2)

implicit none

integer :: iter,ier2
complex(8) :: Rmat(iter,iter)
complex(8) :: Smat(iter,iter)
complex(8) :: alpha(iter),betav(iter)
complex(8) :: evec(iter,iter)
character :: JOBVL*1, JOBVR*1
integer :: INFO, LDA, LDB, LDVL, LDVR, LWORK, N
complex(8),allocatable :: A( :, : ), ALPHA2( : )
complex(8),allocatable :: B( :, : ), BETA( : ), VL( :, : )
complex(8),allocatable :: VR( :, : ), WORK( : )

ier2=0

JOBVL='N'
JOBVR='V'
N=iter
LDA=iter
LDB=iter

LDVL=iter
LDVR=iter

LWORK=8*iter

allocate ( A( LDA, N ), ALPHA2( N ) )
allocate ( B( LDB, N ), BETA( N ), VL( LDVL, N ) )
allocate ( VR( LDVR, N ), WORK( LWORK ) )

A=Rmat
B=Smat

call ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA2,      &
       BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

alpha=ALPHA2
betav=BETA
evec=VR

deallocate (A,ALPHA2,B,BETA,VL,VR,WORK)

return

end subroutine C_gen_eigen

end module rmmdiis_sub
