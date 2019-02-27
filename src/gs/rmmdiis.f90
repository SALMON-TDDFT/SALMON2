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
module rmmdiis_sub
  implicit none

contains

!=======================================================================
!============================================================== RMM-DIIS
! This routine is RMM-DIIS
! J. Soc. Mat. Sci., Japan, vol.52 (3), p.260-265. (in Japanese)

subroutine rmmdiis(mg,nspin,info,stencil,spsi,itotmst,mst,num_kpoints_rd,hvol,iflag_diisjump,elp3,esp,norm_diff_psi_stock,   &
                   info_ob,bnmat,cnmat,hgs,ppg,vlocal,iparaway_ob)
  use inputoutput, only: ncg,ispin,lambda1_diis,lambda2_diis
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use calc_allob_sub
  use diis_core_sub
  !$ use omp_lib
  implicit none
  
  type(s_rgrid),intent(in) :: mg
  integer,intent(in)   :: nspin
  type(s_wf_info) :: info
  type(s_wavefunction) :: spsi
  type(s_stencil) :: stencil
  type(s_pp_grid) :: ppg
  integer,intent(in)    :: itotmst
  integer,intent(in)    :: mst(2)
  integer,intent(in)    :: num_kpoints_rd
  real(8),intent(in)    :: hvol
  integer,intent(out)   :: iflag_diisjump
  real(8),intent(out)   :: elp3(3000)
  real(8),intent(in)    :: esp(itotmst,num_kpoints_rd)
  real(8),intent(out)   :: norm_diff_psi_stock(itotmst,1)
  type(s_wf_info)       :: info_ob
  real(8),intent(in)    :: cnmat(0:12,0:12),bnmat(0:12,0:12)
  real(8),intent(in)    :: hgs(3)
  real(8),intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),nspin)
  integer,intent(in)    :: iparaway_ob
  integer,parameter :: nd=4
  integer :: j,ind
  integer :: iob,iob_allob,iter,ix,iy,iz
  integer :: nspin_1
  type(s_wavefunction)  :: stpsi
  type(s_wavefunction)  :: shtpsi
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
  
  allocate(stpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  nspin_1=1
  allocate(v(nspin_1))
  allocate(v(nspin_1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

  allocate (htphi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (phibox(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (Rbox(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0:ncg))
  allocate (R1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0:ncg))
  allocate (phibar(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0:ncg))
  allocate (Rbar(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0:ncg))
  allocate (psi_stock(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:nspin,  &
                      1:nspin*info%numo,info%ik_s:info%ik_e,1))
  
  allocate (iobcheck(1:itotmst,0:ncg))
  iobcheck=0
  
  !$OMP parallel do private(iz,iy,ix)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    stpsi%rwf(ix,iy,iz,1,1,1,1)=0.d0
  end do
  end do
  end do
  
  if(nspin*info%numo>=1)then
    allocate (iflagdiis(1:nspin*info%numo))
    allocate (epsdiis(1:nspin*info%numo,0:ncg))
    allocate (Rnorm(1:nspin*info%numo,0:ncg))
  end if 
  
  ! Flag for convergence
  if(nspin*info%numo >= 1) iflagdiis=1
  
  if(nspin*info%numo >= 1) then
    phi=0.d0
    psi_stock(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:nspin,  &
                      1:info%numo,info%ik_s:info%ik_e,1)=   &
      spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:nspin,  &
                      1:info%numo,info%ik_s:info%ik_e,1)
  end if
  
  iflag_diisjump=0
  
  do iob=1,nspin*info%numo
    call calc_allob(iob,iob_allob,iparaway_ob,itotmst,mst,nspin*info%numo)
    if(iob>info%numo)then
      is=2
    else
      is=1
    end if
  
    call setv(mg,vlocal,v,iob_allob,mst)

    Iteration : do iter=1,ncg
  
    if(iter == 1) then
  ! Obtain residual vector R_0
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
        phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,0)=   &
          spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,is,iob-(is-1)*info%numo,1,1)
      end do

  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        stpsi%rwf(ix,iy,iz,1,1,1,1)=phi(ix,iy,iz,0)
      end do
      end do
      end do

      call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,ppg)

  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        htphi(ix,iy,iz)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
      end do
      end do
      end do

      call inner_product3(mg,phi(mg%is(1),mg%is(2),mg%is(3),0),htphi(mg%is(1),mg%is(2),mg%is(3)),rbox1,elp3)
  
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
        R1(:,:,iz,0)=htphi(:,:,iz)-rbox1*hvol*phi(:,:,iz,0)
      end do 
      epsdiis(iob,0)=rbox1*hvol
      call inner_product3(mg,R1(mg%is(1),mg%is(2),mg%is(3),0),R1(mg%is(1),mg%is(2),mg%is(3),0),rbox1,elp3)
      Rnorm(iob,0)=rbox1*hvol
  
    else
  ! Solve by Lagrange's method of undetermined multipliers, and obtain 
  ! Rbar from previous combinations of phi and R.
      if(iflagdiis(iob) == 1)then
        call diis_core(mg,itotmst,hvol,phi,R1,phibar,Rbar,iob,iter,iobcheck)
      end if
    end if
  
    if(iflagdiis(iob) == 1)then
  
      if(iter == 1) then
  !$OMP parallel do
        do iz=mg%is(3),mg%ie(3)
          phi(:,:,iz,iter)=phi(:,:,iz,0)-lambda1_diis*R1(:,:,iz,0)
        end do
      else
  !$OMP parallel do
        do iz=mg%is(3),mg%ie(3)
          phi(:,:,iz,iter)=phibar(:,:,iz,iter-1)-lambda2_diis*Rbar(:,:,iz,iter-1)
        end do
      end if
  
  ! normalization
      call inner_product3(mg,phi(mg%is(1),mg%is(2),mg%is(3),iter),phi(mg%is(1),mg%is(2),mg%is(3),iter),rbox1,elp3)
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
        phi(:,:,iz,iter)=phi(:,:,iz,iter)/sqrt(rbox1*hvol)
      end do
  
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        stpsi%rwf(ix,iy,iz,1,1,1,1)=phi(ix,iy,iz,iter)
      end do
      end do
      end do

      call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,ppg)

  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        htphi(ix,iy,iz)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
      end do
      end do
      end do
  
      call inner_product3(mg,phi(mg%is(1),mg%is(2),mg%is(3),iter),htphi(mg%is(1),mg%is(2),mg%is(3)),rbox1,elp3)
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
        R1(:,:,iz,iter)=htphi(:,:,iz)-rbox1*hvol*phi(:,:,iz,iter)
      end do
  
      call inner_product3(mg,phi(mg%is(1),mg%is(2),mg%is(3),iter),htphi(mg%is(1),mg%is(2),mg%is(3)),rbox1,elp3)
      epsdiis(iob,iter)=rbox1*hvol
  
      call inner_product3(mg,R1(mg%is(1),mg%is(2),mg%is(3),iter),R1(mg%is(1),mg%is(2),mg%is(3),iter),rbox1,elp3)
      Rnorm(iob,iter)=rbox1*hvol
  
  ! judgement for closing loop.
  ! The ratio of Rnorm is set to 0.3 as well as Kresse-Furthmuller.
      if(iter >= 2) then
        if(iter >= 3 .and. epsdiis(iob,iter) > epsdiis(iob,iter-1)) then
  !$OMP parallel do
          do iz=mg%is(3),mg%ie(3)
            spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,is,iob-(is-1)*info%numo,1,1) =   &
              phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,iter-1)
          end do
          iflagdiis(iob)=0
        else if(-(epsdiis(iob,iter)-epsdiis(iob,iter-1)) <= 1.0d-8 .or.      &
                Rnorm(iob,iter)/Rnorm(iob,0) <= 0.3d0) then
  !$OMP parallel do
          do iz=mg%is(3),mg%ie(3)
            spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,is,iob-(is-1)*info%numo,1,1) =   &
              phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,iter)
          end do
          iflagdiis(iob)=0
        end if
      end if
  
      if(iter == ncg) then
  !$OMP parallel do
        do iz=mg%is(3),mg%ie(3)
          spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,is,iob-(is-1)*info%numo,1,1) =   &
            phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,ncg)
        end do
       end if
      if(iter == 1 .and. iflag_diisjump == 1) then
  !$OMP parallel do
        do iz=mg%is(3),mg%ie(3)
          spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,is,iob-(is-1)*info%numo,1,1) =   &
            phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,1)
        end do
      end if
  
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        stpsi%rwf(ix,iy,iz,1,1,1,1)=phi(ix,iy,iz,iter)
      end do
      end do
      end do

      call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,ppg)

  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        htphi(ix,iy,iz)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
      end do
      end do
      end do

      call inner_product3(mg,phi(mg%is(1),mg%is(2),mg%is(3),iter),htphi(mg%is(1),mg%is(2),mg%is(3)),rbox1,elp3)
      
      end if
  
    end do Iteration
  
  end do        ! loop for iob
  
  iflag_diisjump=0
  do iob=1,nspin*info%numo
    call calc_allob(iob,iob_allob,iparaway_ob,itotmst,mst,nspin*info%numo)
    if(iob>info%numo)then
      is=2
    else
      is=1
    end if

    call setv(mg,vlocal,v,iob_allob,mst)

  !$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      stpsi%rwf(ix,iy,iz,1,1,1,1)=spsi%rwf(ix,iy,iz,is,iob-(is-1)*info%numo,1,1)
    end do
    end do
    end do

    call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,ppg)

  !$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      htphi(ix,iy,iz)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
    end do
    end do
    end do

    rbox1=sum(spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob-(is-1)*info%numo,1,1)* &
              htphi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))*hvol

    if(rbox1-esp(iob,1)>5.d0) iflag_diisjump=1
  end do
  
  if(iflag_diisjump==0)then
    continue
  else if(iflag_diisjump==1)then
    spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:nspin,  &
                      1:info%numo,info%ik_s:info%ik_e,1)=   &
      psi_stock(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:nspin,  &
                      1:info%numo,info%ik_s:info%ik_e,1)
    do iob=1,nspin*info%numo
      call calc_allob(iob,iob_allob,iparaway_ob,itotmst,mst,nspin*info%numo)
      if(iob>info%numo)then
        is=2
      else
        is=1
      end if
    
      call setv(mg,vlocal,v,iob_allob,mst)

  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
        phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,0)=   &
          spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,is,iob-(is-1)*info%numo,1,1)
      end do
  
  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        stpsi%rwf(ix,iy,iz,1,1,1,1)=phi(ix,iy,iz,0)
      end do
      end do
      end do

      call hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,ppg)

  !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        htphi(ix,iy,iz)=shtpsi%rwf(ix,iy,iz,1,1,1,1)
      end do
      end do
      end do
  
      call inner_product3(mg,phi(mg%is(1),mg%is(2),mg%is(3),0),htphi(mg%is(1),mg%is(2),mg%is(3)),rbox1,elp3)
  
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
        R1(:,:,iz,0)=htphi(:,:,iz)-rbox1*hvol*phi(:,:,iz,0)
        spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,is,iob-(is-1)*info%numo,1,1)=   &
          phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,0)+  &
            lambda2_diis*R1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),iz,0)
      end do 
  
    end do
  
    rnorm_diff_psi=0.d0
    do iob=1,nspin*info%numo
      if(iob>info%numo)then
        is=2
      else
        is=1
      end if
      phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0)=    &
         abs(spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob-(is-1)*info%numo,1,1)-   &
             psi_stock(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob-(is-1)*info%numo,1,1))
      rbox1=sum(phi(:,:,:,0)*phi(:,:,:,0))*hvol
      rnorm_diff_psi(iob,1)=rbox1
    end do
    call comm_summation(rnorm_diff_psi,norm_diff_psi_stock,itotmst,nproc_group_global)
  end if
  
  
  deallocate(htphi)
  deallocate(phibox,Rbox,phi,R1,phibar,Rbar)
  
  if(nspin*info%numo>=1)then
    deallocate (iflagdiis,epsdiis,Rnorm)
  end if 
  deallocate(iobcheck) 
  
  deallocate(stpsi%rwf,shtpsi%rwf)
  deallocate(v(nspin_1)%f)
  deallocate(v)

  return
  
end subroutine rmmdiis

subroutine setv(mg,vlocal,v,iob_allob,mst)
  use inputoutput, only: ispin
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_scalar)        :: v(1)
  real(8),intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  integer,intent(in)    :: iob_allob
  integer,intent(in)    :: mst(2)
  integer :: ix,iy,iz

  if(iob_allob<=mst(1))then
  !$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,1)
    end do
    end do
    end do
  else
  !$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,2)
    end do
    end do
    end do
  end if

end subroutine

subroutine hpsi_test_diis(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,ppg)
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  use hpsi_sub, only: hpsi
  implicit none
  type(s_wavefunction)  :: stpsi
  type(s_wavefunction)  :: shtpsi
  type(s_wf_info)       :: info_ob
  type(s_rgrid),intent(in) :: mg
  type(s_scalar)        :: v(1)
  integer :: nspin_1
  type(s_stencil) :: stencil
  type(s_pp_grid) :: ppg

  call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin_1,stencil,ppg)

end subroutine hpsi_test_diis

end module rmmdiis_sub
