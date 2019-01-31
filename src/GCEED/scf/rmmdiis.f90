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
!=======================================================================
!============================================================== RMM-DIIS
! This routine is RMM-DIIS
! J. Soc. Mat. Sci., Japan, vol.52 (3), p.260-265. (in Japanese)

SUBROUTINE rmmdiis(mg,info,spsi)
use structures, only: s_rgrid,s_wf_info,s_wavefunction
use salmon_parallel, only: nproc_group_global
use salmon_communication, only: comm_summation
use scf_data
use hpsi2_sub
use new_world_sub
!$ use omp_lib
implicit none

type(s_rgrid),intent(in) :: mg
type(s_wf_info) :: info
type(s_wavefunction) :: spsi
integer :: iob,iter,ix,iy,iz
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
real(8) :: rnorm_diff_psi(itotMST,1)
real(8) :: tpsi(mg%is_array(1):mg%ie_array(1),  &
                mg%is_array(2):mg%ie_array(2),  &
                mg%is_array(3):mg%ie_array(3))

allocate (htphi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
allocate (phibox(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
allocate (Rbox(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
allocate (phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0:Ncg))
allocate (R1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0:Ncg))
allocate (phibar(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0:Ncg))
allocate (Rbar(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),0:Ncg))
allocate (psi_stock(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,  &
                    1:info%numo,info%ik_s:info%ik_e,1))

allocate (iobcheck(1:itotMST,0:Ncg))
iobcheck=0

iwk_size=2
call make_iwksta_iwkend


!$OMP parallel do private(iz,iy,ix)
do iz=mg%is_array(3),mg%ie_array(3)
do iy=mg%is_array(2),mg%ie_array(2)
do ix=mg%is_array(1),mg%ie_array(1)
  tpsi(ix,iy,iz)=0.d0
end do
end do
end do

if(info%numo>=1)then
  allocate (iflagdiis(1:info%numo))
  allocate (epsdiis(1:info%numo,0:Ncg))
  allocate (Rnorm(1:info%numo,0:Ncg))
end if 

! Flag for convergence
if(info%numo >= 1) iflagdiis=1

if(info%numo >= 1) then
  phi=0.d0
  psi_stock=spsi%rwf
end if

iflag_diisjump=0

do iob=1,info%numo

  Iteration : do iter=1,Ncg

  if(iter == 1) then
! Obtain residual vector R_0
!$OMP parallel do
    do iz=mg%is(3),mg%ie(3)
      phi(:,:,iz,0)=spsi%rwf(:,:,iz,1,iob,1,1)
    end do
!$OMP parallel do private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    tpsi(ix,iy,iz)=phi(ix,iy,iz,0)
  end do
  end do
  end do

    call hpsi2(tpsi,htphi(:,:,:),iob,1,0,0)
    call inner_product3(mg,phi(mg_sta(1),mg_sta(2),mg_sta(3),0),htphi(mg_sta(1),mg_sta(2),mg_sta(3)),rbox1,elp3)

!$OMP parallel do
    do iz=mg%is(3),mg%ie(3)
      R1(:,:,iz,0)=htphi(:,:,iz)-rbox1*Hvol*phi(:,:,iz,0)
    end do 
    epsdiis(iob,0)=rbox1*Hvol
    call inner_product3(mg,R1(mg_sta(1),mg_sta(2),mg_sta(3),0),R1(mg_sta(1),mg_sta(2),mg_sta(3),0),rbox1,elp3)
    Rnorm(iob,0)=rbox1*Hvol

  else
! Solve by Lagrange's method of undetermined multipliers, and obtain 
! Rbar from previous combinations of phi and R.
    if(iflagdiis(iob) == 1)then
      call diis_core(mg,phi,R1,phibar,Rbar,iob,iter,iobcheck)
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
    call inner_product3(mg,phi(mg_sta(1),mg_sta(2),mg_sta(3),iter),phi(mg_sta(1),mg_sta(2),mg_sta(3),iter),rbox1,elp3)
!$OMP parallel do
    do iz=mg%is(3),mg%ie(3)
      phi(:,:,iz,iter)=phi(:,:,iz,iter)/sqrt(rbox1*Hvol)
    end do

!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      tpsi(ix,iy,iz)=phi(ix,iy,iz,iter)
    end do
    end do
    end do

    call hpsi2(tpsi,htphi(:,:,:),iob,1,0,0)
    call inner_product3(mg,phi(mg_sta(1),mg_sta(2),mg_sta(3),iter),htphi(mg_sta(1),mg_sta(2),mg_sta(3)),rbox1,elp3)
!$OMP parallel do
    do iz=mg%is(3),mg%ie(3)
      R1(:,:,iz,iter)=htphi(:,:,iz)-rbox1*Hvol*phi(:,:,iz,iter)
    end do

    call inner_product3(mg,phi(mg_sta(1),mg_sta(2),mg_sta(3),iter),htphi(mg_sta(1),mg_sta(2),mg_sta(3)),rbox1,elp3)
    epsdiis(iob,iter)=rbox1*Hvol

    call inner_product3(mg,R1(mg_sta(1),mg_sta(2),mg_sta(3),iter),R1(mg_sta(1),mg_sta(2),mg_sta(3),iter),rbox1,elp3)
    Rnorm(iob,iter)=rbox1*Hvol

! judgement for closing loop.
! The ratio of Rnorm is set to 0.3 as well as Kresse-Furthmuller.
    if(iter >= 2) then
      if(iter >= 3 .and. epsdiis(iob,iter) > epsdiis(iob,iter-1)) then
!$OMP parallel do
        do iz=mg%is(3),mg%ie(3)
          spsi%rwf(:,:,iz,1,iob,1,1) = phi(:,:,iz,iter-1)
        end do
        iflagdiis(iob)=0
      else if(-(epsdiis(iob,iter)-epsdiis(iob,iter-1)) <= 1.0d-8 .or.      &
              Rnorm(iob,iter)/Rnorm(iob,0) <= 0.3d0) then
!$OMP parallel do
        do iz=mg%is(3),mg%ie(3)
          spsi%rwf(:,:,iz,1,iob,1,1) = phi(:,:,iz,iter)
        end do
        iflagdiis(iob)=0
      end if
    end if

    if(iter == Ncg) then
!$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
        spsi%rwf(:,:,iz,1,iob,1,1) = phi(:,:,iz,Ncg)
      end do
     end if
    if(iter == 1 .and. iflag_diisjump == 1) then
!$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
        spsi%rwf(:,:,iz,1,iob,1,1) = phi(:,:,iz,1)
      end do
    end if

!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      tpsi(ix,iy,iz)=phi(ix,iy,iz,iter)
    end do
    end do
    end do

    call hpsi2(tpsi,htphi(:,:,:),iob,1,0,0)
    call inner_product3(mg,phi(mg_sta(1),mg_sta(2),mg_sta(3),iter),htphi(mg_sta(1),mg_sta(2),mg_sta(3)),rbox1,elp3)
    
    end if

  end do Iteration

end do        ! loop for iob

iflag_diisjump=0
do iob=1,info%numo
!$OMP parallel do private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    tpsi(ix,iy,iz)=spsi%rwf(ix,iy,iz,1,iob,1,1)
  end do
  end do
  end do

  call hpsi2(tpsi,htphi(:,:,:),iob,1,0,0)
  call inner_product3(mg,spsi%rwf(mg_sta(1),mg_sta(2),mg_sta(3),1,iob,1,1),htphi(mg_sta(1),mg_sta(2),mg_sta(3)),rbox1,elp3)
  rbox1=sum(spsi%rwf(:,:,:,1,iob,1,1)*htphi(:,:,:))*Hvol
  if(rbox1-esp(iob,1)>5.d0) iflag_diisjump=1
end do

if(iflag_diisjump==0)then
  continue
else if(iflag_diisjump==1)then
  spsi%rwf=psi_stock
  do iob=1,info%numo
  
!$OMP parallel do
    do iz=mg%is(3),mg%ie(3)
      phi(:,:,iz,0)=spsi%rwf(:,:,iz,1,iob,1,1)
    end do
!$OMP parallel do private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      tpsi(ix,iy,iz)=phi(ix,iy,iz,0)
    end do
    end do
    end do

    call hpsi2(tpsi,htphi(:,:,:),iob,1,0,0)

    call inner_product3(mg,phi(mg_sta(1),mg_sta(2),mg_sta(3),0),htphi(mg_sta(1),mg_sta(2),mg_sta(3)),rbox1,elp3)

!$OMP parallel do
    do iz=mg%is(3),mg%ie(3)
      R1(:,:,iz,0)=htphi(:,:,iz)-rbox1*Hvol*phi(:,:,iz,0)
      spsi%rwf(:,:,iz,1,iob,1,1)=phi(:,:,iz,0)+lambda2_diis*R1(:,:,iz,0)
    end do 

  end do

  rnorm_diff_psi=0.d0
  do iob=1,info%numo
    phi(:,:,:,0)=abs(spsi%rwf(:,:,:,1,iob,1,1)-psi_stock(:,:,:,1,iob,1,1))
    rbox1=sum(phi(:,:,:,0)*phi(:,:,:,0))*Hvol
    rnorm_diff_psi(iob,1)=rbox1
  end do
  call comm_summation(rnorm_diff_psi,norm_diff_psi_stock,itotMST,nproc_group_global)
end if


deallocate(htphi)
deallocate(phibox,Rbox,phi,R1,phibar,Rbar)

if(info%numo>=1)then
  deallocate (iflagdiis,epsdiis,Rnorm)
end if 
deallocate(iobcheck) 

return

END SUBROUTINE rmmdiis

