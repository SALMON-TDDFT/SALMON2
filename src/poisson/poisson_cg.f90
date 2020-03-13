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
!=======================================================================
module poisson_cg_sub
  implicit none

contains

!============================ Hartree potential (Solve Poisson equation)
subroutine poisson_cg(lg,mg,ng,info,system,poisson,trho,tVh,srg_ng,stencil)
  use inputoutput, only: threshold_cg
  use structures, only: s_rgrid,s_parallel_info,s_dft_system,s_poisson,s_sendrecv_grid,s_stencil
  use communication, only: comm_is_root, comm_summation
  use math_constants, only : pi
  use sendrecv_grid, only: update_overlap_real8
  use poisson_boundary_sub
  
  implicit none
  integer,parameter :: ndh=4
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_parallel_info),intent(in) :: info
  type(s_dft_system),intent(in) :: system
  type(s_poisson),intent(inout) :: poisson
  real(8) :: trho(mg%is(1):mg%ie(1),    &
                  mg%is(2):mg%ie(2),      &
                  mg%is(3):mg%ie(3))
  real(8) :: tVh(mg%is(1):mg%ie(1),    &
                 mg%is(2):mg%ie(2),      &
                 mg%is(3):mg%ie(3))
  type(s_sendrecv_grid),intent(inout) :: srg_ng
  type(s_stencil),intent(in) :: stencil
  
  integer,parameter :: maxiter=1000
  integer :: ix,iy,iz,iter
  real(8) :: sum1,sum2,ak,ck
  real(8) :: tottmp
  real(8) :: totbox
  real(8) :: rlap_wk(ng%is_array(1):ng%ie_array(1),    &
                     ng%is_array(2):ng%ie_array(2),      &
                     ng%is_array(3):ng%ie_array(3))
  real(8) :: zk(ng%is_array(1):ng%ie_array(1),   &
                ng%is_array(2):ng%ie_array(2),   &
                ng%is_array(3):ng%ie_array(3))
  real(8) :: pk(ng%is_array(1):ng%ie_array(1),   &
                ng%is_array(2):ng%ie_array(2),   &
                ng%is_array(3):ng%ie_array(3))
  
  call poisson_boundary(lg,mg,ng,info,system,poisson,trho,pk)
  
!------------------------- C-G minimization
  
!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=ng%is_array(3),ng%ie_array(3)
  do iy=ng%is_array(2),ng%ie_array(2)
  do ix=ng%is_array(1),ng%ie_array(1)
    zk(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    pk(ix,iy,iz)=tVh(ix,iy,iz)
    zk(ix,iy,iz)=-4.d0*pi*trho(ix,iy,iz)
  end do
  end do
  end do
  call update_overlap_real8(srg_ng, ng, pk)
  call laplacian_poisson(ng,pk,rlap_wk,stencil%coef_lap0,stencil%coef_lap)
  
!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=ng%is_array(3),ng%ie_array(3)
  do iy=ng%is_array(2),ng%ie_array(2)
  do ix=ng%is_array(1),ng%ie_array(1)
    pk(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    zk(ix,iy,iz)=zk(ix,iy,iz)-rlap_wk(ix,iy,iz)
    pk(ix,iy,iz)=zk(ix,iy,iz)
  end do
  end do
  end do
  
  sum1=0.d0
!$omp parallel do reduction(+ : sum1) private(iz,iy,ix) collapse(2)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    sum1=sum1+zk(ix,iy,iz)**2*system%hvol
  end do
  end do
  end do
  
  if(info%isize_r==1)then
  else
    call comm_summation(sum1,sum2,info%icomm_r)
    sum1=sum2
  end if
  
  iteration : do iter=1,maxiter
  
    call update_overlap_real8(srg_ng, ng, pk)
    call laplacian_poisson(ng,pk,rlap_wk,stencil%coef_lap0,stencil%coef_lap)
  
    totbox=0d0
!$omp parallel do reduction(+ : totbox) private(iz,iy,ix) collapse(2)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      totbox=totbox+(zk(ix,iy,iz)*rlap_wk(ix,iy,iz))
    end do
    end do
    end do
  
    if(info%isize_r==1)then
      tottmp=totbox
    else
      call comm_summation(totbox,tottmp,info%icomm_r)
    end if
  
    ak=sum1/tottmp/system%hvol
  
!$omp parallel do private(iz,iy,ix) firstprivate(ak) collapse(2)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
       tVh(ix,iy,iz)=tVh(ix,iy,iz)+ak*pk(ix,iy,iz)
       zk(ix,iy,iz)=zk(ix,iy,iz)-ak*rlap_wk(ix,iy,iz)
    end do
    end do
    end do
  
    totbox=0d0
!$omp parallel do reduction(+ : totbox) private(iz,iy,ix) collapse(2)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      totbox=totbox+zk(ix,iy,iz)**2
    end do
    end do
    end do
  
    if(info%isize_r==1)then
      tottmp=totbox
    else
      call comm_summation(totbox,tottmp,info%icomm_r)
    end if
  
    sum2=tottmp*system%hvol
  
    if ( abs(sum2) < threshold_cg*dble(lg%num(1)*lg%num(2)*lg%num(3)) ) exit
  
    ck=sum2/sum1 ; sum1=sum2
  
!$omp parallel do private(iz,iy,ix) firstprivate(ck) collapse(2)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      pk(ix,iy,iz)=zk(ix,iy,iz)+ck*pk(ix,iy,iz)
    end do
    end do
    end do
     
  end do iteration
  
  poisson%iterVh=iter
  if ( poisson%iterVh>maxiter .and. comm_is_root(info%id_r)) then
     write(*,*) "Warning:Vh iteration is not converged"
     write(*,'("||tVh(i)-tVh(i-1)||**2/(# of grids) = ",e15.8)') &
                                sum2/dble(lg%num(1)*lg%num(2)*lg%num(3))
  end if
  
  return

end subroutine poisson_cg

subroutine laplacian_poisson(ng,pk,rlap_wk,lap0,lapt)
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: ng
  real(8),intent(in) :: pk(ng%is_array(1):ng%ie_array(1),  &
                           ng%is_array(2):ng%ie_array(2),  &
                           ng%is_array(3):ng%ie_array(3))
  real(8),intent(out) :: rlap_wk(ng%is_array(1):ng%ie_array(1),  &
                                 ng%is_array(2):ng%ie_array(2),  &
                                 ng%is_array(3):ng%ie_array(3))
  real(8),intent(in)  :: lap0,lapt(4,3)
  integer :: ix,iy,iz

!$omp parallel do private(iz,iy,ix) collapse(2)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rlap_wk(ix,iy,iz)=-2.d0*lap0*pk(ix,iy,iz)+(  &
                      lapt(1,1)*(pk(ix+1,iy,iz) + pk(ix-1,iy,iz)) &
                     +lapt(2,1)*(pk(ix+2,iy,iz) + pk(ix-2,iy,iz)) &
                     +lapt(3,1)*(pk(ix+3,iy,iz) + pk(ix-3,iy,iz)) &
                     +lapt(4,1)*(pk(ix+4,iy,iz) + pk(ix-4,iy,iz)) &
                     +lapt(1,2)*(pk(ix,iy+1,iz) + pk(ix,iy-1,iz)) &
                     +lapt(2,2)*(pk(ix,iy+2,iz) + pk(ix,iy-2,iz)) &
                     +lapt(3,2)*(pk(ix,iy+3,iz) + pk(ix,iy-3,iz)) &
                     +lapt(4,2)*(pk(ix,iy+4,iz) + pk(ix,iy-4,iz)) &
                     +lapt(1,3)*(pk(ix,iy,iz+1) + pk(ix,iy,iz-1)) &
                     +lapt(2,3)*(pk(ix,iy,iz+2) + pk(ix,iy,iz-2)) &
                     +lapt(3,3)*(pk(ix,iy,iz+3) + pk(ix,iy,iz-3)) &
                     +lapt(4,3)*(pk(ix,iy,iz+4) + pk(ix,iy,iz-4)))
  end do
  end do
  end do

  return 

end subroutine laplacian_poisson

end module poisson_cg_sub
