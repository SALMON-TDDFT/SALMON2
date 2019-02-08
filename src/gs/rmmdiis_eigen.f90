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
subroutine diis_core(mg,itotmst,mst,hvol,phi,R1,phibar,Rbar,iob,iter,pcheck)
  use inputoutput, only: ncg
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_korbital
  use inner_product_sub
  use eigen_sub
  !$ use omp_lib
  implicit none
  
  type(s_rgrid),intent(in) :: mg
  integer,intent(in) :: itotmst
  integer,intent(in) :: mst(2)
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
  character(30) :: commname
  
  commname='nproc_group_korbital'
  
  allocate(Rmat(iter,iter))
  allocate(Smat(iter,iter))
  allocate(alpha(iter),betav(iter),eval(iter))
  allocate(evec(iter,iter))
        
  do ii=0,iter-1
    do jj=0,iter-1
      call inner_product(mg,R1(:,:,:,ii),R1(:,:,:,jj),rbox,commname)
      Rmat(ii+1,jj+1)=rbox*hvol
  
      call inner_product(mg,phi(:,:,:,ii),phi(:,:,:,jj),rbox,commname)
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
      
    call inner_product(mg,phibar(:,:,:,iter-1),phibar(:,:,:,iter-1),rbox,commname)
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
      call inner_product(mg,phibar(:,:,:,iter-1),phibar(:,:,:,iter-1),rbox,commname)
   
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
  
      call inner_product(mg,phibar(:,:,:,iter-1),phibar(:,:,:,iter-1),rbox,commname)
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
