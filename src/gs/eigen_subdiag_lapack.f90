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
module eigen_subdiag_sub
  implicit none

contains

subroutine eigen_subdiag(Rmat,evec,iter,ier2)
  implicit none
  
  integer :: iter,ier2
  real(8) :: Rmat(iter,iter)
  real(8) :: evec(iter,iter)
  
  character(1) :: JOBZ,UPLO
  integer :: N
  real(8) :: A(iter,iter)
  integer :: LDA
  real(8) :: W(iter)
  real(8) :: WORK(3*iter-1)
  integer :: LWORK
  
  ier2=0
  
  JOBZ='V'
  UPLO='L'
  N=iter
  A=Rmat
  LDA=iter
  LWORK=3*iter-1
  call DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,ier2)
  evec=A

end subroutine eigen_subdiag

subroutine eigen_subdiag_periodic(Rmat,evec,iter,ier2)
  implicit none
  character :: JOBZ, UPLO
  integer :: LWORK
  integer :: iter,ier2
  real(8),allocatable :: RWORK(:)
  real(8) :: W(iter)
  complex(8) :: Rmat(iter,iter)
  complex(8),allocatable :: WORK(:)
  complex(8) :: evec(iter,iter)

  ier2=0

  JOBZ='V'
  UPLO='U'

  LWORK=2*iter-1
  allocate(WORK(LWORK))
  allocate(RWORK(3*iter-2))

  evec(:,:)=Rmat(:,:)

  call ZHEEV(JOBZ,UPLO,iter,evec,iter,W,WORK,LWORK,RWORK,ier2)

  deallocate(WORK,RWORK)

end subroutine eigen_subdiag_periodic
 
end module eigen_subdiag_sub
