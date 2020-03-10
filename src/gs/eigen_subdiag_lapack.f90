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
module eigen_lapack
  implicit none

  public :: eigen_dsyev, eigen_zheev

contains

  subroutine eigen_dsyev(h,e,v)
    implicit none
    real(8), intent(in) :: h(:,:)
    real(8), intent(out) :: e(:)
    real(8), intent(out) :: v(:,:)
    real(8), allocatable :: work(:)
    integer :: n, lwork, info

    n = ubound(h,1)
    lwork = 3*n-1
    allocate(work(lwork))
    v=h
    call dsyev('V', 'U', n, v, n, e, work, lwork, info)
    deallocate(work)
    return
  end subroutine eigen_dsyev

  subroutine eigen_zheev(h,e,v)
    implicit none
    complex(8), intent(in)  :: h(:,:)
    real(8), intent(out) :: e(:)
    complex(8), intent(out) :: v(:,:)
    complex(8), allocatable :: work(:)
    real(8), allocatable :: rwork(:)
    integer :: n,lwork,info

    n = ubound(h,1)
    lwork = 2*n-1
    allocate(work(lwork), rwork(3*n-2))
    v=h
    call ZHEEV('V', 'U', n, v, n, e, work, lwork, rwork, info)
    deallocate(work, rwork)
    return
  end subroutine eigen_zheev

end module eigen_lapack
