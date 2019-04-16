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
module gram_schmidt_orth
  implicit none

contains

  subroutine gram_schmidt()
    implicit none
    ! We need a few types of implementetaion for the following cases:
    ! 1. GS for isolated system
    ! 2. GS for large-scale isolated system (domain parallel)
    ! 3. GS for regular periodic system 
    !  3-1. Small number of k-points (based on ARTED Gram-Schmidt type-1)
    !  3-2. Large number of k-points (based on ARTED Gram-Schmidt type-2)
    ! 4. GS for large-scale periodic system (domain parallel)

    return
  end subroutine

  ! Gram_Schmitd for orbital(b) omp
  subroutine gram_schmidt_omp_b()
    implicit none
    return
  end subroutine gram_schmidt_omp_b

  ! Gram_Schmitd for k-space omp
  subroutine gram_schmidt_omp_k()
    implicit none
    return
  end subroutine gram_schmidt_omp_k

  ! Gram_Schmitd for realspace(r) omp
  subroutine gram_schmidt_omp_r()
    implicit none
    return
  end subroutine gram_schmidt_omp_r

end module gram_schmidt_orth
