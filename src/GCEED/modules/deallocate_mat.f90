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
MODULE deallocate_mat_sub

use allocate_mat_sub
implicit none

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE deallocate_mat
  implicit none

deallocate (matbox_m,matbox_m2)
deallocate (matbox_l,matbox_l2)
deallocate (cmatbox_m,cmatbox_m2)
deallocate (cmatbox_l,cmatbox_l2)

deallocate (rho_tmp)

END SUBROUTINE deallocate_mat

!======================================================================

END MODULE deallocate_mat_sub
