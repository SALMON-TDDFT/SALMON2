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
module allocate_psl_sub

use scf_data
use read_pslfile_sub

real(8),allocatable :: dVloc_G(:,:)
real(8),allocatable :: dVloc_G_tmp(:,:)
complex(8),allocatable :: rhoion_G(:)
complex(8),allocatable :: rhoion_G_tmp(:)

contains
!==================================================================================================

subroutine allocate_psl(lg)
use structures, only: s_rgrid
use scf_data
implicit none
type(s_rgrid),intent(in) :: lg

if(iperiodic==3)then
  allocate(dVloc_G(lg%num(1)*lg%num(2)*lg%num(3),MKI))
  allocate(dVloc_G_tmp(lg%num(1)*lg%num(2)*lg%num(3),MKI))
  allocate(rhoion_G(lg%num(1)*lg%num(2)*lg%num(3)))
  allocate(rhoion_G_tmp(lg%num(1)*lg%num(2)*lg%num(3)))
end if

end subroutine allocate_psl
!==================================================================================================
subroutine deallocate_psl

if(iperiodic==3)then
  deallocate(dVloc_G,rhoion_G)
  deallocate(dVloc_G_tmp,rhoion_G_tmp)
end if

end subroutine deallocate_psl
!==================================================================================================
end module allocate_psl_sub
