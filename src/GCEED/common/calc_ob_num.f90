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
module calc_iobnum_sub
  implicit none

contains

subroutine calc_iobnum(tmst,info,tiobnum,nproc_ob)
  use inputoutput, only : ispin
  use structures, only : s_orbital_parallel
  implicit none
  integer :: tmst,tiobnum
  type(s_orbital_parallel),intent(in) :: info
  integer :: ttmst
  integer :: nproc_ob
  
  if(ispin==0)then
    ttmst=tmst
  else if(ispin==1)then
    ttmst=tmst/2
  end if
  
  tiobnum=(info%id_o+1)*ttmst/nproc_ob-info%id_o*ttmst/nproc_ob
  
  if(ispin==1)then
    tiobnum=tiobnum*2
  end if
  
end subroutine calc_iobnum

end module calc_iobnum_sub
