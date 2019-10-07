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
subroutine calc_pmax(iobmax,info)
use structures, only: s_orbital_parallel
use calc_iobnum_sub
use scf_data
implicit none
integer :: iobmax
type(s_orbital_parallel),intent(in) :: info

if(iSCFRT==1)then
  iobmax=iobnum
else if(iSCFRT==2)then
  if(ilsda==0)then
    call calc_iobnum(ifMST(1),info,iobmax,nproc_ob)
  else if(ilsda==1)then
    iobmax=iobnum
  end if
end if

end subroutine calc_pmax
