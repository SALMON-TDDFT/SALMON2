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
module calc_myob_sub
  implicit none

contains

subroutine calc_myob(iob,iob_myob,ilsda,nproc_ob,itotmst,mst,iobnum)
  use calc_iquotient_sub
  implicit none
  integer,intent(in)  :: iob
  integer,intent(out) :: iob_myob
  integer,intent(in)  :: ilsda,nproc_ob,itotmst,mst(2),iobnum
  integer :: iquotient,iob_min
  integer :: iob_tmp
  
  if(ilsda==0)then
    call calc_iquotient(iob,nproc_ob,itotmst,iquotient)
    iob_min=itotmst*iquotient/nproc_ob
    iob_myob=iob-iob_min
  else
    if(iob<=mst(1))then
      call calc_iquotient(iob,nproc_ob,mst(1),iquotient)
      iob_min=mst(1)*iquotient/nproc_ob
      iob_myob=iob-iob_min
    end if
  end if

end subroutine calc_myob

end module calc_myob_sub
