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
module calc_allob_sub
  implicit none

contains

subroutine calc_allob(iob,info,iob_allob,itotmst,mst,iobnum)
  use inputoutput, only: ispin,nproc_ob
  use structures, only: s_orbital_parallel
  implicit none
  integer,intent(in) :: iob
  type(s_orbital_parallel),intent(in) :: info
  integer,intent(out) ::  iob_allob
  integer,intent(in) :: itotmst,mst(2),iobnum
  integer :: iob_tmp
  
  if(ispin==0)then
    iob_allob=info%id_o*itotmst/nproc_ob+iob
  else
    if(iob<=iobnum/2)then
      iob_allob=info%id_o*mst(1)/nproc_ob+iob
    else
      iob_tmp=iob-iobnum/2
      iob_allob=info%id_o*mst(2)/nproc_ob+iob_tmp+mst(1)
    end if
  end if
  
end subroutine calc_allob

end module calc_allob_sub
