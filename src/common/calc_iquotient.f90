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
module calc_iquotient_sub
  implicit none

contains

subroutine calc_iquotient(iob,nproc_ob,itotmst,iquotient)
  implicit none
  integer,intent(in) :: iob,nproc_ob,itotmst
  integer,intent(out) :: iquotient
  
  if(mod(iob*nproc_ob,itotmst)==0)then
    iquotient=iob*nproc_ob/itotmst-1
  else
    iquotient=iob*nproc_ob/itotmst
  end if

end subroutine calc_iquotient

end module calc_iquotient_sub
