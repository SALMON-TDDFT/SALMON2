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
module calc_iobnum_sub
  implicit none

contains

subroutine calc_iobnum(tmst,tproc,trank,tiobnum,nproc_ob,iparaway_ob)
  use inputoutput, only : ispin
  implicit none
  integer :: tmst,tproc,trank,tiobnum
  integer :: ttmst
  integer :: nproc_ob,iparaway_ob
  
  if(ispin==0)then
    ttmst=tmst
  else if(ispin==1)then
    ttmst=tmst/2
  end if
  
  if(iparaway_ob==1)then
    tiobnum=(trank+1)*ttmst/nproc_ob-trank*ttmst/nproc_ob
  else if(iparaway_ob==2)then
    if(mod(ttmst,tproc)==0)then
      tiobnum=ttmst/tproc
    else
      if(trank<mod(ttmst,tproc))then
        tiobnum=ttmst/tproc+1
      else
        tiobnum=ttmst/tproc
      end if
    end if
  end if
  
  if(ispin==1)then
    tiobnum=tiobnum*2
  end if
  
end subroutine calc_iobnum

end module calc_iobnum_sub
