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
subroutine calc_myob(iob,iob_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst,iobnum)
  use salmon_parallel, only: nproc_id_kgrid
  implicit none
  integer,intent(in)  :: iob
  integer,intent(out) :: iob_myob
  integer,intent(in)  :: ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin(2),mst(2),iobnum
  integer :: iquotient,iob_min
  integer :: iob_tmp
  
  if(ilsda==0)then
    if(iparaway_ob==1)then
      call calc_iquotient(iob,nproc_ob,itotmst,iquotient)
      iob_min=itotmst*iquotient/nproc_ob
      iob_myob=iob-iob_min
    else if(iparaway_ob==2)then
      iob_myob=(iob-1)/nproc_ob+1
    end if
  else
    if(iparaway_ob==1)then
      if(iob<=mst(1))then
        call calc_iquotient(iob,nproc_ob,mst(1),iquotient)
        iob_min=mst(1)*iquotient/nproc_ob
        iob_myob=iob-iob_min
      else
        iob_tmp=iob-mst(1)
        call calc_iquotient(iob_tmp,nproc_ob,mst(1),iquotient)
        iob_min=mst(1)*iquotient/nproc_ob
        iob_myob=iob_tmp-iob_min+iobnum/2
      end if
    else if(iparaway_ob==2)then
      if(iob<=mst(1))then
        iob_myob=(iob-1)/nproc_ob+1
      else
        iob_tmp=iob-mst(1)
        iob_myob=(iob_tmp-1)/nproc_ob+1+iobnum/2
      end if
    end if
  end if

end subroutine calc_myob
