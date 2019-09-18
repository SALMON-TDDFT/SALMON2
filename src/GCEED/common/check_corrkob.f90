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
module check_corrkob_sub
  implicit none

contains

subroutine check_corrkob(iob,info,ik,icorr_p,ilsda,nproc_ob,k_sta,k_end,mst)
  use calc_iquotient_sub
  use structures, only: s_orbital_parallel
  implicit none
  type(s_orbital_parallel),intent(in) :: info
  integer,intent(in)  :: iob,ik
  integer,intent(out) :: icorr_p
  integer,intent(in)  :: ilsda,nproc_ob,k_sta,k_end,mst(2)
  integer :: iquotient
  integer :: iob_tmp

  if(ilsda==0.or.iob<=mst(1))then
    iob_tmp=iob
  else
    iob_tmp=iob-mst(1)
  end if 

  call calc_iquotient(iob_tmp,nproc_ob,mst(1),iquotient)
  if(info%id_o==iquotient.and.ik>=k_sta.and.ik<=k_end)then
    icorr_p=1
  else
    icorr_p=0
  end if
  
end subroutine check_corrkob

end module check_corrkob_sub
